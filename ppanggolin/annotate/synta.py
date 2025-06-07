#!/usr/bin/env python3

# default libraries
import logging
import os
import tempfile
from io import TextIOWrapper
from multiprocessing import Value
from subprocess import Popen, PIPE
import ast
from collections import defaultdict
from typing import Dict, List, Optional, Union, Generator, Tuple
from pathlib import Path
import shutil

# install libraries
from pyrodigal import GeneFinder, Sequence

# local libraries
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.utils import (
    is_compressed,
    read_compressed_or_not,
    check_tools_availability,
)

contig_counter: Value = Value("i", 0)


def init_contig_counter(value: Value):
    """Initialize the contig counter for later use"""
    global contig_counter
    contig_counter = value


def reverse_complement(seq: str):
    """reverse complement the given dna sequence

    :param seq: sequence which need to be reversed

    :return: reverse sequence
    """

    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
    }
    # see https://www.bioinformatics.org/sms/iupac.html for the code.
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq


def launch_aragorn(
    fna_file: str, org: Organism, contig_to_length: Dict[str, int]
) -> defaultdict:
    """
    Launches Aragorn to annotate tRNAs.

    :param fna_file: file-like object containing the uncompressed fasta sequences
    :param org: Organism which will be annotated

    :return: Annotated genes in a list of gene objects
    """
    locustag = org.name
    cmd = ["aragorn", "-t", "-gcbact", "-l", "-w", fna_file]
    logging.getLogger("PPanGGOLiN").debug(f"aragorn command : {' '.join(cmd)}")

    if shutil.which("aragorn") is None:
        raise FileNotFoundError(
            f"Command 'aragorn' not found. Please install it and try again."
        )

    p = Popen(cmd, stdout=PIPE)
    # loading the whole thing, reverting it to 'pop' in order.
    file_data = p.communicate()[0].decode().split("\n")[::-1]
    gene_objs = defaultdict(set)
    c = 0
    contig_name = ""
    while len(file_data) != 0:
        line = file_data.pop()
        if line.startswith(">"):
            contig_name = line.replace(">", "").split()[0]
            file_data.pop()  # then next line must be removed too.
        elif len(line) > 0:  # if the line isn't empty, there's data to get.
            line_data = line.split()
            start, stop = map(int, ast.literal_eval(line_data[2].replace("c", "")))
            if start < 1 or stop < 1:
                # In some case aragorn gives negative coordinates. This case is just ignored.
                logging.warning(
                    f"Aragorn gives non valid coordinates for a RNA gene in contig {contig_name}: {line_data}. This RNA is ignored."
                )
                continue
            if (
                start > contig_to_length[contig_name]
                or stop > contig_to_length[contig_name]
            ):
                logging.warning(
                    f"Aragorn gives non valide coordinates for a RNA gene in contig {contig_name}. "
                    f"Gene coordinates exceed contig length ({contig_to_length[contig_name]}): "
                    f"{line_data}. This RNA is ignored."
                )
                continue

            c += 1
            gene = RNA(rna_id=locustag + "_tRNA_" + str(c).zfill(4))
            gene.fill_annotations(
                start=start,
                stop=stop,
                strand="-" if line_data[2].startswith("c") else "+",
                gene_type="tRNA",
                product=line_data[1] + line_data[4],
            )
            gene_objs[contig_name].add(gene)
    return gene_objs


def launch_prodigal(
    contig_sequences: Dict[str, str],
    org: Organism,
    code: int = 11,
    use_meta: bool = False,
) -> defaultdict:
    """
    Launches Prodigal to annotate CDS. Takes a fna file name and a locustag to give an ID to the pred genes.

    :param contig_sequences: Dict containing contig sequences for pyrodigal
    :param org: Organism which will be annotated
    :param code: Translation table (genetic code) to use.
    :param use_meta: use meta procedure in Prodigal

    :return: Annotated genes in a list of gene objects
    """
    gene_objs = defaultdict(set)
    sequences = {
        contig_name: Sequence(sequence)
        for contig_name, sequence in contig_sequences.items()
    }
    gene_finder = GeneFinder(
        meta=use_meta,  # '-p meta' if meta is true else '-p single'
        closed=True,  # -c: Closed ends. Do not allow genes to run off edges.
        mask=True,  # -m: Treat runs of N as masked sequence; don't build genes across them.
        min_gene=120,  # This is to prevent error with mmseqs translatenucs that cut too short sequences
    )

    if not use_meta:
        gene_finder.train(
            *contig_sequences.values(), force_nonsd=False, translation_table=code
        )  # -g: Specify a translation table to use (default 11).
    gene_counter = 1
    for contig_name, sequence in sequences.items():
        for pred in gene_finder.find_genes(sequence):
            gene = Gene(gene_id=f"{org.name}_CDS_{str(gene_counter).zfill(4)}")
            gene.fill_annotations(
                start=pred.begin,
                stop=pred.end,
                strand="-" if pred.strand == -1 else "+",
                gene_type="CDS",
                genetic_code=code,
            )
            gene_counter += 1
            gene_objs[contig_name].add(gene)
    return gene_objs


def launch_infernal(
    fna_file: str, org: Organism, tmpdir: str, kingdom: str = "bacteria"
) -> defaultdict:
    """
    Launches Infernal in hmmer-only mode to annotate rRNAs.

    :param fna_file: file-like object containing the uncompressed fasta sequences
    :param org: Organism which will be annotated
    :param kingdom: Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    :param tmpdir: Path to temporary directory

    :return: Annotated genes in a list of gene objects.
    """
    locustag = org.name
    modelfile = ""
    if kingdom == "bacteria":
        modelfile = (
            os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_bact.cm"
        )
    elif kingdom == "archaea":
        modelfile = (
            os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_arch.cm"
        )

    tmp_file = tempfile.NamedTemporaryFile(mode="r", dir=tmpdir)
    cmd = [
        "cmscan",
        "--tblout",
        tmp_file.name,
        "--hmmonly",
        "--cpu",
        str(1),
        "--noali",
        modelfile,
        fna_file,
    ]
    logging.getLogger("PPanGGOLiN").debug(f"infernal command : {' '.join(cmd)}")
    if shutil.which("cmscan") is None:
        raise FileNotFoundError(
            f"Command 'cmscan' from the tool Infernal is not found. Please install it and try again."
        )

    p = Popen(cmd, stdout=open(os.devnull, "w"), stderr=PIPE)
    err = p.communicate()[1].decode().split()
    if err:
        if err[0] == "Error: ":
            raise Exception(
                f"Infernal (cmscan) failed with error:  '{' '.join(err)}'. If you never used this script,"
                f" you should press the .cm file using cmpress executable from Infernal. "
                f"You should find the file in '{os.path.dirname(os.path.realpath(__file__))}/rRNA_DB/'."
            )
        raise Exception(
            f"An error occurred with Infernal. Error is:  '{' '.join(err)}'."
        )
    # never managed to test what happens if the .cm files are compressed with a 'bad' version of infernal,
    # so if that happens you are on your own.

    gene_objs = defaultdict(set)
    c = 0
    for line in tmp_file:
        if not line.startswith("#"):
            c += 1
            line_data = line.split()
            strand = line_data[9]
            start, stop = map(
                int,
                (
                    (line_data[8], line_data[7])
                    if strand == "-"
                    else (line_data[7], line_data[8])
                ),
            )
            gene = RNA(rna_id=locustag + "_rRNA_" + str(c).zfill(4))
            gene.fill_annotations(
                start=start,
                stop=stop,
                strand=strand,
                gene_type="rRNA",
                product=" ".join(line_data[17:]),
            )
            gene_objs[line_data[2]].add(gene)
    return gene_objs


def check_sequence_tuple(name: str, sequence: str):
    """
    Checks and validates a sequence name and its corresponding sequence.

    :param name: The name (header) of the sequence, typically extracted from the FASTA file header.
    :param sequence: The sequence string corresponding to the name, containing the nucleotide or protein sequence.

    :return: A tuple containing the validated name and sequence.

    :raises ValueError:
        - If the sequence is empty, a ValueError is raised with a message containing the header name.
        - If the name is empty, a ValueError is raised with a message containing a preview of the sequence.
    """
    if not sequence:
        raise ValueError(f"Found an empty sequence with header '{name}'")

    if not name:
        raise ValueError(
            f"Found a sequence with empty name (sequence starts as '{sequence[:60]}')"
        )

    return name, sequence


def parse_fasta(
    fna_file: Union[TextIOWrapper, list]
) -> Generator[Tuple[str, str], None, None]:
    """Yields each sequence name and sequence from a FASTA file or stream as a tuple.

    :param fna_file: Input FASTA file or list of lines as sequences.
    :yield: Tuple with contig header (without '>') and sequence.
    :raises ValueError: If the file does not contain valid FASTA format.
    """
    name = None
    sequence = ""

    for line in fna_file:
        line = line.strip()

        if line.startswith(">"):  # New header
            if name:  # Yield previous header and sequence if available
                yield check_sequence_tuple(name, sequence)

            name = line[1:].split()[
                0
            ]  # Strip '>' and extract the first word as the name
            sequence = ""

        elif line:  # Only append non-empty lines
            sequence += line

        else:
            # You can skip or handle empty lines here if required
            pass

    # Yield the final contig if exists
    if name:
        yield check_sequence_tuple(name, sequence)

    # Check if there was any valid data (at least one header and sequence)
    if not name:
        raise ValueError("The file does not contain any valid FASTA content.")


def get_contigs_from_fasta_file(
    org: Organism, fna_file: Union[TextIOWrapper, list]
) -> Dict[str, str]:
    """Processes contigs from a parsed FASTA generator and stores in a dictionary.

    :param org: Organism instance to update with contig info.
    :param fna_file: Input FASTA file or list of lines as sequences.
    :return: Dictionary with contig names as keys and sequences as values.
    """

    global contig_counter
    contigs = {}

    for contig_name, sequence in parse_fasta(fna_file):

        # Retrieve or create the contig
        try:
            contig = org.get(contig_name)
        except KeyError:
            with contig_counter.get_lock():
                contig = Contig(contig_counter.value, contig_name)
                contig_counter.value += 1
            org.add(contig)

        # Update contig information
        if contig.length is not None and contig.length != len(sequence):
            raise ValueError(
                f"Length mismatch for contig {contig_name}: expected {contig.length}, found {len(sequence)} from the fasta sequence."
            )

        contig.length = len(sequence)
        contigs[contig_name] = sequence.upper()

    return contigs


def write_tmp_fasta(contigs: dict, tmpdir: str) -> tempfile._TemporaryFileWrapper:
    """
     Writes a temporary fna formatted file and returns the file-like object. Useful in case of  compressed input file.
     The file will be deleted when close() is called.

    :param contigs: Contigs sequences of each contig
    :param tmpdir: path to temporary directory

    :return: fasta file
    """

    tmp_file = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir)
    for header in contigs.keys():
        tmp_file.write(f">{header}\n")
        j = 0
        while j < len(contigs[header]):
            tmp_file.write(contigs[header][j : j + 60] + "\n")
            j += 60
    tmp_file.flush()  # force write what remains in the buffer.
    return tmp_file


def syntaxic_annotation(
    org: Organism,
    fasta_file: TextIOWrapper,
    contig_sequences: Dict[str, str],
    tmpdir: str,
    norna: bool = False,
    kingdom: str = "bacteria",
    code: int = 11,
    use_meta: bool = False,
) -> defaultdict:
    """
    Runs the different software for the syntaxic annotation.

    :param org: Organism which will be annotated
    :param fasta_file: file-like object containing the uncompressed fasta sequences
    :param contig_sequences: Dict containing contig sequences for pyrodigal
    :param tmpdir: Path to temporary directory
    :param norna: Use to avoid annotating RNA features.
    :param kingdom: Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    :param code: Translation table (genetic code) to use.
    :param use_meta: Use meta prodigal procedure

    :return: list of genes in the organism
    """

    # launching tools for syntaxic annotation
    genes = defaultdict(list)
    for contig_name, genes_from_contig in launch_prodigal(
        contig_sequences=contig_sequences, org=org, code=code, use_meta=use_meta
    ).items():
        genes[contig_name].extend(genes_from_contig)
    if not norna:
        contig_to_length = {
            contig_name: len(contig_seq)
            for contig_name, contig_seq in contig_sequences.items()
        }

        for contig_name, genes_from_contig in launch_aragorn(
            fna_file=fasta_file.name, org=org, contig_to_length=contig_to_length
        ).items():
            genes[contig_name].extend(genes_from_contig)
        for contig_name, genes_from_contig in launch_infernal(
            fna_file=fasta_file.name, org=org, kingdom=kingdom, tmpdir=tmpdir
        ).items():
            genes[contig_name].extend(genes_from_contig)
    fasta_file.close()  # closing either tmp file or original fasta file.
    return genes


def overlap_filter(all_genes: defaultdict, allow_overlap: bool = False) -> defaultdict:
    """
    Removes the CDS that overlap with RNA genes.

    :param all_genes: Dictionary with complete list of genes
    :param allow_overlap: Use to not remove genes overlapping with RNA features

    :return: Dictionary with genes filtered
    """

    sorted_genes = defaultdict(list)
    for key, genes in all_genes.items():
        tmp_genes = sorted(genes, key=lambda x: x.start)
        rm_genes = set()
        if not allow_overlap:
            for i, gene_i in enumerate(tmp_genes):
                if i + 1 < len(tmp_genes):
                    gene_j = tmp_genes[i + 1]
                    if (
                        gene_i.type != "CDS"
                        and gene_j.type == "CDS"
                        and gene_i.stop > gene_j.start
                    ):
                        rm_genes.add(gene_j)
                    elif (
                        gene_i.type == "CDS"
                        and gene_j.type != "CDS"
                        and gene_i.stop > gene_j.start
                    ):
                        rm_genes.add(gene_i)

        for gene in rm_genes:
            tmp_genes.remove(gene)
        cds_counter = 0
        for gene in tmp_genes:
            if gene.type == "CDS":
                gene.position = cds_counter
                cds_counter += 1
        sorted_genes[key] = tmp_genes
    return sorted_genes


def get_dna_sequence(contig_seq: str, gene: Union[Gene, RNA]) -> str:
    """Return the gene sequence

    :param contig_seq: Contig sequence
    :param gene: Gene

    :return: str
    """

    # check contig coordinate is in scope of contig seq length
    highest_position = max((stop for _, stop in gene.coordinates))
    assert highest_position <= len(
        contig_seq
    ), f"Coordinates of gene {gene} exceed length of the contig. Gene coordinates {gene.coordinates} vs contig length {len(contig_seq)}"

    # Extract gene seq
    seq = "".join([contig_seq[start - 1 : stop] for start, stop in gene.coordinates])

    # check length of extracted seq
    assert len(seq) == len(gene), (
        f"The gene sequence of {gene} extracted from the contig does not have the expected length: "
        f"extracted seq length {len(seq)}nt vs expected length based on gene coordinates ({gene.coordinates}) {len(gene)}nt "
    )

    if gene.strand == "+":
        return seq
    elif gene.strand == "-":
        return reverse_complement(seq)


def annotate_organism(
    org_name: str,
    file_name: Path,
    circular_contigs: List[str],
    tmpdir: str,
    code: int = 11,
    norna: bool = False,
    kingdom: str = "bacteria",
    allow_overlap: bool = False,
    procedure: Optional[str] = None,
) -> Organism:
    """
    Function to annotate a single organism

    :param org_name: Name of the organism / genome
    :param file_name: Path to the fasta file containing organism sequences
    :param circular_contigs: list of contigs
    :param code: Translation table (genetic code) to use.
    :param kingdom: Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    :param norna: Use to avoid annotating RNA features.
    :param tmpdir: Path to temporary directory
    :param allow_overlap: Use to not remove genes overlapping with RNA features
    :param procedure: prodigal procedure used

    :return: Complete organism object for pangenome
    """

    if not norna:
        check_tools_availability({"aragorn": "", "cmscan": " from the tool Infernal"})

    org = Organism(org_name)

    fasta_file = read_compressed_or_not(file_name)

    contig_sequences = get_contigs_from_fasta_file(org, fasta_file)
    if is_compressed(file_name):  # TODO simply copy file with shutil.copyfileobj
        fasta_file = write_tmp_fasta(contig_sequences, tmpdir)
    if procedure is None:  # prodigal procedure is not force by user
        max_contig_len = max(len(contig) for contig in org.contigs)
        if max_contig_len < 20000:  # case of short sequence
            use_meta = True
            logging.getLogger("PPanGGOLiN").info(
                f"Using the metagenomic mode to predict genes for {org_name}, as all its contigs are < 20KB in size."
            )

        else:
            use_meta = False
    else:
        use_meta = True if procedure == "meta" else False
    genes = syntaxic_annotation(
        org, fasta_file, contig_sequences, tmpdir, norna, kingdom, code, use_meta
    )
    genes = overlap_filter(genes, allow_overlap=allow_overlap)

    for contig_name, genes in genes.items():
        contig = org.get(contig_name)
        contig.is_circular = True if contig.name in circular_contigs else False
        for gene in genes:
            gene.add_sequence(get_dna_sequence(contig_sequences[contig.name], gene))
            gene.fill_parents(org, contig)
            if isinstance(gene, Gene):
                contig.add(gene)
            elif isinstance(gene, RNA):
                contig.add_rna(gene)
    return org
