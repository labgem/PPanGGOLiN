#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import os
import tempfile
from io import TextIOWrapper
from multiprocessing import Value
from subprocess import Popen, PIPE
import ast
from collections import defaultdict
from typing import Dict, List, Union
from pathlib import Path

# install libraries
from pyrodigal import GeneFinder, Sequence

# local libraries
from ppanggolin.genome import Organism, Gene, RNA, Contig
from ppanggolin.utils import is_compressed, read_compressed_or_not


contig_counter: Value = Value('i', 0)


def init_contig_counter(value: Value):
    """Initialize the contig counter for later use"""
    global contig_counter
    contig_counter = value


def reverse_complement(seq: str):
    """reverse complement the given dna sequence

    :param seq: sequence which need to be reversed

    :return: reverse sequence
    """

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'R': 'Y', 'Y': 'R',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    # see https://www.bioinformatics.org/sms/iupac.html for the code.
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq


def launch_aragorn(fna_file: str, org: Organism) -> defaultdict:
    """
    Launches Aragorn to annotate tRNAs.

    :param fna_file: file-like object containing the uncompressed fasta sequences
    :param org: Organism which will be annotated

    :return: Annotated genes in a list of gene objects
    """
    locustag = org.name
    cmd = ["aragorn", "-t", "-gcbact", "-l", "-w", fna_file]
    logging.getLogger("PPanGGOLiN").debug(f"aragorn command : {' '.join(cmd)}")
    p = Popen(cmd, stdout=PIPE)
    # loading the whole thing, reverting it to 'pop' in order.
    file_data = p.communicate()[0].decode().split("\n")[:: -1]
    gene_objs = defaultdict(set)
    c = 0
    header = ""
    while len(file_data) != 0:
        line = file_data.pop()
        if line.startswith(">"):
            header = line.replace(">", "").split()[0]
            file_data.pop()  # then next line must be removed too.
        elif len(line) > 0:  # if the line isn't empty, there's data to get.
            line_data = line.split()
            start, stop = map(int, ast.literal_eval(line_data[2].replace("c", "")))
            c += 1
            gene = RNA(rna_id=locustag + '_tRNA_' + str(c).zfill(4))
            gene.fill_annotations(start=start, stop=stop, strand="-" if line_data[2].startswith("c") else "+",
                                  gene_type="tRNA", product=line_data[1] + line_data[4])
            gene_objs[header].add(gene)
    return gene_objs


def launch_prodigal(contig_sequences: Dict[str, str], org: Organism, code: int = 11, use_meta: bool = False) -> defaultdict:
    """
    Launches Prodigal to annotate CDS. Takes a fna file name and a locustag to give an ID to the pred genes.

    :param contig_sequences: Dict containing contig sequences for pyrodigal
    :param org: Organism which will be annotated
    :param code: Translation table (genetic code) to use.
    :param use_meta: use meta procedure in Prodigal

    :return: Annotated genes in a list of gene objects
    """
    gene_objs = defaultdict(set)
    sequences = {contig_name: Sequence(sequence) for contig_name, sequence in contig_sequences.items()}
    gene_finder = GeneFinder(
        meta=use_meta,  # '-p meta' if meta is true else '-p single'
        closed=True,  # -c: Closed ends. Do not allow genes to run off edges.
        mask=True,  # -m: Treat runs of N as masked sequence; don't build genes across them.
        min_gene=120  # This is to prevent erreur with mmseqs translatenucs that cut too short sequences
    )
    gene_finder.train(max(sequences.values(), key=len), force_nonsd=False,
                      translation_table=code)  # -g: Specify a translation table to use (default 11).
    gene_counter = 1
    for contig_name, sequence in sequences.items():
        for pred in gene_finder.find_genes(sequence):
            gene = Gene(gene_id=f"{org.name}_CDS_{str(gene_counter).zfill(4)}")
            gene.fill_annotations(start=pred.begin, stop=pred.end, strand='-' if pred.strand == -1 else '+',
                                  gene_type="CDS", genetic_code=code)
            gene_counter += 1
            gene_objs[contig_name].add(gene)
    return gene_objs


def launch_infernal(fna_file: str, org: Organism, tmpdir: str, kingdom: str = "bacteria") -> defaultdict:
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
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_bact.cm"
    elif kingdom == "archaea":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_arch.cm"

    tmp_file = tempfile.NamedTemporaryFile(mode="r", dir=tmpdir)
    cmd = ["cmscan", "--tblout", tmp_file.name, "--hmmonly", "--cpu", str(1), "--noali", modelfile, fna_file]
    logging.getLogger("PPanGGOLiN").debug(f"infernal command : {' '.join(cmd)}")
    p = Popen(cmd, stdout=open(os.devnull, "w"), stderr=PIPE)
    err = p.communicate()[1].decode().split()
    if err:
        if err[0] == 'Error: ':
            raise Exception(f"Infernal (cmscan) failed with error:  '{' '.join(err)}'. If you never used this script,"
                            f" you should press the .cm file using cmpress executable from Infernal. "
                            f"You should find the file in '{os.path.dirname(os.path.realpath(__file__))}/rRNA_DB/'.")
        raise Exception(f"An error occurred with Infernal. Error is:  '{' '.join(err)}'.")
    # never managed to test what happens if the .cm files are compressed with a 'bad' version of infernal,
    # so if that happens you are on your own.

    gene_objs = defaultdict(set)
    c = 0
    for line in tmp_file:
        if not line.startswith("#"):
            c += 1
            line_data = line.split()
            strand = line_data[9]
            start, stop = map(int, (line_data[8], line_data[7]) if strand == "-" else (line_data[7], line_data[8]))
            gene = RNA(rna_id=locustag + "_rRNA_" + str(c).zfill(4))
            gene.fill_annotations(start=start, stop=stop, strand=strand, gene_type="rRNA",
                                  product=" ".join(line_data[17:]))
            gene_objs[line_data[2]].add(gene)
    return gene_objs


def read_fasta(org: Organism, fna_file: Union[TextIOWrapper, list]) -> Dict[str, str]:
    """ Reads a fna file (or stream, or string) and stores it in a dictionary with contigs as key and sequence as value.

    :param org: Organism corresponding to fasta file
    :param fna_file: Input fasta file with sequences or list of each line as sequence

    :return: Dictionnary with contig_name as keys and contig sequence in values
    """
    global contig_counter

    try:
        contigs = {}
        contig_seq = ""
        contig = None
        for line in fna_file:
            if line.startswith('>'):
                if len(contig_seq) >= 1:  # contig filter = 1
                    contigs[contig.name] = contig_seq.upper()
                    contig.length = len(contig_seq)
                contig_seq = ""
                try:
                    contig = org.get(line.split()[0][1:])
                except KeyError:
                    with contig_counter.get_lock():
                        contig = Contig(contig_counter.value, line.split()[0][1:])
                        contig_counter.value += 1
                    org.add(contig)
            else:
                contig_seq += line.strip()
        if len(contig_seq) >= 1:  # processing the last contig
            contigs[contig.name] = contig_seq.upper()
            contig.length = len(contig_seq)

    except AttributeError as e:
        raise AttributeError(f"{e}\nAn error was raised when reading file: '{fna_file.name}'. "
                             f"One possibility for this error is that the file did not start with a '>' "
                             f"as it would be expected from a fna file.")
    except Exception:  # To manage other exception which can occur
        raise Exception("Unexpected error. Please check your input file and if everything looks fine, "
                        "please post an issue on our github")
    return contigs


def write_tmp_fasta(contigs: dict, tmpdir: str) -> tempfile._TemporaryFileWrapper:
    """
     Writes a temporary fna formated file and returns the file-like object. Useful in case of  compressed input file.
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
            tmp_file.write(contigs[header][j: j + 60] + "\n")
            j += 60
    tmp_file.flush()  # force write what remains in the buffer.
    return tmp_file


def syntaxic_annotation(org: Organism, fasta_file: TextIOWrapper, contig_sequences: Dict[str, str],
                        tmpdir: str, norna: bool = False, kingdom: str = "bacteria",
                        code: int = 11, use_meta: bool = False) -> defaultdict:
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
    for key, items in launch_prodigal(contig_sequences=contig_sequences, org=org, code=code, use_meta=use_meta).items():
        genes[key].extend(items)
    if not norna:
        for key, items in launch_aragorn(fna_file=fasta_file.name, org=org).items():
            genes[key].extend(items)
        for key, items in launch_infernal(fna_file=fasta_file.name, org=org, kingdom=kingdom, tmpdir=tmpdir).items():
            genes[key].extend(items)
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
                    if gene_i.type != "CDS" and gene_j.type == "CDS" and gene_i.stop > gene_j.start:
                        rm_genes.add(gene_j)
                    elif gene_i.type == "CDS" and gene_j.type != "CDS" and gene_i.stop > gene_j.start:
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


def get_dna_sequence(contig_seq: str, gene: Gene) -> str:
    """Return the gene sequence

    :param contig_seq: Contig sequence
    :param gene: Gene

    :return: str
    """
    if gene.strand == "+":
        return contig_seq[gene.start - 1:gene.stop]
    elif gene.strand == "-":
        return reverse_complement(contig_seq[gene.start - 1:gene.stop])


def annotate_organism(org_name: str, file_name: Path, circular_contigs: List[str], tmpdir: str,
                      code: int = 11, norna: bool = False, kingdom: str = "bacteria",
                      allow_overlap: bool = False, procedure: str = None) -> Organism:
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
    org = Organism(org_name)

    fasta_file = read_compressed_or_not(file_name)

    contig_sequences = read_fasta(org, fasta_file)
    if is_compressed(file_name):  # TODO simply copy file with shutil.copyfileobj
        fasta_file = write_tmp_fasta(contig_sequences, tmpdir)
    if procedure is None:  # prodigal procedure is not force by user
        max_contig_len = max(len(contig) for contig in org.contigs)
        if max_contig_len < 20000:  # case of short sequence
            use_meta = True
        else:
            use_meta = False
    else:
        use_meta = True if procedure == "meta" else False
    genes = syntaxic_annotation(org, fasta_file, contig_sequences, tmpdir, norna, kingdom, code, use_meta)
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
