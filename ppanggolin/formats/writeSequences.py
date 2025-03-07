#!/usr/bin/env python3

# default libraries
import argparse
import logging
import re
from pathlib import Path
from typing import Dict, Iterable, Union
import tempfile
import shutil

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene, Organism

from ppanggolin.utils import (
    write_compressed_or_not,
    mk_outdir,
    create_tmpdir,
    read_compressed_or_not,
    restricted_float,
    detect_filetype,
    run_subprocess,
    check_tools_availability,
)
from ppanggolin.formats.readBinaries import (
    check_pangenome_info,
    write_genes_from_pangenome_file,
    write_fasta_gene_fam_from_pangenome_file,
    write_fasta_prot_fam_from_pangenome_file,
)

module_regex = re.compile(r"^module_\d+")  # \d == [0-9]
poss_values = [
    "all",
    "persistent",
    "shell",
    "cloud",
    "rgp",
    "softcore",
    "core",
    module_regex,
]
poss_values_log = f"Possible values are {', '.join(poss_values[:-1])}, module_X with X being a module id."


def check_write_sequences_args(args: argparse.Namespace) -> None:
    """Check arguments compatibility in CLI

    :param args: argparse namespace arguments from CLI

    :raises argparse.ArgumentTypeError: if region is given but neither fasta nor anno is given
    """
    if args.regions is not None and args.fasta is None and args.anno is None:
        raise argparse.ArgumentError(
            argument=None,
            message="The --regions options requires the use of --anno or --fasta "
            "(You need to provide the same file used to compute the pangenome)",
        )


def write_gene_sequences_from_annotations(
    genes_to_write: Iterable[Gene],
    output: Path,
    add: str = "",
    compress: bool = False,
    disable_bar: bool = False,
):
    """
    Writes the CDS sequences to a File object,
    and adds the string provided through `add` in front of it.
    Loads the sequences from previously computed or loaded annotations.

    :param genes_to_write: Genes to write.
    :param output: Path to output file to write sequences.
    :param add: Add prefix to gene ID.
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar.
    """
    logging.getLogger("PPanGGOLiN").info(
        f"Writing all CDS sequences in {output.absolute()}"
    )
    with write_compressed_or_not(output, compress) as file_obj:
        for gene in tqdm(genes_to_write, unit="gene", disable=disable_bar):
            if gene.type == "CDS":
                file_obj.write(f">{add}{gene.ID}\n")
                file_obj.write(f"{gene.dna}\n")


def create_mmseqs_db(
    sequences: Iterable[Path],
    db_name: str,
    tmpdir: Path,
    db_mode: int = 0,
    db_type: int = 0,
) -> Path:
    """Create a MMseqs2 database from a sequences file.

    :param sequences: File with the sequences
    :param db_name: name of the database
    :param tmpdir: Temporary directory to save the MMSeqs2 files
    :param db_mode: Createdb mode 0: copy data, 1: soft link data and write new index (works only with single line fasta/q)
    :param db_type: Database type 0: auto, 1: amino acid 2: nucleotides

    :return: Path to the MMSeqs2 database
    """
    assert db_mode in [0, 1], f"Createdb mode must be 0 or 1, given {db_mode}"
    assert db_type in [0, 1, 2], f"dbtype must be 0, 1 or 2, given {db_type}"

    seq_nucdb = tmpdir / db_name
    cmd = ["mmseqs", "createdb", "--createdb-mode", db_mode, "--dbtype", db_type]
    cmd += [seq.absolute().as_posix() for seq in sequences] + [seq_nucdb.absolute()]
    cmd = list(map(str, cmd))
    logging.getLogger("PPanGGOLiN").info("Creating sequence database...")
    run_subprocess(cmd, msg="MMSeqs createdb failed with the following error:\n")
    return seq_nucdb


def translate_genes(
    sequences: Union[Path, Iterable[Path]],
    tmpdir: Path,
    cpu: int = 1,
    is_single_line_fasta: bool = False,
    code: int = 11,
) -> Path:
    """Translate nucleotide sequences into MMSeqs2 amino acid sequences database

    :param sequences: File with the nucleotide sequences
    :param tmpdir: Temporary directory to save the MMSeqs2 files
    :param cpu: Number of available threads to use
    :param is_single_line_fasta: Allow to use soft link in MMSeqs2 database
    :param code: Translation code to use

    :return: Path to the MMSeqs2 database
    """
    check_tools_availability(["mmseqs"])

    seq_nucdb = create_mmseqs_db(
        [sequences] if isinstance(sequences, Path) else sequences,
        "nucleotides_db",
        tmpdir,
        db_mode=1 if is_single_line_fasta else 0,
        db_type=2,
    )
    logging.getLogger("PPanGGOLiN").debug("Translate sequence ...")
    seqdb = tmpdir / "translate_db"
    cmd = list(
        map(
            str,
            [
                "mmseqs",
                "translatenucs",
                seq_nucdb,
                seqdb,
                "--threads",
                cpu,
                "--translation-table",
                code,
            ],
        )
    )
    run_subprocess(cmd, msg="MMSeqs translatenucs failed with the following error:\n")
    return seqdb


def write_gene_protein_sequences(
    pangenome_filename: str,
    output: Path,
    gene_filter: str,
    soft_core: float = 0.95,
    compress: bool = False,
    keep_tmp: bool = False,
    tmp: Path = None,
    cpu: int = 1,
    code: int = 11,
    disable_bar: bool = False,
):
    """Write all amino acid sequences from given genes in pangenome

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param proteins: Selected partition of gene
    :param soft_core: Soft core threshold to use
    :param compress: Compress the file in .gz
    :param keep_tmp: Keep temporary directory
    :param tmp: Path to temporary directory
    :param cpu: Number of threads available
    :param code: Genetic code use to translate nucleotide sequences to protein sequences
    :param disable_bar: Disable progress bar
    """

    check_tools_availability(["mmseqs"])

    with create_tmpdir(
        tmp if tmp is not None else Path(tempfile.gettempdir()),
        basename="translateGenes",
        keep_tmp=keep_tmp,
    ) as tmpdir:

        write_genes_from_pangenome_file(
            pangenome_filename=pangenome_filename,
            gene_filter=gene_filter,
            output=tmpdir,
            compress=compress,
            disable_bar=disable_bar,
        )

        genes_sequence_tmp_file = (
            tmpdir / f"{gene_filter}_genes.fna{'.gz' if compress else ''}"
        )
        translate_db = translate_genes(
            sequences=genes_sequence_tmp_file,
            tmpdir=tmpdir,
            cpu=cpu,
            is_single_line_fasta=True,
            code=code,
        )

        outpath = output / f"{gene_filter}_protein_genes.faa"

        logging.getLogger("PPanGGOLiN").info(
            "Translating nucleotide gene sequence in protein sequences with mmseqs convert2fasta"
        )

        cmd = list(map(str, ["mmseqs", "convert2fasta", translate_db, outpath]))
        run_subprocess(
            cmd, msg="MMSeqs convert2fasta failed with the following error:\n"
        )
        if compress:
            with write_compressed_or_not(outpath, compress) as compress_file:
                with open(outpath) as sequence_file:
                    shutil.copyfileobj(sequence_file, compress_file)
            outpath.unlink()
            logging.getLogger("PPanGGOLiN").info(
                f"Done writing the gene protein sequences : '{outpath}.gz'"
            )
        else:
            logging.getLogger("PPanGGOLiN").info(
                f"Done writing the gene protein sequences : '{outpath}'"
            )


def read_fasta_or_gff(file_path: Path) -> Dict[str, str]:
    """
    Read the genome file in fasta or gbff format

    :param file_path: Path to genome file

    :return: Dictionary with all sequences associated to contig
    """
    sequence_dict = {}
    seqname = ""
    seq = ""
    in_fasta_part = False
    with read_compressed_or_not(file_path) as f:
        for line in f:
            if line.startswith(">"):
                in_fasta_part = True
            if in_fasta_part:
                if line.startswith(">"):
                    if seq != "":
                        sequence_dict[seqname] = seq
                        seq = ""
                    seqname = line[1:].strip().split()[0]
                else:
                    seq += line.strip()
        if seq != "":
            sequence_dict[seqname] = seq
    return sequence_dict


def read_fasta_gbk(file_path: Path) -> Dict[str, str]:
    """
    Read the genome file in gbk format

    :param file_path: Path to genome file

    :return: Dictionary with all sequences associated to contig
    """
    # line.startswith("ORIGIN"):
    sequence_dict = {}
    lines = read_compressed_or_not(file_path).readlines()[::-1]
    contig_id, contig_locus_id = ("", "")
    while len(lines) != 0:
        line = lines.pop()
        # beginning of contig
        if line.startswith("LOCUS"):
            contig_locus_id = line.split()[1]
            # If contig_id is not specified in VERSION afterward like with Prokka,
            # in that case we use the one in LOCUS.
            while not line.startswith("FEATURES"):
                if line.startswith("VERSION"):
                    contig_id = line[12:].strip()
                line = lines.pop()
        if contig_id == "":
            contig_id = contig_locus_id
        while not line.startswith("ORIGIN"):
            line = lines.pop()  # stuff
        line = lines.pop()  # first sequence line.
        sequence = ""
        while not line.startswith("//"):
            sequence += line[10:].replace(" ", "").strip().upper()
            line = lines.pop()
        # get each gene's sequence.
        sequence_dict[contig_id] = sequence
        # end of contig
    return sequence_dict


def read_genome_file(genome_file: Path, organism: Organism) -> Dict[str, str]:
    """
    Read the genome file associated to organism to extract sequences

    :param genome_file: Path to a fasta file or gbff/gff file
    :param organism: organism object

    :return: Dictionary with all sequences associated to contig

    :raises TypeError: If the file containing sequences is not recognized
    :raises KeyError: If their inconsistency between pangenome contigs and the given contigs
    """
    filetype = detect_filetype(genome_file)
    if filetype in ["fasta", "gff"]:
        contig_to_sequence = read_fasta_or_gff(genome_file)
    elif filetype == "gbff":
        contig_to_sequence = read_fasta_gbk(genome_file)
    else:
        raise TypeError(f"Unknown filetype detected: '{genome_file}'")

    # check_contig_names
    if set(contig_to_sequence) != {contig.name for contig in organism.contigs}:
        raise KeyError(
            f"Contig name inconsistency detected in genome '{organism.name}' between the "
            f"information stored in the pangenome file and the contigs found in '{genome_file}'."
        )

    return contig_to_sequence


def write_spaced_fasta(sequence: str, space: int = 60) -> str:
    """
    Write a maximum of element per line

    :param sequence: sequence to write
    :param space: maximum of size for one line

    :return: a sequence of maximum space character
    """
    seq = ""
    j = 0
    while j < len(sequence):
        seq += sequence[j : j + space] + "\n"
        j += space
    return seq


def write_regions_sequences(
    pangenome: Pangenome,
    output: Path,
    regions: str,
    fasta: Path = None,
    anno: Path = None,
    compress: bool = False,
    disable_bar: bool = False,
):
    """
    Write representative amino acid sequences of gene families.

    :param pangenome: Pangenome object with gene families sequences
    :param output: Path to output directory
    :param regions: Write the RGP nucleotide sequences
    :param fasta: A tab-separated file listing the organism names, fasta filepath of its genomic sequences
    :param anno: A tab-separated file listing the organism names, and the gff/gbff filepath of its annotations
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar

    :raises SyntaxError: if no tabulation are found in list genomes file
    """
    assert (
        fasta is not None or anno is not None
    ), "Write regions requires to use anno or fasta, not any provided"

    organisms_file = fasta if fasta is not None else anno
    org_dict = {}
    for line in read_compressed_or_not(organisms_file):
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            raise ValueError(
                f"No tabulation separator found in given --fasta or --anno file: '{organisms_file}'"
            )
        org_dict[elements[0]] = Path(elements[1])
        if not org_dict[
            elements[0]
        ].exists():  # Check tsv sanity test if it's not one it's the other
            org_dict[elements[0]] = organisms_file.parent.joinpath(
                org_dict[elements[0]]
            )

    logging.getLogger("PPanGGOLiN").info(f"Writing {regions} rgp genomic sequences...")

    if regions == "complete":
        regions_to_write = (
            region for region in pangenome.regions if not region.is_contig_border
        )
    else:
        regions_to_write = pangenome.regions

    regions_to_write = sorted(regions_to_write, key=lambda x: x.organism.name)
    # order regions by organism, so that we only have to read one genome at the time

    outname = output / f"{regions}_rgp_genomic_sequences.fasta"
    with write_compressed_or_not(outname, compress) as fasta:
        loaded_genome = ""
        for region in tqdm(regions_to_write, unit="rgp", disable=disable_bar):
            if region.organism.name != loaded_genome:
                organism = region.organism
                genome_sequence = read_genome_file(org_dict[organism.name], organism)
            fasta.write(f">{region.name}\n")
            fasta.write(
                write_spaced_fasta(
                    genome_sequence[region.contig.name][region.start : region.stop], 60
                )
            )
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing the regions nucleotide sequences: "
        f"'{outname}{'.gz' if compress else ''}'"
    )


def write_sequence_files(
    pangenome: Pangenome,
    output: Path,
    fasta: Path = None,
    anno: Path = None,
    soft_core: float = 0.95,
    regions: str = None,
    genes: str = None,
    proteins: str = None,
    gene_families: str = None,
    prot_families: str = None,
    compress: bool = False,
    disable_bar: bool = False,
    **translate_kwgs,
):
    """
    Main function to write sequence file from pangenome

    :param pangenome: Pangenome object containing sequences
    :param output: Path to output directory
    :param fasta: A tab-separated file listing the organism names, fasta filepath of its genomic sequences
    :param anno: A tab-separated file listing the organism names, and the gff/gbff filepath of its annotations
    :param soft_core: Soft core threshold to use
    :param regions: Write the RGP nucleotide sequences
    :param genes: Write all nucleotide CDS sequences
    :param proteins: Write amino acid CDS sequences.
    :param gene_families: Write representative nucleotide sequences of gene families.
    :param prot_families: Write representative amino acid sequences of gene families.
    :param compress: Compress the file in .gz
    :param disable_bar: Disable progress bar
    """

    if gene_families is not None:

        logging.getLogger("PPanGGOLiN").info(
            "Writing the representative nucleotide sequences "
            "of the gene families by reading the pangenome file directly."
        )

        write_fasta_gene_fam_from_pangenome_file(
            pangenome_filename=pangenome.file,
            family_filter=gene_families,
            soft_core=soft_core,
            output=output,
            compress=compress,
            disable_bar=disable_bar,
        )
        gene_families = None

    if prot_families is not None:

        logging.getLogger("PPanGGOLiN").info(
            "Writing the representative protein sequences "
            "of the gene families by reading the pangenome file directly."
        )
        write_fasta_prot_fam_from_pangenome_file(
            pangenome_filename=pangenome.file,
            family_filter=prot_families,
            soft_core=soft_core,
            output=output,
            compress=compress,
            disable_bar=disable_bar,
        )

        prot_families = None

    if genes is not None:

        logging.getLogger("PPanGGOLiN").info(
            "Writing gene nucleotide sequences by reading the pangenome file directly."
        )
        write_genes_from_pangenome_file(
            pangenome_filename=pangenome.file,
            gene_filter=genes,
            soft_core=soft_core,
            output=output,
            compress=compress,
            disable_bar=disable_bar,
        )

        genes = None

    if proteins is not None:

        logging.getLogger("PPanGGOLiN").info(
            "Writing gene protein sequences by reading the pangenome file directly."
        )
        write_gene_protein_sequences(
            pangenome_filename=pangenome.file,
            output=output,
            gene_filter=proteins,
            soft_core=soft_core,
            compress=compress,
            disable_bar=disable_bar,
            **translate_kwgs,
        )

        proteins = None

    if regions is not None:
        # load pangenome when writing region sequence
        check_pangenome_info(
            pangenome,
            need_annotations=True,
            need_families=True,
            need_rgp=True,
            disable_bar=disable_bar,
        )

        write_regions_sequences(
            pangenome, output, regions, fasta, anno, compress, disable_bar
        )


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    check_write_sequences_args(args)
    translate_kwgs = {
        "code": args.translation_table,
        "cpu": args.cpu,
        "tmp": args.tmpdir,
        "keep_tmp": args.keep_tmp,
    }
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    write_sequence_files(
        pangenome,
        args.output,
        fasta=args.fasta,
        anno=args.anno,
        soft_core=args.soft_core,
        regions=args.regions,
        genes=args.genes,
        proteins=args.proteins,
        gene_families=args.gene_families,
        prot_families=args.prot_families,
        compress=args.compress,
        disable_bar=args.disable_prog_bar,
        **translate_kwgs,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "fasta", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_seq(parser)
    return parser


def filter_values(arg_value: str):
    """
    Check filter value to ensure they are in the expected format.

    :param arg_value: Argument value that is being tested.

    :return: The same argument if it is valid.

    :raises argparse.ArgumentTypeError: If the argument value is not in the expected format.
    """
    if arg_value in poss_values or module_regex.match(arg_value):
        return arg_value
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid choice '{arg_value}'. {poss_values_log}"
        )


def parser_seq(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of fasta command

    :param parser: parser for align argument
    """

    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory where the file(s) will be written",
    )

    context = parser.add_argument_group(
        title="Contextually required arguments",
        description="With --regions, the following arguments are required:",
    )
    context.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the fasta filepath of its genomic "
        "sequence(s) (the fastas can be compressed with gzip). One line per genome.",
    )
    context.add_argument(
        "--anno",
        required=False,
        type=Path,
        help="A tab-separated file listing the genome names, and the gff/gbff filepath of its "
        "annotations (the files can be compressed with gzip). One line per genome. "
        "If this is provided, those annotations will be used.",
    )

    onereq = parser.add_argument_group(
        title="Output file",
        description="At least one of the following argument is required. "
        "Indicating 'all' writes all elements. Writing a partition "
        "('persistent', 'shell' or 'cloud') write the elements associated "
        "to said partition. Writing 'rgp' writes elements associated to RGPs",
    )
    onereq.add_argument(
        "--genes",
        required=False,
        type=filter_values,
        help=f"Write all nucleotide CDS sequences. {poss_values_log}",
    )
    onereq.add_argument(
        "--proteins",
        required=False,
        type=filter_values,
        help=f"Write representative amino acid sequences of genes. {poss_values_log}",
    )
    onereq.add_argument(
        "--prot_families",
        required=False,
        type=filter_values,
        help=f"Write representative amino acid sequences of gene families. {poss_values_log}",
    )
    onereq.add_argument(
        "--gene_families",
        required=False,
        type=filter_values,
        help=f"Write representative nucleotide sequences of gene families. {poss_values_log}",
    )

    optional = parser.add_argument_group(title="Optional arguments")
    # could make choice to allow customization
    optional.add_argument(
        "--regions",
        required=False,
        type=str,
        choices=["all", "complete"],
        help="Write the RGP nucleotide sequences (requires --anno or --fasta used to compute "
        "the pangenome to be given)",
    )
    optional.add_argument(
        "--soft_core",
        required=False,
        type=restricted_float,
        default=0.95,
        help="Soft core threshold to use if 'softcore' partition is chosen",
    )
    optional.add_argument(
        "--compress",
        required=False,
        action="store_true",
        help="Compress the files in .gz",
    )
    optional.add_argument(
        "--translation_table",
        required=False,
        default="11",
        help="Translation table (genetic code) to use.",
    )
    optional.add_argument(
        "--cpu", required=False, default=1, type=int, help="Number of available threads"
    )
    optional.add_argument(
        "--tmpdir",
        required=False,
        type=Path,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
    )
    optional.add_argument(
        "--keep_tmp",
        required=False,
        default=False,
        action="store_true",
        help="Keeping temporary files (useful for debugging).",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_seq(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
