#!/usr/bin/env python3
# coding:utf-8

# default libraries
from _io import TextIOWrapper
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict
from typing import List, Tuple, Set, Dict, IO
from pathlib import Path

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.utils import mk_outdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Spot
from ppanggolin.figures.draw_spot import draw_selected_spots, subgraph


def createdb(file_obj: TextIOWrapper, tmpdir: tempfile.TemporaryDirectory) -> IO:
    """
    Create a MMseqs2 sequence database with the given fasta file

    :param file_obj: Fasta file
    :param tmpdir: temporary directory

    :return: DB file
    """
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name)
    cmd = ["mmseqs", "createdb", file_obj.name, seqdb.name, '--dbtype', '0']
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return seqdb


def align_seq_to_pang(pang_file: IO, seq_file: TextIOWrapper, output: Path,
                      tmpdir: tempfile.TemporaryDirectory, cpu: int = 1, no_defrag: bool = False,
                      identity: float = 0.8, coverage: float = 0.8) -> Path:
    """
    Align pangenome sequences against fasta sequence

    :param pang_file: File with sequences in pangenome
    :param seq_file: File with sequences from input file
    :param output: Path of the output directory
    :param tmpdir: Temporary directory to align sequences
    :param cpu: Number of available cpu
    :param no_defrag: Allow to pass the defragmentation step
    :param identity: minimal identity threshold for the alignment
    :param coverage: minimal identity threshold for the alignment

    :return: Alignement result file
    """

    pang_db = createdb(pang_file, tmpdir)
    seq_db = createdb(seq_file, tmpdir)
    cov_mode = "0"  # coverage of query and target
    if not no_defrag:
        cov_mode = "1"  # coverage of target
    aln_db = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name)
    cmd = ["mmseqs", "search", seq_db.name, pang_db.name, aln_db.name, tmpdir.name, "-a", "--min-seq-id", str(identity),
           "-c", str(coverage), "--cov-mode", cov_mode, "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning sequences to cluster representatives...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    outfile = output.absolute()/"input_to_pangenome_associations.blast-tab_tmp"  # write a tmp file of the results
    cmd = ["mmseqs", "convertalis", seq_db.name, pang_db.name, aln_db.name, outfile.as_posix(), "--format-mode", "2"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    pang_db.close()
    seq_db.close()
    aln_db.close()

    return outfile


def read_alignments(aln_res: Path, pangenome: Pangenome) -> Tuple[Dict[str, GeneFamily], str]:
    """
    Read alignment result to link input sequence to pangenome

    :param aln_res: Alignement result file
    :param pangenome: Input pangenome

    :return: Dictionnary with sequence link to pangenome gene families and actual name of resulting alignment file
    """
    seq2pang = {}
    outname = open(aln_res.absolute().as_posix().replace("_tmp", ""), "w")  # write the actual result file
    with open(aln_res, "r") as alnFile:
        for line in alnFile:
            line = line.replace("ppanggolin_", "")  # remove the 'ppanggolin_' bit of the id
            outname.write(line)
            line = line.split()
            if seq2pang.get(line[0]) is None:  # if no results were found yet
                seq2pang[line[0]] = pangenome.get_gene_family(line[1])  # then the best hit is the first one we see.
    outname.close()
    return seq2pang, outname.name


def get_seq(seq_file: TextIOWrapper) -> Set[str]:
    """
    get sequence from sequence input file

    :param seq_file: file containing sequences

    :return: set of sequences
    """
    seqset = set()
    for line in seq_file:
        if line.startswith(">"):
            seqset.add(line[1:])
    return seqset


def write_gene_fam_sequences(pangenome: Pangenome, file_obj: IO, add: str = ""):
    """
    Export the sequence of genes in families

    :param pangenome: Pangenome containing families
    :param file_obj: Temporary file where sequences will be written
    :param add: Add prefix to sequence name
    """
    for fam in pangenome.gene_families:
        file_obj.write(">" + add + fam.name + "\n")
        file_obj.write(fam.sequence + "\n")
    file_obj.flush()


def project_partition(seq_to_pang: Dict[str, GeneFamily], seq_set: Set[str], output: Path) -> str:
    """
    Project the partition of each sequence from the input file

    :param seq_to_pang: dictionnary which link sequence and pangenome
    :param seq_set: input sequences
    :param output: Path of the output directory

    :return: Path to file which contain partition projection
    """

    partition_proj = output.absolute()/"sequences_partition_projection.tsv"
    with open(partition_proj, "w") as partProjFile:
        for key, pang_fam in seq_to_pang.items():
            partProjFile.write(key + "\t" + pang_fam.named_partition + "\n")
        for remaining_seq in (seq_to_pang.keys() & seq_set):
            partProjFile.write(remaining_seq + "\tcloud\n")  # if there is no hit, it's going to be cloud genes.
    return partition_proj


def get_fam_to_rgp(pangenome, multigenics: set) -> dict:
    """
    Associate families to the RGP they belong to, and those they are bordering

    :param pangenome: Input pangenome
    :param multigenics: multigenics families

    :return: Dictionnary link families to RGP
    """
    fam2rgp = defaultdict(list)
    for rgp in pangenome.regions:
        for fam in rgp.families:
            fam2rgp[fam].append(rgp.name)
        for fam in [gene.family for border in rgp.get_bordering_genes(pangenome.parameters["spots"]["set_size"],
                                                                      multigenics) for gene in border]:
            fam2rgp[fam].append(rgp.name)
    return fam2rgp


def get_fam_to_spot(pangenome: Pangenome, multigenics: Set[GeneFamily]) -> Tuple[Dict[str, List[Spot]],
                                                                                 Dict[str, List[Spot]]]:
    """
    Reads a pangenome object to link families and spots and indicate where each family is.

    :param pangenome: Input pangenome
    :param multigenics: multigenics families

    :return: Dictionary of family to RGP and family to spot
    """
    # those are to be replaced as spots should be stored in the pangenome, and in the h5.
    fam2spot = defaultdict(list)
    fam2border = defaultdict(list)
    for spot in pangenome.spots:
        fams = set()
        fams_border = set()
        for rgp in spot.regions:
            fams |= rgp.families
            fams_border |= set([gene.family for border in  # Set of families in border of spot
                                rgp.get_bordering_genes(pangenome.parameters["spots"]["set_size"], multigenics)
                                for gene in border])
        for fam in fams:
            fam2spot[fam].append(spot)
        for fam in fams_border:
            fam2border[fam].append(spot)
    return fam2spot, fam2border


def add_spot_str(spot: Spot) -> str:
    # TODO define as self.__str__ in spot
    """
    allow to map spot set

    :param spot: spot which will be return

    :return: Str with spot ID
    """
    return "spot_" + str(spot.ID)


def draw_spot_gexf(spots: set, output: Path, multigenics: set, fam_to_mod: dict, set_size: int = 3):
    """
    Draw a gexf graph of the spot

    :param spots: spot find in the alignment between pangenome and input sequences
    :param output: Path of the output directory
    :param multigenics: multigenics families
    :param fam_to_mod: dictionnary which link families and modules
    :param set_size:
    """
    for spot in spots:
        fname = output/f"spot_{str(spot.ID)}.gexf"
        subgraph(spot, fname, set_size=set_size, multigenics=multigenics, fam_to_mod=fam_to_mod)


def get_seq_info(seq_to_pang: dict, pangenome: Pangenome, output: Path, draw_related: bool = False, disable_bar=False):
    """
    Get sequences information after alignment

    :param seq_to_pang: Alignment result
    :param pangenome: Pangenome which contain information
    :param output: Path of the output directory
    :param draw_related: Draw figures and graphs in a gexf format of spots associated to the input sequences
    :param disable_bar: disable progress bar
    :return:
    """
    logging.getLogger().info("Writing RGP and spot information related to hits in the pangenome")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    finfo = open(output/"info_input_seq.tsv", "w")
    finfo.write("input\tfamily\tpartition\tspot_list_as_member\tspot_list_as_border\trgp_list\n")
    fam2rgp = get_fam_to_rgp(pangenome, multigenics)
    fam2spot, fam2border = get_fam_to_spot(pangenome, multigenics)
    spot_list = set()
    for seq, panfam in seq_to_pang.items():
        finfo.write(seq + '\t' + panfam.name + "\t" + panfam.named_partition + "\t" + ",".join(
            map(add_spot_str, fam2spot[panfam])) + "\t" + ",".join(
            map(add_spot_str, fam2border[panfam])) + "\t" + ','.join(fam2rgp[panfam]) + "\n")
        spot_list |= set(fam2spot[panfam])
        spot_list |= set(fam2border[panfam])
    finfo.close()
    if draw_related:
        drawn_spots = set()
        for spot in spot_list:
            if len(spot.get_uniq_ordered_set()) > 1:
                drawn_spots.add(spot)
        logging.getLogger().info(f"Drawing the {len(drawn_spots)} spots with more than 1 organization "
                                 f"related to hits of the input sequences...")
        draw_selected_spots(drawn_spots, pangenome, output, pangenome.parameters["spots"]["overlapping_match"],
                            pangenome.parameters["spots"]["exact_match"], pangenome.parameters["spots"]["set_size"],
                            disable_bar=disable_bar)

        fam2mod = {}  # fam2module
        if pangenome.status["modules"] != "No":
            for mod in pangenome.modules:
                for fam in mod.families:
                    fam2mod[fam] = f"module_{mod.ID}"

        draw_spot_gexf(drawn_spots, output, multigenics=multigenics, fam_to_mod=fam2mod)

    logging.getLogger().info(f"File listing RGP and spots where sequences of interest are located : "
                             f"{output/'info_input_seq.tsv'}")


def get_seq2pang(pangenome: Pangenome, sequence_file: Path, output: Path, tmpdir: tempfile.TemporaryDirectory,
                 cpu: int = 1, no_defrag: bool = False, identity: float = 0.8,
                 coverage: float = 0.8) -> Tuple[set, str, dict]:
    """
    Assign a pangenome gene family to the input sequences.

    :param pangenome: Pangenome with gene families to align with the given input sequences
    :param sequence_file: Sequences in a .fasta file to align with the given Pangenome
    :param output: Path of the output directory
    :param tmpdir: Temporary directory
    :param cpu: number of CPU cores to use
    :param no_defrag: do not use the defrag workflow if true
    :param identity: minimal identity threshold for the alignment
    :param coverage: minimal identity threshold for the alignment

    :return: sequence set, blast-tab result file string, and sequences aligned with families
    """
    tmp_pang_file = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name)

    write_gene_fam_sequences(pangenome, tmp_pang_file, add="ppanggolin_")

    with read_compressed_or_not(sequence_file) as seqFileObj:
        seq_set = get_seq(seqFileObj)
        align_file = align_seq_to_pang(tmp_pang_file, seqFileObj, output, tmpdir, cpu, no_defrag, identity, coverage)

    seq2pang, align_file = read_alignments(align_file, pangenome)

    tmp_pang_file.close()

    return seq_set, align_file, seq2pang


def align(pangenome: Pangenome, sequence_file: Path, output: Path, tmpdir: Path, identity: float = 0.8,
          coverage: float = 0.8, no_defrag: bool = False, cpu: int = 1, getinfo: bool = False,
          draw_related: bool = False, disable_bar: bool = False):
    """
    Main function to align pangenome sequences with fasta file using MMSeqs2

    :param pangenome: Pangenome with gene families to align with the given input sequences
    :param sequence_file: Sequences in a .fasta file to align with the given Pangenome
    :param output: Path of the output directory
    :param tmpdir: Temporary directory
    :param identity: minimal identity threshold for the alignment
    :param coverage: minimal coverage threshold for the alignment
    :param no_defrag: do not use the defrag workflow if true
    :param cpu: number of CPU cores to use
    :param getinfo: Extract info related to the best hit of each query, such as the RGP it is in, or the spots.
    :param draw_related: Draw figures and graphs in a gexf format of spots associated to the input sequences
    :param disable_bar: Disable the progresse bar
    """

    if pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    # could be possible either by picking a representative somehow, or by aligning on genes rather than on
    # families, if they are in the pangenome.

    if getinfo or draw_related:
        need_mod = False
        if pangenome.status["modules"] != "No":
            # modules are not required to be loaded, but if they have been computed we load them.
            need_mod = True
        check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_partitions=True, need_rgp=True,
                             need_spots=True, need_modules=need_mod, disable_bar=disable_bar)
    else:
        check_pangenome_info(pangenome, need_families=True, disable_bar=disable_bar)

    new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)

    seq_set, align_file, seq2pang = get_seq2pang(pangenome, sequence_file, output, new_tmpdir, cpu, no_defrag, identity,
                                                 coverage)

    if getinfo or draw_related:  # TODO Add getinfo to function and remove if
        get_seq_info(seq2pang, pangenome, output, draw_related, disable_bar=disable_bar)
    part_proj = project_partition(seq2pang, seq_set, output)  # write the partition assignation only
    logging.getLogger().info(f"sequences partition projection : '{part_proj}'")
    logging.getLogger().info(f"{len(seq2pang)} sequences over {len(seq_set)} have at least one hit in the pangenome.")
    logging.getLogger().info(f"Blast-tab file of the alignment : '{align_file}'")

    new_tmpdir.cleanup()


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    align(pangenome=pangenome, sequence_file=args.sequences, output=args.output, tmpdir=args.tmpdir,
          identity=args.identity, coverage=args.coverage, no_defrag=args.no_defrag, cpu=args.cpu, getinfo=args.getinfo,
          draw_related=args.draw_related, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("align", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_align(parser)
    return parser


def parser_align(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of align command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-S', '--sequences', required=True, type=Path,
                          help="sequences (nucleotides or amino acids) to align on the pangenome gene families")

    required.add_argument('-p', '--pangenome', required=True, type=Path, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family. (default: False)")
    optional.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("--translation_table", required=False, default="11",
                          help="Translation table (genetic code) to use.")
    optional.add_argument("--getinfo", required=False, action="store_true",
                          help="Use this option to extract info related to the best hit of each query, "
                               "such as the RGP it is in, or the spots.")
    optional.add_argument("--draw_related", required=False, action="store_true",
                          help="Draw figures and provide graphs in a gexf format of the eventual spots"
                               " associated to the input sequences")
    # but does not use the option
    optional.add_argument("--use_pseudo", required=False, action="store_true",
                          help="In the context of provided annotation, use this option to read pseudogenes. "
                               "(Default behavior is to ignore them)")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_align(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--tmpdir", required=False, type=Path, default=Path(tempfile.gettempdir()),
                        help="directory for storing temporary files")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
