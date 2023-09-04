#!/usr/bin/env python3
# coding:utf-8

# default libraries
from _io import TextIOWrapper
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict
from typing import List, Tuple, Set, Dict, IO, Iterator
from pathlib import Path
import time

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.utils import mk_outdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Spot
from ppanggolin.figures.draw_spot import draw_selected_spots, subgraph
from ppanggolin.formats.readBinaries import get_gene_sequences_from_file

def createdb(file_obj: TextIOWrapper, tmpdir: Path, delete_tmp_file: bool = True) -> IO:
    """
    Create a MMseqs2 sequence database with the given fasta file

    :param file_obj: Fasta file
    :param tmpdir: temporary directory

    :return: DB file
    """
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=False)
    cmd = ["mmseqs", "createdb", file_obj.name, seqdb.name, '--dbtype', '0']
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return seqdb


def align_seq_to_pang(pang_file: IO, seq_file: TextIOWrapper, output: Path,
                      tmpdir: Path, cpu: int = 1, no_defrag: bool = False,
                      identity: float = 0.8, coverage: float = 0.8, is_nucleotid:bool = False, translation_table: int = None) -> Path:
    """
    Align fasta sequence to pangenome sequences. 

    :param pang_file: File with sequences in pangenome
    :param seq_file: File with sequences from input file
    :param output: Path of the output directory
    :param tmpdir: Temporary directory to align sequences
    :param cpu: Number of available cpu
    :param no_defrag: Allow to pass the defragmentation step
    :param identity: minimal identity threshold for the alignment
    :param coverage: minimal identity threshold for the alignment
    :param is_nucleotid: Is the sequence file are nucleotide sequences. If True, sequences are translated by mmseqs
    :param translation_table: Translation table to use, if sequences are nucleotide and need to be translated.

    :return: Alignement result file
    """
    translate_first = True
    pang_file_is_nt = True

    pang_db = createdb(pang_file, tmpdir)
    seq_db = createdb(seq_file, tmpdir)

    if pang_file_is_nt:
        logging.getLogger().debug(f"Pangenomes sequences will be translated by mmseqs with translation table {translation_table}")

        pangdb_aa = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=False, prefix="pangenome_seq", suffix=".aa.DB")

        cmd = ["mmseqs", "translatenucs", pang_db.name, pangdb_aa.name, "--translation-table", 
               f"{translation_table}", "--threads", str(cpu)]
        logging.getLogger().debug(" ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        pang_db = pangdb_aa

    if translate_first and is_nucleotid:
        logging.getLogger().debug(f"Input sequences will be translated by mmseqs with translation table {translation_table}")

        seqdb_aa = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=False, prefix="input_seq", suffix=".aa.DB")

        cmd = ["mmseqs", "translatenucs", seq_db.name, seqdb_aa.name, "--translation-table",
               f"{translation_table}", "--threads", str(cpu)]
        
        logging.getLogger().debug(" ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        seq_db = seqdb_aa

    cov_mode = "1"  # coverage of target
    if no_defrag:    
       cov_mode = "0"  # coverage of query and target

    with tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.as_posix(), prefix="aln_result_db_file", delete=False) as aln_db:
        cmd = ["mmseqs", "search", seq_db.name, pang_db.name, aln_db.name, tmpdir.as_posix(), "-a", "--min-seq-id", str(identity),
            "-c", str(coverage), "--cov-mode", cov_mode, "--threads", str(cpu),  "--max-seqs", str(1)]
        if not translate_first and is_nucleotid:
            logging.getLogger().debug(f"Input sequences will be translated by mmseqs search with translation table {translation_table}")
            cmd += ["--translation-table", f"{translation_table}", "--translate", "0" ]
    
        
        logging.getLogger().info("Aligning sequences")# to cluster representatives...")
        logging.getLogger().debug(" ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)

        with tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, prefix="aln_result_db_file", suffix = ".tsv", delete=False) as outfile:
            cmd = ["mmseqs", "convertalis", seq_db.name, pang_db.name, aln_db.name, outfile.name, "--format-mode", "2"]

            logging.getLogger().info("Extracting alignments...")
            logging.getLogger().debug(" ".join(cmd))
            subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)


    pang_db.close()
    seq_db.close()

    return outfile.name


def associate_input_seq_to_gene_family_from_aln_all(aln_res: Path, outdir:Path, pangenome: Pangenome) -> Tuple[Dict[str, GeneFamily], str]:
    """
    Read alignment result to link input sequences to pangenome gene family. 
    Alignment have been made against all genes of the pangenome.

    :param aln_res: Alignement result file
    :param outdir: Output directory
    :param pangenome: Input pangenome

    :return: Dictionnary with sequence link to pangenome gene families and actual path to the cleaned alignment file
    """

    seq2pang = {}
    result_file = outdir / f"alignment_input_seqs_to_pangenome_gene_families.tsv"  # write the actual result file 
    logging.getLogger(f'Get write alignment file in {result_file}')

    with open(aln_res, "r") as alnFile, open(result_file, "w") as outfile :
        for line in alnFile:
            line_splitted = line.split()
            
            line_splitted[1] = line_splitted[1].replace("ppanggolin_", "")  # remove the 'ppanggolin_' bit of the id

            outfile.write("\t".join(line_splitted) + "\n")

            input_seq_id, gene_id = line_splitted[0:2]

            if seq2pang.get(input_seq_id) is None:  # if no results were found yet
                seq2pang[input_seq_id] = pangenome.get_gene(gene_id).family  # then the best hit is the first one we see.

    return seq2pang, outfile


def associate_input_seq_to_gene_family_from_aln_rep(aln_res: Path, outdir:Path, pangenome: Pangenome) -> Tuple[Dict[str, GeneFamily], str]:
    """
    Read alignment result to link input sequences to pangenome gene family.
    Alignment have been made against representative sequence of gene families of the pangenome.

    :param aln_res: Alignement result file
    :param outdir: Output directory
    :param pangenome: Input pangenome

    :return: Dictionnary with sequence link to pangenome gene families and actual path to the cleaned alignment file
    """
    seq2pang = {}
    result_file = outdir / f"alignment_input_seqs_to_pangenome_gene_families.tsv"  # write the actual result file 
    logging.getLogger(f'Get write alignment file in {result_file}')

    with open(aln_res, "r") as alnFile, open(result_file, "w") as outfile :
        for line in alnFile:
            line_splitted = line.split()
            
            line_splitted[1] = line_splitted[1].replace("ppanggolin_", "")  # remove the 'ppanggolin_' bit of the id

            outfile.write("\t".join(line_splitted) + "\n")

            input_seq_id, gene_family_id = line_splitted[0:2]

            if seq2pang.get(input_seq_id) is None:  # if no results were found yet
                seq2pang[input_seq_id] = pangenome.get_gene_family(gene_family_id)  # then the best hit is the first one we see.

    return seq2pang, outfile


def get_seq(seq_file: TextIOWrapper) -> Set[str]:
    """
    get sequence if from sequence input file in fasta format

    :param seq_file: file containing sequences

    :return: set of sequences
    """
    seqset = set()
    for line in seq_file:
        if line.startswith(">"):
            seqset.add(line[1:].split()[0].strip())
    return seqset


def write_gene_fam_sequences(pangenome: Pangenome, file_obj: IO, add: str = ""):
    """
    Export the sequence of gene families

    :param pangenome: Pangenome containing families
    :param file_obj: Temporary file where sequences will be written
    :param add: Add prefix to sequence name
    """
    for fam in pangenome.gene_families:
        file_obj.write(">" + add + fam.name + "\n")
        file_obj.write(fam.sequence + "\n")
    # file_obj.flush()

def write_all_gene_sequences(pangenome: Pangenome, file_obj: IO, add: str = "", disable_bar:bool = False):
    """
    Export the sequence of pangenome genes

    :param pangenome: Pangenome containing genes
    :param file_obj: Temporary file where sequences will be written
    :param add: Add prefix to sequence name
    """
    gene_ids_to_write = {gene.ID for fam in pangenome.gene_families for gene in fam.genes}
    get_gene_sequences_from_file(pangenome.file, file_obj, gene_ids_to_write,
                                         disable_bar=disable_bar)
    
    # for gene in pangenome.genes:
    #     file_obj.write(">" + add + gene.ID + "\n")
    #     file_obj.write(gene.protein + "\n")
    # file_obj.flush()

def project_and_write_partition(seqid_to_gene_family: Dict[str, GeneFamily], seq_set: Set[str], output: Path) -> Path:
    """
    Project the partition of each sequence from the input file and write them in a file

    :param seqid_to_gene_family: dictionnary which link sequence and pangenome gene family
    :param seq_set: input sequences
    :param output: Path of the output directory

    :return: Path to file which contain partition projection
    """

    partition_proj = output.absolute() / "sequences_partition_projection.tsv"
    with open(partition_proj, "w") as partProjFile:
        for input_seq, pangFam in seqid_to_gene_family.items():
            partProjFile.write(input_seq + "\t" + pangFam.named_partition + "\n")
        for remainingSeq in  seq_set - seqid_to_gene_family.keys():
            partProjFile.write(remainingSeq + "\tcloud\n")  # if there is no hit, it's going to be cloud genes.
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
        for fam in [gene.family for border in rgp.get_bordering_genes(pangenome.parameters["spot"]["set_size"],
                                                                      multigenics) for gene in border]:
            fam2rgp[fam].append(rgp.name)
    return fam2rgp


def get_fam_to_spot(pangenome: Pangenome, multigenics: Set[GeneFamily]) \
        -> Tuple[Dict[str, List[Spot]], Dict[str, List[Spot]]]:
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
                                rgp.get_bordering_genes(pangenome.parameters["spot"]["set_size"], multigenics)
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
        fname = output / f"spot_{str(spot.ID)}.gexf"
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
    logging.getLogger("PPanGGOLiN").info("Writing RGP and spot information related to hits in the pangenome")
    multigenics = pangenome.get_multigenics(pangenome.parameters["rgp"]["dup_margin"])

    finfo = open(output / "info_input_seq.tsv", "w")
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
        logging.getLogger("PPanGGOLiN").info(f"Drawing the {len(drawn_spots)} spots with more than 1 organization "
                                             f"related to hits of the input sequences...")
        draw_selected_spots(drawn_spots, pangenome, output, pangenome.parameters["spot"]["overlapping_match"],
                            pangenome.parameters["spot"]["exact_match_size"], pangenome.parameters["spot"]["set_size"],
                            disable_bar=disable_bar)

        fam2mod = {}  # fam2module
        if pangenome.status["modules"] != "No":
            for mod in pangenome.modules:
                for fam in mod.families:
                    fam2mod[fam] = f"module_{mod.ID}"

        draw_spot_gexf(drawn_spots, output, multigenics=multigenics, fam_to_mod=fam2mod)

    logging.getLogger("PPanGGOLiN").info(f"File listing RGP and spots where sequences of interest are located : "
                                         f"{output / 'info_input_seq.tsv'}")


def get_seq2pang(pangenome: Pangenome, sequence_file: Path, output: Path, tmpdir: Path,
                 cpu: int = 1, no_defrag: bool = False, identity: float = 0.8,
                 coverage: float = 0.8, is_nucleotide: bool = False, translation_table: int = 11,
                 target_type: str = "representative") -> Tuple[set, str, dict]:
    """
    Assign gene families from a pangenome to input sequences.

    This function aligns input sequences to gene families in a pangenome using MMseqs2 and assigns them
    to appropriate gene families based on alignment results.

    :param pangenome: Annotated pangenome containing gene families.
    :param sequence_file: Path to a FASTA file containing input sequences to align.
    :param output: Path to the output directory where alignment results will be stored.
    :param tmpdir: Temporary directory for intermediate files.
    :param cpu: Number of CPU cores to use for the alignment (default: 1).
    :param no_defrag: If True, the defragmentation workflow is skipped (default: False).
    :param identity: Minimum identity threshold for the alignment (default: 0.8).
    :param coverage: Minimum coverage threshold for the alignment (default: 0.8).
    :param is_nucleotide: Set to True if the sequence file contains nucleotide sequences to be translated.
                          If True, sequences will be translated using the specified translation table (default: False).
    :param translation_table: Translation table to use if sequences need to be translated (default: 11).
    :param target_type: Type of target sequences to align input sequences with, either 'all' or 'representative' 
                        (default: 'representative').

    :return: A tuple containing the set of input sequences, the path to the alignment result file, 
             and a dictionary mapping input sequences to gene families.
    :raises ValueError: If the `target_type` is not 'all' or 'representative'.

    """

    with tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.as_posix(), delete=False, suffix=".faa") as tmp_pang_file:
        if target_type == "representative":
            write_gene_fam_sequences(pangenome, tmp_pang_file, add="ppanggolin_")
        elif target_type == "all":
            write_all_gene_sequences(pangenome, tmp_pang_file, add="ppanggolin_")
        else:
            raise ValueError('Invalid value for target_type. It should be "all" or "representative".')

        with read_compressed_or_not(sequence_file) as seqFileObj:
            seq_set = get_seq(seqFileObj)
            align_file = align_seq_to_pang(tmp_pang_file, seqFileObj, output, tmpdir, cpu, 
                                           no_defrag, identity, coverage, is_nucleotide, translation_table )
        if target_type == "representative":
            seq2pang, align_file = associate_input_seq_to_gene_family_from_aln_rep(align_file, output, pangenome)
        else:
            seq2pang, align_file = associate_input_seq_to_gene_family_from_aln_all(align_file, output, pangenome)

    return seq_set, align_file, seq2pang


def align(pangenome: Pangenome, sequence_file: Path, output: Path, identity: float = 0.8,
          coverage: float = 0.8, no_defrag: bool = False, cpu: int = 1, getinfo: bool = False,
          draw_related: bool = False, tmpdir: Path = None, disable_bar: bool = False, keep_tmp=False):
    """
    Aligns pangenome sequences with sequences in a FASTA file using MMSeqs2.

    :param pangenome: Pangenome containing gene families to align with the input sequences.
    :param sequence_file: Path to a FASTA file containing sequences to align with the pangenome.
    :param output: Path to the output directory.
    :param identity: Minimum identity threshold for the alignment.
    :param coverage: Minimum coverage threshold for the alignment.
    :param no_defrag: If True, the defrag workflow will not be used.
    :param cpu: Number of CPU cores to use.
    :param getinfo: If True, extract information related to the best hit of each query, such as the RGP it is in or the spots.
    :param draw_related: If True, draw figures and graphs in a gexf format of spots associated with the input sequences.
    :param tmpdir: Temporary directory for intermediate files.
    :param disable_bar: If True, disable the progress bar.
    :param keep_tmp: If True, keep temporary files.
    """


    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
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

    if keep_tmp:
        
        dir_name = 'align_tmpdir' +  time.strftime("_%Y-%m-%d_%H.%M.%S",time.localtime())
        tmp_path = Path(tmpdir) / dir_name
        mk_outdir(tmp_path, force=True)
        logging.getLogger().info(f'Temporary files will be written {tmp_path} and kept for reference.')

    else:
        # if keep tmp is false, TemporaryDirectory created and then removed 
        new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
        tmp_path = Path(new_tmpdir.name)
        print(tmp_path)

    seq_set, align_file, seq2pang = get_seq2pang(pangenome, sequence_file, output, tmp_path, cpu, no_defrag, identity,
                                                 coverage)

    if getinfo or draw_related:  # TODO Add getinfo to function and remove if
        get_seq_info(seq2pang, pangenome, output, draw_related, disable_bar=disable_bar)

    part_proj = project_and_write_partition(seq2pang, seq_set, output)  # write the partition assignation only
    logging.getLogger().info(f"sequences partition projection : '{part_proj}'")
    logging.getLogger().info(f"{len(seq2pang)} sequences over {len(seq_set)} have at least one hit in the pangenome.")
    logging.getLogger().info(f"Blast-tab file of the alignment : '{align_file.name}'")

    # new_tmpdir.cleanup()


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
          draw_related=args.draw_related, disable_bar=args.disable_prog_bar, keep_tmp=args.keep_tmp)


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
    required.add_argument('-S', '--sequences', required=False, type=Path,
                          help="sequences (nucleotides or amino acids) to align on the pangenome gene families")

    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome .h5 file")
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
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    optional.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                        help="Keeping temporary files (useful for debugging).")

if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)
    parser_align(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
