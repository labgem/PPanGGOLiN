#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import tempfile
import subprocess
from collections import defaultdict
import os
import argparse
from typing import TextIO, Tuple, Dict, Set
from pathlib import Path

# installed libraries
from networkx import Graph
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.utils import read_compressed_or_not, restricted_float
from ppanggolin.formats.writeBinaries import write_pangenome, erase_pangenome
from ppanggolin.formats.readBinaries import check_pangenome_info, get_gene_sequences_from_file
from ppanggolin.formats.writeSequences import write_gene_sequences_from_annotations


# Global functions
def check_pangenome_former_clustering(pangenome: Pangenome, force: bool = False):
    """
    Checks pangenome status and .h5 files for former clusterings, delete them if allowed or raise an error

    :param pangenome: Annotated Pangenome
    :param force: Force to write on existing pangenome information
    """
    if pangenome.status["genesClustered"] == "inFile" and not force:
        raise Exception("You are trying to cluster genes that are already clustered together. If you REALLY want to "
                        "do that, use --force (it will erase everything except annotation data in your HDF5 file!)")
    elif pangenome.status["genesClustered"] == "inFile" and force:
        erase_pangenome(pangenome, gene_families=True)


# Clustering functions
def check_pangenome_for_clustering(pangenome: Pangenome, tmp_file: TextIO, force: bool = False,
                                   disable_bar: bool = False):
    """
    Check the pangenome statuses and write the gene sequences in the provided tmpFile.
    (whether they are written in the .h5 file or currently in memory)

    :param pangenome: Annotated Pangenome
    :param tmp_file: Temporary file
    :param force: Force to write on existing pangenome information
    :param disable_bar: Allow to disable progress bar
    """
    check_pangenome_former_clustering(pangenome, force)
    if pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
        # we append the gene ids by 'ppanggolin' to avoid crashes from mmseqs when sequence IDs are only numeric.
        write_gene_sequences_from_annotations(pangenome, tmp_file, add="ppanggolin_", disable_bar=disable_bar)
    elif pangenome.status["geneSequences"] == "inFile":
        get_gene_sequences_from_file(pangenome.file, tmp_file, add="ppanggolin_",
                                     disable_bar=disable_bar)  # write CDS sequences to the tmpFile
    else:
        tmp_file.close()  # closing the tmp file since an exception will be raised.
        raise Exception("The pangenome does not include gene sequences, thus it is impossible to cluster "
                        "the genes in gene families. Either provide clustering results (see --clusters), "
                        "or provide a way to access the gene sequence during the annotation step "
                        "(having the fasta in the gff files, or providing the fasta files through the --fasta option)")


def first_clustering(sequences: TextIO, tmpdir: Path, cpu: int = 1, code: int = 11, coverage: float = 0.8,
                     identity: float = 0.8, mode: int = 1) -> Tuple[Path, Path]:
    """
    Make a first clustering of all sequences in pangenome

    :param sequences: Sequence from pangenome
    :param tmpdir: Temporary directory
    :param cpu: number of CPU cores to use
    :param code: Genetic code used
    :param coverage: minimal coverage threshold for the alignment
    :param identity: minimal identity threshold for the alignment
    :param mode: MMseqs2 clustering mode

    :return: path to representative sequence file and path to tsv clustering result
    """
    seq_nucdb = tmpdir/'nucleotid_sequences_db'
    cmd = list(map(str, ["mmseqs", "createdb", sequences.name, seq_nucdb]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    logging.getLogger("PPanGGOLiN").info("Creating sequence database...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger("PPanGGOLiN").debug("Translate sequence ...")
    seqdb = tmpdir/'aa_db'
    cmd = list(map(str, ["mmseqs", "translatenucs", seq_nucdb, seqdb, "--threads", cpu, "--translation-table", code]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger("PPanGGOLiN").info("Clustering sequences...")
    cludb = tmpdir/'cluster_db'
    cmd = list(map(str, ["mmseqs", "cluster", seqdb, cludb, tmpdir, "--cluster-mode", mode, "--min-seq-id",
                         identity, "-c", coverage, "--threads", cpu, "--kmer-per-seq", 80, "--max-seqs", 300]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger("PPanGGOLiN").info("Extracting cluster representatives...")
    repdb = tmpdir/'representative_db'
    cmd = list(map(str, ["mmseqs", "result2repseq", seqdb, cludb, repdb]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    reprfa = tmpdir/'representative_sequences.fasta'
    cmd = list(map(str, ["mmseqs", "result2flat", seqdb, seqdb, repdb, reprfa, "--use-fasta-header"]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger("PPanGGOLiN").info("Writing gene to family informations")
    outtsv = tmpdir/'families_tsv'
    cmd = list(map(str, ["mmseqs", "createtsv", seqdb, seqdb, cludb, outtsv, "--threads", cpu, "--full-header"]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    return reprfa, outtsv


def read_faa(faa_file_name: Path) -> Dict[str, str]:
    """
    Read a faa file to link pangenome families to sequences.

    :param faa_file_name: path to the faa file

    :return: dictionary with families ID as key and sequence as value
    """
    fam2seq = {}
    head = ""
    with open(faa_file_name, "r") as faaFile:
        for line in faaFile:
            if line.startswith('>'):
                head = line[1:].strip().replace("ppanggolin_", "")  # remove the eventual addition
            else:
                fam2seq[head] = line.strip()
    return fam2seq


def align_rep(faa_file: Path, tmpdir: Path, cpu: int = 1, coverage: float = 0.8, identity: float = 0.8) -> Path:
    """
    Align representative sequence

    :param faa_file: sequence of representative family
    :param tmpdir: Temporary directory
    :param cpu: number of CPU cores to use
    :param coverage: minimal coverage threshold for the alignment
    :param identity: minimal identity threshold for the alignment

    :return: Result of alignment
    """
    logging.getLogger("PPanGGOLiN").debug("Create database")
    seqdb = tmpdir/'rep_sequence_db'
    cmd = list(map(str, ["mmseqs", "createdb", faa_file, seqdb]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger("PPanGGOLiN").info("Aligning cluster representatives...")
    alndb = tmpdir/'rep_alignment_db'
    cmd = list(map(str, ["mmseqs", "search", seqdb, seqdb, alndb, tmpdir, "-a", "--min-seq-id", identity,
                         "-c", coverage, "--cov-mode", 1, "--threads", cpu]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    logging.getLogger("PPanGGOLiN").info("Extracting alignments...")
    outfile = tmpdir/'rep_families.tsv'
    cmd = list(map(str, ["mmseqs", "convertalis", seqdb, seqdb, alndb, outfile,
                         "--format-output", "query,target,qlen,tlen,bits"]))
    logging.getLogger("PPanGGOLiN").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    return outfile


def read_tsv(tsv_file_name: Path) -> Tuple[Dict[str, Tuple[str, bool]], Dict[str, Set[str]]]:
    """Reading tsv file

    :param tsv_file_name: path to the tsv

    :return: two dictionnary which link genes and families
    """
    genes2fam = {}
    fam2genes = defaultdict(set)
    with open(tsv_file_name, "r") as tsvfile:
        for line in tsvfile:
            line = line.replace('"', '').replace("ppanggolin_", "").split()
            # remove the '"' char which protects the fields, and the eventual addition
            genes2fam[line[1]] = (line[0], False)  # fam id, and it's a gene (and not a fragment)
            fam2genes[line[0]].add(line[1])
    return genes2fam, fam2genes


def refine_clustering(tsv: str, aln_file: str, fam_to_seq: dict) -> Tuple[Dict[str, Tuple[str, bool]], Dict[str, str]]:
    """
    Refine clustering by removing fragment

    :param tsv: First clusterin result
    :param aln_file: Reprensentative alignment result
    :param fam_to_seq: Dictionary which link families to sequence

    :return: Two dictionary which link genes and families
    """
    simgraph = Graph()
    genes2fam, fam2genes = read_tsv(tsv)
    logging.getLogger("PPanGGOLiN").info(f"Starting with {len(fam_to_seq)} families")
    # create the nodes
    for fam, genes in fam2genes.items():
        simgraph.add_node(fam, nbgenes=len(genes))

    # add the edges
    with open(aln_file, "r") as alnfile:
        for line in alnfile:
            line = line.replace('"', '').replace("ppanggolin_", "").split()  # remove the eventual addition

            if line[0] != line[1]:
                simgraph.add_edge(line[0], line[1], score=float(line[4]))
                simgraph.nodes[line[0]]["length"] = int(line[2])
                simgraph.nodes[line[1]]["length"] = int(line[3])

    for node, nodedata in simgraph.nodes(data=True):
        choice = (None, 0, 0, 0)
        for neighbor in simgraph.neighbors(node):
            nei = simgraph.nodes[neighbor]
            score = simgraph[neighbor][node]["score"]
            if nei["length"] > nodedata["length"] and nei["nbgenes"] >= nodedata["nbgenes"] and choice[3] < score:
                choice = (genes2fam[neighbor][0], nei["length"], nei["nbgenes"], score)
                # `genes2fam[neighbor]` instead of just neighbor in case that family has been assigned already
                # (this is for smaller fragments that are closer to other fragments than the actual gene family)

        if choice[0] is not None:
            genestochange = fam2genes[node]
            for gene in genestochange:
                genes2fam[gene] = (str(choice[0]), True)
                fam2genes[choice[0]].add(gene)
            del fam2genes[node]
    new_fam_to_seq = {}
    for fam in fam2genes:
        new_fam_to_seq[fam] = fam_to_seq[fam]
    logging.getLogger("PPanGGOLiN").info(f"Ending with {len(new_fam_to_seq)} gene families")
    return genes2fam, new_fam_to_seq


def read_fam2seq(pangenome: Pangenome, fam_to_seq: Dict[str, str]):
    """
    Add gene family to pangenome and sequences to gene families

    :param pangenome: Annotated pangenome
    :param fam_to_seq: Dictionary which link families and sequences
    """
    logging.getLogger("PPanGGOLiN").info("Adding protein sequences to the gene families")
    for family, protein in fam_to_seq.items():
        fam = GeneFamily(pangenome.max_fam_id, family)
        fam.add_sequence(protein)
        pangenome.add_gene_family(fam)


def read_gene2fam(pangenome: Pangenome, gene_to_fam: dict, disable_bar: bool = False):
    """
    Add gene to pangenome families

    :param pangenome: Annotated Pangenome
    :param gene_to_fam: Dictionary which link gene to families
    :param disable_bar: Allow to disable progress bar
    """
    logging.getLogger("PPanGGOLiN").info(f"Adding {len(gene_to_fam)} genes to the gene families")

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False
    if link and len(gene_to_fam) != pangenome.number_of_genes:  # then maybe there are genes with identical IDs
        raise Exception("Something unexpected happened during clustering (have less genes clustered than genes "
                        "in the pangenome). A probable reason is that two genes in two different organisms have "
                        "the same IDs; If you are sure that all of your genes have non identical IDs,  please post an "
                        "issue at https://github.com/labgem/PPanGGOLiN/")
    for gene, (family, is_frag) in tqdm(gene_to_fam.items(), unit="gene", total=len(gene_to_fam), disable=disable_bar):
        try:
            fam = pangenome.get_gene_family(family)
        except KeyError:  # Family not found so create and add
            fam = GeneFamily(pangenome.max_fam_id, family)
            pangenome.add_gene_family(fam)
        if link:  # doing the linking if the annotations are loaded.
            gene_obj = pangenome.get_gene(gene)
        else:
            gene_obj = Gene(gene)
        gene_obj.is_fragment = is_frag
        fam.add(gene_obj)


def clustering(pangenome: Pangenome, tmpdir: Path, cpu: int = 1, defrag: bool = True, code: int = 11,
               coverage: float = 0.8, identity: float = 0.8, mode: int = 1, force: bool = False,
               disable_bar: bool = False):
    """
    Main function to cluster pangenome gene sequences into families

    :param pangenome: Annoatated Pangenome
    :param tmpdir: Path to temporary directory
    :param cpu: number of CPU cores to use
    :param defrag: Allow to remove fragment
    :param code: Genetic code used
    :param coverage: minimal coverage threshold for the alignment
    :param identity: minimal identity threshold for the alignment
    :param mode: MMseqs2 clustering mode
    :param force: force to write in the pangenome
    :param disable_bar: Allow to disable progress bar
    """

    newtmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
    tmp_path = Path(newtmpdir.name)
    with open(tmp_path/'nucleotid_sequences', "w") as sequence_file:
        check_pangenome_for_clustering(pangenome, sequence_file, force, disable_bar=disable_bar)
        logging.getLogger("PPanGGOLiN").info("Clustering all of the genes sequences...")
        rep, tsv = first_clustering(sequence_file, tmp_path, cpu, code, coverage, identity, mode)

    fam2seq = read_faa(rep)
    if not defrag:
        logging.getLogger("PPanGGOLiN").debug("No defragmentation")
        genes2fam, _ = read_tsv(tsv)
    else:
        logging.getLogger("PPanGGOLiN").info("Associating fragments to their original gene family...")
        aln = align_rep(rep, tmp_path, cpu, coverage, identity)
        genes2fam, fam2seq = refine_clustering(tsv, aln, fam2seq)
        pangenome.status["defragmented"] = "Computed"
    newtmpdir.cleanup()
    read_fam2seq(pangenome, fam2seq)
    read_gene2fam(pangenome, genes2fam, disable_bar=disable_bar)

    pangenome.status["genesClustered"] = "Computed"
    pangenome.status["geneFamilySequences"] = "Computed"

    pangenome.parameters["cluster"] = {}
    pangenome.parameters["cluster"]["coverage"] = coverage
    pangenome.parameters["cluster"]["identity"] = identity
    pangenome.parameters["cluster"]["defragmentation"] = defrag
    pangenome.parameters["cluster"]["translation_table"] = code
    pangenome.parameters["cluster"]["read_clustering_from_file"] = False


# Read clustering
def mk_local_to_gene(pangenome: Pangenome) -> dict:
    """Creates a dictionary that stores local identifiers, if all local identifiers are unique (and if they exist)

    :param pangenome: Input Pangenome

    :return: Dictionary with local identifiers
    """
    local_dict = {}
    for gene in pangenome.genes:
        old_len = len(local_dict)
        local_dict[gene.local_identifier] = gene
        if len(local_dict) == old_len:
            if pangenome.parameters["annotation"]["read_annotations_from_file"] and not \
                    pangenome.parameters["annotation"]["used_local_identifiers"]:
                raise Exception(f"'{gene.local_identifier}' was found multiple times used as an identifier. "
                                f"The identifier of the genes (locus_tag, protein_id in gbff, ID in gff) were not "
                                f"unique throughout all of the files. It is thus impossible to differentiate the genes."
                                f" To use this function while importing annotation, all identifiers MUST be unique "
                                f"throughout all of your genomes")
            return {}  # local identifiers are not unique.
    return local_dict


def infer_singletons(pangenome: Pangenome):
    """Creates a new family for each gene with no associated family

    :param pangenome: Input pangenome
    """
    singleton_counter = 0
    for gene in pangenome.genes:
        if gene.family is None:
            fam = GeneFamily(family_id=pangenome.max_fam_id, name=gene.ID)
            fam.add(gene)
            pangenome.add_gene_family(fam)
            singleton_counter += 1
    logging.getLogger("PPanGGOLiN").info(f"Inferred {singleton_counter} singleton families")


def read_clustering(pangenome: Pangenome, families_tsv_file: Path, infer_singleton: bool = False, force: bool = False,
                    disable_bar: bool = False):
    """
    Get the pangenome information, the gene families and the genes with an associated gene family.
    Reads a families tsv file from mmseqs2 output and adds the gene families and the genes to the pangenome.

    :param pangenome: Input Pangenome
    :param families_tsv_file: MMseqs2 clustering results
    :param infer_singleton: creates a new family for each gene with no associated family
    :param force: force to write in the pangenome
    :param disable_bar: Allow to disable progress bar
    """
    check_pangenome_former_clustering(pangenome, force)
    check_pangenome_info(pangenome, need_annotations=True, disable_bar=disable_bar)

    logging.getLogger("PPanGGOLiN").info(f"Reading {families_tsv_file.name} the gene families file ...")
    filesize = os.stat(families_tsv_file).st_size
    families_tsv_file = read_compressed_or_not(families_tsv_file)
    frag = False  # the genome annotations are necessarily loaded.
    nb_gene_with_fam = 0
    local_dict = mk_local_to_gene(pangenome)
    bar = tqdm(total=filesize, unit="bytes", disable=disable_bar)
    line_counter = 0
    for line in families_tsv_file:
        line_counter += 1
        bar.update(len(line))
        try:
            elements = [el.strip() for el in line.split()]  # 2 or 3 fields expected
            if len(elements) <= 1:
                raise ValueError("No tabulation separator found in gene families file")
            (fam_id, gene_id, is_frag) = elements if len(elements) == 3 else elements + ["Na"]  # case of 2 fields
            try:
                gene_obj = pangenome.get_gene(gene_id)
            except KeyError:
                gene_obj = local_dict.get(gene_id)
            if gene_obj is not None:
                nb_gene_with_fam += 1
                try:
                    fam = pangenome.get_gene_family(fam_id)
                except KeyError:  # Family not found so create and add
                    fam = GeneFamily(pangenome.max_fam_id, fam_id)
                    pangenome.add_gene_family(fam)
                gene_obj.is_fragment = True if is_frag == "F" else False  # F for Fragment
                fam.add(gene_obj)
            if is_frag == "F":
                frag = True
        except Exception:
            raise Exception(f"line {line_counter} of the file '{families_tsv_file.name}' raised an error.")
    bar.close()
    families_tsv_file.close()
    if nb_gene_with_fam < pangenome.number_of_genes:  # not all genes have an associated cluster
        if nb_gene_with_fam == 0:
            raise Exception("No gene ID in the cluster file matched any gene ID from the annotation step."
                            " Please ensure that the annotations that you loaded previously and the clustering results "
                            "that you have used the same gene IDs. If you use .gff files it is the identifier stored in"
                            " the field 'ID'. If you use .gbff files it is the identifier stored in 'locus_tag'.")
        else:
            if infer_singleton:
                infer_singletons(pangenome)
            else:
                raise Exception(f"Some genes ({pangenome.number_of_genes - nb_gene_with_fam}) did not have an associated "
                                f"cluster. Either change your cluster file so that each gene has a cluster, "
                                f"or use the --infer_singletons option to infer a cluster for each non-clustered gene.")
    pangenome.status["genesClustered"] = "Computed"
    if frag:  # if there was fragment information in the file.
        pangenome.status["defragmented"] = "Computed"
    pangenome.parameters["cluster"] = {}
    pangenome.parameters["cluster"]["read_clustering_from_file"] = True
    pangenome.parameters["cluster"]["infer_singletons"] = infer_singleton


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.clusters is None:
        if args.infer_singletons is True:
            logging.getLogger("PPanGGOLiN").warning("--infer_singletons option is not compatible with clustering "
                                                    "creation. To infer singleton you should give a clustering")
        clustering(pangenome, args.tmpdir, args.cpu, defrag=not args.no_defrag, code=args.translation_table,
                   coverage=args.coverage, identity=args.identity, mode=args.mode, force=args.force,
                   disable_bar=args.disable_prog_bar)
        logging.getLogger("PPanGGOLiN").info("Done with the clustering")
    else:
        if None in [args.tmpdir, args.cpu, args.no_defrag, args.translation_table,
                    args.coverage, args.identity, args.mode]:
            logging.getLogger("PPanGGOLiN").warning("You are using an option compatible only with clustering creation.")
        read_clustering(pangenome, args.clusters, args.infer_singletons, args.force, disable_bar=args.disable_prog_bar)
        logging.getLogger("PPanGGOLiN").info("Done reading the cluster file")
    write_pangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("cluster", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_clust(parser)
    return parser


def parser_clust(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of cluster command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome .h5 file")
    clust = parser.add_argument_group(title="Clustering arguments")
    clust.add_argument("--identity", required=False, type=restricted_float, default=0.8,
                       help="Minimal identity percent for two proteins to be in the same cluster")
    clust.add_argument("--coverage", required=False, type=restricted_float, default=0.8,
                       help="Minimal coverage of the alignment for two proteins to be in the same cluster")
    clust.add_argument("--mode", required=False, default="1", choices=["0", "1", "2", "3"],
                       help="the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component),"
                            " 2: CD-HIT-like, 3: CD-HIT-like (lowmem)")
    clust.add_argument('--no_defrag', required=False, default=False, action="store_true",
                       help="DO NOT Use the defragmentation strategy to link potential fragments "
                            "with their original gene family.")
    clust.add_argument("--translation_table", required=False, default="11",
                       help="Translation table (genetic code) to use.")

    clust.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")

    read = parser.add_argument_group(title="Read clustering arguments")
    read.add_argument('--clusters', required=False, type=Path,
                      help="A tab-separated list containing the result of a clustering. One line per gene. "
                           "First column is cluster ID, and second is gene ID")
    read.add_argument("--infer_singletons", required=False, action="store_true",
                      help="When reading a clustering result with --clusters, if a gene is not in the provided file"
                           " it will be placed in a cluster where the gene is the only member.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)
    parser_clust(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
