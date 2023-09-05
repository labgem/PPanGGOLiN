#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import tempfile
import time
from pathlib import Path

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.genome import Gene, Contig
from ppanggolin.utils import mk_outdir, restricted_float, add_gene, connected_components, create_tmpdir
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import  project_and_write_partition, get_input_seq_to_family_with_rep, get_input_seq_to_family_with_all
from ppanggolin.region import GeneContext


def search_gene_context_in_pangenome(pangenome: Pangenome, output: Path, tmpdir: Path, sequence_file: Path = None,
                                     families: Path = None, transitive: int = 4, identity: float = 0.5,
                                     coverage: float = 0.8, use_representatives: bool = False, jaccard: float = 0.85, no_defrag: bool = False,
                                     cpu: int = 1, disable_bar=True, translation_table:int=11, keep_tmp:bool = False):
    """
    Main function to search common gene contexts between sequence set and pangenome families

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :param sequence_file: Path to file containing the sequences
    :param families: Path to file containing families name
    :param output: Path to output directory
    :param tmpdir: Path to temporary directory
    :param transitive: number of genes to check on both sides of a family aligned with an input sequence
    :param identity: minimum identity threshold between sequences and gene families for the alignment
    :param coverage: minimum coverage threshold between sequences and gene families for the alignment
    :param use_representatives: Use representative sequences of gene families rather than all sequences to align input genes
    :param jaccard: Jaccard index to filter edges in graph
    :param no_defrag: do not use the defrag workflow if true
    :param cpu: Number of core used to process
    :param disable_bar: Allow preventing bar progress print
    :param translation_table: Translation table ID for nucleotide sequences.
    :param keep_tmp: If True, keep temporary files.

    """
    # check statuses and load info
    if sequence_file is not None and pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)
    gene_families = {}
    fam_2_seq = None
    if sequence_file is not None:
        # Alignment of sequences on pangenome families
        with create_tmpdir(main_dir=tmpdir, basename="align_input_seq_tmpdir", keep_tmp=keep_tmp) as new_tmpdir:
        
            if use_representatives:
                seq_set, _, seq2pan = get_input_seq_to_family_with_rep(pangenome, sequence_file, output, new_tmpdir,
                                                            cpu, no_defrag, identity=identity, coverage=coverage,
                                                            translation_table=translation_table)
            else:
                seq_set, _, seq2pan = get_input_seq_to_family_with_all(pangenome=pangenome, sequence_file=sequence_file, 
                                                                                    output=output, tmpdir=new_tmpdir,
                                                                                    cpu=cpu, no_defrag=no_defrag,
                                                                                    identity=identity, coverage=coverage,
                                                                                    translation_table=translation_table)
        
        project_and_write_partition(seq2pan, seq_set, output)

        for k, v in seq2pan.items():
            gene_families[v.name] = v
        fam_2_seq = fam2seq(seq2pan)

    if families is not None:
        with open(families, 'r') as f:
            for fam_name in f.read().splitlines():
                gene_families[fam_name] = pangenome.get_gene_family(fam_name)

    # Compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger("PPanGGOLiN").info("Building the graph...")
    g = compute_gene_context_graph(families=gene_families, t=transitive, disable_bar=disable_bar)
    logging.getLogger("PPanGGOLiN").info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common gene contexts")
    logging.getLogger("PPanGGOLiN").debug(f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges")

    # extract the modules from the graph
    common_components = compute_gene_context(g, jaccard)

    families = set()
    for gene_context in common_components:
        families |= gene_context.families

    if len(families) != 0:
        export_to_dataframe(families, common_components, fam_2_seq, output)
    else:
        logging.getLogger("PPanGGOLiN").info("No gene contexts were found")

    logging.getLogger("PPanGGOLiN").info(f"Computing gene contexts took {round(time.time() - start_time, 2)} seconds")


def compute_gene_context_graph(families: dict, t: int = 4, disable_bar: bool = False) -> nx.Graph:
    """
    Construct the graph of gene contexts between families of the pangenome

    :param families: Gene families of interest
    :param t: transitive value
    :param disable_bar: Prevents progress bar printing

    :return: Graph of gene contexts between interesting gene families of the pangenome
    """

    g = nx.Graph()
    for family in tqdm(families.values(), unit="families", disable=disable_bar):
        for gene in family.genes:
            contig = gene.contig.genes
            pos_left, in_context_left, pos_right, in_context_right = extract_gene_context(gene, contig, families, t)
            if in_context_left or in_context_right:
                for env_gene in contig[pos_left:pos_right + 1]:
                    _compute_gene_context_graph(g, env_gene, contig, pos_right)
    return g


def _compute_gene_context_graph(g: nx.Graph, env_gene: Gene, contig: Contig, pos_r: int):
    """
    Compute graph of gene contexts between one gene and the other part of the contig

    :param: Graph of gene contexts between interesting gene families of the pangenome
    :param env_gene: Gene of the current position
    :param contig: Current contig to search a gene context
    :param pos_r: Gene to search a gene context
    """

    g.add_node(env_gene.family)
    add_gene(g.nodes[env_gene.family], env_gene, fam_split=False)
    pos = env_gene.position + 1
    while pos <= pos_r:
        if env_gene.family != contig[pos].family:
            g.add_edge(env_gene.family, contig[pos].family)
            edge = g[env_gene.family][contig[pos].family]
            add_gene(edge, env_gene)
            add_gene(edge, contig[pos])
        pos += 1


def extract_gene_context(gene: Gene, contig: list, families: dict, t: int = 4) -> (int, bool, int, bool):
    """
    Extract gene context and whether said gene context exists

    :param gene: Gene of interest
    :param contig: list of genes in contig
    :param families: Alignment results
    :param t: transitive value

    :return: Position of the context and if it exists for each side ('left' and 'right')
    """

    pos_left, pos_right = (max(0, gene.position - t),
                           min(gene.position + t, len(contig) - 1))  # Gene positions to compare family
    in_context_left, in_context_right = (False, False)
    while pos_left < gene.position and not in_context_left:
        if contig[pos_left].family in families.values():
            in_context_left = True
        else:
            pos_left += 1

    while pos_right > gene.position and not in_context_right:
        if contig[pos_right].family in families.values():
            in_context_right = True
        else:
            pos_right -= 1

    return pos_left, in_context_left, pos_right, in_context_right


def compute_gene_context(g: nx.Graph, jaccard: float = 0.85) -> set:
    """
    Compute the gene contexts in the graph

    :param g: Graph of gene contexts between interesting gene families of the pangenome
    :param jaccard: Jaccard index

    :return: Set of gene contexts find in graph
    """

    gene_contexts = set()
    c = 1
    for comp in connected_components(g, removed=set(), weight=jaccard):
        gene_contexts.add(GeneContext(gc_id=c, families=comp))
        c += 1
    return gene_contexts


def fam2seq(seq_to_pan: dict) -> dict:
    """
    Create a dictionary with gene families as keys and list of sequences id as values

    :param seq_to_pan: Dictionary storing the sequence ids as keys and the gene families
                       to which they are assigned as values

    :return: Dictionary reversed
    """

    fam_2_seq = {}
    for sequence, family in seq_to_pan.items():
        if family.ID in fam_2_seq.keys():
            fam_2_seq[family.ID].append(sequence)
        else:
            fam_2_seq[family.ID] = [sequence]
    return fam_2_seq


def export_to_dataframe(families: set, gene_contexts: set, fam_to_seq: dict, output: str):
    """ Export the results into dataFrame

    :param families: Families related to the connected components
    :param gene_contexts: connected components found in the pangenome
    :param fam_to_seq: Dictionary with gene families as keys and list of sequence ids as values
    :param output: output path
    """

    logging.getLogger("PPanGGOLiN").debug(f"There are {len(families)} families among {len(gene_contexts)} gene contexts")

    lines = []
    for gene_context in gene_contexts:
        for family in gene_context.families:
            line = [gene_context.ID]
            if fam_to_seq is None or fam_to_seq.get(family.ID) is None:
                line += [family.name, None, len(family.organisms), family.named_partition]
            else:
                line += [family.name, ','.join(fam_to_seq.get(family.ID)),
                         len(family.organisms), family.named_partition]
            lines.append(line)
    df = pd.DataFrame(lines,
                      columns=["GeneContext ID", "Gene family name", "Sequence ID", "Nb Genomes", "Partition"]
                      ).set_index("GeneContext ID")
    df.sort_values(["GeneContext ID", "Sequence ID"], na_position='last').to_csv(
        path_or_buf=f"{output}/gene_contexts.tsv", sep="\t", na_rep='NA')
    logging.getLogger("PPanGGOLiN").info(f"detected gene context(s) are listed in: '{output}/gene_contexts.tsv'")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    if not any([args.sequences, args.family]):
        raise Exception("At least one of --sequences or --family option must be given")
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    search_gene_context_in_pangenome(pangenome=pangenome, output=args.output, tmpdir=args.tmpdir,
                                     sequence_file=args.sequences, families=args.family, transitive=args.transitive,
                                     identity=args.identity, coverage=args.coverage, 
                                     use_representatives=args.fast, jaccard=args.jaccard,
                                     no_defrag=args.no_defrag, cpu=args.cpu, disable_bar=args.disable_prog_bar,
                                     translation_table=args.translation_table, keep_tmp=args.keep_tmp)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """

    parser = sub_parser.add_parser("context", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_context(parser)
    return parser


def parser_context(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of context command

    :param parser: parser for align argument
    """

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome.h5 file")
    required.add_argument('-o', '--output', required=False, type=Path,
                          help="Output directory where the file(s) will be written")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('-S', '--sequences', required=False, type=Path,
                        help="Fasta file with the sequences of interest")
    onereq.add_argument('-F', '--family', required=False, type=Path,
                        help="List of family IDs of interest from the pangenome")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family.")
    optional.add_argument("--fast", required=False, action="store_true",
                            help="Use representative sequences of gene families for input gene alignment. "
                                "This option is recommended for faster processing but may be less sensitive. "
                                "By default, all pangenome genes are used for alignment. "
                                "This argument makes sense only when --sequence is provided.")
    optional.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("--translation_table", required=False, default="11",
                          help="The translation table (genetic code) to use when the input sequences are nucleotide sequences. ")
    optional.add_argument("-t", "--transitive", required=False, type=int, default=4,
                          help="Size of the transitive closure used to build the graph. This indicates the number of "
                               "non related genes allowed in-between two related genes. Increasing it will improve "
                               "precision but lower sensitivity a little.")
    optional.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                          help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
                               "will improve precision but lower sensitivity a lot.")
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    optional.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                        help="Keeping temporary files (useful for debugging).")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_context(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
