#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tempfile
import time
import logging
import os
from typing import List, Dict, Tuple


# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.genome import Gene, Contig
from ppanggolin.utils import mk_outdir, restricted_float, add_gene, connected_components
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import get_seq2pang, project_partition
from ppanggolin.region import GeneContext


def search_gene_context_in_pangenome(pangenome: Pangenome, output: str, tmpdir: str, sequences: str = None,
                                     families: str = None, transitive: int = 4, identity: float = 0.5,
                                     coverage: float = 0.8, jaccard: float = 0.85, window_size: int = 1, no_defrag: bool = False,
                                     cpu: int = 1, disable_bar=True):
    """
    Main function to search common gene contexts between sequence set and pangenome families

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :param sequences: Path to file containing the sequences
    :param families: Path to file containing families name
    :param output: Path to output directory
    :param tmpdir: Path to temporary directory
    :param transitive: number of genes to check on both sides of a family aligned with an input sequence
    :param identity: minimum identity threshold between sequences and gene families for the alignment
    :param coverage: minimum coverage threshold between sequences and gene families for the alignment
    :param jaccard: Jaccard index to filter edges in graph
    :param window_size: Number of genes to consider in the gene context.
    :param no_defrag: do not use the defrag workflow if true
    :param cpu: Number of core used to process
    :param disable_bar: Allow preventing bar progress print
    """

    # check statuses and load info
    if sequences is not None and pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering has been made by PPanGGOLiN.")

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)
    gene_families = {}
    fam_2_seq = None
    if sequences is not None:
        # Alignment of sequences on pangenome families
        new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
        seq_set, _, seq2pan = get_seq2pang(pangenome, sequences, output, new_tmpdir, cpu, no_defrag, identity,
                                           coverage)
        project_partition(seq2pan, seq_set, output)
        
        new_tmpdir.cleanup()

        for pan_family in seq2pan.values():
            gene_families[pan_family.name] = pan_family
        
        fam_2_seq = fam2seq(seq2pan)

    if families is not None:
        with open(families, 'r') as f:
            for fam_name in f.read().splitlines():
                gene_families[fam_name] = pangenome.get_gene_family(fam_name)

    half_window = round((window_size-1)/2)
    logging.info(f'Window size of {half_window*2 + 1}. Gene context will include {half_window} genes on each side of the target gene.')

    # Compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger().info("Building the graph...")
    gene_context_graph = compute_gene_context_graph(families=gene_families, t=transitive, half_window=half_window, disable_bar=disable_bar)
    logging.getLogger().info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common gene contexts")
    logging.getLogger().debug(f"There are {nx.number_of_nodes(gene_context_graph)} nodes and {nx.number_of_edges(gene_context_graph)} edges")
    
    
    # extract the modules from the graph
    common_components = compute_gene_context(gene_context_graph, jaccard)

    families = set()
    for gene_context in common_components:
        families |= gene_context.families

    if len(families) != 0:
        export_to_dataframe(families, common_components, fam_2_seq, output)
    else:
        logging.getLogger().info(f"No gene contexts were found")

    logging.getLogger().info(f"Computing gene contexts took {round(time.time() - start_time, 2)} seconds")



    # for e, data in gene_context_graph(data=True):

    # nx.write_graphml_lxml(gene_context_graph, os.path.join(output, "context.graphml"))
    

def compute_gene_context_graph(families: dict, t: int = 4, half_window: int = 0, disable_bar: bool = False) -> nx.Graph:
    """
    Construct the graph of gene contexts between families of the pan

    :param families: Gene families of interest
    :param t: transitive value
    :param half_window: An integer specifying the number of genes to include in the context on each side of the gene of interest.
    :param disable_bar: Prevents progress bar printing

    :return: Graph of gene contexts between interesting gene families of the pan
    """

    g = nx.Graph()
    for family in tqdm(families.values(), unit="families", disable=disable_bar):
        for gene in family.genes:
            contig = gene.contig.genes
            pos_left, in_context_left, pos_right, in_context_right = extract_gene_context(gene, contig, families, t, half_window)
            if in_context_left or in_context_right:
                for env_gene in contig[pos_left:pos_right + 1]:
                    _compute_gene_context_graph(g, env_gene, contig, pos_right)
    return g


def _compute_gene_context_graph(g: nx.Graph, env_gene: Gene, contig: Contig, pos_r: int):
    """
    Compute graph of gene contexts between one gene and the other part of the contig

    :param: Graph of gene contexts between interesting gene families of the pan
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



def extract_gene_context(gene: Gene, contig: List[Gene], families: Dict[str, str], t: int = 4, half_window: int = 0) -> Tuple[int, bool, int, bool]:
    """
    Determine the left and rigth position of the gene context and whether said gene context exists. 

    :param gene: Gene of interest
    :param contig: list of genes in contig
    :param families: Alignment results
    :param t: transitive value
    :param half_window: An integer specifying the number of genes to include in the context on each side of the gene of interest.

    :return: Position of the context and if it exists for each side ('left' and 'right')
    """

    search_window = max(t, half_window)

    pos_left, pos_right = (max(0, gene.position - search_window),
                           min(gene.position + search_window, len(contig) - 1))  # Gene positions to compare family
    
    in_context_left, in_context_right = (False, False)
    while pos_left < gene.position and not in_context_left:
        if gene.position - pos_left <= half_window:
            # position is in the window 
            in_context_left = True

        elif contig[pos_left].family in families.values():
            in_context_left = True
        else:
            pos_left += 1

    while pos_right > gene.position and not in_context_right:
        if pos_right - gene.position <= half_window:
            in_context_right = True
        elif contig[pos_right].family in families.values():
            in_context_right = True
        else:
            pos_right -= 1

    return pos_left, in_context_left, pos_right, in_context_right


def compute_gene_context(g: nx.Graph, jaccard: float = 0.85) -> set:
    """
    Compute the gene contexts in the graph

    :param g: Graph of gene contexts between interesting gene families of the pan
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
    :param gene_contexts: connected components found in the pan
    :param fam_to_seq: Dictionary with gene families as keys and list of sequence ids as values
    :param output: output path
    """

    logging.getLogger().debug(f"There are {len(families)} families among {len(gene_contexts)} gene contexts")

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
    logging.getLogger(f"detected gene context(s) are listed in: '{output}/gene_contexts.tsv'")


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
                                     sequences=args.sequences, families=args.family, transitive=args.transitive,
                                     identity=args.identity, coverage=args.coverage, jaccard=args.jaccard, window_size=args.window_size,
                                     no_defrag=args.no_defrag, cpu=args.cpu, disable_bar=args.disable_prog_bar)


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
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome.h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('-S', '--sequences', required=False, type=str,
                        help="Fasta file with the sequences of interest")
    onereq.add_argument('-F', '--family', required=False, type=str,
                        help="List of family IDs of interest from the pan")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family.")
    optional.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("-t", "--transitive", required=False, type=int, default=4,
                          help="Size of the transitive closure used to build the graph. This indicates the number of "
                               "non related genes allowed in-between two related genes. Increasing it will improve "
                               "precision but lower sensitivity a little.")
    optional.add_argument("-w", "--window_size", required=False, type=int, default=1,
                        help="Number of genes adjacent to a gene of interest to consider in the gene context even if they are non related genes.")
    
    optional.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                          help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
                               "will improve precision but lower sensitivity a lot.")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_context(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
                        help="directory for storing temporary files")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
