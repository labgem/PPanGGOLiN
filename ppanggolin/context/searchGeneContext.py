#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tempfile
import time
import logging
import os
from typing import List, Dict, Tuple, Iterable, Hashable, Iterator, Set
from itertools import zip_longest, chain
from collections import defaultdict

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.genome import Gene, Contig
from ppanggolin.utils import mk_outdir, restricted_float, connected_components
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import get_seq2pang, project_partition
from ppanggolin.region import GeneContext
from ppanggolin.geneFamily import GeneFamily


def search_gene_context_in_pangenome(pangenome: Pangenome, output: str, tmpdir: str, sequences: str = None,
                                     families: str = None, transitive: int = 4, identity: float = 0.5,
                                     coverage: float = 0.8, jaccard_threshold: float = 0.85, window_size: int = 1, no_defrag: bool = False,
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
    :param jaccard_threshold: Jaccard index threshold to filter edges in graph
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

    gene_families = set()
    fam_2_seq = None
    
    if sequences is not None:
        # Alignment of sequences on pangenome families
        new_tmpdir = tempfile.TemporaryDirectory(dir=tmpdir)
        seq_set, _, seq2pan = get_seq2pang(pangenome, sequences, output, new_tmpdir, cpu, no_defrag, identity,
                                           coverage)
        project_partition(seq2pan, seq_set, output)
        
        new_tmpdir.cleanup()

        for pan_family in seq2pan.values():
            gene_families.add(pan_family)
        
        fam_2_seq = fam2seq(seq2pan)

    if families is not None:
        with open(families, 'r') as f:
            for fam_name in f.read().splitlines():
                gene_families.add(pangenome.get_gene_family(fam_name))

    # half_window = round((window_size-1)/2)
    # logging.info(f'Window size of {half_window*2 + 1}. Gene context will include {half_window} genes on each side of the target gene.')

    # Compute the graph with transitive closure size provided as parameter
    start_time = time.time()

    logging.getLogger().info("Building the graph...")
    
    gene_context_graph = compute_gene_context_graph(families=gene_families, transitive=transitive, window_size=window_size, disable_bar=disable_bar)
    
    logging.getLogger().info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common gene contexts")
    
    logging.getLogger().debug(f"There are {nx.number_of_nodes(gene_context_graph)} nodes and {nx.number_of_edges(gene_context_graph)} edges")

    compute_edge_metrics(gene_context_graph, jaccard_threshold)


    write_graph(gene_context_graph, output, gene_families)
    
    # extract the modules from the graph
    # common_components = compute_gene_context(gene_context_graph, jaccard)

    # families = set()
    # for gene_context in common_components:
    #     families |= gene_context.families

    # if len(families) != 0:
    #     export_to_dataframe(families, common_components, fam_2_seq, output)
    # else:
    #     logging.getLogger().info(f"No gene contexts were found")

    logging.getLogger().info(f"Computing gene contexts took {round(time.time() - start_time, 2)} seconds")

def write_graph(context_graph:nx.Graph, output_dir:str, famillies_of_interest):

    def filter_edge_attribute(data):
        return {k:v for k, v in data.items() if type(v) != set}
        
    G = nx.Graph()

    G.add_edges_from(((f1.name, f2.name) for f1,f2 in context_graph.edges()))

    edges_with_attributes = {(f1.name, f2.name):filter_edge_attribute(d) for f1,f2,d in context_graph.edges(data=True)}
    
    nx.set_edge_attributes(G, edges_with_attributes)

    nodes_data = {f.name:{"organisms":len(f.organisms), "genes":len(f.genes), "famillies_of_interest": f in famillies_of_interest} for f in context_graph.nodes()}
    
    for f, d in G.nodes(data=True):
        d.update(nodes_data[f])

    nx.write_graphml_lxml(G, os.path.join(output_dir, "graph_context.graphml"))

def compute_edge_metrics(context_graph: nx.Graph, gene_proportion_cutoff:float ):
    # compute jaccard on organism
    for f1, f2, data in context_graph.edges(data=True):
        
        data['jaccard_organism'] = len(data['organisms'])/len(f1.organisms | f2.organisms)
        
        data['min_jaccard_organism'] = len(data['organisms'])/min(len(f1.organisms), len(f2.organisms))
        data['max_jaccard_organism'] = len(data['organisms'])/max(len(f1.organisms), len(f2.organisms))
    

        # print("f1", f1)
        # print("len data[f1]", len(data[f1]))
        
        # print('len(f1.genes)', len(f1.genes))
        
        f1_gene_proportion = len(data[f1])/len(f1.genes)
        f2_gene_proportion = len(data[f2])/len(f2.genes)
        
        data[f'f1'] = f1.name
        data[f'f2'] = f2.name
        data[f'f1_gene_proportion'] = f1_gene_proportion
        data[f'f2_gene_proportion'] = f2_gene_proportion
                        
        data[f'is_gene_proportion_>_{gene_proportion_cutoff}'] = (f1_gene_proportion >= gene_proportion_cutoff) and (f2_gene_proportion >= gene_proportion_cutoff)
        
        # for k, v in data.items():
        #     if type(v) != set:
        #         print(k, v)
        # print('===')  
                        


def add_edges_to_context_graph(context_graph: nx.Graph,
                               contig_genes: Iterable[Gene],
                               contig_windows: List[Tuple[int, int]],
                               t: int,
                               is_circular: bool):
    """
    Add edges to the context graph based on contig genes and windows.

    :param context_graph: The context graph to which edges will be added.
    :param contig_genes: An iterable of genes in the contig.
    :param contig_windows: A list of tuples representing the start and end positions of contig windows.
    :param t: The number of next genes to consider when adding edges.
    :param is_circular: A boolean indicating if the contig is circular.

    """
    for window_start, window_end in contig_windows:
        for gene_index in range(window_start, window_end + 1):
            gene = contig_genes[gene_index]
            next_genes = get_n_next_genes_index(gene_index, next_genes_count=t, 
                                                contig_size=len(contig_genes), is_circular=is_circular)
            next_genes = list(next_genes)

            for next_gene_index in next_genes:
                # Check if the next gene is within the contig windows
                if not any(lower <= next_gene_index <= upper for (lower, upper) in contig_windows):
                    # next_gene_index is not in any range of genes in the context
                    # so it is ignored along with all following genes
                    break
                
                next_gene = contig_genes[next_gene_index]
                if next_gene.family == gene.family:
                    # If the next gene has the same family, the two genes refer to the same node
                    # so they are ignored
                    continue
                
                context_graph.add_edge(gene.family, next_gene.family)

                # Add edge attributes
                edge_dict = context_graph[gene.family][next_gene.family]
                add_val_to_edge_attribute(edge_dict, gene.family, gene)
                add_val_to_edge_attribute(edge_dict, next_gene.family, next_gene)

                add_val_to_edge_attribute(edge_dict, "organisms", gene.organism)

                update_edge_attribute_counter(edge_dict, "gene_pairs")
                
                assert gene.organism == next_gene.organism 


def add_val_to_edge_attribute(edge_dict: dict, attribute_key, attribute_value):
    """
    Add an edge attribute value to the edge dictionary set.

    :param edge_dict: The dictionary containing the edge attributes.
    :param attribute_key: The key of the attribute.
    :param attribute_value: The value of the attribute to be added.

    """

    try:
        edge_dict[attribute_key].add(attribute_value)
    except KeyError:
        edge_dict[attribute_key] = {attribute_value}


def update_edge_attribute_counter(edge_dict: dict, key:Hashable):
    """
    Update the counter for an edge attribute in the edge dictionary.

    :param edge_dict: The dictionary containing the edge attributes.
    :param key: The key of the attribute.

    """

    try:
        edge_dict[key] += 1
    except KeyError:
        edge_dict[key] = 1


def get_n_next_genes_index(current_index: int, next_genes_count: int, contig_size: int, is_circular: bool = False) -> Iterator[int]:
    """
    Generate the indices of the next genes based on the current index and contig properties.

    :param current_index: The index of the current gene.
    :param next_genes_count: The number of next genes to consider.
    :param contig_size: The total number of genes in the contig.
    :param is_circular: Flag indicating whether the contig is circular (default: False).
    :return: An iterator yielding the indices of the next genes.

    Raises:
    - IndexError: If the current index is out of range for the given contig size.

    """

    # Check if the current index is out of range
    if current_index >= contig_size:
        raise IndexError(f'current gene index is out of range. '
                         f"Contig has {contig_size} genes while the given gene index is {current_index}")
    if is_circular:
        next_genes = chain(range(current_index+1, contig_size), range(0, current_index))
    else:
        next_genes = range(current_index+1, contig_size)
    
    for i, next_gene_index in enumerate(next_genes):
        if i == next_genes_count:
            break
        yield next_gene_index
        
def extract_contig_window(contig_size: int, positions_of_interest: Iterable[int], window_size: int, is_circular:bool = False):
    """
    Extracts contiguous windows around positions of interest within a contig.

    :param contig_size: Number of genes in contig.
    :param positions_of_interest: An iterable containing the positions of interest.
    :param window_size: The size of the window to extract around each position of interest.
    :param is_circular: Indicates if the contig is circular.
    :return: Yields tuples representing the start and end positions of each contiguous window.
    """
    windows_coordinates = []

    # Sort the positions of interest
    sorted_positions = sorted(positions_of_interest)

    # Check if any position of interest is out of range
    if sorted_positions[0] <0 or sorted_positions[-1] >= contig_size:
        raise IndexError(f'Positions of interest are out of range. '
                         f"Contig has {contig_size} genes while given min={sorted_positions[0]} & max={sorted_positions[-1]} positions")

    if is_circular:
        first_position = sorted_positions[0]
        last_position = sorted_positions[-1]
        # in a circular contig, if the window of a gene of interest overlaps the end/start of the contig
        # an out of scope position is added to the sorted positions to take into account those positions
        # the returned window are always checked that its positions are not out of range... 
        # so there's no chance to find an out of scope position in final list
        if first_position - window_size < 0:
            out_of_scope_position = (contig_size ) + first_position
            sorted_positions.append(out_of_scope_position)
    
        if last_position + window_size >= contig_size :
            out_of_scope_position = last_position - contig_size
            sorted_positions.insert(0, out_of_scope_position)
            
    start_po = max(sorted_positions[0] - window_size, 0)
    
    for position, next_po in zip_longest(sorted_positions, sorted_positions[1:]):
        
        if next_po is None:
            # If there are no more positions, add the final window
            end_po = min(position + window_size, contig_size-1)
            windows_coordinates.append((start_po, end_po))
            
        elif position + window_size +1 < next_po - window_size:
            # If there is a gap between positions, add the current window 
            # and update the start position for the next window
            end_po = min(position + window_size, contig_size-1)
            
            windows_coordinates.append((start_po, end_po))
            
            start_po = max(next_po - window_size, 0)
            
    return windows_coordinates

def get_contig_to_genes(gene_families: Iterable[GeneFamily]) -> Dict[Contig, Set[Gene]]:
    """
    Group genes from specified gene families by contig.

    :param gene_families: An iterable of gene families object.
    
    :return: A dictionary mapping contigs to sets of genes.
    """
    
    contig_to_genes_of_interest = defaultdict(set)
    for gene_family in gene_families:
        for gene in gene_family.genes:
            contig = gene.contig
            contig_to_genes_of_interest[contig].add(gene)
    return contig_to_genes_of_interest


def compute_gene_context_graph(families: Iterable[GeneFamily], transitive: int = 4, window_size: int = 0, disable_bar: bool = False) -> nx.Graph:
    """
    Construct the graph of gene contexts between families of the pangenome.

    :param families: An iterable of gene families.
    :param transitive: Size of the transitive closure used to build the graph.
    :param window_size: Size of the window for extracting gene contexts (default: 0).
    :param disable_bar: Flag to disable the progress bar (default: False).

    :return: The constructed gene context graph.
    """

    context_graph = nx.Graph()
    
    contig_to_genes_of_interest = get_contig_to_genes(families)
    
    for contig, genes_of_interest in  tqdm(contig_to_genes_of_interest.items(), unit="contig", total=len(contig_to_genes_of_interest), disable=disable_bar):
        logging.debug(f'Processing {len(genes_of_interest)} genes of interest in contig {contig}')
        
        genes_count = len(contig.genes)
        
        genes_of_interest_positions = [g.position for g in genes_of_interest]

        contig_windows = extract_contig_window(genes_count, genes_of_interest_positions, 
                                                            window_size=window_size, is_circular=contig.is_circular)
        
        add_edges_to_context_graph(context_graph,
                                contig.genes,
                                contig_windows,
                                transitive,
                                contig.is_circular)
    return context_graph


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
    """
    Export the results into dataFrame

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
                                     identity=args.identity, coverage=args.coverage, jaccard_threshold=args.jaccard, window_size=args.window_size,
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
    optional.add_argument("-w", "--window_size", required=False, type=int, default=5,
                        help="Number of neighboring genes that are considered on each side of "
                        "a gene of interest when searching for conserved genomic contexts.")
    
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
