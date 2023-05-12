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
                                     cpu: int = 1, write_context_graph:bool = False, disable_bar=True):
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
    :param write_context_graph: Write graph of the contexts
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
    
    logging.getLogger().debug(f"Context graph made of {nx.number_of_nodes(gene_context_graph)} nodes and {nx.number_of_edges(gene_context_graph)} edges")

    compute_edge_metrics(gene_context_graph, jaccard_threshold)

    # Filter graph 
    filter_flag = f'is_jaccard_gene_>_{jaccard_threshold}'
    
    filtered_graph = nx.subgraph_view(gene_context_graph, filter_edge=lambda n1, n2: gene_context_graph[n1][n2][filter_flag] )

    logging.getLogger().debug(f"Filtering context graph on {filter_flag}")
    logging.getLogger().debug(f"Context graph made of {nx.number_of_nodes(filtered_graph)} nodes and {nx.number_of_edges(filtered_graph)} edges")
    connected_components = nx.connected_components(filtered_graph)

    # Connected component graph Filtering

    # remove singleton famillies
    connected_components = (component for component in connected_components if len(component) > 1)  

    # remove component made only of famillies not initially requested
    connected_components = (component for component in connected_components if component & gene_families)
    
    gene_contexts = {GeneContext(gc_id=i, families=component) for i, component in enumerate(connected_components) }

    families_in_contexts = {family for gene_context in gene_contexts for family in gene_context.families}
    
    graph_with_final_contexts = nx.subgraph_view(gene_context_graph, filter_node=lambda n: n in families_in_contexts)

    if write_context_graph:
        write_graph(graph_with_final_contexts, output, gene_families, gene_contexts)

    if len(families_in_contexts) != 0:
        logging.getLogger().debug(f"There are {len(families_in_contexts)} families among {len(gene_contexts)} gene contexts")
        
        output_file = os.path.join(output, "gene_contexts.tsv")

        export_context_to_dataframe(gene_contexts, fam_2_seq, output_file)

    else:
        logging.getLogger().info(f"No gene contexts were found")

    logging.getLogger().info(f"Computing gene contexts took {round(time.time() - start_time, 2)} seconds")


    # # Finding connected components with panmodule functions 
    # # extract the modules from the graph

    # logging.getLogger().debug(f"panmodule style:")
    # gene_contexts_pandmodule_way = compute_gene_context(gene_context_graph, jaccard_threshold)

    # # remove singleton famillies
    # gene_contexts_pandmodule_way = (context for context in gene_contexts_pandmodule_way if len(context.families) >  1 )  

    # # remove component made only of famillies not initially requested
    # gene_contexts_pandmodule_way = [context for context in gene_contexts_pandmodule_way if context.families & gene_families]

    # families_in_contexts = {family for gene_context in gene_contexts_pandmodule_way for family in gene_context.families}
    
    # logging.getLogger().debug(f"There are {len(families_in_contexts)} families among {len(gene_contexts_pandmodule_way)} gene contexts")


    # output_file = os.path.join(output, "gene_contexts_panmodule_style.tsv")
    # export_context_to_dataframe(gene_contexts_pandmodule_way, fam_2_seq, output_file)



def write_graph(context_graph: nx.Graph, output_dir: str, famillies_of_interest: Set[GeneFamily], gene_contexts:List[GeneContext]):
    """
    Write a graph to file with node and edge attributes. 
    
    This function writes a graph to a file in the GraphML format or in GEXF format. The original context graph contains 
    ppanggolin objects as nodes and lists and dictionaries in edge attributes. Since these objects 
    cannot be written to the output graph, this function creates a new graph that contains only 
    writable objects.

    :param context_graph: A NetworkX Graph object representing the graph.
    :param output_dir: The output directory where the graph file will be written.
    :param famillies_of_interest: A list of node objects that are of interest.
    :param gene_contexts: List of gene context, used to add context id to node of the graph

    """
    def filter_attribute(data:dict):
        """
        Helper function to filter the edge attributes.

        :param data: The edge attribute data.
        :return: A filtered dictionary containing only non-collection attributes.
        """
        return {k:v for k, v in data.items() if type(v) not in [set, dict, list]}
    
    G = nx.Graph()

    G.add_edges_from(((f1.name, f2.name) for f1,f2 in context_graph.edges()))

    edges_with_attributes = {(f1.name, f2.name):filter_attribute(d) for f1,f2,d in context_graph.edges(data=True)}
    
    nx.set_edge_attributes(G, edges_with_attributes)
    
    nodes_attributes_filtered = {f.name:filter_attribute(d) for f,d in context_graph.nodes(data=True)}

    # on top of attributes already contained in node of context graph
    # add organisms and genes count that have the family, the partition and if the family was in initially requested 
    nodes_family_data = {f.name:{"organisms":len(f.organisms), 
                                 "partition":f.partition,
                          "genes":len(f.genes), 
                          "famillies_of_interest": f in famillies_of_interest} for f in context_graph.nodes()}
    
    family_name_to_context_id = {family.name:context.ID for context in gene_contexts for family in context.families}

    for f, d in G.nodes(data=True):
        d.update(nodes_family_data[f])
        d.update(nodes_attributes_filtered[f])
        d['context_id'] = family_name_to_context_id[f]
        
    
    graphml_file = os.path.join(output_dir, "graph_context.graphml")
    logging.info(f'Writting context graph in {graphml_file}')
    nx.write_graphml_lxml(G, graphml_file)

    gexf_file = os.path.join(output_dir, "graph_context.gexf")
    logging.info(f'Writting context graph in {gexf_file}')
    nx.readwrite.gexf.write_gexf(G, gexf_file)


def compute_edge_metrics(context_graph: nx.Graph, gene_proportion_cutoff: float) -> None:
    """
    Compute various metrics on the edges of the context graph.

    :param context_graph: The context graph.
    :param gene_proportion_cutoff: The minimum proportion of shared genes between two features for their edge to be considered significant.
    """
    # compute jaccard on organism and on genes
    for f1, f2, data in context_graph.edges(data=True):
        
        data['jaccard_organism'] = len(data['organisms'])/len(f1.organisms | f2.organisms)
        
        f1_gene_proportion = len(data['genes'][f1])/len(f1.genes)
        f2_gene_proportion = len(data['genes'][f2])/len(f2.genes)
        
        data[f'f1'] = f1.name
        data[f'f2'] = f2.name
        data[f'f1_jaccard_gene'] = f1_gene_proportion
        data[f'f2_jaccard_gene'] = f2_gene_proportion
                        
        data[f'is_jaccard_gene_>_{gene_proportion_cutoff}'] = (f1_gene_proportion >= gene_proportion_cutoff) and (f2_gene_proportion >= gene_proportion_cutoff)
        
        # the following commented out lines are additional metrics that could be used

        # data['min_jaccard_organism'] = len(data['organisms'])/min(len(f1.organisms), len(f2.organisms))
        # data['max_jaccard_organism'] = len(data['organisms'])/max(len(f1.organisms), len(f2.organisms))
        # f1_gene_proportion_partial = len(data['genes'][f1])/len(context_graph.nodes[f1]['genes'])
        # f2_gene_proportion_partial = len(data['genes'][f2])/len(context_graph.nodes[f2]['genes'])
        # data[f'f1_jaccard_gene_partital'] = f1_gene_proportion_partial
        # data[f'f2_jaccard_gene_partital'] = f2_gene_proportion_partial          


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

            for i, next_gene_index in enumerate(next_genes):
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

                if i == 0:
                    context_graph[gene.family][next_gene.family]['adjacent_family'] = True


                # Add node attributes
                node_gene_dict = context_graph.nodes[gene.family]
                next_gene_gene_dict = context_graph.nodes[next_gene.family]

                increment_attribute_counter(node_gene_dict, "genes_count")
                increment_attribute_counter(next_gene_gene_dict, "genes_count")

                add_val_to_dict_attribute(node_gene_dict, "genes", gene)
                add_val_to_dict_attribute(next_gene_gene_dict, "genes", next_gene)


                # Add edge attributes
                edge_dict = context_graph[gene.family][next_gene.family]
                try:
                    genes_edge_dict = edge_dict['genes']
                except:
                    genes_edge_dict = {}
                    edge_dict['genes'] = genes_edge_dict
                
                add_val_to_dict_attribute(genes_edge_dict, gene.family, gene)
                add_val_to_dict_attribute(genes_edge_dict, next_gene.family, next_gene)

                add_val_to_dict_attribute(edge_dict, "organisms", gene.organism)

                increment_attribute_counter(edge_dict, "gene_pairs")
                
                assert gene.organism == next_gene.organism


def add_val_to_dict_attribute(attr_dict: dict, attribute_key, attribute_value):
    """
    Add an attribute value to a edge or node dictionary set.

    :param attr_dict: The dictionary containing the edge/node attributes.
    :param attribute_key: The key of the attribute.
    :param attribute_value: The value of the attribute to be added.

    """

    try:
        attr_dict[attribute_key].add(attribute_value)
    except KeyError:
        attr_dict[attribute_key] = {attribute_value}


def increment_attribute_counter(edge_dict: dict, key:Hashable):
    """
    Increment the counter for an edge/node attribute in the edge/node dictionary.

    :param edge_dict: The dictionary containing the attributes.
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


# def compute_gene_context(g: nx.Graph, jaccard: float = 0.85) -> set:
#     """
#     Compute the gene contexts in the graph

#     :param g: Graph of gene contexts between interesting gene families of the pan
#     :param jaccard: Jaccard index

#     :return: Set of gene contexts find in graph
#     """

#     gene_contexts = set()
#     c = 1
#     for comp in connected_components(g, removed=set(), weight=jaccard):
#         gene_contexts.add(GeneContext(gc_id=c, families=comp))
#         c += 1
#     return gene_contexts


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


def export_context_to_dataframe(gene_contexts: set, fam_to_seq: dict, output: str):
    """
    Export the results into dataFrame

    :param gene_contexts: connected components found in the pan
    :param fam_to_seq: Dictionary with gene families as keys and list of sequence ids as values
    :param output: output path
    """

    lines = []
    for gene_context in gene_contexts:
        for family in gene_context.families:

            family_info = {"GeneContext ID":gene_context.ID,
                           "Gene family name": family.name,
                           "Sequence ID":None if fam_to_seq is None else ','.join(fam_to_seq.get(family.ID)),
                           "Nb Genomes":len(family.organisms),
                           "Partition": family.named_partition  }
            lines.append(family_info)
            
    df = pd.DataFrame(lines).set_index("GeneContext ID")
    
    df = df.sort_values(["GeneContext ID", "Sequence ID"], na_position='last')

    df.to_csv(output, sep="\t", na_rep='NA')

    logging.getLogger().debug(f"detected gene context(s) are listed in: '{output}")


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
                                     no_defrag=args.no_defrag, cpu=args.cpu, write_context_graph=args.write_graph, disable_bar=args.disable_prog_bar)


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
    optional.add_argument('--write_graph',  action="store_true",
                    help="Write context graph in GEXF format.")

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
