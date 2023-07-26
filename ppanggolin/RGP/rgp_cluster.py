#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
import os
from itertools import combinations
from collections.abc import Callable
from collections import defaultdict
from multiprocessing.pool import Pool
from typing import Dict, List, Tuple, Set


# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd


# local libraries
from ppanggolin.genome import Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region
from ppanggolin.formats import check_pangenome_info
from ppanggolin.utils import restricted_float, mk_outdir


def compute_grr(rgp_a_families: set, rgp_b_families: set, mode: Callable) -> float:
    """
    Compute gene repertoire relatedness (GRR) between two rgp.
    mode can be the function min to compute min GRR or max to compute max_grr

    :param rgp_a: rgp A
    :param rgp_b: rgp B
    :param mode: min or max function

    :return : grr value between 0 and 1
    """

    grr = len((rgp_a_families & rgp_b_families)) / \
        mode(len(rgp_a_families), len(rgp_b_families))

    return grr


def compute_jaccard_index(rgp_a_families: set, rgp_b_families: set) -> float:
    """
    Compute jaccard index between two rgp based on their famillies.

    :param rgp_a: rgp A
    :param rgp_b: rgp B

    :return : jaccard index
    """

    jaccard_index = len((rgp_a_families & rgp_b_families)) / \
        len(rgp_a_families | rgp_b_families)

    return jaccard_index


def get_rgp_info_dict(regions: List[Region], region_to_spot: dict) -> Dict[int, dict]:
    """
    Format RGP information into a dictionary for adding to the graph.

    This function takes a list of RGPs and a dictionary mapping each RGP to its corresponding spot ID,
    and formats the RGP information into a dictionary for further processing or addition to a graph.

    :param regions: A list of RGPs.
    :param region_to_spot: A dictionary mapping each RGP to its corresponding spot ID.
    :return: A dictionary with RGP id as key and a dictionaries containing information on the corresponding RGP as value.
    """

    region_attributes = {}
    for region in regions:

        region_info = {"contig": region.contig.name,
                       'organism': region.organism.name,
                       "name": region.name,
                       "genes_count": len(region.genes),
                       "is_contig_border": region.is_contig_border,
                       "is_whole_contig": region.is_whole_contig,
                       "spot_id": str(region_to_spot.get(region, "No spot"))}
        
        region_info['families_count'] = len(region.families)

        region_attributes[region.ID] = region_info

    return region_attributes

# def add_rgp_metadata_to_graph(graph, pangenome):
#     """
#     """

#     source_fields = {m.source: m.fields for rgp in pangenome.regions if len(list(rgp.metadata)) > 0 for m in rgp.metadata}


#     for rgp in pangenome.regions:
        
#         for source_metadata_rgps in pangenome.metadata_sources("RGPs"):
#             metadata_source_count = 0 # counter to count how many metadata for the current source a rgp has. 
#             print(source_metadata_rgps)
#             to_concat = defaultdict(list)

#             for rgp_metadata in rgp.metadata:
#                 if rgp_metadata.source == source_metadata_rgps:
#                     metadata_source_count += 1

#                     for field in rgp_metadata.fields:
#                         to_concat[field].append(str(rgp_metadata.get(field)))
                        
#                 for field in source_fields[source_metadata_rgps]:
#                     concatenated_fields = '|'.join(to_concat[field])
#                     if concatenated_fields != "":
#                         graph.nodes[rgp.ID][f"{source_metadata_rgps}_{field}"] = concatenated_fields
            

def add_rgp_metadata_to_graph(graph, pangenome):
    """
    """
    
    metadata_sources = pangenome.metadata_sources("RGPs")

    for rgp in pangenome.regions:
        
        source_field_2_value = defaultdict(list)

        source_metadata_counter = defaultdict(int)
        for rgp_metadata in rgp.metadata:
            source = rgp_metadata.source
            source_metadata_counter[source] += 1

            for field in rgp_metadata.fields:
                source_field_2_value[f"{source}_{field}"].append(str(rgp_metadata.get(field)))
    
        # for source in metadata_sources:
        #     metadata_count = source_metadata_counter[source] # number of metadata associated to the current metadata for source
        #     graph.nodes[rgp.ID][f"{source}_metadata_count"] = metadata_count
            
        for source_field, values in source_field_2_value.items():
            concatenated_values = '|'.join(values)
            if concatenated_values != '':
                graph.nodes[rgp.ID][source_field] = concatenated_values

            
            
def add_identical_rgps_info(rgp_graph: nx.Graph, rgp_to_identical_rgps: Dict[Region, list]):
    """
    Add identical rgps info in the graph as node attributes.

    :params rgp_graph: Graph with rgp id as node and grr value as edges
    :params rgp_to_identical_rgps: dict with uniq RGP as key and set of identical rgps as value  
    """

    for rgp, identical_rgps in rgp_to_identical_rgps.items():
        rgp_graph.add_node(rgp.ID,
                           identical_rgp_count=len(identical_rgps),
                           identical_rgp_names=';'.join([rgp.name for rgp in identical_rgps]))


def differentiate_spot_in_identical_rgps(rgp_graph: nx.Graph, rgp_to_identical_rgps: Dict[Region, List[Region]], rgp_to_spot: Dict[Region, int]):
    """
    Adds nodes and edges to the `rgp_graph` for RGPs that have identical families but different spot ID.

    :param rgp_graph: a NetworkX graph representing RGP relationships
    :param rgp_to_identical_rgps: a dict mapping RGPs to sets of RGPs that are identical to them except for their spot ID
    :param rgp_to_spot: a dict mapping RGPs to their spot ID
    """
    identical_rgp_diff_spot_count = 0
    for rgp, identical_rgps in rgp_to_identical_rgps.items():
        # Create a defaultdict that maps spot IDs to sets of RGPs
        spot_to_rgps = defaultdict(set)
        for identical_rgp in identical_rgps:
            # For each identical RGP, add it to the set corresponding to its spot ID
            spot = rgp_to_spot.get(identical_rgp, None)
            spot_to_rgps[spot].add(identical_rgp)

        # For each set of strictly identical RGPs (i.e., all RGPs in the set have the same spot ID)
        for spot, strictly_identical_rgps in spot_to_rgps.items():
            # is the spot identical as the main rgp that have been used to compute grr?
            if spot == rgp_to_spot.get(rgp, None):
                rgp_graph.add_node(rgp.ID, identical_rgp_fam_and_spot=len(
                    strictly_identical_rgps)+1)
            else:
                # If the spot ID is different from the main RGP, add a node for one of the strictly identical RGPs
                # and set the `identical_rgp_fam_and_spot` and `identical_rgp` attributes.
                identical_rgp_diff_spot_count += 1
                strictly_identical_rgp = strictly_identical_rgps.pop()
                rgp_graph.add_node(strictly_identical_rgp.ID, identical_rgp_fam_and_spot=len(
                    strictly_identical_rgps)+1, identical_rgp=True)
                rgp_graph.add_edge(rgp.ID, strictly_identical_rgp.ID, grr=1.0, min_grr=1.0,
                                   max_grr=1.0, identical_famillies=True, different_spot=True)

    logging.info(f'{identical_rgp_diff_spot_count} RGPs with identical families but with different spot id have been added to the graph for visualisation purpose.')


def add_edges_to_identical_rgps(rgp_graph: nx.Graph, rgp_to_identical_rgps: Dict[Region, Set[Region]]):
    """
    Adds edges to the RGP graph between identical RGPs.

    :param rgp_graph: The RGP graph to add edges to.
    :param rgp_to_identical_rgps: A dictionary mapping RGPs to sets of identical RGPs.
    """

    identical_edge_data = {'grr': 1.0, 'max_grr': 1.0,
                           'min_grr': 1.0,
                           "identical_famillies": True}

    for rgp, identical_rgps in rgp_to_identical_rgps.items():

        rgp_graph.add_nodes_from([ident_rgp.ID for ident_rgp in identical_rgps if ident_rgp.ID != rgp.ID ], identical_rgp = True)

        # add edge between identical rgp with metrics at 1 (perfect score)
        edges_to_add = [(rgp_a.ID, rgp_b.ID, identical_edge_data)
                        for rgp_a, rgp_b in combinations(identical_rgps | {rgp}, 2)]

        # replicate all edges that connect main rgp to all identical rgps
        for connected_rgp in rgp_graph.neighbors(rgp.ID):
            edge_data = rgp_graph[rgp.ID][connected_rgp]
            edges_to_add += [(identical_rgp.ID, connected_rgp, edge_data)
                             for identical_rgp in identical_rgps]

        rgp_graph.add_edges_from(edges_to_add)


def dereplicate_rgp(rgps: list, disable_bar: bool = False) -> Dict[Region, set]:
    """
    Dereplicate RGPs that have the same families.

    :params rgps: list of rgps
    :param disable_bar: Disable progress bar

    :return : dict with uniq RGP as key and set of identical rgps as value  
    """
    logging.info(f'Dereplicating {len(rgps)} RGPs')
    families_to_rgps = defaultdict(list)

    for rgp in tqdm(rgps, total=len(rgps), unit="RGP", disable=disable_bar):
        families_to_rgps[tuple(
            sorted((f.ID for f in rgp.families)))].append(rgp)

    uniq_rgps = {}
    for rgps in families_to_rgps.values():

        uniq_rgps[rgps[0]] = set(rgps)

    logging.info(f'{len(uniq_rgps)} uniq RGPs')
    return uniq_rgps


def compute_rgp_metric(rgp_pair: Tuple[int, int],
                       rgp_to_families: Dict[int, set],
                       rgp_to_contigborder: Dict[int, bool],
                       grr_cutoff: float,
                       grr_metric: str) -> Tuple[int, int, dict]:
    """
    Compute GRR metric between two RGPs.

    :param rgp_pair: Pair of RGP IDs to compute the metric for
    :param rgp_to_families: Mapping of RGP ID to set of families it contains
    :param rgp_to_contigborder: Mapping of RGP ID to whether it is at a contig border or not
    :param grr_cutoff: Minimum GRR value required for the RGP pair to be considered
    :param grr_metric: grr mode between min_grr, max_grr and incomplete_aware_grr

    :returns: Tuple containing the IDs of the two RGPs and the computed metrics as a dictionary
    """

    # Unpack RGP pair
    rgp_a, rgp_b = rgp_pair

    # Get families for each RGP
    rgp_a_fam = rgp_to_families[rgp_a]
    rgp_b_fam = rgp_to_families[rgp_b]

    # Compute metrics between 2 rgp if they share at least one family
    if rgp_a_fam & rgp_b_fam:
        edge_metrics = {}

        # RGP at a contig border are seen as incomplete and min GRR is used instead of max GRR
        if rgp_to_contigborder[rgp_a] or rgp_to_contigborder[rgp_b]:
            edge_metrics["incomplete_aware_grr"] = compute_grr(
                rgp_a_fam, rgp_b_fam, min)
        else:
            edge_metrics["incomplete_aware_grr"] = compute_grr(
                rgp_a_fam, rgp_b_fam, max)

        # Compute max and min GRR metrics
        edge_metrics['max_grr'] = compute_grr(rgp_a_fam, rgp_b_fam, max)
        edge_metrics['min_grr'] = compute_grr(rgp_a_fam, rgp_b_fam, min)

        # Number of shared families can be useful when visualising the graph
        edge_metrics['shared_family'] = len(rgp_a_fam & rgp_b_fam)

        # Only return the metrics if the GRR value is above the cutoff
        if edge_metrics[grr_metric] >= grr_cutoff:
            return (rgp_a, rgp_b, edge_metrics)


def launch_rgp_metric(pack: tuple) -> tuple:
    """ 
    Allow to launch in multiprocessing the rgp metric

    :param pack: Pack of argument for rgp metrics

    :return: edge metrics 
    """
    return compute_rgp_metric(*pack)


def get_rgp_family_ids(rgp: Region) -> dict:
    """
    Get the set of family IDs contained within an RGP.

    :param rgp: The RGP object to get family IDs from.
    :return: A set of family IDs contained within the RGP.
    """
    return {f.ID for f in rgp.families}


def cluster_rgp_on_grr(G: nx.Graph, clustering_attribute: str = "grr"):
    """
    Cluster rgp based on grr using louvain communities clustering. 

    :param G: NetworkX graph object representing the RGPs and their relationships
    :param clustering_attribute: Attribute of the graph to use for clustering (default is "grr")
    """

    partitions = nx.algorithms.community.louvain_communities(
        G, weight=clustering_attribute)

    # Add partition index in node attributes
    for i, cluster_nodes in enumerate(partitions):
        nx.set_node_attributes(
            G, {node: f"cluster_{i}" for node in cluster_nodes}, name=f"{clustering_attribute}_cluster")

    logging.info(
        f"Graph has {len(partitions)} clusters using {clustering_attribute}")


def write_rgp_cluster_table(outfile: str, grr_graph: nx.Graph,
                            rgp_to_identical_rgps: Dict[Region, Set[Region]],
                            grr_metric: str,
                            rgp_to_spot: Dict[Region, int]) -> None:
    """
    Writes RGP cluster info to a TSV file using pandas.

    :param outfile: name of the tsv file
    :param grr_graph: The GRR graph.
    :param rgp_to_identical_rgps: A dictionary mapping an RGP to a set of identical RGPs.
    :param grr_metric: The GRR metric used for clustering.
    :param rgp_to_spot: A dictionary mapping an RGP to its spot.
    :return: None
    """

    def get_spot_id(rgp:Region, rgp_to_spot:Dict[Region, int]) -> str:
        """
        Return Spot ID associated to an RGP. 
        It adds the prefix "spot_" to the spot ID.
        When no spot is associated to the RGP, then the string "No spot" is return 

        :params rgp: RGP id
        :params rgp_to_spot: A dictionary mapping an RGP to its spot .

        :return: Spot ID of the given RGP with the prefix spot_ or "No spot". 
        """
        if rgp in rgp_to_spot:
            return f"spot_{rgp_to_spot[rgp]}"
        else:
            return "No spot"

    all_rgps_infos = []
    for rgp, identical_rgps in rgp_to_identical_rgps.items():

        cluster = grr_graph.nodes[rgp.ID][f'{grr_metric}_cluster']
        all_rgps_infos += [{"RGP": r.name, "cluster": cluster,
                            "spot_id": get_spot_id(rgp, rgp_to_spot)} for r in identical_rgps]

    df = pd.DataFrame(all_rgps_infos)
    df.to_csv(outfile, sep='\t', index=False)


def cluster_rgp(pangenome, grr_cutoff, output, basename, cpu, ignore_incomplete_rgp, unmerge_identical_rgps, grr_metric, disable_bar):
    """
    Main function to cluster regions of genomic plasticity based on their GRR

    :param pangenome: pangenome object
    :param output: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """

    # check statuses and load info
    check_pangenome_info(pangenome, need_families=True, need_annotations=True,
                         disable_bar=disable_bar, need_rgp=True, need_spots=True, need_metadata=True, metatype="RGPs")

    if pangenome.regions == 0:
        raise Exception(
            "The pangenome has no RGPs. The clustering of RGP is then not possible.")

    grr_graph = nx.Graph()

    # add all rgp as node
    if ignore_incomplete_rgp:
        valid_rgps = [
            rgp for rgp in pangenome.regions if not rgp.is_contig_border]
        ignored_rgp_count = len(pangenome.regions) - len(valid_rgps)
        total_rgp_count = len(pangenome.regions)
        logging.info(f'Ignoring {ignored_rgp_count}/{total_rgp_count} ({100*(ignored_rgp_count)/total_rgp_count:.2f}%) '
                     'RGPs that are located at a contig border and are likely incomplete.')
        if len(valid_rgps) == 0:
            raise Exception(
                "The pangenome has no complete RGPs. The clustering of RGP is then not possible.")
    else:
        valid_rgps = pangenome.regions

    rgp_to_identical_rgps = dereplicate_rgp(
        valid_rgps, disable_bar=disable_bar)

    uniq_rgps = list(rgp_to_identical_rgps)
    rgp_count = len(uniq_rgps)

    grr_graph.add_nodes_from(
        (rgp.ID for rgp in uniq_rgps), identical_rgp=False)

    # Creating dictonnaries paring rgp ID with their families ids
    rgp_to_families = {rgp.ID: get_rgp_family_ids(
        rgp) for rgp in uniq_rgps}
    
    rgp_to_iscontigborder = {
        rgp.ID: rgp.is_contig_border for rgp in uniq_rgps}

    # compute grr for all possible pair of rgp
    pairs_count = int((rgp_count**2 - rgp_count)/2)
    logging.info(
        f'Computing GRR metric for {pairs_count:,} pairs of RGP using {cpu} cpus...')

    # use the ID of the rgp to make pair rather than their name to save memory
    rgp_pairs = combinations(rgp_to_families, 2)

    ideal_chunk_size = 50000
    chunk_count = (pairs_count/ideal_chunk_size) + cpu
    chunk_size = int(pairs_count/chunk_count)+1
    logging.debug(
        f'Computing GRR metric in ~{chunk_count:.2f} chunks of {chunk_size} pairs')

    # create the argument iterator to use in parallell
    arg_iter = ((rgp_pair, rgp_to_families, rgp_to_iscontigborder,
                grr_cutoff, grr_metric) for rgp_pair in rgp_pairs)

    pairs_of_rgps_metrics = []
    with Pool(processes=cpu) as p:
        for pair_metrics in tqdm(p.imap_unordered(launch_rgp_metric, arg_iter, chunksize=chunk_size),
                                 unit="pair of RGPs", total=pairs_count, disable=disable_bar):
            if pair_metrics:
                pairs_of_rgps_metrics.append(pair_metrics)

    grr_graph.add_edges_from(pairs_of_rgps_metrics)

    if unmerge_identical_rgps:
        add_edges_to_identical_rgps(grr_graph, rgp_to_identical_rgps)

    # cluster rgp based on grr value
    logging.info(
        f"Louvain_communities clustering of RGP  based on {grr_metric} on {grr_graph}.")

    cluster_rgp_on_grr(grr_graph, grr_metric)

    rgp_to_spot = {region: int(spot.ID)
                   for spot in pangenome.spots for region in spot.regions}

    if not unmerge_identical_rgps:
        logging.info(f"Add info in the graph on identical RGPs merged")
        add_identical_rgps_info(grr_graph, rgp_to_identical_rgps)
        differentiate_spot_in_identical_rgps(
            grr_graph, rgp_to_identical_rgps, rgp_to_spot)

    # add some attribute to the graph nodes.
    logging.info(f"Add RGP information to the graph")
    region_infos = get_rgp_info_dict(pangenome.regions, rgp_to_spot)

    nx.set_node_attributes(grr_graph, region_infos)

    add_rgp_metadata_to_graph(grr_graph, pangenome)

    # writting graph in gexf format
    graph_file_name = os.path.join(output, f"{basename}.gexf")
    logging.info(f"Writting graph in gexf format in {graph_file_name}.")
    nx.readwrite.gexf.write_gexf(grr_graph, graph_file_name)

    nx.readwrite.graphml.write_graphml(
        grr_graph, os.path.join(output, f"{basename}.graphml"))

    outfile = os.path.join(output, f"{basename}.tsv")
    logging.info(f"Writting rgp clusters in tsv format in {outfile}")
    write_rgp_cluster_table(
        outfile, grr_graph, rgp_to_identical_rgps, grr_metric, rgp_to_spot)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()

    mk_outdir(args.output, args.force)

    pangenome.add_file(args.pangenome)

    cluster_rgp(pangenome, grr_cutoff=args.grr_cutoff, output=args.output,
                basename=args.basename, cpu=args.cpu, ignore_incomplete_rgp=args.ignore_incomplete_rgp,
                unmerge_identical_rgps=args.no_identical_rgp_merging,
                grr_metric=args.grr_metric, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for cluster_rgp command

    :return : parser arguments for cluster_rgp command
    """
    parser = sub_parser.add_parser(
        "rgp_cluster", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_cluster_rgp(parser)
    return parser


def parser_cluster_rgp(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of rgp command

    :param parser: parser for cluster_rgp argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True,
                          type=str, help="The pangenome .h5 file")

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument('--grr_cutoff', required=False, type=restricted_float, default=0.8,
                          help="Min gene repertoire relatedness metric used in the rgp clustering")
    optional.add_argument('--grr_metric', required=False, type=str, default="incomplete_aware_grr",
                          help="The grr (Gene Repertoire Relatedness) is used to assess the similarity between two RGPs based on their gene families. "
                          "There are three different modes for calculating the grr value: 'min_grr', 'max_grr' or  'incomplete_aware_grr'."
                          " 'min_grr': Computes the number of gene families shared between the two RGPs and divides it by the smaller number of gene families among the two RGPs. "
                          " 'max_grr': Calculates the number of gene families shared between the two RGPs and divides it by the larger number of gene families among the two RGPs. "
                          " 'incomplete_aware_grr' (default): If at least one RGP is considered incomplete, which occurs when it is located at the border of a contig, "
                          "the 'min_grr' mode is used. Otherwise, the 'max_grr' mode is applied.",
                          choices=['incomplete_aware_grr', "min_grr", "max_grr"])

    optional.add_argument('--ignore_incomplete_rgp', required=False, action="store_true",
                          help="Do not cluster RGPs located on a contig border which are likely incomplete.")

    optional.add_argument('--no_identical_rgp_merging', required=False, action="store_true",
                          help="Do not merge in one node identical RGP (i.e. having the same family content) before clustering.")

    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")

    optional.add_argument("--basename", required=False,
                          default="rgp_cluster", help="basename for the output file")

    optional.add_argument('-o', '--output', required=False, type=str,
                          default="rgp_clustering", help="Output directory")
