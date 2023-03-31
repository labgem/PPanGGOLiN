#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
import time
import os
from itertools import combinations, product
from collections.abc import Callable
from collections import defaultdict 
from multiprocessing.pool import Pool
from itertools import islice
import time

# installed libraries
from tqdm import tqdm
import networkx as nx

# local libraries
from ppanggolin.genome import Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import restricted_float, mk_outdir

def compute_grr(rgp_a_families:set, rgp_b_families:set, mode:Callable) -> float:
    """
    Compute gene repertoire relatedness (GRR) between two rgp.
    mode can be the function min to compute min GRR or max to compute max_grr

    :param rgp_a: rgp A
    :param rgp_b: rgp B
    :param mode: min or max function

    :return : grr value between 0 and 1
    """
    
    max_grr = len((rgp_a_families & rgp_b_families))/mode(len(rgp_a_families), len(rgp_b_families))

    return max_grr

def compute_jaccard_index(rgp_a_families:set, rgp_b_families:set) -> float:
    """
    Compute jaccard index between two rgp based on their famillies.

    :param rgp_a: rgp A
    :param rgp_b: rgp B

    :return : jaccard index
    """
    
    jaccard_index = len((rgp_a_families & rgp_b_families))/len(rgp_a_families | rgp_b_families)

    return jaccard_index


def get_rgp_info_dict(regions:list, region_to_spot:set) -> nx.Graph :
    """

    """

    region_attributes = {}
    for region in regions:

        region_info = {"contig":region.contig.name,
                       'organism':region.organism.name,
                        "name":region.name,
                        "genes_count":len(region.genes),
                        "is_contig_border":region.is_contig_border, 
                        "is_whole_contig":region.is_whole_contig,
                        "spot_id":region_to_spot.get(region, "No spot")}
        
        region_info['famillies'] =';'.join([str(f.ID)  for f in region.families])
        region_info['families_count'] = len(region.families)

        region_attributes[region.ID] = region_info

    return region_attributes


def add_identical_rgps_info(rgp_graph, rgp_to_identical_rgps):
    """
    """

    for rgp, identical_rgps in rgp_to_identical_rgps.items():
        rgp_graph.add_node(rgp.ID, 
                            identical_rgp_count=len(identical_rgps)+1, # +1 is to count the main rgp use as node in the graph
                            identical_rgp_names= ';'.join([rgp.name for rgp in identical_rgps])) 


def differentiate_spot_in_identical_rgps(rgp_graph, rgp_to_identical_rgps, rgp_to_spot):
    """
    """
    identical_rgp_diff_spot_count = 0
    for rgp, identical_rgps in rgp_to_identical_rgps.items():
        # print(rgp.name, rgp.ID, rgp_to_spot.get(rgp, None))
        # print('IDENTICAL CONTENT', [(rgp.name, get_rgp_family_ids(rgp), rgp_to_spot.get(rgp, None)) for rgp in identical_rgps])
        spot_to_rgps = defaultdict(set)
        for identical_rgp in identical_rgps:
            spot = rgp_to_spot.get(identical_rgp, None)
            spot_to_rgps[spot].add(identical_rgp)

        # print()
        for spot, strictly_identical_rgps in spot_to_rgps.items():
            if spot == rgp_to_spot.get(rgp, None): # is the spot identical as the main rgp that have been used to compute grr?
                rgp_graph.add_node(rgp.ID, identical_rgp_fam_and_spot=len(strictly_identical_rgps))
            else:
                identical_rgp_diff_spot_count += 1
                strictly_identical_rgp = strictly_identical_rgps.pop()
                rgp_graph.add_node(strictly_identical_rgp.ID, identical_rgp_fam_and_spot=len(strictly_identical_rgps)+1, identical_rgp=True)
                rgp_graph.add_edge(rgp.ID, strictly_identical_rgp.ID, grr=1.0, min_grr=1.0, max_grr=1.0, jaccard_index=1.0, identical_famillies=True, different_spot=True)
        
    logging.info(f'{identical_rgp_diff_spot_count} identical RGPs but with different spot id have been added to the graph for visualisation purpose.')

def add_edges_to_identical_rgps(rgp_graph, rgp_to_identical_rgps):
    """
    """

    identical_edge_data = {'grr': 1.0, 'max_grr': 1.0, 
                           'min_grr': 1.0, 'jaccard_index': 1.0,
                           "identical_famillies":True}

    for rgp, identical_rgps in rgp_to_identical_rgps.items():
        # print(rgp.name)
        # print("identical", [rgp.name for rgp in identical_rgps])


        if not identical_rgps:
            continue
        
        # add edge between identical rgp with metrics at 1 (perfect score)
        edges_to_add = [(rgp_a.ID, rgp_b.ID, identical_edge_data) for rgp_a, rgp_b in combinations(identical_rgps | {rgp}, 2)]

        # replicate all edges that connect main rgp to all identical rgps

        for connected_rgp in rgp_graph.neighbors(rgp.ID):
            edge_data = rgp_graph[rgp.ID][connected_rgp]
            edges_to_add += [(identical_rgp.ID, connected_rgp, edge_data) for identical_rgp in identical_rgps]

        rgp_graph.add_edges_from(edges_to_add)





def dereplicate_rgp(rgps: list, disable_bar=False) -> dict:
    """
    Dereplicate RGPs.

    :params rgps:
    :return : dict with uniq RGP as key and set of identical rgps as value  
    """
    logging.info(f'Dereplicating {len(rgps)} RGPs')
    families_to_rgps = defaultdict(set)

    for rgp in tqdm(rgps, total=len(rgps), unit="RGP", disable=disable_bar):
        families_to_rgps[tuple(sorted((f.ID for f in rgp.families)))].add(rgp)

    uniq_rgps = {}
    for rgps in families_to_rgps.values():

        uniq_rgps[rgps.pop()] = rgps

    logging.info(f'{len(uniq_rgps)} uniq RGPs')
    return uniq_rgps

def compute_rgp_metric(rgp_pair, rgp_to_families, rgp_to_contigborder, grr_cutoff):
    """
    """
    # rgp_pair, rgp_to_families = rgp_pair_and_rgp_to_families
    rgp_a, rgp_b = rgp_pair
    
    rgp_a_fam = rgp_to_families[rgp_a]
    rgp_b_fam = rgp_to_families[rgp_b]

    # compute metrics between 2 rgp if they share at least one familly
    if rgp_a_fam & rgp_b_fam:
        edge_metrics = {}
        
        if rgp_to_contigborder[rgp_a] or rgp_to_contigborder[rgp_b]:
            # RGP at a contig border are seen as incomplete and min GRR is used instead of max GRR
            edge_metrics['grr']  = compute_grr(rgp_a_fam, rgp_b_fam, min)
        else:
            edge_metrics['grr']  = compute_grr(rgp_a_fam, rgp_b_fam, max)

        edge_metrics['max_grr']  = compute_grr(rgp_a_fam, rgp_b_fam, max)
        edge_metrics['min_grr'] = compute_grr(rgp_a_fam, rgp_b_fam, min)
        edge_metrics['jaccard_index']  = compute_jaccard_index(rgp_a_fam, rgp_b_fam)

        # number of shared fam can be useful when visualising the graph
        edge_metrics['shared_family']  = len(rgp_a_fam & rgp_b_fam)

        if edge_metrics['grr'] >= grr_cutoff:
            return (rgp_a, rgp_b, edge_metrics)


def launch_rgp_metric(pack: tuple) -> tuple:
    """ Allow to launch in multiprocessing the rgp metric

    :param pack: Pack of argument for rgp metrics

    :return: edge metrics 
    """
    return compute_rgp_metric(*pack)


def get_rgp_family_ids(rgp):
    return {f.ID for f in rgp.families}

def cluster_nodes(G, clustering_attributes, prefix=""):
            
        for weight in clustering_attributes:
            partitions = nx.algorithms.community.louvain_communities(G, weight=weight)
        
            # Add partition index in node attributes
            for i, cluster_nodes in enumerate(partitions):
                nx.set_node_attributes(G, {node:f"cluster_{i}" for node in cluster_nodes}, name=f"{prefix}_{weight}_cluster")

            logging.info(f"{prefix}: Graph has {len(partitions)} clusters using {weight}")


def cluster_rgp(pangenome, grr_cutoff, output, basename, cpu, ignore_incomplete_rgp, unmerge_identical_rgps, disable_bar):
    """
    Main function to cluster regions of genomic plasticity based on their GRR

    :param pangenome: pangenome object
    :param output: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """

    # check statuses and load info
    check_pangenome_info(pangenome, need_families=True, need_annotations=True,
                        disable_bar=disable_bar, need_rgp=True, need_spots=True)


    if pangenome.regions == 0:
        raise Exception("The pangenome has no RGPs. The clustering of RGP is then not possible.")
    
    grr_graph = nx.Graph()

    # add all rgp as node
    if ignore_incomplete_rgp:
        valid_rgps = [rgp for rgp in pangenome.regions if not rgp.is_contig_border]
        ignored_rgp_count = len(pangenome.regions) - len(valid_rgps)
        total_rgp_count = len(pangenome.regions)
        logging.info(f'Ignoring {ignored_rgp_count}/{total_rgp_count} ({100*(ignored_rgp_count)/total_rgp_count:.2f}%) '
                      'RGPs that are located at a contig border and are likely incomplete.')
        if len(valid_rgps) == 0:
            raise Exception("The pangenome has no complete RGPs. The clustering of RGP is then not possible.")
    else:
        valid_rgps = pangenome.regions
        

    rgp_to_identical_rgps = dereplicate_rgp(valid_rgps, disable_bar=disable_bar)
    
    uniq_rgps = list(rgp_to_identical_rgps)
    rgp_count = len(uniq_rgps)


    grr_graph.add_nodes_from((rgp.ID for rgp in uniq_rgps), identical_rgp=False)


    # Creating dictonnaries paring rgp ID with their families ids
    rgp_to_families = {rgp.ID:get_rgp_family_ids(rgp) for rgp_index, rgp in enumerate(uniq_rgps)}
    rgp_to_iscontigborder = {rgp.ID:rgp.is_contig_border for rgp_index, rgp in enumerate(uniq_rgps)}

    # compute grr for all possible pair of rgp
    pairs_count = int((rgp_count**2 - rgp_count)/2)
    logging.info(f'Computing GRR metric for {pairs_count:,} pairs of RGP using {cpu} cpus...')

     # use the ID of the rgp to make pair rather than their name to save memory 
    rgp_pairs = combinations(rgp_to_families, 2)

    ideal_chunk_size = 50000
    chunk_count = (pairs_count/ideal_chunk_size) + cpu
    chunk_size = int(pairs_count/chunk_count )+1 
    logging.debug(f'Processing RGP pairs in  ~{chunk_count:.2f} chunks of {chunk_size} pairs')
    

    arg_iter = ((rgp_pair, rgp_to_families, rgp_to_iscontigborder, grr_cutoff) for rgp_pair in rgp_pairs)

    pairs_of_rgps_metrics = []
    with Pool(processes=cpu) as p:
        for pair_metrics in tqdm(p.imap_unordered(launch_rgp_metric, arg_iter, chunksize=chunk_size),
                                 unit="pair of RGPs", total=pairs_count, disable=disable_bar):
            if pair_metrics:
                pairs_of_rgps_metrics.append(pair_metrics)

    grr_graph.add_edges_from(pairs_of_rgps_metrics)

    # cluster rgp based on grr
    logging.info(f"Louvain_communities clustering of RGP on {grr_graph}.")
    clustering_attributes = ["grr"]
    cluster_nodes(grr_graph, clustering_attributes, prefix=f"merged_rgp") 

    if unmerge_identical_rgps:
        add_edges_to_identical_rgps(grr_graph, rgp_to_identical_rgps)

        # cluster rgp based on grr
        logging.info(f"Louvain_communities clustering of RGP on {grr_graph}.")
        clustering_attributes = ["grr"]
        cluster_nodes(grr_graph, clustering_attributes, prefix=f"unmerged_rgp") 


    # logging.info(f"Clustering RGPs on filtered graph")
    # thresolds = [0.6, 0.7, 0.8, 0.9]
    # for edge_metric, threshold in product(clustering_attributes, thresolds):
    #     edges_to_rm = [(u,v) for u,v,e in grr_graph.edges(data=True) if e[edge_metric] < threshold]
    #     grr_graph_filtered = nx.restricted_view(grr_graph, nodes=[], edges=edges_to_rm )
    #     logging.info(f"Filtering graph edges with {edge_metric}<{threshold}: {grr_graph_filtered}")

    #     # cluster filtered graph
    #     cluster_nodes(grr_graph_filtered, [edge_metric], prefix=f"{edge_metric}>{threshold}")


    logging.info(f"Manage identical RGP in the graph")
    rgp_to_spot =  {region:spot.ID  for spot in pangenome.spots  for region in spot.regions}
    add_identical_rgps_info(grr_graph, rgp_to_identical_rgps) 

    differentiate_spot_in_identical_rgps(grr_graph, rgp_to_identical_rgps, rgp_to_spot)
    
    

    # add some attribute to the graph nodes.
    logging.info(f"Add RGP information to the graph")
    region_infos = get_rgp_info_dict(pangenome.regions, rgp_to_spot)
    nx.set_node_attributes(grr_graph, region_infos)
    
    # logging.info(f"Writting graph in graphml format with multiple thresholds.")

    # for edge_metric, threshold in product(clustering_attributes, thresolds):
    #     edges_to_rm = [(u,v) for u,v,e in grr_graph.edges(data=True) if e[edge_metric] < threshold]
    #     grr_graph_filtered = nx.restricted_view(grr_graph, nodes=[], edges=edges_to_rm )

    #     nx.readwrite.graphml.write_graphml(grr_graph_filtered, os.path.join(output + f"{basename}_{edge_metric}-{threshold}.graphml"))


    # writting graph in gexf format
    logging.info(f"Writting graph in gexf and graphml format.")
    # nx.readwrite.gexf.write_gexf(grr_graph, os.path.join(output + f"{basename}.gexf"))
    nx.readwrite.graphml.write_graphml(grr_graph, os.path.join(output + f"{basename}.graphml"))



def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    
    mk_outdir(args.output, args.force)

    pangenome.add_file(args.pangenome)

    cluster_rgp(pangenome, args.grr_cutoff, args.output, args.basename, 
                args.cpu, args.ignore_incomplete_rgp, args.no_identical_rgp_merging,  args.disable_prog_bar)



def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for cluster_rgp command

    :return : parser arguments for cluster_rgp command
    """
    parser = sub_parser.add_parser("rgp_cluster", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_cluster_rgp(parser)
    return parser


def parser_cluster_rgp(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of rgp command

    :param parser: parser for cluster_rgp argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")


    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument('--grr_cutoff', required=False, type=restricted_float, default=0.5,
                          help="Min gene repertoire relatedness score used in the rgp clustering")
    
    optional.add_argument('--ignore_incomplete_rgp', required=False, action="store_true",
                          help="Do not cluster RGPs located on a contig border which are likely incomplete.")
    
    optional.add_argument('--no_identical_rgp_merging', required=False, action="store_true",
                          help="Do not merge in one node identical RGP (i.e. having the same family content) before clustering.")
    
    # optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")

    optional.add_argument("--basename", required=False, default="rgp_cluster", help="basename for the output file")

    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")