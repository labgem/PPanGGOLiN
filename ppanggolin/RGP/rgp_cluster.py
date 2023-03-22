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
from multiprocessing import get_context
from itertools import islice


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

        region_attributes[region.name] = region_info

    return region_attributes




def manage_identical_rgps(rgp_graph, uniq_rgps, rgp_to_spot):
    """
    """
    

    for rgp, identical_rgps in uniq_rgps.items():

        spot_to_rgps = defaultdict(set)
        for identical_rgp in identical_rgps:
            spot = rgp_to_spot.get(identical_rgp, None)
            spot_to_rgps[spot].add(identical_rgp)

        for spot, strictly_identical_rgps in spot_to_rgps.items():

            if spot == rgp_to_spot.get(rgp, None): # is spot same as the main rgp that have been used to compute grr?
                rgp_graph.add_node(rgp.name, strict_identical_rgp=len(strictly_identical_rgps))
            else:
                strictly_identical_rgp = strictly_identical_rgps.pop()
                rgp_graph.add_node(strictly_identical_rgp.name, strict_identical_rgp=len(strictly_identical_rgps), identical_rgp=True)
                rgp_graph.add_edge(rgp.name, identical_rgp.name, max_grr=1.0, min_grr=1.0, jaccard_index=1.0, identical_famillies=True)
                
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

def compute_rgp_metric(pairs_of_rgps):
    """
    """
    logging.debug(f'in compute_rgp_metric: computing metrics for {len(pairs_of_rgps)=}')
    pairs_of_rgps_metrics = []
    for (rgp_a, rgp_a_fam), (rgp_b, rgp_b_fam) in pairs_of_rgps:
        # compute metrics between 2 rgp if they share at least one familly
        if rgp_a_fam & rgp_b_fam:
            edge_metrics = {}
            edge_metrics['min_grr'] = compute_grr(rgp_a_fam, rgp_b_fam, min)
            edge_metrics['max_grr']  = compute_grr(rgp_a_fam, rgp_b_fam, max)
            edge_metrics['jaccard_index']  = compute_jaccard_index(rgp_a_fam, rgp_b_fam)
            edge_metrics['shared_family']  = len(rgp_a_fam & rgp_b_fam)

            pairs_of_rgps_metrics.append((rgp_a, rgp_b, edge_metrics))

    return pairs_of_rgps_metrics

def simplify_rgp_object(rgp):
    return rgp.name, {f.ID for f in rgp.families}


def make_chunks(iterable, size):
    """
    """

    i = iter(iterable)
    piece = [(simplify_rgp_object(rgp1),simplify_rgp_object(rgp2)) for rgp1, rgp2 in islice(i, size)]
    while piece:
        yield piece
        piece =  [(simplify_rgp_object(rgp1),simplify_rgp_object(rgp2)) for rgp1, rgp2 in islice(i, size)]

def cluster_nodes(G, clustering_attributes = ["min_grr", "max_grr", "jaccard_index"], prefix=""):
            
        for weight in clustering_attributes:
            partitions = nx.algorithms.community.louvain_communities(G, weight=weight)
        
            # Add partition index in node attributes
            for i, cluster_nodes in enumerate(partitions):
                nx.set_node_attributes(G, {node:f"cluster_{i}" for node in cluster_nodes}, name=f"{prefix}_{weight}_cluster")

            logging.info(f"{prefix}: Graph has {len(partitions)} clusters using {weight}")

def cluster_rgp(pangenome, output, basename, cpu, disable_bar):
    """
    Main function to cluster regions of genomic plasticity based on their GRR

    :param pangenome: pangenome object
    :param output: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """

    # check statuses and load info
    check_pangenome_info(pangenome, need_families=True, need_annotations=True,
                        disable_bar=disable_bar, need_rgp=True, need_spots=True)

    grr_graph = nx.Graph()

    # add all rgp as node

    uniq_rgps = dereplicate_rgp(pangenome.regions, disable_bar=disable_bar )

    grr_graph.add_nodes_from((rgp.name for rgp in uniq_rgps), identical_rgp=False)

    # compute grr for all possible pair of rgp
    rgp_count = len(uniq_rgps)
    pairs_count = int((rgp_count**2 - rgp_count)/2)
    logging.info(f'Computing GRR metric for {pairs_count:,} pairs of RGP using {cpu} cpus...')
    
    rgp_pairs = combinations(uniq_rgps, 2)
    chunk_size = int(pairs_count/(cpu*200)) +1
    logging.debug(f'Spliting RGP pairs in {pairs_count/(cpu*200)} chunks of {chunk_size} pairs')
    chunks_of_rgp_pairs = make_chunks(rgp_pairs, chunk_size)

    # for rgp_a, rgp_b in tqdm(combinations(uniq_rgps, 2), total=pairs_count, unit="RGP pairs", disable=disable_bar) :
    with get_context('fork').Pool(processes=cpu) as p:
        for pairs_of_rgps_metrics in tqdm(p.imap_unordered(compute_rgp_metric, chunks_of_rgp_pairs), unit=f"pairs of RGPs", unit_scale=round(pairs_count/(cpu*10)),
                             total=cpu*10, disable=disable_bar):
            logging.debug(f'adding {len(pairs_of_rgps_metrics)} edges to the graph...')
            grr_graph.add_edges_from(pairs_of_rgps_metrics)

        p.close()
        p.join()
    

    # cluster rgp based on grr
    logging.info(f"Louvain_communities clustering of RGP.")
    clustering_attributes = ["min_grr", "max_grr", "jaccard_index"]
    cluster_nodes(grr_graph, clustering_attributes, prefix="unfiltered") 


    logging.info(f"Clustering RGPs on filtered graph")
    thresolds = [0.1, 0.2, 0.5, 0.8]
    for edge_metric, threshold in product(clustering_attributes, thresolds):
        edges_to_rm = [(u,v) for u,v,e in grr_graph.edges(data=True) if e[edge_metric] < threshold]
        grr_graph_filtered = nx.restricted_view(grr_graph, nodes=[], edges=edges_to_rm )
        logging.info(f"Filtering graph edges with {edge_metric}<{threshold}: {grr_graph_filtered}")

        # cluster filtered graph
        cluster_nodes(grr_graph_filtered, [edge_metric], prefix=f"{edge_metric}>{threshold}")


    logging.info(f"Manage identical RGP in the graph")
    rgp_to_spot =  {region:spot.ID  for spot in pangenome.spots  for region in spot.regions}
    manage_identical_rgps(grr_graph, uniq_rgps, rgp_to_spot)
    
    # add some attribute to the graph nodes.
    logging.info(f"Add RGP information to the graph")
    region_infos = get_rgp_info_dict(pangenome.regions, rgp_to_spot)
    nx.set_node_attributes(grr_graph, region_infos)
    
    logging.info(f"Writting graph in graphml format with multiple thresholds.")
    thresolds = [0.1, 0.2, 0.5, 0.8]
    for edge_metric, threshold in product(clustering_attributes, thresolds):
        edges_to_rm = [(u,v) for u,v,e in grr_graph.edges(data=True) if e[edge_metric] < threshold]
        grr_graph_filtered = nx.restricted_view(grr_graph, nodes=[], edges=edges_to_rm )

        nx.readwrite.graphml.write_graphml(grr_graph_filtered, os.path.join(output + f"{basename}_{edge_metric}-{threshold}.graphml"))


    # writting graph in gexf format
    logging.info(f"Writting graph in gexf and graphml format.")
    nx.readwrite.gexf.write_gexf(grr_graph, os.path.join(output + f"{basename}.gexf"))
    nx.readwrite.graphml.write_graphml(grr_graph, os.path.join(output + f"{basename}.graphml"))


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    
    mk_outdir(args.output, args.force)

    pangenome.add_file(args.pangenome)

    cluster_rgp(pangenome, args.output, args.basename, args.cpu, args.disable_prog_bar)



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

    optional.add_argument('--grr_cutoff', required=False, type=restricted_float, default=0.8,
                          help="Gene repertoire relatedness score cutoff used in the rgp clustering")
    
    # optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")

    optional.add_argument("--basename", required=False, default="rgp_cluster", help="basename for the output file")

    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")