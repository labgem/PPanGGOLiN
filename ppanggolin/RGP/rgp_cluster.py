#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
import time
import os
from itertools import combinations
from collections.abc import Callable

# installed libraries
from tqdm import tqdm
import networkx as nx

# local libraries
from ppanggolin.genome import Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import restricted_float, mk_outdir

def compute_grr(rgp_a:Region, rgp_b:Region, mode:Callable):
    """
    Compute gene repertoire relatedness (GRR) between two rgp.
    mode can be the function min to compute min GRR or max to compute max_grr

    :param rgp_a: rgp A
    :param rgp_b: rgp B
    :param mode: min or max function

    """
    
    max_grr = len((rgp_a.families & rgp_b.families))/mode(len(rgp_a.families), len(rgp_b.families))

    return max_grr



def cluster_rgp(pangenome, output, disable_bar):
    """
    Main function to cluster regions of genomic plasticity based on their GRR

    :param pangenome: pangenome object
    :param output: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """

    # check statuses and load info
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=False, need_partitions=True,
                        disable_bar=disable_bar, need_rgp=True)

    grr_graph = nx.Graph()

    # add all rgp as node
    grr_graph.add_nodes_from(pangenome.regions)

    # add some attribute to the graph. not really important..
    region_attributes = {}
    for region in pangenome.regions:
        
        
        region_info = {"contig":region.contig.name,
                       'organism':region.organism.name}
        
        region_info['famillies'] =';'.join([str(f.ID)  for f in region.families])

        region_attributes[region] = region_info


    nx.set_node_attributes(grr_graph, region_attributes)

    # get all possible pair of rgp
    rgp_pairs = combinations(pangenome.regions, 2)
    
    for rgp_a, rgp_b in rgp_pairs:
        min_grr = compute_grr(rgp_a, rgp_b, min)
        max_grr = compute_grr(rgp_a, rgp_b, max)

        # add edge if at least one of grr is > 0
        if min_grr or max_grr:
            grr_graph.add_edge(rgp_a, rgp_b, max_grr=max_grr, min_grr=min_grr)
    
    # cluster rgp based on grr
    min_grr_partitions = nx.algorithms.community.louvain_communities(grr_graph, weight='min_grr')
    max_grr_partitions = nx.algorithms.community.louvain_communities(grr_graph, weight='max_grr')
    
    # Add partition index in node attributes
    for i, part in enumerate(min_grr_partitions):
        for node in part:
            grr_graph.add_node(node, min_grr_cluster=i)

    for i, part in enumerate(max_grr_partitions):
        for node in part:
            grr_graph.add_node(node, max_grr_cluster=i)
        
    logging.info(f"GRR graph has {len(min_grr_partitions)} partitions using min_grr")
    logging.info(f"GRR graph has {len(max_grr_partitions)} partitions using max_grr")

    # writting graph in gexf format 
    nx.readwrite.gexf.write_gexf(grr_graph, output + "/grrGraph.gexf")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    
    mk_outdir(args.output, args.force)

    pangenome.add_file(args.pangenome)

    cluster_rgp(pangenome, args.output, args.disable_prog_bar)



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

    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")