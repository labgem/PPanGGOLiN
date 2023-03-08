#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
import time
import os

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.genome import Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import restricted_float, mk_outdir

def cluster_rgp(pangenome, disable_bar):
    """
    Main function to cluster regions of genomic plasticity

    :param pangenome: pangenome object
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=False, need_partitions=True,
                        disable_bar=disable_bar, need_rgp=True)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    print(args.pangenome)
    
    mk_outdir(args.output, args.force)

    pangenome.add_file(args.pangenome)

    cluster_rgp(pangenome, args.disable_prog_bar)



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