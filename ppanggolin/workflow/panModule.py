#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import time
import argparse

# local libraries
from ppanggolin.utils import add_step_specific_args
from ppanggolin.workflow.all import launch_workflow


"""a global workflow that does everything in one go."""


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    launch_workflow(args, subparser, panrgp=False, panmodule=True)

def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("panmodule", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Input arguments", description="The possible input arguments :")
    required.add_argument('--fasta', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). "
                               "One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the gff filepath of its annotations (the gffs can be compressed). "
                               "One line per organism. This option can be used alone "
                               "IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument("--clusters", required=False, type=str,
                          help="a tab-separated file listing the cluster names, the gene IDs, "
                               "and optionally whether they are a fragment or not.")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")
    optional.add_argument("--rarefaction", required=False, action="store_true",
                          help="Use to compute the rarefaction curves (WARNING: can be time consuming)")
    optional.add_argument("--only_pangenome", required=False, action="store_true",
                          help="Only generate the HDF5 pangenome file")

    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    
    add_step_specific_args(optional)

    return parser
