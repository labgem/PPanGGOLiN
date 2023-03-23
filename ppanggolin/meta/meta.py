#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
from pathlib import Path
import csv

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.meta.metafamilies import metadata_to_families


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.assign == "families":
        logging.getLogger().debug("Begin gene families metadata assignment...")
        metadata_to_families(pangenome, args.metadata, args.source, args.omit, args.force, args.disable_bar)
    elif args.assign == "genes":
        raise NotImplementedError("Option not implemented yet !")
    elif args.assign == "genomes":
        raise NotImplementedError("Option not implemented yet !")
    elif args.assign == "RGPs":
        raise NotImplementedError("Option not implemented yet !")
    elif args.assign == "spots":
        raise NotImplementedError("Option not implemented yet !")
    elif args.assign == "modules":
        raise NotImplementedError("Option not implemented yet !")



def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("metadata", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_meta(parser)
    return parser


def parser_meta(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of graph command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-m', '--metadata', type=Path, nargs='?',
                          help='Gene families annotation in TSV file. See our github for more detail about format')
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    required.add_argument("-a", "--assign", required=True, type=str, nargs="?",
                          choices=["families", "genomes", "genes", "RGPs", "spots", "modules"],
                          help="Select to which pangenome element metadata will be assigned.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--omit", required=False, action="store_true",
                          help="Allow to pass if a key in metadata is not find in pangenome")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_meta(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
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
