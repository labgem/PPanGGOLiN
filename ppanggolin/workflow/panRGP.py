#!/usr/bin/env python3

# default libraries
import argparse

# local libraries
from ppanggolin.workflow.all import launch_workflow, add_workflow_args


"""a global workflow that does everything in one go."""


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    launch_workflow(args, panrgp=True, panmodule=False)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line
    :param sub_parser : sub_parser for all command
    :return : parser arguments for all command
    """
    parser = sub_parser.add_parser(
        "panrgp", formatter_class=argparse.RawTextHelpFormatter
    )

    add_workflow_args(parser)

    return parser
