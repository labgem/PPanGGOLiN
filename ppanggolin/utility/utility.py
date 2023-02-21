#!/usr/bin/env python3

import argparse

from ppanggolin.annotate.annotate import parser_annot
from ppanggolin.cluster.cluster import parser_clust
from ppanggolin.graph.makeGraph import parser_graph
from ppanggolin.nem.rarefaction import parser_rarefaction
from ppanggolin.nem.partition import parser_partition
from ppanggolin.formats.writeFlat import  parser_flat
from ppanggolin.RGP.genomicIsland import parser_rgp
from ppanggolin.RGP.spot import parser_spot
from ppanggolin.mod.module import parser_module



""" Utility scripts to help formating input files of PPanggolin."""


def get_default_argument_lines(parser_fct, 
                        general_params: list = ['help', 'fasta', 'clusters', 'anno', 'cpu', "output"],
                        comment: bool = True, 
                        indentation: str ='    ') -> dict:
    """
    Get default arguments for a specific command using the argument parser of this command and format them for the yaml output.

    Help and possible values of the argument is added as comment line. 

    :param parser_fct: parser function of the command.
    :param general_params: General parameters to remove from the expected arguments.
    :param indentation: Indentation to use for starting a line.

    :return: default arguments for the given command 
    """

    parser = argparse.ArgumentParser()  

    parser_fct(parser)

    # remove required arguments. Config file ca only contained optional arguments
    parser._actions = [p_action for p_action in parser._actions if p_action.required == False]

    # remove general arguments to only expect arguments specific to the step.
    parser._actions = [p_action for p_action in parser._actions if p_action.dest not in general_params]

    arg_default_lines = []
    for action in parser._actions:
        if action.help == "==SUPPRESS==":
            # arg with suppressed help are ignored.
            continue

        # Add the help as comment
        arg_default_lines.append(f"{indentation}# {action.help}")

        if action.choices:
            arg_default_lines.append(f"{indentation}# Choices: {', '.join(action.choices)}")

        # When default is None, it is replaced by False to omit the arg and get the None value as expected.
        default = action.default if action.default is not None else False
        arg_default_lines.append(f"{indentation}{action.dest}: {default}")

    return arg_default_lines


def write_yaml_default_config(output: str, step_to_arg_parser: dict, comment:bool):
    """
    Write default config in yaml format.


    :param output: output file name. 
    :param step_to_arg_parser: key = name of the step and value = specific parser function for the step.
    :param comment: if true, help comments describing the arguments are not added to the yaml file.
    """
    
    step_arg_lines = []
    for step_name, parser_fct in step_to_arg_parser.items():

        arg_lines = [f"{step_name}:"]
        arg_lines += get_default_argument_lines(parser_fct, comment = comment)

        step_arg_lines.append('\n'.join(arg_lines) +'\n')        

    with open(output, 'w') as fl:
        fl.write('\n'.join(step_arg_lines))


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if args.utility == "default_config":
        print(args.output)

        step_to_arg_parser = {
            "annotate":parser_annot, 
            "cluster": parser_clust,
            "graph": parser_graph,
            "partition": parser_partition,
            "rarefaction": parser_rarefaction,
            "write": parser_flat,
            "rgp":parser_rgp,
            "spot":parser_spot, 
            "module":parser_module
        }

        add_comment = not args.no_comment

        write_yaml_default_config(args.output, step_to_arg_parser, add_comment)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for info command

    :return : parser arguments for info command
    """
    parser = sub_parser.add_parser("utility", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_utility(parser)
    return parser


def parser_utility(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of info command

    :param parser: parser for info argument
    """

    subparsers = parser.add_subparsers(help='help for subcommand', dest="utility")

    # create the parser for the "default_config" command
    parser_config = subparsers.add_parser('default_config', help='Generate default config file to use as a template in workflow commands.')
    parser_config.add_argument('-o','--output', type=str, default='default_config.yaml',
                         help='output config file with default parameters written in yaml.')
    parser_config.add_argument("--no_comment", required=False, action="store_true",
                                help="Does not add help for each argument as a comment in the yaml file.")

    # create the parser for the "command_2" command
    # parser_b = subparsers.add_parser('command_2', help='help for command_2')
    # parser_b.add_argument('-b', type=str, help='help for b')
    # parser_b.add_argument('-c', type=str, action='store', default='', help='test')


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_utility(main_parser)
    launch(main_parser.parse_args())