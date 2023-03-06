#!/usr/bin/env python3

import argparse

import ppanggolin.nem.rarefaction
import ppanggolin.graph
import ppanggolin.annotate
import ppanggolin.cluster
import ppanggolin.figures
import ppanggolin.formats
import ppanggolin.info
import ppanggolin.metrics
import ppanggolin.align
import ppanggolin.RGP
import ppanggolin.mod
import ppanggolin.context
import ppanggolin.workflow
import ppanggolin.utility


""" Utility scripts to help formating input files of PPanggolin."""


def get_default_argument_lines(parser_fct, 
                        general_params: list,
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
        if comment:
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
    general_params = ['help', 'fasta', 'clusters', 'anno', "output", 'basename']
    step_arg_lines = []
    for step_name, parser_fct in step_to_arg_parser.items():

        arg_lines = [f"{step_name}:"]
        arg_lines += get_default_argument_lines(parser_fct, general_params= general_params, comment = comment)

        step_arg_lines.append('\n'.join(arg_lines) +'\n')        

    with open(output, 'w') as fl:
        fl.write('\n'.join(step_arg_lines))


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    subcommand_to_subparser = {
            "annotate":ppanggolin.annotate.subparser,
            "cluster":ppanggolin.cluster.subparser,
            "graph":ppanggolin.graph.subparser,
            "partition":ppanggolin.nem.partition.subparser,
            "rarefaction":ppanggolin.nem.rarefaction.subparser,
            "workflow":ppanggolin.workflow.workflow.subparser,
            "panrgp":ppanggolin.workflow.panRGP.subparser,
            "panModule":ppanggolin.workflow.panModule.subparser,
            "all":ppanggolin.workflow.all.subparser,
            "draw":ppanggolin.figures.subparser,
            "write":ppanggolin.formats.writeFlat.subparser,
            "fasta":ppanggolin.formats.writeSequences.subparser,
            "msa":ppanggolin.formats.writeMSA.subparser,
            "metrics":ppanggolin.metrics.metrics.subparser,
            "align":ppanggolin.align.subparser,
            "rgp":ppanggolin.RGP.genomicIsland.subparser,
            "spot":ppanggolin.RGP.spot.subparser,
            "module":ppanggolin.mod.subparser,
            "context":ppanggolin.context.subparser,
            "info":ppanggolin.info.subparser,
            "utility":ppanggolin.utility.subparser}

    add_comment = not args.no_comment

    write_yaml_default_config(args.output, subcommand_to_subparser, add_comment)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for info command

    :return : parser arguments for info command
    """
    parser = sub_parser.add_parser("default_config", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_default_config(parser)
    return parser


def parser_default_config(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of info command

    :param parser: parser for info argument
    """
    subcommands = ['annotate', 'cluster', 'graph', 'partition', 'rarefaction', 'workflow', 'panrgp', 'panModule',
                  'all', 'draw', 'write', 'fasta', 'msa', 'metrics', 'align', 'rgp', 'spot', 'module', 'context']

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    # create the parser for the "default_config" command

    required.add_argument('--subcommand', required=True, type=str, help="The subcommand forwhich to generate a config file with default values", choices=subcommands)
    
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument('-o','--output', type=str, default='default_config.yaml',
        help='output config file with default parameters written in yaml.')
    
    optional.add_argument("--no_comment", action="store_true",
        help="Does not add help for each argument as a comment in the yaml file."), 

if __name__ == '__main__':
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
                        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
                        formatter_class=argparse.RawTextHelpFormatter)

    parser_default_config(main_parser)
    launch(main_parser.parse_args())