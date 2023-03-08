#!/usr/bin/env python3

import argparse
import logging
import os

from ppanggolin.utils import get_subcommand_parser, check_log

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

def split(list_object:list, chunk_count:int) -> list([list]):
    """
    Split list into n chunk. 

    :params list_object: list to split
    :params chunk_count: number of final chunk

    :return : list of chunk of the initial list. 
    """
    quotient, remainder = divmod(len(list_object), chunk_count)

    return [list_object[index*quotient+min(index, remainder):(index+1)*quotient+min(index+1, remainder)] for index in range(chunk_count)]


def split_comment_string(comment_string:str, max_word_count:int=20, prefix="\n    # "):
    """
    """
    splitted_comment = comment_string.split()
    word_count = len(splitted_comment)
    line_count = round(word_count/max_word_count) + 1

    comment_lines = [' '.join(words) for words in split(splitted_comment, line_count)]
    
    return prefix.join(comment_lines)

def get_input_argument_lines(argument_actions:list([argparse._SubParsersAction])) -> dict:
    """
    Manage input argument from a specific list of parser actions and format them for the yaml output.

    Input arguments are commented in the config file: as no default is valid.
    Help and possible values of the argument is added as comment line. 

    :param argument_actions: list of parser action for input arguments. 

    :return: default arguments for the given command 
    """

    arg_default_lines = []
    for action in argument_actions:

        # Add the help as comment
        arg_default_lines.append(f"    # {split_comment_string(action.help)}")

        arg_default_lines.append(f"  # {action.dest}: <{action.dest} file>")

    return arg_default_lines

def get_default_argument_lines(argument_actions:list([argparse._SubParsersAction])) -> dict:
    """
    Get default arguments for a specific list of parser actions and format them for the yaml output.

    Help and possible values of the argument is added as comment line. 

    :param argument_actions: list of parser action arguments. 

    :return: default arguments for the given command 
    """

    arg_default_lines = []
    for action in argument_actions:
        if action.help == "==SUPPRESS==":
            # arg with suppressed help are ignored.
            continue

        # Add the help as comment
        arg_default_lines.append(f"    # {action.help}")

        if action.choices:
            arg_default_lines.append(f"    # Choices: {', '.join([str(choice) for choice in action.choices])}")
        
        # When default is None, it is replaced by False to omit the arg and get the None value as expected.
        default = action.default if action.default is not None else False
        arg_default_lines.append(f"    {action.dest}: {default}")

    return arg_default_lines

def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if os.path.exists(args.output) and not args.force: 
        raise FileExistsError(f"{args.output} already exists. Use -f if you want to overwrite it.")

    
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
            "context":ppanggolin.context.subparser}
    
    input_params = ['fasta', 'anno', 'clusters', 'pangenome']
    general_params = ['output', 'basename', 'rarefaction', 'only_pangenome', 'tmpdir', 'verbose', 'log', 'disable_prog_bar', 'force']
    ignored_params = ['config', 'help']
    unspecific_params = input_params + general_params + ignored_params

    workflow_subcommands = ['all', 'workflow', 'panrgp', 'panmodule']
    workflow_dependencies = ["annotate", "cluster", "graph", "partition", "write"] #, "rgp", "spot", "module" ]

    write_flag_default_in_wf = ["csv", "Rtab", "gexf", "light_gexf",
                        'projection', 'stats', 'json', 'partitions']
    
    if args.command in ['panrgp', 'all']:
        workflow_dependencies += ["rgp", "spot"]
        write_flag_default_in_wf += ['regions', 'spots', 'borders']

    if args.command in ['panmodule', 'all']:
        workflow_dependencies.append('module')
        write_flag_default_in_wf.append('modules')

    if args.command == 'all':
        write_flag_default_in_wf.append('spot_modules')

    parser_fct = subcommand_to_subparser[args.command]

    _, sub = get_subcommand_parser(parser_fct, args.command)
    
    inputs_actions = []
    general_actions = []
    specific_actions = []

    for parser_action in sub._actions:
        if parser_action.dest in ignored_params:
            continue

        if parser_action.dest in input_params:
            inputs_actions.append(parser_action)

        elif parser_action.dest in general_params:
            general_actions.append(parser_action)

        else:
            specific_actions.append(parser_action)

    arg_lines = ['input_parameters:']
    arg_lines += get_input_argument_lines(inputs_actions)

    arg_lines.append('\ngeneral_parameters:')
    arg_lines += get_default_argument_lines(general_actions)

    if args.command not in workflow_subcommands:
        arg_lines.append(f"\n{args.command}:")
        arg_lines += get_default_argument_lines(specific_actions)
    else:
        for wf_subcmd in workflow_dependencies:
            _, sub = get_subcommand_parser(subcommand_to_subparser[wf_subcmd], wf_subcmd )
            specific_subcmd_actions = [sub_action for sub_action in sub._actions if sub_action.dest not in unspecific_params]

            # overwrite some default value for write cmd in a workflow context
            if wf_subcmd == 'write':
                for sub_action in specific_subcmd_actions:
                    if sub_action.dest in write_flag_default_in_wf:
                        sub_action.default = True

            arg_lines.append(f"\n{wf_subcmd}:")
            arg_lines += get_default_argument_lines(specific_subcmd_actions)

    logging.info(f'Writting default config in {args.output}')
    with open(args.output, 'w') as fl:
        fl.write('\n'.join(arg_lines) + '\n')
        


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

    required.add_argument('-c', '--command', required=True, type=str, help="The subcommand forwhich to generate a config file with default values", choices=subcommands)
    
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument('-o','--output', type=str, default='default_config.yaml',
        help='output config file with default parameters written in yaml.')
    optional.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    optional.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")

    optional.add_argument('-f', '--force', action="store_true",
                        help="Overwrite the given output file if it exists.")
        
if __name__ == '__main__':
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
                        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
                        formatter_class=argparse.RawTextHelpFormatter)

    parser_default_config(main_parser)
    launch(main_parser.parse_args())