#!/usr/bin/env python3

# default libraries
import argparse
import logging
import os
from pathlib import Path
from typing import List

# local libraries
from ppanggolin.utils import (
    get_subcommand_parser,
    check_log,
    ALL_INPUT_PARAMS,
    ALL_GENERAL_PARAMS,
    WORKFLOW_SUBCOMMANDS,
    ALL_WORKFLOW_DEPENDENCIES,
    WRITE_PAN_FLAG_DEFAULT_IN_WF,
    WRITE_GENOME_FLAG_DEFAULT_IN_WF,
    DRAW_FLAG_DEFAULT_IN_WF,
)
from ppanggolin import SUBCOMMAND_TO_SUBPARSER

""" Utility scripts to help formatting input files of PPanggolin."""


def split(list_object: list, chunk_count: int) -> List[List[int]]:
    """
    Split list into n chunk.

    :params list_object: list to split
    :params chunk_count: number of final chunk

    :return : list of chunk of the initial list.
    """
    quotient, remainder = divmod(len(list_object), chunk_count)

    return [
        list_object[
            index * quotient
            + min(index, remainder) : (index + 1) * quotient
            + min(index + 1, remainder)
        ]
        for index in range(chunk_count)
    ]


def split_comment_string(
    comment_string: str, max_word_count: int = 20, prefix: str = "\n    # "
) -> str:
    """
    Split a line of comment into multiple line.

    :params comment_string: comment string to split
    :params max_word_count: maximum number of word per line
    :params prefix: prefix used to start a new comment line

    :return : the split comment line.
    """

    splitted_comment = comment_string.split()
    word_count = len(splitted_comment)
    line_count = round(word_count / max_word_count) + 1

    comment_lines = [" ".join(words) for words in split(splitted_comment, line_count)]

    return prefix.join(comment_lines)


def get_input_argument_lines(
    argument_actions: List[argparse._SubParsersAction],
) -> List[str]:
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


def get_default_argument_lines(
    argument_actions: List[argparse._SubParsersAction],
) -> List[str]:
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
            arg_default_lines.append(
                f"    # Choices: {', '.join([str(choice) for choice in action.choices])}"
            )

        # When default is None, it is replaced by False to omit the arg and get the None value as expected.
        default = action.default if action.default is not None else False
        arg_default_lines.append(f"    {action.dest}: {default}")

    return arg_default_lines


def deduplicate_actions(
    actions: List[argparse._SubParsersAction],
) -> List[argparse._SubParsersAction]:
    """
    Deduplicate duplicate actions based on their dest.

    When two actions with the same dest attribute, only the first one is kept in the returned list.

    :param actions: list of parser action arguments.

    :return: list of parser action arguments deduplicated.

    """
    dedup_actions = []
    dedup_names = []
    for action in actions:
        if action.dest in dedup_names:
            # action has been already seen.
            continue
        dedup_names.append(action.dest)
        dedup_actions.append(action)
    return dedup_actions


def launch_default_config(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    initial_command = args.default_config

    if args.output.exists() and not args.force:
        raise FileExistsError(
            f"{args.output} already exists. Use -f if you want to overwrite it."
        )

    ignored_params = ["config", "help"]

    workflow_dependencies = {
        sub_cmd
        for sub_cmd in ALL_WORKFLOW_DEPENDENCIES
        if sub_cmd not in ["rgp", "spot", "module"]
    }

    if initial_command in ["panrgp", "all"]:
        workflow_dependencies |= {"rgp", "spot"}

    if initial_command in ["panmodule", "all"]:
        workflow_dependencies.add("module")

    if initial_command in WORKFLOW_SUBCOMMANDS:
        # it is clearer if the order of the subcommand is conserved in wf config file
        commands = [initial_command] + [
            sub_cmd
            for sub_cmd in ALL_WORKFLOW_DEPENDENCIES
            if sub_cmd in workflow_dependencies
        ]
    elif initial_command == "projection":
        commands = [initial_command] + ["annotate"]

    else:
        commands = [initial_command]

    sub_cmd_to_actions = {}
    inputs_actions = []
    general_actions = []

    for sub_command in commands:

        parser_fct = SUBCOMMAND_TO_SUBPARSER[sub_command]

        _, sub = get_subcommand_parser(parser_fct, sub_command)

        specific_actions = []

        # overwrite some default value for write cmd in a workflow context
        if initial_command in WORKFLOW_SUBCOMMANDS and sub_command in [
            "write_pangenome",
            "write_genomes",
        ]:
            for sub_action in sub._actions:
                if (
                    sub_action.dest
                    in WRITE_PAN_FLAG_DEFAULT_IN_WF + WRITE_GENOME_FLAG_DEFAULT_IN_WF
                ):
                    sub_action.default = True
        # overwrite some default value for draw cmd in a workflow context
        if initial_command in WORKFLOW_SUBCOMMANDS and sub_command == "draw":
            for sub_action in sub._actions:
                if sub_action.dest in DRAW_FLAG_DEFAULT_IN_WF:
                    sub_action.default = True

        for parser_action in sub._actions:
            if parser_action.dest in ignored_params:
                continue

            if parser_action.dest in ALL_INPUT_PARAMS:
                if sub_command == initial_command:
                    # with workflow dependencies, we do not use their input params
                    # as input params are given by workflow cmds
                    inputs_actions.append(parser_action)

            elif parser_action.dest in ALL_GENERAL_PARAMS:
                general_actions.append(parser_action)

            else:
                specific_actions.append(parser_action)

        sub_cmd_to_actions[sub_command] = specific_actions

    inputs_actions = deduplicate_actions(inputs_actions)
    general_actions = deduplicate_actions(general_actions)

    arg_lines = ["input_parameters:"]
    arg_lines += get_input_argument_lines(inputs_actions)

    arg_lines.append("\ngeneral_parameters:")
    arg_lines += get_default_argument_lines(general_actions)

    for sub_command, specific_actions in sub_cmd_to_actions.items():

        if sub_command in WORKFLOW_SUBCOMMANDS:
            continue

        arg_lines.append(f"\n{sub_command}:")
        arg_lines += get_default_argument_lines(specific_actions)

    logging.getLogger("PPanGGOLiN").info(f"Writing default config in {args.output}")
    with open(args.output, "w") as fl:
        fl.write("\n".join(arg_lines) + "\n")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    if args.default_config is not None:
        launch_default_config(args)

    # elif args.another_util_args is not None:
    #     launch_another_utils()


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for info command

    :return : parser arguments for info command
    """
    parser = sub_parser.add_parser(
        "utils", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_default_config(parser)
    return parser


def parser_default_config(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of utils command

    :param parser: parser for utils argument
    """

    subcommands = list(SUBCOMMAND_TO_SUBPARSER.keys())

    required = parser.add_argument_group(
        title="Required arguments",
        description="All of the following arguments are required :",
    )

    required.add_argument(
        "--default_config",
        required=False,
        type=str,
        default=None,  # nargs="*",,
        help="Generate a config file with default values for the given subcommand.",
        choices=subcommands,
    )

    optional = parser.add_argument_group(title="Config arguments")

    optional.add_argument(
        "-o",
        "--output",
        type=Path,
        default="default_config.yaml",
        help="name and path of the config file with default parameters written in yaml.",
    )

    optional.add_argument(
        "--verbose",
        required=False,
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)",
    )

    optional.add_argument(
        "--log",
        required=False,
        type=check_log,
        default="stdout",
        help="log output file",
    )

    optional.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite the given output file if it exists.",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_default_config(main_parser)

    launch(main_parser.parse_args())
