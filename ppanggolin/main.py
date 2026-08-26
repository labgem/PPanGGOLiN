#!/usr/bin/env python3

# default libraries
import sys

if sys.version_info < (3, 9):  # minimum is python3.9
    raise AssertionError(
        "Minimum python version to run PPanGGOLiN is 3.9. Your current python version is "
        + ".".join(map(str, sys.version_info))
    )
import argparse

# local modules
import ppanggolin.pangenome
from ppanggolin.utils import (
    RawTextHelpFormatterWithDefaults,
    check_input_files,
    set_verbosity_level,
    add_common_arguments,
    manage_cli_and_config_args,
)
import ppanggolin.nem.partition
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
import ppanggolin.meta
import ppanggolin.utility

from ppanggolin import (
    SUBCOMMAND_TO_SUBPARSER,
    epilog,
    pan_epilog,
    rgp_epilog,
    mod_epilog,
    version,
)


def _build_command_summary(parser: argparse.ArgumentParser) -> str:
    """Build the human-readable command list from each subparser's own metadata."""
    category_order = [
        "Basic",
        "Expert",
        "Output",
        "Regions of Genomic Plasticity",
        "Analysis using reference pangenomes",
        "Utility command",
    ]

    subparser_action = next(
        (
            action
            for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        ),
        None,
    )
    if subparser_action is None:
        return ""

    grouped = {category: [] for category in category_order}
    for command_name, command_parser in subparser_action.choices.items():
        category = getattr(command_parser, "category", "General")
        description = (getattr(command_parser, "description", "") or "").strip()
        grouped.setdefault(category, []).append((command_name, description))

    lines = [
        "All of the following subcommands have their own set of options. To see them for a given subcommand,",
        " use it with -h or --help, as such:",
        "  ppanggolin <subcommand> -h",
        "",
    ]

    for category in category_order:
        if category not in grouped or not grouped[category]:
            continue
        lines.append(f"  {category}:")
        for command_name, description in sorted(
            grouped[category], key=lambda item: item[0]
        ):
            lines.append(f"    {command_name:<15} {description}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def build_parser() -> argparse.ArgumentParser:
    """Build the PPanGGOLiN command-line parser without parsing user input."""

    argparse.RawTextHelpFormatter = RawTextHelpFormatterWithDefaults

    parser = argparse.ArgumentParser(
        prog="ppanggolin",
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=RawTextHelpFormatterWithDefaults,
        epilog=epilog + pan_epilog + rgp_epilog + mod_epilog,
    )

    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + version
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", title="subcommands", metavar=""
    )
    subparsers.required = True  # because python3 sent subcommands to hell apparently

    for sub_cmd, sub_fct in SUBCOMMAND_TO_SUBPARSER.items():
        sub = sub_fct(subparsers)
        add_common_arguments(sub)
        sub.epilog = epilog
        if sub_cmd not in ["rgp", "spot", "module", "rgp_cluster"]:
            sub.epilog += pan_epilog
        if sub_cmd not in [
            "annotate",
            "cluster",
            "graph",
            "partition",
            "rarefaction",
            "workflow",
        ]:
            if sub_cmd not in ["module", "panmodule"]:
                sub.epilog += rgp_epilog
            if sub_cmd not in ["rgp", "spot", "rgp_cluster", "panrgp"]:
                sub.epilog += mod_epilog

        # if getattr(sub, "description", None):
        #     sub.help = sub.description

    parser.epilog = _build_command_summary(parser) + parser.epilog
    return parser


def cmd_line() -> argparse.Namespace:
    """Manage the command line argument given by user

    :return: arguments given and readable by PPanGGOLiN
    """
    parser = build_parser()

    # print help if no subcommand is specified
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # print help if only the command is given
    subcommands = parser._subparsers._group_actions[0].choices
    for sub_name, sub in subcommands.items():
        if len(sys.argv) == 2 and sub_name == sys.argv[1]:
            sub.print_help()
            exit(0)

    # First parse args to check that nothing is missing or not expected in cli and throw help when requested
    args = parser.parse_args()
    if hasattr(args, "config"):
        # the two subcommand with no common args does not have config parameter. so we can skip this part for them.
        args = manage_cli_and_config_args(
            args.subcommand, args.config, SUBCOMMAND_TO_SUBPARSER
        )
    else:
        set_verbosity_level(args)

    if args.subcommand == "annotate" and args.fasta is None and args.anno is None:
        parser.error(
            "Please provide either a sequence file using the --fasta option or "
            "an annotation file using the --anno option to enable annotation. "
            "Use the command line or the config file."
        )

    cmds_pangenome_required = [
        "cluster",
        "info",
        "module",
        "graph",
        "align",
        "context",
        "write_pangenome",
        "write_genomes",
        "write_metadata",
        "msa",
        "draw",
        "partition",
        "rarefaction",
        "spot",
        "fasta",
        "metrics",
        "rgp",
        "projection",
        "metadata",
    ]
    if args.subcommand in cmds_pangenome_required and args.pangenome is None:
        parser.error(
            "Please specify a pangenome file using the --pangenome argument, "
            "either through the command line or the config file."
        )

    if args.subcommand == "align" and args.sequences is None:
        parser.error(
            "Please provide sequences (nucleotides or amino acids) for alignment "
            "with the pangenome gene families using the --sequences argument, "
            "either through the command line or the config file."
        )

    if args.subcommand == "projection":
        # check argument correctness and determine input mode (single or multiple files) and add it to args.
        input_mode = ppanggolin.projection.projection.check_projection_arguments(
            args, parser
        )
        setattr(args, "input_mode", input_mode)

    if args.subcommand == "metadata":
        # check argument correctness and determine input mode (single or multiple files) and add it to args.
        input_mode = ppanggolin.meta.meta.check_metadata_arguments(args, parser)
        setattr(args, "input_mode", input_mode)

    return args


def main():
    """
    Run the command given by user and set / check some things.

    :return:
    """
    args = cmd_line()

    if hasattr(args, "pangenome") and args.pangenome is not None:
        check_input_files(args.pangenome)

    if args.subcommand == "annotate":
        ppanggolin.annotate.launch(args)
    elif args.subcommand == "cluster":
        ppanggolin.cluster.launch(args)
    elif args.subcommand == "graph":
        ppanggolin.graph.launch(args)
    elif args.subcommand == "partition":
        ppanggolin.nem.partition.launch(args)
    elif args.subcommand == "workflow":
        ppanggolin.workflow.workflow.launch(args)
    elif args.subcommand == "rarefaction":
        ppanggolin.nem.rarefaction.launch(args)
    elif args.subcommand == "draw":
        ppanggolin.figures.launch(args)
    elif args.subcommand == "write_pangenome":
        ppanggolin.formats.writeFlatPangenome.launch(args)
    elif args.subcommand == "write_genomes":
        ppanggolin.formats.writeFlatGenomes.launch(args)
    elif args.subcommand == "write_metadata":
        ppanggolin.formats.writeFlatMetadata.launch(args)
    elif args.subcommand == "fasta":
        ppanggolin.formats.writeSequences.launch(args)
    elif args.subcommand == "msa":
        ppanggolin.formats.writeMSA.launch(args)
    elif args.subcommand == "info":
        ppanggolin.info.launch(args)
    elif args.subcommand == "metrics":
        ppanggolin.metrics.metrics.launch(args)
    elif args.subcommand == "align":
        ppanggolin.align.launch(args)
    elif args.subcommand == "projection":
        ppanggolin.projection.projection.launch(args)
    elif args.subcommand == "rgp":
        ppanggolin.RGP.genomicIsland.launch(args)
    elif args.subcommand == "spot":
        ppanggolin.RGP.spot.launch(args)
    elif args.subcommand == "rgp_cluster":
        ppanggolin.RGP.rgp_cluster.launch(args)
    elif args.subcommand == "panrgp":
        ppanggolin.workflow.panRGP.launch(args)
    elif args.subcommand == "module":
        ppanggolin.mod.launch(args)
    elif args.subcommand == "panmodule":
        ppanggolin.workflow.panModule.launch(args)
    elif args.subcommand == "all":
        ppanggolin.workflow.all.launch(args)
    elif args.subcommand == "context":
        ppanggolin.context.launch(args)
    elif args.subcommand == "metadata":
        ppanggolin.meta.launch(args)
    elif args.subcommand == "utils":
        ppanggolin.utility.launch(args)


if __name__ == "__main__":
    main()
