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


def cmd_line() -> argparse.Namespace:
    """Manage the command line argument given by user

    :return: arguments given and readable by PPanGGOLiN
    """
    # need to manually write the description so that it's displayed into groups of subcommands ....
    desc = "\n"
    desc += (
        "All of the following subcommands have their own set of options. To see them for a given subcommand,"
        " use it with -h or --help, as such:\n"
    )
    desc += "  ppanggolin <subcommand> -h\n"
    desc += "\n"
    desc += "  Basic:\n"
    desc += "    all           Easy workflow to run all possible analysis\n"
    desc += "    workflow      Easy workflow to run a pangenome analysis in one go\n"
    desc += (
        "    panrgp        Easy workflow to run a pangenome analysis with genomic islands and spots of"
        " insertion detection\n"
    )
    desc += "    panmodule     Easy workflow to run a pangenome analysis with module prediction\n"
    desc += "  \n"
    desc += "  Expert:\n"
    desc += "    annotate      Annotate genomes\n"
    desc += "    cluster       Cluster genes into gene families\n"
    desc += "    graph         Create the pangenome graph\n"
    desc += "    partition     Partition the pangenome graph\n"
    desc += "    rarefaction   Compute the rarefaction curve of the pangenome\n"
    desc += "    metadata      Add metadata to elements in yout pangenome\n"
    desc += "  \n"
    desc += "  Output:\n"
    desc += "    draw              Draw figures representing the pangenome through different aspects\n"
    desc += "    write_pangenome   Writes 'flat' files that represent the pangenome and its elements for use with other software.\n"
    desc += "    write_genomes     Writes 'flat' files that represent the genomes along with their associated pangenome elements.\n"
    desc += "    write_metadata    Writes 'TSV' files that represent the metadata associated with elements of the pangenome.\n"
    desc += "    fasta             Writes fasta files for different elements of the pangenome.\n"
    desc += (
        "    info              Prints information about a given pangenome graph file.\n"
    )
    desc += "    metrics           Compute several metrics on a given pangenome.\n"
    desc += "  \n"
    desc += "  Regions of Genomic Plasticity:\n"
    desc += "    rgp           Predicts Regions of Genomic Plasticity in the genomes of your pangenome.\n"
    desc += "    spot          Predicts spots in your pangenome.\n"
    desc += "    module        Predicts functional modules in your pangenome.\n"
    desc += "    rgp_cluster   Cluster RGPs based on their gene families.\n"
    desc += "  \n"
    desc += "  Analysis using reference pangenomes:\n"
    desc += "    msa          Compute Multiple Sequence Alignments for pangenome gene families.\n"
    desc += (
        "    align        Aligns a genome or a set of proteins to the pangenome gene families and "
        "predicts information from it.\n"
    )
    desc += "    context      Local genomic context analysis.\n"
    desc += "    projection   Annotates external genomes with an existing pangenome.\n"
    desc += "  \n"
    desc += "  Utility command:\n"
    desc += "    utils      Helper side commands."

    parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog + pan_epilog + rgp_epilog + mod_epilog,
    )

    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + version
    )

    subparsers = parser.add_subparsers(
        metavar="", dest="subcommand", title="subcommands", description=desc
    )
    subparsers.required = True  # because python3 sent subcommands to hell apparently

    # print help if no subcommand is specified
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # manage command parser to use command arguments
    subs = []
    for sub_cmd, sub_fct in SUBCOMMAND_TO_SUBPARSER.items():
        sub = sub_fct(subparsers)
        # add options common to all subcommands
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
        subs.append(sub)

    # manage command without common arguments
    sub_info = ppanggolin.info.subparser(subparsers)
    sub_utils = ppanggolin.utility.utils.subparser(subparsers)
    subs += [sub_info, sub_utils]

    # print help if only the command is given
    for sub in subs:
        if len(sys.argv) == 2 and sub.prog.split()[1] == sys.argv[1]:
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
