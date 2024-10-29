#!/usr/bin/env python3

# default libraries
import argparse
import tables
import logging

from pathlib import Path
import yaml

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.readBinaries import read_info
from ppanggolin.metrics.fluidity import compute_genomes_fluidity, fam_fluidity


def check_already_computed_metric(
    pangenome: Pangenome,
    genomes_fluidity: bool = False,
    print_metric: bool = True,
    recompute: bool = False,
):
    """
    Check if one of the asked metrics is not already computed

    :param pangenome: pangenome object
    :param genomes_fluidity: Ask to compute genome fluidity
    :param print_metric: Print metrics if already computed
    :param recompute: Are metrics going to be recompute
    """
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        if genomes_fluidity and "genomes_fluidity" in info_group._v_attrs._f_list():
            logging.getLogger("PPanGGOLiN").warning(
                "Genome fluidity has been already computed. "
                "Use --force if you want to compute it again"
            )
            if print_metric and not recompute:
                print_computed_metric(info_group._v_attrs["genomes_fluidity"])
            return True
    return False


def compute_metrics(
    pangenome: Pangenome,
    genomes_fluidity: bool = False,
    families_fluidity: bool = False,
    disable_bar: bool = False,
) -> dict:
    """Compute the metrics

    :param pangenome: pangenome which will be used to compute the genomes' fluidity
    :param genomes_fluidity: Ask to compute genome fluidity
    :param families_fluidity: Ask to compute family fluidity
    :param disable_bar: Disable the progress bar

    :return: dictionary with all the metrics computed
    """

    metrics_dict = {}
    if genomes_fluidity:
        metrics_dict["genomes_fluidity"] = compute_genomes_fluidity(
            pangenome, disable_bar
        )
    if families_fluidity:
        metrics_dict["families_fluidity"] = fam_fluidity(pangenome, disable_bar)

    return metrics_dict


def write_metrics(
    pangenome: Pangenome, metrics_dict: dict, print_metrics: bool = False
):
    """
    Write the metrics computed in the pangenome

    :param pangenome: pangenome which will be used to compute the genomes' fluidity
    :param metrics_dict: dictionary with all the metrics computed
    :param print_metrics: enable print of metrics
    """
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        logging.getLogger("PPanGGOLiN").debug("H5f open")
        if "genomes_fluidity" in metrics_dict.keys():
            logging.getLogger("PPanGGOLiN").info(
                "Writing genome fluidity of the pangenome."
            )
            info_group._v_attrs.genomes_fluidity = metrics_dict["genomes_fluidity"]

        # After all metrics have been written
        if print_metrics:
            print_computed_metric(metrics_dict["genomes_fluidity"])


def print_computed_metric(metrics_dict: dict):
    """
    Print metrics in yaml format

    :params metrics_dict: Dict of computed metrics
    """
    metric_dict = {
        "Genomes_fluidity": {key: round(val, 3) for key, val in metrics_dict.items()}
    }
    metric_yaml = yaml.dump(
        metric_dict, default_flow_style=False, sort_keys=False, indent=4
    )
    print(metric_yaml)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if not any(x for x in [args.genome_fluidity]):
        raise Exception("You did not indicate which metric you want to compute.")

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    print_metrics = not args.no_print_info

    logging.getLogger("PPanGGOLiN").debug(
        "Check if one of the metrics was already computed"
    )
    is_metric_already_computed = check_already_computed_metric(
        pangenome,
        genomes_fluidity=args.genome_fluidity,
        print_metric=print_metrics,
        recompute=args.recompute_metrics,
    )

    if not is_metric_already_computed or args.recompute_metrics:
        logging.getLogger("PPanGGOLiN").info("Metrics computation begin")
        metrics_dictionary = compute_metrics(
            pangenome,
            disable_bar=args.disable_prog_bar,
            genomes_fluidity=args.genome_fluidity,
        )
        logging.getLogger("PPanGGOLiN").info("Metrics computation done")

        write_metrics(pangenome, metrics_dictionary, print_metrics=print_metrics)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "metrics", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_metrics(parser)
    return parser


def parser_metrics(parser: argparse.ArgumentParser):
    """
    Argument parser for the 'metrics' command.

    :param parser: Argument parser for the 'metrics' command.
    """
    required = parser.add_argument_group(
        title="Required arguments", description="Specify the required argument:"
    )
    required.add_argument(
        "-p",
        "--pangenome",
        required=False,
        type=Path,
        help="Path to the pangenome .h5 file",
    )

    onereq = parser.add_argument_group(
        title="Input file", description="Choose one of the following arguments:"
    )
    onereq.add_argument(
        "--genome_fluidity",
        required=False,
        action="store_true",
        default=False,
        help="Compute the pangenome genomic fluidity.",
    )

    optional = parser.add_argument_group(
        title="Optional arguments",
        description="Specify optional arguments with default values:",
    )
    optional.add_argument(
        "--no_print_info",
        required=False,
        action="store_true",
        default=False,
        help="Suppress printing the metrics result. "
        "Metrics are saved in the pangenome and viewable using 'ppanggolin info'.",
    )
    optional.add_argument(
        "--recompute_metrics",
        action="store_true",
        help="Force re-computation of metrics if already computed.",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_metrics(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
