#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tables
import logging

# local libraries
from ppanggolin.pangenome import Pangenome

# metrics libraries
from ppanggolin.metrics.fluidity import genomes_fluidity


def compute_metrics(pangenome, fluidity, disable_bar=False):
    """Compute the metrics
    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :type pangenome: Pangenome
    :param fluidity: Ask to compute fluidity
    :type fluidity: bool
    :param disable_bar: Disable the progress bar
    :type disable_bar: bool


    :return: dictionary with all the metrics computed
    :rtype: dict
    """

    metrics_dict = {}
    if fluidity:
        metrics_dict['fluidity'] = genomes_fluidity(pangenome, disable_bar)

    return metrics_dict


def write_metrics(pangenome, metrics_dict):
    """Write the metrics computed in the pangenome
    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :type pangenome: Pangenome
    :param metrics_dict: dictionary with all the metrics computed
    :type metrics_dict: dict
    """
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        logging.getLogger().debug("Fluidity computation")
        if 'fluidity' in metrics_dict.keys():
            info_group._v_attrs.fluidity = metrics_dict['fluidity']


def launch(args):
    if not any(x for x in [args.fluidity]):
        raise Exception("You did not indicate which metric you want to compute.")
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)

    logging.getLogger().info("Metrics computation begin")
    metrics_dictionary = compute_metrics(pangenome, fluidity=args.fluidity, disable_bar=args.disable_prog_bar)
    logging.getLogger().info("Metrics computation done")

    write_metrics(pangenome, metrics_dictionary)


def metricsSubparser(sub_parser):
    """
    Parser arguments specific to metrics command

    :param sub_parser : sub_parser for align command
    :type sub_parser : argparse._SubParsersAction

    :return : parser arguments for align command
    :rtype : argparse.ArgumentParser
    """
    parser = sub_parser.add_parser("metrics", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('--fluidity', required=False, action="store_true",
                        help="Compute the pangenome genomic fluidity")
    return parser
