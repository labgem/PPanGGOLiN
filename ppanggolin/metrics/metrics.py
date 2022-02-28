#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tables
import logging

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.readBinaries import checkPangenomeInfo, readInfo
from ppanggolin.formats.writeBinaries import writeInfoModules

# metrics libraries
from ppanggolin.metrics.fluidity import genomes_fluidity, fam_fluidity


def check_metric(pangenome, all=False, genome_fluidity=False, family_fluidity=False, force=False):
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        if genome_fluidity or all:
            if 'genome_fluidity' in info_group._v_attrs._f_list() and not force:
                raise Exception("Genome fluidity was already compute. "
                                "Please use -f option if you REALLY want to compute again")

        if family_fluidity or all:
            if 'family_fluidity' in info_group._v_attrs._f_list() and not force:
                raise Exception("Family fluidity was already compute. "
                                "Please use -f option if you REALLY want to compute again")



def compute_metrics(pangenome, all=False, genome_fluidity=False, family_fluidity=False, disable_bar=False):
    """Compute the metrics
    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :type pangenome: Pangenome
    :param all: compute all the metrics
    :type all: bool
    :param genome_fluidity: Ask to compute genome fluidity
    :type genome_fluidity: bool
    :param family_fluidity: Ask to compute family fluidity
    :type family_fluidity: bool
    :param disable_bar: Disable the progress bar
    :type disable_bar: bool


    :return: dictionary with all the metrics computed
    :rtype: dict
    """

    metrics_dict = {}
    if genome_fluidity or all:
        metrics_dict['genome_fluidity'] = genomes_fluidity(pangenome, disable_bar)
    if family_fluidity or all:
        metrics_dict['family_fluidity'] = fam_fluidity(pangenome, disable_bar)
    if info_modules:
        checkPangenomeInfo(pangenome, needFamilies=True, needModules=True)
        metrics_dict['info_modules'] = True
    return metrics_dict


def write_metrics(pangenome, metrics_dict, no_print_info=False):
    """Write the metrics computed in the pangenome
    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :type pangenome: Pangenome
    :param metrics_dict: dictionary with all the metrics computed
    :type metrics_dict: dict
    """
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        logging.getLogger().debug("H5f open")
        if 'genome_fluidity' in metrics_dict.keys():
            logging.getLogger().info("Writing genome fluidity in pangenome")
            info_group._v_attrs.genome_fluidity = metrics_dict['genome_fluidity']

        if 'family_fluidity' in metrics_dict.keys():
            logging.getLogger().info("Writing family fluidity in pangenome")
            info_group._v_attrs.family_fluidity = metrics_dict['family_fluidity']
           
        if 'info_modules' in metrics_dict.keys():
            writeInfoModules(pangenome, h5f)

        # After all metrics was written
        if not no_print_info:
            readInfo(h5f)

        readInfo(h5f)


def launch(args):
    if not any(x for x in [args.genome_fluidity, args.family_fluidity, args.info_modules, args.all]):
        raise Exception("You did not indicate which metric you want to compute.")
        
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)

    logging.getLogger().debug("Check if one of the metrics was already compute")
    check_metric(pangenome, all=args.all, genome_fluidity=args.genome_fluidity, family_fluidity=args.family_fluidity,
                 force=args.force)
    logging.getLogger().info("Metrics computation begin")
    metrics_dictionary = compute_metrics(pangenome, all=args.all, genome_fluidity=args.genome_fluidity, info_modules=args.info_modules,
                                         family_fluidity=args.family_fluidity, disable_bar=args.disable_prog_bar)
    logging.getLogger().info("Metrics computation done")

    write_metrics(pangenome, metrics_dictionary, no_print_info=args.no_print_info)


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
    onereq.add_argument('--genome_fluidity', required=False, action="store_true", default=False,
                        help="Compute the pangenome genomic fluidity.")
    # help="Compute the pangenome genomic and/or family fluidity.")
    onereq.add_argument('--info_modules', required=False, action='store_true', default=False,
                        help='Compute more information about modules')
    onereq.add_argument('--family_fluidity', required=False, action="store_true", default=False,
                        help=argparse.SUPPRESS)
    onereq.add_argument('--all', required=False, action="store_true", default=False,
                        help="Compute all the metrics")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument('--no_print_info', required=False, action="store_true", default=False,
                          help="Don't show the metrics result. "
                               "All the metric are saved in your pangenome and visible with ppanggolin info.")
    # optional.add_argument('--genome_only', required=False, action="store_true", default=False,
    #                       help="Compute the genome fluidity only")
    # optional.add_argument('--family_only', required=False, action="store_true", default=False,
    #                       help="Compute the genome fluidity only")
    return parser
