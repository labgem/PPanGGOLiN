#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tables
import logging
from pathlib import Path
# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.readBinaries import check_pangenome_info, read_info
from ppanggolin.formats.writeBinaries import write_info_modules

# metrics libraries
from ppanggolin.metrics.fluidity import gen_fluidity, fam_fluidity


def check_metric(pangenome: Pangenome, genomes_fluidity: bool = False, families_fluidity: bool = False,
                 info_modules: bool = False):
    """
    Check if one of the asked metrics is not already computed

    :param pangenome: pangenome object
    :param genomes_fluidity: Ask to compute genome fluidity
    :param families_fluidity: Ask to compute family fluidity
    :param info_modules: Ask to compute more information about module
    """
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        if genomes_fluidity:
            if 'genomes_fluidity' in info_group._v_attrs._f_list():
                raise Exception("Genome fluidity was already compute. "
                                "Please use -f option if you REALLY want to compute again")

        if families_fluidity:
            if 'families_fluidity' in info_group._v_attrs._f_list():
                raise Exception("Family fluidity was already compute. "
                                "Please use -f option if you REALLY want to compute again")

        if info_modules:
            if any(x in info_group._v_attrs._f_list() for x in ['CloudSpecInModules', 'PersistentSpecInModules',
                                                                'ShellSpecInModules', 'numberOfFamiliesInModules',
                                                                'StatOfFamiliesInModules']):
                raise Exception("Supplementary information on modules was already compute. "
                                "Please use -f option if you REALLY want to compute again")


def compute_metrics(pangenome: Pangenome, genomes_fluidity: bool = False, families_fluidity: bool = False,
                    info_modules: bool = False, disable_bar: bool = False) -> dict:
    """Compute the metrics

    :param pangenome: pangenome which will be used to compute the genomes' fluidity
    :param genomes_fluidity: Ask to compute genome fluidity
    :param families_fluidity: Ask to compute family fluidity
    :param info_modules: Ask to compute more information about module
    :param disable_bar: Disable the progress bar

    :return: dictionary with all the metrics computed
    """

    metrics_dict = {}
    if genomes_fluidity:
        metrics_dict['genomes_fluidity'] = gen_fluidity(pangenome, disable_bar)
    if families_fluidity:
        metrics_dict['families_fluidity'] = fam_fluidity(pangenome, disable_bar)
    if info_modules:
        check_pangenome_info(pangenome, need_families=True, need_modules=True)
        metrics_dict['info_modules'] = True
    return metrics_dict


def write_metrics(pangenome: Pangenome, metrics_dict: dict, no_print_info: bool = False):
    """
    Write the metrics computed in the pangenome

    :param pangenome: pangenome which will be used to compute the genomes' fluidity
    :param metrics_dict: dictionary with all the metrics computed
    :param no_print_info: disable print of information
    """
    with tables.open_file(pangenome.file, "a") as h5f:
        info_group = h5f.root.info
        logging.getLogger("PPanGGOLiN").debug("H5f open")
        if 'genomes_fluidity' in metrics_dict.keys():
            logging.getLogger("PPanGGOLiN").info("Writing genome fluidity in pangenome")
            info_group._v_attrs.genomes_fluidity = metrics_dict['genomes_fluidity']

        if 'info_modules' in metrics_dict.keys():
            logging.getLogger("PPanGGOLiN").info("Writing modules information in pangenome")
            write_info_modules(pangenome, h5f)

        # After all metrics was written
        if not no_print_info:
            read_info(h5f)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if not any(x for x in [args.genome_fluidity, args.info_modules, args.all]):
        raise Exception("You did not indicate which metric you want to compute.")
    args_dict = {'genomes_fluidity': args.genome_fluidity,
                 'info_modules': args.info_modules}
    if args.all:
        for arg in args_dict.keys():
            args_dict[arg] = True

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    logging.getLogger("PPanGGOLiN").debug("Check if one of the metrics was already compute")
    if not args.force:
        check_metric(pangenome, **args_dict)
    logging.getLogger("PPanGGOLiN").info("Metrics computation begin")
    metrics_dictionary = compute_metrics(pangenome, disable_bar=args.disable_prog_bar, **args_dict)
    logging.getLogger("PPanGGOLiN").info("Metrics computation done")

    write_metrics(pangenome, metrics_dictionary, no_print_info=args.no_print_info)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("metrics", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_metrics(parser)
    return parser


def parser_metrics(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of metrics command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome .h5 file")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('--genome_fluidity', required=False, action="store_true", default=False,
                        help="Compute the pangenome genomic fluidity.")
    # help="Compute the pangenome genomic and/or family fluidity.")
    onereq.add_argument('--info_modules', required=False, action='store_true', default=False,
                        help='Compute more information about modules')
    onereq.add_argument('--all', required=False, action="store_true", default=False,
                        help="Compute all the metrics")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument('--no_print_info', required=False, action="store_true", default=False,
                          help="Don't show the metrics result. "
                               "All the metric are saved in your pangenome and visible with ppanggolin info.")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_metrics(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
