#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import time
import os
from pathlib import Path

# Installed libraries

# local libraries
from ppanggolin.utils import mk_outdir, restricted_float
from ppanggolin.pangenome import Pangenome
from ppanggolin.figures.draw_spot import draw_spots
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.figures.ucurve import draw_ucurve


def check_spot_args(args: argparse.Namespace):
    """
    Check whether the draw_spots and spots arguments are valid.

    :param args: The parsed command line arguments.
    :type args: argparse.Namespace
    :raises argparse.ArgumentError: If args.spots is specified but args.draw_spots is False.
    """
    default_arg_spots = 'all'
    if not args.draw_spots and args.spots != default_arg_spots:
        raise argparse.ArgumentError(None, "The --spots argument cannot be used when --draw_spots is not specified.")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    check_spot_args(args)

    mk_outdir(args.output, args.force)

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.tile_plot:
        draw_tile_plot(pangenome, args.output, args.nocloud, disable_bar=args.disable_prog_bar)
    if args.ucurve:
        draw_ucurve(pangenome, args.output, soft_core=args.soft_core, disable_bar=args.disable_prog_bar)
    if args.draw_spots:
        draw_spots(pangenome=pangenome, output=args.output, spot_list=args.spots, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """

    parser = sub_parser.add_parser("draw", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_draw(parser)
    return parser


def parser_draw(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of draw command

    :param parser: parser for align argument
    """
    date = time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome.h5 file")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=Path,
                          default=Path(f'ppanggolin_output{date}_PID{str(os.getpid())}'),
                          help="Output directory")
    optional.add_argument("--tile_plot", required=False, default=False, action="store_true",
                          help="draw the tile plot of the pangenome")
    optional.add_argument("--nocloud", required=False, default=False, action="store_true",
                          help="Do not draw the cloud in the tile plot")
    optional.add_argument("--soft_core", required=False, default=0.95, type=restricted_float,
                          help="Soft core threshold to use")
    optional.add_argument("--ucurve", required=False, default=False, action="store_true",
                          help="draw the U-curve of the pangenome")
    optional.add_argument("--draw_spots", required=False, default=False, action="store_true",
                          help="draw plots for spots of the pangenome")
    optional.add_argument("--spots", required=False, default='all', nargs='+',
                          help="a comma-separated list of spots to draw (or 'all' to draw all spots, or 'synteny' to draw spots with different RGP syntenies).")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_draw(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
