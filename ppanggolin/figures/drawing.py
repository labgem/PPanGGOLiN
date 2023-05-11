#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import time
import os

# Installed libraries

# local libraries
from ppanggolin.utils import mk_outdir
from ppanggolin.pangenome import Pangenome
from ppanggolin.figures.draw_spot import draw_spots
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.figures.ucurve import draw_ucurve


# TODO check args
# def check_draw_args(args):
#     if ...:
#         raise argparse.ArgumentError("Incompatible arguments")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.tile_plot:
        draw_tile_plot(pangenome, args.output, args.nocloud, disable_bar=args.disable_prog_bar)
    if args.ucurve:
        draw_ucurve(pangenome, args.output, soft_core=args.soft_core, disable_bar=args.disable_prog_bar)
    if args.spots != '':
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

    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=str, help="The pangenome.h5 file")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--tile_plot", required=False, default=False, action="store_true",
                          help="draw the tile plot of the pan")
    optional.add_argument("--nocloud", required=False, default=False, action="store_true",
                          help="Do not draw the cloud in the tile plot")
    optional.add_argument("--soft_core", required=False, default=0.95, help="Soft core threshold to use")
    optional.add_argument("--ucurve", required=False, default=False, action="store_true",
                          help="draw the U-curve of the pangenome")
    optional.add_argument("--spots", required=False, type=str, default='',
                          help="a comma-separated list of spots to draw (or 'all' to draw all spots)")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_draw(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
