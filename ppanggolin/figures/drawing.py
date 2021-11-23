#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import time
import os

# local libraries
from ppanggolin.utils import mkOutdir
from ppanggolin.pangenome import Pangenome
from ppanggolin.figures.draw_spot import drawSpots
from ppanggolin.figures.tile_plot import drawTilePlot
from ppanggolin.figures.ucurve import drawUCurve


def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.tile_plot:
        drawTilePlot(pangenome, args.output, args.nocloud, disable_bar=args.disable_prog_bar)
    if args.ucurve:
        drawUCurve(pangenome, args.output, soft_core=args.soft_core, disable_bar=args.disable_prog_bar)
    if args.spots != '':
        drawSpots(pangenome=pangenome, output=args.output, spot_list=args.spots, disable_bar=args.disable_prog_bar)


def figureSubparser(subparser):
    parser = subparser.add_parser("draw", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--tile_plot", required=False, default=False, action="store_true",
                          help="draw the tile plot of the pangenome")
    optional.add_argument("--nocloud", required=False, default=False, action="store_true",
                          help="Do not draw the cloud in the tile plot")
    optional.add_argument("--soft_core", required=False, default=0.95, help="Soft core threshold to use")
    optional.add_argument("--ucurve", required=False, default=False, action="store_true",
                          help="draw the U-curve of the pangenome")
    optional.add_argument("--spots", required=False, type=str, default='',
                          help="a comma-separated list of spots to draw (or 'all' to draw all spots)")

    return parser
