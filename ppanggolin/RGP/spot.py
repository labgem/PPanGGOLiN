#!/usr/bin/env python3


# default libraries
import logging
import argparse
import time
import os

# installed libraries
import networkx as nx

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Spot
from ppanggolin.formats import checkPangenomeInfo, writePangenome, ErasePangenome
from ppanggolin.utils import mkOutdir


def compBorder(border1, border2, overlapping_match, exact_match, set_size):
    if border1[0:exact_match] == border2[0:exact_match]:
        return True
    elif len(border1) == set_size and len(border2) == set_size:
        for ikb in range(1, set_size - overlapping_match + 1):
            if border1[0:len(border2[ikb:])] == border2[ikb:]:
                return True
        for ib in range(1, set_size - overlapping_match + 1):
            if border1[ib:] == border2[0:len(border1[ib:])]:
                return True
    return False


def checkSim(pairBorder1, pairBorder2, overlapping_match, exact_match, set_size):
    """
    Checks if the two pairs of 'exact_match' first gene families are identical,
    or eventually if they overlap in an ordered way at least 'overlapping_match'
    """
    b1pair = [False, False]
    b2pair = [False, False]
    for countb1, b1 in enumerate(pairBorder1):
        for countb2, b2 in enumerate(pairBorder2):
            if compBorder(b2, b1, overlapping_match, exact_match, set_size):
                b1pair[countb1] = True
                b2pair[countb2] = True

    if b1pair[0] and b1pair[1] and b2pair[0] and b2pair[1]:
        return True
    return False


def makeSpotGraph(rgps, multigenics, output, spot_graph=False, overlapping_match=2, set_size=3, exact_match=1):
    def addNewNode(g, rgp, borders):
        blocks = str(sorted([[gene.family.ID for gene in borders[0]], [gene.family.ID for gene in borders[1]]],
                            key=lambda x: x[0]))
        g.add_node(blocks)
        try:
            g.nodes[blocks]["nb_rgp"] += 1
            g.nodes[blocks]["rgp"].add(rgp)
        except KeyError:
            g.nodes[blocks]["nb_rgp"] = 1
            g.nodes[blocks]["border1"] = [gene.family for gene in borders[1]]
            g.nodes[blocks]["border0"] = [gene.family for gene in borders[0]]
            g.nodes[blocks]["rgp"] = {rgp}

    spotGraph = nx.Graph()
    lost = 0
    used = 0
    for rgp in rgps:
        border = rgp.getBorderingGenes(set_size, multigenics)
        if len(border[0]) < set_size or len(border[1]) < set_size:
            lost += 1
        else:
            used += 1
            addNewNode(spotGraph, rgp, border)
    logging.getLogger().info(f"{lost} RGPs were not used as they are on a contig border (or have less than {set_size} "
                             f"persistent gene families until the contig border)")
    logging.getLogger().info(f"{used} RGPs are being used to predict spots of insertion")
    nodeList = list(spotGraph.nodes)
    logging.getLogger().info(f"{len(nodeList)} number of different pairs of flanking gene families")
    for i, nodei in enumerate(nodeList[:-1]):
        for nodej in nodeList[i + 1:]:
            nodeObji = spotGraph.nodes[nodei]
            nodeObjj = spotGraph.nodes[nodej]
            if checkSim([nodeObji["border0"], nodeObji["border1"]], [nodeObjj["border0"], nodeObjj["border1"]],
                        overlapping_match, exact_match, set_size):
                spotGraph.add_edge(nodei, nodej)

    spots = []
    spot_id = 0
    for comp in nx.algorithms.components.connected_components(spotGraph):
        curr_spot = Spot(spot_id)
        spots.append(curr_spot)
        for node in comp:
            curr_spot.addRegions(spotGraph.nodes[node]["rgp"])
        spot_id += 1

    if spot_graph:
        for node in spotGraph.nodes:
            del spotGraph.nodes[node]["border0"]
            del spotGraph.nodes[node]["border1"]
            del spotGraph.nodes[node]["rgp"]

        nx.readwrite.gexf.write_gexf(spotGraph, output + "/spotGraph.gexf")
    return spots


def checkPangenomeFormerSpots(pangenome, force):
    """ checks pangenome status and .h5 files for former spots, delete them if allowed or raise an error """
    if pangenome.status["spots"] == "inFile" and not force:
        raise Exception("You are trying to detect spots on a pangenome which already has predicted spots. "
                        "If you REALLY want to do that, use --force (it will erase spots previously predicted).")
    elif pangenome.status["spots"] == "inFile" and force:
        ErasePangenome(pangenome, spots=True)


def predictHotspots(pangenome, output, force=False, cpu=1, spot_graph=False, overlapping_match=2, set_size=3,
                    exact_match=1, disable_bar=False):
    # check that given parameters for hotspot computation make sense
    checkParameterLogic(overlapping_match, set_size, exact_match)
    # check for formerly computed stuff, and erase if allowed
    checkPangenomeFormerSpots(pangenome, force)
    # check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions=True,
                       needRGP=True, disable_bar=disable_bar)

    # get multigenic gene families
    logging.getLogger().info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    logging.getLogger().info("Detecting hotspots in the pangenome...")

    # predict spots
    spots = makeSpotGraph(pangenome.regions, multigenics, output, spot_graph, overlapping_match, set_size, exact_match)

    if len(spots) == 0:
        logging.getLogger().warning("No spots were detected.")
    else:
        logging.getLogger().info(f"{len(spots)} spots were detected")

    pangenome.addSpots(spots)
    pangenome.status["spots"] = "Computed"
    pangenome.parameters["spots"] = {}
    pangenome.parameters["spots"]["set_size"] = set_size
    pangenome.parameters["spots"]["overlapping_match"] = overlapping_match
    pangenome.parameters["spots"]["exact_match"] = exact_match

def checkParameterLogic(overlapping_match, set_size, exact_match):
    if overlapping_match >= set_size:
        raise Exception(f'--overlapping_match_hotspot ({overlapping_match}) cannot be bigger than (or equal to) '
                        f'--set_size_hotspot ({set_size})')
    if exact_match > set_size:
        raise Exception(f'--exact_match_size_hotspot ({exact_match}) cannot be bigger than '
                        f'--set_size_hotspot ({set_size})')

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.spot_graph:
        mkOutdir(args.output, args.force)
    if args.draw_hotspots or args.interest or args.fig_margin or args.priority:
        logging.getLogger().warning(
            "Options to draw the spots with the 'ppanggolin spot' subcommand have been deprecated, "
            "and are now dealt with in a dedicated subcommand 'ppanggolin drawspot'.")
    predictHotspots(pangenome, args.output, force=args.force, cpu=args.cpu, spot_graph=args.spot_graph,
                    overlapping_match=args.overlapping_match, set_size=args.set_size, exact_match=args.exact_match_size,
                    disable_bar=args.disable_prog_bar)
    writePangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def spotSubparser(subparser):
    parser = subparser.add_parser("spot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--spot_graph", required=False, action="store_true",
                          help="Writes a graph in .gexf format of pairs of blocks of single copy markers flanking RGPs,"
                               " supposedly belonging to the same hotspot")
    optional.add_argument("--overlapping_match", required=False, type=int, default=2,
                          help="The number of 'missing' persistent genes allowed when comparing flanking genes during "
                               "hotspot computations")
    optional.add_argument("--set_size", required=False, type=int, default=3,
                          help="Number of single copy markers to use as flanking genes for a RGP during "
                               "hotspot computation")
    optional.add_argument("--exact_match_size", required=False, type=int, default=1,
                          help="Number of perfectly matching flanking single copy markers required to associate RGPs "
                               "during hotspot computation (Ex: If set to 1, two RGPs are in the same hotspot "
                               "if both their 1st flanking genes are the same)")
    optional.add_argument("--draw_hotspots", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    optional.add_argument("--interest", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    optional.add_argument("--fig_margin", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    optional.add_argument("--priority", required=False, action="store_true",
                          help=argparse.SUPPRESS)  # This ensures compatibility with the old API
    # but does not use the option
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    return parser
