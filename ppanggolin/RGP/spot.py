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
from ppanggolin.region import Region, Spot
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import mk_outdir


def comp_border(border1: list, border2: list, overlapping_match: int = 2,
                set_size: int = 3, exact_match: int = 1) -> bool:
    """
    Compare two border

    :param border1:
    :param border2:
    :param overlapping_match:
    :param set_size:
    :param exact_match:

    :return:
    """
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


def check_sim(pair_border1: list, pair_border2: list, overlapping_match: int = 2,
              set_size: int = 3, exact_match: int = 1) -> bool:
    """
    Checks if the two pairs of exact_match first gene families are identical,
    or eventually if they overlap in an ordered way at least 'overlapping_match'

    :param pair_border1: First flanking gene families pair
    :param pair_border2: Second flanking gene families pair
    :param overlapping_match: Number of missing persistent genes allowed when comparing flanking genes
    :param set_size: Number of single copy markers to use as flanking genes for RGP during hotspot computation
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs

    :return: Whether identical gene families or not
    """
    b1pair = [False, False]
    b2pair = [False, False]
    for countb1, b1 in enumerate(pair_border1):
        for countb2, b2 in enumerate(pair_border2):
            if comp_border(b2, b1, overlapping_match, set_size, exact_match):
                b1pair[countb1] = True
                b2pair[countb2] = True

    if b1pair[0] and b1pair[1] and b2pair[0] and b2pair[1]:
        return True
    return False


def make_spot_graph(rgps: list, multigenics: set, output: str, spot_graph: bool = False, overlapping_match: int = 2,
                    set_size: int = 3, exact_match: int = 1) -> list:
    """
    Create a spot graph from pangenome RGP

    :param rgps: list of pangenome RGP
    :param multigenics: pangenome graph multigenic persistent families
    :param output: Output directory to save the spot graph
    :param spot_graph: Writes gexf graph of pairs of blocks of single copy markers flanking RGPs from same hotspot
    :param overlapping_match: Number of missing persistent genes allowed when comparing flanking genes
    :param set_size: Number of single copy markers to use as flanking genes for RGP during hotspot computation
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs

    :return: list of computed spot
    """

    def add_new_node(g: nx.Graph, region: Region, borders: list):
        """
        Add bordering region as node to graph

        :param g: spot graph
        :param region: region in spot
        :param borders: bordering families in spot
        """
        blocks = str(sorted([[gene.family.ID for gene in borders[0]], [gene.family.ID for gene in borders[1]]],
                            key=lambda x: x[0]))
        g.add_node(blocks)
        try:
            g.nodes[blocks]["nb_rgp"] += 1
            g.nodes[blocks]["rgp"].add(region)
        except KeyError:
            g.nodes[blocks]["nb_rgp"] = 1
            g.nodes[blocks]["border1"] = [gene.family for gene in borders[1]]
            g.nodes[blocks]["border0"] = [gene.family for gene in borders[0]]
            g.nodes[blocks]["rgp"] = {region}

    graph_spot = nx.Graph()
    lost = 0
    used = 0
    for rgp in rgps:
        border = rgp.get_bordering_genes(set_size, multigenics)
        if len(border[0]) < set_size or len(border[1]) < set_size:
            lost += 1
        else:
            used += 1
            add_new_node(graph_spot, rgp, border)
    logging.getLogger().info(f"{lost} RGPs were not used as they are on a contig border (or have less than {set_size} "
                             f"persistent gene families until the contig border)")
    logging.getLogger().info(f"{used} RGPs are being used to predict spots of insertion")
    node_list = list(graph_spot.nodes)
    logging.getLogger().info(f"{len(node_list)} number of different pairs of flanking gene families")
    for i, nodei in enumerate(node_list[:-1]):
        for nodej in node_list[i + 1:]:
            node_obj_i = graph_spot.nodes[nodei]
            node_obj_j = graph_spot.nodes[nodej]
            if check_sim([node_obj_i["border0"], node_obj_i["border1"]], [node_obj_j["border0"], node_obj_j["border1"]],
                         overlapping_match, set_size, exact_match):
                graph_spot.add_edge(nodei, nodej)
    spots = []
    spot_id = 0
    for comp in nx.algorithms.components.connected_components(graph_spot):
        curr_spot = Spot(spot_id)
        spots.append(curr_spot)
        for node in comp:
            curr_spot.add_regions(graph_spot.nodes[node]["rgp"])
        spot_id += 1

    if spot_graph:
        for node in graph_spot.nodes:
            del graph_spot.nodes[node]["border0"]
            del graph_spot.nodes[node]["border1"]
            del graph_spot.nodes[node]["rgp"]

        nx.readwrite.gexf.write_gexf(graph_spot, output + "/spotGraph.gexf")
        nx.readwrite.graphml.write_graphml(graph_spot, output + "/spotGraph.graphml")
    return spots


def check_pangenome_former_spots(pangenome: Pangenome, force: bool = False):
    """
    checks pangenome status and .h5 files for former spots, delete them if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    """
    if pangenome.status["spots"] == "inFile" and not force:
        raise Exception("You are trying to detect spots on a pangenome which already has predicted spots. "
                        "If you REALLY want to do that, use --force (it will erase spots previously predicted).")
    elif pangenome.status["spots"] == "inFile" and force:
        erase_pangenome(pangenome, spots=True)


def predict_hotspots(pangenome: Pangenome, output: str, spot_graph: bool = False, overlapping_match: int = 2,
                     set_size: int = 3, exact_match: int = 1, force: bool = False, disable_bar: bool = False):
    """
    Main function to predict hotspot

    :param pangenome: Blank pangenome object
    :param output: Output directory to save the spot graph
    :param spot_graph: Writes gexf graph of pairs of blocks of single copy markers flanking RGPs from same hotspot
    :param overlapping_match: Number of missing persistent genes allowed when comparing flanking genes
    :param set_size: Number of single copy markers to use as flanking genes for RGP during hotspot computation
    :param exact_match: Number of perfectly matching flanking single copy markers required to associate RGPs
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    # check that given parameters for hotspot computation make sense
    if overlapping_match >= set_size:
        raise Exception(f'--overlapping_match_hotspot ({overlapping_match}) cannot be bigger than (or equal to) '
                        f'--set_size_hotspot ({set_size})')
    if exact_match > set_size:
        raise Exception(f'--exact_match_size_hotspot ({exact_match}) cannot be bigger than '
                        f'--set_size_hotspot ({set_size})')
    # check for formerly computed stuff, and erase if allowed
    check_pangenome_former_spots(pangenome, force)
    # check statuses and load info
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_partitions=True,
                         need_rgp=True, disable_bar=disable_bar)

    # get multigenic gene families
    logging.getLogger().info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])

    logging.getLogger().info("Detecting hotspots in the pan...")

    # predict spots
    spots = make_spot_graph(pangenome.regions, multigenics, output, spot_graph, overlapping_match, set_size,
                            exact_match)

    if len(spots) == 0:
        logging.getLogger().warning("No spots were detected.")
    else:
        logging.getLogger().info(f"{len(spots)} spots were detected")

    pangenome.add_spots(spots)
    pangenome.status["spots"] = "Computed"
    pangenome.parameters["spots"] = {}
    pangenome.parameters["spots"]["set_size"] = set_size
    pangenome.parameters["spots"]["overlapping_match"] = overlapping_match
    pangenome.parameters["spots"]["exact_match"] = exact_match


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    if args.spot_graph:
        mk_outdir(args.output, args.force)
    if args.draw_hotspots or args.interest or args.fig_margin or args.priority:
        logging.getLogger().warning(
            "Options to draw the spots with the 'ppanggolin spot' subcommand have been deprecated, "
            "and are now dealt with in a dedicated subcommand 'ppanggolin drawspot'.")
    predict_hotspots(pangenome, args.output, force=args.force, spot_graph=args.spot_graph,
                     overlapping_match=args.overlapping_match, set_size=args.set_size,
                     exact_match=args.exact_match_size, disable_bar=args.disable_prog_bar)
    write_pangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("spot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_spot(parser)
    return parser


def parser_spot(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of spot command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=False, type=str, help="The pangenome .h5 file")
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


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_spot(main_parser)
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
