#!/usr/bin/env python3

# default libraries
import logging
import argparse
import time
from pathlib import Path

# installed libraries
from tqdm import tqdm
import networkx as nx

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Module
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import restricted_float, add_gene, connected_components


def check_pangenome_former_modules(pangenome: Pangenome, force: bool = False):
    """
    Checks pangenome status and .h5 files for former modules, delete them if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on pangenome by erasing already present modules
    """
    if pangenome.status["modules"] == "inFile" and not force:
        raise Exception(
            "You are trying to detect modules on a pangenome which already has predicted modules. "
            "If you REALLY want to do that, use --force (it will erase modules previously predicted)."
        )
    elif pangenome.status["modules"] == "inFile" and force:
        erase_pangenome(pangenome, modules=True)


def compute_mod_graph(pangenome: Pangenome, t: int = 1, disable_bar: bool = False):
    """
    Computes a graph using all provided genomes with a transitive closure of size t

    :param pangenome: pangenome with organisms to compute the graph
    :param t: the size of the transitive closure
    :param disable_bar: whether to show a progress bar or not
    """

    g = nx.Graph()
    for org in tqdm(
        pangenome.organisms,
        total=pangenome.number_of_organisms,
        unit="genome",
        disable=disable_bar,
    ):
        for contig in org.contigs:
            if contig.number_of_genes > 0:
                start_gene = contig[0]
                g.add_node(start_gene.family)
                add_gene(g.nodes[start_gene.family], start_gene, fam_split=False)
                for i, gene in enumerate(contig.genes):
                    for j, a_gene in enumerate(
                        contig.get_genes(i + 1, i + t + 2, outrange_ok=True),
                        start=i + 1,
                    ):
                        g.add_edge(gene.family, a_gene.family)
                        edge = g[gene.family][a_gene.family]
                        add_gene(edge, gene)
                        add_gene(edge, a_gene)
                        if (
                            j == i + t + 1 or i == 0
                        ):  # if it's the last gene of the series, or the first series
                            add_gene(g.nodes[a_gene.family], a_gene, fam_split=False)

    return g


def compute_modules(
    g: nx.Graph, multi: set, weight: float = 0.85, min_fam: int = 2, size: int = 3
):
    """
    Computes modules using a graph built by :func:`ppanggolin.mod.module.compute_mod_graph` and different parameters
    defining how restrictive the modules will be.

    :param g: The networkx graph from :func:`ppanggolin.mod.module.compute_mod_graph`
    :param multi: a set of families :class:`ppanggolin.geneFamily.GeneFamily` considered multigenic
    :param weight: the minimal jaccard under which edges are not considered
    :param min_fam: the minimal number of presence under which the family is not considered
    :param size: Minimal number of gene family in a module
    """

    # removing families with low presence
    removed = {fam for fam in g.nodes if fam.number_of_organisms < min_fam}

    modules = set()
    c = 0
    for comp in connected_components(g, removed, weight):
        if len(comp) >= size and not any(
            fam.named_partition == "persistent" and fam not in multi for fam in comp
        ):
            # keep only the modules with at least 'size' non-multigenic genes and
            # remove 'persistent' and non-multigenic modules
            modules.add(Module(module_id=c, families=comp))
            c += 1
    return modules


def predict_modules(
    pangenome: Pangenome,
    dup_margin: float = 0.05,
    size: int = 3,
    min_presence: int = 2,
    transitive: int = 4,
    jaccard: float = 0.85,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Main function to predict module

    :param pangenome: Pangenome object with Gene Families, Annotation and Partition
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    :param size: Minimal number of gene family in a module
    :param min_presence: Minimum number of times the module needs to be present in the pangenome to be reported.
    :param transitive: Size of the transitive closure used to build the graph.
    :param jaccard: minimum jaccard similarity used to filter edges between gene families.
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    # check statuses and load info
    check_pangenome_former_modules(pangenome, force)
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_partitions=True,
        disable_bar=disable_bar,
    )

    # compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger("PPanGGOLiN").info("Building the graph...")
    g = compute_mod_graph(pangenome, t=transitive, disable_bar=disable_bar)
    logging.getLogger("PPanGGOLiN").info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find modules in"
    )
    logging.getLogger("PPanGGOLiN").info(
        f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges"
    )

    start_time = time.time()
    # get all multigenic gene families
    multi = pangenome.get_multigenics(dup_margin, persistent=False)

    # extract the modules from the graph
    modules = compute_modules(g, multi, jaccard, min_presence, size=size)

    fams = set()
    for mod in modules:
        fams |= set(mod.families)
        pangenome.add_module(mod)

    logging.getLogger("PPanGGOLiN").info(
        f"There are {len(fams)} families among {len(modules)} modules"
    )
    logging.getLogger("PPanGGOLiN").info(
        f"Computing modules took {round(time.time() - start_time, 2)} seconds"
    )

    pangenome.status["modules"] = "Computed"
    pangenome.parameters["module"] = {}
    pangenome.parameters["module"]["size"] = size
    pangenome.parameters["module"]["min_presence"] = min_presence
    pangenome.parameters["module"]["transitive"] = transitive
    pangenome.parameters["module"]["jaccard"] = jaccard
    pangenome.parameters["module"]["dup_margin"] = dup_margin


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    predict_modules(
        pangenome=pangenome,
        dup_margin=args.dup_margin,
        size=args.size,
        min_presence=args.min_presence,
        transitive=args.transitive,
        jaccard=args.jaccard,
        force=args.force,
        disable_bar=args.disable_prog_bar,
    )
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "module", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_module(parser)
    return parser


def parser_module(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of module command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "--size",
        required=False,
        type=int,
        default=3,
        help="Minimal number of gene family in a module",
    )
    optional.add_argument(
        "--dup_margin",
        required=False,
        type=restricted_float,
        default=0.05,
        help="minimum ratio of genomes in which the family must have multiple genes"
        " for it to be considered 'duplicated'",
    )
    optional.add_argument(
        "-m",
        "--min_presence",
        required=False,
        type=int,
        default=2,
        help="Minimum number of times the module needs to be present in the pangenome to be reported."
        " Increasing it will improve precision but lower sensitivity.",
    )
    optional.add_argument(
        "-t",
        "--transitive",
        required=False,
        type=int,
        default=4,
        help="Size of the transitive closure used to build the graph. "
        "This indicates the number of non related genes allowed in-between two related genes. "
        "Increasing it will improve precision but lower sensitivity a little.",
    )
    optional.add_argument(
        "-s",
        "--jaccard",
        required=False,
        type=restricted_float,
        default=0.85,
        help="minimum jaccard similarity used to filter edges between gene families. "
        "Increasing it will improve precision but lower sensitivity a lot.",
    )

    optional.add_argument(
        "-c",
        "--cpu",
        required=False,
        default=1,
        type=int,
        help="Number of available cpus",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_module(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
