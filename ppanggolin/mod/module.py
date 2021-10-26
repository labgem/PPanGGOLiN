#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse
import time
import os
import tempfile
import subprocess
from itertools import combinations
from statistics import mean, median
from collections import defaultdict

# installed libraries
from tqdm import tqdm
import networkx as nx
from gmpy2 import xmpz, popcount  # pylint: disable=no-name-in-module

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Module
from ppanggolin.formats import checkPangenomeInfo, writePangenome, ErasePangenome
from ppanggolin.utils import mkOutdir, restricted_float, add_gene, connected_components


def checkPangenomeFormerModules(pangenome, force):
    """ checks pangenome status and .h5 files for former modules, delete them if allowed or raise an error """
    if pangenome.status["modules"] == "inFile" and not force:
        raise Exception("You are trying to detect modules on a pangenome which already has predicted modules. "
                        "If you REALLY want to do that, use --force (it will erase modules previously predicted).")
    elif pangenome.status["modules"] == "inFile" and force:
        ErasePangenome(pangenome, modules=True)


def predictModules(pangenome, cpu, tmpdir, force=False, dup_margin=0.05, size=3, min_presence=2, transitive=4,
                   jaccard=0.85, disable_bar=False):
    # check statuses and load info
    checkPangenomeFormerModules(pangenome, force)
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needPartitions=True, disable_bar=disable_bar)

    # compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger().info("Building the graph...")
    g = compute_mod_graph(pangenome.organisms, t=transitive, disable_bar=disable_bar)
    logging.getLogger().info(f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find modules in")
    logging.getLogger().info(f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges")

    start_time = time.time()
    # get all multigenic gene families
    multi = pangenome.get_multigenics(dup_margin, persistent=False)

    # extract the modules from the graph
    modules = compute_modules(g, multi, jaccard, min_presence, size=size)

    fams = set()
    for mod in modules:
        fams |= mod.families

    logging.getLogger().info(f"There are {len(fams)} families among {len(modules)} modules")
    logging.getLogger().info(f"Computing modules took {round(time.time() - start_time, 2)} seconds")

    pangenome.addModules(modules)

    pangenome.status["modules"] = "Computed"
    pangenome.parameters["modules"] = {}
    pangenome.parameters["modules"]["size"] = size
    pangenome.parameters["modules"]["min_presence"] = min_presence
    pangenome.parameters["modules"]["transitive"] = transitive
    pangenome.parameters["modules"]["jaccard"] = jaccard
    pangenome.parameters["modules"]["dup_margin"] = dup_margin


def compute_mod_graph(organisms, t=1, disable_bar=False):
    """
    Computes a graph using all provided genomes with a transitive closure of size t
    
    :param organisms: the list of organisms to compute the graph with
    :type list: list[:class:`ppanggolin.genome.Organism`]
    :param t: the size of the transitive closure
    :type t: int
    :param disable_bar: whether to show a progress bar or not
    :type disable_bar: bool
    """

    g = nx.Graph()
    for org in tqdm(organisms, unit="genome", disable=disable_bar):
        for contig in org.contigs:
            if len(contig.genes) > 0:
                start_gene = contig.genes[0]
                g.add_node(start_gene.family)
                add_gene(g.nodes[start_gene.family], start_gene, fam_split=False)
                for i, gene in enumerate(contig.genes):
                    for j, a_gene in enumerate(contig.genes[i + 1:i + t + 2], start=i + 1):
                        g.add_edge(gene.family, a_gene.family)
                        edge = g[gene.family][a_gene.family]
                        add_gene(edge, gene)
                        add_gene(edge, a_gene)
                        if j == i + t + 1 or i == 0:  # if it's the last gene of the serie, or the first serie
                            add_gene(g.nodes[a_gene.family], a_gene, fam_split=False)

    return g


def compute_modules(g, multi, weight, min_fam, size):
    """
    Computes modules using a graph built by :func:`ppanggolin.mod.module.compute_mod_graph` and different parameters
    defining how restrictive the modules will be.

    :param g: The networkx graph from :func:`ppanggolin.mod.module.compute_mod_graph`
    :type g: :class:`networkx.Graph`
    :param multi: a set of families :class:`ppanggolin.geneFamily.GeneFamily` considered multigenic
    :type multi: set
    :param weight: the minimal jaccard under which edges are not considered
    :type weight: float
    :param min_fam: the minimal number of presence under which the family is not considered
    :type min_fam: int
    """
    # removing families with low presence
    removed = set([fam for fam in g.nodes if len(fam.organisms) < min_fam])

    modules = set()
    c = 0
    for comp in connected_components(g, removed, weight):
        if len(comp) >= size and not any(fam.namedPartition == "persistent" and
                                         fam not in multi for fam in comp):
            # keep only the modules with at least 'size' non-multigenic genes and
            # remove 'persistent' and non-multigenic modules
            modules.add(Module(ID=c, families=comp))
            c += 1
    return modules


def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    predictModules(pangenome=pangenome, cpu=args.cpu, tmpdir=args.tmpdir, force=args.force, dup_margin=args.dup_margin,
                   size=args.size, min_presence=args.min_presence, transitive=args.transitive, jaccard=args.jaccard,
                   disable_bar=args.disable_prog_bar)
    writePangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def moduleSubparser(subparser):
    parser = subparser.add_parser("module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--size", required=False, type=int, default=3,
                          help="Minimal number of gene family in a module")
    optional.add_argument("--dup_margin", required=False, type=restricted_float, default=0.05,
                          help="minimum ratio of organisms in which the family must have multiple genes"
                               " for it to be considered 'duplicated'")
    optional.add_argument("-m", "--min_presence", required=False, type=int, default=2,
                          help="Minimum number of times the module needs to be present in the pangenome to be reported."
                               " Increasing it will improve precision but lower sensitivity.")
    optional.add_argument("-t", "--transitive", required=False, type=int, default=4,
                          help="Size of the transitive closure used to build the graph. "
                               "This indicates the number of non related genes allowed in-between two related genes. "
                               "Increasing it will improve precision but lower sensitivity a little.")
    optional.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                          help="minimum jaccard similarity used to filter edges between gene families. "
                               "Increasing it will improve precision but lower sensitivity a lot.")

    required = parser.add_argument_group(title="Required arguments",
                                         description="One of the following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    return parser
