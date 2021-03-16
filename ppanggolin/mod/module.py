#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import argparse
import time
import os
import tempfile
import subprocess
from itertools import combinations
from statistics import mean, median
from collections import defaultdict
#installed libraries
from tqdm import tqdm
import networkx as nx
from gmpy2 import xmpz, popcount#pylint: disable=no-name-in-module
#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Module
from ppanggolin.formats import checkPangenomeInfo, writePangenome, ErasePangenome
from ppanggolin.utils import mkOutdir, restricted_float

def connected_components(g, removed, weight):
    """
        Yields subgraphs of each connected component you get when filtering edges based on the given weight.
    """
    for v in g.nodes:
        if v not in removed:
            c = set(_plain_bfs(g, v, removed, weight))
            yield c
            removed.update(c)

def _plain_bfs(g, source, removed, weight):
    """A fast BFS node generator, copied from networkx then adapted to the current use case"""
    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in removed:
                yield v
                removed.add(v)

                for n in g.neighbors(v):
                    if n not in removed:
                        edge_genes_v = g[v][n]["genes"][v]
                        edge_genes_n = g[v][n]["genes"][n]
                        #if the edge is indeed existent for most genes of both families, we use it
                        if len(edge_genes_n) / len(g.nodes[n]["genes"]) >= weight and len(edge_genes_v) / len(g.nodes[v]["genes"]) >= weight:
                            nextlevel.add(n)

def checkPangenomeFormerModules(pangenome, force):
    """ checks pangenome status and .h5 files for former modules, delete them if allowed or raise an error """
    if pangenome.status["modules"] == "inFile" and force == False:
        raise Exception("You are trying to detect modules on a pangenome which already has predicted modules. If you REALLY want to do that, use --force (it will erase modules previously predicted).")
    elif pangenome.status["modules"] == "inFile" and force == True:
        ErasePangenome(pangenome, modules = True)

def predictModules(pangenome, output, cpu, tmpdir, force=False, dup_margin=0.05, size=3, min_presence=3, transitive=5, jaccard=0.85, show_bar=True):
    #check statuses and load info
    checkPangenomeFormerModules(pangenome, force)
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needPartitions = True)

    ##compute the graph with transitive closure size provided as parameter
    start_time = time.time()
    logging.getLogger().info("Building the graph...")
    g = compute_mod_graph(pangenome.organisms, t=transitive, show_bar=show_bar)
    logging.getLogger().info(f"Took {round(time.time() - start_time,2)} seconds to build the graph to find modules in")
    logging.getLogger().info(f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges")

    start_time = time.time()
    #get all multigenic gene families
    multi = pangenome.get_multigenics(dup_margin, persistent=False)

    #extract the modules from the graph
    modules = compute_modules(g, multi, jaccard, min_presence, size=size)

    fams = set()
    for mod in modules:
        fams |= mod.families

    logging.getLogger().info(f"There are {len(fams)} families among {len(modules)} modules")
    logging.getLogger().info(f"Computing modules took {round(time.time() - start_time,2)} seconds")

    pangenome.addModules(modules)

    pangenome.status["modules"] = "Computed"
    pangenome.parameters["modules"] = {}
    pangenome.parameters["modules"]["size"] = size
    pangenome.parameters["modules"]["min_presence"] = min_presence
    pangenome.parameters["modules"]["transitive"] = transitive
    pangenome.parameters["modules"]["jaccard"] = jaccard
    pangenome.parameters["modules"]["dup_margin"] = dup_margin

def add_gene(obj, gene, fam_split=True):
    if fam_split:
        try:
            obj["genes"][gene.family].add(gene)
        except KeyError:
            try:
                obj["genes"][gene.family] = set([gene])
            except KeyError:
                obj["genes"] = {gene.family:set([gene])}
    else:
        try:
            obj["genes"].add(gene)
        except KeyError:
            obj["genes"] = set([gene])

def compute_mod_graph(organisms, t=1, show_bar=True):
    """
    Computes a graph using all provided genomes with a transitive closure of size t
    
    :param organisms: the list of organisms to compute the graph with
    :type list: list[:class:`ppanggolin.genome.Organism`]
    :param t: the size of the transitive closure
    :type t: int
    :param show_bar: whether to show a progress bar or not
    :type show_bar: bool
    """

    g = nx.Graph()
    for org in tqdm(organisms, unit="genome", disable=not show_bar):
        for contig in org.contigs:
            if len(contig.genes) > 0:
                start_gene = contig.genes[0]
                g.add_node(start_gene.family)
                add_gene(g.nodes[start_gene.family], start_gene, fam_split=False)
                for i, gene in enumerate(contig.genes):
                    for j, a_gene in enumerate(contig.genes[i+1:i+t+2], start=i+1):
                        g.add_edge(gene.family, a_gene.family)
                        edge = g[gene.family][a_gene.family]
                        add_gene(edge, gene)
                        add_gene(edge, a_gene)
                        if j == i+t+1 or i == 0:#if it's the last gene of the serie, or the first serie
                            add_gene(g.nodes[a_gene.family], a_gene, fam_split=False)

    return g

def lookForMultiNeighbors(g, multi, mod, classified, weight):
    linked_multi = set()
    for v in mod.core:
        for n in g.neighbors(v):
            if n in multi and n not in classified:
                edge_genes_v = g[v][n]["genes"][v]
                #if the edge is indeed existent for most genes of the module's families, we use it
                if len(edge_genes_v) / len(g.nodes[v]["genes"]) >= weight:
                    linked_multi.add(n)
    return linked_multi

def compute_modules(g, multi, weight, min_fam, size):
    """
    Computes modules using a graph built by :func:`ppanggolin.mod.module.compute_mod_graph` and differents parameters defining how restrictive the modules will be.

    :param g: The networkx graph from :func:`ppanggolin.mod.module.compute_mod_graph`
    :type g: :class:`networkx.Graph`
    :param multi: a set of families :class:`ppanggolin.geneFamily.GeneFamily` considered multigenic
    :type multi: set
    :param weight: the minimal jaccard under which edges are not considered
    :type weight: float
    :param min_fam: the minimal number of presence under which the family is not considered
    :type min_fam: int
    """
    #removing families with low presence
    removed = set([fam for fam in g.nodes if len(fam.organisms) < min_fam])

    modules = set()
    #define core modules
    classified_families = set(removed)
    c = 0
    for comp in connected_components(g, removed, weight):
        if len(comp) >= size:#keep only the modules with at least 'size' non-multigenic genes.
            if not any(fam.namedPartition == "persistent" and fam not in multi for fam in comp):#remove 'persistent' and non-multigenic modules
                modules.add(Module(ID = c, core=comp))
                classified_families |= comp
                c += 1

    ##parse the graph to get multigenics near modules that are not classified, to associate them to modules
    nb_added = 0
    nb_mod = 0
    for mod in modules:
        curr_multi = lookForMultiNeighbors(g, multi, mod, classified_families, weight)

        if len(curr_multi) > 0:
            nb_added += len(curr_multi)
            nb_mod +=1
            mod.associate_families(curr_multi)
    logging.getLogger().info(f"Associated {nb_added} multigenic families to {nb_mod} modules")
    return modules

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    mkOutdir(args.output, args.force)
    predictModules(pangenome = pangenome, output=args.output, cpu = args.cpu, tmpdir = args.tmpdir, force=args.force, dup_margin=args.dup_margin, size=args.size, min_presence=args.min_presence, transitive=args.transitive, jaccard=args.jaccard, show_bar=args.show_prog_bars)
    writePangenome(pangenome, pangenome.file, args.force, show_bar=args.show_prog_bars)


def moduleSubparser(subparser):
    parser = subparser.add_parser("module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument("--size", required=False, type=int, default=3, help = "Minimal number of gene family in a module")
    optional.add_argument("--dup_margin", required=False, type=restricted_float, default=0.05, help = "minimum ratio of organisms in which the family must have multiple genes for it to be considered 'duplicated'")
    optional.add_argument("--min_presence", required=False, type=int, default=5, help = "Minimum number of times the module needs to be present in the pangenome to be reported. Increasing it will improve precision but lower sensitivity.")
    optional.add_argument("--transitive",required=False, type=int, default = 5, help = "Size of the transitive closure used to build the graph. This indicates the number of non related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little.")
    optional.add_argument("--jaccard", required=False, type=restricted_float, default=0.85, help = "minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot." )
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser