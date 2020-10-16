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

#installed libraries
from tqdm import tqdm
import networkx as nx
from gmpy2 import xmpz, popcount#pylint: disable=no-name-in-module
#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.region import Region
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir


class Module:
    def __init__(self, geneFamilies):
        if not all(isinstance(fam, GeneFamily) for fam in geneFamilies):
            raise Exception(f"You provided elements that were not GeneFamily objetcs. Modules are only made of GeneFamily")
        self.families = set(geneFamilies)

    @property
    def distance(self):
        if hasattr(self, "_distance"):
            return self._distance
        else:
            self._get_dist()
            return self._distance

    def _get_dist(self):
        union = xmpz(0)#infinitely long vector of 0-bit
        inter = xmpz(-1)#infinitely long vector of 1-bit
        for fam in self.families:
            union = union | fam.bitarray
            inter = inter & fam.bitarray
        self._distance = float(popcount(inter) / popcount(union))


def getFam2Mod(modules):
    fam2mod = {}
    for mod in modules:
        for fam in mod.families:
            fam2mod[fam] = mod
    return fam2mod


def connected_components(g, removed, weight):
    """
        Yields subgraphs of each connected component.
    """
    for v in g.nodes:
        if v not in removed:
            c = set(_plain_bfs(g, v, removed, weight))
            yield c
            removed.update(c)

def _plain_bfs(g, source, removed, weight):
    """A fast BFS node generator, copied from networkx"""
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
                        # edge_rgps = set([rgp.organism for rgp in g[v][n]["rgp"]])
                        # if len(edge_rgps) / popcount( v.bitarray | n.bitarray) >= weight:
                        #     nextlevel.add(n)


def set_mod_identifiers(fam2mod, modules):
    fam2id = {}
    mod2id = {}
    for i, mod in enumerate(modules):
        mod2id[mod] = i
        for fam in mod.families:
            fam2id[fam] = i
    return fam2id, mod2id

def predictModules(pangenome, output, cpu, tmpdir):
    #check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions = True, needRGP = True)
    ##do the module thing
    pangenome.computeFamilyBitarrays()#might need that

    start_time = time.time()
    logging.getLogger().info("Building the graph...")
    g = compute_rgp_graph(pangenome.organisms, t=5)
    logging.getLogger().info(f"Took {round(time.time() - start_time,2)} seconds to build the graph")
    logging.getLogger().info(f"There are {nx.number_of_nodes(g)} nodes and {nx.number_of_edges(g)} edges")

    start_time = time.time()
    modules = compute_modules(g, 0.9, 5)
    fam2mod = getFam2Mod(modules)
    logging.getLogger().info(f"There are {len(fam2mod)} families among {len(modules)} modules")
    logging.getLogger().info(f"Computing modules took {round(time.time() - start_time,2)} seconds")

    fam2id, mod2id = set_mod_identifiers(fam2mod, modules)
    # write_cgview(pangenome._orgGetter["GCF_000005845.2_ASM584v2"], output, fam2id)


def write_cgview(genome, output, fam2id):
    fgenes = open(output + "/" + genome.name+"_cgview_genes.tab","w")
    fparts = open(output+"/" + genome.name+"_cgview_partitions.tab","w")
    fmods = open(output+"/" + genome.name + "_cgview_modules.tab","w")
    fgenes.write("name\ttype\tstart\tstop\tstrand\n")
    fparts.write("name\ttype\tstart\tstop\tstrand\n")
    fmods.write("name\ttype\tstart\tstop\tstrand\n")
    prev_size = 0
    for contig in genome.contigs:
        for gene in contig.genes + list(contig.RNAs):
            fgenes.write(f"{gene.name}\t{gene.type}\t{gene.start + prev_size}\t{gene.stop + prev_size}\t{gene.strand}\n")
            if gene.type == "CDS":
                fparts.write(f"{gene.family.name}\t{gene.family.namedPartition}\t{gene.start + prev_size}\t{gene.stop + prev_size}\t{gene.strand}\n")
                mod = fam2id.get(gene.family)
                if mod is not None:
                    #then the family has an assigned module
                    fmods.write(f"{gene.family.name}\tmodule{mod}\t{gene.start + prev_size}\t{gene.stop + prev_size}\t{gene.strand}\n")
    fgenes.close()
    fparts.close()
    fmods.close()

def add_rgp(obj, rgp):
    try:
        obj["rgp"].add(rgp)
    except KeyError:
        obj["rgp"] = set([rgp])

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

def compute_rgp_graph(organisms, t=1):
    """
    Computes a graph using the RGPs only with a transitive closure of size t
    
    :param organisms: the list of organisms to compute the graph with
    :type list: list[:class:`ppanggolin.genome.Organism`]
    :param t: the size of the transitive closure
    :type t: int
    """
    # g = nx.Graph()
    # for region in pangenome.regions:
    #     start_gene = region.startGene
        
    #     g.add_node(start_gene.family)
    #     # add_gene(g.nodes[start_gene.family], start_gene)
    #     # add_rgp(g.nodes[start_gene.family], region)
    #     for i, gene in enumerate(region.genes):
    #         for j, a_gene in enumerate(region.genes[i+1:i+t+2], start=i+1):
    #             g.add_edge(gene.family, a_gene.family)
    #             edge = g[gene.family][a_gene.family]
    #             #add_rgp(edge, region)#add rgp to the edge
    #             add_gene(edge, gene)
    #             add_gene(edge, a_gene)
    #             #if j == i+t+1 or i == 0:#if it's the last gene of the serie, or the first serie
    #                 #add_rgp(g.nodes[a_gene.family], region)#add rgp to the eventual new node
    #                 #add_gene(g.nodes[a_gene.family], a_gene)

    g = nx.Graph()
    for org in tqdm(organisms, unit="genome", disable = True):
        for contig in org.contigs:
            start_gene = contig.genes[0]
            g.add_node(start_gene.family)
            add_gene(g.nodes[start_gene.family], start_gene, fam_split=False)
            # add_rgp(g.nodes[start_gene.family], region)
            for i, gene in enumerate(contig.genes):
                for j, a_gene in enumerate(contig.genes[i+1:i+t+2], start=i+1):
                    g.add_edge(gene.family, a_gene.family)
                    edge = g[gene.family][a_gene.family]
                    #add_rgp(edge, region)#add rgp to the edge
                    add_gene(edge, gene)
                    add_gene(edge, a_gene)
                    if j == i+t+1 or i == 0:#if it's the last gene of the serie, or the first serie
                        #add_rgp(g.nodes[a_gene.family], region)#add rgp to the eventual new node
                        add_gene(g.nodes[a_gene.family], a_gene, fam_split=False)

    return g

def compute_modules(g, weight, min_fam):
    """

    :param weight: the minimal weight under which edges are not considered
    :type weight: float
    :param min_fam: the minimal number of presence under which the family is not considered
    :type min_fam: int
    """
    #removing families with low presence
    removed = set([fam for fam in g.nodes if len(fam.organisms) < min_fam])
    max_size = 0
    modules = set()
    for comp in connected_components(g, removed, weight):
        if len(comp) >= 3:
            if len(comp) > max_size:
                max_size = len(comp)
            modules.add(Module(geneFamilies=comp))

    return modules

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    mkOutdir(args.output, args.force)
    predictModules(pangenome = pangenome, output=args.output, cpu = args.cpu, tmpdir = args.tmpdir)
    #write modules to the hdf5 file

def moduleSubparser(subparser):
    parser = subparser.add_parser("module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument("--size", required=False, type=int, default=3, help = "Minimal number of gene family in a module")
    optional.add_argument("--min_presence", required=False, type=int, default=5, help = "Minimum number of times the module needs to be present in the pangenome to be reported")

    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser