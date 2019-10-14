#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import argparse
import time
import os
from collections import defaultdict
import random

#installed libraries
from tqdm import tqdm
import networkx as nx

#local libraries
from ppanggolin.pangenome import Pangenome, Region
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir

class MatriceNode:
    def __init__(self, state, score, prev, gene):
        self.state = state  # state of the node. 1 for RGP and 0 for not RGP.
        self.score = score if score > 0 else 0  # current score of the node
        self.prev = prev  # previous matriceNode
        self.gene = gene  # gene this node corresponds to

    def changes(self, state, score):
        # state of the node. 1 for RGP and 0 for not RGP.
        self.state = 1 if score >= 0 else 0
        # current score of the node. If the given score is negative, set to 0.
        self.score = score if score >= 0 else 0

def extractRGP(contig, node, ID):
    """
        Extract the region from the given starting node
    """
    new_region = Region(contig.name + "#" + str(ID))
    while node.state:
        new_region.append(node.gene)
        node.state = 0
        node.score = 0
        node = node.prev
        if node is None:#it's the end of the contig and the end of the region.
            break
    return new_region

def rewriteMatrix(contig, matrix, index, persistent, continuity, multi):
    """
        ReWrite the matrice from the given index of the node that started a region.
    """
    prev = matrix[index]
    index += 1
    if index > len(matrix) and contig.is_circular:
        index = 0
    # else the node was the last one of the contig, and there is nothing to do
    if index < len(matrix):
        nextNode = matrix[index]
        nbPerc = 0
        while nextNode.state:  # while the old state is not 0, recompute the scores.
            if nextNode.gene.family.namedPartition == "persistent" and nextNode.gene.family not in multi:
                modif = -pow(persistent, nbPerc)
                nbPerc += 1
            else:
                modif = continuity
                nbPerc = 0

            curr_score = modif + prev.score
            curr_state = 1 if curr_score >= 0 else 0
            # scores can't be negative. If they are, they'll be set to 0.
            matrix[index].changes(curr_state, curr_score)
            index += 1
            if index >= len(matrix):
                if contig.is_circular:
                    index = 0
                else:
                    # else we're at the end of the contig, so there are no more computations. Get out of the loop
                    break

            prev = nextNode
            nextNode = matrix[index]

def initMatrices(contig, persistent_penalty, variable_gain, multi ):
    """initialize the vector of score/state nodes"""
    mat = []
    prev = None
    nbPerc = 0
    for gene in contig.genes:
        if gene.family.namedPartition == "persistent" and gene.family not in multi:
            modif = -pow(persistent_penalty, nbPerc)
            nbPerc += 1
        else:
            modif = variable_gain
            nbPerc = 0

        curr_score = modif + prev.score if prev is not None else modif
        curr_state = 1 if curr_score >= 0 else 0
        prev = MatriceNode(curr_state, curr_score, prev, gene)
        mat.append(prev)

    # if the contig is circular and we're in a rgp state, we need to continue from the "starting" gene until we leave rgp state.
    if contig.is_circular and curr_state:
        # the previous node of the first processed gene is the last node.
        mat[0].prev = prev
        lastNode = prev  # saving the last node that was inserted.
        curr_score = prev.score
        c = 0
        nbPerc = 0
        while curr_state:  # while state is rgp.
            matNode = mat[c]
            if matNode == lastNode:  # then we've parsed the entire contig twice. The whole sequence is a rgp so we're stopping the iteration now, otherwise we'll loop indefinitely
                break

            if matNode.gene.family.namedPartition == "persistent" and matNode.gene.family not in multi:
                modif = -pow(persistent_penalty, nbPerc)
                nbPerc += 1
            else:
                modif = variable_gain
                nbPerc = 0

            curr_score = modif + prev.score
            curr_state = 1 if curr_score >= 0 else 0
            matNode.changes(curr_state, curr_score)
            c += 1
    return mat

def mkRegions(contig, matrix, min_length, min_score, persistent, continuity, multi):
    # processing matrix and 'emptying' it to get the regions.
    def maxIndexNode(lst):
        """gets the last node with the highest score from a list of matriceNode"""
        if isinstance(lst, list):
            # init with the first element of the list
            maxScore = lst[0].score
            maxIndex = 0
            for index, node in enumerate(lst):
                if node.score >= maxScore:
                    maxScore = node.score
                    maxIndex = index
            return (maxScore, maxIndex)
        else:
            raise TypeError("List of matriceNode is expected. The detected type was " + type(lst))

    contigRegions = set()
    val, index = maxIndexNode(matrix)
    c = 0
    while val >= min_score:
        new_region = extractRGP(contig, matrix[index], len(contigRegions))
        new_region.score = val
        if (new_region[0].stop - new_region[-1].start) > min_length:
            contigRegions.add(new_region)
        rewriteMatrix(contig, matrix, index, persistent, continuity, multi)
        val, index = maxIndexNode(matrix)
        c += 1
    return contigRegions

def compute_org_rgp(organism, persistent_penalty, variable_gain, min_length, min_score, multigenics):
    orgRegions = set()
    for contig in organism.contigs:
        ## can definitely multiprocess this part, as not THAT much information is needed...
        matrix = initMatrices(contig, persistent_penalty, variable_gain, multigenics)
        orgRegions |= mkRegions(contig, matrix, min_length, min_score, persistent_penalty, variable_gain, multigenics)

    return orgRegions

def get_multigenics(pangenome, dup_margin):
    """
        Returns the multigenic persistent families of the pangenome graph. A family will be considered multigenic if it is duplicated in more than 5% of the genomes where it is present.
    """
    multigenics = set()
    for fam in pangenome.geneFamilies:
        if fam.namedPartition == "persistent":
            dup=len([genes for org, genes in fam.getOrgDict().items() if len(genes) > 1])
            if (dup / len(fam.organisms)) >= dup_margin:#tot / nborgs >= 1.05
                multigenics.add(fam)
    logging.getLogger().info(f"{len(multigenics)} gene families are defined as being multigenic. (duplicated in more than {dup_margin} of the genomes)")
    return multigenics


def spot_distribution(spots, output):
    """takes in spots are a list of sets of rgps"""
    fdistrib = open(output + "/spot_rgp_distribution.tsv","w")
    for c, rgps in enumerate(spots):
        fdistrib.write(str(c) + "\t" + str(len(rgps)) + "\t" + str(len(getUniqRGP(rgps))) + "\n")
    fdistrib.close()

def getBorderingGenes(rgp, multigenics):
        border = [None, None]
        pos = rgp.genes[-1].position
        while border[0] is None:
            curr_fam = None
            if pos == 0:
                if rgp.contig.is_circular:
                    curr_fam = rgp.contig.genes[-1].family
                else:
                    border[0] = ""
            else:
                curr_fam = rgp.contig.genes[pos -1].family
            if curr_fam is not None and curr_fam not in multigenics and curr_fam.namedPartition == "persistent":
                border[0] = curr_fam
            else:
                pos -= 1

        pos = rgp.genes[0].position
        while border[1] is None:
            curr_fam = None
            if pos == len(rgp.contig.genes)-1:
                if rgp.contig.is_circular:
                    curr_fam = rgp.contig.genes[0].family
                else:
                    border[1] = ""
            else:
                curr_fam = rgp.contig.genes[pos+1].family
            if curr_fam is not None and curr_fam not in multigenics:
                border[1] = curr_fam
            else:
                pos+=1
        return border

def getUniqRGP(rgpList):
        uniqRGP = set()
        for rgp in rgpList:
            z = True
            for seenRgp in uniqRGP:
                if rgp == seenRgp:
                    z = False
                    break
            if z:
                uniqRGP.add(rgp)
        return uniqRGP

def checkSim(pairKnownBorder, pairBorderGenes, overlapping_match, exact_match, set_size):
    """ Checks if the two pairs of 'exact_match' first gene families are identical, or eventually if they overlap in an ordered way at least 'overlapping_match'"""
    kbpair = [False, False]
    bpair = [False, False]
    for countkb, kb in enumerate(pairKnownBorder):
        for countb, b in enumerate(pairBorderGenes):
            match = False
            if b[0:exact_match] == kb[0:exact_match]:
                match = True
            elif len(b) == set_size and len(kb) == set_size:
                if not match:
                    for ikb in range(1, set_size-overlapping_match+1):
                        if b[0:len(kb[ikb:])] == kb[ikb:]:
                            match = True
                            break
                if not match:
                    for ib in range(1, set_size-overlapping_match+1):
                        if b[ib:] == kb[0:len(b[ib:])]:
                            match = True
                            break
            if match:
                kbpair[countkb] = True
                bpair[countb] = True

    if kbpair[0] and kbpair[1] and bpair[0] and bpair[1]:
        return True
    return False

def makeSpotGraph(rgps, multigenics, output, spot_graph, overlapping_match, set_size, exact_match):
    def addNewNode(g, rgp, borders):
        blocks = str(sorted([[fam.ID for fam in borders[0]],[fam.ID for fam in borders[1]]], key = lambda x : x[0]))
        g.add_node(blocks)
        try:
            g.nodes[blocks]["nb_rgp"]+=1
            g.nodes[blocks]["rgp"].add(rgp)
        except KeyError:
            g.nodes[blocks]["nb_rgp"] = 1
            g.nodes[blocks]["border1"] = borders[1]
            g.nodes[blocks]["border0"] = borders[0]
            g.nodes[blocks]["rgp"] = set([rgp])

    spotGraph = nx.Graph()
    lost = 0
    for rgp in rgps:
        border = rgp.getBorderingGenes(set_size, multigenics)
        if len(border[0]) < set_size or len(border[1]) < set_size:
            lost+=1
        else:
            addNewNode(spotGraph, rgp, border)
    logging.getLogger().info(f"{lost} RGPs were not used as they are on a contig border (or have less than {set_size} persistent gene families until the contig border)")
    nodeList = list(spotGraph.nodes)
    logging.getLogger().info(f"{len(nodeList)} number of different pairs of flanking gene families")
    for i, nodei in enumerate(nodeList[:-1]):
        for nodej in nodeList[i+1:]:
            nodeObji = spotGraph.nodes[nodei]
            nodeObjj = spotGraph.nodes[nodej]
            if checkSim([nodeObji["border0"], nodeObji["border1"]], [nodeObjj["border0"], nodeObjj["border1"]], overlapping_match, exact_match, set_size):
                spotGraph.add_edge(nodei, nodej)

    spots = []
    for comp in nx.algorithms.components.connected_components(spotGraph):
        spots.append([ [], set() ])
        for node in comp:
            spots[-1][1] |= spotGraph.nodes[node]["rgp"]
            spots[-1][0].append([spotGraph.nodes[node]["border1"], spotGraph.nodes[node]["border0"]])

    if spot_graph:
        for node in spotGraph.nodes:
            del spotGraph.nodes[node]["border0"]
            del spotGraph.nodes[node]["border1"]
            del spotGraph.nodes[node]["rgp"]

        nx.readwrite.gexf.write_gexf(spotGraph, output + "/spotGraph.gexf")
    return spots

def makeFlanking(spots, output):
    flankGraph = nx.Graph()
    c = 0
    for borders, rgps in spots:
        flankGraph.add_node(c)
        flankGraph.nodes[c]["nb_rgp"] = len(rgps)
        flankGraph.nodes[c]["nb_organisations"] = len(getUniqRGP(rgps))
        bords = set()
        for border in borders:
            bords.add(frozenset(border[0]))
            bords.add(frozenset(border[1]))
        flankGraph.nodes[c]["borders"] = frozenset(bords)
        c+=1

    nodeList = list(flankGraph.nodes)
    for i, nodei in enumerate(nodeList[:-1]):
        for nodej in nodeList[i+1:]:
            inter = len(flankGraph.nodes[nodei]["borders"] & flankGraph.nodes[nodej]["borders"])
            if inter != 0:
                flankGraph.add_edge(nodei, nodej, nb_bord=inter)
    for node in nodeList:
        del flankGraph.nodes[node]["borders"]
    nx.readwrite.gexf.write_gexf(flankGraph, output + "/flankGraph.gexf")

def detect_hotspots(pangenome, multigenics, output, spot_graph = False, flanking_graph = False, overlapping_match = 2, set_size = 3, exact_match = 1):
    logging.getLogger().info("Detecting hotspots in the pangenome...")

    spots = makeSpotGraph(pangenome.regions, multigenics, output, spot_graph, overlapping_match, set_size, exact_match)

    nb_above1perc = len([ rgps for _, rgps in spots if len(rgps) > len(pangenome.organisms) * 0.05])
    logging.getLogger().info(f"There are {len(spots)} spots in this pangenome, and {nb_above1perc} of them have RGPs in more than 5% of the organisms.")

    if flanking_graph:
        makeFlanking(spots, output)

    spot_rgps = [ set(rgps) for _, rgps in spots]#list of rgps, regrouped by spots.
    spot_distribution(spot_rgps, output)

    return spots

def predictRGP(pangenome, output, persistent_penalty = 3, variable_gain = 1, min_length = 3000, min_score = 4, dup_margin = 0.05, spot_graph = False,flanking_graph = False,overlapping_match = 2, set_size = 3, exact_match = 1, cpu = 1):
    #check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions = True)

    logging.getLogger().info("Detecting multigenic families...")
    multigenics = get_multigenics(pangenome, dup_margin)
    pangenomeRGP = set()
    logging.getLogger().info("Compute Regions of Genomic Plasticity ...")
    bar = tqdm(pangenome.organisms, unit = "genomes")
    for org in bar:
        pangenomeRGP |= compute_org_rgp(org, persistent_penalty, variable_gain, min_length, min_score, multigenics)
    logging.getLogger().info(f"Predicted {len(pangenomeRGP)} RGP")
    pangenome.addRegions(pangenomeRGP)

    detect_hotspots(pangenome, multigenics, output, spot_graph, flanking_graph,overlapping_match, set_size, exact_match)

    #save parameters and save status
    pangenome.parameters["RGP"] = {}
    pangenome.parameters["RGP"]["persistent_penalty"] = persistent_penalty
    pangenome.parameters["RGP"]["variable_gain"] = variable_gain
    pangenome.parameters["RGP"]["min_length"] = min_length
    pangenome.parameters["RGP"]["min_score"] = min_score
    pangenome.parameters["RGP"]["dup_margin"] = dup_margin
    pangenome.status['predictedRGP'] = "Computed"

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    # if args.flanking_graph or args.spot_graph:
    mkOutdir(args.output, args.force)
    predictRGP(pangenome, args.output, args.persistent_penalty, args.variable_gain, args.min_length, args.min_score, args.dup_margin, args.spot_graph, args.flanking_graph, args.overlapping_match_size_hotspot, args.set_size_hotspot, args.exact_match_size_hotspot, args.cpu)

def rgpSubparser(subparser):
    parser = subparser.add_parser("rgp", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--persistent_penalty', required=False, type=int, default=3, help="Penalty score to apply to persistent genes")
    optional.add_argument('--variable_gain', required=False, type=int, default=1, help="Gain score to apply to variable genes")
    optional.add_argument('--min_score', required=False, type=int, default=4, help="Minimal score wanted for considering a region as being a RGP")
    optional.add_argument('--min_length', required=False, type=int, default=3000, help="Minimum length (bp) of a region to be considered a RGP")
    optional.add_argument("--dup_margin", required = False, type=int, default=0.05, help="Minimum ratio of organisms where the family is present in which the family must have multiple genes for it to be considered 'duplicated'" )
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--spot_graph", required = False, action="store_true", help = "Writes a graph in a .gexf format of pairs of blocks of single copy markers flanking RGPs, supposedly belonging to the same hotspot")
    optional.add_argument("--flanking_graph", required = False, action="store_true", help = "Writes a graph in a .gexf format of common blocks of single copy markers flanking RGPs, supposedly with some common origin")
    optional.add_argument("--overlapping_match_size_hotspot", required=False, type=int, default = 2, help="The number of 'missing' persistent genes allowed when comparing flanking genes during hotspot computations")
    optional.add_argument("--set_size_hotspot", required = False, type = int, default = 3, help = "Number of single copy markers to use as flanking genes for a RGP during hotspot computation")
    optional.add_argument("--exact_match_size_hotspot", required = False, type= int, default = 1, help = "Number of perfecty matching flanking single copy markers required to associate RGPs during hotspot computation (Ex: If set to 1, two RGPs are in the same hotspot if both their 1st flanking genes are the same)")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser