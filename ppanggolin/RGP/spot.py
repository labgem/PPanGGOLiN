#!/usr/bin/env python3


#default libraries
import logging
import argparse
import time
import os
from collections import defaultdict, Counter
import random
from operator import attrgetter
from statistics import mean, stdev
from multiprocessing import Pool
from random import randint, shuffle

#installed libraries
from tqdm import tqdm
import networkx as nx

import rpy2
from rpy2 import robjects
from rpy2.robjects.packages import importr

from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram
from numpy import percentile

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region, Spot
from ppanggolin.formats import checkPangenomeInfo, writePangenome, ErasePangenome
from ppanggolin.utils import mkOutdir, jaccard_similarities

def compBorder(border1, border2, overlapping_match, exact_match, set_size):
    if border1[0:exact_match] == border2[0:exact_match]:
        return True
    elif len(border1) == set_size and len(border2) == set_size:
        for ikb in range(1, set_size-overlapping_match+1):
            if border1[0:len(border2[ikb:])] == border2[ikb:]:
                return True
        for ib in range(1, set_size-overlapping_match+1):
            if border1[ib:] == border2[0:len(border1[ib:])]:
                return True
    return False

def checkSim(pairBorder1, pairBorder2, overlapping_match, exact_match, set_size):
    """ Checks if the two pairs of 'exact_match' first gene families are identical, or eventually if they overlap in an ordered way at least 'overlapping_match'"""
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
        blocks = str(sorted([[gene.family.ID for gene in borders[0]],[gene.family.ID for gene in borders[1]]], key = lambda x : x[0]))
        g.add_node(blocks)
        try:
            g.nodes[blocks]["nb_rgp"]+=1
            g.nodes[blocks]["rgp"].add(rgp)
        except KeyError:
            g.nodes[blocks]["nb_rgp"] = 1
            g.nodes[blocks]["border1"] = [gene.family for gene in borders[1]]
            g.nodes[blocks]["border0"] = [gene.family for gene in borders[0]]
            g.nodes[blocks]["rgp"] = set([rgp])

    spotGraph = nx.Graph()
    lost = 0
    used = 0
    for rgp in rgps:
        border = rgp.getBorderingGenes(set_size, multigenics)
        if len(border[0]) < set_size or len(border[1]) < set_size:
            lost+=1
        else:
            used+=1
            addNewNode(spotGraph, rgp, border)
    logging.getLogger().info(f"{lost} RGPs were not used as they are on a contig border (or have less than {set_size} persistent gene families until the contig border)")
    logging.getLogger().info(f"{used} RGPs are being used to predict spots of insertion")
    nodeList = list(spotGraph.nodes)
    logging.getLogger().info(f"{len(nodeList)} number of different pairs of flanking gene families")
    for i, nodei in enumerate(nodeList[:-1]):
        for nodej in nodeList[i+1:]:
            nodeObji = spotGraph.nodes[nodei]
            nodeObjj = spotGraph.nodes[nodej]
            if checkSim([nodeObji["border0"], nodeObji["border1"]], [nodeObjj["border0"], nodeObjj["border1"]], overlapping_match, exact_match, set_size):
                spotGraph.add_edge(nodei, nodej)

    spots = []
    spot_id = 0
    for comp in nx.algorithms.components.connected_components(spotGraph):
        curr_spot = Spot(spot_id)
        spots.append(curr_spot)
        for node in comp:
            curr_spot.addRegions(spotGraph.nodes[node]["rgp"])
        spot_id+=1

    if spot_graph:
        for node in spotGraph.nodes:
            del spotGraph.nodes[node]["border0"]
            del spotGraph.nodes[node]["border1"]
            del spotGraph.nodes[node]["rgp"]

        nx.readwrite.gexf.write_gexf(spotGraph, output + "/spotGraph.gexf")
    return spots

def checkPangenomeFormerSpots(pangenome, force):
    """ checks pangenome status and .h5 files for former spots, delete them if allowed or raise an error """
    if pangenome.status["spots"] == "inFile" and force == False:
        raise Exception("You are trying to detect spots on a pangenome which already has predicted spots. If you REALLY want to do that, use --force (it will erase spots previously predicted).")
    elif pangenome.status["spots"] == "inFile" and force == True:
        ErasePangenome(pangenome, spots = True)


def predictHotspots(pangenome, output, force=False, cpu = 1, spot_graph = False, overlapping_match = 2, set_size = 3, exact_match = 1, draw_hotspot = False, interest = "", priority="name,ID", show_bar=True):
    #check that given parameters for hotspot computation make sense
    checkParameterLogic(overlapping_match, set_size, exact_match, priority)
    #check for formerly computed stuff, and erase if allowed
    checkPangenomeFormerSpots(pangenome, force)
    #check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions = True, needRGP=True, show_bar=show_bar)

    #get multigenic gene families
    logging.getLogger().info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])
    
    logging.getLogger().info("Detecting hotspots in the pangenome...")

    #predict spots
    spots = makeSpotGraph(pangenome.regions, multigenics, output, spot_graph, overlapping_match, set_size, exact_match)

    if len(spots) == 0:
        logging.getLogger().warning("No spots were detected.")
    else:
        logging.getLogger().info(f"{len(spots)} spots were detected")
    #define elements of interest (e.g. gene name, product substring) to search in gene annotations
    if interest != "":
        elements = [ el.strip() for el in interest.split(',') ]
    else:
        elements = []
    #draw spots of interest
    if draw_hotspot:
        drawn_spots = select_spots(pangenome, spots, elements)
        if len(drawn_spots)>0:
            logging.getLogger().info(f"Drawing {len(drawn_spots)} spots...")
            draw_spots(drawn_spots, output, cpu, overlapping_match, exact_match, set_size, multigenics, elements, priority, show_bar=show_bar)
        else:
            logging.getLogger().warning("No spot has enough variability to be drawn in a meaningful figure")
    pangenome.addSpots(spots)
    pangenome.status["spots"] = "Computed"
    pangenome.parameters["spots"] = {}
    pangenome.parameters["spots"]["set_size"] = set_size
    pangenome.parameters["spots"]["overlapping_match"] = overlapping_match
    pangenome.parameters["spots"]["exact_match"] = exact_match

def subgraph(spot, output, filename, with_border=True, set_size=3, multigenics = None):
        """ write a pangenome subgraph of the gene families of a spot in gexf format"""
        g = nx.Graph()
        for rgp in spot.regions:
            if with_border:
                borders = rgp.getBorderingGenes(set_size, multigenics)
                minpos = min([ gene.position for border in borders for gene in border ])
                maxpos = max([ gene.position for border in borders for gene in border ])
            else:
                minpos = rgp.startGene.position
                maxpos = rgp.stopGene.position
            GeneList = rgp.contig.genes[minpos:maxpos+1]
            prev = None
            for gene in GeneList:
                g.add_node(gene.family.name, partition = gene.family.namedPartition)
                try:
                    g.nodes[gene.family.name]["occurrence"] += 1
                except KeyError:
                    g.nodes[gene.family.name]["occurrence"] = 1
                if gene.name != "":
                    if "name" in g.nodes[gene.family.name]:
                        try:
                            g.nodes[gene.family.name]["name"][gene.name] +=1
                        except KeyError:
                            g.nodes[gene.family.name]["name"][gene.name] = 1
                    else:
                        g.nodes[gene.family.name]["name"] = Counter([gene.name])
                if prev is not None:
                    g.add_edge(gene.family.name, prev)
                    try:
                        g[gene.family.name][prev]["rgp"].add(rgp)
                    except KeyError:
                        g[gene.family.name][prev]["rgp"] = set(rgp)
                prev = gene.family.name
        for node1, node2 in g.edges:
            g[node1][node2]["weight"] = len(g[node1][node2]["rgp"]) / len(spot.regions)
            del g[node1][node2]["rgp"]
        for node in g.nodes:
            if "name" in g.nodes[node]:
                g.nodes[node]["name"] = g.nodes[node]["name"].most_common(1)[0][0]

        nx.write_gexf(g, output+"/"+filename)

#useless atm, but could be useful for the futur should we want to study spots cross-species.
def writeBorders_spots(spots, pangenome, output):
    fout = open(output + "/spots_families.faa","w")
    n_spot = 0

    for spot in spots:
        n_border_group = 0
        for borders in spot[0]:
            n_border = 0
            prevalence = len(spot.regions)
            organisations = len(spot.getUniqSynteny())
            for border in borders:
                order=0
                for fam in border:
                    fout.write(f">{fam.name}_prevalence-{prevalence}_synt-{organisations}_spot-{n_spot}_group-{n_border_group}_border-{n_border}_order-{order}\n")
                    fout.write(f"{fam.sequence}\n")
                    order+=1
                n_border+=1
            n_border_group+=1
        n_spot+=1
    fout.close()

def select_spots(pangenome, spots, elements, min_uniq_orders = 2):
    if min_uniq_orders <= 1:
        raise Exception("Can't draw spots with only a single (or less) gene order")
    to_draw= []
    for spot in spots:
        nb_uniq = len(spot.getUniqOrderedSet())
        if nb_uniq >= min_uniq_orders:
            to_draw.append(spot)
    return to_draw

def checkParameterLogic(overlapping_match, set_size, exact_match, priority):
    if overlapping_match >= set_size:
        raise Exception(f'--overlapping_match_hotspot ({overlapping_match}) cannot be bigger than (or equal to) --set_size_hotspot ({set_size})')
    if exact_match > set_size:
        raise Exception(f'--exact_match_size_hotspot ({exact_match}) cannot be bigger than --set_size_hotspot ({set_size})')

    for p in priority.split(','):
        if p.lower() not in ["name","id","family"]:
            raise Exception(f"You have indicated a label which is not supported with --label_priority. You indicated '{p}'. Supported labels are name, id and family")

def makeColorsForFams(fams):
    famcol = {}
    for fam in fams:
        col =  list(random.choices(range(256), k=3))
        famcol[fam] = '#%02x%02x%02x' % (col[0], col[1], col[2])
    return famcol

def orderGeneLists(geneLists, ordered_counts, overlapping_match, exact_match, set_size):
    geneLists = lineOrderGeneLists(geneLists, overlapping_match, exact_match, set_size)
    return rowOrderGeneLists(geneLists, ordered_counts)

def rowOrderGeneLists(geneLists, ordered_counts):
    famDict = defaultdict(set)

    for index, genelist in enumerate([genelist[0] for genelist in geneLists]):
        for gene in genelist:
            if hasattr(gene,"family"):
                famDict[gene.family].add(index)
    all_indexes = []
    all_columns = []
    data = []
    for famIndex, RGPindexes in enumerate(famDict.values()):
        all_indexes.extend([famIndex] * len(RGPindexes))
        all_columns.extend(RGPindexes)
        data.extend([1.0]*len(RGPindexes))

    mat_p_a = csc_matrix((data, (all_indexes,all_columns)), shape = (len(famDict),len(geneLists)), dtype='float')
    dist    = pdist(1 - jaccard_similarities(mat_p_a,0).todense())
    hc      = linkage(dist, 'single')

    dendro = dendrogram(hc,no_plot=True)

    new_geneLists = [ geneLists[index] for index in dendro["leaves"]]
    new_ordered_counts = [ ordered_counts[index] for index in dendro["leaves"] ]
    return new_geneLists, new_ordered_counts

def lineOrderGeneLists(geneLists, overlapping_match, exact_match, set_size):
    classified = set([0])#first gene list has the right order
    new_classify = set()
    
    nbloop = 0
    to_classify = set(range(1, len(geneLists)))#the others may (or may not) have it
    while len(to_classify) != 0:
        for classIndex in classified:
            base_border1 = [ gene.family for gene in geneLists[classIndex][1][0] ]
            base_border2 = [ gene.family for gene in geneLists[classIndex][1][1] ]
            for unclassIndex in list(to_classify):
                border1 = [ gene.family for gene in geneLists[unclassIndex][1][0] ]
                border2 = [ gene.family for gene in geneLists[unclassIndex][1][1] ]
                if compBorder(base_border1, border1, overlapping_match, exact_match, set_size) or compBorder(base_border2, border2, overlapping_match, exact_match, set_size):
                    to_classify.discard(unclassIndex)
                    new_classify.add(unclassIndex)
                elif compBorder(base_border2, border1, overlapping_match, exact_match, set_size) or compBorder(base_border1, border2, overlapping_match, exact_match, set_size):
                    geneLists[unclassIndex][0] = geneLists[unclassIndex][0][::-1]#reverse the order of the genes to match the 'reference'
                    to_classify.discard(unclassIndex)
                    new_classify.add(unclassIndex)
        classified = new_classify#the newly classified will help to check the unclassified, the formerly classified are not useful for what remains (if sth remains)
        new_classify = set()
        nbloop+=1
        if nbloop>= 10:
            print("infinit loop !")
            print( [ gene.family.ID for border in geneLists[0][1] for gene in border ])
            print("#######")
            for iden in to_classify:
                print( [ gene.family.ID for border in geneLists[iden][1] for gene in border  ])
            print("####")
            for iden in classified:
                print( [ gene.family.ID for border in geneLists[iden][1] for gene in border  ])
            raise Exception()
    return geneLists

def defineElementsOfInterest(genelist, elements):
    present_EOI = set()
    for gene in genelist:
        if gene.name in elements:
            present_EOI.add(gene.name)
        for el in elements:
            if el in gene.product:
                present_EOI.add(el)
    return sorted(present_EOI)#sort EOI so that they are in the same order between different runs of the same analysis

def drawCurrSpot(genelists, ordered_counts, elements, famCol, filename, priority):
    rdframes = []
    annotList = []
    partitionColors = {"shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF"}
    importr("genoPlotR")
    dna_seg = robjects.r["dna_seg"]
    annotation = robjects.r["annotation"]
    middle = robjects.r["middle"]

    longest_gene_list = 0

    present_EOI = set()

    for index, GeneList in enumerate(genelists):
        genelist = GeneList[0]
        present_EOI |= defineElementsOfInterest(genelist, elements)
        if len(genelist) > longest_gene_list:
            longest_gene_list = len(genelist)
        if genelist[0].start < genelist[1].start:
            ordered=True
            start = genelist[0].start
        else:
            ordered = False
            start = genelist[0].stop
        df = {'name':[],'start':[],'end':[],'strand':[],'col':[], 'fill':[], "gene_type":[]}
        gene_names = []
        
        for gene in genelist:
            if 'RNA' in gene.type:
                gene_names.append(' ' + gene.product)
            else:
                for p in priority.split(','):
                    if p == 'name':
                        if gene.name != "":
                            gene_names.append(' ' + gene.name)
                            break
                    elif p == "ID":
                        if gene.local_identifier != "":
                            gene_names.append(' ' + gene.local_identifier)
                        else:
                            gene_names.append(' ' + str(gene.ID))
                        break#necessarily has an ID
                    elif p == 'family':
                        gene_names.append(' ' + str(gene.family.name))
                        break#necessarily has a family
            df['name'].append(gene.ID)
            if ordered:
                if gene.strand == "+":
                    df['start'].append(gene.start - start)
                    df['end'].append(gene.stop - start)
                    df['strand'].append(gene.strand)
                else:
                    df["end"].append(gene.start - start)
                    df["start"].append(gene.stop - start)
                    df["strand"].append(gene.strand)
            else:
                if gene.strand == "+":
                    df["start"].append(abs(gene.stop - start))
                    df['end'].append(abs(gene.start - start))
                    df["strand"].append("-")#we invert it because we reversed the order !
                else:
                    df["start"].append(abs(gene.start - start))
                    df["end"].append(abs(gene.stop - start))
                    df["strand"].append("+")#we invert it because we reversed the order !
            if gene.type == "CDS":
                df['col'].append(partitionColors[gene.family.namedPartition])
                df['fill'].append(famCol[gene.family])
                df["gene_type"].append("side_blocks")
            elif gene.type == "tRNA":
                df['col'].append('#000000')
                df['gene_type'].append("bars")
                df['fill'].append('#000000')
            else:
                df['col'].append('#000000')
                df['gene_type'].append("headless_arrows")
                df['fill'].append('#000000')
        df["name"] = robjects.StrVector(df["name"])
        df["start"] = robjects.IntVector(df["start"])
        df["end"] = robjects.IntVector(df["end"])
        df["strand"] = robjects.StrVector(df["strand"])
        df["col"] = robjects.StrVector(df["col"])
        df["fill"] = robjects.StrVector(df["fill"])
        df["gene_type"] = robjects.StrVector(df["gene_type"])
        dnasegObj = dna_seg(robjects.DataFrame(df))
        annot = annotation(x1 = middle(dnasegObj), text=robjects.StrVector(gene_names), rot = 20)
        annotList.append((f'{GeneList[2].organism.name}, x'+str(ordered_counts[index]), annot))
        rdframes.append((f'{GeneList[2].organism.name}, x'+str(ordered_counts[index]), dnasegObj))
    filename = filename +('_' + "_".join(present_EOI) if len(present_EOI) > 0 else "") + ".png"
    Rannot = robjects.ListVector(annotList)
    Rdna_segs = robjects.ListVector(rdframes)
    return Rdna_segs, Rannot, rdframes, longest_gene_list, filename

def _spotDrawing(Rdna_segs, Rannot, rdframes, longest_gene_list, filename):
    try:
        plot_gene_map = robjects.r["plot_gene_map"]
        grdevices = importr('grDevices')
        grdevices.png(file=filename, width = longest_gene_list * 70, height= len(rdframes) * 60)#pylint: disable=no-member
        plot_gene_map(dna_segs = Rdna_segs, annotations=Rannot, lwd = 4)
        grdevices.dev_off()#pylint: disable=no-member
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        logging.getLogger().warning(f"{os.path.basename(filename)} cannot be drawn as the spot is probably too large to be rendered")

def draw_spots(spots, output, cpu, overlapping_match, exact_match, set_size, multigenics, elements, priority, show_bar=False):
    logging.getLogger().info("Selecting and ordering genes among regions...")
    bar = tqdm(range(len(spots)), unit = "spot", disable = not show_bar)
    spots_to_draw = []
    for spot in spots:
        uniqRGPS = frozenset(spot.getUniqOrderedSet())
        Fams = set()
        GeneLists = []

        countUniq = spot.countUniqOrderedSet()

        #order unique rgps by occurrences
        sortedUniqRGPs = sorted(uniqRGPS, key = lambda x : countUniq[x], reverse=True)
        for rgp in sortedUniqRGPs:
            borders = rgp.getBorderingGenes(set_size, multigenics)
            minpos = min([ gene.position for border in borders for gene in border ])
            maxpos = max([ gene.position for border in borders for gene in border ])
            GeneList = rgp.contig.genes[minpos:maxpos+1]
            minstart = min([ gene.start for border in borders for gene in border ])
            maxstop = max([ gene.stop for border in borders for gene in border ])
            RNAstoadd = set()
            for rna in rgp.contig.RNAs:
                if rna.start > minstart and rna.start < maxstop:
                    RNAstoadd.add(rna)
            GeneList.extend(RNAstoadd)
            GeneList = sorted(GeneList, key = lambda x : x.start)

            Fams |= { gene.family for gene in GeneList if gene.type == "CDS"}

            GeneLists.append([GeneList, borders, rgp])
        famcol = makeColorsForFams(Fams)
        ordered_counts = sorted(countUniq.values(), reverse = True)
        GeneLists, ordered_counts = orderGeneLists(GeneLists, ordered_counts, overlapping_match, exact_match, set_size)
        fname = output + '/spot_' + str(spot.ID)
        # spots_to_draw.append((GeneLists, ordered_counts, elements, famcol, fname))
        spots_to_draw.append(drawCurrSpot(GeneLists, ordered_counts, elements, famcol, fname, priority))#make R dataframes, and plot them using genoPlotR.
        bar.update()
    logging.getLogger().info("Drawing spots...")
    bar = tqdm(spots_to_draw, unit = "spot drawn", disable=not show_bar)
    with Pool(cpu) as p:
        p.starmap(_spotDrawing, bar)
    bar.close()

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.spot_graph or args.draw_hotspots:
        mkOutdir(args.output, args.force)
    predictHotspots(pangenome, args.output, force=args.force, cpu = args.cpu, spot_graph=args.spot_graph, overlapping_match=args.overlapping_match, set_size=args.set_size, exact_match=args.exact_match_size, draw_hotspot=args.draw_hotspots, interest=args.interest, priority=args.label_priority, show_bar=args.show_prog_bars)
    writePangenome(pangenome, pangenome.file, args.force, show_bar=args.show_prog_bars)


def spotSubparser(subparser):
    parser = subparser.add_parser("spot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--spot_graph", required = False, action="store_true", help = "Writes a graph in a .gexf format of pairs of blocks of single copy markers flanking RGPs, supposedly belonging to the same hotspot")
    optional.add_argument("--draw_hotspots", required=False, action="store_true", help = "Draws a figure representing all of the hotspots syntenies")
    optional.add_argument("--overlapping_match", required=False, type=int, default = 2, help="The number of 'missing' persistent genes allowed when comparing flanking genes during hotspot computations")
    optional.add_argument("--set_size", required = False, type = int, default = 3, help = "Number of single copy markers to use as flanking genes for a RGP during hotspot computation")
    optional.add_argument("--exact_match_size", required = False, type= int, default = 1, help = "Number of perfecty matching flanking single copy markers required to associate RGPs during hotspot computation (Ex: If set to 1, two RGPs are in the same hotspot if both their 1st flanking genes are the same)")
    optional.add_argument("--interest",required=False, type=str, default="tRNA,integrase",help = "Comma separated list of elements to indicate in figure file names if present whenD drawing hotspots (any text which can be found as 'gene name' or 'product' if annotations were provided). Defaults are 'tRNA' and 'integrase' which are often found near the most dynamic spots of insertion")
    optional.add_argument("--label_priority", required=False, type=str, default='name,ID', help = "Option to use with --draw_hotspots. Will indicate what to write in the figure labels as a comma-separated list. Order gives priority. Possible values are: name (for gene names), ID (for the gene IDs), family (for the family IDs)")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser