#!/usr/bin/env python3


#default libraries
import logging
import argparse
import time
import os
from collections import defaultdict, Counter
import random
from operator import attrgetter
from statistics import mean, variance
from random import randint, shuffle

#installed libraries
from tqdm import tqdm, trange
import networkx as nx
from rpy2 import robjects
from rpy2.robjects.packages import importr

from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram
import colorlover as cl
from numpy import percentile

#local libraries
from ppanggolin.pangenome import Pangenome, Region
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir, jaccard_similarities


def spot_distribution(spots, pangenome, output):
    """takes in spots are a list of sets of rgps"""
    fdistrib = open(output + "/spot_rgp_distribution.tsv","w")
    for rgps in spots:
        fdistrib.write(str(len(rgps)) + "\t" + str(len(getUniqRGP(rgps))) + "\t" + str(round(len(rgps) / len(pangenome.organisms),2)) +"\n")
    fdistrib.close()

def countUniqContent(rgpList):
    uniqRGP = Counter()
    for rgp in rgpList:
        z = True
        for seenRgp in uniqRGP:
            if rgp.families == seenRgp.families:
                z = False
                uniqRGP[seenRgp] +=1
        if z:
            uniqRGP[rgp]+=1
    return uniqRGP

def countUniqRGP(rgpList):
    uniqRGP = Counter()
    for rgp in rgpList:
        z = True
        for seenRgp in uniqRGP:
            if rgp == seenRgp:
                z = False
                uniqRGP[seenRgp]+=1
                break
        if z:
            uniqRGP[rgp] +=1
    return uniqRGP

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

def checkSim(pairKnownBorder, pairBorderGenes, overlapping_match, exact_match, set_size):
    """ Checks if the two pairs of 'exact_match' first gene families are identical, or eventually if they overlap in an ordered way at least 'overlapping_match'"""
    kbpair = [False, False]
    bpair = [False, False]
    for countkb, kb in enumerate(pairKnownBorder):
        for countb, b in enumerate(pairBorderGenes):
            if compBorder(b, kb, overlapping_match, exact_match, set_size):
                kbpair[countkb] = True
                bpair[countb] = True

    if kbpair[0] and kbpair[1] and bpair[0] and bpair[1]:
        return True
    return False

def makeSpotGraph(rgps, multigenics, output, spot_graph, overlapping_match, set_size, exact_match):
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

def predictHotspots(pangenome, output, spot_graph = False, flanking_graph = False, overlapping_match = 2, set_size = 3, exact_match = 1, draw_hotspot = False, interest = "", write_spots = True):
    
    #check that given parameters for hotspot computation make sense
    checkParameterLogic(overlapping_match, set_size, exact_match)

    #check statuses and load info
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=False, needPartitions = True, needRGP=True)

    logging.getLogger().info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(pangenome.parameters["RGP"]["dup_margin"])
    
    logging.getLogger().info("Detecting hotspots in the pangenome...")

    spots = makeSpotGraph(pangenome.regions, multigenics, output, spot_graph, overlapping_match, set_size, exact_match)

    if flanking_graph:
        makeFlanking(spots, output)

    # spot_rgps = [ set(rgps) for _, rgps in spots]#list of rgps, regrouped by spots.
    #spot_distribution(spot_rgps, pangenome, output)
    #ADD : Associate rgps on contig borders to known spots. (using both gene content and bordering gene families?)
    if interest != "":
        elements = [ el.strip() for el in interest.split(',') ]
    else:
        elements = []

    if draw_hotspot:
        drawn_spots = select_spots(pangenome, spots, elements)
        if len(drawn_spots)>0:
            draw_spots(drawn_spots, output, overlapping_match, exact_match, set_size, multigenics, elements)#TODO: add a parameter to control how much presence is needed for a 'hotspot'

    all_spot_fams = set()
    for spot in spots:
        for rgp in spot[1]:
            all_spot_fams |= rgp.families
    maxHTg = uniform_spots(len(spots),len(all_spot_fams) )

    if write_spots:
        writeSpots(spots, output, elements)
        summarize_spots(spots, output, maxHTg)
        #writeBorders_spots(spots,pangenome, output)
    return spots

def uniform_spots(S, N):
    logging.getLogger().info(f"There are {N} gene families among {S} spots")
    el1 = int(N / 3)#number of 1 gene elements
    el3 = int(((2*N) / 3)/3)#number of 3 gene elements
    logging.getLogger().info(f"There are {el1} elements of 1 gene and {el3} elements of 3 genes")
    # p_el = 1 / S#proba of an element to be placed in a given spot
    maxHTdist = []

    els = [1] * el1 + [3] * el3
    for _ in trange(1000, unit = "sample"):
        th_spots = Counter()
        shuffle(els)
        for i in els:
            th_spots[randint(1,S)]+=i
        maxHTdist.append(th_spots.most_common(1)[0][1])
    return percentile(maxHTdist, 95)

def summarize_spots(spots, output, nbFamLimit):
    fout = open(output + "/summarize_spots.tsv","w")
    fout.write("spot\tsorensen\tturnover\tnestedness\tnb_rgp\tnb_families\tnb_organisations\tnb_content\tmean_spot_fluidity\tvar_spot_fluidity\tmean_nb_genes\tvar_nb_genes\tstatus\n")
    logging.getLogger().info("Computing sorensen, turnover and nestedness indexes for spots with more than 1 rgp...")
    n_spot = 0
    bar = tqdm(spots, unit = "spot")#can multi
    for spot in bar:
        if len(spot[1]) > 1:
            tot_fams = set()
            summin = 0
            spot_fluidity=[]
            summax = 0
            rgp_list = list(spot[1])
            len_uniq_content = len(countUniqContent(spot[1]))
            uniq_dic = countUniqRGP(spot[1])
            uniq_list = list(uniq_dic.keys())
            nbuniq_organizations = len(uniq_list)
            size_list = []
            for i, rgp in enumerate(uniq_list[:-1]):
                tot_fams |= rgp.families
                size_list.append(len(rgp.genes))
                spot_fluidity.extend([0] * uniq_dic[rgp] )
                for rgp2 in uniq_list[i+1:]:
                    interfams = set(rgp.families) & set(rgp2.families)
                    spot_fluidity.extend([(((len(rgp.families) - len(interfams)) + (len(rgp2.families) - len(interfams))) / (len(rgp.families) + len(rgp2.families)))] * (uniq_dic[rgp] * uniq_dic[rgp2]))
                    summin += min(len(rgp.families) - len(interfams), len(rgp2.families) - len(interfams)) * uniq_dic[rgp] * uniq_dic[rgp2]
                    summax += max(len(rgp.families) - len(interfams), len(rgp2.families) - len(interfams)) * uniq_dic[rgp] * uniq_dic[rgp2]
            tot_fams |= rgp_list[-1].families
            size_list.extend([len(rgp_list[-1].genes)] * uniq_dic[uniq_list[-1]] )
            spot_fluidity.extend([0] * uniq_dic[uniq_list[-1]] )
            mean_size = mean(size_list)
            var_size = variance(size_list)
            mean_spot_fluidity = mean(spot_fluidity)
            var_spot_fluidity = variance(spot_fluidity)
            sumSiSt = sum([ len(rgp.families) for rgp in rgp_list ])-len(tot_fams)
            sorensen = (summin + summax) / (2*sumSiSt + summin + summax )
            turnover = summin / (sumSiSt + summin)
            nestedness = sorensen - turnover
            status = "hotspot" if nbFamLimit <= len(tot_fams) else "coldspot"
            fout.write("\t".join(map(str,[f"spot_{n_spot}", sorensen, turnover, nestedness, len(rgp_list), len(tot_fams), nbuniq_organizations, len_uniq_content, mean_spot_fluidity,var_spot_fluidity, mean_size,var_size, status])) + "\n")
        n_spot+=1
    bar.update()
    fout.close()

def writeSpots(spots, output, elements):
    fout = open(output + "/spots.tsv","w")
    fout.write("spot_id\trgp_id\tinterest\n")
    n_spot = 0
    for spot in spots:
        for rgp in spot[1]:
            curr_intest = defineElementsOfInterest(rgp.genes, elements)

            fout.write(f"spot_{n_spot}\t{rgp.name}\t")
            fout.write("-\n" if len(curr_intest) == 0 else ','.join(curr_intest)+"\n")
        n_spot+=1
    fout.close()

def writeBorders_spots(spots, pangenome, output):
    fout = open(output + "/spots_families.faa","w")
    n_spot = 0
    
    for spot in spots:
        n_border_group = 0
        for borders in spot[0]:
            n_border = 0
            prevalence = len(spot[1])
            organisations = len(getUniqRGP(spot[1]))
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

def select_spots(pangenome, spots, elements, min_presence_ratio=0.05, min_org_ratio=0.01):
    to_draw= []
    z=False
    for spot in spots:
        nb_uniq = len(getUniqRGP(spot[1]))
        if nb_uniq > 1:
            for rgp in spot[1]:
                for gene in rgp.genes:
                    if gene.name in elements or any(x in gene.product for x in elements):
                        z=True
        if len(spot[1]) > len(pangenome.organisms)*min_presence_ratio and nb_uniq > max(len(pangenome.organisms) * min_org_ratio, 2):
            z=True
        if z:
            to_draw.append(spot)
            z=False
    return to_draw

def checkParameterLogic(overlapping_match, set_size, exact_match):
    if overlapping_match >= set_size:
        raise Exception(f'--overlapping_match_hotspot ({overlapping_match}) cannot be bigger than (or equal to) --set_size_hotspot ({set_size})')
    if exact_match > set_size:
        raise Exception(f'--exact_match_size_hotspot ({exact_match}) cannot be bigger than --set_size_hotspot ({set_size})')

def makeColorsForFams(fams):
    # potentialColors = cl.to_numeric([ col for val in cl.scales['8'].values() for val2 in val.values() for col in val2 ])
    # random.shuffle(potentialColors)
    famcol = {}
    # if len(fams) < len(potentialColors):#can't display if not (families would have the same color)
    for fam in fams:
        col =  list(random.choices(range(256), k=3))
        famcol[fam] = '#%02x%02x%02x' % (col[0], col[1], col[2])
    return famcol

def countRGPoccurrence(uniqRGPS, currRGP):
    countRGPs = Counter()
    for rgp in currRGP:
        for uniqrgp in uniqRGPS:
            if rgp == uniqrgp:
                countRGPs[uniqrgp] += 1
                break
    return countRGPs

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
        if 'RNA' in gene.type:
            present_EOI.add(gene.type)
        if 'integrase' in gene.product.lower():
            present_EOI.add('integrase')
        if gene.name in elements or any(x in gene.product for x in elements):
            present_EOI.add(gene.name)
    return present_EOI

def drawCurrSpot(genelists, ordered_counts, elements, famCol, filename):
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
                gene_names.append(' ' + gene.name)
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
    filename = filename +('_' + "_".join(present_EOI) if len(present_EOI) > 0 else "")
    Rannot = robjects.ListVector(annotList)
    Rdna_segs = robjects.ListVector(rdframes)
    plot_gene_map = robjects.r["plot_gene_map"]
    grdevices = importr('grDevices')
    grdevices.png(file=filename +".png", width = longest_gene_list * 70, height= len(rdframes) * 60)
    plot_gene_map(dna_segs = Rdna_segs, annotations=Rannot, lwd = 4)
    grdevices.dev_off()

def draw_spots(spots, output, overlapping_match, exact_match, set_size, multigenics, elements):
    logging.getLogger().info("Drawing the hotspots of the pangenome")
    bar = tqdm(range(len(spots)), unit = "region")
    for i, spot in enumerate(spots):
        uniqRGPS = frozenset(getUniqRGP(spot[1]))
        Fams = set()
        GeneLists = []

        countUniq = countRGPoccurrence(uniqRGPS, spot[1])

        #order unique rgps by occurrences
        sortedUniqRGPs = sorted(uniqRGPS, key = lambda x : countUniq[x], reverse=True)
        for rgp in sortedUniqRGPs:
            GeneList = []
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
        fname = output + '/hotspot_' + str(i)
        drawCurrSpot(GeneLists, ordered_counts, elements, famcol, fname)#make R dataframes, and plot them using genoPlotR.
        bar.update()
    bar.close()

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    mkOutdir(args.output, args.force)
    predictHotspots(pangenome, args.output, spot_graph=args.spot_graph, flanking_graph=args.flanking_graph, overlapping_match=args.overlapping_match, set_size=args.set_size, exact_match=args.exact_match_size, draw_hotspot=args.draw_hotspots, interest=args.interest)

def hotspotSubparser(subparser):
    parser = subparser.add_parser("hotspot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--spot_graph", required = False, action="store_true", help = "Writes a graph in a .gexf format of pairs of blocks of single copy markers flanking RGPs, supposedly belonging to the same hotspot")
    optional.add_argument("--flanking_graph", required = False, action="store_true", help = "Writes a graph in a .gexf format of common blocks of single copy markers flanking RGPs, supposedly with some common origin")
    optional.add_argument("--draw_hotspots", required=False, action="store_true", help = "Draws a figure representing all of the hotspots syntenies")
    optional.add_argument("--overlapping_match", required=False, type=int, default = 2, help="The number of 'missing' persistent genes allowed when comparing flanking genes during hotspot computations")
    optional.add_argument("--set_size", required = False, type = int, default = 3, help = "Number of single copy markers to use as flanking genes for a RGP during hotspot computation")
    optional.add_argument("--exact_match_size", required = False, type= int, default = 1, help = "Number of perfecty matching flanking single copy markers required to associate RGPs during hotspot computation (Ex: If set to 1, two RGPs are in the same hotspot if both their 1st flanking genes are the same)")
    optional.add_argument("--interest",required=False, type=str, default="",help = "Comma separated list of elements to flag when drawing hotspots")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    return parser