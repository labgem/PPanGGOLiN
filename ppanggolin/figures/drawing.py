#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import time
import os
from collections import defaultdict
#installed libraries
from tqdm import tqdm
import plotly.graph_objs as go
import plotly.offline as out_plotly
import numpy
from scipy.spatial.distance import squareform, pdist
from scipy.sparse import csr_matrix, bsr_matrix, csc_matrix
from scipy.stats import iqr, linregress, pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
import colorlover as cl

#local libraries
from ppanggolin.utils import mkOutdir
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import readPangenome, getNumberOfOrganisms

def checkPangenomeTilePlot(pangenome):
    if pangenome.status["partitionned"] == "No":
        logging.getLogger().warning("Your pangenome has not been partitionned ! So no tile_plot can be drawn.")
        return False
    if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"]:
       pass
    elif pangenome.status["genomesAnnotated"] == "inFile":
        readPangenome(pangenome, annotation = True)
    else:
        raise Exception("You want to partition an unannotated pangenome")
    if pangenome.status["genesClustered"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["genesClustered"] == "inFile":
        readPangenome(pangenome, geneFamilies= True)
    else:
        raise Exception("You want to partition a pangenome whose genes have not been clustered")
    if pangenome.status["neighborsGraph"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["neighborsGraph"] == "inFile":
        readPangenome(pangenome, graph=True)#whether it is faster to compute it or to load it will have to be checked on bigger graphs.
    else:
        raise Exception("You want to partition a pangenome whose neighbors graph has not been computed.")
    return True

def checkPangenomeUCurve(pangenome):
    if pangenome.status["genomesAnnotated"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["genomesAnnotated"] == "inFile":
        readPangenome(pangenome, annotation = True)
    else:
        raise Exception("You want to partition an unannotated pangenome")
    if pangenome.status["genesClustered"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["genesClustered"] == "inFile":
        readPangenome(pangenome, geneFamilies= True)
    else:
        raise Exception("You want to partition a pangenome whose genes have not been clustered")
    if pangenome.status["neighborsGraph"] in ["Computed","Loaded"]:
        pass
    elif pangenome.status["neighborsGraph"] == "inFile":
        readPangenome(pangenome, graph=True)#whether it is faster to compute it or to load it will have to be checked on bigger graphs.
    else:
        raise Exception("You want to partition a pangenome whose neighbors graph has not been computed.")

def jaccard_similarities(mat,jaccard_similarity_th):
    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat
    # for rows
    aa = numpy.repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]
    similarities = ab.copy()
    similarities.data /= (aa + bb - ab.data)
    similarities.data[similarities.data<jaccard_similarity_th] = 0
    similarities.eliminate_zeros()
    return similarities

def drawTilePlot(pangenome, output, nocloud = False):
    can_draw = checkPangenomeTilePlot(pangenome)
    if not can_draw:
        return
    if len(pangenome.organisms) > 500 and nocloud is False:
        logging.getLogger().warning("You asked to draw a tile plot for a lot of organisms (>500). Your browser will probably not be able to open it.")
    logging.getLogger().info("Drawing the tile plot...")
    data        = []
    all_indexes = []
    all_columns = []
    fam2index = {}
    index2fam = {}
    if nocloud:
        families = { fam for fam in pangenome.geneFamilies if not fam.partition.startswith("C")}
    else:
        families =  set(pangenome.geneFamilies)
    org_index = pangenome.getIndex()
    index2org = {}
    for org, index  in org_index.items():
        index2org[index] = org
    COLORS = {"pangenome":"black", "exact_accessory":"#EB37ED", "exact_core" :"#FF2828", "soft_core":"#c7c938", "soft_accessory":"#996633","shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF", "undefined":"#828282"}

    logging.getLogger().info("start with matrice")

    for row, fam in enumerate(families):
        new_col=[ org_index[org] for org in fam.organisms ]
        all_indexes.extend([row]*len(new_col))
        all_columns.extend(new_col)
        data.extend([1.0]*len(new_col))
        index2fam[row]=fam.name
        fam2index[fam.name] = row

    mat_p_a = csc_matrix((data, (all_indexes,all_columns)), shape = (len(families),len(pangenome.organisms)), dtype='float')
    dist    = pdist(1 - jaccard_similarities(mat_p_a,0).todense())
    hc      = linkage(dist, 'single')

    dendro = dendrogram(hc,no_plot=True)
    logging.getLogger().info("done with making the dendrogram to order the organisms on the plot")

    order_organisms = [ index2org[index] for index in dendro["leaves"]]

    binary_data = []
    text_data   = []
    fam_order   = []
    partitions_dict = defaultdict(list)
    shell_subs = set()
    for fam in families:
        partitions_dict[fam.partition].append(fam)
        if fam.partition.startswith("S"):
            shell_subs.add(fam.partition)#number of elements will tell the number of subpartitions
    ordered_nodes = []
    ordored_nodes_p = sorted(partitions_dict["P"], key=lambda n:len(n.organisms), reverse=True)
    ordored_nodes_c = sorted(partitions_dict["C"], key=lambda n:len(n.organisms), reverse=True)
    sep_p = len(ordored_nodes_p)-0.5
    separators = [sep_p]
    shell_NA = None
    if len(shell_subs)==1:
        ordored_nodes_s = sorted(partitions_dict[shell_subs.pop()], key=lambda n:len(n.organisms), reverse=True)
        ordered_nodes = ordored_nodes_p+ordored_nodes_s+ordored_nodes_c
        separators.append(separators[len(separators)-1]+len(ordored_nodes_s))
        separators.append(separators[len(separators)-1]+len(ordored_nodes_c))
    else:
        ordered_nodes = ordored_nodes_p
        for subpartition in sorted(shell_subs):
            if subpartition=="S_":
                shell_NA=len(separators)-1
            ordored_nodes_s = sorted(partitions_dict[subpartition], key=lambda n:len(n.organisms), reverse=True)
            ordered_nodes+= ordored_nodes_s
            separators.append(separators[len(separators)-1]+len(ordored_nodes_s))
        ordered_nodes+=ordored_nodes_c
        separators.append(separators[len(separators)-1]+len(ordored_nodes_c))
    
    logging.getLogger().info("Getting the gene name(s) and the number for each tile of the plot ...")
    for node in ordered_nodes:
        fam_order.append('\u200c' + node.name)
        data = node.organisms
        binary_data.append([len(node.getGenesPerOrg(org)) if org in data else numpy.nan for org in order_organisms])
        text_data.append([("\n".join(map(str,node.getGenesPerOrg(org)))) if org in data else numpy.nan for org in order_organisms])
    
    xaxis_values = [ '\u200c'+org.name for org in order_organisms ]

    logging.getLogger().info("Done extracting names and numbers. Making the heatmap ...")

    heatmap = go.Heatmap(z              = binary_data,
                            x              = xaxis_values,
                            y              = fam_order,
                            text           = text_data,
                            zauto          = False,
                            zmin           = 1,
                            zmax           = 2,
                            autocolorscale = False,
                            colorscale     = [[0.50, 'rgb(100, 15, 78)'],[1, 'rgb(59, 157, 50)']],
                            colorbar       = dict(title     = 'Presence/Absence',
                                                titleside = 'top',
                                                tickmode  = 'array',
                                                tickvals  = [1,2],
                                                ticktext  = ['Presence','Multicopy'],
                                                ticks     = 'outside'))
    shell_color=None
    if len(shell_subs)>1:
        if "S_" not in shell_subs:
            shell_color = cl.interp(cl.flipper()['seq']['9']['Greens'][1:7],len(shell_subs))
        else:
            shell_color = cl.interp(cl.flipper()['seq']['9']['Greens'][1:7],len(shell_subs)-1)
    shapes = []
    sep_prec=0
    for nb, sep in enumerate(separators):
        color = None
        if nb==0:
            color=COLORS["persistent"]
        elif nb==(len(separators)-1):
            color=COLORS["cloud"]
        elif len(shell_subs)>1:
            if shell_NA is not None and nb==shell_NA:
                color=COLORS["shell"]
            else:
                color=shell_color.pop()
        else:
            color=COLORS["shell"]
        shapes.append(dict(type='line', x0=-1, x1=-1, y0=sep_prec, y1=sep, line = dict(dict(width=10, color=color))))
        shapes.append(dict(type='line', x0=len(pangenome.organisms), x1=len(pangenome.organisms), y0=sep_prec, y1=sep, line = dict(dict(width=10, color=color))))
        shapes.append(dict(type='line', x0=-1, x1=len(pangenome.organisms), y0=sep, y1=sep, line = dict(dict(width=1, color=color))))
        sep_prec=sep

    layout = go.Layout(title  = "presence/absence matrix",
                        xaxis  = go.layout.XAxis(ticktext=xaxis_values,
                                                title='organisms',
                                                tickvals=xaxis_values,
                                                automargin=True,
                                                tickfont=dict(size=10)),
                        yaxis  = go.layout.YAxis(ticktext=fam_order,
                                                tickvals=fam_order,
                                                title='gene families',
                                                automargin=True,
                                                tickfont=dict(size=10)),
                        shapes = shapes)
    logging.getLogger().info("Drawing the figure itself...")
    out_plotly.plot(go.Figure(data=[heatmap], layout=layout), filename = output+"/tile_plot.html", auto_open=False)
    logging.getLogger().info(f"Done with the tile plot : '{output+'/tile_plot.html'}' ")

def drawUCurve(pangenome, output, soft_core = 0.95):
    checkPangenomeUCurve(pangenome)
    logging.getLogger().info("Drawing the U-shaped curve...")
    max_bar = 0
    count = defaultdict(lambda : defaultdict(int))
    is_partitionned = False
    has_undefined = False
    for fam in pangenome.geneFamilies:
        nb_org  = len(fam.organisms)
        if fam.partition != "":
            is_partitionned = True
            if fam.partition == "U":
                has_undefined = True
            count[nb_org][fam.namedPartition]+=1
        count[nb_org]["pangenome"]+=1
        max_bar = count[nb_org]["pangenome"] if count[nb_org]["pangenome"] > max_bar else max_bar
    data_plot = []
    chao = "NA"
    if count[1]["pangenome"] > 0:
        chao = round(len(pangenome.geneFamilies) + ((count[0]["pangenome"]^2)/(count[1]["pangenome"]*2)),2)
    COLORS = {"pangenome":"black", "exact_accessory":"#EB37ED", "exact_core" :"#FF2828", "soft_core":"#c7c938", "soft_accessory":"#996633","shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF", "undefined":"#828282"}
    if is_partitionned and not has_undefined:
        persistent_values = []
        shell_values      = []
        cloud_values      = []
        for nb_org in range(1,len(pangenome.organisms)+1):
            persistent_values.append(count[nb_org]["persistent"])
            shell_values.append(count[nb_org]["shell"])
            cloud_values.append(count[nb_org]["cloud"])
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=persistent_values,name='persistent', marker=dict(color = COLORS["persistent"])))
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=shell_values,name='shell', marker=dict(color = COLORS["shell"])))
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=cloud_values,name='cloud', marker=dict(color = COLORS["cloud"])))
    else:
        text = 'undefined' if has_undefined else "pangenome"
        undefined_values = []
        for nb_org in range(1,len(pangenome.organisms)+1):
            undefined_values.append(count[nb_org][text])
        data_plot.append(go.Bar(x=list(range(1,len(pangenome.organisms)+1)),y=undefined_values,name=text, marker=dict(color = COLORS[text])))
    layout = None
    x = len(pangenome.organisms)*soft_core
    layout =  go.Layout(title = "Gene families frequency distribution (U shape), chao="+str(chao),
                        xaxis = dict(title='Occurring in x genomes'),
                        yaxis = dict(title='# of gene families (F)'),
                        barmode='stack', shapes=[dict(type='line', x0=x, x1=x, y0=0, y1=max_bar, line = dict(dict(width=5, dash='dashdot', color="grey")))])

    fig = go.Figure(data=data_plot, layout=layout)
    out_plotly.plot(fig, filename = output+"/Ushaped_plot.html", auto_open=False)
    logging.getLogger().info(f"Done drawing the U-shaped curve : '{output+'/Ushaped_plot.html'}'")

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.tile_plot:
        drawTilePlot(pangenome, args.output, args.nocloud)
    if args.ucurve:
        drawUCurve(pangenome, args.output, soft_core = args.soft_core)

def figureSubparser(subparser):
    parser = subparser.add_parser("draw",help = "Draw figures representing the pangenome through different aspects")
    optional = parser.add_argument_group(title = "Optional arguments")

    optional.add_argument('-o','--output', required=False, type=str, default="ppanggolin_output"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory")
    optional.add_argument("--tile_plot",required = False, default = False, action = "store_true",help = "draw the tile plot")
    optional.add_argument("--nocloud", required=False, default = False, action = "store_true", help = "Do not draw the cloud in the tile plot")
    optional.add_argument("--soft_core",required=False, default = 0.95, help = "Soft core threshold to use")
    optional.add_argument("--ucurve",required = False, default = False, action = "store_true",help = "draw the U-curve")

    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    
    return parser