#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
from collections import defaultdict

#installed libraries
import numpy
from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram
import plotly.graph_objs as go
import plotly.offline as out_plotly
import colorlover as cl
#local libraries
from ppanggolin.formats import checkPangenomeInfo

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
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, needGraph=True)
    if pangenome.status["partitionned"] == "No":
        raise Exception("Cannot draw the tile plot as your pangenome has not been partitionned")
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
                        shapes = shapes,
                        plot_bgcolor='#ffffff')
    logging.getLogger().info("Drawing the figure itself...")
    out_plotly.plot(go.Figure(data=[heatmap], layout=layout), filename = output+"/tile_plot.html", auto_open=False)
    logging.getLogger().info(f"Done with the tile plot : '{output+'/tile_plot.html'}' ")
