#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from collections import defaultdict
from pathlib import Path

# installed libraries
import numpy
from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram
import plotly.graph_objs as go
import plotly.offline as out_plotly
import colorlover as cl

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import jaccard_similarities


def draw_tile_plot(pangenome: Pangenome, output: Path, nocloud: bool = False, disable_bar: bool = False):
    """
    Draw a tile plot from a partitioned pangenome

    :param pangenome: Partitioned pangenome
    :param output: Path to output directory
    :param nocloud: Do not draw the cloud partition
    :param disable_bar: Allow to disable progress bar
    """

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True, disable_bar=disable_bar)
    if pangenome.status["partitioned"] == "No":
        raise Exception("Cannot draw the tile plot as your pangenome has not been partitioned")
    if pangenome.number_of_organisms > 500 and nocloud is False:
        logging.getLogger("PPanGGOLiN").warning("You asked to draw a tile plot for a lot of genomes (>500). "
                                                "Your browser will probably not be able to open it.")
    logging.getLogger("PPanGGOLiN").info("Drawing the tile plot...")
    data = []
    all_indexes = []
    all_columns = []
    fam2index = {}
    index2fam = {}
    if nocloud:
        families = {fam for fam in pangenome.gene_families if not fam.partition.startswith("C")}
    else:
        families = set(pangenome.gene_families)
    org_index = pangenome.get_org_index()
    index2org = {}
    for org, index in org_index.items():
        index2org[index] = org
    colors = {"pangenome": "black", "exact_accessory": "#EB37ED", "exact_core": "#FF2828", "soft_core": "#c7c938",
              "soft_accessory": "#996633", "shell": "#00D860", "persistent": "#F7A507", "cloud": "#79DEFF",
              "undefined": "#828282"}

    logging.getLogger("PPanGGOLiN").info("start with matrice")

    for row, fam in enumerate(families):
        new_col = [org_index[org] for org in fam.organisms]
        all_indexes.extend([row] * len(new_col))
        all_columns.extend(new_col)
        data.extend([1.0] * len(new_col))
        index2fam[row] = fam.name
        fam2index[fam.name] = row

    mat_p_a = csc_matrix((data, (all_indexes, all_columns)), shape=(len(families), pangenome.number_of_organisms),
                         dtype='float')
    dist = pdist(1 - jaccard_similarities(mat_p_a, 0).todense())
    hc = linkage(dist, 'single')

    dendro_org = dendrogram(hc, no_plot=True)
    logging.getLogger("PPanGGOLiN").info("done with making the dendrogram to order the genomes on the plot")

    order_organisms = [index2org[index] for index in dendro_org["leaves"]]

    binary_data = []
    text_data = []
    fam_order = []
    partitions_dict = defaultdict(list)
    shell_subs = set()
    for fam in families:
        partitions_dict[fam.partition].append(fam)
        if fam.partition.startswith("S"):
            shell_subs.add(fam.partition)  # number of elements will tell the number of subpartitions
    ordered_nodes_p = sorted(partitions_dict["P"], key=lambda n: n.number_of_organisms, reverse=True)
    ordered_nodes_c = sorted(partitions_dict["C"], key=lambda n: n.number_of_organisms, reverse=True)
    sep_p = len(ordered_nodes_p) - 0.5
    separators = [sep_p]
    shell_na = None
    if len(shell_subs) == 1:
        ordered_nodes_s = sorted(partitions_dict[shell_subs.pop()], key=lambda n: n.number_of_organisms, reverse=True)
        ordered_nodes = ordered_nodes_p + ordered_nodes_s + ordered_nodes_c
        separators.append(separators[len(separators) - 1] + len(ordered_nodes_s))
        separators.append(separators[len(separators) - 1] + len(ordered_nodes_c))
    else:
        ordered_nodes = ordered_nodes_p
        for subpartition in sorted(shell_subs):
            if subpartition == "S_":
                shell_na = len(separators) - 1
            ordered_nodes_s = sorted(partitions_dict[subpartition], key=lambda n: n.number_of_organisms, reverse=True)
            ordered_nodes += ordered_nodes_s
            separators.append(separators[len(separators) - 1] + len(ordered_nodes_s))
        ordered_nodes += ordered_nodes_c
        separators.append(separators[len(separators) - 1] + len(ordered_nodes_c))

    logging.getLogger("PPanGGOLiN").info("Getting the gene name(s) and the number for each tile of the plot ...")
    for node in ordered_nodes:
        fam_order.append(node.name)
        data = set(node.organisms)
        binary_data.append([len(list(node.get_genes_per_org(org))) if org in data else numpy.nan for org in order_organisms])
        text_data.append([("\n".join(map(str, node.get_genes_per_org(org))))
                          if org in data else numpy.nan for org in order_organisms])

    xaxis_values = [org.name for org in order_organisms]

    logging.getLogger("PPanGGOLiN").info("Done extracting names and numbers. Making the heatmap ...")

    heatmap = go.Heatmap(z=binary_data,
                         x=xaxis_values,
                         y=fam_order,
                         text=text_data,
                         zauto=False,
                         zmin=1,
                         zmax=2,
                         autocolorscale=False,
                         colorscale=[[0.50, 'rgb(100, 15, 78)'], [1, 'rgb(59, 157, 50)']],
                         colorbar=dict(title='Presence/Absence',
                                       titleside='top',
                                       tickmode='array',
                                       tickvals=[1, 2],
                                       ticktext=['Presence', 'Multicopy'],
                                       ticks='outside'))
    shell_color = None
    if len(shell_subs) > 1:
        if "S_" not in shell_subs:
            shell_color = cl.interp(cl.flipper()['seq']['9']['Greens'][1:7], len(shell_subs))
        else:
            shell_color = cl.interp(cl.flipper()['seq']['9']['Greens'][1:7], len(shell_subs) - 1)
    shapes = []
    sep_prec = 0
    for nb, sep in enumerate(separators):
        if nb == 0:
            color = colors["persistent"]
        elif nb == (len(separators) - 1):
            color = colors["cloud"]
        elif len(shell_subs) > 1:
            if shell_na is not None and nb == shell_na:
                color = colors["shell"]
            else:
                color = shell_color.pop()
        else:
            color = colors["shell"]
        shapes.append(dict(type='line', x0=-1, x1=-1, y0=sep_prec, y1=sep, line=dict(dict(width=10, color=color))))
        shapes.append(dict(type='line', x0=pangenome.number_of_organisms, x1=pangenome.number_of_organisms, y0=sep_prec, y1=sep,
                           line=dict(dict(width=10, color=color))))
        shapes.append(dict(type='line', x0=-1, x1=pangenome.number_of_organisms, y0=sep, y1=sep,
                           line=dict(dict(width=1, color=color))))
        sep_prec = sep

    layout = go.Layout(title="presence/absence matrix",
                       xaxis=go.layout.XAxis(ticktext=xaxis_values,
                                             title='genomes',
                                             tickvals=xaxis_values,
                                             automargin=True,
                                             tickfont=dict(size=10)),
                       yaxis=go.layout.YAxis(ticktext=fam_order,
                                             tickvals=fam_order,
                                             title='gene families',
                                             automargin=True,
                                             tickfont=dict(size=10)),
                       shapes=shapes,
                       plot_bgcolor='#ffffff')
    logging.getLogger("PPanGGOLiN").info("Drawing the figure itself...")

    fig = go.Figure(data=[heatmap])
    fig.update_layout(layout)
    out_plotly.plot(fig, filename=output.as_posix() + "/tile_plot.html", auto_open=False)
    logging.getLogger("PPanGGOLiN").info(f"Done with the tile plot : '{output / 'tile_plot.html'}' ")
