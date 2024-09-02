#!/usr/bin/env python3

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
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import numpy as np


# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import jaccard_similarities


def draw_tile_plot(pangenome: Pangenome,
                   output: Path,
                   nocloud: bool = False,
                   draw_dendrogram:bool=False,
                   disable_bar: bool = False,):
    """
    Draw a tile plot from a partitioned pangenome.

    :param pangenome: Partitioned pangenome
    :param output: Path to output directory
    :param nocloud: Do not draw the cloud partition
    :param disable_bar: Allow to disable progress bar
    """

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True, disable_bar=disable_bar)
    if pangenome.status["partitioned"] == "No":
        raise Exception("Cannot draw the tile plot as your pangenome has not been partitioned")

    if pangenome.number_of_organisms > 500 and not nocloud:
        logging.getLogger("PPanGGOLiN").warning(
            "You asked to draw a tile plot for a lot of genomes (>500). Your browser will probably not be able to open it."
        )


    logging.getLogger("PPanGGOLiN").info("Drawing the tile plot...")

    families, org_index, index2org = prepare_data_structures(pangenome, nocloud)
    
    mat_p_a, index2fam = build_presence_absence_matrix(families, org_index)
    order_organisms, dendrogram_fig = generate_dendrogram(mat_p_a, org_index)

    binary_data, text_data, fam_order, separators = process_tile_data(families, order_organisms)

    fig, dendro = create_tile_plot(binary_data, text_data, fam_order, separators, order_organisms, dendrogram_fig, draw_dendrogram)
    #save_plot(fig, output)
    filename = output / "tile_plot.html"
    fig.write_html(filename)

    logging.getLogger("PPanGGOLiN").info(f"Done with the tile plot: '{filename}'")

    
    return fig, dendro



def warn_large_genomes(pangenome: Pangenome, nocloud: bool):
    """Warn if the number of genomes is large, potentially causing performance issues."""


def prepare_data_structures(pangenome: Pangenome, nocloud: bool):
    """Prepare data structures for plotting."""
    if nocloud:
        families = {fam for fam in pangenome.gene_families if not fam.partition.startswith("C")}
    else:
        families = set(pangenome.gene_families)
    
    org_index = pangenome.get_org_index()
    index2org = {index: org for org, index in org_index.items()}
    
    return families, org_index, index2org


def build_presence_absence_matrix(families: set, org_index: dict):
    """Build the presence-absence matrix for gene families."""
    data, all_indexes, all_columns = [], [], []
    fam2index, index2fam = {}, {}
    
    for row, fam in enumerate(families):
        new_col = [org_index[org] for org in fam.organisms]
        all_indexes.extend([row] * len(new_col))
        all_columns.extend(new_col)
        data.extend([1.0] * len(new_col))
        index2fam[row] = fam.name
        fam2index[fam.name] = row
    
    mat_p_a = csc_matrix((data, (all_indexes, all_columns)), shape=(len(families), len(org_index)), dtype='float')
    return mat_p_a, index2fam


def generate_dendrogram(mat_p_a, org_index):
    """Generate the order of organisms based on a dendrogram."""



    # dist = pdist(1 - jaccard_similarities(mat_p_a, 0).todense())
    # hc = linkage(dist, 'single')
    # dendro_org = dendrogram(hc, no_plot=True)
    # return [index2org[index] for index in dendro_org["leaves"]]

    genom_names = [org.name for org in org_index]

    name_to_org =  {org.name:org for org in org_index}

    distance_matrice = 1 - jaccard_similarities(mat_p_a, 0).todense()

    dendrogram_fig = ff.create_dendrogram(distance_matrice, labels=genom_names, orientation='bottom')


    for i in range(len(dendrogram_fig['data'])):
        dendrogram_fig['data'][i]['yaxis'] = 'y2'
        # dendrogram_fig['data'][i]['showlegend'] = False



    order_organisms = [name_to_org[org_name] for org_name in dendrogram_fig['layout']['xaxis']['ticktext']]

    # dendrogram_fig.update_layout(width=800, height=800)
    return order_organisms, dendrogram_fig


def process_tile_data(families, order_organisms):
    """Process data for each tile in the plot."""
    binary_data, text_data, fam_order = [], [], []
    partitions_dict = defaultdict(list)
    shell_subs = set()
    
    for fam in families:
        partitions_dict[fam.partition].append(fam)
        if fam.partition.startswith("S"):
            shell_subs.add(fam.partition)
    
    ordered_nodes, separators = order_nodes(partitions_dict, shell_subs)
    
    for node in ordered_nodes:
        fam_order.append(node.name)
        data = set(node.organisms)
        binary_data.append([len(list(node.get_genes_per_org(org))) if org in data else np.nan for org in order_organisms])
        text_data.append([("\n".join(map(str, node.get_genes_per_org(org))) if org in data else np.nan) for org in order_organisms])

    text_data = get_heatmap_hover_text(ordered_nodes, order_organisms)
    
    return binary_data, text_data, fam_order, separators


def order_nodes(partitions_dict, shell_subs):
    """Order gene families based on their partitions."""
    ordered_nodes_p = sorted(partitions_dict["P"], key=lambda n: n.number_of_organisms, reverse=True)
    ordered_nodes_c = sorted(partitions_dict["C"], key=lambda n: n.number_of_organisms, reverse=True)
    separators = [len(ordered_nodes_p) - 0.5]
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
    
    return ordered_nodes, separators


def create_partition_shapes(separators, xval_max, heatmap_row, colors):
    """Create the shapes for plot separators."""


    shapes = []
    sep_prec = 0
    xref='x1'
    yref=f'y{heatmap_row}'

    for nb, sep in enumerate(separators):

        if nb == 0:
            color = colors["persistent"]
        elif nb == (len(separators) - 1):
            color = colors["cloud"]
        else:
            color = colors["shell"]
        shapes.append(dict(type='line', x0=-1, x1=-1, y0=sep_prec, y1=sep, line=dict(width=10, color=color), xref= xref, yref=yref, name="Persistent", showlegend=True))
        shapes.append(dict(type='line', x0=xval_max, x1=xval_max, y0=sep_prec, y1=sep, line=dict(width=10, color=color), xref= xref, yref=yref, name="Persistent", showlegend=True))
        shapes.append(dict(type='line', x0=-1, x1=xval_max, y0=sep, y1=sep, line=dict(width=1, color=color), xref= xref, yref=yref, name="Persistent", showlegend=True))
        sep_prec = sep

    return shapes


def metadata_stringify(gene):
    metadata_str = ''
    if gene.has_metadata():
        
        metadata_str = f'<br><br>{gene.ID} metadata'
        for metadata in gene.metadata:
            metadata_str += f"<br>metadata source: {metadata.source}"
            metadata_dict = metadata.to_dict()
            metadata_str += '<br>'.join((f"{key}: {value}" for key,value in metadata_dict.items()))
            
    return metadata_str
    
def get_heatmap_hover_text(ordered_families, order_organisms):


    text_data = []
    for family in ordered_families:
        text_per_family = []
        for org in order_organisms:
            if org in family.organisms:
                gene_count = len(list(family.get_genes_per_org(org)))
                genes = "<br>- ".join(map(str, family.get_genes_per_org(org)))
    
                names = ";".join( (gene.name for gene in family.get_genes_per_org(org) if gene.name))
                if names:
                    print(names)
                extra_gene_info = f"genes:<br>- {genes}"
                if names:
                    extra_gene_info += f'<br>names:{names}'
                    
                metadata = "<br>".join( (metadata_stringify(gene) for gene in family.get_genes_per_org(org) if gene.has_metadata())) 
                extra_gene_info += metadata
                
            else:
                gene_count = 0
                extra_gene_info = numpy.nan
    
            gene_info = f"genome:{org.name}<br>family:{family.name}<br>gene_count:{gene_count}<br>{extra_gene_info}" 
            text_per_family.append(gene_info)
                               
        text_data.append(text_per_family)
    return text_data
    



def create_tile_plot(binary_data, text_data, fam_order, separators, order_organisms, dendrogram_fig, draw_dendrogram):
    """Create the heatmap tile plot using Plotly."""
  
    xaxis_values = [org.name for org in order_organisms]
    colors =  {
        "pangenome": "black", "exact_accessory": "#EB37ED", "exact_core": "#FF2828",
        "soft_core": "#c7c938", "soft_accessory": "#996633", "shell": "#00D860",
        "persistent": "#F7A507", "cloud": "#79DEFF", "undefined": "#828282",
        "presence":"#005AB5", # blue
        "multicopy":'#DC3220' # red
    }


    heatmap = [go.Heatmap(z=binary_data,
                 x=xaxis_values,
                 y=fam_order,
                 text=text_data,
                 zauto=False,
                 zmin=0,
                 zmax=2,
                 autocolorscale=False,
                 hoverinfo="text",
                 colorscale=[[0, '#ffffff'],[0.33, '#ffffff'],
                             [0.33, colors['presence']],[0.66, colors['presence']],
                             [0.66, colors['multicopy']], [1, colors['multicopy']]],
                 colorbar=dict(title='Presence/Absence',
                               titleside='top',
                               tickmode='array',
                               tickvals=[0.33, 1, 1.66],
                               ticktext=['Absence', 'Presence', 'Multicopy'],
                               len=0.27,
                               outlinecolor='black',
                               outlinewidth=0.5,
                               ticks=None, orientation='v'))]



    if draw_dendrogram:
        heatmap_row = 2
        fig = make_subplots(rows=2, cols=1,
                    shared_xaxes=True,
                    vertical_spacing=0.01,
                row_heights=[0.2, 0.8])
    
        for data in dendrogram_fig['data']:
            fig.add_trace(data,  row=1, col=1)

    else:
        heatmap_row = 1
        fig = make_subplots(rows=1, cols=1)


    heatmap[0]['x'] = dendrogram_fig['layout']['xaxis']['tickvals']

    for data in heatmap:
        
        fig.add_trace(data, row=heatmap_row, col=1)

    layout = go.Layout(title="Presence-Absence Matrix",
                plot_bgcolor='#ffffff')

    if draw_dendrogram:
        layout.xaxis2 = dendrogram_fig.layout.xaxis
    else:
        layout.xaxis = dendrogram_fig.layout.xaxis

    fig.update_layout(layout)



    fig.update_xaxes(
        ticklen=0,
        title="Genomes"
    )
    fig.update_yaxes(
        ticklen=0,
        title='Gene Families',
        tickfont=dict(size=10),
        automargin=True,
    )
    fig.layout.yaxis.title = None

    xmax = dendrogram_fig['layout']['xaxis']['tickvals'][-1] +  dendrogram_fig['layout']['xaxis']['tickvals'][0] + 0.5
    shapes = create_partition_shapes(separators, xmax, heatmap_row, colors)
    
    fig.update_layout(go.Layout(shapes=shapes, 
                                showlegend=False))

    fig.write_html("file.html")

    fig.update_layout({'width':1000, 'height':1000,

                        })
    return fig, dendrogram_fig

    