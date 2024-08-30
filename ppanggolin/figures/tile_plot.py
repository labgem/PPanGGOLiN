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

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import jaccard_similarities


def draw_tile_plot(pangenome: Pangenome, output: Path=None, nocloud: bool = False, disable_bar: bool = False):
    """
    Draw a tile plot from a partitioned pangenome.

    :param pangenome: Partitioned pangenome
    :param output: Path to output directory
    :param nocloud: Do not draw the cloud partition
    :param disable_bar: Allow to disable progress bar
    """
    validate_pangenome(pangenome, disable_bar)
    warn_large_genomes(pangenome, nocloud)
    logging.getLogger("PPanGGOLiN").info("Drawing the tile plot...")

    families, org_index, index2org = prepare_data_structures(pangenome, nocloud)
    mat_p_a, index2fam = build_presence_absence_matrix(families, org_index)

    order_organisms = generate_dendrogram_order(mat_p_a, index2org)
    binary_data, text_data, fam_order, separators = process_tile_data(families, order_organisms)

    fig = create_tile_plot(binary_data, text_data, fam_order, separators, pangenome, order_organisms)
    #save_plot(fig, output)
    #logging.getLogger("PPanGGOLiN").info(f"Done with the tile plot: '{output / 'tile_plot.html'}'")
    return fig

def validate_pangenome(pangenome: Pangenome, disable_bar: bool):
    """Check if the pangenome is ready for plotting."""
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True, disable_bar=disable_bar)
    if pangenome.status["partitioned"] == "No":
        raise Exception("Cannot draw the tile plot as your pangenome has not been partitioned")


def warn_large_genomes(pangenome: Pangenome, nocloud: bool):
    """Warn if the number of genomes is large, potentially causing performance issues."""
    if pangenome.number_of_organisms > 500 and not nocloud:
        logging.getLogger("PPanGGOLiN").warning(
            "You asked to draw a tile plot for a lot of genomes (>500). Your browser will probably not be able to open it."
        )


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


def generate_dendrogram_order(mat_p_a, index2org):
    """Generate the order of organisms based on a dendrogram."""
    dist = pdist(1 - jaccard_similarities(mat_p_a, 0).todense())
    hc = linkage(dist, 'single')
    dendro_org = dendrogram(hc, no_plot=True)
    return [index2org[index] for index in dendro_org["leaves"]]


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


def get_color_scheme():
    """Get the color scheme for the plot."""
    return {
        "pangenome": "black", "exact_accessory": "#EB37ED", "exact_core": "#FF2828",
        "soft_core": "#c7c938", "soft_accessory": "#996633", "shell": "#00D860",
        "persistent": "#F7A507", "cloud": "#79DEFF", "undefined": "#828282"
    }


def create_shapes(separators, pangenome, colors):
    """Create the shapes for plot separators."""
    shapes = []
    sep_prec = 0
    
    for nb, sep in enumerate(separators):
        print(sep)
        if nb == 0:
            color = colors["persistent"]
        elif nb == (len(separators) - 1):
            color = colors["cloud"]
        else:
            color = colors["shell"]
        shapes.append(dict(type='line', x0=-1, x1=-1, y0=sep_prec, y1=sep, line=dict(width=10, color=color)))
        shapes.append(dict(type='line', x0=pangenome.number_of_organisms, x1=pangenome.number_of_organisms, y0=sep_prec, y1=sep, line=dict(width=10, color=color)))
        shapes.append(dict(type='line', x0=-1, x1=pangenome.number_of_organisms, y0=sep, y1=sep, line=dict(width=1, color=color)))
        sep_prec = sep
    
    return shapes


def save_plot(fig, output: Path):
    """Save the plot to an HTML file."""
    out_plotly(fig, filename=output / "tile_plot.html", auto_open=False)


def metadata_stringify(gene):
    metadata_str = ''
    if gene.has_metadata():
        
        metadata_str = f'<br><br>{gene.ID} metadata'
        for metadata in gene.metadata:
            metadata_str += f"<br>metadata source: {metadata.source}"
            metadata_dict = metadata.to_dict()
            metadata_str += '<br>'.join((f"{key}: {value}" for key,value in metadata_dict.items()))
            print('=='*10)
            print(metadata_str)
            
    return metadata_str
    
def get_heatmap_hover_text(ordered_families, order_organisms):


    text_data = []
    f = False
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
    