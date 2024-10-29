#!/usr/bin/env python3

# default libraries
import logging
from collections import defaultdict
from pathlib import Path
from itertools import cycle
from typing import List, Tuple, Dict, Set, Optional

# installed libraries
from scipy.sparse import csc_matrix
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import numpy as np

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.genome import Organism
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import jaccard_similarities


def draw_tile_plot(
    pangenome: Pangenome,
    output: Path,
    nocloud: bool = False,
    draw_dendrogram: bool = False,
    add_metadata: bool = False,
    metadata_sources: Optional[Set[str]] = None,
    disable_bar: bool = False,
):
    """
    Draw a tile plot from a partitioned pangenome.

    :param pangenome: Partitioned pangenome.
    :param output: Path to the output directory where the tile plot will be saved.
    :param nocloud: If True, exclude the cloud partition from the plot.
    :param draw_dendrogram: If True, include a dendrogram in the tile plot.
    :param disable_bar: If True, disable the progress bar during processing.
    """

    # Check if the pangenome has the required information and is partitioned
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=disable_bar,
        need_metadata=add_metadata,
        sources=metadata_sources,
    )
    if pangenome.status["partitioned"] == "No":
        raise Exception(
            "Cannot draw the tile plot as the pangenome has not been partitioned."
        )

    # Warn if there are more than 32767 genomes, as the output might be too large for browsers to handle
    if pangenome.number_of_organisms > 32767:
        logging.getLogger("PPanGGOLiN").warning(
            "You requested to draw a tile plot for a large number of genomes (>32k). "
            "This may result in a plot that is too large for web browsers to render."
        )
    if pangenome.number_of_gene_families > 32767 and not nocloud:
        logging.getLogger("PPanGGOLiN").warning(
            "You requested to draw a tile plot for a pangenome with a large number of families (>32k). "
            "This may result in a plot that is too large for web browsers to render."
            "You can use the --nocloud flag to exclude cloud families from the plot. "
        )

    logging.getLogger("PPanGGOLiN").info(
        "Starting the process of drawing the tile plot..."
    )

    # Prepare the data structures required for generating the tile plot
    families, org_index = prepare_data_structures(pangenome, nocloud)

    # Build the presence-absence matrix for the families and generate the dendrogram if required
    mat_p_a = build_presence_absence_matrix(families, org_index)
    order_organisms, dendrogram_fig = generate_dendrogram(mat_p_a, org_index)

    # Process the data to be displayed in the tile plot
    binary_data, text_data, fam_order, separators = process_tile_data(
        families, order_organisms
    )

    # Create the tile plot figure with or without the dendrogram
    fig = create_tile_plot(
        binary_data,
        text_data,
        fam_order,
        separators,
        order_organisms,
        dendrogram_fig,
        draw_dendrogram,
    )

    # Save the plot to the specified output directory
    filename = output / "tile_plot.html"
    fig.write_html(filename)

    logging.getLogger("PPanGGOLiN").info(
        f"Tile plot successfully created and saved to: '{filename}'."
    )

    return fig


def prepare_data_structures(pangenome: Pangenome, nocloud: bool) -> Tuple[set, dict]:
    """
    Prepare data structures required for generating the tile plot.

    :param pangenome: Partitioned pangenome containing gene families and organism data.
    :param nocloud: If True, exclude gene families belonging to the cloud partition.
    :return: A tuple containing a set of gene families to be plotted and a dictionary mapping organisms to their indices.
    """

    # Exclude gene families in the cloud partition if 'nocloud' is True; otherwise, include all gene families
    if nocloud:
        families = {
            fam for fam in pangenome.gene_families if not fam.partition.startswith("C")
        }
    else:
        families = set(pangenome.gene_families)

    # Get the organism index mapping from the pangenome
    org_index = pangenome.get_org_index()
    return families, org_index


def build_presence_absence_matrix(families: set, org_index: dict) -> csc_matrix:
    """
    Build the presence-absence matrix for gene families.

    This matrix indicates the presence (1) or absence (0) of each gene family across different organisms.

    :param families: A set of gene families to be included in the matrix.
    :param org_index: A dictionary mapping each organism to its respective index in the matrix.
    :return: A sparse matrix (Compressed Sparse Column format) representing the presence-absence of gene families.
    """

    # Initialize lists to store matrix data in a sparse format
    data, all_indexes, all_columns = [], [], []

    # Iterate through each gene family to populate the presence-absence matrix
    for row, fam in enumerate(families):
        # Find the indices of organisms that have the current gene family
        new_col = [org_index[org] for org in fam.organisms]
        all_indexes.extend([row] * len(new_col))  # Row index repeated for each presence
        all_columns.extend(new_col)  # Corresponding column indices for the organisms
        data.extend([1.0] * len(new_col))  # All presences are marked with 1.0

    # Construct the presence-absence matrix using Compressed Sparse Column format
    mat_p_a = csc_matrix(
        (data, (all_indexes, all_columns)),
        shape=(len(families), len(org_index)),
        dtype="float",
    )

    return mat_p_a


def generate_dendrogram(mat_p_a: csc_matrix, org_index: dict) -> Tuple[List, go.Figure]:
    """
    Generate the order of organisms based on a dendrogram.

    :param mat_p_a: Sparse matrix representing the presence-absence of gene families.
    :param org_index: Dictionary mapping organism names to their respective indices in the matrix.
    :return: A tuple containing the ordered list of organisms and the dendrogram figure.
    """

    # Extract organism names from the org_index dictionary
    genom_names = [org.name for org in org_index]

    # Create a mapping from organism names to organism objects
    name_to_org = {org.name: org for org in org_index}

    # Compute the distance matrix using Jaccard similarity
    distance_matrice = 1 - jaccard_similarities(mat_p_a, 0).todense()

    # Create a dendrogram figure using the computed distance matrix
    dendrogram_fig = ff.create_dendrogram(
        distance_matrice, labels=genom_names, orientation="bottom"
    )

    # Adjust the dendrogram figure to make it match with the heatmap later on
    for i in range(len(dendrogram_fig["data"])):
        dendrogram_fig["data"][i][
            "yaxis"
        ] = "y2"  # Aligns dendrogram data on a secondary y-axis
        dendrogram_fig["data"][i][
            "showlegend"
        ] = False  # Hides legends in the dendrogram

    # Extract the ordered list of organisms from the dendrogram tick labels
    order_organisms = [
        name_to_org[org_name]
        for org_name in dendrogram_fig["layout"]["xaxis"]["ticktext"]
    ]

    return order_organisms, dendrogram_fig


def process_tile_data(
    families: set, order_organisms: List
) -> Tuple[List[List[float]], List[List[str]], List[str], List[Tuple[str, float]]]:
    """
    Process data for each tile in the plot.

    :param families: A set of gene families to be processed.
    :param order_organisms: The ordered list of organisms for the tile plot.
    :return: A tuple containing binary data, text data, family order, and separators for the plot.
    """
    binary_data, text_data, fam_order = [], [], []
    partitions_dict = defaultdict(list)
    shell_subs = set()

    # Group families by partition and identify shell subpartitions
    for fam in families:
        partitions_dict[fam.partition].append(fam)
        if fam.partition.startswith("S"):
            shell_subs.add(fam.partition)

    ordered_nodes, separators = order_nodes(partitions_dict, shell_subs)

    # Populate binary and text data for each family
    for node in ordered_nodes:
        fam_order.append(node.name)
        data = set(node.organisms)
        binary_data.append(
            [
                len(list(node.get_genes_per_org(org))) if org in data else np.nan
                for org in order_organisms
            ]
        )
        text_data.append(
            [
                (
                    "\n".join(map(str, node.get_genes_per_org(org)))
                    if org in data
                    else np.nan
                )
                for org in order_organisms
            ]
        )

    # Generate hover text for the heatmap
    text_data = get_heatmap_hover_text(ordered_nodes, order_organisms)

    return binary_data, text_data, fam_order, separators


def order_nodes(
    partitions_dict: dict, shell_subs: set
) -> Tuple[List, List[Tuple[str, float]]]:
    """
    Order gene families based on their partitions.

    :param partitions_dict: A dictionary where keys are partition names and values are lists of gene families in each partition.
    :param shell_subs: A set of shell subpartition names.
    :return: A tuple containing the ordered list of gene families and a list of partition separators.
    """

    # Sort persistent and cloud partitions by the number of organisms in descending order
    ordered_nodes_p = sorted(
        partitions_dict["P"], key=lambda n: n.number_of_organisms, reverse=True
    )
    ordered_nodes_c = sorted(
        partitions_dict["C"], key=lambda n: n.number_of_organisms, reverse=True
    )

    partition_separators = [("Persistent", len(ordered_nodes_p) - 0.5)]
    ordered_nodes = ordered_nodes_p

    # Sort shell subpartitions and add them to the ordered nodes list
    for subpartition in sorted(shell_subs):
        partition_name = "Shell" if len(shell_subs) == 1 else f"Shell_{subpartition}"
        ordered_nodes_s = sorted(
            partitions_dict[subpartition],
            key=lambda n: n.number_of_organisms,
            reverse=True,
        )
        ordered_nodes += ordered_nodes_s
        partition_separators.append(
            (partition_name, partition_separators[-1][1] + len(ordered_nodes_s))
        )

    # Append cloud partition to the ordered nodes list
    ordered_nodes += ordered_nodes_c
    partition_separators.append(
        ("Cloud", partition_separators[-1][1] + len(ordered_nodes_c))
    )

    return ordered_nodes, partition_separators


def create_partition_shapes(
    separators: List[Tuple[str, float]],
    xval_max: float,
    heatmap_row: int,
    partition_to_color: Dict[str, str],
) -> List[dict]:
    """
    Create the shapes for plot separators to visually distinguish partitions in the plot.

    :param separators: A list of tuples containing partition names and their corresponding separator positions.
    :param xval_max: The maximum x-value for the plot.
    :param heatmap_row: The row number of the heatmap.
    :param partition_to_color: A dictionary mapping partition names to their corresponding colors.
    :return: A list of shape dictionaries for Plotly to use in the plot.
    """

    shapes = []
    sep_prec = 0
    xref = "x1"
    yref = f"y{heatmap_row}"

    for partition_name, sep in separators:
        color = partition_to_color[partition_name]

        # Left vertical line for partition separator
        shapes.append(
            dict(
                type="line",
                x0=-1,
                x1=-1,
                y0=sep_prec,
                y1=sep,
                line=dict(width=10, color=color),
                xref=xref,
                yref=yref,
                name=partition_name,
                showlegend=True,
                legendgroup=partition_name,
            )
        )

        # Right vertical line for partition separator
        shapes.append(
            dict(
                type="line",
                x0=xval_max,
                x1=xval_max,
                y0=sep_prec,
                y1=sep,
                line=dict(width=10, color=color),
                xref=xref,
                yref=yref,
                name=partition_name,
                showlegend=False,
                legendgroup=partition_name,
            )
        )

        # Horizontal line across the partition boundary
        shapes.append(
            dict(
                type="line",
                x0=-1,
                x1=xval_max,
                y0=sep,
                y1=sep,
                line=dict(width=1, color=color),
                xref=xref,
                yref=yref,
                name=partition_name,
                showlegend=False,
                legendgroup=partition_name,
            )
        )

        sep_prec = sep

    return shapes


def metadata_stringify(gene) -> str:
    """
    Convert gene metadata to a formatted string.

    :param gene: The gene object with potential metadata.
    :return: A formatted string containing gene metadata information.
    """
    metadata_str = ""
    if gene.has_metadata():
        metadata_str = f"<br><br>{gene.ID} metadata"
        for metadata in gene.metadata:
            metadata_str += f"<br>metadata source: {metadata.source}<br>"
            metadata_dict = metadata.to_dict()
            metadata_str += "<br>".join(
                (f"{key}: {value}" for key, value in metadata_dict.items())
            )

    return metadata_str


def get_heatmap_hover_text(
    ordered_families: List, order_organisms: List
) -> List[List[str]]:
    """
    Generate hover text for the heatmap cells.

    :param ordered_families: The list of ordered gene families.
    :param order_organisms: The list of ordered organisms.
    :return: A 2D list of strings representing hover text for each heatmap cell.
    """
    text_data = []

    for family in ordered_families:
        text_per_family = []

        for org in order_organisms:
            if org in family.organisms:
                # gene_count = len(list(family.get_genes_per_org(org)))
                genes = ";".join(map(str, family.get_genes_per_org(org)))
                names = ";".join(
                    (gene.name for gene in family.get_genes_per_org(org) if gene.name)
                )

                # Compile additional information about genes
                extra_gene_info = f"genes:{genes}"
                if names:
                    extra_gene_info += f"<br>names:{names}"

                metadata = "<br>".join(
                    (
                        metadata_stringify(gene)
                        for gene in family.get_genes_per_org(org)
                        if gene.has_metadata()
                    )
                )
                extra_gene_info += metadata
            else:
                # gene_count = 0
                extra_gene_info = (
                    np.nan
                )  # Using np.nan instead of numpy.nan for consistency with numpy import

            # To get a more explicit hover. But is quite heavy on the finam html
            # gene_info = f"genome:{org.name}<br>family:{family.name}<br>gene_count:{gene_count}<br>{extra_gene_info}"
            # Light version:
            gene_info = extra_gene_info
            text_per_family.append(gene_info)

        text_data.append(text_per_family)

    return text_data


def create_tile_plot(
    binary_data: List[List[float]],
    text_data: List[List[str]],
    fam_order: List[str],
    partition_separator: List[tuple],
    order_organisms: List[
        Organism
    ],  # Replace 'Any' with the appropriate type if available
    dendrogram_fig: go.Figure,
    draw_dendrogram: bool,
) -> go.Figure:
    """
    Create the heatmap tile plot using Plotly.

    :param binary_data: The binary presence-absence matrix data.
    :param text_data: Hover text data for each cell in the heatmap.
    :param fam_order: List of gene family names in the desired order.
    :param partition_separator: List of tuples containing partition names and their separator positions.
    :param order_organisms: List of organisms in the desired order.
    :param dendrogram_fig: Plotly figure object for the dendrogram.
    :param draw_dendrogram: Flag indicating whether to draw the dendrogram.
    :return: A Plotly Figure object representing the tile plot.
    """

    xaxis_values = [org.name for org in order_organisms]

    heatmap_color = {"presence": "#005AB5", "multicopy": "#DC3220"}  # blue  # red

    green_colors = [
        "rgb(35,139,69)",
        "rgb(65,171,93)",
        "rgb(116,196,118)",
        "rgb(161,217,155)",
        "rgb(199,233,192)",
        "rgb(229,245,224)",
    ]

    shell_color_generator = cycle(green_colors)

    partition_to_color = {
        "Persistent": "#F7A507",
        "Cloud": "#79DEFF",
        "Shell_S1": "#00D860",
    }
    partition_to_color.update(
        {
            partition: next(shell_color_generator)
            for partition, _ in partition_separator
            if partition not in partition_to_color
        }
    )

    heatmap = [
        go.Heatmap(
            z=binary_data,
            x=xaxis_values,
            y=fam_order,
            text=text_data,
            zauto=False,
            zmin=0,
            zmax=2,
            autocolorscale=False,
            hovertemplate="genome: %{x}<br>family: %{y}<br>gene_count: %{z}<br>%{text}   <extra></extra>",
            colorscale=[
                [0, "#ffffff"],
                [0.33, "#ffffff"],
                [0.33, heatmap_color["presence"]],
                [0.66, heatmap_color["presence"]],
                [0.66, heatmap_color["multicopy"]],
                [1, heatmap_color["multicopy"]],
            ],
            colorbar=dict(
                title="Presence/Absence",
                titleside="top",
                tickmode="array",
                tickvals=[0.33, 1, 1.66],
                ticktext=["Absence", "Presence", "Multicopy"],
                len=0.27,
                outlinecolor="black",
                outlinewidth=0.5,
                ticks=None,
                orientation="v",
            ),
        )
    ]

    if draw_dendrogram:
        heatmap_row = 2
        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.01,
            row_heights=[0.2, 0.8],
        )

        for data in dendrogram_fig["data"]:
            fig.add_trace(data, row=1, col=1)

    else:
        heatmap_row = 1
        fig = make_subplots(rows=1, cols=1)

    heatmap[0]["x"] = dendrogram_fig["layout"]["xaxis"]["tickvals"]

    for data in heatmap:

        fig.add_trace(data, row=heatmap_row, col=1)

    layout = go.Layout(title="Presence-Absence Matrix", plot_bgcolor="#ffffff")

    if draw_dendrogram:
        layout.xaxis2 = dendrogram_fig.layout.xaxis
    else:
        layout.xaxis = dendrogram_fig.layout.xaxis

    fig.update_layout(layout)

    fig.update_xaxes(ticklen=0, title="Genomes")
    fig.update_yaxes(
        ticklen=0,
        title="Gene Families",
        tickfont=dict(size=10),
        automargin=True,
    )
    if draw_dendrogram:

        fig.layout.yaxis.title = None
        fig.layout.yaxis2.title = dict(text="Gene Families")
        fig.layout.xaxis.title = None

    xmax = (
        dendrogram_fig["layout"]["xaxis"]["tickvals"][-1]
        + dendrogram_fig["layout"]["xaxis"]["tickvals"][0]
        + 0.5
    )
    shapes = create_partition_shapes(
        partition_separator, xmax, heatmap_row, partition_to_color
    )

    fig.update_layout(
        go.Layout(
            shapes=shapes,
            showlegend=True,
        )
    )

    fig.update_layout(
        legend=dict(
            title="Family Partition",
            traceorder="reversed",
        )
    )

    fig.update_layout(
        {
            "width": 1000,
            "height": 1000,
        }
    )

    return fig
