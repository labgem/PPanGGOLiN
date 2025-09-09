# default libraries
import logging
from collections import defaultdict
from pathlib import Path
from typing import Tuple, Dict, Union

# installed libraries
import plotly.graph_objects as go

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.pangenome import Pangenome


def compute_family_counts(
    pangenome: Pangenome,
) -> Tuple[Dict[int, Dict[str, int]], int, bool, bool, Union[float, str]]:
    """
    Compute gene family distribution across genomes, and estimate the total family richness using Chao1.

    This function computes:
      - The number of gene families found in exactly `n` genomes (for each `n`)
      - Their counts by partition (e.g. persistent, shell, cloud, etc.)
      - Whether the pangenome is partitioned
      - Whether any gene family has an undefined partition
      - A Chao1 estimate of the total number of gene families (to account for unobserved families)

    :param pangenome: A Pangenome object with annotated and partitioned gene families
    :return: A tuple containing:
        - family_count_by_org_and_part: dict[int][str] → count of gene families by number of organisms and partition
        - max_family_count: the highest number of gene families for any given organism count (for scaling plots)
        - is_partitioned: whether at least one gene family has a defined partition
        - has_undefined: whether any gene family has the undefined ("U") partition
        - chao: Chao1 estimate for total gene family richness (float or "NA" if not computable)
    """

    family_count_by_org_and_part = defaultdict(lambda: defaultdict(int))
    is_partitioned = False
    has_undefined = False
    max_family_count = 0

    for fam in pangenome.gene_families:
        nb_org = fam.number_of_organisms

        if fam.partition:
            is_partitioned = True
            if fam.partition == "U":
                has_undefined = True
            family_count_by_org_and_part[nb_org][fam.named_partition.capitalize()] += 1

        family_count_by_org_and_part[nb_org]["pangenome"] += 1
        max_family_count = max(
            max_family_count, family_count_by_org_and_part[nb_org]["pangenome"]
        )

    singleton = family_count_by_org_and_part[1]["pangenome"]
    doubleton = family_count_by_org_and_part[2]["pangenome"]
    chao: Union[float, str] = "NA"

    if singleton > 0 and doubleton > 0:
        chao = round(
            pangenome.number_of_gene_families + (singleton**2 / (2 * doubleton)),
            2,
        )

    return (
        family_count_by_org_and_part,
        max_family_count,
        is_partitioned,
        has_undefined,
        chao,
    )


def build_ucurve_plot(
    number_of_organisms: int,
    count: Dict[int, Dict[str, int]],
    max_bar: int,
    is_partitioned: bool,
    has_undefined: bool,
    chao: float,
    soft_core: float,
    number_of_gene_families: int,
) -> go.Figure:
    """
    Build a U-curve bar plot of gene family frequencies across genomes.

    This function visualizes the distribution of gene families (e.g., persistent, shell, cloud)
    across genomes in the pangenome using a stacked bar chart. It supports both partitioned
    and non-partitioned pangenomes, and highlights the soft core threshold as a vertical dashed line.

    :param number_of_organism: Number of organisms in the pangenome.
    :param count: Nested dictionary mapping genome counts to gene family counts by partition.
    :param max_bar: Maximum y-axis value to set the plot height.
    :param is_partitioned: Flag indicating if the pangenome is partitioned into persistent, shell, cloud.
    :param has_undefined: Flag indicating presence of undefined gene families.
    :param chao: Chao1 richness estimator value, shown in the plot title.
    :param soft_core: Threshold fraction defining the soft core genome, used to draw a vertical line.
    :param number_of_gene_families: Total number of gene families in the pangenome.

    :return: A Plotly Figure object containing the U-curve bar plot.
    """
    # Color definitions for each gene family category
    color_palette = {
        "pangenome": "black",
        "exact_accessory": "#EB37ED",
        "exact_core": "#FF2828",
        "soft_core": "#c7c938",
        "soft_accessory": "#996633",
        "Shell": "#00D860",
        "Persistent": "#F7A507",
        "Cloud": "#79DEFF",
        "undefined": "#828282",
    }

    data_plot = []

    if is_partitioned and not has_undefined:
        # Build stacked bars for persistent, shell, and cloud partitions
        for part in ["Persistent", "Shell", "Cloud"]:
            values = [
                count[nb_org][part] for nb_org in range(1, number_of_organisms + 1)
            ]
            data_plot.append(
                go.Bar(
                    x=list(range(1, number_of_organisms + 1)),
                    y=values,
                    name=part,
                    marker=dict(color=color_palette.get(part, "grey")),
                )
            )
    else:
        # For non-partitioned or with undefined families, plot a single category bar
        text = "undefined" if has_undefined else "pangenome"
        values = [count[nb_org][text] for nb_org in range(1, number_of_organisms + 1)]
        data_plot.append(
            go.Bar(
                x=list(range(1, number_of_organisms + 1)),
                y=values,
                name=text,
                marker=dict(color=color_palette.get(text, "grey")),
            )
        )

    # Calculate position of the soft core threshold line (vertical)
    soft_core_pos = number_of_organisms * soft_core

    layout = go.Layout(
        title=dict(
            text="Gene families frequency distribution (U shape)",
            # font=dict(size=16),
            subtitle=dict(
                text=f"Gene Families: {number_of_gene_families} | Chao1: {chao}",
                font=dict(color="gray", size=14),
            ),
        ),
        xaxis=dict(title="Occurring in x genomes"),
        yaxis=dict(title="# of gene families (F)"),
        barmode="stack",
        shapes=[
            dict(
                type="line",
                x0=soft_core_pos,
                x1=soft_core_pos,
                y0=0,
                y1=max_bar,
                line=dict(dash="dashdot", color="grey"),
                # This line marks the soft core genome threshold
            )
        ],
        plot_bgcolor="#ffffff",
    )

    return go.Figure(data=data_plot, layout=layout)


def draw_ucurve(
    pangenome: Pangenome,
    output: Path,
    soft_core: float = 0.95,
    disable_bar: bool = False,
):
    """
    Draws the U-shaped curve of gene family frequency distribution.

    :param pangenome: Partitioned pangenome
    :param output: Path to output directory
    :param soft_core: Soft core threshold to use
    :param disable_bar: Allow to disable progress bar
    """
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=False,
        disable_bar=disable_bar,
    )

    logging.getLogger("PPanGGOLiN").info("Drawing the U-shaped curve...")

    count, max_bar, is_partitioned, has_undefined, chao = compute_family_counts(
        pangenome
    )
    fig = build_ucurve_plot(
        pangenome.number_of_organisms,
        count,
        max_bar,
        is_partitioned,
        has_undefined,
        chao,
        soft_core,
        pangenome.number_of_gene_families,
    )

    output_file = output / "Ushaped_plot.html"
    fig.write_html(output_file)

    logging.getLogger("PPanGGOLiN").info(
        f"Done drawing the U-shaped curve: '{output_file}'"
    )
