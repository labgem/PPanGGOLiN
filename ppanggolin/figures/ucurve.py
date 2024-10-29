# default libraries
import logging
from collections import defaultdict
from pathlib import Path

# installed libraries
import plotly.graph_objs as go
import plotly.offline as out_plotly

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.pangenome import Pangenome


def draw_ucurve(
    pangenome: Pangenome,
    output: Path,
    soft_core: float = 0.95,
    disable_bar: bool = False,
):
    """

    :param pangenome: Partitioned pangenome
    :param output: Path to output directory
    :param soft_core: Soft core threshold to use
    :param disable_bar: Allow to disable progress bar
    :return:
    """
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=disable_bar,
    )
    logging.getLogger("PPanGGOLiN").info("Drawing the U-shaped curve...")
    max_bar = 0
    count = defaultdict(lambda: defaultdict(int))
    is_partitioned = False
    has_undefined = False
    for fam in pangenome.gene_families:
        nb_org = fam.number_of_organisms
        if fam.partition != "":
            is_partitioned = True
            if fam.partition == "U":
                has_undefined = True
            count[nb_org][fam.named_partition] += 1
        count[nb_org]["pangenome"] += 1
        max_bar = (
            count[nb_org]["pangenome"]
            if count[nb_org]["pangenome"] > max_bar
            else max_bar
        )
    data_plot = []
    chao = "NA"
    if count[1]["pangenome"] > 0:
        chao = round(
            pangenome.number_of_gene_families
            + (pow(count[1]["pangenome"], 2) / (count[2]["pangenome"] * 2)),
            2,
        )
    colors = {
        "pangenome": "black",
        "exact_accessory": "#EB37ED",
        "exact_core": "#FF2828",
        "soft_core": "#c7c938",
        "soft_accessory": "#996633",
        "shell": "#00D860",
        "persistent": "#F7A507",
        "cloud": "#79DEFF",
        "undefined": "#828282",
    }

    if is_partitioned and not has_undefined:
        persistent_values = []
        shell_values = []
        cloud_values = []
        for nb_org in range(1, pangenome.number_of_organisms + 1):
            persistent_values.append(count[nb_org]["persistent"])
            shell_values.append(count[nb_org]["shell"])
            cloud_values.append(count[nb_org]["cloud"])
        data_plot.append(
            go.Bar(
                x=list(range(1, pangenome.number_of_organisms + 1)),
                y=persistent_values,
                name="persistent",
                marker=dict(color=colors["persistent"]),
            )
        )
        data_plot.append(
            go.Bar(
                x=list(range(1, pangenome.number_of_organisms + 1)),
                y=shell_values,
                name="shell",
                marker=dict(color=colors["shell"]),
            )
        )
        data_plot.append(
            go.Bar(
                x=list(range(1, pangenome.number_of_organisms + 1)),
                y=cloud_values,
                name="cloud",
                marker=dict(color=colors["cloud"]),
            )
        )
    else:
        text = "undefined" if has_undefined else "pangenome"
        undefined_values = [
            count[nb_org][text]
            for nb_org in range(1, pangenome.number_of_organisms + 1)
        ]

        data_plot.append(
            go.Bar(
                x=list(range(1, pangenome.number_of_organisms + 1)),
                y=undefined_values,
                name=text,
                marker=dict(color=colors[text]),
            )
        )
    x = pangenome.number_of_organisms * soft_core
    layout = go.Layout(
        title="Gene families frequency distribution (U shape), chao=" + str(chao),
        xaxis=dict(title="Occurring in x genomes"),
        yaxis=dict(title="# of gene families (F)"),
        barmode="stack",
        shapes=[
            dict(
                type="line",
                x0=x,
                x1=x,
                y0=0,
                y1=max_bar,
                line=dict(dict(width=5, dash="dashdot", color="grey")),
            )
        ],
        plot_bgcolor="#ffffff",
    )

    fig = go.Figure(data=data_plot, layout=layout)
    out_plotly.plot(
        fig, filename=output.as_posix() + "/Ushaped_plot.html", auto_open=False
    )
    logging.getLogger("PPanGGOLiN").info(
        f"Done drawing the U-shaped curve : '{output.as_posix() + '/Ushaped_plot.html'}'"
    )
