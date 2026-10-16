#!/usr/bin/env python3

# default libraries
import argparse
import logging
from collections import Counter
import random
import time
from multiprocessing import get_context
import os
import warnings
from pathlib import Path
from typing import Tuple, Dict

# installed libraries
from tqdm import tqdm
import gmpy2
import numpy as np
from pandas import Series, read_csv
import plotly.offline as out_plotly
import plotly.graph_objs as go
import scipy.optimize as optimization

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir
from ppanggolin.formats import check_pangenome_info
from ppanggolin.nem.model import NemResult, solve
from ppanggolin.nem.partition import evaluate_nb_partitions
from ppanggolin.nem.nem_input import (
    NemIndex,
    build_nem_index,
    build_nem_input,
    subset_nem_index,
)


def partition_sample(
    sample_id: int,
    sample: NemIndex,
    beta: float = 2.5,
    max_degree_smoothing: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    k: int = -1,
    krange: list = None,
    seed: int = 42,
) -> Tuple[Dict[str, int], int]:
    """

    :param sample_id: Index of the sample group organisms
    :param sample: NEM index of the sample group organisms
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param max_degree_smoothing:  Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param chunk_size: Size of the chunks when performing partitioning using chunks of organisms.
    :param k: Number of partitions to use
    :param krange: Range of K values to test when detecting K automatically.
    :param seed: seed used to generate random numbers

    :return: Count of each partition and parameters for the given sample index
    """
    genomes = set(sample.org_names)

    krange = [3, 20] if krange is None else krange

    if k < 3:
        k = evaluate_nb_partitions(
            index=sample,
            max_degree_smoothing=max_degree_smoothing,
            free_dispersion=free_dispersion,
            chunk_size=chunk_size,
            krange=krange,
            seed=seed,
        )

    if len(genomes) <= chunk_size:  # fits in one chunk
        nem_input = build_nem_input(sample, max_degree_smoothing)
        votes = solve(
            nem_input, nem_input.scaled_beta(beta), free_dispersion, k
        ).partitions
    else:  # going to need multiple partitioning for this sample...
        votes = {fam: {"P": 0, "S": 0, "C": 0, "U": 0} for fam in sample.fam_names}
        validated = set()

        def count_votes(
            result: NemResult,
        ):
            """
            Validate partition assignation to families

            :param result: Partitioning results
            """
            for family, label in result.partitions.items():
                votes[family][label[0]] += 1
                nb_votes = sum(votes[family].values())
                if (
                    nb_votes > len(genomes) / chunk_size
                    and max(votes[family].values()) >= nb_votes * 0.5
                ) or (nb_votes > len(genomes)):
                    if family not in validated:
                        if max(votes[family].values()) < nb_votes * 0.5:
                            votes[family]["U"] = len(genomes)
                        validated.add(family)

        rng = random.Random(f"{seed}-{sample_id}")  # workers are reseeded at fork
        times_drawn = Counter()
        for org in genomes:
            times_drawn[org] = 0
        min_draws = len(genomes) / chunk_size

        while len(validated) < len(votes):
            chunks = []

            while not all(val >= min_draws for val in times_drawn.values()):
                # each family must be tested at least len(select_organisms)/chunk_size times.
                shuffled = sorted(genomes)
                rng.shuffle(shuffled)
                while len(shuffled) > chunk_size:
                    chunks.append(set(shuffled[:chunk_size]))
                    for org in chunks[-1]:
                        times_drawn[org] += 1
                    shuffled = shuffled[chunk_size:]
            for chunk in chunks:
                nem_input = build_nem_input(
                    subset_nem_index(sample, chunk), max_degree_smoothing
                )
                count_votes(
                    solve(nem_input, nem_input.scaled_beta(beta), free_dispersion, k)
                )
            min_draws += 1  # if families are left to validate, resample more
    if len(votes) == 0:
        counts = {
            "persistent": "NA",
            "shell": "NA",
            "cloud": "NA",
            "undefined": "NA",
            "K": k,
        }
    else:
        counts = {"persistent": 0, "shell": 0, "cloud": 0, "undefined": 0, "K": k}

        for val in votes.values():
            if isinstance(val, str):
                label = val
            else:
                label = max(val, key=val.get)
            if label.startswith("P"):
                counts["persistent"] += 1
            elif label.startswith("C"):
                counts["cloud"] += 1
            elif label.startswith("S"):
                counts["shell"] += 1
            else:
                counts["undefined"] += 1
    return counts, sample_id


def _partition_sample_worker(
    args: Tuple[int, NemIndex, float, int, bool, int, int, list, int]
) -> Tuple[Tuple[Dict[str, int], int]]:
    """
    Launch partition_sample in multiprocessing

    :param args: {index: int, sample: NemIndex, beta: float, max_degree_smoothing: int, free_dispersion: bool,
                  chunk_size: int, k: int, krange: list, seed: int}
    :return: Count of each partition and parameters for the given sample index
    """
    return partition_sample(*args)


def draw_curve(output: Path, data: list, max_sampling: int = 10):
    """
    Draw the rarefaction curve and associated data

    :param output: output directory path to draw the rarefaction curve and associated data
    :param max_sampling: Maximum number of organisms in a sample
    :param data:
    """
    logging.getLogger("PPanGGOLiN").info("Drawing the rarefaction curve ...")
    raref_name = output / "rarefaction.csv"
    raref = open(raref_name, "w")
    raref.write(
        ",".join(
            [
                "genomes_count",
                "persistent",
                "shell",
                "cloud",
                "undefined",
                "exact_core",
                "exact_accessory",
                "soft_core",
                "soft_accessory",
                "pangenome",
                "K",
            ]
        )
        + "\n"
    )
    for part in data:
        raref.write(
            ",".join(
                map(
                    str,
                    [
                        part["nborgs"],
                        part["persistent"],
                        part["shell"],
                        part["cloud"],
                        part["undefined"],
                        part["exact_core"],
                        part["exact_accessory"],
                        part["soft_core"],
                        part["soft_accessory"],
                        part["exact_core"] + part["exact_accessory"],
                        part["K"],
                    ],
                )
            )
            + "\n"
        )
    raref.close()

    def heap_law(n, p_kappa, p_gamma) -> float:
        return p_kappa * n**p_gamma

    def poly_area(p_x: list, p_y: list) -> float:
        return 0.5 * np.abs(np.dot(p_x, np.roll(p_y, 1)) - np.dot(p_y, np.roll(p_x, 1)))

    annotations = []
    traces = []
    data_raref = read_csv(raref_name, index_col=False)
    params_file = open(output / "rarefaction_parameters.csv", "w")
    params_file.write(
        "partition,kappa,gamma,kappa_std_error,gamma_std_error,IQR_area\n"
    )
    for partition in [
        "persistent",
        "shell",
        "cloud",
        "undefined",
        "exact_core",
        "exact_accessory",
        "soft_core",
        "soft_accessory",
        "pangenome",
    ]:
        percentiles_75 = Series(
            {
                i: np.nanpercentile(
                    data_raref[data_raref["genomes_count"] == i][partition], 75
                )
                for i in range(1, max_sampling + 1)
            }
        ).dropna()
        percentiles_25 = Series(
            {
                i: np.nanpercentile(
                    data_raref[data_raref["genomes_count"] == i][partition], 25
                )
                for i in range(1, max_sampling + 1)
            }
        ).dropna()
        mins = Series(
            {
                i: np.min(data_raref[data_raref["genomes_count"] == i][partition])
                for i in range(1, max_sampling + 1)
            }
        ).dropna()
        maxs = Series(
            {
                i: np.max(data_raref[data_raref["genomes_count"] == i][partition])
                for i in range(1, max_sampling + 1)
            }
        ).dropna()
        medians = Series(
            {
                i: np.median(data_raref[data_raref["genomes_count"] == i][partition])
                for i in range(1, max_sampling + 1)
            }
        ).dropna()
        means = Series(
            {
                i: np.mean(data_raref[data_raref["genomes_count"] == i][partition])
                for i in range(1, max_sampling + 1)
            }
        ).dropna()
        initial_kappa_gamma = np.array([0.0, 0.0])
        x = percentiles_25.index.tolist()
        x += list(reversed(percentiles_25.index.tolist()))
        area_iqr = poly_area(x, percentiles_25.tolist() + percentiles_75.tolist())
        nb_org_min_fitting = 15
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
        try:
            all_values = data_raref[data_raref["genomes_count"] > nb_org_min_fitting][
                partition
            ].dropna()
            res = optimization.curve_fit(
                heap_law,
                data_raref.loc[all_values.index]["genomes_count"],
                all_values,
                initial_kappa_gamma,
            )
            kappa, gamma = res[0]
            error_k, error_g = np.sqrt(
                np.diag(res[1])
            )  # to calculate the fitting error.
            # The variance of parameters are the diagonal elements of the variance-co variance matrix,
            # and the standard error is the square root of it. source :
            # https://stackoverflow.com/questions/25234996/getting-standard-error-associated-with-parameter-estimates-from-scipy-optimize-c
            if np.isinf(error_k) and np.isinf(error_g):
                params_file.write(
                    ",".join([partition, "NA", "NA", "NA", "NA", str(area_iqr)]) + "\n"
                )
            else:
                params_file.write(
                    ",".join(
                        [
                            partition,
                            str(kappa),
                            str(gamma),
                            str(error_k),
                            str(error_g),
                            str(area_iqr),
                        ]
                    )
                    + "\n"
                )
                regression = np.apply_along_axis(
                    heap_law,
                    0,
                    range(nb_org_min_fitting + 1, max_sampling + 1),
                    kappa,
                    gamma,
                )
                regression_sd_top = np.apply_along_axis(
                    heap_law,
                    0,
                    range(nb_org_min_fitting + 1, max_sampling + 1),
                    kappa - error_k,
                    gamma + error_g,
                )
                regression_sd_bottom = np.apply_along_axis(
                    heap_law,
                    0,
                    range(nb_org_min_fitting + 1, max_sampling + 1),
                    kappa + error_k,
                    gamma - error_g,
                )
                traces.append(
                    go.Scatter(
                        x=list(range(nb_org_min_fitting + 1, max_sampling + 1)),
                        y=regression,
                        name=partition + ": Heaps' law",
                        line=dict(color=colors[partition], width=4, dash="dash"),
                        visible="legendonly" if partition == "undefined" else True,
                    )
                )
                traces.append(
                    go.Scatter(
                        x=list(range(nb_org_min_fitting + 1, max_sampling + 1)),
                        y=regression_sd_top,
                        name=partition + ": Heaps' law error +",
                        line=dict(color=colors[partition], width=1, dash="dash"),
                        visible="legendonly" if partition == "undefined" else True,
                    )
                )
                traces.append(
                    go.Scatter(
                        x=list(range(nb_org_min_fitting + 1, max_sampling + 1)),
                        y=regression_sd_bottom,
                        name=partition + ": Heaps' law error -",
                        line=dict(color=colors[partition], width=1, dash="dash"),
                        visible="legendonly" if partition == "undefined" else True,
                    )
                )
                annotations.append(
                    dict(
                        x=max_sampling,
                        y=heap_law(max_sampling, kappa, gamma),
                        ay=0,
                        ax=50,
                        text="F="
                        + str(round(kappa, 0))
                        + "N"
                        + "<sup>"
                        + str(round(gamma, 5))
                        + "</sup><br>IQRarea="
                        + str(round(area_iqr, 2)),
                        showarrow=True,
                        arrowhead=7,
                        font=dict(size=10, color="white"),
                        align="center",
                        arrowcolor=colors[partition],
                        bordercolor="#c7c7c7",
                        borderwidth=2,
                        borderpad=4,
                        bgcolor=colors[partition],
                        opacity=0.8,
                    )
                )
        except (TypeError, RuntimeError, ValueError):  # if fitting doesn't work
            params_file.write(
                ",".join([partition, "NA", "NA", "NA", "NA", str(area_iqr)]) + "\n"
            )

        traces.append(
            go.Scatter(
                x=medians.index,
                y=medians,
                name=partition + " : medians",
                mode="lines+markers",
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=maxs.subtract(medians),
                    arrayminus=medians.subtract(mins),
                    visible=True,
                    color=colors[partition],
                    thickness=0.5,
                ),
                line=dict(color=colors[partition], width=1),
                marker=dict(color=colors[partition], symbol=3, size=8, opacity=0.5),
                visible="legendonly" if partition == "undefined" else True,
            )
        )
        traces.append(
            go.Scatter(
                x=means.index,
                y=means,
                name=partition + " : means",
                mode="markers",
                marker=dict(color=colors[partition], symbol=4, size=8, opacity=0.5),
                visible="legendonly" if partition == "undefined" else True,
            )
        )
        # up = percentiles_75
        # down = percentiles_25
        # IQR_area = up.append(down[::-1])
        # traces.append(go.Scatter(x=IQR_area.index,
        #                          y=IQR_area,
        #                          name = "IQR",
        #                          fill='toself',
        #                          mode="lines",
        #                          hoveron="points",
        #                          #hovertext=[str(round(e)) for e in half_stds.multiply(2)],
        #                          line=dict(color=COLORS[partition]),
        #                          marker=dict(color = COLORS[partition]),
        #                          visible = "legendonly" if partition == "undefined" else True))
        traces.append(
            go.Scatter(
                x=percentiles_75.index,
                y=percentiles_75,
                name=partition + " : 3rd quartile",
                mode="lines",
                hoveron="points",
                # hovertext=[str(round(e)) for e in half_stds.multiply(2)],
                line=dict(color=colors[partition]),
                marker=dict(color=colors[partition]),
                visible="legendonly" if partition == "undefined" else True,
            )
        )
        traces.append(
            go.Scatter(
                x=percentiles_25.index,
                y=percentiles_25,
                name=partition + " : 1st quartile",
                fill="tonexty",
                mode="lines",
                hoveron="points",
                # hovertext=[str(round(e)) for e in half_stds.multiply(2)],
                line=dict(color=colors[partition]),
                marker=dict(color=colors[partition]),
                visible="legendonly" if partition == "undefined" else True,
            )
        )
    layout = go.Layout(
        title="Rarefaction curve ",
        titlefont=dict(size=20),
        xaxis=dict(title="size of genome subsets (N)"),
        yaxis=dict(title="# of gene families (F)"),
        annotations=annotations,
        plot_bgcolor="#ffffff",
    )
    fig = go.Figure(data=traces, layout=layout)
    out_plotly.plot(
        fig, filename=output.as_posix() + "/rarefaction_curve.html", auto_open=False
    )
    params_file.close()


def make_rarefaction_curve(
    pangenome: Pangenome,
    output: Path,
    beta: float = 2.5,
    depth: int = 30,
    min_sampling: int = 1,
    max_sampling: int = 100,
    max_degree_smoothing: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    k: int = -1,
    krange: list = None,
    cpu: int = 1,
    seed: int = 42,
    reestimate_k: bool = False,
    soft_core: float = 0.95,
    disable_bar: bool = False,
):
    """
    Main function to make the rarefaction curve

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :param output: output directory path to draw the rarefaction curve and associated data
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param depth: Number of samplings at each sampling point
    :param min_sampling: Minimum number of organisms in a sample
    :param max_sampling: Maximum number of organisms in a sample
    :param max_degree_smoothing: Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param chunk_size: Size of the chunks when performing partitioning using chunks of organisms.
    :param k: Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically.
    :param krange: Range of K values to test when detecting K automatically.
    :param cpu: Number of available core
    :param seed: seed used to generate random numbers
    :param reestimate_k: recompute the number of partitions for each sample between the values provided by krange
    :param soft_core: Soft core threshold
    :param disable_bar: Disable progress bar
    """
    krange = [3, -1] if krange is None else list(krange)
    try:
        krange[0] = (
            pangenome.parameters["partition"]["# final nb of partitions"]
            if krange[0] < 0
            else krange[0]
        )
        krange[1] = (
            pangenome.parameters["partition"]["# final nb of partitions"]
            if krange[1] < 0
            else krange[1]
        )
    except KeyError:
        krange = [3, 20]
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=disable_bar,
    )
    nem_index = build_nem_index(pangenome)

    if float(pangenome.number_of_organisms) < max_sampling:
        max_sampling = pangenome.number_of_organisms
    else:
        max_sampling = int(max_sampling)

    if k < 3 and reestimate_k is False:  # estimate K once and for all.
        try:
            k = pangenome.parameters["partition"]["# final nb of partitions"]
            logging.getLogger("PPanGGOLiN").info(f"Reuse the number of partitions {k}")
        except KeyError:
            logging.getLogger("PPanGGOLiN").info(
                "Estimating the number of partitions..."
            )
            k = evaluate_nb_partitions(
                index=nem_index,
                max_degree_smoothing=max_degree_smoothing,
                free_dispersion=free_dispersion,
                chunk_size=chunk_size,
                krange=krange,
                cpu=cpu,
                seed=seed,
            )
            logging.getLogger("PPanGGOLiN").info(
                f"The number of partitions has been evaluated at {k}"
            )

    logging.getLogger("PPanGGOLiN").info("Extracting samples ...")
    random.seed(seed)
    samples = []
    for i in range(min_sampling, max_sampling):  # each point
        for _ in range(depth):  # number of samples per points
            samples.append(set(random.sample(list(pangenome.organisms), i + 1)))
    logging.getLogger("PPanGGOLiN").info(
        f"Done sampling genomes in the pan, there are {len(samples)} samples"
    )
    rows = []

    logging.getLogger("PPanGGOLiN").info("Computing bitarrays for each family...")
    org_bit = pangenome.compute_family_bitarrays()
    logging.getLogger("PPanGGOLiN").info(
        "Done computing bitarrays. Comparing them to get exact and soft core stats "
        f"for {len(samples)} samples..."
    )
    bar = tqdm(
        range(len(samples) * pangenome.number_of_gene_families),
        unit="gene family",
        disable=disable_bar,
    )
    for sample in samples:
        # make the sample's organism bitarray.
        sample_bitarray = gmpy2.xmpz()  # pylint: disable=no-member
        for org in sample:
            sample_bitarray[org_bit[org]] = 1

        counts = Counter()
        counts["soft_core"] = 0
        counts["exact_core"] = 0
        counts["exact_accessory"] = 0
        counts["soft_accessory"] = 0
        for fam in pangenome.gene_families:
            nb_common = gmpy2.popcount(
                fam.bitarray & sample_bitarray
            )  # pylint: disable=no-member
            counts["nborgs"] = len(sample)
            if nb_common != 0:  # in that case the node 'does not exist'
                if nb_common == len(sample):
                    counts["exact_core"] += 1
                else:
                    counts["exact_accessory"] += 1

                if float(nb_common) >= len(sample) * soft_core:
                    counts["soft_core"] += 1
                else:
                    counts["soft_accessory"] += 1
            bar.update()
        rows.append(counts)
    bar.close()
    # done with frequency of each family for each sample.

    with get_context("fork").Pool(processes=cpu) as p:
        # launch partitioning
        logging.getLogger("PPanGGOLiN").info(" Partitioning all samples...")
        bar = tqdm(range(len(samples)), unit="samples partitioned", disable=disable_bar)
        order = list(range(len(samples)))
        random.shuffle(
            order
        )  # shuffling the processing so that the progress bar is closer to reality.
        args = (
            (
                sample_id,
                subset_nem_index(nem_index, [org.name for org in samples[sample_id]]),
                beta,
                max_degree_smoothing,
                free_dispersion,
                chunk_size,
                k,
                krange,
                seed,
            )
            for sample_id in order
        )
        for result in p.imap_unordered(_partition_sample_worker, args):
            rows[result[1]] = {**result[0], **rows[result[1]]}
            bar.update()
    bar.close()

    logging.getLogger("PPanGGOLiN").info("Done  partitioning everything")
    warnings.filterwarnings("ignore")
    draw_curve(output, rows, max_sampling)
    warnings.resetwarnings()
    logging.getLogger("PPanGGOLiN").info("Done making the rarefaction curves")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    make_rarefaction_curve(
        pangenome=pangenome,
        output=args.output,
        beta=args.beta,
        depth=args.depth,
        min_sampling=args.min,
        max_sampling=args.max,
        max_degree_smoothing=args.max_degree_smoothing,
        free_dispersion=args.free_dispersion,
        chunk_size=args.chunk_size,
        k=args.nb_of_partitions,
        krange=args.krange,
        cpu=args.cpu,
        seed=args.seed,
        reestimate_k=args.reestimate_K,
        soft_core=args.soft_core,
        disable_bar=args.disable_prog_bar,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "rarefaction",
        description="Compute the rarefaction curve of the pangenome",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.description = "Compute the rarefaction curve of the pangenome"
    parser.category = "Expert"
    parser_rarefaction(parser)
    return parser


def parser_rarefaction(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of graph command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p",
        "--pangenome",
        required=False,
        type=Path,
        help="Path to the pangenome .h5 file.",
    )

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "-b",
        "--beta",
        required=False,
        default=2.5,
        type=float,
        help="beta is the strength of the smoothing using the graph topology during  partitioning. "
        "0 will deactivate spatial smoothing.",
    )
    optional.add_argument(
        "--depth",
        required=False,
        default=30,
        type=int,
        help="Number of samplings at each sampling point",
    )
    optional.add_argument(
        "--min",
        required=False,
        default=1,
        type=int,
        help="Minimum number of genomes in a sample",
    )
    optional.add_argument(
        "--max",
        required=False,
        type=float,
        default=100,
        help="Maximum number of genomes in a sample (if above the number of provided genomes, "
        "the provided genomes will be the maximum)",
    )

    optional.add_argument(
        "-ms",
        "--max_degree_smoothing",
        required=False,
        default=10,
        type=float,
        help="max. degree of the nodes to be included in the smoothing process.",
    )
    optional.add_argument(
        "-o",
        "--output",
        required=False,
        type=Path,
        default=Path(
            f"ppanggolin_output{time.strftime('DATE%Y-%m-%d_HOUR%H.%M.%S', time.localtime())}"
            f"_PID{str(os.getpid())}"
        ),
        help="Output directory.",
    )
    optional.add_argument(
        "-fd",
        "--free_dispersion",
        required=False,
        default=False,
        action="store_true",
        help="use if the dispersion around the centroid vector of each partition during must be free."
        " It will be the same for all genomes by default.",
    )
    optional.add_argument(
        "-ck",
        "--chunk_size",
        required=False,
        default=500,
        type=int,
        help="Size of the chunks when performing partitioning using chunks of genomes. "
        "Chunk partitioning will be used automatically "
        "if the number of genomes is above this number.",
    )
    optional.add_argument(
        "-K",
        "--nb_of_partitions",
        required=False,
        default=-1,
        type=int,
        help="Number of partitions to use. Must be at least 2. "
        "By default reuse K if it exists else compute it.",
    )
    optional.add_argument(
        "--reestimate_K",
        required=False,
        action="store_true",
        help=" Will recompute the number of partitions for each sample "
        "(between the values provided by --krange) (VERY intensive. Can take a long time.)",
    )
    optional.add_argument(
        "-Kmm",
        "--krange",
        nargs=2,
        required=False,
        type=int,
        default=[3, -1],
        help="Range of K values to test when detecting K automatically. "
        "Default between 3 and the K previously computed "
        "if there is one, or 20 if there are none.",
    )
    optional.add_argument(
        "--soft_core",
        required=False,
        type=float,
        default=0.95,
        help="Soft core threshold",
    )
    optional.add_argument(
        "-se",
        "--seed",
        type=int,
        default=42,
        help="seed used to generate random numbers",
    )
    optional.add_argument(
        "-c",
        "--cpu",
        required=False,
        default=1,
        type=int,
        help="Number of available CPUs.",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_rarefaction(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
