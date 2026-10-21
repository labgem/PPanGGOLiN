#!/usr/bin/env python3

# default libraries
import argparse
import logging
import os
import random
import time
import warnings
from collections import Counter
from dataclasses import dataclass, replace
from multiprocessing import get_context
from pathlib import Path

# installed libraries
import gmpy2
import numpy as np
import plotly.graph_objs as go
import plotly.offline as out_plotly
import scipy.optimize as optimization
from pandas import DataFrame, Series, read_csv
from tqdm import tqdm

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.nem.nem_input import NemIndex, build_nem_index, subset_nem_index
from ppanggolin.nem.partition import (
    ChunkVotes,
    NemSettings,
    draw_chunks,
    _keep_index,
    evaluate_nb_partitions,
    inherited_index,
    partition_chunk,
)
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir

logger = logging.getLogger("PPanGGOLiN")

CURVE_COLUMNS = [
    "persistent",
    "shell",
    "cloud",
    "undefined",
    "exact_core",
    "exact_accessory",
    "soft_core",
    "soft_accessory",
    "pangenome",
]
COLORS = {
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
NB_GENOMES_MIN_FITTING = 15


@dataclass
class RarefactionOptions(NemSettings):
    """Parameters of the rarefaction command.

    :param chunk_size: Partition a sample in chunks of this many genomes when it has more.
    :param k: Number of partitions. Under 3, it is reused from the partition or estimated.
    :param krange: Range of k values tested when estimating k. Negative bounds take the pangenome's k.
    :param seed: Seed used to generate random numbers.
    :param depth: Number of samples drawn at each sample size.
    :param min_sampling: Smallest sample size (minus one for the range).
    :param max_sampling: Largest sample size.
    :param reestimate_k: Estimate k again for each sample.
    :param soft_core: Soft core threshold.
    """

    chunk_size: int
    k: int
    krange: list[int] | None
    seed: int
    depth: int
    min_sampling: int
    max_sampling: int
    reestimate_k: bool
    soft_core: float

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "RarefactionOptions":
        """Read the parameters from the arguments of the rarefaction command.

        :param args: Arguments of the rarefaction command.
        :return: The parameters.
        """
        return cls(
            beta=args.beta,
            max_degree_smoothing=args.max_degree_smoothing,
            free_dispersion=args.free_dispersion,
            chunk_size=args.chunk_size,
            k=args.nb_of_partitions,
            krange=args.krange,
            seed=args.seed,
            depth=args.depth,
            min_sampling=args.min,
            max_sampling=args.max,
            reestimate_k=args.reestimate_K,
            soft_core=args.soft_core,
        )


def _count_partitions(labels: dict[str, str], k: int) -> dict[str, int | str]:
    """Count the families of each partition.

    :param labels: Label of each family.
    :param k: Number of partitions.
    :return: Number of persistent, shell, cloud and undefined families, and k.
    """
    if len(labels) == 0:
        return {
            "persistent": "NA",
            "shell": "NA",
            "cloud": "NA",
            "undefined": "NA",
            "K": k,
        }
    counts = {"persistent": 0, "shell": 0, "cloud": 0, "undefined": 0, "K": k}
    for label in labels.values():
        if label.startswith("P"):
            counts["persistent"] += 1
        elif label.startswith("C"):
            counts["cloud"] += 1
        elif label.startswith("S"):
            counts["shell"] += 1
        else:
            counts["undefined"] += 1
    return counts


def partition_sample(
    sample_id: int, sample: NemIndex, k: int, options: RarefactionOptions
) -> tuple[dict[str, int | str], int]:
    """Partition the gene families of one sample of genomes.

    :param sample_id: Position of the sample.
    :param sample: NEM index of the sample.
    :param k: Number of partitions. Under 3, it is estimated on the sample.
    :param options: Rarefaction parameters, with the krange already resolved.
    :return: Number of families in each partition, and the sample's position.
    """
    genomes = sample.org_names
    if k < 3:
        k = evaluate_nb_partitions(
            index=sample,
            output=None,
            max_degree_smoothing=options.max_degree_smoothing,
            free_dispersion=options.free_dispersion,
            chunk_size=options.chunk_size,
            krange=options.krange,
            draw_icl=False,
            cpu=1,  # already in a pool worker
            seed=options.seed,
            disable_bar=True,
        )

    if len(genomes) <= options.chunk_size:
        labels = partition_chunk(sample, k, options).partitions
    else:
        votes = ChunkVotes(sample.fam_names, len(genomes), options.chunk_size)
        # forked workers are reseeded, so each sample has its own generator
        rng = random.Random(f"{options.seed}-{sample_id}")
        times_drawn = Counter({genome: 0 for genome in genomes})
        min_draws = len(genomes) / options.chunk_size
        while not votes.done:
            for chunk in draw_chunks(
                genomes, options.chunk_size, times_drawn, min_draws, rng
            ):
                result = partition_chunk(subset_nem_index(sample, chunk), k, options)
                votes.add(result.partitions)
            min_draws += 1  # if families are left to validate, draw more chunks
        labels = votes.labels()
    return _count_partitions(labels, k), sample_id


def _partition_sample_worker(
    args: tuple[int, list[str], int, RarefactionOptions],
) -> tuple[dict[str, int | str], int]:
    """Partition one sample of the index the pool inherited.

    :param args: (position of the sample, its genomes, k, options).
    :return: Number of families in each partition, and the sample's position.
    """
    sample_id, genomes, k, options = args
    return partition_sample(
        sample_id, subset_nem_index(inherited_index(), genomes), k, options
    )


def write_rarefaction_table(output: Path, rows: list[dict]) -> Path:
    """Write the family counts of each sample.

    :param output: Output directory.
    :param rows: Family counts of each sample.
    :return: Path of the table.
    """
    table = output / "rarefaction.csv"
    with open(table, "w") as out:
        out.write(
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
        for row in rows:
            out.write(
                ",".join(
                    map(
                        str,
                        [
                            row["nborgs"],
                            row["persistent"],
                            row["shell"],
                            row["cloud"],
                            row["undefined"],
                            row["exact_core"],
                            row["exact_accessory"],
                            row["soft_core"],
                            row["soft_accessory"],
                            row["exact_core"] + row["exact_accessory"],
                            row["K"],
                        ],
                    )
                )
                + "\n"
            )
    return table


def heap_law(n, p_kappa, p_gamma) -> float:
    """Heaps' law: number of families in n genomes.

    :param n: Number of genomes.
    :param p_kappa: Scale.
    :param p_gamma: Exponent.
    :return: The number of families.
    """
    return p_kappa * n**p_gamma


def _poly_area(p_x: list, p_y: list) -> float:
    """Area of a polygon, by the shoelace formula.

    :param p_x: Abscissas of the vertices.
    :param p_y: Ordinates of the vertices.
    :return: The area.
    """
    return 0.5 * np.abs(np.dot(p_x, np.roll(p_y, 1)) - np.dot(p_y, np.roll(p_x, 1)))


def _summary_series(
    data: DataFrame, column: str, max_sampling: int
) -> dict[str, Series]:
    """Summarize a column at each sample size.

    :param data: Rarefaction table.
    :param column: Column to summarize.
    :param max_sampling: Largest sample size.
    :return: First and third quartiles, minimum, maximum, median and mean at each size.
    """
    sizes = range(1, max_sampling + 1)
    values = {i: data[data["genomes_count"] == i][column] for i in sizes}
    return {
        "q3": Series({i: np.nanpercentile(values[i], 75) for i in sizes}).dropna(),
        "q1": Series({i: np.nanpercentile(values[i], 25) for i in sizes}).dropna(),
        "min": Series({i: np.min(values[i]) for i in sizes}).dropna(),
        "max": Series({i: np.max(values[i]) for i in sizes}).dropna(),
        "median": Series({i: np.median(values[i]) for i in sizes}).dropna(),
        "mean": Series({i: np.mean(values[i]) for i in sizes}).dropna(),
    }


def _fit_heaps_law(data: DataFrame, column: str) -> tuple[float, float, float, float]:
    """Fit Heaps' law on the samples larger than `NB_GENOMES_MIN_FITTING`.

    :param data: Rarefaction table.
    :param column: Column to fit.
    :raises RuntimeError: If the fit does not converge.
    :return: kappa, gamma and their standard errors.
    """
    values = data[data["genomes_count"] > NB_GENOMES_MIN_FITTING][column].dropna()
    params, covariance = optimization.curve_fit(
        heap_law,
        data.loc[values.index]["genomes_count"],
        values,
        np.array([0.0, 0.0]),
    )
    kappa, gamma = params
    error_k, error_g = np.sqrt(np.diag(covariance))
    return kappa, gamma, error_k, error_g


def _heaps_law_traces(
    column: str, fit: tuple, area_iqr: float, max_sampling: int
) -> tuple[list[go.Scatter], dict]:
    """Plot a Heaps' law fit with its error band.

    :param column: Fitted column.
    :param fit: kappa, gamma and their standard errors.
    :param area_iqr: Area between the quartile curves.
    :param max_sampling: Largest sample size.
    :return: The fit traces and their annotation.
    """
    kappa, gamma, error_k, error_g = fit
    sizes = range(NB_GENOMES_MIN_FITTING + 1, max_sampling + 1)
    visible = "legendonly" if column == "undefined" else True
    traces = []
    for name, width, p_kappa, p_gamma in (
        (": Heaps' law", 4, kappa, gamma),
        (": Heaps' law error +", 1, kappa - error_k, gamma + error_g),
        (": Heaps' law error -", 1, kappa + error_k, gamma - error_g),
    ):
        traces.append(
            go.Scatter(
                x=list(sizes),
                y=np.apply_along_axis(heap_law, 0, sizes, p_kappa, p_gamma),
                name=column + name,
                line=dict(color=COLORS[column], width=width, dash="dash"),
                visible=visible,
            )
        )
    annotation = dict(
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
        arrowcolor=COLORS[column],
        bordercolor="#c7c7c7",
        borderwidth=2,
        borderpad=4,
        bgcolor=COLORS[column],
        opacity=0.8,
    )
    return traces, annotation


def _summary_traces(column: str, summary: dict[str, Series]) -> list[go.Scatter]:
    """Plot the medians, means and quartiles of a column.

    :param column: Summarized column.
    :param summary: Output of `_summary_series`.
    :return: The traces.
    """
    color = COLORS[column]
    visible = "legendonly" if column == "undefined" else True
    medians = summary["median"]
    return [
        go.Scatter(
            x=medians.index,
            y=medians,
            name=column + " : medians",
            mode="lines+markers",
            error_y=dict(
                type="data",
                symmetric=False,
                array=summary["max"].subtract(medians),
                arrayminus=medians.subtract(summary["min"]),
                visible=True,
                color=color,
                thickness=0.5,
            ),
            line=dict(color=color, width=1),
            marker=dict(color=color, symbol=3, size=8, opacity=0.5),
            visible=visible,
        ),
        go.Scatter(
            x=summary["mean"].index,
            y=summary["mean"],
            name=column + " : means",
            mode="markers",
            marker=dict(color=color, symbol=4, size=8, opacity=0.5),
            visible=visible,
        ),
        go.Scatter(
            x=summary["q3"].index,
            y=summary["q3"],
            name=column + " : 3rd quartile",
            mode="lines",
            hoveron="points",
            line=dict(color=color),
            marker=dict(color=color),
            visible=visible,
        ),
        go.Scatter(
            x=summary["q1"].index,
            y=summary["q1"],
            name=column + " : 1st quartile",
            fill="tonexty",
            mode="lines",
            hoveron="points",
            line=dict(color=color),
            marker=dict(color=color),
            visible=visible,
        ),
    ]


def draw_curve(output: Path, rows: list[dict], max_sampling: int):
    """Write the rarefaction table, fit Heaps' law on each column and draw the curves.

    :param output: Output directory.
    :param rows: Family counts of each sample.
    :param max_sampling: Largest sample size.
    """
    logger.info("Drawing the rarefaction curve ...")
    data = read_csv(write_rarefaction_table(output, rows), index_col=False)
    annotations = []
    traces = []
    with open(output / "rarefaction_parameters.csv", "w") as params_file:
        params_file.write(
            "partition,kappa,gamma,kappa_std_error,gamma_std_error,IQR_area\n"
        )
        for column in CURVE_COLUMNS:
            summary = _summary_series(data, column, max_sampling)
            q1 = summary["q1"]
            x = q1.index.tolist() + list(reversed(q1.index.tolist()))
            area_iqr = _poly_area(x, q1.tolist() + summary["q3"].tolist())
            try:
                fit = _fit_heaps_law(data, column)
            except (TypeError, RuntimeError, ValueError):
                fit = None
            if fit is None or (np.isinf(fit[2]) and np.isinf(fit[3])):
                params_file.write(
                    ",".join([column, "NA", "NA", "NA", "NA", str(area_iqr)]) + "\n"
                )
            else:
                params_file.write(
                    ",".join([column, *map(str, fit), str(area_iqr)]) + "\n"
                )
                fit_traces, annotation = _heaps_law_traces(
                    column, fit, area_iqr, max_sampling
                )
                traces += fit_traces
                annotations.append(annotation)
            traces += _summary_traces(column, summary)

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


def _resolve_krange(pangenome: Pangenome, krange: list[int] | None) -> list[int]:
    """Replace negative krange bounds with the pangenome's number of partitions.

    :param pangenome: Pangenome, possibly already partitioned.
    :param krange: Range of k values, None for [3, k of the pangenome].
    :return: The range, [3, 20] if a bound is negative and the pangenome has no k.
    """
    krange = [3, -1] if krange is None else list(krange)
    try:
        for i in (0, 1):
            if krange[i] < 0:
                krange[i] = pangenome.parameters["partition"][
                    "# final nb of partitions"
                ]
    except KeyError:
        krange = [3, 20]
    return krange


def _draw_samples(
    organisms: list,
    min_sampling: int,
    max_sampling: int,
    depth: int,
    rng: random.Random,
) -> list[set]:
    """Draw `depth` random samples of genomes at each sample size.

    :param organisms: Genomes to draw from.
    :param min_sampling: Smallest sample size, minus one.
    :param max_sampling: Largest sample size.
    :param depth: Number of samples at each size.
    :param rng: Random generator.
    :return: The samples.
    """
    return [
        set(rng.sample(organisms, i + 1))
        for i in range(min_sampling, max_sampling)
        for _ in range(depth)
    ]


def _exact_and_soft_counts(
    pangenome: Pangenome, samples: list[set], soft_core: float, disable_bar: bool
) -> list[Counter]:
    """Count the exact and soft core and accessory families of each sample.

    :param pangenome: Pangenome.
    :param samples: Samples of genomes.
    :param soft_core: Soft core threshold.
    :param disable_bar: Disable the progress bar.
    :return: The counts of each sample.
    """
    logger.info("Computing bitarrays for each family...")
    org_bit = pangenome.compute_family_bitarrays()
    logger.info(
        "Done computing bitarrays. Comparing them to get exact and soft core stats "
        f"for {len(samples)} samples..."
    )
    rows = []
    with tqdm(
        total=len(samples) * pangenome.number_of_gene_families,
        unit="gene family",
        disable=disable_bar,
    ) as bar:
        for sample in samples:
            sample_bitarray = gmpy2.xmpz()  # pylint: disable=no-member
            for org in sample:
                sample_bitarray[org_bit[org]] = 1

            counts = Counter(
                soft_core=0, exact_core=0, exact_accessory=0, soft_accessory=0
            )
            counts["nborgs"] = len(sample)
            for fam in pangenome.gene_families:
                nb_common = gmpy2.popcount(
                    fam.bitarray & sample_bitarray
                )  # pylint: disable=no-member
                if nb_common != 0:  # otherwise the family is not in the sample
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
    return rows


def _partition_samples(
    nem_index: NemIndex,
    samples: list[set],
    k: int,
    options: RarefactionOptions,
    rng: random.Random,
    rows: list[dict],
    cpu: int,
    disable_bar: bool,
):
    """Partition every sample and add its partition counts to its row.

    :param nem_index: NEM index of all the genomes.
    :param samples: Samples of genomes.
    :param k: Number of partitions, under 3 to estimate it on each sample.
    :param options: Rarefaction parameters, with the krange already resolved.
    :param rng: Random generator, used to shuffle the processing order.
    :param rows: Family counts of each sample, updated in place.
    :param cpu: Number of available CPUs.
    :param disable_bar: Disable the progress bar.
    """
    with get_context("fork").Pool(
        cpu, initializer=_keep_index, initargs=(nem_index,)
    ) as pool:
        logger.info(" Partitioning all samples...")
        with tqdm(
            total=len(samples), unit="samples partitioned", disable=disable_bar
        ) as bar:
            order = list(range(len(samples)))
            # shuffled so that the progress bar is closer to reality
            rng.shuffle(order)
            args = (
                (
                    sample_id,
                    [org.name for org in samples[sample_id]],
                    k,
                    options,
                )
                for sample_id in order
            )
            for counts, sample_id in pool.imap_unordered(
                _partition_sample_worker, args
            ):
                rows[sample_id] = {**counts, **rows[sample_id]}
                bar.update()


def make_rarefaction_curve(
    pangenome: Pangenome,
    options: RarefactionOptions,
    *,
    output: Path,
    cpu: int,
    disable_bar: bool,
):
    """Partition random samples of genomes of growing size and draw the rarefaction curves.

    :param pangenome: Pangenome with gene families and neighbors graph.
    :param options: Rarefaction parameters.
    :param output: Output directory.
    :param cpu: Number of available CPUs.
    :param disable_bar: Disable the progress bars.
    """
    krange = _resolve_krange(pangenome, options.krange)
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=disable_bar,
    )
    nem_index = build_nem_index(pangenome)

    if float(pangenome.number_of_organisms) < options.max_sampling:
        max_sampling = pangenome.number_of_organisms
    else:
        max_sampling = int(options.max_sampling)

    k = options.k
    if k < 3 and not options.reestimate_k:  # estimate k once and for all
        try:
            k = pangenome.parameters["partition"]["# final nb of partitions"]
            logger.info(f"Reuse the number of partitions {k}")
        except KeyError:
            logger.info("Estimating the number of partitions...")
            k = evaluate_nb_partitions(
                index=nem_index,
                output=None,
                max_degree_smoothing=options.max_degree_smoothing,
                free_dispersion=options.free_dispersion,
                chunk_size=options.chunk_size,
                krange=krange,
                draw_icl=False,
                cpu=cpu,
                seed=options.seed,
                disable_bar=disable_bar,
            )
            logger.info(f"The number of partitions has been evaluated at {k}")

    logger.info("Extracting samples ...")
    rng = random.Random(options.seed)
    samples = _draw_samples(
        list(pangenome.organisms),
        options.min_sampling,
        max_sampling,
        options.depth,
        rng,
    )
    logger.info(f"Done sampling genomes in the pan, there are {len(samples)} samples")

    rows = _exact_and_soft_counts(pangenome, samples, options.soft_core, disable_bar)
    _partition_samples(
        nem_index,
        samples,
        k,
        replace(options, krange=krange),
        rng,
        rows,
        cpu,
        disable_bar,
    )

    logger.info("Done  partitioning everything")
    warnings.filterwarnings("ignore")
    draw_curve(output, rows, max_sampling)
    warnings.resetwarnings()
    logger.info("Done making the rarefaction curves")


def launch(args: argparse.Namespace):
    """Command launcher.

    :param args: All arguments provided by the user.
    """
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    make_rarefaction_curve(
        pangenome,
        RarefactionOptions.from_args(args),
        output=args.output,
        cpu=args.cpu,
        disable_bar=args.disable_prog_bar,
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create the subparser of the rarefaction command.

    :param sub_parser: Subparsers of the main parser.
    :return: The rarefaction parser.
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
    """Add the arguments of the rarefaction command.

    :param parser: Parser to add them to.
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
