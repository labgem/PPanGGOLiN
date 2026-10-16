#!/usr/bin/env python3

# default libraries
import logging
import random
import time
from multiprocessing import get_context
import os
import argparse
from collections import defaultdict, Counter
import math
from pathlib import Path

# installed libraries
from tqdm import tqdm
import plotly.offline as out_plotly
import plotly.graph_objs as go

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.nem.model import KScore, NemResult, score, solve
from ppanggolin.nem.nem_input import (
    NemIndex,
    build_nem_index,
    build_nem_input,
    subset_nem_index,
)

logger = logging.getLogger("PPanGGOLiN")


class ChunkVotes:
    """Labels each gene family got over the chunks of genomes it was partitioned in."""

    def __init__(self, families: list[str], nb_genomes: int, chunk_size: int):
        """Start with no votes.

        :param families: Names of the gene families to vote on.
        :param nb_genomes: Number of genomes the chunks are drawn from.
        :param chunk_size: Number of genomes in a chunk.
        """
        self.votes = {fam: {"P": 0, "S": 0, "C": 0, "U": 0} for fam in families}
        self.validated = set()
        self.nb_genomes = nb_genomes
        self.chunk_size = chunk_size

    @property
    def done(self) -> bool:
        """Whether every family has a validated label."""
        return len(self.validated) >= len(self.votes)

    def add(self, partitions: dict[str, str]):
        """Count the labels of one chunk and validate the families with a clear majority.

        A family left without an absolute majority after more votes than genomes is undefined.

        :param partitions: Label of each family of the chunk.
        """
        for family, label in partitions.items():
            votes = self.votes[family]
            votes[label[0]] += 1
            nb_votes = sum(votes.values())
            majority = max(votes.values()) >= nb_votes * 0.5
            if (
                nb_votes > self.nb_genomes / self.chunk_size and majority
            ) or nb_votes > self.nb_genomes:
                if family not in self.validated:
                    if not majority:
                        votes["U"] = self.nb_genomes
                    self.validated.add(family)

    def labels(self) -> dict[str, str]:
        """Give each family its most voted label.

        :return: Label of each family.
        """
        return {
            family: max(votes, key=votes.get) for family, votes in self.votes.items()
        }


def draw_chunks(
    genomes: list[str],
    chunk_size: int,
    times_drawn: Counter,
    min_draws: float,
    rng: random.Random,
) -> list[list[str]]:
    """Shuffle the genomes into chunks until each was drawn at least `min_draws` times.

    :param genomes: Genome names, in a stable order.
    :param chunk_size: Number of genomes in a chunk.
    :param times_drawn: Number of chunks each genome is already in, updated in place.
    :param min_draws: Number of chunks each genome must be in.
    :param rng: Random generator used to shuffle.
    :return: The new chunks.
    """
    chunks = []
    while not all(n >= min_draws for n in times_drawn.values()):
        shuffled = list(genomes)
        rng.shuffle(shuffled)
        while len(shuffled) > chunk_size:
            chunks.append(shuffled[:chunk_size])
            for genome in chunks[-1]:
                times_drawn[genome] += 1
            shuffled = shuffled[chunk_size:]
    return chunks


def partition_chunk(
    chunk: NemIndex,
    k: int,
    beta: float = 2.5,
    max_degree_smoothing: int = 10,
    free_dispersion: bool = False,
) -> NemResult:
    """

    :param chunk: NEM index of the sample group
    :param k: Number of partitions to use
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param max_degree_smoothing:  Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :return:
    """
    nem_input = build_nem_input(chunk, max_degree_smoothing)
    return solve(nem_input, nem_input.scaled_beta(beta), free_dispersion, k)


def _partition_chunk_worker(pack: tuple) -> NemResult:
    """run partitioning
    :param pack: {chunk: NemIndex, k: int, beta: float, max_degree_smoothing: int, free_dispersion: bool}

    :return:
    """
    return partition_chunk(*pack)


def _score_k_worker(args: tuple[NemIndex, int, bool, int]) -> KScore:
    """Score one k in a pool.

    :param args: (index, max_degree_smoothing, free_dispersion, k).
    :return: The score of k.
    """
    index, max_degree_smoothing, free_dispersion, k = args
    nem_input = build_nem_input(index, max_degree_smoothing)
    return score(nem_input, 0, free_dispersion, k, max_iter=10)


def _score_ks(
    index: NemIndex,
    ks: range,
    max_degree_smoothing: int,
    free_dispersion: bool,
    cpu: int,
    disable_bar: bool,
) -> list[KScore]:
    """Fit NEM without smoothing for each k.

    :param index: Index of the genomes to score k on.
    :param ks: Values of k to score.
    :param max_degree_smoothing: Maximum degree of a family included in the smoothing.
    :param free_dispersion: Free dispersion per genome instead of shared.
    :param cpu: Number of available CPUs.
    :param disable_bar: Disable the progress bar.
    :return: One score per k.
    """
    if cpu <= 1:  # also the case inside rarefaction's daemonic workers
        nem_input = build_nem_input(index, max_degree_smoothing)
        return [score(nem_input, 0, free_dispersion, k, max_iter=10) for k in ks]

    k_scores = []
    args = ((index, max_degree_smoothing, free_dispersion, k) for k in ks)
    with get_context("fork").Pool(processes=cpu) as pool:
        with tqdm(
            total=len(ks), unit="Number of partitions", disable=disable_bar
        ) as bar:
            for k_score in pool.imap_unordered(_score_k_worker, args):
                k_scores.append(k_score)
                bar.update()
    return k_scores


def _information_criteria(
    k_scores: list[KScore], nb_genomes: int, nb_families: int, free_dispersion: bool
) -> tuple[defaultdict, defaultdict, defaultdict]:
    """Compute the BIC and ICL of each k.

    :param k_scores: Score of each k.
    :param nb_genomes: Number of genomes the scores were computed on.
    :param nb_families: Number of families the scores were computed on.
    :param free_dispersion: Free dispersion per genome instead of shared.
    :return: BIC, ICL and log-likelihood of each k.
    """
    bics, icls, lls = defaultdict(float), defaultdict(float), defaultdict(float)
    for k_score in k_scores:
        k = k_score.k
        nb_params = k * (nb_genomes + 1 + (nb_genomes if free_dispersion else 1))
        bics[k] = k_score.log_likelihood - 0.5 * (math.log(nb_params) * nb_families)
        icls[k] = bics[k] - k_score.entropy
        lls[k] = k_score.log_likelihood
    return bics, icls, lls


def _select_k(icls: defaultdict, icl_margin: float) -> tuple[int, int, int]:
    """Choose the lowest k whose ICL is within a margin of the best one.

    :param icls: ICL of each k.
    :param icl_margin: Margin, as a fraction of the ICL range.
    :return: The chosen k (at least 3), the selected k and the k of maximal ICL.
    """
    chosen_k = best_k = 3
    max_icl_k = 0
    if len(icls) > 3:
        max_icl_k = max(sorted(icls), key=icls.get)
        delta_icl = (icls[max_icl_k] - min(icls.values())) * icl_margin
        best_k = min(
            k
            for k, icl in icls.items()
            if icl >= icls[max_icl_k] - delta_icl and k <= max_icl_k
        )
        chosen_k = best_k if best_k >= 3 else chosen_k
    return chosen_k, best_k, max_icl_k


def _draw_icl_curve(
    output: Path,
    bics: defaultdict,
    icls: defaultdict,
    lls: defaultdict,
    best_k: int,
    max_icl_k: int,
    icl_margin: float,
    max_k: int,
):
    """Plot the BIC, ICL and log-likelihood of each k.

    :param output: Output directory.
    :param bics: BIC of each k.
    :param icls: ICL of each k.
    :param lls: Log-likelihood of each k.
    :param best_k: Selected k.
    :param max_icl_k: k of maximal ICL.
    :param icl_margin: Margin used to select k.
    :param max_k: Largest k tested.
    """
    traces = [
        go.Scatter(
            x=sorted(values.keys()),
            y=[values[key] for key in sorted(values.keys())],
            name=name,
            mode="lines+markers",
        )
        for name, values in (("BIC", bics), ("ICL", icls), ("log likelihood", lls))
    ]
    line = dict(width=1, dash="dashdot", color="black")
    layout = go.Layout(
        title=f"ICL curve (best K is {best_k}, ICL_th= is {icl_margin})",
        titlefont=dict(size=20),
        xaxis=dict(title="number of partitions"),
        yaxis=dict(title="ICL, BIC, log likelihood"),
        plot_bgcolor="#ffffff",
        shapes=[
            dict(type="line", x0=best_k, x1=best_k, y0=0, y1=icls[best_k], line=line),
            dict(
                type="line",
                x0=max_icl_k,
                x1=max_icl_k,
                y0=0,
                y1=icls[max_icl_k],
                line=line,
            ),
            dict(
                type="line",
                x0=best_k,
                x1=max_icl_k,
                y0=icls[max_icl_k],
                y1=icls[max_icl_k],
                line=line,
            ),
            dict(
                type="line",
                x0=2,
                x1=max_k,
                y0=icls[best_k],
                y1=icls[best_k],
                line=line,
            ),
        ],
    )
    fig = go.Figure(data=traces, layout=layout)
    out_plot = output / f"ICL_curve_K{best_k}.html"
    out_plotly.plot(fig, filename=out_plot.as_posix(), auto_open=False)


def evaluate_nb_partitions(
    index: NemIndex,
    output: Path | None = None,
    max_degree_smoothing: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    krange: list[int] | None = None,
    icl_margin: float = 0.05,
    draw_icl: bool = False,
    cpu: int = 1,
    seed: int = 42,
    disable_bar: bool = False,
) -> int:
    """Estimate the number of partitions of the pangenome.

    :param index: NEM index of the genomes to estimate k on.
    :param output: Output directory, needed to draw the ICL curve.
    :param max_degree_smoothing: Maximum degree of a family included in the smoothing.
    :param free_dispersion: Free dispersion per genome instead of shared.
    :param chunk_size: Estimate k on a random chunk of this many genomes when there are more.
    :param krange: Range of k values to test.
    :param icl_margin: Margin used to select the lowest k maximizing the ICL.
    :param draw_icl: Draw the ICL curve of the tested k values.
    :param cpu: Number of available CPUs.
    :param seed: Seed used to generate random numbers.
    :param disable_bar: Disable the progress bar.
    :return: The estimated number of partitions.
    """
    krange = [3, 20] if krange is None else krange
    genomes = index.org_names
    if len(genomes) > chunk_size:
        genomes = random.Random(seed).sample(genomes, chunk_size)
        index = subset_nem_index(index, genomes)

    ks = range(krange[0] - 1, krange[1] + 1)
    k_scores = _score_ks(
        index, ks, max_degree_smoothing, free_dispersion, cpu, disable_bar
    )
    bics, icls, lls = _information_criteria(
        k_scores, len(genomes), len(index.fam_names), free_dispersion
    )
    chosen_k, best_k, max_icl_k = _select_k(icls, icl_margin)

    if len(bics) > 0 and draw_icl:
        _draw_icl_curve(
            output, bics, icls, lls, best_k, max_icl_k, icl_margin, krange[1]
        )
    return chosen_k


def check_pangenome_former_partition(pangenome: Pangenome, force: bool = False):
    """checks pangenome status and .h5 files for former partitions, delete them if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    """
    if pangenome.status["partitioned"] == "inFile" and not force:
        raise Exception(
            "You are trying to partition a pangenome already partitioned."
            " If you REALLY want to do that, "
            "use --force (it will erase partitions and every feature computed from them."
        )
    elif pangenome.status["partitioned"] == "inFile" and force:
        erase_pangenome(pangenome, partition=True)


def partition(
    pangenome: Pangenome,
    output: Path = None,
    beta: float = 2.5,
    max_degree_smoothing: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    k: int = -1,
    krange: list = None,
    icl_margin: float = 0.05,
    draw_icl: bool = False,
    cpu: int = 1,
    seed: int = 42,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Partitioning the pangenome

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :param output: output directory path to draw ICL
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param max_degree_smoothing: Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param chunk_size: Size of the chunks when performing partitioning using chunks of genomes.
    :param k: Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically.
    :param krange: Range of K values to test when detecting K automatically.
    :param icl_margin: margin use to select the lowest K in maximizing ICL
    :param draw_icl: draw the ICL curve for all the tested K values.
    :param cpu: Number of available core
    :param seed: seed used to generate random numbers
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    krange = [3, 20] if krange is None else krange
    if draw_icl and output is None:
        raise ValueError(
            "Combination of option impossible: "
            "You asked to draw the ICL curves but did not provide an output directory!"
        )
    check_pangenome_former_partition(pangenome, force)
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=disable_bar,
    )
    index = build_nem_index(pangenome)
    genomes = index.org_names

    if len(genomes) <= 10:
        logger.warning(
            f"The number of selected genomes is too low ({len(genomes)} "
            f"genomes used) to robustly partition the graph"
        )

    pangenome.parameters["partition"] = {}
    pangenome.parameters["partition"]["beta"] = beta
    pangenome.parameters["partition"]["max_degree_smoothing"] = max_degree_smoothing
    pangenome.parameters["partition"]["free_dispersion"] = free_dispersion
    pangenome.parameters["partition"]["ICL_margin"] = icl_margin
    pangenome.parameters["partition"]["seed"] = seed
    if len(genomes) > chunk_size:
        pangenome.parameters["partition"]["chunk_size"] = chunk_size
    pangenome.parameters["partition"]["# computed nb of partitions"] = False

    # the K value initially given by the user
    pangenome.parameters["partition"]["nb_of_partitions"] = k
    if k < 2:
        pangenome.parameters["partition"]["# computed nb of partitions"] = True
        logger.info("Estimating the optimal number of partitions...")
        k = evaluate_nb_partitions(
            index,
            output,
            max_degree_smoothing,
            free_dispersion,
            chunk_size,
            krange,
            icl_margin,
            draw_icl,
            cpu,
            seed,
            disable_bar,
        )
        logger.info(f"The number of partitions has been evaluated at {k}")

    pangenome.parameters["partition"]["# final nb of partitions"] = k
    pangenome.parameters["partition"]["krange"] = krange

    start_partitioning = time.time()
    logger.info("Partitioning...")
    if chunk_size < len(genomes):
        families = [fam.name for fam in pangenome.gene_families]
        rng = random.Random(seed)
        votes = ChunkVotes(families, len(genomes), chunk_size)
        times_drawn = Counter({genome: 0 for genome in genomes})
        min_draws = len(genomes) / chunk_size
        nb_chunks = 0
        while not votes.done:
            chunks = draw_chunks(genomes, chunk_size, times_drawn, min_draws, rng)
            nb_chunks += len(chunks)
            # sliced as the pool asks for them, so one chunk is in transit at a time
            args = (
                (
                    subset_nem_index(index, chunk),
                    k,
                    beta,
                    max_degree_smoothing,
                    free_dispersion,
                )
                for chunk in chunks
            )

            logger.info("Launching NEM")
            with get_context("fork").Pool(processes=cpu) as p:
                # launch partitioning
                bar = tqdm(
                    range(len(chunks)),
                    unit=" samples partitioned",
                    disable=disable_bar,
                )
                for result in p.imap_unordered(_partition_chunk_worker, args):
                    votes.add(result.partitions)
                    bar.update()

                bar.close()
                min_draws += 1  # if families are left to validate, draw more chunks
                logger.debug(
                    f"There are {len(votes.validated)} validated families out of {len(families)} families."
                )
                p.close()
                p.join()
        labels = votes.labels()

        logger.info(
            f"Did {nb_chunks} partitioning with chunks of size {chunk_size} among "
            f"{len(genomes)} genomes in {round(time.time() - start_partitioning, 2)} seconds."
        )
    else:
        nem_input = build_nem_input(index, max_degree_smoothing)
        labels = solve(
            nem_input, nem_input.scaled_beta(beta), free_dispersion, k
        ).partitions
        if not labels:
            raise Exception(
                "Statistical partitioning does not work on your data. "
                "This usually happens because you used very few (<15) genomes."
            )
        logger.info(
            f"Partitioned {len(genomes)} genomes in "
            f"{round(time.time() - start_partitioning, 2)} seconds."
        )

    for fam_name, part in labels.items():
        pangenome.get_gene_family(fam_name).partition = part

    pangenome.status["partitioned"] = "Computed"


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if args.draw_ICL:
        mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    partition(
        pangenome,
        args.output,
        args.beta,
        args.max_degree_smoothing,
        args.free_dispersion,
        args.chunk_size,
        args.nb_of_partitions,
        args.krange,
        args.ICL_margin,
        args.draw_ICL,
        args.cpu,
        args.seed,
        args.force,
        disable_bar=args.disable_prog_bar,
    )
    logger.debug("Write partition in pangenome")
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )
    logger.debug("Partitioning is finished")


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "partition", formatter_class=argparse.RawTextHelpFormatter
    )
    parser.description = "Partition the pangenome graph"
    parser.category = "Expert"
    parser_partition(parser)
    return parser


def parser_partition(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of partition command

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
        help="beta is the strength of the smoothing using the graph topology during partitioning. "
        "0 will deactivate spatial smoothing.",
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
        "If under 2, it will be detected automatically.",
    )
    optional.add_argument(
        "-Kmm",
        "--krange",
        nargs=2,
        required=False,
        type=int,
        default=[3, 20],
        help="Range of K values to test when detecting K automatically.",
    )
    optional.add_argument(
        "-im",
        "--ICL_margin",
        required=False,
        type=float,
        default=0.05,
        help="K is detected automatically by maximizing ICL. However at some point the ICL "
        "reaches a plateau. Therefore we are looking for the minimal value of K without "
        "significant gain from the larger values of K measured by ICL. For that we take the "
        "lowest K that is found within a given 'margin' of the maximal ICL value. Basically, "
        "change this option only if you truly understand it, otherwise just leave it be.",
    )
    optional.add_argument(
        "--draw_ICL",
        required=False,
        default=False,
        action="store_true",
        help="Use if you want to draw the ICL curve for all the tested K values. "
        "Will not be done if K is given.",
    )
    optional.add_argument(
        "-se",
        "--seed",
        type=int,
        default=42,
        help="Seed used to generate random numbers.",
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

    parser_partition(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
