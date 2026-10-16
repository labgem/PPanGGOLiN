#!/usr/bin/env python3

# default libraries
import argparse
import logging
import math
import os
import random
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from multiprocessing import get_context
from collections import defaultdict, Counter
from pathlib import Path

# installed libraries
import plotly.graph_objs as go
import plotly.offline as out_plotly
from tqdm import tqdm

# local libraries
from ppanggolin.formats import check_pangenome_info, erase_pangenome, write_pangenome
from ppanggolin.nem.model import KScore, NemResult, score, solve
from ppanggolin.nem.nem_input import (
    NemIndex,
    build_nem_index,
    build_nem_input,
    subset_nem_index,
)
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir

logger = logging.getLogger("PPanGGOLiN")


@dataclass
class PartitionOptions:
    """Parameters of the partition command.

    :param output: Output directory, needed to draw the ICL curve.
    :param beta: Strength of the smoothing using the graph topology. 0 deactivates it.
    :param max_degree_smoothing: Maximum degree of a family included in the smoothing.
    :param free_dispersion: Free dispersion around each partition's centroid per genome.
    :param chunk_size: Partition in chunks of this many genomes when there are more.
    :param k: Number of partitions. Under 2, it is estimated.
    :param krange: Range of k values tested when estimating k.
    :param icl_margin: Margin used to select the lowest k maximizing the ICL.
    :param draw_icl: Draw the ICL curve of the tested k values.
    :param cpu: Number of available CPUs.
    :param seed: Seed used to generate random numbers.
    :param force: Overwrite a former partition.
    :param disable_bar: Disable the progress bars.
    """

    output: Path | None = None
    beta: float = 2.5
    max_degree_smoothing: int = 10
    free_dispersion: bool = False
    chunk_size: int = 500
    k: int = -1
    krange: list[int] = field(default_factory=lambda: [3, 20])
    icl_margin: float = 0.05
    draw_icl: bool = False
    cpu: int = 1
    seed: int = 42
    force: bool = False
    disable_bar: bool = False


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
    """Partition the gene families of one chunk of genomes.

    :param chunk: NEM index of the chunk.
    :param k: Number of partitions.
    :param beta: Strength of smoothing using graph topology during partitioning. 0 deactivates spatial smoothing.
    :param max_degree_smoothing: Maximum degree of a family included in the smoothing process.
    :param free_dispersion: Free dispersion per genome rather than shared.
    :return: The partitions of the chunk's families.
    """
    nem_input = build_nem_input(chunk, max_degree_smoothing)
    return solve(nem_input, nem_input.scaled_beta(beta), free_dispersion, k)


def _partition_chunk_worker(pack: tuple[NemIndex, int, float, int, bool]) -> NemResult:
    """Run `partition_chunk` in a pool."""
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
    """Check the pangenome for a former partition and erase it if allowed.

    :param pangenome: Pangenome object.
    :param force: Allow to force write on Pangenome file.
    :raises Exception: If the pangenome is already partitioned and force is False.
    """
    if pangenome.status["partitioned"] == "inFile" and not force:
        raise Exception(
            "You are trying to partition a pangenome already partitioned."
            " If you REALLY want to do that, "
            "use --force (it will erase partitions and every feature computed from them."
        )
    elif pangenome.status["partitioned"] == "inFile" and force:
        erase_pangenome(pangenome, partition=True)


def _record_parameters(pangenome: Pangenome, options: PartitionOptions):
    """Store the partition parameters in the pangenome.

    :param pangenome: Pangenome being partitioned.
    :param options: Partition parameters.
    """
    parameters = {
        "beta": options.beta,
        "max_degree_smoothing": options.max_degree_smoothing,
        "free_dispersion": options.free_dispersion,
        "ICL_margin": options.icl_margin,
        "seed": options.seed,
    }
    if pangenome.number_of_organisms > options.chunk_size:
        parameters["chunk_size"] = options.chunk_size
    parameters["# computed nb of partitions"] = options.k < 2
    # the k initially given by the user
    parameters["nb_of_partitions"] = options.k
    pangenome.parameters["partition"] = parameters


def _partition_in_chunks(
    index: NemIndex, families: list[str], k: int, options: PartitionOptions
) -> tuple[dict[str, str], int]:
    """Partition chunks of genomes until every family has a validated label.

    :param index: NEM index of all the genomes.
    :param families: Names of all the gene families.
    :param k: Number of partitions.
    :param options: Partition parameters.
    :return: Label of each family, and the number of chunks partitioned.
    """
    genomes = index.org_names
    rng = random.Random(options.seed)
    votes = ChunkVotes(families, len(genomes), options.chunk_size)
    times_drawn = Counter({genome: 0 for genome in genomes})
    min_draws = len(genomes) / options.chunk_size
    nb_chunks = 0
    while not votes.done:
        chunks = draw_chunks(genomes, options.chunk_size, times_drawn, min_draws, rng)
        nb_chunks += len(chunks)
        # sliced as the pool asks for them, so one chunk is in transit at a time
        args = (
            (
                subset_nem_index(index, chunk),
                k,
                options.beta,
                options.max_degree_smoothing,
                options.free_dispersion,
            )
            for chunk in chunks
        )
        logger.info("Launching NEM")
        with get_context("fork").Pool(processes=options.cpu) as pool:
            with tqdm(
                total=len(chunks),
                unit=" samples partitioned",
                disable=options.disable_bar,
            ) as bar:
                for result in pool.imap_unordered(_partition_chunk_worker, args):
                    votes.add(result.partitions)
                    bar.update()
        min_draws += 1  # if families are left to validate, draw more chunks
        logger.debug(
            f"There are {len(votes.validated)} validated families out of {len(families)} families."
        )
    return votes.labels(), nb_chunks


def _partition_whole(
    index: NemIndex, k: int, options: PartitionOptions
) -> dict[str, str]:
    """Partition all the genomes at once.

    :param index: NEM index of all the genomes.
    :param k: Number of partitions.
    :param options: Partition parameters.
    :raises Exception: If NEM labels no family.
    :return: Label of each family.
    """
    result = partition_chunk(
        index, k, options.beta, options.max_degree_smoothing, options.free_dispersion
    )
    if not result.partitions:
        raise Exception(
            "Statistical partitioning does not work on your data. "
            "This usually happens because you used very few (<15) genomes."
        )
    return result.partitions


def partition(pangenome: Pangenome, options: PartitionOptions | None = None):
    """Partition the gene families of a pangenome.

    :param pangenome: Pangenome with gene families and neighbors graph.
    :param options: Partition parameters.
    :raises ValueError: If the ICL curve is asked for without an output directory.
    """
    options = PartitionOptions() if options is None else options
    if options.draw_icl and options.output is None:
        raise ValueError(
            "Combination of option impossible: "
            "You asked to draw the ICL curves but did not provide an output directory!"
        )
    check_pangenome_former_partition(pangenome, options.force)
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=options.disable_bar,
    )
    index = build_nem_index(pangenome)
    nb_genomes = len(index.org_names)

    if nb_genomes <= 10:
        logger.warning(
            f"The number of selected genomes is too low ({nb_genomes} "
            f"genomes used) to robustly partition the graph"
        )

    _record_parameters(pangenome, options)
    k = options.k
    if k < 2:
        logger.info("Estimating the optimal number of partitions...")
        k = evaluate_nb_partitions(
            index,
            options.output,
            options.max_degree_smoothing,
            options.free_dispersion,
            options.chunk_size,
            options.krange,
            options.icl_margin,
            options.draw_icl,
            options.cpu,
            options.seed,
            options.disable_bar,
        )
        logger.info(f"The number of partitions has been evaluated at {k}")
    pangenome.parameters["partition"]["# final nb of partitions"] = k
    pangenome.parameters["partition"]["krange"] = options.krange

    start = time.time()
    logger.info("Partitioning...")
    if options.chunk_size < nb_genomes:
        families = [fam.name for fam in pangenome.gene_families]
        labels, nb_chunks = _partition_in_chunks(index, families, k, options)
        logger.info(
            f"Did {nb_chunks} partitioning with chunks of size {options.chunk_size} among "
            f"{nb_genomes} genomes in {round(time.time() - start, 2)} seconds."
        )
    else:
        labels = _partition_whole(index, k, options)
        logger.info(
            f"Partitioned {nb_genomes} genomes in "
            f"{round(time.time() - start, 2)} seconds."
        )

    for fam_name, label in labels.items():
        pangenome.get_gene_family(fam_name).partition = label
    pangenome.status["partitioned"] = "Computed"


def launch(args: argparse.Namespace):
    """Command launcher.

    :param args: All arguments provided by the user.
    """
    if args.draw_ICL:
        mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    partition(
        pangenome,
        PartitionOptions(
            output=args.output,
            beta=args.beta,
            max_degree_smoothing=args.max_degree_smoothing,
            free_dispersion=args.free_dispersion,
            chunk_size=args.chunk_size,
            k=args.nb_of_partitions,
            krange=args.krange,
            icl_margin=args.ICL_margin,
            draw_icl=args.draw_ICL,
            cpu=args.cpu,
            seed=args.seed,
            force=args.force,
            disable_bar=args.disable_prog_bar,
        ),
    )
    logger.debug("Write partition in pangenome")
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )
    logger.debug("Partitioning is finished")


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create the subparser of the partition command.

    :param sub_parser: Subparsers of the main parser.
    :return: The partition parser.
    """
    parser = sub_parser.add_parser(
        "partition", formatter_class=argparse.RawTextHelpFormatter
    )
    parser.description = "Partition the pangenome graph"
    parser.category = "Expert"
    parser_partition(parser)
    return parser


def parser_partition(parser: argparse.ArgumentParser):
    """Add the arguments of the partition command.

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
