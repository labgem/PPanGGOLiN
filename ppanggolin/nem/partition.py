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
from typing import Union, Tuple

from tqdm import tqdm
import plotly.offline as out_plotly
import plotly.graph_objs as go

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.nem.model import set_threads, solve
from ppanggolin.nem.nem_input import (
    NemIndex,
    build_nem_index,
    build_nem_input,
    subset_nem_index,
)

logger = logging.getLogger("PPanGGOLiN")


def run_partitioning(
    nem_input: tuple,
    nb_org: int,
    beta: float = 2.5,
    free_dispersion: bool = False,
    kval: int = 3,
    itermax: int = 100,
    just_log_likelihood: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """
    Run NEM on inputs built by `build_nem_input`

    :param nem_input: (presence, graph, index_fam, edges_weight, nb_fam)
    :param nb_org: Number of organisms
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param kval: Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically.
    :param itermax: Maximum iteration to compute partitioning
    :param just_log_likelihood: Return only nem parameter result

    :return: Nem parameters and if not just log likelihood the families associated to partition
    """
    logger.debug("run_partitioning...")
    presence, graph, index_fam, _, _ = nem_input
    return solve(
        presence,
        graph,
        index_fam,
        nb_org,
        beta=beta,
        free_dispersion=free_dispersion,
        kval=kval,
        itermax=itermax,
        just_log_likelihood=just_log_likelihood,
    )


def nem_single(
    args: Tuple[NemIndex, int, int, float, bool, int, int, bool]
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """
    Allow to run partitioning in multiprocessing to evaluate partition number

    :param args: {chunk: NemIndex, sm_degree: int, nb_org: int, beta: float, free_dispersion: bool,
                  kval: int, itermax: int, just_log_likelihood: bool}
    :return: Result of run partitioning
    """
    chunk, sm_degree, *nem_args = args
    return run_partitioning(build_nem_input(chunk, sm_degree), *nem_args)


def partition_nem(
    chunk: NemIndex,
    kval: int,
    beta: float = 2.5,
    sm_degree: int = 10,
    free_dispersion: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """

    :param chunk: NEM index of the sample group
    :param kval: Number of partitions to use
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param sm_degree:  Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :return:
    """
    *_, edges_weight, nb_fam = nem_input = build_nem_input(chunk, sm_degree)
    return run_partitioning(
        nem_input,
        len(chunk.org_names),
        beta * (nb_fam / edges_weight),
        free_dispersion,
        kval=kval,
    )


def nem_samples(
    pack: tuple,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """run partitioning
    :param pack: {chunk: NemIndex, kval: int, beta: float, sm_degree: int, free_dispersion: bool}

    :return:
    """
    return partition_nem(*pack)


def evaluate_nb_partitions(
    index: NemIndex,
    output: Path = None,
    sm_degree: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    krange: list = None,
    icl_margin: float = 0.05,
    draw_icl: bool = False,
    cpu: int = 1,
    seed: int = 42,
    disable_bar: bool = False,
) -> int:
    """
    Evaluate the optimal number of partition for the pangenome

    :param index: NEM index of the genomes to evaluate K on
    :param output: output directory path to draw ICL
    :param sm_degree: Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param chunk_size: Size of the chunks when performing partitioning using chunks of organisms.
    :param krange: Range of K values to test when detecting K automatically.
    :param icl_margin: margin use to select the lowest K in maximizing ICL
    :param draw_icl: draw the ICL curve for all the tested K values.
    :param cpu: Number of available core
    :param seed: seed used to generate random numbers
    :param disable_bar: Disable progress bar

    :return: Ideal number of partition computed
    """

    random.seed(seed)  # Seeded before forking
    select_organisms = index.org_names
    if len(select_organisms) > chunk_size:
        select_organisms = random.sample(select_organisms, chunk_size)
        index = subset_nem_index(index, select_organisms)

    nb_fam = len(index.fam_names)
    max_icl_k = 0
    nem_args = [
        (len(select_organisms), 0, free_dispersion, k, 10, True)
        for k in range(krange[0] - 1, krange[1] + 1)
    ]
    all_log_likelihood = []

    if cpu > 1:
        bar = tqdm(
            range(len(nem_args)),
            unit="Number of partitions",
            disable=disable_bar,
        )
        with get_context("fork").Pool(
            processes=cpu, initializer=set_threads, initargs=(1,)
        ) as p:
            for result in p.imap_unordered(
                nem_single, ((index, sm_degree, *arguments) for arguments in nem_args)
            ):
                all_log_likelihood.append(result)
                bar.update()
            p.close()
            p.join()
        bar.close()
    else:  # for the case where it is called in a daemonic subprocess with a single cpu
        nem_input = build_nem_input(index, sm_degree)
        for arguments in nem_args:
            all_log_likelihood.append(run_partitioning(nem_input, *arguments))

    all_bics = defaultdict(float)
    all_icls = defaultdict(float)
    all_lls = defaultdict(float)
    for k_candidate, log_likelihood, entropy in all_log_likelihood:
        if log_likelihood is not None:
            nb_params = k_candidate * (
                len(select_organisms)
                + 1
                + (len(select_organisms) if free_dispersion else 1)
            )
            all_bics[k_candidate] = log_likelihood - 0.5 * (
                math.log(nb_params) * nb_fam
            )  # Calculate BIC
            all_icls[k_candidate] = all_bics[k_candidate] - entropy
            all_lls[k_candidate] = log_likelihood

    chosen_k = 3
    best_k = chosen_k

    if len(all_bics) > 3:
        max_icl_k = max(sorted(all_icls), key=all_icls.get)
        delta_icl = (all_icls[max_icl_k] - min(all_icls.values())) * icl_margin
        best_k = min(
            {
                k
                for k, icl in all_icls.items()
                if icl >= all_icls[max_icl_k] - delta_icl and k <= max_icl_k
            }
        )
        chosen_k = best_k if best_k >= 3 else chosen_k

    if len(all_bics) > 0 and draw_icl:
        traces = [
            go.Scatter(
                x=sorted(all_bics.keys()),
                y=[all_bics[key] for key in sorted(all_bics.keys())],
                name="BIC",
                mode="lines+markers",
            ),
            go.Scatter(
                x=sorted(all_icls.keys()),
                y=[all_icls[key] for key in sorted(all_icls.keys())],
                name="ICL",
                mode="lines+markers",
            ),
            go.Scatter(
                x=sorted(all_lls.keys()),
                y=[all_lls[key] for key in sorted(all_lls.keys())],
                name="log likelihood",
                mode="lines+markers",
            ),
        ]
        layout = go.Layout(
            title="ICL curve (best K is "
            + str(best_k)
            + ", ICL_th= is "
            + str(icl_margin)
            + ")",
            titlefont=dict(size=20),
            xaxis=dict(title="number of partitions"),
            yaxis=dict(title="ICL, BIC, log likelihood"),
            plot_bgcolor="#ffffff",
            shapes=[
                dict(
                    type="line",
                    x0=best_k,
                    x1=best_k,
                    y0=0,
                    y1=all_icls[best_k],
                    line=dict(dict(width=1, dash="dashdot", color="black")),
                ),
                dict(
                    type="line",
                    x0=max_icl_k,
                    x1=max_icl_k,
                    y0=0,
                    y1=all_icls[max_icl_k],
                    line=dict(dict(width=1, dash="dashdot", color="black")),
                ),
                dict(
                    type="line",
                    x0=best_k,
                    x1=max_icl_k,
                    y0=all_icls[max_icl_k],
                    y1=all_icls[max_icl_k],
                    line=dict(dict(width=1, dash="dashdot", color="black")),
                ),
                dict(
                    type="line",
                    x0=2,
                    x1=krange[1],
                    y0=all_icls[best_k],
                    y1=all_icls[best_k],
                    line=dict(dict(width=1, dash="dashdot", color="black")),
                ),
            ],
        )
        fig = go.Figure(data=traces, layout=layout)
        out_plot = output / f"ICL_curve_K{str(best_k)}.html"
        out_plotly.plot(fig, filename=out_plot.as_posix(), auto_open=False)
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
    sm_degree: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    kval: int = -1,
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
    :param sm_degree: Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param chunk_size: Size of the chunks when performing partitioning using chunks of organisms.
    :param kval: Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically.
    :param krange: Range of K values to test when detecting K automatically.
    :param icl_margin: margin use to select the lowest K in maximizing ICL
    :param draw_icl: draw the ICL curve for all the tested K values.
    :param cpu: Number of available core
    :param seed: seed used to generate random numbers
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    kmm = [3, 20] if krange is None else krange
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
    organisms = index.org_names

    if len(organisms) <= 10:
        logger.warning(
            f"The number of selected genomes is too low ({len(organisms)} "
            f"genomes used) to robustly partition the graph"
        )

    pangenome.parameters["partition"] = {}
    pangenome.parameters["partition"]["beta"] = beta
    pangenome.parameters["partition"]["max_degree_smoothing"] = sm_degree
    pangenome.parameters["partition"]["free_dispersion"] = free_dispersion
    pangenome.parameters["partition"]["ICL_margin"] = icl_margin
    pangenome.parameters["partition"]["seed"] = seed
    if len(organisms) > chunk_size:
        pangenome.parameters["partition"]["chunk_size"] = chunk_size
    pangenome.parameters["partition"]["# computed nb of partitions"] = False

    # the K value initially given by the user
    pangenome.parameters["partition"]["nb_of_partitions"] = kval
    if kval < 2:
        pangenome.parameters["partition"]["# computed nb of partitions"] = True
        logger.info("Estimating the optimal number of partitions...")
        kval = evaluate_nb_partitions(
            index,
            output,
            sm_degree,
            free_dispersion,
            chunk_size,
            kmm,
            icl_margin,
            draw_icl,
            cpu,
            seed,
            disable_bar,
        )
        logger.info(f"The number of partitions has been evaluated at {kval}")

    pangenome.parameters["partition"]["# final nb of partitions"] = kval
    pangenome.parameters["partition"]["krange"] = kmm

    partitioning_results = {}

    families = set()
    cpt_partition = {}
    random.seed(seed)

    for fam in pangenome.gene_families:
        families.add(fam)
        if chunk_size < len(organisms):
            cpt_partition[fam.name] = {"P": 0, "S": 0, "C": 0, "U": 0}

    start_partitioning = time.time()
    logger.info("Partitioning...")
    pansize = len(families)
    if chunk_size < len(organisms):
        validated = set()

        def validate_family(res):
            """
            Validate partition assignation to families

            :param res: Partitioning results
            """
            for node, nem_class in res[0].items():
                cpt_partition[node][nem_class[0]] += 1
                sum_partionning = sum(cpt_partition[node].values())
                if (
                    sum_partionning > len(organisms) / chunk_size
                    and max(cpt_partition[node].values()) >= sum_partionning * 0.5
                ) or (sum_partionning > len(organisms)):
                    if node not in validated:
                        if max(cpt_partition[node].values()) < sum_partionning * 0.5:
                            cpt_partition[node]["U"] = len(organisms)
                            # if despite len(select_organisms) partionning,
                            # an abosolute majority is not found then the families is set to undefined
                        validated.add(node)

        samples = []
        org_nb_sample = Counter()
        for org in organisms:
            org_nb_sample[org] = 0
        condition = len(organisms) / chunk_size
        while len(validated) < pansize:
            prev = len(samples)  # if we've been sampling already, samples is not empty.
            while not all(val >= condition for val in org_nb_sample.values()):
                # each family must be tested at least len(select_organisms)/chunk_size times.
                shuffled_orgs = list(organisms)
                random.shuffle(shuffled_orgs)  # seeded shuffle
                while len(shuffled_orgs) > chunk_size:
                    samples.append(shuffled_orgs[:chunk_size])
                    for org in samples[-1]:
                        org_nb_sample[org] += 1
                    shuffled_orgs = shuffled_orgs[chunk_size:]
            args = (
                (
                    subset_nem_index(index, samp),
                    kval,
                    beta,
                    sm_degree,
                    free_dispersion,
                )
                for samp in samples[prev:]
            )

            logger.info("Launching NEM")
            with get_context("fork").Pool(
                processes=cpu, initializer=set_threads, initargs=(1,)
            ) as p:
                # launch partitioning
                bar = tqdm(
                    range(len(samples) - prev),
                    unit=" samples partitioned",
                    disable=disable_bar,
                )
                for result in p.imap_unordered(nem_samples, args):
                    validate_family(result)
                    bar.update()

                bar.close()
                condition += (
                    1  # if len(validated) < pan_size, we will want to resample more.
                )
                logger.debug(
                    f"There are {len(validated)} validated families out of {pansize} families."
                )
                p.close()
                p.join()
        for fam, data in cpt_partition.items():
            partitioning_results[fam] = max(data, key=data.get)

        # need to compute the median vectors of each partition ???
        partitioning_results = [partitioning_results, []]  # introduces a 'non feature'.

        logger.info(
            f"Did {len(samples)} partitioning with chunks of size {chunk_size} among "
            f"{len(organisms)} genomes in {round(time.time() - start_partitioning, 2)} seconds."
        )
    else:
        *_, edges_weight, nb_fam = nem_input = build_nem_input(index, sm_degree)
        partitioning_results = run_partitioning(
            nem_input,
            len(organisms),
            beta * (nb_fam / edges_weight),
            free_dispersion,
            kval=kval,
        )
        if not partitioning_results[0]:
            raise Exception(
                "Statistical partitioning does not work on your data. "
                "This usually happens because you used very few (<15) genomes."
            )
        logger.info(
            f"Partitioned {len(organisms)} genomes in "
            f"{round(time.time() - start_partitioning, 2)} seconds."
        )

    # pangenome.savePartitionParameters(K, beta, free_dispersion, sm_degree, partitioning_results[1], chunk_size)

    for fam_name, part in partitioning_results[0].items():
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
