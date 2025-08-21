#!/usr/bin/env python3

# default libraries
import logging
import random
import tempfile
import time
from multiprocessing import get_context
import os
import argparse
from collections import defaultdict, Counter
import math
from shutil import copytree
from pathlib import Path

# installed libraries
from typing import Union, Tuple, List

from tqdm import tqdm
import plotly.offline as out_plotly
import plotly.graph_objs as go

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome

# cython library (local)
import nem_stats

pan = Pangenome()
samples = []


def run_partitioning(
    nem_dir_path: Path,
    nb_org: int,
    beta: float = 2.5,
    free_dispersion: bool = False,
    kval: int = 3,
    seed: int = 42,
    init: str = "param_file",
    keep_files: bool = False,
    itermax: int = 100,
    just_log_likelihood: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """
    Main function to make partitioning

    :param nem_dir_path: Path to directory with nem files
    :param nb_org: Number of organisms
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param kval: Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically.
    :param seed: seed used to generate random numbers
    :param init: Initiate nem parameters with pangenome parameters or randomly
    :param keep_files: True if you want to keep the NEM files
    :param itermax: Maximum iteration to compute partitioning
    :param just_log_likelihood: Return only nem parameter result

    :return: Nem parameters and if not just log likelihood the families associated to partition
    """
    logging.getLogger("PPanGGOLiN").debug("run_partitioning...")
    if init == "param_file":
        with open(nem_dir_path / f"nem_file_init_{str(kval)}.m", "w") as m_file:
            m_file.write("1 ")  # 1 to initialize parameter,
            m_file.write(" ".join([str(round(1 / float(kval), 2))] * (kval - 1)) + " ")
            # 1/K give the initial proportion to each class
            # (the last proportion is automatically determined by subtraction in nem)
            mu = []
            epsilon = []
            step = 0.5 / (math.ceil(kval / 2))
            pichenette = 0.1 if kval == 2 else 0
            for k in range(1, kval + 1):
                if k <= kval / 2:
                    mu += ["1"] * nb_org
                    epsilon += [str((step * k) - pichenette)] * nb_org
                else:
                    mu += ["0"] * nb_org
                    epsilon += [str((step * (kval - k + 1)) - pichenette)] * nb_org

            m_file.write(" ".join(mu) + " " + " ".join(epsilon))

    algo = b"nem"  # fuzzy classification by mean field approximation
    model = b"bern"  # multivariate Bernoulli mixture model
    proportion = b"pk"  # equal proportion :  "p_"     varying proportion : "pk"

    variance_model = b"skd" if free_dispersion else b"sk_"
    # one variance per partition and organism : "sdk"      one variance per partition,
    # same in all organisms : "sd_"   one variance per organism,
    # same in all partition : "s_d"    same variance in organisms and partitions : "s__"

    convergence = b"clas"
    convergence_th = 0.01
    # (INIT_SORT, init_random, init_param_file, INIT_FILE, INIT_LABEL, INIT_NB) = range(0,6)
    init_random, init_param_file = range(1, 3)
    logging.getLogger("PPanGGOLiN").debug("Running NEM...")
    logging.getLogger("PPanGGOLiN").debug(
        [
            nem_dir_path.as_posix().encode("ascii") + b"/nem_file",
            kval,
            algo,
            beta,
            convergence,
            convergence_th,
            b"fuzzy",
            itermax,
            True,
            model,
            proportion,
            variance_model,
            init_param_file if init in ["param_file", "init_from_old"] else init_random,
            nem_dir_path.as_posix().encode("ascii")
            + b"/nem_file_init_"
            + str(kval).encode("ascii")
            + b".m",
            nem_dir_path.as_posix().encode("ascii")
            + b"/nem_file_"
            + str(kval).encode("ascii"),
            seed,
        ]
    )
    nem_stats.nem(
        Fname=nem_dir_path.as_posix().encode("ascii") + b"/nem_file",
        nk=kval,
        algo=algo,
        beta=beta,
        convergence=convergence,
        convergence_th=convergence_th,
        format=b"fuzzy",
        it_max=itermax,
        dolog=True,
        model_family=model,
        proportion=proportion,
        dispersion=variance_model,
        init_mode=(
            init_param_file if init in ["param_file", "init_from_old"] else init_random
        ),
        init_file=nem_dir_path.as_posix().encode("ascii")
        + b"/nem_file_init_"
        + str(kval).encode("ascii")
        + b".m",
        out_file_prefix=nem_dir_path.as_posix().encode("ascii")
        + b"/nem_file_"
        + str(kval).encode("ascii"),
        seed=seed,
    )

    logging.getLogger("PPanGGOLiN").debug("After running NEM...")

    no_nem = False
    nem_out_path = nem_dir_path / f"nem_file_{str(kval)}.uf"
    if nem_out_path.is_file():
        logging.getLogger("PPanGGOLiN").debug("Reading NEM results...")
    elif not just_log_likelihood:
        # logging.getLogger("PPanGGOLiN").warning("No NEM output file found: "+ nem_dir_path+"/nem_file_"+str(K)+".uf")
        no_nem = True
    else:
        logging.getLogger("PPanGGOLiN").debug(
            f"No NEM output file found: {nem_out_path.absolute().as_posix()}"
        )
        no_nem = True
    index_fam = []

    with open(nem_dir_path / "nem_file.index") as index_nem_file:
        for line in index_nem_file:
            index_fam.append(line.split("\t")[1].strip())

    partitions_list = ["U"] * len(index_fam)
    all_parameters = {}
    log_likelihood = None
    entropy = None
    try:
        with (
            open(nem_dir_path / f"nem_file_{str(kval)}.uf") as partitions_nem_file,
            open(nem_dir_path / f"nem_file_{str(kval)}.mf") as parameters_nem_file,
        ):
            parameters = parameters_nem_file.readlines()
            log_likelihood = float(parameters[2].split()[3])

            sum_mu_k = []
            sum_epsilon_k = []

            for k, line in enumerate(parameters[-kval:]):
                vector = line.split()
                mu_k = [bool(float(mu_kj)) for mu_kj in vector[0:nb_org]]
                epsilon_k = [float(epsilon_kj) for epsilon_kj in vector[nb_org + 1 :]]
                proportion = float(vector[nb_org])
                sum_mu_k.append(sum(mu_k))
                sum_epsilon_k.append(sum(epsilon_k))
                if k == 0:
                    all_parameters["persistent"] = (mu_k, epsilon_k, proportion)
                elif k == kval - 1:
                    all_parameters["cloud"] = (mu_k, epsilon_k, proportion)
                else:
                    all_parameters["shell_" + str(k)] = (mu_k, epsilon_k, proportion)

            parti = {0: "P", kval - 1: "C"}

            for i in range(1, kval - 1):
                parti[i] = "S" + str(i)
            entropy = 0

            for i, line in enumerate(partitions_nem_file):
                elements = [float(el) for el in line.split()]
                if just_log_likelihood:
                    entropy += sum(
                        [
                            math.log(float(el)) * float(el) if float(el) > 0 else 0
                            for el in elements
                        ]
                    )
                else:
                    max_prob = max([float(el) for el in elements])
                    positions_max_prob = [
                        pos for pos, prob in enumerate(elements) if prob == max_prob
                    ]
                    if len(positions_max_prob) > 1 or max_prob < 0.5:
                        partitions_list[i] = (
                            "S_"  # SHELL in case of doubt gene families is attributed to shell
                        )
                    else:
                        partitions_list[i] = parti[positions_max_prob.pop()]
    except OSError:
        logging.getLogger("PPanGGOLiN").warning(
            "Partitioning did not work (the number of genomes used is probably too low), "
            f"see logs here to obtain more details {nem_dir_path.as_posix()}"
        )
        return {}, None, None  # return empty objects

    except ValueError:
        # return the default partitions_list which correspond to undefined
        pass

    if just_log_likelihood:
        return kval, log_likelihood, entropy
    else:
        return dict(zip(index_fam, partitions_list)), all_parameters, log_likelihood


def nem_single(
    args: List[Tuple[Path, int, float, bool, int, int, str, bool, int, bool]]
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """
    Allow to run partitioning in multiprocessing to evaluate partition number

    :param args: {nem_dir_path: str, nb_org: int, beta: float, free_dispersion: bool, kval: int, seed: int,
                  init: str, keep_files: bool, itermax: int, just_log_likelihood: bool}
    :return: Result of run partitioning
    """
    return run_partitioning(*args)


def partition_nem(
    index: int,
    kval: int,
    beta: float = 2.5,
    sm_degree: int = 10,
    free_dispersion: bool = False,
    seed: int = 42,
    init: str = "param_file",
    tmpdir: Path = None,
    keep_tmp_files: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """

    :param index: Index of the sample group
    :param tmpdir: temporary directory path
    :param kval: Number of partitions to use
    :param beta: strength of the smoothing using the graph topology during partitioning. 0 deactivate spatial smoothing
    :param sm_degree:  Maximum degree of the nodes to be included in the smoothing process.
    :param free_dispersion: use if the dispersion around the centroid vector of each partition during must be free.
    :param seed: seed used to generate random numbers
    :param init: Initiate nem parameters with pangenome parameters or randomly
    :param keep_tmp_files: True if you want to keep the temporary NEM files

    :return:
    """
    currtmpdir = tmpdir / f"{str(index)}"  # unique directory name
    samp = samples[index]  # org_samples accessible because it is a global variable.

    edges_weight, nb_fam = write_nem_input_files(
        tmpdir=currtmpdir, organisms=samp, sm_degree=sm_degree
    )
    return run_partitioning(
        currtmpdir,
        len(samp),
        beta * (nb_fam / edges_weight),
        free_dispersion,
        kval=kval,
        seed=seed,
        init=init,
        keep_files=keep_tmp_files,
    )


def nem_samples(
    pack: tuple,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """run partitioning
    :param pack: {index: int, tmpdir: str, beta: float, sm_degree: int, free_dispersion: bool, kval: int, seed: int, init: str, keep_tmp_files: bool}

    :return:
    """
    return partition_nem(*pack)


def write_nem_input_files(
    tmpdir: Path, organisms: set, sm_degree: int = 10
) -> Tuple[float, int]:
    """
    Create and format input files for partitioning with NEM

    :param tmpdir: temporary directory path
    :param organisms: Set of organism from pangenome
    :param sm_degree: Maximum degree of the nodes to be included in the smoothing process.

    :return: total edge weight to ponderate beta and number of families
    """
    mk_outdir(tmpdir, force=False)
    total_edges_weight = 0

    with open(tmpdir / "column_org_file", "w") as org_file:
        org_file.write(" ".join([f'"{org.name}"' for org in organisms]) + "\n")

    logging.getLogger("PPanGGOLiN").debug(
        "Writing nem_file.str nem_file.index nem_file.nei and nem_file.dat files"
    )
    with (
        open(tmpdir / "nem_file.str", "w") as str_file,
        open(tmpdir / "nem_file.index", "w") as index_file,
        open(tmpdir / "nem_file.nei", "w") as nei_file,
        open(tmpdir / "nem_file.dat", "w") as dat_file,
    ):

        nei_file.write("1\n")
        index_fam = {}

        index_org = {}
        default_dat = []
        for index, org in enumerate(organisms):
            default_dat.append("0")
            index_org[org] = index

        for fam in pan.gene_families:
            fam_organisms = set(fam.organisms)
            # could use bitarrays if this part is limiting?
            if not organisms.isdisjoint(fam_organisms):
                curr_dat = list(default_dat)
                curr_orgs = fam_organisms & organisms
                for org in curr_orgs:
                    curr_dat[index_org[org]] = "1"
                dat_file.write("\t".join(curr_dat) + "\n")
                index_fam[fam] = len(index_fam) + 1
                index_file.write(f"{len(index_fam)}\t{fam.name}\n")

        for fam in index_fam.keys():
            row_fam = []
            row_dist_score = []
            neighbor_number = 0
            sum_dist_score = 0
            for edge in fam.edges:  # iter on the family's edges.
                coverage = sum(
                    [
                        len(gene_list)
                        for org, gene_list in edge.get_organisms_dict().items()
                        if org in organisms
                    ]
                )
                if coverage == 0:
                    continue  # nothing interesting to write, this edge does not exist with this subset of organisms.
                distance_score = coverage / len(organisms)
                sum_dist_score += distance_score
                row_fam.append(
                    str(index_fam[edge.target if fam == edge.source else edge.source])
                )
                row_dist_score.append(str(round(distance_score, 4)))
                neighbor_number += 1
            if neighbor_number > 0 and float(neighbor_number) < sm_degree:
                total_edges_weight += sum_dist_score
                nei_file.write(
                    "\t".join(
                        [
                            str(item)
                            for sublist in [
                                [index_fam[fam]],
                                [neighbor_number],
                                row_fam,
                                row_dist_score,
                            ]
                            for item in sublist
                        ]
                    )
                    + "\n"
                )
            else:
                nei_file.write(str(index_fam[fam]) + "\t0\n")

        str_file.write("S\t" + str(len(index_fam)) + "\t" + str(len(organisms)) + "\n")
    return total_edges_weight / 2, len(index_fam)


def evaluate_nb_partitions(
    organisms: set,
    output: Path = None,
    sm_degree: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    krange: list = None,
    icl_margin: float = 0.05,
    draw_icl: bool = False,
    cpu: int = 1,
    seed: int = 42,
    tmpdir: Path = None,
    disable_bar: bool = False,
) -> int:
    """
    Evaluate the optimal number of partition for the pangenome

    :param organisms: Set of organisms from pangenome
    :param tmpdir: temporary directory path
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
    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    newtmpdir = tmpdir / "eval_partitions"

    if len(organisms) > chunk_size:
        select_organisms = set(random.sample(list(organisms), chunk_size))
    else:
        select_organisms = set(organisms)

    _, nb_fam = write_nem_input_files(newtmpdir, select_organisms, sm_degree)
    max_icl_k = 0
    args_partitionning = []
    for k in range(krange[0] - 1, krange[1] + 1):
        args_partitionning.append(
            (
                newtmpdir,
                len(select_organisms),
                0,
                free_dispersion,
                k,
                seed,
                "param_file",
                True,
                10,
                True,
            )
        )  # follow order run_partitionning args
    all_log_likelihood = []

    if cpu > 1:
        bar = tqdm(
            range(len(args_partitionning)),
            unit="Number of partitions",
            disable=disable_bar,
        )
        with get_context("fork").Pool(processes=cpu) as p:
            for result in p.imap_unordered(nem_single, args_partitionning):
                all_log_likelihood.append(result)
                bar.update()
            p.close()
            p.join()
        bar.close()
    else:  # for the case where it is called in a daemonic subprocess with a single cpu
        for arguments in args_partitionning:
            all_log_likelihood.append(nem_single(arguments))

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
        max_icl_k = max(all_icls, key=all_icls.get)
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
    tmpdir: Path = None,
    keep_tmp_files: bool = False,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Partitioning the pangenome

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :param tmpdir: temporary directory path
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
    :param keep_tmp_files: True if you want to keep the temporary NEM files
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    kmm = [3, 20] if krange is None else krange
    global samples
    global pan

    pan = pangenome
    if draw_icl and output is None:
        raise Exception(
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
    organisms = set(pangenome.organisms)

    if keep_tmp_files:
        # Create a temporary directory without auto-cleanup
        tmp_dir = tempfile.mkdtemp(dir=tmpdir)
        tmp_path = Path(tmp_dir)
    else:
        # Create a temporary directory with auto-cleanup
        tmp_dir = tempfile.TemporaryDirectory(dir=tmpdir)
        tmp_path = Path(tmp_dir.name)

    if len(organisms) <= 10:
        logging.getLogger("PPanGGOLiN").warning(
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
        logging.getLogger("PPanGGOLiN").info(
            "Estimating the optimal number of partitions..."
        )
        kval = evaluate_nb_partitions(
            organisms,
            output,
            sm_degree,
            free_dispersion,
            chunk_size,
            kmm,
            icl_margin,
            draw_icl,
            cpu,
            seed,
            tmp_path,
            disable_bar,
        )
        logging.getLogger("PPanGGOLiN").info(
            f"The number of partitions has been evaluated at {kval}"
        )

    pangenome.parameters["partition"]["# final nb of partitions"] = kval
    pangenome.parameters["partition"]["krange"] = kmm
    init = "param_file"

    partitioning_results = {}

    families = set()
    cpt = 0
    cpt_partition = {}
    random.seed(seed)

    for fam in pangenome.gene_families:
        families.add(fam)
        if chunk_size < len(organisms):
            cpt_partition[fam.name] = {"P": 0, "S": 0, "C": 0, "U": 0}

    start_partitioning = time.time()
    logging.getLogger("PPanGGOLiN").info("Partitioning...")
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

        org_nb_sample = Counter()
        for org in organisms:
            org_nb_sample[org] = 0
        condition = len(organisms) / chunk_size
        while len(validated) < pansize:
            prev = len(samples)  # if we've been sampling already, samples is not empty.
            while not all(val >= condition for val in org_nb_sample.values()):
                # each family must be tested at least len(select_organisms)/chunk_size times.
                shuffled_orgs = list(organisms)  # copy select_organisms
                random.shuffle(shuffled_orgs)  # shuffle the copied list
                while len(shuffled_orgs) > chunk_size:
                    samples.append(set(shuffled_orgs[:chunk_size]))
                    for org in samples[-1]:
                        org_nb_sample[org] += 1
                    shuffled_orgs = shuffled_orgs[chunk_size:]
            args = []
            # tmpdir, beta, sm_degree, free_dispersion, K, seed
            for i, _ in enumerate(samples[prev:], start=prev):
                args.append(
                    (
                        i,
                        kval,
                        beta,
                        sm_degree,
                        free_dispersion,
                        seed,
                        init,
                        tmp_path,
                        keep_tmp_files,
                    )
                )

            logging.getLogger("PPanGGOLiN").info("Launching NEM")
            with get_context("fork").Pool(processes=cpu) as p:
                # launch partitioning
                bar = tqdm(
                    range(len(args)), unit=" samples partitioned", disable=disable_bar
                )
                for result in p.imap_unordered(nem_samples, args):
                    validate_family(result)
                    bar.update()

                bar.close()
                condition += (
                    1  # if len(validated) < pan_size, we will want to resample more.
                )
                logging.getLogger("PPanGGOLiN").debug(
                    f"There are {len(validated)} validated families out of {pansize} families."
                )
                p.close()
                p.join()
        for fam, data in cpt_partition.items():
            partitioning_results[fam] = max(data, key=data.get)

        # need to compute the median vectors of each partition ???
        partitioning_results = [partitioning_results, []]  # introduces a 'non feature'.

        logging.getLogger("PPanGGOLiN").info(
            f"Did {len(samples)} partitioning with chunks of size {chunk_size} among "
            f"{len(organisms)} genomes in {round(time.time() - start_partitioning, 2)} seconds."
        )
    else:
        edges_weight, nb_fam = write_nem_input_files(
            tmp_path / f"{str(cpt)}", organisms, sm_degree=sm_degree
        )
        partitioning_results = run_partitioning(
            tmp_path / f"{str(cpt)}",
            len(organisms),
            beta * (nb_fam / edges_weight),
            free_dispersion,
            kval=kval,
            seed=seed,
            init=init,
            keep_files=keep_tmp_files,
        )
        if partitioning_results == [{}, None, None]:
            raise Exception(
                "Statistical partitioning does not work on your data. "
                "This usually happens because you used very few (<15) genomes."
            )
        cpt += 1
        logging.getLogger("PPanGGOLiN").info(
            f"Partitioned {len(organisms)} genomes in "
            f"{round(time.time() - start_partitioning, 2)} seconds."
        )

    # pangenome.savePartitionParameters(K, beta, free_dispersion, sm_degree, partitioning_results[1], chunk_size)

    for fam_name, part in partitioning_results[0].items():
        pangenome.get_gene_family(fam_name).partition = part

    pangenome.status["partitioned"] = "Computed"
    if not keep_tmp_files:
        tmp_dir.cleanup()
    else:
        copytree(tmp_path, output / "NEM_files/")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    if args.draw_ICL or args.keep_tmp_files:
        mk_outdir(args.output, args.force)
    global pan
    pan.add_file(args.pangenome)
    partition(
        pan,
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
        args.tmpdir,
        args.keep_tmp_files,
        args.force,
        disable_bar=args.disable_prog_bar,
    )
    logging.getLogger("PPanGGOLiN").debug("Write partition in pangenome")
    write_pangenome(pan, pan.file, args.force, disable_bar=args.disable_prog_bar)
    logging.getLogger("PPanGGOLiN").debug("Partitioning is finished")


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser(
        "partition", formatter_class=argparse.RawTextHelpFormatter
    )
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
        "-p", "--pangenome", required=False, type=Path, help="The pangenome.h5 file"
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
        help="Output directory",
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
        "--keep_tmp_files",
        required=False,
        default=False,
        action="store_true",
        help="Use if you want to keep the temporary NEM files",
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
        help="Number of available cpus",
    )
    optional.add_argument(
        "--tmpdir",
        required=False,
        type=str,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
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
