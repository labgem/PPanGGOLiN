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
from typing import Dict, List, Optional, Tuple, Union

from tqdm import tqdm
import plotly.offline as out_plotly
import plotly.graph_objs as go

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome

# cython library (local)
import nem_stats

# Module-level state shared with forked worker processes.
# These are written in the main process *before* any Pool is created so that
# forked child processes inherit the objects via copy-on-write without any
# serialisation (pickling a large Pangenome for every task would be slow).
_pangenome: Pangenome = Pangenome()
_samples: List[set] = []


# ─────────────────────────────────────────────────────────────────────────────
# NEM file helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_family_index(nem_dir: Path) -> List[str]:
    """Return gene-family names in NEM row order."""
    with open(nem_dir / "nem_file.index") as fh:
        return [line.split("\t")[1].strip() for line in fh]


def _partition_label_map(kval: int) -> Dict[int, str]:
    """Map 0-based cluster index to a partition label (P / S1…Sk-2 / C)."""
    labels: Dict[int, str] = {0: "P", kval - 1: "C"}
    for i in range(1, kval - 1):
        labels[i] = f"S{i}"
    return labels


# ─────────────────────────────────────────────────────────────────────────────
# C NEM backend
# ─────────────────────────────────────────────────────────────────────────────

def _write_cnem_init_file(path: Path, kval: int, nb_org: int) -> None:
    """Write the NEM parameter initialisation file (.m)."""
    step = 0.5 / math.ceil(kval / 2)
    pichenette = 0.1 if kval == 2 else 0.0
    mu: List[str] = []
    epsilon: List[str] = []
    for k in range(1, kval + 1):
        if k <= kval / 2:
            mu += ["1"] * nb_org
            epsilon += [str((step * k) - pichenette)] * nb_org
        else:
            mu += ["0"] * nb_org
            epsilon += [str((step * (kval - k + 1)) - pichenette)] * nb_org
    with open(path, "w") as fh:
        fh.write("1 ")
        # Keep enough precision so the (k-1) proportions never sum above 1.
        base_prop = format(1.0 / float(kval), ".12g")
        fh.write(" ".join([base_prop] * (kval - 1)) + " ")
        fh.write(" ".join(mu) + " " + " ".join(epsilon))


def _parse_cnem_output(
    nem_dir: Path,
    kval: int,
    nb_org: int,
    just_log_likelihood: bool,
) -> Union[Tuple[dict, dict, float], Tuple[int, float, float]]:
    """Read .uf and .mf output files produced by C NEM."""
    index_fam = _read_family_index(nem_dir)
    labels = _partition_label_map(kval)

    with (
        open(nem_dir / f"nem_file_{kval}.uf") as uf,
        open(nem_dir / f"nem_file_{kval}.mf") as mf,
    ):
        params_lines = mf.readlines()
        log_likelihood = float(params_lines[2].split()[3])

        all_parameters: Dict[str, tuple] = {}
        for k, line in enumerate(params_lines[-kval:]):
            vector = line.split()
            mu_k = [bool(float(v)) for v in vector[:nb_org]]
            epsilon_k = [float(v) for v in vector[nb_org + 1:]]
            proportion = float(vector[nb_org])
            if k == 0:
                all_parameters["persistent"] = (mu_k, epsilon_k, proportion)
            elif k == kval - 1:
                all_parameters["cloud"] = (mu_k, epsilon_k, proportion)
            else:
                all_parameters[f"shell_{k}"] = (mu_k, epsilon_k, proportion)

        partitions: Dict[str, str] = {}
        entropy = 0.0
        for i, line in enumerate(uf):
            probs = [float(x) for x in line.split()]
            if just_log_likelihood:
                entropy += sum(math.log(p) * p for p in probs if p > 0)
            else:
                max_p = max(probs)
                best = [j for j, p in enumerate(probs) if p == max_p]
                partitions[index_fam[i]] = (
                    "S_" if (len(best) > 1 or max_p < 0.5) else labels[best[0]]
                )

    if just_log_likelihood:
        return kval, log_likelihood, entropy
    return partitions, all_parameters, log_likelihood


def _log_cnem_failure(
    logger,
    err,
    nem_dir: Path,
    kval: int,
    nb_org: int,
    beta: float,
    seed: int,
    free_dispersion: bool,
    init: str,
    itermax: int,
    just_log_likelihood: bool,
) -> None:
    if just_log_likelihood:
        logger.warning(
            "A C NEM run failed while estimating the optimal number of partitions "
            "(candidate K=%d skipped). See: %s", kval, nem_dir
        )
    else:
        logger.warning(
            "A C NEM run failed for this dataset/chunk and was skipped. See: %s",
            nem_dir,
        )
    logger.debug(
        "C NEM failure: error=%r, kval=%d, nb_org=%d, beta=%.4f, seed=%d, "
        "free_dispersion=%s, init=%s, itermax=%d, "
        "index_exists=%s, uf_exists=%s, mf_exists=%s",
        err, kval, nb_org, beta, seed, free_dispersion, init, itermax,
        (nem_dir / "nem_file.index").is_file(),
        (nem_dir / f"nem_file_{kval}.uf").is_file(),
        (nem_dir / f"nem_file_{kval}.mf").is_file(),
    )


def _run_cnem(
    nem_dir: Path,
    nb_org: int,
    beta: float,
    free_dispersion: bool,
    kval: int,
    seed: int,
    init: str,
    itermax: int,
    just_log_likelihood: bool,
) -> Union[Tuple[dict, dict, float], Tuple[int, float, float], Tuple[dict, None, None]]:
    """Run the C NEM extension and return partitioning results."""
    logger = logging.getLogger("PPanGGOLiN")

    if init in ("param_file", "init_from_old"):
        _write_cnem_init_file(nem_dir / f"nem_file_init_{kval}.m", kval, nb_org)

    init_random, init_param_file = 1, 2
    use_param_init = init in ("param_file", "init_from_old")

    logger.debug("Running C NEM…")
    nem_stats.nem(
        Fname=str(nem_dir / "nem_file").encode("ascii"),
        nk=kval,
        algo=b"nem",
        beta=beta,
        convergence=b"clas",
        convergence_th=0.01,
        format=b"fuzzy",
        it_max=itermax,
        dolog=True,
        model_family=b"bern",
        proportion=b"pk",
        dispersion=b"skd" if free_dispersion else b"sk_",
        init_mode=init_param_file if use_param_init else init_random,
        init_file=(str(nem_dir / f"nem_file_init_{kval}") + ".m").encode("ascii"),
        out_file_prefix=str(nem_dir / f"nem_file_{kval}").encode("ascii"),
        seed=seed,
    )

    uf_path = nem_dir / f"nem_file_{kval}.uf"
    if not uf_path.is_file():
        logger.debug(
            "C NEM produced no output (run may have failed): "
            "kval=%d, nb_org=%d, beta=%.4f, seed=%d", kval, nb_org, beta, seed
        )
        return ({}, None, None) if not just_log_likelihood else (kval, None, None)

    try:
        return _parse_cnem_output(nem_dir, kval, nb_org, just_log_likelihood)
    except (OSError, ValueError) as err:
        _log_cnem_failure(logger, err, nem_dir, kval, nb_org, beta, seed,
                          free_dispersion, init, itermax, just_log_likelihood)
        return ({}, None, None) if not just_log_likelihood else (kval, None, None)


# ─────────────────────────────────────────────────────────────────────────────
# Python NEM (pynem) backend
# ─────────────────────────────────────────────────────────────────────────────

def _run_pynem(
    nem_dir: Path,
    nb_org: int,
    beta: float,
    free_dispersion: bool,
    kval: int,
    seed: int,
    init: str,
    itermax: int,
    just_log_likelihood: bool,
) -> Union[Tuple[dict, dict, float], Tuple[int, float, float], Tuple[dict, None, None]]:
    """Run the pure-Python pynem backend and return partitioning results."""
    import numpy as np
    import pynem as _pynem
    from ppanggolin.nem.pynem_adapter import PPanGGOLiNEM

    logger = logging.getLogger("PPanGGOLiN")
    logger.debug("Running Python NEM (pynem)…")

    G = _pynem.io.read_graph(nem_dir / "nem_file")

    init_mode = "random"
    init_kwargs: dict = {}
    if init in ("param_file", "init_from_old"):
        step = 0.5 / math.ceil(kval / 2)
        pichenette = 0.1 if kval == 2 else 0.0
        centers = np.zeros((kval, nb_org))
        dispersions = np.zeros((kval, nb_org))
        for k in range(1, kval + 1):
            if k <= kval / 2:
                centers[k - 1] = 1.0
                dispersions[k - 1] = (step * k) - pichenette
            else:
                centers[k - 1] = 0.0
                dispersions[k - 1] = (step * (kval - k + 1)) - pichenette
        init_mode = "param_file"
        init_kwargs = dict(
            init_centers=centers,
            init_dispersions=dispersions,
            init_proportions=np.full(kval, 1.0 / kval),
        )

    try:
        model = PPanGGOLiNEM(
            n_clusters=kval,
            beta=beta,
            family="bernoulli",
            dispersion="skd" if free_dispersion else "sk_",
            proportion="pk",
            init=init_mode,
            max_iter=itermax,
            tol=0.01,
            convergence="classification",
            random_state=seed,
            verbose=0,
            **init_kwargs,
        )
        model.fit(G)
    except Exception as exc:
        logger.warning("Python NEM (pynem) run failed: %r", exc)
        return ({}, None, None) if not just_log_likelihood else (kval, None, None)

    index_fam = _read_family_index(nem_dir)

    if just_log_likelihood:
        ll = float(model.criteria_["L"])
        C = model.membership_
        entropy = float(np.sum(C * np.log(np.maximum(C, 1e-300))))
        return kval, ll, entropy

    labels = _partition_label_map(kval)
    C = model.membership_
    partitions: Dict[str, str] = {}
    for i, name in enumerate(index_fam):
        row = C[i]
        max_p = float(row.max())
        best = [j for j, p in enumerate(row) if p == max_p]
        partitions[name] = (
            "S_" if (len(best) > 1 or max_p < 0.5) else labels[best[0]]
        )

    all_parameters: Dict[str, tuple] = {}
    for k in range(kval):
        mu_k = [bool(round(v)) for v in model.centers_[k]]
        epsilon_k = list(map(float, model.dispersions_[k]))
        prop_k = float(model.proportions_[k])
        if k == 0:
            all_parameters["persistent"] = (mu_k, epsilon_k, prop_k)
        elif k == kval - 1:
            all_parameters["cloud"] = (mu_k, epsilon_k, prop_k)
        else:
            all_parameters[f"shell_{k}"] = (mu_k, epsilon_k, prop_k)

    return partitions, all_parameters, float(model.criteria_["L"])


# ─────────────────────────────────────────────────────────────────────────────
# Unified dispatcher
# ─────────────────────────────────────────────────────────────────────────────

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
    backend: str = "cnem",
    # Deprecated: use backend="pynem" instead
    use_pynem: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """Run NEM partitioning using the requested backend.

    :param nem_dir_path: directory containing NEM input files
    :param nb_org: number of organisms
    :param beta: spatial-smoothing strength (0 → standard EM)
    :param free_dispersion: per-organism dispersion per cluster
    :param kval: number of clusters
    :param seed: random seed
    :param init: initialisation mode ('param_file' or 'random')
    :param keep_files: retain intermediate NEM files
    :param itermax: maximum EM iterations
    :param just_log_likelihood: return only log-likelihood statistics (for K selection)
    :param backend: 'cnem' (C extension, default) or 'pynem' (pure Python)
    :param use_pynem: deprecated alias for backend='pynem'
    :return:
        * ``just_log_likelihood=False``: ``(partitions_dict, parameters_dict, log_likelihood)``
        * ``just_log_likelihood=True``:  ``(kval, log_likelihood, entropy)``
    """
    if use_pynem and backend == "cnem":
        backend = "pynem"
    runner = _run_pynem if backend == "pynem" else _run_cnem
    return runner(
        nem_dir_path, nb_org, beta, free_dispersion,
        kval, seed, init, itermax, just_log_likelihood,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Multiprocessing workers
# ─────────────────────────────────────────────────────────────────────────────

def _partition_sample_worker(args: tuple):
    """Pool worker: write NEM files for one organism sub-sample and partition it.

    Uses the module-level ``_samples`` list and ``_pangenome`` object inherited
    by forked child processes at fork() time.
    """
    (
        index, kval, beta, sm_degree, free_dispersion,
        seed, init, tmpdir, keep_tmp_files, backend,
    ) = args
    currtmpdir = tmpdir / str(index)
    samp = _samples[index]
    edges_weight, nb_fam = write_nem_input_files(currtmpdir, samp, sm_degree)
    return run_partitioning(
        currtmpdir,
        len(samp),
        beta * (nb_fam / edges_weight),
        free_dispersion,
        kval=kval,
        seed=seed,
        init=init,
        keep_files=keep_tmp_files,
        backend=backend,
    )


def nem_single(args) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """Legacy wrapper: call run_partitioning with positional args tuple.

    :param args: positional arguments for run_partitioning
    :return: Result of run_partitioning
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
    use_pynem: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """Legacy wrapper; prefer _partition_sample_worker for new code.

    :param index: Index of the sample group
    :param tmpdir: temporary directory path
    :param kval: Number of partitions to use
    :param beta: spatial-smoothing strength
    :param sm_degree: maximum adjacency degree for smoothing
    :param free_dispersion: allow per-organism dispersion per partition
    :param seed: random seed
    :param init: initialisation mode
    :param keep_tmp_files: retain temporary NEM files
    :param use_pynem: deprecated; sets backend to 'pynem'
    :return: partitioning result tuple
    """
    backend = "pynem" if use_pynem else "cnem"
    return _partition_sample_worker(
        (index, kval, beta, sm_degree, free_dispersion,
         seed, init, tmpdir, keep_tmp_files, backend)
    )


def nem_samples(pack: tuple) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """Legacy wrapper for _partition_sample_worker.

    :param pack: argument tuple for _partition_sample_worker
    :return: partitioning result tuple
    """
    return _partition_sample_worker(pack)


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

        for fam in _pangenome.gene_families:
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


def _draw_icl_curve(
    all_bics: dict,
    all_icls: dict,
    all_lls: dict,
    best_k: int,
    max_icl_k: int,
    icl_margin: float,
    krange: list,
    output: Path,
) -> None:
    """Write a Plotly ICL/BIC/log-likelihood curve to *output* as HTML."""
    x_sorted = sorted(all_bics.keys())
    traces = [
        go.Scatter(x=x_sorted, y=[all_bics[k] for k in x_sorted],
                   name="BIC", mode="lines+markers"),
        go.Scatter(x=x_sorted, y=[all_icls[k] for k in x_sorted],
                   name="ICL", mode="lines+markers"),
        go.Scatter(x=x_sorted, y=[all_lls[k] for k in x_sorted],
                   name="log likelihood", mode="lines+markers"),
    ]
    layout = go.Layout(
        title=f"ICL curve (best K={best_k}, ICL_margin={icl_margin})",
        titlefont=dict(size=20),
        xaxis=dict(title="number of partitions"),
        yaxis=dict(title="ICL / BIC / log likelihood"),
        plot_bgcolor="#ffffff",
        shapes=[
            dict(type="line", x0=best_k, x1=best_k, y0=0,
                 y1=all_icls[best_k],
                 line=dict(width=1, dash="dashdot", color="black")),
            dict(type="line", x0=max_icl_k, x1=max_icl_k, y0=0,
                 y1=all_icls[max_icl_k],
                 line=dict(width=1, dash="dashdot", color="black")),
            dict(type="line", x0=best_k, x1=max_icl_k,
                 y0=all_icls[max_icl_k], y1=all_icls[max_icl_k],
                 line=dict(width=1, dash="dashdot", color="black")),
            dict(type="line", x0=2, x1=krange[1],
                 y0=all_icls[best_k], y1=all_icls[best_k],
                 line=dict(width=1, dash="dashdot", color="black")),
        ],
    )
    fig = go.Figure(data=traces, layout=layout)
    out_plotly.plot(fig, filename=str(output / f"ICL_curve_K{best_k}.html"),
                    auto_open=False)


def evaluate_nb_partitions(
    organisms: set,
    output: Optional[Path] = None,
    sm_degree: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    krange: Optional[list] = None,
    icl_margin: float = 0.05,
    draw_icl: bool = False,
    cpu: int = 1,
    seed: int = 42,
    tmpdir: Optional[Path] = None,
    disable_bar: bool = False,
    backend: str = "cnem",
    # Deprecated: use backend="pynem" instead
    use_pynem: bool = False,
) -> int:
    """Evaluate the optimal number of partitions for the pangenome.

    :param organisms: organisms to include in the evaluation
    :param output: directory for the optional ICL plot
    :param sm_degree: maximum adjacency degree for smoothing
    :param free_dispersion: allow per-organism dispersion per partition
    :param chunk_size: if more organisms than this, sub-sample for K evaluation
    :param krange: [K_min, K_max] search range
    :param icl_margin: fraction of max-ICL range used to prefer lower K
    :param draw_icl: write the ICL curve HTML
    :param cpu: parallel workers
    :param seed: random seed
    :param tmpdir: base temporary directory
    :param disable_bar: suppress tqdm bar
    :param backend: 'cnem' or 'pynem'
    :param use_pynem: deprecated alias for backend='pynem'
    :return: optimal K (≥ 3)
    """
    if use_pynem and backend == "cnem":
        backend = "pynem"

    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    krange = krange or [3, 20]
    newtmpdir = tmpdir / "eval_partitions"

    if len(organisms) > chunk_size:
        # Sort by name for reproducible sampling
        select_organisms = set(
            random.sample(sorted(organisms, key=lambda o: o.name), chunk_size)
        )
    else:
        select_organisms = set(organisms)

    _, nb_fam = write_nem_input_files(newtmpdir, select_organisms, sm_degree)

    # Build one args tuple per candidate K (positional → run_partitioning)
    args_partionning = [
        (
            newtmpdir,
            len(select_organisms),
            0,            # beta=0: no spatial smoothing during K scan
            free_dispersion,
            k,
            seed,
            "param_file",
            True,         # keep_files
            10,           # itermax (coarse scan)
            True,         # just_log_likelihood
            backend,
        )
        for k in range(krange[0] - 1, krange[1] + 1)
    ]

    all_log_likelihood = []
    if cpu > 1:
        bar = tqdm(
            range(len(args_partionning)),
            unit="Number of partitions",
            disable=disable_bar,
        )
        with get_context("fork").Pool(processes=cpu) as p:
            for result in p.imap_unordered(nem_single, args_partionning):
                all_log_likelihood.append(result)
                bar.update()
            p.close()
            p.join()
        bar.close()
    else:
        for arguments in args_partionning:
            all_log_likelihood.append(nem_single(arguments))

    all_bics: Dict[int, float] = defaultdict(float)
    all_icls: Dict[int, float] = defaultdict(float)
    all_lls: Dict[int, float] = defaultdict(float)

    for k_candidate, log_likelihood, entropy in all_log_likelihood:
        if log_likelihood is not None:
            nb_params = k_candidate * (
                len(select_organisms)
                + 1
                + (len(select_organisms) if free_dispersion else 1)
            )
            bic = log_likelihood - 0.5 * math.log(nb_params) * nb_fam
            all_bics[k_candidate] = bic
            all_icls[k_candidate] = bic - entropy
            all_lls[k_candidate] = log_likelihood

    chosen_k = 3
    best_k = chosen_k
    max_icl_k = chosen_k

    if len(all_bics) > 3:
        max_icl_k = max(all_icls, key=all_icls.get)
        delta_icl = (all_icls[max_icl_k] - min(all_icls.values())) * icl_margin
        best_k = min(
            k for k, icl in all_icls.items()
            if icl >= all_icls[max_icl_k] - delta_icl and k <= max_icl_k
        )
        chosen_k = best_k if best_k >= 3 else chosen_k

    if all_bics and draw_icl and output is not None:
        _draw_icl_curve(all_bics, all_icls, all_lls, best_k, max_icl_k,
                        icl_margin, krange, output)

    return chosen_k


def check_pangenome_former_partition(pangenome: Pangenome, force: bool = False) -> None:
    """Check pangenome status for former partitions; erase them if *force* is True.

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


# ─────────────────────────────────────────────────────────────────────────────
# Chunked partitioning helper
# ─────────────────────────────────────────────────────────────────────────────

def _run_chunked_partitioning(
    organisms: set,
    tmp_path: Path,
    kval: int,
    beta: float,
    sm_degree: int,
    free_dispersion: bool,
    seed: int,
    init: str,
    keep_tmp_files: bool,
    backend: str,
    chunk_size: int,
    cpu: int,
    disable_bar: bool,
) -> Dict[str, str]:
    """Partition by repeated random sub-sampling and majority vote.

    Organisms are shuffled into chunks of *chunk_size*.  Each chunk is
    partitioned independently.  A family's final label is the one that
    receives ≥ 50 % of votes after sufficient sampling.

    :return: mapping family_name → partition_label
    """
    logger = logging.getLogger("PPanGGOLiN")

    # Vote tallies: first character of label ('P', 'S', 'C', or 'U')
    cpt_partition: Dict[str, Dict[str, int]] = {}
    families = set()
    for fam in _pangenome.gene_families:
        families.add(fam)
        cpt_partition[fam.name] = {"P": 0, "S": 0, "C": 0, "U": 0}

    pansize = len(families)
    validated: set = set()
    org_nb_sample: Counter = Counter({org: 0 for org in organisms})
    condition = len(organisms) / chunk_size
    start = time.time()

    def _tally_result(res) -> None:
        """Record votes from one chunk result and mark families validated."""
        for node, nem_class in res[0].items():
            cpt_partition[node][nem_class[0]] += 1
            sum_partionning = sum(cpt_partition[node].values())
            if (
                sum_partionning > len(organisms) / chunk_size
                and max(cpt_partition[node].values()) >= sum_partionning * 0.5
            ) or (sum_partionning > len(organisms)):
                if node not in validated:
                    if max(cpt_partition[node].values()) < sum_partionning * 0.5:
                        # No majority even after exhaustive sampling → undefined
                        cpt_partition[node]["U"] = len(organisms)
                    validated.add(node)

    while len(validated) < pansize:
        prev = len(_samples)
        # Ensure every organism appears in at least `condition` samples
        while not all(val >= condition for val in org_nb_sample.values()):
            shuffled_orgs = list(organisms)
            random.shuffle(shuffled_orgs)
            while len(shuffled_orgs) > chunk_size:
                _samples.append(set(shuffled_orgs[:chunk_size]))
                for org in _samples[-1]:
                    org_nb_sample[org] += 1
                shuffled_orgs = shuffled_orgs[chunk_size:]

        args = [
            (i, kval, beta, sm_degree, free_dispersion,
             seed, init, tmp_path, keep_tmp_files, backend)
            for i in range(prev, len(_samples))
        ]

        logger.info("Launching NEM on %d chunk samples (backend=%s, cpu=%d)…",
                    len(args), backend, cpu)
        with get_context("fork").Pool(processes=cpu) as p:
            bar = tqdm(range(len(args)), unit=" samples partitioned",
                       disable=disable_bar)
            for result in p.imap_unordered(nem_samples, args):
                _tally_result(result)
                bar.update()
            bar.close()
            condition += 1
            logger.debug(
                "%d validated families out of %d.", len(validated), pansize
            )
            p.close()
            p.join()

    logger.info(
        "Did %d partitionings with chunks of size %d among %d genomes in %.1f s.",
        len(_samples), chunk_size, len(organisms), time.time() - start,
    )
    return {fam: max(data, key=data.get) for fam, data in cpt_partition.items()}


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def partition(
    pangenome: Pangenome,
    output: Optional[Path] = None,
    beta: float = 2.5,
    sm_degree: int = 10,
    free_dispersion: bool = False,
    chunk_size: int = 500,
    kval: int = -1,
    krange: Optional[list] = None,
    icl_margin: float = 0.05,
    draw_icl: bool = False,
    cpu: int = 1,
    seed: int = 42,
    tmpdir: Optional[Path] = None,
    keep_tmp_files: bool = False,
    force: bool = False,
    disable_bar: bool = False,
    backend: str = "cnem",
) -> None:
    """Assign every gene family in *pangenome* to a partition.

    Partitions are Persistent (P), Shell (S1…Sk-2), and Cloud (C) according
    to the Neighborhood EM model.

    :param pangenome: pangenome to partition
    :param output: output directory (required when *draw_icl* or *keep_tmp_files*)
    :param beta: spatial-smoothing strength; 0 disables smoothing
    :param sm_degree: maximum adjacency degree included in the smoothing
    :param free_dispersion: allow per-organism dispersion per partition
    :param chunk_size: max organisms per NEM run; excess triggers chunked mode
    :param kval: number of partitions (< 2 → auto-detected)
    :param krange: [K_min, K_max] search range for auto K
    :param icl_margin: tolerance margin for choosing the lowest adequate K
    :param draw_icl: write an ICL curve HTML to *output*
    :param cpu: parallel worker count
    :param seed: random seed
    :param tmpdir: base directory for temporary files
    :param keep_tmp_files: copy NEM temporary files to *output*
    :param force: overwrite an existing partition in the pangenome
    :param disable_bar: suppress tqdm progress bars
    :param backend: 'cnem' (C extension, default) or 'pynem' (pure Python)
    :param use_pynem: deprecated alias for backend='pynem'
    """
    global _pangenome, _samples

    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    kmm = [3, 20] if krange is None else krange

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
    if len(organisms) <= 10:
        logging.getLogger("PPanGGOLiN").warning(
            f"The number of selected genomes is too low ({len(organisms)} "
            f"genomes used) to robustly partition the graph"
        )

    # Expose pangenome and sample list via module globals so forked workers
    # can access them without serialisation.
    _pangenome = pangenome
    _samples = []
    random.seed(seed)

    if keep_tmp_files:
        tmp_path = Path(tempfile.mkdtemp(dir=tmpdir))
        tmp_cleanup = None
    else:
        tmp_dir_obj = tempfile.TemporaryDirectory(dir=tmpdir)
        tmp_path = Path(tmp_dir_obj.name)
        tmp_cleanup = tmp_dir_obj

    pangenome.parameters["partition"] = {
        "beta": beta,
        "max_degree_smoothing": sm_degree,
        "free_dispersion": free_dispersion,
        "ICL_margin": icl_margin,
        "seed": seed,
        "nb_of_partitions": kval,
        "# computed nb of partitions": False,
        "krange": kmm,
    }
    if len(organisms) > chunk_size:
        pangenome.parameters["partition"]["chunk_size"] = chunk_size

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
            backend=backend,
        )
        logging.getLogger("PPanGGOLiN").info(
            f"The number of partitions has been evaluated at {kval}"
        )

    pangenome.parameters["partition"]["# final nb of partitions"] = kval
    init = "param_file"

    logging.getLogger("PPanGGOLiN").info("Partitioning...")
    start_partitioning = time.time()

    if chunk_size < len(organisms):
        # ── Chunked mode ───────────────────────────────────────────────────
        family_partitions = _run_chunked_partitioning(
            organisms=organisms,
            tmp_path=tmp_path,
            kval=kval,
            beta=beta,
            sm_degree=sm_degree,
            free_dispersion=free_dispersion,
            seed=seed,
            init=init,
            keep_tmp_files=keep_tmp_files,
            backend=backend,
            chunk_size=chunk_size,
            cpu=cpu,
            disable_bar=disable_bar,
        )
    else:
        # ── Single run ─────────────────────────────────────────────────────
        nem_dir = tmp_path / "0"
        edges_weight, nb_fam = write_nem_input_files(nem_dir, organisms, sm_degree)
        result = run_partitioning(
            nem_dir,
            len(organisms),
            beta * (nb_fam / edges_weight),
            free_dispersion,
            kval=kval,
            seed=seed,
            init=init,
            keep_files=keep_tmp_files,
            backend=backend,
        )
        if not result[0]:
            raise Exception(
                "Statistical partitioning does not work on your data. "
                "This usually happens because you used very few (<15) genomes."
            )
        family_partitions = result[0]
        logging.getLogger("PPanGGOLiN").info(
            f"Partitioned {len(organisms)} genomes in "
            f"{round(time.time() - start_partitioning, 2)} seconds."
        )

    for fam_name, part in family_partitions.items():
        pangenome.get_gene_family(fam_name).partition = part

    pangenome.status["partitioned"] = "Computed"

    if keep_tmp_files and output is not None:
        copytree(tmp_path, output / "NEM_files/")
    elif tmp_cleanup is not None:
        tmp_cleanup.cleanup()


def launch(args: argparse.Namespace):
    """Command launcher.

    :param args: All arguments provided by user
    """
    if args.draw_ICL or args.keep_tmp_files:
        mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    backend = "pynem" if getattr(args, "use_pynem", False) else "cnem"
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
        args.tmpdir,
        args.keep_tmp_files,
        args.force,
        disable_bar=args.disable_prog_bar,
        backend=backend,
    )
    logging.getLogger("PPanGGOLiN").debug("Write partition in pangenome")
    write_pangenome(pangenome, pangenome.file, args.force,
                    disable_bar=args.disable_prog_bar)
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
    optional.add_argument(
        "--use_pynem",
        required=False,
        default=False,
        action="store_true",
        help="Use the pure-Python pynem backend instead of the default C NEM implementation. "
        "Requires the pynem package to be installed.",
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
