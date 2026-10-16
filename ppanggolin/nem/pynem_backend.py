#!/usr/bin/env python3

# default libraries
import os
from typing import Tuple, Union

# installed libraries
import networkx as nx
import numpy as np
from pynem import NEM
from threadpoolctl import threadpool_limits

try:
    import numba
except ImportError:
    numba = None

# pynem 0.6.1's parallel kernels run on OpenMP, which kills workers forked after them
if numba is not None and "NUMBA_THREADING_LAYER" not in os.environ:
    numba.config.THREADING_LAYER = "forksafe"


def _build_param_init(K, n_org):
    """
    Structured parameter initialization.

    :param K: number of partitions
    :param n_org: number of genomes (feature dimension D)

    :return: initial modes (K, n_org), dispersions (K, n_org) and proportions (K,)
    """
    p = float(format(1.0 / K, ".12g"))
    proportions = np.array([p] * (K - 1) + [1.0 - p * (K - 1)], dtype=float)

    centers = np.zeros((K, n_org))
    dispersions = np.zeros((K, n_org))
    step = 0.5 / ((K + 1) // 2)
    # Nudge dispersions away from the 0.5 boundary when there are only two classes.
    offset = 0.1 if K == 2 else 0.0
    for k in range(1, K + 1):
        if k <= K / 2:
            mu, eps = 1.0, step * k - offset
        else:
            mu, eps = 0.0, step * (K - k + 1) - offset
        centers[k - 1, :] = mu
        dispersions[k - 1, :] = eps
    return centers, dispersions, proportions


def set_threads(nb_threads: int):
    """Let pynem's numba kernels and BLAS use at most `nb_threads` threads in this process.

    :param nb_threads: Number of threads.
    """
    if numba is not None:
        numba.set_num_threads(max(1, min(nb_threads, numba.config.NUMBA_NUM_THREADS)))
    threadpool_limits(nb_threads)


def solve(
    presence: np.ndarray,
    graph: nx.DiGraph,
    index_fam: list,
    nb_org: int,
    beta: float = 2.5,
    free_dispersion: bool = False,
    kval: int = 3,
    itermax: int = 100,
    just_log_likelihood: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """
    Run NEM on pre-built inputs

    :param presence: (nb_fam, nb_org) presence/absence matrix
    :param graph: neighbourhood graph, nodes 0..nb_fam-1, edges carrying `weight`
    :param index_fam: gene-family names in row order of `presence`
    :param nb_org: number of genomes
    :param beta: spatial smoothing coefficient
    :param free_dispersion: whether dispersion is free per genome (`skd`) or shared (`sk_`)
    :param kval: number of partitions
    :param itermax: maximum number of iterations
    :param just_log_likelihood: return only the log-likelihood diagnostics

    :return: partition per family, per-class parameters, log-likelihood --
        or, with `just_log_likelihood`, (kval, log-likelihood, entropy)
    """
    model = NEM(
        n_clusters=kval,
        beta=beta,
        family="bernoulli",
        dispersion="skd" if free_dispersion else "sk_",
        proportion="pk",
        algorithm="nem",
        init="param",
        init_params=_build_param_init(kval, nb_org),
        site_update="seq",
        convergence="classification",
        tol=0.01,
        missing="ignore",
        max_iter=itermax,
    )
    model.fit(presence, graph=graph)

    membership = model.membership_
    log_likelihood = float(model.criteria_["M"])

    if just_log_likelihood:
        positive = membership[membership > 0]
        entropy = float(np.sum(positive * np.log(positive)))
        return kval, log_likelihood, entropy

    all_parameters = {}
    for k in range(kval):
        mu_k = [bool(round(float(mu_kj))) for mu_kj in model.centers_[k]]
        epsilon_k = [float(epsilon_kj) for epsilon_kj in model.dispersions_[k]]
        proportion = float(model.proportions_[k])
        if k == 0:
            all_parameters["persistent"] = (mu_k, epsilon_k, proportion)
        elif k == kval - 1:
            all_parameters["cloud"] = (mu_k, epsilon_k, proportion)
        else:
            all_parameters["shell_" + str(k)] = (mu_k, epsilon_k, proportion)

    parti = {0: "P", kval - 1: "C"}
    for i in range(1, kval - 1):
        parti[i] = "S" + str(i)
    partitions_list = []
    for row in membership:
        max_prob = float(row.max())
        positions_max_prob = [pos for pos, prob in enumerate(row) if prob == max_prob]
        if len(positions_max_prob) > 1 or max_prob < 0.5:
            partitions_list.append("S_")
        else:
            partitions_list.append(parti[positions_max_prob[0]])

    return dict(zip(index_fam, partitions_list)), all_parameters, log_likelihood
