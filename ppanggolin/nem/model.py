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


def _initial_parameters(k, nb_genomes):
    """
    Structured parameter initialization.

    :param k: number of partitions
    :param nb_genomes: number of genomes (feature dimension D)

    :return: initial modes (k, nb_genomes), dispersions (k, nb_genomes) and proportions (k,)
    """
    p = float(format(1.0 / k, ".12g"))
    proportions = np.array([p] * (k - 1) + [1.0 - p * (k - 1)], dtype=float)

    centers = np.zeros((k, nb_genomes))
    dispersions = np.zeros((k, nb_genomes))
    step = 0.5 / ((k + 1) // 2)
    # Nudge dispersions away from the 0.5 boundary when there are only two classes.
    offset = 0.1 if k == 2 else 0.0
    for c in range(1, k + 1):
        if c <= k / 2:
            mu, eps = 1.0, step * c - offset
        else:
            mu, eps = 0.0, step * (k - c + 1) - offset
        centers[c - 1, :] = mu
        dispersions[c - 1, :] = eps
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
    fam_names: list,
    nb_genomes: int,
    beta: float = 2.5,
    free_dispersion: bool = False,
    k: int = 3,
    max_iter: int = 100,
    scores_only: bool = False,
) -> Union[Tuple[dict, None, None], Tuple[int, float, float], Tuple[dict, dict, float]]:
    """
    Run NEM on pre-built inputs

    :param presence: (nb_families, nb_genomes) presence/absence matrix
    :param graph: neighbourhood graph, nodes 0..nb_families-1, edges carrying `weight`
    :param fam_names: gene-family names in row order of `presence`
    :param nb_genomes: number of genomes
    :param beta: spatial smoothing coefficient
    :param free_dispersion: whether dispersion is free per genome (`skd`) or shared (`sk_`)
    :param k: number of partitions
    :param max_iter: maximum number of iterations
    :param scores_only: return only the log-likelihood diagnostics

    :return: partition per family, per-class parameters, log-likelihood --
        or, with `scores_only`, (k, log-likelihood, entropy)
    """
    model = NEM(
        n_clusters=k,
        beta=beta,
        family="bernoulli",
        dispersion="skd" if free_dispersion else "sk_",
        proportion="pk",
        algorithm="nem",
        init="param",
        init_params=_initial_parameters(k, nb_genomes),
        site_update="seq",
        convergence="classification",
        tol=0.01,
        missing="ignore",
        max_iter=max_iter,
    )
    model.fit(presence, graph=graph)

    membership = model.membership_
    log_likelihood = float(model.criteria_["M"])

    if scores_only:
        positive = membership[membership > 0]
        entropy = float(np.sum(positive * np.log(positive)))
        return k, log_likelihood, entropy

    all_parameters = {}
    for c in range(k):
        mu = [bool(round(float(mu_j))) for mu_j in model.centers_[c]]
        epsilon = [float(epsilon_j) for epsilon_j in model.dispersions_[c]]
        proportion = float(model.proportions_[c])
        if c == 0:
            all_parameters["persistent"] = (mu, epsilon, proportion)
        elif c == k - 1:
            all_parameters["cloud"] = (mu, epsilon, proportion)
        else:
            all_parameters["shell_" + str(c)] = (mu, epsilon, proportion)

    labels = {0: "P", k - 1: "C"}
    for i in range(1, k - 1):
        labels[i] = "S" + str(i)
    family_labels = []
    for row in membership:
        max_prob = float(row.max())
        best = [c for c, prob in enumerate(row) if prob == max_prob]
        if len(best) > 1 or max_prob < 0.5:
            family_labels.append("S_")
        else:
            family_labels.append(labels[best[0]])

    return dict(zip(fam_names, family_labels)), all_parameters, log_likelihood
