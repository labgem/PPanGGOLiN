#!/usr/bin/env python3

# default libraries
from dataclasses import dataclass

# installed libraries
import numpy as np
from pynem import NEM

# local libraries
from ppanggolin.nem.nem_input import NemInput


@dataclass
class NemResult:
    """Outcome of one NEM partitioning.

    :param partitions: Label of each gene family: P, S<k>, S_ or C.
    :param parameters: Centers, dispersions and proportion of each class.
    :param log_likelihood: Log-likelihood of the fit.
    """

    partitions: dict[str, str]
    parameters: dict[str, tuple[list[bool], list[float], float]]
    log_likelihood: float


@dataclass
class KScore:
    """Fit of NEM with k partitions, used to choose k.

    :param k: Number of partitions.
    :param log_likelihood: Log-likelihood of the fit.
    :param entropy: Entropy of the fuzzy classification.
    """

    k: int
    log_likelihood: float
    entropy: float


def _initial_parameters(k, nb_genomes):
    """
    Structured parameter initialization.

    :param k: number of partitions
    :param nb_genomes: number of genomes (feature dimension D)

    :return: initial modes (k, nb_genomes), dispersions (k, nb_genomes) and proportions (k,)
    """
    p = round(1.0 / k, 2)
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


def _fit(
    nem_input: NemInput,
    beta: float,
    free_dispersion: bool,
    k: int,
    max_iter: int
) -> NEM:
    """Fit NEM on an input.

    :param nem_input: Input to partition.
    :param beta: Smoothing strength, already scaled.
    :param free_dispersion: Free dispersion per genome (`skd`) instead of shared (`sk_`).
    :param k: Number of partitions.
    :param max_iter: Maximum number of iterations.
    :return: The fitted model.
    """
    model = NEM(
        n_clusters=k,
        beta=beta,
        family="bernoulli",
        dispersion="skd" if free_dispersion else "sk_",
        proportion="pk",
        algorithm="nem",
        init="param",
        init_params=_initial_parameters(k, nem_input.nb_genomes),
        site_update="seq",
        convergence="classification",
        tol=0.01,
        missing="ignore",
        max_iter=max_iter,
    )
    model.fit(nem_input.presence, graph=nem_input.graph)
    return model


def solve(
    nem_input: NemInput,
    beta: float = 2.5,
    free_dispersion: bool = False,
    k: int = 3,
    max_iter: int = 100,
) -> NemResult:
    """Partition the gene families of an input.

    :param nem_input: Input to partition.
    :param beta: Smoothing strength, already scaled.
    :param free_dispersion: Free dispersion per genome instead of shared.
    :param k: Number of partitions.
    :param max_iter: Maximum number of iterations.
    :return: The partitions and the fitted parameters.
    """
    model = _fit(nem_input, beta, free_dispersion, k, max_iter)

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
    for row in model.membership_:
        max_prob = float(row.max())
        best = [c for c, prob in enumerate(row) if prob == max_prob]
        if len(best) > 1 or max_prob < 0.5:
            family_labels.append("S_")  # shell in case of doubt
        else:
            family_labels.append(labels[best[0]])

    return NemResult(
        partitions=dict(zip(nem_input.fam_names, family_labels)),
        parameters=all_parameters,
        log_likelihood=float(model.criteria_["M"]),
    )


def score(
    nem_input: NemInput,
    beta: float = 0.0,
    free_dispersion: bool = False,
    k: int = 3,
    max_iter: int = 100,
) -> KScore:
    """Fit NEM with k partitions and report how well it fits.

    :param nem_input: Input to partition.
    :param beta: Smoothing strength, already scaled.
    :param free_dispersion: Free dispersion per genome instead of shared.
    :param k: Number of partitions.
    :param max_iter: Maximum number of iterations.
    :return: The log-likelihood and classification entropy.
    """
    model = _fit(nem_input, beta, free_dispersion, k, max_iter)
    membership = model.membership_
    positive = membership[membership > 0]
    return KScore(
        k=k,
        log_likelihood=float(model.criteria_["M"]),
        entropy=float(np.sum(positive * np.log(positive))),
    )
