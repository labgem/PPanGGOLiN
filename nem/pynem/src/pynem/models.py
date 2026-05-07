"""Distribution families and parameter estimation for NEM."""

import numpy as np
from enum import Enum

EPSILON = 1e-20


class Family(Enum):
    NORMAL = "normal"
    LAPLACE = "laplace"
    BERNOULLI = "bernoulli"


class Dispersion(Enum):
    S__ = "s__"   # single dispersion for all
    SK_ = "sk_"   # per-class
    S_D = "s_d"   # per-variable
    SKD = "skd"   # full (per-class, per-variable)


class Proportion(Enum):
    EQUAL = "p_"  # equal proportions = 1/K
    FREE = "pk"   # free proportions


def compute_log_density(X, centers, dispersions, proportions, family):
    """Compute log(p_k * f_k(x_i)) for all i and k.

    Parameters
    ----------
    X : (N, D) array
    centers : (K, D) array
    dispersions : (K, D) array
    proportions : (K,) array
    family : Family enum

    Returns
    -------
    log_pkfki : (N, K) array
        log(p_k * f_k(x_i)) for each observation i and class k.

    Notes
    -----
    All three families are fully vectorised over (N, K, D) using numpy
    broadcasting — no Python loops over K or D.  NaN entries in *X* are
    treated as missing observations and contribute 0 to the log-density.
    """
    obs = ~np.isnan(X)  # (N, D) — True where observed
    X_filled = np.where(obs, X, 0.0)  # (N, D) — NaN replaced by 0 for arithmetic
    obs3d = obs[:, np.newaxis, :]  # (N, 1, D) → broadcasts to (N, K, D)

    log_pk = np.log(np.maximum(proportions, EPSILON))  # (K,)

    if family == Family.NORMAL:
        v = np.clip(dispersions, EPSILON, None)  # (K, D)
        diff = X_filled[:, np.newaxis, :] - centers[np.newaxis, :, :]  # (N, K, D)
        contrib = np.where(
            obs3d,
            -0.5 * (np.log(2.0 * np.pi * v[np.newaxis]) + diff**2 / v[np.newaxis]),
            0.0,
        )  # (N, K, D)

    elif family == Family.LAPLACE:
        v = np.clip(dispersions, EPSILON, None)  # (K, D)
        absdif = np.abs(X_filled[:, np.newaxis, :] - centers[np.newaxis, :, :])
        contrib = np.where(
            obs3d,
            -(np.log(2.0 * v[np.newaxis]) + absdif / v[np.newaxis]),
            0.0,
        )  # (N, K, D)

    elif family == Family.BERNOULLI:
        eps = np.clip(dispersions, EPSILON, 1.0 - EPSILON)  # (K, D)
        log1me = np.log1p(-eps)  # (K, D) = log(1-eps)
        log_ratio = np.log((1.0 - eps) / eps)  # (K, D)
        absdif = np.abs(X_filled[:, np.newaxis, :] - centers[np.newaxis, :, :])
        contrib = np.where(
            obs3d,
            log1me[np.newaxis] - absdif * log_ratio[np.newaxis],
            0.0,
        )  # (N, K, D)

    else:
        raise ValueError(f"Unknown family: {family}")

    log_fki = contrib.sum(axis=2)  # (N, K)
    return log_pk[np.newaxis, :] + log_fki  # (N, K)


def estimate_parameters(
    X,
    C,
    family,
    dispersion_model,
    proportion_model,
    miss_mode="replace",
    old_centers=None,
    old_dispersions=None,
):
    """M-step: estimate model parameters from soft classification.

    Parameters
    ----------
    X : (N, D) array
    C : (N, K) array — soft classification
    family : Family enum
    dispersion_model : Dispersion enum
    proportion_model : Proportion enum
    miss_mode : str — 'replace' or 'ignore'
    old_centers : (K, D) array or None
    old_dispersions : (K, D) array or None

    Returns
    -------
    dict with keys 'centers', 'dispersions', 'proportions'.
    """
    N, D = X.shape
    K = C.shape[1]
    observed = ~np.isnan(X)  # (N, D)
    obs3d = observed[:, np.newaxis, :]  # (N, 1, D) → broadcasts to (N, K, D)
    X_filled = np.where(observed, X, 0.0)  # (N, D)

    # Class sizes
    N_K = np.maximum(C.sum(axis=0), EPSILON)  # (K,)

    # Per-class, per-variable observed sizes — vectorised
    N_KD = np.maximum(
        (C[:, :, np.newaxis] * obs3d).sum(axis=0),  # (K, D)
        EPSILON,
    )

    # --- Centers ---
    if family == Family.LAPLACE:
        centers = _estimate_laplace_centers(X, C, observed, K, D, N_K)
    else:
        centers = _estimate_mean_centers(
            X_filled, observed, obs3d, C, K, D, N_K, N_KD, miss_mode, old_centers
        )

    # --- Inertia (vectorised over K) ---
    diff = X_filled[:, np.newaxis, :] - centers[np.newaxis, :, :]  # (N, K, D)
    diff_masked = diff * obs3d  # zero out unobserved positions
    if family == Family.NORMAL:
        Iner_KD = (C[:, :, np.newaxis] * diff_masked**2).sum(axis=0)  # (K, D)
        # Missing-data correction (Normal + replace mode only)
        if miss_mode == "replace" and old_dispersions is not None:
            n_miss = np.maximum(N_K[:, np.newaxis] - N_KD, 0.0)  # (K, D)
            Iner_KD += n_miss * old_dispersions
    else:  # Laplace or Bernoulli
        Iner_KD = (C[:, :, np.newaxis] * np.abs(diff_masked)).sum(axis=0)  # (K, D)

    # --- Dispersions ---
    dispersions = np.maximum(
        _inertia_to_dispersions(
            Iner_KD, N_K, N_KD, N, D, K, dispersion_model, miss_mode
        ),
        EPSILON,
    )

    # --- Proportions ---
    if proportion_model == Proportion.EQUAL:
        proportions = np.full(K, 1.0 / K)
    else:
        proportions = N_K / N
        proportions = np.maximum(proportions, EPSILON)
        proportions /= proportions.sum()

    return {
        "centers": centers,
        "dispersions": dispersions,
        "proportions": proportions,
    }


def _estimate_mean_centers(
    X_filled, observed, obs3d, C, K, D, N_K, N_KD, miss_mode, old_centers
):
    """Weighted mean center estimation (Normal, Bernoulli) — vectorised."""
    # weighted_sum[k, d] = sum_i c_ik * X_filled[i, d] * obs[i, d]
    weighted_sum = (C[:, :, np.newaxis] * X_filled[:, np.newaxis, :] * obs3d).sum(
        axis=0
    )  # (K, D)

    if miss_mode == "replace" and old_centers is not None:
        n_miss = np.maximum(N_K[:, np.newaxis] - N_KD, 0.0)  # (K, D)
        weighted_sum += n_miss * old_centers
        return weighted_sum / N_K[:, np.newaxis]
    else:
        return weighted_sum / N_KD


def _estimate_laplace_centers(X, C, observed, K, D, N_K):
    """Weighted median center estimation for Laplace family."""
    centers = np.zeros((K, D))
    for k in range(K):
        weights = C[:, k]
        for d in range(D):
            obs_mask = observed[:, d]
            if obs_mask.sum() == 0:
                centers[k, d] = 0.0
                continue
            vals = X[obs_mask, d]
            w = weights[obs_mask]
            idx = np.argsort(vals)
            cumw = np.cumsum(w[idx])
            half = cumw[-1] / 2.0
            median_idx = min(np.searchsorted(cumw, half), len(vals) - 1)
            centers[k, d] = vals[idx[median_idx]]
    return centers


def _inertia_to_dispersions(Iner_KD, N_K, N_KD, N, D, K, model, miss_mode):
    """Convert inertia matrix to dispersions according to model."""
    dispersions = np.zeros((K, D))

    if model == Dispersion.S__:
        if miss_mode == "replace":
            v = Iner_KD.sum() / (N * D)
        else:
            v = Iner_KD.sum() / N_KD.sum()
        dispersions[:] = v

    elif model == Dispersion.SK_:
        for k in range(K):
            if miss_mode == "replace":
                vk = Iner_KD[k].sum() / (D * N_K[k])
            else:
                vk = Iner_KD[k].sum() / N_KD[k].sum()
            dispersions[k, :] = vk

    elif model == Dispersion.S_D:
        for d in range(D):
            if miss_mode == "replace":
                vd = Iner_KD[:, d].sum() / N
            else:
                vd = Iner_KD[:, d].sum() / N_KD[:, d].sum()
            dispersions[:, d] = vd

    elif model == Dispersion.SKD:
        if miss_mode == "replace":
            dispersions = Iner_KD / N_K[:, np.newaxis]
        else:
            dispersions = Iner_KD / N_KD

    return dispersions


def estimate_parameters(X, C, family, dispersion_model, proportion_model,
                        miss_mode="replace", old_centers=None,
                        old_dispersions=None):
    """M-step: estimate model parameters from soft classification.

    Parameters
    ----------
    X : (N, D) array
    C : (N, K) array — soft classification
    family : Family enum
    dispersion_model : Dispersion enum
    proportion_model : Proportion enum
    miss_mode : str — 'replace' or 'ignore'
    old_centers : (K, D) array or None
    old_dispersions : (K, D) array or None

    Returns
    -------
    dict with keys 'centers', 'dispersions', 'proportions'.
    """
    N, D = X.shape
    K = C.shape[1]
    observed = ~np.isnan(X)  # (N, D)

    # Class sizes
    N_K = C.sum(axis=0)  # (K,)
    N_K = np.maximum(N_K, EPSILON)

    # Per-class, per-variable observed sizes
    N_KD = np.zeros((K, D))
    for k in range(K):
        N_KD[k] = (C[:, k:k+1] * observed).sum(axis=0)
    N_KD = np.maximum(N_KD, EPSILON)

    # --- Centers ---
    if family == Family.LAPLACE:
        centers = _estimate_laplace_centers(X, C, observed, K, D, N_K)
    else:
        centers = _estimate_mean_centers(X, C, observed, K, D, N, N_K, N_KD,
                                         miss_mode, old_centers)

    # --- Inertia ---
    Iner_KD = np.zeros((K, D))
    X_filled = X.copy()
    X_filled[~observed] = 0.0

    for k in range(K):
        diff = X_filled - centers[k]
        diff[~observed] = 0.0
        if family == Family.NORMAL:
            Iner_KD[k] = (C[:, k:k+1] * diff ** 2 * observed).sum(axis=0)
        else:  # Laplace or Bernoulli
            Iner_KD[k] = (C[:, k:k+1] * np.abs(diff) * observed).sum(axis=0)

        # Missing data correction for REPLACE mode (Normal only)
        if miss_mode == "replace" and family == Family.NORMAL and old_dispersions is not None:
            n_miss_kd = N_K[k] - N_KD[k]
            Iner_KD[k] += np.maximum(n_miss_kd, 0) * old_dispersions[k]

    # --- Dispersions ---
    dispersions = _inertia_to_dispersions(
        Iner_KD, N_K, N_KD, N, D, K, dispersion_model, miss_mode
    )

    # Clamp dispersions from below
    dispersions = np.maximum(dispersions, EPSILON)

    # --- Proportions ---
    if proportion_model == Proportion.EQUAL:
        proportions = np.full(K, 1.0 / K)
    else:
        proportions = N_K / N
        proportions = np.maximum(proportions, EPSILON)
        proportions /= proportions.sum()

    return {
        "centers": centers,
        "dispersions": dispersions,
        "proportions": proportions,
    }


def _estimate_mean_centers(X, C, observed, K, D, N, N_K, N_KD, miss_mode,
                           old_centers):
    """Weighted mean center estimation (Normal, Bernoulli)."""
    centers = np.zeros((K, D))
    X_filled = X.copy()
    X_filled[~observed] = 0.0

    for k in range(K):
        weighted_sum = (C[:, k:k+1] * X_filled * observed).sum(axis=0)
        if miss_mode == "replace" and old_centers is not None:
            n_miss_kd = N_K[k] - N_KD[k]
            weighted_sum += np.maximum(n_miss_kd, 0) * old_centers[k]
            centers[k] = weighted_sum / N_K[k]
        else:
            centers[k] = weighted_sum / N_KD[k]

    return centers


def _estimate_laplace_centers(X, C, observed, K, D, N_K):
    """Weighted median center estimation for Laplace family."""
    centers = np.zeros((K, D))
    for k in range(K):
        weights = C[:, k]
        for d in range(D):
            obs_mask = observed[:, d]
            if obs_mask.sum() == 0:
                centers[k, d] = 0.0
                continue
            vals = X[obs_mask, d]
            w = weights[obs_mask]
            # Weighted median
            idx = np.argsort(vals)
            vals_sorted = vals[idx]
            w_sorted = w[idx]
            cumw = np.cumsum(w_sorted)
            half = cumw[-1] / 2.0
            median_idx = np.searchsorted(cumw, half)
            median_idx = min(median_idx, len(vals_sorted) - 1)
            centers[k, d] = vals_sorted[median_idx]
    return centers


def _inertia_to_dispersions(Iner_KD, N_K, N_KD, N, D, K, model, miss_mode):
    """Convert inertia matrix to dispersions according to model."""
    dispersions = np.zeros((K, D))

    if model == Dispersion.S__:
        if miss_mode == "replace":
            v = Iner_KD.sum() / (N * D)
        else:
            v = Iner_KD.sum() / N_KD.sum()
        dispersions[:] = v

    elif model == Dispersion.SK_:
        for k in range(K):
            if miss_mode == "replace":
                vk = Iner_KD[k].sum() / (D * N_K[k])
            else:
                vk = Iner_KD[k].sum() / N_KD[k].sum()
            dispersions[k, :] = vk

    elif model == Dispersion.S_D:
        for d in range(D):
            if miss_mode == "replace":
                vd = Iner_KD[:, d].sum() / N
            else:
                vd = Iner_KD[:, d].sum() / N_KD[:, d].sum()
            dispersions[:, d] = vd

    elif model == Dispersion.SKD:
        if miss_mode == "replace":
            for k in range(K):
                dispersions[k] = Iner_KD[k] / N_K[k]
        else:
            dispersions = Iner_KD / N_KD

    return dispersions
