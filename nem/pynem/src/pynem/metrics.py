"""Evaluation metrics for comparing clusterings."""

import numpy as np
from scipy.special import comb


def adjusted_rand_index(labels_true, labels_pred):
    """Compute the Adjusted Rand Index (ARI) between two clusterings.

    The ARI is a measure of agreement between two partitions, adjusted for
    chance. It ranges from -1 to 1, where 1 means perfect agreement and 0
    means agreement no better than random.

    Parameters
    ----------
    labels_true : array-like of shape (N,)
        Ground truth labels.
    labels_pred : array-like of shape (N,)
        Predicted labels.

    Returns
    -------
    ari : float
        Adjusted Rand Index.
    """
    labels_true = np.asarray(labels_true).ravel()
    labels_pred = np.asarray(labels_pred).ravel()
    assert len(labels_true) == len(labels_pred), "Label arrays must have same length"

    # Build contingency table
    classes_true = np.unique(labels_true)
    classes_pred = np.unique(labels_pred)
    n = len(labels_true)

    # Map labels to indices
    true_map = {c: i for i, c in enumerate(classes_true)}
    pred_map = {c: i for i, c in enumerate(classes_pred)}

    contingency = np.zeros((len(classes_true), len(classes_pred)), dtype=np.int64)
    for t, p in zip(labels_true, labels_pred):
        contingency[true_map[t], pred_map[p]] += 1

    # Sum of C(n_ij, 2) over all cells
    sum_comb_c = sum(comb(contingency[i, j], 2)
                     for i in range(contingency.shape[0])
                     for j in range(contingency.shape[1]))

    # Row and column sums
    row_sums = contingency.sum(axis=1)
    col_sums = contingency.sum(axis=0)

    sum_comb_rows = sum(comb(r, 2) for r in row_sums)
    sum_comb_cols = sum(comb(c, 2) for c in col_sums)

    total_comb = comb(n, 2)

    # ARI formula
    expected = sum_comb_rows * sum_comb_cols / total_comb
    max_index = 0.5 * (sum_comb_rows + sum_comb_cols)

    if max_index == expected:
        return 1.0

    ari = (sum_comb_c - expected) / (max_index - expected)
    return float(ari)
