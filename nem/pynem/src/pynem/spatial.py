"""Graph/neighborhood handling for NEM."""

import numpy as np
import networkx as nx
import scipy.sparse


class NeighborhoodSystem:
    """Adjacency structure extracted from a networkx graph.

    Parameters
    ----------
    G : nx.Graph
        Graph with optional 'weight' edge attribute (default 1.0).

    Notes
    -----
    A scipy CSR sparse matrix is built in ``__init__`` so that
    ``compute_all_contexts`` reduces to a single sparse BLAS call
    instead of a Python-level double loop.
    """

    def __init__(self, G):
        self._n = G.number_of_nodes()
        # Build adjacency lists (kept for individual neighbor lookups)
        self._neighbors = []
        for i in range(self._n):
            neighs = [(j, G.edges[i, j].get("weight", 1.0)) for j in G.neighbors(i)]
            self._neighbors.append(neighs)

        # Pre-build sparse adjacency matrix for fast batch operations (W[i,j] = w_ij)
        rows, cols, data = [], [], []
        for i, neighbors in enumerate(self._neighbors):
            for j, w in neighbors:
                rows.append(i)
                cols.append(j)
                data.append(w)
        self._W = scipy.sparse.csr_matrix(
            (data, (rows, cols)), shape=(self._n, self._n), dtype=np.float64
        )

    @property
    def n_nodes(self):
        return self._n

    @property
    def max_neighbors(self):
        return int(self._W.getnnz(axis=1).max()) if self._W.nnz > 0 else 0

    def neighbors(self, i):
        """Return list of (index, weight) for neighbors of node i."""
        return self._neighbors[i]

    def spatial_context(self, i, C, K):
        """Compute spatial context for node i.

        Returns
        -------
        context : array of shape (K,)
            context[k] = sum_{j in N(i)} w_ij * c_jk
        """
        context = np.zeros(K)
        for j, w in self._neighbors[i]:
            context += w * C[j]
        return context

    def compute_all_contexts(self, C):
        """Compute spatial context for all nodes at once via sparse matmul.

        Equivalent to ``contexts[i, k] = sum_{j in N(i)} w_ij * c_jk`` but
        implemented as ``W @ C`` using a pre-built CSR sparse matrix so that
        no Python-level loop is needed.

        Parameters
        ----------
        C : (N, K) array — classification matrix.

        Returns
        -------
        contexts : (N, K) array
        """
        return self._W @ C

    def compute_G(self, C):
        """Compute geographic cohesion: G = sum_i c_i · context_i."""
        contexts = self._W @ C
        return float((C * contexts).sum())
