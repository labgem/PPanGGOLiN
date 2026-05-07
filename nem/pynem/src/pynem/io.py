"""Read/write NEM file formats and convert to/from networkx graphs."""

import numpy as np
import networkx as nx
from pathlib import Path


def read_str(path):
    """Parse a .str file.

    Returns
    -------
    type_char : str
        'S' (spatial graph), 'I' (image grid), or 'N' (non-spatial).
    n : int
        Number of nodes (nl*nc for images).
    d : int
        Number of variables per node.
    """
    with open(path) as f:
        parts = f.read().split()
    type_char = parts[0]
    if type_char == "I":
        # Image format: I nl nc d
        nl, nc, d = int(parts[1]), int(parts[2]), int(parts[3])
        return type_char, nl * nc, d
    else:
        # Spatial/non-spatial: S n d or N n d
        return type_char, int(parts[1]), int(parts[2])


def read_dat(path, n, d):
    """Parse a .dat file.

    Returns
    -------
    X : np.ndarray of shape (n, d)
        Feature matrix. Missing values are NaN.
    """
    X = np.loadtxt(path).reshape(n, d)
    return X


def read_nei(path, n):
    """Parse a .nei file.

    Parameters
    ----------
    path : str or Path
        Path to the .nei file.
    n : int
        Number of nodes (from .str file).

    Returns
    -------
    G : nx.Graph
        Graph with nodes 0..n-1. Edges have 'weight' attribute (default 1.0).
    """
    G = nx.Graph()
    G.add_nodes_from(range(n))

    with open(path) as f:
        # First line: 0 = unweighted, 1 = weighted
        weighted = int(f.readline().strip())

        for line in f:
            tokens = line.split()
            if len(tokens) < 2:
                continue
            node_id = int(tokens[0]) - 1  # C code uses 1-based indexing
            nb_neigh = int(tokens[1])
            neighbors = [int(tokens[2 + j]) - 1 for j in range(nb_neigh)]

            if weighted:
                # Weights follow the neighbor list on the same line
                weights = [float(tokens[2 + nb_neigh + j]) for j in range(nb_neigh)]
                for j, neigh in enumerate(neighbors):
                    if 0 <= neigh < n and weights[j] != 0.0:
                        G.add_edge(node_id, neigh, weight=weights[j])
            else:
                for neigh in neighbors:
                    if 0 <= neigh < n:
                        G.add_edge(node_id, neigh, weight=1.0)

    return G


def read_graph(basename):
    """Load a complete NEM dataset from basename.{str,dat,nei}.

    Parameters
    ----------
    basename : str or Path
        Base path (e.g. 'examples/essai').

    Returns
    -------
    G : nx.Graph
        Graph where ``G.nodes[i]['features']`` is a numpy array of shape (D,).
    """
    basename = Path(basename)
    str_path = str(basename) + ".str"
    type_char, n, d = read_str(str_path)
    X = read_dat(str(basename) + ".dat", n, d)
    G = read_nei(str(basename) + ".nei", n)

    for i in range(n):
        G.nodes[i]["features"] = X[i]

    G.graph["type"] = type_char
    G.graph["n"] = n
    G.graph["d"] = d

    # Store image dimensions for visualization
    if type_char == "I":
        with open(str_path) as f:
            parts = f.read().split()
        G.graph["nl"] = int(parts[1])
        G.graph["nc"] = int(parts[2])

    return G


def write_classification(path, labels):
    """Write hard classification to a .cf file.

    Parameters
    ----------
    path : str or Path
        Output file path.
    labels : array-like of shape (N,)
        Class labels (1-based).
    """
    labels = np.asarray(labels)
    with open(path, "w") as f:
        for label in labels:
            f.write(f"{label}\n")


def write_fuzzy(path, membership):
    """Write fuzzy classification to a .uf file.

    Parameters
    ----------
    path : str or Path
        Output file path.
    membership : np.ndarray of shape (N, K)
        Membership probabilities.
    """
    with open(path, "w") as f:
        for row in membership:
            f.write(" ".join(f"{v:.6f}" for v in row) + "\n")
