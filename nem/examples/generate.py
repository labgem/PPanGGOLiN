"""Synthetic data generators for testing NEM.

Provides two data classes:
- SBMData: Stochastic Block Model graph with Gaussian emissions
- PottsImageData: Potts/Strauss model on a grid with Gaussian emissions

Each class generates labeled graph data and can export to NEM file formats
(.str, .dat, .nei).

Usage:
    python generate.py          # generates example datasets in examples/
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path


# ---------------------------------------------------------------------------
# Export functions
# ---------------------------------------------------------------------------

def write_str(path, type_char, *args):
    """Write a .str file.

    Parameters
    ----------
    path : str or Path
    type_char : str
        'S' for spatial graph (args: n, d), 'I' for image grid (args: nl, nc, d).
    *args : int
        For 'S': (n, d). For 'I': (nl, nc, d).
    """
    with open(path, "w") as f:
        f.write(type_char + " " + " ".join(str(a) for a in args) + "\n")


def write_dat(path, X):
    """Write a .dat file (N x D feature matrix).

    Parameters
    ----------
    path : str or Path
    X : np.ndarray of shape (N, D)
    """
    np.savetxt(path, X, fmt="%.6f", delimiter="\t")


def write_labels(path, labels):
    """Write true labels to a .cf file (1-based, space-separated).

    Parameters
    ----------
    path : str or Path
    labels : np.ndarray of shape (N,)
        0-based class labels.
    """
    labels_1based = np.asarray(labels) + 1
    with open(path, "w") as f:
        f.write(" ".join(str(l) for l in labels_1based) + "\n")


def write_nei(path, G):
    """Write a .nei file from a networkx Graph.

    Uses 1-based node indexing and unweighted format, matching the C code.

    Parameters
    ----------
    path : str or Path
    G : nx.Graph
        Nodes must be integers 0..N-1.
    """
    n = G.number_of_nodes()
    with open(path, "w") as f:
        f.write("0\n")  # unweighted
        for i in range(n):
            neighbors = sorted(G.neighbors(i))
            nb = len(neighbors)
            neigh_str = " ".join(str(j + 1) for j in neighbors)
            f.write(f"{i + 1} {nb} {neigh_str}\n")


# ---------------------------------------------------------------------------
# SBMData
# ---------------------------------------------------------------------------

class SBMData:
    """Stochastic Block Model graph with Gaussian emissions.

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Number of classes.
    d : int
        Dimension of feature vectors.
    p_in : float
        Intra-class edge probability.
    p_out : float
        Inter-class edge probability.
    centers : array-like of shape (k, d)
        Class centers for Gaussian emissions.
    sigma : float
        Standard deviation for Gaussian emissions (isotropic).
    seed : int or None
        Random seed.
    """

    def __init__(self, n, k, d, p_in, p_out, centers, sigma=1.0, seed=None):
        self.rng = np.random.default_rng(seed)
        centers = np.asarray(centers, dtype=float)
        assert centers.shape == (k, d), f"centers must be ({k}, {d})"

        # Balanced partition
        sizes = [n // k] * k
        for i in range(n % k):
            sizes[i] += 1

        # Edge probability matrix
        probs = np.full((k, k), p_out)
        np.fill_diagonal(probs, p_in)

        # Generate SBM graph
        self.graph = nx.stochastic_block_model(
            sizes, probs.tolist(), seed=int(self.rng.integers(2**31))
        )
        # Remove the 'block' attribute added by networkx, keep only our labels
        # Build ground-truth labels (0-based)
        self.labels = np.array(
            [self.graph.nodes[i]["block"] for i in range(n)]
        )

        # Generate features
        self.features = np.zeros((n, d))
        for i in range(n):
            c = self.labels[i]
            self.features[i] = self.rng.normal(centers[c], sigma, size=d)
            self.graph.nodes[i]["features"] = self.features[i]

        self.centers = centers.copy()
        self.n = n
        self.k = k
        self.d = d

    def plot(self, ax=None):
        """Draw the graph colored by true class labels."""
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 6))
        pos = nx.spring_layout(self.graph, seed=42)
        cmap = plt.cm.get_cmap("tab10", self.k)
        nx.draw_networkx_edges(self.graph, pos, ax=ax, alpha=0.15)
        nx.draw_networkx_nodes(
            self.graph, pos, ax=ax,
            node_color=self.labels,
            cmap=cmap, vmin=0, vmax=self.k - 1,
            node_size=60, alpha=0.9,
        )
        ax.set_title(f"SBM — {self.n} nodes, {self.k} classes")
        ax.axis("off")
        return ax

    def export(self, basename):
        """Write basename.{str,dat,nei,true.cf} files."""
        write_str(f"{basename}.str", "S", self.n, self.d)
        write_dat(f"{basename}.dat", self.features)
        write_nei(f"{basename}.nei", self.graph)
        write_labels(f"{basename}.true.cf", self.labels)
        print(f"Exported: {basename}.{{str,dat,nei,true.cf}}")


# ---------------------------------------------------------------------------
# PottsImageData
# ---------------------------------------------------------------------------

class PottsImageData:
    """Potts model on a 2D grid with Gaussian emissions.

    Generates spatially coherent labels via Gibbs sampling on a Potts MRF,
    then draws Gaussian features for each pixel.

    Parameters
    ----------
    nl : int
        Number of rows (lines).
    nc : int
        Number of columns.
    k : int
        Number of classes.
    beta : float
        Potts interaction parameter (higher = more spatial coherence).
    centers : array-like of shape (k, d)
        Class centers for Gaussian emissions.
    sigma : float
        Standard deviation for Gaussian emissions.
    seed : int or None
        Random seed.
    n_gibbs : int
        Number of Gibbs sweeps for label simulation.
    """

    def __init__(self, nl, nc, k, beta, centers, sigma=1.0, seed=None,
                 n_gibbs=200):
        self.rng = np.random.default_rng(seed)
        centers = np.asarray(centers, dtype=float)
        d = centers.shape[1]
        assert centers.shape == (k, d), f"centers must be ({k}, {d})"

        self.nl = nl
        self.nc = nc
        self.k = k
        self.d = d
        self.centers = centers.copy()
        n = nl * nc

        # Initialize labels randomly
        labels = self.rng.integers(0, k, size=(nl, nc))

        # Gibbs sampling on Potts model (4-connectivity)
        for sweep in range(n_gibbs):
            for r in range(nl):
                for c in range(nc):
                    # Count neighbor labels
                    counts = np.zeros(k)
                    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        rr, cc = r + dr, c + dc
                        if 0 <= rr < nl and 0 <= cc < nc:
                            counts[labels[rr, cc]] += 1
                    # Potts energy: exp(beta * count_of_same_class)
                    logprob = beta * counts
                    logprob -= logprob.max()  # numerical stability
                    prob = np.exp(logprob)
                    prob /= prob.sum()
                    labels[r, c] = self.rng.choice(k, p=prob)

        self.labels_image = labels
        self.labels = labels.ravel()

        # Generate features
        self.features = np.zeros((n, d))
        for i in range(n):
            c = self.labels[i]
            self.features[i] = self.rng.normal(centers[c], sigma, size=d)

        # Build grid graph
        self.graph = nx.grid_2d_graph(nl, nc)
        # Relabel to integer nodes 0..N-1 (row-major order)
        mapping = {(r, c): r * nc + c for r in range(nl) for c in range(nc)}
        self.graph = nx.relabel_nodes(self.graph, mapping)
        for i in range(n):
            self.graph.nodes[i]["features"] = self.features[i]

    def plot(self, figsize=(12, 5)):
        """Display label map and first feature channel (or RGB if D=3)."""
        fig, axes = plt.subplots(1, 2, figsize=figsize)

        # Label map
        cmap = plt.cm.get_cmap("tab10", self.k)
        axes[0].imshow(self.labels_image, cmap=cmap, vmin=0,
                       vmax=self.k - 1, interpolation="nearest")
        axes[0].set_title(f"Labels — {self.nl}x{self.nc}, {self.k} classes")
        axes[0].axis("off")

        # Feature image
        feat_img = self.features.reshape(self.nl, self.nc, self.d)
        if self.d == 3:
            # Normalize to [0, 1] for RGB display
            fmin = feat_img.min(axis=(0, 1), keepdims=True)
            fmax = feat_img.max(axis=(0, 1), keepdims=True)
            rgb = (feat_img - fmin) / (fmax - fmin + 1e-8)
            axes[1].imshow(rgb, interpolation="nearest")
            axes[1].set_title("Features (RGB)")
        else:
            axes[1].imshow(feat_img[:, :, 0], cmap="viridis",
                           interpolation="nearest")
            axes[1].set_title("Features (channel 0)")
        axes[1].axis("off")

        plt.tight_layout()
        return fig, axes

    def export(self, basename):
        """Write basename.{str,dat,nei,true.cf} files.

        Uses image format ('I nl nc d') for .str and writes the 4-connectivity
        neighborhood for compatibility with spatial mode.
        """
        write_str(f"{basename}.str", "I", self.nl, self.nc, self.d)
        write_dat(f"{basename}.dat", self.features)
        write_nei(f"{basename}.nei", self.graph)
        write_labels(f"{basename}.true.cf", self.labels)
        print(f"Exported: {basename}.{{str,dat,nei,true.cf}}")


# ---------------------------------------------------------------------------
# Main: generate example datasets
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    script_dir = Path(__file__).parent

    # --- SBM example: 100 nodes, 3 classes, 3D features ---
    print("Generating SBM dataset (100 nodes, 3 classes)...")
    sbm = SBMData(
        n=100, k=3, d=3,
        p_in=0.3, p_out=0.02,
        centers=[[0, 0, 0], [4, 4, 4], [0, 4, 0]],
        sigma=1.0, seed=42,
    )
    sbm.export(str(script_dir / "sbm_100_3"))

    # --- Potts image: 20x25 grid, 3 classes, 3D features ---
    print("Generating Potts image dataset (20x25, 3 classes)...")
    potts = PottsImageData(
        nl=20, nc=25, k=3, beta=1.5,
        centers=[[0, 0, 0], [4, 4, 4], [0, 4, 0]],
        sigma=1.0, seed=42,
    )
    potts.export(str(script_dir / "potts_20x25_3"))

    print("Done.")
