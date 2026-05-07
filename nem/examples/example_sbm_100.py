"""Example: SBM graph with 100 nodes, 2D features.

Generates synthetic data, runs NEM, computes Adjusted Rand Index,
and saves visualization figures for the README.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Add parent dir so pynem is importable without install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pynem" / "src"))


import pynem
from generate import SBMData

# ── Generate synthetic data ──────────────────────────────────────────────

sbm = SBMData(
    n=100, k=3, d=2,
    p_in=0.2, p_out=0.02,
    centers=[[0, 0], [3, 3], [0, 3]],
    sigma=1.0, seed=42,
)
sbm.export(str(Path(__file__).parent / "sbm_100_2"))

# ── Load and fit NEM ─────────────────────────────────────────────────────

G = pynem.io.read_graph(str(Path(__file__).parent / "sbm_100_2"))
model = pynem.NEM(n_clusters=3, beta=1.0, family="normal")
model.fit(G)

# ── Compute Adjusted Rand Index ──────────────────────────────────────────

true_labels_1based = sbm.labels + 1  # convert 0-based to 1-based
ari = pynem.metrics.adjusted_rand_index(true_labels_1based, model.labels_)
print(f"Adjusted Rand Index: {ari:.4f}")

# ── Visualization: results overview ──────────────────────────────────────

fig = pynem.viz.plot_results(
    G, model,
    true_labels=true_labels_1based,
    var_names=["Variable 1", "Variable 2"],
)
fig.suptitle(f"NEM on SBM (100 nodes, 3 classes) — ARI = {ari:.2f}",
             fontsize=14, y=1.02)
fig.savefig(str(Path(__file__).parent / "sbm_100_2_results.png"),
            dpi=150, bbox_inches="tight")
print(f"Saved: sbm_100_2_results.png")

# ── Visualization: convergence ───────────────────────────────────────────

fig2, ax2 = plt.subplots(figsize=(6, 4))
pynem.viz.plot_convergence(model.history_, ax=ax2)
fig2.savefig(str(Path(__file__).parent / "sbm_100_2_convergence.png"),
             dpi=150, bbox_inches="tight")
print(f"Saved: sbm_100_2_convergence.png")

plt.show()
