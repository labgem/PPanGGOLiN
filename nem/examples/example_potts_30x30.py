"""Example: Potts image with 30x30 grid, 2D features.

Generates synthetic data on a 2D grid via Gibbs sampling (Potts model),
runs NEM, computes Adjusted Rand Index, and saves visualization figures
for the README.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Add parent dir so pynem is importable without install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pynem" / "src"))

import pynem
from generate import PottsImageData

# ── Generate synthetic data ──────────────────────────────────────────────

potts = PottsImageData(
    nl=30,
    nc=30,
    k=3,
    beta=0.8,
    centers=[[0, 0], [3, 3], [0, 3]],
    sigma=1.0,
    seed=42,
)
potts.export(str(Path(__file__).parent / "potts_30x30_3"))

# ── Load and fit NEM ─────────────────────────────────────────────────────

G = pynem.io.read_graph(str(Path(__file__).parent / "potts_30x30_3"))
model = pynem.NEM(n_clusters=3, beta=1.0, family="normal")
model.fit(G)

# ── Compute Adjusted Rand Index ──────────────────────────────────────────

true_labels_1based = potts.labels + 1  # convert 0-based to 1-based
ari = pynem.metrics.adjusted_rand_index(true_labels_1based, model.labels_)
print(f"Adjusted Rand Index: {ari:.4f}")

# ── Visualization: results overview ──────────────────────────────────────

fig = pynem.viz.plot_results(
    G, model,
    true_labels=true_labels_1based,
    var_names=["Variable 1", "Variable 2"],
)
fig.suptitle(f"NEM on Potts (30×30 grid, 3 classes) — ARI = {ari:.2f}",
             fontsize=14, y=1.02)
fig.savefig(str(Path(__file__).parent / "potts_30x30_3_results.png"),
            dpi=150, bbox_inches="tight")
print("Saved: potts_30x30_3_results.png")

# ── Visualization: convergence ───────────────────────────────────────────

fig2, ax2 = plt.subplots(figsize=(6, 4))
pynem.viz.plot_convergence(model.history_, ax=ax2)
fig2.savefig(str(Path(__file__).parent / "potts_30x30_3_convergence.png"),
             dpi=150, bbox_inches="tight")
print("Saved: potts_30x30_3_convergence.png")

plt.show()
