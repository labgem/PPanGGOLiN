#!/usr/bin/env python3
"""
Test: compare C NEM (nem_stats) vs Python NEM (pynem) on PPanGGOLiN-format input files.

Steps:
  1. Create synthetic binary input files (Bernoulli, similar to ppanggolin's gene-family matrix)
  2. Run C NEM via nem_stats (current ppanggolin backend)
  3. Run Python NEM via pynem
  4. Compare partition assignments with true labels and with each other

Multiple scenarios are tested, progressively harder:
  - Easy: wide separation between classes (easy to cluster)
  - Realistic: overlap between classes, graph from ppanggolin-like adjacency
  - Real ppanggolin NEM files (if a pangenome HDF5 is available)

Usage:
    conda activate ppanggo_dev
    python test_pynem_vs_cnem.py
"""

import math
import sys
import tempfile
import numpy as np
from pathlib import Path

# ── make pynem importable without installing ──────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "nem" / "pynem" / "src"))

import pynem
import nem_stats

# ppanggolin imports for real-data test
import ppanggolin.nem.partition as _partition_mod
from ppanggolin.nem.pynem_adapter import PPanGGOLiNEM
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import check_pangenome_info
from ppanggolin.formats.readBinaries import get_status


# ─────────────────────────────────────────────────────────────────────────────
# Step 0 – Real pangenome test (via ppanggolin's own reader + write_nem_input_files)
# ─────────────────────────────────────────────────────────────────────────────

def run_real_pangenome_test(h5_path: Path, beta: float = 2.5, kval: int = 3,
                            seed: int = 42, sm_degree: int = 10):
    """Load a ppanggolin HDF5, write NEM input files exactly as partition.py
    does, then compare C NEM vs Python NEM.
    """
    print(f"\n  Loading pangenome: {h5_path.name}")
    pangenome = Pangenome()
    pangenome.file = str(h5_path)
    get_status(pangenome, h5_path)
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=True,
        disable_bar=True,
    )
    organisms = set(pangenome.organisms)
    n_orgs = len(organisms)
    print(f"  {n_orgs} organisms | {pangenome.number_of_gene_families} gene families")

    # Set the global 'pan' that write_nem_input_files reads
    _partition_mod.pan = pangenome

    with tempfile.TemporaryDirectory() as tmpdir_str:
        tmpdir = Path(tmpdir_str) / "nem_files"  # subdir so mk_outdir can create it
        edges_weight, nb_fam = _partition_mod.write_nem_input_files(
            tmpdir=tmpdir, organisms=organisms, sm_degree=sm_degree
        )
        effective_beta = beta * (nb_fam / edges_weight) if edges_weight > 0 else beta
        print(f"  {nb_fam} families in NEM | effective beta={effective_beta:.4f}")

        print("\n  [C NEM]")
        c_labels = run_c_nem(tmpdir, n_orgs=n_orgs, kval=kval,
                             beta=effective_beta, seed=seed)
        if c_labels is None:
            print("  C NEM FAILED.")
            return
        c_labels = np.array(c_labels)
        print(f"  {label_counts(c_labels, kval)}")

        print("\n  [Python NEM]")
        py_labels = run_python_nem(tmpdir, kval=kval, beta=effective_beta,
                                   seed=seed, use_param_file_init=True)
        py_labels = np.array(py_labels)
        print(f"  {label_counts(py_labels, kval)}")

        print("\n  [Comparison]")
        ari_cp = pynem.metrics.adjusted_rand_index(c_labels, py_labels)
        print(f"  ARI(C NEM vs Python NEM) = {ari_cp:.4f}")

        from scipy.optimize import linear_sum_assignment
        confusion = np.zeros((kval, kval), dtype=int)
        for cl, pl in zip(c_labels, py_labels):
            confusion[cl - 1, pl - 1] += 1
        row_ind, col_ind = linear_sum_assignment(-confusion)
        n_agree = confusion[row_ind, col_ind].sum()
        print(f"  Label agreement (best perm): {n_agree}/{nb_fam} = "
              f"{n_agree / nb_fam:.1%}")

        # Also compare ppanggolin partition name assignments
        with open(tmpdir / "nem_file.index") as f:
            index_fam = [line.split("\t")[1].strip() for line in f]
        parti_map = {0: "P", kval - 1: "C"}
        for i in range(1, kval - 1):
            parti_map[i] = f"S{i}"

        # Read the C NEM .uf file directly to get ppanggolin-style assignments
        uf_path = tmpdir / f"nem_file_{kval}.uf"
        c_partitions = {}
        with open(uf_path) as f:
            for i, line in enumerate(f):
                probs = [float(v) for v in line.split()]
                max_p = max(probs)
                best = [j for j, p in enumerate(probs) if p == max_p]
                if len(best) > 1 or max_p < 0.5:
                    c_partitions[index_fam[i]] = "S_"
                else:
                    c_partitions[index_fam[i]] = parti_map[best[0]]

        # Python NEM assignments
        py_partitions = {}
        for i, label in enumerate(py_labels):
            k_idx = int(label) - 1
            py_partitions[index_fam[i]] = parti_map.get(k_idx, "S_")

        # Count where they agree on partition name
        name_agree = sum(1 for fam in index_fam
                         if c_partitions[fam] == py_partitions[fam])
        print(f"  Partition-name agreement: {name_agree}/{nb_fam} = "
              f"{name_agree / nb_fam:.1%}")

        # Show distribution per partition name
        from collections import Counter
        c_cnt = Counter(c_partitions.values())
        py_cnt = Counter(py_partitions.values())
        print(f"  C NEM    partitions: {dict(sorted(c_cnt.items()))}")
        print(f"  Python NEM partitions: {dict(sorted(py_cnt.items()))}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 – Create synthetic NEM input files (same format as ppanggolin)
# ─────────────────────────────────────────────────────────────────────────────

def create_synthetic_nem_files(
    tmpdir: Path,
    n_families: int = 150,
    n_orgs: int = 20,
    k: int = 3,
    probs: list = None,
    graph_type: str = "chain",
    seed: int = 42,
):
    """Write .str / .dat / .nei / .index files with k-class Bernoulli structure.

    Parameters
    ----------
    probs       : list of length k — presence probability per class.
                  Defaults to [0.90, 0.50, 0.10] for k=3.
    graph_type  : 'chain' (simple sequential) or 'random' (Erdős–Rényi-like).

    Returns
    -------
    true_labels : np.ndarray, shape (n_families,)  — 0-based true class
    X           : np.ndarray, shape (n_families, n_orgs)  — binary matrix
    """
    rng = np.random.default_rng(seed)
    if probs is None:
        # evenly spaced from high to low presence
        probs = [0.90 - 0.80 * i / (k - 1) for i in range(k)]

    per_class = n_families // k
    labels = np.array([i // per_class for i in range(n_families)])
    labels[per_class * k :] = k - 1  # remainder → last class

    X = np.zeros((n_families, n_orgs), dtype=int)
    for i in range(n_families):
        X[i] = rng.binomial(1, probs[labels[i]], n_orgs)

    tmpdir.mkdir(parents=True, exist_ok=True)

    # .str
    with open(tmpdir / "nem_file.str", "w") as f:
        f.write(f"S\t{n_families}\t{n_orgs}\n")

    # .dat  (tab-separated, same as ppanggolin)
    with open(tmpdir / "nem_file.dat", "w") as f:
        for row in X:
            f.write("\t".join(map(str, row)) + "\n")

    # .nei  (weighted, 1-based, same format as ppanggolin)
    with open(tmpdir / "nem_file.nei", "w") as f:
        f.write("1\n")  # weighted flag

        if graph_type == "chain":
            # Simple chain: each node connected to previous and next
            for i in range(n_families):
                neighbors, weights = [], []
                if i > 0:
                    neighbors.append(i)       # 1-based
                    weights.append(1.0)
                if i < n_families - 1:
                    neighbors.append(i + 2)   # 1-based
                    weights.append(1.0)
                if neighbors:
                    f.write(
                        f"{i + 1}\t{len(neighbors)}\t"
                        + "\t".join(map(str, neighbors))
                        + "\t"
                        + "\t".join(map(str, weights))
                        + "\n"
                    )
                else:
                    f.write(f"{i + 1}\t0\n")
        else:
            # Random sparse graph (Erdős–Rényi, p ≈ 4/n_families)
            p_edge = 4.0 / max(n_families - 1, 1)
            adj = {i: [] for i in range(n_families)}
            for i in range(n_families):
                for j in range(i + 1, n_families):
                    if rng.random() < p_edge:
                        w = round(float(rng.uniform(0.1, 1.0)), 4)
                        adj[i].append((j, w))
                        adj[j].append((i, w))
            for i in range(n_families):
                neighbors = [j + 1 for j, _ in adj[i]]
                weights   = [w     for _, w in adj[i]]
                if neighbors:
                    f.write(
                        f"{i + 1}\t{len(neighbors)}\t"
                        + "\t".join(map(str, neighbors))
                        + "\t"
                        + "\t".join(map(str, weights))
                        + "\n"
                    )
                else:
                    f.write(f"{i + 1}\t0\n")

    # .index
    with open(tmpdir / "nem_file.index", "w") as f:
        for i in range(n_families):
            f.write(f"{i + 1}\tfam_{i + 1}\n")

    print(f"  {n_families} families | {n_orgs} organisms | {k} true classes | graph={graph_type}")
    for ci in range(k):
        print(f"    class {ci}: {(labels == ci).sum()} families  "
              f"(p_presence ≈ {probs[ci]:.0%})")

    return labels, X


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 – Run C NEM (nem_stats, current ppanggolin backend)
# ─────────────────────────────────────────────────────────────────────────────

def write_init_file(nem_dir: Path, kval: int, n_orgs: int):
    """Write the .m parameter-initialisation file, exactly as partition.py does."""
    with open(nem_dir / f"nem_file_init_{kval}.m", "w") as f:
        f.write("1 ")
        base_prop = format(1.0 / float(kval), ".12g")
        f.write(" ".join([base_prop] * (kval - 1)) + " ")

        mu, epsilon = [], []
        step = 0.5 / math.ceil(kval / 2)
        pichenette = 0.1 if kval == 2 else 0.0
        for ki in range(1, kval + 1):
            if ki <= kval / 2:
                mu += ["1"] * n_orgs
                epsilon += [str((step * ki) - pichenette)] * n_orgs
            else:
                mu += ["0"] * n_orgs
                epsilon += [str((step * (kval - ki + 1)) - pichenette)] * n_orgs

        f.write(" ".join(mu) + " " + " ".join(epsilon))


def run_c_nem(
    nem_dir: Path,
    n_orgs: int,
    kval: int = 3,
    beta: float = 2.5,
    seed: int = 42,
    free_dispersion: bool = False,
):
    """Run C NEM and return 1-based class labels (list of int), or None on failure."""
    write_init_file(nem_dir, kval, n_orgs)

    _init_random, init_param_file = range(1, 3)
    variance_model = b"skd" if free_dispersion else b"sk_"

    nem_stats.nem(
        Fname=(nem_dir / "nem_file").as_posix().encode("ascii"),
        nk=kval,
        algo=b"nem",
        beta=beta,
        convergence=b"clas",
        convergence_th=0.01,
        format=b"fuzzy",
        it_max=100,
        dolog=True,
        model_family=b"bern",
        proportion=b"pk",
        dispersion=variance_model,
        init_mode=init_param_file,
        init_file=(nem_dir / f"nem_file_init_{kval}.m").as_posix().encode("ascii"),
        out_file_prefix=(nem_dir / f"nem_file_{kval}").as_posix().encode("ascii"),
        seed=seed,
    )

    uf_path = nem_dir / f"nem_file_{kval}.uf"
    if not uf_path.exists():
        print(f"  ⚠  C NEM did not produce {uf_path.name}")
        return None

    labels = []
    with open(uf_path) as f:
        for line in f:
            probs = [float(v) for v in line.split()]
            labels.append(int(np.argmax(probs)) + 1)  # 1-based
    return labels


# ─────────────────────────────────────────────────────────────────────────────
# Step 3 – Run Python NEM (pynem)
# ─────────────────────────────────────────────────────────────────────────────

def ppanggolin_init_params(kval: int, n_orgs: int, free_dispersion: bool = False):
    """Compute initial centers and dispersions exactly as ppanggolin's partition.py does.

    Returns
    -------
    centers     : np.ndarray (K, N_org)
    dispersions : np.ndarray (K, N_org)
    proportions : np.ndarray (K,)
    """
    step = 0.5 / math.ceil(kval / 2)
    pichenette = 0.1 if kval == 2 else 0.0

    centers = np.zeros((kval, n_orgs))
    dispersions = np.zeros((kval, n_orgs))
    for k in range(1, kval + 1):
        if k <= kval / 2:
            centers[k - 1] = 1.0
            dispersions[k - 1] = (step * k) - pichenette
        else:
            centers[k - 1] = 0.0
            dispersions[k - 1] = (step * (kval - k + 1)) - pichenette

    proportions = np.full(kval, 1.0 / kval)
    return centers, dispersions, proportions


def run_python_nem(
    nem_dir: Path,
    kval: int = 3,
    beta: float = 2.5,
    seed: int = 42,
    free_dispersion: bool = False,
    use_param_file_init: bool = True,
):
    """Run pynem on the .str/.dat/.nei files and return 1-based labels."""
    G = pynem.io.read_graph(nem_dir / "nem_file")
    n_orgs = G.graph["d"]

    dispersion = "skd" if free_dispersion else "sk_"

    if use_param_file_init:
        centers, dispersions, proportions = ppanggolin_init_params(kval, n_orgs, free_dispersion)
        init_mode = "param_file"
        init_kwargs = dict(
            init_centers=centers,
            init_dispersions=dispersions,
            init_proportions=proportions,
        )
    else:
        init_mode = "sort"
        init_kwargs = {}

    model = PPanGGOLiNEM(
        n_clusters=kval,
        beta=beta,
        family="bernoulli",
        dispersion=dispersion,
        proportion="pk",
        init=init_mode,
        max_iter=100,
        tol=0.01,
        convergence="classification",
        random_state=seed,
        verbose=0,
        **init_kwargs,
    )
    model.fit(G)

    print(f"  converged in {model.n_iter_} iterations")
    print(f"  U={model.criteria_['U']:.4f}  D={model.criteria_['D']:.4f}  "
          f"G={model.criteria_['G']:.4f}")

    return model.labels_  # 1-based numpy array


# ─────────────────────────────────────────────────────────────────────────────
# Step 4 – Compare
# ─────────────────────────────────────────────────────────────────────────────

def label_counts(labels, k):
    counts = np.bincount(np.asarray(labels), minlength=k + 1)[1:]  # skip 0
    return "  ".join(f"class {i+1}: {counts[i]}" for i in range(k))


def main():
    K = 3
    SEED = 42

    scenarios = [
        {
            "name": "Easy — well-separated, chain graph",
            "n_fam": 150, "n_org": 20, "beta": 2.5,
            "probs": [0.90, 0.50, 0.10], "graph": "chain",
        },
        {
            "name": "Harder — overlapping classes, chain graph",
            "n_fam": 150, "n_org": 20, "beta": 2.5,
            "probs": [0.80, 0.60, 0.30], "graph": "chain",
        },
        {
            "name": "Realistic — many organisms, random graph",
            "n_fam": 300, "n_org": 50, "beta": 2.5,
            "probs": [0.90, 0.50, 0.10], "graph": "random",
        },
        {
            "name": "No spatial smoothing (beta=0) — pure EM",
            "n_fam": 150, "n_org": 20, "beta": 0.0,
            "probs": [0.90, 0.50, 0.10], "graph": "chain",
        },
    ]

    for sc in scenarios:
        print("\n" + "═" * 65)
        print(f"Scenario: {sc['name']}")
        print("═" * 65)

        with tempfile.TemporaryDirectory() as tmpdir_str:
            tmpdir = Path(tmpdir_str)

            print("\n[1] Create synthetic NEM files")
            true_labels, _ = create_synthetic_nem_files(
                tmpdir,
                n_families=sc["n_fam"],
                n_orgs=sc["n_org"],
                k=K,
                probs=sc["probs"],
                graph_type=sc["graph"],
                seed=SEED,
            )

            print("\n[2] C NEM  (nem_stats)")
            c_labels = run_c_nem(
                tmpdir, n_orgs=sc["n_org"], kval=K, beta=sc["beta"], seed=SEED
            )
            if c_labels is None:
                print("  C NEM FAILED — skipping scenario.")
                continue
            c_labels = np.array(c_labels)
            print(f"  {label_counts(c_labels, K)}")

            print("\n[3] Python NEM  (pynem, param_file init)")
            py_labels = run_python_nem(
                tmpdir, kval=K, beta=sc["beta"], seed=SEED, use_param_file_init=True
            )
            py_labels = np.array(py_labels)
            print(f"  {label_counts(py_labels, K)}")

            print("\n[4] Comparison")
            true_1based = true_labels + 1
            ari_c  = pynem.metrics.adjusted_rand_index(true_1based, c_labels)
            ari_py = pynem.metrics.adjusted_rand_index(true_1based, py_labels)
            ari_cp = pynem.metrics.adjusted_rand_index(c_labels, py_labels)

            print(f"  ARI(C NEM      vs true) = {ari_c:.4f}")
            print(f"  ARI(Python NEM vs true) = {ari_py:.4f}")
            print(f"  ARI(C NEM vs Python NEM) = {ari_cp:.4f}")

            from scipy.optimize import linear_sum_assignment
            confusion = np.zeros((K, K), dtype=int)
            for cl, pl in zip(c_labels, py_labels):
                confusion[cl - 1, pl - 1] += 1
            row_ind, col_ind = linear_sum_assignment(-confusion)
            n_agree = confusion[row_ind, col_ind].sum()
            print(f"  Label agreement (best perm): {n_agree}/{sc['n_fam']} = "
                  f"{n_agree / sc['n_fam']:.1%}")

    print("\n" + "═" * 65)
    print("All scenarios complete.")

    # ── Real pangenome test ──────────────────────────────────────────────────
    real_h5 = Path(__file__).parent / "test_pangenome" / \
              "GTDB_refseq_s__Acinetobacter_schindleri_id54.h5"
    if real_h5.exists():
        print("\n" + "═" * 65)
        print("Scenario: Real pangenome (Acinetobacter schindleri, GTDB_refseq)")
        print("═" * 65)
        run_real_pangenome_test(real_h5, beta=2.5, kval=3, seed=42)
    else:
        print(f"\n[SKIP] Real pangenome file not found: {real_h5}")
        print("       Download it with:")
        print("       pangbank search-pangenomes --outdir test_pangenome "
              "-c GTDB_refseq -t 's__Acinetobacter schindleri' --download")


if __name__ == "__main__":
    main()
