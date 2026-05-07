"""Diagnose why C NEM and pynem disagree on ~175 families."""
import sys, math, tempfile, shutil
from pathlib import Path
import numpy as np
from collections import Counter

sys.path.insert(0, ".")
import pynem as _pynem
from ppanggolin.nem.pynem_adapter import PPanGGOLiNEM
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.readBinaries import get_status
from ppanggolin.formats import check_pangenome_info
import ppanggolin.nem.partition as _part
from ppanggolin.nem.partition import write_nem_input_files, run_partitioning

h5 = Path("test_pangenome/GTDB_refseq_s__Acinetobacter_schindleri_id54.h5")
pan = Pangenome()
pan.file = str(h5)
get_status(pan, h5)
check_pangenome_info(pan, need_families=True, need_annotations=True,
                     need_graph=True, disable_bar=True)
_part.pan = pan

kval = 3
nb_org = pan.number_of_organisms
beta_raw = 2.5

td_obj = tempfile.TemporaryDirectory()
td = td_obj.name
try:
    nem_dir = Path(td) / "nem_files"
    edges_weight, nb_fam = write_nem_input_files(nem_dir, set(pan.organisms), 10)
    beta = beta_raw * (nb_fam / edges_weight)

    # --- C NEM ---
    parts_c, params_c, _ = run_partitioning(
        nem_dir, nb_org, beta=beta, kval=kval, seed=42, use_pynem=False)

    # --- pynem ---
    G = _pynem.io.read_graph(nem_dir / "nem_file")
    step = 0.5 / math.ceil(kval / 2)
    centers = np.zeros((kval, nb_org))
    dispersions = np.zeros((kval, nb_org))
    for k in range(1, kval + 1):
        if k <= kval / 2:
            centers[k-1] = 1.0
            dispersions[k-1] = step * k
        else:
            centers[k-1] = 0.0
            dispersions[k-1] = step * (kval - k + 1)

    model = PPanGGOLiNEM(
        n_clusters=kval, beta=beta, family="bernoulli", dispersion="sk_",
        proportion="pk", init="param_file", max_iter=100, tol=0.01,
        convergence="classification", random_state=42, verbose=0,
        init_centers=centers, init_dispersions=dispersions,
        init_proportions=np.full(kval, 1.0 / kval),
    )
    model.fit(G)

    index_fam = []
    with open(nem_dir / "nem_file.index") as f:
        for line in f:
            index_fam.append(line.split("\t")[1].strip())

    parti = {0: "P", kval - 1: "C"}
    for i in range(1, kval - 1):
        parti[i] = "S" + str(i)

    C_membership = model.membership_
    parts_py = {}
    max_probs = {}
    for i, name in enumerate(index_fam):
        row = C_membership[i]
        mp = float(max(row))
        pos = [j for j, p in enumerate(row) if p == mp]
        parts_py[name] = "S_" if (len(pos) > 1 or mp < 0.5) else parti[pos[0]]
        max_probs[name] = mp

    # --- Read C NEM .mf header for iteration info ---
    mf_lines = (nem_dir / f"nem_file_{kval}.mf").read_text().splitlines()
    print("=== C NEM .mf header ===")
    for line in mf_lines[:6]:
        print(" ", line)

    # --- Compare final cluster parameters ---
    cnem_params = {}
    for k, line in enumerate(mf_lines[-kval:]):
        vector = line.split()
        proportions_c = float(vector[nb_org])
        epsilon_k = [float(v) for v in vector[nb_org + 1:]]
        label = ["P", "S1", "C"][k]
        cnem_params[label] = {"proportion": proportions_c,
                               "mean_epsilon": np.mean(epsilon_k)}
    pynem_params = {}
    for k_idx, label in [(0, "P"), (1, "S1"), (2, "C")]:
        pynem_params[label] = {"proportion": model.proportions_[k_idx],
                                "mean_epsilon": np.mean(model.dispersions_[k_idx])}

    print("\n=== Final cluster parameters ===")
    print(f"  {'Label':<5} {'C NEM prop':>12} {'pynem prop':>12} {'C NEM eps':>12} {'pynem eps':>12}")
    for label in ["P", "S1", "C"]:
        cp = cnem_params.get(label, {})
        pp = pynem_params.get(label, {})
        print(f"  {label:<5} {cp.get('proportion', 0):>12.4f} {pp.get('proportion', 0):>12.4f}"
              f" {cp.get('mean_epsilon', 0):>12.4f} {pp.get('mean_epsilon', 0):>12.4f}")

finally:
    td_obj.cleanup()

print(f"\npynem ran {model.n_iter_} EM iterations")

agree = [n for n in parts_c if parts_c[n] == parts_py.get(n)]
disagree = [n for n in parts_c if parts_c[n] != parts_py.get(n)]
print(f"Agreeing: {len(agree)}, Disagreeing: {len(disagree)}")

conf_agree = [max_probs[n] for n in agree if n in max_probs]
conf_disagree = [max_probs[n] for n in disagree if n in max_probs]
print(f"\nMax membership confidence (pynem side):")
print(f"  Agreeing families    mean={np.mean(conf_agree):.4f}  min={np.min(conf_agree):.4f}")
print(f"  Disagreeing families mean={np.mean(conf_disagree):.4f}  min={np.min(conf_disagree):.4f}")

transitions = Counter((parts_c[n], parts_py[n]) for n in disagree)
print(f"\nTransitions (C NEM -> pynem):")
for (a, b), cnt in sorted(transitions.items(), key=lambda x: -x[1]):
    print(f"  {a} -> {b}: {cnt}")

buckets = Counter()
for n in disagree:
    if n in max_probs:
        mp = max_probs[n]
        if mp < 0.55:
            buckets["[0.50,0.55)"] += 1
        elif mp < 0.60:
            buckets["[0.55,0.60)"] += 1
        elif mp < 0.70:
            buckets["[0.60,0.70)"] += 1
        elif mp < 0.80:
            buckets["[0.70,0.80)"] += 1
        else:
            buckets["[0.80,1.00]"] += 1

print(f"\nMax-membership bins for disagreeing families (pynem):")
for label in ["[0.50,0.55)", "[0.55,0.60)", "[0.60,0.70)", "[0.70,0.80)", "[0.80,1.00]"]:
    print(f"  {label}: {buckets.get(label, 0)}")
