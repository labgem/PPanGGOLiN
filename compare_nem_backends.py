#!/usr/bin/env python3
"""compare_nem_backends.py
~~~~~~~~~~~~~~~~~~~~~~~~~
Run ``ppanggolin partition`` with both the C NEM backend (default) and the
pynem backend (``--use_pynem``), then compare wall time, peak memory, and
partition assignments on the same pangenome.

Usage
-----
    conda activate ppanggo_dev
    python compare_nem_backends.py path/to/pangenome.h5 [-K 3] [-b 2.5] [--cpu 2]
"""

import argparse
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np


# ── partition-running ─────────────────────────────────────────────────────────

def run_partition(h5_copy: Path, extra_args: list, label: str):
    """Run ``ppanggolin partition`` and return ``(wall_s, peak_rss_mb, timing_breakdown)``.

    Uses ``/usr/bin/time -v`` when available to capture peak RSS.
    Prints ppanggolin's own log output to stdout in real time via Popen.
    Exits the process on non-zero return code.
    ``timing_breakdown`` is ``{step_name: total_elapsed_s}`` parsed from
    ``[TIMING backend] step=<name> elapsed_s=<value>`` lines in stderr.
    """
    has_time_cmd = Path("/usr/bin/time").exists()

    cmd = (
        ["/usr/bin/time", "-v"] if has_time_cmd else []
    ) + [
        "ppanggolin", "partition",
        "-p", str(h5_copy),
        "--force",
        *extra_args,
    ]

    t0 = time.perf_counter()
    # ppanggolin writes tqdm progress bars to stderr; INFO-level log lines go to stdout.
    # /usr/bin/time -v appends its report to stderr after the process ends.
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout_lines, stderr_lines = [], []
    # Drain stderr (tqdm + /usr/bin/time output) in real time.
    for line in proc.stderr:
        stderr_lines.append(line)
        sys.stdout.write("  | " + line)
        sys.stdout.flush()
    # Drain stdout after stderr is exhausted (INFO logs + [TIMING] lines land here).
    for line in proc.stdout:
        if "TIMING" not in line:
            print(line.rstrip())
        stdout_lines.append(line)
    proc.wait()
    wall_s = time.perf_counter() - t0

    # Parse peak RSS from /usr/bin/time -v output
    peak_rss_mb = None
    if has_time_cmd:
        for line in stderr_lines:
            m = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", line)
            if m:
                peak_rss_mb = int(m.group(1)) / 1024  # KB → MB

    # Parse per-step timing: "[TIMING backend] step=<name> elapsed_s=<value>"
    # These INFO-level lines go to stdout (ppanggolin default log target).
    timing_breakdown: dict = {}
    timing_re = re.compile(r"\[TIMING \w+\] step=(\S+) elapsed_s=([\d.]+)")
    for line in stdout_lines + stderr_lines:
        m = timing_re.search(line)
        if m:
            step, elapsed = m.group(1), float(m.group(2))
            timing_breakdown[step] = timing_breakdown.get(step, 0.0) + elapsed

    if proc.returncode != 0:
        print(f"\n[{label}] FAILED (rc={proc.returncode})")
        sys.exit(1)

    return wall_s, peak_rss_mb, timing_breakdown


# ── partition-reading ─────────────────────────────────────────────────────────

def read_partitions(h5_path: Path) -> dict:
    """Return ``{family_name: partition_label}`` by loading the h5 via ppanggolin."""
    from ppanggolin.pangenome import Pangenome
    from ppanggolin.formats.readBinaries import get_status
    from ppanggolin.formats import check_pangenome_info

    pan = Pangenome()
    pan.file = str(h5_path)
    get_status(pan, h5_path)
    check_pangenome_info(pan, need_families=True, need_partitions=True, disable_bar=True)
    return {fam.name: fam.partition for fam in pan.gene_families}


# ── statistics ───────────────────────────────────────────────────────────────

def adjusted_rand_index(labels_a, labels_b) -> float:
    """Compute ARI without sklearn, using ``scipy.special.comb``."""
    from scipy.special import comb

    a = np.asarray(labels_a)
    b = np.asarray(labels_b)
    n = len(a)

    classes_a = np.unique(a)
    classes_b = np.unique(b)
    ia = {v: i for i, v in enumerate(classes_a)}
    ib = {v: i for i, v in enumerate(classes_b)}
    cont = np.zeros((len(classes_a), len(classes_b)), dtype=np.int64)
    for x, y in zip(a, b):
        cont[ia[x], ib[y]] += 1

    sum_comb_c = sum(comb(n_ij, 2, exact=True) for n_ij in cont.ravel())
    sum_comb_a = sum(comb(row.sum(), 2, exact=True) for row in cont)
    sum_comb_b = sum(comb(col.sum(), 2, exact=True) for col in cont.T)
    comb_n = comb(n, 2, exact=True)

    if comb_n == 0:
        return 1.0
    expected = sum_comb_a * sum_comb_b / comb_n
    max_val = (sum_comb_a + sum_comb_b) / 2
    if max_val == expected:
        return 1.0 if sum_comb_c == expected else 0.0
    return (sum_comb_c - expected) / (max_val - expected)


def _partition_group(label: str) -> str:
    """Map a label to its top-level group: 'C', 'P', or 'S'."""
    if label == "C":
        return "C"
    if label == "P":
        return "P"
    return "S"  # S1, S2, S_, …


def compare_partitions(p_c: dict, p_py: dict) -> dict:
    """Compare two ``{family_name: label}`` dicts and return agreement stats."""
    families = sorted(p_c.keys())
    labels_c  = [p_c[f]  for f in families]
    labels_py = [p_py[f] for f in families]
    n = len(families)

    direct_agree = sum(a == b for a, b in zip(labels_c, labels_py))

    all_labels = sorted(set(labels_c) | set(labels_py))
    enc = {lbl: i for i, lbl in enumerate(all_labels)}
    int_c  = [enc[l] for l in labels_c]
    int_py = [enc[l] for l in labels_py]
    ari = adjusted_rand_index(int_c, int_py)

    # --- per-label agreement (families assigned to that label by C NEM) ---
    # For each exact label (C, P, S1, S2, …) and for the merged S group,
    # count how many families that C NEM called that label were also called
    # the same label (or the same S group) by pynem.
    all_shell_labels = sorted(
        {lbl for lbl in all_labels if _partition_group(lbl) == "S"}
    )
    per_label: dict = {}   # exact label → (n_in_cnem, n_agreed)
    per_group: dict = {}   # 'C', 'P', 'S' → (n_in_cnem, n_agreed)

    for lbl in all_labels:
        group = _partition_group(lbl)
        pairs = [(lc, lp) for lc, lp in zip(labels_c, labels_py) if lc == lbl]
        n_in = len(pairs)
        n_ok = sum(1 for _, lp in pairs if lp == lbl)
        per_label[lbl] = (n_in, n_ok)
        # accumulate group totals
        prev_in, prev_ok = per_group.get(group, (0, 0))
        per_group[group] = (prev_in + n_in, prev_ok + n_ok)

    return {
        "n_families": n,
        "direct_agree": direct_agree,
        "direct_agree_pct": direct_agree / n * 100,
        "ari": ari,
        "cnt_c":  Counter(labels_c),
        "cnt_py": Counter(labels_py),
        "per_label": per_label,        # {label: (n_in_cnem, n_agreed_exact)}
        "per_group": per_group,        # {'C','P','S': (n_in_cnem, n_agreed_same_group)}
        "all_shell_labels": all_shell_labels,
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Compare C NEM vs pynem backends for ppanggolin partition.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("pangenome", type=Path, help="Input pangenome .h5 file")
    ap.add_argument("-K", "--nb_of_partitions", type=int, default=-1,
                    metavar="K",
                    help="Number of partitions (-1 = auto-detect)")
    ap.add_argument("-b", "--beta", type=float, default=2.5,
                    help="Beta (spatial smoothing strength)")
    ap.add_argument("-ms", "--max_degree_smoothing", type=float, default=10,
                    help="Max. degree for smoothing")
    ap.add_argument("-c", "--cpu", type=int, default=1, help="CPU count")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--free_dispersion", action="store_true",
                    help="Use free dispersion")
    ap.add_argument("--chunk_size", type=int, default=500,
                    help="Chunk size for large pangenomes")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO)

    if not args.pangenome.exists():
        sys.exit(f"Error: {args.pangenome} does not exist.")

    shared_args = [
        "-K", str(args.nb_of_partitions),
        "-b", str(args.beta),
        "-ms", str(args.max_degree_smoothing),
        "-c", str(args.cpu),
        "--seed", str(args.seed),
        "--chunk_size", str(args.chunk_size),
        "--disable_prog_bar"
    ]
    if args.free_dispersion:
        shared_args.append("--free_dispersion")

    SEP = "=" * 65
    print(SEP)
    print(f"Pangenome     : {args.pangenome}")
    print(f"Shared options: {' '.join(shared_args)}")
    print(SEP)

    with tempfile.TemporaryDirectory(prefix="ppang_compare_") as tmpdir:
        h5_c  = Path(tmpdir) / "pangenome_cnem.h5"
        h5_py = Path(tmpdir) / "pangenome_pynem.h5"
        shutil.copy2(args.pangenome, h5_c)
        shutil.copy2(args.pangenome, h5_py)

        # ── C NEM ──────────────────────────────────────────────────────
        print(f"\n{'─'*65}")
        print("[1/2] Running C NEM backend …")
        print(f"{'─'*65}")
        wall_c, rss_c, timing_c = run_partition(h5_c, shared_args, "C NEM")

        # ── pynem ──────────────────────────────────────────────────────
        print(f"\n{'─'*65}")
        print("[2/2] Running pynem backend (--use_pynem) …")
        print(f"{'─'*65}")
        wall_py, rss_py, timing_py = run_partition(h5_py, shared_args + ["--use_pynem"], "pynem")

        # ── read partitions ────────────────────────────────────────────
        print("\nReading partition results from both runs …")
        parts_c  = read_partitions(h5_c)
        parts_py = read_partitions(h5_py)
    # temp dir cleaned up; run comparison on in-memory dicts

    stats = compare_partitions(parts_c, parts_py)
    n         = stats["n_families"]
    cnt_c     = stats["cnt_c"]
    cnt_py    = stats["cnt_py"]
    all_lbl   = sorted(set(cnt_c) | set(cnt_py))
    speedup   = wall_c / wall_py if wall_py > 0 else float("inf")

    # ── print results ─────────────────────────────────────────────────
    print(f"\n{SEP}")
    print("RESULTS")
    print(SEP)

    print(f"\n  {'Metric':<35} {'C NEM':>10} {'pynem':>10}")
    print(f"  {'-'*57}")
    print(f"  {'Wall time (s)':<35} {wall_c:>10.2f} {wall_py:>10.2f}")
    if speedup >= 1:
        print(f"  {'Speedup (C NEM / pynem)':<35} {speedup:>10.2f}x {'':>10}  ← pynem faster")
    else:
        print(f"  {'Speedup (pynem / C NEM)':<35} {1/speedup:>10.2f}x {'':>10}  ← C NEM faster")

    if rss_c is not None or rss_py is not None:
        rss_c_str  = f"{rss_c:.1f}"  if rss_c  is not None else "N/A"
        rss_py_str = f"{rss_py:.1f}" if rss_py is not None else "N/A"
        print(f"  {'Peak RSS (MB)':<35} {rss_c_str:>10} {rss_py_str:>10}")

    print(f"\n  Partition counts ({n} gene families total):")
    print(f"  {'Label':<12} {'C NEM':>10} {'pynem':>10}")
    print(f"  {'-'*34}")
    for lbl in all_lbl:
        print(f"  {lbl:<12} {cnt_c.get(lbl, 0):>10} {cnt_py.get(lbl, 0):>10}")
    print(f"  {'-'*34}")
    print(f"  {'Total':<12} {n:>10} {n:>10}")

    print(f"\n  Partition agreement:")
    print(f"  {'Direct agreement':<35} {stats['direct_agree']:>7} / {n}"
          f"  ({stats['direct_agree_pct']:.2f} %)")
    print(f"  {'Adjusted Rand Index (ARI)':<35} {stats['ari']:.4f}")

    # per-group agreement (C NEM families as reference)
    per_label = stats["per_label"]
    per_group = stats["per_group"]
    all_shell  = stats["all_shell_labels"]

    print(f"\n  Per-group agreement (reference = C NEM assignment):")
    hdr = f"  {'Group/Label':<14} {'in C NEM':>10} {'agreed':>10} {'%':>8}"
    print(hdr)
    print(f"  {'-'*44}")
    for group in ("P", "C", "S"):
        n_in, n_ok = per_group.get(group, (0, 0))
        pct = n_ok / n_in * 100 if n_in else float("nan")
        label_str = f"{group} (all shell)" if group == "S" and len(all_shell) > 1 else group
        print(f"  {label_str:<14} {n_in:>10} {n_ok:>10} {pct:>7.2f}%")
        # per sub-shell breakdown when there is more than one shell label
        if group == "S":
            for slbl in all_shell:
                n_in_s, n_ok_s = per_label.get(slbl, (0, 0))
                pct_s = n_ok_s / n_in_s * 100 if n_in_s else float("nan")
                print(f"  {'  '+slbl:<14} {n_in_s:>10} {n_ok_s:>10} {pct_s:>7.2f}%")
    print(f"  {'-'*44}")

    # ── per-step timing breakdown ──────────────────────────────────────
    if timing_c or timing_py:
        all_steps = sorted(set(timing_c) | set(timing_py))
        print(f"\n  Per-step timing breakdown (summed across all chunks):")
        print(f"  {'Step':<20} {'C NEM':>10} {'pynem':>10} {'ratio':>8}")
        print(f"  {'-'*50}")
        for step in all_steps:
            tc = timing_c.get(step)
            tp = timing_py.get(step)
            tc_str  = f"{tc:.4f}s"  if tc  is not None else "  N/A   "
            tp_str  = f"{tp:.4f}s"  if tp  is not None else "  N/A   "
            if tc and tp:
                ratio = tp / tc
                ratio_str = f"{ratio:.2f}x"
                indicator = " ← slower" if ratio > 1.1 else (" ← faster" if ratio < 0.9 else "")
            else:
                ratio_str = "  N/A"
                indicator = ""
            print(f"  {step:<20} {tc_str:>10} {tp_str:>10} {ratio_str:>8}{indicator}")

    print(f"\n{SEP}\n")


if __name__ == "__main__":
    main()
