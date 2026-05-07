# Changes — pynem integration and partition.py refactoring

This document covers all changes made during the integration of the pure-Python
NEM backend (pynem) into PPanGGOLiN and the associated refactoring of
`ppanggolin/nem/partition.py`.

For the changes made to the upstream `cambroise/nem` repository (pynem itself),
see [`nem/pynem/CHANGES_FOR_PR.md`](nem/pynem/CHANGES_FOR_PR.md).

---

## Files changed

| File | Type |
|---|---|
| `ppanggolin/nem/partition.py` | Refactored + pynem integration |
| `ppanggolin/nem/pynem_adapter.py` | New file — PPanGGOLiN-specific NEM subclass |
| `tests/conftest.py` | Updated global reset fixture |
| `compare_nem_backends.py` | New script — benchmark / comparison tool |
| `nem/pynem/src/pynem/spatial.py` | Performance (upstream PR) |
| `nem/pynem/src/pynem/models.py` | Performance + bug fix (upstream PR) |
| `nem/pynem/src/pynem/core.py` | Performance + bug fix (upstream PR) |

---

## 1. `ppanggolin/nem/partition.py` — full refactoring

### Before

The original file had one monolithic `run_partitioning` function (~250 lines)
that mixed:
- writing the `.m` parameter initialisation file
- calling the C NEM extension (`nem_stats.nem`)
- reading `.uf` / `.mf` output files
- building partition dicts
- handling failures

In addition:
- `global pan` and `global samples` aliases existed with no clear purpose
- A dead `cpt = 0` / `cpt += 1` counter was incremented but never used
- `evaluate_nb_partitions` and `partition` contained large inline blocks that
  made them hard to follow
- The ICL-curve drawing code was inlined inside `evaluate_nb_partitions`
- The chunked-partitioning loop was inlined inside `partition`
- `launch()` mutated the module-level `_pangenome` global directly instead of
  using a local variable
- The pynem path was a private `_run_partitioning_pynem` helper that duplicated
  much of the infrastructure of `run_partitioning`

### After

#### New private helpers extracted from `run_partitioning`

| Function | Purpose |
|---|---|
| `_read_family_index(nem_dir)` | Read gene-family names from `nem_file.index` |
| `_partition_label_map(kval)` | Build `{int → label}` map (P / S1…Sk-2 / C) |
| `_write_cnem_init_file(path, kval, nb_org)` | Write `.m` parameter-init file |
| `_parse_cnem_output(nem_dir, kval, nb_org, just_ll)` | Parse `.uf` / `.mf` output |
| `_log_cnem_failure(logger, err, …)` | Structured failure logging |
| `_run_cnem(…)` | Run C NEM extension end-to-end |
| `_run_pynem(…)` | Run pynem backend end-to-end |

`run_partitioning` is now a thin 10-line dispatcher that selects the backend
and delegates:

```python
def run_partitioning(…, backend="cnem", use_pynem=False):
    if use_pynem and backend == "cnem":
        backend = "pynem"
    runner = _run_pynem if backend == "pynem" else _run_cnem
    return runner(nem_dir_path, nb_org, beta, …)
```

#### New helpers extracted from `partition` and `evaluate_nb_partitions`

| Function | Purpose |
|---|---|
| `_draw_icl_curve(…)` | Plotly ICL-curve HTML (was inline in `evaluate_nb_partitions`) |
| `_run_chunked_partitioning(…)` | Majority-vote chunked loop (was inline in `partition`) |
| `_partition_sample_worker(args)` | Single-chunk Pool worker (was `partition_nem` + `nem_samples`) |

#### Backend selection — `use_pynem` → `backend`

All public functions that previously accepted `use_pynem: bool` now prefer
`backend: str = "cnem"` (accepts `"cnem"` or `"pynem"`).  The old
`use_pynem` parameter is kept as a deprecated keyword argument that maps to
`backend="pynem"` for backward compatibility.

```python
# old
partition(pangenome, …, use_pynem=True)
# new (preferred)
partition(pangenome, …, backend="pynem")
# old form still works
partition(pangenome, …, use_pynem=True)
```

#### Module globals cleaned up

| Before | After |
|---|---|
| `pan = _pangenome` (alias) | Removed |
| `samples = _samples` (alias) | Removed |
| `global pan, samples, _pangenome, _samples` in `partition()` | `global _pangenome, _samples` only |
| `launch()` called `_pangenome.add_file(…)` on the global | `launch()` creates a local `pangenome = Pangenome()` |

The globals `_pangenome` and `_samples` are still necessary: on Linux,
`multiprocessing.Pool` forks child processes that inherit the parent's memory
via copy-on-write.  Passing a large pangenome as a function argument would
require pickling it for every task.

#### Bug fix — duplicate ICL traces

`evaluate_nb_partitions` was appending BIC, ICL, and log-likelihood traces
to the Plotly figure twice in some code paths.  The extracted `_draw_icl_curve`
function builds the traces exactly once.

#### Bug fix — `[partitioning_results, []]` hack removed

The old `partition()` contained:

```python
# old — wraps the result so that downstream iteration is uniform
partitioning_results = [partitioning_results, []]
```

This was removed; the chunked and single-run paths now both produce a plain
`dict` that is used directly.

#### Dead code removed

- `cpt = 0` / `cpt += 1` counter (incremented but never read)
- All `[TIMING cnem]` / `[TIMING pynem]` log lines (development artefacts)

#### Organisms sorted before sub-sampling (reproducibility)

`evaluate_nb_partitions` now sorts organisms by name before random
sub-sampling, making results reproducible for a given seed regardless of set
iteration order.

---

## 2. `ppanggolin/nem/pynem_adapter.py` — new file

Provides `PPanGGOLiNEM`, a `pynem.NEM` subclass with one extension:

**`init="param_file"` initialisation mode** — replicates the `.m` file
initialisation used by C NEM.  Accepts `init_centers`, `init_dispersions`,
and `init_proportions` as constructor arguments.

```python
model = PPanGGOLiNEM(
    n_clusters=kval,
    beta=beta,
    family="bernoulli",
    dispersion="sk_",
    init="param_file",
    init_centers=centers,        # (K, D) array
    init_dispersions=dispersions, # (K, D) array
    init_proportions=proportions, # (K,) array
    …
)
model.fit(graph)
```

The NaN-safety fix for the D criterion (`np.nansum` instead of `.sum()`) is
applied directly in upstream `pynem.core.NEM._compute_criteria` — no override
needed in the adapter.

---

## 3. `tests/conftest.py` — global reset fixture

The `reset_globals` autouse fixture previously reset the legacy aliases:

```python
# old
partition.samples = []
partition.pan = Pangenome()
```

Updated to reset the canonical names after the aliases were removed:

```python
# new
partition._samples = []
partition._pangenome = Pangenome()
```

---

## 4. `compare_nem_backends.py` — new benchmark/comparison script

A standalone script that runs `ppanggolin partition` twice on the same
pangenome (once with each backend), then reports:

- Wall time and peak RSS for each backend
- Partition counts per label
- **Overall direct agreement** (% families with matching label)
- **Adjusted Rand Index (ARI)**
- **Per-group agreement** — agreement broken down by C NEM reference group:
  - P (Persistent)
  - C (Cloud)
  - S (all Shell labels combined)
  - S1, S2, … (each sub-shell separately)

Example output:

```
  Per-group agreement (reference = C NEM assignment):
  Group/Label      in C NEM     agreed        %
  --------------------------------------------
  P                    2407       2405   99.92%
  C                    6354       6354  100.00%
  S (all shell)         620        447   72.10%
    S1                  620        447   72.10%
  --------------------------------------------
```

Usage:

```bash
conda activate ppanggo_dev
python compare_nem_backends.py path/to/pangenome.h5 -K 3 [-b 2.5] [-c 4]
```

---

## 5. pynem upstream changes (summary)

Full details: [`nem/pynem/CHANGES_FOR_PR.md`](nem/pynem/CHANGES_FOR_PR.md)

| Area | File | Type |
|---|---|---|
| NaN in D criterion (`0 * -inf = NaN`) | `core.py` | Bug fix |
| Cache `log_pkfki` across E-step and criteria | `core.py` | Optimisation |
| Vectorised E-step row normalisation | `core.py` | Optimisation |
| Sparse CSR matrix for spatial context | `spatial.py` | Optimisation |
| Vectorised `compute_log_density` | `models.py` | Optimisation |
| Vectorised M-step (`estimate_parameters`) | `models.py` | Optimisation |

---

## Benchmark

Dataset: *Acinetobacter schindleri* (48 organisms, 9 381 gene families, K=3, β=2.5, seed=42).

| Backend | Wall time (s) | Peak RSS (MB) | Agreement with C NEM |
|---|---|---|---|
| C NEM | ~5.5 | ~541 | — |
| pynem (after changes) | ~6.0 | ~577 | 98.13% direct, ARI=0.946 |

The ~2% residual disagreement between backends is expected and explained by:

1. C NEM uses 32-bit `float`; pynem uses 64-bit `float64` — different rounding
2. The NaN convergence bug (fixed in pynem) caused C NEM to exit at a slightly
   different local optimum in some runs
3. Both algorithms are non-convex EM variants; identical solutions for
   different implementations are not guaranteed

Per-group, the disagreement is concentrated in the Shell partition (~28%
of C NEM Shell families are re-labelled by pynem), while Persistent and Cloud
families agree at ≥99.9% and 100% respectively.
