"""Layer-specific expression features in the pooled full tree (Track E), step E0.

E0 fixes the labels and the per-round criterion before any feature is ranked
(see ``plans/track_e_layer_features.md``). Decisions agreed with the author on
2026-10-15:

* **Cohort and rounds.** The pooled cosine context
  (``pooled_analysis.build_context``) and its contraction rounds, counted
  bottom-up (round 1 holds the 491 measured leaves). Rounds 5-8 hold only
  71 cells and are pooled into one stratum, so the strata are 1, 2, 3, 4 and
  5-8.
* **Label.** A cell's Schema D mixture: the mean one-hot Schema D class over
  its terminal descendants in the complete lineage (``data/cell_lineage.json``).
  Cohort leaves cut off by the tracking time limit therefore get the mixture of
  the cells they will produce; biological terminals are one-hot.
  ``programmed_death`` is a class by default (``--exclude-death`` drops it).
* **Criterion: equal-relatedness separation.** Over pairs of cells in one
  stratum, grouped by lineage path length (sisters 2, first cousins 4, ...),
  the Spearman correlation between cosine expression distance and Bray-Curtis
  distance between mixtures, averaged over path lengths with equal weight, so
  close relatives count as much as the far more numerous distant pairs (the
  pair-weighted mean is reported beside it). Comparing only equally related
  pairs keeps features that merely mark ancestry from scoring well. Ancestor-descendant pairs (possible only in the
  pooled stratum) are excluded. A path length enters only if it has at least
  ``MIN_PATH_PAIRS`` pairs and its mixture distances vary; which path lengths
  enter depends on the labels alone, so every feature set is scored on the same
  pairs. A feature set whose distances are constant within a path length
  contributes zero there. The sisters-only correlation is reported beside it
  (the Track A analogue), with a bootstrap over sister pairs.
* **Lineage memory.** Among non-ancestral pairs of a stratum that share their
  majority Schema D class and have near-identical mixtures (Bray-Curtis at most
  ``MEMORY_BC_MAX``), the Spearman correlation between expression distance and
  path length, against profiles shuffled among the stratum's cells of the same
  majority class. All same-majority pairs are reported as a sensitivity check.
* **Coverage.** Per round and per stratum, how much of each feature set carries
  reporter signal in raw intensities (``s3.csv``).

Feature sets in E0 are the global ones from Track A: the published top-20, the
HV pool (105 proteins) and the whole panel (210), with random 20-sets from the
HV pool as reference. Nothing here uses travel or a front.

Run from the repository root in the dev environment::

    python -m full_tree_pareto.layer_features --run-id e0_20261015
    python -m full_tree_pareto.layer_features --run-id e0_20261015 --verify-only

An existing completed run is never overwritten. Hash-pinned modules
(``pooled_analysis``, ``representation_criteria``) are imported, never edited.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata

from full_tree_pareto import pooled_analysis as pa
from terminal_pareto import representation_criteria as rc

ROOT = pa.ROOT
VERSION = "layer-features-e0-2"
DEFAULT_RUN = "e0_groupweight_20261015"
RUN_ROOT = ROOT / "full_tree_pareto/output/runs/layer_features_v1"
POOLED_RUN = pa.DEFAULT_RUN
PROTEIN_ZSCORE = "data/protein/aggregated_all/s3_zscore.csv"
PROTEIN_RAW = "data/protein/aggregated_all/s3.csv"
TRACK_A = "terminal_pareto/output/runs/representation_criteria_v1/criteria_metric_20261014/analysis"
CODE_SOURCES = ("full_tree_pareto/layer_features.py", "full_tree_pareto/pooled_analysis.py",
                "terminal_pareto/representation_criteria.py", rc.RANKING_SCRIPT)
DATA_SOURCES = ("data/cell_lineage.json", rc.CELL_KEY, PROTEIN_ZSCORE, PROTEIN_RAW,
                f"{TRACK_A}/candidates.csv", f"{TRACK_A}/cell_types.csv")

STRATA = {"1": (1,), "2": (2,), "3": (3,), "4": (4,), "5-8": (5, 6, 7, 8)}
CANDIDATES = ("top-20", "HV pool", "whole panel")
DEATH = "programmed_death"
SEED = 20261015
N_BOOT = 2000          # bootstrap resamples of sister pairs
N_PERM = 1000          # majority-class shuffles for lineage memory
N_RANDOM = 200         # random reference sets from the HV pool
RANDOM_K = 20
MIN_PATH_PAIRS = 10    # a path length enters the criterion with at least this many pairs
MEMORY_BC_MAX = 0.2    # memory pairs: mixtures at most this Bray-Curtis apart
COVERAGE_MIN = 0.10    # a protein is expressed in a set of cells if raw > 0 in at least this share
COVERAGE_MAX = 0.90    # ... and variable if also at most this share


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_path(run_id):
    if not run_id or Path(run_id).name != run_id or run_id in (".", ".."):
        raise ValueError("Run identity must be a single directory name")
    return RUN_ROOT / run_id


# ---- labels -------------------------------------------------------------------------

def schema_d_mixtures(cells, *, include_death=True):
    """Mean one-hot Schema D class over each cell's terminal descendants in the complete lineage.

    Returns the mixture table (cells without any counted descendant are dropped)
    and the number of terminal descendants per cell.
    """
    ranking = rc.load_ranking_module()
    terminals, descendants = ranking.lineage_terminals()
    classes = ranking.terminal_cell_types(terminals)
    known = set(terminals) | set(descendants)
    missing = [c for c in cells if c not in known]
    if missing:
        raise ValueError(f"Cells absent from the complete lineage: {missing[:5]}")
    keys = sorted(set(classes.values()) - (set() if include_death else {DEATH}))
    rows, counts = {}, {}
    for cell in cells:
        tally = Counter(classes[t] for t in (descendants.get(cell) or [cell]))
        counts[cell] = sum(tally.values())
        total = sum(tally[k] for k in keys)
        if total:
            rows[cell] = [tally[k] / total for k in keys]
    return pd.DataFrame.from_dict(rows, orient="index", columns=keys), counts


def majority_class(mixtures):
    """Largest share; ties go to the first class in alphabetical order."""
    return mixtures.columns[np.argmax(mixtures.to_numpy(), axis=1)].to_numpy()


# ---- design ---------------------------------------------------------------------------

@dataclass
class Stratum:
    name: str
    rounds: tuple
    cells: list            # cell names, in pooled canonical order
    a: np.ndarray          # pair index arrays into ``cells`` (non-ancestral pairs only)
    b: np.ndarray
    path: np.ndarray       # lineage path length per pair
    bc: np.ndarray         # Bray-Curtis distance between mixtures per pair
    paths: list            # (path length, pair indices, ranks of bc) entering the criterion
    sisters: np.ndarray    # pair indices of sisters
    majority: np.ndarray   # majority class per cell
    groups: list           # index arrays of cells sharing a majority class (two or more)
    memory: dict           # "near-identical" / "same majority" -> pair indices


def build_stratum(name, rounds, cells, mixtures, lineage, *, min_path_pairs=MIN_PATH_PAIRS,
                  memory_bc_max=MEMORY_BC_MAX):
    cells = [c for c in cells if c in mixtures.index]
    a, b = np.triu_indices(len(cells), 1)
    ancestral = np.array([cells[i] in lineage.ancestors[cells[j]] or cells[j] in lineage.ancestors[cells[i]]
                          for i, j in zip(a, b)], dtype=bool)
    a, b = a[~ancestral], b[~ancestral]
    path = np.array([lineage.path_length(cells[i], cells[j]) for i, j in zip(a, b)], dtype=int)
    bc = squareform(pdist(mixtures.loc[cells].to_numpy(float), "braycurtis"))[a, b]
    paths = []
    for value in np.unique(path):
        index = np.flatnonzero(path == value)
        if len(index) >= min_path_pairs and np.ptp(bc[index]) > 0:
            paths.append((int(value), index, rankdata(bc[index])))
    majority = majority_class(mixtures.loc[cells])
    by_class = {}
    for i, label in enumerate(majority):
        by_class.setdefault(label, []).append(i)
    groups = [np.array(g) for g in by_class.values() if len(g) > 1]
    same = majority[a] == majority[b]
    memory = {"near-identical": np.flatnonzero(same & (bc <= memory_bc_max)), "same majority": np.flatnonzero(same)}
    return Stratum(name, tuple(rounds), cells, a, b, path, bc, paths, np.flatnonzero(path == 2),
                   majority, groups, memory)


def spearman(x, y_ranks):
    """Spearman correlation against precomputed ranks; zero when ``x`` is constant."""
    if np.ptp(x) == 0:
        return 0.0
    return float(np.corrcoef(rankdata(x), y_ranks)[0, 1])


def pair_distances(D, stratum, perm=None):
    """Pair values of a square distance matrix; with ``perm`` cell i carries the profile of cell perm[i]."""
    if perm is not None:
        return D[perm[stratum.a], perm[stratum.b]]
    return D[stratum.a, stratum.b]


def cosine_matrix(X):
    return squareform(pdist(X, "cosine"))


def separation(d, stratum):
    """Equal-relatedness separation (path lengths weighted equally, and by pair count), per-path,
    sisters-only and all-pairs correlations."""
    per_path = [(value, len(index), spearman(d[index], ranks)) for value, index, ranks in stratum.paths]
    weights = np.array([n for _, n, _ in per_path], float)
    rhos = np.array([r for *_, r in per_path], float)
    equal = float(rhos.mean()) if len(per_path) else float("nan")
    weighted = float(np.dot(weights, rhos) / weights.sum()) if len(per_path) else float("nan")
    sis = stratum.sisters
    sisters = spearman(d[sis], rankdata(stratum.bc[sis])) if len(sis) > 2 and np.ptp(stratum.bc[sis]) > 0 else float("nan")
    return dict(separation=equal, separation_pair_weighted=weighted, sisters=sisters,
                all_pairs=spearman(d, rankdata(stratum.bc))), per_path


def sister_bootstrap(d, stratum, rng, n_boot):
    sis = stratum.sisters
    values = np.empty(n_boot)
    for k in range(n_boot):
        pick = sis[rng.integers(len(sis), size=len(sis))]
        bc = stratum.bc[pick]
        values[k] = spearman(d[pick], rankdata(bc)) if np.ptp(bc) > 0 else np.nan
    return [float(q) for q in np.nanquantile(values, [.025, .975])]


def majority_permutation(stratum, rng):
    perm = np.arange(len(stratum.cells))
    for group in stratum.groups:
        perm[group] = rng.permutation(group)
    return perm


def lineage_memory(X, stratum, rng, n_perm):
    """Spearman of expression distance with path length among same-class pairs, against majority-class shuffles."""
    rows = []
    D = cosine_matrix(X)
    perms = [majority_permutation(stratum, rng) for _ in range(n_perm)]
    for version, index in stratum.memory.items():
        row = dict(memory_pairs=version, n_pairs=len(index), n_sister_pairs=int(np.sum(stratum.path[index] == 2)))
        if len(index) < 3 or np.ptp(stratum.path[index]) == 0:
            rows.append(row)
            continue
        ranks = rankdata(stratum.path[index])
        observed = spearman(pair_distances(D, stratum)[index], ranks)
        null = np.array([spearman(pair_distances(D, stratum, p)[index], ranks) for p in perms])
        row.update(rho=observed, null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                   null_lo=float(np.quantile(null, .025)), null_hi=float(np.quantile(null, .975)),
                   excess=observed - float(null.mean()),
                   z=float((observed - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) > 0 else float("nan"),
                   p=float((1 + np.sum(null >= observed)) / (1 + n_perm)))
        rows.append(row)
    return rows


# ---- context ------------------------------------------------------------------------------

@dataclass
class Context:
    pooled: object        # pooled_analysis.Context
    zscore: pd.DataFrame  # cells x 210 z-scored proteins (NaN -> 0, as in the pooled run)
    raw: pd.DataFrame     # cells x 210 raw intensities (NaN -> 0)
    mixtures: pd.DataFrame
    descendants: dict     # cell -> number of terminal descendants
    rounds: dict          # cell -> contraction round
    strata: dict          # name -> Stratum
    candidates: dict      # name -> features
    lineage: object
    include_death: bool
    sources: dict


def load_context(*, include_death=True):
    pooled = pa.build_context()
    zscore = pd.read_csv(ROOT / PROTEIN_ZSCORE, index_col=0).T.fillna(0)
    raw = pd.read_csv(ROOT / PROTEIN_RAW, index_col=0).fillna(0)[zscore.columns]
    track_a = pd.read_csv(ROOT / TRACK_A / "candidates.csv")
    track_a = track_a[track_a.dataset == "CE protein"].set_index("candidate").features.str.split(";")
    candidates = {name: list(track_a[name]) for name in CANDIDATES}
    ranking = rc.load_ranking_module()
    lineage = rc.build_lineage(rc.dl.load_json(str(ROOT / "data/cell_lineage.json")), ranking.map_names)
    rounds = {pooled.names[c]: r + 1 for r, layer in enumerate(pooled.layers) for c in layer}
    cells = [n for n in pooled.names if n in rounds]
    mixtures, descendants = schema_d_mixtures(cells, include_death=include_death)
    strata = {name: build_stratum(name, members, [c for c in cells if rounds[c] in members], mixtures, lineage)
              for name, members in STRATA.items()}
    sources = {name: file_hash(ROOT / name) for name in (*CODE_SOURCES, *DATA_SOURCES)}
    sources[str((POOLED_RUN / "analysis/context.json").relative_to(ROOT))] = file_hash(POOLED_RUN / "analysis/context.json")
    return Context(pooled, zscore, raw, mixtures, descendants, rounds, strata, candidates, lineage, include_death, sources)


def matrix(ctx, stratum, features):
    return ctx.zscore.loc[stratum.cells, list(features)].to_numpy(float)


def check_reproduction(ctx):
    """The context is the pooled cosine run's; labels and feature sets agree with Track A."""
    failures = []
    saved = json.loads((POOLED_RUN / "analysis/context.json").read_text())
    if saved["context_key"] != ctx.pooled.identity["context_key"]:
        failures.append("pooled context differs from the cosine run")
    if ctx.candidates["top-20"] != ctx.pooled.identity["proteins"]:
        failures.append("Track A top-20 differs from the pooled run's proteins")
    delta = float(np.max(np.abs(ctx.zscore.loc[ctx.pooled.names, ctx.candidates["top-20"]].to_numpy(float)
                                - ctx.pooled.expression)))
    if delta > 0:
        failures.append(f"z-scored top-20 differs from the pooled expression by {delta}")
    if sorted(ctx.candidates["whole panel"]) != sorted(ctx.zscore.columns):
        failures.append("whole panel differs from the z-scored table")
    types = pd.read_csv(ROOT / TRACK_A / "cell_types.csv").set_index("cell").schema_d
    shared = [c for c in types.index if c in ctx.mixtures.index]
    wrong = [c for c in shared if ctx.mixtures.loc[c].idxmax() != types[c] or ctx.mixtures.loc[c].max() != 1]
    if len(shared) < 200 or wrong:
        failures.append(f"Schema D labels disagree with Track A on {len(wrong)} of {len(shared)} terminal cells")
    unlabeled = [c for c in ctx.rounds if c not in ctx.mixtures.index]
    if ctx.include_death and unlabeled:
        failures.append(f"{len(unlabeled)} cohort cells have no Schema D mixture")
    return dict(context_key=ctx.pooled.identity["context_key"], top20_max_abs_delta=delta,
                track_a_terminals_checked=len(shared), cells_without_mixture=unlabeled, failures=failures)


# ---- computation ----------------------------------------------------------------------------

def label_table(ctx):
    rows = []
    for name, stratum in ctx.strata.items():
        for cell, majority in zip(stratum.cells, stratum.majority):
            rows.append(dict(cell=cell, round=ctx.rounds[cell], stratum=name, terminal_descendants=ctx.descendants[cell],
                             majority=majority, **ctx.mixtures.loc[cell].to_dict()))
    return pd.DataFrame(rows)


def design_table(ctx):
    rows = []
    for name, s in ctx.strata.items():
        entering = {value for value, _, _ in s.paths}
        for value in np.unique(s.path):
            index = s.path == value
            rows.append(dict(stratum=name, path_length=int(value), n_pairs=int(index.sum()),
                             mean_bc=float(s.bc[index].mean()), in_criterion=int(value) in entering,
                             memory_near_identical=int(np.isin(np.flatnonzero(index), s.memory["near-identical"]).sum()),
                             memory_same_majority=int(np.isin(np.flatnonzero(index), s.memory["same majority"]).sum())))
    return pd.DataFrame(rows)


def coverage_tables(ctx):
    """Reporter signal per round and stratum: summary per feature set, and per protein."""
    scopes = [(f"round {r}", [c for c, k in ctx.rounds.items() if k == r]) for r in range(1, len(ctx.pooled.layers) + 1)]
    scopes += [(f"stratum {name}", s.cells) for name, s in ctx.strata.items()]
    summary, proteins = [], []
    for scope, cells in scopes:
        R = ctx.raw.loc[cells]
        expressed = (R > 0).mean()
        log_var = np.log1p(R).var(ddof=1)
        for protein in R.columns:
            proteins.append(dict(scope=scope, protein=protein, expressed_fraction=float(expressed[protein]),
                                 log1p_variance=float(log_var[protein])))
        for name, features in ctx.candidates.items():
            e = expressed[features]
            summary.append(dict(scope=scope, n_cells=len(cells), candidate=name, n_features=len(features),
                                expressed=int((e >= COVERAGE_MIN).sum()),
                                variable=int(((e >= COVERAGE_MIN) & (e <= COVERAGE_MAX)).sum()),
                                nonzero_fraction=float((R[features] > 0).to_numpy().mean()),
                                silent_cells=int((R[features] == 0).all(axis=1).sum())))
    return pd.DataFrame(summary), pd.DataFrame(proteins)


def draw_random_sets(ctx, *, n=N_RANDOM, k=RANDOM_K, seed=SEED, max_rejections=10_000):
    """Random k-sets from the Track A HV pool; sets leaving a zero profile in any stratum are redrawn."""
    pool = ctx.candidates["HV pool"]
    rng = np.random.default_rng([seed, 99, k])
    rows = []
    for draw in range(n):
        for rejected in range(max_rejections + 1):
            chosen = [pool[i] for i in rng.choice(len(pool), size=k, replace=False)]
            if all((np.linalg.norm(matrix(ctx, s, chosen), axis=1) > 0).all() for s in ctx.strata.values()):
                break
        else:
            raise ValueError("Could not draw a valid random set")
        rows.append(dict(draw=draw, rejected_draws=rejected, features=";".join(chosen)))
    return pd.DataFrame(rows)


def _candidate_job(ctx, name, stratum_name, seed, n_boot, n_perm):
    stratum = ctx.strata[stratum_name]
    X = matrix(ctx, stratum, ctx.candidates[name])
    if not (np.linalg.norm(X, axis=1) > 0).all():
        raise ValueError(f"{name} leaves a zero profile in stratum {stratum_name}")
    d = pair_distances(cosine_matrix(X), stratum)
    row, per_path = separation(d, stratum)
    rng = np.random.default_rng(seed)
    lo, hi = sister_bootstrap(d, stratum, rng, n_boot) if n_boot else (float("nan"), float("nan"))
    ident = dict(stratum=stratum_name, candidate=name)
    criteria = dict(**ident, n_cells=len(stratum.cells), n_pairs=len(d), n_criterion_pairs=sum(n for _, n, _ in per_path),
                    n_sister_pairs=len(stratum.sisters), **row, sisters_lo=lo, sisters_hi=hi)
    paths = [dict(**ident, path_length=v, n_pairs=n, rho=r) for v, n, r in per_path]
    memory = [dict(**ident, **m) for m in lineage_memory(X, stratum, rng, n_perm)]
    return criteria, paths, memory


def _random_job(ctx, chunk):
    rows = []
    for draw, features in chunk:
        for name, stratum in ctx.strata.items():
            row, _ = separation(pair_distances(cosine_matrix(matrix(ctx, stratum, features.split(";"))), stratum), stratum)
            rows.append(dict(draw=draw, stratum=name, **row))
    return rows


def compute(ctx, *, n_jobs, n_boot=N_BOOT, n_perm=N_PERM, n_random=N_RANDOM, seed=SEED):
    from joblib import Parallel, delayed

    jobs = [(c, s, [seed, ci, si]) for ci, c in enumerate(CANDIDATES) for si, s in enumerate(STRATA)]
    results = Parallel(n_jobs=min(n_jobs, len(jobs)))(delayed(_candidate_job)(ctx, c, s, sd, n_boot, n_perm)
                                                      for c, s, sd in jobs)
    random_sets = draw_random_sets(ctx, n=n_random, seed=seed) if n_random else pd.DataFrame()
    work = [] if random_sets.empty else list(random_sets[["draw", "features"]].itertuples(index=False, name=None))
    chunks = [work[i:i + 10] for i in range(0, len(work), 10)]
    random_metrics = pd.DataFrame([r for rows in Parallel(n_jobs=n_jobs)(delayed(_random_job)(ctx, ch) for ch in chunks)
                                   for r in rows])
    criteria = pd.DataFrame([r for r, _, _ in results])
    if not random_metrics.empty:
        criteria = criteria.merge(placement(criteria, random_metrics), on=["stratum", "candidate"])
    return dict(criteria=criteria, path_criteria=pd.DataFrame([p for _, ps, _ in results for p in ps]),
                memory=pd.DataFrame([m for _, _, ms in results for m in ms]),
                random_sets=random_sets, random_metrics=random_metrics)


def placement(criteria, random_metrics):
    """Share of random sets below each candidate (ties count half), and the random median and 5-95% range."""
    rows = []
    for _, row in criteria.iterrows():
        ref = random_metrics[random_metrics.stratum == row.stratum]
        out = dict(stratum=row.stratum, candidate=row.candidate)
        for column in ("separation", "separation_pair_weighted", "sisters"):
            values = ref[column].to_numpy(float)
            out[f"{column}_random_median"] = float(np.median(values))
            out[f"{column}_random_q05"], out[f"{column}_random_q95"] = (float(q) for q in np.quantile(values, [.05, .95]))
            out[f"{column}_percentile"] = float(np.mean(values < row[column]) + .5 * np.mean(values == row[column]))
        rows.append(out)
    return pd.DataFrame(rows)


def settings(*, include_death=True, n_boot=N_BOOT, n_perm=N_PERM, n_random=N_RANDOM, seed=SEED):
    return dict(seed=seed, n_boot=n_boot, n_perm=n_perm, n_random=n_random, random_k=RANDOM_K,
                strata={k: list(v) for k, v in STRATA.items()}, include_death=include_death,
                label="schema_d_mixture_of_terminal_descendants_complete_lineage",
                criterion="spearman_cosine_vs_braycurtis_within_path_length_equal_weight_per_path_length",
                min_path_pairs=MIN_PATH_PAIRS, memory_bc_max=MEMORY_BC_MAX,
                memory_null="profiles_shuffled_within_majority_schema_d_class",
                coverage_min=COVERAGE_MIN, coverage_max=COVERAGE_MAX, candidates=list(CANDIDATES),
                random_pool="track_a_protein_hv_pool", numpy_version=np.__version__)


TABLES = ("labels", "design", "coverage", "protein_coverage", "criteria", "path_criteria", "memory",
          "random_sets", "random_metrics")


def tables(ctx, results):
    coverage, protein_coverage = coverage_tables(ctx)
    return dict(labels=label_table(ctx), design=design_table(ctx), coverage=coverage,
                protein_coverage=protein_coverage, **results)


def write_run(out, *, n_jobs, include_death=True, **overrides):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    ctx = load_context(include_death=include_death)
    report = check_reproduction(ctx)
    if report["failures"]:
        raise ValueError("Reproduction gates failed:\n  " + "\n  ".join(report["failures"]))
    print("Pooled context, top-20 and Schema D labels reproduced.", flush=True)
    output = tables(ctx, compute(ctx, n_jobs=n_jobs, **overrides))
    analysis.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        output[name].to_csv(analysis / f"{name}.csv", index=False)
    (analysis / "reproduction.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    run_settings = settings(include_death=include_death, **overrides)
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=run_settings, sources=ctx.sources),
                                         sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=run_settings, sources=ctx.sources,
                  elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def _compare(fresh, saved, label, tol=1e-12):
    fresh, saved = fresh.reset_index(drop=True), saved.reset_index(drop=True)
    if len(fresh) != len(saved) or list(fresh.columns) != list(saved.columns):
        raise ValueError(f"Replay mismatch {label}: shape or columns")
    for column in saved.columns:
        if pd.api.types.is_numeric_dtype(saved[column]) and not pd.api.types.is_bool_dtype(saved[column]):
            if not np.allclose(fresh[column].astype(float), saved[column].astype(float), rtol=tol, atol=tol, equal_nan=True):
                raise ValueError(f"Replay mismatch {label} column {column}")
        elif not (fresh[column].astype(str) == saved[column].astype(str)).all():
            raise ValueError(f"Replay mismatch {label} column {column}")


def validate_run(out, *, n_jobs=1, check_sources=True):
    """Check file and source hashes, then recompute every table and compare at 1e-12."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported layer-features cache")
    for name, expected in record["files"].items():
        if file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    s = record["settings"]
    ctx = load_context(include_death=s["include_death"])
    if check_reproduction(ctx)["failures"]:
        raise ValueError("Reproduction gates fail on current inputs")
    fresh = tables(ctx, compute(ctx, n_jobs=n_jobs, n_boot=s["n_boot"], n_perm=s["n_perm"],
                                n_random=s["n_random"], seed=s["seed"]))
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        _compare(fresh[name], saved, name)
    return dict(analysis_id=record["analysis_id"], tables=len(TABLES),
                verified_source_hashes=len(record["sources"]) if check_sources else 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--jobs", type=int, default=min(24, os.cpu_count() or 1))
    parser.add_argument("--exclude-death", action="store_true", help="mixtures over living Schema D classes only")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = run_path(args.run_id)
    if not args.verify_only:
        record = write_run(out, n_jobs=args.jobs, include_death=not args.exclude_death)
        print(json.dumps({k: record[k] for k in ("analysis_id", "elapsed_seconds")}, indent=2))
    print(json.dumps(validate_run(out, n_jobs=args.jobs), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
