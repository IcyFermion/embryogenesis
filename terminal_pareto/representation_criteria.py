"""Choosing the cell-state representation without using the front (Track A1b).

The Pareto front is the outcome of the terminal analysis, so the expression
representation behind its cell-state cost must be chosen by criteria computed
without travel or the front. This add-on scores candidate representations and
measures how much of the natural lineage's cell-state advantage comes from
lineage memory. Methods are written up in ``METHODS_cell_state_representation.md``.

Datasets (fixed cohorts, cosine distance throughout):

* ``CE protein``: z-scored C. elegans protein, primary 275 pooled edges
  (Figures 2-6). The gold standard for C. elegans.
* ``CE RNA`` and ``CB RNA``: TF expression of each species on the 187
  cross-species edges (Figures 8-9). RNA serves only this terminal
  cross-species comparison, so one TF set is used for both species.

Candidates per modality: the published classifier ``top-20``, the ``HV pool``
(top half by log1p variance over the cohort cells; for RNA ranked jointly over
both species) and the ``whole panel``.

Criteria, none of which uses travel or the front:

* **Sibling separation** (fate discrimination): Cohen's d of cosine distance
  between terminal sisters of different type and sisters of the same type,
  for Schema D types (the classifier labels) and for fine types (WormBase
  cell descriptions), with bootstrap intervals over sister pairs.
* **Lineage memory**: among cohort cells of the same fine type, the Spearman
  correlation between expression distance and lineage path length, against a
  null that shuffles profiles within each fine type (type composition kept,
  lineage alignment broken). Also by path-length bin and by cell type.
* A reference distribution from random 20-feature sets of each HV pool.

Control on the front: **fate-only** cell-state cost, in which every child is
represented by the mean profile of its fine type (Schema D class when it has
no description) while parents keep their measured profiles. The fronts for
measured and fate-only costs are computed with the published sweep
(``cross_species_analysis.analyze`` through ``feature_sensitivity.evaluate``).

Before any computation, the measured top-20 fronts must reproduce the
published caches (primary P0 row, all six 187-edge configurations including
C. briggsae) and the sibling separation of the protein top-20 must reproduce
``elegans_protein_linear_baseline.ipynb`` (cosine 0.684, Euclidean 0.637).

Run from the repository root in the dev environment::

    python terminal_pareto/representation_criteria.py --run-id criteria_20261014
    python terminal_pareto/representation_criteria.py --run-id criteria_20261014 --verify-only

An existing completed run is never overwritten.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import feature_sensitivity as fs

VERSION = "terminal-representation-criteria-v1"
DEFAULT_RUN = "criteria_20261014"
RUN_ROOT = ROOT / "terminal_pareto/output/runs/representation_criteria_v1"
RANKING_SCRIPT = "expression_embedding/elegans_protein_feature_ranking.py"
CELL_KEY = "data/2023-06-29_entropy_cell_key_V2.csv"
CB_RNA = "data/c_briggsae/science.adu8249/c_briggsae_tf.csv"
CODE_SOURCES = ("terminal_pareto/representation_criteria.py", RANKING_SCRIPT)

SEED = 20261014
N_BOOT = 2000           # bootstrap resamples of sister pairs
N_PERM = 1000           # within-type shuffles for lineage memory
N_PERM_TYPE = 500       # shuffles within one type, memory by type
N_RANDOM = 200          # random reference sets per modality
RANDOM_K = 20
MIN_TYPE_CELLS = 6      # memory by type: types with at least this many cohort cells
BIN_MIN_PAIRS = 30      # path-length bins merge until they hold this many pairs
INTERVALS = fs.INTERVALS
CANDIDATES = ("top-20", "HV pool", "whole panel")
# dataset -> (modality, expression frame key, cohort for criteria)
DATASETS = {"CE protein": ("protein", "ce_protein", "primary"),
            "CE RNA": ("rna", "ce_rna", "cross187"),
            "CB RNA": ("rna", "cb_rna", "cross187")}
# Fronts: dataset -> (cohort with that dataset's travel, geometries)
FRONTS = {"CE protein": ("primary", ("raw3d",)),
          "CE RNA": ("cross187", ("raw3d", "xy")),
          "CB RNA": ("cross187_cb", ("raw3d", "xy"))}
# Notebook outputs (elegans_protein_linear_baseline.ipynb, cell 27), printed to 3 decimals.
NOTEBOOK_SIBLING_D = {"cosine": 0.684, "euclidean": 0.637}
NOTEBOOK_TOL = 5e-4
REPRODUCTION_TOL = fs.REPRODUCTION_TOL


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_path(run_id):
    if not run_id or Path(run_id).name != run_id or run_id in (".", ".."):
        raise ValueError("Run identity must be a single directory name")
    return RUN_ROOT / run_id


def load_ranking_module():
    """The protein-ranking script holds the Schema D cell-type assignment used by the classifiers."""
    spec = importlib.util.spec_from_file_location("elegans_protein_feature_ranking", ROOT / RANKING_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---- lineage and cell types ---------------------------------------------------

@dataclass
class Lineage:
    ancestors: dict  # cell -> ancestors from the root, nearest last
    order: dict      # cell -> depth-first position

    def path_length(self, a, b):
        """Divisions between two cells through their last common ancestor (sisters 2, first cousins 4)."""
        A, B = self.ancestors[a] + [a], self.ancestors[b] + [b]
        n = 0
        while n < min(len(A), len(B)) and A[n] == B[n]:
            n += 1
        return len(A) + len(B) - 2 * n


def build_lineage(tree, rename):
    ancestors, order = {}, {}

    def visit(node, path):
        name = rename(node["did"]); ancestors[name] = path; order[name] = len(order)
        for child in node.get("children", []):
            visit(child, path + [name])

    visit(tree, [])
    return Lineage(ancestors, order)


def cell_types(ranking):
    """Schema D class and fine type (WormBase description) for every terminal cell."""
    terminals, _ = ranking.lineage_terminals()
    coarse = ranking.terminal_cell_types(terminals)
    key = pd.read_csv(ROOT / CELL_KEY).drop_duplicates("wormweb.lineage").set_index("wormweb.lineage")
    fine = {}
    for cell in terminals:
        if cell in key.index and isinstance(key.at[cell, "wormweb.desc"], str):
            fine[cell] = key.at[cell, "wormweb.desc"]
    return coarse, fine


# ---- criteria -------------------------------------------------------------------

@dataclass
class Design:
    cells: list
    sisters: dict      # "coarse"/"fine" -> (same-type index pairs, different-type index pairs)
    groups: list       # index arrays of cells sharing a fine type (two or more cells)
    pairs: np.ndarray  # same-fine-type index pairs
    lineage_distance: np.ndarray
    lineage_rank: np.ndarray
    bins: list


def path_bins(distance, min_pairs=BIN_MIN_PAIRS):
    """Consecutive path lengths merged until each bin holds at least ``min_pairs`` pairs."""
    values, bins, low = sorted(set(distance.astype(int))), [], None
    for value in values:
        low = value if low is None else low
        if np.sum((distance >= low) & (distance <= value)) >= min_pairs:
            bins.append((low, value)); low = None
    if low is not None:
        bins[-1] = (bins[-1][0], values[-1]) if bins else (low, values[-1])
    return bins


def build_design(terms, coarse, fine, lineage, *, min_pairs=BIN_MIN_PAIRS):
    cells = [child for child, _ in terms]
    position = {cell: i for i, cell in enumerate(cells)}
    families = {}
    for child, parent in terms:
        families.setdefault(parent, []).append(child)
    sisters_all = [(position[a], position[b]) for kids in families.values()
                   for i, a in enumerate(kids) for b in kids[i + 1:]]
    sisters = {}
    for name, types in (("coarse", coarse), ("fine", fine)):
        typed = [(i, j) for i, j in sisters_all if cells[i] in types and cells[j] in types]
        sisters[name] = tuple(np.array([(i, j) for i, j in typed if (types[cells[i]] == types[cells[j]]) == same],
                                       dtype=int).reshape(-1, 2) for same in (True, False))
    by_type = {}
    for i, cell in enumerate(cells):
        if cell in fine:
            by_type.setdefault(fine[cell], []).append(i)
    groups = [np.array(g) for g in by_type.values() if len(g) > 1]
    pairs = np.array([(a, b) for g in groups for k, a in enumerate(g) for b in g[k + 1:]], dtype=int).reshape(-1, 2)
    distance = np.array([lineage.path_length(cells[a], cells[b]) for a, b in pairs], float)
    return Design(cells, sisters, groups, pairs, distance, rankdata(distance), path_bins(distance, min_pairs))


def cohen_d(D, same, diff):
    if len(same) < 2 or len(diff) < 2:
        return float("nan")
    s, t = D[same[:, 0], same[:, 1]], D[diff[:, 0], diff[:, 1]]
    pooled = np.sqrt((s.var(ddof=1) + t.var(ddof=1)) / 2)
    return float((t.mean() - s.mean()) / pooled) if pooled > 0 else 0.0


def within_type_permutation(design, rng):
    """Shuffle cells within each fine type; cells without a fine type stay in place."""
    perm = np.arange(len(design.cells))
    for group in design.groups:
        perm[group] = rng.permutation(group)
    return perm


def pair_distances(D, design, perm=None):
    if perm is not None:
        D = D[np.ix_(perm, perm)]
    return D[design.pairs[:, 0], design.pairs[:, 1]]


def spearman_with_lineage(distances, design):
    return float(np.corrcoef(rankdata(distances), design.lineage_rank)[0, 1])


def bin_means(distances, design):
    d, x = distances, design.lineage_distance
    return [float(np.mean(d[(x >= lo) & (x <= hi)])) for lo, hi in design.bins]


def criteria(X, design, *, seed, n_boot=N_BOOT, n_perm=N_PERM):
    """Sibling separation with bootstrap intervals; lineage memory with its within-type null."""
    D = squareform(pdist(X, "cosine"))
    rng = np.random.default_rng(seed)
    row = dict(n_cells=len(X), n_distinct_profiles=int(len(np.unique(X.round(10), axis=0))))
    for name in ("coarse", "fine"):
        same, diff = design.sisters[name]
        row[f"sib_{name}"] = cohen_d(D, same, diff)
        row[f"n_sisters_{name}_same"], row[f"n_sisters_{name}_diff"] = len(same), len(diff)
        boot = [cohen_d(D, same[rng.integers(len(same), size=len(same))], diff[rng.integers(len(diff), size=len(diff))])
                for _ in range(n_boot)]
        row[f"sib_{name}_lo"], row[f"sib_{name}_hi"] = (float(q) for q in np.quantile(boot, [.025, .975]))
    observed = pair_distances(D, design)
    row["memory_rho"] = spearman_with_lineage(observed, design)
    null_rho, null_bins = np.empty(n_perm), np.empty((n_perm, len(design.bins)))
    for k in range(n_perm):
        shuffled = pair_distances(D, design, within_type_permutation(design, rng))
        null_rho[k] = spearman_with_lineage(shuffled, design)
        null_bins[k] = bin_means(shuffled, design)
    row.update(n_same_type_pairs=len(design.pairs), memory_null_mean=float(null_rho.mean()),
               memory_null_sd=float(null_rho.std(ddof=1)),
               memory_null_lo=float(np.quantile(null_rho, .025)), memory_null_hi=float(np.quantile(null_rho, .975)),
               memory_excess=row["memory_rho"] - float(null_rho.mean()),
               memory_z=float((row["memory_rho"] - null_rho.mean()) / null_rho.std(ddof=1)),
               memory_p=float((1 + np.sum(null_rho >= row["memory_rho"])) / (1 + n_perm)))
    observed_bins = bin_means(observed, design)
    curve = [dict(bin_lo=int(lo), bin_hi=int(hi),
                  n_pairs=int(np.sum((design.lineage_distance >= lo) & (design.lineage_distance <= hi))),
                  observed_mean=observed_bins[i], null_median=float(np.median(null_bins[:, i])),
                  null_q05=float(np.quantile(null_bins[:, i], .05)), null_q95=float(np.quantile(null_bins[:, i], .95)))
             for i, (lo, hi) in enumerate(design.bins)]
    return row, curve


def point_criteria(X, design):
    """Criteria without intervals or nulls, for the random reference sets."""
    D = squareform(pdist(X, "cosine"))
    return dict(sib_coarse=cohen_d(D, *design.sisters["coarse"]), sib_fine=cohen_d(D, *design.sisters["fine"]),
                memory_rho=spearman_with_lineage(pair_distances(D, design), design),
                n_distinct_profiles=int(len(np.unique(X.round(10), axis=0))))


def memory_by_type(X, design, fine, lineage, *, seed, n_perm=N_PERM_TYPE, min_cells=MIN_TYPE_CELLS):
    """Within one fine type: Spearman of expression distance and path length, shuffled within the type."""
    D = squareform(pdist(X, "cosine"))
    rng = np.random.default_rng(seed)
    rows = []
    for group in design.groups:
        if len(group) < min_cells:
            continue
        cells = [design.cells[i] for i in group]
        a, b = np.triu_indices(len(group), 1)
        path = rankdata([lineage.path_length(cells[i], cells[j]) for i, j in zip(a, b)])
        sub = D[np.ix_(group, group)]
        observed = sub[a, b]
        if np.ptp(observed) == 0:
            rho, p = np.nan, np.nan
        else:
            rho = float(np.corrcoef(rankdata(observed), path)[0, 1])
            null = np.empty(n_perm)
            for k in range(n_perm):
                perm = rng.permutation(len(group))
                null[k] = np.corrcoef(rankdata(sub[np.ix_(perm, perm)][a, b]), path)[0, 1]
            p = float((1 + np.sum(null >= rho)) / (1 + n_perm))
        founders = sorted({founder(cell) for cell in cells})
        rows.append(dict(fine_type=fine[cells[0]], n_cells=len(group), founders=";".join(founders), rho=rho, p=p,
                         n_distinct_profiles=int(len(np.unique(X[group].round(10), axis=0)))))
    return rows


def founder(cell):
    for name in ("AB", "MS", "E", "C", "D"):
        if cell.startswith(name):
            return name
    return "P4"


# ---- fronts ----------------------------------------------------------------------

def state_matrix(frame, features, terms, fate_groups=None):
    """Cosine cell-state cost; with ``fate_groups`` each child is its group's mean profile."""
    if fate_groups is None:
        return fs.state_matrix(frame, features, terms)
    children, parents = zip(*terms)
    P = frame.loc[list(parents), list(features)].to_numpy(float)
    C = frame.loc[list(children), list(features)]
    C = C.groupby(np.array([fate_groups[c] for c in children])).transform("mean")
    return cdist(P, C.to_numpy(float), "cosine")


CURVE_COLUMNS = ("alpha_travel", "travel", "cell_state", "D1", "D2", "u", "edge_retention")
FRONT_METRICS = ("u_L", "d_LP", "d_NP", "d_lineage_sd", "natural_state_z", "natural_travel_z", "max_edge_retention",
                 "retention_at_nearest", "exactly_on_front", "travel_sigma", "state_sigma", "natural_travel",
                 "natural_cell_state")


def front(inputs, cohort, geometry, state):
    row, curve = fs.evaluate(inputs, None, cohort, geometry, "representation", state=state)
    return {k: row[k] for k in FRONT_METRICS}, curve[list(CURVE_COLUMNS)]


# ---- inputs ----------------------------------------------------------------------

@dataclass
class Context:
    inputs: object
    frames: dict       # frame key -> expression table
    candidates: dict   # modality -> {candidate: features}
    designs: dict      # dataset -> Design
    coarse: dict
    fine: dict
    fate_groups: dict  # terminal cell -> fine type, or "Schema D: <class>" without a description
    lineage: Lineage
    sources: dict


def load_context():
    inputs = fs.load_inputs()
    ranking = load_ranking_module()
    coarse, fine = cell_types(ranking)
    lineage = build_lineage(dl.load_json(str(ROOT / "data/cell_lineage.json")), ranking.map_names)
    cross = inputs.cohorts["cross187"]
    with np.load(fs.CROSS_ANALYSIS / "assignments_and_costs.npz", allow_pickle=False) as arrays:  # hash-checked by fs
        cb_travel = {g: np.asarray(arrays[f"{g}__cb_rna__travel"]) for g in ("raw3d", "xy")}
        inputs.published_states.update({("cross187", g, "cb_rna"): np.asarray(arrays[f"{g}__cb_rna__state"])
                                        for g in ("raw3d", "xy")})
    inputs.cohorts["cross187_cb"] = fs.Cohort(cross.terms, cross.parents, cross.groups, cb_travel)
    frames = {"ce_protein": inputs.frames["ce_protein"], "ce_rna": inputs.frames["ce_rna"],
              "cb_rna": dl.load_cb_rna(str(ROOT / CB_RNA))}
    cross_cells = sorted({name for edge in cross.terms for name in edge})
    rna_pool = fs.variability_scores({"ce": frames["ce_rna"], "cb": frames["cb_rna"]}, inputs.rankings["rna"], cross_cells)
    protein_pool = inputs.hv_scores["protein"]
    candidates = {
        "protein": {"top-20": inputs.rankings["protein"][:20],
                    "HV pool": protein_pool.loc[protein_pool.in_hv_pool, "feature"].tolist(),
                    "whole panel": list(inputs.rankings["protein"])},
        "rna": {"top-20": inputs.rankings["rna"][:20],
                "HV pool": rna_pool.loc[rna_pool.in_hv_pool, "feature"].tolist(),
                "whole panel": list(inputs.rankings["rna"])},
    }
    for modality, published in inputs.published_features.items():
        if candidates[modality]["top-20"] != published:
            raise ValueError(f"{modality} top-20 differs from the published set")
    designs = {name: build_design(inputs.cohorts[cohort].terms, coarse, fine, lineage)
               for name, (_, _, cohort) in DATASETS.items()}
    fate_groups = {cell: fine.get(cell, f"Schema D: {coarse[cell]}") for cell in coarse}
    sources = dict(inputs.sources)
    sources.update({name: file_hash(ROOT / name) for name in (*CODE_SOURCES, CELL_KEY, CB_RNA)})
    return Context(inputs, frames, candidates, designs, coarse, fine, fate_groups, lineage, sources)


def matrix(ctx, dataset, features):
    _, frame, _ = DATASETS[dataset]
    return ctx.frames[frame].loc[ctx.designs[dataset].cells, list(features)].to_numpy(float)


# ---- reproduction gates --------------------------------------------------------------

def check_reproduction(ctx):
    """Measured top-20 fronts reproduce the published caches; protein sibling separation reproduces the notebook."""
    published, _ = fs.check_published(ctx.inputs)
    report = dict(feature_sensitivity=published, cb_rna={}, notebook={})
    failures = list(published["failures"])
    cached = pd.read_csv(fs.CROSS_ANALYSIS / "metrics.csv")
    cached = cached[(cached.cohort == "base") & (cached.config == "cb_rna")].set_index("geometry")
    for geometry in ("raw3d", "xy"):
        state = state_matrix(ctx.frames["cb_rna"], ctx.candidates["rna"]["top-20"], ctx.inputs.cohorts["cross187_cb"].terms)
        delta = float(np.max(np.abs(state - ctx.inputs.published_states["cross187", geometry, "cb_rna"])))
        report["cb_rna"][f"{geometry}/state_matrix"] = delta
        if delta > 1e-12:
            failures.append(f"CB RNA state matrix {geometry}: {delta}")
        row, _ = front(ctx.inputs, "cross187_cb", geometry, state)
        for name in ("u_L", "d_LP", "d_NP", "max_edge_retention", "retention_at_nearest", "travel_sigma", "state_sigma"):
            gap = float(abs(row[name] - cached.at[geometry, name]))
            report["cb_rna"][f"{geometry}/{name}"] = dict(actual=float(row[name]), expected=float(cached.at[geometry, name]),
                                                         abs_delta=gap)
            if not gap <= REPRODUCTION_TOL:
                failures.append(f"CB RNA {geometry} {name}: {row[name]!r} vs {cached.at[geometry, name]!r}")
    # Notebook representation: per-cell z-score across proteins, then per-protein standardization, all terminals.
    from sklearn.preprocessing import StandardScaler
    raw = pd.read_csv(ROOT / fs.PROTEIN_RAW, index_col=0).T.fillna(0)
    raw = raw[~(raw == 0).all(axis=1)]
    z = raw.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    X = pd.DataFrame(StandardScaler().fit_transform(z.values.T.astype(float)),
                     index=[load_ranking_module().map_names(n) for n in z.columns], columns=z.index)
    terminals = [c for c in X.index if c in ctx.coarse]
    families = {}
    for cell in terminals:
        families.setdefault(ctx.lineage.ancestors[cell][-1], []).append(cell)
    pairs = [(a, b) for kids in families.values() for i, a in enumerate(kids) for b in kids[i + 1:]]
    top = ctx.candidates["protein"]["top-20"]
    for metric, expected in NOTEBOOK_SIBLING_D.items():
        values = {same: np.array([cdist(X.loc[[a], top], X.loc[[b], top], metric)[0, 0] for a, b in pairs
                                  if (ctx.coarse[a] == ctx.coarse[b]) == same]) for same in (True, False)}
        pooled = np.sqrt((values[True].var(ddof=1) + values[False].var(ddof=1)) / 2)
        d = float((values[False].mean() - values[True].mean()) / pooled)
        report["notebook"][metric] = dict(actual=d, expected=expected, n_same=len(values[True]), n_diff=len(values[False]))
        if not abs(d - expected) <= NOTEBOOK_TOL:
            failures.append(f"notebook sibling separation {metric}: {d} vs {expected}")
    report["failures"] = failures
    return report


# ---- computation -----------------------------------------------------------------------

def candidate_table(ctx):
    rows = []
    for dataset, (modality, _, cohort) in DATASETS.items():
        for c_index, (candidate, features) in enumerate(ctx.candidates[modality].items()):
            rows.append(dict(dataset=dataset, modality=modality, cohort=cohort, candidate=candidate,
                             n_features=len(features), features=";".join(features)))
    return pd.DataFrame(rows)


def draw_random_sets(ctx, *, n=N_RANDOM, k=RANDOM_K, seed=SEED, max_rejections=10_000):
    """Random k-sets from each modality's HV pool; sets leaving a zero vector in any dataset are redrawn."""
    rows = []
    for m_index, modality in enumerate(("protein", "rna")):
        pool = ctx.candidates[modality]["HV pool"]
        datasets = [d for d, (m, _, _) in DATASETS.items() if m == modality]
        rng = np.random.default_rng([seed, 99, m_index, k])
        for draw in range(n):
            for rejected in range(max_rejections + 1):
                chosen = [pool[i] for i in rng.choice(len(pool), size=k, replace=False)]
                if all((np.linalg.norm(matrix(ctx, d, chosen), axis=1) > 0).all() for d in datasets):
                    break
            else:
                raise ValueError(f"Could not draw a valid {modality} set")
            rows.append(dict(modality=modality, draw=draw, rejected_draws=rejected, features=";".join(chosen)))
    return pd.DataFrame(rows)


def _criteria_job(ctx, dataset, candidate, seed, n_boot, n_perm):
    modality = DATASETS[dataset][0]
    row, curve = criteria(matrix(ctx, dataset, ctx.candidates[modality][candidate]), ctx.designs[dataset],
                          seed=seed, n_boot=n_boot, n_perm=n_perm)
    return dict(dataset=dataset, candidate=candidate, **row), [dict(dataset=dataset, candidate=candidate, **c) for c in curve]


def _random_job(ctx, chunk):
    rows = []
    for modality, draw, features in chunk:
        for dataset in (d for d, (m, _, _) in DATASETS.items() if m == modality):
            rows.append(dict(modality=modality, draw=draw, dataset=dataset,
                             **point_criteria(matrix(ctx, dataset, features.split(";")), ctx.designs[dataset])))
    return rows


def _front_job(ctx, dataset, geometry, candidate, cost):
    modality, frame, _ = DATASETS[dataset]
    cohort = FRONTS[dataset][0]
    state = state_matrix(ctx.frames[frame], ctx.candidates[modality][candidate], ctx.inputs.cohorts[cohort].terms,
                         ctx.fate_groups if cost == "fate-only" else None)
    row, curve = front(ctx.inputs, cohort, geometry, state)
    return dict(dataset=dataset, cohort=cohort, geometry=geometry, candidate=candidate, cost=cost, **row), curve


def compute(ctx, *, n_jobs, n_boot=N_BOOT, n_perm=N_PERM, n_perm_type=N_PERM_TYPE, n_random=N_RANDOM, seed=SEED):
    from joblib import Parallel, delayed

    jobs = [(d, c, [seed, di, ci]) for di, d in enumerate(DATASETS) for ci, c in enumerate(CANDIDATES)]
    results = Parallel(n_jobs=min(n_jobs, len(jobs)))(delayed(_criteria_job)(ctx, d, c, s, n_boot, n_perm) for d, c, s in jobs)
    table = pd.DataFrame([r for r, _ in results])
    curves = pd.DataFrame([c for _, cs in results for c in cs])

    random_sets = draw_random_sets(ctx, n=n_random, seed=seed) if n_random else pd.DataFrame()
    work = [] if random_sets.empty else list(random_sets[["modality", "draw", "features"]].itertuples(index=False, name=None))
    chunks = [work[i:i + 10] for i in range(0, len(work), 10)]
    random_metrics = pd.DataFrame([r for rows in Parallel(n_jobs=n_jobs)(delayed(_random_job)(ctx, ch) for ch in chunks)
                                   for r in rows])

    by_type = []
    for di, (dataset, (modality, _, _)) in enumerate(DATASETS.items()):
        rows = memory_by_type(matrix(ctx, dataset, ctx.candidates[modality]["top-20"]), ctx.designs[dataset], ctx.fine,
                              ctx.lineage, seed=[seed, 7, di], n_perm=n_perm_type)
        by_type += [dict(dataset=dataset, candidate="top-20", **r) for r in rows]

    front_jobs = [(d, g, c, cost) for d, (_, geometries) in FRONTS.items() for g in geometries
                  for c in CANDIDATES for cost in ("measured", "fate-only")]
    fronts = Parallel(n_jobs=n_jobs)(delayed(_front_job)(ctx, *job) for job in front_jobs)
    front_table = pd.DataFrame([r for r, _ in fronts])
    front_curves = {f"{d}__{g}__{c}__{cost}__{column}".replace(" ", "_"): curve[column].to_numpy(float)
                    for (d, g, c, cost), (_, curve) in zip(front_jobs, fronts) for column in curve.columns}
    return dict(criteria=table, memory_curves=curves, random_sets=random_sets, random_metrics=random_metrics,
                memory_by_type=pd.DataFrame(by_type), fronts=front_table, front_curves=front_curves)


def cell_type_table(ctx):
    cells = sorted({c for d in ctx.designs.values() for c in d.cells}, key=ctx.lineage.order.get)
    return pd.DataFrame(dict(cell=cells, founder=[founder(c) for c in cells], schema_d=[ctx.coarse[c] for c in cells],
                             fine_type=[ctx.fine.get(c, "") for c in cells], fate_group=[ctx.fate_groups[c] for c in cells],
                             in_protein_cohort=[c in set(ctx.designs["CE protein"].cells) for c in cells],
                             in_cross_species_cohort=[c in set(ctx.designs["CE RNA"].cells) for c in cells]))


def settings(n_boot=N_BOOT, n_perm=N_PERM, n_perm_type=N_PERM_TYPE, n_random=N_RANDOM, seed=SEED):
    return dict(seed=seed, n_boot=n_boot, n_perm=n_perm, n_perm_type=n_perm_type, n_random=n_random, random_k=RANDOM_K,
                min_type_cells=MIN_TYPE_CELLS, bin_min_pairs=BIN_MIN_PAIRS, intervals=INTERVALS,
                hv_fraction=fs.HV_FRACTION, rna_hv_pool="log1p_variance_mean_rank_ce_cb_cross187",
                distance="cosine", candidates=list(CANDIDATES), datasets={k: list(v) for k, v in DATASETS.items()},
                fate_only="children_as_fine_type_mean_parents_measured", numpy_version=np.__version__)


TABLES = ("candidates", "cell_types", "criteria", "memory_curves", "memory_by_type", "random_sets", "random_metrics", "fronts")


def write_run(out, *, n_jobs, **overrides):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    ctx = load_context()
    report = check_reproduction(ctx)
    if report["failures"]:
        raise ValueError("Reproduction gates failed:\n  " + "\n  ".join(report["failures"]))
    print("Published top-20 fronts and notebook sibling separation reproduced.", flush=True)
    results = compute(ctx, n_jobs=n_jobs, **overrides)
    analysis.mkdir(parents=True, exist_ok=True)
    tables = dict(candidates=candidate_table(ctx), cell_types=cell_type_table(ctx), **{
        k: results[k] for k in TABLES if k in results})
    for name in TABLES:
        tables[name].to_csv(analysis / f"{name}.csv", index=False)
    np.savez_compressed(analysis / "front_curves.npz", **results["front_curves"])
    (analysis / "reproduction.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    run_settings = settings(**overrides)
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=run_settings, sources=ctx.sources),
                                         sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=run_settings, sources=ctx.sources,
                  elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def _compare(fresh, saved, label, keys, tol=1e-12):
    fresh, saved = fresh.reset_index(drop=True), saved.reset_index(drop=True)
    if len(fresh) != len(saved):
        raise ValueError(f"Replay mismatch {label}: {len(fresh)} vs {len(saved)} rows")
    for column in saved.columns:
        if column in keys:
            if not (fresh[column].astype(str) == saved[column].astype(str)).all():
                raise ValueError(f"Replay mismatch {label} column {column}")
        elif not np.allclose(fresh[column].astype(float), saved[column].astype(float), rtol=tol, atol=tol, equal_nan=True):
            raise ValueError(f"Replay mismatch {label} column {column}")


def validate_run(out, *, n_jobs=1, check_sources=True, replay_random=20):
    """Check hashes, then recompute every table (random metrics on a seeded sample) and compare."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported representation-criteria cache")
    for name, expected in record["files"].items():
        if file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    s = record["settings"]
    overrides = dict(n_boot=s["n_boot"], n_perm=s["n_perm"], n_perm_type=s["n_perm_type"], n_random=s["n_random"],
                     seed=s["seed"])
    ctx = load_context()
    if check_reproduction(ctx)["failures"]:
        raise ValueError("Reproduction gates fail on current inputs")
    saved = {name: pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""]) for name in TABLES}
    keys = {"dataset", "modality", "cohort", "candidate", "features", "cell", "founder", "schema_d", "fine_type",
            "fate_group", "founders", "geometry", "cost", "in_protein_cohort", "in_cross_species_cohort", "exactly_on_front"}
    _compare(candidate_table(ctx), saved["candidates"], "candidates", keys)
    _compare(cell_type_table(ctx), saved["cell_types"].fillna(""), "cell_types", keys)
    _compare(draw_random_sets(ctx, n=s["n_random"], seed=s["seed"]), saved["random_sets"], "random_sets", keys | {"draw"})
    sample = saved["random_sets"].sample(min(replay_random, len(saved["random_sets"])), random_state=0)
    fresh = pd.DataFrame(_random_job(ctx, list(sample[["modality", "draw", "features"]].itertuples(index=False, name=None))))
    stored = saved["random_metrics"].merge(fresh[["modality", "draw", "dataset"]], on=["modality", "draw", "dataset"])
    _compare(fresh.sort_values(["modality", "draw", "dataset"]), stored.sort_values(["modality", "draw", "dataset"]),
             "random_metrics", keys)
    results = compute(ctx, n_jobs=n_jobs, **{**overrides, "n_random": 0})
    for name in ("criteria", "memory_curves", "memory_by_type", "fronts"):
        _compare(results[name], saved[name], name, keys)
    with np.load(analysis / "front_curves.npz") as arrays:
        if set(arrays.files) != set(results["front_curves"]):
            raise ValueError("Replay mismatch front curve keys")
        for key, values in results["front_curves"].items():
            if not np.allclose(arrays[key], values, rtol=1e-12, atol=1e-12):
                raise ValueError(f"Replay mismatch front curve {key}")
    return dict(analysis_id=record["analysis_id"], replayed_random=len(sample), criteria=len(results["criteria"]),
                fronts=len(results["fronts"]), verified_source_hashes=len(record["sources"]) if check_sources else 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--jobs", type=int, default=min(24, os.cpu_count() or 1))  # 24 physical cores
    parser.add_argument("--check-only", action="store_true", help="only run the reproduction gates")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.check_only:
        report = check_reproduction(load_context())
        print(json.dumps(report, indent=2, sort_keys=True, default=float))
        return 1 if report["failures"] else 0
    out = run_path(args.run_id)
    if not args.verify_only:
        record = write_run(out, n_jobs=args.jobs)
        print(json.dumps({k: record[k] for k in ("analysis_id", "elapsed_seconds")}, indent=2))
    print(json.dumps(validate_run(out, n_jobs=args.jobs), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
