"""Expression-feature calibration of the terminal Pareto geometry (Track A1).

Is the front geometry (position u, distances d_LP and d_NP, maximum edge
retention) a property of the selected markers or of the cost structure? For
k in {5, 10, 20, 30} this compares the classifier's top-k features with random
feature sets of the same size, drawn from the full panel and from a highly
variable pool, on five terminal configurations. The classifier's bottom-k
features, over the whole ranking and within the highly variable pool, are
negative controls placed in the same random distributions:

* the primary C. elegans protein cohort (275 pooled edges), and
* the C. elegans RNA-matched cohort: the 217 of those edges whose cells all
  have RNA, with CE protein and CE RNA, each in existing-3D and XY travel.

C. briggsae is excluded at the author's request (2026-10-07): its tracking has
larger z-axis uncertainty, which would add variation unrelated to the features.
Without the C. briggsae match the RNA cohort grows from the 187 edges of
Figures 8-9 to 217, the maximum: every three-embryo tracked edge with CE RNA is
already among the 275.

Pooled 3D travel is read from the published primary cache (hash-checked) and
restricted, keeping its fixed 275-edge denominators. XY travel is recomputed
from the three tracking files as in ``cross_species_analysis`` (XY totals on
the 275 edges) and must equal the cached 187-edge XY matrix on those edges.
Only the cosine cell-state matrix changes per feature set. Every sweep is
``cross_species_analysis.analyze`` (exact first-cousin scaling, 301 weights,
checked endpoint ties). Before the survey, the published top-20 sets must
reproduce the primary and the 187-edge cached metric tables, and the primary
cohort is also checked against the Figure 5 landmark implementation.

Run from the repository root in the dev environment::

    python terminal_pareto/feature_sensitivity.py --run-id survey_20261014
    python terminal_pareto/feature_sensitivity.py --run-id survey_20261014 --verify-only

An existing completed run is never overwritten.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.analysis_context import _sha256_array
from terminal_pareto.cross_species_analysis import analyze, pool_matrices, spatial_matrix

VERSION = "terminal-feature-sensitivity-v1"
DEFAULT_RUN = "survey_20261014"
RUN_ROOT = ROOT / "terminal_pareto/output/runs/feature_sensitivity_v1"
PRIMARY_ANALYSIS = ROOT / "terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920/analysis"
CROSS_ANALYSIS = ROOT / "terminal_pareto/output/runs/cross_species_terminal_v1/pooled_comparison_20260930/analysis"

PROTEIN_ZSCORE = "data/protein/aggregated_all/s3_zscore.csv"
PROTEIN_RAW = "data/protein/aggregated_all/s3.csv"
PROTEIN_RANKING = "expression_embedding/results/elegans_protein_linear_baseline/protein_feature_ranking.csv"
PROTEIN_TOP20 = "expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv"
CE_RNA = "data/c_briggsae/science.adu8249/c_elegans_tf.csv"
RNA_COEFFICIENTS = "expression_embedding/results/cross_species_rna_linear/linear_coefficients.csv"
RNA_TOP20 = "expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv"
CODE_SOURCES = tuple(f"terminal_pareto/{name}" for name in (
    "feature_sensitivity.py", "cross_species_analysis.py", "subtree_analysis.py", "front_coordinates.py",
    "pareto_engine.py", "lineage_metrics.py", "data_loader.py", "analysis_context.py",
    "fig5_table1_ce_canonical_metrics.py"))

KS = (5, 10, 20, 30)
DRAWS = 200
HV_FRACTION = 0.5
SEED = 20261007
INTERVALS = 300
MODALITIES = ("protein", "rna")
# Fixed sets from the classifier ranking: the top k, the bottom k, and the
# bottom k among highly variable features (low rank without silent panels).
REFERENCE_KINDS = ("top", "bottom", "bottom_hv")
KINDS = (*REFERENCE_KINDS, "random_all", "random_hv")
# (cohort, geometry, config) -> modality. "primary" is the 275-edge Figure 2-6
# cohort, "rna_matched" its 217 RNA-covered edges, "cross187" the Figure 8-9
# cohort, used only to reproduce the published top-20 metrics.
CONFIGS = {
    ("primary", "raw3d", "ce_protein"): "protein",
    ("rna_matched", "raw3d", "ce_protein"): "protein",
    ("rna_matched", "xy", "ce_protein"): "protein",
    ("rna_matched", "raw3d", "ce_rna"): "rna",
    ("rna_matched", "xy", "ce_rna"): "rna",
}
REPRODUCED = {("primary", "raw3d", "ce_protein"): "protein",
              **{("cross187", g, c): m for g in ("raw3d", "xy") for c, m in (("ce_protein", "protein"), ("ce_rna", "rna"))}}
N_RNA_MATCHED = 217
ON_FRONT_TOL = 1e-6  # null-SD units, as D_LINEAGE_TOL in the Figure 5 analysis
REPRODUCTION_TOL = 1e-9


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_path(run_id):
    if not run_id or Path(run_id).name != run_id or run_id in (".", ".."):
        raise ValueError("Run identity must be a single directory name")
    return RUN_ROOT / run_id


class ZeroVectorError(ValueError):
    """A selected molecular vector is zero, so cosine distance is undefined."""


def state_matrix(frame, features, terms):
    """Cosine cell-state cost, parent slots in rows and children in columns."""
    if len(set(features)) != len(features) or not frame.index.is_unique:
        raise ValueError("Expected unique identities and unique molecular features")
    children, parents = zip(*terms)
    p = frame.loc[list(parents), list(features)].to_numpy(float)
    c = frame.loc[list(children), list(features)].to_numpy(float)
    for values in (p, c):
        if not np.isfinite(values).all():
            raise ValueError("Selected molecular vectors must be finite")
        if not (np.linalg.norm(values, axis=1) > 0).all():
            raise ZeroVectorError("Selected molecular vectors must be nonzero")
    return cdist(p, c, metric="cosine")


@dataclass
class Cohort:
    terms: list
    parents: np.ndarray
    groups: list
    travel: dict  # geometry -> pooled travel matrix

    @property
    def cells(self):
        return sorted({name for edge in self.terms for name in edge})


@dataclass
class Inputs:
    cohorts: dict
    frames: dict  # config -> expression table used for its state matrix
    rankings: dict  # modality -> features in classifier order
    hv_scores: dict  # modality -> per-feature variability table
    published_states: dict  # (cohort, geometry, config) -> cached top-20 state matrix
    published_features: dict  # modality -> published top-20 list
    sources: dict


def load_primary(tree_index):
    manifest = json.loads((PRIMARY_ANALYSIS / "analysis_manifest.json").read_text())
    with np.load(PRIMARY_ANALYSIS / "analysis_context.npz", allow_pickle=False) as arrays:
        travel = np.asarray(arrays["travel_matrix"])
        state = np.asarray(arrays["expression_matrix"])
        children, parents = arrays["children"].tolist(), arrays["natural_parents"].tolist()
    if _sha256_array(travel) != manifest["matrix_hashes"]["travel"] or \
            _sha256_array(state) != manifest["matrix_hashes"]["expression"]:
        raise ValueError("Primary context arrays do not match their manifest")
    for name, expected in manifest["source_hashes"].items():
        if file_hash(ROOT / name) != expected:
            raise ValueError(f"Primary context input changed: {name}")
    if manifest["profile"] != "pooled_tracking_v1" or len(children) != 275:
        raise ValueError("Unexpected primary context")
    terms = list(zip(children, parents))
    cohort = Cohort(terms, np.asarray(parents), pe.build_ancestor_groups(children, tree_index, 2), {"raw3d": travel})
    return cohort, state


def primary_xy_travel(primary):
    """Pooled XY travel on the 275 edges, as built by cross_species_analysis.load_inputs."""
    with np.load(PRIMARY_ANALYSIS / "analysis_context.npz", allow_pickle=False) as arrays:
        labels = arrays["component_labels"].tolist()
        cached = {label: np.asarray(arrays[f"travel_component_{label}"]) for label in labels}
    components_3d, components_xy = {}, {}
    for label, path, cutoff in dl.CE_REPLICATES:
        xyz, _ = dl.load_elegans_tracking(cutoff, path)
        xyz = {name: np.asarray(value, dtype=float) for name, value in xyz.items()}
        components_3d[label] = spatial_matrix(xyz, primary.terms, "raw3d")
        components_xy[label] = spatial_matrix(xyz, primary.terms, "xy")
    if set(labels) != set(components_3d):
        raise ValueError("Tracking replicates differ from the primary context")
    for label in labels:  # the coordinates read here are those behind the cached 3D matrices
        np.testing.assert_allclose(components_3d[label], cached[label], rtol=1e-12, atol=1e-12)
    totals = {label: float(np.trace(matrix)) for label, matrix in components_xy.items()}
    return pool_matrices(components_xy, totals)


def restrict(cohort, terms, tree_index):
    lookup = {edge: i for i, edge in enumerate(cohort.terms)}
    index = np.asarray([lookup[edge] for edge in terms])
    children = [child for child, _ in terms]
    return Cohort(list(terms), cohort.parents[index], pe.build_ancestor_groups(children, tree_index, 2),
                  {geometry: matrix[np.ix_(index, index)] for geometry, matrix in cohort.travel.items()})


def load_cross187(tree_index):
    record = json.loads((CROSS_ANALYSIS / "provenance.json").read_text())
    path = CROSS_ANALYSIS / "assignments_and_costs.npz"
    if record["version"] != "pooled-terminal-cross-species-v1" or file_hash(path) != record["files"][path.name]:
        raise ValueError("Cross-species cache does not match its provenance")
    for name in (PROTEIN_ZSCORE, PROTEIN_TOP20, CE_RNA, RNA_TOP20, "data/cell_lineage.json"):
        if file_hash(ROOT / name) != record["sources"][name]:
            raise ValueError(f"Cross-species input changed: {name}")
    with np.load(path, allow_pickle=False) as arrays:
        children, parents = arrays["terminals"].tolist(), arrays["parents"].tolist()
        travel, states = {}, {}
        for geometry in ("raw3d", "xy"):
            # CE protein and RNA share pooled travel (checked by the cache's validate_run).
            travel[geometry] = np.asarray(arrays[f"{geometry}__ce_protein__travel"])
            for config in ("ce_protein", "ce_rna"):
                states["cross187", geometry, config] = np.asarray(arrays[f"{geometry}__{config}__state"])
    if [list(edge) for edge in zip(children, parents)] != record["ordered_edges"] or len(children) != 187:
        raise ValueError("Unexpected cross-species cohort")
    cohort = Cohort(list(zip(children, parents)), np.asarray(parents),
                    pe.build_ancestor_groups(children, tree_index, 2), travel)
    return cohort, states


def variability_scores(frames, features, cells):
    """Rank features by variance of log1p values over cohort cells (mean rank if several frames).

    Variance on the z-scored protein table is one by construction, so the raw
    (non-negative) intensities and stored TF values are used. Ties break by name.
    """
    table = pd.DataFrame(index=pd.Index(features, name="feature"))
    for label, frame in frames.items():
        values = frame.loc[cells, features].to_numpy(float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"{label}: expected finite non-negative values")
        table[f"{label}_log1p_variance"] = np.log1p(values).var(axis=0, ddof=1)
        table[f"{label}_rank"] = table[f"{label}_log1p_variance"].rank(ascending=False, method="average")
    table["mean_rank"] = table[[f"{label}_rank" for label in frames]].mean(axis=1)
    table = table.reset_index().sort_values(["mean_rank", "feature"]).reset_index(drop=True)
    table["hv_rank"] = np.arange(1, len(table) + 1)
    table["in_hv_pool"] = table["hv_rank"] <= int(round(HV_FRACTION * len(table)))
    return table


def load_inputs():
    lineage = dl.load_json(str(ROOT / "data/cell_lineage.json"))
    tree_index = lm.build_lineage_tree_index(lineage)
    primary, primary_state = load_primary(tree_index)
    cross187, cross_states = load_cross187(tree_index)
    protein = dl.load_protein_expression(str(ROOT / PROTEIN_ZSCORE))
    ce_rna = dl.load_ce_rna(str(ROOT / CE_RNA))
    primary.travel["xy"] = primary_xy_travel(primary)
    np.testing.assert_allclose(restrict(primary, cross187.terms, tree_index).travel["xy"], cross187.travel["xy"],
                               rtol=1e-14, atol=1e-15)
    matched = restrict(primary, [edge for edge in primary.terms if all(name in ce_rna.index for name in edge)],
                       tree_index)
    if len(matched.terms) != N_RNA_MATCHED or not set(cross187.terms) <= set(matched.terms):
        raise ValueError(f"Expected {N_RNA_MATCHED} RNA-covered edges containing the 187, found {len(matched.terms)}")

    protein_ranking = pd.read_csv(ROOT / PROTEIN_RANKING)["protein"].tolist()
    coefficients = pd.read_csv(ROOT / RNA_COEFFICIENTS, index_col=0)
    importance = coefficients.abs().max(axis=1)
    rna_ranking = sorted(importance.index, key=lambda name: (-importance[name], name))
    published = {"protein": dl.load_prot_sel(str(ROOT / PROTEIN_TOP20)), "rna": dl.load_rna_sel(str(ROOT / RNA_TOP20))}
    for modality, ranking in (("protein", protein_ranking), ("rna", rna_ranking)):
        if ranking[:20] != published[modality]:
            raise ValueError(f"{modality} ranking does not reproduce the published top 20")
    if set(protein_ranking) != set(protein.columns) or len(protein_ranking) != 210:
        raise ValueError("Protein ranking must cover the 210 z-scored features")
    # The ranking keeps the 237 TFs shared with C. briggsae used to select the panel.
    if not set(rna_ranking) <= set(ce_rna.columns) or len(rna_ranking) != 237:
        raise ValueError("RNA ranking must cover the 237 shared TFs")

    raw_protein = pd.read_csv(ROOT / PROTEIN_RAW, index_col=0).fillna(0)
    hv_scores = {
        # One protein pool for both cohorts, defined on the larger primary cohort.
        "protein": variability_scores({"s3": raw_protein}, protein_ranking, primary.cells),
        "rna": variability_scores({"ce": ce_rna}, rna_ranking, matched.cells),
    }
    frames = {"ce_protein": protein, "ce_rna": ce_rna}
    published_states = {("primary", "raw3d", "ce_protein"): primary_state, **cross_states}
    sources = {name: file_hash(ROOT / name) for name in sorted({
        PROTEIN_ZSCORE, PROTEIN_RAW, PROTEIN_RANKING, PROTEIN_TOP20, CE_RNA, RNA_COEFFICIENTS,
        RNA_TOP20, "data/cell_lineage.json", *CODE_SOURCES, *(str(Path(path).relative_to(ROOT))
                                                              for _, path, _ in dl.CE_REPLICATES)})}
    for path in (PRIMARY_ANALYSIS / "analysis_manifest.json", PRIMARY_ANALYSIS / "analysis_context.npz",
                 PRIMARY_ANALYSIS / "ce_subtree_canonical_metrics.csv", CROSS_ANALYSIS / "provenance.json",
                 CROSS_ANALYSIS / "assignments_and_costs.npz", CROSS_ANALYSIS / "metrics.csv"):
        sources[str(path.relative_to(ROOT))] = file_hash(path)
    return Inputs({"primary": primary, "rna_matched": matched, "cross187": cross187}, frames,
                  {"protein": protein_ranking, "rna": rna_ranking}, hv_scores,
                  published_states, published, sources)


def feature_set_is_valid(inputs, modality, features):
    """Cosine needs a nonzero vector for every cell of every configuration."""
    for (cohort, geometry, config), kind in CONFIGS.items():
        if kind != modality or geometry != "raw3d":
            continue
        frame = inputs.frames[config]
        values = frame.loc[inputs.cohorts[cohort].cells, list(features)].to_numpy(float)
        if not (np.linalg.norm(values, axis=1) > 0).all():
            return False
    return True


def draw_feature_sets(inputs, *, ks=KS, draws=DRAWS, seed=SEED, max_rejections=10_000):
    """Reference sets plus seeded random draws; draws leaving a zero vector are redrawn and counted."""
    rows = []
    for m_index, modality in enumerate(MODALITIES):
        ranking = inputs.rankings[modality]
        position = {name: i for i, name in enumerate(ranking)}
        scores = inputs.hv_scores[modality]
        hv_pool = set(scores.loc[scores.in_hv_pool, "feature"])
        pools = {"random_all": ranking, "random_hv": [name for name in ranking if name in hv_pool]}
        for k in ks:
            top = ranking[:k]
            candidates = [("top", 0, None, 0, top), ("bottom", 0, None, 0, ranking[-k:]),
                          ("bottom_hv", 0, None, 0, pools["random_hv"][-k:])]
            for p_index, kind in enumerate(("random_all", "random_hv"), start=1):
                pool = pools[kind]
                if k >= len(pool):
                    raise ValueError(f"{modality} {kind} pool is too small for k={k}")
                generator_seed = [seed, m_index, p_index, k]
                rng = np.random.default_rng(generator_seed)
                for draw in range(draws):
                    rejected = 0
                    while True:
                        chosen = [pool[i] for i in rng.choice(len(pool), size=k, replace=False)]
                        if feature_set_is_valid(inputs, modality, chosen):
                            break
                        rejected += 1
                        if rejected > max_rejections:
                            raise ValueError(f"Could not draw a valid {modality} set of size {k}")
                    candidates.append((kind, draw, generator_seed, rejected, chosen))
            for kind, draw, generator_seed, rejected, features in candidates:
                if kind in REFERENCE_KINDS and not feature_set_is_valid(inputs, modality, features):
                    raise ValueError(f"{kind} {k} {modality} set has a zero vector")
                features = sorted(features, key=position.get)
                rows.append(dict(
                    set_id=f"{modality}__{kind}__k{k:02d}__{draw:03d}", modality=modality, kind=kind, k=k,
                    draw=draw, generator_seed="" if generator_seed is None else "-".join(map(str, generator_seed)),
                    rejected_draws=rejected, n_in_hv_pool=sum(name in hv_pool for name in features),
                    n_in_top_k=len(set(features) & set(top)),
                    best_classifier_rank=min(position[name] for name in features) + 1,
                    features=";".join(features)))
    return pd.DataFrame(rows)


def evaluate(inputs, features, cohort, geometry, config, *, intervals=INTERVALS, state=None):
    """Canonical metrics for one feature set and configuration (plus the sweep curve)."""
    data = inputs.cohorts[cohort]
    travel = data.travel[geometry]
    if state is None:
        state = state_matrix(inputs.frames[config], features, data.terms)
    result = analyze(travel, state, data.parents, data.groups, intervals=intervals,
                     identity=f"{VERSION}:{cohort}:{geometry}:{config}")
    m, curve = result["metrics"], result["curve"]
    sx, ss = m["travel_sigma"], m["state_sigma"]
    t, s = curve["travel"].to_numpy(), curve["cell_state"].to_numpy()
    d_lineage_sd = float(np.min(np.hypot((t - m["natural_travel"]) / sx, (s - m["natural_cell_state"]) / ss)))
    travel_span, state_span = float(t[0] - t[-1]), float(s[-1] - s[0])
    keep = ("n", "u_L", "d_LP", "d_NP", "max_edge_retention", "retention_at_nearest", "retention_travel",
            "retention_state", "natural_travel", "natural_cell_state", "natural_D1", "natural_D2",
            "null_mean_travel", "null_mean_state", "null_D1", "null_D2", "travel_sigma", "state_sigma",
            "null_covariance", "nearest_index", "maximum_index")
    row = {key: m[key] for key in keep}
    row.update(
        u_max_retention=float(curve["u"].iloc[m["maximum_index"]]),
        n_max_retention_ties=int((curve["edge_retention"] == m["max_edge_retention"]).sum()),
        natural_travel_z=(m["natural_travel"] - m["null_mean_travel"]) / sx,
        natural_state_z=(m["natural_cell_state"] - m["null_mean_state"]) / ss,
        travel_span=travel_span, state_span=state_span,
        travel_span_sd=travel_span / sx, state_span_sd=state_span / ss,
        d_lineage_sd=d_lineage_sd, exactly_on_front=bool(d_lineage_sd <= ON_FRONT_TOL),
        n_distinct_front_points=int(len(np.unique(np.round(np.c_[t, s], 12), axis=0))),
    )
    return row, curve


def check_published(inputs):
    """Published top-20 sets must reproduce both cached metric tables before any survey."""
    from terminal_pareto.fig5_table1_ce_canonical_metrics import compute_front_landmarks
    from terminal_pareto.subtree_analysis import first_cousin_null_summary

    report, failures = dict(state_matrices={}, primary={}, cross187={}), []

    def compare(label, actual, expected, tol=REPRODUCTION_TOL):
        delta = float(abs(actual - expected))
        if not delta <= tol:
            failures.append(f"{label}: {actual!r} vs {expected!r}")
        return dict(actual=float(actual), expected=float(expected), abs_delta=delta)

    rows = {}
    for key, modality in REPRODUCED.items():
        cohort, geometry, config = key
        state = state_matrix(inputs.frames[config], inputs.published_features[modality], inputs.cohorts[cohort].terms)
        delta = float(np.max(np.abs(state - inputs.published_states[key])))
        report["state_matrices"]["/".join(key)] = delta
        if delta > 1e-12:
            failures.append(f"state matrix {key}: max |delta| {delta}")
        rows[key] = evaluate(inputs, None, *key, state=state)[0]

    # Primary cohort: Figure 5 canonical P0 row and the Figure 5 implementation.
    canonical = pd.read_csv(PRIMARY_ANALYSIS / "ce_subtree_canonical_metrics.csv").set_index("subtree").loc["P0"]
    row = rows["primary", "raw3d", "ce_protein"]
    data = inputs.cohorts["primary"]
    travel, state = data.travel["raw3d"], inputs.published_states["primary", "raw3d", "ce_protein"]
    nodes = [child for child, _ in data.terms]
    null = first_cousin_null_summary(travel, state, nodes, pe.build_grandparent_map(
        dl.load_json(str(ROOT / "data/cell_lineage.json"))), seed=42)
    xa, ea, edge, _ = pe.compute_std_scaled_pareto(travel, state, list(data.parents), null["_raw"], iteration=INTERVALS)
    landmarks = compute_front_landmarks(xa, ea, edge, nx=null["null_xyz_mean_sigma"], ny=null["null_exp_mean_sigma"])
    for name, ours, cached_name in (("u", "u_L", "u_lineage_lp"), ("d_LP", "d_LP", "d_lp"),
                                    ("d_NP", "d_NP", "d_np"), ("max_retention", "max_edge_retention", "max_er"),
                                    ("natural_D1", "natural_D1", "D1_lineage"), ("natural_D2", "natural_D2", "D2_lineage")):
        report["primary"][f"{name}_vs_cache"] = compare(f"primary {name} vs cache", row[ours], canonical[cached_name])
        report["primary"][f"{name}_vs_fig5"] = compare(f"primary {name} vs fig5", row[ours], landmarks[cached_name])
    report["primary"]["exactly_on_front_vs_cache"] = dict(actual=row["exactly_on_front"],
                                                         expected=bool(canonical["exactly_on_front"]))
    if row["exactly_on_front"] != bool(canonical["exactly_on_front"]):
        failures.append("primary exactly_on_front")

    cached = pd.read_csv(CROSS_ANALYSIS / "metrics.csv")
    cached = cached[cached.cohort == "base"].set_index(["geometry", "config"])
    for (cohort, geometry, config), row in rows.items():
        if cohort != "cross187":
            continue
        expected = cached.loc[(geometry, config)]
        for name in ("u_L", "d_LP", "d_NP", "max_edge_retention", "retention_at_nearest",
                     "travel_sigma", "state_sigma"):
            report["cross187"][f"{geometry}/{config}/{name}"] = compare(
                f"cross187 {geometry}/{config} {name}", row[name], expected[name])
    report["failures"] = failures
    return report, rows


def _evaluate_chunk(inputs, chunk):
    rows = []
    for set_id, modality, features in chunk:
        for key, kind in CONFIGS.items():
            if kind != modality:
                continue
            row = dict(set_id=set_id, cohort=key[0], geometry=key[1], config=key[2], status="ok")
            try:
                row.update(evaluate(inputs, features.split(";"), *key)[0])
            except ValueError as error:
                # Mostly silent panels give (nearly) identical vectors for every cell, so the
                # cell-state cost has no first-cousin variance and the front is undefined.
                degenerate = "null variances must be positive" in str(error)
                row["status"] = "degenerate_state" if degenerate else f"failed: {error}"
            rows.append(row)
    return rows


def run_survey(inputs, feature_sets, *, n_jobs, chunk_size=10):
    from joblib import Parallel, delayed

    work = list(feature_sets[["set_id", "modality", "features"]].itertuples(index=False, name=None))
    chunks = [work[i:i + chunk_size] for i in range(0, len(work), chunk_size)]
    results = Parallel(n_jobs=n_jobs, verbose=5)(delayed(_evaluate_chunk)(inputs, chunk) for chunk in chunks)
    return pd.DataFrame([row for chunk in results for row in chunk])


SUMMARY_METRICS = ("u_L", "d_LP", "d_NP", "max_edge_retention", "retention_at_nearest", "u_max_retention",
                   "natural_travel_z", "natural_state_z", "travel_span_sd", "state_span_sd")


def summarize(feature_sets, metrics):
    """Per configuration, k, reference set and metric: its value against both random distributions.

    ``*_percentile`` is the share of random sets below the reference value plus
    half the ties; degenerate sets are counted, not summarized (a degenerate
    reference has value NaN).
    """
    frame = metrics.merge(feature_sets[["set_id", "kind", "k"]], on="set_id")
    rows = []
    for (cohort, geometry, config, k), group in frame.groupby(["cohort", "geometry", "config", "k"], sort=False):
        ok = group[group.status == "ok"]
        for reference in REFERENCE_KINDS:
            if not (group.kind == reference).any():
                continue
            chosen = ok[ok.kind == reference]
            for metric in SUMMARY_METRICS:
                value = float(chosen[metric].iloc[0]) if len(chosen) else np.nan
                row = dict(cohort=cohort, geometry=geometry, config=config, k=k, reference=reference, metric=metric,
                           value=value)
                for kind in ("random_all", "random_hv"):
                    values = ok.loc[ok.kind == kind, metric].to_numpy(float)
                    row[f"{kind}_n"] = len(values)
                    row[f"{kind}_degenerate"] = int(((group.kind == kind) & (group.status == "degenerate_state")).sum())
                    for name, q in (("q05", .05), ("q25", .25), ("median", .5), ("q75", .75), ("q95", .95)):
                        row[f"{kind}_{name}"] = float(np.quantile(values, q)) if len(values) else np.nan
                    row[f"{kind}_percentile"] = (100 * (np.mean(values < value) + .5 * np.mean(values == value))
                                                 if len(values) and np.isfinite(value) else np.nan)
                rows.append(row)
    return pd.DataFrame(rows)


def settings(ks=KS, draws=DRAWS, seed=SEED):
    return dict(ks=list(ks), draws=draws, seed=seed, kinds=list(KINDS), hv_fraction=HV_FRACTION, hv_score="log1p_variance_mean_rank",
                intervals=INTERVALS, configs=["/".join(key) for key in CONFIGS], on_front_tol=ON_FRONT_TOL,
                zero_vector_policy="redraw_and_count", numpy_version=np.__version__)


def write_run(out, *, n_jobs, ks=KS, draws=DRAWS, seed=SEED):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    inputs = load_inputs()
    report, _ = check_published(inputs)
    if report["failures"]:
        raise ValueError("Published top-20 sets do not reproduce the caches:\n  " + "\n  ".join(report["failures"]))
    print(f"Published top-20 metrics reproduced for all {len(REPRODUCED)} cached configurations.", flush=True)
    feature_sets = draw_feature_sets(inputs, ks=ks, draws=draws, seed=seed)
    metrics = run_survey(inputs, feature_sets, n_jobs=n_jobs)
    analysis.mkdir(parents=True, exist_ok=True)
    feature_sets.to_csv(analysis / "feature_sets.csv", index=False)
    metrics.to_csv(analysis / "metrics.csv", index=False)
    pools = pd.concat([table.assign(modality=modality) for modality, table in inputs.hv_scores.items()],
                      ignore_index=True)
    pools.to_csv(analysis / "hv_pools.csv", index=False)
    summarize(feature_sets, metrics).to_csv(analysis / "summary.csv", index=False)
    (analysis / "reproduction.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    run_settings = settings(ks, draws, seed)
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=run_settings, sources=inputs.sources),
                                         sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=run_settings, sources=inputs.sources,
                  n_feature_sets=len(feature_sets), n_evaluations=len(metrics),
                  n_degenerate=int((metrics.status == "degenerate_state").sum()),
                  n_failed=int(metrics.status.str.startswith("failed").sum()), elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def validate_run(out, *, replay=25, check_sources=True):
    """Check hashes, then replay every reference row and a seeded sample of random rows."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported feature-sensitivity cache")
    for name, expected in record["files"].items():
        if file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    inputs = load_inputs()
    run_settings = record["settings"]
    feature_sets = pd.read_csv(analysis / "feature_sets.csv", keep_default_na=False)
    redrawn = draw_feature_sets(inputs, ks=run_settings["ks"], draws=run_settings["draws"], seed=run_settings["seed"])
    pd.testing.assert_frame_equal(feature_sets, redrawn, check_dtype=False)
    metrics = pd.read_csv(analysis / "metrics.csv")
    pd.testing.assert_frame_equal(pd.read_csv(analysis / "summary.csv"), summarize(feature_sets, metrics))
    ok = metrics[metrics.status == "ok"]
    sample = ok[ok.set_id.str.split("__").str[1].isin(REFERENCE_KINDS)]
    sample = pd.concat([sample, ok[~ok.index.isin(sample.index)].sample(replay, random_state=0)])
    features = feature_sets.set_index("set_id")["features"]
    numeric = [c for c in metrics.columns if c not in ("set_id", "cohort", "geometry", "config", "status")]
    for _, row in sample.iterrows():
        fresh, _ = evaluate(inputs, features[row.set_id].split(";"), row.cohort, row.geometry, row.config)
        for name in numeric:
            if not np.isclose(float(fresh[name]), float(row[name]), rtol=1e-12, atol=1e-12):
                raise ValueError(f"Replay mismatch {row.set_id} {row.cohort}/{row.geometry}/{row.config} {name}")
    return dict(analysis_id=record["analysis_id"], feature_sets=len(feature_sets), evaluations=len(metrics),
                replayed=len(sample), verified_source_hashes=len(record["sources"]) if check_sources else 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--jobs", type=int, default=min(24, os.cpu_count() or 1))  # 24 physical cores
    parser.add_argument("--draws", type=int, default=DRAWS)
    parser.add_argument("--check-only", action="store_true", help="only reproduce the published top-20 metrics")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.check_only:
        report, _ = check_published(load_inputs())
        print(json.dumps(report, indent=2, sort_keys=True))
        return 1 if report["failures"] else 0
    out = run_path(args.run_id)
    if not args.verify_only:
        record = write_run(out, n_jobs=args.jobs, draws=args.draws)
        print(json.dumps({k: record[k] for k in ("analysis_id", "n_feature_sets", "n_evaluations", "n_degenerate", "n_failed",
                                                 "elapsed_seconds")}, indent=2))
    print(json.dumps(validate_run(out), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
