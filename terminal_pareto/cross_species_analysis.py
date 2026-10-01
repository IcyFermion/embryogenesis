"""Pooled terminal comparison of CE protein, CE RNA, and CB AF16 RNA.

Run from the repository root in the dev environment. The published CE pooled
context supplies fixed global normalization totals; CB uses a separate fixed
two-AF16 reference cohort. Neither the publication nor archived pilot is written.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.analysis_context import build_analysis_context, validate_existing_context_manifest
from terminal_pareto.front_coordinates import EndpointTransform
from terminal_pareto.subtree_analysis import build_type_map, exact_cousin_stats

VERSION = "pooled-terminal-cross-species-v1"
DEFAULT_RUN = "pooled_comparison_20260930"
RUN_ROOT = ROOT / "terminal_pareto/output/runs/cross_species_terminal_v1"
PRIMARY = ("ce_protein", "ce_rna", "cb_rna")
LABELS = {"ce_protein": "C. elegans protein", "ce_rna": "C. elegans RNA", "cb_rna": "C. briggsae RNA"}
GEOMETRIES = ("raw3d", "xy")


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_path(run_id):
    if not run_id or Path(run_id).name != run_id or run_id in (".", ".."):
        raise ValueError("Run identity must be a single directory name")
    return RUN_ROOT / run_id


def pool_matrices(components, denominators):
    """Mean normalized candidate edge costs, with denominators fixed upstream."""
    if not components or set(components) != set(denominators):
        raise ValueError("Each component requires one fixed normalization total")
    values = []
    shape = next(iter(components.values())).shape
    for key, matrix in components.items():
        matrix = np.asarray(matrix, float)
        total = float(denominators[key])
        if matrix.shape != shape or not np.isfinite(matrix).all() or not np.isfinite(total) or total <= 0:
            raise ValueError("Invalid pooled matrix or normalization total")
        values.append(matrix / total)
    return np.mean(values, axis=0)


def endpoint_assignment(primary, secondary):
    """Checked secondary-objective tie handling from the archived pilot.

    Diminish perturbations after scaling matrices to unit maximum. Require the
    primary optimum within floating summation tolerance and secondary agreement
    at consecutive accepted perturbations; this is not exhaustive enumeration.
    """
    rows, ordinary = linear_sum_assignment(primary)
    p = primary / max(float(np.max(np.abs(primary))), 1e-15)
    s = secondary / max(float(np.max(np.abs(secondary))), 1e-15)
    best = np.sum(p[rows, ordinary], dtype=np.longdouble)
    tolerance = 32 * np.finfo(float).eps * len(primary)
    previous = None
    for power in range(4, 14):
        epsilon = 10.0 ** (-power)
        _, chosen = linear_sum_assignment(p + epsilon * s)
        gap = float(np.sum(p[rows, chosen], dtype=np.longdouble) - best)
        value = float(s[rows, chosen].sum())
        if abs(gap) <= tolerance:
            if previous is not None and abs(value - previous) <= tolerance:
                return chosen, dict(
                    primary_gap=float(primary[rows, chosen].sum() - primary[rows, ordinary].sum()),
                    secondary_saving=float(secondary[rows, ordinary].sum() - secondary[rows, chosen].sum()),
                    perturbation=epsilon, primary_tolerance_scaled=tolerance,
                )
            previous = value
        else:
            previous = None
    raise ValueError("Endpoint tie handling failed to preserve the primary optimum")


def spatial_matrix(xyz, terms, geometry):
    children, parents = zip(*terms)
    p = np.asarray([xyz[name] for name in parents], dtype=float)
    c = np.asarray([xyz[name] for name in children], dtype=float)
    if geometry == "xy":
        p, c = p[:, 1:], c[:, 1:]
    elif geometry != "raw3d":
        raise ValueError(f"Unknown geometry: {geometry}")
    if not np.isfinite(p).all() or not np.isfinite(c).all():
        raise ValueError("Nonfinite spatial coordinates")
    return cdist(p, c)


def expression_matrix(frame, features, terms):
    if len(features) != 20 or len(set(features)) != 20 or not frame.index.is_unique:
        raise ValueError("Expected unique identities and 20 unique molecular features")
    children, parents = zip(*terms)
    p = frame.loc[list(parents), features].to_numpy(float)
    c = frame.loc[list(children), features].to_numpy(float)
    for values in (p, c):
        if not np.isfinite(values).all() or not (np.linalg.norm(values, axis=1) > 0).all():
            raise ValueError("Selected molecular vectors must be finite and nonzero")
    return cdist(p, c, metric="cosine")


@dataclass
class Inputs:
    terms: list
    lineage: dict
    matrices: dict
    metadata: dict
    manifest: pd.DataFrame
    references: pd.DataFrame
    coverage: pd.DataFrame


def load_inputs():
    # Read and validate the adopted context without writing into its run.
    ce = build_analysis_context("pooled_tracking_v1", run_id="migration_candidate_20260920")
    validate_existing_context_manifest(ce)
    ce_rna, cb_rna, features = dl.load_ce_rna(), dl.load_cb_rna(), dl.load_rna_sel()
    metadata_path = ROOT / "data/c_briggsae/yiming/metadata.csv"
    metadata = pd.read_csv(metadata_path).set_index("ID")
    cb_xyz, cb_edges, cb_tracking = {}, {}, []
    for replicate in ("210519ZZY0874p1", "210519ZZY0874p5"):
        label, path, cutoff = next(row for row in dl.CB_REPLICATES_NEW if row[0].startswith(replicate))
        if metadata.loc[replicate, "Type"] != "AF16":
            raise ValueError(f"Expected AF16 tracking: {label}")
        xyz, valid = dl.load_briggsae_tracking(cutoff, path)
        children, parents = dl.collect_terminals(ce.lineage, set(valid) & set(cb_rna.index))
        cb_xyz[replicate] = xyz
        cb_edges[replicate] = set(zip(children, parents))
        cb_tracking.append(dict(label=replicate, path=str(Path(path).relative_to(ROOT)), cutoff=cutoff))
    cb_common = set.intersection(*cb_edges.values())
    ce_rna_common = {edge for edge in ce.ordered_edges if all(name in ce_rna.index for name in edge)}
    common = set(ce.ordered_edges) & ce_rna_common & cb_common
    terms = [edge for edge in ce.ordered_edges if edge in common]
    if len(terms) != 187:
        raise ValueError(f"Inputs changed: expected 187 common edges, found {len(terms)}")
    if len({child for child, _ in terms}) != len(terms):
        raise ValueError("Repeated terminal identity")
    # CB reference has the two-embryo/RNA edge intersection before CE matching.
    children, parents = dl.collect_terminals(ce.lineage, set(lm.build_lineage_tree_index(ce.lineage)))
    cb_terms = [edge for edge in zip(children, parents) if edge in cb_common]
    ce_indices = ce.edge_indices(terms)
    cb_lookup = {edge: index for index, edge in enumerate(cb_terms)}
    cb_indices = np.asarray([cb_lookup[edge] for edge in terms])
    selectors = {"ce": np.ix_(ce_indices, ce_indices), "cb": np.ix_(cb_indices, cb_indices)}
    state = {
        "ce_protein": ce.cost_matrices(terms)[1],
        "ce_rna": expression_matrix(ce_rna, features, terms),
        "cb_rna": expression_matrix(cb_rna, features, terms),
    }
    matrices, reference_rows, totals = {}, [], {}
    for geometry in GEOMETRIES:
        ce_components = (ce.component_travel_matrices if geometry == "raw3d" else {
            label: spatial_matrix(xyz, ce.ordered_edges, geometry)
            for label, xyz in ce.component_xyz_maps.items()
        })
        cb_components = {label: spatial_matrix(xyz, cb_terms, geometry) for label, xyz in cb_xyz.items()}
        ce_totals = (ce.component_natural_totals if geometry == "raw3d" else {
            label: float(np.trace(matrix)) for label, matrix in ce_components.items()
        })
        cb_totals = {label: float(np.trace(matrix)) for label, matrix in cb_components.items()}
        totals[geometry] = dict(ce=ce_totals, cb=cb_totals)
        ce_pooled = pool_matrices(ce_components, ce_totals)[selectors["ce"]]
        cb_pooled = pool_matrices(cb_components, cb_totals)[selectors["cb"]]
        if geometry == "raw3d":
            np.testing.assert_array_equal(ce_pooled, ce.cost_matrices(terms)[0])
        for config in PRIMARY:
            matrices[geometry, config] = (cb_pooled if config == "cb_rna" else ce_pooled, state[config])
        for config in ("ce_protein", "ce_rna"):
            for label, matrix in ce_components.items():
                matrices[geometry, f"{config}__{label}"] = (matrix[selectors["ce"]] / ce_totals[label], state[config])
        for label, matrix in cb_components.items():
            matrices[geometry, f"cb_rna__{label}"] = (matrix[selectors["cb"]] / cb_totals[label], state["cb_rna"])
        for species, ref_terms, components, denoms in (
            ("ce", ce.ordered_edges, ce_components, ce_totals),
            ("cb", cb_terms, cb_components, cb_totals),
        ):
            for label in components:
                reference_rows.extend(dict(species=species, geometry=geometry, replicate=label,
                    terminal=child, natural_parent=parent, normalization_total=denoms[label],
                    in_comparison=(child, parent) in common) for child, parent in ref_terms)
    type_map = build_type_map([child for child, _ in terms])
    annotation_path = ROOT / "data/c_briggsae/science.adu8249/shared/CellTable_20241027.txt"
    annotation = pd.read_csv(annotation_path, sep="\t", keep_default_na=False).set_index("Lineage")
    manifest = []
    for child, parent in terms:
        if child in annotation.index and annotation.loc[child, "Parent"] != parent:
            raise ValueError(f"Annotation/tree parent mismatch for {child}")
        manifest.append(dict(terminal=child, natural_parent=parent, cell_type=type_map[child],
            child_annotation_missing=annotation.loc[child, "Missing"] if child in annotation.index else "absent",
            parent_annotation_missing=annotation.loc[parent, "Missing"] if parent in annotation.index else "absent"))
    source_paths = list(ce.source_hashes) + [
        "data/c_briggsae/science.adu8249/c_elegans_tf.csv",
        "data/c_briggsae/science.adu8249/c_briggsae_tf.csv",
        "expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv",
        str(metadata_path.relative_to(ROOT)), str(annotation_path.relative_to(ROOT)),
    ] + [row["path"] for row in cb_tracking]
    source_paths += [f"terminal_pareto/{name}" for name in (
        "cross_species_analysis.py", "data_loader.py", "analysis_context.py", "front_coordinates.py",
        "pareto_engine.py", "subtree_analysis.py", "lineage_metrics.py")]
    info = dict(version=VERSION, ce_context_cache_key=ce.cache_key,
        cohort_size=len(terms), ce_reference_size=len(ce.ordered_edges), cb_reference_size=len(cb_terms),
        ordered_edges=terms, tracking=dict(ce=[dict(label=x.label, path=str(Path(x.path).relative_to(ROOT)),
            cutoff=x.cutoff) for x in ce.spec.tracking], cb=cb_tracking),
        normalization_totals=totals, features=dict(ce_protein=ce.prot_sel, ce_rna=features, cb_rna=features),
        sources={name: file_hash(ROOT / name) for name in sorted(set(source_paths))},
        limitations=["CB relative axial calibration remains unverified in existing 3D",
            "Cutoff developmental alignment and RNA measured/aggregated/imputed provenance remain unresolved",
            "Protein and RNA use different selected panels; modality is not isolated",
            "Tracking replicates share molecular matrices; their spread is not molecular replicate uncertainty",
            "Shared lineage annotation checks joins, not independent cross-species parentage",
            "Weighted sweeps are sampled optima; interior tied optima are not enumerated"])
    coverage = pd.DataFrame([dict(config=config, available=n, matched=len(terms)) for config, n in
                            (("ce_protein", len(ce.ordered_edges)), ("ce_rna", len(ce_rna_common)),
                             ("cb_rna", len(cb_terms)))])
    return Inputs(terms, ce.lineage, matrices, info, pd.DataFrame(manifest), pd.DataFrame(reference_rows), coverage)


def analyze(travel, state, parents, groups, *, intervals, identity, null_draws=1000, seed=42):
    mx, ms, vx, vs, covariance = exact_cousin_stats(travel, state, groups)
    if vx <= 0 or vs <= 0:
        raise ValueError("Both exact first-cousin null variances must be positive")
    sx, ss = np.sqrt(vx), np.sqrt(vs)
    rows = np.arange(len(parents))
    state_endpoint, state_tie = endpoint_assignment(state, travel)
    travel_endpoint, travel_tie = endpoint_assignment(travel, state)
    permutations = []
    for step in range(intervals + 1):
        alpha = step / intervals
        _, perm = linear_sum_assignment(alpha * travel / sx + (1 - alpha) * state / ss)
        if step == 0:
            perm = state_endpoint
        elif step == intervals:
            perm = travel_endpoint
        permutations.append(perm)
    permutations = np.asarray(permutations)
    t = travel[rows, permutations].sum(axis=1)
    s = state[rows, permutations].sum(axis=1)
    retention = (parents[None, :] == parents[permutations]).mean(axis=1)
    if np.any(np.diff(t) > 1e-8) or np.any(np.diff(s) < -1e-8):
        raise ValueError("Sweep cost ordering failed")
    transform = EndpointTransform.from_endpoints(reference_analysis_id=identity,
        travel_optimum_assignment_id=f"{identity}:sweep:{intervals}",
        state_optimum_assignment_id=f"{identity}:sweep:0",
        travel_optimum_costs=(t[-1], s[-1]), state_optimum_costs=(t[0], s[0]))
    x, y = transform.transform(t, s)
    natural = np.array([np.trace(travel), np.trace(state)])
    nx, ny = map(float, transform.transform(*natural))
    null_x, null_y = map(float, transform.transform(mx, ms))
    nearest = int(np.argmin(np.hypot(x - nx, y - ny)))
    length = np.r_[0.0, np.cumsum(np.hypot(np.diff(x[::-1]), np.diff(y[::-1])))]
    u = (length / length[-1])[::-1]
    candidates = np.flatnonzero(retention == retention.max())
    maximum = int(candidates[np.argmin(np.abs(u[candidates] - u[nearest]))])
    curve = pd.DataFrame(dict(sweep_index=rows if intervals + 1 == len(rows) else np.arange(intervals + 1),
        alpha_travel=np.linspace(0, 1, intervals + 1), travel=t, cell_state=s,
        D1=x, D2=y, u=u, edge_retention=retention,
        nearest_assignment=np.arange(intervals + 1) == nearest,
        maximum_retention=np.arange(intervals + 1) == maximum))
    metrics = dict(n=len(parents), intervals=intervals, d_LP=float(np.hypot(x[nearest] - nx, y[nearest] - ny)),
        d_NP=float(np.hypot(x[nearest] - null_x, y[nearest] - null_y)), u_L=float(u[nearest]),
        max_edge_retention=float(retention[maximum]), retention_at_nearest=float(retention[nearest]),
        retention_travel=float(retention[-1]), retention_state=float(retention[0]),
        natural_travel=float(natural[0]), natural_cell_state=float(natural[1]), natural_D1=nx, natural_D2=ny,
        null_mean_travel=mx, null_mean_state=ms, null_variance_travel=vx, null_variance_state=vs,
        null_covariance=covariance, travel_sigma=sx, state_sigma=ss,
        null_D1=null_x, null_D2=null_y, nearest_index=nearest, maximum_index=maximum,
        **{f"travel_endpoint_{key}": value for key, value in travel_tie.items()},
        **{f"state_endpoint_{key}": value for key, value in state_tie.items()})
    return dict(curve=curve, metrics=metrics, permutations=permutations, transform=transform)


def write_analysis(inputs, out, *, intervals=300, dense_intervals=1200):
    if intervals <= 0 or dense_intervals <= intervals or dense_intervals % intervals:
        raise ValueError("The dense grid must be a larger nested multiple of the base grid")
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use --render-only or a new --run-id")
    analysis.mkdir(parents=True, exist_ok=True)
    settings = dict(intervals=intervals, dense_intervals=dense_intervals, null_draws=1000, full_random_draws=2000, seed=42)
    identity = hashlib.sha256(json.dumps(dict(inputs=inputs.metadata, settings=settings), sort_keys=True).encode()).hexdigest()
    arrays = dict(terminals=np.asarray([c for c, _ in inputs.terms]), parents=np.asarray([p for _, p in inputs.terms]))
    parents = arrays["parents"]
    tree_index = lm.build_lineage_tree_index(inputs.lineage)
    groups = pe.build_ancestor_groups(arrays["terminals"], tree_index, 2)
    fronts, metrics, clouds, convergence, transforms = [], [], [], [], {}
    for (geometry, config), (travel, state) in inputs.matrices.items():
        key = f"{geometry}__{config}"
        result = analyze(travel, state, parents, groups, intervals=intervals, identity=f"{identity}:{key}")
        ident = dict(geometry=geometry, config=config, cohort="base")
        fronts.append(result["curve"].assign(**ident))
        metrics.append(dict(**ident, **result["metrics"]))
        for name, value in (("travel", travel), ("state", state), ("permutations", result["permutations"])):
            arrays[f"{key}__{name}"] = value
        transforms[key] = result["transform"].metadata()
        if config in PRIMARY:
            for family, steps, seed, draws in (("first_cousin", 2, 42, 1000),
                ("second_cousin", 3, 43, 1000), ("third_cousin", 4, 44, 1000), ("full_random", None, 45, 2000)):
                null_groups = ([list(range(len(parents)))] if steps is None else
                               pe.build_ancestor_groups(arrays["terminals"], tree_index, steps))
                rt, rs = pe.compute_group_shuffle_costs(travel, state, null_groups, n_random=draws, seed=seed)
                dx, dy = result["transform"].transform(rt, rs)
                clouds.append(pd.DataFrame(dict(travel=rt, cell_state=rs, D1=dx, D2=dy,
                    geometry=geometry, config=config, family=family)))
            dense = analyze(travel, state, parents, groups, intervals=dense_intervals, identity=f"{identity}:{key}:dense")
            stride = dense_intervals // intervals
            for column in ("travel", "cell_state"):
                np.testing.assert_allclose(result["curve"][column], dense["curve"][column].to_numpy()[::stride], rtol=1e-12, atol=1e-10)
            fronts.append(dense["curve"].assign(geometry=geometry, config=config, cohort="dense"))
            metrics.append(dict(geometry=geometry, config=config, cohort="dense", **dense["metrics"]))
            arrays[f"{key}__dense_permutations"] = dense["permutations"]
            convergence.append(dict(geometry=geometry, config=config, **{
                f"delta_{name}": dense["metrics"][name] - result["metrics"][name]
                for name in ("d_LP", "d_NP", "u_L", "max_edge_retention")}))
        print(key, {name: round(result["metrics"][name], 5) for name in ("d_LP", "d_NP", "max_edge_retention")}, flush=True)
    np.savez_compressed(analysis / "assignments_and_costs.npz", **arrays)
    pd.concat(fronts, ignore_index=True).to_csv(analysis / "fronts.csv", index=False)
    pd.DataFrame(metrics).to_csv(analysis / "metrics.csv", index=False)
    pd.concat(clouds, ignore_index=True).to_csv(analysis / "null_clouds.csv", index=False)
    pd.DataFrame(convergence).to_csv(analysis / "convergence.csv", index=False)
    inputs.manifest.to_csv(analysis / "cell_manifest.csv", index=False)
    inputs.references.to_csv(analysis / "normalization_reference.csv", index=False)
    inputs.coverage.to_csv(analysis / "coverage.csv", index=False)
    record = dict(**inputs.metadata, settings=settings, analysis_id=identity, endpoint_transforms=transforms,
        assignment_convention="parent slot i receives child permutation[i]; child j receives parents[argsort(permutation)[j]]",
        files={p.name: file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def validate_run(out, *, check_sources=True):
    """Replay every cached assignment and verify transforms, null draws and identity."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported comparison cache")
    for name, expected in record["files"].items():
        if file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    fronts, clouds = pd.read_csv(analysis / "fronts.csv"), pd.read_csv(analysis / "null_clouds.csv")
    metrics = pd.read_csv(analysis / "metrics.csv")
    count = 0
    with np.load(analysis / "assignments_and_costs.npz", allow_pickle=False) as arrays:
        parents, children = arrays["parents"], arrays["terminals"]
        rows = np.arange(len(parents))
        if len(parents) != record["cohort_size"] or len(set(children)) != len(children):
            raise ValueError("Invalid cached cohort")
        for (geometry, config, cohort), frame in fronts.groupby(["geometry", "config", "cohort"]):
            key = f"{geometry}__{config}"
            suffix = "dense_permutations" if cohort == "dense" else "permutations"
            perms = arrays[f"{key}__{suffix}"]
            frame = frame.sort_values("sweep_index")
            np.testing.assert_array_equal(np.sort(perms, axis=1), np.broadcast_to(rows, perms.shape))
            t, s = (arrays[f"{key}__{axis}"][rows, perms].sum(axis=1) for axis in ("travel", "state"))
            np.testing.assert_allclose(t, frame.travel, rtol=1e-12, atol=1e-10)
            np.testing.assert_allclose(s, frame.cell_state, rtol=1e-12, atol=1e-10)
            np.testing.assert_allclose((parents[None, :] == parents[perms]).mean(axis=1), frame.edge_retention, atol=1e-14, rtol=0)
            metadata = record["endpoint_transforms"][key]
            transform = EndpointTransform(**{k: v for k, v in metadata.items() if k not in ("display_mode", "clipping", "axis_limits")})
            x, y = transform.transform(t, s)
            np.testing.assert_allclose([x, y], [frame.D1, frame.D2], atol=1e-12, rtol=1e-12)
            np.testing.assert_allclose([[x[-1], y[-1]], [x[0], y[0]]], [[0, 1], [1, 0]], atol=1e-12)
            count += len(perms)
            if config in PRIMARY and cohort == "base":
                null = clouds[(clouds.geometry == geometry) & (clouds.config == config)]
                qx, qy = transform.transform(null.travel.to_numpy(), null.cell_state.to_numpy())
                np.testing.assert_allclose([qx, qy], [null.D1, null.D2], atol=1e-12, rtol=1e-12)
                m = metrics[(metrics.geometry == geometry) & (metrics.config == config) & (metrics.cohort == cohort)].iloc[0]
                nearest = int(np.argmin(np.hypot(x - m.natural_D1, y - m.natural_D2)))
                np.testing.assert_allclose(np.hypot(x[nearest] - m.natural_D1, y[nearest] - m.natural_D2), m.d_LP, atol=1e-12)
        for geometry in GEOMETRIES:
            np.testing.assert_array_equal(arrays[f"{geometry}__ce_protein__travel"], arrays[f"{geometry}__ce_rna__travel"])
            for config, replicas in (("ce_protein", [x["label"] for x in record["tracking"]["ce"]]),
                                    ("cb_rna", [x["label"] for x in record["tracking"]["cb"]])):
                pooled = np.mean([arrays[f"{geometry}__{config}__{r}__travel"] for r in replicas], axis=0)
                np.testing.assert_allclose(pooled, arrays[f"{geometry}__{config}__travel"], rtol=1e-14, atol=1e-14)
    return dict(verified_assignments=count, verified_source_hashes=len(record["sources"]),
                cohort_size=record["cohort_size"], analysis_id=record["analysis_id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--intervals", type=int, default=300)
    parser.add_argument("--dense-intervals", type=int, default=1200)
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    out = run_path(args.run_id)
    if not args.render_only and not args.verify_only:
        write_analysis(load_inputs(), out, intervals=args.intervals, dense_intervals=args.dense_intervals)
    report = validate_run(out)
    print(json.dumps(report, indent=2), flush=True)
    if not args.verify_only:
        from terminal_pareto.fig_terminal_cross_species import render
        render(out)
        validation = out / "validation"
        validation.mkdir(exist_ok=True)
        (validation / "replay.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"Comparison ready: {out / 'figures/terminal_cross_species_existing_3d.png'}")


if __name__ == "__main__":
    main()
