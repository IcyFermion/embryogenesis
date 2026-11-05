"""Terminal-pipeline cache adapters (read-only; never run solvers or nulls)."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput
from publication.provenance import sha256
from terminal_pareto import cross_species_analysis as tcs
from terminal_pareto.front_coordinates import EndpointTransform

CROSS_SPECIES_FAMILY = "terminal-cross-species"
CROSS_SPECIES_INPUTS = ("provenance.json", "fronts.csv", "metrics.csv", "null_clouds.csv",
                        "assignments_and_costs.npz")
CROSS_SPECIES_REFERENCES = dict(
    main=("first_cousin", "second_cousin", "third_cousin"), inset="full_random",
    inset_title="Full-random shuffle", inset_legend="Full-random shuffle (insets)", display_draws=None)


def _run(run) -> Path:
    return Path(run) if isinstance(run, Path) else tcs.run_path(run)


def transform(metadata: dict) -> EndpointTransform:
    return EndpointTransform(**{k: v for k, v in metadata.items()
                                if k not in ("display_mode", "clipping", "axis_limits")})


def tracking_transfers(record: dict, analysis: Path) -> dict:
    """Individual-embryo optima re-scored under pooled 3D costs and pooled anchors (replay only)."""
    children = np.arange(record["cohort_size"])
    transfers = {}
    with np.load(analysis / "assignments_and_costs.npz", allow_pickle=False) as arrays:
        for config in tcs.PRIMARY:
            key = f"raw3d__{config}"
            endpoint = transform(record["endpoint_transforms"][key])
            species = "cb" if config == "cb_rna" else "ce"
            curves = []
            for j, replica in enumerate(record["tracking"][species]):
                perms = arrays[f"{key}__{replica['label']}__permutations"]
                t = arrays[f"{key}__travel"][children, perms].sum(axis=1)
                s = arrays[f"{key}__state"][children, perms].sum(axis=1)
                x, y = endpoint.transform(t, s)
                label = f"AF16 {'p1' if j == 0 else 'p5'}" if species == "cb" else f"Embryo {j + 1}"
                curves.append((label, np.asarray(x), np.asarray(y)))
            transfers[config] = curves
    return transfers


def cross_species(run=tcs.DEFAULT_RUN) -> FrontComparisonInput:
    """Pooled 187-edge terminal comparison; cosine molecular distance."""
    run = _run(run)
    if not (run / "analysis/provenance.json").is_file():
        raise FileNotFoundError(
            f"No terminal comparison cache at {run}; create it with "
            f"`python terminal_pareto/cross_species_analysis.py --run-id {run.name}`")
    report = tcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    primary = list(tcs.PRIMARY)
    fronts, metrics, clouds = (pd.read_csv(analysis / name) for name in ("fronts.csv", "metrics.csv", "null_clouds.csv"))
    fronts = fronts[(fronts.cohort == "base") & fronts.config.isin(primary)].reset_index(drop=True)
    metrics = (metrics[(metrics.cohort == "base") & metrics.config.isin(primary)]
               .rename(columns=CANONICAL_FIELDS).reset_index(drop=True))
    clouds = clouds[clouds.config.isin(primary)].reset_index(drop=True)
    return FrontComparisonInput(
        family=CROSS_SPECIES_FAMILY, run=run, analysis_id=report["analysis_id"],
        edges=int(record["cohort_size"]), configs=tuple(tcs.PRIMARY), geometries=tuple(tcs.GEOMETRIES),
        fronts=fronts, metrics=metrics, null_mean="analytic first-cousin",
        caption=dict(edges=int(record["cohort_size"]), cb_reference_size=int(record["cb_reference_size"])),
        validation=report,
        input_files={f"analysis/{name}": sha256(analysis / name) for name in CROSS_SPECIES_INPUTS},
        clouds=clouds, references=CROSS_SPECIES_REFERENCES,
        extras=dict(tracking_transfers=tracking_transfers(record, analysis)))


PRIMARY_FAMILY = "terminal-primary"
PRIMARY_PROFILE = "pooled_tracking_v1"
PRIMARY_RUN = "migration_candidate_20260920"
PRIMARY_MIN_CELLS = 12  # floor of the cached subtree analysis (its manifests and validators)
SUBTREE_FLOOR = 24  # subtrees shown in the paper: see _resolution (Track D6)
PRIMARY_INTERVALS = 300


def _within_type(analysis: Path, context):
    """Cached Figure 6B-C tables; endpoint landmarks re-derived from the cached aggregate."""
    fronts = pd.read_csv(analysis / "within_type_fronts.csv")
    summary = pd.read_csv(analysis / "within_type_summary.csv")
    aggregate = pd.read_csv(analysis / "type_preserving_aggregate_front.csv")
    gx, ge = aggregate.global_travel_sigma.to_numpy(), aggregate.global_state_sigma.to_numpy()
    rx, re = aggregate.restricted_travel_sigma.to_numpy(), aggregate.restricted_state_sigma.to_numpy()
    last = len(aggregate) - 1
    shared = EndpointTransform.from_endpoints(
        reference_analysis_id=context.cache_key,
        travel_optimum_assignment_id=f"unrestricted:sweep:{last}",
        state_optimum_assignment_id="unrestricted:sweep:0",
        travel_optimum_costs=(gx[-1], ge[-1]), state_optimum_costs=(gx[0], ge[0]))
    np.testing.assert_allclose(np.column_stack(shared.transform(gx, ge)),
                               aggregate[["global_endpoint_travel", "global_endpoint_state"]], atol=1e-12)
    np.testing.assert_allclose(np.column_stack(shared.transform(rx, re)),
                               aggregate[["restricted_endpoint_travel", "restricted_endpoint_state"]], atol=1e-12)
    metadata = shared.metadata(clipping=False)
    natural = shared.transform(np.asarray([0.0]), np.asarray([0.0]))
    endpoints = dict(
        travel_penalty_sigma=float(rx[-1] - gx[-1]), state_penalty_sigma=float(re[0] - ge[0]),
        global_travel_reduction_sigma=float(-gx[-1]), restricted_travel_reduction_sigma=float(-rx[-1]),
        global_state_reduction_sigma=float(-ge[0]), restricted_state_reduction_sigma=float(-re[0]),
        endpoint_reference=shared.reference_analysis_id,
        endpoint_natural_travel=float(natural[0][0]), endpoint_natural_state=float(natural[1][0]),
        endpoint_metadata=metadata)
    return dict(fronts=fronts, summary=summary, aggregate=aggregate, endpoints=endpoints)


def _subtree_map(context, analysis: Path, min_cells: int, floor: int):
    from publication.figures.terminal_subtree_map import collect_nodes, terminal_inorder
    from terminal_pareto import data_loader as dl
    from terminal_pareto import lineage_metrics as lm
    from terminal_pareto.subtree_analysis import build_type_map, load_validated_subtree_summary
    summary = load_validated_subtree_summary(analysis, min_cells, context)
    summary = summary[summary.n >= floor].reset_index(drop=True)
    on_front = set(summary[summary["natural_on_front"].astype(bool)]["subtree"])
    by_name = summary.set_index("subtree")
    nodes, edges, *_ = collect_nodes(
        context.lineage, lm.build_lineage_tree_index(context.lineage), set(summary["subtree"]),
        terminal_inorder(context.lineage, dl.map_names), build_type_map(context.terminal_nodes),
        set(context.terminal_nodes), dl.map_names)
    for name, node in nodes.items():
        node["on_front"] = name in on_front
        node["n"] = sum(node["composition"].values())
        node["max_er"] = float(by_name.loc[name, "max_er"])
    return summary, nodes, edges


S3_DIR = "s3_cross_geometry"


def _s3(context, analysis: Path):
    """Validated S3 projection table and references (computed by figS3_cross_geometry.py)."""
    directory = analysis / S3_DIR
    provenance = json.loads((directory / "provenance.json").read_text())
    table = directory / "cross_geometry_coordinates.csv"
    if provenance.get("projection_table_sha256") != sha256(table):
        raise ValueError("S3 projection table differs from its provenance record")
    if provenance.get("context_cache_key") != context.cache_key:
        raise ValueError("S3 projections belong to a different analysis context")
    return pd.read_csv(table), provenance["references"], provenance


# Track D runs behind the subtree floor and the Figure 5 caption note (Track D6).
TRACK_D_RUNS = dict(compare="compare_20261015", control="apoptosis_subtrees_20261021",
                    single="apoptosis_single_20261022")
TRACK_D_TABLES = dict(compare=("subtree_metrics",), control=("decision",), single=("single",))
TRACK_D_LINKS = dict(compare_run="compare_20261015", apoptosis_run="apoptosis_20261019")
RESOLVED_NULL_CELLS = 12  # a resolved subtree's first-cousin null shuffles at least this many cells (D5a)
ZERO = 1e-9  # the D4b tie tolerance: d_SD at or below it is on the front


def _track_d(run: Path, analysis: Path):
    """Hash-checked Track D tables keyed (stage, name), and the input files keyed relative to ``run``."""
    from terminal_pareto import subtree_symmetry as ss
    tables, files = {}, []
    for stage, run_id in TRACK_D_RUNS.items():
        folder = ss.run_path(run_id) / "analysis"
        record = json.loads((folder / "provenance.json").read_text())
        files.append(folder / "provenance.json")
        for name in TRACK_D_TABLES[stage]:
            path = folder / f"{name}.csv"
            if sha256(path) != record["table_hashes"][name]:
                raise ValueError(f"{path} differs from its provenance record")
            tables[stage, name] = pd.read_csv(path)
            files.append(path)
        linked = {k: v for k, v in record.items() if k in TRACK_D_LINKS and v != TRACK_D_LINKS[k]}
        if linked:
            raise ValueError(f"{run_id} was built from other Track D runs: {linked}")
        if stage == "compare":
            if Path(ss.ROOT / record["canonical_analysis"]).resolve() != analysis.resolve():
                raise ValueError(f"{run_id} compares against another canonical analysis")
            changed = [n for n, digest in record["canonical_hashes"].items() if sha256(analysis / n) != digest]
            if changed:
                raise ValueError(f"Canonical files changed since {run_id}: {changed}")
    inputs = {os.path.relpath(path.resolve(), Path(run).resolve()): sha256(path) for path in files}
    return tables, inputs


def _resolution(tables, canonical: pd.DataFrame, floor: int = SUBTREE_FLOOR) -> dict:
    """Check that ``floor`` is the size above which every subtree is resolved (Track D6 rule).

    Resolved: the first-cousin null shuffles at least ``RESOLVED_NULL_CELLS``
    cells and the front spans at least one null SD on both axes. Every subtree
    at or above the floor must be resolved, and no subtree below it may be
    larger than the largest unresolved one (so the floor drops nothing the
    rule would keep).
    """
    from terminal_pareto import subtree_symmetry as ss
    metrics = tables["compare", "subtree_metrics"].query("travel == 'pooled'").set_index("subtree")
    if set(metrics.index) != set(canonical.subtree):
        raise ValueError("Track D subtree table and canonical table cover different subtrees")
    shown = canonical.set_index("subtree")
    for key, column in (("u_L", "u_lineage_lp"), ("d_LP", "d_lp")):
        delta = np.abs(metrics[key] - shown.loc[metrics.index, column]).max()
        if delta > 1e-9:
            raise ValueError(f"Track D {key} does not reproduce the canonical table ({delta:.2e})")
    resolved = ((metrics.null_support >= RESOLVED_NULL_CELLS)
                & (metrics[["travel_span_sd", "state_span_sd"]].min(axis=1) >= ss.TINY_SPAN_SD))
    largest = metrics[~resolved].n.idxmax()
    kept = metrics.n >= floor
    if not resolved[kept].all() or (metrics.n[~kept] > metrics.at[largest, "n"]).any():
        raise ValueError(f"A floor of {floor} cells no longer follows the resolution rule")
    return dict(floor=floor, kept=int(kept.sum()), resolved=int(resolved.sum()),
                largest_unresolved=(largest, int(metrics.at[largest, "n"])))


def _note_evidence(tables, shown: pd.DataFrame) -> dict:
    """What the Figure 5 caption note on ABplp cites (Track D4b-D4c)."""
    single = tables["single", "single"].query("subtree == 'ABplp' and variant == 'birth_endpoint'")
    alone = single[single.treatment == "only"]
    control = tables["control", "decision"].set_index(["subtree", "variant"]).loc[("ABplp", "birth_endpoint")]
    return dict(
        farthest=shown.loc[shown.d_lp.idxmax(), "subtree"],
        deaths=sorted(alone.death),
        alone_on_front={r.death: r.status == "ok" and r.d_SD <= ZERO for r in alone.itertuples()},
        control_draws=int(control.n_ok),
        control_on_front=int(round(control.d_SD_p * (control.n_ok + 1))) - 1)


def _global_display(result, context):
    """Endpoint coordinates for the cached global front and its null clouds."""
    from terminal_pareto.front_coordinates import build_endpoint_transform
    keys = {1: "first_cousin", 2: "second_cousin", 3: "third_cousin", "full": "full_random"}
    transform = build_endpoint_transform(
        result.raw_front_travel, result.raw_front_state,
        travel_optimum_index=context.spec.sweep_intervals, state_optimum_index=0,
        reference_analysis_id=result.analysis_cache_key,
        assignment_ids=[f"sweep:{index}" for index in range(context.spec.sweep_intervals + 1)])
    x, y = transform.transform(result.raw_front_travel, result.raw_front_state)
    natural_x, natural_y = transform.transform(np.asarray([result.natural_costs[0]]),
                                               np.asarray([result.natural_costs[1]]))
    return dict(x=x, y=y, nulls={key: transform.transform(*result.null_raw[raw]) for key, raw in keys.items()},
                natural=(float(natural_x[0]), float(natural_y[0])),
                xlabel="Travel distance\n(fraction of endpoint cost span)",
                ylabel="Cell-state distance\n(fraction of endpoint cost span)",
                metadata=transform.metadata(clipping=False))


def primary(run_id: str = PRIMARY_RUN, *, output_root: Path | None = None):
    """Primary pooled 275-edge terminal analysis (cosine molecular distance)."""
    from terminal_pareto import fig6a_figs2_ce_cell_types as cell_types
    from terminal_pareto.analysis_context import (DEFAULT_OUTPUT_ROOT, build_analysis_context,
                                                  validate_existing_context_manifest)
    from terminal_pareto.fig5_table1_ce_canonical_metrics import load_validated_canonical_metrics
    from terminal_pareto.front_coordinates import null_sd_coordinates
    from terminal_pareto.global_analysis import load_global_analysis
    from terminal_pareto.subtree_explore import _merge_small
    from publication.contracts import TerminalPrimaryInput

    context = build_analysis_context(PRIMARY_PROFILE, run_id=run_id, output_root=output_root or DEFAULT_OUTPUT_ROOT,
                                     sweep_intervals=PRIMARY_INTERVALS)
    analysis = context.run_paths.analysis
    if not (analysis / "analysis_manifest.json").is_file():
        raise FileNotFoundError(f"No pooled terminal analysis at {analysis}; run "
                                "`python terminal_pareto/analysis_pipeline.py` to create it")
    validate_existing_context_manifest(context)
    result = load_global_analysis(context)
    keys = {1: "first_cousin", 2: "second_cousin", 3: "third_cousin", "full": "full_random"}
    nulls = {k: null_sd_coordinates(*result.null_raw[name], natural_costs=result.natural_costs,
                                    null_stds=result.null_stds) for k, name in keys.items()}
    display = _global_display(result, context)
    canonical = load_validated_canonical_metrics(analysis / "ce_subtree_canonical_metrics.csv",
                                                 min_cells=PRIMARY_MIN_CELLS, iteration=PRIMARY_INTERVALS,
                                                 context=context)
    run = context.run_paths.root
    track_d, track_d_inputs = _track_d(run, analysis)
    resolution = _resolution(track_d, canonical)
    canonical = canonical[canonical.n >= SUBTREE_FLOOR].reset_index(drop=True)
    canonical_display = canonical[canonical["endpoint_ok"].astype(bool) & canonical["u_monotone"].astype(bool)].copy()
    summary, nodes, edges = _subtree_map(context, analysis, PRIMARY_MIN_CELLS, SUBTREE_FLOOR)
    retention, keypoints = cell_types.load_validated_cell_type_caches(analysis, PRIMARY_INTERVALS, context)
    within = _within_type(analysis, context)
    s3_frame, s3_references, s3_provenance = _s3(context, analysis)
    files = ["analysis_manifest.json", "global_terminal_analysis.json", "global_terminal_analysis.npz",
             "ce_subtree_canonical_metrics.csv", f"subtree_summary_min{PRIMARY_MIN_CELLS}.csv",
             "global_front_retention_by_cell_type.csv", "global_front_by_cell_type.csv",
             "within_type_fronts.csv", "within_type_summary.csv", "type_preserving_aggregate_front.csv",
             f"{S3_DIR}/cross_geometry_coordinates.csv", f"{S3_DIR}/provenance.json"]
    inputs = {f"analysis/{name}": sha256(analysis / name) for name in files}
    inputs.update(track_d_inputs)
    return TerminalPrimaryInput(
        run=run, analysis_id=context.cache_key, edges=len(context.terminal_nodes),
        subtrees=len(canonical), min_cells=SUBTREE_FLOOR,
        global_twr=result.twr, global_nulls=nulls, global_display=display,
        canonical=canonical, canonical_display=canonical_display, subtree_summary=summary,
        subtree_nodes=nodes, subtree_edges=edges,
        cell_type_retention=cell_types.merge_small_retention(retention), cell_type_keypoints=_merge_small(keypoints),
        within_type=within, s3_projections=s3_frame, s3_references=s3_references,
        notes=_note_evidence(track_d, canonical_display),
        validation=dict(context_cache_key=context.cache_key, global_analysis=result.analysis_cache_key,
                        s3_projection_count=s3_provenance["projection_count"],
                        canonical_rows=len(canonical), subtree_map_nodes=len(nodes), resolution=resolution),
        input_files=inputs)
