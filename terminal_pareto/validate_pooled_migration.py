"""Scientific and cache validation for the pooled-travel migration candidate.

Numerical only. Figure organization, wrapper wording and page checks moved to
``publication/tests`` when figures moved to the publication front end.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto import pareto_engine as pe
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
)
from terminal_pareto.front_coordinates import build_endpoint_transform
from terminal_pareto.global_analysis import (
    compute_global_analysis,
    load_global_analysis,
)


MIGRATION_VERSION = "pooled-travel-validation-v1"
S3_SOURCES = ("embryo1", "embryo2", "embryo3", "pooled")


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _check(checks: list[dict], name: str, condition: bool, detail: str = ""):
    checks.append({"name": name, "passed": bool(condition), "detail": detail})
    print(f"[{'PASS' if condition else 'FAIL'}] {name} {detail}")


def _checkpoint_mismatches(checkpoint: Path, manifest_name: str) -> list[str]:
    """Compare a pre-organization hash snapshot with its original files."""
    record = json.loads((checkpoint / manifest_name).read_text())
    root = Path(record["root"])
    # Promotion freezes the previous publication at this location. Historical
    # checkpoints continue to verify those exact bytes, not the new release.
    legacy_publication = ROOT / "terminal_pareto/output/legacy/embryo1/publication"
    if manifest_name == "accepted_publication_sha256.json" and legacy_publication.exists():
        root = legacy_publication
    return [
        relative for relative, expected in record["files"].items()
        if not (root / relative).exists()
        or _hash(root / relative) != expected
    ]


def _validate_s3_projection(pooled, result, projections: Path) -> tuple[bool, str]:
    """Independently replay every source/target/sweep projection."""
    table_path = projections / "cross_geometry_coordinates.csv"
    provenance_path = projections / "provenance.json"
    if not table_path.exists() or not provenance_path.exists():
        return False, "missing S3 projection table or provenance"
    frame = pd.read_csv(table_path)
    provenance = json.loads(provenance_path.read_text())
    if len(frame) != 4 * 4 * (pooled.spec.sweep_intervals + 1):
        return False, f"unexpected projection count={len(frame)}"
    if (provenance.get("projection_count") != len(frame)
            or provenance.get("projection_table_sha256") != _hash(table_path)
            or provenance.get("context_cache_key") != pooled.cache_key):
        return False, "projection provenance identity mismatch"
    references = provenance.get("references", {})
    display = provenance.get("display", {})
    offsets = display.get("inset_offsets_from_natural", {})
    if offsets != {"travel": [-0.026, 0.036], "state": [-0.09, 0.09]}:
        return False, "inset offsets differ from declared layout"
    input_hashes = provenance.get("input_sha256", {})
    if not input_hashes or any(
        not Path(path).exists() or _hash(Path(path)) != value
        for path, value in input_hashes.items()
    ):
        return False, "S3 source hashes are missing or have changed"

    historical = (ROOT / "terminal_pareto" / "output"
                  / "tracking_geometry_sensitivity" / "assignments_and_costs.npz")
    with np.load(historical, allow_pickle=False) as saved:
        assignments = {
            source: np.argsort(saved[f"{source}_assignments"], axis=1)
            for source in S3_SOURCES[:3]
        }
        matrices = {
            source: saved[f"{source}_travel_matrix"].copy()
            for source in S3_SOURCES[:3]
        }
        state_matrix = saved["expression_matrix"].copy()
        alpha = saved["alpha"].copy()
        if not np.array_equal(saved["children"], pooled.terminal_nodes):
            return False, "historical source cohort differs"
    assignments["pooled"] = result.assignments
    matrices["pooled"] = pooled.travel_matrix
    rows = np.arange(len(pooled.terminal_nodes))
    errors = []
    for target in S3_SOURCES:
        reference = references.get(target)
        if not reference or reference.get("display_mode") != "endpoint":
            return False, f"missing endpoint reference for {target}"
        native = assignments[target]
        native_travel = matrices[target][native, rows].sum(axis=1)
        native_state = state_matrix[native, rows].sum(axis=1)
        expected_anchor = (
            native_travel[-1], native_state[-1],
            native_travel[0], native_state[0],
        )
        actual_anchor = (
            reference["travel_at_travel_optimum"],
            reference["state_at_travel_optimum"],
            reference["travel_at_state_optimum"],
            reference["state_at_state_optimum"],
        )
        if not np.allclose(actual_anchor, expected_anchor, atol=1e-10, rtol=0):
            return False, f"endpoint anchors differ for {target}"
        natural = reference["natural_display"]
        expected_limits = {
            "travel": [natural[0] + offset for offset in offsets["travel"]],
            "state": [natural[1] + offset for offset in offsets["state"]],
        }
        if not np.allclose(
            display["inset_limits_by_target"][target]["travel"],
            expected_limits["travel"], atol=1e-12, rtol=0,
        ) or not np.allclose(
            display["inset_limits_by_target"][target]["state"],
            expected_limits["state"], atol=1e-12, rtol=0,
        ):
            return False, f"inset limits differ for {target}"
        for source in S3_SOURCES:
            block = frame[(frame.target == target) & (frame.source == source)]
            block = block.sort_values("step")
            if len(block) != len(alpha) or not np.array_equal(
                block.step.to_numpy(), np.arange(len(alpha))
            ) or not np.allclose(block.alpha, alpha, atol=1e-15, rtol=0):
                return False, f"sweep order differs for {target}/{source}"
            assignment = assignments[source]
            travel = matrices[target][assignment, rows].sum(axis=1)
            state = state_matrix[assignment, rows].sum(axis=1)
            x = ((travel - reference["travel_at_travel_optimum"])
                 / reference["travel_span"])
            y = ((state - reference["state_at_state_optimum"])
                 / reference["state_span"])
            errors.append(float(np.max(np.abs(
                block[["travel", "state", "x", "y"]].to_numpy()
                - np.column_stack((travel, state, x, y))))))
            if not (block.in_place.to_numpy(dtype=bool) == (target == source)).all():
                return False, f"native/transfer flag differs for {target}/{source}"
            expected_hash = hashlib.sha256(assignment.tobytes()).hexdigest()
            if provenance.get("source_assignments_sha256", {}).get(source) != expected_hash:
                return False, f"assignment hash differs for {source}"
    maximum = max(errors)
    return maximum < 1e-10, f"projections={len(frame)} max error={maximum:.3g}"


def validate(run_id: str, output_root: Path) -> dict:
    checks: list[dict] = []
    pooled = build_analysis_context(
        "pooled_tracking_v1", run_id=run_id, output_root=output_root)
    matched = build_analysis_context(
        "embryo1_matched", run_id="validation", output_root=output_root)
    legacy = build_analysis_context(
        "embryo1_legacy", run_id="validation", output_root=output_root)

    # Legacy reconstruction uses the unchanged single-embryo matrix builder.
    rebuilt, expression, _ = pe.build_cost_matrices(
        legacy.terminal_nodes, legacy.terminal_parents,
        legacy.component_xyz_maps["embryo1"], legacy.protein_exp,
        legacy.prot_sel)
    _check(checks, "legacy travel matrix exact",
           np.array_equal(rebuilt, legacy.travel_matrix))
    _check(checks, "legacy expression matrix exact",
           np.array_equal(expression, legacy.expression_matrix))
    from terminal_pareto.fig2_fig3_ce_terminal_pareto import _load_analysis
    legacy_twr, legacy_nulls = _load_analysis()
    legacy_profile_result = compute_global_analysis(legacy)
    twr_keys = ("xyz_arr", "exp_arr", "traditional_er",
                "lineage_mean_dist", "per_cell_dist")
    front_errors = [
        float(np.max(np.abs(np.asarray(legacy_profile_result.twr[key])
                            - np.asarray(legacy_twr[key]))))
        for key in twr_keys
    ]
    _check(checks, "legacy global front and structure reproduce",
           max(front_errors) == 0.0
           and legacy_profile_result.twr["kp"] == legacy_twr["kp"],
           f"max absolute error={max(front_errors):.3g}")
    null_errors = []
    for key, raw_key in ((1, "first_cousin"), (2, "second_cousin"),
                         (3, "third_cousin"), ("full", "full_random")):
        actual_x = ((legacy_profile_result.null_raw[raw_key][0]
                     - legacy_profile_result.natural_costs[0])
                    / legacy_profile_result.null_stds[0])
        actual_y = ((legacy_profile_result.null_raw[raw_key][1]
                     - legacy_profile_result.natural_costs[1])
                    / legacy_profile_result.null_stds[1])
        null_errors.extend([
            float(np.max(np.abs(actual_x - legacy_nulls[key][0]))),
            float(np.max(np.abs(actual_y - legacy_nulls[key][1]))),
        ])
    _check(checks, "legacy null draws reproduce to roundoff",
           max(null_errors) < 1e-12,
           f"max absolute error={max(null_errors):.3g}")

    fixture_path = (ROOT / "terminal_pareto" / "output"
                    / "tracking_metric_comparison"
                    / "candidate_travel_matrices.npz")
    with np.load(fixture_path, allow_pickle=False) as fixture:
        _check(checks, "pooled ordered children match audited fixture",
               np.array_equal(np.asarray(pooled.terminal_nodes),
                              fixture["children"]))
        _check(checks, "pooled ordered parents match audited fixture",
               np.array_equal(np.asarray(pooled.terminal_parents),
                              fixture["natural_parents"]))
        _check(checks, "pooled matrix exact reconstruction",
               np.array_equal(pooled.travel_matrix,
                              fixture["normalized_mean"]))
        _check(checks, "component natural totals match audited fixture",
               np.allclose(
                   [pooled.component_natural_totals[label]
                    for label in ("embryo1", "embryo2", "embryo3")],
                   fixture["natural_travel"], atol=1e-12, rtol=0))

    _check(checks, "profile edge counts",
           (len(legacy.terminal_nodes), len(matched.terminal_nodes),
            len(pooled.terminal_nodes)) == (299, 275, 275),
           "legacy/matched/pooled = 299/275/275")
    _check(checks, "profile subtree counts",
           (len(legacy.subtrees(12)), len(matched.subtrees(12)),
            len(pooled.subtrees(12))) == (43, 42, 42),
           "legacy/matched/pooled = 43/42/42")
    _check(checks, "pooled natural total is one",
           np.isclose(np.trace(pooled.travel_matrix), 1.0,
                      atol=1e-12, rtol=0))

    # Every pooled assignment total must equal the component-average ratio.
    rng = np.random.default_rng(428)
    rows = np.arange(len(pooled.terminal_nodes))
    identity = np.arange(len(rows))
    assignments = [identity] + [rng.permutation(len(rows)) for _ in range(8)]
    pooling_errors = []
    for assignment in assignments:
        pooled_total = pooled.travel_matrix[rows, assignment].sum()
        component_mean = np.mean([
            matrix[rows, assignment].sum()
            / pooled.component_natural_totals[label]
            for label, matrix in pooled.component_travel_matrices.items()
        ])
        pooling_errors.append(abs(pooled_total - component_mean))
    _check(checks, "pooled assignment identity",
           max(pooling_errors) < 1e-12,
           f"max absolute error={max(pooling_errors):.3g}")

    result = load_global_analysis(pooled)
    transform = build_endpoint_transform(
        result.raw_front_travel, result.raw_front_state,
        travel_optimum_index=pooled.spec.sweep_intervals,
        state_optimum_index=0,
        reference_analysis_id=result.analysis_cache_key,
    )
    x, y = transform.transform(
        result.raw_front_travel, result.raw_front_state)
    inverse_travel, inverse_state = transform.inverse(x, y)
    _check(checks, "endpoint anchors",
           np.allclose([x[-1], y[-1], x[0], y[0]], [0, 1, 1, 0],
                       atol=1e-12, rtol=0))
    _check(checks, "endpoint inverse",
           np.allclose(inverse_travel, result.raw_front_travel,
                       atol=1e-13, rtol=1e-13)
           and np.allclose(inverse_state, result.raw_front_state,
                           atol=1e-13, rtol=1e-13))
    raw_dominance = ((result.raw_front_travel[:, None]
                      <= result.raw_front_travel[None, :])
                     & (result.raw_front_state[:, None]
                        <= result.raw_front_state[None, :]))
    display_dominance = ((x[:, None] <= x[None, :])
                         & (y[:, None] <= y[None, :]))
    _check(checks, "display preserves dominance",
           np.array_equal(raw_dominance, display_dominance))

    projections = pooled.run_paths.analysis / "s3_cross_geometry"
    s3_ok, s3_detail = _validate_s3_projection(pooled, result, projections)
    _check(checks, "S3 saved assignments and 4,816 projections replay",
           s3_ok, s3_detail)
    subtree_manifest_path = (pooled.run_paths.analysis
                             / "subtree_summary_min12_manifest.json")
    subtree_manifest_ok = False
    if subtree_manifest_path.exists():
        subtree_manifest = json.loads(subtree_manifest_path.read_text())
        subtree_manifest_ok = (
            subtree_manifest.get("profile") == pooled.profile
            and subtree_manifest.get("context_cache_key") == pooled.cache_key
            and subtree_manifest.get("min_cells") == 12
            and subtree_manifest.get("terminal_edges") == 275
        )
    _check(checks, "Figure 4 subtree cache identity matches pooled run",
           subtree_manifest_ok)

    canonical_manifest_path = (pooled.run_paths.analysis
                               / "fig5_canonical_cache_manifest.json")
    canonical_manifest_ok = False
    if canonical_manifest_path.exists():
        canonical_manifest = json.loads(canonical_manifest_path.read_text())
        canonical_manifest_ok = (
            canonical_manifest.get("profile") == pooled.profile
            and canonical_manifest.get("context_cache_key") == pooled.cache_key
            and canonical_manifest.get("min_cells") == 12
            and canonical_manifest.get("iteration") == 300
            and canonical_manifest.get("terminal_edges") == 275
        )
    _check(checks, "Figure 5 canonical cache identity matches pooled run",
           canonical_manifest_ok)

    fig6a_manifest_path = (pooled.run_paths.analysis
                           / "fig6a_cell_type_cache_manifest.json")
    fig6a_manifest_ok = False
    if fig6a_manifest_path.exists():
        fig6a_manifest = json.loads(fig6a_manifest_path.read_text())
        fig6a_manifest_ok = (
            fig6a_manifest.get("profile") == pooled.profile
            and fig6a_manifest.get("context_cache_key") == pooled.cache_key
            and fig6a_manifest.get("iteration") == 300
            and fig6a_manifest.get("terminal_edges") == 275
        )
    _check(checks, "Figure 6A cache identity manifest matches pooled run",
           fig6a_manifest_ok)
    legacy_metrics = pd.read_csv(
        ROOT / "terminal_pareto" / "output" / "runs"
        / "baseline_legacy_20260920" / "archive"
        / "ce_subtree_canonical_metrics.csv")
    comparison = pd.read_csv(
        ROOT / "terminal_pareto" / "output"
        / "tracking_metric_comparison" / "matched_subtree_metrics.csv")
    matched_metrics = comparison[comparison["metric"] == "embryo1_matched"]
    pooled_metrics = pd.read_csv(
        pooled.run_paths.analysis / "ce_subtree_canonical_metrics.csv")
    profile_rows = []
    for label, frame, edge_count in (
        ("legacy embryo 1", legacy_metrics, 299),
        ("matched embryo 1", matched_metrics, 275),
        ("pooled", pooled_metrics, 275),
    ):
        on_column = ("exactly_on_front" if "exactly_on_front" in frame
                     else "natural_on_grid")
        profile_rows.append({
            "profile": label,
            "edges": edge_count,
            "eligible_subtrees": len(frame),
            "natural_on_sampled_front": int(frame[on_column].astype(bool).sum()),
            "P0_d_LP": float(frame.loc[frame["subtree"] == "P0", "d_lp"].iloc[0]),
            "P0_max_edge_retention": float(
                frame.loc[frame["subtree"] == "P0", "max_er"].iloc[0]),
        })
    profile_table = pd.DataFrame(profile_rows)

    baseline_manifest = (ROOT / "terminal_pareto" / "output" / "runs"
                         / "baseline_legacy_20260920"
                         / "archive_sha256.json")
    _check(checks, "baseline archive manifest present",
           baseline_manifest.exists())
    baseline_mismatches: list[str] = []
    if baseline_manifest.exists():
        archived_hashes = json.loads(baseline_manifest.read_text())["files"]
        # The immutable baseline is the legacy regression target after promotion.
        accepted_root = baseline_manifest.parent / "archive"
        for relative, expected_hash in archived_hashes.items():
            current = accepted_root / relative
            if not current.exists() or _hash(current) != expected_hash:
                baseline_mismatches.append(relative)
    _check(checks, "archived legacy baseline unchanged",
           not baseline_mismatches,
           f"mismatches={baseline_mismatches[:5]}")

    checkpoint = (pooled.run_paths.validation
                  / "organization_checkpoint_20260922")
    for label, manifest_name in (
        ("analysis caches", "analysis_sha256.json"),
        ("accepted publication assets", "accepted_publication_sha256.json"),
        ("tracking sensitivity inputs", "tracking_sensitivity_sha256.json"),
    ):
        mismatches = (_checkpoint_mismatches(checkpoint, manifest_name)
                      if (checkpoint / manifest_name).exists()
                      else ["missing pre-organization checkpoint"])
        _check(checks, f"pre-organization {label} unchanged",
               not mismatches, f"mismatches={mismatches[:5]}")

    report = {
        "version": MIGRATION_VERSION,
        "profile": pooled.profile,
        "run_id": run_id,
        "context_cache_key": pooled.cache_key,
        "analysis_cache_key": result.analysis_cache_key,
        "checks": checks,
        "profiles": profile_rows,
        "artifact_hashes": {
            "analysis_manifest": _hash(
                pooled.run_paths.analysis / "analysis_manifest.json"),
            "global_cache": _hash(
                pooled.run_paths.analysis / "global_terminal_analysis.npz"),
            "s3_projection_table": _hash(projections / "cross_geometry_coordinates.csv"),
            "baseline_manifest": _hash(baseline_manifest),
        },
    }
    validation = pooled.run_paths.validation
    validation.mkdir(parents=True, exist_ok=True)
    (validation / "migration_validation.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = list(profile_table.columns)
    table_lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in profile_table.itertuples(index=False, name=None):
        rendered = [f"{value:.4f}" if isinstance(value, float) else str(value)
                    for value in row]
        table_lines.append("| " + " | ".join(rendered) + " |")
    table = "\n".join(table_lines)
    lines = [
        "# Pooled-travel migration validation",
        "",
        "The three profiles separate the legacy cohort, the coverage change, "
        "and the pooled travel aggregation. Canonical subtree metrics use "
        "exact first-cousin moments for the matched and pooled profiles.",
        "",
        table,
        "",
        "## Checks",
        "",
    ]
    lines.extend(
        f"- [{'x' if row['passed'] else ' '}] {row['name']}"
        + (f" — {row['detail']}" if row["detail"] else "")
        for row in checks
    )
    (validation / "migration_validation.md").write_text("\n".join(lines) + "\n")
    if not all(row["passed"] for row in checks):
        raise AssertionError("Pooled migration validation failed")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="migration_candidate_20260920")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args(argv)
    validate(args.run_id, args.output_root)


if __name__ == "__main__":
    main()
