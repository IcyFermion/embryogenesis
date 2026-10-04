"""Regression and identity tests for pooled terminal analysis inputs."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto.analysis_context import build_analysis_context
from terminal_pareto.front_coordinates import (
    DegenerateEndpointSpan,
    EndpointTransform,
)
from terminal_pareto import pareto_engine as pe
from terminal_pareto.fig5_table1_ce_canonical_metrics import (
    CANONICAL_CACHE_MANIFEST,
    load_validated_canonical_metrics,
    write_canonical_cache_manifest,
)
from terminal_pareto.fig6a_figs2_ce_cell_types import (
    KEYPOINT_CACHE,
    FULL_RETENTION_CACHE,
    load_validated_cell_type_caches,
    write_cell_type_cache_manifest,
)
from terminal_pareto.global_analysis import (
    compute_global_analysis,
    load_global_analysis,
    save_global_analysis,
)
from terminal_pareto.subtree_analysis import (
    exact_cousin_stats,
    load_validated_subtree_summary,
    subtree_summary_manifest_path,
    write_subtree_summary_cache_manifest,
)


def test_profile_edge_and_subtree_regressions(tmp_path):
    legacy = build_analysis_context("embryo1_legacy", output_root=tmp_path)
    matched = build_analysis_context("embryo1_matched", output_root=tmp_path)
    pooled = build_analysis_context("pooled_tracking_v1", output_root=tmp_path)
    assert len(legacy.terminal_nodes) == 299
    assert len(legacy.subtrees(12)) == 43
    assert len(matched.terminal_nodes) == 275
    assert len(matched.subtrees(12)) == 42
    assert len(pooled.terminal_nodes) == 275
    assert len(pooled.subtrees(12)) == 42
    assert pooled.ordered_edges == matched.ordered_edges


def test_legacy_and_pooled_matrix_reconstruction(tmp_path):
    legacy = build_analysis_context("embryo1_legacy", output_root=tmp_path)
    rebuilt, expression, _ = pe.build_cost_matrices(
        legacy.terminal_nodes, legacy.terminal_parents,
        legacy.component_xyz_maps["embryo1"], legacy.protein_exp,
        legacy.prot_sel)
    np.testing.assert_array_equal(legacy.travel_matrix, rebuilt)
    np.testing.assert_array_equal(legacy.expression_matrix, expression)

    pooled = build_analysis_context("pooled_tracking_v1", output_root=tmp_path)
    expected = np.mean([
        matrix / pooled.component_natural_totals[label]
        for label, matrix in pooled.component_travel_matrices.items()
    ], axis=0)
    np.testing.assert_allclose(pooled.travel_matrix, expected, atol=0, rtol=0)
    np.testing.assert_allclose(np.trace(pooled.travel_matrix), 1.0,
                               atol=1e-12, rtol=0)


def test_fixed_global_normalization_for_subtrees(tmp_path):
    pooled = build_analysis_context("pooled_tracking_v1", output_root=tmp_path)
    name, terms = pooled.subtrees(12)[1]
    travel, expression = pooled.cost_matrices(terms)
    indices = pooled.edge_indices(terms)
    expected = np.mean([
        matrix[np.ix_(indices, indices)]
        / pooled.component_natural_totals[label]
        for label, matrix in pooled.component_travel_matrices.items()
    ], axis=0)
    np.testing.assert_allclose(travel, expected, atol=0, rtol=0)
    np.testing.assert_array_equal(
        expression, pooled.expression_matrix[np.ix_(indices, indices)])
    assert not np.isclose(np.trace(travel), 1.0)


def test_pooling_identities_synthetic():
    rng = np.random.default_rng(187)
    matrix = rng.uniform(0.1, 2.0, size=(7, 7))
    total = float(np.trace(matrix))
    identical = np.mean([matrix / total] * 3, axis=0)
    np.testing.assert_allclose(identical, matrix / total)

    matrices = [rng.uniform(0.1, 2.0, size=(7, 7)) for _ in range(3)]
    totals = [float(np.trace(value)) for value in matrices]
    reference = np.mean([value / total for value, total
                         in zip(matrices, totals)], axis=0)
    order = [2, 0, 1]
    permuted = np.mean([matrices[i] / totals[i] for i in order], axis=0)
    np.testing.assert_allclose(reference, permuted)
    scaled = matrices.copy()
    scaled[1] = scaled[1] * 17.0
    scaled_totals = totals.copy()
    scaled_totals[1] *= 17.0
    invariant = np.mean([value / total for value, total
                         in zip(scaled, scaled_totals)], axis=0)
    np.testing.assert_allclose(reference, invariant)


def test_endpoint_transform_roundtrip_and_outside_points():
    transform = EndpointTransform.from_endpoints(
        reference_analysis_id="fixture",
        travel_optimum_assignment_id="A",
        state_optimum_assignment_id="B",
        travel_optimum_costs=(2.0, 9.0),
        state_optimum_costs=(8.0, 3.0),
    )
    travel = np.array([2.0, 8.0, 5.0, 10.0])
    state = np.array([9.0, 3.0, 6.0, 1.0])
    x, y = transform.transform(travel, state)
    np.testing.assert_allclose([x[0], y[0]], [0.0, 1.0])
    np.testing.assert_allclose([x[1], y[1]], [1.0, 0.0])
    assert x[-1] > 1.0 and y[-1] < 0.0
    inverse = transform.inverse(x, y)
    np.testing.assert_allclose(inverse[0], travel)
    np.testing.assert_allclose(inverse[1], state)


def test_endpoint_transform_rejects_degenerate_span():
    try:
        EndpointTransform.from_endpoints(
            reference_analysis_id="fixture",
            travel_optimum_assignment_id="A",
            state_optimum_assignment_id="B",
            travel_optimum_costs=(2.0, 4.0),
            state_optimum_costs=(2.0, 3.0),
        )
    except DegenerateEndpointSpan:
        pass
    else:
        raise AssertionError("Degenerate travel span was accepted")


def test_exact_pooled_null_matches_shared_permutation_enumeration():
    rng = np.random.default_rng(912)
    components = rng.uniform(0.1, 3.0, size=(3, 4, 4))
    totals = np.trace(components, axis1=1, axis2=2)
    pooled = np.mean(components / totals[:, None, None], axis=0)
    state = rng.uniform(0.1, 2.0, size=(4, 4))
    groups = [[0, 1], [2, 3]]
    permutations = [
        np.array([a, b, c, d])
        for a, b in ((0, 1), (1, 0))
        for c, d in ((2, 3), (3, 2))
    ]
    rows = np.arange(4)
    travel = np.asarray([pooled[rows, perm].sum() for perm in permutations])
    expression = np.asarray([state[rows, perm].sum()
                             for perm in permutations])
    mx, me, vx, ve, cov = exact_cousin_stats(pooled, state, groups)
    np.testing.assert_allclose(
        [mx, me, vx, ve, cov],
        [travel.mean(), expression.mean(), travel.var(), expression.var(),
         np.mean((travel - travel.mean())
                 * (expression - expression.mean()))],
        atol=1e-14, rtol=1e-14,
    )
    component_ratios = np.asarray([
        [matrix[rows, perm].sum() / total
         for matrix, total in zip(components, totals)]
        for perm in permutations
    ])
    np.testing.assert_allclose(travel, component_ratios.mean(axis=1))


def test_global_cache_rejects_cross_profile_reuse(tmp_path):
    pooled = build_analysis_context(
        "pooled_tracking_v1", output_root=tmp_path,
        run_id="pooled_fixture", sweep_intervals=4, null_draws=8)
    result = compute_global_analysis(pooled, n_lineage_null=3)
    # The production loader fixes n_lineage_null=100 in its cache identity;
    # write a deliberately non-production cache and ensure it is rejected.
    save_global_analysis(result, pooled)
    try:
        load_global_analysis(pooled)
    except ValueError:
        pass
    else:
        raise AssertionError("A cache with different analysis settings was accepted")

    matched = build_analysis_context(
        "embryo1_matched", output_root=tmp_path,
        run_id="pooled_fixture", sweep_intervals=4, null_draws=8)
    # Point the matched profile at its own empty profile-specific directory.
    try:
        load_global_analysis(matched)
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("A pooled cache leaked into a matched-profile run")


def test_cell_type_cache_mismatch_rejected_before_loading(tmp_path):
    analysis_out = tmp_path / "analysis"
    analysis_out.mkdir()
    full = analysis_out / FULL_RETENTION_CACHE.name
    keypoint = analysis_out / KEYPOINT_CACHE.name
    # Deliberately invalid CSV proves a configuration mismatch is rejected
    # before pandas attempts to load either table.
    full.write_bytes(b"not a csv cache")
    keypoint.write_bytes(b"also not a csv cache")
    context_300 = SimpleNamespace(
        profile="pooled_tracking_v1", cache_key="context-300",
        terminal_nodes=list(range(275)))
    write_cell_type_cache_manifest(analysis_out, 300, context_300)
    before = {path.name: path.read_bytes()
              for path in (full, keypoint,
                           analysis_out / "fig6a_cell_type_cache_manifest.json")}
    context_4 = SimpleNamespace(
        profile="pooled_tracking_v1", cache_key="context-4",
        terminal_nodes=list(range(275)))
    try:
        load_validated_cell_type_caches(analysis_out, 4, context_4)
    except ValueError as error:
        assert "configuration mismatch" in str(error)
    else:
        raise AssertionError("Mismatched cell-type cache was accepted")
    after = {path.name: path.read_bytes()
             for path in (full, keypoint,
                          analysis_out / "fig6a_cell_type_cache_manifest.json")}
    assert after == before


def test_figure4_context_mismatch_rejected_before_loading(tmp_path):
    context = build_analysis_context(
        "pooled_tracking_v1", output_root=tmp_path,
        run_id="figure4-cache-fixture")
    context.write()
    analysis_out = context.run_paths.analysis
    summary = analysis_out / "subtree_summary_min12.csv"
    # Invalid bytes prove rejection occurs before pandas reads the table.
    summary.write_bytes(b"not a subtree summary")
    write_subtree_summary_cache_manifest(analysis_out, 12, context)
    context_manifest = analysis_out / "analysis_manifest.json"
    payload = json.loads(context_manifest.read_text())
    payload["cache_key"] = "deliberately-mismatched-context"
    context_manifest.write_text(json.dumps(payload, indent=2) + "\n")
    cache_manifest = subtree_summary_manifest_path(analysis_out, 12)
    before = {path.name: path.read_bytes()
              for path in (summary, cache_manifest, context_manifest)}
    try:
        load_validated_subtree_summary(analysis_out, 12, context)
    except ValueError as error:
        assert "context manifest mismatch" in str(error)
    else:
        raise AssertionError("Figure 4 accepted a mismatched context manifest")
    after = {path.name: path.read_bytes()
             for path in (summary, cache_manifest, context_manifest)}
    assert after == before


def test_figure5_cache_identity_rejected_before_loading(tmp_path):
    context = build_analysis_context(
        "pooled_tracking_v1", output_root=tmp_path,
        run_id="figure5-cache-fixture")
    context.write()
    analysis_out = context.run_paths.analysis
    metrics = analysis_out / "ce_subtree_canonical_metrics.csv"
    curves = analysis_out / "ce_subtree_canonical_curves.npz"
    summary = analysis_out / "subtree_summary_min12.csv"
    metrics.write_bytes(b"not canonical metrics")
    curves.write_bytes(b"not canonical curves")
    summary.write_bytes(b"not a subtree summary")
    write_canonical_cache_manifest(analysis_out, 12, 300, context)
    manifest = analysis_out / CANONICAL_CACHE_MANIFEST
    payload = json.loads(manifest.read_text())
    payload["iteration"] = 4
    manifest.write_text(json.dumps(payload, indent=2) + "\n")
    before = {path.name: path.read_bytes()
              for path in (metrics, curves, summary, manifest,
                           analysis_out / "analysis_manifest.json")}
    try:
        load_validated_canonical_metrics(
            metrics, min_cells=12, iteration=300, context=context)
    except ValueError as error:
        assert "configuration mismatch" in str(error)
    else:
        raise AssertionError("Figure 5 accepted mismatched cache identity")
    after = {path.name: path.read_bytes()
             for path in (metrics, curves, summary, manifest,
                          analysis_out / "analysis_manifest.json")}
    assert after == before


def main() -> None:
    """Run the focused migration suite without requiring pytest."""
    with tempfile.TemporaryDirectory(prefix="pooled-migration-tests-") as tmp:
        tmp_path = Path(tmp)
        test_profile_edge_and_subtree_regressions(tmp_path)
        test_legacy_and_pooled_matrix_reconstruction(tmp_path)
        test_fixed_global_normalization_for_subtrees(tmp_path)
        test_global_cache_rejects_cross_profile_reuse(tmp_path)
        test_cell_type_cache_mismatch_rejected_before_loading(tmp_path)
        test_figure4_context_mismatch_rejected_before_loading(tmp_path)
        test_figure5_cache_identity_rejected_before_loading(tmp_path)
    test_pooling_identities_synthetic()
    test_endpoint_transform_roundtrip_and_outside_points()
    test_endpoint_transform_rejects_degenerate_span()
    test_exact_pooled_null_matches_shared_permutation_enumeration()
    print("11 pooled-migration tests passed")


if __name__ == "__main__":
    main()
