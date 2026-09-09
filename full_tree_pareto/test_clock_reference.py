"""Scientific regression checks using the measured publication configuration.

Run: conda run -n dev python -m unittest full_tree_pareto.test_clock_reference
"""

import copy
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from full_tree_pareto import publication_analysis as analysis


class ClockReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = analysis.load_ce_protein_context()
        cls.reference = analysis.separate_clock_reference(
            cls.context, n_samples=300, seed=242, capture_edge_norms=True
        )

    def test_scope_fit_and_recomposed_costs(self):
        c, ref = self.context, self.reference
        mapping = np.asarray(c.tree.lineage_id_mapping)
        ids = ref.evaluated_tree_ids
        parents = np.asarray(c.tree.parent_list)[ids]
        delta = c.optimization.exp_mat[mapping[ids]] - c.optimization.exp_mat[mapping[parents]]
        self.assertEqual(len(set(ids)), 1000)
        self.assertEqual(len(ref.root_tree_ids), 4)
        np.testing.assert_allclose(ref.covariance[3:, 3:], delta.T @ delta / 1000)
        np.testing.assert_array_equal(ref.covariance[:3, 3:], np.zeros((3, 20)))
        self.assertAlmostEqual(np.linalg.norm(delta, axis=1).sum(), c.optimization.lineage_exp_cost)
        np.testing.assert_allclose(ref.edge_travel.sum(axis=1), ref.travel_cost)
        np.testing.assert_allclose(ref.edge_cell_state.sum(axis=1), ref.cell_state_cost)

    def test_duration_perturbation_cannot_change_protein_fit(self):
        altered = copy.copy(self.context)
        altered.tree = copy.copy(self.context.tree)
        altered.tree.branch_time_length = np.asarray(self.context.tree.branch_time_length) * 4
        other = analysis.separate_clock_reference(altered, n_samples=3, seed=242)
        np.testing.assert_allclose(other.covariance[3:, 3:], self.reference.covariance[3:, 3:])
        np.testing.assert_allclose(other.covariance[:3, :3] * 4, self.reference.covariance[:3, :3])

    def test_gaussian_second_moments_for_both_clocks(self):
        ref = self.reference
        durations = np.asarray(self.context.tree.branch_time_length)[ref.evaluated_tree_ids]
        for norms, cov, clock in [
            (ref.edge_travel, ref.covariance[:3, :3], durations),
            (ref.edge_cell_state, ref.covariance[3:, 3:], np.ones(1000)),
        ]:
            # Independent Gaussian quadratic forms have mean tr(C) and
            # variance 2 tr(C^2). This checks the sampler against that identity.
            squared_totals = (norms**2).sum(axis=1)
            expected = clock.sum() * np.trace(cov)
            se = np.sqrt(2 * (clock**2).sum() * np.trace(cov @ cov) / len(norms))
            self.assertLess(abs(squared_totals.mean() - expected), 6 * se)

    def test_determinism_and_legacy_regression(self):
        again = analysis.separate_clock_reference(self.context, n_samples=300, seed=242)
        np.testing.assert_array_equal(again.travel_cost, self.reference.travel_cost)
        np.testing.assert_array_equal(again.cell_state_cost, self.reference.cell_state_cost)
        old = analysis.parametric_brownian_reference(self.context, n_samples=300, seed=242)
        # Values independently recorded before the shared sampler refactor.
        self.assertAlmostEqual(old.travel_cost.mean(), 4827.897001883872, places=8)
        self.assertAlmostEqual(old.cell_state_cost.mean(), 5240.695830423984, places=8)

    def test_invalid_draw_count(self):
        with self.assertRaises(ValueError):
            analysis.separate_clock_reference(self.context, n_samples=0)

    def test_publication_selects_new_reference_and_rejects_old_provenance(self):
        fronts = pd.read_csv(analysis.CACHE_ROOT / "ce_full_tree_layerwise_fronts.csv")
        spanning = pd.read_csv(analysis.CACHE_ROOT / "ce_full_tree_spanning_forest.csv")
        with (
            patch.object(analysis, "N_PARAMETRIC_REFERENCE", 300),
            patch.object(analysis, "separate_clock_reference", return_value=self.reference) as sampler,
            patch.object(analysis, "parametric_brownian_reference", side_effect=AssertionError("Historical sampler selected")),
        ):
            heuristics, nulls, summary, ref = analysis.collective_analysis(self.context, fronts, spanning)
            sampler.assert_called_once_with(self.context, n_samples=300, seed=242)
            self.assertIs(ref, self.reference)
            selected = nulls["null_model"].eq(analysis.PUBLICATION_REFERENCE_LABEL)
            self.assertEqual(selected.sum(), 300)
            self.assertNotIn("Parametric Brownian reference", set(nulls["null_model"]))
            np.testing.assert_array_equal(nulls.loc[selected, "cell_state_cost"], ref.cell_state_cost)
            stale = nulls.copy()
            stale.loc[selected, "source"] = "all-edge fixed-topology parametric bootstrap"
            with self.assertRaisesRegex(AssertionError, "separate-clock"):
                analysis.validate_collective(heuristics, stale, summary)


if __name__ == "__main__":
    unittest.main()
