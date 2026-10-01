"""Focused scientific checks for pooling and endpoint-tied RNA comparisons."""

import itertools
import unittest

import numpy as np

from terminal_pareto.cross_species_analysis import analyze, endpoint_assignment, pool_matrices
from terminal_pareto.subtree_analysis import exact_cousin_stats


class CrossSpeciesTests(unittest.TestCase):
    def test_pooling_averages_distances_and_uses_fixed_reference_totals(self):
        a = np.array([[3., 9.], [6., 7.]])
        b = np.array([[8., 4.], [12., 12.]])
        pooled = pool_matrices({"a": a, "b": b}, {"a": 10., "b": 20.})
        self.assertAlmostEqual(np.trace(pooled), 1.)
        # Selecting a cohort must preserve global denominators, including a
        # natural total below one. Renormalizing the subset changes weights.
        restricted = pool_matrices({"a": a[:1, :1], "b": b[:1, :1]}, {"a": 10., "b": 20.})
        self.assertAlmostEqual(restricted[0, 0], .35)
        np.testing.assert_array_equal(restricted, pooled[:1, :1])

    def test_pooling_rejects_missing_or_nonpositive_denominator(self):
        for totals in ({}, {"a": 0}, {"a": np.nan}):
            with self.assertRaises(ValueError):
                pool_matrices({"a": np.ones((2, 2))}, totals)

    def test_endpoint_tie_selects_secondary_optimum(self):
        primary = np.zeros((2, 2))
        secondary = np.array([[5., 1.], [1., 5.]])
        chosen, record = endpoint_assignment(primary, secondary)
        np.testing.assert_array_equal(chosen, [1, 0])
        self.assertEqual(record["primary_gap"], 0)

    def test_secondary_cannot_replace_unique_primary_optimum(self):
        primary = np.array([[0., 1e-8, 1.], [1e-8, 0., 1.], [1., 1., 0.]])
        secondary = np.array([[10., 0., 10.], [0., 10., 10.], [10., 10., 0.]])
        chosen, record = endpoint_assignment(primary, secondary)
        np.testing.assert_array_equal(chosen, [0, 1, 2])
        self.assertEqual(record["primary_gap"], 0)

    def test_pooled_null_variance_includes_shared_assignment_covariance(self):
        a = np.array([[1., 4.], [3., 2.]])
        b = 2 * a
        state = np.array([[4., 1.], [2., 3.]])
        pooled = pool_matrices({"a": a, "b": b}, {"a": 3., "b": 6.})
        moments = exact_cousin_stats(pooled, state, [[0, 1]])
        draws = np.array([[pooled[np.arange(2), p].sum(), state[np.arange(2), p].sum()]
                          for p in itertools.permutations(range(2))])
        np.testing.assert_allclose(moments[:2], draws.mean(axis=0))
        np.testing.assert_allclose(moments[2:4], draws.var(axis=0))
        self.assertAlmostEqual(moments[4], np.cov(draws.T, ddof=0)[0, 1])
        independent_average_variance = sum(exact_cousin_stats(x / total, state, [[0, 1]])[2]
                                           for x, total in ((a, 3.), (b, 6.))) / 4
        self.assertGreater(moments[2], independent_average_variance)

    def test_retention_counts_biological_parents_and_metrics_attained_points(self):
        travel = np.array([[2., 1., 4.], [1., 2., 4.], [4., 4., 1.]])
        state = np.array([[1., 3., 2.], [3., 1., 2.], [2., 2., 1.]])
        parents = np.array(["same", "same", "other"])
        result = analyze(travel, state, parents, [[0, 1, 2]], intervals=12, identity="toy")
        curve, metric = result["curve"], result["metrics"]
        self.assertAlmostEqual(curve.iloc[-1].edge_retention, 1.)
        np.testing.assert_allclose(curve.iloc[[-1, 0]][["D1", "D2"]].to_numpy(float), [[0, 1], [1, 0]])
        nearest = curve.iloc[metric["nearest_index"]]
        self.assertAlmostEqual(np.hypot(nearest.D1 - metric["natural_D1"], nearest.D2 - metric["natural_D2"]), metric["d_LP"])


if __name__ == "__main__":
    unittest.main()
