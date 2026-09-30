"""Run with python -m unittest full_tree_pareto.test_pooled_analysis."""
import itertools
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.optimize import linear_sum_assignment

from full_tree_pareto import pooled_analysis as pa
from terminal_pareto.front_coordinates import EndpointTransform


class PooledTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c = pa.build_context()

    def test_cohort_and_round_partition(self):
        c = self.c
        self.assertEqual((len(c.names), len(c.internal), len(c.leaves), len(c.edges)), (978, 487, 491, 974))
        self.assertEqual(len(c.roots), 4)
        self.assertEqual(set(np.concatenate(c.layers)), set(c.edges))
        self.assertEqual(sum(map(len, c.layers)), len(c.edges))
        self.assertTrue(np.all(c.depths[c.edges]-c.depths[c.parents[c.edges]] == 1))

    def test_pooled_geometry_replay(self):
        c = self.c
        for alpha in (0., .5, 1.):
            parent = pa.layerwise(c, c.combined(alpha))
            totals = [np.linalg.norm(x[c.edges]-x[parent[c.edges]], axis=1).sum() for x in c.xyz]
            self.assertAlmostEqual(np.mean(totals/c.denominators), c.score(parent)[0], places=12)
        self.assertAlmostEqual(c.natural[0], 1., places=13)
        # Do not substitute distances in an averaged coordinate cloud.
        avg = (c.xyz/c.denominators[:, None, None]).mean(axis=0)
        self.assertGreater(abs(np.linalg.norm(avg[c.edges]-avg[c.parents[c.edges]], axis=1).sum()-1), 1e-4)

    def test_layerwise_sum_and_endpoints(self):
        c = self.c
        arrays = [pa.layerwise(c, c.combined(a)) for a in (0, .5, 1)]
        for parent in arrays:
            pa.validate_assignment(c, parent)
            subtotal = np.sum([[c.travel[ch, parent[ch]].sum(), c.state[ch, parent[ch]].sum()] for ch in c.layers], axis=0)
            np.testing.assert_allclose(subtotal, c.score(parent))
        for children in c.layers+[c.edges]:
            costs = np.array([[c.travel[children, p[children]].sum(), c.state[children, p[children]].sum()] for p in arrays])
            transform = EndpointTransform.from_endpoints(reference_analysis_id="test", travel_optimum_assignment_id="2", state_optimum_assignment_id="0", travel_optimum_costs=costs[2], state_optimum_costs=costs[0])
            x, y = transform.transform(costs[:, 0], costs[:, 1])
            np.testing.assert_allclose([x[2], y[2], x[0], y[0]], [0, 1, 1, 0])
            np.testing.assert_allclose(np.column_stack(transform.inverse(x, y)), costs)
            # Aggregate feasibility is the product of round-wise assignments,
            # not an unrestricted all-edge bipartite assignment.
            if len(children) == len(c.edges):
                continue
            for j, objective in enumerate((c.travel, c.state)):
                matrix = objective[np.ix_(children, c.parents[children])]
                r, col = linear_sum_assignment(matrix)
                self.assertAlmostEqual(matrix[r, col].sum(), costs[2 if j == 0 else 0, j], places=10)

    def test_rebuilding_constraints(self):
        c = self.c
        for solver in (pa.spanning, pa.top_down, pa.bottom_up):
            for alpha in (0, 1):
                parent = solver(c, c.combined(alpha))
                pa.validate_assignment(c, parent, fixed_roots=solver is not pa.bottom_up)
        invalid = c.parents.copy()
        invalid[c.edges[0]] = c.edges[0]
        with self.assertRaises(ValueError):
            pa.validate_assignment(c, invalid)

    def test_cousin_exact_moments(self):
        c = self.c
        variances = np.zeros(2)
        for group in c.cousin_groups:
            leaves = c.leaves[group]
            slots = c.parents[leaves]
            totals = [[c.travel[leaves, list(p)].sum(), c.state[leaves, list(p)].sum()]
                      for p in itertools.permutations(slots)]
            variances += np.var(totals, axis=0)
        np.testing.assert_allclose(variances, c.scales**2, rtol=1e-10)

    def test_clamping_is_not_conditioning(self):
        # Root -> internal -> leaf: free internal mean stays at root, not the
        # endpoint bridge mean. Root and leaf remain exact in every sample.
        values = np.array([[0.], [1.], [4.]])
        parents = np.array([-1, 0, 1])
        totals, covariance, states = pa.forward_clamped(values, parents, np.ones(3), np.array([2]), 40000, np.random.default_rng(5), return_states=True)
        self.assertAlmostEqual(covariance[0, 0], 5.)
        self.assertLess(abs(states[:, 1].mean()), .04)
        self.assertLess(abs(states[:, 1].var()-5), .1)
        self.assertTrue(np.all(states[:, 0] == 0))
        self.assertTrue(np.all(states[:, 2] == 4))
        np.testing.assert_allclose(totals, np.abs(states[:, 1, 0])+np.abs(4-states[:, 1, 0]))

    def test_gaussian_cohort_and_determinism(self):
        c = self.c
        a, components, _ = pa.gaussian(c, 12, 242)
        b, _, _ = pa.gaussian(c, 12, 242)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_allclose(a[:, 0], (components[:, :3]/c.denominators).mean(axis=1))
        np.testing.assert_array_equal(a[:, 1], components[:, 3])
        # Check every actual leaf and root in a full 20-dimensional draw.
        _, _, states = pa.forward_clamped(c.expression, c.parents, np.ones(len(c.names)), c.leaves, 3, np.random.default_rng(10), return_states=True)
        fixed = np.r_[c.roots, c.leaves]
        np.testing.assert_array_equal(states[:, fixed], np.broadcast_to(c.expression[fixed], states[:, fixed].shape))

    def test_cache_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"cache.npz"
            identity = pa.cache_identity(self.c, 2, 10, 42)
            pa.save_cache(path, identity, data=np.arange(3))
            pa.verify_cache(path, identity)
            with self.assertRaises(ValueError):
                pa.verify_cache(path, identity | {"draws": 11})
            with path.open("ab") as f:
                f.write(b"changed")
            with self.assertRaises(ValueError):
                pa.verify_cache(path, identity)

    def test_paired_reconstruction_toy(self):
        parents = np.array([-1, 0, 0, 1, 1, 2, 2])
        x = np.array([0., 2., 5., 1., 3., 4., 6.])
        cost = np.abs(x[:, None]-x[None, :])
        c = replace(self.c, names=list("abcdefg"), parents=parents,
                    internal=np.array([0, 1, 2]), leaves=np.array([3, 4, 5, 6]),
                    edges=np.arange(1, 7), roots=np.array([0]), travel=cost,
                    state=cost**2, scales=np.ones(2))
        for alpha in (0, .5, 1):
            p = pa.paired(c, c.combined(alpha))
            pa.validate_assignment(c, p, fixed_roots=False)
            self.assertEqual(np.count_nonzero(p >= 0), 6)


if __name__ == "__main__":
    unittest.main()
