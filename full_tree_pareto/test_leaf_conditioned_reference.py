"""Exact small-tree and publication-scope checks for Gaussian conditioning."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from full_tree_pareto import publication_analysis as analysis
from full_tree_pareto.leaf_conditioned_reference import ConditionalBlock, prepare, compare, plot


class LeafConditioningTests(unittest.TestCase):
    def test_chain_bridge_mean_and_joint_covariance(self):
        block = ConditionalBlock([0, 1, 2], [1, 2, 3], [1., 1., 1.],
                                 np.array([[0., 0.], [90., -80.], [40., 30.], [3., 6.]]),
                                 [True, False, False, True], [[2., .5], [.5, 1.]])
        np.testing.assert_allclose(block.mean, [[0, 0], [1, 2], [2, 4], [3, 6]])
        np.testing.assert_allclose(block.node_covariance, [[2/3, 1/3], [1/3, 2/3]])
        np.testing.assert_allclose(block.edge_variance, [2/3]*3)
        samples = block.sample(30_000, np.random.default_rng(456))
        flattened = samples[:, block.free].reshape(30_000, -1)
        np.testing.assert_allclose(flattened.mean(axis=0), block.mean[block.free].ravel(), atol=.025)
        np.testing.assert_allclose(np.cov(flattened, rowvar=False),
                                   np.kron(block.node_covariance, block.covariance), atol=.035)
        squared = ((samples[:, block.children]-samples[:, block.parents])**2).sum(axis=(1, 2))
        self.assertLess(abs(squared.mean()-block.analytic_squared_total()),
                        6*squared.std()/np.sqrt(len(squared)))

    def test_branching_conditions_on_both_daughters(self):
        # One free parent, three fixed neighbors with conductances 1/2, 1, 1/4.
        block = ConditionalBlock([0, 1, 1], [1, 2, 3], [2., 1., 4.],
                                 np.array([[0.], [999.], [3.], [-2.]]),
                                 [True, False, True, True], [[1.]])
        self.assertAlmostEqual(block.mean[1, 0], (3-2/4)/(1/2+1+1/4))
        self.assertAlmostEqual(block.node_covariance[0, 0], 1/(1/2+1+1/4))

    def test_root_only_recovers_independent_increment_law(self):
        block = ConditionalBlock([0, 1, 1], [1, 2, 3], [2., 1., 4.],
                                 np.array([[5.], [8.], [9.], [10.]]),
                                 [True, False, False, False], [[3.]])
        np.testing.assert_allclose(block.mean, 5)
        np.testing.assert_allclose(block.edge_variance, [2, 1, 4])
        samples = block.sample(30_000, np.random.default_rng(6))
        increments = (samples[:, block.children]-samples[:, block.parents])[:, :, 0]
        np.testing.assert_allclose(np.cov(increments, rowvar=False), np.diag([6, 3, 12]), atol=.2)

    def test_free_observations_do_not_enter_conditional_mean(self):
        args = ([0, 1], [1, 2], [1., 2.])
        fixed = [True, False, True]
        first = ConditionalBlock(*args, np.array([[2.], [99.], [8.]]), fixed, [[1.]])
        second = ConditionalBlock(*args, np.array([[2.], [-999.], [8.]]), fixed, [[1.]])
        np.testing.assert_array_equal(first.mean, second.mean)
        np.testing.assert_array_equal(first.sample(10, np.random.default_rng(8)),
                                      second.sample(10, np.random.default_rng(8)))

    def test_invalid_clock_and_draw_count(self):
        with self.assertRaises(ValueError):
            ConditionalBlock([0, 1], [1, 2], [0., 1.], np.zeros((3, 1)),
                             [True, False, True], [[1.]])
        block = ConditionalBlock([0, 1], [1, 2], [1., 1.], np.zeros((3, 1)),
                                 [True, False, True], [[1.]])
        with self.assertRaises(ValueError):
            block.sample(0, np.random.default_rng(1))

    def test_publication_scope_and_exact_both_block_boundaries(self):
        context = analysis.load_ce_protein_context(max_workers=1)
        reference = analysis.separate_clock_reference(context, n_samples=2, seed=242)
        edges, nodes, blocks = prepare(context, reference.covariance)
        self.assertEqual(len(edges), 1000)
        self.assertEqual(int(nodes.fixed_boundary.sum()), 508)
        self.assertEqual(int(nodes.terminal_in_measured_tree.sum()), 504)
        for name, block in blocks.items():
            self.assertEqual(block.free.sum(), 496)
            states = block.sample(20, np.random.default_rng(19))
            for draw in states:
                np.testing.assert_array_equal(draw[block.fixed], block.values[block.fixed])
            observed = np.linalg.norm(block.values[block.children]-block.values[block.parents], axis=1)
            np.testing.assert_allclose(observed, edges[f"{name}_norm"])
            # Leaf conditioning cannot add marginal node uncertainty relative
            # to the same process conditioned only on roots.
            root_only = ConditionalBlock(block.parents, block.children, block.clock, block.values,
                                         nodes.optimization_root.to_numpy(bool), block.covariance)
            prior_indices = np.flatnonzero(root_only.free)
            selected = np.flatnonzero(np.isin(prior_indices, np.flatnonzero(block.free)))
            prior = root_only.node_covariance[np.ix_(selected, selected)]
            self.assertGreater(np.linalg.eigvalsh(prior-block.node_covariance).min(), -1e-8)

    def test_end_to_end_outputs_and_plot(self):
        context = analysis.load_ce_protein_context(max_workers=1)
        outputs = compare(context, n_draws=20, seed=5)
        self.assertEqual(len(outputs["summary.csv"]), 4)
        self.assertEqual(len(outputs["draws.csv"]), 40)
        self.assertEqual(int(outputs["conditional_mean_states.csv"].fixed_boundary.sum()), 508)
        with TemporaryDirectory() as directory:
            plot(outputs, Path(directory))
            self.assertGreater((Path(directory)/"leaf_conditioned_comparison.png").stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
