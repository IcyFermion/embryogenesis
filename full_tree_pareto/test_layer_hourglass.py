import unittest

import numpy as np

from full_tree_pareto import layer_hourglass as lh
from full_tree_pareto import round_references as rr
from terminal_pareto.subtree_analysis import exact_cousin_stats


class UnitTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.blocks = rng.random((7, 7)), rng.random((7, 7))
        self.groups = [np.array([0, 1, 2]), np.array([4, 5])]
        self.units = self.groups + [np.array([3]), np.array([6])]

    def test_unit_sums_equal_exact_shuffle_mean(self):
        table = lh.unit_costs(self.blocks, self.units)
        exact = exact_cousin_stats(*self.blocks, self.groups)
        np.testing.assert_allclose(table.sum(axis=0), [np.trace(self.blocks[0]), np.trace(self.blocks[1]), *exact[:2]],
                                   rtol=1e-12)

    def test_singletons_carry_no_advantage(self):
        table = lh.unit_costs(self.blocks, [np.array([3])])
        np.testing.assert_allclose(lh.advantage(table[0]), [0, 0], atol=1e-15)

    def test_identity_resample_is_the_point_estimate(self):
        table = lh.unit_costs(self.blocks, self.units)
        identity = np.arange(len(self.units))[None, :]
        np.testing.assert_allclose(lh.bootstrap(table, identity)[0], lh.advantage(table.sum(axis=0)), rtol=1e-12)

    def test_layer_numbering(self):
        self.assertEqual([lh.layer_of(r) for r in (1, 2, 3, 8)], [8, 7, 6, 1])


class CacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trees, _ = lh.load()

    def test_units_partition_every_round(self):
        index = rr.tree_index()
        for tree in self.trees.values():
            for children in tree["layers"]:
                for family in lh.FAMILIES:
                    members = np.concatenate(lh.units(tree["names"], children, index, family))
                    np.testing.assert_array_equal(np.sort(members), np.arange(len(children)))

    def test_reproduces_round_references(self):
        results, _ = lh.compute(self.trees, draws=20)
        report = lh.check_reproduction(results["layers"], self.trees)
        self.assertEqual(report["rows_checked"], 8 + 3 * len(self.trees["partial"]["layers"]))


if __name__ == "__main__":
    unittest.main()
