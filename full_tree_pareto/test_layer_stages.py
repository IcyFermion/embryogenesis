import unittest

import numpy as np
import pandas as pd

from full_tree_pareto import layer_stages as ls


class StageTests(unittest.TestCase):
    def test_alive_counts_births_and_divisions(self):
        birth = pd.Series([0., 10., 10., 20.])
        death = pd.Series([10., np.nan, 20., np.nan])
        np.testing.assert_array_equal(ls.alive_at(np.array([0., 10., 20., np.nan]), birth, death), [1, 2, 2, -1])

    def test_stage_bins(self):
        self.assertEqual([ls.stage_of(t) for t in (10., 56., 180.9, 181., 400.)], ["8", "15", "100", "200", "600"])
        self.assertIsNone(ls.stage_of(np.nan))

    def test_layer_birth_stages(self):
        _, layers = ls.birth_tables()
        pooled = layers[layers.tree == "pooled"].set_index("layer")
        np.testing.assert_array_equal(pooled.alive_median.loc[[8, 7, 6, 5, 4]], [370, 192, 126, 73, 38])
        self.assertTrue((layers.cells == layers.edges).all())


if __name__ == "__main__":
    unittest.main()
