"""Checks for the apoptotic-cell sensitivity analysis (Track D4)."""

import unittest

import numpy as np
import pandas as pd

from terminal_pareto import apoptosis_sensitivity as ap
from terminal_pareto import subtree_symmetry as ss


class ControlTests(unittest.TestCase):
    def test_draws_are_distinct_and_reproducible(self):
        pools = {0: [10, 11, 12], 1: [10, 11, 12], 2: [12, 13, 14]}
        first = ap.draw_controls(None, pools, np.random.default_rng(1))
        again = ap.draw_controls(None, pools, np.random.default_rng(1))
        self.assertEqual(first, again)
        for seed in range(50):
            picked = ap.draw_controls(None, pools, np.random.default_rng(seed))
            self.assertEqual(len(set(picked)), 3)
            self.assertTrue(set(picked) <= {10, 11, 12, 13, 14})

    def test_shared_neighbours_fall_back_to_the_next_nearest(self):
        ranked = list(range(10, 10 + ap.POOL_SIZE + 3))
        pools = {d: ranked for d in range(ap.POOL_SIZE + 2)}
        picked = ap.draw_controls(None, pools, np.random.default_rng(0))
        self.assertEqual(len(set(picked)), ap.POOL_SIZE + 2)

    def test_no_candidate_left_raises(self):
        with self.assertRaises(AssertionError):
            ap.draw_controls(None, {0: [5], 1: [5]}, np.random.default_rng(0))

    def test_p_value_counts_ties_and_adds_one(self):
        control = pd.Series([1.0, 2.0, 3.0])
        self.assertEqual(ap.p_value(control, 2.0, lower_is_better=True), 3 / 4)
        self.assertEqual(ap.p_value(control, 0.5, lower_is_better=True), 1 / 4)
        self.assertEqual(ap.p_value(control, 3.0, lower_is_better=False), 2 / 4)


class StoredRunTests(unittest.TestCase):
    def test_stored_run_replays(self):
        out = ss.run_path(ap.DEFAULT_RUN)
        if not (out / "analysis/provenance.json").exists():
            self.skipTest("no stored apoptosis run")
        tables = ap.verify_run(out)
        whole = tables["whole_cohort"].query("variant == 'drop_deaths' and travel == 'pooled'")
        self.assertEqual(set(whole.n), {245})
        self.assertEqual(set(whole.n_deaths), {0})
        self.assertEqual(len(tables["controls"]), ap.N_DRAWS * 2)


if __name__ == "__main__":
    unittest.main()
