"""Checks for the subtree-level matched control (Track D4b)."""

import unittest

import pandas as pd

from terminal_pareto import apoptosis_subtree_control as sc
from terminal_pareto import subtree_symmetry as ss


class RuleTests(unittest.TestCase):
    def test_float_noise_counts_as_a_tie(self):
        control = pd.Series([0.0, 4.5e-15, 0.3])
        self.assertEqual(sc.p_value(control, 0.0, lower_is_better=True), 3 / 4)
        self.assertEqual(sc.p_value(pd.Series([1 - 1e-15, 0.9]), 1.0, lower_is_better=False), 2 / 3)

    def test_roles(self):
        self.assertEqual(sc.role("ABplp"), "primary")
        self.assertEqual(sc.role("ABplppp"), "supporting")
        self.assertEqual(sc.role("ABarp"), "descriptive")


class StoredRunTests(unittest.TestCase):
    def test_stored_run_replays(self):
        out = ss.run_path(sc.DEFAULT_RUN)
        if not (out / "analysis/provenance.json").exists():
            self.skipTest("no stored subtree control run")
        tables = sc.verify_run(out)
        tested = set(tables["observed"].query("tested").subtree)
        self.assertIn(sc.PRIMARY, tested)
        self.assertNotIn("ABalpp", tested)
        self.assertEqual(len(tables["controls"]), sc.N_DRAWS * 2 * len(tested))


if __name__ == "__main__":
    unittest.main()
