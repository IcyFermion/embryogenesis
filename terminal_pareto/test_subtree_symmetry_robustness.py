"""Checks for the robustness and classification of paired subtrees (Track D2)."""

import unittest

from terminal_pareto import subtree_symmetry as ss
from terminal_pareto import subtree_symmetry_robustness as sr


class StoredRunTests(unittest.TestCase):
    def test_stored_run_replays(self):
        out = ss.run_path(sr.DEFAULT_RUN)
        if not (out / "analysis/provenance.json").exists():
            self.skipTest("no stored robustness run")
        tables = sr.verify_run(out)
        labels = tables["classification"].set_index(["subtree_x", "subtree_y"]).label
        self.assertEqual(labels[("ABplpa", "ABprpa")], "numerical")   # u only, tiny travel span
        self.assertEqual(labels[("MSaa", "MSpa")], "data")            # one-group null, late-born cells
        variants = tables["variant_metrics"]
        self.assertEqual(set(variants.status), {"ok", "degenerate_null", "degenerate_span"})


if __name__ == "__main__":
    unittest.main()
