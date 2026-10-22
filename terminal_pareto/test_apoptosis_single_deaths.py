"""Checks for one-death-at-a-time treatments in the posterior ABpl family (Track D4c)."""

import json
import unittest

import pandas as pd

from terminal_pareto import apoptosis_single_deaths as sd
from terminal_pareto import subtree_symmetry as ss


class MirrorTests(unittest.TestCase):
    def test_mirror_swaps_the_fourth_letter(self):
        self.assertEqual(sd.mirror("ABplpppapp"), "ABprpppapp")
        self.assertEqual(sd.mirror("ABprp"), "ABplp")
        self.assertEqual(sd.side("ABprpp"), "right")
        with self.assertRaises(ValueError):
            sd.mirror("MSaa")

    def test_subtrees_pair_left_with_right(self):
        self.assertEqual(sd.SUBTREES[4:], tuple(sd.mirror(s) for s in sd.SUBTREES[:4]))


class StoredRunTests(unittest.TestCase):
    def setUp(self):
        self.analysis = ss.run_path(sd.DEFAULT_RUN) / "analysis"
        if not (self.analysis / "provenance.json").exists():
            self.skipTest("no stored single-death run")

    def test_stored_tables_match_their_hashes_and_readouts_hold(self):
        record = json.loads((self.analysis / "provenance.json").read_text())
        for name in sd.TABLES:
            self.assertEqual(ss.file_hash(self.analysis / f"{name}.csv"), record["table_hashes"][name], name)
        readouts = pd.read_csv(self.analysis / "readouts.csv").set_index("readout").holds
        self.assertTrue(readouts.all())

    def test_stored_run_replays(self):
        tables = sd.verify_run(self.analysis.parent)
        self.assertEqual(set(tables["single"].subtree), set(sd.SUBTREES))


if __name__ == "__main__":
    unittest.main()
