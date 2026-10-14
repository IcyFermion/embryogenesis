"""Checks for the cell-state representation criteria (Track A1b)."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from terminal_pareto import representation_criteria as rc

# A small lineage: root -> a, b; a -> aa, ap; b -> ba, bp; aa -> aaa, aap; ap -> apa, app.
TREE = {"did": "r", "children": [
    {"did": "a", "children": [{"did": "aa", "children": [{"did": "aaa"}, {"did": "aap"}]},
                              {"did": "ap", "children": [{"did": "apa"}, {"did": "app"}]}]},
    {"did": "b", "children": [{"did": "ba"}, {"did": "bp"}]}]}
LINEAGE = rc.build_lineage(TREE, lambda name: name)
TERMS = [("aaa", "aa"), ("aap", "aa"), ("apa", "ap"), ("app", "ap"), ("ba", "b"), ("bp", "b")]


class LineageTests(unittest.TestCase):
    def test_path_length_counts_divisions(self):
        self.assertEqual(LINEAGE.path_length("aaa", "aap"), 2)  # sisters
        self.assertEqual(LINEAGE.path_length("aaa", "apa"), 4)  # first cousins
        self.assertEqual(LINEAGE.path_length("aaa", "aa"), 1)   # parent and child
        self.assertEqual(LINEAGE.path_length("aaa", "ba"), 5)

    def test_bins_merge_until_enough_pairs(self):
        self.assertEqual(rc.path_bins(np.array([2, 2, 4, 4, 4, 6]), min_pairs=2), [(2, 2), (4, 6)])  # short tail joins the last bin
        self.assertEqual(rc.path_bins(np.array([2, 4, 4, 6, 6, 6]), min_pairs=3), [(2, 4), (6, 6)])


class CriteriaTests(unittest.TestCase):
    def setUp(self):
        self.coarse = {c: "x" for c, _ in TERMS}
        self.fine = {"aaa": "m", "aap": "m", "apa": "m", "app": "n", "ba": "n", "bp": "n"}
        self.design = rc.build_design(TERMS, self.coarse, self.fine, LINEAGE, min_pairs=1)

    def test_design_pairs_and_sisters(self):
        same, diff = self.design.sisters["fine"]
        self.assertEqual(sorted(map(tuple, same)), [(0, 1), (4, 5)])
        self.assertEqual(sorted(map(tuple, diff)), [(2, 3)])
        self.assertEqual(len(self.design.pairs), 3 + 3)  # two types of three cells

    def test_within_type_permutation_keeps_type_composition(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            perm = rc.within_type_permutation(self.design, rng)
            for group in self.design.groups:
                self.assertEqual(sorted(perm[group]), sorted(group))

    def test_cohen_d(self):
        D = np.zeros((4, 4)); D[0, 1], D[2, 3] = 1.0, 3.0
        same, diff = np.array([[0, 1], [0, 1]]), np.array([[2, 3], [2, 3]])
        self.assertEqual(rc.cohen_d(D, same, diff), 0.0)  # zero pooled variance
        D[1, 0] = 2.0
        self.assertAlmostEqual(rc.cohen_d(D, np.array([[0, 1], [1, 0]]), np.array([[2, 3], [2, 3]])),
                               (3 - 1.5) / np.sqrt(0.5 / 2))

    def test_memory_detects_lineage_inherited_state(self):
        rng = np.random.default_rng(1)
        # Balanced binary tree of 6 levels; one fine type; profiles inherit a random walk along the lineage.
        tree = {"did": "r", "children": []}
        def grow(node, depth):
            if depth == 6:
                return
            node["children"] = [{"did": node["did"] + "a"}, {"did": node["did"] + "p"}]
            for child in node["children"]:
                grow(child, depth + 1)
        grow(tree, 0)
        lineage = rc.build_lineage(tree, lambda name: name)
        leaves = [n for n in lineage.order if len(n) == 7]
        terms = [(leaf, leaf[:-1]) for leaf in leaves]
        fine = {leaf: "t" for leaf in leaves}
        design = rc.build_design(terms, fine, fine, lineage, min_pairs=10)
        state = {"r": rng.normal(size=8) + 5}
        for name in sorted(lineage.order, key=len)[1:]:
            state[name] = state[name[:-1]] + rng.normal(scale=1.0, size=8)
        inherited = np.array([state[leaf] for leaf in leaves])
        row, curve = rc.criteria(inherited, design, seed=0, n_boot=10, n_perm=200)
        self.assertGreater(row["memory_excess"], 0.2)
        self.assertLess(row["memory_p"], 0.01)
        self.assertLess(curve[0]["observed_mean"], curve[0]["null_q05"])  # sisters closer than any shuffle
        exchangeable = rng.normal(size=inherited.shape) + 5
        row, _ = rc.criteria(exchangeable, design, seed=0, n_boot=10, n_perm=200)
        self.assertLess(abs(row["memory_excess"]), 0.1)

    def test_fate_only_state_replaces_children_by_group_mean(self):
        frame = pd.DataFrame(np.arange(1, 25, dtype=float).reshape(12, 2), columns=["f", "g"],
                             index=["aaa", "aap", "apa", "app", "ba", "bp", "aa", "ap", "b", "a", "r", "x"])
        groups = {"aaa": "m", "aap": "m", "apa": "m", "app": "n", "ba": "n", "bp": "n"}
        S = rc.state_matrix(frame, ["f", "g"], TERMS, groups)
        m = frame.loc[["aaa", "aap", "apa"]].mean().to_numpy()
        expected = 1 - frame.loc["aa"].to_numpy() @ m / (np.linalg.norm(frame.loc["aa"]) * np.linalg.norm(m))
        self.assertAlmostEqual(S[0, 0], expected)
        np.testing.assert_allclose(S[:, 0], S[:, 1])  # same group, same column
        measured = rc.state_matrix(frame, ["f", "g"], TERMS)
        self.assertFalse(np.allclose(measured[:, 0], measured[:, 1]))


class RunTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ctx = rc.load_context()

    def test_candidates_and_cohorts(self):
        c = self.ctx.candidates
        self.assertEqual(c["protein"]["top-20"], self.ctx.inputs.published_features["protein"])
        self.assertEqual(c["rna"]["top-20"], self.ctx.inputs.published_features["rna"])
        self.assertEqual((len(c["protein"]["HV pool"]), len(c["rna"]["HV pool"])), (105, 118))
        self.assertEqual((len(c["protein"]["whole panel"]), len(c["rna"]["whole panel"])), (210, 237))
        self.assertEqual(len(self.ctx.designs["CE protein"].cells), 275)
        self.assertEqual(self.ctx.designs["CE RNA"].cells, self.ctx.designs["CB RNA"].cells)
        self.assertEqual(len(self.ctx.designs["CE RNA"].cells), 187)

    def test_reproduction_gates(self):
        report = rc.check_reproduction(self.ctx)
        self.assertEqual(report["failures"], [])
        self.assertEqual(report["notebook"]["cosine"]["n_same"], 106)

    def test_small_run_writes_once_and_replays(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run"
            small = dict(n_boot=20, n_perm=20, n_perm_type=10, n_random=3)
            original = rc.load_context
            rc.load_context = lambda: self.ctx  # reuse the loaded inputs
            try:
                rc.write_run(out, n_jobs=4, **small)
                report = rc.validate_run(out, n_jobs=4, replay_random=4)
                with self.assertRaises(FileExistsError):
                    rc.write_run(out, n_jobs=1, **small)
            finally:
                rc.load_context = original
            self.assertEqual(report["criteria"], 9)
            self.assertEqual(report["fronts"], 30)
            fronts = pd.read_csv(out / "analysis/fronts.csv")
            top = fronts[(fronts.candidate == "top-20") & (fronts.cost == "measured")].set_index(["dataset", "geometry"])
            self.assertAlmostEqual(top.at[("CE protein", "raw3d"), "d_LP"], 0.015, places=3)


if __name__ == "__main__":
    unittest.main()
