"""Checks for the left-right subtree pairing (Track D0)."""

import unittest

from terminal_pareto import subtree_symmetry as ss

# Two mirrored subtrees Xl and Xr, plus one unrelated subtree Y.
CELLS = {
    "Xla": "ASEL", "Xlpa": "MU BOD", "Xlpp": "DA1",
    "Xra": "ASER", "Xrpa": "MU BOD", "Xrpp": "EXC CELL",
    "Ya": "HYP2DL", "Yp": "XRAY",
}
DESCS = {c: "described" for c in CELLS}
DESCS["Ya"] = "Postembryonic blast cell; stochastic fate (aka P1/2 L)"
LINEAGE = ss.Lineage(
    terminals=sorted(CELLS),
    descendants={"Xl": ["Xla", "Xlpa", "Xlpp"], "Xr": ["Xra", "Xrpa", "Xrpp"],
                 "Xlp": ["Xlpa", "Xlpp"], "Xrp": ["Xrpa", "Xrpp"], "Y": ["Ya", "Yp"]})


def classes(deaths=()):
    return ss.classify(sorted(CELLS), CELLS, DESCS, set(deaths))


class ClassTests(unittest.TestCase):
    def test_kinds(self):
        cls = classes(deaths=["Yp"])
        self.assertEqual((cls["Xla"]["fate_class"], cls["Xla"]["side"], cls["Xla"]["kind"]), ("ASE", "L", "sided"))
        self.assertEqual(cls["Xlpp"]["fate_class"], "DA")
        self.assertEqual(cls["Xlpp"]["kind"], "serial")
        self.assertEqual(cls["Xrpp"]["kind"], "unsided")
        self.assertEqual((cls["Ya"]["fate_class"], cls["Ya"]["side"]), ("P1/2", "L"))  # alias beats the name
        self.assertEqual(cls["Yp"]["kind"], "death")

    def test_side_needs_the_swapped_name(self):
        cls = ss.classify(["a", "b"], {"a": "HYP2DL", "b": "AVL"}, {"a": "x", "b": "x"}, set())
        self.assertEqual(cls["a"]["side"], None)
        self.assertEqual(cls["b"]["side"], None)

    def test_unannotated(self):
        cls = ss.classify(["ABplapapaa"], {"ABplapapaa": "ABPLAPAPAA"}, {"ABplapapaa": float("nan")}, set())
        self.assertEqual(cls["ABplapapaa"]["kind"], "unannotated")

    def test_compatibility(self):
        cls = classes()
        self.assertTrue(ss.compatible(cls["Xla"], cls["Xra"]))
        self.assertFalse(ss.compatible(cls["Xla"], cls["Xla"]))       # same side
        self.assertTrue(ss.compatible(cls["Xlpa"], cls["Xrpa"]))      # unsided, same class
        self.assertFalse(ss.compatible(cls["Xlpp"], cls["Xrpp"]))


class MatchingTests(unittest.TestCase):
    def setUp(self):
        self.cls = classes()

    def test_fate_matching(self):
        self.assertEqual(ss.fate_matching(LINEAGE.descendants["Xl"], LINEAGE.descendants["Xr"], self.cls), (2, 1))
        self.assertEqual(ss.fate_matching(["Xla", "Xra"], ["Xra"], self.cls), (1, 1))

    def test_position_match_and_leftovers(self):
        mode, frac = ss.position_match("Xl", "Xr", LINEAGE, self.cls)
        self.assertEqual(mode, "same")
        self.assertAlmostEqual(frac, 4 / 6)
        pos, leftover = ss.mirror_pairs("Xl", "Xr", LINEAGE, self.cls, mode)
        self.assertEqual(sorted(pos), [("Xla", "Xra"), ("Xlpa", "Xrpa")])
        self.assertEqual(leftover, {"Xlpp", "Xrpp"})

    def test_swap_phi(self):
        lineage = ss.Lineage(["Lal", "Lar", "Ral", "Rar"],
                             {"La": ["Lal", "Lar"], "Ra": ["Ral", "Rar"]})
        names = {"Lal": "AL", "Lar": "BL", "Ral": "BR", "Rar": "AR"}
        cls = ss.classify(sorted(names), names, {c: "x" for c in names}, set())
        self.assertEqual(ss.position_match("La", "Ra", lineage, cls), ("swap_lr", 1.0))

    def test_pair_metrics(self):
        m = ss.pair_metrics("Xl", "Xr", LINEAGE, self.cls)
        self.assertAlmostEqual(m["fate_match"], 4 / 6)
        self.assertEqual(m["sided_match"], 1.0)
        self.assertEqual((m["laterality_x"], m["laterality_y"]), (1.0, -1.0))
        self.assertEqual(m["self_symmetry_x"], 0.0)


class LineageTests(unittest.TestCase):
    def test_parent_of_founders(self):
        self.assertEqual(ss.parent_of("MS"), "EMS")
        self.assertEqual(ss.parent_of("Z2"), "P4")
        self.assertEqual(ss.parent_of("ABplpa"), "ABplp")
        self.assertIsNone(ss.parent_of("P0"))


class CompareTests(unittest.TestCase):
    def test_null_support_ignores_sister_only_groups(self):
        parents = ["p", "p", "q", "r", "r"]
        # Group 1: sisters only (one parent). Group 2: cousins (two parents).
        self.assertEqual(ss.null_support(None, parents, [[0, 1], [2, 3, 4]]), 3)

    def test_depth(self):
        self.assertEqual(ss.depth("P0"), 0)
        self.assertEqual(ss.depth("ABplpa"), 5)  # AB, ABp, ABpl, ABplp, ABplpa
        self.assertEqual(ss.depth("MSaa"), 5)


class StoredRunTests(unittest.TestCase):
    def check(self, run_id):
        out = ss.run_path(run_id)
        if not (out / "analysis/provenance.json").exists():
            self.skipTest(f"no stored run {run_id}")
        return ss.verify_run(out)

    def test_pairing_run_replays(self):
        self.check(ss.DEFAULT_RUN)

    def test_compare_run_replays(self):
        tables = self.check(ss.DEFAULT_COMPARE_RUN)
        tiers = tables["pair_tiers"].set_index(["subtree_x", "subtree_y"]).tier
        self.assertEqual(tiers[("ABplpa", "ABprpa")], "priority")
        self.assertEqual(tiers[("MSaa", "MSpa")], "priority")
        self.assertEqual(tiers[("ABpl", "ABpr")], "excluded")   # sided match 0.83 < 0.85


if __name__ == "__main__":
    unittest.main()
