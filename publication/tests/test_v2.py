import unittest

from publication.v2.data import tree_distance


class V2Tests(unittest.TestCase):
    def test_tree_distance_counts_divisions_through_the_common_ancestor(self):
        index = {"r": dict(ancestors=[]), "a": dict(ancestors=["r"]), "b": dict(ancestors=["r"]),
                 "a1": dict(ancestors=["a", "r"]), "b1": dict(ancestors=["b", "r"])}
        self.assertEqual(tree_distance(index, "a1", "a1"), 0)
        self.assertEqual(tree_distance(index, "a1", "a"), 1)
        self.assertEqual(tree_distance(index, "a1", "b"), 3)
        self.assertEqual(tree_distance(index, "a1", "b1"), 4)

    def test_palette_holds_one_valid_colour_per_concept(self):
        from matplotlib.colors import is_color_like
        from publication import palette
        for group, colors in palette.GROUPS.items():
            for name, color in colors.items():
                self.assertTrue(is_color_like(color), f"{group}.{name}")
        # Each measure has its own colour.
        self.assertEqual(len(set(palette.MEASURE.values())), len(palette.MEASURE))
        # The pie legend and the subtree figure need every fate class of the cohort.
        from publication.figures.terminal_subtree_map import TYPE_COLORS
        self.assertEqual(set(palette.FATE), set(TYPE_COLORS))

    def test_subtree_order_keeps_counterparts_adjacent(self):
        from publication.v2.panels import SUBTREE_ORDER, SUBTREE_PAIRS
        for a, b in SUBTREE_PAIRS:
            group = next(names for names in SUBTREE_ORDER.values() if a in names)
            self.assertIn(b, group)
            self.assertEqual(abs(group.index(a) - group.index(b)), 1, (a, b))
        order = [n for names in SUBTREE_ORDER.values() for n in names]
        self.assertEqual(len(order), len(set(order)))

    def test_subtree_distances_share_one_axis_without_r(self):
        import pandas as pd
        from publication.v2.panels import subtree_metrics
        canonical = pd.DataFrame(dict(u_lineage_lp=[.2, 0.], d_lp=[.06, 0.], d_np=[.85, .3],
                                      relative_distance=[.24, 0.]))
        columns = subtree_metrics(canonical)
        fields = [[field for field, _ in column.series] for column in columns]
        # Author, 2026-11-10: d_LP and d_CP on the same axis, r no longer drawn.
        self.assertEqual(fields, [["u_lineage_lp"], ["d_lp", "d_np"]])
        self.assertEqual(columns[1].xmax, 1.)

    def test_every_figure_of_the_manuscript_set_has_a_caption(self):
        from publication.v2 import captions, manuscript
        captioned = [*captions.AGREED, *captions.DRAFT, *captions.RELEASED]
        # Every figure of the free-standing set comes from the v2 build with exactly one caption (2026-11-13).
        self.assertEqual(sorted(manuscript.MAIN), sorted(captioned))
        self.assertEqual({n.replace("--", "-") for n in manuscript.BACKLOG_V2}, set(captions.BACKLOG))
        from publication.v2.__main__ import figure_order
        self.assertEqual(sorted(["figS2", "fig7", "fig1"], key=figure_order), ["fig1", "fig7", "figS2"])

    def test_presentation_sources_cover_build_captions_and_set(self):
        from publication.v2.__main__ import presentation_sources
        sources = presentation_sources()
        for name in ("__main__.py", "panels.py", "data.py", "captions.py", "manuscript.py", "assets/null_models.tex"):
            self.assertIn(f"publication/v2/{name}", sources)
        self.assertIn("publication/palette.py", sources)
        self.assertNotIn("publication/tests/test_v2.py", sources)

    def test_set_verification_rejects_missing_extra_and_modified_files(self):
        import tempfile
        from pathlib import Path
        from publication.provenance import sha256
        from publication.v2.manuscript import check_files
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "backlog").mkdir()
            (root / "fig1.tex").write_text("a")
            (root / "backlog/fig8.tex").write_text("b")
            (root / "manifest.json").write_text("{}")
            files = {name: sha256(root / name) for name in ("fig1.tex", "backlog/fig8.tex")}
            check_files(root, files)
            (root / "fig1.tex").write_text("edited")
            with self.assertRaisesRegex(ValueError, "modified"):
                check_files(root, files)
            (root / "fig1.tex").write_text("a")
            (root / "fig1.aux").write_text("")
            with self.assertRaisesRegex(ValueError, "unexpected"):
                check_files(root, files)
            (root / "fig1.aux").unlink()
            (root / "backlog/fig8.tex").unlink()
            with self.assertRaisesRegex(ValueError, "missing"):
                check_files(root, files)

    def test_promotion_repoints_current_and_keeps_history(self):
        import json
        import os
        import tempfile
        from pathlib import Path
        from publication.v2.manuscript import point_current
        with tempfile.TemporaryDirectory() as tmp:
            sets = Path(tmp)
            for name in ("set_a", "set_b"):
                (sets / name).mkdir()
            point_current(sets, sets / "set_a", dict(commit="a"))
            entry = point_current(sets, sets / "set_b", dict(commit="b"))
            self.assertEqual(os.readlink(sets / "current"), "set_b")
            self.assertEqual(entry["previous"], "set_a")
            history = json.loads((sets / "history.json").read_text())
            self.assertEqual([h["set"] for h in history], ["set_a", "set_b"])
            with self.assertRaises(ValueError):
                point_current(sets / "set_a", sets / "set_b", {})


if __name__ == "__main__":
    unittest.main()
