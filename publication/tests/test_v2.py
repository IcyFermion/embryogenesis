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

    def test_every_figure_of_the_release_set_has_a_caption(self):
        from publication.v2 import captions, manuscript
        captioned = {*captions.AGREED, *captions.DRAFT}
        self.assertEqual(set(manuscript.MAIN), captioned)
        self.assertFalse(set(manuscript.ADOPTED) & captioned)
        self.assertEqual({n.replace("--", "-") for n in manuscript.BACKLOG_V2}, set(captions.BACKLOG))
        from publication.v2.__main__ import figure_order
        self.assertEqual(sorted(["figS2", "fig7", "fig1"], key=figure_order), ["fig1", "fig7", "figS2"])

    def test_release_set_renumbers_released_wrappers(self):
        from publication.v2.manuscript import renumber
        tex = ("\\begin{document}\n\\renewcommand{\\thefigure}{S3}\n\\begin{figure}[p]\n"
               "\\includegraphics[width=\\textwidth]{figS3_four_geometries.pdf}\n"
               "\\includegraphics{figS3A_panel.pdf}\\input{figS30_other.tex}\n\\end{figure}\n")
        new, names = renumber(tex, "figS3_ce_tracking_robustness", "figS1_ce_tracking_robustness", "S1")
        self.assertIn("\\renewcommand{\\thefigure}{S1}\n", new)
        self.assertNotIn("{S3}", new)
        self.assertEqual(names, {"figS3_four_geometries.pdf": "figS1_four_geometries.pdf",
                                 "figS3A_panel.pdf": "figS1A_panel.pdf", "figS30_other.tex": "figS30_other.tex"})
        self.assertIn("{figS1_four_geometries.pdf}", new)
        # A figure that keeps its number keeps its files.
        tex = "\\setcounter{figure}{5}\n\\includegraphics{fig6A_ce_retention_heatmap.pdf}\n"
        new, names = renumber(tex, "fig6_ce_cell_types", "fig6_ce_cell_types", "6")
        self.assertIn("\\setcounter{figure}{5}\n", new)
        self.assertEqual(names, {"fig6A_ce_retention_heatmap.pdf": "fig6A_ce_retention_heatmap.pdf"})
        with self.assertRaises(ValueError):
            renumber("\\includegraphics{figS3_x.pdf}\n", "figS3_a", "figS1_a", "S1")


if __name__ == "__main__":
    unittest.main()
