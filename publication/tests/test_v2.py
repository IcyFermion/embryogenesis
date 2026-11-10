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


if __name__ == "__main__":
    unittest.main()
