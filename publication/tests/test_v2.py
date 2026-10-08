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


if __name__ == "__main__":
    unittest.main()
