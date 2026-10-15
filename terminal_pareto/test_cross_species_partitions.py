import unittest

import numpy as np

from terminal_pareto import cross_species_partitions as tcp


class PartitionTests(unittest.TestCase):
    def test_child_retention_follows_slot_convention(self):
        parents = np.array(["A", "A", "B", "B"])
        identity = np.arange(4)[None, :]
        np.testing.assert_array_equal(tcp.child_retention(parents, identity), np.ones((1, 4), bool))
        # Slot 0 (parent A) receives child 2 (natural parent B) and slot 2 receives child 0.
        swapped = np.array([[2, 1, 0, 3], [1, 0, 3, 2]])
        np.testing.assert_array_equal(tcp.child_retention(parents, swapped),
                                      [[False, True, False, True], [True, True, True, True]])

    def test_subtrees_and_type_partition(self):
        index = {"t1": dict(ancestors=["p", "g", "r"], depth=3), "t2": dict(ancestors=["p", "g", "r"], depth=3),
                 "t3": dict(ancestors=["q", "g", "r"], depth=3), "t4": dict(ancestors=["s", "h", "r"], depth=3),
                 "g": dict(depth=1), "r": dict(depth=0), "p": dict(depth=2), "q": dict(depth=2), "s": dict(depth=2),
                 "h": dict(depth=1)}
        blocks = tcp.subtrees(np.array(["t1", "t2", "t3", "t4"]), index, min_cells=2)
        self.assertEqual(list(blocks), ["r", "g", "p"])
        np.testing.assert_array_equal(blocks["g"], [0, 1, 2])
        types = tcp.type_partition(np.array(["neuron"] * 3 + ["glial", "other"]), min_cells=3)
        self.assertEqual(list(types), ["neuron"] * 3 + ["other", "other"])

    def test_saved_partitions_load(self):
        try:
            loaded = tcp.load()
        except FileNotFoundError:
            self.skipTest("Partitions not built")
        whole = loaded["metrics"][(loaded["metrics"].kind == "subtree") & (loaded["metrics"].name == "P0")]
        self.assertEqual(set(whole.n), {187})
        self.assertFalse(whole.degenerate.any())


if __name__ == "__main__":
    unittest.main()
