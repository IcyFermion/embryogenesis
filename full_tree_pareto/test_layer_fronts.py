import unittest

import numpy as np

from full_tree_pareto import layer_fronts as lfr
from full_tree_pareto import pooled_analysis as pa


class LayerFrontsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ctx, cls.zscore, cls.sets, cls.permutations, cls.references, _ = lfr.load()

    def test_sets_follow_strata(self):
        e1, _ = lfr.load_sets()
        for r in range(1, 9):
            self.assertEqual(self.sets["top-20", r], self.ctx.identity["proteins"])
            self.assertEqual(self.sets["per-round", r], e1["5-8" if r >= 5 else str(r)])
            self.assertEqual(self.sets["shared", r], e1["shared"])

    def test_round_solver_matches_pooled_layerwise(self):
        state, scales = self.ctx.state, lfr.leaf_scales(self.ctx, self.ctx.state)
        np.testing.assert_allclose(scales, self.ctx.scales, rtol=1e-12)
        sweeps = [lfr.solve_round(self.ctx, state, scales, children) for children in self.ctx.layers]
        for index in (0, 150, 300):
            parent = pa.layerwise(self.ctx, self.ctx.combined(index / lfr.INTERVALS))
            for children, sweep in zip(self.ctx.layers, sweeps):
                np.testing.assert_array_equal(sweep[index], parent[children])

    def test_top20_reproduces_cached_rounds(self):
        sets = {k: v for k, v in self.sets.items() if k[0] == "top-20"}
        results, matrices = lfr.compute(self.ctx, self.zscore, sets, self.permutations)
        report = lfr.check_reproduction(self.ctx, results, matrices, self.references, sets)
        self.assertEqual(report["rounds"], len(self.ctx.layers))


if __name__ == "__main__":
    unittest.main()
