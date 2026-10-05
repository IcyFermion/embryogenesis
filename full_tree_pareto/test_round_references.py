import json
import unittest

import numpy as np
import pandas as pd

from full_tree_pareto import cross_species_analysis as fcs
from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import round_references as rr


class RoundReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ctx = pa.build_context()
        cls.index = rr.tree_index()

    def test_round_draws_keep_slots_and_score_only_the_round(self):
        ctx = self.ctx
        for r, children in enumerate(ctx.layers, 1):
            blocks = rr.round_blocks(ctx.travel, ctx.state, ctx.parents, children)
            natural = np.array([np.trace(blocks[0]), np.trace(blocks[1])])
            identity = np.arange(len(children))[None, :]
            np.testing.assert_allclose(rr.score_permutations(blocks, identity)[0], natural)
            for family in rr.FAMILIES:
                groups = rr.round_groups(ctx.names, children, self.index, family)
                perms = rr.sample_permutations(len(children), groups, 5, r)
                rr.validate_permutations(perms, len(children), groups)
                for perm, cost in zip(perms, rr.score_permutations(blocks, perms)):
                    parent = ctx.parents.copy()
                    parent[children[perm]] = ctx.parents[children]
                    pa.validate_assignment(ctx, parent)
                    np.testing.assert_allclose(ctx.score(parent) - ctx.natural, cost - natural, atol=1e-9)

    def test_exact_moments_match_draws(self):
        ctx = self.ctx
        children = ctx.layers[1]
        blocks = rr.round_blocks(ctx.travel, ctx.state, ctx.parents, children)
        groups = rr.round_groups(ctx.names, children, self.index, "first_cousin")
        stats = rr.exact_moments(blocks, groups)
        draws = rr.score_permutations(blocks, rr.sample_permutations(len(children), groups, 20000, 7))
        self.assertTrue(np.all(np.abs(draws.mean(axis=0) - stats[:2]) < 5*np.sqrt(stats[2:4]/len(draws))))
        np.testing.assert_allclose(draws.var(axis=0), stats[2:4], rtol=0.05)

    def test_validator_rejects_crossing_groups(self):
        groups = [np.array([0, 1]), np.array([2, 3])]
        with self.assertRaises(ValueError):
            rr.validate_permutations(np.array([[2, 1, 0, 3]]), 4, groups)

    def test_landmarks_reproduce_cached_partial_forest_metrics(self):
        run = fcs.run_path(fcs.DEFAULT_RUN)
        analysis = run / "analysis"
        record = json.loads((analysis / "provenance.json").read_text())
        metrics = pd.read_csv(analysis / "metrics.csv")
        with np.load(analysis / "matrices.npz") as data:
            arrays = {k: data[k] for k in data.files}
        for geometry in fcs.GEOMETRIES:
            for config in fcs.PRIMARY:
                key = f"{geometry}__{config}"
                ctx = fcs.forest(record["names"], arrays["parents"], arrays[f"{key}__travel"], arrays[f"{key}__state"])
                cached = fcs.load_cache(analysis / f"{key}.npz", dict(analysis_id=record["analysis_id"], geometry=geometry,
                                        config=config, settings=record["settings"]))
                saved = metrics[(metrics.geometry == geometry) & (metrics.config == config) & (metrics.cohort == "base")].iloc[0]
                retention = (cached["parents"][:, ctx.edges] == ctx.parents[ctx.edges]).mean(axis=1)
                null = [saved.null_mean_travel, saved.null_mean_state]
                _, mine, _ = rr.landmarks(cached["costs"], ctx.natural, null, retention, key)
                for field in ("u_L", "d_LP", "d_NP", "max_edge_retention", "retention_at_nearest"):
                    self.assertAlmostEqual(mine[field], saved[field], places=10, msg=f"{key} {field}")

    def test_saved_caches_load(self):
        for run, target in ((pa.DEFAULT_RUN, "pooled"), (fcs.run_path(fcs.DEFAULT_RUN), "cross_species")):
            if not (run / rr.SUBDIR / "validation.json").exists():
                self.skipTest("Round references not built")
            loaded = rr.load(run, target)
            self.assertTrue(np.isfinite(np.concatenate(list(loaded["costs"].values()))).all())
            self.assertEqual(sorted(loaded["metrics"]["round"].unique()), list(range(1, len(loaded["identity"]["rounds"]) + 1)))


if __name__ == "__main__":
    unittest.main()
