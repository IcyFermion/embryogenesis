import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from full_tree_pareto import layer_features as lf
from terminal_pareto import representation_criteria as rc


class LayerFeatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ctx = lf.load_context()
        ranking = rc.load_ranking_module()
        cls.terminals, cls.descendants = ranking.lineage_terminals()
        cls.classes = ranking.terminal_cell_types(cls.terminals)

    def test_reproduction_gates_pass(self):
        self.assertEqual(lf.check_reproduction(self.ctx)["failures"], [])

    def test_mixtures_follow_complete_lineage(self):
        m = self.ctx.mixtures
        np.testing.assert_allclose(m.sum(axis=1), 1, atol=1e-12)
        biological = [c for c in m.index if c in set(self.terminals)]
        self.assertTrue((m.loc[biological].max(axis=1) == 1).all())
        # A leaf cut off by the tracking limit carries the mixture of the cells it will produce.
        cutoff = [c for c in self.ctx.strata["1"].cells if c not in set(self.terminals)]
        self.assertEqual(len(cutoff), 216)
        cell = cutoff[0]
        counts = pd.Series([self.classes[t] for t in self.descendants[cell]]).value_counts()
        np.testing.assert_allclose(m.loc[cell, counts.index], counts / counts.sum(), atol=1e-12)
        self.assertEqual(self.ctx.descendants[cell], len(self.descendants[cell]))

    def test_exclude_death_renormalizes_over_living_classes(self):
        cells = self.ctx.strata["2"].cells
        living, _ = lf.schema_d_mixtures(cells, include_death=False)
        self.assertNotIn(lf.DEATH, living.columns)
        np.testing.assert_allclose(living.sum(axis=1), 1, atol=1e-12)
        full = self.ctx.mixtures.loc[living.index, living.columns]
        np.testing.assert_allclose(living, full.div(full.sum(axis=1), axis=0), atol=1e-12)

    def test_strata_partition_rounds_and_pairs_are_valid(self):
        ctx = self.ctx
        cells = [c for s in ctx.strata.values() for c in s.cells]
        self.assertEqual(len(cells), len(set(cells)))
        self.assertEqual(set(cells), set(ctx.rounds))
        self.assertEqual([len(s.cells) for s in ctx.strata.values()], [491, 223, 123, 66, 71])
        parent = dict(zip(ctx.pooled.names, (ctx.pooled.names[p] if p >= 0 else None for p in ctx.pooled.parents)))
        for s in ctx.strata.values():
            self.assertTrue(all(ctx.rounds[c] in s.rounds for c in s.cells))
            for i, j in zip(s.a, s.b):
                self.assertNotIn(s.cells[i], ctx.lineage.ancestors[s.cells[j]])
                self.assertNotIn(s.cells[j], ctx.lineage.ancestors[s.cells[i]])
            for k in s.sisters:
                self.assertEqual(parent[s.cells[s.a[k]]], parent[s.cells[s.b[k]]])
            # Path lengths entering the criterion depend on labels only.
            for value, index, ranks in s.paths:
                self.assertGreaterEqual(len(index), lf.MIN_PATH_PAIRS)
                self.assertTrue((s.path[index] == value).all())
                self.assertGreater(np.ptp(s.bc[index]), 0)
            for index in s.memory.values():
                self.assertTrue((s.majority[s.a[index]] == s.majority[s.b[index]]).all())
            self.assertTrue((s.bc[s.memory["near-identical"]] <= lf.MEMORY_BC_MAX).all())
        # Only the pooled stratum can hold ancestor-descendant pairs; they are removed.
        pooled = ctx.strata["5-8"]
        self.assertEqual(len(pooled.a), len(pooled.cells) * (len(pooled.cells) - 1) // 2 - 143)
        for name in ("1", "2", "3", "4"):
            s = ctx.strata[name]
            self.assertEqual(len(s.a), len(s.cells) * (len(s.cells) - 1) // 2)

    def test_separation_on_constructed_distances(self):
        s = self.ctx.strata["3"]
        row, per_path = lf.separation(s.bc.copy(), s)
        self.assertAlmostEqual(row["separation"], 1.0)
        self.assertAlmostEqual(row["separation_pair_weighted"], 1.0)
        self.assertAlmostEqual(row["sisters"], 1.0)
        self.assertAlmostEqual(row["all_pairs"], 1.0)
        row, per_path = lf.separation(np.ones(len(s.a)), s)
        self.assertEqual(row["separation"], 0.0)
        # Each path length counts once, however many pairs it holds.
        d = np.zeros(len(s.a))
        value, index, ranks = s.paths[-1]
        d[index] = ranks
        row, per_path = lf.separation(d, s)
        self.assertAlmostEqual(row["separation"], 1 / len(s.paths))
        self.assertAlmostEqual(row["separation_pair_weighted"], len(index) / sum(n for _, n, _ in per_path))
        self.assertEqual(sum(n for _, n, _ in per_path), sum(len(i) for _, i, _ in s.paths))

    def test_majority_shuffle_stays_within_class(self):
        s = self.ctx.strata["2"]
        rng = np.random.default_rng(0)
        for _ in range(20):
            perm = lf.majority_permutation(s, rng)
            np.testing.assert_array_equal(np.sort(perm), np.arange(len(s.cells)))
            np.testing.assert_array_equal(s.majority[perm], s.majority)
        D = lf.cosine_matrix(lf.matrix(self.ctx, s, self.ctx.candidates["top-20"]))
        perm = lf.majority_permutation(s, rng)
        expected = D[np.ix_(perm, perm)][s.a, s.b]
        np.testing.assert_array_equal(lf.pair_distances(D, s, perm), expected)

    def test_run_write_verify_and_tamper_guards(self):
        small = dict(n_boot=20, n_perm=20, n_random=6)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run"
            lf.write_run(out, n_jobs=2, **small)
            with self.assertRaises(FileExistsError):
                lf.write_run(out, n_jobs=2, **small)
            self.assertEqual(lf.validate_run(out, n_jobs=2)["tables"], len(lf.TABLES))
            analysis = out / "analysis"
            criteria = pd.read_csv(analysis / "criteria.csv")
            criteria.loc[0, "separation"] += .1
            criteria.to_csv(analysis / "criteria.csv", index=False)
            with self.assertRaises(ValueError):
                lf.validate_run(out, n_jobs=2)
            # Even a rehashed incorrect value must fail the replay.
            record = json.loads((analysis / "provenance.json").read_text())
            record["files"]["criteria.csv"] = lf.file_hash(analysis / "criteria.csv")
            (analysis / "provenance.json").write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "Replay mismatch criteria"):
                lf.validate_run(out, n_jobs=2)


if __name__ == "__main__":
    unittest.main()
