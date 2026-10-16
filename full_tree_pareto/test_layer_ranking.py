import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from full_tree_pareto import layer_features as lf
from full_tree_pareto import layer_ranking as lr


class LayerRankingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ctx, cls.pools, cls.halves, _, _ = lr.load(n_splits=3)

    def test_protein_score_is_univariate_e0_criterion(self):
        s = self.ctx.strata["3"]
        Z = lf.matrix(self.ctx, s, list(self.ctx.zscore.columns))
        scores = lr.protein_scores(Z, s)
        for j in (0, 17, 101, 209):
            d = np.abs(Z[s.a, j] - Z[s.b, j])
            self.assertAlmostEqual(scores[j], lf.separation(d, s)[0]["separation"], places=12)

    def test_restrict_keeps_within_half_pairs_and_e0_rules(self):
        s = self.ctx.strata["2"]
        keep = np.zeros(len(s.cells), dtype=bool)
        keep[::2] = True
        r = lr.restrict(s, keep)
        self.assertTrue(keep[r.a].all() and keep[r.b].all())
        self.assertEqual(len(r.a), int((keep[s.a] & keep[s.b]).sum()))
        for value, index, _ in r.paths:
            self.assertGreaterEqual(len(index), lf.MIN_PATH_PAIRS)
            self.assertTrue((r.path[index] == value).all())
        for name, index in r.memory.items():
            full = {(s.a[i], s.b[i]) for i in s.memory[name] if keep[s.a[i]] and keep[s.b[i]]}
            self.assertEqual({(r.a[i], r.b[i]) for i in index}, full)
        np.testing.assert_array_equal(r.path[r.sisters], 2)

    def test_splits_keep_families_and_balance_cells(self):
        for name, halves in self.halves.items():
            labels = lr.families(self.ctx, self.ctx.strata[name])
            largest = np.bincount(labels).max()
            for split in range(3):
                (_, _, train_a, _, test_a), (_, _, train_b, _, _) = [h for h in halves if h[0] == split]
                np.testing.assert_array_equal(train_a, ~train_b)
                for f in np.unique(labels):
                    self.assertEqual(len(set(train_a[labels == f])), 1)
                self.assertLessEqual(abs(int(train_a.sum()) - int((~train_a).sum())), largest)
                self.assertTrue(np.isin(test_a.a, np.flatnonzero(~train_a)).all())

    def test_top_k_breaks_ties_by_name(self):
        self.assertEqual(lr.top_k(np.array([1., 2., 2., 0.]), ["d", "c", "a", "b"], 3), ["a", "c", "d"])

    def test_gate_rule(self):
        def heldout(stratum, wins, pct):
            n = len(wins)
            return pd.DataFrame(dict(stratum=stratum, per_round=np.where(wins, 1., 0.), top20=.5, shared=.5,
                                     per_round_percentile=pct, top20_percentile=.5, shared_percentile=.5,
                                     win=wins, win_shared=wins))
        passing = heldout("2", [True] * 8 + [False] * 2, .95)
        table = pd.concat([passing, heldout("3", [True] * 8 + [False] * 2, .95),
                           heldout("4", [True] * 7 + [False] * 3, .99)])
        gate = lr.gate_table(table)
        self.assertEqual(gate.passes.tolist(), [True, True, False])
        self.assertFalse(lr.decision(gate)["proceed_to_e2"])
        table = pd.concat([passing, heldout("3", [True] * 10, .95), heldout("4", [True] * 10, .9)])
        self.assertTrue(lr.decision(lr.gate_table(table))["proceed_to_e2"])

    def test_run_write_verify_and_tamper_guards(self):
        small = dict(n_splits=2, n_random=4, n_perm=5)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run"
            record = lr.write_run(out, n_jobs=4, **small)
            self.assertEqual(record["settings"]["k"], lr.K)
            with self.assertRaises(FileExistsError):
                lr.write_run(out, n_jobs=4, **small)
            self.assertEqual(lr.validate_run(out, n_jobs=4)["tables"], len(lr.TABLES))
            analysis = out / "analysis"
            sets = pd.read_csv(analysis / "sets.csv")
            self.assertEqual(sets.features.str.split(";").str.len().tolist(), [lr.K] * len(sets))
            heldout = pd.read_csv(analysis / "heldout.csv")
            heldout.loc[0, "per_round"] += .1
            heldout.to_csv(analysis / "heldout.csv", index=False)
            with self.assertRaises(ValueError):
                lr.validate_run(out, n_jobs=4)
            provenance = json.loads((analysis / "provenance.json").read_text())
            provenance["files"]["heldout.csv"] = lf.file_hash(analysis / "heldout.csv")
            (analysis / "provenance.json").write_text(json.dumps(provenance))
            with self.assertRaisesRegex(ValueError, "Replay mismatch heldout"):
                lr.validate_run(out, n_jobs=4)


if __name__ == "__main__":
    unittest.main()
