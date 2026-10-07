"""Checks for the expression-feature calibration survey (Track A1)."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from terminal_pareto import feature_sensitivity as fs


class StateMatrixTests(unittest.TestCase):
    def test_cosine_parent_rows_and_child_columns(self):
        frame = pd.DataFrame({"a": [1., 0., 1.], "b": [0., 1., 1.]}, index=["p", "c1", "c2"])
        state = fs.state_matrix(frame, ["a", "b"], [("c1", "p"), ("c2", "p")])
        np.testing.assert_allclose(state, [[1., 1 - 1 / np.sqrt(2)], [1., 1 - 1 / np.sqrt(2)]])

    def test_zero_vector_is_rejected(self):
        frame = pd.DataFrame({"a": [1., 0.], "b": [1., 0.]}, index=["p", "c"])
        with self.assertRaises(fs.ZeroVectorError):
            fs.state_matrix(frame, ["a", "b"], [("c", "p")])

    def test_variability_pool_is_top_fraction_by_log1p_variance(self):
        cells = [f"x{i}" for i in range(6)]
        frame = pd.DataFrame({"flat": [1.] * 6, "low": [0, 0, 0, 0, 0, 1.], "high": [0, 9, 0, 9, 0, 9.],
                              "mid": [0, 2, 0, 2, 0, 2.]}, index=cells)
        table = fs.variability_scores({"raw": frame}, ["flat", "low", "high", "mid"], cells)
        self.assertEqual(table.feature.tolist(), ["high", "mid", "low", "flat"])
        self.assertEqual(table.loc[table.in_hv_pool, "feature"].tolist(), ["high", "mid"])
        with self.assertRaises(ValueError):
            fs.variability_scores({"raw": -frame}, ["low"], cells)


    def test_summary_counts_degenerate_sets_without_summarizing_them(self):
        sets = pd.DataFrame(dict(set_id=["t", "a1", "a2", "a3", "h1"], kind=["top", "random_all", "random_all",
                                 "random_all", "random_hv"], k=5))
        metrics = pd.DataFrame(dict(set_id=sets.set_id, cohort="primary", geometry="raw3d", config="ce_protein",
                                    status=["ok", "ok", "ok", "degenerate_state", "ok"]))
        for metric in fs.SUMMARY_METRICS:
            metrics[metric] = [0.5, 0.2, 0.5, np.nan, 0.9]
        row = fs.summarize(sets, metrics).iloc[0]
        self.assertEqual((row.random_all_n, row.random_all_degenerate), (2, 1))
        self.assertEqual(row.random_all_percentile, 75.0)  # one below, one tied
        self.assertEqual(row.random_hv_percentile, 0.0)


class SurveyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = fs.load_inputs()

    def test_configurations_exclude_c_briggsae(self):
        self.assertFalse(any("cb" in config for _, _, config in fs.CONFIGS))
        self.assertEqual(set(self.inputs.frames), {"ce_protein", "ce_rna"})

    def test_rna_matched_cohort_restricts_primary_travel(self):
        primary, matched = self.inputs.cohorts["primary"], self.inputs.cohorts["rna_matched"]
        self.assertEqual(len(matched.terms), 217)
        rna = self.inputs.frames["ce_rna"]
        self.assertEqual(matched.terms, [e for e in primary.terms if all(n in rna.index for n in e)])
        index = [primary.terms.index(edge) for edge in matched.terms]
        for geometry in ("raw3d", "xy"):
            np.testing.assert_array_equal(matched.travel[geometry], primary.travel[geometry][np.ix_(index, index)])
            # Fixed 275-edge denominators: the pooled natural total over all 275 edges is one.
            self.assertAlmostEqual(np.trace(primary.travel[geometry]), 1.0)

    def test_published_top20_reproduces_both_caches(self):
        report, _ = fs.check_published(self.inputs)
        self.assertEqual(report["failures"], [])
        self.assertEqual(len(report["state_matrices"]), len(fs.REPRODUCED))

    def test_pools_and_rankings(self):
        for modality, size in (("protein", 210), ("rna", 237)):
            ranking = self.inputs.rankings[modality]
            self.assertEqual(ranking[:20], self.inputs.published_features[modality])
            pool = self.inputs.hv_scores[modality]
            self.assertEqual(len(pool), size)
            self.assertEqual(int(pool.in_hv_pool.sum()), round(size * fs.HV_FRACTION))

    def test_draws_are_seeded_valid_and_recorded(self):
        first = fs.draw_feature_sets(self.inputs, ks=(5,), draws=4)
        second = fs.draw_feature_sets(self.inputs, ks=(5,), draws=4)
        pd.testing.assert_frame_equal(first, second)
        self.assertEqual(len(first), 2 * (1 + 2 * 4))
        for row in first.itertuples():
            features = row.features.split(";")
            self.assertEqual(len(set(features)), 5)
            self.assertTrue(fs.feature_set_is_valid(self.inputs, row.modality, features))
            if row.kind == "top":
                self.assertEqual(features, self.inputs.rankings[row.modality][:5])
            if row.kind == "random_hv":
                self.assertEqual(row.n_in_hv_pool, 5)
        other = fs.draw_feature_sets(self.inputs, ks=(5,), draws=4, seed=fs.SEED + 1)
        self.assertFalse(first.features.equals(other.features))

    def test_small_run_writes_once_and_replays(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run"
            record = fs.write_run(out, n_jobs=2, ks=(10,), draws=2)
            self.assertEqual(record["n_evaluations"], 3 * 5 + 2 * 5)
            self.assertEqual(record["n_failed"], 0)
            report = fs.validate_run(out, replay=3)
            self.assertEqual(report["evaluations"], 25)
            with self.assertRaises(FileExistsError):
                fs.write_run(out, n_jobs=1, ks=(10,), draws=2)


if __name__ == "__main__":
    unittest.main()
