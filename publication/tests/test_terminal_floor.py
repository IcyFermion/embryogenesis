"""Subtree floor of Figures 4-5, S1 and Table 1 (Track D6): resolution rule, caption note and counts."""

from __future__ import annotations

import copy
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import pandas as pd

from publication.adapters import terminal as adapter
from publication.captions import terminal as captions
from publication.figures import terminal_canonical
from publication.tests.support import no_experiments
from terminal_pareto import subtree_symmetry as ss

RUN = ss.ROOT / "terminal_pareto/output/runs" / adapter.PRIMARY_PROFILE / adapter.PRIMARY_RUN
ANALYSIS = RUN / "analysis"


def load(analysis=ANALYSIS):
    with no_experiments():
        return adapter._track_d(RUN, analysis)


class Floor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables, cls.inputs = load()
        cls.canonical = pd.read_csv(ANALYSIS / "ce_subtree_canonical_metrics.csv")

    def test_the_floor_follows_the_resolution_rule(self):
        resolution = adapter._resolution(self.tables, self.canonical)
        self.assertEqual(resolution, dict(floor=24, kept=18, resolved=23, largest_unresolved=("ABala", 22)))

    def test_a_floor_that_drops_a_resolved_subtree_is_refused(self):
        with self.assertRaisesRegex(ValueError, "resolution rule"):
            adapter._resolution(self.tables, self.canonical, floor=26)  # would drop MSp (25 cells)

    def test_an_unresolved_subtree_above_the_floor_is_refused(self):
        tables = copy.deepcopy(self.tables)
        metrics = tables["compare", "subtree_metrics"]
        metrics.loc[metrics.subtree == "MSp", "null_support"] = 11
        with self.assertRaisesRegex(ValueError, "resolution rule"):
            adapter._resolution(tables, self.canonical)

    def test_inputs_are_recorded_relative_to_the_primary_run(self):
        self.assertEqual(len(self.inputs), 6)
        for name, digest in self.inputs.items():
            self.assertTrue(name.startswith("../../subtree_symmetry_v1/"), name)
            self.assertEqual(adapter.sha256(RUN / name), digest)

    def test_a_changed_table_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for run_id in adapter.TRACK_D_RUNS.values():
                shutil.copytree(ss.run_path(run_id) / "analysis", root / run_id / "analysis",
                                ignore=shutil.ignore_patterns("controls.csv", "pools.csv"))
            path = root / adapter.TRACK_D_RUNS["single"] / "analysis/single.csv"
            path.write_text(path.read_text() + "\n")
            with mock.patch.object(ss, "RUN_ROOT", root), self.assertRaisesRegex(ValueError, "provenance record"):
                load()

    def test_caption_counts_match_the_cached_subtrees(self):
        summary = pd.read_csv(ANALYSIS / f"subtree_summary_min{adapter.PRIMARY_MIN_CELLS}.csv")
        shown = summary[summary.n >= adapter.SUBTREE_FLOOR]
        self.assertEqual((captions.SUBTREES, captions.FLOOR, captions.ON_FRONT),
                         (len(shown), adapter.SUBTREE_FLOOR, int(shown.natural_on_front.sum())))
        self.assertEqual(int((self.canonical.n >= adapter.SUBTREE_FLOOR).sum()), captions.SUBTREES)


class Note(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tables, _ = load()
        canonical = pd.read_csv(ANALYSIS / "ce_subtree_canonical_metrics.csv")
        cls.shown = canonical[canonical.n >= adapter.SUBTREE_FLOOR]
        cls.evidence = adapter._note_evidence(tables, cls.shown)

    def test_note_holds_and_fails_when_the_evidence_changes(self):
        terminal_canonical.check_notes(self.evidence)
        changed = copy.deepcopy(self.evidence)
        changed["control_on_front"] = 3
        with self.assertRaisesRegex(ValueError, "0 of 1000"):
            terminal_canonical.check_notes(changed)
        changed = copy.deepcopy(self.evidence)
        changed["alone_on_front"]["ABplppaaap"] = True
        with self.assertRaisesRegex(ValueError, "only ABplpppapp"):
            terminal_canonical.check_notes(changed)

    def test_rows_render_without_computation(self):
        with tempfile.TemporaryDirectory() as tmp, no_experiments():
            terminal_canonical.plot_primary_rows(self.shown, tmp, notes=self.evidence)
            self.assertTrue((Path(tmp) / "fig5B_ce_canonical_all_subtrees.pdf").is_file())
            self.assertFalse((Path(tmp) / "fig5C_ce_mirrored_pairs.pdf").exists())

    def test_caption_carries_the_note_and_no_pair_panel(self):
        text = captions.fig5_ce_canonical_summary()
        for piece in (r"$\dagger$", "ABplpppapp", "0 of 1000", "fig5B_ce_canonical_all_subtrees.pdf"):
            self.assertIn(piece, text)
        for piece in ("fig5C", "mirrored", r"$\ddagger$", r"$\S$"):
            self.assertNotIn(piece, text)


if __name__ == "__main__":
    unittest.main()
