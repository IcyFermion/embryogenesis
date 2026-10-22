"""Mirrored-pair layer of Figure 5: Track D loader, caption-note checks, pair numbers and panel C."""

from __future__ import annotations

import copy
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from publication import notation
from publication.adapters import terminal as adapter
from publication.captions.terminal import CAPTIONED_FIGURES
from publication.figures import terminal_symmetry as sym
from publication.tests.support import no_experiments
from terminal_pareto import subtree_symmetry as ss

RUN = ss.ROOT / "terminal_pareto/output/runs" / adapter.PRIMARY_PROFILE / adapter.PRIMARY_RUN
ANALYSIS = RUN / "analysis"
MARKER = r"\Delta_{\mathrm{test}}"


def load(analysis=ANALYSIS):
    canonical = pd.read_csv(analysis / "ce_subtree_canonical_metrics.csv")
    with no_experiments():
        return adapter._symmetry(RUN, analysis, canonical)


class Loader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.symmetry, cls.inputs = load()

    def test_nine_pairs_with_both_partners(self):
        pairs = self.symmetry["pairs"]
        self.assertEqual(len(pairs), 9)
        self.assertEqual(pairs.tier.value_counts().to_dict(), dict(strict=5, priority=2, exploratory=2))
        self.assertEqual(set(pairs.left), {"ABalp", "ABalpa", "ABplp", "ABplpa", "ABplpp", "ABplppp", "MSa", "MSaa", "MSap"})

    def test_inputs_are_recorded_relative_to_the_primary_run(self):
        self.assertEqual(len(self.inputs), 12)
        for name, digest in self.inputs.items():
            self.assertTrue(name.startswith("../../subtree_symmetry_v1/"), name)
            self.assertEqual(adapter.sha256(RUN / name), digest)

    def test_a_changed_table_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for run_id in adapter.SYMMETRY_RUNS.values():
                shutil.copytree(ss.run_path(run_id) / "analysis", root / run_id / "analysis",
                                ignore=shutil.ignore_patterns("controls.csv", "pools.csv"))
            path = root / adapter.SYMMETRY_RUNS["classification"] / "analysis/classification.csv"
            path.write_text(path.read_text() + "\n")
            with mock.patch.object(ss, "RUN_ROOT", root), self.assertRaisesRegex(ValueError, "provenance record"):
                load()

    def test_notes_fail_when_the_evidence_changes(self):
        sym.check_notes(self.symmetry)
        changed = copy.deepcopy(self.symmetry)
        changed["evidence"]["abplp_control_on_front"] = 3
        with self.assertRaisesRegex(ValueError, "0 of 1000"):
            sym.check_notes(changed)

    def test_pair_numbers_follow_the_display_order(self):
        tags = sym.row_tags(self.symmetry["pairs"])
        self.assertEqual(len(tags), 18)
        self.assertEqual(tags["ABplpa"], dict(number=2, side="left", partner="ABprpa", mark="†", metrics=("u",)))
        self.assertEqual((tags["ABprpa"]["number"], tags["ABprpa"]["side"]), (2, "right"))
        self.assertEqual((tags["MSpa"]["number"], tags["MSpa"]["mark"]), (6, "‡"))
        self.assertEqual((tags["ABarap"]["number"], tags["ABarap"]["mark"]), (9, None))

    def test_panel_c_renders_without_computation_and_carries_notation(self):
        with tempfile.TemporaryDirectory() as tmp, no_experiments(), \
                notation.override(natural_front_distance_sd=dict(symbol=MARKER)):
            fig = sym.plot_mirrored_pairs(self.symmetry, tmp)
            try:
                titles = [ax.get_title() for ax in fig.axes]
                self.assertTrue(any(MARKER in t for t in titles))
                self.assertTrue((Path(tmp) / "fig5C_ce_mirrored_pairs.pdf").is_file())
            finally:
                plt.close(fig)


class Caption(unittest.TestCase):
    def test_caption_carries_panel_c_symbols_and_notes(self):
        text = CAPTIONED_FIGURES["fig5_ce_canonical_summary"][1]()
        for piece in (r"$\tilde{d}_{LP}$", "$z_T$", "$z_S$", "fig5C_ce_mirrored_pairs.pdf", r"$\dagger$",
                      r"$\ddagger$", r"$\S$", "ABplpppapp", "0 of 1000"):
            self.assertIn(piece, text)
        with notation.override(natural_state_position=dict(symbol=MARKER)):
            self.assertIn(f"${MARKER}$", CAPTIONED_FIGURES["fig5_ce_canonical_summary"][1]())


if __name__ == "__main__":
    unittest.main()
