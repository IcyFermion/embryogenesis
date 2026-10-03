"""Pooled full-tree Figure 7/S4: cache-only build against the approved working layout."""

from __future__ import annotations

from pathlib import Path
import shutil
import tempfile
import unittest

from publication import provenance, registry
from publication.parity import compare
from publication.tests.support import no_experiments

APPROVED = provenance.ROOT / "full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_20260927/publication"
HAS_TEX = shutil.which("tectonic") is not None or shutil.which("pdflatex") is not None


@unittest.skipUnless(HAS_TEX, "tectonic or pdflatex required")
class CacheOnlyPooledBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="publication-pooled-"))
        with no_experiments():
            cls.record = registry.build("full-tree-pooled", cls.tmp / "build")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_inventory_matches_approved_working_build(self):
        approved = {p.name for p in APPROVED.iterdir() if p.suffix in (".pdf", ".png", ".tex")}
        self.assertEqual(set(self.record["files"]), approved)
        self.assertEqual(self.record["figure_numbers"],
                         {"fig7_ce_full_tree_layerwise": "7", "figs_ce_full_tree_heuristics": "S4"})

    def test_panels_match_approved_layout(self):
        names = [n for n in self.record["files"] if n.endswith(".png")]
        for name, result in compare(self.tmp / "build", APPROVED, names).items():
            self.assertTrue(result["identical"], name)

    def test_only_the_reviewed_figure7_warning_is_tolerated(self):
        compiled = {stem: figure["compiled"] for stem, figure in self.record["figures"].items()}
        self.assertEqual(compiled["figs_ce_full_tree_heuristics"]["warnings"], [])
        self.assertTrue(all("Float too large for page by 8.2687pt" in w
                            for w in compiled["fig7_ce_full_tree_layerwise"]["warnings"]))
        self.assertEqual({c["pages"] for c in compiled.values()}, {1})


if __name__ == "__main__":
    unittest.main()
