"""Terminal primary family: cache-only build, production inventory and notation propagation."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from publication import notation, provenance, registry
from publication.adapters import terminal as adapter
from publication.captions.terminal import CAPTIONED_FIGURES
from publication.figures import terminal_canonical
from publication.tables import subtree_statistics_tex
from publication.tests.support import no_experiments

PRODUCTION = provenance.ROOT / "terminal_pareto/output/publication"
HAS_TEX = shutil.which("tectonic") is not None or shutil.which("pdflatex") is not None
MARKER = r"\Delta_{\mathrm{test}}"


@unittest.skipUnless(HAS_TEX, "tectonic or pdflatex required")
class CacheOnlyTerminalBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="publication-terminal-"))
        with no_experiments():
            cls.record = registry.build("terminal-primary", cls.tmp / "build")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_every_wrapper_compiles_to_one_labelled_page(self):
        self.assertEqual(set(self.record["figures"]), set(CAPTIONED_FIGURES))
        for stem, figure in self.record["figures"].items():
            self.assertEqual(figure["compiled"]["pages"], 1, stem)
            self.assertEqual(figure["compiled"]["warnings"], [], stem)

    def test_rendered_inventory_matches_production_figures(self):
        released = set(json.loads((PRODUCTION / "release_manifest.json").read_text())["files"])
        figure_assets = {name for name in released if Path(name).suffix in (".pdf", ".png", ".svg", ".tex")
                         and not name.startswith(("fig8_", "fig9_"))}
        self.assertEqual(set(self.record["files"]), figure_assets)
        self.assertTrue(provenance.verify_build(self.tmp / "build")["ready"])

    def test_numerical_inputs_are_recorded_and_unchanged(self):
        run = Path(self.record["run"])
        for name, digest in self.record["input_files"].items():
            self.assertEqual(provenance.sha256(run / name), digest, name)


class TerminalNotation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with no_experiments():
            cls.data = adapter.primary()

    def test_override_reaches_captions_table_and_canonical_artists(self):
        with notation.override(null_front_distance=dict(symbol=MARKER, plain="Delta_test")):
            for stem in ("fig1_ce_endpoint_normalization_amendment", "fig5_ce_canonical_summary"):
                self.assertIn(f"${MARKER}$", CAPTIONED_FIGURES[stem][1](), stem)
            table = subtree_statistics_tex(self.data.canonical, self.data.subtree_summary, self.data.min_cells)
            self.assertIn(rf"\({MARKER}\)", table)
            fig, ax = plt.subplots()
            try:
                terminal_canonical.plot_definition_dlp(ax)
                self.assertTrue(any(MARKER in text.get_text() for text in ax.texts))
            finally:
                plt.close(fig)
        self.assertIn(r"\(d_{CP}\)", subtree_statistics_tex(self.data.canonical, self.data.subtree_summary, 12))
        self.assertNotIn("d_{NP}", CAPTIONED_FIGURES["fig5_ce_canonical_summary"][1]())

    def test_table_rows_cover_every_endpoint_valid_subtree(self):
        table = subtree_statistics_tex(self.data.canonical, self.data.subtree_summary, self.data.min_cells)
        valid = self.data.canonical[self.data.canonical["endpoint_ok"].astype(bool)]
        rows = [line for line in table.splitlines() if line.split(" & ")[0] in set(valid.subtree)]
        self.assertEqual(sorted(line.split(" & ")[0] for line in rows), sorted(valid.subtree))


if __name__ == "__main__":
    unittest.main()
