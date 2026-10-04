"""Figures 9/11 vertical slice: cache-only builds, artists, notation and provenance.

Uses the published numerical runs read-only. Run from the repository root:

    python -m unittest discover -s publication/tests -t .
"""

from __future__ import annotations

import json
import io
from contextlib import ExitStack, redirect_stdout
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from publication import notation, provenance, registry, style
from publication.adapters import full_tree as full_tree_adapter, terminal as terminal_adapter
from publication.canonical_panels import canonical_records
from publication.figures.cross_species import compose_overlay
from publication.tests.support import no_experiments

HAS_TEX = shutil.which("tectonic") is not None or shutil.which("pdflatex") is not None
CROSS_SPECIES = ("terminal-cross-species", "full-tree-cross-species")
FAMILY_TITLES = {key: next(spec for spec in registry.FAMILIES[key].figures if spec.stem.endswith("_overlays"))
                 for key in CROSS_SPECIES}


class CrossSpeciesSlice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        style.configure()
        with no_experiments():
            cls.inputs = {terminal_adapter.CROSS_SPECIES_FAMILY: terminal_adapter.cross_species(),
                          full_tree_adapter.CROSS_SPECIES_FAMILY: full_tree_adapter.cross_species()}

    def compose(self, family):
        spec = FAMILY_TITLES[family]
        return compose_overlay(self.inputs[family], **spec.render.keywords)

    def test_scopes_are_explicit(self):
        terminal = self.inputs["terminal-cross-species"]
        partial = self.inputs["full-tree-cross-species"]
        self.assertEqual((terminal.edges, partial.edges), (187, 454))
        self.assertEqual(partial.caption["nodes"], 485)
        self.assertEqual(partial.caption["roots"], 31)
        for data in self.inputs.values():
            self.assertEqual(len(data.metrics), 6)

    def test_canonical_metrics_match_cache_and_use_closest_point(self):
        for data in self.inputs.values():
            cached = pd.read_csv(data.run / "analysis/metrics.csv")
            for record in canonical_records(data):
                row = cached[(cached.geometry == record["geometry"]) & (cached.config == record["config"])
                             & (cached.cohort == "base")].iloc[0]
                self.assertAlmostEqual(record["canonical_position"], row.u_L, places=11)
                self.assertAlmostEqual(record["natural_front_distance"], row.d_LP, places=11)
                self.assertAlmostEqual(record["null_front_distance"], row.d_NP, places=11)
                self.assertEqual(record["nearest_index"], int(row.nearest_index))

    def test_artists_show_every_value_and_connector(self):
        for family, data in self.inputs.items():
            fig, _ = self.compose(family)
            try:
                records = canonical_records(data)
                artists = {child.get_gid(): child for ax in fig.axes for child in ax.get_children() if child.get_gid()}
                for key in ("canonical_position", "natural_front_distance", "null_front_distance"):
                    offsets = artists[f"canonical:{key}"].get_offsets()
                    np.testing.assert_allclose(offsets[:, 0], [r[key] for r in records], rtol=0, atol=1e-12)
                every_axis = [a for ax in fig.axes for a in (ax, *ax.child_axes)]
                connectors = [child for ax in every_axis for child in ax.get_lines()
                              if (child.get_gid() or "").startswith("natural_to_closest:")]
                self.assertEqual(len(connectors), 2 * len(records))  # main axes and inset
                for record in records:
                    _, metric = data.select(record["geometry"], record["config"])
                    gid = f"natural_to_closest:{record['geometry']}:{record['config']}"
                    for line in (c for c in connectors if c.get_gid() == gid):
                        np.testing.assert_allclose(np.column_stack(line.get_data()),
                                                   [[metric.natural_D1, metric.natural_D2],
                                                    [record["closest_D1"], record["closest_D2"]]], atol=1e-12)
            finally:
                plt.close(fig)

    def test_notation_override_propagates_to_both_families(self):
        marker = r"\Delta_{\mathrm{test}}"
        with notation.override(null_front_distance=dict(symbol=marker, plain="Delta_test")):
            for family, data in self.inputs.items():
                fig, _ = self.compose(family)
                titles = [ax.get_title() for ax in fig.axes]
                plt.close(fig)
                self.assertIn(f"Null-mean distance, ${marker}$", titles)
                caption = FAMILY_TITLES[family].caption(data.caption)
                self.assertIn(rf"\({marker}\)", caption)
                if family == "terminal-cross-species":
                    self.assertIn(rf"\({marker}=\lVert", caption)
                self.assertNotIn("d_{CP}", caption)
        caption = FAMILY_TITLES["terminal-cross-species"].caption(self.inputs["terminal-cross-species"].caption)
        self.assertIn(r"\(d_{CP}\)", caption)
        self.assertNotIn("d_{NP}", caption)

    def test_override_does_not_touch_numerical_inputs(self):
        before = {f: dict(d.input_files) for f, d in self.inputs.items()}
        with notation.override(null_front_distance=dict(symbol="x")):
            for family, data in self.inputs.items():
                after = {name: provenance.sha256(data.run / name) for name in data.input_files}
                self.assertEqual(after, before[family])

    def test_shared_layout_imports_no_back_end(self):
        for module in ("artists", "canonical_panels", "figures/cross_species", "captions/cross_species",
                       "notation", "style", "assembly", "contracts"):
            text = (provenance.PACKAGE / f"{module}.py").read_text()
            self.assertNotIn("terminal_pareto", text.replace("``terminal_pareto/plot_style.py``", ""))
            self.assertNotIn("full_tree_pareto", text)


@unittest.skipUnless(HAS_TEX, "tectonic or pdflatex required")
class CacheOnlyBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="publication-test-"))
        cls.records = {}
        with no_experiments():
            for family in CROSS_SPECIES:
                cls.records[family] = registry.build(family, cls.tmp / family)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_build_inventory_and_verification(self):
        for family, record in self.records.items():
            specs = registry.FAMILIES[family].figures
            self.assertEqual(set(record["files"]), {name for spec in specs for name in spec.files})
            status = provenance.verify_build(self.tmp / family)
            self.assertTrue(status["integrity"] and status["ready"])
            for spec in specs:
                self.assertIn(rf"\label{{{spec.tex_label}}}", (self.tmp / family / f"{spec.stem}.tex").read_text())

    def test_changed_presentation_source_is_stale_not_blessed(self):
        family = "terminal-cross-species"
        live = provenance.presentation_sources()
        changed = dict(live, **{"publication/notation.py": "0" * 64})
        with mock.patch.object(provenance, "presentation_sources", return_value=changed):
            status = provenance.verify_build(self.tmp / family)
        self.assertTrue(status["integrity"])
        self.assertFalse(status["ready"])
        self.assertEqual(status["stale_sources"], ["publication/notation.py"])

    def test_notation_change_alone_marks_build_stale(self):
        with notation.override(null_front_distance=dict(symbol="x")):
            self.assertFalse(provenance.verify_build(self.tmp / "full-tree-cross-species")["ready"])

    def test_altered_artifact_or_input_fails_integrity(self):
        copy = self.tmp / "altered"
        shutil.copytree(self.tmp / "full-tree-cross-species", copy)
        try:
            manifest = json.loads((copy / provenance.MANIFEST).read_text())
            manifest["input_files"]["analysis/metrics.csv"] = "0" * 64
            (copy / provenance.MANIFEST).write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "Numerical inputs changed"):
                provenance.verify_build(copy)
            (copy / "fig11_full_tree_cross_species_overlays.tex").write_text("tampered")
            with self.assertRaisesRegex(ValueError, "Changed or missing build artifacts"):
                provenance.verify_build(copy, check_inputs=False)
        finally:
            shutil.rmtree(copy)

    def test_destination_guards(self):
        with self.assertRaisesRegex(ValueError, "inside back-end output"):
            registry.build("terminal-cross-species", provenance.ROOT / "terminal_pareto/output/publication/x")
        with self.assertRaisesRegex(ValueError, "new or empty"):
            registry.build("terminal-cross-species", self.tmp / "terminal-cross-species")


class RetiredBackEndDrawing(unittest.TestCase):
    """Back ends compute only; their retired render hooks point to the front end."""

    def test_render_only_hooks_redirect(self):
        from full_tree_pareto import fig_cross_species
        from terminal_pareto import fig_terminal_cross_species
        for module in (fig_terminal_cross_species, fig_cross_species):
            output = io.StringIO()
            with redirect_stdout(output), self.assertRaises(SystemExit) as stopped:
                module.render(Path("run"))
            self.assertEqual(stopped.exception.code, 0)
            self.assertIn("python -m publication build", output.getvalue())

    def test_successful_numerical_commands_exit_successfully(self):
        from terminal_pareto import cross_species_analysis as terminal
        from full_tree_pareto import cross_species_analysis as full_tree
        for module, operations in ((terminal, ("load_inputs", "write_analysis", "validate_run")),
                                   (full_tree, ("build",))):
            with self.subTest(module=module.__name__), ExitStack() as stack:
                stack.enter_context(mock.patch("sys.argv", ["analysis", "--run-id", "test-completed-run"]))
                calls = [stack.enter_context(mock.patch.object(module, name, return_value={}))
                         for name in operations]
                output = stack.enter_context(redirect_stdout(io.StringIO()))
                with self.assertRaises(SystemExit) as stopped:
                    module.main()
                self.assertEqual(stopped.exception.code, 0)
                for call in calls:
                    call.assert_called_once()
                self.assertIn("Numerical results are ready", output.getvalue())

    def test_numerical_failures_still_propagate(self):
        from terminal_pareto import cross_species_analysis as terminal
        from full_tree_pareto import cross_species_analysis as full_tree
        for module in (terminal, full_tree):
            with self.subTest(module=module.__name__), \
                    mock.patch("sys.argv", ["analysis", "--render-only"]), \
                    mock.patch.object(module, "validate_run", side_effect=ValueError("invalid cache")), \
                    self.assertRaisesRegex(ValueError, "invalid cache"):
                module.main()

    def test_plot_style_shim_reexports_shared_style(self):
        from terminal_pareto import plot_style as legacy_style
        for name in ("COLORS", "SPECIES_COLORS", "NULL_MODEL_COLORS", "SEMANTIC_COLORS", "configure"):
            self.assertIs(getattr(legacy_style, name), getattr(style, name))


if __name__ == "__main__":
    unittest.main()
