"""Numbering, caption and archive boundaries for Figures 8 and 9."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from terminal_pareto.cross_species_analysis import file_hash
from terminal_pareto.cross_species_publication import archive_previous, captions, write_wrappers
from terminal_pareto.fig_terminal_cross_species import (
    PRIMARY, canonical_point, canonical_records, canonical_summary, comparison, overlay_panel,
)


def figure_fixture():
    fronts, metrics, clouds = [], [], []
    for geometry in ("raw3d", "xy"):
        for config in PRIMARY:
            fronts.extend(dict(geometry=geometry, config=config, cohort="base", sweep_index=i,
                               D1=x, D2=y, edge_retention=er) for i, (x, y, er) in enumerate(
                                   ((1., 0., .2), (.3, .3, .5), (0., 1., .8))))
            metrics.append(dict(geometry=geometry, config=config, cohort="base", nearest_index=1,
                                maximum_index=2, natural_D1=.32, natural_D2=.33,
                                null_D1=.5, null_D2=.5, u_L=.5,
                                d_LP=np.hypot(.02, .03), d_NP=np.hypot(.2, .2)))
            clouds.extend(dict(geometry=geometry, config=config, family=family, D1=x, D2=y)
                          for family in ("first_cousin", "second_cousin", "third_cousin", "full_random")
                          for x, y in ((.5, .6), (.52, .61), (.54, .63)))
    return tuple(pd.DataFrame(rows) for rows in (fronts, metrics, clouds))


class PublicationTests(unittest.TestCase):
    def test_wrappers_match_publication_style_and_numbering(self):
        with tempfile.TemporaryDirectory() as folder:
            wrappers = write_wrappers(Path(folder), dict(cohort_size=187, cb_reference_size=196))
            self.assertEqual(len(wrappers), 2)
            for number, path in zip((8, 9), wrappers):
                text = path.read_text()
                self.assertIn(r"\documentclass[10pt]{article}", text)
                self.assertIn(rf"\setcounter{{figure}}{{{number - 1}}}", text)
                self.assertIn(f"{path.stem}_panel.pdf", text)
                self.assertIn("187", text)
                self.assertIn("traditional embryo-tracking techniques", text)
                self.assertIn("provisional", text)
                self.assertNotIn("axial scale unverified", text)
            self.assertIn("196-edge reference", captions(dict(cohort_size=187, cb_reference_size=196))[8])
            caption9 = captions(dict(cohort_size=187, cb_reference_size=196))[9]
            self.assertIn(r"\textbf{(C)}", caption9)
            self.assertIn(r"d_{LP}", caption9)
            self.assertIn(r"d_{NP}", caption9)

    def test_canonical_metrics_use_closest_attained_point_not_maximum_retention(self):
        fronts, metrics, _ = figure_fixture()
        curve = fronts[(fronts.geometry == "raw3d") & (fronts.config == PRIMARY[0])]
        metric = metrics.iloc[0]
        closest, values = canonical_point(curve, metric)
        self.assertEqual(closest.sweep_index, 1)
        self.assertNotEqual(closest.sweep_index, metric.maximum_index)
        self.assertAlmostEqual(values["u_L"], .5)
        self.assertAlmostEqual(values["d_LP"], np.hypot(.02, .03))
        self.assertAlmostEqual(values["d_NP"], np.hypot(.2, .2))
        metric = metric.copy()
        metric.d_NP = 99.
        with self.assertRaises(AssertionError):
            canonical_point(curve, metric)

    def test_metric_panel_displays_all_six_rows_and_eighteen_values(self):
        fronts, metrics, _ = figure_fixture()
        records = canonical_records(fronts, metrics)
        self.assertEqual([(r["geometry"], r["config"]) for r in records],
                         [(g, c) for g in ("raw3d", "xy") for c in PRIMARY])
        fig = plt.figure()
        try:
            axes = canonical_summary(fig, fronts, metrics)
            for ax, field in zip(axes, ("u_L", "d_LP", "d_NP")):
                points = next(c for c in ax.collections if c.get_gid() == f"canonical:{field}")
                np.testing.assert_allclose(points.get_offsets(),
                                           np.column_stack(([r[field] for r in records], [1, 2, 3, 6, 7, 8])))
        finally:
            plt.close(fig)

    def test_figure8_has_no_natural_to_maximum_connector(self):
        fronts, metrics, clouds = figure_fixture()
        captured = []
        with patch("terminal_pareto.fig_terminal_cross_species.save", side_effect=lambda fig, *_: captured.append(fig)):
            comparison(fronts, metrics, clouds, None, Path("unused"), 187)
        try:
            self.assertEqual(len(captured[0].axes), 7)  # Six panels and the color bar.
            for ax in captured[0].axes[:6]:
                self.assertEqual(len(ax.lines), 3)  # Two zero guides and the front only.
        finally:
            plt.close(captured[0])

    def test_overlay_distance_connectors_end_at_closest_point(self):
        fronts, metrics, _ = figure_fixture()
        fig, ax = plt.subplots()
        try:
            overlay_panel(ax, fronts, metrics, "raw3d", legend=False)
            for region in (ax, ax.child_axes[0]):
                connections = [line for line in region.lines if (line.get_gid() or "").startswith("natural_to_closest:")]
                self.assertEqual(len(connections), 3)
                for line in connections:
                    np.testing.assert_allclose(line.get_xdata(), [.32, .3])
                    np.testing.assert_allclose(line.get_ydata(), [.33, .3])
        finally:
            plt.close(fig)

    def test_rebuild_preserves_previous_publication_with_hashes(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)
            source = run / "publication"
            source.mkdir()
            figure = source / "old.pdf"
            figure.write_bytes(b"previous figure fixture")
            expected = file_hash(figure)
            archive = archive_previous(run)
            self.assertEqual(file_hash(archive / "publication/old.pdf"), expected)
            self.assertEqual(json.loads((archive / "manifest.json").read_text())["files"]["old.pdf"], expected)
            self.assertEqual(file_hash(figure), expected)

    def test_empty_run_has_no_spurious_archive(self):
        with tempfile.TemporaryDirectory() as folder:
            self.assertIsNone(archive_previous(Path(folder)))


if __name__ == "__main__":
    unittest.main()
