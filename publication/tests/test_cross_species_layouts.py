"""Artist-level checks for Figures 8-11 on small fixtures (ported from the back-end test suites)."""

from __future__ import annotations

from pathlib import Path
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from publication.adapters.full_tree import CROSS_SPECIES_REFERENCES as FULL_TREE_REFERENCES
from publication.adapters.terminal import CROSS_SPECIES_REFERENCES as TERMINAL_REFERENCES
from publication.canonical_panels import canonical_point, canonical_records, canonical_summary
from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput
from publication.figures import cross_species as layouts
from publication.registry import FAMILIES

PRIMARY = ("ce_protein", "ce_rna", "cb_rna")


def fixture(references, family="fixture"):
    fronts, metrics, clouds = [], [], []
    for geometry in ("raw3d", "xy"):
        for config in PRIMARY:
            fronts.extend(dict(geometry=geometry, config=config, sweep_index=i, D1=x, D2=y, edge_retention=er)
                          for i, (x, y, er) in enumerate(((1., 0., .2), (.3, .3, .5), (0., 1., .8))))
            metrics.append(dict(geometry=geometry, config=config, nearest_index=1, maximum_index=2,
                                natural_D1=.32, natural_D2=.33, null_D1=.5, null_D2=.5, u_L=.5,
                                d_LP=np.hypot(.02, .03), d_NP=np.hypot(.2, .2)))
            clouds.extend(dict(geometry=geometry, config=config, family=name, D1=x, D2=y)
                          for name in (*references["main"], references["inset"])
                          for x, y in ((.5, .6), (.52, .61), (.54, .63)))
    return FrontComparisonInput(
        family=family, run=Path("unused"), analysis_id="fixture", edges=454, configs=PRIMARY,
        geometries=("raw3d", "xy"), fronts=pd.DataFrame(fronts),
        metrics=pd.DataFrame(metrics).rename(columns=CANONICAL_FIELDS), null_mean="analytic first-cousin",
        clouds=pd.DataFrame(clouds), references=references)


def options(family, index):
    return FAMILIES[family].figures[index].render.keywords


class CanonicalMetrics(unittest.TestCase):
    def test_closest_attained_point_not_maximum_retention(self):
        data = fixture(TERMINAL_REFERENCES)
        curve, metric = data.select("raw3d", PRIMARY[0])
        closest, values = canonical_point(curve, metric)
        self.assertEqual(closest.sweep_index, 1)
        self.assertNotEqual(closest.sweep_index, metric.maximum_index)
        self.assertAlmostEqual(values["canonical_position"], .5)
        self.assertAlmostEqual(values["natural_front_distance"], np.hypot(.02, .03))
        self.assertAlmostEqual(values["null_front_distance"], np.hypot(.2, .2))
        metric = metric.copy()
        metric.null_front_distance = 99.
        with self.assertRaises(AssertionError):
            canonical_point(curve, metric)

    def test_metric_panel_displays_all_six_rows_and_eighteen_values(self):
        data = fixture(TERMINAL_REFERENCES)
        records = canonical_records(data)
        self.assertEqual([(r["geometry"], r["config"]) for r in records],
                         [(g, c) for g in ("raw3d", "xy") for c in PRIMARY])
        fig = plt.figure()
        try:
            axes = canonical_summary(fig, data)
            for ax, key in zip(axes, ("canonical_position", "natural_front_distance", "null_front_distance")):
                points = next(c for c in ax.collections if c.get_gid() == f"canonical:{key}")
                np.testing.assert_allclose(points.get_offsets(),
                                           np.column_stack(([r[key] for r in records], [1, 2, 3, 6, 7, 8])))
        finally:
            plt.close(fig)


class Comparisons(unittest.TestCase):
    def check_comparison(self, family, references, inset_title):
        fig, limits = layouts.compose_comparison(fixture(references), **options(family, 0))
        try:
            self.assertEqual(len(fig.axes), 7)  # Six panels and the color bar.
            self.assertEqual(len(limits["reference_insets"]), 6)
            for i, ax in enumerate(fig.axes[:6]):
                self.assertEqual(len(ax.lines), 3)  # Two zero guides and the front; no metric connectors.
                self.assertEqual(len(ax.child_axes), 1)
                self.assertEqual(ax.child_axes[0].get_title(), inset_title)
                self.assertEqual(ax.get_title(), ("C. elegans protein", "C. elegans RNA", "C. briggsae RNA")[i % 3])
                label = next(t for t in ax.texts if (t.get_gid() or "").startswith("comparison_geometry:"))
                self.assertEqual(label.get_text(), "3D tracking" if i < 3 else "2D (XY) tracking")
                self.assertEqual(label.get_fontweight(), "bold")
                self.assertNotIn("embryo", ax.get_title())
        finally:
            plt.close(fig)

    def test_terminal_figure8(self):
        self.check_comparison("terminal-cross-species", TERMINAL_REFERENCES, "Full-random shuffle")

    def test_partial_forest_figure10(self):
        self.check_comparison("full-tree-cross-species", FULL_TREE_REFERENCES, "Random rebuild")


class Overlays(unittest.TestCase):
    def test_connectors_end_at_closest_point_and_metrics_are_drawn(self):
        for family, references in (("terminal-cross-species", TERMINAL_REFERENCES),
                                   ("full-tree-cross-species", FULL_TREE_REFERENCES)):
            fig, _ = layouts.compose_overlay(fixture(references), **options(family, 1))
            try:
                self.assertEqual(len(fig.axes), 5)
                for ax in fig.axes[:2]:
                    for region in (ax, ax.child_axes[0]):
                        lines = [l for l in region.lines if (l.get_gid() or "").startswith("natural_to_closest:")]
                        self.assertEqual(len(lines), 3)
                        for line in lines:
                            np.testing.assert_allclose(line.get_xdata(), [.32, .3])
                            np.testing.assert_allclose(line.get_ydata(), [.33, .3])
                for ax, key in zip(fig.axes[2:], ("canonical_position", "natural_front_distance", "null_front_distance")):
                    points = next(c for c in ax.collections if c.get_gid() == f"canonical:{key}")
                    self.assertEqual(len(points.get_offsets()), 6)
            finally:
                plt.close(fig)


if __name__ == "__main__":
    unittest.main()
