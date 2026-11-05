"""Track E: Figure 6B-C panels (revised layout) and supplement S5 from cache-only inputs.

Uses the published numerical runs read-only. Run from the repository root:

    python -m unittest discover -s publication/tests -t .
"""

from __future__ import annotations

from pathlib import Path
import shutil
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from publication import provenance, registry, style
from publication.adapters import full_tree as adapter
from publication.figures.full_tree_layers import compose_cost, compose_position, compose_travel_decomposition
from publication.tests.support import no_experiments

HAS_TEX = shutil.which("tectonic") is not None or shutil.which("pdflatex") is not None
FAMILY = "full-tree-layers"


class LayerInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        style.configure()
        with no_experiments():
            cls.data = adapter.layers()
        from full_tree_pareto import layer_features as lf
        cls.saved = pd.read_csv(lf.run_path(adapter.LAYER_RUNS["hourglass"]) / "analysis/layers.csv")

    def test_layers_three_to_eight(self):
        t = self.data.layers
        self.assertEqual(list(t.layer), [3, 4, 5, 6, 7, 8])
        self.assertEqual(list(t.edges), [18, 38, 66, 123, 223, 491])
        self.assertEqual(self.data.caption["left_out"], [5, 10])
        self.assertEqual((self.data.caption["first_groups"], self.data.caption["first_points"]), (6, 5))
        self.assertEqual(self.data.caption["cells"], 959)

    def test_relative_cost_is_one_minus_the_saved_advantage(self):
        rows = self.saved.query("series == 'pooled_ce_protein' and family == 'first_cousin' and layer >= 3")
        rows = rows.sort_values("layer")
        t = self.data.layers
        for column in ("travel", "state"):
            np.testing.assert_allclose(t[column], 1 - rows[f"{column}_advantage"], rtol=0, atol=1e-12)
            np.testing.assert_allclose(t[f"{column}_low"], 1 - rows[f"{column}_high"], rtol=0, atol=1e-12)
            np.testing.assert_allclose(t[f"{column}_high"], 1 - rows[f"{column}_low"], rtol=0, atol=1e-12)
        # Lineage beats the shuffle on travel at every drawn layer, by the least in layers 6-7.
        self.assertTrue((t.travel_high < 1).all())
        self.assertEqual(set(t.nlargest(2, "travel").layer), {6, 7})

    def test_fewest_rewired_is_one_minus_the_maximum_similarity(self):
        from full_tree_pareto import pooled_analysis as pa
        from full_tree_pareto import round_references as rr
        with no_experiments():
            metrics = rr.load(pa.DEFAULT_RUN, "pooled")["metrics"]
        metrics = metrics[metrics.config == "ce_protein"].assign(layer=lambda f: 9 - f["round"]).set_index("layer")
        t = self.data.layers.set_index("layer")
        np.testing.assert_allclose(t.fewest_rewired, 1 - metrics.loc[t.index, "max_edge_retention"], rtol=0, atol=1e-12)
        self.assertEqual(t.fewest_rewired.idxmax(), 6)

    def test_decomposition_signs_match_the_caption(self):
        changes = self.data.changes.set_index(["tree", "contrast"])
        for tree in ("pooled", "partial"):
            for contrast in ("L6 vs L8", "L7 vs L8", "L6 vs L4-5"):
                row = changes.loc[(tree, contrast)]
                self.assertGreater(row.aunt, -5, (tree, contrast))
                self.assertGreater(row.natural, row.aunt, (tree, contrast))
            row = changes.loc[(tree, "L7 vs L4-5")]
            self.assertLess(row.aunt, row.natural)
            self.assertLess(row.aunt, 0)
        random = self.data.per_edge.query("tree == 'pooled'").random * 1e3
        self.assertLess(random.max() - random.min(), .2)

    def test_figures_draw_the_input_values(self):
        t = self.data.layers
        for compose, columns in ((compose_position, ("front_distance", "fewest_rewired")), (compose_cost, ("travel", "state"))):
            fig, layout = compose(self.data)
            try:
                points = {child.get_gid(): child.get_offsets() for ax in fig.axes for child in ax.get_children()
                          if (child.get_gid() or "").startswith("layer:")}
                self.assertEqual(set(points), {f"layer:{c}" for c in columns})
                for column in columns:
                    scale = 100 if column == "fewest_rewired" else 1  # drawn as a percentage
                    np.testing.assert_allclose(points[f"layer:{column}"][:, 1], t[column] * scale, atol=1e-10)
                self.assertEqual(layout["layers"], [3, 4, 5, 6, 7, 8])
            finally:
                plt.close(fig)
        low, high = layout["relative_limits"]
        self.assertLessEqual(low, min(t.travel_low.min(), t.state_low.min()))
        self.assertGreaterEqual(high, max(t.travel_high.max(), t.state_high.max()))
        fig, layout = compose_travel_decomposition(self.data)
        plt.close(fig)
        self.assertEqual(layout["contrasts"], ["L6 vs L8", "L7 vs L8", "L6 vs L4-5", "L7 vs L4-5"])

    def test_dependencies_exclude_other_families_layouts(self):
        found = set(provenance.presentation_sources(FAMILY))
        for name in ("figures/full_tree_layers.py", "captions/full_tree_layers.py", "adapters/full_tree.py"):
            self.assertIn(f"publication/{name}", found)
        self.assertFalse({name for name in found if "/terminal" in name or "cross_species" in name
                          or name.endswith("figures/full_tree.py")})


@unittest.skipUnless(HAS_TEX, "tectonic or pdflatex required")
class CacheOnlyLayerBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="publication-layers-"))
        with no_experiments():
            cls.record = registry.build(FAMILY, cls.tmp / "build")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_inventory_and_numbers(self):
        self.assertEqual(set(self.record["files"]), set(registry.FAMILIES[FAMILY].owned_assets))
        self.assertEqual(self.record["figure_numbers"], {"figS5_full_tree_layer_travel": "S5"})
        self.assertEqual(provenance.verify_build(self.tmp / "build")["ready"], True)

    def test_compiles_on_one_page_without_warnings(self):
        for figure in self.record["figures"].values():
            self.assertEqual(figure["compiled"]["pages"], 1)
            self.assertEqual(figure["compiled"]["warnings"], [])


if __name__ == "__main__":
    unittest.main()
