"""Pooled full-tree Figure 7/S4 (backlog figures since 2026-11-13): cache-only build and its inventory."""

from __future__ import annotations

from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from publication import registry
from publication.adapters import full_tree as full_tree_adapter
from publication.tests.support import no_experiments

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

    def test_inventory_matches_the_family_assets(self):
        self.assertEqual(set(self.record["files"]), set(registry.FAMILIES["full-tree-pooled"].owned_assets))
        self.assertEqual(self.record["figure_numbers"],
                         {"fig7_ce_full_tree_layerwise": "7", "figs_ce_full_tree_heuristics": "S4"})

    def test_only_vertical_overflow_is_tolerated(self):
        compiled = {stem: figure["compiled"] for stem, figure in self.record["figures"].items()}
        for result in compiled.values():
            self.assertTrue(all("Float too large for page" in w for w in result["warnings"]))
        self.assertEqual({c["pages"] for c in compiled.values()}, {1})


class EndpointDisplay(unittest.TestCase):
    def test_endpoint_transform_does_not_change_raw_scores(self):
        frame = pd.DataFrame(dict(weight_index=[0, 1, 2], travel=[3., 2., 1.], state=[1., 2., 4.]))
        before = frame.copy(deep=True)
        transform = full_tree_adapter.endpoint(frame, "aggregate")
        x, y = transform.transform(frame.travel.to_numpy(), frame.state.to_numpy())
        pd.testing.assert_frame_equal(frame, before)
        np.testing.assert_allclose([x[2], y[2], x[0], y[0]], [0, 1, 1, 0])
        # Reference values outside the endpoints are not clipped.
        nx, ny = transform.transform(5., 7.)
        self.assertGreater(nx, 1)
        self.assertGreater(ny, 1)



class Ownership(unittest.TestCase):
    def test_families_own_disjoint_assets(self):
        owned = [set(family.owned_assets) for family in registry.FAMILIES.values()]
        self.assertTrue(all(owned))
        self.assertEqual(sum(map(len, owned)), len(set().union(*owned)))
