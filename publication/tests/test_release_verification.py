"""Frozen-release integrity versus live-source readiness in the back-end release verifiers."""

from __future__ import annotations

import unittest
from unittest import mock

from full_tree_pareto import cross_species_publication as full_tree_release
from terminal_pareto import cross_species_analysis as tcs
from terminal_pareto import publication_release as terminal_release

STALE = ["terminal_pareto/fig_terminal_cross_species.py"]


class StaleSourceIsNotCorruption(unittest.TestCase):
    def test_full_tree_production_verifies_with_stale_presentation_source(self):
        with mock.patch.object(full_tree_release, "stale_files", return_value=STALE):
            record = full_tree_release.verify_release()
        self.assertEqual(record["stale_presentation_sources"], STALE)

    def test_full_tree_readiness_rejects_stale_presentation_source(self):
        run = full_tree_release.run_path(full_tree_release.DEFAULT_RUN)
        with mock.patch.object(full_tree_release, "stale_files", return_value=STALE), \
                self.assertRaisesRegex(ValueError, "Changed presentation source"):
            full_tree_release.verify(run, current_source=True)

    def test_terminal_comparison_source_verifies_with_stale_presentation_source(self):
        run = tcs.run_path(tcs.DEFAULT_RUN)
        with mock.patch.object(terminal_release, "stale_files", return_value=STALE):
            record = terminal_release.validate_cross_species_source(run)
        self.assertEqual(record["figure_numbers"], {"fig8_terminal_cross_species_comparison": 8,
                                                    "fig9_terminal_cross_species_overlays": 9})

    def test_frozen_artifact_change_still_fails(self):
        run = full_tree_release.run_path(full_tree_release.DEFAULT_RUN)
        real = full_tree_release.check_hashes

        def tampered(base, hashes):
            if str(base).endswith("figures"):
                raise ValueError("Hash mismatch: simulated")
            return real(base, hashes)
        with mock.patch.object(full_tree_release, "check_hashes", tampered), \
                self.assertRaisesRegex(ValueError, "simulated"):
            full_tree_release.verify(run)


if __name__ == "__main__":
    unittest.main()
