"""Per-family presentation dependencies: only a family's own sources make its builds stale."""

from __future__ import annotations

import unittest
from unittest import mock

from publication import provenance
from publication.registry import FAMILIES

SHARED = {"publication/notation.py", "publication/style.py", "publication/assembly.py",
          "publication/provenance.py", "publication/registry.py", "publication/contracts.py"}


def deps(family):
    return set(provenance.presentation_sources(family))


class FamilyDependencies(unittest.TestCase):
    def test_every_family_tracks_the_shared_core(self):
        for key in FAMILIES:
            self.assertLessEqual(SHARED, deps(key), key)

    def test_terminal_primary_tracks_its_layouts_and_asset(self):
        terminal = deps("terminal-primary")
        for name in ("figures/terminal_within_type.py", "figures/terminal_subtree_map.py", "captions/terminal.py",
                     "tables.py", "adapters/terminal.py", "assets/fig3A_ce_null_models.tex"):
            self.assertIn(f"publication/{name}", terminal)
        self.assertNotIn("publication/figures/cross_species.py", terminal)

    def test_cross_species_families_exclude_other_families_layouts(self):
        for key in ("terminal-cross-species", "full-tree-cross-species"):
            found = deps(key)
            self.assertIn("publication/figures/cross_species.py", found)
            self.assertIn("publication/captions/cross_species.py", found)
            # The terminal adapter module serves two families; only the called loader's imports count.
            for name in ("figures/terminal_within_type.py", "figures/terminal_subtree_map.py",
                         "figures/full_tree.py", "captions/terminal.py", "assets/fig3A_ce_null_models.tex"):
                self.assertNotIn(f"publication/{name}", found, key)

    def test_full_tree_pooled_excludes_terminal_layouts(self):
        found = deps("full-tree-pooled")
        self.assertIn("publication/figures/full_tree.py", found)
        self.assertFalse({name for name in found if "/terminal" in name or "cross_species" in name})

    def test_maintenance_tools_and_tests_never_count(self):
        everything = set(provenance.presentation_sources())
        for key in FAMILIES:
            for name in deps(key) | everything:
                self.assertNotIn("/tests/", name)
                self.assertNotIn(name, {"publication/migration.py", "publication/parity.py"})

    def test_changes_are_judged_against_the_familys_own_dependencies(self):
        recorded = {key: provenance.presentation_sources(key) for key in ("terminal-cross-species", "terminal-primary")}
        added = dict(recorded["terminal-cross-species"], **{"publication/new_dependency.py": "0" * 64})
        self.assertEqual(provenance.changed_sources(recorded["terminal-cross-species"], added),
                         ["publication/new_dependency.py"])
        real = provenance.sha256

        def edited(path):  # simulate an edit to a Figure 6 module only
            return "0" * 64 if str(path).endswith("figures/terminal_within_type.py") else real(path)
        with mock.patch.object(provenance, "sha256", edited):
            live = {key: provenance.presentation_sources(key) for key in recorded}
        self.assertEqual(provenance.changed_sources(recorded["terminal-cross-species"], live["terminal-cross-species"]), [])
        self.assertEqual(provenance.changed_sources(recorded["terminal-primary"], live["terminal-primary"]),
                         ["publication/figures/terminal_within_type.py"])

if __name__ == "__main__":
    unittest.main()
