"""Mixed full-tree release on temporary fixtures: Figure 7/S4 in, Figures 10/11 preserved.

Uses the real pooled and cross-species release verifiers; only the scientific
replay of the cross-species run is stubbed. Production is never touched.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

from full_tree_pareto import cross_species_publication as cross_species
from full_tree_pareto import publication_build as pooled
from publication import provenance, release

HISTORICAL = ("brownian_covariance_mle_derivation.pdf", "parametric_brownian_bootstrap.pdf",
              "separate_clock_reference.pdf", "fig7A_ce_full_tree_layerwise_rounds.pdf",
              "fig7A_ce_full_tree_layerwise_rounds.png", "fig7_ce_full_tree_layerwise.pdf",
              "fig7_ce_full_tree_layerwise.tex", "figs_ce_full_tree_heuristics.pdf",
              "figs_ce_full_tree_heuristics_panel.pdf", "figs_ce_full_tree_heuristics_panel.png",
              *release.FULL_TREE_POOLED_RETIRES)
METHODS_PDFS = HISTORICAL[:3]
POOLED = ("fig7A_ce_full_tree_layerwise_rounds.pdf", "fig7A_ce_full_tree_layerwise_rounds.png",
          "fig7B_ce_full_tree_collective.pdf", "fig7B_ce_full_tree_collective.png",
          "fig7_ce_full_tree_layerwise.pdf", "fig7_ce_full_tree_layerwise.tex",
          "figs_ce_full_tree_heuristics.pdf", "figs_ce_full_tree_heuristics.tex",
          "figs_ce_full_tree_heuristics_panel.pdf", "figs_ce_full_tree_heuristics_panel.png")


class MixedFullTreeRelease(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="publication-release-"))
        self.production = self.root / "output/publication"
        self.production.mkdir(parents=True)
        for name in HISTORICAL:
            (self.production / name).write_bytes(f"historical {name}".encode())
        for name in cross_species.publication_files():
            (self.production / name).write_bytes(f"released {name}".encode())
        assembly = json.dumps(dict(analysis_id="toy")).encode()
        (self.production / cross_species.ASSEMBLY_COPY).write_bytes(assembly)
        cs_run = self.root / "cs_run"
        (cs_run / "publication").mkdir(parents=True)
        (cs_run / "publication/publication_manifest.json").write_bytes(assembly)
        files = {name: provenance.sha256(self.production / name)
                 for name in sorted(cross_species.publication_files() | {cross_species.ASSEMBLY_COPY})}
        preserved = {name: provenance.sha256(self.production / name) for name in HISTORICAL}
        (self.production / cross_species.RELEASE_MANIFEST).write_text(json.dumps(dict(
            version="full-tree-cross-species-release-1", run=str(cs_run), analysis_id="toy",
            figure_numbers={stem: number for stem, (number, _) in cross_species.FIGURES.items()},
            files=files, preserved_files=preserved, previous_publication=None)))
        self.cross_species_bytes = {name: (self.production / name).read_bytes() for name in files}
        self.before = release.inventory(self.production)
        run = self.root / "pooled_run"
        (run / "analysis").mkdir(parents=True)
        (run / "analysis/fronts.csv").write_text("method,scope\n")
        self.build = self.root / "build"
        self.build.mkdir()
        for name in POOLED:
            (self.build / name).write_bytes(f"pooled {name}".encode())
        self.write_build_manifest(POOLED, run)
        self.archive_root = self.root / "output/legacy/releases"
        self.verify_patch = mock.patch.object(cross_species, "verify", return_value=dict(analysis_id="toy"))
        self.verify_patch.start()

    def tearDown(self):
        self.verify_patch.stop()
        shutil.rmtree(self.root)

    def write_build_manifest(self, names, run):
        provenance.write_manifest(
            self.build, family="full-tree-pooled", artifacts=[self.build / name for name in names],
            inputs=dict(run=str(run), analysis_id="pooled",
                        input_files={"analysis/fronts.csv": provenance.sha256(run / "analysis/fronts.csv")}),
            extra=dict(figure_numbers={"fig7_ce_full_tree_layerwise": "7", "figs_ce_full_tree_heuristics": "S4"}))

    def run_release(self, **options):
        return release.release_full_tree_pooled(self.build, self.production, archive_root=self.archive_root,
                                                when="20261003T000000000000Z", **options)

    def test_release_preserves_figures_10_11_and_both_verifiers_agree(self):
        result = self.run_release()
        for name, data in self.cross_species_bytes.items():
            self.assertEqual((self.production / name).read_bytes(), data, name)
        for name in POOLED:
            self.assertEqual((self.production / name).read_bytes(), (self.build / name).read_bytes(), name)
        for name in release.FULL_TREE_POOLED_RETIRES:
            self.assertFalse((self.production / name).exists(), name)
        for name in METHODS_PDFS:
            self.assertEqual((self.production / name).read_bytes(), f"historical {name}".encode())
        pooled.verify_release(self.production)
        cross = cross_species.verify_release(self.production, check_archive=False)
        self.assertEqual(set(cross["preserved_files"]) & set(POOLED), set(POOLED))
        self.assertIn(release.FULL_TREE_RELEASE_MANIFEST, cross["preserved_files"])
        self.assertEqual(len(cross["preserved_files_history"]), 1)
        archived = self.archive_root / "20261003T000000000000Z"
        release.check(archived / "publication", self.before)
        self.assertEqual(json.loads((archived / "manifest.json").read_text())["files"], self.before)
        self.assertEqual(sorted(result["removed"]), sorted(release.FULL_TREE_POOLED_RETIRES))

    def test_release_refuses_assets_owned_by_another_family(self):
        intruder = "fig10_full_tree_cross_species_comparison.pdf"
        (self.build / intruder).write_bytes(b"not mine")
        self.write_build_manifest((*POOLED, intruder), Path(json.loads(
            (self.build / provenance.MANIFEST).read_text())["run"]))
        with self.assertRaisesRegex(ValueError, "does not own"):
            self.run_release()
        self.assertEqual(release.inventory(self.production), self.before)
        self.assertFalse(self.archive_root.exists())

    def test_failed_final_verification_restores_production(self):
        calls = []

        def flaky(directory):
            calls.append(directory)
            if len(calls) == 2:
                raise ValueError("injected final failure")
            return pooled.verify_release(directory)
        with self.assertRaisesRegex(ValueError, "injected"):
            self.run_release(verifiers=(flaky, cross_species.verify_release))
        self.assertEqual(release.inventory(self.production), self.before)
        self.assertFalse((self.archive_root / "20261003T000000000000Z").exists())

    def test_failed_staged_verification_leaves_production_untouched(self):
        def failing(directory):
            raise ValueError("staged verification failed")
        with self.assertRaisesRegex(ValueError, "staged"):
            self.run_release(verifiers=(failing, cross_species.verify_release))
        self.assertEqual(release.inventory(self.production), self.before)
        self.assertFalse(self.archive_root.exists())
        self.assertEqual([p.name for p in self.production.parent.iterdir() if p.name.startswith(".release")], [])

    def test_altered_build_is_refused_before_staging(self):
        (self.build / POOLED[0]).write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "Changed or missing build artifacts"):
            self.run_release()
        self.assertEqual(release.inventory(self.production), self.before)


if __name__ == "__main__":
    unittest.main()
