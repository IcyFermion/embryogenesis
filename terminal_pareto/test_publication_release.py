"""Release-boundary tests; numerical science is covered by migration validation."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from terminal_pareto import publication_release as release


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        run = self.root / "runs/pooled_tracking_v1/test"
        self.source = run / "publication/endpoint"
        self.source.mkdir(parents=True)
        (run / "analysis").mkdir()
        (run / "analysis/cache.npz").write_bytes(b"unchanged scientific data")
        for stem in release.WRAPPERS:
            (self.source / (stem + ".pdf")).write_bytes(b"new pdf")
            (self.source / (stem + ".tex")).write_text("new wrapper")
        self.published = self.root / "publication"
        self.published.mkdir()
        (self.published / "old.pdf").write_bytes(b"legacy pdf")
        self.validation = patch("terminal_pareto.validate_pooled_migration.validate",
                                return_value={"context_cache_key": "test", "checks": [True]})
        self.pdfinfo = patch.object(release.subprocess, "check_output", return_value="Pages: 1\n")
        self.validation.start()
        self.pdfinfo.start()
        self.addCleanup(self.validation.stop)
        self.addCleanup(self.pdfinfo.stop)

    def test_promotion_freezes_legacy_and_detects_tampering(self):
        release.promote("test", self.root)
        legacy = self.root / "legacy/embryo1/publication/old.pdf"
        self.assertEqual(legacy.read_bytes(), b"legacy pdf")
        release.verify(self.root)
        (self.published / (release.WRAPPERS[0] + ".pdf")).write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "mismatch"):
            release.verify(self.root)

    def test_failed_final_verification_restores_previous_publication(self):
        with patch.object(release, "verify", side_effect=ValueError("verification failed")):
            with self.assertRaisesRegex(ValueError, "verification failed"):
                release.promote("test", self.root)
        self.assertEqual((self.published / "old.pdf").read_bytes(), b"legacy pdf")
        self.assertFalse((self.published / "release_manifest.json").exists())
        self.assertFalse(list(self.root.glob(".publication-stage-*")))

    def test_incomplete_build_does_not_replace_publication(self):
        (self.source / (release.WRAPPERS[-1] + ".tex")).unlink()
        with self.assertRaises(FileNotFoundError):
            release.promote("test", self.root)
        self.assertEqual((self.published / "old.pdf").read_bytes(), b"legacy pdf")

    def comparison_fixture(self):
        run = self.root / "runs/cross_species_terminal_v1/comparison"
        source = run / "publication"
        source.mkdir(parents=True)
        (run / "analysis").mkdir()
        (run / "analysis/cache.npz").write_bytes(b"unchanged comparison data")
        for name in release.cross_species_files() - {release.CROSS_SPECIES_MANIFEST}:
            (source / name).write_bytes(f"comparison {name}".encode())
        assembled = dict(analysis_id="comparison", figure_numbers=dict(zip(release.CROSS_SPECIES_STEMS, (8, 9))),
                         validation=dict(verified_assignments=10), files=release.file_inventory(source),
                         analysis_files=release.file_inventory(run / "analysis"))
        (source / "publication_manifest.json").write_text(json.dumps(assembled))
        return run, assembled

    def test_comparison_is_additive_archived_and_hash_verified(self):
        release.promote("test", self.root)
        before = release.file_inventory(self.published)
        run, assembled = self.comparison_fixture()
        source_before = release.file_inventory(run)
        with patch.object(release, "validate_cross_species_source", return_value=assembled):
            release.promote_cross_species("comparison", self.root)
        record = release.verify(self.root)
        for name, sha in before.items():
            if name != "release_manifest.json":
                self.assertEqual(release.digest(self.published / name), sha)
        previous = self.root / record["cross_species"]["previous_publication"]
        self.assertEqual(release.file_inventory(previous), before)
        self.assertEqual(json.loads((previous.parent / "manifest.json").read_text())["files"], before)
        self.assertEqual(release.file_inventory(run), source_before)
        self.assertEqual(set(record["cross_species"]["files"]), release.cross_species_files())
        (self.published / (release.CROSS_SPECIES_STEMS[0] + "_panel.png")).write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "mismatch"):
            release.verify(self.root)

    def test_future_pooled_release_retains_comparison_assets_and_provenance(self):
        release.promote("test", self.root)
        _, assembled = self.comparison_fixture()
        with patch.object(release, "validate_cross_species_source", return_value=assembled):
            release.promote_cross_species("comparison", self.root)
        previous = release.verify(self.root)["cross_species"]
        release.promote("test", self.root)
        self.assertEqual(release.verify(self.root)["cross_species"], previous)

    def test_comparison_final_failure_rolls_back_without_touching_source(self):
        release.promote("test", self.root)
        before = release.file_inventory(self.published)
        run, assembled = self.comparison_fixture()
        source_before = release.file_inventory(run)
        original_verify = release.verify

        def fail_final(output_root, **kwargs):
            if "publication" not in kwargs and "cross_species" in json.loads(
                    (Path(output_root) / "publication/release_manifest.json").read_text()):
                raise ValueError("injected final failure")
            return original_verify(output_root, **kwargs)

        with patch.object(release, "validate_cross_species_source", return_value=assembled), \
                patch.object(release, "verify", side_effect=fail_final), \
                self.assertRaisesRegex(ValueError, "injected"):
            release.promote_cross_species("comparison", self.root)
        self.assertEqual(release.file_inventory(self.published), before)
        self.assertEqual(release.file_inventory(run), source_before)
        self.assertFalse(list(self.root.glob(".comparison-stage-*")))
        self.assertFalse(list((self.root / "legacy/releases").iterdir()))

    def test_invalid_comparison_source_leaves_release_untouched(self):
        release.promote("test", self.root)
        before = release.file_inventory(self.published)
        with patch.object(release, "validate_cross_species_source", side_effect=ValueError("changed source")), \
                self.assertRaisesRegex(ValueError, "changed source"):
            release.promote_cross_species("comparison", self.root)
        self.assertEqual(release.file_inventory(self.published), before)


if __name__ == "__main__":
    unittest.main()
