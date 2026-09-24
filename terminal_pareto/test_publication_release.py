"""Release-boundary tests; numerical science is covered by migration validation."""
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


if __name__ == "__main__":
    unittest.main()
