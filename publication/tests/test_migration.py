"""Figure-material migration: archive, verify and guarded purge on a temporary tree."""

from __future__ import annotations

from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

from publication import migration


class MigrationArchive(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="publication-migration-"))
        (self.root / "runs/a/publication").mkdir(parents=True)
        (self.root / "runs/a/publication/fig.pdf").write_bytes(b"figure")
        (self.root / "runs/a/analysis").mkdir()
        (self.root / "runs/a/analysis/metrics.csv").write_text("numbers\n")
        (self.root / "runs/a/build.log").write_text("render log\n")
        batches = {migration.DEFAULT_BATCH: ("runs/a/publication", "runs/a/build.log", "runs/missing")}
        self.patches = [mock.patch.object(migration, "ROOT", self.root),
                        mock.patch.object(migration, "ARCHIVE_ROOT", self.root / "archive"),
                        mock.patch.object(migration, "BATCHES", batches)]
        for patch in self.patches:
            patch.start()

    def tearDown(self):
        for patch in self.patches:
            patch.stop()
        shutil.rmtree(self.root)

    def test_archive_then_purge_removes_only_figure_material(self):
        self.assertEqual(migration.archive(), dict(paths=2, files=2))
        self.assertEqual(migration.purge()["removed"], ["runs/a/publication", "runs/a/build.log"])
        self.assertFalse((self.root / "runs/a/publication").exists())
        self.assertTrue((self.root / "runs/a/analysis/metrics.csv").exists())
        self.assertEqual((self.root / "archive" / migration.DEFAULT_BATCH / "runs/a/publication/fig.pdf").read_bytes(), b"figure")

    def test_purge_refuses_changed_originals_and_corrupt_archives(self):
        migration.archive()
        (self.root / "runs/a/publication/new.pdf").write_bytes(b"added later")
        with self.assertRaisesRegex(ValueError, "changed since archiving"):
            migration.purge()
        (self.root / "runs/a/publication/new.pdf").unlink()
        (self.root / "archive" / migration.DEFAULT_BATCH / "runs/a/build.log").write_text("corrupted")
        with self.assertRaisesRegex(ValueError, "archive mismatch"):
            migration.purge()
        self.assertTrue((self.root / "runs/a/publication/fig.pdf").exists())

    def test_archive_is_written_once(self):
        migration.archive()
        with self.assertRaises(FileExistsError):
            migration.archive()


CHECKPOINT = ("terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920/validation/"
              "organization_checkpoint_20260922/accepted_publication_sha256.json")
EMBRYO1 = "terminal_pareto/output/legacy/embryo1/publication"


class ArchivedEmbryo1Figures(unittest.TestCase):
    """The accepted pre-pooled figures (checkpoint 2026-09-22) survive intact in the legacy archive."""

    def test_archive_matches_pre_organization_checkpoint(self):
        from publication.provenance import ROOT, sha256
        archived = migration.ARCHIVE_ROOT / "legacy_20261004" / EMBRYO1
        if not archived.exists():
            self.skipTest("legacy figure archive not created on this machine")
        record = __import__("json").loads((ROOT / CHECKPOINT).read_text())
        mismatched = [name for name, digest in record["files"].items()
                      if not (archived / name).is_file() or sha256(archived / name) != digest]
        self.assertEqual(mismatched, [])


if __name__ == "__main__":
    unittest.main()
