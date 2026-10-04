"""Unified production releases on temporary fixtures (production itself is never touched)."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

from publication import provenance, release
from publication.registry import BUILD_SET, FAMILIES

TERMINAL, PARTIAL = "terminal-cross-species", "full-tree-cross-species"


class ProductionRelease(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="publication-release-"))
        self.run = self.root / "run"
        (self.run / "analysis").mkdir(parents=True)
        (self.run / "analysis/metrics.csv").write_text("x\n1\n")
        self.production = self.root / "production"
        self.archive = self.root / "archive"
        self.set_dir = self.root / "set"
        self.make_set("v1")

    def tearDown(self):
        shutil.rmtree(self.root)

    def make_family(self, directory: Path, key: str, tag: str, names=None):
        directory.mkdir(parents=True)
        names = list(names or FAMILIES[key].owned_assets)
        for name in names:
            (directory / name).write_bytes(f"{tag} {name}".encode())
        provenance.write_manifest(
            directory, family=key, artifacts=[directory / name for name in names],
            inputs=dict(run=str(self.run), analysis_id=f"{key}-id",
                        input_files={"analysis/metrics.csv": provenance.sha256(self.run / "analysis/metrics.csv")}),
            extra=dict(figure_numbers=FAMILIES[key].numbers))

    def make_set(self, tag, names=None):
        if self.set_dir.exists():
            shutil.rmtree(self.set_dir)
        for key in (TERMINAL, PARTIAL):
            self.make_family(self.set_dir / key, key, tag, names.get(key) if names else None)
        (self.set_dir / BUILD_SET).write_text(json.dumps(dict(families={
            key: dict(directory=key) for key in (TERMINAL, PARTIAL)})))

    def apply(self, **options):
        return release.release(self.set_dir, production=self.production, archive_root=self.archive, **options)

    def test_first_release_creates_flat_production_with_family_sections(self):
        result = self.apply(when="T1")
        status = release.verify_production(self.production)
        self.assertEqual(set(status["families"]), {TERMINAL, PARTIAL})
        manifest = json.loads((self.production / release.PRODUCTION_MANIFEST).read_text())
        self.assertIsNone(manifest["history"][0]["previous"])
        for key in (TERMINAL, PARTIAL):
            for name in FAMILIES[key].owned_assets:
                self.assertEqual((self.production / name).read_bytes(), f"v1 {name}".encode())
        self.assertEqual(sorted(result["families"]), sorted([TERMINAL, PARTIAL]))

    def test_family_rerelease_changes_only_its_assets_and_retires_dropped_ones(self):
        self.apply(when="T1")
        other = {name: (self.production / name).read_bytes() for name in FAMILIES[PARTIAL].owned_assets}
        dropped = FAMILIES[TERMINAL].owned_assets[0]
        kept = [n for n in FAMILIES[TERMINAL].owned_assets if n != dropped]
        self.make_set("v2", names={TERMINAL: kept})
        before = release.inventory(self.production)
        result = self.apply(families=[TERMINAL], when="T2")
        self.assertEqual(result["removed"], [dropped])
        for name, data in other.items():
            self.assertEqual((self.production / name).read_bytes(), data)
        for name in kept:
            self.assertEqual((self.production / name).read_bytes(), f"v2 {name}".encode())
        release.check(self.archive / "T2/publication", before)
        release.verify_production(self.production)

    def test_build_listing_another_familys_asset_is_refused(self):
        intruder = FAMILIES[PARTIAL].owned_assets[0]
        shutil.copy2(self.set_dir / PARTIAL / intruder, self.set_dir / TERMINAL / intruder)
        names = [*FAMILIES[TERMINAL].owned_assets, intruder]
        shutil.rmtree(self.set_dir / TERMINAL)
        self.make_family(self.set_dir / TERMINAL, TERMINAL, "v1", names)
        with self.assertRaisesRegex(ValueError, "does not own"):
            self.apply()
        self.assertFalse(self.archive.exists())

    def test_failed_final_verification_restores_previous_production(self):
        self.apply(when="T1")
        before = release.inventory(self.production)
        self.make_set("v2")
        calls = []
        real = release.verify_production

        def flaky(directory, **kwargs):
            calls.append(directory)
            if len(calls) == 2:
                raise ValueError("injected final failure")
            return real(directory, **kwargs)
        with mock.patch.object(release, "verify_production", flaky), \
                self.assertRaisesRegex(ValueError, "injected"):
            self.apply(when="T2")
        self.assertEqual(release.inventory(self.production), before)
        self.assertFalse((self.archive / "T2").exists())

    def test_stale_build_is_refused(self):
        live = provenance.presentation_sources()
        with mock.patch.object(provenance, "presentation_sources",
                               return_value=dict(live, **{"publication/style.py": "0" * 64})), \
                self.assertRaisesRegex(ValueError, "stale"):
            self.apply()

    def test_verify_detects_tampering_and_strays(self):
        self.apply(when="T1")
        (self.production / "stray.pdf").write_bytes(b"?")
        with self.assertRaisesRegex(ValueError, "Unexpected or missing"):
            release.verify_production(self.production)
        (self.production / "stray.pdf").unlink()
        name = FAMILIES[TERMINAL].owned_assets[0]
        (self.production / name).write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "Hash mismatch"):
            release.verify_production(self.production)

    def test_rehearsal_leaves_production_untouched(self):
        self.apply(when="T1")
        before = release.inventory(self.production)
        self.make_set("v2")
        result = release.rehearse(self.set_dir, production=self.production)
        self.assertTrue(result["production_unchanged"])
        self.assertEqual(release.inventory(self.production), before)


if __name__ == "__main__":
    unittest.main()
