"""Numbering, caption and archive boundaries for Figures 8 and 9.

Artist-level checks for these figures are in publication/tests/test_cross_species_layouts.py.
"""

import json
from pathlib import Path
import tempfile
import unittest

from terminal_pareto.cross_species_analysis import file_hash
from terminal_pareto.cross_species_publication import archive_previous, captions, write_wrappers


class PublicationTests(unittest.TestCase):
    def test_wrappers_match_publication_style_and_numbering(self):
        with tempfile.TemporaryDirectory() as folder:
            wrappers = write_wrappers(Path(folder), dict(cohort_size=187, cb_reference_size=196))
            self.assertEqual(len(wrappers), 2)
            for number, path in zip((8, 9), wrappers):
                text = path.read_text()
                self.assertIn(r"\documentclass[10pt]{article}", text)
                self.assertIn(rf"\setcounter{{figure}}{{{number - 1}}}", text)
                self.assertIn(f"{path.stem}_panel.pdf", text)
                self.assertIn("187", text)
                self.assertIn("traditional embryo-tracking techniques", text)
                self.assertIn("provisional", text)
                self.assertNotIn("axial scale unverified", text)
            self.assertIn("196-edge reference", captions(dict(cohort_size=187, cb_reference_size=196))[8])
            caption9 = captions(dict(cohort_size=187, cb_reference_size=196))[9]
            self.assertIn(r"\textbf{(C)}", caption9)
            self.assertIn(r"d_{LP}", caption9)
            self.assertIn(r"d_{CP}", caption9)

    def test_rebuild_preserves_previous_publication_with_hashes(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)
            source = run / "publication"
            source.mkdir()
            figure = source / "old.pdf"
            figure.write_bytes(b"previous figure fixture")
            expected = file_hash(figure)
            archive = archive_previous(run)
            self.assertEqual(file_hash(archive / "publication/old.pdf"), expected)
            self.assertEqual(json.loads((archive / "manifest.json").read_text())["files"]["old.pdf"], expected)
            self.assertEqual(file_hash(figure), expected)

    def test_empty_run_has_no_spurious_archive(self):
        with tempfile.TemporaryDirectory() as folder:
            self.assertIsNone(archive_previous(Path(folder)))


if __name__ == "__main__":
    unittest.main()
