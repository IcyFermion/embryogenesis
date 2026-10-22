"""Reserve notes: cache-only build, one page, and statement checks."""

from __future__ import annotations

import copy
from pathlib import Path
import shutil
import tempfile
import unittest

from publication.notes import apoptosis
from publication.tests.support import no_experiments

HAS_TEX = shutil.which("tectonic") is not None or shutil.which("pdflatex") is not None


class ApoptosisNote(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with no_experiments():
            tables, cls.inputs = apoptosis.load()
        cls.facts = apoptosis.facts(tables)

    def test_statements_hold_and_fail_when_the_data_change(self):
        apoptosis.check(self.facts)
        changed = copy.deepcopy(self.facts)
        changed["travel"][1] = (3.0, 3.4)
        with self.assertRaisesRegex(ValueError, "farther in every embryo"):
            apoptosis.check(changed)
        changed = copy.deepcopy(self.facts)
        changed["rows"][0]["right_max"] = .4
        with self.assertRaisesRegex(ValueError, "right partner off its front"):
            apoptosis.check(changed)

    def test_inputs_are_the_four_track_d_runs(self):
        self.assertEqual({name.split("/")[0] for name in self.inputs}, set(apoptosis.RUNS.values()))

    @unittest.skipUnless(HAS_TEX, "tectonic or pdflatex required")
    def test_builds_to_one_page_without_computation(self):
        with tempfile.TemporaryDirectory() as tmp, no_experiments():
            record = apoptosis.build(Path(tmp))
            self.assertEqual(record["compiled"]["pages"], 1)
            self.assertEqual(record["compiled"]["warnings"], [])
            self.assertTrue((Path(tmp) / record["file"]).is_file())


if __name__ == "__main__":
    unittest.main()
