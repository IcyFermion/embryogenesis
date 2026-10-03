import json
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from full_tree_pareto import publication_build as pb
from full_tree_pareto import pooled_analysis as pa


class BuildTests(unittest.TestCase):
    def test_shared_worker_defaults_and_overrides(self):
        from full_tree_pareto import resume_paired as rp
        self.assertEqual(pb.DEFAULT_WORKERS, 24)
        self.assertEqual(rp.run.__defaults__[1], pb.DEFAULT_WORKERS)
        for parse in (pb.parse_args, rp.parse_args):
            self.assertEqual(parse([]).workers, 24)
            for count in (1, 6, 12, 32):
                self.assertEqual(parse(["--workers", str(count)]).workers, count)
            for invalid in ("0", "-1", "1.5", "abc"):
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    parse(["--workers", invalid])

    def test_endpoint_transform_does_not_change_raw_scores(self):
        frame = pd.DataFrame(dict(weight_index=[0, 1, 2], travel=[3., 2., 1.], state=[1., 2., 4.]))
        before = frame.copy(deep=True)
        transform = pb.endpoint(frame, "aggregate")
        x, y = transform.transform(frame.travel.to_numpy(), frame.state.to_numpy())
        pd.testing.assert_frame_equal(frame, before)
        np.testing.assert_allclose([x[2], y[2], x[0], y[0]], [0, 1, 1, 0])
        # Reference values outside the endpoints are not clipped.
        nx, ny = transform.transform(5., 7.)
        self.assertGreater(nx, 1)
        self.assertGreater(ny, 1)

    def test_release_hashes_check_both_figures_and_caches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run, publication = root / "run", root / "publication"
            (run / "analysis").mkdir(parents=True)
            publication.mkdir()
            figure = publication / "figure.pdf"
            (run / "analysis/cousin_shuffles").mkdir()
            cache = run / "analysis/cousin_shuffles/result.npz"
            figure.write_bytes(b"test figure fixture")
            cache.write_bytes(b"test cache fixture")
            record = dict(run=str(run), files={figure.name: pa.digest(figure)},
                          analysis_files={str(cache.relative_to(run / "analysis")): pa.digest(cache)})
            pa.write_json(publication / "release_manifest.json", record)
            self.assertEqual(pb.verify_release(publication), record)
            cache.write_bytes(b"tampered")
            with self.assertRaises(ValueError):
                pb.verify_release(publication)
            cache.write_bytes(b"test cache fixture")
            figure.write_bytes(b"tampered figure")
            with self.assertRaises(ValueError):
                pb.verify_release(publication)

    def test_working_layout_archive_retains_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "publication").mkdir()
            figure = run / "publication/figure.pdf"
            figure.write_bytes(b"preceding layout")
            archive = pb.archive_working_publication(run)
            self.assertEqual((archive / "publication/figure.pdf").read_bytes(), figure.read_bytes())
            record = json.loads((archive / "manifest.json").read_text())
            self.assertEqual(record["files"]["figure.pdf"], pa.digest(figure))
            figure.write_bytes(b"new layout")
            self.assertEqual((archive / "publication/figure.pdf").read_bytes(), b"preceding layout")

    def test_wrapper_captions_follow_sweep_settings(self):
        ctx = pa.build_context()
        with tempfile.TemporaryDirectory() as tmp:
            pb.write_wrappers(ctx, Path(tmp))
            for stem in pb.WRAPPERS:
                text = (Path(tmp) / f"{stem}.tex").read_text()
                for token in ("INTERNAL", "LEAVES", "EDGES", "ROUNDS", "NWEIGHTS", "NSHOWN", "NDRAWS"):
                    self.assertNotIn(token, text)
                self.assertIn(f"{pb.INTERVALS+1:,}", text)
                self.assertIn(f"{pb.DRAWS:,}", text)
            main = (Path(tmp) / f"{pb.WRAPPERS[0]}.tex").read_text()
            supplement = (Path(tmp) / f"{pb.WRAPPERS[1]}.tex").read_text()
            self.assertNotIn("Gaussian", main)
            self.assertIn("third-cousin", main)
            self.assertIn("retained only here", supplement)
        self.assertNotIn(pa.REFERENCE, pb.MAIN_REFERENCES)
        self.assertIn(pa.REFERENCE, pb.SUPPLEMENT_REFERENCES)

    def test_failed_promotion_restores_previous_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "full_tree_pareto/output"
            run = output / "runs/test"
            (run / "analysis").mkdir(parents=True)
            (run / "publication").mkdir()
            (run / "analysis" / "cache.npz").write_bytes(b"cache")
            (run / "publication" / "figure.pdf").write_bytes(b"new")
            old = output / "publication"
            old.mkdir()
            (old / "figure.pdf").write_bytes(b"old")
            original = pb.verify_release
            calls = []

            def fail_on_target(path):
                calls.append(path)
                if Path(path) == old:
                    raise ValueError("simulated post-swap failure")
                return original(path)

            with mock.patch.object(pa, "ROOT", root), \
                    mock.patch.object(pa, "build", return_value=(object(), None, None, None)), \
                    mock.patch.object(pb.cr, "build"), \
                    mock.patch.object(pb.subprocess, "check_output", return_value="Pages: 1\n"), \
                    mock.patch.object(pb, "verify_release", side_effect=fail_on_target):
                with self.assertRaises(ValueError):
                    pb.promote(run)
            self.assertEqual(len(calls), 2)
            self.assertEqual((old / "figure.pdf").read_bytes(), b"old")
            self.assertEqual(sorted(p.name for p in old.iterdir()), ["figure.pdf"])
            self.assertFalse((output / "legacy").exists() and any((output / "legacy/releases").iterdir()))
            self.assertEqual([p for p in output.iterdir() if p.name.startswith(".pooled-stage-")], [])

    def test_pooled_promotion_refuses_to_replace_cross_species_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            publication = root / "full_tree_pareto/output/publication"
            publication.mkdir(parents=True)
            manifest = publication / "cross_species_release_manifest.json"
            manifest.write_text('{"release_promoted": true}\n')
            figure = publication / "fig10_full_tree_cross_species_comparison.pdf"
            figure.write_bytes(b"accepted comparison")
            before = {p.name: p.read_bytes() for p in publication.iterdir()}
            with mock.patch.object(pa, "ROOT", root), mock.patch.object(pa, "build") as build:
                with self.assertRaisesRegex(ValueError, "Figures 10/11"):
                    pb.promote(root / "unbuilt-run")
                build.assert_not_called()
            self.assertEqual({p.name: p.read_bytes() for p in publication.iterdir()}, before)
            self.assertEqual(list(publication.parent.iterdir()), [publication])


if __name__ == "__main__":
    unittest.main()
