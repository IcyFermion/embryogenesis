"""Pooled full-tree pipeline settings (figure tests live in publication/tests)."""
import contextlib
import io
import unittest

from full_tree_pareto import pooled_pipeline as pp


class PipelineTests(unittest.TestCase):
    def test_shared_worker_defaults_and_overrides(self):
        from full_tree_pareto import resume_paired as rp
        self.assertEqual(pp.DEFAULT_WORKERS, 24)
        self.assertEqual(rp.run.__defaults__[1], pp.DEFAULT_WORKERS)
        for parse in (pp.parse_args, rp.parse_args):
            self.assertEqual(parse([]).workers, 24)
            for count in (1, 6, 12, 32):
                self.assertEqual(parse(["--workers", str(count)]).workers, count)
            for invalid in ("0", "-1", "1.5", "abc"):
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    parse(["--workers", invalid])

    def test_sweep_settings_are_the_cached_identity_settings(self):
        self.assertEqual(pp.sweep_settings(), dict(intervals=300, draws=10000, seed=42))


if __name__ == "__main__":
    unittest.main()
