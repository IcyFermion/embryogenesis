"""Compute the pooled terminal analysis caches behind Figures 1-6, S1-S3 and Table 1.

python terminal_pareto/analysis_pipeline.py [--profile pooled_tracking_v1] [--run-id ID]

Numerical only: every step writes into the run's ``analysis/`` folder (logs in
``validation/``). Figures, captions and Table 1 are produced by
``python -m publication build --family terminal-primary`` and released with
``python -m publication release``. Validate the result with
``python terminal_pareto/validate_pooled_migration.py``.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto.analysis_context import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
    validate_existing_context_manifest,
)

S3_DIR = "s3_cross_geometry"
S3_FILES = ("cross_geometry_coordinates.csv", "provenance.json")
TRACKING_CACHE = ROOT / "terminal_pareto/output/tracking_geometry_sensitivity/assignments_and_costs.npz"


def _run(command: list[str], *, log: Path, env: dict[str, str]) -> None:
    result = subprocess.run(command, cwd=ROOT, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(result.stdout)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout)


def _script(name: str, *args: str) -> list[str]:
    return [sys.executable, str(ROOT / "terminal_pareto" / name), *args]


def build(profile: str, run_id: str, output_root: Path) -> Path:
    context = build_analysis_context(profile, run_id=run_id, output_root=output_root)
    if (context.run_paths.analysis / "analysis_manifest.json").exists():
        validate_existing_context_manifest(context)
    else:
        context.write()
    paths = context.run_paths.create()
    env = dict(os.environ)
    env.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="terminal-mpl-"))
    common = ["--profile", profile, "--run-id", run_id, "--output-root", str(output_root)]
    # fig5_table1 and figS3_cross_geometry are validator-pinned and also write
    # presentation files; those go to a discarded folder inside the repository
    # (fig5_table1 reports its output path relative to the repository root).
    with tempfile.TemporaryDirectory(prefix=".discard-", dir=ROOT / "terminal_pareto" / "output") as discard:
        steps = [
            ("global_analysis", _script("fig2_fig3_ce_terminal_pareto.py", *common)),
            ("subtree_analysis", _script("subtree_analysis.py", *common, "--min-cells", "12")),
            ("canonical_metrics", _script("fig5_table1_ce_canonical_metrics.py", *common, "--min-cells", "12",
                                          "--iteration", "300", "--out", discard)),
            ("cell_type_caches", _script("fig6a_figs2_ce_cell_types.py", *common, "--iteration", "300")),
            ("within_type", _script("fig6bc_ce_within_type.py", *common, "--iteration", "300")),
            ("s3_projections", _script("figS3_cross_geometry.py", "--run-id", run_id,
                                       "--output-root", str(output_root), "--tracking-cache", str(TRACKING_CACHE),
                                       "--out", discard, "--publication")),
        ]
        for label, command in steps:
            _run(command, log=paths.validation / f"{label}.log", env=env)
        target = paths.analysis / S3_DIR
        target.mkdir(parents=True, exist_ok=True)
        for name in S3_FILES:
            shutil.copy2(Path(discard) / name, target / name)
    # The (hash-pinned) context helper creates display folders; figures live in publication/.
    for folder in sorted((paths.root / "publication").glob("**/"), reverse=True):
        if folder.is_dir() and not any(folder.iterdir()):
            folder.rmdir()
    return paths.root


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default="pooled_tracking_v1", choices=("pooled_tracking_v1",))
    parser.add_argument("--run-id", default="migration_candidate_20260920")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args(argv)
    root = build(args.profile, args.run_id, args.output_root)
    print(f"Pooled terminal analysis caches complete: {root / 'analysis'}")


if __name__ == "__main__":
    main()
