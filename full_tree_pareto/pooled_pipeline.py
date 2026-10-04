"""Build or replay-validate the pooled full-tree caches (Figure 7 and supplement inputs).

python -m full_tree_pareto.pooled_pipeline [--run PATH] [--workers N] [--layout-only]

Numerical only: figures are drawn by ``python -m publication build --family
full-tree-pooled`` and released with ``python -m publication release``.
``--layout-only`` replays existing caches without solving (the former
layout-only build's validation step).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import cousin_references as cr

# Single source for the sweep/draw settings used by builds and checkpoints.
# pooled_analysis.py is hash-pinned in cache identities, so these live here
# rather than there; changing them requires a new run directory.
INTERVALS = 300
DRAWS = 10000
SEED = 42
DEFAULT_WORKERS = 24


def sweep_settings():
    return dict(intervals=INTERVALS, draws=DRAWS, seed=SEED)


def positive_workers(value):
    """Validate the shared CLI worker option without affecting cache identity."""
    try:
        workers = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("workers must be a positive integer") from exc
    if workers < 1:
        raise argparse.ArgumentTypeError("workers must be a positive integer")
    return workers


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", type=Path, default=pa.DEFAULT_RUN)
    parser.add_argument("--workers", type=positive_workers, default=DEFAULT_WORKERS,
                        help=f"Parallel solver processes (default: {DEFAULT_WORKERS}; use 1 for serial execution).")
    parser.add_argument("--layout-only", action="store_true",
                        help="Replay-validate existing caches; never solve or draw references")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    ctx, fronts, _, nulls = pa.build(run=args.run, workers=args.workers, layout_only=args.layout_only,
                                     **sweep_settings())
    nulls.update(cr.build(ctx, args.run, draws=DRAWS, seed=SEED, layout_only=args.layout_only))
    print(f"Validated pooled full-tree caches in {args.run}: {len(fronts)} front rows, "
          f"{len(nulls)} reference families. Draw figures with "
          "`python -m publication build --family full-tree-pooled --output-dir PATH`.", flush=True)


if __name__ == "__main__":
    main()
