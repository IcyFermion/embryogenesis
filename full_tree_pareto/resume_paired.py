"""Resume expensive paired matching with verified, per-weight checkpoints.

Run before pooled_pipeline when the paired aggregate cache is missing.
Uses the unchanged production solver and cache identity; never relabels caches.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto.pooled_pipeline import DEFAULT_WORKERS, positive_workers, sweep_settings

SETTINGS = sweep_settings()
INTERVALS = SETTINGS["intervals"]
WEIGHTS = INTERVALS + 1


def stamp(message):
    print(f"{datetime.now(timezone.utc).isoformat()} {message}", flush=True)


def run(run_dir=pa.DEFAULT_RUN, workers=DEFAULT_WORKERS):
    workers = positive_workers(workers)
    ctx = pa.build_context()
    identity = pa.cache_identity(ctx, INTERVALS, SETTINGS["draws"], SETTINGS["seed"])
    target = Path(run_dir) / "analysis" / "paired_bottom-up.npz"
    if target.exists():
        pa.verify_cache(target, identity)
        stamp("Complete paired cache already exists")
        return
    checkpoint = Path(run_dir) / "checkpoints" / "paired_bottom-up"
    checkpoint.mkdir(parents=True, exist_ok=True)
    parents = {}
    for index in range(WEIGHTS):
        path = checkpoint / f"weight_{index:03d}.npz"
        if not path.exists():
            continue
        # A crash between the two atomic renames leaves an uncommitted file;
        # retain it as evidence and recompute only that weight.
        if not path.with_suffix(".json").exists():
            stamp(f"Uncommitted checkpoint {index}; will recompute")
            continue
        pa.verify_cache(path, identity | {"weight_index": index})
        with np.load(path) as cache:
            p = cache["parents"]
            costs = cache["costs"]
        replay = pa.validate_assignment(ctx, p, fixed_roots=False)
        np.testing.assert_allclose(replay, costs, rtol=1e-12, atol=1e-12)
        parents[index] = p
    missing = [i for i in range(WEIGHTS) if i not in parents]
    stamp(f"Verified {len(parents)}/{WEIGHTS} checkpoints; starting {len(missing)} weights with {workers} workers")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(pa._solve, (ctx, pa.METHODS[-1], i, INTERVALS)): i for i in missing}
        for future in as_completed(futures):
            index = futures[future]
            p = future.result()
            costs = pa.validate_assignment(ctx, p, fixed_roots=False)
            final = checkpoint / f"weight_{index:03d}.npz"
            temporary = checkpoint / f"pending_{index:03d}.npz"
            pa.save_cache(temporary, identity | {"weight_index": index}, parents=p, costs=costs)
            temporary.replace(final)
            temporary.with_suffix(".json").replace(final.with_suffix(".json"))
            parents[index] = p
            stamp(f"Saved weight {index}/{INTERVALS}; {len(parents)}/{WEIGHTS} complete")
    ordered = np.array([parents[i] for i in range(WEIGHTS)])
    costs = np.array([pa.validate_assignment(ctx, p, fixed_roots=False) for p in ordered])
    temporary = target.with_name("paired_pending.npz")
    pa.save_cache(temporary, identity, parents=ordered, costs=costs)
    temporary.replace(target)
    temporary.with_suffix(".json").replace(target.with_suffix(".json"))
    stamp(f"All {WEIGHTS} paired weights verified and aggregate cache saved")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=pa.DEFAULT_RUN)
    parser.add_argument("--workers", type=positive_workers, default=DEFAULT_WORKERS,
                        help=f"Parallel solver processes (default: {DEFAULT_WORKERS}).")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    run(args.run, args.workers)
