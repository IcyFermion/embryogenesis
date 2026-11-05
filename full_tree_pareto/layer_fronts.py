"""Round-by-round layerwise fronts under the E1 feature sets (Track E, E1 front check).

Feature selection in E0/E1 never looks at a front. This add-on asks the
follow-up question afterwards: given the sets E1 chose, how do the per-round
Pareto fronts compare with those under the global top-20?

Each configuration swaps only the cell-state matrix of the pooled full tree
(``pooled_analysis.build_context``); travel, the cohort, the rounds and the
sweep are unchanged:

* ``top-20``: the published global set; must reproduce the cached layerwise
  round fronts and the per-round canonical metrics of ``round_references``.
* ``per-round``: round r is solved and scored with the E1 set of its stratum
  (rounds 5-8 share the 5-8 set).
* ``shared``: the E1 shared-set control, one set for every round.

For each configuration and round: cosine distance between z-scored profiles
of the set's proteins (NaN -> 0, as in the pooled context); optimizer scaling
by that matrix's global leaf first-cousin null SDs (the pooled rule); a
layerwise sweep of ``INTERVALS + 1`` weights solving only that round's
assignment; the round's exact first-cousin null mean; and the canonical
metrics of ``round_references.landmarks`` in the round's own endpoint
coordinates. The cached round-reference first-cousin permutations are
re-scored for display clouds and for the cell-state z of the natural lineage.

Run from the repository root in the dev environment::

    python -m full_tree_pareto.layer_fronts --run-id e1_fronts_20261015
    python -m full_tree_pareto.layer_fronts --run-id e1_fronts_20261015 --verify-only
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from full_tree_pareto import layer_features as lf
from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import round_references as rr
from terminal_pareto.subtree_analysis import exact_cousin_stats

ROOT = lf.ROOT
VERSION = "layer-features-e1-fronts-1"
DEFAULT_RUN = "e1_fronts_20261015"
E1_RUN = lf.RUN_ROOT / "e1_20261015"
INTERVALS = 300
CONFIGS = ("top-20", "per-round", "shared")
CODE_SOURCES = ("full_tree_pareto/layer_fronts.py", "full_tree_pareto/layer_features.py",
                "full_tree_pareto/pooled_analysis.py", "full_tree_pareto/round_references.py",
                "terminal_pareto/front_coordinates.py", "terminal_pareto/subtree_analysis.py")
DATA_SOURCES = (lf.PROTEIN_ZSCORE,)
TABLES = ("fronts", "metrics", "clouds", "sets")
CLOUD_DISPLAY = 300        # first-cousin draws per round and configuration kept for display


def stratum_of(round_number):
    return next(name for name, members in lf.STRATA.items() if round_number in members)


def state_matrix(zscore, names, features):
    X = zscore.loc[names, list(features)].to_numpy(float)
    if not (np.linalg.norm(X, axis=1) > 0).all():
        raise ValueError("Cosine cell-state cost needs a nonzero profile for every cell")
    return cdist(X, X, metric="cosine")


def leaf_scales(ctx, state):
    """Global leaf first-cousin null SDs, as ``pa.build_context`` computes them for the top-20."""
    leaves = ctx.leaves
    groups = [g for g in ctx.cousin_groups if len(g) > 1]
    stats = exact_cousin_stats(ctx.travel[np.ix_(leaves, ctx.parents[leaves])],
                               state[np.ix_(leaves, ctx.parents[leaves])], groups)
    return np.sqrt(np.asarray(stats[2:4], dtype=float))


def solve_round(ctx, state, scales, children, intervals=INTERVALS):
    """Layerwise sweep of one round: (intervals+1, len(children)) parent indices."""
    slots = ctx.parents[children]
    travel, cell = ctx.travel[np.ix_(children, slots)] / scales[0], state[np.ix_(children, slots)] / scales[1]
    sweep = np.empty((intervals + 1, len(children)), dtype=int)
    for i in range(intervals + 1):
        alpha = i / intervals
        r, c = linear_sum_assignment(alpha * travel + (1 - alpha) * cell)
        sweep[i, r] = slots[c]
    return sweep


def load_sets():
    """E1 sets, after checking the file against the E1 run's recorded hash."""
    record = json.loads((E1_RUN / "analysis/provenance.json").read_text())
    path = E1_RUN / "analysis/sets.csv"
    if lf.file_hash(path) != record["files"]["sets.csv"]:
        raise ValueError("E1 sets.csv does not match its provenance")
    sets = pd.read_csv(path).set_index("stratum").features.str.split(";")
    return {name: list(features) for name, features in sets.items()}, record["analysis_id"]


def load():
    ctx = pa.build_context()
    zscore = pd.read_csv(ROOT / lf.PROTEIN_ZSCORE, index_col=0).T.fillna(0)
    e1_sets, e1_id = load_sets()
    references = rr.load(pa.DEFAULT_RUN, "pooled")
    if references["identity"]["context_key"] != ctx.identity["context_key"]:
        raise ValueError("Round references were drawn on a different pooled context")
    with np.load(pa.DEFAULT_RUN / rr.SUBDIR / "draws.npz") as data:
        permutations = {r: data[f"round_{r}__first_cousin"] for r in range(1, len(ctx.layers) + 1)}
    sets = {("top-20", r): ctx.identity["proteins"] for r in range(1, len(ctx.layers) + 1)}
    sets.update({("per-round", r): e1_sets[stratum_of(r)] for r in range(1, len(ctx.layers) + 1)})
    sets.update({("shared", r): e1_sets["shared"] for r in range(1, len(ctx.layers) + 1)})
    return ctx, zscore, sets, permutations, references, dict(e1_analysis_id=e1_id,
                                                             round_reference_identity=references["identity_sha"])


def source_hashes():
    paths = (*CODE_SOURCES, *DATA_SOURCES, str((E1_RUN / "analysis/sets.csv").relative_to(ROOT)),
             str((pa.DEFAULT_RUN / "analysis/context.json").relative_to(ROOT)))
    return {name: lf.file_hash(ROOT / name) for name in paths}


def compute(ctx, zscore, sets, permutations, *, intervals=INTERVALS):
    groups_index = rr.tree_index()
    matrices, fronts, metrics, clouds = {}, [], [], []
    for (config, r), features in sets.items():
        key = tuple(features)
        if key not in matrices:
            state = state_matrix(zscore, ctx.names, features)
            matrices[key] = state, leaf_scales(ctx, state)
        state, scales = matrices[key]
        children = ctx.layers[r - 1]
        sweep = solve_round(ctx, state, scales, children, intervals)
        blocks = rr.round_blocks(ctx.travel, state, ctx.parents, children)
        natural = np.array([np.trace(blocks[0]), np.trace(blocks[1])])
        groups = rr.round_groups(ctx.names, children, groups_index, "first_cousin")
        stats = rr.exact_moments(blocks, groups) if groups else np.r_[natural, 0., 0., 0.]
        costs = np.column_stack((ctx.travel[children[None, :], sweep].sum(axis=1),
                                 state[children[None, :], sweep].sum(axis=1)))
        retention = (sweep == ctx.parents[children]).mean(axis=1)
        frame, metric, transform = rr.landmarks(costs, natural, stats[:2], retention,
                                                f"{VERSION}:{config}:round_{r}")
        cloud = rr.score_permutations(blocks, permutations[r])
        sd_state = np.sqrt(stats[3])
        metrics.append(dict(config=config, round=r, stratum=stratum_of(r), edges=len(children),
                            scale_travel=scales[0], scale_state=scales[1],
                            natural_travel=natural[0], natural_state=natural[1],
                            null_mean_travel=stats[0], null_mean_state=stats[1],
                            null_sd_travel=np.sqrt(stats[2]), null_sd_state=sd_state,
                            state_advantage_z=(stats[1] - natural[1]) / sd_state if sd_state > 0 else np.nan,
                            **metric))
        if frame is None:
            raise ValueError(f"Degenerate endpoint span: {config} round {r}")
        fronts.append(frame.assign(config=config, round=r))
        shown = cloud[:CLOUD_DISPLAY]
        D1, D2 = transform.transform(shown[:, 0], shown[:, 1])
        clouds.append(pd.DataFrame(dict(config=config, round=r, draw=np.arange(len(shown)), travel=shown[:, 0],
                                        cell_state=shown[:, 1], D1=D1, D2=D2)))
    set_rows = pd.DataFrame([dict(config=c, round=r, features=";".join(f)) for (c, r), f in sets.items()])
    return dict(fronts=pd.concat(fronts, ignore_index=True), metrics=pd.DataFrame(metrics),
                clouds=pd.concat(clouds, ignore_index=True), sets=set_rows), matrices


def check_reproduction(ctx, results, matrices, references, sets):
    """The top-20 configuration must reproduce the pooled state, cached fronts and round metrics."""
    state, scales = matrices[tuple(sets["top-20", 1])]
    np.testing.assert_allclose(state, ctx.state, rtol=0, atol=1e-12)
    np.testing.assert_allclose(scales, ctx.scales, rtol=1e-12)
    cached = pd.read_csv(pa.DEFAULT_RUN / "analysis/fronts.csv")
    cached = cached[cached.method == pa.METHODS[0]]
    mine = results["fronts"][results["fronts"].config == "top-20"]
    rounds = sorted(mine["round"].unique())
    for r in rounds:
        saved = cached[cached.scope == f"round_{r}"].sort_values("weight_index")
        np.testing.assert_allclose(mine[mine["round"] == r][["travel", "cell_state", "edge_retention"]].to_numpy(),
                                   saved[["travel", "state", "retention"]].to_numpy(), rtol=1e-12, atol=1e-12)
    columns = ["natural_D1", "natural_D2", "null_D1", "null_D2", "u_L", "d_LP", "d_NP", "max_edge_retention",
               "retention_at_nearest", "nearest_index", "maximum_index"]
    metrics = results["metrics"].query("config == 'top-20'").set_index("round")[columns]
    saved = references["metrics"].set_index("round")[columns]
    np.testing.assert_allclose(metrics.to_numpy(float), saved.loc[metrics.index].to_numpy(float),
                               rtol=1e-10, atol=1e-12)
    return dict(rounds=len(rounds), fronts_points=int(len(mine)), metrics=columns,
                state_matrix="identical to pooled context", scales=ctx.scales.tolist())


def settings(*, intervals=INTERVALS):
    return dict(intervals=intervals, configs=list(CONFIGS), cloud_display=CLOUD_DISPLAY,
                scaling="global leaf first-cousin null SDs of each configuration's state matrix",
                null="exact first-cousin moments per round; cached round-reference permutations for clouds",
                numpy_version=np.__version__)


def run_path(run_id):
    return lf.run_path(run_id)


def write_run(out, *, intervals=INTERVALS):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    run_settings = settings(intervals=intervals)
    ctx, zscore, sets, permutations, references, upstream = load()
    sources = source_hashes()
    results, matrices = compute(ctx, zscore, sets, permutations, intervals=intervals)
    report = check_reproduction(ctx, results, matrices, references, sets)
    print("Top-20 reproduces the cached round fronts and metrics.", flush=True)
    analysis.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        results[name].to_csv(analysis / f"{name}.csv", index=False)
    (analysis / "reproduction.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=run_settings, upstream=upstream,
                                              context_key=ctx.identity["context_key"], sources=sources),
                                         sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=run_settings, upstream=upstream,
                  context_key=ctx.identity["context_key"], sources=sources,
                  elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: lf.file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def validate_run(out, *, check_sources=True):
    """Check file and source hashes, then re-solve every sweep and compare at 1e-12."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported layer-fronts cache")
    for name, expected in record["files"].items():
        if lf.file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if lf.file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    ctx, zscore, sets, permutations, references, upstream = load()
    if upstream != record["upstream"] or ctx.identity["context_key"] != record["context_key"]:
        raise ValueError("Upstream E1 run or pooled context changed")
    fresh, matrices = compute(ctx, zscore, sets, permutations, intervals=record["settings"]["intervals"])
    check_reproduction(ctx, fresh, matrices, references, sets)
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        lf._compare(fresh[name], saved, name)
    return dict(analysis_id=record["analysis_id"], tables=len(TABLES),
                verified_source_hashes=len(record["sources"]) if check_sources else 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = run_path(args.run_id)
    if not args.verify_only:
        record = write_run(out)
        print(json.dumps({k: record[k] for k in ("analysis_id", "elapsed_seconds")}, indent=2))
    print(json.dumps(validate_run(out), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
