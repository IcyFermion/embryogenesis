"""Per-round feature ranking for the pooled full tree (Track E, step E1).

Builds on the E0 labels and criterion in ``layer_features.py`` (imported, not
edited, so the E0 run keeps verifying). Settings agreed with the author on
2026-10-15:

* **One k, one method.** In every stratum (rounds 1, 2, 3, 4, 5-8) the top
  ``K = 20`` proteins of the 210, ranked one protein at a time on the E0
  criterion: within each path-length group of the stratum's pairs, the
  Spearman correlation between the absolute difference of the protein's
  z-scored values and the Bray-Curtis distance between descendant mixtures,
  averaged over path lengths with equal weight. Ties break by name.
* **Held-out evaluation.** Ranking and scoring on the same pairs is circular.
  Each stratum is split in two ``N_SPLITS`` times, keeping first-cousin
  families (cells sharing their grandmother in the complete lineage) together
  and balancing cell counts. A set is ranked on one half and scored with the
  E0 criterion (cosine, equal weight per path length) on pairs within the other
  half, in both directions.
* **Reference.** ``N_RANDOM`` random 20-sets from the stratum's own highly
  variable pool (top half of the panel by log1p variance of raw intensities
  over the stratum's cells), scored on the same held-out halves.
* **Gate.** A per-round set passes in a stratum if, on held-out halves, it
  beats the global top-20 in at least ``WIN_FRACTION`` of halves and its median
  percentile among the random sets is at least ``PERCENTILE_MIN``. E2 goes
  ahead only if it passes in each of strata 2, 3 and 4.
* **Shared-set control** (reported, not part of the gate). One set for all
  strata ranked by the same criterion: protein scores averaged over strata
  with equal weight. Split ``i`` side ``A`` trains on the training halves of
  every stratum, so the control is held out exactly like the per-round sets.
  The global top-20 was chosen by a terminal classifier, not by the criterion;
  this control separates "ranked on the criterion" from "ranked per round".
* **Reported beside it.** The per-round set ranked on the whole stratum, its
  overlap with the global top-20, selection frequency over the training halves
  (subsampling stability), its coverage, and its lineage memory next to the
  global top-20's (paired shuffles).

Nothing here uses travel or a front.

Run from the repository root in the dev environment::

    python -m full_tree_pareto.layer_ranking --run-id e1_20261015
    python -m full_tree_pareto.layer_ranking --run-id e1_20261015 --verify-only
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from full_tree_pareto import layer_features as lf
from terminal_pareto import feature_sensitivity as fs

ROOT = lf.ROOT
VERSION = "layer-features-e1-2"
DEFAULT_RUN = "e1_20261015"
RUN_ROOT = lf.RUN_ROOT
CODE_SOURCES = ("full_tree_pareto/layer_ranking.py", "terminal_pareto/feature_sensitivity.py")

SEED = 20261016
K = 20
N_SPLITS = 50          # two-way splits per stratum; each gives two held-out halves
N_RANDOM = 200         # random K-sets from each stratum's HV pool
N_PERM = lf.N_PERM     # majority-class shuffles for lineage memory
HV_FRACTION = fs.HV_FRACTION
WIN_FRACTION = 0.8
PERCENTILE_MIN = 0.9
GATE_STRATA = ("2", "3", "4")


def run_path(run_id):
    if not run_id or Path(run_id).name != run_id or run_id in (".", ".."):
        raise ValueError("Run identity must be a single directory name")
    return RUN_ROOT / run_id


# ---- ranking -------------------------------------------------------------------------

def restrict(stratum, keep):
    """The stratum's pairs with both cells in ``keep`` (a boolean mask over its cells).

    Cell indices are unchanged; path-length groups are re-selected on the kept
    pairs with the E0 rules.
    """
    mask = keep[stratum.a] & keep[stratum.b]
    a, b, path, bc = stratum.a[mask], stratum.b[mask], stratum.path[mask], stratum.bc[mask]
    paths = []
    for value in np.unique(path):
        index = np.flatnonzero(path == value)
        if len(index) >= lf.MIN_PATH_PAIRS and np.ptp(bc[index]) > 0:
            paths.append((int(value), index, rankdata(bc[index])))
    memory = {k: np.flatnonzero(np.isin(np.flatnonzero(mask), v)) for k, v in stratum.memory.items()}
    return replace(stratum, a=a, b=b, path=path, bc=bc, paths=paths, sisters=np.flatnonzero(path == 2), memory=memory)


def protein_scores(Z, stratum):
    """Univariate E0 criterion per protein (columns of ``Z``, rows = stratum cells)."""
    total = np.zeros(Z.shape[1])
    for _, index, ranks in stratum.paths:
        diff = np.abs(Z[stratum.a[index]] - Z[stratum.b[index]])
        x = rankdata(diff, axis=0)
        x -= x.mean(axis=0)
        y = ranks - ranks.mean()
        norm = np.sqrt((x ** 2).sum(axis=0) * (y ** 2).sum())
        total += np.divide(x.T @ y, norm, out=np.zeros(Z.shape[1]), where=norm > 0)
    return total / len(stratum.paths) if stratum.paths else np.full(Z.shape[1], np.nan)


def top_k(scores, features, k=K):
    order = sorted(range(len(features)), key=lambda i: (-scores[i], features[i]))
    return [features[i] for i in order[:k]]


def held_out_score(ctx, stratum, features):
    d = lf.pair_distances(lf.cosine_matrix(lf.matrix(ctx, stratum, features)), stratum)
    return lf.separation(d, stratum)[0]["separation"]


# ---- splits ----------------------------------------------------------------------------

def families(ctx, stratum):
    """First-cousin family of each stratum cell: its grandmother in the complete lineage."""
    keys = [(ctx.lineage.ancestors[c][-2] if len(ctx.lineage.ancestors[c]) >= 2 else c) for c in stratum.cells]
    _, labels = np.unique(keys, return_inverse=True)
    return labels


def split_halves(labels, rng):
    """Families in random order, each to the half with fewer cells so far."""
    sizes = np.bincount(labels)
    side = np.zeros(len(sizes), dtype=bool)
    count = [0, 0]
    for f in rng.permutation(len(sizes)):
        side[f] = count[1] < count[0]
        count[int(side[f])] += sizes[f]
    return side[labels]


def make_halves(ctx, stratum, n_splits, seed):
    """(split, side, train mask, train stratum, test stratum) for both directions of every split."""
    labels = families(ctx, stratum)
    rng = np.random.default_rng(seed)
    halves = []
    for split in range(n_splits):
        in_b = split_halves(labels, rng)
        for side, test in (("A", in_b), ("B", ~in_b)):
            halves.append((split, side, ~test, restrict(stratum, ~test), restrict(stratum, test)))
    return halves


# ---- context ---------------------------------------------------------------------------

def hv_pool(ctx, stratum):
    table = fs.variability_scores({"s3": ctx.raw}, list(ctx.zscore.columns), stratum.cells)
    if not np.isclose(HV_FRACTION * len(table), round(HV_FRACTION * len(table))):
        raise ValueError("HV pool size is not an integer")
    return table


def load(*, n_splits=N_SPLITS, seed=SEED):
    ctx = lf.load_context()
    report = lf.check_reproduction(ctx)
    if report["failures"]:
        raise ValueError("E0 reproduction gates failed:\n  " + "\n  ".join(report["failures"]))
    pools = {name: hv_pool(ctx, s) for name, s in ctx.strata.items()}
    halves = {name: make_halves(ctx, s, n_splits, [seed, 1, i]) for i, (name, s) in enumerate(ctx.strata.items())}
    sources = dict(ctx.sources)
    sources.update({name: lf.file_hash(ROOT / name) for name in CODE_SOURCES})
    return ctx, pools, halves, report, sources


def draw_random_sets(ctx, pools, *, n=N_RANDOM, k=K, seed=SEED, max_rejections=10_000):
    rows = []
    for i, (name, stratum) in enumerate(ctx.strata.items()):
        pool = pools[name].loc[pools[name].in_hv_pool, "feature"].tolist()
        rng = np.random.default_rng([seed, 99, i, k])
        for draw in range(n):
            for rejected in range(max_rejections + 1):
                chosen = [pool[j] for j in rng.choice(len(pool), size=k, replace=False)]
                if (np.linalg.norm(lf.matrix(ctx, stratum, chosen), axis=1) > 0).all():
                    break
            else:
                raise ValueError(f"Could not draw a valid set for stratum {name}")
            rows.append(dict(stratum=name, draw=draw, rejected_draws=rejected, features=";".join(chosen)))
    return pd.DataFrame(rows)


# ---- computation -------------------------------------------------------------------------

def _score_job(ctx, name, halves):
    """Protein scores on each training half."""
    Z = lf.matrix(ctx, ctx.strata[name], list(ctx.zscore.columns))
    return [(name, split, side, protein_scores(Z, train)) for split, side, _, train, _ in halves]


def _evaluate_job(ctx, name, halves, sets):
    top20 = ctx.candidates["top-20"]
    rows = []
    for split, side, train_mask, _, test in halves:
        selected, shared = sets[name, split, side], sets["shared", split, side]
        rows.append(dict(stratum=name, split=split, side=side, train_cells=int(train_mask.sum()),
                         test_cells=int((~train_mask).sum()), test_pairs=len(test.a),
                         test_path_lengths=len(test.paths), per_round=held_out_score(ctx, test, selected),
                         top20=held_out_score(ctx, test, top20), shared=held_out_score(ctx, test, shared),
                         overlap_top20=len(set(selected) & set(top20)), overlap_shared=len(set(selected) & set(shared))))
    return rows


def _random_job(ctx, name, halves, draws):
    rows = []
    for draw, features in draws:
        D = lf.cosine_matrix(lf.matrix(ctx, ctx.strata[name], features.split(";")))
        for split, side, _, _, test in halves:
            rows.append(dict(stratum=name, draw=draw, split=split, side=side,
                             separation=lf.separation(lf.pair_distances(D, test), test)[0]["separation"]))
    return rows


def _memory_job(ctx, name, features_by_set, seed, n_perm):
    stratum = ctx.strata[name]
    rows = []
    for label, features in features_by_set.items():
        # Same seed for every set: shuffles are paired across sets.
        rng = np.random.default_rng(seed)
        for m in lf.lineage_memory(lf.matrix(ctx, stratum, features), stratum, rng, n_perm):
            rows.append(dict(stratum=name, feature_set=label, **m))
    return rows


def compute(ctx, pools, halves, *, n_jobs, n_random=N_RANDOM, n_perm=N_PERM, seed=SEED):
    from joblib import Parallel, delayed

    features = list(ctx.zscore.columns)
    top20 = ctx.candidates["top-20"]
    # Full-stratum ranking: the per-round sets themselves, and the shared control.
    ranking, sets = [], {}
    full_scores = {name: protein_scores(lf.matrix(ctx, s, features), s) for name, s in ctx.strata.items()}
    shared_scores = np.mean(list(full_scores.values()), axis=0)
    sets["shared"] = top_k(shared_scores, features)
    ranking += [dict(stratum="shared", protein=f, score=float(s), in_global_top20=f in top20)
                for f, s in zip(features, shared_scores)]
    for name, stratum in ctx.strata.items():
        scores = full_scores[name]
        sets[name] = top_k(scores, features)
        pool = set(pools[name].loc[pools[name].in_hv_pool, "feature"])
        ranking += [dict(stratum=name, protein=f, score=float(s), in_global_top20=f in top20,
                         in_stratum_hv_pool=f in pool) for f, s in zip(features, scores)]
    ranking = pd.DataFrame(ranking)
    ranking = ranking.sort_values(["stratum", "score", "protein"], ascending=[True, False, True]).reset_index(drop=True)
    ranking["rank"] = ranking.groupby("stratum").cumcount() + 1
    ranking["selected"] = ranking["rank"] <= K

    chunks = [(name, h[i:i + 10]) for name, h in halves.items() for i in range(0, len(h), 10)]
    trained = [t for ts in Parallel(n_jobs=n_jobs)(delayed(_score_job)(ctx, name, h) for name, h in chunks) for t in ts]
    half_sets = {(name, split, side): top_k(scores, features) for name, split, side, scores in trained}
    by_half = {}
    for name, split, side, scores in trained:
        by_half.setdefault((split, side), []).append(scores)
    if any(len(v) != len(ctx.strata) for v in by_half.values()):
        raise ValueError("Every split must train in every stratum")
    half_sets.update({("shared", split, side): top_k(np.mean(v, axis=0), features) for (split, side), v in by_half.items()})
    heldout = pd.DataFrame([r for rows in Parallel(n_jobs=n_jobs)(
        delayed(_evaluate_job)(ctx, name, h, half_sets) for name, h in chunks) for r in rows])
    chosen = pd.DataFrame([dict(stratum=name, split=split, side=side, features=";".join(f))
                           for (name, split, side), f in half_sets.items()])

    random_sets = draw_random_sets(ctx, pools, n=n_random, seed=seed)
    work = [(name, list(g[["draw", "features"]].itertuples(index=False, name=None))) for name, g in random_sets.groupby("stratum", sort=False)]
    jobs = [(name, draws[i:i + 10]) for name, draws in work for i in range(0, len(draws), 10)]
    random_heldout = pd.DataFrame([r for rows in Parallel(n_jobs=n_jobs)(
        delayed(_random_job)(ctx, name, halves[name], d) for name, d in jobs) for r in rows])
    heldout = add_percentiles(heldout, random_heldout)

    memory = pd.DataFrame([r for rows in Parallel(n_jobs=min(n_jobs, len(sets)))(
        delayed(_memory_job)(ctx, name, {"per-round": sets[name], "shared": sets["shared"], "top-20": top20},
                             [seed, 7, i], n_perm)
        for i, name in enumerate(ctx.strata)) for r in rows])
    stability = selection_stability(ranking, chosen, sets)
    return dict(ranking=ranking.merge(stability, on=["stratum", "protein"]), heldout=heldout, half_sets=chosen,
                random_sets=random_sets, random_heldout=random_heldout, memory=memory,
                sets=set_table(ctx, sets, chosen, heldout, memory), gate=gate_table(heldout))


def add_percentiles(heldout, random_heldout):
    """Percentile of the per-round set and of the global top-20 among random sets on the same held-out half."""
    keyed = random_heldout.groupby(["stratum", "split", "side"]).separation.apply(np.asarray)
    out = heldout.copy()
    for column in ("per_round", "top20", "shared"):
        out[f"{column}_percentile"] = [float(np.mean(keyed[(r.stratum, r.split, r.side)] < r[column])
                                             + .5 * np.mean(keyed[(r.stratum, r.split, r.side)] == r[column]))
                                       for _, r in out.iterrows()]
    out["random_median"] = [float(np.median(keyed[(r.stratum, r.split, r.side)])) for _, r in out.iterrows()]
    out["win"] = out.per_round > out.top20
    out["win_shared"] = out.per_round > out.shared
    return out


def selection_stability(ranking, chosen, sets):
    rows = []
    for name, group in chosen.groupby("stratum", sort=False):
        counts = pd.Series([f for s in group.features for f in s.split(";")]).value_counts()
        for protein in ranking.loc[ranking.stratum == name, "protein"]:
            rows.append(dict(stratum=name, protein=protein, selection_frequency=float(counts.get(protein, 0) / len(group))))
    return pd.DataFrame(rows)


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b)


def set_table(ctx, sets, chosen, heldout, memory):
    top20 = ctx.candidates["top-20"]
    rows = []
    for name, features in ((n, f) for n, f in sets.items() if n != "shared"):
        stratum = ctx.strata[name]
        expressed = (ctx.raw.loc[stratum.cells, features] > 0).mean()
        in_sample = lf.separation(lf.pair_distances(lf.cosine_matrix(lf.matrix(ctx, stratum, features)), stratum), stratum)[0]
        halves = chosen[chosen.stratum == name].features.str.split(";")
        h = heldout[heldout.stratum == name]
        mem = memory[(memory.stratum == name) & (memory.memory_pairs == "near-identical")].set_index("feature_set")
        rows.append(dict(stratum=name, features=";".join(features), overlap_top20=len(set(features) & set(top20)),
                         jaccard_top20=jaccard(features, top20),
                         mean_jaccard_half_sets=float(np.mean([jaccard(features, s) for s in halves])),
                         variable=int(((expressed >= lf.COVERAGE_MIN) & (expressed <= lf.COVERAGE_MAX)).sum()),
                         in_sample_separation=in_sample["separation"], in_sample_sisters=in_sample["sisters"],
                         overlap_shared=len(set(features) & set(sets["shared"])),
                         heldout_per_round_median=float(h.per_round.median()), heldout_top20_median=float(h.top20.median()),
                         heldout_shared_median=float(h.shared.median()),
                         memory_excess=float(mem.at["per-round", "excess"]), memory_excess_shared=float(mem.at["shared", "excess"]),
                         memory_excess_top20=float(mem.at["top-20", "excess"])))
    shared = sets["shared"]
    rows.append(dict(stratum="shared", features=";".join(shared), overlap_top20=len(set(shared) & set(top20)),
                     jaccard_top20=jaccard(shared, top20),
                     mean_jaccard_half_sets=float(np.mean([jaccard(shared, s) for s in
                                                           chosen[chosen.stratum == "shared"].features.str.split(";")]))))
    return pd.DataFrame(rows)


def gate_table(heldout):
    rows = []
    for name, h in heldout.groupby("stratum", sort=False):
        win = float(h.win.mean())
        pct = float(h.per_round_percentile.median())
        rows.append(dict(stratum=name, halves=len(h), win_fraction=win, median_percentile=pct,
                         top20_median_percentile=float(h.top20_percentile.median()),
                         median_gain=float((h.per_round - h.top20).median()),
                         shared_median_percentile=float(h.shared_percentile.median()),
                         win_fraction_vs_shared=float(h.win_shared.mean()),
                         median_gain_vs_shared=float((h.per_round - h.shared).median()),
                         passes=win >= WIN_FRACTION and pct >= PERCENTILE_MIN, gate_stratum=name in GATE_STRATA))
    return pd.DataFrame(rows)


def decision(gate):
    passed = gate.set_index("stratum").passes
    proceed = bool(all(passed[s] for s in GATE_STRATA))
    return dict(proceed_to_e2=proceed, passing_strata=[s for s in passed.index if passed[s]],
                rule=f"per-round set passes in each of strata {', '.join(GATE_STRATA)}: held-out win fraction >= "
                     f"{WIN_FRACTION} against the global top-20 and median random percentile >= {PERCENTILE_MIN}")


def settings(*, n_splits=N_SPLITS, n_random=N_RANDOM, n_perm=N_PERM, seed=SEED):
    return dict(seed=seed, k=K, n_splits=n_splits, n_random=n_random, n_perm=n_perm, hv_fraction=HV_FRACTION,
                ranking="univariate_e0_criterion_abs_difference_equal_weight_per_path_length",
                split="two_way_first_cousin_families_balanced_cells", reference="random_k_sets_from_stratum_hv_pool",
                win_fraction=WIN_FRACTION, percentile_min=PERCENTILE_MIN, gate_strata=list(GATE_STRATA),
                e0=lf.settings(), numpy_version=np.__version__)


TABLES = ("ranking", "sets", "gate", "heldout", "half_sets", "random_sets", "random_heldout", "memory")


def write_run(out, *, n_jobs, **overrides):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    run_settings = settings(**overrides)
    ctx, pools, halves, report, sources = load(n_splits=run_settings["n_splits"], seed=run_settings["seed"])
    print("E0 gates reproduced.", flush=True)
    results = compute(ctx, pools, halves, n_jobs=n_jobs, n_random=run_settings["n_random"],
                      n_perm=run_settings["n_perm"], seed=run_settings["seed"])
    analysis.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        results[name].to_csv(analysis / f"{name}.csv", index=False)
    pd.concat([p.assign(stratum=name) for name, p in pools.items()]).to_csv(analysis / "hv_pools.csv", index=False)
    lf.pa.write_json(analysis / "decision.json", decision(results["gate"]))
    (analysis / "reproduction.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=run_settings, sources=sources),
                                         sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=run_settings, sources=sources,
                  elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: lf.file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def validate_run(out, *, n_jobs=1, check_sources=True):
    """Check file and source hashes, then recompute every table and compare at 1e-12."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported layer-ranking cache")
    for name, expected in record["files"].items():
        if lf.file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if lf.file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    s = record["settings"]
    ctx, pools, halves, _, _ = load(n_splits=s["n_splits"], seed=s["seed"])
    fresh = compute(ctx, pools, halves, n_jobs=n_jobs, n_random=s["n_random"], n_perm=s["n_perm"], seed=s["seed"])
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        lf._compare(fresh[name], saved, name)
    if json.loads((analysis / "decision.json").read_text()) != decision(fresh["gate"]):
        raise ValueError("Replay mismatch decision")
    return dict(analysis_id=record["analysis_id"], tables=len(TABLES), decision=decision(fresh["gate"]),
                verified_source_hashes=len(record["sources"]) if check_sources else 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--jobs", type=int, default=min(24, os.cpu_count() or 1))
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = run_path(args.run_id)
    if not args.verify_only:
        record = write_run(out, n_jobs=args.jobs)
        print(json.dumps({k: record[k] for k in ("analysis_id", "elapsed_seconds")}, indent=2))
    print(json.dumps(validate_run(out, n_jobs=args.jobs), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
