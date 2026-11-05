"""Per-layer lineage advantage over cousin shuffles, with bootstrap intervals (Track E, E3).

Each contraction round of the full tree is its own assignment problem. Per
round, the natural lineage's travel and cell-state costs are compared with the
exact mean of a within-group shuffle (first or second cousins):

    advantage = 1 - natural cost / shuffle mean

Intervals come from resampling the round's shuffle groups with replacement
(children outside every group are singleton units whose natural cost equals
their shuffle mean). Travel and cell state are resampled together, and the
three partial-forest configurations share resample indices, since they share
the forest and its groups.

Text and figures number layers top-down: layer = 9 - round. The contrast
below was fixed in ``plans/track_e_layer_features.md`` before any interval was
computed:

* mid = mean advantage of layers 6 and 7 (rounds 3 and 2);
* late flank = layer 8 (round 1); early flank = mean of layers 4 and 5
  (rounds 5 and 4);
* dip = mid - flank, for each flank.

The dip is claimed only if both 95% intervals exclude zero for pooled
C. elegans travel against first cousins. Everything else is sensitivity.

Inputs are the hash-verified caches of the pooled cosine run and the
cross-species partial forest (raw3d geometry); nothing is solved. With no
resampling, every per-round natural cost and first-cousin mean equals the
cached ``round_references`` metrics.

Run from the repository root in the dev environment::

    python -m full_tree_pareto.layer_hourglass --run-id hourglass_20261019
    python -m full_tree_pareto.layer_hourglass --run-id hourglass_20261019 --verify-only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time

import numpy as np
import pandas as pd

from full_tree_pareto import cross_species_analysis as fcs
from full_tree_pareto import layer_features as lf
from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import round_references as rr

ROOT = lf.ROOT
VERSION = "layer-hourglass-e3-1"
DEFAULT_RUN = "hourglass_20261019"
CROSS_RUN = fcs.run_path(fcs.DEFAULT_RUN)
GEOMETRY = "raw3d"
FAMILIES = ("first_cousin", "second_cousin")
DRAWS = 2000
SEED = 4300
LEVEL = 0.95
TOP = 9                                   # layer = TOP - round
MID, LATE, EARLY = (6, 7), (8,), (4, 5)   # layers
PRIMARY = dict(series="pooled_ce_protein", family="first_cousin", cost="travel")
CODE_SOURCES = ("full_tree_pareto/layer_hourglass.py", "full_tree_pareto/round_references.py",
                "full_tree_pareto/pooled_analysis.py", "full_tree_pareto/cross_species_analysis.py",
                "terminal_pareto/subtree_analysis.py", "terminal_pareto/pareto_engine.py")
TABLES = ("layers", "contrasts")


def layer_of(round_number):
    return TOP - round_number


def units(names, children, index, family):
    """Shuffle groups as positions within ``children``, plus singletons for unshuffled children."""
    groups = rr.round_groups(names, children, index, family)
    covered = np.zeros(len(children), dtype=bool)
    for g in groups:
        covered[g] = True
    return groups + [np.array([i]) for i in np.flatnonzero(~covered)]


def unit_costs(blocks, unit_list):
    """(units, 4): natural travel, natural state, shuffle-mean travel, shuffle-mean state."""
    rows = []
    for g in unit_list:
        travel, state = (b[np.ix_(g, g)] for b in blocks)
        rows.append([np.trace(travel), np.trace(state), travel.mean(axis=1).sum(), state.mean(axis=1).sum()])
    return np.array(rows, dtype=float)


def advantage(costs):
    """1 - natural / shuffle mean for travel and state; ``costs`` is (..., 4) summed over units."""
    return 1 - costs[..., :2] / costs[..., 2:]


def resample_indices(n_units, draws, rng):
    return rng.integers(0, n_units, size=(draws, n_units))


def bootstrap(table, indices):
    """(draws, 2) advantages from resampled unit sums."""
    return advantage(table[indices].sum(axis=1))


def interval(values, level=LEVEL):
    tail = (1 - level) / 2 * 100
    return np.percentile(values, [tail, 100 - tail], axis=0)


# ---- inputs -------------------------------------------------------------------------

def load_pooled():
    ctx = pa.build_context()
    references = rr.load(pa.DEFAULT_RUN, "pooled")
    if references["identity"]["context_key"] != ctx.identity["context_key"]:
        raise ValueError("Round references were drawn on a different pooled context")
    forests = {"pooled_ce_protein": dict(travel=ctx.travel, state=ctx.state, config="ce_protein")}
    return dict(names=ctx.names, parents=ctx.parents, layers=ctx.layers, forests=forests,
                references=references["metrics"], identity=references["identity_sha"],
                context_key=ctx.identity["context_key"])


def load_partial():
    report = fcs.validate_run(CROSS_RUN)
    record = json.loads((CROSS_RUN / "analysis/provenance.json").read_text())
    references = rr.load(CROSS_RUN, "cross_species")
    if references["identity"]["analysis_id"] != report["analysis_id"]:
        raise ValueError("Round references were drawn on a different cross-species analysis")
    with np.load(CROSS_RUN / "analysis/matrices.npz", allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    names, forests, layers = record["names"], {}, None
    for config in fcs.PRIMARY:
        key = f"{GEOMETRY}__{config}"
        ctx = fcs.forest(names, arrays["parents"], arrays[f"{key}__travel"], arrays[f"{key}__state"])
        layers = ctx.layers
        forests[f"partial_{config}"] = dict(travel=ctx.travel, state=ctx.state, config=config)
        species = "cb" if config == "cb_rna" else "ce"
        if config == "ce_rna":
            continue    # same embryos and travel as ce_protein
        for replicate in record["tracking"][species]:
            label = replicate["label"]
            forests[f"embryo_{species}_{label}"] = dict(travel=arrays[f"{GEOMETRY}__{species}__{label}"],
                                                        state=ctx.state, config=config)
    metrics = references["metrics"]
    return dict(names=names, parents=np.asarray(arrays["parents"], dtype=int), layers=layers, forests=forests,
                references=metrics[metrics.geometry == GEOMETRY], identity=references["identity_sha"],
                analysis_id=report["analysis_id"])


# ---- computation --------------------------------------------------------------------

def analyze_tree(tree, index, *, draws, seed, tree_seed):
    """Per-round point estimates and bootstrap draws for every forest of one tree."""
    rows, samples = [], {}
    for r, children in enumerate(tree["layers"], 1):
        for f, family in enumerate(FAMILIES):
            unit_list = units(tree["names"], children, index, family)
            rng = np.random.default_rng([seed, tree_seed, r, f])
            indices = resample_indices(len(unit_list), draws, rng)
            for series, forest in tree["forests"].items():
                blocks = rr.round_blocks(forest["travel"], forest["state"], tree["parents"], children)
                table = unit_costs(blocks, unit_list)
                point = advantage(table.sum(axis=0))
                draws_ = bootstrap(table, indices)
                samples[series, family, r] = draws_
                (t_lo, s_lo), (t_hi, s_hi) = interval(draws_)
                total = table.sum(axis=0)
                rows.append(dict(series=series, family=family, round=r, layer=layer_of(r), edges=len(children),
                                 units=len(unit_list), shuffled=int(sum(len(g) for g in unit_list if len(g) > 1)),
                                 natural_travel=total[0], natural_state=total[1],
                                 shuffle_travel=total[2], shuffle_state=total[3],
                                 travel_advantage=point[0], travel_low=t_lo, travel_high=t_hi,
                                 state_advantage=point[1], state_low=s_lo, state_high=s_hi))
    return rows, samples


def contrasts(samples, layers_table):
    """Pre-specified dips (mid - flank) per series, family and cost."""
    rows = []
    rounds = lambda group: [TOP - layer for layer in group]
    for (series, family), frame in layers_table.groupby(["series", "family"], sort=False):
        present = set(frame["round"])
        if not all(r in present for r in rounds(MID + LATE + EARLY)):
            continue
        point = frame.set_index("round")
        side = lambda group, col: np.mean([point.loc[r, col] for r in rounds(group)])
        drawn = lambda group: np.mean([samples[series, family, r] for r in rounds(group)], axis=0)
        mid = drawn(MID)
        for flank_name, flank in (("late", LATE), ("early", EARLY)):
            diff = mid - drawn(flank)
            low, high = interval(diff)
            for c, cost in enumerate(("travel", "state")):
                estimate = side(MID, f"{cost}_advantage") - side(flank, f"{cost}_advantage")
                rows.append(dict(series=series, family=family, cost=cost, flank=flank_name,
                                 mid_layers="6-7", flank_layers="-".join(map(str, sorted(flank))),
                                 mid=side(MID, f"{cost}_advantage"), flank_value=side(flank, f"{cost}_advantage"),
                                 dip=estimate, low=low[c], high=high[c], excludes_zero=bool(high[c] < 0 or low[c] > 0)))
    return pd.DataFrame(rows)


def decision(contrast_table):
    mask = np.all([contrast_table[k] == v for k, v in PRIMARY.items()], axis=0)
    primary = contrast_table[mask]
    passed = bool(len(primary) == 2 and (primary.high < 0).all())
    return dict(rule="both 95% intervals of mid - flank below zero (pooled C. elegans travel, first cousins)",
                passed=passed, dips={row.flank: [row.dip, row.low, row.high] for row in primary.itertuples()})


def check_reproduction(layers_table, trees):
    """Natural costs and first-cousin means must equal the cached round references."""
    checked = 0
    for tree, prefix in ((trees["pooled"], "pooled_"), (trees["partial"], "partial_")):
        refs = tree["references"]
        mine = layers_table[layers_table.series.str.startswith(prefix) & (layers_table.family == "first_cousin")]
        for row in mine.itertuples():
            ref = refs[(refs.config == tree["forests"][row.series]["config"]) & (refs["round"] == row.round)]
            if len(ref) != 1:
                raise ValueError(f"No cached metrics for {row.series} round {row.round}")
            ref = ref.iloc[0]
            np.testing.assert_allclose([row.natural_travel, row.natural_state, row.shuffle_travel, row.shuffle_state],
                                       [ref.natural_travel, ref.natural_state, ref.null_mean_travel, ref.null_mean_state],
                                       rtol=1e-12, atol=1e-12)
            checked += 1
    return dict(rows_checked=checked, tolerance=1e-12,
                compared="natural and exact first-cousin mean costs against round_references metrics")


def compute(trees, *, draws=DRAWS, seed=SEED):
    index = rr.tree_index()
    rows, samples = [], {}
    for tree_seed, name in enumerate(("pooled", "partial")):
        r, s = analyze_tree(trees[name], index, draws=draws, seed=seed, tree_seed=tree_seed)
        rows += r
        samples.update(s)
    layers_table = pd.DataFrame(rows)
    contrast_table = contrasts(samples, layers_table)
    return dict(layers=layers_table, contrasts=contrast_table), samples


# ---- run folder ---------------------------------------------------------------------

def source_hashes():
    paths = (*CODE_SOURCES,
             str((pa.DEFAULT_RUN / rr.SUBDIR / "validation.json").relative_to(ROOT)),
             str((CROSS_RUN / rr.SUBDIR / "validation.json").relative_to(ROOT)))
    return {name: lf.file_hash(ROOT / name) for name in paths}


def settings(*, draws=DRAWS, seed=SEED):
    return dict(draws=draws, seed=seed, level=LEVEL, families=list(FAMILIES), geometry=GEOMETRY,
                layer="9 - round", mid=list(MID), late=list(LATE), early=list(EARLY), primary=PRIMARY,
                unit="shuffle groups within a round (unshuffled children as singletons)",
                numpy_version=np.__version__)


def load():
    trees = dict(pooled=load_pooled(), partial=load_partial())
    upstream = dict(pooled_round_references=trees["pooled"]["identity"],
                    pooled_context_key=trees["pooled"]["context_key"],
                    cross_species_round_references=trees["partial"]["identity"],
                    cross_species_analysis_id=trees["partial"]["analysis_id"])
    return trees, upstream


def write_run(out, *, draws=DRAWS, seed=SEED):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    run_settings = settings(draws=draws, seed=seed)
    trees, upstream = load()
    sources = source_hashes()
    results, samples = compute(trees, draws=draws, seed=seed)
    report = check_reproduction(results["layers"], trees)
    verdict = decision(results["contrasts"])
    analysis.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        results[name].to_csv(analysis / f"{name}.csv", index=False)
    np.savez_compressed(analysis / "draws.npz", **{f"{s}__{f}__round_{r}": v for (s, f, r), v in samples.items()})
    (analysis / "reproduction.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (analysis / "decision.json").write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=run_settings, upstream=upstream,
                                              sources=sources), sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=run_settings, upstream=upstream, sources=sources,
                  elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: lf.file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def validate_run(out, *, check_sources=True):
    """Check file and source hashes, then recompute every table and draw at 1e-12."""
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported layer-hourglass cache")
    for name, expected in record["files"].items():
        if lf.file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if lf.file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    trees, upstream = load()
    if upstream != record["upstream"]:
        raise ValueError("Upstream round references changed")
    run_settings = record["settings"]
    fresh, samples = compute(trees, draws=run_settings["draws"], seed=run_settings["seed"])
    check_reproduction(fresh["layers"], trees)
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        lf._compare(fresh[name], saved, name)
    with np.load(analysis / "draws.npz") as data:
        if set(data.files) != {f"{s}__{f}__round_{r}" for s, f, r in samples}:
            raise ValueError("Replay mismatch: draw keys")
        for (s, f, r), value in samples.items():
            np.testing.assert_allclose(value, data[f"{s}__{f}__round_{r}"], rtol=1e-12, atol=1e-12)
    if decision(fresh["contrasts"]) != json.loads((analysis / "decision.json").read_text()):
        raise ValueError("Replay mismatch: decision")
    return dict(analysis_id=record["analysis_id"], tables=len(TABLES), draw_sets=len(samples),
                verified_source_hashes=len(record["sources"]) if check_sources else 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = lf.run_path(args.run_id)
    if not args.verify_only:
        record = write_run(out)
        print(json.dumps({k: record[k] for k in ("analysis_id", "elapsed_seconds")}, indent=2))
    print(json.dumps(validate_run(out), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
