"""Apoptotic cells in the terminal analysis (Track D4).

Thirty of the 275 cohort children undergo apoptosis (the ``death`` kind of the
D0 cell table, from ``data/apoptotic_cells.txt``). Their place in the lineage
is exact, but their travel endpoint is wherever the track stopped before the
cutoff. That can be a corpse, a track lost at death or the cutoff itself, and
it differs between embryos. There are no death times in the data, so an
endpoint at death is not available.

Variants of the canonical pooled context:

- ``as_is``: the canonical matrices;
- ``drop_deaths``: apoptotic children removed;
- ``birth_endpoint``: apoptotic children stay, with their first tracked
  position in each embryo as the endpoint (they do not migrate). Pooled travel
  keeps the canonical per-embryo natural totals as denominators.

Each variant is computed for the whole cohort and for every canonical subtree,
with the measured and the fate-only cell-state cost, and D1's pair comparison
is redone so that a change in any within-pair difference shows up.

**Control.** Each draw picks one living cohort child per death, matched on
depth and birth frame: deaths are visited in random order, and each picks at
random among the ``POOL_SIZE`` nearest living children (same depth first, then
birth frame) not yet taken in that draw. Deaths born close together share
neighbours, so a fixed pool would run out. The same draw is used for both
operations: removing the picked cells, and moving their endpoints to birth.

**Decision rule** (fixed in the Track D brief before the run): if dropping the
deaths brings the natural lineage closer to the front than the matched control
on both d_SD (lower) and retention at the nearest front point (higher), each
one-sided at ``ALPHA``, the departure from the front is concentrated on
apoptotic cells and is reported as a finding. Otherwise deaths stay in the
canonical analysis and ``drop_deaths`` goes to the supplement.

Run: ``python -m terminal_pareto.apoptosis_sensitivity [--run-id ID]``;
``--verify-only`` replays every table.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from terminal_pareto import data_loader as dl
from terminal_pareto import subtree_symmetry as ss
from terminal_pareto import subtree_symmetry_robustness as sr
from terminal_pareto.front_coordinates import DegenerateEndpointSpan

VERSION = "terminal-apoptosis-sensitivity-v1"
DEFAULT_RUN = "apoptosis_20261019"
CODE_SOURCES = ("terminal_pareto/apoptosis_sensitivity.py",) + sr.CODE_SOURCES
VARIANTS = ("as_is", "drop_deaths", "birth_endpoint")
COSTS = ("measured", "fate_only")
FRONT_METRICS = ("d_SD", "natural_travel_z", "natural_state_z", "u_L", "d_LP", "d_NP", "max_edge_retention",
                 "retention_at_nearest", "travel_span_sd", "state_span_sd", "null_support", "tiny_span")
GATE_COLUMNS = (("d_SD", "d_lineage_sd"), ("natural_state_z", "natural_state_z"),
                ("natural_travel_z", "natural_travel_z"), ("u_L", "u_L"), ("d_LP", "d_LP"),
                ("max_edge_retention", "max_edge_retention"), ("retention_at_nearest", "retention_at_nearest"))
N_DRAWS = 1000
SEED = 20261019
POOL_SIZE = 8              # matched control: candidates per death, nearest in birth frame at the same depth
ALPHA = 0.05
UM = 0.1625                # tracks are in pixels, as in data_loader.load_elegans_tracking


@dataclass
class Context:
    compare: ss.CompareContext
    cells: pd.DataFrame        # D0 cell table, indexed by cell
    tiers: pd.DataFrame
    pair_comparison: pd.DataFrame
    fate_state: np.ndarray
    dead: np.ndarray           # bool over cohort children
    components: dict           # embryo -> canonical travel component
    parent_xyz: dict           # embryo -> parent endpoint per row
    child_xyz: dict            # embryo -> child endpoint per column
    birth_xyz: dict            # embryo -> child first position per column
    tracks: dict               # embryo -> per-child track summary
    gates: dict
    runs: dict


def _tracks(path, cut, names):
    frame = pd.read_csv(path, sep="\t")
    frame = frame[frame.name.isin(names)].sort_values(["name", "t"])
    first = frame.groupby("name").first()
    last = frame.groupby("name").last()
    # Same [z, y, x] order as the endpoints from data_loader.
    return pd.DataFrame(dict(born=first.t, end=last.t,
                             birth_0=first.z * UM, birth_1=first.y * UM, birth_2=first.x * UM)), cut


def load_context(pairing_run=ss.DEFAULT_RUN, compare_run=ss.DEFAULT_COMPARE_RUN):
    base = sr.load_context(pairing_run, compare_run)
    compare = base.compare
    children, parents = compare.children, compare.parents
    gates = dict(base.gates)
    dead = np.array([base.cells.at[c, "kind"] == "death" for c in children])
    components, parent_xyz, child_xyz, birth_xyz, tracks = {}, {}, {}, {}, {}
    for label, path, cut in dl.CE_REPLICATES:
        xyz, _ = dl.load_elegans_tracking(cut, path=path)
        p = np.array([xyz[q] for q in parents], dtype=float)
        c = np.array([xyz[q] for q in children], dtype=float)
        rebuilt = cdist(p, c)
        gates[f"{label}_component_max_abs_delta"] = float(np.max(np.abs(rebuilt - compare.travel[label])))
        if gates[f"{label}_component_max_abs_delta"] > 1e-12:
            raise AssertionError(f"{label}: rebuilt travel does not match the canonical component")
        summary, _ = _tracks(path, cut, set(children))
        summary = summary.loc[children]
        components[label], parent_xyz[label], child_xyz[label] = compare.travel[label], p, c
        birth_xyz[label] = summary[["birth_0", "birth_1", "birth_2"]].to_numpy()
        tracks[label] = summary.assign(cut=cut)
    ctx = Context(compare, base.cells, base.tiers, base.pair_comparison, base.variants["fate_only"][3], dead,
                  components, parent_xyz, child_xyz, birth_xyz, tracks, gates,
                  dict(pairing=pairing_run, compare=compare_run))
    # Whole-cohort gate: the Track A criteria front (top-20 cosine, measured and fate-only).
    criteria = sr._check_tables(sr.CRITERIA, ("fronts",))["fronts"].query(
        "dataset == 'CE protein' and geometry == 'raw3d' and candidate == 'top-20' and metric == 'cosine'"
    ).set_index("cost")
    travel = pooled_travel(ctx)
    for cost, theirs in (("measured", "measured"), ("fate_only", "fate-only")):
        row, _ = ss.front_metrics(travel, state(ctx, cost), parents, children, compare.tree_index,
                                  f"gate:whole:{cost}")
        for ours, column in GATE_COLUMNS:
            sr._close(f"whole_{cost}_{ours}", row[ours], criteria.at[theirs, column], gates)
    return ctx


# --- variant matrices ---------------------------------------------------------

def state(ctx, cost):
    return ctx.compare.state if cost == "measured" else ctx.fate_state


def embryo_travel(ctx, label, moved=()):
    """One embryo's travel component, with the ``moved`` children's endpoints at birth."""
    if not len(moved):
        return ctx.components[label]
    child = ctx.child_xyz[label].copy()
    child[list(moved)] = ctx.birth_xyz[label][list(moved)]
    return cdist(ctx.parent_xyz[label], child)


def pooled_travel(ctx, moved=()):
    """Mean of the embryo components over their canonical natural totals."""
    return np.mean([embryo_travel(ctx, e, moved) / np.trace(ctx.components[e]) for e in ss.TRAVELS[1:]], axis=0)


def variant_matrices(ctx, variant, cells=None):
    """(kept indices, travel by key) for a variant; ``cells`` defaults to the deaths."""
    n = len(ctx.compare.children)
    cells = np.flatnonzero(ctx.dead) if cells is None else np.asarray(cells)
    if variant == "as_is":
        return np.arange(n), dict(ctx.compare.travel)
    if variant == "drop_deaths":
        return np.setdiff1d(np.arange(n), cells), dict(ctx.compare.travel)
    if variant == "birth_endpoint":
        travel = {e: embryo_travel(ctx, e, cells) for e in ss.TRAVELS[1:]}
        travel["pooled"] = pooled_travel(ctx, cells)
        return np.arange(n), travel
    raise ValueError(variant)


def _front(ctx, idx, travel, cost, identity):
    children = [ctx.compare.children[i] for i in idx]
    try:
        row, _ = ss.front_metrics(travel[np.ix_(idx, idx)], state(ctx, cost)[np.ix_(idx, idx)],
                                  ctx.compare.parents[idx], children, ctx.compare.tree_index, identity)
    except DegenerateEndpointSpan:
        return dict(n=len(idx), status="degenerate_span")
    except ValueError as error:
        if "variances must be positive" not in str(error):
            raise
        return dict(n=len(idx), status="degenerate_null")
    return dict(n=len(idx), status="ok", **{k: row[k] for k in FRONT_METRICS})


# --- tables -------------------------------------------------------------------

def whole_cohort(ctx):
    rows = []
    for variant in VARIANTS:
        idx, travel = variant_matrices(ctx, variant)
        for cost in COSTS:
            for key in ss.TRAVELS:
                row = _front(ctx, idx, travel[key], cost, f"{VERSION}:whole:{variant}:{cost}:{key}")
                rows.append(dict(variant=variant, cost=cost, travel=key, n_deaths=int(ctx.dead[idx].sum()), **row))
    return pd.DataFrame(rows)


def subtree_metrics(ctx):
    position = {c: i for i, c in enumerate(ctx.compare.children)}
    rows = []
    for variant in VARIANTS:
        kept, travel = variant_matrices(ctx, variant)
        kept = set(kept.tolist())
        for node in ctx.compare.canonical_metrics.index:
            idx = np.array(sorted(position[c] for c in ctx.compare.lineage.descendants[node]
                                  if c in position and position[c] in kept))
            for cost in COSTS:
                for key in (ss.TRAVELS if cost == "measured" else ("pooled",)):
                    head = dict(variant=variant, cost=cost, travel=key, subtree=node, depth=ss.depth(node),
                                n_deaths=int(ctx.dead[idx].sum()) if len(idx) else 0)
                    if len(idx) < ss.MIN_EDGES:
                        rows.append(dict(head, n=len(idx), status="below_min_edges"))
                        continue
                    rows.append(dict(head, **_front(ctx, idx, travel[key], cost,
                                                    f"{VERSION}:subtree:{variant}:{cost}:{key}:{node}")))
    return pd.DataFrame(rows)


def pair_comparison(ctx, subtrees):
    """D1's comparison of kept pairs in each variant, with D2's 'differs' flag."""
    out = []
    for variant in VARIANTS:
        table = subtrees[(subtrees.variant == variant) & (subtrees.cost == "measured")]
        ok = table.groupby("subtree").status.agg(lambda s: (s == "ok").all())
        metrics = table[table.subtree.isin(ok.index[ok])]
        pairs = ctx.tiers[ctx.tiers.subtree_x.isin(metrics.subtree) & ctx.tiers.subtree_y.isin(metrics.subtree)]
        refs = ss.reference_pairs(ctx.compare, metrics, pairs)
        comparison = ss.pair_comparison(metrics, pairs, refs)
        comparison["differs"] = ((comparison.reference_percentile >= sr.DIFFERS_PERCENTILE)
                                 & (comparison.abs_difference > comparison.embryo_floor))
        out.append(comparison.assign(variant=variant))
    out = pd.concat(out, ignore_index=True)
    stored = ctx.pair_comparison.set_index(["subtree_x", "subtree_y", "metric"])
    mine = out[out.variant == "as_is"].set_index(["subtree_x", "subtree_y", "metric"])
    for column in ("abs_difference", "reference_percentile", "embryo_floor"):
        delta = float(np.nanmax(np.abs(mine[column] - stored.loc[mine.index, column])))
        ctx.gates[f"pair_comparison_as_is_{column}_max_abs_delta"] = delta
        if delta > ss.REPRODUCTION_TOL:
            raise AssertionError(f"as_is pair comparison does not reproduce D1 ({column})")
    base = mine.differs
    out["differs_as_is"] = [bool(base.get((r.subtree_x, r.subtree_y, r.metric), False)) for r in out.itertuples()]
    out["changed"] = out.differs != out.differs_as_is
    return out


def _birth_frame(ctx, cell):
    return float(np.median([ctx.tracks[e].at[cell, "born"] for e in ss.TRAVELS[1:]]))


def control_pools(ctx):
    """For each death, every living cohort child ranked by depth difference, then birth-frame gap."""
    children = ctx.compare.children
    living = [i for i in range(len(children)) if not ctx.dead[i]]
    info = {i: (ss.depth(children[i]), _birth_frame(ctx, children[i])) for i in range(len(children))}
    pools = {}
    for d in np.flatnonzero(ctx.dead):
        depth, born = info[d]
        ranked = sorted(living, key=lambda i: (abs(info[i][0] - depth), abs(info[i][1] - born), children[i]))
        pools[int(d)] = ranked
    return pools, info


def draw_controls(ctx, pools, rng):
    """One matched draw: a distinct living child for each death, from its nearest free candidates."""
    taken = []
    for d in rng.permutation(sorted(pools)):
        free = [i for i in pools[int(d)] if i not in taken][:POOL_SIZE]
        if not free:
            raise AssertionError("No living candidate left for a death")
        taken.append(int(rng.choice(free)))
    return sorted(taken)


def controls(ctx):
    pools, _ = control_pools(ctx)
    rng = np.random.default_rng(SEED)
    rows = []
    for k in range(N_DRAWS):
        picked = draw_controls(ctx, pools, rng)
        names = ";".join(ctx.compare.children[i] for i in picked)
        for variant in VARIANTS[1:]:
            idx, travel = variant_matrices(ctx, variant, picked)
            row = _front(ctx, idx, travel["pooled"], "measured", f"{VERSION}:control:{variant}:{k}")
            rows.append(dict(draw=k, variant=variant, cells=names, **row))
    return pd.DataFrame(rows)


def p_value(control, observed, lower_is_better):
    hits = (control <= observed) if lower_is_better else (control >= observed)
    return float((1 + hits.sum()) / (len(control) + 1))


def decision(whole, control):
    rows = []
    observed = whole[(whole.cost == "measured") & (whole.travel == "pooled")].set_index("variant")
    for variant in VARIANTS[1:]:
        c = control[(control.variant == variant) & (control.status == "ok")]
        row = dict(variant=variant, n_draws=len(c))
        for metric, lower in (("d_SD", True), ("retention_at_nearest", False),
                              ("natural_travel_z", True), ("natural_state_z", True)):
            obs = observed.at[variant, metric]
            row.update({f"{metric}_observed": obs, f"{metric}_as_is": observed.at["as_is", metric],
                        f"{metric}_control_median": float(c[metric].median()),
                        f"{metric}_control_q05": float(c[metric].quantile(0.05)),
                        f"{metric}_control_q95": float(c[metric].quantile(0.95)),
                        f"{metric}_p": p_value(c[metric], obs, lower)})
        row["finding"] = bool(row["d_SD_p"] < ALPHA and row["retention_at_nearest_p"] < ALPHA)
        rows.append(row)
    return pd.DataFrame(rows)


def deaths(ctx):
    """One row per apoptotic cohort child and embryo: track, costs, null role and control pool."""
    from terminal_pareto import pareto_engine as pe
    children, parents = ctx.compare.children, ctx.compare.parents
    groups = pe.build_ancestor_groups(children, ctx.compare.tree_index, 2)
    in_null = {i for g in groups if len({parents[j] for j in g}) > 1 for i in g}
    pools, info = control_pools(ctx)
    subtrees = ctx.compare.canonical_metrics.index
    rows = []
    for d in np.flatnonzero(ctx.dead):
        c = children[d]
        pool = pools[int(d)][:POOL_SIZE]
        for e in ss.TRAVELS[1:]:
            t = ctx.tracks[e].loc[c]
            rows.append(dict(
                cell=c, embryo=e, depth=info[d][0], born=int(t.born), track_end=int(t.end), cut=int(t.cut),
                ends_before_cut=bool(t.end < t.cut), frames_to_endpoint=int(min(t.end, t.cut) - t.born),
                birth_to_endpoint=float(np.linalg.norm(ctx.child_xyz[e][d] - ctx.birth_xyz[e][d])),
                natural_travel=float(ctx.components[e][d, d]), natural_state=float(ctx.compare.state[d, d]),
                in_cousin_null=d in in_null,
                n_subtrees=int(sum(c in ctx.compare.lineage.descendants[s] for s in subtrees)),
                pool_depth_matched=int(sum(info[i][0] == info[d][0] for i in pool)),
                pool_median_birth_gap=float(np.median([abs(info[i][1] - info[d][1]) for i in pool]))))
    return pd.DataFrame(rows)


TABLES = ("deaths", "whole_cohort", "subtree_metrics", "pair_comparison", "controls", "decision")


def compute(ctx):
    whole = whole_cohort(ctx)
    subtrees = subtree_metrics(ctx)
    control = controls(ctx)
    return dict(deaths=deaths(ctx), whole_cohort=whole, subtree_metrics=subtrees,
                pair_comparison=pair_comparison(ctx, subtrees), controls=control,
                decision=decision(whole, control))


def provenance(tables, ctx):
    sources = CODE_SOURCES + (ss.CELL_KEY, ss.LINEAGE, ss.PROTEIN, ss.SELECTED) + tuple(
        str(Path(path).resolve().relative_to(ss.ROOT)) for _, path, _ in dl.CE_REPLICATES)
    caches = {str((folder / name).relative_to(ss.ROOT)): ss.file_hash(folder / name) for folder, name in (
        (ss.CANONICAL, "analysis_context.npz"), (ss.CANONICAL, "analysis_manifest.json"),
        (sr.CRITERIA, "provenance.json"), (sr.CROSS, "provenance.json"))}
    return dict(
        version=VERSION, pairing_run=ctx.runs["pairing"], compare_run=ctx.runs["compare"],
        cache_hashes=caches, source_hashes={s: ss.file_hash(ss.ROOT / s) for s in sources},
        settings=dict(variants=list(VARIANTS), costs=list(COSTS), n_draws=N_DRAWS, seed=SEED, pool_size=POOL_SIZE,
                      alpha=ALPHA, differs_percentile=sr.DIFFERS_PERCENTILE, min_edges=ss.MIN_EDGES,
                      intervals=ss.INTERVALS),
        gates=ctx.gates,
        table_hashes={t: hashlib.sha256(tables[t].to_csv(index=False).encode()).hexdigest() for t in TABLES})


def write_run(out, **runs):
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists; runs are write-once")
    ctx = load_context(**runs)
    tables = compute(ctx)
    analysis = out / "analysis"
    analysis.mkdir(parents=True)
    for name in TABLES:
        tables[name].to_csv(analysis / f"{name}.csv", index=False)
    (analysis / "provenance.json").write_text(json.dumps(provenance(tables, ctx), indent=2, sort_keys=True) + "\n")
    return tables


def verify_run(out):
    analysis = Path(out) / "analysis"
    stored = json.loads((analysis / "provenance.json").read_text())
    if stored["version"] != VERSION:
        raise AssertionError(f"Run version {stored['version']} != {VERSION}")
    ctx = load_context(stored["pairing_run"], stored["compare_run"])
    tables = compute(ctx)
    fresh = provenance(tables, ctx)
    for group in ("cache_hashes", "source_hashes"):
        changed = [k for k in stored[group] if stored[group][k] != fresh[group].get(k) and k not in CODE_SOURCES]
        if changed:
            raise AssertionError(f"{group} changed: {changed}")
    code_changed = [k for k in CODE_SOURCES if stored["source_hashes"].get(k) != fresh["source_hashes"][k]]
    if code_changed:
        print(f"note: code changed since the run ({', '.join(code_changed)}); tables replayed")
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        now = pd.read_csv(io.StringIO(tables[name].to_csv(index=False)), keep_default_na=False, na_values=[""])
        pd.testing.assert_frame_equal(saved, now, check_exact=False, rtol=0, atol=1e-12)
    return tables


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--pairing-run", default=ss.DEFAULT_RUN)
    parser.add_argument("--compare-run", default=ss.DEFAULT_COMPARE_RUN)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = ss.run_path(args.run_id)
    if args.verify_only:
        verify_run(out)
        print(f"verified {out}")
    else:
        write_run(out, pairing_run=args.pairing_run, compare_run=args.compare_run)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
