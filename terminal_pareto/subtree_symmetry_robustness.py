"""Robustness and classification of left-right subtree differences (Track D2).

Reads the D0 pairing and D1 comparison runs of ``subtree_symmetry`` and asks,
for every kept pair, whether a within-pair difference survives:

- **XY travel** instead of 3D (pooled over the three embryos);
- the **fate-only cell-state cost** (children replaced by the mean profile of
  their fine type, as in ``representation_criteria``);
- the **187-edge cross-species cohort**: C. elegans protein and RNA, and
  C. briggsae RNA, in 3D and XY, where both partners keep at least 12 edges.
  RNA is resolved by cell type, so it tests fate-level agreement only;
- a **matched cohort**: only cells whose mirror counterpart is also in the
  cohort, which removes coverage differences between partners.

Every variant except the matched cohort is computed for all subtrees with at
least 12 edges, so each has its own matched random-pair reference. Each
variant is gated on an existing cache before use.

Each pair is then classified from recorded evidence:

- **numerical**: the difference is in an endpoint metric (u, d_LP, d_CP) and a
  partner's endpoint span is under one null SD;
- **data**: a partner's first-cousin null rests on under half its cells, the
  cohorts are unbalanced (ratio < 0.8), late-born cells are missing from one
  side only, or the difference disappears in the matched cohort;
- **biological**: the difference survives the data checks and the partners
  differ in fate (cells without a compatible counterpart, or known
  left-right asymmetric fates).

Run: ``python -m terminal_pareto.subtree_symmetry_robustness [--run-id ID]``;
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

from terminal_pareto import data_loader as dl
from terminal_pareto import subtree_symmetry as ss
from terminal_pareto.front_coordinates import DegenerateEndpointSpan

VERSION = "terminal-subtree-symmetry-robustness-v1"
DEFAULT_RUN = "robustness_20261015"
CRITERIA = ss.ROOT / "terminal_pareto/output/runs/representation_criteria_v1/criteria_metric_20261014/analysis"
CROSS = ss.ROOT / "terminal_pareto/output/runs/cross_species_terminal_v1/pooled_comparison_20260930/analysis"
CODE_SOURCES = ("terminal_pareto/subtree_symmetry_robustness.py", "terminal_pareto/subtree_symmetry.py",
                "terminal_pareto/feature_sensitivity.py", "terminal_pareto/representation_criteria.py")
CROSS_CONFIGS = tuple((g, c) for g in ("raw3d", "xy") for c in ("ce_protein", "ce_rna", "cb_rna"))
METRICS = ("d_SD", "natural_state_z", "natural_travel_z", "u_L", "d_LP", "d_NP")
ENDPOINT_METRICS = ("u_L", "d_LP", "d_NP")
DIFFERS_PERCENTILE = 0.9   # a difference is notable at or above this reference percentile ...
LOW_NULL_SUPPORT = 0.5     # ... a partner's null rests on under this share of its cells ...
COHORT_IMBALANCE = 0.8     # ... or the smaller cohort is under this share of the larger
LATE_BORN_GAP = 3          # untracked cells on one side minus the other, in any embryo


@dataclass
class Context:
    compare: ss.CompareContext
    tiers: pd.DataFrame
    pair_comparison: pd.DataFrame
    pairs: pd.DataFrame        # D0 candidate pairs (coverage columns)
    mismatches: pd.DataFrame
    cells: pd.DataFrame
    variants: dict             # variant -> (children, parents, travel, state)
    gates: dict
    runs: dict


def _check_tables(folder, names):
    record = json.loads((folder / "provenance.json").read_text())
    # Track D runs hash table text; the criteria run hashes files.
    stored = record.get("table_hashes") or {k.removesuffix(".csv"): v for k, v in record["files"].items()}
    for name in names:
        if ss.file_hash(folder / f"{name}.csv") != stored[name]:
            raise ValueError(f"{folder / name}.csv does not match its provenance")
    return {name: pd.read_csv(folder / f"{name}.csv") for name in names}


def _close(label, actual, expected, gates, tol=ss.REPRODUCTION_TOL):
    delta = float(abs(actual - expected))
    gates[label] = delta
    if not delta <= tol:
        raise AssertionError(f"{label}: {actual!r} vs {expected!r}")


def load_context(pairing_run=ss.DEFAULT_RUN, compare_run=ss.DEFAULT_COMPARE_RUN):
    from terminal_pareto import feature_sensitivity as fs
    from terminal_pareto import representation_criteria as rc

    compare = ss.load_compare_context(pairing_run)
    d0 = _check_tables(ss.run_path(pairing_run) / "analysis", ("cells", "candidate_pairs", "pair_mismatches"))
    d1 = _check_tables(ss.run_path(compare_run) / "analysis", ("pair_tiers", "pair_comparison"))
    gates, variants = {}, {}
    children, parents = compare.children, compare.parents
    terms = list(zip(children, parents.tolist()))

    # XY: pooled exactly as the cross-species cache builds it (3D components are checked inside).
    xy = fs.primary_xy_travel(type("Primary", (), {"terms": terms})())
    variants["xy"] = (children, parents, xy, compare.state)

    # Fate-only: children as the mean profile of their fate group over the 275 cohort children.
    criteria = _check_tables(CRITERIA, ("cell_types", "fronts"))
    groups = criteria["cell_types"].set_index("cell").fate_group.to_dict()
    frame = dl.load_protein_expression(str(ss.ROOT / ss.PROTEIN))
    features = dl.load_prot_sel(str(ss.ROOT / ss.SELECTED))
    measured = rc.state_matrix(frame, features, terms)
    gates["measured_state_max_abs_delta"] = float(np.max(np.abs(measured - compare.state)))
    if gates["measured_state_max_abs_delta"] > 1e-12:
        raise AssertionError("Measured state does not reproduce the canonical matrix")
    fate = rc.state_matrix(frame, features, terms, fate_groups=groups)
    variants["fate_only"] = (children, parents, compare.travel["pooled"], fate)
    expected = criteria["fronts"].query(
        "dataset == 'CE protein' and geometry == 'raw3d' and candidate == 'top-20' "
        "and cost == 'fate-only' and metric == 'cosine'").iloc[0]
    row, _ = ss.front_metrics(compare.travel["pooled"], fate, parents, children, compare.tree_index, "gate:fate")
    for ours, theirs in (("d_SD", "d_lineage_sd"), ("natural_state_z", "natural_state_z"), ("u_L", "u_L"),
                         ("d_LP", "d_LP")):
        _close(f"fate_only_P0_{ours}", row[ours], expected[theirs], gates)

    # Cross-species cohort: published 187-edge matrices.
    record = json.loads((CROSS / "provenance.json").read_text())
    path = CROSS / "assignments_and_costs.npz"
    if ss.file_hash(path) != record["files"][path.name]:
        raise ValueError("Cross-species cache does not match its provenance")
    cached = pd.read_csv(CROSS / "metrics.csv").query("cohort == 'base'").set_index(["geometry", "config"])
    with np.load(path, allow_pickle=False) as arrays:
        cross_children, cross_parents = arrays["terminals"].tolist(), np.asarray(arrays["parents"].tolist())
        for geometry, config in CROSS_CONFIGS:
            key = f"{geometry}__{config}"
            variants[f"cross187_{key}"] = (cross_children, cross_parents, np.asarray(arrays[f"{key}__travel"]),
                                           np.asarray(arrays[f"{key}__state"]))
    for geometry, config in CROSS_CONFIGS:
        c, p, t, s = variants[f"cross187_{geometry}__{config}"]
        row, _ = ss.front_metrics(t, s, p, c, compare.tree_index, f"gate:{geometry}:{config}")
        for metric in ("u_L", "d_LP", "d_NP", "max_edge_retention"):
            _close(f"cross187_{geometry}_{config}_P0_{metric}", row[metric],
                   cached.loc[(geometry, config), metric], gates)
    # XY gate: the 275-edge XY matrix restricted to the 187 edges is the cached one.
    lookup = {edge: i for i, edge in enumerate(terms)}
    index = np.asarray([lookup[edge] for edge in zip(cross_children, cross_parents.tolist())])
    gates["xy_vs_cross187_max_abs_delta"] = float(np.max(np.abs(
        xy[np.ix_(index, index)] - variants["cross187_xy__ce_protein"][2])))
    if gates["xy_vs_cross187_max_abs_delta"] > 1e-12:
        raise AssertionError("XY travel does not match the cross-species cache")
    return Context(compare, d1["pair_tiers"], d1["pair_comparison"], d0["candidate_pairs"],
                   d0["pair_mismatches"], d0["cells"].set_index("cell"), variants, gates,
                   dict(pairing=pairing_run, compare=compare_run))


# --- variant fronts ---------------------------------------------------------

def variant_metrics(ctx):
    """Every subtree with at least 12 edges in each variant's cohort."""
    rows = []
    nodes = list(ctx.compare.lineage.descendants)
    for variant, (children, parents, travel, state) in ctx.variants.items():
        position = {c: i for i, c in enumerate(children)}
        for node in sorted(nodes):
            idx = np.array(sorted(position[c] for c in ctx.compare.lineage.descendants[node] if c in position))
            if len(idx) < ss.MIN_EDGES:
                continue
            cells = [children[i] for i in idx]
            try:
                row, _ = ss.front_metrics(travel[np.ix_(idx, idx)], state[np.ix_(idx, idx)], parents[idx], cells,
                                          ctx.compare.tree_index, f"{VERSION}:{variant}:{node}")
            except (ValueError, DegenerateEndpointSpan) as error:
                # Type-level states (fate-only, RNA) can leave a subtree with no cousin
                # variance, or with coinciding travel and state optima.
                if isinstance(error, DegenerateEndpointSpan):
                    row = dict(n=len(cells), status="degenerate_span")
                elif "variances must be positive" in str(error):
                    row = dict(n=len(cells), status="degenerate_null")
                else:
                    raise
            rows.append({"variant": variant, "subtree": node, "depth": ss.depth(node), "status": "ok", **row})
    out = pd.DataFrame(rows)
    return out


def matched_metrics(ctx):
    """Each kept pair on the 3D pooled front, keeping only cells whose counterpart is in the cohort."""
    cls = {c: dict(fate_class=r.fate_class, side=r.side if isinstance(r.side, str) else None)
           for c, r in ctx.cells.iterrows()}
    position = {c: i for i, c in enumerate(ctx.compare.children)}
    rows = []
    for r in ctx.tiers[ctx.tiers.tier != "excluded"].itertuples():
        pos, _ = ss.mirror_pairs(r.subtree_x, r.subtree_y, ctx.compare.lineage, cls, r.position_phi)
        kept = [(a, b) for a, b in pos if a in position and b in position]
        for side, cells in ((r.subtree_x, [a for a, _ in kept]), (r.subtree_y, [b for _, b in kept])):
            row = dict(subtree_x=r.subtree_x, subtree_y=r.subtree_y, subtree=side, n=len(cells))
            if len(cells) >= 4:
                idx = np.array(sorted(position[c] for c in cells))
                try:
                    front, _ = ss.front_metrics(
                        ctx.compare.travel["pooled"][np.ix_(idx, idx)], ctx.compare.state[np.ix_(idx, idx)],
                        ctx.compare.parents[idx], [ctx.compare.children[i] for i in idx],
                        ctx.compare.tree_index, f"{VERSION}:matched:{side}")
                    row.update({k: v for k, v in front.items() if k != "n"}, status="ok")
                except (ValueError, DegenerateEndpointSpan) as error:
                    row["status"] = ("degenerate_span" if isinstance(error, DegenerateEndpointSpan)
                                     else "degenerate_null")
            else:
                row["status"] = "too_few_cells"
            rows.append(row)
    return pd.DataFrame(rows)


def _reference_diffs(table, lineage, x, y, metric, mirror):
    """|difference| over unrelated, non-mirror pairs matched for size and depth (as in D1)."""
    table = table.set_index("subtree")
    if x not in table.index or y not in table.index:
        return None
    lo, hi = sorted((table.at[x, "n"], table.at[y, "n"]))
    d = (ss.depth(x) + ss.depth(y)) / 2
    nodes = list(table.index)
    out = []
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            if lineage.related(a, b) or tuple(sorted((a, b))) in mirror:
                continue
            na, nb = sorted((table.at[a, "n"], table.at[b, "n"]))
            if (abs(ss.depth(a) - d) <= ss.REF_DEPTH and abs(ss.depth(b) - d) <= ss.REF_DEPTH
                    and na >= ss.REF_SIZE[0] * lo and nb <= ss.REF_SIZE[1] * hi):
                out.append(abs(table.at[a, metric] - table.at[b, metric]))
    return np.asarray(out)


def pair_robustness(ctx, metrics, matched):
    """Within-pair differences in every variant, with matched random references."""
    mirror = {tuple(sorted(p)) for p in zip(ctx.tiers.subtree_x, ctx.tiers.subtree_y)}
    primary = ctx.pair_comparison.set_index(["subtree_x", "subtree_y", "metric"])
    rows = []
    for r in ctx.tiers[ctx.tiers.tier != "excluded"].itertuples():
        x, y = r.subtree_x, r.subtree_y
        for metric in METRICS:
            base = primary.loc[(x, y, metric)]
            sign = np.sign(base.value_x - base.value_y)
            rows.append(dict(subtree_x=x, subtree_y=y, tier=r.tier, variant="primary_3d", metric=metric,
                             value_x=base.value_x, value_y=base.value_y, abs_difference=base.abs_difference,
                             n_reference=base.n_reference, reference_percentile=base.reference_percentile,
                             same_sign=True))
            for variant in ctx.variants:
                table = metrics[(metrics.variant == variant) & (metrics.status == "ok")]
                vx = table.loc[table.subtree == x, metric]
                vy = table.loc[table.subtree == y, metric]
                if vx.empty or vy.empty:
                    continue
                ref = _reference_diffs(table, ctx.compare.lineage, x, y, metric, mirror)
                diff = abs(vx.iloc[0] - vy.iloc[0])
                rows.append(dict(subtree_x=x, subtree_y=y, tier=r.tier, variant=variant, metric=metric,
                                 value_x=vx.iloc[0], value_y=vy.iloc[0], abs_difference=diff,
                                 n_reference=len(ref),
                                 reference_percentile=(float((ref < diff).mean() + 0.5 * (ref == diff).mean())
                                                       if len(ref) else np.nan),
                                 same_sign=bool(np.sign(vx.iloc[0] - vy.iloc[0]) == sign or sign == 0)))
            m = matched[(matched.subtree_x == x) & (matched.subtree_y == y) & (matched.status == "ok")]
            if len(m) == 2 and metric in m:
                vx, vy = m.set_index("subtree").loc[[x, y], metric]
                rows.append(dict(subtree_x=x, subtree_y=y, tier=r.tier, variant="matched_3d", metric=metric,
                                 value_x=vx, value_y=vy, abs_difference=abs(vx - vy), n_reference=0,
                                 reference_percentile=np.nan,
                                 same_sign=bool(np.sign(vx - vy) == sign or sign == 0)))
    return pd.DataFrame(rows)


# --- classification ---------------------------------------------------------

def classify(ctx, robustness, matched):
    """One row per kept pair: evidence flags, the differing readouts and a class."""
    comparison = ctx.pair_comparison
    pairs = ctx.pairs.set_index(["subtree_x", "subtree_y"])
    cohort = set(ctx.cells.index[ctx.cells.in_cohort])
    rows = []
    for r in ctx.tiers[ctx.tiers.tier != "excluded"].itertuples():
        x, y = r.subtree_x, r.subtree_y
        c = comparison[(comparison.subtree_x == x) & (comparison.subtree_y == y)].set_index("metric")
        p = pairs.loc[(x, y)]
        rob = robustness[(robustness.subtree_x == x) & (robustness.subtree_y == y)]
        differs = [m for m in METRICS if c.at[m, "reference_percentile"] >= DIFFERS_PERCENTILE
                   and c.at[m, "abs_difference"] > c.at[m, "embryo_floor"]]
        support = min(c.at["d_SD", "null_support_x"] / p.cohort_x, c.at["d_SD", "null_support_y"] / p.cohort_y)
        late_gap = max(abs(p[f"not_tracked_{e}_x"] - p[f"not_tracked_{e}_y"]) for e in ss.TRAVELS[1:])
        mism = ctx.mismatches[(ctx.mismatches.subtree_x == x) & (ctx.mismatches.subtree_y == y)]
        fate_cells = mism[mism.fate_unmatched & ~mism.kind.isin(["death", "unannotated"])]
        known = [c_ for c_ in ctx.compare.lineage.descendants[x] + ctx.compare.lineage.descendants[y]
                 if isinstance(ctx.cells.at[c_, "known_asymmetry"], str) and ctx.cells.at[c_, "known_asymmetry"]]
        flags = dict(
            tiny_span=bool(c.at["d_SD", "tiny_span_x"] or c.at["d_SD", "tiny_span_y"]),
            low_null_support=bool(support < LOW_NULL_SUPPORT),
            cohort_imbalance=bool(p.cohort_ratio < COHORT_IMBALANCE),
            late_born_one_side=bool(late_gap >= LATE_BORN_GAP))
        survives = {}
        for m in differs:
            sub = rob[rob.metric == m].set_index("variant")
            survives[m] = {v: bool(sub.at[v, "same_sign"] and (v == "matched_3d" or
                                                               sub.at[v, "reference_percentile"] >= DIFFERS_PERCENTILE))
                           for v in sub.index if v != "primary_3d"}
        matched_gone = any(not s.get("matched_3d", True) for s in survives.values())
        reasons = []
        if differs:
            if all(m in ENDPOINT_METRICS for m in differs) and flags["tiny_span"]:
                label = "numerical"
                reasons.append("only endpoint metrics differ and a span is under one null SD")
            elif flags["low_null_support"] or flags["cohort_imbalance"] or flags["late_born_one_side"] or matched_gone:
                label = "data"
                reasons += [k for k in ("low_null_support", "cohort_imbalance", "late_born_one_side") if flags[k]]
                if matched_gone:
                    reasons.append("difference gone in the matched cohort")
            elif len(fate_cells) or known:
                label = "biological"
                reasons.append(f"{len(fate_cells)} cells without a fate counterpart")
            else:
                label = "unexplained"
        else:
            label = "symmetric"
        rows.append(dict(
            subtree_x=x, subtree_y=y, tier=r.tier, label=label, differs=";".join(differs),
            reasons="; ".join(reasons),
            survives=";".join(f"{m}:{','.join(v for v, ok in s.items() if ok) or 'none'}"
                              for m, s in survives.items()),
            min_null_support_share=support, cohort_ratio=p.cohort_ratio, late_born_gap=int(late_gap),
            n_fate_unmatched=len(fate_cells),
            n_fate_unmatched_in_cohort=int(fate_cells.cell.isin(cohort).sum()),
            fate_unmatched_names=";".join(sorted(set(fate_cells.name))),
            known_asymmetry=";".join(sorted({ctx.cells.at[k, "name"] for k in known})),
            **flags))
    return pd.DataFrame(rows)


# --- run --------------------------------------------------------------------

TABLES = ("variant_metrics", "matched_metrics", "pair_robustness", "classification")


def compute(ctx):
    metrics = variant_metrics(ctx)
    matched = matched_metrics(ctx)
    robustness = pair_robustness(ctx, metrics, matched)
    return dict(variant_metrics=metrics, matched_metrics=matched, pair_robustness=robustness,
                classification=classify(ctx, robustness, matched))


def provenance(tables, ctx):
    sources = CODE_SOURCES + (ss.CELL_KEY, ss.LINEAGE, ss.PROTEIN, ss.SELECTED)
    caches = {str((folder / name).relative_to(ss.ROOT)): ss.file_hash(folder / name) for folder, name in (
        (ss.CANONICAL, "analysis_context.npz"), (ss.CANONICAL, "analysis_manifest.json"),
        (CRITERIA, "provenance.json"), (CROSS, "provenance.json"), (CROSS, "assignments_and_costs.npz"),
        (CROSS, "metrics.csv"))}
    return dict(
        version=VERSION, pairing_run=ctx.runs["pairing"], compare_run=ctx.runs["compare"],
        cache_hashes=caches, source_hashes={s: ss.file_hash(ss.ROOT / s) for s in sources},
        settings=dict(variants=list(ctx.variants), metrics=list(METRICS), differs_percentile=DIFFERS_PERCENTILE,
                      low_null_support=LOW_NULL_SUPPORT, cohort_imbalance=COHORT_IMBALANCE,
                      late_born_gap=LATE_BORN_GAP, min_edges=ss.MIN_EDGES, intervals=ss.INTERVALS),
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
