"""Matched control for apoptotic cells inside subtrees (Track D4b).

D4 found that, inside the posterior ABpl subtrees, d_SD goes to 0 when the
deaths are dropped or held at their birth position. This asks whether that is
specific to deaths: each draw replaces the subtree's deaths with as many
living cells of the same subtree, matched on depth and birth frame with the
D4 sampler, and applies the same operation to them.

Subtrees tested are every canonical subtree with a death and a pooled as-is
d_SD above ``ZERO``, except P0 (the whole cohort, done in D4). ABplp is primary,
the rest of the posterior ABpl family supports it, and the others are
descriptive.

**Decision rule** (fixed in the Track D brief before the run): if, in ABplp,
dropping the deaths beats the matched control on d_SD (lower) and retention at
the nearest front point (higher), each one-sided at ``ALPHA``, the ABplp/ABprp
difference is caused by the deaths' endpoints. Otherwise it is a small-subtree
effect that any matched removal reproduces.

Run: ``python -m terminal_pareto.apoptosis_subtree_control [--run-id ID]``;
``--verify-only`` replays every table.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from terminal_pareto import apoptosis_sensitivity as ap
from terminal_pareto import subtree_symmetry as ss

VERSION = "terminal-apoptosis-subtree-control-v1"
DEFAULT_RUN = "apoptosis_subtrees_20261021"
CODE_SOURCES = ("terminal_pareto/apoptosis_subtree_control.py",) + ap.CODE_SOURCES
VARIANTS = ap.VARIANTS
PRIMARY = "ABplp"
SUPPORTING = ("ABp", "ABpl", "ABplpp", "ABplppp")
METRICS = (("d_SD", True), ("retention_at_nearest", False), ("natural_travel_z", True), ("natural_state_z", True))
N_DRAWS = 1000
SEED = 20261021
ZERO = 1e-9                # d_SD at or below this is on the front (ABalpp's as-is 4.5e-15 is float noise)


def subtree_index(ctx, node):
    position = {c: i for i, c in enumerate(ctx.compare.children)}
    return np.array(sorted(position[c] for c in ctx.compare.lineage.descendants[node] if c in position))


def front(ctx, node, variant, cells, identity):
    """Pooled measured-cost front of ``node`` after applying ``variant`` to ``cells``."""
    kept, travel = ap.variant_matrices(ctx, variant, cells)
    idx = np.intersect1d(subtree_index(ctx, node), kept)
    if len(idx) < ss.MIN_EDGES:
        return dict(n=len(idx), status="below_min_edges")
    return ap._front(ctx, idx, travel["pooled"], "measured", identity)


def role(node):
    return "primary" if node == PRIMARY else "supporting" if node in SUPPORTING else "descriptive"


def observed(ctx, stored):
    """Every variant on every subtree with a death; gated against D4's subtree table."""
    rows = []
    for node in ctx.compare.canonical_metrics.index:
        idx = subtree_index(ctx, node)
        deaths = idx[ctx.dead[idx]]
        if node == "P0" or not len(deaths):
            continue
        for variant in VARIANTS:
            rows.append(dict(subtree=node, variant=variant, n_deaths=len(deaths),
                             **front(ctx, node, variant, deaths, f"{VERSION}:observed:{variant}:{node}")))
    out = pd.DataFrame(rows)
    theirs = stored.query("cost == 'measured' and travel == 'pooled'").set_index(["subtree", "variant"])
    mine = out.set_index(["subtree", "variant"])
    ok = mine.status == "ok"
    if not (theirs.loc[mine.index, "status"] == mine.status).all():
        raise AssertionError("Subtree statuses do not reproduce D4")
    delta = float(np.max(np.abs(mine.loc[ok, "d_SD"] - theirs.loc[mine.index[ok], "d_SD"])))
    ctx.gates["observed_d_SD_max_abs_delta"] = delta
    if delta > 1e-12:
        raise AssertionError("Subtree d_SD does not reproduce D4")
    as_is = out[out.variant == "as_is"].set_index("subtree").d_SD
    tested = [s for s in as_is.index if as_is[s] > ZERO]
    if PRIMARY not in tested or not set(SUPPORTING) <= set(tested):
        raise AssertionError("The ABpl family is not among the tested subtrees")
    out["tested"] = out.subtree.isin(tested)
    out["role"] = out.subtree.map(role)
    return out


def control_pools(ctx, node, info):
    """For each death in ``node``, the subtree's living cells ranked by depth gap, then birth-frame gap."""
    idx = subtree_index(ctx, node)
    living = [int(i) for i in idx if not ctx.dead[i]]
    names = ctx.compare.children
    return {int(d): sorted(living, key=lambda i: (abs(info[i][0] - info[d][0]), abs(info[i][1] - info[d][1]),
                                                  names[i]))
            for d in idx[ctx.dead[idx]]}


def _info(ctx):
    return {i: (ss.depth(c), ap._birth_frame(ctx, c)) for i, c in enumerate(ctx.compare.children)}


def pools_table(ctx, nodes):
    info = _info(ctx)
    names = ctx.compare.children
    rows = []
    for node in nodes:
        for d, ranked in control_pools(ctx, node, info).items():
            for rank, i in enumerate(ranked[:ap.POOL_SIZE]):
                rows.append(dict(subtree=node, death=names[d], rank=rank, candidate=names[i],
                                 depth_gap=abs(info[i][0] - info[d][0]), birth_gap=abs(info[i][1] - info[d][1])))
    return pd.DataFrame(rows)


def controls(ctx, nodes):
    info = _info(ctx)
    rows = []
    for node in nodes:
        pools = control_pools(ctx, node, info)
        rng = np.random.default_rng([SEED, zlib.crc32(node.encode())])
        for k in range(N_DRAWS):
            picked = ap.draw_controls(ctx, pools, rng)
            names = ";".join(ctx.compare.children[i] for i in picked)
            for variant in VARIANTS[1:]:
                row = front(ctx, node, variant, picked, f"{VERSION}:control:{variant}:{node}:{k}")
                rows.append(dict(subtree=node, draw=k, variant=variant, cells=names, **row))
    return pd.DataFrame(rows)


def p_value(control, observed, lower_is_better):
    """D4's one-sided p, with ties counted within ``ZERO``."""
    shifted = observed + ZERO if lower_is_better else observed - ZERO
    return ap.p_value(control, shifted, lower_is_better)


def decision(observed_table, control):
    obs = observed_table.set_index(["subtree", "variant"])
    rows = []
    for (node, variant), c in control.groupby(["subtree", "variant"], sort=False):
        o = obs.loc[(node, variant)]
        row = dict(subtree=node, role=role(node), variant=variant, status=o.status,
                   n_draws=len(c), n_ok=int((c.status == "ok").sum()), n_distinct=c.cells.nunique())
        c = c[c.status == "ok"]
        if o.status == "ok" and len(c):
            for metric, lower in METRICS:
                row.update({f"{metric}_as_is": obs.at[(node, "as_is"), metric], f"{metric}_observed": o[metric],
                            f"{metric}_control_median": float(c[metric].median()),
                            f"{metric}_control_q05": float(c[metric].quantile(0.05)),
                            f"{metric}_control_q95": float(c[metric].quantile(0.95)),
                            f"{metric}_p": p_value(c[metric], o[metric], lower)})
            row["finding"] = bool(row["d_SD_p"] < ap.ALPHA and row["retention_at_nearest_p"] < ap.ALPHA)
        rows.append(row)
    return pd.DataFrame(rows)


TABLES = ("observed", "pools", "controls", "decision")


def compute(ctx, stored):
    seen = observed(ctx, stored)
    nodes = list(dict.fromkeys(seen[seen.tested].subtree))
    control = controls(ctx, nodes)
    return dict(observed=seen, pools=pools_table(ctx, nodes), controls=control,
                decision=decision(seen, control))


def _stored(apoptosis_run):
    return pd.read_csv(ss.run_path(apoptosis_run) / "analysis/subtree_metrics.csv")


def provenance(tables, ctx, apoptosis_run):
    analysis = ss.run_path(apoptosis_run) / "analysis"
    return dict(
        version=VERSION, pairing_run=ctx.runs["pairing"], compare_run=ctx.runs["compare"],
        apoptosis_run=apoptosis_run,
        cache_hashes={str((analysis / n).relative_to(ss.ROOT)): ss.file_hash(analysis / n)
                      for n in ("subtree_metrics.csv", "provenance.json")},
        source_hashes={s: ss.file_hash(ss.ROOT / s) for s in CODE_SOURCES},
        settings=dict(variants=list(VARIANTS), primary=PRIMARY, supporting=list(SUPPORTING), n_draws=N_DRAWS,
                      seed=SEED, zero=ZERO, pool_size=ap.POOL_SIZE, alpha=ap.ALPHA, min_edges=ss.MIN_EDGES),
        gates=ctx.gates,
        table_hashes={t: hashlib.sha256(tables[t].to_csv(index=False).encode()).hexdigest() for t in TABLES})


def write_run(out, apoptosis_run=ap.DEFAULT_RUN, **runs):
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists; runs are write-once")
    ctx = ap.load_context(**runs)
    tables = compute(ctx, _stored(apoptosis_run))
    analysis = out / "analysis"
    analysis.mkdir(parents=True)
    for name in TABLES:
        tables[name].to_csv(analysis / f"{name}.csv", index=False)
    (analysis / "provenance.json").write_text(
        json.dumps(provenance(tables, ctx, apoptosis_run), indent=2, sort_keys=True) + "\n")
    return tables


def verify_run(out):
    analysis = Path(out) / "analysis"
    stored = json.loads((analysis / "provenance.json").read_text())
    if stored["version"] != VERSION:
        raise AssertionError(f"Run version {stored['version']} != {VERSION}")
    ctx = ap.load_context(stored["pairing_run"], stored["compare_run"])
    tables = compute(ctx, _stored(stored["apoptosis_run"]))
    fresh = provenance(tables, ctx, stored["apoptosis_run"])
    changed = [k for k in stored["cache_hashes"] if stored["cache_hashes"][k] != fresh["cache_hashes"].get(k)]
    if changed:
        raise AssertionError(f"cache_hashes changed: {changed}")
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
    parser.add_argument("--apoptosis-run", default=ap.DEFAULT_RUN)
    parser.add_argument("--pairing-run", default=ss.DEFAULT_RUN)
    parser.add_argument("--compare-run", default=ss.DEFAULT_COMPARE_RUN)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = ss.run_path(args.run_id)
    if args.verify_only:
        verify_run(out)
        print(f"verified {out}")
    else:
        write_run(out, apoptosis_run=args.apoptosis_run, pairing_run=args.pairing_run, compare_run=args.compare_run)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
