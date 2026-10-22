"""One apoptotic cell at a time in the posterior ABpl family (Track D4c).

D4b showed that the deaths' endpoints cause the ABplp/ABprp difference. This
stage stores which death carries it and whether the same treatment ever moves
a right partner off its front.

For each subtree of ``SUBTREES`` and each of its deaths d, under
``drop_deaths`` and ``birth_endpoint``:

- ``only``: d alone is treated;
- ``all_but``: every other death of the subtree is treated (skipped when d is
  the only one).

``none`` (no death treated) and ``all`` (every death) must reproduce D4's
subtree table. The endpoint table pairs every death in the family with its
mirror cell (l <-> r at the fourth letter), embryo by embryo.

Readouts (fixed in the Track D brief before the run; descriptive):
``single_cell`` (ABplpppapp alone puts ABplp on its front and the other deaths
do not, under both variants), ``right_unmoved`` (no treatment moves a right
subtree off its front) and ``endpoint_offset`` (ABplpppapp's travel exceeds
its mirror's in every embryo).

Run: ``python -m terminal_pareto.apoptosis_single_deaths [--run-id ID]``;
``--verify-only`` replays every table.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

from terminal_pareto import apoptosis_sensitivity as ap
from terminal_pareto import apoptosis_subtree_control as sc
from terminal_pareto import subtree_symmetry as ss

VERSION = "terminal-apoptosis-single-deaths-v1"
DEFAULT_RUN = "apoptosis_single_20261022"
CODE_SOURCES = ("terminal_pareto/apoptosis_single_deaths.py", "terminal_pareto/apoptosis_subtree_control.py") \
    + ap.CODE_SOURCES
LEFT = ("ABpl", "ABplp", "ABplpp", "ABplppp")
VARIANTS = ap.VARIANTS[1:]
CELL = "ABplpppapp"
SUBTREE = "ABplp"
ZERO = sc.ZERO


def mirror(name):
    """Left-right counterpart in the posterior AB lineage: ABpl... <-> ABpr..."""
    if not name.startswith(("ABpl", "ABpr")):
        raise ValueError(name)
    return name[:3] + {"l": "r", "r": "l"}[name[3]] + name[4:]


SUBTREES = LEFT + tuple(mirror(s) for s in LEFT)


def side(name):
    return "left" if name[3] == "l" else "right"


def deaths_in(ctx, node):
    idx = sc.subtree_index(ctx, node)
    return [int(i) for i in idx[ctx.dead[idx]]]


def single(ctx):
    """Front metrics of each subtree with no, all, one or all-but-one of its deaths treated."""
    names = ctx.compare.children
    rows = []
    for node in SUBTREES:
        dead = deaths_in(ctx, node)
        sets = [("none", None, [])] + [("all", None, dead)]
        for d in dead:
            sets.append(("only", names[d], [d]))
            if len(dead) > 1:
                sets.append(("all_but", names[d], [i for i in dead if i != d]))
        for treatment, death, cells in sets:
            variants = ("as_is",) if treatment == "none" else VARIANTS
            for variant in variants:
                row = sc.front(ctx, node, variant, cells, f"{VERSION}:{node}:{treatment}:{death}:{variant}")
                rows.append(dict(subtree=node, side=side(node), mirror=mirror(node), n_deaths=len(dead),
                                 treatment=treatment, death=death, variant=variant, n_treated=len(cells), **row))
    return pd.DataFrame(rows)


def check_d4(ctx, table, stored):
    """``none`` and ``all`` rows reproduce D4's pooled measured subtree table."""
    theirs = stored.query("cost == 'measured' and travel == 'pooled'").set_index(["subtree", "variant"])
    mine = table[table.treatment.isin(["none", "all"])].set_index(["subtree", "variant"])
    if not (theirs.loc[mine.index, "status"] == mine.status).all():
        raise AssertionError("Subtree statuses do not reproduce D4")
    ok = mine.status == "ok"
    delta = float(np.max(np.abs(mine.loc[ok, "d_SD"] - theirs.loc[mine.index[ok], "d_SD"])))
    ctx.gates["d4_d_SD_max_abs_delta"] = delta
    if delta > 1e-12:
        raise AssertionError("Subtree d_SD does not reproduce D4")


def endpoints(ctx):
    """Each death of the family next to its mirror cell, per embryo."""
    names = list(ctx.compare.children)
    index = {c: i for i, c in enumerate(names)}
    family = sorted({names[d] for node in SUBTREES for d in deaths_in(ctx, node)})
    rows = []
    for cell in family:
        other = mirror(cell)
        j = index.get(other)
        for e in ss.TRAVELS[1:]:
            row = dict(cell=cell, side=side(cell), mirror=other, embryo=e,
                       mirror_in_cohort=j is not None, mirror_dies=bool(j is not None and ctx.dead[j]))
            for prefix, k in (("", index[cell]), ("mirror_", j)):
                if k is None:
                    continue
                t = ctx.tracks[e].loc[names[k]]
                row.update({f"{prefix}born": int(t.born),
                            f"{prefix}frames_to_endpoint": int(min(t.end, t.cut) - t.born),
                            f"{prefix}travel": float(ctx.components[e][k, k]),
                            f"{prefix}birth_to_endpoint":
                                float(np.linalg.norm(ctx.child_xyz[e][k] - ctx.birth_xyz[e][k]))})
            rows.append(row)
    return pd.DataFrame(rows)


def readouts(table, ends):
    """The three readouts fixed in the brief."""
    t = table.set_index(["subtree", "treatment", "death", "variant"])
    on = lambda key: t.loc[key, "status"] == "ok" and t.loc[key, "d_SD"] <= ZERO
    single_cell = all(on((SUBTREE, "only", CELL, v)) and not on((SUBTREE, "all_but", CELL, v)) for v in VARIANTS)
    right = table[(table.side == "right") & (table.status == "ok")]
    right_unmoved = bool((right.d_SD <= ZERO).all())
    cell = ends[ends.cell == CELL]
    offset = bool((cell.travel > cell.mirror_travel).all())
    return pd.DataFrame([
        dict(readout="single_cell", holds=single_cell,
             detail=f"{SUBTREE}: {CELL} only -> on front; all other deaths -> off front (both variants)"),
        dict(readout="right_unmoved", holds=right_unmoved,
             detail=f"{len(right)} right-subtree fronts, max d_SD {right.d_SD.max():.3g}"),
        dict(readout="endpoint_offset", holds=offset,
             detail="travel " + ", ".join(f"{a:.1f}/{b:.1f}" for a, b in zip(cell.travel, cell.mirror_travel))
                    + " um (left/right) in embryos 1-3"),
    ])


TABLES = ("single", "endpoints", "readouts")


def compute(ctx, stored):
    table = single(ctx)
    check_d4(ctx, table, stored)
    ends = endpoints(ctx)
    return dict(single=table, endpoints=ends, readouts=readouts(table, ends))


def provenance(tables, ctx, apoptosis_run):
    analysis = ss.run_path(apoptosis_run) / "analysis"
    return dict(
        version=VERSION, pairing_run=ctx.runs["pairing"], compare_run=ctx.runs["compare"],
        apoptosis_run=apoptosis_run,
        cache_hashes={str((analysis / n).relative_to(ss.ROOT)): ss.file_hash(analysis / n)
                      for n in ("subtree_metrics.csv", "provenance.json")},
        source_hashes={s: ss.file_hash(ss.ROOT / s) for s in CODE_SOURCES},
        settings=dict(subtrees=list(SUBTREES), variants=list(VARIANTS), cell=CELL, subtree=SUBTREE, zero=ZERO,
                      min_edges=ss.MIN_EDGES),
        gates=ctx.gates,
        table_hashes={t: hashlib.sha256(tables[t].to_csv(index=False).encode()).hexdigest() for t in TABLES})


def write_run(out, apoptosis_run=ap.DEFAULT_RUN, **runs):
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists; runs are write-once")
    ctx = ap.load_context(**runs)
    tables = compute(ctx, sc._stored(apoptosis_run))
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
    tables = compute(ctx, sc._stored(stored["apoptosis_run"]))
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
