"""Left-right subtree symmetry in the terminal Pareto analysis (Track D).

D0, pairing. Every embryonic terminal cell (the WormGUIDES lineage leaves) is
given a fate class and a side from its WormBase name:

- **sided**: a name whose L/R-swapped form is also an embryonic terminal name
  (ASEL/ASER, IL2DL/IL2DR). The P cells are sided by their "aka P1/2 L" alias.
- **serial**: ventral-cord motor neurons (DA, DB, DD), one class per letter
  pair, without a side.
- **unsided**: every other name, shared (MU BOD, HYP7) or unique (RIH,
  EXC CELL); programmed deaths (class ``death``) and cells without a WormBase
  name (``unannotated``).

Two cells are mirror-compatible when they share a class and, if both are
sided, sit on opposite sides. A sided cell is never compatible with an unsided
one. Two subtrees X and Y are compared in two ways:

- **fate match**: the maximum number of mirror-compatible pairs between the
  terminal cells of X and Y, ignoring lineage position, as a fraction
  2 m / (|X| + |Y|). This is what pairs the anterior AB branches, whose
  counterparts do not correspond by name.
- **position match**: X + s is paired with Y + phi(s), where phi keeps the
  lineage suffix or swaps l and r. The fraction of cells whose mirror exists
  and is compatible, under the better phi.

The sided share of the fate match is the left-right evidence. Unsided matches
(body-wall muscle, deaths) show fate agreement but not laterality, so a pair
built from them alone is not a left-right pair.

Coverage reads the canonical pooled context: which counterpart cells are
lineage leaves of the analysis, tracked in each embryo, and in the 275-edge
cohort. Nothing here recomputes a front.

Run: ``python -m terminal_pareto.subtree_symmetry [--run-id ID]``; add
``--verify-only`` to recompute an existing run and compare every table.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import io
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto import data_loader as dl  # noqa: E402

VERSION = "terminal-subtree-symmetry-v1"
DEFAULT_RUN = "pairing_20261014"
RUN_ROOT = ROOT / "terminal_pareto/output/runs/subtree_symmetry_v1"
CANONICAL = ROOT / ("terminal_pareto/output/runs/pooled_tracking_v1/"
                    "migration_candidate_20260920/analysis")
CELL_KEY = "data/2023-06-29_entropy_cell_key_V2.csv"
APOPTOTIC = "data/apoptotic_cells.txt"
LINEAGE = "data/cell_lineage.json"
PROTEIN = "data/protein/aggregated_all/s3_zscore.csv"
CODE_SOURCES = ("terminal_pareto/subtree_symmetry.py", "terminal_pareto/data_loader.py")
CANONICAL_FILES = ("analysis_manifest.json", "cell_manifest.csv", "subtree_manifest.csv")

SERIAL = ("DA", "DB", "DD")
# Sister-level founder structure (WormGUIDES names EMS's daughters MS and E).
FOUNDER_PARENT = {"AB": "P0", "P1": "P0", "EMS": "P1", "P2": "P1", "MS": "EMS",
                  "E": "EMS", "C": "P2", "P3": "P2", "D": "P3", "P4": "P3",
                  "Z2": "P4", "Z3": "P4"}
# Left-right asymmetric fates among sided pairs, from the literature.
KNOWN_ASYMMETRY = {
    "ASE": "ASEL/ASER: left-right asymmetric chemoreceptor expression and function",
    "AWC": "AWCL/AWCR: stochastic AWC-ON/AWC-OFF choice",
}
# Pairs named in the Track D brief.
BRIEF_PAIRS = (("ABpl", "ABpr"), ("ABpla", "ABpra"), ("ABplp", "ABprp"),
               ("ABplpa", "ABprpa"), ("ABplpp", "ABprpp"), ("ABplppp", "ABprppp"),
               ("MSa", "MSp"), ("MSaa", "MSpa"), ("MSap", "MSpp"))
# Author decisions (2026-10-15): strict pairing threshold, two priority pairs
# suggested by embryologists, anterior pairs exploratory, d_SD primary.
MIN_SIDED_MATCH = 0.85
MIN_EDGES = 12
PRIORITY_PAIRS = (("ABplpa", "ABprpa"), ("MSaa", "MSpa"))
EXPLORATORY_PAIRS = (("ABalp", "ABara"), ("ABalpa", "ABarap"))
PHI = {"same": str.maketrans("", ""), "swap_lr": str.maketrans("lr", "rl")}


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_path(run_id):
    return RUN_ROOT / run_id


# --- lineage --------------------------------------------------------------

def parent_of(cell):
    if cell == "P0":
        return None
    return FOUNDER_PARENT.get(cell, cell[:-1])


@dataclass
class Lineage:
    terminals: list            # embryonic terminal cells, sorted
    descendants: dict          # internal node -> terminal cells under it

    def related(self, x, y):
        return bool(set(self.descendants[x]) & set(self.descendants[y]))


def build_lineage(key):
    nodes = set(key["wormguides.lineage"].dropna())
    # WormGUIDES lists EMS, P2, P3 and P4 without lineage-suffixed daughters.
    nodes |= set(FOUNDER_PARENT) | {"P0"}
    children = collections.defaultdict(list)
    for node in nodes:
        parent = parent_of(node)
        if parent is not None:
            if parent not in nodes:
                raise ValueError(f"{node} has no parent {parent} in the lineage")
            children[parent].append(node)
    terminals = sorted(n for n in nodes if n not in children)
    descendants = {}

    def collect(node):
        if node not in children:
            return [node]
        out = [c for child in sorted(children[node]) for c in collect(child)]
        descendants[node] = sorted(out)
        return out

    collect("P0")
    if descendants["P0"] != terminals:
        raise AssertionError("Lineage is not a single tree rooted at P0")
    return Lineage(terminals, descendants)


# --- cell classes ---------------------------------------------------------

def _swap(name):
    return name[:-1] + {"L": "R", "R": "L"}[name[-1]]


def classify(cells, names, descs, deaths):
    """Class, side and kind of every embryonic terminal cell."""
    alive = {names[c] for c in cells if c not in deaths}
    rows = {}
    for cell in cells:
        name, desc = names[cell], descs[cell]
        alias = re.search(r"aka (P\d+/\d+) ([LR])\)?", desc) if isinstance(desc, str) else None
        if cell in deaths:
            row = ("death", None, "death")
        elif name.upper() == cell.upper():
            row = ("unannotated", None, "unannotated")
        elif alias:
            row = (alias.group(1), alias.group(2), "sided")
        elif name[-1] in "LR" and _swap(name) in alive:
            row = (name[:-1], name[-1], "sided")
        elif re.fullmatch(r"(%s)\d+" % "|".join(SERIAL), name):
            row = (re.match(r"[A-Z]+", name).group(), None, "serial")
        else:
            row = (name, None, "unsided")
        rows[cell] = dict(cell=cell, name=name, fate_class=row[0], side=row[1], kind=row[2])
    return rows


def compatible(a, b):
    if a["fate_class"] != b["fate_class"]:
        return False
    if a["side"] is None and b["side"] is None:
        return True
    return a["side"] is not None and b["side"] is not None and a["side"] != b["side"]


# --- matching -------------------------------------------------------------

def fate_matching(cells_x, cells_y, cls):
    """Maximum mirror-compatible matching; returns (total, sided) pair counts."""
    count_x = collections.Counter((cls[c]["fate_class"], cls[c]["side"]) for c in cells_x)
    count_y = collections.Counter((cls[c]["fate_class"], cls[c]["side"]) for c in cells_y)
    total = sided = 0
    for (fate, side), n in count_x.items():
        if side is None:
            total += min(n, count_y.get((fate, None), 0))
        else:
            m = min(n, count_y.get((fate, {"L": "R", "R": "L"}[side]), 0))
            total += m
            sided += m
    return total, sided


def position_pairs(x, y, lineage, mode):
    """Cells of X paired with Y + phi(suffix); unmatched cells map to None."""
    in_y = set(lineage.descendants.get(y, [y]))
    return {c: (y + c[len(x):].translate(PHI[mode])
                if y + c[len(x):].translate(PHI[mode]) in in_y else None)
            for c in lineage.descendants.get(x, [x])}


def position_match(x, y, lineage, cls):
    """Best phi and the fraction of X and Y cells with a compatible mirror."""
    best = None
    for mode in PHI:
        ok = 0
        for a, b in ((x, y), (y, x)):
            for c, m in position_pairs(a, b, lineage, mode).items():
                ok += m is not None and compatible(cls[c], cls[m])
        n = len(lineage.descendants[x]) + len(lineage.descendants[y])
        if best is None or ok / n > best[1]:
            best = (mode, ok / n)
    return best


def laterality(cells, cls):
    sides = collections.Counter(cls[c]["side"] for c in cells if cls[c]["side"])
    n = sides["L"] + sides["R"]
    return (sides["L"] - sides["R"]) / n if n else np.nan


def self_symmetry(cells, cls):
    """Share of the sided cells whose counterpart lies in the same subtree."""
    count = collections.Counter((cls[c]["fate_class"], cls[c]["side"]) for c in cells)
    sided = [c for c in cells if cls[c]["side"]]
    hits = sum(count.get((cls[c]["fate_class"], {"L": "R", "R": "L"}[cls[c]["side"]]), 0) > 0
               for c in sided)
    return hits / len(sided) if sided else np.nan


def pair_metrics(x, y, lineage, cls):
    cx, cy = lineage.descendants[x], lineage.descendants[y]
    total, sided = fate_matching(cx, cy, cls)
    n_sided = sum(cls[c]["side"] is not None for c in cx + cy)
    mode, pos = position_match(x, y, lineage, cls)
    return dict(
        n_terminal_x=len(cx), n_terminal_y=len(cy),
        fate_match=2 * total / (len(cx) + len(cy)),
        sided_cells=n_sided,
        sided_match=2 * sided / n_sided if n_sided else np.nan,
        sided_share=sided / total if total else np.nan,
        position_match=pos, position_phi=mode,
        laterality_x=laterality(cx, cls), laterality_y=laterality(cy, cls),
        self_symmetry_x=self_symmetry(cx, cls), self_symmetry_y=self_symmetry(cy, cls),
    )


# --- inputs ---------------------------------------------------------------

@dataclass
class Context:
    lineage: Lineage
    cells: dict                # cell -> class row
    descs: dict
    canonical: dict            # analysis manifest
    cell_manifest: pd.DataFrame
    subtree_manifest: pd.DataFrame
    leaves: set                # analysis (cell_lineage.json) leaves
    valid: dict                # embryo -> tracked names at the cutoff
    birth: dict                # embryo -> cell -> first time point (no cutoff)
    protein: set


def _json_leaves(node, out):
    children = node.get("children", [])
    if not children:
        out.add(dl.map_names(node["did"]))
    for child in children:
        _json_leaves(child, out)
    return out


def load_context():
    key = pd.read_csv(ROOT / CELL_KEY)
    lineage = build_lineage(key)
    named = key.drop_duplicates("wormweb.lineage").set_index("wormweb.lineage")
    missing = [c for c in lineage.terminals if c not in named.index]
    if missing:
        raise ValueError(f"Cells without a WormBase entry: {missing[:5]}")
    names = {c: named.at[c, "wormweb.cell"] for c in lineage.terminals}
    descs = {c: named.at[c, "wormweb.desc"] for c in lineage.terminals}
    deaths = set((ROOT / APOPTOTIC).read_text().split())
    cells = classify(lineage.terminals, names, descs, deaths)

    canonical = json.loads((CANONICAL / "analysis_manifest.json").read_text())
    for rel, digest in canonical["source_hashes"].items():
        if file_hash(ROOT / rel) != digest:
            raise ValueError(f"{rel} differs from the canonical context")
    valid, birth = {}, {}
    for spec in canonical["spec"]["tracking"]:
        path = ROOT / f"data/{spec['label']}/tracks.txt"
        _, names_ok = dl.load_elegans_tracking(spec["cutoff"], path=path)
        valid[spec["label"]] = set(names_ok)
        tracks = pd.read_csv(path, sep="\t", usecols=["t", "name"])
        birth[spec["label"]] = tracks.groupby("name")["t"].min().to_dict()
    protein = set(pd.read_csv(ROOT / PROTEIN, index_col=0).columns)
    return Context(
        lineage=lineage, cells=cells, descs=descs, canonical=canonical,
        cell_manifest=pd.read_csv(CANONICAL / "cell_manifest.csv"),
        subtree_manifest=pd.read_csv(CANONICAL / "subtree_manifest.csv").set_index("subtree"),
        leaves=_json_leaves(dl.load_json(ROOT / LINEAGE), set()),
        valid=valid, birth=birth, protein=protein)


# --- tables ---------------------------------------------------------------

def coverage_status(ctx, cell, embryo):
    if cell not in ctx.leaves:
        return "not_analysis_leaf"
    if cell not in ctx.protein:
        return "no_protein_profile"
    if cell not in ctx.valid[embryo]:
        return "not_tracked"
    if parent_of(cell) not in ctx.valid[embryo]:
        return "parent_not_tracked"
    return "edge"


def cell_table(ctx):
    cls = ctx.cells
    by_class = collections.defaultdict(list)
    for c, row in cls.items():
        by_class[(row["fate_class"], row["side"])].append(c)
    cohort = set(ctx.cell_manifest.loc[ctx.cell_manifest.included_profile, "terminal_cell"])
    embryos = list(ctx.valid)
    rows = []
    for c in ctx.lineage.terminals:
        row = dict(cls[c])
        row["description"] = ctx.descs[c] if isinstance(ctx.descs[c], str) else ""
        partners = (by_class[(row["fate_class"], {"L": "R", "R": "L"}[row["side"]])]
                    if row["side"] else [])
        row["counterparts"] = ";".join(partners)
        row["n_counterparts"] = len(partners)
        row["known_asymmetry"] = KNOWN_ASYMMETRY.get(row["fate_class"], "")
        for e in embryos:
            row[f"status_{e}"] = coverage_status(ctx, c, e)
            row[f"birth_{e}"] = ctx.birth[e].get(c, np.nan)
        row["in_cohort"] = c in cohort
        rows.append(row)
    return pd.DataFrame(rows)


def check_coverage_replay(ctx, cells):
    """Edge flags rebuilt here must equal the canonical cell manifest."""
    manifest = ctx.cell_manifest.set_index("terminal_cell")
    for e in ctx.valid:
        rebuilt = set(cells.loc[cells[f"status_{e}"] == "edge", "cell"])
        canonical = set(manifest.index[manifest[f"present_{e}"]])
        if rebuilt != canonical:
            raise AssertionError(f"{e}: coverage differs from the canonical manifest: "
                                 f"{sorted(rebuilt ^ canonical)[:10]}")
    if not set(cells.loc[cells.in_cohort, "cell"]) <= set(manifest.index):
        raise AssertionError("Cohort cells outside the canonical manifest")


def partner_table(ctx):
    """Best fate-matched partner of every eligible subtree among unrelated nodes."""
    eligible = list(ctx.subtree_manifest.index[ctx.subtree_manifest.eligible_min12])
    nodes = [n for n, cells in ctx.lineage.descendants.items() if len(cells) >= 2]
    rows = []
    for x in eligible:
        scores = []
        for y in nodes:
            if y == x or ctx.lineage.related(x, y):
                continue
            total, _ = fate_matching(ctx.lineage.descendants[x], ctx.lineage.descendants[y], ctx.cells)
            n = len(ctx.lineage.descendants[x]) + len(ctx.lineage.descendants[y])
            scores.append((2 * total / n, -abs(len(ctx.lineage.descendants[y])
                                               - len(ctx.lineage.descendants[x])), y))
        if not scores:
            continue
        scores.sort(reverse=True)
        best, runner = scores[0], scores[1] if len(scores) > 1 else (np.nan, 0, "")
        rows.append(dict(subtree=x, partner=best[2], runner_up=runner[2],
                         runner_up_fate_match=runner[0],
                         **pair_metrics(x, best[2], ctx.lineage, ctx.cells)))
    out = pd.DataFrame(rows)
    best = dict(zip(out.subtree, out.partner))
    out["mutual"] = [best.get(p) == s for s, p in zip(out.subtree, out.partner)]
    return out


def candidate_pairs(ctx, partners):
    pairs = {tuple(sorted(p)) for p in BRIEF_PAIRS}
    for r in partners.itertuples():
        if r.mutual:
            pairs.add(tuple(sorted((r.subtree, r.partner))))
    cohort = set(ctx.cell_manifest.loc[ctx.cell_manifest.included_profile, "terminal_cell"])
    manifest = ctx.cell_manifest.set_index("terminal_cell")
    rows = []
    for x, y in sorted(pairs):
        row = dict(subtree_x=x, subtree_y=y,
                   source="brief" if (x, y) in {tuple(sorted(p)) for p in BRIEF_PAIRS} else "mutual_best",
                   **pair_metrics(x, y, ctx.lineage, ctx.cells))
        cx, cy = ctx.lineage.descendants[x], ctx.lineage.descendants[y]
        for tag, cells, node in (("x", cx, x), ("y", cy, y)):
            row[f"analysis_leaves_{tag}"] = sum(c in ctx.leaves for c in cells)
            for e in ctx.valid:
                row[f"edges_{e}_{tag}"] = int(manifest.loc[manifest.index.intersection(cells), f"present_{e}"].sum())
            row[f"cohort_{tag}"] = sum(c in cohort for c in cells)
            row[f"eligible_{tag}"] = bool(ctx.subtree_manifest.eligible_min12.get(node, False))
        for tag, cells in (("x", cx), ("y", cy)):
            for e in ctx.valid:
                status = collections.Counter(coverage_status(ctx, c, e) for c in cells)
                for reason in ("not_tracked", "parent_not_tracked"):
                    row[f"{reason}_{e}_{tag}"] = status[reason]
            row[f"no_protein_{tag}"] = sum(c in ctx.leaves and c not in ctx.protein for c in cells)
        pos, _ = mirror_pairs(x, y, ctx.lineage, ctx.cells, row["position_phi"])
        row["position_pairs"] = len(pos)
        row["cohort_position_pairs"] = sum(a in cohort and b in cohort for a, b in pos)
        row["discordant_position_pairs"] = sum((a in cohort) != (b in cohort) for a, b in pos)
        for e in ctx.valid:
            # Birth of X's cell minus its mirror's, in time points; > 0: X side born later.
            lag = [ctx.birth[e][a] - ctx.birth[e][b] for a, b in pos
                   if a in ctx.birth[e] and b in ctx.birth[e]]
            row[f"birth_lag_median_{e}"] = float(np.median(lag)) if lag else np.nan
        total, sided = fate_matching([c for c in cx if c in cohort], [c for c in cy if c in cohort], ctx.cells)
        n = row["cohort_x"] + row["cohort_y"]
        row["cohort_fate_match"] = 2 * total / n if n else np.nan
        row["cohort_ratio"] = (min(row["cohort_x"], row["cohort_y"]) / max(row["cohort_x"], row["cohort_y"])
                               if max(row["cohort_x"], row["cohort_y"]) else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def mirror_pairs(x, y, lineage, cls, mode):
    """Position-compatible pairs under phi, then the fate matching of the rest.

    Returns (position pairs, cells left without any compatible partner).
    Pairing by position first never lowers the matching: both cells of a
    position pair share a class, so every class count drops on both sides.
    """
    pos = [(c, m) for c, m in position_pairs(x, y, lineage, mode).items()
           if m is not None and compatible(cls[c], cls[m])]
    paired = {c for pair in pos for c in pair}
    rest_x = [c for c in lineage.descendants[x] if c not in paired]
    rest_y = [c for c in lineage.descendants[y] if c not in paired]
    return pos, _fate_leftover(rest_x, rest_y, cls) | _fate_leftover(rest_y, rest_x, cls)


def mismatch_table(ctx, pairs):
    """Cells of each candidate pair whose lineage mirror is not compatible.

    ``mirror`` is the cell at the mirrored lineage position under the pair's
    phi. ``fate_unmatched`` marks cells with no compatible partner anywhere in
    the other subtree once position pairs are taken.
    """
    cls = ctx.cells
    rows = []
    for r in pairs.itertuples():
        x, y = r.subtree_x, r.subtree_y
        pos, leftover = mirror_pairs(x, y, ctx.lineage, cls, r.position_phi)
        paired = {c for pair in pos for c in pair}
        for a, b in ((x, y), (y, x)):
            for c, m in position_pairs(a, b, ctx.lineage, r.position_phi).items():
                if c in paired:
                    continue
                rows.append(dict(subtree_x=x, subtree_y=y, side_of_pair=a, cell=c,
                                 name=cls[c]["name"], fate_class=cls[c]["fate_class"], kind=cls[c]["kind"],
                                 mirror=m or "", mirror_name=cls[m]["name"] if m else "",
                                 mirror_kind=cls[m]["kind"] if m else "",
                                 fate_unmatched=c in leftover))
    return pd.DataFrame(rows)


def _fate_leftover(cells_x, cells_y, cls):
    """Cells of X beyond the matched count of their class (deterministic order)."""
    count_y = collections.Counter((cls[c]["fate_class"], cls[c]["side"]) for c in cells_y)
    used = collections.Counter()
    out = set()
    for c in sorted(cells_x):
        fate, side = cls[c]["fate_class"], cls[c]["side"]
        target = (fate, {"L": "R", "R": "L"}[side] if side else None)
        if used[(fate, side)] < count_y.get(target, 0):
            used[(fate, side)] += 1
        else:
            out.add(c)
    return out


# --- run ------------------------------------------------------------------

TABLES = ("cells", "subtree_partners", "candidate_pairs", "pair_mismatches")


def compute(ctx):
    cells = cell_table(ctx)
    check_coverage_replay(ctx, cells)
    partners = partner_table(ctx)
    pairs = candidate_pairs(ctx, partners)
    return dict(cells=cells, subtree_partners=partners, candidate_pairs=pairs,
                pair_mismatches=mismatch_table(ctx, pairs))


def provenance(tables):
    sources = list(CODE_SOURCES) + [CELL_KEY, APOPTOTIC, LINEAGE, PROTEIN]
    sources += [f"data/{e}/tracks.txt" for e in ("embryo1", "embryo2", "embryo3")]
    return dict(
        version=VERSION,
        canonical_analysis=str(CANONICAL.relative_to(ROOT)),
        canonical_hashes={f: file_hash(CANONICAL / f) for f in CANONICAL_FILES},
        source_hashes={s: file_hash(ROOT / s) for s in sources},
        settings=dict(serial_classes=list(SERIAL), position_phi=list(PHI),
                      known_asymmetry=KNOWN_ASYMMETRY, brief_pairs=[list(p) for p in BRIEF_PAIRS]),
        table_hashes={t: hashlib.sha256(tables[t].to_csv(index=False).encode()).hexdigest()
                      for t in TABLES},
    )


def write_run(out):
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists; runs are write-once")
    tables = compute(load_context())
    analysis = out / "analysis"
    analysis.mkdir(parents=True)
    for name in TABLES:
        tables[name].to_csv(analysis / f"{name}.csv", index=False)
    (analysis / "provenance.json").write_text(json.dumps(provenance(tables), indent=2, sort_keys=True) + "\n")
    return tables


def verify_run(out):
    """Recompute every table and compare with the stored run."""
    analysis = Path(out) / "analysis"
    stored = json.loads((analysis / "provenance.json").read_text())
    if stored["version"] != VERSION:
        raise AssertionError(f"Run version {stored['version']} != {VERSION}")
    tables = compute(load_context())
    fresh = provenance(tables)
    for group in ("canonical_hashes", "source_hashes"):
        changed = [k for k in stored[group] if stored[group][k] != fresh[group].get(k)]
        if changed:
            raise AssertionError(f"{group} changed: {changed}")
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        pd.testing.assert_frame_equal(saved, _roundtrip(tables[name]),
                                      check_exact=False, rtol=0, atol=1e-12)
    return tables


def _roundtrip(frame):
    return pd.read_csv(io.StringIO(frame.to_csv(index=False)), keep_default_na=False, na_values=[""])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    out = run_path(args.run_id)
    if args.verify_only:
        verify_run(out)
        print(f"verified {out}")
    else:
        write_run(out)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
