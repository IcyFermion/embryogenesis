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


# --- D1: compare paired subtrees ------------------------------------------

DEFAULT_COMPARE_RUN = "compare_20261015"
TRAVELS = ("pooled", "embryo1", "embryo2", "embryo3")
INTERVALS = 300
REPRODUCTION_TOL = 1e-9
TINY_SPAN_SD = 1.0        # an endpoint span under one null SD: u and d_CP are unstable
REF_DEPTH = 1             # reference pairs: each member within this many generations
REF_SIZE = (0.67, 1.5)    # ... and edge counts within these factors of the pair's
COMPARE_METRICS = ("d_SD", "natural_state_z", "natural_travel_z", "u_L", "d_LP", "d_NP",
                   "max_edge_retention", "travel_span_sd", "state_span_sd")
CANONICAL_COLUMNS = (("u_L", "u_lineage_lp"), ("d_LP", "d_lp"), ("d_NP", "d_np"),
                     ("max_edge_retention", "max_er"))


def depth(node):
    d = 0
    while node != "P0":
        node, d = parent_of(node), d + 1
    return d


@dataclass
class CompareContext:
    lineage: Lineage
    children: list
    parents: np.ndarray
    state: np.ndarray
    travel: dict               # TRAVELS -> matrix
    tree_index: dict
    canonical_metrics: pd.DataFrame
    pairing: dict              # pairing tables
    pairing_run: str


def load_compare_context(pairing_run=DEFAULT_RUN):
    from terminal_pareto import lineage_metrics as lm
    from terminal_pareto.analysis_context import _sha256_array

    manifest = json.loads((CANONICAL / "analysis_manifest.json").read_text())
    for rel, digest in manifest["source_hashes"].items():
        if file_hash(ROOT / rel) != digest:
            raise ValueError(f"{rel} differs from the canonical context")
    with np.load(CANONICAL / "analysis_context.npz", allow_pickle=False) as arrays:
        travel = {"pooled": np.asarray(arrays["travel_matrix"])}
        for e in TRAVELS[1:]:
            travel[e] = np.asarray(arrays[f"travel_component_{e}"])
        state = np.asarray(arrays["expression_matrix"])
        children = arrays["children"].tolist()
        parents = np.asarray(arrays["natural_parents"].tolist())
    for key, value in (("travel", travel["pooled"]), ("expression", state),
                       *((f"travel_component_{e}", travel[e]) for e in TRAVELS[1:])):
        if _sha256_array(value) != manifest["matrix_hashes"][key]:
            raise ValueError(f"Canonical {key} matrix does not match its manifest")
    analysis = run_path(pairing_run) / "analysis"
    stored = json.loads((analysis / "provenance.json").read_text())
    pairing = {}
    for name in TABLES:
        text = (analysis / f"{name}.csv").read_text()
        if hashlib.sha256(text.encode()).hexdigest() != stored["table_hashes"][name]:
            raise ValueError(f"Pairing table {name} does not match its provenance")
        pairing[name] = pd.read_csv(io.StringIO(text))
    return CompareContext(
        lineage=build_lineage(pd.read_csv(ROOT / CELL_KEY)), children=children, parents=parents,
        state=state, travel=travel, tree_index=lm.build_lineage_tree_index(dl.load_json(ROOT / LINEAGE)),
        canonical_metrics=pd.read_csv(CANONICAL / "ce_subtree_canonical_metrics.csv").set_index("subtree"),
        pairing=pairing, pairing_run=pairing_run)


def subtree_edges(ctx, node):
    """Cohort edge indices under ``node``, in canonical edge order."""
    if node == "P0":
        return np.arange(len(ctx.children))
    under = set(ctx.lineage.descendants[node])
    return np.array([i for i, c in enumerate(ctx.children) if c in under])


def null_support(cells, parents, groups):
    """Cells whose first-cousin group spans at least two natural parents.

    A group of sisters alone has one parent, so shuffling it never changes a
    cost; only the other groups give the null its variance.
    """
    return sum(len(g) for g in groups if len({parents[i] for i in g}) > 1)


def subtree_front(ctx, node, travel_key):
    """Front metrics for one subtree, through ``cross_species_analysis.analyze``."""
    from terminal_pareto import pareto_engine as pe
    from terminal_pareto.cross_species_analysis import analyze

    idx = subtree_edges(ctx, node)
    cells = [ctx.children[i] for i in idx]
    groups = pe.build_ancestor_groups(cells, ctx.tree_index, 2)
    travel = ctx.travel[travel_key][np.ix_(idx, idx)]
    state = ctx.state[np.ix_(idx, idx)]
    result = analyze(travel, state, ctx.parents[idx], groups, intervals=INTERVALS,
                     identity=f"{VERSION}:{node}:{travel_key}")
    m, curve = result["metrics"], result["curve"]
    sx, ss_ = m["travel_sigma"], m["state_sigma"]
    t, s = curve["travel"].to_numpy(), curve["cell_state"].to_numpy()
    row = dict(subtree=node, travel=travel_key, n=len(idx), depth=depth(node),
               null_support=null_support(cells, ctx.parents[idx], groups),
               **{k: m[k] for k in ("u_L", "d_LP", "d_NP", "max_edge_retention", "retention_at_nearest",
                                    "natural_travel", "natural_cell_state", "null_mean_travel",
                                    "null_mean_state", "travel_sigma", "state_sigma")})
    row.update(
        d_SD=float(np.min(np.hypot((t - m["natural_travel"]) / sx, (s - m["natural_cell_state"]) / ss_))),
        natural_travel_z=(m["natural_travel"] - m["null_mean_travel"]) / sx,
        natural_state_z=(m["natural_cell_state"] - m["null_mean_state"]) / ss_,
        travel_optimum=float(t[-1]), state_optimum=float(s[0]),
        travel_span_sd=float(t[0] - t[-1]) / sx, state_span_sd=float(s[-1] - s[0]) / ss_)
    row["tiny_span"] = min(row["travel_span_sd"], row["state_span_sd"]) < TINY_SPAN_SD
    curve = curve.assign(subtree=node, travel_key=travel_key)
    return row, curve


def subtree_metrics(ctx):
    eligible = [n for n in ctx.canonical_metrics.index]
    rows, curves = [], []
    for node in eligible:
        for key in TRAVELS:
            row, curve = subtree_front(ctx, node, key)
            rows.append(row)
            curves.append(curve)
    return pd.DataFrame(rows), pd.concat(curves, ignore_index=True)


def check_reproduction(ctx, metrics):
    """Pooled fronts must reproduce the canonical subtree table."""
    pooled = metrics[metrics.travel == "pooled"].set_index("subtree")
    worst = 0.0
    for node, row in ctx.canonical_metrics.iterrows():
        if pooled.at[node, "n"] != row["n"]:
            raise AssertionError(f"{node}: {pooled.at[node, 'n']} edges, canonical {row['n']}")
        for ours, theirs in CANONICAL_COLUMNS:
            delta = abs(pooled.at[node, ours] - row[theirs])
            worst = max(worst, delta)
            if not delta <= REPRODUCTION_TOL:
                raise AssertionError(f"{node} {ours}: {pooled.at[node, ours]} vs canonical {row[theirs]}")
    return worst


def pair_tiers(ctx):
    """Candidate pairs with the author's tiers (decisions of 2026-10-15)."""
    pairs = ctx.pairing["candidate_pairs"].copy()
    partners = ctx.pairing["subtree_partners"].set_index("subtree")
    key = lambda x, y: tuple(sorted((x, y)))
    priority = {key(*p) for p in PRIORITY_PAIRS}
    exploratory = {key(*p) for p in EXPLORATORY_PAIRS}
    tiers = []
    for r in pairs.itertuples():
        mutual = (partners.partner.get(r.subtree_x) == r.subtree_y
                  and partners.partner.get(r.subtree_y) == r.subtree_x)
        strict = (mutual and r.sided_match >= MIN_SIDED_MATCH
                  and min(r.cohort_x, r.cohort_y) >= MIN_EDGES)
        pair = key(r.subtree_x, r.subtree_y)
        tiers.append("priority" if pair in priority and strict else
                     "strict" if strict else "exploratory" if pair in exploratory else "excluded")
    pairs["tier"] = tiers
    return pairs


def reference_pairs(ctx, metrics, pairs):
    """Unrelated, non-mirror pairs of eligible subtrees matched to each kept pair."""
    pooled = metrics[metrics.travel == "pooled"].set_index("subtree")
    nodes = list(pooled.index)
    mirror = {tuple(sorted(p)) for p in zip(pairs.subtree_x, pairs.subtree_y)}
    pool = [(a, b) for i, a in enumerate(nodes) for b in nodes[i + 1:]
            if not ctx.lineage.related(a, b) and tuple(sorted((a, b))) not in mirror]
    rows = []
    for r in pairs[pairs.tier != "excluded"].itertuples():
        lo, hi = sorted((pooled.at[r.subtree_x, "n"], pooled.at[r.subtree_y, "n"]))
        d = (depth(r.subtree_x) + depth(r.subtree_y)) / 2
        for a, b in pool:
            na, nb = sorted((pooled.at[a, "n"], pooled.at[b, "n"]))
            if (abs(depth(a) - d) <= REF_DEPTH and abs(depth(b) - d) <= REF_DEPTH
                    and na >= REF_SIZE[0] * lo and nb <= REF_SIZE[1] * hi):
                rows.append(dict(subtree_x=r.subtree_x, subtree_y=r.subtree_y, ref_a=a, ref_b=b))
    return pd.DataFrame(rows)


def pair_comparison(metrics, pairs, references):
    """Within-pair differences against matched random pairs and the embryo floor.

    ``reference_percentile``: share of matched random pairs with a smaller
    absolute difference, plus half the ties. ``embryo_floor``: the largest
    difference between two tracking embryos for either partner.
    """
    pooled = metrics[metrics.travel == "pooled"].set_index("subtree")
    per_embryo = metrics[metrics.travel != "pooled"]
    rows = []
    for r in pairs[pairs.tier != "excluded"].itertuples():
        ref = references[(references.subtree_x == r.subtree_x) & (references.subtree_y == r.subtree_y)]
        for metric in COMPARE_METRICS:
            x, y = pooled.at[r.subtree_x, metric], pooled.at[r.subtree_y, metric]
            diff = abs(x - y)
            ref_diff = np.abs(pooled.loc[ref.ref_a, metric].to_numpy() - pooled.loc[ref.ref_b, metric].to_numpy())
            floor = 0.0
            for node in (r.subtree_x, r.subtree_y):
                values = per_embryo.loc[per_embryo.subtree == node, metric].to_numpy()
                floor = max(floor, float(values.max() - values.min()))
            rows.append(dict(
                subtree_x=r.subtree_x, subtree_y=r.subtree_y, tier=r.tier, metric=metric,
                value_x=x, value_y=y, abs_difference=diff, n_reference=len(ref_diff),
                reference_median=float(np.median(ref_diff)) if len(ref_diff) else np.nan,
                reference_percentile=(float((ref_diff < diff).mean() + 0.5 * (ref_diff == diff).mean())
                                      if len(ref_diff) else np.nan),
                embryo_floor=floor,
                tiny_span_x=bool(pooled.at[r.subtree_x, "tiny_span"]),
                tiny_span_y=bool(pooled.at[r.subtree_y, "tiny_span"]),
                null_support_x=int(pooled.at[r.subtree_x, "null_support"]),
                null_support_y=int(pooled.at[r.subtree_y, "null_support"])))
    return pd.DataFrame(rows)


def compute_compare(ctx):
    metrics, curves = subtree_metrics(ctx)
    worst = check_reproduction(ctx, metrics)
    pairs = pair_tiers(ctx)
    references = reference_pairs(ctx, metrics, pairs)
    tables = dict(pair_tiers=pairs, subtree_metrics=metrics, reference_pairs=references,
                  pair_comparison=pair_comparison(metrics, pairs, references), front_curves=curves)
    return tables, dict(reproduction_max_abs_delta=worst)


# --- D1: cell-level checks -------------------------------------------------

DEFAULT_CELLS_RUN = "cells_20261015"
SELECTED = "expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv"


@dataclass
class CellContext:
    lineage: Lineage
    cells: pd.DataFrame        # pairing cells table, indexed by cell
    cls: dict
    tiers: pd.DataFrame
    profiles: pd.DataFrame     # top-20 protein profiles, cells in rows
    positions: dict            # embryo -> cell -> xyz at its last time <= cutoff
    canonical: dict
    runs: dict


def load_cell_context(pairing_run=DEFAULT_RUN, compare_run=DEFAULT_COMPARE_RUN):
    from scipy.spatial.distance import cosine

    pairing = run_path(pairing_run) / "analysis"
    compare = run_path(compare_run) / "analysis"
    for folder, names in ((pairing, TABLES), (compare, ("pair_tiers",))):
        stored = json.loads((folder / "provenance.json").read_text())["table_hashes"]
        for name in names:
            if file_hash(folder / f"{name}.csv") != stored[name]:
                raise ValueError(f"{folder / name}.csv does not match its provenance")
    cells = pd.read_csv(pairing / "cells.csv").set_index("cell")
    cls = {c: dict(cell=c, name=r["name"], fate_class=r["fate_class"],
                   side=r["side"] if isinstance(r["side"], str) else None, kind=r["kind"])
           for c, r in cells.iterrows()}
    canonical = json.loads((CANONICAL / "analysis_manifest.json").read_text())
    profiles = dl.load_protein_expression()[dl.load_prot_sel()].astype(float)
    # Gate: profiles reproduce the canonical natural-edge cell-state costs.
    with np.load(CANONICAL / "analysis_context.npz", allow_pickle=False) as arrays:
        state = np.asarray(arrays["expression_matrix"]).diagonal()
        edges = zip(arrays["children"].tolist(), arrays["natural_parents"].tolist())
        delta = max(abs(cosine(profiles.loc[p], profiles.loc[c]) - s) for (c, p), s in zip(edges, state))
    if delta > REPRODUCTION_TOL:
        raise AssertionError(f"Profiles do not reproduce the canonical state costs ({delta})")
    positions = {}
    for spec in canonical["spec"]["tracking"]:
        xyz, _ = dl.load_elegans_tracking(spec["cutoff"], path=ROOT / f"data/{spec['label']}/tracks.txt")
        positions[spec["label"]] = {k: np.asarray(v, dtype=float) for k, v in xyz.items()}
    return CellContext(build_lineage(pd.read_csv(ROOT / CELL_KEY)), cells, cls,
                       pd.read_csv(compare / "pair_tiers.csv"), profiles, positions, canonical,
                       dict(pairing=pairing_run, compare=compare_run, state_gate=delta))


def path_length(a, b):
    """Divisions separating two cells through their last common ancestor."""
    chain = lambda c: [c] + ([] if c == "P0" else chain(parent_of(c)))
    ca, cb = chain(a), chain(b)
    common = next(x for x in ca if x in cb)
    return ca.index(common) + cb.index(common)


def counterpart_pairs(ctx, x, y, phi):
    """Position-compatible pairs, then sided name counterparts across X and Y."""
    pos, _ = mirror_pairs(x, y, ctx.lineage, ctx.cls, phi)
    out = [(a, b, "position") for a, b in pos]
    seen = {c for pair in pos for c in pair}
    in_y = set(ctx.lineage.descendants[y])
    for a in ctx.lineage.descendants[x]:
        partners = str(ctx.cells.at[a, "counterparts"]).split(";")
        if a in seen or ctx.cells.at[a, "n_counterparts"] != 1 or partners[0] not in in_y or partners[0] in seen:
            continue
        out.append((a, partners[0], "name"))
    return out


def counterpart_expression(ctx):
    """Cosine distance of counterparts against other cells of the same fine type.

    ``*_percentile``: share of same-type cells closer to cell a than its
    counterpart. ``far`` keeps only same-type cells at least as distant in the
    lineage as the counterpart, which removes lineage memory from the reference.
    """
    from scipy.spatial.distance import cdist

    desc = ctx.cells["description"].fillna("")
    have = [c for c in ctx.profiles.index if c in ctx.cells.index]
    rows = []
    for r in ctx.tiers[ctx.tiers.tier != "excluded"].itertuples():
        for a, b, basis in counterpart_pairs(ctx, r.subtree_x, r.subtree_y, r.position_phi):
            if a not in ctx.profiles.index or b not in ctx.profiles.index or ctx.cls[a]["kind"] in ("death", "unannotated"):
                continue
            same = [c for c in have if desc[a] and desc[c] == desc[a] and c not in (a, b)]
            far = [c for c in same if path_length(a, c) >= path_length(a, b)]
            d_all = cdist(ctx.profiles.loc[[a]], ctx.profiles.loc[[b] + same], "cosine")[0]
            d, ref = d_all[0], d_all[1:]
            ref_far = ref[[same.index(c) for c in far]]
            rows.append(dict(subtree_x=r.subtree_x, subtree_y=r.subtree_y, tier=r.tier, basis=basis,
                             cell_a=a, name_a=ctx.cls[a]["name"], cell_b=b, name_b=ctx.cls[b]["name"],
                             path_length=path_length(a, b), distance=d,
                             n_same_type=len(same), same_type_percentile=float((ref < d).mean()) if len(ref) else np.nan,
                             same_type_median=float(np.median(ref)) if len(ref) else np.nan,
                             n_far=len(far), far_percentile=float((ref_far < d).mean()) if len(ref_far) else np.nan))
    return pd.DataFrame(rows)


def mirror_planes(ctx):
    """Per embryo, the plane that best reflects left cells onto right cells.

    Normal: mean right-minus-left vector over unique sided counterpart pairs in
    the cohort; the plane passes through the mean midpoint.
    """
    cohort = ctx.cells.index[ctx.cells.in_cohort]
    pairs = []
    for c in cohort:
        partner = str(ctx.cells.at[c, "counterparts"])
        if ctx.cls[c]["side"] == "L" and ctx.cells.at[c, "n_counterparts"] == 1 and partner in cohort:
            pairs.append((c, partner))
    rows, planes = [], {}
    for e, xyz in ctx.positions.items():
        left = np.array([xyz[a] for a, _ in pairs])
        right = np.array([xyz[b] for _, b in pairs])
        v = right - left
        normal = v.mean(axis=0) / np.linalg.norm(v.mean(axis=0))
        point = ((left + right) / 2).mean(axis=0)
        planes[e] = (normal, point)
        residual = np.linalg.norm(reflect(left, normal, point) - right, axis=1)
        rows.append(dict(embryo=e, n_pairs=len(pairs), normal_x=normal[0], normal_y=normal[1], normal_z=normal[2],
                         median_separation=float(np.median(np.linalg.norm(v, axis=1))),
                         median_axis_cosine=float(np.median(v @ normal / np.linalg.norm(v, axis=1))),
                         median_residual=float(np.median(residual)), max_residual=float(residual.max())))
    return pd.DataFrame(rows), planes


def reflect(points, normal, point):
    points = np.asarray(points, dtype=float)
    return points - 2 * ((points - point) @ normal)[..., None] * normal


def _angle(u, v):
    return float(np.degrees(np.arccos(np.clip(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))


def counterpart_travel(ctx, planes):
    """Natural-edge travel of counterparts in each embryo's mirrored frame.

    ``mirror_angle``: angle between a's reflected displacement (parent to
    child) and b's. ``reference_angle``: median of the same angle over every
    non-counterpart pairing of a with cohort cells of the partner subtree.
    """
    cohort = set(ctx.cells.index[ctx.cells.in_cohort])
    rows = []
    for r in ctx.tiers[ctx.tiers.tier != "excluded"].itertuples():
        pairs = [(a, b, basis) for a, b, basis in counterpart_pairs(ctx, r.subtree_x, r.subtree_y, r.position_phi)
                 if a in cohort and b in cohort]
        others = [c for c in ctx.lineage.descendants[r.subtree_y] if c in cohort]
        for e, xyz in ctx.positions.items():
            normal, point = planes[e]
            step = lambda c: xyz[c] - xyz[parent_of(c)]
            mirrored = lambda c: reflect(step(c)[None, :] + point, normal, point)[0] - point
            for a, b, basis in pairs:
                ref = [_angle(mirrored(a), step(c)) for c in others if c != b]
                rows.append(dict(subtree_x=r.subtree_x, subtree_y=r.subtree_y, tier=r.tier, embryo=e, basis=basis,
                                 cell_a=a, name_a=ctx.cls[a]["name"], cell_b=b, name_b=ctx.cls[b]["name"],
                                 length_a=float(np.linalg.norm(step(a))), length_b=float(np.linalg.norm(step(b))),
                                 mirror_angle=_angle(mirrored(a), step(b)),
                                 reference_angle=float(np.median(ref)) if ref else np.nan,
                                 position_residual=float(np.linalg.norm(reflect(xyz[a], normal, point) - xyz[b]))))
    out = pd.DataFrame(rows)
    # Measurement floor: the same edge's length across embryos.
    span = lambda col, cell: out.groupby(cell)[col].transform(lambda v: v.max() - v.min())
    out["length_floor_a"] = span("length_a", "cell_a")
    out["length_floor_b"] = span("length_b", "cell_b")
    return out


def cell_summary(expression, travel):
    rows = []
    for (x, y, tier), g in expression.groupby(["subtree_x", "subtree_y", "tier"], sort=False):
        t = travel[(travel.subtree_x == x) & (travel.subtree_y == y)]
        rows.append(dict(
            subtree_x=x, subtree_y=y, tier=tier, n_expression_pairs=len(g),
            median_distance=g.distance.median(), median_same_type_percentile=g.same_type_percentile.median(),
            share_closest_of_type=float((g.same_type_percentile[g.n_same_type > 0] == 0).mean()),
            median_far_percentile=g.far_percentile.median(),
            n_travel_pairs=t.cell_a.nunique(),
            median_abs_length_difference=float((t.length_a - t.length_b).abs().median()) if len(t) else np.nan,
            median_length_floor=float(pd.concat([t.length_floor_a, t.length_floor_b]).median()) if len(t) else np.nan,
            median_mirror_angle=t.mirror_angle.median(), median_reference_angle=t.reference_angle.median(),
            median_position_residual=t.position_residual.median()))
    return pd.DataFrame(rows)


def compute_cells(ctx):
    expression = counterpart_expression(ctx)
    planes, fitted = mirror_planes(ctx)
    travel = counterpart_travel(ctx, fitted)
    return dict(counterpart_expression=expression, mirror_planes=planes, counterpart_travel=travel,
                cell_summary=cell_summary(expression, travel))


# --- run ------------------------------------------------------------------

TABLES = ("cells", "subtree_partners", "candidate_pairs", "pair_mismatches")
COMPARE_TABLES = ("pair_tiers", "subtree_metrics", "reference_pairs", "pair_comparison", "front_curves")
CELL_TABLES = ("counterpart_expression", "mirror_planes", "counterpart_travel", "cell_summary")


def compute(ctx):
    cells = cell_table(ctx)
    check_coverage_replay(ctx, cells)
    partners = partner_table(ctx)
    pairs = candidate_pairs(ctx, partners)
    return dict(cells=cells, subtree_partners=partners, candidate_pairs=pairs,
                pair_mismatches=mismatch_table(ctx, pairs))


def _table_hashes(tables, names):
    return {t: hashlib.sha256(tables[t].to_csv(index=False).encode()).hexdigest() for t in names}


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
        table_hashes=_table_hashes(tables, TABLES),
    )


def compare_provenance(tables, ctx, gates):
    files = CANONICAL_FILES + ("analysis_context.npz", "ce_subtree_canonical_metrics.csv")
    return dict(
        version=VERSION, stage="compare", pairing_run=ctx.pairing_run,
        canonical_analysis=str(CANONICAL.relative_to(ROOT)),
        canonical_hashes={f: file_hash(CANONICAL / f) for f in files},
        source_hashes={s: file_hash(ROOT / s) for s in CODE_SOURCES + (LINEAGE, CELL_KEY)},
        settings=dict(intervals=INTERVALS, tiny_span_sd=TINY_SPAN_SD, reference_depth=REF_DEPTH,
                      reference_size=list(REF_SIZE), min_sided_match=MIN_SIDED_MATCH, min_edges=MIN_EDGES,
                      priority_pairs=[list(p) for p in PRIORITY_PAIRS],
                      exploratory_pairs=[list(p) for p in EXPLORATORY_PAIRS],
                      reproduction_tol=REPRODUCTION_TOL),
        gates=gates,
        table_hashes=_table_hashes(tables, COMPARE_TABLES),
    )


def cells_provenance(tables, ctx):
    sources = CODE_SOURCES + (CELL_KEY, LINEAGE, PROTEIN, SELECTED)
    sources += tuple(f"data/{e}/tracks.txt" for e in ctx.positions)
    return dict(
        version=VERSION, stage="cells", pairing_run=ctx.runs["pairing"], compare_run=ctx.runs["compare"],
        canonical_hashes={f: file_hash(CANONICAL / f) for f in CANONICAL_FILES + ("analysis_context.npz",)},
        source_hashes={s: file_hash(ROOT / s) for s in sources},
        settings=dict(expression="top-20 protein, cosine", tiers="all but excluded",
                      cutoffs={t["label"]: t["cutoff"] for t in ctx.canonical["spec"]["tracking"]}),
        gates=dict(state_cost_max_abs_delta=ctx.runs["state_gate"]),
        table_hashes=_table_hashes(tables, CELL_TABLES),
    )


def _stage(stage, pairing_run=DEFAULT_RUN, compare_run=DEFAULT_COMPARE_RUN):
    """Compute one stage; returns (tables, table names, provenance record)."""
    if stage == "pairing":
        tables = compute(load_context())
        return tables, TABLES, provenance(tables)
    if stage == "compare":
        ctx = load_compare_context(pairing_run)
        tables, gates = compute_compare(ctx)
        return tables, COMPARE_TABLES, compare_provenance(tables, ctx, gates)
    ctx = load_cell_context(pairing_run, compare_run)
    tables = compute_cells(ctx)
    return tables, CELL_TABLES, cells_provenance(tables, ctx)


def write_run(out, stage="pairing", pairing_run=DEFAULT_RUN, compare_run=DEFAULT_COMPARE_RUN):
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists; runs are write-once")
    tables, names, record = _stage(stage, pairing_run, compare_run)
    analysis = out / "analysis"
    analysis.mkdir(parents=True)
    for name in names:
        tables[name].to_csv(analysis / f"{name}.csv", index=False)
    (analysis / "provenance.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return tables


def verify_run(out):
    """Recompute every table and compare with the stored run."""
    analysis = Path(out) / "analysis"
    stored = json.loads((analysis / "provenance.json").read_text())
    if stored["version"] != VERSION:
        raise AssertionError(f"Run version {stored['version']} != {VERSION}")
    tables, names, fresh = _stage(stored.get("stage", "pairing"), stored.get("pairing_run", DEFAULT_RUN),
                                  stored.get("compare_run", DEFAULT_COMPARE_RUN))
    # Inputs must be identical. Code may grow (later stages share this module),
    # so code hashes are only reported; the exact table replay below is the check.
    for group in ("canonical_hashes", "source_hashes"):
        changed = [k for k in stored[group] if stored[group][k] != fresh[group].get(k)
                   and k not in CODE_SOURCES]
        if changed:
            raise AssertionError(f"{group} changed: {changed}")
    code_changed = [k for k in CODE_SOURCES if stored["source_hashes"].get(k) != fresh["source_hashes"][k]]
    if code_changed:
        print(f"note: code changed since the run ({', '.join(code_changed)}); tables replayed")
    for name in names:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        pd.testing.assert_frame_equal(saved, _roundtrip(tables[name]),
                                      check_exact=False, rtol=0, atol=1e-12)
    return tables


def _roundtrip(frame):
    return pd.read_csv(io.StringIO(frame.to_csv(index=False)), keep_default_na=False, na_values=[""])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--stage", choices=("pairing", "compare", "cells"), default="pairing")
    parser.add_argument("--run-id")
    parser.add_argument("--pairing-run", default=DEFAULT_RUN)
    parser.add_argument("--compare-run", default=DEFAULT_COMPARE_RUN)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    run_id = args.run_id or {"pairing": DEFAULT_RUN, "compare": DEFAULT_COMPARE_RUN,
                             "cells": DEFAULT_CELLS_RUN}[args.stage]
    out = run_path(run_id)
    if args.verify_only:
        verify_run(out)
        print(f"verified {out}")
    else:
        write_run(out, args.stage, args.pairing_run, args.compare_run)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
