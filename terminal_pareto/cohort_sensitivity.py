"""Terminal cohort sensitivity: cells that end embryogenesis undivided but are internal in the reference tree.

Track D, stage D7 (author, 2026-11-08). The terminal cohort keeps a cell only if
it is a leaf of ``data/cell_lineage.json`` and it and its parent are tracked in
all three embryos at the stage cutoffs and have protein data. The reference
tree continues past embryogenesis (for example the L1 intestinal nuclear
divisions), so cells that never divide again in the embryo can be internal
nodes and fall out of the cohort. At the cutoffs the E lineage has 16 cells,
and only Ealpa and Earpa are reference leaves.

This add-on adds every such cell and reruns the published global terminal
sweep. A cell is added when it is usable in all three embryos (tracked at the
cutoff, with protein data), none of its reference children is usable, and none
of its reference children appears anywhere in the three tracking files (they
divide only after tracking ends). Cells whose children appear after the cutoff
are left out, as the stage-matched design intends.

Everything else follows the primary analysis: pooled 3D travel (each embryo's
matrix divided by its natural-lineage total on the cohort), top-20 protein
cosine cell state, exact first-cousin scaling, 301 weights
(``cross_species_analysis.analyze``). A reproduction gate first rebuilds the
primary 275-edge matrices and metrics from the same code.

Run from the repository root in the dev environment::

    python -m terminal_pareto.cohort_sensitivity --run-id embryonic_terminals_20261108
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from terminal_pareto import data_loader as dl
from terminal_pareto import feature_sensitivity as fs
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.cross_species_analysis import analyze, pool_matrices, spatial_matrix

VERSION = "cohort-sensitivity-d7-1"
ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "terminal_pareto/output/runs/cohort_sensitivity_v1"
PRIMARY_METRICS = fs.PRIMARY_ANALYSIS / "ce_subtree_canonical_metrics.csv"
FOUNDERS = ("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "D", "Z2", "Z3")
SUBTREE_FLOOR = 24
TOL = 1e-9
SOURCES = ("data/cell_lineage.json", fs.PROTEIN_ZSCORE, fs.PROTEIN_TOP20,
           *(str(Path(path).relative_to(ROOT)) for _, path, _ in dl.CE_REPLICATES),
           "terminal_pareto/cohort_sensitivity.py", "terminal_pareto/cross_species_analysis.py")


def founder(cell):
    return next((f for f in FOUNDERS if cell.startswith(f)), "other")


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cohorts():
    """Primary 275 edges and the edges added by the embryonic-terminal rule, with the per-embryo coordinates."""
    lineage = dl.load_json(str(ROOT / "data/cell_lineage.json"))
    tree_index = lm.build_lineage_tree_index(lineage)
    protein = dl.load_protein_expression(str(ROOT / fs.PROTEIN_ZSCORE))
    children = {}
    for name, entry in tree_index.items():
        if entry["ancestors"]:
            children.setdefault(entry["ancestors"][0], []).append(name)
    xyz, usable, tracked = {}, None, set()
    for label, path, cutoff in dl.CE_REPLICATES:
        coordinates, valid = dl.load_elegans_tracking(cutoff, path)
        xyz[label] = {name: np.asarray(value, float) for name, value in coordinates.items()}
        usable = set(valid) if usable is None else usable & set(valid)
        tracked |= set(pd.read_csv(path, sep="\t")["name"].dropna())
    usable &= set(protein.index)
    primary, _ = fs.load_primary(tree_index)
    frontier = sorted(n for n in usable if not any(k in usable for k in children.get(n, ())))
    leaves = {n for n in frontier if not children.get(n)}
    if leaves != {child for child, _ in primary.terms}:
        raise ValueError("The reference-leaf rule no longer reproduces the primary cohort")
    rows = []
    for cell in frontier:
        kids = children.get(cell, [])
        parent = tree_index[cell]["ancestors"][0]
        if parent not in usable:
            raise ValueError(f"Parent of {cell} is not usable")
        status = ("reference leaf" if not kids else
                  "divides after tracking ends" if not any(k in tracked for k in kids) else
                  "divides after the cutoff")
        rows.append(dict(cell=cell, parent=parent, founder=founder(cell), status=status,
                         reference_children=";".join(sorted(kids))))
    table = pd.DataFrame(rows)
    added = [tuple(r) for r in table.loc[table.status == "divides after tracking ends", ["cell", "parent"]].values]
    return dict(tree_index=tree_index, protein=protein, xyz=xyz, primary=primary, table=table,
                terms={"primary": list(primary.terms), "embryonic": list(primary.terms) + added})


def matrices(data, terms):
    """Pooled 3D travel (natural-lineage totals of these edges) and top-20 cosine cell state."""
    components = {label: spatial_matrix(xyz, terms, "raw3d") for label, xyz in data["xyz"].items()}
    travel = pool_matrices(components, {label: float(np.trace(m)) for label, m in components.items()})
    features = dl.load_prot_sel(str(ROOT / fs.PROTEIN_TOP20))
    state = fs.state_matrix(data["protein"], list(features), terms)
    return travel, state


def evaluate(data, name, terms):
    travel, state = matrices(data, terms)
    parents = np.asarray([p for _, p in terms])
    groups = pe.build_ancestor_groups([c for c, _ in terms], data["tree_index"], 2)
    result = analyze(travel, state, parents, groups, intervals=fs.INTERVALS, identity=f"{VERSION}:{name}")
    m, curve = result["metrics"], result["curve"]
    sx, ss = m["travel_sigma"], m["state_sigma"]
    t, s = curve["travel"].to_numpy(), curve["cell_state"].to_numpy()
    row = dict(cohort=name, n=m["n"], u_L=m["u_L"], d_LP=m["d_LP"], d_CP=m["d_NP"],
               max_lineage_similarity=m["max_edge_retention"], similarity_travel=m["retention_travel"],
               similarity_state=m["retention_state"],
               natural_travel_z=(m["natural_travel"] - m["null_mean_travel"]) / sx,
               natural_state_z=(m["natural_cell_state"] - m["null_mean_state"]) / ss,
               d_LP_null_sd=float(np.min(np.hypot((t - m["natural_travel"]) / sx, (s - m["natural_cell_state"]) / ss))))
    row.update({f"n_{f}": sum(founder(c) == f for c, _ in terms) for f in FOUNDERS})
    return row, travel, state


def subtree_counts(data):
    """Terminal edges per named subtree in both cohorts (floor of the subtree figure: 24)."""
    counts = {}
    for name, terms in data["terms"].items():
        for child, _ in terms:
            for ancestor in data["tree_index"][child]["ancestors"]:
                counts.setdefault(ancestor, {}).setdefault(name, 0)
                counts[ancestor][name] += 1
    table = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    table = table[(table >= SUBTREE_FLOOR).any(axis=1)].sort_values("embryonic", ascending=False)
    return table.rename_axis("subtree").reset_index()


def run(run_id):
    if not run_id or Path(run_id).name != run_id:
        raise ValueError("Run identity must be a single directory name")
    out = RUN_ROOT / run_id / "analysis"
    if (out / "provenance.json").exists():
        raise FileExistsError(f"Run {run_id} exists; use a new --run-id")
    data = cohorts()
    rows, gate = [], {}
    for name, terms in data["terms"].items():
        row, travel, state = evaluate(data, name, terms)
        rows.append(row)
        if name == "primary":  # reproduction gate against the published primary context and Figure 5 metrics
            gate["travel"] = float(np.max(np.abs(travel - data["primary"].travel["raw3d"])))
            _, cached_state = fs.load_primary(data["tree_index"])
            gate["state"] = float(np.max(np.abs(state - cached_state)))
            p0 = pd.read_csv(PRIMARY_METRICS).set_index("subtree").loc["P0"]
            gate["d_LP"] = abs(row["d_LP"] - p0.d_lp)
            gate["u_L"] = abs(row["u_L"] - p0.u_lineage_lp)
            if max(gate.values()) > TOL:
                raise ValueError(f"Primary cohort does not reproduce: {gate}")
    out.mkdir(parents=True)
    files = {"metrics.csv": pd.DataFrame(rows), "cells.csv": data["table"], "subtrees.csv": subtree_counts(data)}
    for name, frame in files.items():
        frame.to_csv(out / name, index=False)
    record = dict(version=VERSION, reproduction=gate, settings=dict(intervals=fs.INTERVALS, floor=SUBTREE_FLOOR),
                  sources={name: file_hash(ROOT / name) for name in SOURCES},
                  files={name: file_hash(out / name) for name in files})
    (out / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")
    return out, rows, gate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)
    out, rows, gate = run(args.run_id)
    print(json.dumps(dict(out=str(out), reproduction=gate), indent=2))
    print(pd.DataFrame(rows).set_index("cohort").T.to_string())


if __name__ == "__main__":
    main()
