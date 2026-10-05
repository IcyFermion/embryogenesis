"""Cross-species terminal comparison by subtree and by cell type (add-on to the 187-edge run).

python -m terminal_pareto.cross_species_partitions [--run-id ID]

Reads the frozen cross-species terminal cache (``cross_species_analysis``) and
writes ``<run>/analysis/partitions/``. Nothing frozen is modified.

* Subtrees: every lineage node with at least ``MIN_SUBTREE`` matched terminal
  edges (the Figure 4/5 rule, applied to the 187 matched edges). Each subtree
  is its own assignment problem over its terminals and their parent slots,
  analysed exactly as the whole cohort (``cross_species_analysis.analyze``:
  exact first-cousin scaling within the subtree, 301 weights, endpoint ties).
* Cell types: the C. elegans annotation joined on lineage names, assumed for
  the matched C. briggsae cells (author decision 2026-10-05; caption caveat).
  Types with at least ``MIN_TYPE`` cells get a within-type front (assignments
  never cross types); smaller types are merged into ``other`` for the
  retention decomposition only.
* Retention by type along each cached global front (no new solves): a child
  is retained when its assigned parent slot belongs to its natural parent.

Deterministic: no random draws (null means are exact).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from terminal_pareto import cross_species_analysis as tcs
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.front_coordinates import DegenerateEndpointSpan

VERSION = "terminal-cross-species-partitions-1"
MIN_SUBTREE = 12
MIN_TYPE = 15
SUBDIR = "analysis/partitions"


def subtrees(terminals, tree_index, min_cells=MIN_SUBTREE):
    """{node: terminal positions} for lineage nodes with >= min_cells matched terminals, largest first."""
    members = {}
    for i, cell in enumerate(terminals):
        for ancestor in tree_index[cell]["ancestors"]:
            members.setdefault(ancestor, []).append(i)
    kept = {k: np.asarray(v) for k, v in members.items() if len(v) >= min_cells}
    return dict(sorted(kept.items(), key=lambda kv: (-len(kv[1]), tree_index[kv[0]]["depth"], kv[0])))


def type_partition(types, min_cells=MIN_TYPE):
    """Display types: large types kept, the rest merged into 'other'."""
    counts = pd.Series(types).value_counts()
    return np.array([t if counts[t] >= min_cells and t != "other" else "other" for t in types])


def child_retention(parents, permutations):
    """(weights, children): child j retained when its assigned slot is a natural-parent slot of j."""
    # Slot i receives child permutation[i], so child j sits in slot argsort(permutation)[j].
    slots = np.argsort(permutations, axis=1)
    return parents[slots] == parents[None, :]


NULL_FIELDS = ("d_NP", "null_D1", "null_D2", "null_mean_travel", "null_mean_state")


def analyze_block(travel, state, parents, terminals, positions, tree_index, *, intervals, identity,
                  full_random_fallback=False):
    """``tcs.analyze`` on one block; returns (curve, metrics) or (None, reason).

    With ``full_random_fallback`` a block whose first-cousin null has zero
    variance is scaled by its full-random null instead (as the CE within-type
    analysis falls back from degenerate cousin nulls); its first-cousin null
    mean and d_NP are then not reported.
    """
    ix = np.asarray(positions)
    groups = pe.build_ancestor_groups(terminals[ix], tree_index, 2)
    blocks = travel[np.ix_(ix, ix)], state[np.ix_(ix, ix)], parents[ix]
    null_model = "first_cousin"
    try:
        try:
            result = tcs.analyze(*blocks, groups, intervals=intervals, identity=identity)
        except DegenerateEndpointSpan:
            raise
        except ValueError:
            if not full_random_fallback:
                raise
            null_model = "full_random_fallback"
            result = tcs.analyze(*blocks, [list(range(len(ix)))], intervals=intervals, identity=identity)
    except DegenerateEndpointSpan:
        return None, "degenerate endpoint span"
    except ValueError as error:
        return None, str(error)
    metrics = {k: v for k, v in result["metrics"].items() if np.isscalar(v)}
    if null_model != "first_cousin":
        metrics.update({k: np.nan for k in NULL_FIELDS})
    metrics.update(null_model=null_model, first_cousin_groups=len(groups), first_cousin_shuffled=int(sum(map(len, groups))))
    return result["curve"], metrics


def identity_for(record):
    sources = [Path(__file__), tcs.ROOT / "terminal_pareto/cross_species_analysis.py",
               tcs.ROOT / "terminal_pareto/pareto_engine.py", tcs.ROOT / "terminal_pareto/lineage_metrics.py",
               tcs.ROOT / "terminal_pareto/front_coordinates.py", tcs.ROOT / "terminal_pareto/subtree_analysis.py"]
    return dict(version=VERSION, analysis_id=record["analysis_id"], intervals=record["settings"]["intervals"],
                min_subtree=MIN_SUBTREE, min_type=MIN_TYPE, numpy_version=np.__version__,
                source_hashes={str(p.relative_to(tcs.ROOT)): tcs.file_hash(p) for p in sources})


def compute(run):
    report = tcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    intervals = record["settings"]["intervals"]
    manifest = pd.read_csv(analysis / "cell_manifest.csv")
    with np.load(analysis / "assignments_and_costs.npz", allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    terminals, parents = arrays["terminals"], arrays["parents"]
    if not np.array_equal(manifest.terminal.to_numpy(), terminals):
        raise ValueError("Cell manifest and cached terminal order differ")
    tree_index = lm.build_lineage_tree_index(tcs.dl.load_json(tcs.ROOT / "data/cell_lineage.json"))
    blocks = subtrees(terminals, tree_index)
    types = type_partition(manifest.cell_type.to_numpy())
    type_blocks = {t: np.flatnonzero(types == t) for t in sorted(set(types)) if t != "other"}
    fronts, metrics, retention = [], [], []
    saved_metrics = pd.read_csv(analysis / "metrics.csv")
    for geometry in tcs.GEOMETRIES:
        for config in tcs.PRIMARY:
            key = f"{geometry}__{config}"
            travel, state = arrays[f"{key}__travel"], arrays[f"{key}__state"]
            for kind, partition in (("subtree", blocks), ("cell_type", type_blocks)):
                for name, ix in partition.items():
                    curve, result = analyze_block(travel, state, parents, terminals, ix, tree_index,
                                                  intervals=intervals, identity=f"{record['analysis_id']}:{key}:{kind}:{name}",
                                                  full_random_fallback=kind == "cell_type")
                    row = dict(geometry=geometry, config=config, kind=kind, name=name, n=len(ix),
                               depth=tree_index[name]["depth"] if kind == "subtree" else np.nan)
                    if curve is None:
                        metrics.append(dict(row, degenerate=True, reason=result))
                        continue
                    metrics.append(dict(row, degenerate=False, reason="", **result))
                    fronts.append(curve.assign(geometry=geometry, config=config, kind=kind, name=name))
            # The whole-cohort block must reproduce the cached comparison exactly.
            root = next(m for m in reversed(metrics) if m["kind"] == "subtree" and m["n"] == len(terminals))
            saved = saved_metrics[(saved_metrics.geometry == geometry) & (saved_metrics.config == config)
                                  & (saved_metrics.cohort == "base")].iloc[0]
            np.testing.assert_allclose([root[k] for k in ("u_L", "d_LP", "d_NP")], [saved[k] for k in ("u_L", "d_LP", "d_NP")],
                                       rtol=1e-10, atol=1e-12)
            # Retention by type along the cached global sweep.
            front = pd.read_csv(analysis / "fronts.csv")
            front = front[(front.geometry == geometry) & (front.config == config) & (front.cohort == "base")].sort_values("sweep_index")
            kept = child_retention(parents, arrays[f"{key}__permutations"])
            np.testing.assert_allclose(kept.mean(axis=1), front.edge_retention.to_numpy(), atol=1e-14)
            for t in sorted(set(types)):
                cols = types == t
                retention.append(pd.DataFrame(dict(geometry=geometry, config=config, cell_type=t, n=int(cols.sum()),
                    sweep_index=front.sweep_index.to_numpy(), u=front.u.to_numpy(),
                    retention=kept[:, cols].mean(axis=1))))
    return (report, record, pd.concat(fronts, ignore_index=True), pd.DataFrame(metrics),
            pd.concat(retention, ignore_index=True))


def build(run_id=tcs.DEFAULT_RUN):
    run = tcs.run_path(run_id)
    out = run / SUBDIR
    _, record, fronts, metrics, retention = compute(run)
    identity = identity_for(record)
    files = dict(zip(("fronts.csv", "metrics.csv", "type_retention.csv"), (fronts, metrics, retention)))
    if (out / "validation.json").exists():
        saved = json.loads((out / "validation.json").read_text())
        if saved["identity"] != identity:
            raise ValueError("Partition cache has another identity; use a new run")
        for name, frame in files.items():
            if tcs.file_hash(out / name) != saved["hashes"][name]:
                raise ValueError(f"Changed partition file: {name}")
            previous = pd.read_csv(out / name)
            numeric = frame.select_dtypes(include="number").columns
            np.testing.assert_allclose(frame[numeric].to_numpy(float), previous[numeric].to_numpy(float),
                                       rtol=1e-11, atol=1e-12, equal_nan=True)
        return out
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in files.items():
        frame.to_csv(out / name, index=False)
    (out / "validation.json").write_text(json.dumps(dict(identity=identity, hashes={n: tcs.file_hash(out / n) for n in files}),
                                                    indent=2, sort_keys=True) + "\n")
    return out


def load(run_id=tcs.DEFAULT_RUN):
    """Read-only: verify identity, sources and file hashes."""
    run = tcs.run_path(run_id)
    out = run / SUBDIR
    saved = json.loads((out / "validation.json").read_text())
    record = json.loads((run / "analysis/provenance.json").read_text())
    if saved["identity"]["version"] != VERSION or saved["identity"]["analysis_id"] != record["analysis_id"]:
        raise ValueError("Partition cache belongs to another analysis")
    for path, sha in saved["identity"]["source_hashes"].items():
        if tcs.file_hash(tcs.ROOT / path) != sha:
            raise ValueError(f"Changed source since partitions were computed: {path}")
    for name, sha in saved["hashes"].items():
        if tcs.file_hash(out / name) != sha:
            raise ValueError(f"Changed partition file: {name}")
    frames = {name.removesuffix(".csv"): pd.read_csv(out / name) for name in saved["hashes"]}
    return dict(**frames, identity=saved["identity"], record=record,
                identity_sha=hashlib.sha256(json.dumps(saved["identity"], sort_keys=True).encode()).hexdigest())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-id", default=tcs.DEFAULT_RUN)
    args = parser.parse_args(argv)
    out = build(args.run_id)
    metrics = pd.read_csv(out / "metrics.csv")
    print(metrics.groupby(["kind", "geometry", "config"]).agg(blocks=("name", "size"),
          degenerate=("degenerate", "sum")).to_string(), flush=True)
    print(f"Partitions ready: {out}", flush=True)


if __name__ == "__main__":
    main()
