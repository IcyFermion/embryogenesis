"""Developmental stage of each layer, and the Large et al. (2025) expression hourglass by stage (Track E, E4).

Two read-only tables for the E4 figures, so the figures need no ad hoc data:

* **Birth stage per layer.** For every cell of the pooled full tree and of the
  cross-species partial forest: its contraction round, layer (9 - round), birth
  time and the number of cells alive at that time, from the Large et al. cell
  table (``CellTable_20241027.txt``, C. elegans). Summaries per tree and round:
  median and quartiles of cells alive at birth.
* **Between-species JSD by stage.** Recomputed from the released pseudobulk
  tables on all shared gene names, as in Large et al. Fig 5E: progenitors are
  binned by birth time (AB8 ... AB512 boundaries), terminal types form the
  600-cell stage (mean over embryo-time bins). Each stage is placed at the
  number of cells alive at its start. The recomputed medians sit about
  0.02-0.06 below the published ones; the shape, with its minimum at the
  200-cell stage, is reproduced.

Run from the repository root in the dev environment::

    python -m full_tree_pareto.layer_stages --run-id stages_20261020
    python -m full_tree_pareto.layer_stages --run-id stages_20261020 --verify-only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from full_tree_pareto import cross_species_analysis as fcs
from full_tree_pareto import layer_features as lf
from full_tree_pareto import pooled_analysis as pa

ROOT = lf.ROOT
VERSION = "layer-stages-e4-1"
DEFAULT_RUN = "stages_20261020"
LARGE = "data/c_briggsae/science.adu8249/shared"
CELL_TABLE = f"{LARGE}/CellTable_20241027.txt"
TPM = {(kind, species): f"{LARGE}/tpm_matrix{kind}_c_{species}.txt"
       for kind in ("", "_time") for species in ("elegans", "briggsae")}
CROSS_RUN = fcs.run_path(fcs.DEFAULT_RUN)
STAGE_EDGES = (56, 75, 103, 137, 181, 232, 330)   # AB8 ... AB512 birth times (min) in the cell table
STAGES = ("15", "28", "50", "100", "200", "350", "600")
TOP = 9
CODE_SOURCES = ("full_tree_pareto/layer_stages.py", "full_tree_pareto/pooled_analysis.py",
                "full_tree_pareto/cross_species_analysis.py")
TABLES = ("cells", "layers", "jsd_cells", "jsd_stages")


def cell_table():
    table = pd.read_csv(ROOT / CELL_TABLE, sep="\t").drop_duplicates("Lineage").set_index("Lineage")
    birth = pd.to_numeric(table.br_time, errors="coerce")
    death = pd.to_numeric(table.d_time, errors="coerce")
    return table, birth, death


def alive_at(times, birth, death):
    """Cells of the cell table alive at each time (born at or before, not yet divided or dead)."""
    b, d = birth.to_numpy(float), death.to_numpy(float)
    return np.array([int(((b <= t) & ((d > t) | np.isnan(d))).sum()) if np.isfinite(t) else -1 for t in times])


def stage_of(birth_time):
    if not np.isfinite(birth_time):
        return None
    k = int(np.searchsorted(STAGE_EDGES, birth_time, side="right")) - 1
    return "8" if k < 0 else STAGES[k]


def tree_layers():
    """(tree, names, layers) for the pooled full tree and the cross-species partial forest."""
    ctx = pa.build_context()
    record = json.loads((CROSS_RUN / "analysis/provenance.json").read_text())
    with np.load(CROSS_RUN / "analysis/matrices.npz", allow_pickle=False) as data:
        parents = np.asarray(data["parents"], dtype=int)
    return [("pooled", ctx.names, ctx.layers), ("partial", record["names"], pa.contraction_layers(parents))]


def birth_tables():
    _, birth, death = cell_table()
    rows = []
    for tree, names, layers in tree_layers():
        for r, children in enumerate(layers, 1):
            for i in children:
                name = names[i]
                rows.append(dict(tree=tree, round=r, layer=TOP - r, cell=name,
                                 birth_time=float(birth.get(name, np.nan))))
    cells = pd.DataFrame(rows)
    cells["alive_at_birth"] = alive_at(cells.birth_time.to_numpy(), birth, death)
    cells.loc[cells.alive_at_birth < 0, "alive_at_birth"] = np.nan
    timed = cells.dropna(subset=["alive_at_birth"])
    layers = timed.groupby(["tree", "round", "layer"], sort=False).alive_at_birth.agg(
        cells="size", alive_median="median", alive_q25=lambda x: x.quantile(.25),
        alive_q75=lambda x: x.quantile(.75)).reset_index()
    edges = cells.groupby(["tree", "round"], sort=False).size().rename("edges")
    layers = layers.merge(edges, on=["tree", "round"])
    return cells, layers


def jsd_tables():
    """Between-species JSD per progenitor and terminal type, and by stage."""
    table, birth, death = cell_table()

    def per_column(kind):
        ce = pd.read_csv(ROOT / TPM[kind, "elegans"], sep="\t", index_col=0)
        cb = pd.read_csv(ROOT / TPM[kind, "briggsae"], sep="\t", index_col=0)
        genes, columns = ce.index.intersection(cb.index), ce.columns.intersection(cb.columns)
        ce, cb = ce.loc[genes, columns].fillna(0), cb.loc[genes, columns].fillna(0)
        return pd.Series({c: jensenshannon(ce[c].to_numpy(), cb[c].to_numpy(), base=2) for c in columns}), len(genes)

    pseudo, n_genes = per_column("")
    timed, _ = per_column("_time")
    timed = timed[timed.index.str.contains(":")]
    terminal = timed.groupby(timed.index.str.split("_").str[0]).mean()
    rows = [dict(cell=c, kind="progenitor", stage=stage_of(float(birth.get(c, np.nan))), jsd=float(v))
            for c, v in pseudo.items() if c in table.index and table.Cell.get(c) == "progenitor"]
    rows += [dict(cell=c, kind="terminal type", stage="600", jsd=float(v)) for c, v in terminal.items()]
    cells = pd.DataFrame(rows).dropna(subset=["stage"])
    order = ["8", *STAGES]
    starts = dict(zip(STAGES, alive_at(np.array(STAGE_EDGES, dtype=float), birth, death)))
    stages = cells.groupby("stage").jsd.agg(n="size", median="median", q25=lambda x: x.quantile(.25),
                                            q75=lambda x: x.quantile(.75)).reindex(order).dropna(how="all").reset_index()
    stages["alive_at_start"] = stages.stage.map(starts)
    stages["shared_genes"] = n_genes
    return cells.reset_index(drop=True), stages


def compute():
    cells, layers = birth_tables()
    jsd_cells, jsd_stages = jsd_tables()
    return dict(cells=cells, layers=layers, jsd_cells=jsd_cells, jsd_stages=jsd_stages)


# ---- run folder ---------------------------------------------------------------------

def source_hashes():
    paths = (*CODE_SOURCES, CELL_TABLE, *TPM.values(),
             str((CROSS_RUN / "analysis/provenance.json").relative_to(ROOT)))
    return {name: lf.file_hash(ROOT / name) for name in paths}


def write_run(out):
    analysis = out / "analysis"
    record_path = analysis / "provenance.json"
    if record_path.exists():
        raise FileExistsError(f"Numerical run already exists at {out}; use a new --run-id")
    started = time.time()
    sources = source_hashes()
    results = compute()
    analysis.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        results[name].to_csv(analysis / f"{name}.csv", index=False)
    settings = dict(stage_edges=list(STAGE_EDGES), stages=list(STAGES), layer="9 - round")
    identity = hashlib.sha256(json.dumps(dict(version=VERSION, settings=settings, sources=sources),
                                         sort_keys=True).encode()).hexdigest()
    record = dict(version=VERSION, analysis_id=identity, settings=settings, sources=sources,
                  elapsed_seconds=round(time.time() - started, 1),
                  files={p.name: lf.file_hash(p) for p in sorted(analysis.iterdir()) if p.is_file()})
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def validate_run(out, *, check_sources=True):
    analysis = out / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported layer-stages cache")
    for name, expected in record["files"].items():
        if lf.file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch: {name}")
    if check_sources:
        for name, expected in record["sources"].items():
            if lf.file_hash(ROOT / name) != expected:
                raise ValueError(f"Input/code changed since analysis: {name}; create a new numerical run")
    fresh = compute()
    for name in TABLES:
        saved = pd.read_csv(analysis / f"{name}.csv", keep_default_na=False, na_values=[""])
        lf._compare(fresh[name], saved, name)
    return dict(analysis_id=record["analysis_id"], tables=len(TABLES),
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
