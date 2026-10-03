"""Full-tree-pipeline cache adapters (read-only; never run solvers or nulls)."""

from __future__ import annotations

import json

import pandas as pd

from full_tree_pareto import cross_species_analysis as fcs
from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput
from publication.provenance import sha256

CROSS_SPECIES_FAMILY = "full-tree-cross-species"
CROSS_SPECIES_INPUTS = ("provenance.json", "fronts.csv", "metrics.csv", "null_clouds.csv")


def cross_species(run_id: str = fcs.DEFAULT_RUN) -> FrontComparisonInput:
    """Terminal-anchored partial forest (485 cells, 454 edges); Euclidean molecular distance."""
    run = fcs.run_path(run_id)
    if not (run / "analysis/provenance.json").is_file():
        raise FileNotFoundError(
            f"No partial-forest comparison cache at {run}; create it with "
            f"`python -m full_tree_pareto.cross_species_analysis --run-id {run_id}`")
    report = fcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    fronts = pd.read_csv(analysis / "fronts.csv")
    metrics = pd.read_csv(analysis / "metrics.csv")
    fronts = fronts[fronts.cohort == "base"].reset_index(drop=True)
    metrics = metrics[metrics.cohort == "base"].rename(columns=CANONICAL_FIELDS).reset_index(drop=True)
    settings = record["settings"]
    return FrontComparisonInput(
        family=CROSS_SPECIES_FAMILY, run=run, analysis_id=report["analysis_id"],
        edges=int(record["edges"]), configs=tuple(fcs.PRIMARY), geometries=tuple(fcs.GEOMETRIES),
        fronts=fronts, metrics=metrics, null_mean="analytic first-cousin",
        caption=dict(edges=int(record["edges"]), nodes=int(record["cohort_size"]),
                     terminals=int(record["terminal_count"]), roots=int(record["roots"]),
                     round_edges=list(record["round_edges"]), weights=settings["intervals"] + 1,
                     dense_weights=settings["dense_intervals"] + 1, draws=settings["draws"]),
        validation=report,
        input_files={f"analysis/{name}": sha256(analysis / name) for name in CROSS_SPECIES_INPUTS})
