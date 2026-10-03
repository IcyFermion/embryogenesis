"""Terminal-pipeline cache adapters (read-only; never run solvers or nulls)."""

from __future__ import annotations

import json

import pandas as pd

from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput
from publication.provenance import sha256
from terminal_pareto import cross_species_analysis as tcs

CROSS_SPECIES_FAMILY = "terminal-cross-species"
CROSS_SPECIES_INPUTS = ("provenance.json", "fronts.csv", "metrics.csv", "null_clouds.csv")


def cross_species(run_id: str = tcs.DEFAULT_RUN) -> FrontComparisonInput:
    """Pooled 187-edge terminal comparison; cosine molecular distance."""
    run = tcs.run_path(run_id)
    if not (run / "analysis/provenance.json").is_file():
        raise FileNotFoundError(
            f"No terminal comparison cache at {run}; create it with "
            f"`python terminal_pareto/cross_species_analysis.py --run-id {run_id}`")
    report = tcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    fronts = pd.read_csv(analysis / "fronts.csv")
    metrics = pd.read_csv(analysis / "metrics.csv")
    primary = list(tcs.PRIMARY)
    fronts = fronts[(fronts.cohort == "base") & fronts.config.isin(primary)].reset_index(drop=True)
    metrics = (metrics[(metrics.cohort == "base") & metrics.config.isin(primary)]
               .rename(columns=CANONICAL_FIELDS).reset_index(drop=True))
    return FrontComparisonInput(
        family=CROSS_SPECIES_FAMILY, run=run, analysis_id=report["analysis_id"],
        edges=int(record["cohort_size"]), configs=tuple(tcs.PRIMARY), geometries=tuple(tcs.GEOMETRIES),
        fronts=fronts, metrics=metrics, null_mean="analytic first-cousin",
        caption=dict(edges=int(record["cohort_size"]), cb_reference_size=int(record["cb_reference_size"])),
        validation=report,
        input_files={f"analysis/{name}": sha256(analysis / name) for name in CROSS_SPECIES_INPUTS})
