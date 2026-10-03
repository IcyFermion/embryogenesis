"""Terminal-pipeline cache adapters (read-only; never run solvers or nulls)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput
from publication.provenance import sha256
from terminal_pareto import cross_species_analysis as tcs
from terminal_pareto.front_coordinates import EndpointTransform

CROSS_SPECIES_FAMILY = "terminal-cross-species"
CROSS_SPECIES_INPUTS = ("provenance.json", "fronts.csv", "metrics.csv", "null_clouds.csv",
                        "assignments_and_costs.npz")
CROSS_SPECIES_REFERENCES = dict(
    main=("first_cousin", "second_cousin", "third_cousin"), inset="full_random",
    inset_title="Full-random shuffle", inset_legend="Full-random shuffle (insets)", display_draws=None)


def _run(run) -> Path:
    return Path(run) if isinstance(run, Path) else tcs.run_path(run)


def transform(metadata: dict) -> EndpointTransform:
    return EndpointTransform(**{k: v for k, v in metadata.items()
                                if k not in ("display_mode", "clipping", "axis_limits")})


def tracking_transfers(record: dict, analysis: Path) -> dict:
    """Individual-embryo optima re-scored under pooled 3D costs and pooled anchors (replay only)."""
    children = np.arange(record["cohort_size"])
    transfers = {}
    with np.load(analysis / "assignments_and_costs.npz", allow_pickle=False) as arrays:
        for config in tcs.PRIMARY:
            key = f"raw3d__{config}"
            endpoint = transform(record["endpoint_transforms"][key])
            species = "cb" if config == "cb_rna" else "ce"
            curves = []
            for j, replica in enumerate(record["tracking"][species]):
                perms = arrays[f"{key}__{replica['label']}__permutations"]
                t = arrays[f"{key}__travel"][children, perms].sum(axis=1)
                s = arrays[f"{key}__state"][children, perms].sum(axis=1)
                x, y = endpoint.transform(t, s)
                label = f"AF16 {'p1' if j == 0 else 'p5'}" if species == "cb" else f"Embryo {j + 1}"
                curves.append((label, np.asarray(x), np.asarray(y)))
            transfers[config] = curves
    return transfers


def cross_species(run=tcs.DEFAULT_RUN) -> FrontComparisonInput:
    """Pooled 187-edge terminal comparison; cosine molecular distance."""
    run = _run(run)
    if not (run / "analysis/provenance.json").is_file():
        raise FileNotFoundError(
            f"No terminal comparison cache at {run}; create it with "
            f"`python terminal_pareto/cross_species_analysis.py --run-id {run.name}`")
    report = tcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    primary = list(tcs.PRIMARY)
    fronts, metrics, clouds = (pd.read_csv(analysis / name) for name in ("fronts.csv", "metrics.csv", "null_clouds.csv"))
    fronts = fronts[(fronts.cohort == "base") & fronts.config.isin(primary)].reset_index(drop=True)
    metrics = (metrics[(metrics.cohort == "base") & metrics.config.isin(primary)]
               .rename(columns=CANONICAL_FIELDS).reset_index(drop=True))
    clouds = clouds[clouds.config.isin(primary)].reset_index(drop=True)
    return FrontComparisonInput(
        family=CROSS_SPECIES_FAMILY, run=run, analysis_id=report["analysis_id"],
        edges=int(record["cohort_size"]), configs=tuple(tcs.PRIMARY), geometries=tuple(tcs.GEOMETRIES),
        fronts=fronts, metrics=metrics, null_mean="analytic first-cousin",
        caption=dict(edges=int(record["cohort_size"]), cb_reference_size=int(record["cb_reference_size"])),
        validation=report,
        input_files={f"analysis/{name}": sha256(analysis / name) for name in CROSS_SPECIES_INPUTS},
        clouds=clouds, references=CROSS_SPECIES_REFERENCES,
        extras=dict(tracking_transfers=tracking_transfers(record, analysis)))
