"""Full-tree-pipeline cache adapters (read-only; never run solvers or nulls)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from full_tree_pareto import cross_species_analysis as fcs
from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput
from publication.provenance import sha256

CROSS_SPECIES_FAMILY = "full-tree-cross-species"
CROSS_SPECIES_INPUTS = ("provenance.json", "fronts.csv", "metrics.csv", "null_clouds.csv")
CROSS_SPECIES_REFERENCES = dict(
    main=("first_cousin", "second_cousin", "third_cousin"), inset="random_rebuild",
    inset_title="Random rebuild", inset_legend="Random rebuild (insets)", display_draws=1000)


def _run(run) -> Path:
    return Path(run) if isinstance(run, Path) else fcs.run_path(run)


def cross_species(run=fcs.DEFAULT_RUN) -> FrontComparisonInput:
    """Terminal-anchored partial forest (485 cells, 454 edges); Euclidean molecular distance."""
    run = _run(run)
    if not (run / "analysis/provenance.json").is_file():
        raise FileNotFoundError(
            f"No partial-forest comparison cache at {run}; create it with "
            f"`python -m full_tree_pareto.cross_species_analysis --run-id {run.name}`")
    report = fcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    fronts, metrics, clouds = (pd.read_csv(analysis / name) for name in ("fronts.csv", "metrics.csv", "null_clouds.csv"))
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
        input_files={f"analysis/{name}": sha256(analysis / name) for name in CROSS_SPECIES_INPUTS},
        clouds=clouds, references=CROSS_SPECIES_REFERENCES)


POOLED_FAMILY = "full-tree-pooled"
POOLED_DISPLAY_DRAWS = 1000


def endpoint(frame, scope):
    """Endpoint display transform from one layerwise front (travel optimum at the largest weight index)."""
    from full_tree_pareto import pooled_analysis as pa
    from terminal_pareto.front_coordinates import EndpointTransform
    a = frame.loc[frame.weight_index.idxmax()]
    b = frame.loc[frame.weight_index.idxmin()]
    return EndpointTransform.from_endpoints(
        reference_analysis_id=f"{pa.PROFILE}:layerwise:{scope}",
        travel_optimum_assignment_id=f"layerwise:{scope}:{int(a.weight_index)}",
        state_optimum_assignment_id=f"layerwise:{scope}:{int(b.weight_index)}",
        travel_optimum_costs=(a.travel, a.state),
        state_optimum_costs=(b.travel, b.state))


def reference_sets():
    """Figure 7B shows the three cousin shuffles; S4 adds every other reference family."""
    from full_tree_pareto import cousin_references as cr
    from full_tree_pareto import pooled_analysis as pa
    main = (pa.NULLS[0], *cr.REFERENCES)
    return main, (*main, *pa.NULLS[1:])


def caption_meta(data) -> dict:
    return dict(edges=data.edges, leaves=data.leaves, internal=data.internal, round_edges=list(data.round_edges),
                **data.settings)


def pooled_from_results(run, ctx, fronts, layers, nulls):
    """Contract from replay-validated ``pooled_analysis.build`` / ``cousin_references.build`` results."""
    from full_tree_pareto import pooled_analysis as pa
    from full_tree_pareto import pooled_pipeline as settings
    from publication.contracts import FullTreePooledInput

    run = Path(run)
    rounds = [f"round_{i}" for i in range(1, len(ctx.layers) + 1)]
    layerwise = fronts[fronts.method == pa.METHODS[0]]
    analysis = run / "analysis"
    validation = json.loads((analysis / "validation.json").read_text())
    files = sorted(p for p in analysis.rglob("*") if p.is_file())
    return FullTreePooledInput(
        run=run, analysis_id=hashlib.sha256(json.dumps(validation["identity"], sort_keys=True).encode()).hexdigest(),
        edges=len(ctx.edges), leaves=len(ctx.leaves), internal=len(ctx.internal),
        round_edges=tuple(len(layer) for layer in ctx.layers), natural=tuple(map(float, ctx.natural)),
        fronts=fronts, layers=layers, references=dict(nulls), methods=tuple(pa.METHODS),
        main_references=reference_sets()[0], supplement_references=reference_sets()[1],
        inset_reference="Random rebuild",
        round_transforms={scope: endpoint(layerwise[layerwise.scope == scope], scope) for scope in rounds},
        aggregate_transform=endpoint(layerwise[layerwise.scope == "aggregate"], "aggregate"),
        settings=dict(weights=settings.INTERVALS + 1, draws=settings.DRAWS, display_draws=POOLED_DISPLAY_DRAWS),
        validation=dict(assignments_replayed=validation.get("assignments_replayed"), nodes=validation.get("nodes"),
                        edges=validation.get("edges"), rounds=validation.get("rounds")),
        input_files={str(p.relative_to(run)): sha256(p) for p in files})


def pooled(run=None):
    """Pooled full-tree Figure 7/S4 inputs; replays all saved forests and reference draws."""
    from full_tree_pareto import cousin_references as cr
    from full_tree_pareto import pooled_analysis as pa
    from full_tree_pareto import pooled_pipeline as settings

    run = Path(run) if run is not None else Path(pa.DEFAULT_RUN)
    if not (run / "analysis/validation.json").is_file():
        raise FileNotFoundError(f"No pooled full-tree caches at {run}; create them with "
                                "`python -m full_tree_pareto.resume_paired` then `python -m full_tree_pareto.pooled_pipeline`")
    ctx, fronts, layers, nulls = pa.build(run=run, layout_only=True, **settings.sweep_settings())
    nulls = dict(nulls)
    nulls.update(cr.build(ctx, run, draws=settings.DRAWS, seed=settings.SEED, layout_only=True))
    return pooled_from_results(run, ctx, fronts, layers, nulls)
