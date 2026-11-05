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
    """Terminal-anchored partial forest (485 cells, 454 edges); cosine molecular distance."""
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


LAYERS_FAMILY = "full-tree-layers"
LAYER_RUNS = dict(hourglass="hourglass_20261019", stages="stages_20261020")
FIRST_LAYER = 3  # author decision 2026-11-04: layers 3-8 are drawn and claimed; layers 1-2 are left out
FLANKS = {"L8": (8,), "L4-5": (4, 5)}
MID = (6, 7)


def _checked_layer_run(module, run_id):
    """File and recorded source hashes of one Track E run; the tables are not recomputed (the bootstrap draws)."""
    from full_tree_pareto import layer_features as lf
    analysis = lf.run_path(run_id) / "analysis"
    if not (analysis / "provenance.json").is_file():
        raise FileNotFoundError(f"No Track E cache at {analysis.parent}; create it with "
                                f"`python -m {module.__name__} --run-id {run_id}`")
    record = json.loads((analysis / "provenance.json").read_text())
    if record["version"] != module.VERSION:
        raise ValueError(f"Unsupported {module.__name__} cache: {record['version']}")
    for name, expected in record["files"].items():
        if lf.file_hash(analysis / name) != expected:
            raise ValueError(f"Cache hash mismatch in {run_id}: {name}")
    for name, expected in record["sources"].items():
        if lf.file_hash(module.ROOT / name) != expected:
            raise ValueError(f"Input/code changed since {run_id}: {name}; create a new numerical run")
    return analysis, record


def _per_edge(references, stages, tree, config):
    rows = references[(references.config == config) & ~references.degenerate.astype(bool)].copy()
    rows["layer"] = 9 - rows["round"]
    rows = rows[rows.layer >= FIRST_LAYER].merge(stages[stages.tree == tree][["layer", "alive_median"]], on="layer")
    return pd.DataFrame(dict(tree=tree, config=config, layer=rows.layer, edges=rows.edges, born=rows.alive_median,
                             natural=rows.natural_travel / rows.edges,
                             aunt=rows.first_cousin_travel_mean / rows.edges,
                             random=rows.full_random_travel_mean / rows.edges)).sort_values("layer")


def _changes(per_edge):
    """Percent change in per-edge travel from each flank to each mid layer (flank layers averaged)."""
    rows = []
    for tree, costs in per_edge.groupby("tree", sort=False):
        costs = costs.set_index("layer")
        for flank, flank_layers in FLANKS.items():
            base = costs.loc[list(flank_layers), ["natural", "aunt"]].mean()
            for mid in MID:
                rows.append(dict(tree=tree, contrast=f"L{mid} vs {flank}",
                                 natural=100 * (costs.loc[mid, "natural"] / base["natural"] - 1),
                                 aunt=100 * (costs.loc[mid, "aunt"] / base["aunt"] - 1)))
    return pd.DataFrame(rows)


def layers(runs=None):
    """Track E inputs (Figure 6B-C of the revised layout, S5): hash-checked runs, reproduced against the cached
    round references."""
    from full_tree_pareto import layer_features as lf
    from full_tree_pareto import layer_hourglass as lh
    from full_tree_pareto import layer_stages as ls
    from publication.contracts import FullTreeLayersInput

    runs = dict(LAYER_RUNS, **(runs or {}))
    hourglass, h_record = _checked_layer_run(lh, runs["hourglass"])
    stages, s_record = _checked_layer_run(ls, runs["stages"])
    saved = pd.read_csv(hourglass / "layers.csv")
    trees, upstream = lh.load()
    if upstream != h_record["upstream"]:
        raise ValueError("Round references changed since the hourglass run")
    reproduction = lh.check_reproduction(saved, trees)

    stage_table = pd.read_csv(stages / "layers.csv")
    pooled_refs = trees["pooled"]["references"]
    partial_refs = trees["partial"]["references"]
    partial_refs = partial_refs[partial_refs.geometry == lh.GEOMETRY]

    rows = saved[(saved.series == "pooled_ce_protein") & (saved.family == "first_cousin")]
    fronts = pooled_refs[pooled_refs.config == "ce_protein"].assign(layer=lambda f: 9 - f["round"])
    table = (rows.merge(stage_table[stage_table.tree == "pooled"][["layer", "alive_median", "alive_q25", "alive_q75"]],
                        on="layer")
             .merge(fronts[["layer", "d_LP", "distinct_front_points", "max_edge_retention", "degenerate"]],
                    on="layer"))
    table = table[table.layer >= FIRST_LAYER].sort_values("layer")
    if table.degenerate.astype(bool).any():
        raise ValueError("A drawn layer has a degenerate front")
    shown = pd.DataFrame(dict(
        layer=table.layer, round=table["round"], edges=table.edges, groups=table.units,
        born=table.alive_median, born_q25=table.alive_q25, born_q75=table.alive_q75,
        front_distance=table.d_LP, front_points=table.distinct_front_points, fewest_rewired=1 - table.max_edge_retention,
        travel=1 - table.travel_advantage, travel_low=1 - table.travel_high, travel_high=1 - table.travel_low,
        state=1 - table.state_advantage, state_low=1 - table.state_high, state_high=1 - table.state_low)).reset_index(drop=True)
    per_edge = pd.concat([_per_edge(pooled_refs, stage_table, "pooled", "ce_protein"),
                          _per_edge(partial_refs, stage_table, "partial", "cb_rna")], ignore_index=True)

    run = lf.RUN_ROOT
    analysis_ids = dict(hourglass=h_record["analysis_id"], stages=s_record["analysis_id"])
    left_out = saved[(saved.series == "pooled_ce_protein") & (saved.family == "first_cousin")
                     & (saved.layer < FIRST_LAYER)].sort_values("layer")
    settings = h_record["settings"]
    return FullTreeLayersInput(
        run=run, analysis_id=hashlib.sha256(json.dumps(analysis_ids, sort_keys=True).encode()).hexdigest(),
        layers=shown, per_edge=per_edge, changes=_changes(per_edge),
        caption=dict(first_layer=FIRST_LAYER, cells=int(shown.edges.sum()), leaves=int(shown.edges.iloc[-1]),
                     left_out=[int(e) for e in left_out.edges], draws=int(settings["draws"]),
                     level=float(settings["level"]), first_groups=int(shown.groups.iloc[0]),
                     first_points=int(shown.front_points.iloc[0]), first_distance=float(shown.front_distance.iloc[0]),
                     random=[float(per_edge[per_edge.tree == "pooled"].random.min() * 1e3),
                             float(per_edge[per_edge.tree == "pooled"].random.max() * 1e3)]),
        validation=dict(analysis_ids=analysis_ids, reproduction=reproduction,
                        verified_source_hashes=len(h_record["sources"]) + len(s_record["sources"])),
        input_files={str((path).relative_to(run)): sha256(path)
                     for folder in (hourglass, stages) for path in sorted(folder.iterdir()) if path.is_file()})
