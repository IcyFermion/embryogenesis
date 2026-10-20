"""Read-only sketch inputs: per-round fronts, clouds and metrics as comparison contracts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from publication.contracts import CANONICAL_FIELDS, FrontComparisonInput

ROUND_REFERENCES = dict(main=("first_cousin", "second_cousin", "third_cousin"), inset="full_random",
                        inset_title="Full-random shuffle", inset_legend="Full-random shuffle (insets)",
                        display_draws=1000)


def _transform(metadata):
    from terminal_pareto.front_coordinates import EndpointTransform
    return EndpointTransform(**{k: v for k, v in metadata.items() if k not in ("display_mode", "clipping", "axis_limits")})


def _rounds(run, target, *, family, labels):
    """One ``FrontComparisonInput`` per nondegenerate round; ``labels`` maps cache label -> (geometry, config)."""
    from full_tree_pareto import round_references as rr
    loaded = rr.load(run, target)
    fronts, metrics = loaded["fronts"].copy(), loaded["metrics"].copy()
    if "geometry" not in fronts:
        fronts["geometry"], metrics["geometry"] = "raw3d", "raw3d"
    metrics = metrics.rename(columns=CANONICAL_FIELDS)
    geometries = tuple(dict.fromkeys(g for g, _ in labels.values()))
    configs = tuple(dict.fromkeys(c for _, c in labels.values()))
    result, skipped = {}, {}
    for r in sorted(metrics["round"].unique()):
        rows = metrics[metrics["round"] == r]
        if rows.degenerate.astype(bool).any():
            skipped[int(r)] = sorted(f"{g}/{c}" for g, c in rows.loc[rows.degenerate.astype(bool), ["geometry", "config"]].values)
            continue
        clouds = []
        for label, (geometry, config) in labels.items():
            transform = _transform(loaded["transforms"][f"{label}__round_{r}"])
            for fam in rr.FAMILIES:
                costs = loaded["costs"][f"{label}__round_{r}__{fam}"]
                x, y = transform.transform(costs[:, 0], costs[:, 1])
                clouds.append(pd.DataFrame(dict(geometry=geometry, config=config, family=fam, D1=x, D2=y)))
        result[int(r)] = FrontComparisonInput(
            family=family, run=Path(run), analysis_id=loaded["identity_sha"], edges=int(rows.edges.iloc[0]),
            configs=configs, geometries=geometries,
            fronts=fronts[fronts["round"] == r].reset_index(drop=True),
            metrics=rows.reset_index(drop=True), null_mean="exact per-round first-cousin",
            caption=dict(round=int(r), edges=int(rows.edges.iloc[0]), draws=loaded["identity"]["draws"]),
            clouds=pd.concat(clouds, ignore_index=True), references=ROUND_REFERENCES)
    return result, skipped, metrics


def pooled_rounds():
    """CE protein pooled full tree: eight rounds (491 ... 5 edges)."""
    from full_tree_pareto import pooled_analysis as pa
    return _rounds(pa.DEFAULT_RUN, "pooled", family="sketch-pooled-rounds", labels={"ce_protein": ("raw3d", "ce_protein")})


def cross_species_rounds():
    """Partial forest: six rounds (187 ... 2 edges), three configurations x two geometries."""
    from full_tree_pareto import cross_species_analysis as fcs
    labels = {f"{g}__{c}": (g, c) for g in fcs.GEOMETRIES for c in fcs.PRIMARY}
    return _rounds(fcs.run_path(fcs.DEFAULT_RUN), "cross_species", family="sketch-cross-species-rounds", labels=labels)


def terminal_partitions():
    """187-edge cross-species terminal comparison split by subtree and by cell type (read-only)."""
    from terminal_pareto import cross_species_analysis as tcs
    from terminal_pareto import cross_species_partitions as tcp
    loaded = tcp.load()  # Verifies its own identity, sources and files against the validated run.
    analysis = tcs.run_path(tcs.DEFAULT_RUN) / "analysis"
    global_metrics = pd.read_csv(analysis / "metrics.csv")
    global_metrics = global_metrics[global_metrics.cohort == "base"].rename(columns=CANONICAL_FIELDS)
    return dict(metrics=loaded["metrics"].rename(columns=CANONICAL_FIELDS), fronts=loaded["fronts"],
                retention=loaded["type_retention"], manifest=pd.read_csv(analysis / "cell_manifest.csv"),
                global_metrics=global_metrics, record=loaded["record"], identity=loaded["identity"])


def lineage_preorder(names):
    """``names`` in lineage preorder (root first, anterior daughter first)."""
    from terminal_pareto.cross_species_analysis import ROOT
    from utils import load_json, lineage_name_mapping
    wanted, order = set(names), []

    def walk(node):
        name = lineage_name_mapping(node.get("did", node.get("name", "")))
        if name in wanted:
            order.append(name)
        for child in node.get("children", []):
            walk(child)
    walk(load_json(ROOT / "data/cell_lineage.json"))
    return order


def _verified_analysis(out):
    """Read-only: check a layer_features_v1 run's recorded file hashes; returns its analysis folder."""
    import json
    from full_tree_pareto import layer_features as lf
    analysis = Path(out) / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    for name, expected in record["files"].items():
        if lf.file_hash(analysis / name) != expected:
            raise ValueError(f"Changed file in {out}: {name}")
    return analysis, record


def layer_hourglass(hourglass_run="hourglass_20261019", stages_run="stages_20261020"):
    """Track E3/E4 inputs: per-layer advantages and dips, birth stages, Large et al. JSD, per-round d_LP."""
    from full_tree_pareto import cross_species_analysis as fcs
    from full_tree_pareto import layer_features as lf
    from full_tree_pareto import pooled_analysis as pa
    from full_tree_pareto import round_references as rr
    hourglass, h_record = _verified_analysis(lf.run_path(hourglass_run))
    stages, s_record = _verified_analysis(lf.run_path(stages_run))
    pooled = rr.load(pa.DEFAULT_RUN, "pooled")["metrics"].assign(tree="pooled")
    cross = rr.load(fcs.run_path(fcs.DEFAULT_RUN), "cross_species")["metrics"]
    cross = cross[cross.geometry == "raw3d"].assign(tree="partial")
    fronts = pd.concat([pooled, cross], ignore_index=True)[["tree", "config", "round", "edges", "d_LP", "degenerate"]]
    fronts["layer"] = 9 - fronts["round"]
    return dict(layers=pd.read_csv(hourglass / "layers.csv"), contrasts=pd.read_csv(hourglass / "contrasts.csv"),
                stages=pd.read_csv(stages / "layers.csv"), jsd=pd.read_csv(stages / "jsd_stages.csv"), fronts=fronts,
                ids=dict(hourglass=h_record["analysis_id"], stages=s_record["analysis_id"]))
