"""Build exploratory sketches into a new folder: ``python -m publication.sketches --output-dir DIR [--only KEY ...]``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

from publication import style
from publication.figures import cross_species as cs_figures
from publication.provenance import sha256


def _save(fig, out, stem):
    style.save_figure(fig, Path(out) / f"{stem}.png")
    plt.close(fig)
    return f"{stem}.png"


def full_tree_rounds(out):
    """Idea 1: Figure 7 as eight Figure-2-style round panels (aggregate stays in the supplement)."""
    from publication.sketches import data, rounds
    pooled, skipped, metrics = data.pooled_rounds()
    entries = [dict(key="fig7_round_grid", group="Figure 7 by round", file=_save(rounds.round_grid(
        pooled, title="C. elegans full tree: layerwise fronts by contraction round (protein, pooled travel)"), out,
        "fig7_round_grid"), title="All eight rounds, Figure 2 style",
        note="Each round: its own endpoint normalization, per-round cousin clouds (1,000 shown), exact first-cousin "
             "mean C (+), full-random inset, natural lineage (X) joined to the maximum-retention assignment.")]
    for r, d in pooled.items():
        stem = f"fig7_round_{r}"
        rounds.figure2_single(d, out, stem)
        entries.append(dict(key=stem, group="Figure 7 by round", file=f"{stem}.png",
                            title=rounds.round_title(r, d.edges, max(pooled)),
                            note="Exact Figure 2 layout applied to this round."))
    entries.append(dict(key="fig7_round_metrics", group="Figure 7 by round", file=_save(rounds.round_metrics(
        metrics, configs=("ce_protein",), geometries=("raw3d",),
        title="C. elegans full tree: canonical metrics by round"), out, "fig7_round_metrics"),
        title="Canonical metrics by round", note="u, d_LP, d_CP and maximum edge retention per round."))
    return entries


def cross_species_rounds(out):
    """Idea 2c: partial-forest cross-species comparison by round instead of the aggregate (Figures 10-11)."""
    from publication.adapters import full_tree
    from publication.sketches import data, rounds
    per_round, skipped, metrics = data.cross_species_rounds()
    aggregate = full_tree.cross_species()
    note = "; ".join(f"Round {r} omitted where degenerate ({', '.join(v)})" for r, v in skipped.items())
    entries = [dict(key="xs_round_metrics", group="Cross-species partial forest by round", file=_save(
        rounds.round_metrics(metrics[~metrics.degenerate.astype(bool)], configs=aggregate.configs,
                             geometries=aggregate.geometries, aggregate=aggregate.metrics,
                             title="Partial forest (485 cells, 454 edges): canonical metrics by round and configuration",
                             degenerate_note=note + "." if note else None), out, "xs_round_metrics"),
        title="Canonical metrics by round, three configurations",
        note="Replaces the aggregate Figure 11C rows with one row per round; the aggregate row is kept for reference. "
             + (note if note else ""))]
    for r, d in per_round.items():
        stem = f"xs_round_{r}_comparison"
        fig, _ = cs_figures.compose_comparison(d, title=f"Partial forest, round {r} ({{edges}} edges): "
                                               "C. elegans protein / RNA and C. briggsae RNA")
        entries.append(dict(key=stem, group="Cross-species partial forest by round", file=_save(fig, out, stem),
                            title=f"Round {r} · {d.edges} edges, Figure 8/10 layout",
                            note="Per-round cousin clouds and full-random insets."))
    return entries


def cross_species_partitions(out):
    """Idea 2a-b: 187-edge cross-species terminal comparison by subtree and by cell type."""
    from publication.sketches import data, partitions
    d = data.terminal_partitions()
    group_s, group_t = "Cross-species terminal by subtree", "Cross-species terminal by cell type"
    entries = []
    for geometry in ("raw3d", "xy"):
        entries.append(dict(key=f"xs_subtree_map_{geometry}", group=group_s,
                            file=_save(partitions.subtree_map(d, geometry), out, f"xs_subtree_map_{geometry}"),
                            title=f"Subtree map, {geometry}",
                            note="Figure 4 layout; per subtree three dots (CE protein, CE RNA, CB RNA) shaded by d_LP; "
                                 "solid ring = natural exactly on the front."))
        entries.append(dict(key=f"xs_subtree_metrics_{geometry}", group=group_s,
                            file=_save(partitions.subtree_metrics(d, geometry), out, f"xs_subtree_metrics_{geometry}"),
                            title=f"Canonical metrics by subtree, {geometry}",
                            note="Figure 5B-style rows for the 25 subtrees with at least 12 matched terminal edges."))
    entries.append(dict(key="xs_type_retention", group=group_t,
                        file=_save(partitions.type_retention(d), out, "xs_type_retention"),
                        title="Retention by cell type along the whole-cohort fronts",
                        note="Figure 6A analogue for all three configurations and both geometries; no new solves."))
    for geometry in ("raw3d", "xy"):
        entries.append(dict(key=f"xs_within_type_fronts_{geometry}", group=group_t,
                            file=_save(partitions.within_type_fronts(d, geometry), out, f"xs_within_type_fronts_{geometry}"),
                            title=f"Within-type fronts, {geometry}",
                            note="Figure 6B analogue for types with at least 15 cells; configurations overlaid."))
    entries.append(dict(key="xs_type_metrics", group=group_t,
                        file=_save(partitions.type_metrics(d), out, "xs_type_metrics"),
                        title="Within-type canonical metrics", note="Rows are cell types; both geometries."))
    return entries


SKETCHES = {"full-tree-rounds": full_tree_rounds, "cross-species-rounds": cross_species_rounds,
            "cross-species-partitions": cross_species_partitions}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--only", nargs="*", choices=sorted(SKETCHES))
    args = parser.parse_args(argv)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    style.configure()
    index_path = out / "index.json"
    index = json.loads(index_path.read_text()) if index_path.exists() else {}
    for name in args.only or SKETCHES:
        entries = SKETCHES[name](out)
        for entry in entries:
            entry["sha256"] = sha256(out / entry["file"])
        index[name] = entries
        print(f"{name}: {len(entries)} sketches", flush=True)
    index_path.write_text(json.dumps(index, indent=2) + "\n")


if __name__ == "__main__":
    main()
