"""Cross-species terminal sketches by subtree (Figures 4-5 analogues) and by cell type (Figure 6 analogues)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import numpy as np

from publication import notation as nt
from publication import style
from publication.artists import AXIS_LABELS
from publication.sketches.panels import CONFIG_MARKERS, METRIC_COLUMNS, metric_rows

CONFIGS = ("ce_protein", "ce_rna", "cb_rna")
FIGURE5_BOLD = {"P0", "AB", "ABa", "ABp", "P1"}
BRANCHES = {"P0": "Root & AB", "ABa": "ABa branch", "ABp": "ABp branch", "P1": "P1 branch"}
TYPE_ORDER = ("neuron", "muscle", "epithelium", "alimentary", "other")
DISTANCE_CMAP = LinearSegmentedColormap.from_list("d_lp", ["#F7EEF4", style.SEMANTIC_COLORS["lineage_front_distance"], "#5E2350"])
FATE_CAVEAT = ("Cell types: C. elegans annotation joined on lineage names and assumed for the matched C. briggsae cells.")


def _degenerate_note(metrics, kind, geometry):
    rows = metrics[(metrics.kind == kind) & (metrics.geometry == geometry) & metrics.degenerate.astype(bool)]
    if not len(rows):
        return ""
    return "Blank: " + "; ".join(f"{name} ({', '.join(nt.CONFIG_LABELS[c] for c in group.config)})"
                                for name, group in rows.groupby("name")) + " (single optimal assignment)."


def subtree_metrics(data, geometry):
    from publication.sketches.data import lineage_preorder
    metrics = data["metrics"]
    table = metrics[metrics.kind == "subtree"]
    names = lineage_preorder(table.name.unique())
    sizes = table.groupby("name").n.first()
    rows = [(name, f"{name} · {sizes[name]}") for name in names]
    headings = {names.index(k): v for k, v in BRANCHES.items() if k in names}
    note = (f"Each subtree is its own assignment problem. As in Figure 5B, {nt.math('null_front_distance')} is "
            f"omitted: small subtrees have tiny endpoint spans that place {nt.math('null_mean')} far outside the unit "
            "square. The same caveat applies to large d_LP where a front is nearly a point (e.g. C. briggsae P2: "
            "P* differs from the natural lineage by <2% raw travel and <0.1% cell-state cost). "
            + _degenerate_note(metrics, "subtree", geometry))
    columns = tuple(c for c in METRIC_COLUMNS if c[0] != "null_front_distance")
    return metric_rows(table, row_key="name", rows=rows, configs=CONFIGS, geometries=(geometry,), bold=FIGURE5_BOLD,
                       headings=headings, row_height=.26, columns=columns,
                       title=f"Matched terminal edges (187), {nt.GEOMETRY_LABELS[geometry]} tracking: canonical metrics by subtree",
                       note=note)


def subtree_map(data, geometry):
    """Figure 4 layout; each subtree shows three dots (CE protein, CE RNA, CB RNA) shaded by d_LP."""
    from publication.figures import terminal_subtree_map as tsm
    from terminal_pareto import data_loader as dl
    from terminal_pareto import lineage_metrics as lm
    from terminal_pareto.cross_species_analysis import ROOT
    from utils import load_json
    metrics = data["metrics"]
    table = metrics[(metrics.kind == "subtree") & (metrics.geometry == geometry)]
    lineage = load_json(ROOT / "data/cell_lineage.json")
    manifest = data["manifest"]
    type_map = dict(zip(manifest.terminal, manifest.cell_type))
    nodes, edges, *_ = tsm.collect_nodes(lineage, lm.build_lineage_tree_index(lineage), set(table.name),
                                         tsm.terminal_inorder(lineage, dl.map_names), type_map,
                                         set(manifest.terminal), dl.map_names)
    radii = {name: .3 for name in nodes}
    tsm.tidy_display_tree(nodes, edges, radii)
    tsm.stagger_crowded_levels(nodes, radii)
    finite = table.natural_front_distance.dropna()
    norm = Normalize(0, max(.05, float(finite.max())))
    ys = [nd["y"] for nd in nodes.values()]
    fig, ax = plt.subplots(figsize=(15, 2.2 + .72 * (max(ys) - min(ys) + .9)))
    for p, c in edges:
        ax.plot([nodes[p]["x"], nodes[c]["x"]], [nodes[p]["y"], nodes[c]["y"]], color="#B5B5B5", lw=.6, zorder=1)
    dot = .12
    for name, nd in nodes.items():
        for i, config in enumerate(CONFIGS):
            x = nd["x"] + (i - 1) * 2.3 * dot
            row = table[(table.name == name) & (table.config == config)]
            value = row.natural_front_distance.iloc[0] if len(row) else np.nan
            if np.isnan(value):
                ax.add_patch(Circle((x, nd["y"]), dot, facecolor="white", edgecolor="#9A9A9A", lw=.7, ls=":", zorder=3))
                continue
            on_front = value < 1e-9
            ax.add_patch(Circle((x, nd["y"]), dot, facecolor=DISTANCE_CMAP(norm(value)), zorder=3,
                                edgecolor=style.SPECIES_COLORS[config], lw=1.6 if on_front else .9,
                                ls="-" if on_front else (0, (1.5, 1.2))))
        ax.text(nd["x"], nd["y"] - dot - .06, f"{name} {nd['n']}", ha="center", va="top", fontsize=7.5, color="#333333")
    xs = [nd["x"] for nd in nodes.values()]
    ax.set_xlim(min(xs) - .6, max(xs) + .6)
    ax.set_ylim(min(ys) - .5, max(ys) + .35)
    ax.set_aspect("equal")
    ax.axis("off")
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=DISTANCE_CMAP), ax=ax, fraction=.025, pad=.01)
    bar.set_label(f"{nt.label('natural_front_distance')}, {nt.math('natural_front_distance')}")
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=style.SPECIES_COLORS[c],
                      markeredgewidth=1.4, markersize=7, label=nt.CONFIG_LABELS[c]) for c in CONFIGS]
    handles += [Line2D([], [], marker="o", ls="", markerfacecolor="#DDDDDD", markeredgecolor="#555555", markersize=7,
                       label="solid ring: natural exactly on front"),
                Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="#9A9A9A", markersize=7,
                       label="dotted grey: single optimal assignment")]
    ax.legend(handles=handles, loc="upper left", ncol=5, fontsize=7.5, bbox_to_anchor=(0, 0.0), frameon=False)
    ax.set_title(f"Subtrees with ≥12 matched terminal edges, {nt.GEOMETRY_LABELS[geometry]} tracking: "
                 "left to right C. elegans protein, C. elegans RNA, C. briggsae RNA", fontsize=10)
    return fig


def type_retention(data):
    """Figure 6A analogue: natural-edge retention by cell type along each whole-cohort front."""
    retention, metrics = data["retention"], data["global_metrics"]
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 6.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.12, right=.92, top=.86, bottom=.12, hspace=.32, wspace=.06)
    norm = Normalize(0, 1)
    for g, geometry in enumerate(("raw3d", "xy")):
        for c, config in enumerate(CONFIGS):
            ax = axes[g, c]
            frame = retention[(retention.geometry == geometry) & (retention.config == config)]
            order = [t for t in TYPE_ORDER if t in set(frame.cell_type)]
            first = frame[frame.cell_type == order[0]].sort_values("u")
            u = first.u.to_numpy()
            edges = np.r_[0., (u[1:] + u[:-1]) / 2, 1.]
            values = np.array([frame[frame.cell_type == t].sort_values("u").retention.to_numpy() for t in order])
            mesh = ax.pcolormesh(edges, np.arange(len(order) + 1), values, cmap=style.EDGE_RETENTION_CMAP, norm=norm,
                                 shading="flat", rasterized=True)
            metric = metrics[(metrics.geometry == geometry) & (metrics.config == config)].iloc[0]
            ax.axvline(metric.canonical_position, color="white", lw=1.1, ls=(0, (3, 2)))
            ax.set_ylim(len(order), 0)
            counts = frame.groupby("cell_type").n.first()
            ax.set_yticks(np.arange(len(order)) + .5, [f"{t} (n={counts[t]})" for t in order])
            ax.set_title(f"{nt.CONFIG_LABELS[config]} · {nt.GEOMETRY_LABELS[geometry]}", fontsize=9.5)
            if g == 1:
                ax.set_xlabel(f"Position along the front, {nt.math('canonical_position')}")
    bar = fig.colorbar(mesh, cax=fig.add_axes([.935, .12, .012, .74]))
    bar.set_label("Natural edges retained")
    fig.suptitle("Natural-edge retention by cell type along the whole-cohort fronts (187 matched edges)", fontsize=11)
    fig.text(.12, .02, f"Dashed line: {nt.math('closest_point')} of the whole cohort. 'other' merges types with "
             f"<15 cells (glial 6, programmed death 1, ...). {FATE_CAVEAT}", fontsize=7.8, color="#444444")
    return fig


def within_type_fronts(data, geometry):
    """Figure 6B analogue: within-type fronts (assignments never cross types), configurations overlaid."""
    fronts, metrics = data["fronts"], data["metrics"]
    types = [t for t in TYPE_ORDER if t != "other"]
    fig, axes = plt.subplots(1, len(types), figsize=(14.5, 4.1), sharex=True, sharey=True)
    fig.subplots_adjust(left=.07, right=.98, top=.8, bottom=.27, wspace=.1)
    for ax, cell_type in zip(axes, types):
        ax.axhline(0, color="#777777", lw=.6, ls=":")
        ax.axvline(0, color="#777777", lw=.6, ls=":")
        for config, ls in zip(CONFIGS, ("-", "--", "-.")):
            color = style.SPECIES_COLORS[config]
            curve = fronts[(fronts.kind == "cell_type") & (fronts.name == cell_type) & (fronts.geometry == geometry)
                           & (fronts.config == config)].sort_values("sweep_index")
            metric = metrics[(metrics.kind == "cell_type") & (metrics.name == cell_type) & (metrics.geometry == geometry)
                             & (metrics.config == config)].iloc[0]
            if not len(curve):
                continue
            ax.plot(curve.D1, curve.D2, color=color, lw=1.3, ls=ls, alpha=.85)
            ax.scatter(metric.natural_D1, metric.natural_D2, marker=CONFIG_MARKERS[config], s=34, color=color,
                       edgecolor="#222222", lw=.6, zorder=5)
        n = int(metrics[(metrics.kind == "cell_type") & (metrics.name == cell_type)].n.iloc[0])
        ax.set_title(f"{cell_type} (n={n})", fontsize=9.5)
        ax.set_xlabel(AXIS_LABELS[0].replace(" (", "\n("), fontsize=8)
        ax.grid(True, alpha=.3)
    axes[0].set_ylabel(AXIS_LABELS[1], fontsize=8)
    handles = [Line2D([], [], color=style.SPECIES_COLORS[c], ls=ls, marker=CONFIG_MARKERS[c], markeredgecolor="#222222",
                      label=f"{nt.CONFIG_LABELS[c]}: front and natural lineage") for c, ls in zip(CONFIGS, ("-", "--", "-."))]
    fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=8, bbox_to_anchor=(.5, .97))
    fig.suptitle(f"Within-type fronts, {nt.GEOMETRY_LABELS[geometry]} tracking (each in its own endpoint span)",
                 fontsize=11, y=.995)
    fig.text(.07, .02, "Assignments never cross cell types. " + FATE_CAVEAT, fontsize=7.8, color="#444444")
    return fig


def type_metrics(data):
    metrics = data["metrics"]
    table = metrics[metrics.kind == "cell_type"]
    sizes = table.groupby("name").n.first()
    rows = [(t, f"{t} · {sizes[t]}") for t in TYPE_ORDER if t in sizes]
    fallback = sorted(table[table.null_model == "full_random_fallback"].name.unique())
    note = (f"Within-type exact first-cousin null mean {nt.math('null_mean')}"
            + (f"; {', '.join(fallback)}: first-cousin null has zero variance, so the sweep uses full-random "
               f"within-type scaling and {nt.math('null_front_distance')} is not shown" if fallback else "")
            + ". " + FATE_CAVEAT)
    return metric_rows(table, row_key="name", rows=rows, configs=CONFIGS, geometries=("raw3d", "xy"),
                       title="Within-type canonical metrics (187 matched terminal edges)", note=note)
