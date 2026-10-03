"""Figure 7B references in the terminal Figures 8/9 comparison layouts."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from full_tree_pareto.cross_species_analysis import PRIMARY, GEOMETRIES, FAMILIES, ROOT, digest
from terminal_pareto import plot_style as ps
from terminal_pareto import fig_terminal_cross_species as terminal_style
from terminal_pareto.fig2_fig3_ce_terminal_pareto import EDGE_RETENTION_CMAP, proportional_limits

LABELS = terminal_style.LABELS
DISPLAY_DRAWS = 1000
COMPARISON = "partial_forest_cross_species_comparison"
OVERLAY = "partial_forest_cross_species_overlay"


def archive_previous(run, directory):
    source = Path(run) / directory
    if not source.exists() or not any(source.iterdir()):
        return None
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = Path(run) / "layout_history" / stamp / directory
    hashes = {str(p.relative_to(source)): digest(p) for p in source.rglob("*") if p.is_file()}
    shutil.copytree(source, target)
    for name, expected in hashes.items():
        if digest(target / name) != expected:
            raise ValueError(f"Failed layout archive: {name}")
    (target.parent / "manifest.json").write_text(json.dumps(dict(directory=directory, files=hashes), indent=2) + "\n")
    return target.parent


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)


def comparison(fronts, metrics, clouds, out, record):
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 9.8), sharex=True, sharey=True)
    fig.subplots_adjust(left=.06, right=.91, bottom=.16, top=.87, hspace=.32, wspace=.13)
    fig.suptitle(f"Partial-forest layerwise comparison | {record['edges']} shared edges", y=.97, fontsize=12)
    visible_x, visible_y, inset_limits = [], [], {}
    for row, geometry in enumerate(GEOMETRIES):
        for col, config in enumerate(PRIMARY):
            ax = axes[row, col]
            terminal_style.base_axes(ax)
            curve, metric = terminal_style.select(fronts, metrics, geometry, config)
            null = clouds[(clouds.geometry == geometry) & (clouds.config == config)]
            for family in FAMILIES[:3]:
                full = null[null.family == family]
                sample = full.iloc[np.linspace(0, len(full)-1, min(DISPLAY_DRAWS, len(full)), dtype=int)]
                color = ps.NULL_MODEL_COLORS[family]
                ax.scatter(sample.D1, sample.D2, s=7, color=color, alpha=.19, edgecolors="none", rasterized=True, zorder=1)
                mean = (metric.null_D1, metric.null_D2) if family == "first_cousin" else (full.D1.mean(), full.D2.mean())
                ax.scatter(*mean, marker="+", s=34, color=color, lw=1.1, zorder=4)
                visible_x.extend(sample.D1)
                visible_y.extend(sample.D2)
            ax.plot(curve.D1, curve.D2, color=ps.COLORS["blue"], lw=1.3, alpha=.65, zorder=2)
            front = ax.scatter(curve.D1, curve.D2, c=curve.edge_retention, cmap=EDGE_RETENTION_CMAP,
                               norm=Normalize(0, 1), s=16, edgecolors="none", zorder=3)
            ax.scatter(metric.natural_D1, metric.natural_D2, marker="X", s=78, color="#222222", edgecolor="white", lw=.6, zorder=7)
            maximum = curve.iloc[int(metric.maximum_index)]
            ax.scatter(maximum.D1, maximum.D2, s=32, facecolor=EDGE_RETENTION_CMAP(maximum.edge_retention),
                       edgecolor="#222222", lw=1., zorder=7)
            visible_x.extend(curve.D1.to_list() + [metric.natural_D1])
            visible_y.extend(curve.D2.to_list() + [metric.natural_D2])
            full = null[null.family == "random_rebuild"]
            sample = full.iloc[np.linspace(0, len(full)-1, min(DISPLAY_DRAWS, len(full)), dtype=int)]
            inset = ax.inset_axes([.74, .69, .23, .26])
            color = ps.NULL_MODEL_COLORS["full_random"]
            inset.scatter(sample.D1, sample.D2, s=4, color=color, alpha=.18, edgecolors="none", rasterized=True)
            inset.scatter(full.D1.mean(), full.D2.mean(), marker="+", s=26, color=color, lw=1.)
            inset.set_xlim(*proportional_limits(full.D1))
            inset.set_ylim(*proportional_limits(full.D2))
            inset.set_title("Random rebuild", fontsize=7, pad=2)
            inset.set_xlabel("Travel distance", fontsize=6, labelpad=1)
            inset.set_ylabel("Cell-state distance", fontsize=6, labelpad=1)
            inset.tick_params(labelsize=5.5, length=2, pad=1)
            inset.locator_params(axis="both", nbins=3)
            inset.grid(True, alpha=.25)
            inset.spines[["top", "right"]].set_visible(True)
            inset_limits[f"{geometry}__{config}"] = dict(x=inset.get_xlim(), y=inset.get_ylim(), shares_endpoint_transform=True)
            terminal_style.comparison_heading(ax, config, geometry, "ABCDEF"[3*row+col])
            if row == 1:
                ax.set_xlabel("Travel distance\n(fraction of endpoint cost span)")
            if col == 0:
                ax.set_ylabel("Cell-state distance\n(fraction of endpoint cost span)")
    for values, setter in ((visible_x, axes[0, 0].set_xlim), (visible_y, axes[0, 0].set_ylim)):
        low, high = min(values), max(values)
        pad = max(.04*(high-low), .04)
        setter(low-pad, high+pad)
    cb = fig.colorbar(front, cax=fig.add_axes([.934, .16, .013, .71]))
    cb.set_label("Biological-parent edge retention")
    handles = [Line2D([], [], marker="o", ls="", color=ps.NULL_MODEL_COLORS[family], alpha=.7,
                      markersize=5, markeredgecolor="none", label=label) for family, label in (
        ("first_cousin", "First-cousin shuffle"), ("second_cousin", "Second-cousin shuffle"),
        ("third_cousin", "Third-cousin shuffle"), ("full_random", "Random rebuild (insets)"))]
    handles += [Line2D([], [], marker="X", ls="", color="#222222", markersize=7, label="Natural lineage"),
                Line2D([], [], marker="o", ls="", color=ps.COLORS["blue"], markeredgecolor="#222222", markersize=5,
                       label="Maximum edge retention")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.49, .025), ncol=3, fontsize=8)
    limits = dict(x=axes[0, 0].get_xlim(), y=axes[0, 0].get_ylim(), random_rebuild_insets=inset_limits)
    save(fig, out, COMPARISON)
    return limits


def overlay(fronts, metrics, out, record):
    fig = plt.figure(figsize=(12.8, 8.8))
    first = fig.add_axes([.07, .59, .42, .30])
    axes = [first, fig.add_axes([.56, .59, .42, .30], sharex=first, sharey=first)]
    records = terminal_style.canonical_records(fronts, metrics)
    main = metrics[metrics.cohort == "base"]
    close_x = [r["closest_D1"] for r in records] + main.natural_D1.to_list()
    close_y = [r["closest_D2"] for r in records] + main.natural_D2.to_list()
    zoom = ((min(close_x)-.035, max(close_x)+.035), (min(close_y)-.08, max(close_y)+.08))
    for ax, geometry, title in zip(axes, GEOMETRIES, ("A   3D", "B   2D (XY)")):
        handles = terminal_style.overlay_panel(ax, fronts, metrics, geometry, zoom_limits=zoom, legend=False)
        ax.set_title(title, fontsize=10, pad=12)
    # Explicit union, not merely the limits of the last shared subplot.
    base = fronts[fronts.cohort == "base"]
    axes[0].set_xlim(min(base.D1.min(), main.natural_D1.min())-.04, max(base.D1.max(), main.natural_D1.max())+.04)
    axes[0].set_ylim(min(base.D2.min(), main.natural_D2.min())-.04, max(base.D2.max(), main.natural_D2.max())+.04)
    axes[1].set_ylabel("")
    fig.suptitle(f"Partial-forest fronts and canonical metrics | {record['edges']} shared edges", y=.985, fontsize=12)
    fig.text(.5, .947, "Overlay main axes and near-natural zoom limits shared", ha="center", fontsize=8, color="#555555")
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.52, .465), ncol=3, fontsize=8)
    terminal_style.canonical_summary(fig, fronts, metrics)
    limits = dict(x=axes[0].get_xlim(), y=axes[0].get_ylim(), zoom_x=zoom[0], zoom_y=zoom[1])
    save(fig, out, OVERLAY)
    return limits


def render(run):
    ps.configure()
    run = Path(run)
    archive = archive_previous(run, "figures")
    out = run / "figures"
    out.mkdir(exist_ok=True)
    record = json.loads((run / "analysis/provenance.json").read_text())
    fronts, metrics, clouds = (pd.read_csv(run / "analysis" / name) for name in ("fronts.csv", "metrics.csv", "null_clouds.csv"))
    comparison_limits = comparison(fronts, metrics, clouds, out, record)
    overlay_limits = overlay(fronts, metrics, out, record)
    manifest = dict(analysis_id=record["analysis_id"], comparison=COMPARISON, overlay=OVERLAY,
        comparison_rows=list(GEOMETRIES), comparison_columns=list(PRIMARY), retention_scale=[0, 1],
        comparison_geometry_labels="bold second line beneath each configuration title",
        comparison_title_tracking_counts=False,
        comparison_limits=comparison_limits, overlay_limits=overlay_limits,
        comparison_natural_to_maximum_connectors=False, overlay_closest_point_connectors=True,
        canonical_panel="C", canonical_metrics=terminal_style.canonical_records(fronts, metrics),
        previous_layout=str(archive.relative_to(run)) if archive else None,
        plotting_sources={str(p.relative_to(ROOT)): digest(p) for p in (
            Path(__file__), Path(ps.__file__), Path(terminal_style.__file__),
            ROOT / "terminal_pareto/fig2_fig3_ce_terminal_pareto.py")},
        files={p.name: digest(p) for p in out.iterdir() if p.suffix in (".png", ".pdf")})
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return out
