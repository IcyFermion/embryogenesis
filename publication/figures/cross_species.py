"""Cross-species layouts: terminal Figures 8/9 and partial-forest Figures 10/11.

Scope differences (cohort, titles, reference families, display subsets, zoom
rule) are explicit arguments or contract fields supplied by the family
specification, never inferred from the data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np

from publication import notation as nt
from publication import style
from publication.artists import (AXIS_LABELS, base_axes, comparison_heading, overlay_handles,
                                 overlay_panel, proportional_limits)
from publication.canonical_panels import canonical_records, canonical_summary

FAMILY_LABELS = {"first_cousin": "First-cousin shuffle", "second_cousin": "Second-cousin shuffle",
                 "third_cousin": "Third-cousin shuffle"}


def save(fig, out: Path, stem: str):
    style.save_figure(fig, Path(out) / f"{stem}.png")
    plt.close(fig)


def _display(frame, draws):
    """Deterministic evenly spaced display subset; ``None`` shows every row."""
    if draws is None:
        return frame
    return frame.iloc[np.linspace(0, len(frame) - 1, min(draws, len(frame)), dtype=int)]


def compose_comparison(data, *, title: str, geometries=None, front_alpha=.5, front_size=18,
                       colorbar_label="Edge retention"):
    """Figure 2-style fronts with reference clouds: rows are geometries, columns configurations."""
    refs = data.references
    geometries = tuple(geometries or data.geometries)
    combined = len(geometries) == 2
    fig, axes = plt.subplots(len(geometries), len(data.configs), figsize=(14.4, 9.8 if combined else 5.1),
                             sharex=True, sharey=True, squeeze=False)
    if combined:
        fig.subplots_adjust(left=.06, right=.91, bottom=.16, top=.87, hspace=.32, wspace=.13)
    else:
        fig.subplots_adjust(left=.06, right=.91, bottom=.24, top=.80, wspace=.13)
    fig.suptitle(title.format(edges=data.edges), y=.97, fontsize=12)
    visible_x, visible_y, inset_limits = [], [], {}
    for row, geometry in enumerate(geometries):
        for col, config in enumerate(data.configs):
            ax = axes[row, col]
            base_axes(ax)
            curve, metric = data.select(geometry, config)
            null = data.cloud(geometry, config)
            for family in refs["main"]:
                full = null[null.family == family]
                sample = _display(full, refs["display_draws"])
                color = style.NULL_MODEL_COLORS[family]
                ax.scatter(sample.D1, sample.D2, s=7, color=color, alpha=.19, edgecolors="none",
                           rasterized=True, zorder=1)
                # The first-cousin mean is the cached (analytic) canonical null mean.
                mean = ((metric.null_D1, metric.null_D2) if family == "first_cousin"
                        else (full.D1.mean(), full.D2.mean()))
                ax.scatter(*mean, marker="+", s=34, color=color, lw=1.1, zorder=4)
                visible_x.extend(sample.D1)
                visible_y.extend(sample.D2)
            ax.plot(curve.D1, curve.D2, color=style.COLORS["blue"], lw=1.3, alpha=front_alpha, zorder=2)
            front = ax.scatter(curve.D1, curve.D2, c=curve.edge_retention, cmap=style.EDGE_RETENTION_CMAP,
                               norm=Normalize(0, 1), s=front_size, edgecolors="none", zorder=3)
            ax.scatter(metric.natural_D1, metric.natural_D2, marker="X", s=78, color="#222222",
                       edgecolor="white", lw=.6, zorder=7)
            maximum = curve.iloc[int(metric.maximum_index)]
            ax.scatter(maximum.D1, maximum.D2, s=32, facecolor=style.EDGE_RETENTION_CMAP(maximum.edge_retention),
                       edgecolor="#222222", lw=1., zorder=7)
            visible_x.extend(curve.D1.to_list() + [metric.natural_D1])
            visible_y.extend(curve.D2.to_list() + [metric.natural_D2])
            full = null[null.family == refs["inset"]]
            sample = _display(full, refs["display_draws"])
            color = style.NULL_MODEL_COLORS["full_random"]
            inset = ax.inset_axes([.74, .69, .23, .26])
            inset.scatter(sample.D1, sample.D2, s=4, color=color, alpha=.18, edgecolors="none", rasterized=True)
            inset.scatter(full.D1.mean(), full.D2.mean(), marker="+", s=26, color=color, lw=1.)
            inset.set_xlim(*proportional_limits(full.D1))
            inset.set_ylim(*proportional_limits(full.D2))
            inset.set_title(refs["inset_title"], fontsize=7, pad=2)
            inset.set_xlabel("Travel distance", fontsize=6, labelpad=1)
            inset.set_ylabel("Cell-state distance", fontsize=6, labelpad=1)
            inset.tick_params(labelsize=5.5, length=2, pad=1)
            inset.locator_params(axis="both", nbins=3)
            inset.grid(True, alpha=.25)
            inset.spines[["top", "right"]].set_visible(True)
            inset_limits[f"{geometry}__{config}"] = dict(x=inset.get_xlim(), y=inset.get_ylim(),
                                                         shares_endpoint_transform=True)
            comparison_heading(ax, config, geometry, "ABCDEF"[3 * row + col] if combined else "ABC"[col])
            if row == len(geometries) - 1:
                ax.set_xlabel(AXIS_LABELS[0])
            if col == 0:
                ax.set_ylabel(AXIS_LABELS[1])
    for values, setter in ((visible_x, axes[0, 0].set_xlim), (visible_y, axes[0, 0].set_ylim)):
        low, high = min(values), max(values)
        pad = max(.04 * (high - low), .04)
        setter(low - pad, high + pad)
    bar = fig.colorbar(front, cax=fig.add_axes([.934, .16 if combined else .24, .013, .71 if combined else .56]))
    bar.set_label(colorbar_label)
    handles = [Line2D([], [], marker="o", ls="", color=style.NULL_MODEL_COLORS[family], alpha=.7,
                      markersize=5, markeredgecolor="none", label=FAMILY_LABELS[family]) for family in refs["main"]]
    handles.append(Line2D([], [], marker="o", ls="", color=style.NULL_MODEL_COLORS["full_random"], alpha=.7,
                          markersize=5, markeredgecolor="none", label=refs["inset_legend"]))
    handles += [Line2D([], [], marker="X", ls="", color="#222222", markersize=7, label=nt.label("natural_lineage")),
                Line2D([], [], marker="o", ls="", color=style.COLORS["blue"], markeredgecolor="#222222",
                       markersize=5, label="Maximum edge retention")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.49, .025), ncol=3, fontsize=8)
    limits = dict(x=axes[0, 0].get_xlim(), y=axes[0, 0].get_ylim(), reference_insets=inset_limits)
    return fig, limits


def comparison_figure(data, out: Path, stem: str, **options) -> dict:
    fig, limits = compose_comparison(data, **options)
    save(fig, out, stem)
    return limits


def overlay_zoom(data, *, include_closest: bool, geometries=None):
    """Near-natural zoom limits shared by the requested geometries."""
    geometries = geometries or data.geometries
    rows = data.metrics[data.metrics.geometry.isin(geometries)]
    xs, ys = list(rows.natural_D1), list(rows.natural_D2)
    if include_closest:
        records = [r for r in canonical_records(data) if r["geometry"] in geometries]
        xs += [r["closest_D1"] for r in records]
        ys += [r["closest_D2"] for r in records]
    return (min(xs) - .035, max(xs) + .035), (min(ys) - .08, max(ys) + .08)


def _main_limits(ax, data, geometries):
    fronts = data.fronts[data.fronts.geometry.isin(geometries)]
    metrics = data.metrics[data.metrics.geometry.isin(geometries)]
    xs = list(fronts.D1) + list(metrics.natural_D1)
    ys = list(fronts.D2) + list(metrics.natural_D2)
    ax.set_xlim(min(xs) - .04, max(xs) + .04)
    ax.set_ylim(min(ys) - .04, max(ys) + .04)


def compose_overlay(data, *, title: str, zoom_includes_closest: bool):
    """Panels A-B: paired geometry overlays; panel C: canonical metric rows."""
    fig = plt.figure(figsize=(12.8, 8.8))
    first = fig.add_axes([.07, .59, .42, .30])
    axes = [first, fig.add_axes([.56, .59, .42, .30], sharex=first, sharey=first)]
    zoom = overlay_zoom(data, include_closest=zoom_includes_closest)
    for ax, geometry, letter in zip(axes, data.geometries, "AB"):
        overlay_panel(ax, data, geometry, zoom_limits=zoom)
        ax.set_title(f"{letter}   {nt.GEOMETRY_LABELS[geometry]}", fontsize=10, pad=12)
    # Shared main axes: union of every displayed front and natural lineage.
    _main_limits(axes[0], data, data.geometries)
    axes[1].set_ylabel("")
    fig.suptitle(f"{title} | {data.edges} shared edges", y=.985, fontsize=12)
    fig.text(.5, .947, "Overlay main axes and near-natural zoom limits shared", ha="center", fontsize=8, color="#555555")
    fig.legend(handles=overlay_handles(data), loc="lower center", bbox_to_anchor=(.52, .465), ncol=3, fontsize=8)
    canonical_summary(fig, data)
    limits = dict(x=axes[0].get_xlim(), y=axes[0].get_ylim(), zoom_x=zoom[0], zoom_y=zoom[1])
    return fig, limits


def overlay_figure(data, out: Path, stem: str, **options) -> dict:
    fig, limits = compose_overlay(data, **options)
    save(fig, out, stem)
    return limits


def compose_single_overlay(data, geometry: str):
    """One-geometry overlay export with its own legend and zoom (terminal reuse asset)."""
    fig, ax = plt.subplots(figsize=(7.15, 4.85))
    fig.subplots_adjust(left=.12, right=.96, bottom=.17, top=.86)
    overlay_panel(ax, data, geometry, zoom_limits=overlay_zoom(data, include_closest=False, geometries=(geometry,)))
    _main_limits(ax, data, (geometry,))
    ax.legend(handles=overlay_handles(data), loc="upper right", bbox_to_anchor=(.98, .32), fontsize=7)
    fig.suptitle(f"Pooled terminal fronts | {data.edges} shared edges", y=.97, fontsize=10.5)
    subtitle = "3D: pooled distances" if geometry == "raw3d" else "2D (XY): axial coordinates omitted"
    fig.text(.5, .91, subtitle, ha="center", fontsize=7.5, color="#555555")
    return fig


def compose_tracking(data):
    """Individual-embryo optima transferred to pooled 3D costs and pooled anchors."""
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.3), sharex=True, sharey=True)
    fig.subplots_adjust(left=.06, right=.98, top=.77, bottom=.28, wspace=.16)
    colors = [style.COLORS["blue"], style.COLORS["orange"], style.COLORS["green"]]
    all_x, all_y = [], []
    for i, (ax, config) in enumerate(zip(axes, data.configs)):
        base_axes(ax)
        curve, m = data.select("raw3d", config)
        ax.plot(curve.D1, curve.D2, color="#222222", lw=1.8, label="Pooled optimum")
        ax.scatter(m.natural_D1, m.natural_D2, color="#222222", marker="X", s=55, zorder=5)
        all_x.extend(curve.D1.to_list() + [m.natural_D1])
        all_y.extend(curve.D2.to_list() + [m.natural_D2])
        for j, (label, x, y) in enumerate(data.extras["tracking_transfers"][config]):
            all_x.extend(x)
            all_y.extend(y)
            ax.plot(x, y, color=colors[j], lw=1.2, ls="--", alpha=.85, label=label)
        ax.set_title(f"{'ABC'[i]}   {nt.CONFIG_LABELS[config]}", fontsize=9, pad=10)
        ax.set_xlabel("Travel distance\n(pooled endpoint coordinates)")
        ax.legend(loc="upper right", bbox_to_anchor=(.98, .65), fontsize=7)
    axes[0].set_ylabel("Cell-state distance\n(pooled endpoint coordinates)")
    axes[0].set_xlim(min(all_x) - .04, max(all_x) + .04)
    axes[0].set_ylim(min(all_y) - .04, max(all_y) + .04)
    fig.suptitle("Tracking sensitivity: all assignments evaluated using pooled travel", fontsize=11, y=.95)
    fig.text(.5, .87, f"{data.edges} shared edges | 3D pooled distances", ha="center", fontsize=8, color="#555555")
    fig.text(.5, .06, "Solid: pooled optimization. Dashed: individual-embryo optimization transferred to the pooled objective.\n"
             "Every panel uses its pooled front's anchors; molecular matrices and biological-parent capacities stay fixed.",
             ha="center", fontsize=8)
    return fig
