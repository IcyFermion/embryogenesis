"""Pooled full-tree Figure 7 (round fronts A, aggregate comparison B) and supplement S4.

Moved from ``full_tree_pareto/publication_build.py``, which keeps the sweep
settings, cache-building CLI, working-layout archive and release mechanics.
Inputs are replay-validated fronts, reference draws and endpoint transforms.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np

from publication import notation as nt
from publication import style as ps

PANEL_COLUMNS = 4
METHOD_PALETTE = ("#0072B2", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")
REFERENCE_COLORS = {
    "First-cousin shuffle": ps.SEMANTIC_COLORS["first_cousin_null"],
    "Internal-layer shuffle": "#8195A5",
    "Full assignment shuffle": "#9E6BA0",
    "Random rebuild": "#6F3B5C",
    "Terminal-clamped Gaussian reference": "#6B6ECF",
    "Second-cousin shuffle": ps.NULL_MODEL_COLORS["second_cousin"],
    "Third-cousin shuffle": ps.NULL_MODEL_COLORS["third_cousin"],
}
RETENTION_CMAP = LinearSegmentedColormap.from_list(
    "full_tree_retention", ["#17365D", ps.COLORS["blue"], "#72C7EC"])
ROUNDS_STEM = "fig7A_ce_full_tree_layerwise_rounds"
COLLECTIVE_STEM = "fig7B_ce_full_tree_collective"
SUPPLEMENT_STEM = "figs_ce_full_tree_heuristics_panel"


def axis_labels():
    return (f"Pooled travel distance, {nt.math('travel_axis')}",
            f"Cell-state distance, {nt.math('cell_state_axis')}")


def method_colors(data):
    return dict(zip(data.methods, METHOD_PALETTE))


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(Path(out) / f"{stem}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _display(values, draws):
    return values[np.linspace(0, len(values) - 1, min(draws, len(values)), dtype=int)]


def retention_front(ax, frame, transform, natural, color):
    x, y = transform.transform(frame.travel.to_numpy(), frame.state.to_numpy())
    order = np.argsort(x)
    ax.plot(x[order], y[order], color=color, lw=.7, alpha=.65, zorder=3)
    scatter = ax.scatter(x, y, c=frame.retention, cmap=RETENTION_CMAP, norm=Normalize(0, 1), s=9, linewidths=0, zorder=4)
    nx_, ny_ = transform.transform(*np.asarray(natural))
    ax.scatter(nx_, ny_, marker="X", s=42, color=ps.COLORS["black"], edgecolors="white", lw=.5, zorder=7)
    best = np.flatnonzero(frame.retention.to_numpy() == frame.retention.max())
    i = best[len(best)//2]
    ax.scatter(x[i], y[i], marker="o", s=28, c=[frame.retention.iloc[i]], cmap=RETENTION_CMAP, norm=Normalize(0, 1),
               edgecolors=ps.COLORS["black"], lw=.8, zorder=6)
    ax.grid(True, alpha=.3)
    ax.margins(.13)
    return scatter


def plot_rounds(data, out):
    """Figure 7A: one endpoint-normalized layerwise front per contraction round."""
    count = len(data.round_edges)
    rows = math.ceil(count / PANEL_COLUMNS)
    fig, axes = plt.subplots(rows, PANEL_COLUMNS, figsize=(7.15, 1.925 * rows), layout="constrained", squeeze=False)
    for ax in axes.flat[count:]:
        ax.remove()
    axes = axes.flat[:count]
    transforms = {}
    color = method_colors(data)[data.methods[0]]
    for i, ax in enumerate(axes, 1):
        scope = f"round_{i}"
        frame = data.fronts[(data.fronts.method == data.methods[0]) & (data.fronts.scope == scope)]
        transform = data.round_transforms[scope]
        row = data.layers.set_index("scope").loc[scope]
        scatter = retention_front(ax, frame, transform, [row.natural_travel, row.natural_state], color)
        ax.set_title(f"Round {i}" + (" (bottom)" if i == 1 else " (top)" if i == count else "")
                     + f"\n{int(row.edges)} edges", fontsize=8)
        ax.tick_params(labelsize=6.5)
        transforms[scope] = transform.metadata(axis_limits={"x": ax.get_xlim(), "y": ax.get_ylim()})
    cb = fig.colorbar(scatter, ax=list(axes), fraction=.025, pad=.015, aspect=32)
    cb.set_label("Natural edges retained", fontsize=8)
    xlabel, ylabel = axis_labels()
    fig.supxlabel(xlabel, fontsize=8.5)
    fig.supylabel(ylabel, fontsize=8.5)
    fig.text(.003, 1.015, "A", fontweight="bold", fontsize=11)
    save(fig, out, ROUNDS_STEM)
    return transforms


def plot_comparison(data, out, *, supplement=False):
    """Figure 7B (layerwise front with cousin shuffles) or S4 (all methods and references)."""
    transform = data.aggregate_transform
    colors = method_colors(data)
    draws = data.settings["display_draws"]
    fig, ax = plt.subplots(figsize=(7.15, 4.0 if not supplement else 5.2))
    fig.subplots_adjust(left=.10, right=.97, bottom=.14, top=.97)
    models = data.methods if supplement else data.methods[:1]
    references = data.supplement_references if supplement else data.main_references
    handles = []
    for method in models:
        frame = data.fronts[(data.fronts.method == method) & (data.fronts.scope == "aggregate") & data.fronts.nondominated]
        x, y = transform.transform(frame.travel.to_numpy(), frame.state.to_numpy())
        order = np.argsort(x)
        line, = ax.plot(x[order], y[order], color=colors[method], lw=1.8 if method == data.methods[0] else 1.2,
                        ls="-" if method in data.methods[:2] else "--", zorder=4, label=method)
        handles.append(line)
    for name in references:
        values = data.references[name]
        shown = _display(values, draws)
        x, y = transform.transform(shown[:, 0], shown[:, 1])
        ax.scatter(x, y, s=7, alpha=.13, color=REFERENCE_COLORS[name], edgecolors="none", rasterized=True, zorder=1)
        mx, my = transform.transform(*values.mean(axis=0))
        handles.append(ax.scatter(mx, my, marker="+", s=55, lw=1.3, color=REFERENCE_COLORS[name], label=name, zorder=5))
    x, y = transform.transform(*data.natural)
    handles.append(ax.scatter(x, y, marker="X", s=65, c=ps.COLORS["black"], edgecolors="white", lw=.55, zorder=7,
                              label=nt.label("natural_lineage")))
    inset_limits = None
    if not supplement:
        # Broad randomization is contextual: same endpoint transform, separately scaled axes.
        name = data.inset_reference
        values = data.references[name]
        shown = _display(values, draws)
        x, y = transform.transform(shown[:, 0], shown[:, 1])
        inset = ax.inset_axes([.66, .70, .31, .25])
        inset.scatter(x, y, s=6, alpha=.18, color=REFERENCE_COLORS[name], edgecolors="none", rasterized=True)
        inset.scatter(*transform.transform(*values.mean(axis=0)), marker="+", s=35, lw=1.1, color=REFERENCE_COLORS[name])
        inset.set_title(name, fontsize=7, pad=2)
        inset.set_xlabel("Pooled travel distance", fontsize=6, labelpad=1)
        inset.set_ylabel("Cell-state distance", fontsize=6, labelpad=1)
        inset.tick_params(labelsize=5.5, length=2, pad=1)
        inset.xaxis.set_major_locator(plt.MaxNLocator(4))
        inset.yaxis.set_major_locator(plt.MaxNLocator(3))
        inset.grid(True, alpha=.25)
        inset.margins(.10)
        inset_limits = {"x": inset.get_xlim(), "y": inset.get_ylim()}
        handles.append(Line2D([0], [0], marker="+", ls="", color=REFERENCE_COLORS[name], label=f"{name} (inset)"))
    ax.grid(True, alpha=.25)
    xlabel, ylabel = axis_labels()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.margins(.07)
    legend_options = {} if supplement else dict(bbox_to_anchor=(.66, .59), borderaxespad=0)
    ax.legend(handles=handles, loc="upper left", fontsize=6.6, frameon=True,
              facecolor="white", edgecolor="none", framealpha=.9, **legend_options)
    if not supplement:
        ax.text(-.10, 1.015, "B", transform=ax.transAxes, fontweight="bold", fontsize=11)
    metadata = transform.metadata(axis_limits={"x": ax.get_xlim(), "y": ax.get_ylim()})
    metadata["reference_models"] = list(references)
    if inset_limits:
        metadata["random_rebuild_inset"] = dict(axis_limits=inset_limits, shares_endpoint_transform=True)
    save(fig, out, SUPPLEMENT_STEM if supplement else COLLECTIVE_STEM)
    return metadata


def render_panels(data, out) -> dict:
    """All Figure 7/S4 panels; returns the display-transform record."""
    ps.configure()
    transforms = plot_rounds(data, out)
    transforms["main_collective"] = plot_comparison(data, out)
    transforms["supplement"] = plot_comparison(data, out, supplement=True)
    return transforms
