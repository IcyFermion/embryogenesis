"""Terminal Figure 6B-C: within-type fronts and the type-restricted aggregate.

Moved from ``terminal_pareto/fig6bc_ce_within_type.py``, which keeps the
within-type analysis (block sweeps) and its CLI. Inputs are the cached
within-type fronts, summary and aggregate tables plus endpoint metadata.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

from publication import style as ps

TYPE_ORDER = [
    "programmed_death", "epithelium", "neuron", "muscle",
    "alimentary", "glial",
]
TYPE_COLORS = {
    "programmed_death": ps.COLORS["purple"],
    "epithelium": ps.COLORS["orange"],
    "neuron": ps.COLORS["blue"],
    "muscle": ps.COLORS["vermillion"],
    "alimentary": ps.COLORS["green"],
    "glial": ps.COLORS["sky"],
}
TRAVEL_COLOR = ps.SEMANTIC_COLORS["travel"]
CELL_STATE_COLOR = ps.SEMANTIC_COLORS["cell_state"]
UNRESTRICTED_COLOR = "#444444"
TYPE_RESTRICTED_COLOR = ps.COLORS["orange"]


def unique_front(data):
    """Remove consecutive duplicate display coordinates."""
    xy = data[["travel_sigma_per_cell", "state_sigma_per_cell"]].to_numpy()
    keep = np.r_[True, np.any(np.abs(np.diff(xy, axis=0)) > 1e-12, axis=1)]
    return data.loc[keep]


def save_panel_crops(fig, upper_axes, lower_ax, b_artists, c_artists,
                     out_dir):
    """Export Figure 6B and 6C independently from their shared canvas."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_inches = fig.dpi_scale_trans.inverted()
    width, height = fig.get_size_inches()

    upper_box = Bbox.union(
        [axis.get_tightbbox(renderer) for axis in upper_axes]
        + [artist.get_tightbbox(renderer) for artist in b_artists]
    ).transformed(to_inches)
    lower_box = Bbox.union(
        [lower_ax.get_tightbbox(renderer)]
        + [artist.get_tightbbox(renderer) for artist in c_artists]
    ).transformed(to_inches)
    crops = [
        ("fig6B_ce_within_type_fronts", upper_axes, b_artists,
         Bbox.from_extents(0.0, max(0.0, upper_box.y0 - 0.04),
                           width, min(height, upper_box.y1 + 0.04))),
        ("fig6C_ce_type_restricted_aggregate", [lower_ax], c_artists,
         Bbox.from_extents(0.0, max(0.0, lower_box.y0 - 0.04),
                           width, min(height, lower_box.y1 + 0.04))),
    ]

    original_axis_visibility = {axis: axis.get_visible() for axis in fig.axes}
    figure_artists = list(b_artists) + list(c_artists)
    original_artist_visibility = {
        artist: artist.get_visible() for artist in figure_artists
    }
    for stem, visible_axes, visible_artists, crop in crops:
        visible_axes = set(visible_axes)
        visible_artists = set(visible_artists)
        for axis in fig.axes:
            axis.set_visible(axis in visible_axes)
        for artist in figure_artists:
            artist.set_visible(artist in visible_artists)
        fig.savefig(out_dir / f"{stem}.pdf", facecolor="white",
                    bbox_inches=crop, pad_inches=0)
        fig.savefig(out_dir / f"{stem}.png", dpi=300, facecolor="white",
                    bbox_inches=crop, pad_inches=0)

    for axis, was_visible in original_axis_visibility.items():
        axis.set_visible(was_visible)
    for artist, was_visible in original_artist_visibility.items():
        artist.set_visible(was_visible)


def plot_within_type(front_df, summary_df, aggregate, endpoints, *, out_dir,
                     display_mode="null_sd"):
    """Render the within-type panels used below the retention heatmap."""
    if display_mode not in {"null_sd", "endpoint"}:
        raise ValueError(f"Unsupported Figure 6 display mode: {display_mode}")
    endpoint_display = display_mode == "endpoint"
    if endpoint_display:
        local_x, local_y = "endpoint_travel", "endpoint_state"
        natural_x = "endpoint_natural_travel"
        natural_y = "endpoint_natural_state"
        local_xlabel = "Travel distance\n(own endpoint span)"
        local_ylabel = "Cell-state distance\n(own endpoint span)"
        global_x, global_y = "global_endpoint_travel", "global_endpoint_state"
        restricted_x = "restricted_endpoint_travel"
        restricted_y = "restricted_endpoint_state"
        aggregate_natural = (
            endpoints["endpoint_natural_travel"],
            endpoints["endpoint_natural_state"],
        )
        aggregate_xlabel = "Travel distance (unrestricted endpoint span)"
        aggregate_ylabel = "Cell-state distance\n(unrestricted endpoint span)"
    else:
        local_x, local_y = "travel_sigma_per_cell", "state_sigma_per_cell"
        natural_x = natural_y = None
        local_xlabel = "Travel-distance change\n(global null SD per cell)"
        local_ylabel = "Cell-state-distance change\n(global null SD per cell)"
        global_x, global_y = "global_travel_sigma", "global_state_sigma"
        restricted_x = "restricted_travel_sigma"
        restricted_y = "restricted_state_sigma"
        aggregate_natural = (0.0, 0.0)
        aggregate_xlabel = "Travel-distance change (global null SD)"
        aggregate_ylabel = "Cell-state-distance change\n(global null SD)"

    fig = plt.figure(figsize=(7.15, 5.75))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.56)
    upper = outer[0].subgridspec(2, 3, hspace=0.58, wspace=0.34)
    axes = [fig.add_subplot(upper[i, j]) for i in range(2) for j in range(3)]

    for ax, cell_type in zip(axes, TYPE_ORDER):
        data = unique_front(front_df[front_df["type"] == cell_type])
        row = summary_df.set_index("type").loc[cell_type]
        color = TYPE_COLORS[cell_type]
        ax.plot(data[local_x], data[local_y],
                color=color, lw=1.45, zorder=2)
        ax.scatter(data[local_x], data[local_y],
                   s=5, color=color, alpha=0.45, edgecolors="none", zorder=3)
        local_natural_x = (float(data.iloc[0][natural_x])
                           if natural_x is not None else 0.0)
        local_natural_y = (float(data.iloc[0][natural_y])
                           if natural_y is not None else 0.0)
        # Natural lineage stays visible when it coincides with an optimum marker.
        ax.scatter([local_natural_x], [local_natural_y],
                   marker="X", s=34, color="#222222",
                   edgecolor="white", lw=0.5, zorder=6)
        ax.scatter([data.iloc[-1][local_x]],
                   [data.iloc[-1][local_y]],
                   marker="o", s=22, facecolor=TRAVEL_COLOR,
                   edgecolor="#222222", lw=0.45, zorder=5)
        ax.scatter([data.iloc[0][local_x]],
                   [data.iloc[0][local_y]],
                   marker="s", s=22, facecolor=CELL_STATE_COLOR,
                   edgecolor="#222222", lw=0.45, zorder=5)
        ax.axhline(0, color="#999999", lw=0.55, ls=":")
        ax.axvline(0, color="#999999", lw=0.55, ls=":")
        title = cell_type.replace("_", " ")
        if endpoint_display and not bool(data.iloc[0]["endpoint_valid"]):
            title += "\n(null-SD fallback; degenerate endpoints)"
        ax.set_title(f"{title} (n={int(row['n'])})", fontsize=7.4,
                     color="#222222")
        ax.grid(True, alpha=0.22)
        ax.tick_params(labelsize=6.1)

    for ax in axes[3:]:
        ax.set_xlabel(local_xlabel, fontsize=6.7)
    for ax in (axes[0], axes[3]):
        ax.set_ylabel(local_ylabel, fontsize=6.7)
    b_letter = fig.text(0.012, 0.975, "B", fontsize=10, fontweight="bold",
                        ha="left", va="top")
    b_heading = fig.text(
        0.055, 0.972, "Pareto fronts within terminal cell types",
        fontsize=9.2, fontweight="semibold", ha="left", va="top")

    ax = fig.add_subplot(outer[1])
    ax.plot(aggregate[global_x], aggregate[global_y],
            color=UNRESTRICTED_COLOR, lw=1.8,
            label="Unrestricted assignment",
            zorder=3)
    ax.plot(aggregate[restricted_x], aggregate[restricted_y],
            color=TYPE_RESTRICTED_COLOR, lw=1.8,
            label="Assignments restricted within cell type", zorder=3)
    ax.scatter([aggregate_natural[0]], [aggregate_natural[1]],
               marker="X", s=48, color="#222222",
               edgecolor="white", lw=0.55, zorder=6, label="Natural lineage")
    for xcol, ycol, color in [
        (global_x, global_y, UNRESTRICTED_COLOR),
        (restricted_x, restricted_y, TYPE_RESTRICTED_COLOR),
    ]:
        ax.scatter([aggregate.iloc[-1][xcol]], [aggregate.iloc[-1][ycol]],
                   marker="o", s=27, facecolor=color, edgecolor="#222222",
                   lw=0.5, zorder=5)
        ax.scatter([aggregate.iloc[0][xcol]], [aggregate.iloc[0][ycol]],
                   marker="s", s=27, facecolor=color, edgecolor="#222222",
                   lw=0.5, zorder=5)
    ax.axhline(0, color="#888888", lw=0.6, ls=":")
    ax.axvline(0, color="#888888", lw=0.6, ls=":")
    ax.set_xlabel(aggregate_xlabel)
    ax.set_ylabel(aggregate_ylabel)
    ax.set_title("Aggregate consequence of forbidding assignments across cell types",
                 loc="left", pad=7)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=6.4)
    note = (f"Endpoint penalty of type restriction:  "
            f"travel +{endpoints['travel_penalty_sigma']:.2f}σ   ·   "
            f"cell state +{endpoints['state_penalty_sigma']:.2f}σ")
    ax.text(0.015, 0.035, note, transform=ax.transAxes, fontsize=6.5,
            color="#555555", ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="#CCCCCC", lw=0.5, alpha=0.9))
    c_letter = fig.text(0.012, 0.385, "C", fontsize=10, fontweight="bold",
                        ha="left", va="top")

    endpoint_handles = [
        Line2D([0], [0], marker="o", ls="", ms=4.7,
               markerfacecolor=TRAVEL_COLOR, markeredgecolor="#222222",
               label="Travel optimum"),
        Line2D([0], [0], marker="s", ls="", ms=4.7,
               markerfacecolor=CELL_STATE_COLOR, markeredgecolor="#222222",
               label="Cell-state optimum"),
    ]
    endpoint_legend = fig.legend(
        handles=endpoint_handles, loc="upper right", ncol=2,
        bbox_to_anchor=(0.975, 0.975), fontsize=6.2, frameon=False)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.09)
    save_panel_crops(
        fig, axes, ax,
        b_artists=[b_letter, b_heading, endpoint_legend],
        c_artists=[c_letter], out_dir=out_dir,
    )
    plt.close(fig)
