"""Shared publication style for the terminal Pareto figures.

This module contains presentation defaults only; it deliberately has no
dependency on the analysis code or data.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np


# Okabe–Ito-derived palette: distinguishable in common colour-vision deficiencies
# and still separable when printed in greyscale through the accompanying markers.
COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#222222",
    "grey": "#6F6F6F",
    "light_grey": "#D9D9D9",
}

SPECIES_COLORS = {
    "ce_protein": COLORS["blue"],
    "ce_rna": COLORS["orange"],
    "cb_rna": COLORS["green"],
}


def configure() -> None:
    """Apply compact, journal-friendly Matplotlib defaults."""
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.titleweight": "semibold",
            "axes.labelsize": 8.5,
            "axes.labelpad": 3,
            "axes.linewidth": 0.7,
            "axes.edgecolor": COLORS["black"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "legend.handlelength": 1.6,
            "legend.labelspacing": 0.3,
            "lines.linewidth": 1.5,
            "lines.solid_capstyle": "round",
            "grid.color": "#D6D6D6",
            "grid.linewidth": 0.45,
            "grid.alpha": 0.45,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig, png_path) -> None:
    """Export a high-resolution preview and an editable vector PDF."""
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, facecolor="white")
    fig.savefig(png_path.with_suffix(".pdf"), facecolor="white")


def color_ramp(base_color, n, lightness=(0.72, 0.05)):
    """Return ``n`` white-to-base shades that preserve a configuration hue."""
    if n <= 0:
        return []
    base = np.asarray(mcolors.to_rgb(base_color))
    mixes = np.linspace(lightness[0], lightness[1], n)
    return [mcolors.to_hex(base * (1 - mix) + np.ones(3) * mix) for mix in mixes]


def save_panel_crops(
    fig, axes, output_root, config_dirs, filename, extra_axes=None
) -> None:
    """Save each axis of a horizontal comparison figure as its own PNG."""
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    extras = extra_axes if extra_axes is not None else [None] * len(axes)
    for ax, extra_ax, dirname in zip(axes, extras, config_dirs):
        target = Path(output_root) / dirname / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        visible_axes = {ax}
        if extra_ax is not None:
            visible_axes.add(extra_ax)
        original_visibility = {item: item.get_visible() for item in fig.axes}
        for item in fig.axes:
            item.set_visible(item in visible_axes)
        boxes = [ax.get_tightbbox(renderer)]
        if extra_ax is not None:
            boxes.append(extra_ax.get_tightbbox(renderer))
        bbox = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(
            target, dpi=300, bbox_inches=bbox.expanded(1.12, 1.12),
            pad_inches=0.10, facecolor="white",
        )
        for item, visible in original_visibility.items():
            item.set_visible(visible)
