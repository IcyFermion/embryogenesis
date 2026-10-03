"""Cross-species overlay/canonical layouts shared by terminal Figure 9 and full-tree Figure 11.

Scope differences (cohort, title, zoom rule) are explicit arguments supplied
by the family specification, never inferred from the data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from publication import notation as nt
from publication import style
from publication.artists import overlay_handles, overlay_panel
from publication.canonical_panels import canonical_records, canonical_summary


def overlay_zoom(data, *, include_closest: bool):
    """Near-natural zoom limits shared by both geometries."""
    xs, ys = list(data.metrics.natural_D1), list(data.metrics.natural_D2)
    if include_closest:
        records = canonical_records(data)
        xs += [r["closest_D1"] for r in records]
        ys += [r["closest_D2"] for r in records]
    return (min(xs) - .035, max(xs) + .035), (min(ys) - .08, max(ys) + .08)


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
    xs = list(data.fronts.D1) + list(data.metrics.natural_D1)
    ys = list(data.fronts.D2) + list(data.metrics.natural_D2)
    axes[0].set_xlim(min(xs) - .04, max(xs) + .04)
    axes[0].set_ylim(min(ys) - .04, max(ys) + .04)
    axes[1].set_ylabel("")
    fig.suptitle(f"{title} | {data.edges} shared edges", y=.985, fontsize=12)
    fig.text(.5, .947, "Overlay main axes and near-natural zoom limits shared", ha="center", fontsize=8, color="#555555")
    fig.legend(handles=overlay_handles(data), loc="lower center", bbox_to_anchor=(.52, .465), ncol=3, fontsize=8)
    canonical_summary(fig, data)
    limits = dict(x=axes[0].get_xlim(), y=axes[0].get_ylim(), zoom_x=zoom[0], zoom_y=zoom[1])
    return fig, limits


def overlay_figure(data, out: Path, stem: str, **options) -> dict:
    fig, limits = compose_overlay(data, **options)
    style.save_figure(fig, Path(out) / f"{stem}.png")
    plt.close(fig)
    return limits
