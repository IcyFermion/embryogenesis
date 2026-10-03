"""Reusable Matplotlib artists for endpoint-coordinate front figures.

Artists consume ``publication.contracts`` records only; they import neither
back end nor any figure-specific script.
"""

from __future__ import annotations

from matplotlib.lines import Line2D

from publication import notation as nt
from publication import style
from publication.canonical_panels import canonical_point

AXIS_LABELS = ("Travel distance\n(fraction of endpoint cost span)",
               "Cell-state distance\n(fraction of endpoint cost span)")
OVERLAY_LINESTYLES = ("-", "--", "-.")


def base_axes(ax):
    ax.axhline(0, color="#777777", lw=.7, ls=":", zorder=0)
    ax.axvline(0, color="#777777", lw=.7, ls=":", zorder=0)
    ax.grid(True, alpha=.32)
    ax.spines[["top", "right"]].set_visible(True)


def overlay_handles(data):
    handles = [Line2D([], [], color=style.SPECIES_COLORS[config], ls=ls, label=nt.CONFIG_LABELS[config])
               for config, ls in zip(data.configs, OVERLAY_LINESTYLES)]
    handles.extend([
        Line2D([], [], marker="X", ls="", color="#222222", label=nt.label("natural_lineage")),
        Line2D([], [], marker="o", ls="", color="white", markeredgecolor="#222222", label="Maximum edge retention"),
        Line2D([], [], marker="D", ls=(0, (2, 2)), color="#666666", markerfacecolor="white",
               label=nt.heading("closest_point"))])
    return handles


def _closest_overlay(ax, data, geometry, *, main):
    """Fronts, natural lineages, P* diamonds and natural-to-P* connectors."""
    for config, ls in zip(data.configs, OVERLAY_LINESTYLES):
        curve, m = data.select(geometry, config)
        color = style.SPECIES_COLORS[config]
        ax.plot(curve.D1, curve.D2, color=color, lw=1.6 if main else 1.4, ls=ls)
        ax.scatter(m.natural_D1, m.natural_D2, marker="X", s=70 if main else 45, color=color,
                   edgecolor="white", lw=.6 if main else .5, zorder=7 if main else 4)
        if main:
            maximum = curve.iloc[int(m.maximum_index)]
            ax.scatter(maximum.D1, maximum.D2, s=34, color=color, edgecolor="#222222", lw=1, zorder=6)
        closest, _ = canonical_point(curve, m)
        connector, = ax.plot([m.natural_D1, closest.D1], [m.natural_D2, closest.D2], color=color,
                             lw=.9 if main else 1, ls=(0, (2, 2)), zorder=5 if main else 3)
        connector.set_gid(f"natural_to_closest:{geometry}:{config}")
        ax.scatter(closest.D1, closest.D2, marker="D", s=30 if main else 25, facecolor="white",
                   edgecolor=color, lw=1.1 if main else 1, zorder=8 if main else 5)


def overlay_panel(ax, data, geometry, *, zoom_limits):
    """One geometry with every configuration's front and a near-natural zoom inset.

    Main-axis limits are set by the caller so paired panels share them.
    """
    base_axes(ax)
    _closest_overlay(ax, data, geometry, main=True)
    ax.set(xlabel=AXIS_LABELS[0], ylabel=AXIS_LABELS[1])
    # Zoom shares the exact coordinates and anchors of the main axes.
    inset = ax.inset_axes([.48, .40, .47, .47])
    _closest_overlay(inset, data, geometry, main=False)
    inset.set_xlim(*zoom_limits[0])
    inset.set_ylim(*zoom_limits[1])
    inset.set_title("Near natural lineages", fontsize=7.5)
    inset.tick_params(labelsize=6)
    inset.grid(True, alpha=.3)
    inset.spines[["top", "right"]].set_visible(True)
    return inset
