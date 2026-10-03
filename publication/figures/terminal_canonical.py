"""Canonical-metric drawings: Figure 1 amendment, Figure 5 rows and Figure S1.

Moved from ``terminal_pareto/fig5_figs1_ce_canonical_summary.py`` (which keeps
its CLI and validated loader). Inputs are validated canonical-metric tables;
persisted column names (``d_lp``, ``d_np``, ``u_lineage_lp``) are unchanged.
Every displayed symbol comes from ``publication.notation``.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

from publication import notation as nt
from publication import style as ps


def _points():
    """Travel optimum, cell-state optimum, natural lineage, null mean and P* symbols."""
    return tuple(nt.symbol(key) for key in (
        "travel_optimum", "cell_state_optimum", "natural_lineage", "null_mean", "closest_point"))


REGION_MARKERS = {
    "root": "*",
    "AB": "P",
    "ABa": "o",
    "ABp": "^",
    "P1": "s",
}
REGION_LABELS = {
    "root": "Root",
    "AB": "AB",
    "ABa": "ABa lineage",
    "ABp": "ABp lineage",
    "P1": "P1 lineage",
}
REGION_ORDER = ["root", "AB", "ABa", "ABp", "P1"]
MAJOR_SUBTREES = ["P0", "AB", "ABa", "ABp", "P1"]
MAJOR_MARKERS = {"P0": "*", "AB": "P", "ABa": "o", "ABp": "^", "P1": "s"}
U_COLOR = ps.SEMANTIC_COLORS["canonical_u"]
LP_COLOR = ps.SEMANTIC_COLORS["lineage_front_distance"]
NP_COLOR = ps.SEMANTIC_COLORS["first_cousin_null"]
TRAVEL_COLOR = ps.SEMANTIC_COLORS["travel"]
STATE_COLOR = ps.SEMANTIC_COLORS["cell_state"]


def annotate_abplpp_pair(ax, data, u_col, distance_col, color):
    """Label ABplpp/ABplppp jointly only when they truly coincide."""
    nested = data[data["subtree"].isin(["ABplpp", "ABplppp"])]
    if len(nested) != 2:
        return
    coordinates = nested[[u_col, distance_col]].to_numpy(dtype=float)
    # Treat points as one visual anchor only when both coordinates are close
    # on the scale of the full panel. The legacy d_LP values differ slightly
    # but share one marker-sized location; the pooled values are widely
    # separated in u and must be labeled individually.
    u_tolerance = max(0.01, 0.02 * float(np.ptp(data[u_col])))
    distance_scale = max(float(np.abs(data[distance_col]).max()), 0.1)
    distance_tolerance = 0.02 * distance_scale
    visually_coincident = (
        abs(coordinates[0, 0] - coordinates[1, 0]) <= u_tolerance
        and abs(coordinates[0, 1] - coordinates[1, 1]) <= distance_tolerance
    )
    if visually_coincident:
        ax.annotate(
            "ABplpp(p)", tuple(coordinates.mean(axis=0)),
            xytext=(5, 4), textcoords="offset points", fontsize=6.1,
            color=color, ha="left", va="bottom", zorder=8,
        )
        return
    offsets = {
        "ABplpp": (5, -5, "left", "top"),
        "ABplppp": (5, 5, "left", "bottom"),
    }
    for row in nested.itertuples():
        dx, dy, ha, va = offsets[row.subtree]
        ax.annotate(
            row.subtree, (getattr(row, u_col), getattr(row, distance_col)),
            xytext=(dx, dy), textcoords="offset points", fontsize=6.1,
            color=color, ha=ha, va=va, zorder=8,
        )


def add_panel_letter(ax, letter, x=-0.16, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", ha="left", va="top", clip_on=False)


def add_colored_fraction(ax, prefix_x, fraction_x, center_y, prefix,
                         numerator, denominator, color, half_width,
                         fontsize=6.9):
    """Draw a compact fraction with a colored numerator and grey baseline."""
    ax.text(prefix_x, center_y, rf"${prefix}=$", color=color,
            fontsize=fontsize, ha="right", va="center")
    ax.text(fraction_x, center_y + 0.050, rf"${numerator}$", color=color,
            fontsize=fontsize, ha="center", va="center")
    ax.plot([fraction_x - half_width, fraction_x + half_width],
            [center_y + 0.006, center_y + 0.006],
            color="#777777", lw=0.65, solid_capstyle="butt", zorder=7)
    ax.text(fraction_x, center_y - 0.047, rf"${denominator}$",
            color="#777777", fontsize=fontsize, ha="center", va="center")


def plot_definition_r(ax, panel_letter="A"):
    """Supplementary symbolic definition of u and cousin-relative r."""
    T, S, L, C, P = _points()
    t = np.linspace(0.0, 1.0, 15)
    x = 0.13 + 0.72 * t
    y = 0.88 - 0.72 * np.sqrt(t)
    p_idx = 6
    px, py = x[p_idx], y[p_idx]
    lx, ly = px + 0.095, py + 0.095
    nx, ny = px + 0.235, py + 0.235

    # Objective axes and the sampled symbolic Pareto front.
    ax.annotate("", xy=(1.01, 0.0), xytext=(0.055, 0.0),
                annotation_clip=False,
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8))
    ax.annotate("", xy=(0.055, 1.00), xytext=(0.055, 0.0),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8))
    ax.text(1.01, -0.055, "Travel distance", ha="right", va="top",
            fontsize=7, clip_on=False)
    ax.text(0.015, 0.99, "Cell-state\ndistance", ha="right", va="top",
            fontsize=7)
    ax.plot(x, y, color="#A7A7A7", lw=1.3, zorder=1)
    ax.scatter(x, y, s=9, facecolor="white", edgecolor="#777777",
               lw=0.45, zorder=2)

    # u is normalized arc length from the travel optimum to the closest
    # attainable (sampled) front assignment P*.
    ax.plot(x[:p_idx + 1], y[:p_idx + 1], color=U_COLOR,
            lw=2.4, zorder=3)
    ax.scatter([x[0]], [y[0]], s=31, marker="o",
               facecolor=TRAVEL_COLOR, edgecolor="#222222",
               lw=0.5, zorder=5)
    ax.scatter([x[-1]], [y[-1]], s=31, marker="o",
               facecolor=STATE_COLOR, edgecolor="#222222",
               lw=0.5, zorder=5)
    ax.scatter([px], [py], s=29, marker="o", facecolor="#222222",
               edgecolor="white", lw=0.5, zorder=6)
    ax.annotate(f"Travel optimum  {nt.plain('travel_optimum')}", (x[0], y[0]), xytext=(5, 3),
                textcoords="offset points", fontsize=6.2, ha="left",
                va="bottom")
    ax.annotate(nt.plain("cell_state_optimum"), (x[-1], y[-1]), xytext=(4, 1),
                textcoords="offset points", fontsize=6.5, ha="left",
                va="center", fontweight="bold")
    ax.annotate(f"${P}$", (px, py), xytext=(-4, -8),
                textcoords="offset points", fontsize=6.5, ha="right",
                va="top", fontweight="bold")

    # r compares natural-lineage displacement with the displacement of the
    # first-cousin-null mean from the same closest sampled assignment.
    ax.plot([px, nx], [py, ny], color="#888888", lw=0.8,
            ls=(0, (3, 2)), zorder=2)
    ax.plot([px, lx], [py, ly], color=LP_COLOR, lw=2.2,
            zorder=4)
    ax.scatter([lx], [ly], marker="X", s=42, facecolor="#222222",
               edgecolor="white", lw=0.55, zorder=6)
    ax.scatter([nx], [ny], marker="D", s=35,
               facecolor=NP_COLOR, edgecolor="#222222",
               lw=0.5, zorder=6)
    ax.annotate(f"Natural lineage  {nt.plain('natural_lineage')}", (lx, ly), xytext=(6, -3),
                textcoords="offset points", fontsize=6.2, va="center")
    ax.annotate(f"Null mean  {nt.plain('null_mean')}", (nx, ny), xytext=(6, -2),
                textcoords="offset points", fontsize=6.2, va="top")

    add_colored_fraction(
        ax, prefix_x=0.15, fraction_x=0.285, center_y=0.16,
        prefix=nt.symbol("canonical_position"), numerator=rf"\mathrm{{arc}}({T}\rightarrow {P})",
        denominator=rf"\mathrm{{arc}}({T}\rightarrow {S})",
        color=U_COLOR, half_width=0.115,
    )
    add_colored_fraction(
        ax, prefix_x=0.69, fraction_x=0.825, center_y=0.87,
        prefix=nt.symbol("cousin_relative"), numerator=f"|{L}-{P}|", denominator=f"|{C}-{P}|",
        color=LP_COLOR, half_width=0.090,
    )
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    # This is a symbolic construction rather than a metric objective-space
    # panel. Letting it fill the rectangular half-width grid cell materially
    # improves legibility and avoids the unused side margins imposed by an
    # equal data aspect ratio.
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title("Canonical coordinates", loc="left", pad=3)
    add_panel_letter(ax, panel_letter, x=-0.08, y=1.12)


def plot_definition_dlp(ax, panel_letter="A"):
    """Primary symbolic definition of u and both canonical distances."""
    T, S, L, C, P = _points()
    t = np.linspace(0.0, 1.0, 15)
    x = 0.13 + 0.72 * t
    y = 0.88 - 0.72 * np.sqrt(t)
    p_idx = 6
    px, py = x[p_idx], y[p_idx]
    lx, ly = px + 0.115, py + 0.115
    nx, ny = px + 0.255, py + 0.030

    ax.annotate("", xy=(1.01, 0.0), xytext=(0.055, 0.0),
                annotation_clip=False,
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8))
    ax.annotate("", xy=(0.055, 1.00), xytext=(0.055, 0.0),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8))
    ax.text(1.01, -0.055, f"Normalized travel, {nt.math('travel_axis')}", ha="right",
            va="top", fontsize=7, clip_on=False)
    ax.text(0.015, 0.99, f"Normalized\ncell state, {nt.math('cell_state_axis')}", ha="right",
            va="top", fontsize=7)
    ax.plot(x, y, color="#A7A7A7", lw=1.3, zorder=1)
    ax.scatter(x, y, s=9, facecolor="white", edgecolor="#777777",
               lw=0.45, zorder=2)
    ax.plot(x[:p_idx + 1], y[:p_idx + 1], color=U_COLOR,
            lw=2.4, zorder=3)
    ax.plot([px, lx], [py, ly], color=LP_COLOR, lw=2.4,
            zorder=4)
    ax.plot([px, nx], [py, ny], color=NP_COLOR, lw=1.55,
            ls=(0, (2.7, 2.0)), zorder=2)
    ax.scatter([x[0]], [y[0]], s=31, marker="o",
               facecolor=TRAVEL_COLOR, edgecolor="#222222",
               lw=0.5, zorder=5)
    ax.scatter([x[-1]], [y[-1]], s=31, marker="o",
               facecolor=STATE_COLOR, edgecolor="#222222",
               lw=0.5, zorder=5)
    ax.scatter([px], [py], s=29, marker="o", facecolor="#222222",
               edgecolor="white", lw=0.5, zorder=6)
    ax.scatter([lx], [ly], marker="X", s=42, facecolor="#222222",
               edgecolor="white", lw=0.55, zorder=6)
    ax.scatter([nx], [ny], marker="D", s=42, facecolor=NP_COLOR,
               edgecolor="white", lw=0.6, zorder=7)
    ax.annotate(f"${T}=(0,1)$", (x[0], y[0]), xytext=(5, 3),
                textcoords="offset points", fontsize=6.2, ha="left",
                va="bottom")
    ax.annotate(f"${S}=(1,0)$", (x[-1], y[-1]), xytext=(4, 1),
                textcoords="offset points", fontsize=6.2, ha="left",
                va="center")
    ax.annotate(f"Closest assignment ${P}$", (px, py), xytext=(-4, -8),
                textcoords="offset points", fontsize=6.2, ha="right",
                va="top")
    ax.annotate(f"Natural lineage  {nt.plain('natural_lineage')}", (lx, ly), xytext=(6, -2),
                textcoords="offset points", fontsize=6.2, va="center")
    ax.annotate(f"First-cousin null mean, ${C}$", (nx, ny),
                xytext=(5, -6), textcoords="offset points", fontsize=6.0,
                color=NP_COLOR, ha="left", va="top")

    add_colored_fraction(
        ax, prefix_x=0.15, fraction_x=0.285, center_y=0.16,
        prefix=nt.symbol("canonical_position"), numerator=rf"\mathrm{{arc}}({T}\rightarrow {P})",
        denominator=rf"\mathrm{{arc}}({T}\rightarrow {S})",
        color=U_COLOR, half_width=0.115,
    )
    ax.text(0.65, 0.89, f"$D(a)=({nt.symbol('travel_axis')}(a),{nt.symbol('cell_state_axis')}(a))$",
            fontsize=5.8, color="#666666", ha="left", va="center")
    ax.text(0.65, 0.805,
            rf"${nt.symbol('natural_front_distance')}=\|D({L})-D({P})\|_2$",
            fontsize=7.0, color=LP_COLOR, ha="left", va="center")
    ax.text(0.65, 0.720,
            rf"${nt.symbol('null_front_distance')}=\|D({C})-D({P})\|_2$",
            fontsize=6.8, color=NP_COLOR, ha="left", va="center")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title("Endpoint-normalized coordinates", loc="left", pad=3)
    add_panel_letter(ax, panel_letter, x=-0.08, y=1.12)


def scatter_records(ax, data, u_col, distance_col, metric,
                    label_names=(), zoom=False):
    """Plot one canonical scatter using fixed metric color and region shapes."""
    ax.axhline(0.0, color="#AAAAAA", lw=0.65, zorder=1)
    if not zoom and metric == "r":
        ax.axhline(1.0, color="#AAAAAA", lw=0.7,
                   ls=(0, (3, 2)), zorder=1)
        ax.text(1.015, 1.0, "Null reference", fontsize=5.6,
                color="#777777", ha="right", va="bottom")
        ax.text(0.84, 0.012, "On the Pareto front", fontsize=5.6,
                color="#777777", ha="left", va="bottom")

    # Larger subtrees are drawn first so dense clusters do not systematically
    # hide the smaller records. AB is separated from the root star to match
    # the major-subtree panel.
    for row in data.sort_values("n", ascending=False).itertuples():
        display_region = "AB" if row.subtree == "AB" else row.region
        size = 58.0 if row.subtree == "P0" else 43.0
        ax.scatter(
            getattr(row, u_col), getattr(row, distance_col),
            marker=REGION_MARKERS[display_region], s=size,
            facecolor=LP_COLOR, edgecolor="#4B2D43", linewidth=0.55,
            alpha=0.88, zorder=4, clip_on=False,
        )

    placements = {
        "P0": (-4, 6, "right", "bottom"),
        "AB": (4, 5, "left", "bottom"),
        "ABa": (4, -3, "left", "top"),
        "ABp": (-4, 5, "right", "bottom"),
        "P1": (5, -5, "left", "top"),
        "ABala": (-6, 8, "right", "bottom"),
        "ABprppp": (6, 20, "left", "bottom"),
        "ABplpa": (-5, 9, "right", "bottom"),
    }
    for name in label_names:
        rows = data[data["subtree"] == name]
        if rows.empty:
            continue
        row = rows.iloc[0]
        if name == "ABprppp" and row[distance_col] > 0.05:
            dx, dy, ha, va = (6, -4, "left", "top")
        else:
            dx, dy, ha, va = placements.get(
                name, (4, 4, "left", "bottom"))
        ax.annotate(name, (row[u_col], row[distance_col]),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=6.1, ha=ha, va=va, zorder=6)

    if "ABplpp(p)" in label_names:
        annotate_abplpp_pair(
            ax, data, u_col, distance_col, color=LP_COLOR)

    ax.grid(True, axis="y", alpha=0.28)
    ax.set_xlabel(f"Position along front, {nt.plain('canonical_position')}")
    if metric == "dlp":
        ax.set_ylabel(f"Distance to front, {nt.math('natural_front_distance')}")
    else:
        ax.set_ylabel(f"Cousin-relative distance, {nt.math('cousin_relative')}")
    if zoom:
        ux = data[u_col].to_numpy(dtype=float)
        dy = data[distance_col].to_numpy(dtype=float)
        xspan = max(float(np.ptp(ux)), 0.05)
        yspan = max(float(np.ptp(dy)), 0.01)
        ax.set_xlim(float(ux.min() - 0.10 * xspan),
                    float(ux.max() + 0.12 * xspan))
        ax.set_ylim(max(0.0, float(dy.min() - 0.18 * yspan)),
                    float(dy.max() + 0.22 * yspan))
    else:
        ax.set_xlim(-0.025, 1.035)
        ymax = 1.06 if metric == "r" else max(0.10, data[distance_col].max() * 1.15)
        ax.set_ylim(-0.035 if metric == "r" else -0.015, ymax)
        ax.set_xticks([0.0, 0.5, 1.0])
        if metric == "r":
            ax.set_yticks([0.0, 0.5, 1.0])


def plot_major_dlp_dnp(ax, metrics, panel_letter="B"):
    """Show both canonical distances for the five major subtrees without connectors."""
    major = metrics.set_index("subtree").loc[MAJOR_SUBTREES].reset_index()
    for name in ["P1", "ABp", "ABa", "AB", "P0"]:
        row = major[major["subtree"] == name].iloc[0]
        marker = MAJOR_MARKERS[name]
        size = 51 if name == "P0" else 43
        ax.scatter(row["u_lineage_lp"], row["d_lp"], marker=marker, s=size,
                   facecolor=LP_COLOR, edgecolor="#4B2D43", linewidth=0.60,
                   alpha=0.90, zorder=4, clip_on=False)
        ax.scatter(row["u_lineage_lp"], row["d_np"], marker=marker, s=size,
                   facecolor=NP_COLOR, edgecolor="#4A3A1B", linewidth=0.60,
                   alpha=0.92, zorder=4, clip_on=False)

    upper_offsets = {
        "P0": (6, -1, "left", "bottom"),
        "AB": (-6, 5, "right", "bottom"),
        "ABa": (-6, 4, "right", "bottom"),
        "ABp": (6, 4, "left", "bottom"),
        "P1": (6, 1, "left", "bottom"),
    }
    clustered_ab = np.ptp(
        major.loc[major["subtree"].isin(["AB", "ABa", "ABp"]),
                  "u_lineage_lp"].to_numpy(dtype=float)) < 0.03
    lower_offsets = ({
        "P0": (5, 5, "left", "bottom"),
        "AB": (-8, 2, "right", "bottom"),
        "ABa": (7, 2, "left", "bottom"),
        "ABp": (-8, 11, "right", "bottom"),
        "P1": (7, 6, "left", "bottom"),
    } if clustered_ab else {
        "P0": (5, 5, "left", "bottom"),
        "AB": (-7, 3, "right", "bottom"),
        "ABa": (-5, 4, "right", "bottom"),
        "ABp": (-5, 4, "right", "bottom"),
        "P1": (7, 6, "left", "bottom"),
    })
    for row in major.itertuples():
        dx, dy, ha, va = upper_offsets[row.subtree]
        ax.annotate(row.subtree, (row.u_lineage_lp, row.d_np),
                    xytext=(dx, dy), textcoords="offset points", fontsize=5.9,
                    color=NP_COLOR, ha=ha, va=va, zorder=7)
        dx, dy, ha, va = lower_offsets[row.subtree]
        ax.annotate(row.subtree, (row.u_lineage_lp, row.d_lp),
                    xytext=(dx, dy), textcoords="offset points", fontsize=5.9,
                    color=LP_COLOR, ha=ha, va=va, zorder=7)

    # The pooled profile shifts AB, ABa, and ABp left of the legacy panel's
    # former hard-coded lower limit. Derive both ranges from the displayed
    # records so markers and their offset labels stay inside the axes instead
    # of leaking into the component crop.
    major_u = major["u_lineage_lp"].to_numpy(dtype=float)
    u_span = max(float(np.ptp(major_u)), 0.05)
    ax.set_xlim(float(major_u.min() - 0.13 * u_span),
                float(major_u.max() + 0.18 * u_span))
    major_distance_max = float(
        major[["d_lp", "d_np"]].to_numpy(dtype=float).max())
    ax.set_ylim(0.0, max(0.10, 1.17 * major_distance_max))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel(f"Position along front, {nt.math('canonical_position')}")
    ax.set_ylabel("Distance to Pareto front")
    ax.set_title("Major-subtree distances", loc="left", pad=3)
    ax.grid(True, axis="y", alpha=0.28)
    add_panel_letter(ax, panel_letter, x=-0.22, y=1.12)


def plot_all_dlp_with_major_dnp(ax, metrics, panel_letter="C"):
    """Show all natural distances plus five labeled major-subtree null-mean anchors."""
    # Draw smaller subtrees first and give larger subtrees a higher explicit
    # z-order. This keeps major lineages such as AB visible when canonical
    # positions overlap without introducing a second size encoding.
    max_n = float(metrics["n"].max())
    ordered = metrics.sort_values(["n", "subtree"], ascending=[True, True])
    for row in ordered.itertuples():
        marker = MAJOR_MARKERS["AB"] if row.subtree == "AB" else REGION_MARKERS[row.region]
        size = 58.0 if row.subtree == "P0" else 43.0
        subtree_zorder = 3.0 + 2.0 * float(row.n) / max_n
        ax.scatter(row.u_lineage_lp, row.d_lp, marker=marker, s=size,
                   facecolor=LP_COLOR, edgecolor="#4B2D43", linewidth=0.55,
                   alpha=0.86, zorder=subtree_zorder, clip_on=False)

    major = metrics.set_index("subtree").loc[MAJOR_SUBTREES].reset_index()
    null_offsets = {
        "P0": (6, -1, "left", "bottom"),
        "AB": (-7, -2, "right", "top"),
        "ABa": (6, 4, "left", "bottom"),
        "ABp": (-7, 7, "right", "bottom"),
        "P1": (6, 1, "left", "bottom"),
    }
    for row in major.itertuples():
        size = 51 if row.subtree == "P0" else 43
        ax.scatter(row.u_lineage_lp, row.d_np,
                   marker=MAJOR_MARKERS[row.subtree], s=size,
                   facecolor=NP_COLOR, edgecolor="#4A3A1B", linewidth=0.60,
                   alpha=0.92, zorder=7)
        dx, dy, ha, va = null_offsets[row.subtree]
        ax.annotate(row.subtree, (row.u_lineage_lp, row.d_np),
                    xytext=(dx, dy), textcoords="offset points", fontsize=5.8,
                    color=NP_COLOR, ha=ha, va=va, zorder=8)

    placements = {
        "ABala": (-6, 8, "right", "bottom"),
        "ABprppp": (6, 20, "left", "bottom"),
        "ABplpa": (-5, 9, "right", "bottom"),
    }
    for name, (dx, dy, ha, va) in placements.items():
        row = metrics[metrics["subtree"] == name].iloc[0]
        if name == "ABprppp" and row["d_lp"] > 0.05:
            dx, dy, ha, va = (6, -4, "left", "top")
        ax.annotate(name, (row["u_lineage_lp"], row["d_lp"]),
                    xytext=(dx, dy), textcoords="offset points", fontsize=6.1,
                    color=LP_COLOR, ha=ha, va=va, zorder=8)
    annotate_abplpp_pair(
        ax, metrics, "u_lineage_lp", "d_lp", color=LP_COLOR)

    ax.axhline(0.0, color="#AAAAAA", lw=0.65, zorder=1)
    ax.grid(True, axis="y", alpha=0.28)
    ax.set_xlim(-0.025, 1.035)
    # The null-mean distance is an independently displayed series. In the pooled run its P1
    # value exceeds every d_LP value, so a limit based only on d_LP pushes the
    # gold anchors and labels through the title. Retain the legacy d_LP margin
    # when it dominates, but reserve a larger legend band when d_NP dominates.
    dlp_max = float(metrics["d_lp"].max())
    dnp_max = float(major["d_np"].max())
    ymax = max(0.10, 1.15 * dlp_max, 1.35 * dnp_max)
    ax.set_ylim(-0.015, ymax)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xlabel(f"Position along front, {nt.math('canonical_position')}")
    ax.set_ylabel(f"Lineage-to-front distance, {nt.math('natural_front_distance')}")
    ax.set_title("All terminal-cell subtrees", loc="left", pad=5)
    add_panel_letter(ax, panel_letter, x=-0.09, y=1.11)

    shape_specs = [
        ("*", "Root"), ("P", "AB"),
        (REGION_MARKERS["ABa"], REGION_LABELS["ABa"]),
        (REGION_MARKERS["ABp"], REGION_LABELS["ABp"]),
        (REGION_MARKERS["P1"], REGION_LABELS["P1"]),
    ]
    shape_handles = [
        Line2D([0], [0], marker=marker, ls="", ms=6.0,
               markerfacecolor=LP_COLOR, markeredgecolor="#4B2D43",
               markeredgewidth=0.55, label=label)
        for marker, label in shape_specs
    ]
    shape_legend = ax.legend(
        handles=shape_handles, loc="upper center", ncol=5,
        bbox_to_anchor=(0.5, 0.98), fontsize=6.4,
        handletextpad=0.35, columnspacing=0.8, frameon=False,
    )
    ax.add_artist(shape_legend)
    null_handle = Line2D(
        [0], [0], marker="o", ls="", ms=5.3,
        markerfacecolor=NP_COLOR, markeredgecolor="#4A3A1B",
        markeredgewidth=0.55, label=f"Major-subtree null mean, {nt.math('null_front_distance')}",
    )
    ax.legend(handles=[null_handle], loc="upper right",
              bbox_to_anchor=(0.995, 0.88), frameon=False, fontsize=6.2,
              handletextpad=0.35, borderpad=0.1)


def add_metric_correlation_inset(ax, metrics):
    """Show agreement between null-relative r and endpoint d_LP."""
    inset = ax.inset_axes([0.62, 0.42, 0.34, 0.42], zorder=8)
    inset.set_facecolor("white")
    inset.patch.set_alpha(0.96)
    inset.scatter(metrics["d_lp"], metrics["relative_distance"], s=10,
                  facecolor=LP_COLOR, edgecolor="#333333", lw=0.3,
                  alpha=0.85)
    p0 = metrics[metrics["subtree"] == "P0"].iloc[0]
    inset.scatter([p0["d_lp"]], [p0["relative_distance"]], marker="*", s=32,
                  facecolor=LP_COLOR, edgecolor="#222222", lw=0.4,
                  zorder=4)
    inset.annotate("P0", (p0["d_lp"], p0["relative_distance"]),
                   xytext=(3, 2), textcoords="offset points", fontsize=4.8)
    rho = metrics[["d_lp", "relative_distance"]].corr(
        method="spearman").iloc[0, 1]
    inset.text(0.96, 0.06, rf"Spearman $\rho={rho:.2f}$",
               transform=inset.transAxes, ha="right", va="bottom",
               fontsize=5.0, color="#444444")
    inset.set_xlabel(nt.math("natural_front_distance"), fontsize=5.2, labelpad=1)
    inset.set_ylabel(nt.math("cousin_relative"), fontsize=5.2, labelpad=1)
    inset.set_title("Metric agreement", fontsize=5.8, pad=2)
    inset.tick_params(labelsize=4.7, length=2, pad=1)
    inset.grid(True, alpha=0.2, lw=0.35)
    for spine in inset.spines.values():
        spine.set_linewidth(0.5)


def save_component_crops(fig, ax_definition, ax_major, ax_all,
                         component_names, out_dir):
    """Export A/B/C from the shared canvas with matched top-row geometry.

    Panels A and B use the same vertical crop and equal half-canvas widths,
    so their axes remain aligned when LaTeX places them side by side. Panel C
    spans the complete canvas width and includes any child inset.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_inches = fig.dpi_scale_trans.inverted()
    width, height = fig.get_size_inches()
    definition_box = ax_definition.get_tightbbox(renderer).transformed(to_inches)
    major_box = ax_major.get_tightbbox(renderer).transformed(to_inches)
    top = Bbox.union([definition_box, major_box])
    top_y0 = max(0.0, top.y0 - 0.04)
    top_y1 = min(height, top.y1 + 0.04)

    # Use the same crop dimensions for A and B, but center each crop on its
    # complete tight bounding box. This includes B's y-axis label without
    # forcing either component to a different scale in LaTeX.
    top_width = max(definition_box.width, major_box.width) + 0.10

    def centered_x(box):
        x0 = box.x0 - (top_width - box.width) / 2.0
        x0 = min(max(0.0, x0), width - top_width)
        return x0, x0 + top_width

    definition_x0, definition_x1 = centered_x(definition_box)
    major_x0, major_x1 = centered_x(major_box)
    crops = {
        component_names[0]: Bbox.from_extents(
            definition_x0, top_y0, definition_x1, top_y1),
        component_names[1]: Bbox.from_extents(
            major_x0, top_y0, major_x1, top_y1),
    }

    bottom = ax_all.get_tightbbox(renderer).transformed(to_inches)
    crops[component_names[2]] = Bbox.from_extents(
        0.0, max(0.0, bottom.y0 - 0.04),
        width, min(height, bottom.y1 + 0.04),
    )

    original_visibility = {axis: axis.get_visible() for axis in fig.axes}
    for panel_index, (stem, crop) in enumerate(crops.items()):
        if stem is None:
            # The pooled schematic is exported separately as the Figure 1 amendment.
            continue
        if panel_index == 0:
            visible = {ax_definition}
        elif panel_index == 1:
            visible = {ax_major}
        else:
            # Keep panel C and its child correlation inset.
            visible = set(fig.axes) - {ax_definition, ax_major}
        for axis in fig.axes:
            axis.set_visible(axis in visible)
        fig.savefig(out_dir / f"{stem}.pdf", facecolor="white",
                    bbox_inches=crop, pad_inches=0)
        fig.savefig(out_dir / f"{stem}.png", dpi=300, facecolor="white",
                    bbox_inches=crop, pad_inches=0)
    for axis, was_visible in original_visibility.items():
        axis.set_visible(was_visible)


def plot_primary_rows(metrics, out_dir):
    """Separate the null comparison from a fully labelled subtree inventory."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def save(fig, stem):
        # Both components have the same physical width in the wrapper.
        bounds = Bbox.from_bounds(0, 0, *fig.get_size_inches())
        for extension in ("pdf", "png"):
            fig.savefig(out_dir / f"{stem}.{extension}", bbox_inches=bounds,
                        dpi=300, facecolor="white")
        plt.close(fig)

    major = metrics.set_index("subtree").loc[MAJOR_SUBTREES]
    fig = plt.figure(figsize=(7.4, 2.35))
    fig.text(.015, .96, "A", fontsize=12, weight="bold", va="top")
    fig.text(.065, .96, "Five major subtrees: natural lineage versus cousin null",
             fontsize=10.5, weight="bold", va="top")
    handles = [
        Line2D([], [], marker="o", ls="", color=LP_COLOR, ms=6,
               label=f"Natural lineage, {nt.math('natural_front_distance')}"),
        Line2D([], [], marker="o", ls="", color=NP_COLOR, ms=6,
               label=f"First-cousin null mean, {nt.math('null_front_distance')}"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.54, .85),
               ncol=2, frameon=False, fontsize=9, columnspacing=2.5,
               handletextpad=.4, borderaxespad=0)
    ax = fig.add_axes([.13, .22, .83, .48])
    rows = np.arange(len(major))
    ax.hlines(rows, major.d_lp, major.d_np, color="#BBBBBB", lw=1.3, zorder=2)
    for field, color in (("d_lp", LP_COLOR), ("d_np", NP_COLOR)):
        ax.scatter(major[field], rows, color=color, s=35, edgecolor="white",
                   linewidth=.6, zorder=3)
    ax.set_yticks(rows, ["P0 (root)", "AB", "ABa", "ABp", "P1"])
    ax.set_ylim(len(major)-.45, -.55)
    distance_max = max(.25, np.ceil(major[["d_lp", "d_np"]].to_numpy().max()*20)/20)
    ax.set_xlim(-.008*distance_max, 1.035*distance_max)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(.05))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel(f"Distance to the same closest front assignment, {nt.math('closest_point')}", fontsize=9)
    ax.tick_params(axis="y", labelsize=9, length=0, pad=7)
    ax.tick_params(axis="x", labelsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E5E5E5", lw=.6)
    save(fig, "fig5A_ce_canonical_major_subtrees")

    # Rows remove the severe overlap of the many near-zero d_LP records.
    # Root and AB occupy their own section above the ABa branch; the other
    # two columns contain the ABp and P1 branches on identical metric scales.
    root = metrics[metrics.subtree.isin(["P0", "AB"])].set_index("subtree")
    root = root.loc[[name for name in ("P0", "AB") if name in root.index]].reset_index()
    groups = {name: metrics[metrics.region == name].sort_values("subtree")
              for name in ("ABa", "ABp", "P1")}
    p1_order = ["P1", "EMS", "MS", "MSa", "MSaa", "MSap", "MSp",
                "MSpa", "MSpp", "P2", "C", "P3", "D"]
    p1_rank = {name: index for index, name in enumerate(p1_order)}
    groups["P1"] = groups["P1"].sort_values(
        "subtree", key=lambda series: series.map(p1_rank).fillna(len(p1_rank)))
    sections = [
        [("Root & AB", root), ("ABa branch", groups["ABa"])],
        [("ABp branch", groups["ABp"])],
        [("P1 branch", groups["P1"])],
    ]
    displayed = pd.concat([frame for column in sections for _, frame in column])
    if displayed.subtree.duplicated().any() or set(displayed.subtree) != set(metrics.subtree):
        raise ValueError("Figure 5 row groups must include every subtree exactly once")
    row_count = max(sum(len(frame)+2 for _, frame in column)-1 for column in sections)
    fig = plt.figure(figsize=(7.4, 4.5))
    fig.text(.015, .985, "B", fontsize=12, weight="bold", va="top")
    fig.text(.065, .985, f"All {len(metrics)} surveyed subtrees: natural lineage only",
             fontsize=10.5, weight="bold", va="top")
    dmax = max(.2, np.ceil(float(metrics.d_lp.max())*10)/10)
    for left, column in zip((.015, .345, .675), sections):
        axes = [fig.add_axes([left+.105, .115, .077, .68]),
                fig.add_axes([left+.213, .115, .087, .68])]
        for ax, xmax, ticks, title, color in (
            (axes[0], 1., [0., .5, 1.], nt.heading("canonical_position"), U_COLOR),
            (axes[1], dmax, [0., dmax/2, dmax], f"Distance, {nt.math('natural_front_distance')}", LP_COLOR),
        ):
            ax.set_ylim(row_count-.25, -1.1)
            ax.set_xlim(-.07*xmax, 1.07*xmax)
            ax.set_xticks(ticks)
            ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
            ax.tick_params(axis="x", top=True, labeltop=True, bottom=False,
                           labelbottom=False, labelsize=7.5, length=2, pad=3)
            ax.set_yticks([])
            ax.set_title(title, fontsize=8.5, color=color, pad=21)
            ax.spines[["left", "right", "bottom"]].set_visible(False)
            ax.spines["top"].set_color("#CCCCCC")
            ax.grid(axis="x", color="#E8E8E8", lw=.55)
        row_number = 0
        for heading, frame in column:
            # Place group labels across the two metric columns, with a white
            # background so they do not compete with the column gridlines.
            axes[0].text(-1.30, row_number, f"{heading} ({len(frame)})",
                         transform=axes[0].get_yaxis_transform(),
                         fontsize=8.5, weight="bold", va="center", clip_on=False,
                         bbox=dict(facecolor="white", edgecolor="none", pad=1.5),
                         zorder=5)
            for row in frame.itertuples():
                row_number += 1
                axes[0].text(-.12, row_number, row.subtree,
                             transform=axes[0].get_yaxis_transform(), ha="right",
                             va="center", fontsize=8.5, color="#222222",
                             weight="bold" if row.subtree in MAJOR_SUBTREES else "normal")
                for ax, value, color in (
                    (axes[0], row.u_lineage_lp, U_COLOR),
                    (axes[1], row.d_lp, LP_COLOR),
                ):
                    ax.hlines(row_number, 0, value, color=color, lw=.8, alpha=.45)
                    ax.scatter(value, row_number, color=color, s=19,
                               edgecolor="white", linewidth=.35, zorder=3)
            row_number += 2
    fig.text(.02, .045,
             f"${nt.symbol('canonical_position')}=0$: travel optimum; ${nt.symbol('canonical_position')}=1$: cell-state optimum.  Bold names also appear in A.",
             fontsize=8, color="#444444")
    save(fig, "fig5B_ce_canonical_all_subtrees")


def plot_summary(metrics, metric, component_names, out_dir,
                 panel_letters=("A", "B", "C")):
    """Render and split the primary d_LP or supplementary cousin-r summary."""
    if metric == "dlp":
        u_col, distance_col = "u_lineage_lp", "d_lp"
    elif metric == "r":
        u_col, distance_col = "u_lineage", "relative_distance"
    else:
        raise ValueError(f"Unknown metric: {metric}")

    fig = plt.figure(figsize=(7.15, 5.15))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.23],
                            hspace=0.42, wspace=0.34)
    ax_definition = fig.add_subplot(grid[0, 0])
    ax_major = fig.add_subplot(grid[0, 1])
    ax_all = fig.add_subplot(grid[1, :])

    if metric == "dlp":
        plot_definition_dlp(ax_definition, panel_letters[0])
        plot_major_dlp_dnp(ax_major, metrics, panel_letters[1])
        plot_all_dlp_with_major_dnp(ax_all, metrics, panel_letters[2])
    else:
        plot_definition_r(ax_definition, panel_letters[0])
        major = (metrics.set_index("subtree").loc[MAJOR_SUBTREES]
                 .reset_index())
        scatter_records(ax_major, major, u_col, distance_col, metric,
                        label_names=MAJOR_SUBTREES, zoom=True)
        ax_major.set_title("Major subtrees", loc="left", pad=5)
        add_panel_letter(ax_major, panel_letters[1], x=-0.22, y=1.13)

        outlier_labels = ["ABala", "ABprppp", "ABplpa", "ABplpp(p)"]
        scatter_records(ax_all, metrics, u_col, distance_col, metric,
                        label_names=outlier_labels, zoom=False)
        ax_all.set_title("All subtrees (n ≥ 12 terminal cells)", loc="left", pad=5)
        add_panel_letter(ax_all, panel_letters[2], x=-0.09, y=1.11)
        shape_handles = [
            Line2D([0], [0], marker=REGION_MARKERS[region], ls="", ms=6.0,
                   markerfacecolor=LP_COLOR, markeredgecolor="#4B2D43",
                   markeredgewidth=0.55, label=REGION_LABELS[region])
            for region in REGION_ORDER
        ]
        ax_all.legend(handles=shape_handles, loc="upper center", ncol=5,
                      bbox_to_anchor=(0.5, 0.98), fontsize=6.4,
                      handletextpad=0.35, columnspacing=0.8)
        add_metric_correlation_inset(ax_all, metrics)

    fig.subplots_adjust(left=0.10, right=0.96, top=0.95, bottom=0.105)

    save_component_crops(
        fig, ax_definition, ax_major, ax_all,
        component_names=component_names, out_dir=out_dir,
    )
    plt.close(fig)


FIGURE1_STEM = "fig1_endpoint_normalization_amendment"


def render_figure1_amendment(out_dir: Path) -> None:
    """Standalone endpoint/canonical schematic for co-author Figure 1 integration."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ps.configure()
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    fig.subplots_adjust(left=.16, right=.84, bottom=.18, top=.84)
    plot_definition_dlp(ax, panel_letter="")
    for label in ax.texts:
        label.set_fontsize(label.get_fontsize() * 1.35)
    ax.set_title("Endpoint coordinates and canonical metrics",
                 fontsize=13, loc="left", pad=12)
    with plt.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none"}):
        for extension in ("pdf", "png", "svg"):
            fig.savefig(out_dir / f"{FIGURE1_STEM}.{extension}", dpi=300,
                        bbox_inches="tight", facecolor="white")
    plt.close(fig)


SUPPLEMENT_COMPONENTS = ("figS1A_ce_cousin_r_definition", "figS1B_ce_cousin_r_major_subtrees",
                         "figS1C_ce_cousin_r_all_subtrees")
