"""Track E4 figure candidates: the per-layer lineage advantage and the mid-embryo travel dip.

Layers are numbered top-down (layer = 9 - contraction round). Layers 1-3 hold
5-18 edges and carry no claims: advantage panels omit them, and front-distance
panels draw the pooled tree's in grey. The shaded band is the
fate-based cell-sorting window of Schnabel et al. (2006), 50 to about 400 cells;
the dotted line marks the 200-cell stage, where Large et al. (2025) find
between-species expression divergence lowest.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from publication import style

CE, CB = style.SPECIES_COLORS["ce_protein"], style.SPECIES_COLORS["cb_rna"]
TRAVEL, STATE = style.SEMANTIC_COLORS["travel"], style.SEMANTIC_COLORS["cell_state"]
MUTED = "#A8A8A8"
SORTING = (50, 400)          # Schnabel et al. 2006
WAIST = 176                  # cells alive at the start of the Large et al. 200-cell stage
CLAIMED = 4                  # layers >= 4 carry claims
# series -> (label, tree, color, marker, filled)
SERIES = {
    "pooled_ce_protein": ("C. elegans, full tree (protein)", "pooled", CE, "o", True),
    "partial_ce_protein": ("C. elegans, partial forest (protein)", "partial", CE, "o", False),
    "partial_ce_rna": ("C. elegans, partial forest (RNA)", "partial", style.SPECIES_COLORS["ce_rna"], "s", False),
    "partial_cb_rna": ("C. briggsae, partial forest (RNA)", "partial", CB, "^", True),
}
EMBRYOS = {"embryo_ce_embryo1": CE, "embryo_ce_embryo2": CE, "embryo_ce_embryo3": CE,
           "embryo_cb_210519ZZY0874p1": CB, "embryo_cb_210519ZZY0874p5": CB}
ADVANTAGE = "Advantage over\nfirst-cousin shuffle"
STAGE_LABEL = "Embryo size when the layer's cells are born (cells, log scale)"


def series_table(d, series, family="first_cousin"):
    """Per-layer advantages joined with birth stages for one series (claimed layers only)."""
    tree = SERIES[series][1] if series in SERIES else "partial"
    rows = d["layers"].query("series == @series and family == @family")
    stages = d["stages"].query("tree == @tree")[["layer", "alive_median", "alive_q25", "alive_q75"]]
    return rows.merge(stages, on="layer").query("layer >= @CLAIMED").sort_values("layer")


def stage_axis(ax, *, label=True, band_label=False, full=True):
    """Log birth-stage axis; ``full`` also spans layers 1-3 (6-24 cells)."""
    ax.set_xscale("log")
    ticks = [6, 12, 25, 50, 100, 200, 400] if full else [30, 50, 100, 200, 400]
    ax.set_xlim(*((4.5, 700) if full else (27, 650)))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) for x in ticks])
    ax.minorticks_off()
    ax.axvspan(*SORTING, color="#EEEEEE", zorder=0, lw=0)
    ax.axvline(WAIST, color="#777777", lw=.7, ls=":", zorder=1)
    if band_label:
        top = ax.get_ylim()[1]
        ax.text(np.sqrt(SORTING[0] * SORTING[1]), top, "cell sorting\n(Schnabel 2006)", ha="center", va="top",
                fontsize=6.5, color="#555555")
        ax.text(WAIST * 1.06, top, "200-cell\nexpression\nminimum", ha="left", va="top", fontsize=6.5, color="#555555")
    if label:
        ax.set_xlabel(STAGE_LABEL)


def plot_series(ax, t, column, color, marker, filled, *, label=None, ms=4.5, lw=1.2, ribbon=False, offset=1.):
    """Line through the claimed layers, with intervals as whiskers or a ribbon."""
    x = t.alive_median.to_numpy(float) * offset
    y, lo, hi = (t[f"{column}_advantage"].to_numpy(), t[f"{column}_low"].to_numpy(), t[f"{column}_high"].to_numpy())
    ax.plot(x, y, color=color, lw=lw, zorder=3, label=label)
    if ribbon:
        ax.fill_between(x, lo, hi, color=color, alpha=.18, lw=0, zorder=2)
    else:
        ax.vlines(x, lo, hi, color=color, lw=.9, zorder=3)
    ax.scatter(x, y, s=ms**2, marker=marker, facecolor=color if filled else "white", edgecolor=color, lw=1, zorder=4)


def series_handles(keys, *, grey=True):
    out = []
    for k in keys:
        label, _, color, marker, filled = SERIES[k]
        out.append(Line2D([], [], color=color, marker=marker, markerfacecolor=color if filled else "white",
                          markersize=4.5, lw=1.2, label=label))
    if grey:
        out.append(Line2D([], [], color=MUTED, marker="o", markerfacecolor="white", markersize=3.6, ls="",
                          label="Layers 1-3, full tree (5-18 cells; no claims; shown in a only)"))
    return out


def context_handles():
    return [Patch(color="#EEEEEE", label="Fate-based cell sorting, 50-400 cells (Schnabel et al. 2006)"),
            Line2D([], [], color="#777777", lw=.7, ls=":", label="200-cell stage: lowest between-species "
                                                                  "expression divergence (Large et al. 2025)")]


def layer_ticks(ax, t):
    """Top axis naming the layer at each claimed birth stage."""
    top = ax.secondary_xaxis("top")
    top.set_xticks(t.alive_median.to_numpy(float))
    top.set_xticklabels([f"L{int(layer)}" for layer in t.layer])
    top.tick_params(length=2, labelsize=6.5, colors="#555555")
    top.minorticks_off()


# ---- candidates ---------------------------------------------------------------------

def stacked(d):
    """A: four rows on one stage axis: d_LP, travel advantage, cell-state advantage, Large et al. JSD."""
    keys = ("pooled_ce_protein", "partial_ce_protein", "partial_cb_rna")
    fig, axes = plt.subplots(4, 1, figsize=(4.4, 7.2), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.25, 1.25, .9]))
    fig.subplots_adjust(left=.19, right=.97, top=.86, bottom=.08, hspace=.18)
    ax = axes[0]
    for k in keys:
        _, tree, color, marker, filled = SERIES[k]
        config = k.split("_", 1)[1]
        f = d["fronts"].query("tree == @tree and config == @config and not degenerate and edges >= 5")
        st = d["stages"].query("tree == @tree")[["layer", "alive_median"]]
        f = f.merge(st, on="layer").sort_values("layer")
        claimed = f.layer >= CLAIMED
        ax.plot(f.alive_median[claimed], f.d_LP[claimed], color=color, lw=1.2, zorder=3)
        ax.scatter(f.alive_median[claimed], f.d_LP[claimed], s=20, marker=marker, facecolor=color if filled else "white",
                   edgecolor=color, lw=1, zorder=4)
        ax.scatter(f.alive_median[~claimed], f.d_LP[~claimed], s=13, marker=marker, facecolor="white", edgecolor=MUTED,
                   lw=.8, zorder=4)
    ax.set_ylabel("Distance from\nlineage to front ($d_{LP}$)")
    ax.set_ylim(bottom=0)
    stage_axis(ax, label=False)
    layer_ticks(ax, series_table(d, "pooled_ce_protein").query("layer >= @CLAIMED"))
    for ax, column, title in ((axes[1], "travel", "Travel"), (axes[2], "state", "Cell state")):
        for i, k in enumerate(keys):
            _, _, color, marker, filled = SERIES[k]
            plot_series(ax, series_table(d, k), column, color, marker, filled, offset=(.94, 1., 1.06)[i])
        ax.axhline(0, color="#999999", lw=.6, zorder=1)
        ax.set_ylabel(ADVANTAGE)
        ax.set_ylim(-.06, .32)
        ax.text(.02, .95, title, transform=ax.transAxes, ha="left", va="top", fontweight="semibold")
        stage_axis(ax, label=False)
    axes[2].set_ylim(-.06, .3)
    ax = axes[3]
    j = d["jsd"].dropna(subset=["alive_at_start"])
    ax.fill_between(j.alive_at_start, j.q25, j.q75, color="#555555", alpha=.15, lw=0)
    ax.plot(j.alive_at_start, j["median"], color="#333333", lw=1.2, marker="o", ms=3.5)
    ax.set_ylabel("Expression JSD,\nC. elegans vs C. briggsae")
    stage_axis(ax)
    ax.text(.02, .95, "Large et al. 2025, recomputed", transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
            color="#555555")
    for a, letter in zip(axes, "abcd"):
        a.text(-.25, 1.02, letter, transform=a.transAxes, fontweight="bold", fontsize=10, va="bottom")
    fig.legend(handles=series_handles(keys) + context_handles(), loc="upper center", ncol=1, fontsize=6.5,
               bbox_to_anchor=(.55, 1.0))
    return fig


def by_layer(d):
    """B: categorical layer axis; travel and cell state side by side, plus the dip contrasts."""
    keys = ("pooled_ce_protein", "partial_ce_protein", "partial_cb_rna")
    fig = plt.figure(figsize=(7.2, 2.9))
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, .85], left=.08, right=.98, top=.76, bottom=.2, wspace=.34)
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    offsets = np.linspace(-.2, .2, len(keys))
    stage = d["stages"].query("tree == 'pooled'").set_index("layer").alive_median
    layers = range(CLAIMED, 9)
    for ax, column, title in ((axes[0], "travel", "Travel"), (axes[1], "state", "Cell state")):
        ax.axvspan(5.5, 7.5, color="#EEEEEE", lw=0, zorder=0)
        for k, dx in zip(keys, offsets):
            _, _, color, marker, filled = SERIES[k]
            t = series_table(d, k)
            ax.vlines(t.layer + dx, t[f"{column}_low"], t[f"{column}_high"], color=color, lw=.9)
            ax.plot(t.layer + dx, t[f"{column}_advantage"], color=color, lw=.8, alpha=.6)
            ax.scatter(t.layer + dx, t[f"{column}_advantage"], s=18, marker=marker, facecolor=color if filled else "white",
                       edgecolor=color, lw=1, zorder=3)
        ax.axhline(0, color="#999999", lw=.6)
        ax.set_xticks(list(layers))
        ax.set_xticklabels([f"L{layer}\n{int(stage[layer])}" for layer in layers])
        ax.set_xlabel("Layer, and cells at birth (full tree)")
        ax.set_xlim(CLAIMED - .5, 8.5)
        ax.set_ylim(-.05, .32)
        ax.set_title(title, loc="left")
    axes[0].set_ylabel(ADVANTAGE)
    axes[0].text(6.5, .31, "layers 6-7", ha="center", va="top", fontsize=6.5, color="#555555")
    dip_panel(axes[2], d, keys, cost="travel")
    for a, letter in zip(axes, "abc"):
        a.text(-.2, 1.04, letter, transform=a.transAxes, fontweight="bold", fontsize=10, va="bottom")
    fig.legend(handles=series_handles(keys, grey=False), loc="upper center", ncol=3, fontsize=6.5,
               bbox_to_anchor=(.5, 1.0))
    return fig


def dip_panel(ax, d, keys, *, cost, family="first_cousin"):
    c = d["contrasts"].query("cost == @cost and family == @family")
    rows = [(k, flank) for k in keys for flank in ("late", "early")]
    y = np.array([-(i + .5 * (i // 2)) for i in range(len(rows))])   # a half-row gap between series
    for yi, (k, flank) in zip(y, rows):
        _, _, color, marker, filled = SERIES[k]
        r = c.query("series == @k and flank == @flank").iloc[0]
        ax.hlines(yi, r.low, r.high, color=color, lw=1.1)
        ax.scatter(r.dip, yi, s=20, marker=marker, facecolor=color if filled else "white", edgecolor=color, lw=1,
                   zorder=3)
    ax.axvline(0, color="#999999", lw=.6)
    ax.set_yticks(y)
    ax.set_yticklabels(["vs L8" if f == "late" else "vs L4-5" for _, f in rows], fontsize=6.5)
    ax.set_xlabel("Layers 6-7 minus flank\n(travel, 95% interval)")
    ax.set_title("Mid-layer dip", loc="left")
    ax.set_xlim(-.26, .04)


def headline(d):
    """C: one single-column panel: travel and cell-state advantage, C. elegans full tree, C. briggsae overlaid."""
    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    fig.subplots_adjust(left=.18, right=.97, top=.86, bottom=.17)
    pooled, cb = series_table(d, "pooled_ce_protein"), series_table(d, "partial_cb_rna")
    for column, color, name in (("travel", TRAVEL, "travel"), ("state", STATE, "cell state")):
        plot_series(ax, pooled, column, color, "o", True, ribbon=True)
        claimed = cb.layer >= CLAIMED
        ax.plot(cb.alive_median[claimed], cb[f"{column}_advantage"][claimed], color=color, lw=.9, ls=(0, (3, 2)),
                marker="^", ms=3.5, markerfacecolor="white", zorder=3)
        end = pooled.query("layer == 8").iloc[0]
        ax.text(end.alive_median * 1.22, end[f"{column}_advantage"], name, color=color, va="center", fontsize=7.5,
                fontweight="semibold")
    ax.axhline(0, color="#999999", lw=.6)
    ax.set_ylim(-.05, .34)
    stage_axis(ax, full=False)
    ax.set_xlim(27, 900)
    ax.set_ylabel(ADVANTAGE)
    ax.set_xlabel("Embryo size when cells are born (cells, log scale)")
    ax.text(np.sqrt(SORTING[0] * SORTING[1]), .335, "cell sorting", ha="center", va="top", fontsize=6.5, color="#555555")
    ax.text(WAIST * 1.06, .3, "200-cell\nminimum", ha="left", va="top", fontsize=6.5, color="#555555")
    handles = [Line2D([], [], color="#444444", marker="o", ms=4, lw=1.2, label="C. elegans, full tree (95% band)"),
               Line2D([], [], color="#444444", marker="^", ms=3.5, markerfacecolor="white", lw=.9, ls=(0, (3, 2)),
                      label="C. briggsae, partial forest")]
    fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=6.5, bbox_to_anchor=(.56, 1.0))
    return fig


def dip_forest(d):
    """D (supplement): every dip contrast, travel and cell state, both shuffles, every series and embryo."""
    order = list(SERIES) + list(EMBRYOS)
    names = {**{k: v[0] for k, v in SERIES.items()},
             **{k: ("C. elegans embryo " + k[-1]) if "_ce_" in k else ("C. briggsae embryo " + k[-2:]) for k in EMBRYOS}}
    colors = {**{k: v[2] for k, v in SERIES.items()}, **EMBRYOS}
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 3.6), sharey=True)
    fig.subplots_adjust(left=.27, right=.98, top=.82, bottom=.14, wspace=.12)
    panels = [(c, f) for f in ("first_cousin", "second_cousin") for c in ("travel", "state")]
    y = np.arange(len(order))[::-1]
    for ax, (cost, family) in zip(axes, panels):
        c = d["contrasts"].query("cost == @cost and family == @family")
        for yi, k in zip(y, order):
            for flank, dy, marker in (("late", .16, "o"), ("early", -.16, "s")):
                r = c.query("series == @k and flank == @flank")
                if not len(r):
                    continue
                r = r.iloc[0]
                color = colors[k] if not (cost == "travel" and k == "partial_ce_rna") else colors[k]
                ax.hlines(yi + dy, r.low, r.high, color=color, lw=1)
                ax.scatter(r.dip, yi + dy, s=14, marker=marker, facecolor=color if r.excludes_zero else "white",
                           edgecolor=color, lw=.9, zorder=3)
        ax.axvline(0, color="#999999", lw=.6)
        ax.set_title(f"{'Travel' if cost == 'travel' else 'Cell state'}\n{family.split('_')[0]} cousins",
                     fontsize=8)
        ax.set_xlim(-.4, .2)
        ax.axhline(len(order) - len(SERIES) - .5, color="#DDDDDD", lw=.6)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([names[k] for k in order], fontsize=6.5)
    fig.supxlabel("Layers 6-7 minus flank (95% interval; filled = excludes zero)", fontsize=8, y=.02)
    handles = [Line2D([], [], color="#444444", marker="o", ls="", ms=4, label="vs layer 8"),
               Line2D([], [], color="#444444", marker="s", ls="", ms=4, label="vs layers 4-5")]
    fig.legend(handles=handles, loc="upper left", ncol=2, fontsize=6.5, bbox_to_anchor=(.0, 1.0))
    return fig


def decomposition(d):
    """E (supplement): natural and shuffled travel per edge by layer; the dip comes from the natural lineage."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=False)
    fig.subplots_adjust(left=.08, right=.98, top=.8, bottom=.2, wspace=.25)
    for ax, k in zip(axes, ("pooled_ce_protein", "partial_cb_rna")):
        t = series_table(d, k)
        claimed = t.layer >= CLAIMED
        for column, color, name, ls in (("natural_travel", "#222222", "natural lineage", "-"),
                                        ("shuffle_travel", style.SEMANTIC_COLORS["first_cousin_null"],
                                         "first-cousin shuffle (mean)", (0, (3, 2)))):
            per_edge = t[column] / t.edges * 1e3
            ax.plot(t.alive_median[claimed], per_edge[claimed], color=color, ls=ls, marker="o", ms=3.5, lw=1.2,
                    label=name)
            ax.scatter(t.alive_median[~claimed], per_edge[~claimed], s=10, facecolor="white", edgecolor=MUTED, lw=.8)
        stage_axis(ax, full=False)
        ax.set_xlabel("Embryo size when cells are born (cells, log scale)")
        ax.set_title(SERIES[k][0], loc="left", fontsize=8.5)
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel("Travel per edge\n(normalized units × 10³)")
    axes[0].legend(loc="upper left", fontsize=6.5)
    return fig


CANDIDATES = (
    ("hourglass_A_stacked", stacked, "A. Four rows on one stage axis",
     "Front distance, travel advantage, cell-state advantage and the Large et al. JSD share the birth-stage axis. "
     "Column width; the most complete story in one figure."),
    ("hourglass_B_by_layer", by_layer, "B. Layer axis with the dip contrast",
     "Categorical layers with birth stage underneath; travel and cell state side by side, plus the dip intervals. "
     "Double column; closest to how the test was run."),
    ("hourglass_C_headline", headline, "C. One headline panel",
     "Travel and cell-state advantage for the C. elegans full tree with 95% bands; C. briggsae dashed. "
     "Single column; the simplest main-text option."),
    ("hourglass_D_dip_forest", dip_forest, "D. Every dip contrast (supplement)",
     "All series, single embryos, both shuffles, travel and cell state."),
    ("hourglass_E_decomposition", decomposition, "E. Where the travel dip comes from (supplement)",
     "Natural against shuffled travel per edge: mid-layer daughters travel farther, the shuffle cost changes less."),
)
