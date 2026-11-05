"""Per-layer lineage panels (Track E): Figure 6B-C of the revised layout and supplement S5.

Figure 6B-C (half width each, drawn here and placed by ``publication.v2``):
one point per layer at the median embryo size when the layer's cells are
born. B shows the lineage against its front: the front distance, and the
fewest lineage edges any front assignment rewires (1 - maximum lineage
similarity, the point M of Figure 3A). That is another view of the same front,
not independent support. C shows the
natural / first-cousin shuffle cost for travel and cell state. Figure S5
(double column, a family figure): travel per edge to the own mother, an aunt
and a random mother, and which of the first two moves in each mid-layer
contrast. Layers are counted top-down (layer 1 is the top of the tree); the
input holds layers 3-8.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from publication import notation as nt, style

SIZE_LABEL = "Embryo size at birth (cells, log scale)"
RELATIVE_LABEL = "Lineage cost relative to\nfirst-cousin shuffle"
FRONT_COLOR = style.SEMANTIC_COLORS["lineage_front_distance"]
REWIRED_COLOR = style.COLORS["blue"]  # the lineage-similarity color of Figures 3 and 6
TRAVEL, STATE = style.SEMANTIC_COLORS["travel"], style.SEMANTIC_COLORS["cell_state"]
NATURAL, AUNT, RANDOM = "#222222", style.NULL_MODEL_COLORS["first_cousin"], style.NULL_MODEL_COLORS["full_random"]
SIZE_TICKS = (12, 25, 50, 100, 200, 400)
SUMMARY_SIZE = (3.6, 3.4)
OFFSET = 1.05  # travel and cell state sit either side of the layer's birth size in C


def _size_axis(ax, *, limits=(15, 600), ticks=SIZE_TICKS[1:], label=True):
    ax.set_xscale("log")
    ax.set_xlim(*limits)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) for x in ticks])
    ax.minorticks_off()
    if label:
        ax.set_xlabel(SIZE_LABEL)


def _layer_axis(ax, born, layers):
    """Top axis naming the layer at each birth size."""
    top = ax.secondary_xaxis("top")
    top.set_xticks(np.asarray(born, float))
    top.set_xticklabels([f"L{int(layer)}" for layer in layers])
    top.tick_params(length=2, labelsize=6.5, colors="#555555")
    top.minorticks_off()


def _relative_limits(t):
    """Fixed panel range unless an interval would fall outside it."""
    low = min(.5, np.floor(min(t.travel_low.min(), t.state_low.min()) * 20) / 20)
    high = max(1.08, np.ceil(max(t.travel_high.max(), t.state_high.max()) * 20) / 20 + .03)
    return low, high


def compose_position(data):
    """Figure 6B: front distance and the fewest lineage edges rewired by any assignment on the front."""
    t = data.layers
    fig, axes = plt.subplots(2, 1, figsize=SUMMARY_SIZE, sharex=True)
    fig.subplots_adjust(left=.2, right=.96, top=.9, bottom=.14, hspace=.22)
    rows = ((axes[0], "front_distance", FRONT_COLOR, f"Distance to front,\n{nt.math('natural_front_distance')}",
             max(.3, np.ceil(t.front_distance.max() * 20) / 20 + .05)),
            (axes[1], "fewest_rewired", REWIRED_COLOR, "Fewest lineage edges\nrewired on the front", 1.))
    for ax, column, color, label, top in rows:
        _size_axis(ax, label=ax is axes[-1])
        ax.plot(t.born, t[column], color=color, lw=1.2, zorder=3)
        ax.scatter(t.born, t[column], s=18, color=color, zorder=4, clip_on=False, gid=f"layer:{column}")
        ax.set_ylim(0, top)
        ax.set_ylabel(label)
    _layer_axis(axes[0], t.born, t.layer)
    return fig, dict(layers=[int(x) for x in t.layer])


def compose_cost(data):
    """Figure 6C: travel and cell-state cost relative to the first-cousin shuffle, with 95% bootstrap intervals."""
    t = data.layers
    fig, ax = plt.subplots(figsize=SUMMARY_SIZE)
    fig.subplots_adjust(left=.2, right=.96, top=.9, bottom=.14)
    _size_axis(ax)
    ax.axhline(1, color="#999999", lw=.6, zorder=1)
    ax.text(.97, 1.005, "no better than the shuffle", transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=6.5, color="#777777")
    for column, color, label, shift in (("travel", TRAVEL, "Travel", 1 / OFFSET), ("state", STATE, "Cell state", OFFSET)):
        x = t.born * shift
        ax.vlines(x, t[f"{column}_low"], t[f"{column}_high"], color=color, lw=.9, zorder=3)
        ax.plot(x, t[column], color=color, lw=1.2, zorder=3, label=label)
        ax.scatter(x, t[column], s=18, color=color, zorder=4, gid=f"layer:{column}")
    limits = _relative_limits(t)
    ax.set_ylim(*limits)
    ax.set_ylabel(RELATIVE_LABEL)
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    _layer_axis(ax, t.born, t.layer)
    return fig, dict(layers=[int(x) for x in t.layer], relative_limits=[float(x) for x in limits])


def compose_travel_decomposition(data):
    fig = plt.figure(figsize=(7.2, 2.8))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 1.1], left=.08, right=.98, top=.82, bottom=.2, wspace=.3)

    ax = fig.add_subplot(grid[0, 0])
    costs = data.per_edge[data.per_edge.tree == "pooled"]
    for column, color, name, ls in (("random", RANDOM, "random mother (full shuffle)", (0, (1, 1.5))),
                                    ("aunt", AUNT, "an aunt (first-cousin shuffle)", (0, (3, 2))),
                                    ("natural", NATURAL, "own mother (natural lineage)", "-")):
        ax.plot(costs.born, costs[column] * 1e3, color=color, ls=ls, lw=1.2, marker="o", ms=3.2, label=name,
                gid=f"per_edge:{column}")
    _size_axis(ax, limits=(18, 650))
    ax.set_ylim(0, max(4.8, costs.random.max() * 1e3 + .6))
    ax.set_ylabel("Travel per edge, distance to …\n(normalized units × 10³)")
    ax.legend(loc="center right", fontsize=6.5, bbox_to_anchor=(1., .62))
    _layer_axis(ax, costs.born, costs.layer)

    ax = fig.add_subplot(grid[0, 1])
    changes = {tree: frame.set_index("contrast") for tree, frame in data.changes.groupby("tree")}
    order = list(changes["pooled"].index)
    pos, width = np.arange(len(order)), .36
    for shift, column, color, name in ((-width / 2, "natural", NATURAL, "own mother (natural)"),
                                       (width / 2, "aunt", AUNT, "aunt (first-cousin shuffle)")):
        ax.bar(pos + shift, changes["pooled"].loc[order, column], width=width * .92, color=color,
               label=f"{name}, C. elegans full tree", gid=f"change:{column}")
        ax.scatter(pos + shift, changes["partial"].loc[order, column], s=16, marker="^", facecolor="white",
                   edgecolor=color, lw=1, zorder=3)
    ax.scatter([], [], s=16, marker="^", facecolor="white", edgecolor="#555555", lw=1, label="C. briggsae, partial forest")
    ax.axhline(0, color="#555555", lw=.7)
    ax.axvline(1.5, color="#BBBBBB", lw=.6)
    ax.set_xticks(pos)
    ax.set_xticklabels(order)
    ax.set_ylabel("Change in travel per edge,\nmid layer vs flank (%)")
    values = pd.concat([frame[["natural", "aunt"]] for frame in changes.values()]).to_numpy()
    ax.set_ylim(min(-35, np.floor(values.min() / 10) * 10 - 10), max(185, np.ceil(values.max() / 10) * 10 + 15))
    ax.legend(loc="upper right", fontsize=6.5, frameon=False)
    fig.text(.015, .9, "A", fontweight="bold", fontsize=10)
    fig.text(.5, .9, "B", fontweight="bold", fontsize=10)
    return fig, dict(contrasts=order)


def travel_decomposition_figure(data, out: Path, stem: str) -> dict:
    fig, layout = compose_travel_decomposition(data)
    style.save_figure(fig, Path(out) / f"{stem}.png")
    plt.close(fig)
    return layout
