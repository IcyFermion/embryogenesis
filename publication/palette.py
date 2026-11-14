"""Colour palette of the revised manuscript layout: one colour per concept, defined once.

Every figure of ``publication.v2`` takes its colours from here, so a colour is
changed in one place and the next build picks it up everywhere. Concepts are
keyed by meaning, never by hue: a measure keeps its colour in every figure that
shows it (for example ``front_distance`` in Figures 4 and 6B).

Groups:

* ``MEASURE``: quantities and landmarks that recur across figures (lineage,
  front, objectives, lineage similarity, the canonical metrics).
* ``NULL``: the ordered shuffle family (first-, second-, third-cousin, full
  random), warm earth tones.
* ``FATE``: the categorical terminal-fate classes of the pies and cell-type
  panels. Categorical; legends always name them.
* ``SISTERS``, ``FEATURE_SET``: categories of the cell-state feature figure.
* ``NEUTRAL``: text, rules and reference lines.

The figure families (backlog figures, and the released panels of Figures 5
and S1) still draw with ``publication.style``, whose values equal the ones here
today. ``sheet()`` draws the palette for review.
"""

from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap

MEASURE = {
    "natural_lineage": "#D7191C",      # red X: the C. elegans or C. briggsae lineage
    "front": "#222222",                # Pareto front
    "travel": "#009E73",               # travel cost, travel optimum T
    "cell_state": "#D55E00",           # cell-state cost, cell-state optimum S
    "lineage_similarity": "#0072B2",   # lineage similarity, maximum-similarity point M
    "canonical_position": "#6F6F6F",   # position u
    "front_distance": "#CC79A7",       # distance d_LP from the lineage to its front
    "null_front_distance": "#B18F00",  # distance d_CP from the first-cousin null mean to the front (its cloud's hue)
}

NULL = {
    "first_cousin": "#B18F00",
    "second_cousin": "#A6611A",
    "third_cousin": "#7B4F2C",
    "full_random": "#6F3B5C",
}

# Tree distance between the natural and the assigned parent (Figure 4 arcs), light to dark.
TREE_DISTANCE_STOPS = ("#F2D49B", "#D9822B", "#A33B3B", "#4B1D4F")
TREE_DISTANCE = LinearSegmentedColormap.from_list("tree_distance", TREE_DISTANCE_STOPS)

FATE = {
    "neuron": "#0072B2",
    "muscle": "#D55E00",
    "epithelium": "#E69F00",
    "programmed_death": "#6F6F6F",
    "glial": "#56B4E9",
    "alimentary": "#009E73",
    "excretory": "#CC79A7",
    "reproduction": "#F0E442",
    "mesoderm": "#8B5A2B",
    "other": "#D9D9D9",
    "tail": "#A64D79",
}

SISTERS = {
    "same_type": "#E69F00",       # concordant terminal sisters (same fine type)
    "different_type": "#56B4E9",  # discordant terminal sisters
}

FEATURE_SET = {
    "top20": "#222222",           # the published top-20 proteins
    "whole_panel": "#8C8C8C",     # all 210 protein reporters
    "random": "#C8C8C8",          # random 20-protein sets
    "selected": "#0072B2",        # highlight of the selected features in the ranking
}

NEUTRAL = {
    "text": "#222222",
    "muted": "#555555",
    "rule": "#BBBBBB",
    "reference": "#999999",       # reference lines, e.g. ratio 1 or the null reference
    "dendrogram": "#555555",
}

GROUPS = {"MEASURE": MEASURE, "NULL": NULL, "FATE": FATE, "SISTERS": SISTERS, "FEATURE_SET": FEATURE_SET,
          "NEUTRAL": NEUTRAL}

# Hues shared between a measure and a fate class. Pies sit next to the u and d_LP/d_CP columns in the
# subtree figure. Kept for now (author, 2026-11-09); the overall theme is decided when the manuscript wraps up.
SHARED_HUES = sorted((m, f) for m, mc in MEASURE.items() for f, fc in FATE.items() if mc.lower() == fc.lower())


def sheet():
    """A figure of every group's swatches with names and hex values, for the review page.

    A fate colour that is also a measure colour carries an asterisk (``SHARED_HUES``).
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    columns = (("MEASURE", "NULL"), ("FATE",), ("SISTERS", "FEATURE_SET", "NEUTRAL"))
    shared = {f for _, f in SHARED_HUES}
    rows = max(sum(len(GROUPS[g]) + 1 for g in column) for column in columns)
    fig = plt.figure(figsize=(7.4, .2 * rows + .2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, len(columns))
    ax.set_ylim(rows - .4, -.6)
    ax.axis("off")
    for c, column in enumerate(columns):
        row = 0
        for group in column:
            ax.text(c + .03, row, group, fontsize=7.5, weight="bold", va="center")
            row += 1
            for name, color in GROUPS[group].items():
                mark = " *" if group == "FATE" and name in shared else ""
                ax.add_patch(Rectangle((c + .03, row - .36), .12, .72, facecolor=color, edgecolor="#999999", lw=.3))
                ax.text(c + .19, row, name.replace("_", " ") + mark, fontsize=7, va="center")
                ax.text(c + .7, row, color.upper(), fontsize=6.5, va="center", family="monospace", color="#555555")
                row += 1
    ax.text(1.03, rows - .9, "* same colour as a measure", fontsize=6.5, color="#555555", va="center")
    return fig
