"""Panels for the revised layout: plain fronts with direct labels, lineage similarity, rewiring, metric rows."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MPath
import numpy as np

from publication import notation as nt
from publication import palette
from publication import style
from publication import subtrees

FRONT_COLOR = palette.MEASURE["front"]
LINEAGE_RED = palette.MEASURE["natural_lineage"]
SIMILARITY_COLOR = palette.MEASURE["lineage_similarity"]
TRAVEL_COLOR = palette.MEASURE["travel"]
STATE_COLOR = palette.MEASURE["cell_state"]
U_COLOR = palette.MEASURE["canonical_position"]
D_COLOR = palette.MEASURE["front_distance"]
CP_COLOR = palette.MEASURE["null_front_distance"]
CLOUD_LABELS = {"first_cousin": "First-cousin shuffle", "second_cousin": "Second-cousin shuffle",
                "third_cousin": "Third-cousin shuffle"}
DISTANCE_CMAP = palette.TREE_DISTANCE
AXIS = ("Travel distance", "Cell-state distance")


def configure():
    style.configure()
    matplotlib.rcParams.update({"axes.grid": False})


def similarity_label():
    return nt.label("retention")


def _display(values, n=1000):
    values = np.asarray(values)
    return values[np.linspace(0, len(values) - 1, min(n, len(values)), dtype=int)]


INSET_BOX = (.68, .68, .29, .28)  # axes fraction: x, y, width, height


def _place_label(ax, text, anchor, points, color, *, fontsize=8, weight="normal", pad=.025, avoid=(), leader=None):
    """Put ``text`` beside ``anchor`` where it covers the fewest data points (data coordinates)."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    xspan, yspan = np.ptp(ax.get_xlim()), np.ptp(ax.get_ylim())
    best = None
    # ``anchor`` is a point or an extent (xmin, xmax, ymin, ymax): labels sit outside the extent.
    xmin, xmax, ymin, ymax = anchor if len(anchor) == 4 else (anchor[0], anchor[0], anchor[1], anchor[1])
    xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
    directions = ((1, 0, "left", "center"), (1, -1, "left", "top"), (1, 1, "left", "bottom"),
                  (0, -1, "center", "top"), (0, 1, "center", "bottom"), (-1, 0, "right", "center"),
                  (-1, -1, "right", "top"), (-1, 1, "right", "bottom"))
    for reach, (dx, dy, ha, va) in ((r, d) for r in (1, 2.5, 4) for d in directions):
        bx = xmax if dx > 0 else xmin if dx < 0 else xmid
        by = ymax if dy > 0 else ymin if dy < 0 else ymid
        position = (bx + dx * reach * pad * xspan, by + dy * reach * pad * yspan)
        artist = ax.text(*position, text, ha=ha, va=va, color=color, fontsize=fontsize, weight=weight)
        box = artist.get_window_extent(renderer).transformed(ax.transData.inverted())
        inside = ((points[:, 0] >= box.x0) & (points[:, 0] <= box.x1)
                  & (points[:, 1] >= box.y0) & (points[:, 1] <= box.y1)).sum()
        axes_box = artist.get_window_extent(renderer).transformed(ax.transAxes.inverted())
        outside = axes_box.x0 < 0 or axes_box.x1 > 1 or axes_box.y0 < 0 or axes_box.y1 > 1
        # Avoid boxes (x0, y0, x1, y1) in axes fractions: the inset and labels already placed.
        blocked = any(axes_box.x1 > x0 and axes_box.x0 < x1 and axes_box.y1 > y0 and axes_box.y0 < y1
                      for x0, y0, x1, y1 in avoid)
        score = inside + (10_000 if outside else 0) + (5_000 if blocked else 0) + 2 * (reach - 1)
        artist.remove()
        if best is None or score < best[0]:
            best = (score, position, ha, va, reach)
    _, position, ha, va, reach = best
    artist = ax.text(*position, text, ha=ha, va=va, color=color, fontsize=fontsize, weight=weight, zorder=9)
    box = artist.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    if reach > 1 and leader is not None:
        # Far from its points: connect the label to the nearest point it names.
        data_box = artist.get_window_extent(renderer).transformed(ax.transData.inverted())
        corners = np.array([[data_box.x0, data_box.y0], [data_box.x1, data_box.y0],
                            [data_box.x0, data_box.y1], [data_box.x1, data_box.y1],
                            [(data_box.x0 + data_box.x1) / 2, data_box.y0], [(data_box.x0 + data_box.x1) / 2, data_box.y1]])
        scale = np.array([xspan, yspan])
        gaps = np.linalg.norm((corners[:, None, :] - leader[None, :, :]) / scale, axis=2)
        i, j = np.unravel_index(np.argmin(gaps), gaps.shape)
        ax.plot([corners[i, 0], leader[j, 0]], [corners[i, 1], leader[j, 1]], color=color, lw=.6, zorder=8)
    return artist, (box.x0 - .01, box.y0 - .01, box.x1 + .01, box.y1 + .01)


def front_panel(ax, front, natural, clouds, full_random, *, lineage_label, inset=True, direct_labels=True,
                draws=1000, marker_scale=1., label_positions=None):
    """Pareto front (black), natural lineage (red X), cousin clouds with direct labels, full-random inset.

    ``front``: (x, y) in endpoint coordinates; ``clouds``: {family: (x, y)}; ``full_random``: (x, y).
    """
    fx, fy = map(np.asarray, front)
    points = [np.column_stack(front)]
    for family, (x, y) in clouds.items():
        color = palette.NULL[family]
        sx, sy = _display(x, draws), _display(y, draws)
        ax.scatter(sx, sy, s=7 * marker_scale, color=color, alpha=.22, edgecolors="none", rasterized=True, zorder=1)
        points.append(np.column_stack((sx, sy)))
    order = np.argsort(fx)
    ax.plot(fx[order], fy[order], color=FRONT_COLOR, lw=1.4, zorder=3)
    ax.scatter(fx, fy, s=6 * marker_scale, color=FRONT_COLOR, edgecolors="none", zorder=3)
    ax.scatter(*natural, marker="X", s=80 * marker_scale, color=LINEAGE_RED, edgecolor="white", lw=.6, zorder=8)
    xs = np.concatenate([p[:, 0] for p in points] + [[natural[0]]])
    ys = np.concatenate([p[:, 1] for p in points] + [[natural[1]]])
    pad_x, pad_y = .05 * np.ptp(xs), .05 * np.ptp(ys)
    x0, x1 = xs.min() - pad_x, xs.max() + pad_x
    y0, top = ys.min() - pad_y, ys.max() + pad_y
    if inset:
        # Raise the top only as far as needed to keep data under the inset (upper-right 32% x 34%).
        under = xs >= x0 + .64 * (x1 - x0)
        if under.any():
            top = max(top, y0 + (ys[under].max() + pad_y - y0) / .63)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, top)
    if inset:
        x, y = map(np.asarray, full_random)
        color = palette.NULL["full_random"]
        box = ax.inset_axes(list(INSET_BOX))
        box.scatter(x, y, s=4, color=color, alpha=.22, edgecolors="none", rasterized=True)
        box.set_title("Full-random shuffle", fontsize=7, color=color, pad=2)
        box.set_xlabel(AXIS[0], fontsize=6, labelpad=1)
        box.set_ylabel(AXIS[1], fontsize=6, labelpad=1)
        box.tick_params(labelsize=5.5, length=2, pad=1)
        box.locator_params(axis="both", nbins=3)
        for spine in ("top", "right"):
            box.spines[spine].set_visible(True)
    if direct_labels and label_positions:
        # Fixed positions {key: (x, y, ha, va)} in data coordinates; key "lineage" or a cloud family.
        for key, (x, y, ha, va) in label_positions.items():
            lineage = key == "lineage"
            ax.text(x, y, lineage_label if lineage else CLOUD_LABELS[key], ha=ha, va=va, zorder=9,
                    color=LINEAGE_RED if lineage else palette.NULL[key],
                    fontsize=8.5 if lineage else 8, weight="bold" if lineage else "normal")
    elif direct_labels:
        everything = np.vstack(points + [np.array([natural])])
        bx, by, bw, bh = INSET_BOX
        avoid = [(bx - .02, by - .1, bx + bw + .02, by + bh + .06)] if inset else []
        _, box = _place_label(ax, lineage_label, natural, everything, LINEAGE_RED, fontsize=8.5, weight="bold",
                              pad=.03, avoid=avoid)
        avoid.append(box)
        for family, (x, y) in clouds.items():
            x, y = np.asarray(x), np.asarray(y)
            anchor = (np.percentile(x, 2), np.percentile(x, 98), np.percentile(y, 2), np.percentile(y, 98))
            _, box = _place_label(ax, CLOUD_LABELS[family], anchor, everything, palette.NULL[family],
                                  fontsize=8, avoid=avoid, leader=np.column_stack((_display(x, draws), _display(y, draws))))
            avoid.append(box)
    ax.set_xlabel(AXIS[0])
    ax.set_ylabel(AXIS[1])
    return ax


def front_figure(front, natural, clouds, full_random, *, lineage_label, size=(5.6, 4.2), label_positions=None):
    fig, ax = plt.subplots(figsize=size)
    front_panel(ax, front, natural, clouds, full_random, lineage_label=lineage_label, label_positions=label_positions)
    fig.tight_layout()
    return fig


def similarity_curve(ax, x, similarity, natural_x, *, travel_index, state_index, maximum_index, annotate=True):
    """Lineage similarity along the front versus normalized travel distance (Figure 3A)."""
    x, similarity = np.asarray(x), np.asarray(similarity)
    order = np.argsort(x)
    ax.plot(x[order], similarity[order], color=SIMILARITY_COLOR, lw=1.8, zorder=2)
    ax.scatter(natural_x, 1, marker="X", s=70, color=LINEAGE_RED, edgecolor="white", lw=.5, zorder=6)
    marks = ((travel_index, "s", TRAVEL_COLOR, nt.symbol("travel_optimum"), (7, -12)),
             (state_index, "s", STATE_COLOR, nt.symbol("cell_state_optimum"), (-4, 7)),
             (maximum_index, "o", SIMILARITY_COLOR, "M", (7, 2)))
    for index, marker, color, text, offset in marks:
        ax.scatter(x[index], similarity[index], marker=marker, s=46, color=color, edgecolor="#222222", lw=.7, zorder=7)
        if annotate:
            ax.annotate(text, (x[index], similarity[index]), xytext=offset, textcoords="offset points",
                        fontsize=8.5, weight="bold", color=color)
    ax.set_ylim(0, 1.06)
    ax.set_xlabel(AXIS[0])
    ax.set_ylabel(similarity_label())
    return ax


# ---------------------------------------------------------------- rewiring view

def lineage_layout(tree_index, cells, order):
    """Leaf x (lineage order) and internal-node x/depth for the ancestry of ``cells``."""
    leaves = [c for c in order if c in set(cells)]
    x = {c: float(i) for i, c in enumerate(leaves)}
    children = {}
    for c in leaves:
        node = c
        for ancestor in tree_index[c]["ancestors"]:
            children.setdefault(ancestor, set()).add(node)
            node = ancestor
    depth = {n: tree_index[n]["depth"] for n in list(children) + leaves}

    def place(node):
        if node in x:
            return x[node]
        xs = [place(child) for child in children[node]]
        x[node] = (min(xs) + max(xs)) / 2
        return x[node]
    root = min(children, key=lambda n: depth[n])
    place(root)
    return leaves, x, depth, children, root


def draw_dendrogram(ax, x, depth, children, *, baseline, color=palette.NEUTRAL["dendrogram"], lw=.45):
    for node, kids in children.items():
        kx = [x[k] for k in kids]
        ax.plot([min(kx), max(kx)], [-depth[node]] * 2, color=color, lw=lw, solid_capstyle="butt")
        for k in kids:
            bottom = -depth[k] if k in children else baseline
            ax.plot([x[k], x[k]], [-depth[node], bottom], color=color, lw=lw, solid_capstyle="butt")


def rewiring_figure(t, points, *, founders=("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "D")):
    """Natural lineage on top; one row per front point with arcs from each moved terminal to its new parent."""
    from publication.v2.data import edge_changes
    from terminal_pareto import data_loader as dl
    order = subtrees.terminal_inorder(t.lineage, dl.map_names)
    leaves, x, depth, children, root = lineage_layout(t.tree_index, t.nodes, order)
    baseline = -max(depth.values()) - .4
    rows = len(points)
    fig = plt.figure(figsize=(7.4, 1.75 + 1.05 * rows))
    top = fig.add_axes([.12, 1 - 1.45 / fig.get_figheight(), .82, 1.25 / fig.get_figheight()])
    draw_dendrogram(top, x, depth, children, baseline=baseline)
    for name in founders:
        if name in x:
            top.text(x[name], baseline - .6, name, ha="center", va="top", fontsize=7, color="#333333")
    top.set_xlim(-2, len(leaves) + 1)
    top.set_ylim(baseline - 2.2, -depth[root] + .5)
    top.axis("off")
    top.text(-.01, .55, "Natural\nlineage", transform=top.transAxes, ha="right", va="center", fontsize=8.5, weight="bold")
    slot_x = {}
    for c in leaves:
        slot_x.setdefault(t.parents[t.nodes.index(c)], []).append(x[c])
    norm = Normalize(1, 18)
    height = (1.45 + .15) / fig.get_figheight()
    row_h = .95 / fig.get_figheight()
    axes = []
    for i, (letter, title, weight) in enumerate(points):
        ax = fig.add_axes([.12, 1 - height - (i + 1) * row_h - .06 / fig.get_figheight() * i, .82, row_h * .92])
        axes.append(ax)
        changes = edge_changes(t, weight)
        ax.axhline(0, color="#BBBBBB", lw=.6, zorder=0)
        for row in changes.sort_values("tree_distance").itertuples():
            x0 = x[row.child]
            x1 = float(np.mean(slot_x[row.assigned_parent]))
            h = .15 + .85 * min(abs(x1 - x0) / len(leaves) * 2.2, 1)
            path = MPath([(x0, 0), ((x0 + x1) / 2, -2 * h), (x1, 0)], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor=DISTANCE_CMAP(norm(row.tree_distance)),
                                   lw=.75, alpha=.9))
        moved = [x[c] for c in changes.child]
        ax.scatter(moved, np.zeros(len(moved)), s=4, color="#222222", zorder=3, linewidths=0)
        ax.set_xlim(-2, len(leaves) + 1)
        ax.set_ylim(-1.05, .08)
        ax.axis("off")
        similarity = 1 - len(changes) / len(leaves)
        ax.text(-.01, .55, f"{title}\n{similarity_label().lower()} {similarity:.2f}", transform=ax.transAxes,
                ha="right", va="center", fontsize=8)
        ax.text(-.135, 1.0, letter, transform=ax.transAxes, fontsize=11, weight="bold", va="top")
    cax = fig.add_axes([.955, .08, .012, .55])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=DISTANCE_CMAP), cax=cax)
    bar.set_label("Tree distance between natural\nand assigned parent (divisions)", fontsize=7.5)
    bar.ax.tick_params(labelsize=7)
    return fig, axes


# ---------------------------------------------------------------- metric rows

def _pie(ax, center, radius, composition, colors):
    total = sum(composition.values())
    angle = 90.
    for name, count in composition.most_common():
        sweep = 360. * count / total
        ax.add_patch(Wedge(center, radius, angle - sweep, angle, facecolor=colors.get(name, "#D9D9D9"),
                           edgecolor="white", lw=.25))
        angle -= sweep


@dataclass(frozen=True)
class MetricColumn:
    """One metric column of the subtree rows: a heading, a key line, the ``series`` ((field, color), ...) drawn
    on one shared axis, each on a stem from zero, the axis maximum and the column's width in base widths."""
    heading: str
    heading_color: str
    series: tuple
    xmax: float
    width: float = 1.

    def keys(self):
        """The key line: the symbol of each series in its colour, with a dot when the column shows several."""
        dot = "\u25CF " if len(self.series) > 1 else ""
        return [(dot + nt.math(SERIES_SYMBOLS[field]), color) for field, color in self.series]


SERIES_SYMBOLS = {"u_lineage_lp": "canonical_position", "d_lp": "natural_front_distance",
                  "d_np": "null_front_distance"}


def subtree_metrics(canonical):
    """Figure 4 columns: position u, and the distances to the front of the lineage (d_LP) and of the first-cousin
    null mean (d_CP) on one axis (author, 2026-11-10: replaces r, which the two distances let the reader judge)."""
    dmax = max(.5, np.ceil(float(canonical[["d_lp", "d_np"]].to_numpy().max()) * 2) / 2)  # endpoint units: 0-1
    return [MetricColumn(nt.label("canonical_position"), U_COLOR, (("u_lineage_lp", U_COLOR),), 1.),
            MetricColumn("Distance to front", palette.NEUTRAL["text"], (("d_lp", D_COLOR), ("d_np", CP_COLOR)), dmax,
                         width=2.33)]


def mark_value(ax, row, value, mark, xmax, *, fontsize=8):
    """Note mark just right of a plotted readout."""
    ax.text(value + .06 * xmax, row, mark, fontsize=fontsize, color=palette.NEUTRAL["text"], va="center", ha="left",
            clip_on=False, zorder=4)


def metric_columns(fig, left, bottom, height, sections, metrics, *, row_count, pies=None, marks=None, pairs=(),
                   label_width=.1, pie_width=.035, column_width=.085, gap=.028):
    """Subtree rows: [labels | pie | one column per metric] for grouped ``sections``.

    ``sections``: [(heading, [(label, {field: value}, key), ...]), ...]; ``metrics``: ``MetricColumn``s;
    ``pies``: {key: Counter} or None; ``marks``: {key: (note mark, field)} as ``subtrees.NOTES``;
    ``pairs``: (key, key) left-right counterparts, joined by a bracket left of their labels when adjacent.
    Returns {field: axes}; every field of a column maps to that column's axes.
    """
    x = left + label_width
    pie_ax = None
    if pies is not None:
        pie_ax = fig.add_axes([x, bottom, pie_width, height])
        x += pie_width + .008
    axes, columns = {}, []
    for column in metrics:
        width = column_width * column.width
        ax = fig.add_axes([x, bottom, width, height])
        x += width + gap
        columns.append((ax, column))
        for field, _ in column.series:
            axes[field] = ax
        xmax = column.xmax
        ax.set_ylim(row_count - .3, -1.2)
        ax.set_xlim(-.07 * xmax / column.width, xmax * (1 + .07 / column.width))
        ax.set_xticks([0, xmax / 2, xmax])
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, labelsize=7.5,
                       length=2, pad=2)
        ax.set_yticks([])
        ax.spines[["left", "right", "bottom"]].set_visible(False)
        ax.spines["top"].set_color(palette.NEUTRAL["rule"])
        ax.annotate(column.heading, (.5, 1), xycoords="axes fraction", xytext=(0, 27), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8.3, weight="bold", color=column.heading_color)
        keys = column.keys()
        for i, (text, color) in enumerate(keys):
            ax.annotate(text, ((i + .5) / len(keys), 1), xycoords="axes fraction", xytext=(0, 16),
                        textcoords="offset points", ha="center", va="bottom", fontsize=8.3, color=color)
    first = columns[0][0]
    if pie_ax is not None:
        fig_w, fig_h = fig.get_size_inches()
        inch_per_row = height * fig_h / (row_count - .3 + 1.2)
        width_rows = pie_width * fig_w / inch_per_row
        pie_ax.set_xlim(0, width_rows)
        pie_ax.set_ylim(row_count - .3, -1.2)
        pie_ax.axis("off")
    extra = (pie_width + .008) if pie_ax is not None else 0
    label_x = -(extra + .006) / column_width
    heading_x = -(label_width + extra) / column_width
    row, rows = 0, {}
    for heading, entries in sections:
        if heading:
            first.text(heading_x, row, heading, transform=first.get_yaxis_transform(), fontsize=8.5, weight="bold",
                       va="center", clip_on=False)
        for label, values, key in entries:
            row += 1
            rows[key] = row
            first.text(label_x, row, label, transform=first.get_yaxis_transform(), ha="right", va="center",
                       fontsize=8.3, clip_on=False)
            for ax, column in columns:
                # Longest stem first, so each shorter stem stays visible on top of it.
                for field, color in sorted(column.series, key=lambda fc: -values[fc[0]]):
                    ax.hlines(row, 0, values[field], color=color, lw=.8, alpha=.45)
                    ax.scatter(values[field], row, color=color, s=18, edgecolor="white", linewidth=.35, zorder=3)
            if pie_ax is not None and key in pies:
                _pie(pie_ax, (pie_ax.get_xlim()[1] / 2, row), .42, pies[key], PIE_COLORS)
            if key in (marks or {}):
                mark, field = marks[key]
                column = next(c for _, c in columns if field in dict(c.series))  # same offset in inches as width 1
                mark_value(axes[field], row, values[field], mark, column.xmax / column.width)
        row += 2
    bracket_x = heading_x + .3  # axes-fraction units of the first metric column: near the group's left edge
    for a, b in pairs:
        if a in rows and b in rows and abs(rows[a] - rows[b]) == 1:
            y0, y1 = sorted((rows[a], rows[b]))
            first.plot([bracket_x + .08, bracket_x, bracket_x, bracket_x + .08], [y0 - .15, y0 - .15, y1 + .15, y1 + .15],
                       transform=first.get_yaxis_transform(), color=palette.NEUTRAL["rule"], lw=.8, clip_on=False)
    return axes


PIE_COLORS = palette.FATE

# Left-right counterparts among the displayed subtrees (Track D pairing; ABalp/ABara pair by fate).
SUBTREE_PAIRS = (("ABpl", "ABpr"), ("ABplp", "ABprp"), ("ABalp", "ABara"), ("MSa", "MSp"))
# Display order: counterparts adjacent (author, 2026-11-08); ABal/ABar are sisters, kept together above their pair.
SUBTREE_ORDER = {"Root & AB branch": ("P0", "AB"),
                 "ABa branch": ("ABa", "ABal", "ABar", "ABalp", "ABara"),
                 "ABp branch": ("ABp", "ABpl", "ABpr", "ABplp", "ABprp"),
                 "P1 branch": ("P1", "EMS", "MS", "MSa", "MSp", "P2")}


def fate_legend(fig, types, *, y=.02, x=.5, ncol=6, loc="lower center", columnspacing=1.1):
    handles = [Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor=PIE_COLORS[t], markeredgecolor="white",
                      label=t.replace("_", " ")) for t in types]
    return fig.legend(handles=handles, loc=loc, bbox_to_anchor=(x, y), ncol=min(len(types), ncol), fontsize=7.5,
                      frameon=False, handletextpad=.2, columnspacing=columnspacing, title="Terminal fates",
                      title_fontsize=7.5, alignment="left")


def subtree_rows(canonical, nodes, *, notes=None, width=7.4, row_inch=.241):
    """Figure 4: every displayed subtree (at least 24 cells, Track D6) with fate pie, position u and the distances
    to the front of the lineage (d_LP) and of the first-cousin null mean (d_CP), on one axis; left-right
    counterparts adjacent and bracketed.

    Two columns of groups: Root & AB and ABa on the left (the fate legend fills the space below, within the
    left half), ABp and P1 on the right. With ``notes`` (the Track D evidence, checked by
    ``subtrees.check_notes``), ABplp carries the dagger of released Figure 5B.
    """
    table = canonical.set_index("subtree")
    order = [name for names in SUBTREE_ORDER.values() for name in names]
    if sorted(order) != sorted(table.index):
        raise ValueError(f"Subtree order and data differ: {sorted(set(order) ^ set(table.index))}")
    metrics = subtree_metrics(canonical)
    fields = [field for column in metrics for field, _ in column.series]

    def section(heading):
        return heading, [(n, {f: float(table.loc[n, f]) for f in fields}, n) for n in SUBTREE_ORDER[heading]]
    columns = [[section("Root & AB branch"), section("ABa branch")],
               [section("ABp branch"), section("P1 branch")]]
    row_count = max(sum(len(es) + 2 for _, es in column) - 1 for column in columns)
    pies = {name: nodes[name]["composition"] for name in order if name in nodes}
    marks = None
    if notes:
        subtrees.check_notes(notes)
        marks = subtrees.NOTES
    top, bottom = .72, .2  # inches: column heads above, footnote below
    rows_height = row_inch * (row_count + .9)
    height = rows_height + top + bottom
    fig = plt.figure(figsize=(width, height))
    halves = [metric_columns(fig, left, bottom / height, rows_height / height, column, metrics, row_count=row_count,
                             pies=pies, marks=marks, pairs=SUBTREE_PAIRS) for left, column in zip((.0, .5), columns)]
    types = sorted({t for c in pies.values() for t in c}, key=lambda t: -sum(c[t] for c in pies.values()))
    left_rows = sum(len(es) + 2 for _, es in columns[0]) - 1
    # Rows map to the axes as (row + 1.2) / (row_count + .9) from the top; the legend starts below the left column.
    legend_top = (bottom + rows_height * (row_count - .3 - (left_rows + .2)) / (row_count + .9)) / height
    legend = fate_legend(fig, types, x=.012, y=legend_top, ncol=3, loc="upper left", columnspacing=.9)
    # The legend must stay inside the left half: no wider than its rows (labels to the last metric column).
    fig.canvas.draw()
    box = legend.get_window_extent().transformed(fig.transFigure.inverted())
    right = max(ax.get_position().x1 for ax in halves[0].values())
    if box.x1 > right:
        raise ValueError(f"Fate legend reaches x = {box.x1:.3f}, past the left half's rows ({right:.3f})")
    if marks:
        fig.text(.5, .02 / height, "\u2020 refers to the caption.", fontsize=7.5, color="#444444", ha="center")
    return fig


def layer_grid(by_layer, *, columns=4):
    """Figure 6A: each layer's front with cousin clouds and full-random inset (no similarity coloring)."""
    keys = sorted(by_layer)
    rows = int(np.ceil(len(keys) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(7.4, 1.95 * rows + .55), squeeze=False)
    fig.subplots_adjust(left=.08, right=.99, bottom=.15, top=.94, hspace=.48, wspace=.28)
    for ax, layer in zip(axes.flat, keys):
        d = by_layer[layer]
        curve, metric = d.select("raw3d", "ce_protein")
        cloud = d.cloud("raw3d", "ce_protein")
        clouds = {f: (cloud[cloud.family == f].D1.to_numpy(), cloud[cloud.family == f].D2.to_numpy())
                  for f in CLOUD_LABELS}
        full = cloud[cloud.family == "full_random"]
        front_panel(ax, (curve.D1.to_numpy(), curve.D2.to_numpy()), (metric.natural_D1, metric.natural_D2),
                    clouds, (full.D1.to_numpy(), full.D2.to_numpy()), lineage_label="", direct_labels=False,
                    marker_scale=.55, draws=600)
        for child in ax.child_axes:
            child.set_title("")
            child.set_xlabel("")
            child.set_ylabel("")
            child.tick_params(labelsize=4.5)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(labelsize=6.5)
        ax.set_title(f"Layer {layer} · {d.edges} cells", fontsize=8.5)
    for ax in axes.flat[len(keys):]:
        ax.set_visible(False)
    fig.supxlabel(AXIS[0], fontsize=9, y=.075)
    fig.supylabel(AXIS[1], fontsize=9, x=.01)
    handles = [Line2D([], [], marker="o", ls="", color=palette.NULL[f], markersize=5, label=CLOUD_LABELS[f])
               for f in CLOUD_LABELS]
    handles += [Line2D([], [], marker="o", ls="", color=palette.NULL["full_random"], markersize=5,
                       label="Full-random shuffle (insets)"),
                Line2D([], [], color=FRONT_COLOR, lw=1.4, label="Pareto front"),
                Line2D([], [], marker="X", ls="", color=LINEAGE_RED, markersize=7, label="C. elegans lineage")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7, frameon=False, bbox_to_anchor=(.5, 0),
               handletextpad=.3, columnspacing=1.)
    return fig


# ---------------------------------------------------------------- per-layer summaries (Figures 7B-C, 8C)

# Axes of the per-layer panels, as in the full-tree-layers family (copied 2026-11-13 so the set stands alone).
SIZE_LABEL = "Embryo size at birth (cells, log scale)"
RELATIVE_LABEL = "Lineage cost relative to\nfirst-cousin shuffle"
SIZE_TICKS = (12, 25, 50, 100, 200, 400)
SUMMARY_SIZE = (3.6, 3.4)
OFFSET = 1.05  # travel and cell state sit either side of the layer's birth size


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
    top.tick_params(length=2, labelsize=6.5, colors=palette.NEUTRAL["muted"])
    top.minorticks_off()


def _relative_limits(t):
    """Fixed panel range unless an interval would fall outside it."""
    low = min(.5, np.floor(min(t.travel_low.min(), t.state_low.min()) * 20) / 20)
    high = max(1.08, np.ceil(max(t.travel_high.max(), t.state_high.max()) * 20) / 20 + .03)
    return low, high


def layer_position(table, *, dmax=None, size=None):
    """Distance from the lineage to its front and the lowest % of edges differing from it, one point per layer.

    ``table``: layer, born, front_distance, fewest_rewired. Figure 7B (C. elegans) and 8C (C. briggsae) pass
    the same ``dmax`` so their distance axes compare directly.
    """
    t = table
    fig, axes = plt.subplots(2, 1, figsize=size or SUMMARY_SIZE, sharex=True)
    fig.subplots_adjust(left=.2, right=.96, top=.9, bottom=.14, hspace=.22)
    dmax = dmax or max(.3, np.ceil(t.front_distance.max() * 20) / 20 + .05)
    rows = ((axes[0], "front_distance", 1, D_COLOR,
             f"Distance from lineage\nto front, {nt.math('natural_front_distance')}", dmax),
            (axes[1], "fewest_rewired", 100, SIMILARITY_COLOR, "Lowest % of edges\ndiffering, front\nvs lineage", 100))
    for ax, column, scale, color, label, top in rows:
        _size_axis(ax, label=ax is axes[-1])
        ax.plot(t.born, t[column] * scale, color=color, lw=1.2, zorder=3)
        ax.scatter(t.born, t[column] * scale, s=18, color=color, zorder=4, clip_on=False)
        ax.set_ylim(0, top)
        ax.set_ylabel(label)
    _layer_axis(axes[0], t.born, t.layer)
    return fig


def layer_cost(table):
    """Figure 7C: travel and cell-state cost relative to the first-cousin shuffle, with 95% bootstrap intervals."""
    t = table
    fig, ax = plt.subplots(figsize=SUMMARY_SIZE)
    fig.subplots_adjust(left=.2, right=.96, top=.9, bottom=.14)
    _size_axis(ax)
    ax.axhline(1, color=palette.NEUTRAL["reference"], lw=.6, zorder=1)
    ax.text(.97, 1.005, "no better than the shuffle", transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=6.5, color="#777777")
    for column, color, label, shift in (("travel", TRAVEL_COLOR, "Travel", 1 / OFFSET),
                                        ("state", STATE_COLOR, "Cell state", OFFSET)):
        x = t.born * shift
        ax.vlines(x, t[f"{column}_low"], t[f"{column}_high"], color=color, lw=.9, zorder=3)
        ax.plot(x, t[column], color=color, lw=1.2, zorder=3, label=label)
        ax.scatter(x, t[column], s=18, color=color, zorder=4)
    ax.set_ylim(*_relative_limits(t))
    ax.set_ylabel(RELATIVE_LABEL)
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    _layer_axis(ax, t.born, t.layer)
    return fig


def similarity_figure(x, similarity, natural_x, *, travel_index, state_index, maximum_index, size=(3.6, 2.7)):
    """Figure 4A / 8B: lineage similarity along a front against travel distance."""
    fig, ax = plt.subplots(figsize=size)
    similarity_curve(ax, x, similarity, natural_x, travel_index=travel_index, state_index=state_index,
                     maximum_index=maximum_index)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- cell-state features (Figure 2)

def _reporter_labels(proteins):
    """Gene names; the strain is added only where one gene has several reporters in the list."""
    genes = [p.split("_")[0] for p in proteins]
    return [f"{g} ({p.split('_', 1)[1]})" if genes.count(g) > 1 else g for g, p in zip(genes, proteins)]


def _panel_letter(ax, letter, x=-.16, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11, weight="bold", va="bottom", ha="left")


def feature_figure(fs, *, size=(7.4, 6.3)):
    """Figure 2: how the top-20 reporters are chosen and why they suit the cell-state cost.

    A: classifier importance of the top-20 (bars coloured by the fate class of the largest coefficient) and of all
    210 reporters (inset). B: cosine distance of terminal sisters of the same and of different fine type. C: sibling
    separation (Cohen's d, 95% bootstrap intervals) against random 20-reporter sets. D: mean cosine distance of all
    cohort pairs by lineage path length. (A panel E of same-type pairs against the within-type null was dropped by
    the author on 2026-11-09: the measure stays in Methods, the picture was not convincing.)
    """
    fig = plt.figure(figsize=size)
    grid = fig.add_gridspec(2, 6, height_ratios=(1.25, 1), hspace=.62, wspace=1.6, left=.1, right=.98, top=.95,
                            bottom=.09)

    # A: selection.
    ax = fig.add_subplot(grid[0, :3])
    top = fs["ranking"].head(20)
    y = np.arange(len(top))
    ax.barh(y, top.global_importance, color=[palette.FATE.get(c, "#D9D9D9") for c in top.top_class], height=.72)
    ax.set_yticks(y, _reporter_labels(top.protein.tolist()), fontsize=6.6)
    ax.set_ylim(len(top) - .4, -.6)
    ax.set_xlabel("Classifier importance (largest coefficient)", fontsize=7.5)
    ax.set_xlim(0, 1.5)
    ax.set_xticks([0, .2, .4, .6, .8])
    ax.spines["bottom"].set_bounds(0, .9)
    ax.tick_params(axis="y", length=0, pad=2)
    classes = list(dict.fromkeys(top.top_class))
    handles = [Line2D([], [], marker="s", ls="", markersize=5.5, color=palette.FATE[c], label=c.replace("_", " "))
               for c in classes]
    ax.legend(handles=handles, loc="lower right", fontsize=6.3, title="Fate class of the\nlargest coefficient",
              title_fontsize=6.3, handletextpad=.2, borderaxespad=0, alignment="left")
    inset = ax.inset_axes([.68, .66, .3, .3])
    ranks = fs["ranking"]
    inset.plot(ranks["rank"], ranks.global_importance, color=palette.FEATURE_SET["whole_panel"], lw=1)
    inset.plot(ranks["rank"][:20], ranks.global_importance[:20], color=palette.FEATURE_SET["selected"], lw=1.6)
    inset.axvline(20.5, color=palette.NEUTRAL["reference"], lw=.6, ls=(0, (2, 2)))
    inset.set_xlabel(f"Rank of {len(ranks)} reporters", fontsize=6, labelpad=1)
    inset.set_ylabel("Importance", fontsize=6, labelpad=1)
    inset.tick_params(labelsize=5.5, length=2, pad=1)
    inset.text(24, inset.get_ylim()[1] * .92, "top 20", fontsize=6, color=palette.FEATURE_SET["selected"], va="top")
    _panel_letter(ax, "A", x=-.3)

    # B: sister distances, top-20 against the whole panel.
    ax = fig.add_subplot(grid[0, 3:])
    sisters = fs["sisters"]
    sep = fs["separation"].set_index("candidate")
    names = {"top-20": "Top-20", "whole panel": f"All {len(ranks)}"}
    rng = np.random.default_rng(0)
    for k, candidate in enumerate(names):
        for side, kind in ((-.19, "same"), (.19, "different")):
            values = sisters[(sisters.candidate == candidate) & (sisters.kind == kind)].distance.to_numpy()
            x = k + side
            color = palette.SISTERS[f"{kind}_type"]
            ax.boxplot(values, positions=[x], widths=.28, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor="none", edgecolor=color, lw=.9), medianprops=dict(color=color, lw=1.3),
                       whiskerprops=dict(color=color, lw=.8), capprops=dict(color=color, lw=.8))
            ax.scatter(x + rng.uniform(-.09, .09, len(values)), values, s=5, color=color, alpha=.55, lw=0, zorder=3)
        row = sep.loc[candidate]
        ax.text(k, 1.02, f"d = {row.d:.2f}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(names)), list(names.values()))
    ax.set_xlim(-.6, len(names) - .4)
    ax.set_ylabel("Cosine distance between sisters")
    handles = [Line2D([], [], marker="o", ls="", markersize=5, color=palette.SISTERS["same_type"],
                      label=f"same type ({int(sep.loc['top-20'].n_same)} pairs)"),
               Line2D([], [], marker="o", ls="", markersize=5, color=palette.SISTERS["different_type"],
                      label=f"different type ({int(sep.loc['top-20'].n_different)} pairs)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.12), ncol=2, fontsize=6.6,
              title="Terminal sisters", title_fontsize=6.6, handletextpad=.1, columnspacing=.8)
    _panel_letter(ax, "B", x=-.22)

    # C: sibling separation against random sets.
    ax = fig.add_subplot(grid[1, :3])
    order = [("top-20", "Top-20"), ("HV pool", "HV pool"), ("whole panel", f"All {len(ranks)}")]
    random = fs["random"]
    parts = ax.violinplot([random], positions=[len(order)], widths=.7, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(palette.FEATURE_SET["random"]); body.set_edgecolor("none"); body.set_alpha(.9)
    ax.scatter([len(order)], [np.median(random)], s=14, color=palette.FEATURE_SET["whole_panel"], zorder=3)
    for k, (candidate, _) in enumerate(order):
        row = sep.loc[candidate]
        color = palette.FEATURE_SET["top20"] if candidate == "top-20" else palette.FEATURE_SET["whole_panel"]
        ax.vlines(k, row.low, row.high, color=color, lw=1)
        ax.scatter([k], [row.d], s=22, color=color, zorder=3)
    ax.axhline(0, color=palette.NEUTRAL["reference"], lw=.6)
    ax.set_xticks(range(len(order) + 1), [label for _, label in order] + ["Random 20\n(HV pool)"], fontsize=6.6)
    ax.set_ylabel("Sibling separation, d")
    _panel_letter(ax, "C", x=-.22)

    # D: all pairs by lineage path length.
    ax = fig.add_subplot(grid[1, 3:])
    p = fs["by_path"]
    p = p.assign(x=p.path_low + .5)
    ax.fill_between(p.x, p.random_low, p.random_high, color=palette.FEATURE_SET["random"], lw=0,
                    label="random 20 (5–95%)")
    ax.plot(p.x, p["whole panel"], color=palette.FEATURE_SET["whole_panel"], lw=1.2, label=f"all {len(ranks)}")
    ax.plot(p.x, p["top-20"], color=palette.FEATURE_SET["top20"], lw=1.4, marker="o", ms=2.5, label="top-20")
    ax.set_xlabel("Lineage path length (divisions)")
    ax.set_ylabel("Mean cosine distance,\nall cell pairs")
    ax.legend(loc="lower right", fontsize=6.3, handlelength=1.2, borderaxespad=.1)
    ax.set_xlim(0, 22)
    ax.set_xticks([0, 5, 10, 15, 20])
    _panel_letter(ax, "D", x=-.22)

    return fig


# ---------------------------------------------------------------- cell types (Figure 5)

CELL_TYPES = ("programmed_death", "epithelium", "neuron", "muscle", "alimentary", "glial", "other (n<=4)")
WITHIN_TYPES = CELL_TYPES[:6]  # types of at least 15 cells, with their own front
RESTRICTED_COLOR = palette.MEASURE["restricted_front"]
SIMILARITY_BINS, SIMILARITY_BANDWIDTH = 120, .018  # display-only resampling of the heat map (as released Figure 6A)


def _fate(cell_type):
    return "other" if cell_type.startswith("other") else cell_type


def smooth_profiles(u, values, n_bins=SIMILARITY_BINS, bandwidth=SIMILARITY_BANDWIDTH):
    """Resample each row's profile on a uniform u grid and average with a narrow Gaussian (display only).

    Linear resampling keeps densely sampled stretches of the sweep from gaining visual weight; edge padding avoids
    an artificial decay at the ends. Returns the bin edges and the profiles clipped to [0, 1]."""
    u, values = np.asarray(u, float), np.asarray(values, float)
    order = np.argsort(u)
    u, values = u[order], values[:, order]
    edges = np.linspace(0, 1, n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    resampled = np.vstack([np.interp(centers, u, row) for row in values])
    sigma = bandwidth * n_bins
    radius = max(1, int(np.ceil(3 * sigma)))
    kernel = np.exp(-.5 * (np.arange(-radius, radius + 1) / sigma) ** 2)
    kernel /= kernel.sum()
    smoothed = np.vstack([np.convolve(np.pad(row, radius, mode="edge"), kernel, mode="valid") for row in resampled])
    return edges, np.clip(smoothed, 0, 1)


def _unique_front(frame):
    """Drop consecutive duplicate points of a sweep."""
    xy = frame[["travel_sigma_per_cell", "state_sigma_per_cell"]].to_numpy()
    keep = np.r_[True, np.any(np.abs(np.diff(xy, axis=0)) > 1e-12, axis=1)]
    return frame.loc[keep]


def _optima(ax, travel, state, *, edge="#222222", size=30):
    """Travel optimum T and cell-state optimum S as squares in their objective's colour (as in Figure 3)."""
    ax.scatter(*travel, marker="s", s=size, color=TRAVEL_COLOR, edgecolor=edge, lw=.6, zorder=6)
    ax.scatter(*state, marker="s", s=size, color=STATE_COLOR, edgecolor=edge, lw=.6, zorder=6)


def similarity_heatmap(ax, cax, retention):
    """Figure 5A: lineage similarity of each fate class along the unrestricted front (rows), with M marked."""
    types = [c for c in CELL_TYPES if c in set(retention.type)]
    fronts = retention.drop_duplicates("front_index").sort_values("u").reset_index(drop=True)
    values = np.vstack([retention[retention.type == c].set_index("front_index").reindex(fronts.front_index).er
                        .to_numpy() for c in types])
    edges, shown = smooth_profiles(fronts.u.to_numpy(), values)
    top = min(1., max(.1, np.ceil(shown.max() / .05) * .05))
    image = ax.imshow(shown, extent=(edges[0], edges[-1], len(types), 0), aspect="auto", interpolation="none",
                      cmap=palette.SIMILARITY, vmin=0, vmax=top, rasterized=True)
    counts = retention.drop_duplicates("type").set_index("type").n
    ax.set_yticks(np.arange(len(types)) + .5, [f"{_fate(c).replace('_', ' ')} ({counts[c]})" for c in types],
                  fontsize=7.5)
    ax.tick_params(axis="y", length=0, pad=12)
    for row, c in enumerate(types):  # fate swatches, as in the pies of Figure 4
        ax.scatter(-.012, row + .5, marker="s", s=26, color=palette.FATE[_fate(c)], edgecolor="#777777", lw=.3,
                   transform=ax.get_yaxis_transform(), clip_on=False)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, .5, 1], [f"0 ({nt.symbol('travel_optimum')})", "0.5", f"1 ({nt.symbol('cell_state_optimum')})"])
    ax.set_xlabel(f"Position along the front, {nt.symbol('canonical_position')}")
    maximum = fronts[fronts.is_max_er.astype(bool)]
    if len(maximum):
        u = float(maximum.u.iloc[0])
        ax.axvline(u, color=SIMILARITY_COLOR, lw=1., ls=(0, (3, 2)), zorder=4)
        ax.scatter(u, -.02, marker="o", s=34, color=SIMILARITY_COLOR, edgecolor="#222222", lw=.6, clip_on=False,
                   transform=ax.get_xaxis_transform(), zorder=5)
        ax.annotate("M", (u, 1.02), xycoords=("data", "axes fraction"), xytext=(5, 0), textcoords="offset points",
                    fontsize=8.5, weight="bold", color=SIMILARITY_COLOR, va="bottom")
    bar = ax.figure.colorbar(image, cax=cax)
    bar.set_ticks([0, top / 2, top], labels=["0", f"{top / 2:.2f}", f"{top:.2f}"])
    bar.ax.tick_params(labelsize=6.5)
    bar.set_label(similarity_label(), fontsize=7.5)
    return ax


def within_type_front(ax, front, n, cell_type):
    """One panel of Figure 5B: the fate class's own front in its own endpoint coordinates."""
    d = _unique_front(front)
    color = palette.FATE[_fate(cell_type)]
    ax.plot(d.endpoint_travel, d.endpoint_state, color=color, lw=1.5, zorder=2)
    ax.scatter(d.endpoint_travel, d.endpoint_state, s=5, color=color, alpha=.5, edgecolors="none", zorder=3)
    _optima(ax, (d.endpoint_travel.iloc[-1], d.endpoint_state.iloc[-1]),
            (d.endpoint_travel.iloc[0], d.endpoint_state.iloc[0]), size=22)
    ax.scatter(d.endpoint_natural_travel.iloc[0], d.endpoint_natural_state.iloc[0], marker="X", s=48,
               color=LINEAGE_RED, edgecolor="white", lw=.5, zorder=8)
    title = f"{cell_type.replace('_', ' ')} ({n})"
    if not bool(d.endpoint_valid.iloc[0]):
        raise ValueError(f"{cell_type}: degenerate endpoint span; Figure 5B shows only valid spans")
    ax.set_title(title, fontsize=7.5, color=palette.NEUTRAL["text"], pad=3)
    ax.tick_params(labelsize=6.5)
    ax.set_xlim(-.05, 1.05)
    ax.set_ylim(-.05, 1.08)
    ax.set_xticks([0, .5, 1])
    ax.set_yticks([0, .5, 1])
    return ax


def restricted_aggregate(ax, aggregate, endpoints):
    """Figure 5C: the unrestricted front and the front restricted within cell type, in unrestricted coordinates."""
    for prefix, color, label in (("global", FRONT_COLOR, "Unrestricted front"),
                                 ("restricted", RESTRICTED_COLOR, "Assignments restricted within cell type")):
        x, y = aggregate[f"{prefix}_endpoint_travel"], aggregate[f"{prefix}_endpoint_state"]
        ax.plot(x, y, color=color, lw=1.7, label=label, zorder=3)
        _optima(ax, (x.iloc[-1], y.iloc[-1]), (x.iloc[0], y.iloc[0]), edge=color)
    ax.scatter(endpoints["endpoint_natural_travel"], endpoints["endpoint_natural_state"], marker="X", s=70,
               color=LINEAGE_RED, edgecolor="white", lw=.6, zorder=8)
    ax.text(.985, .5, f"Endpoint penalty of type restriction:\ntravel +{endpoints['travel_penalty_sigma']:.2f} SD, "
            f"cell state +{endpoints['state_penalty_sigma']:.2f} SD\n(global first-cousin null SD)",
            transform=ax.transAxes, ha="right", va="center", fontsize=7, color=palette.NEUTRAL["muted"])
    ax.set_xlabel(AXIS[0])
    ax.set_ylabel(AXIS[1])
    return ax


def cell_type_figure(retention, within, *, size=(7.4, 7.9)):
    """Figure 5: lineage similarity by fate class along the front (A), each class's own front (B) and the cost of
    restricting assignments within cell types (C). Same data as released Figure 6, drawn in the v2 style
    (2026-11-13)."""
    fig = plt.figure(figsize=size)
    grid = fig.add_gridspec(3, 1, height_ratios=(1., 1.75, 1.3), hspace=.5, left=.1, right=.97, top=.95,
                            bottom=.12)
    top = grid[0].subgridspec(1, 2, width_ratios=(1, .025), wspace=.03)
    ax = fig.add_subplot(top[0])
    ax_left = .21  # room for the fate labels and swatches
    box = ax.get_position()
    ax.set_position([ax_left, box.y0, box.x1 - ax_left, box.height])
    similarity_heatmap(ax, fig.add_subplot(top[1]), retention)
    _panel_letter(ax, "A", x=-(ax_left - .02) / (box.x1 - ax_left), y=1.04)

    fronts, summary = within["fronts"], within["summary"].set_index("type")
    middle = grid[1].subgridspec(2, 3, hspace=.55, wspace=.3)
    axes = [fig.add_subplot(middle[i, j]) for i in range(2) for j in range(3)]
    for k, (ax, cell_type) in enumerate(zip(axes, WITHIN_TYPES)):
        within_type_front(ax, fronts[fronts.type == cell_type], int(summary.loc[cell_type, "n"]), cell_type)
        if k >= 3:
            ax.set_xlabel(AXIS[0], fontsize=7.5)
        if k % 3 == 0:
            ax.set_ylabel(AXIS[1], fontsize=7.5)
    _panel_letter(axes[0], "B", x=-.42, y=1.08)

    ax = fig.add_subplot(grid[2])
    restricted_aggregate(ax, within["aggregate"], within["endpoints"])
    _panel_letter(ax, "C", x=-.075, y=1.04)
    handles = [Line2D([], [], color=FRONT_COLOR, lw=1.7, label="Unrestricted front"),
               Line2D([], [], color=RESTRICTED_COLOR, lw=1.7, label="Restricted within cell type"),
               Line2D([], [], marker="X", ls="", color=LINEAGE_RED, markersize=7, label="C. elegans lineage"),
               Line2D([], [], marker="s", ls="", color=TRAVEL_COLOR, markeredgecolor="#222222", markersize=5.5,
                      label=f"Travel optimum {nt.symbol('travel_optimum')}"),
               Line2D([], [], marker="s", ls="", color=STATE_COLOR, markeredgecolor="#222222", markersize=5.5,
                      label=f"Cell-state optimum {nt.symbol('cell_state_optimum')}")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7, frameon=False, bbox_to_anchor=(.53, 0),
               handletextpad=.3, columnspacing=1.)
    return fig


# ---------------------------------------------------------------- tracking geometries (Figure S1)

TRACKING_SOURCES = ("embryo1", "embryo2", "embryo3", "pooled")
TRACKING_LABELS = {"embryo1": "Embryo 1", "embryo2": "Embryo 2", "embryo3": "Embryo 3", "pooled": "Pooled (all three)"}
NEAR_NATURAL = {"travel": (-.026, .036), "state": (-.09, .09)}  # inset window around the lineage, as released S3


def tracking_figure(projections, references, *, size=(7.4, 7.)):
    """Figure S1: the four saved sweeps (optimized in embryo 1, 2, 3 or pooled travel) measured in each geometry
    (A-D), each panel in its own geometry's endpoint coordinates, with an inset near the lineage. Same data as
    released S3, drawn in the v2 style (2026-11-13)."""
    fig, axes = plt.subplots(2, 2, figsize=size)
    fig.subplots_adjust(left=.09, right=.98, bottom=.07, top=.85, hspace=.42, wspace=.28)
    sources = [Line2D([], [], color=palette.TRACKING[s], lw=2, label=TRACKING_LABELS[s]) for s in TRACKING_SOURCES]
    fig.legend(handles=sources, title="Assignments optimized using", loc="upper center", bbox_to_anchor=(.53, .995),
               ncol=4, frameon=False, fontsize=7.5, title_fontsize=8, columnspacing=1.6, handletextpad=.5)
    styles = [Line2D([], [], color=palette.NEUTRAL["muted"], lw=1.8, label="Native front"),
              Line2D([], [], color=palette.NEUTRAL["muted"], lw=1.4, ls=(0, (4, 2)), label="Transferred assignments"),
              Line2D([], [], marker="X", ls="", color=LINEAGE_RED, markersize=7, label="C. elegans lineage")]
    fig.legend(handles=styles, loc="upper center", bbox_to_anchor=(.53, .925), ncol=3, frameon=False, fontsize=7.5,
               columnspacing=2., handletextpad=.5)
    x_all, y_all = projections.x.to_numpy(), projections.y.to_numpy()
    for letter, ax, target in zip("ABCD", axes.flat, TRACKING_SOURCES):
        inset = ax.inset_axes([.40, .43, .56, .51])
        for source in [s for s in TRACKING_SOURCES if s != target] + [target]:  # the native front on top
            d = projections[(projections.target == target) & (projections.source == source)].sort_values("step")
            native = source == target
            for where in (ax, inset):
                where.plot(d.x, d.y, color=palette.TRACKING[source], lw=1.8 if native else 1.1,
                           ls="-" if native else (0, (4, 2)), alpha=1 if native else .85)
        nx, ny = references[target]["natural_display"]
        ax.scatter(nx, ny, marker="X", s=55, color=LINEAGE_RED, edgecolor="white", lw=.6, zorder=8)
        inset.scatter(nx, ny, marker="X", s=36, color=LINEAGE_RED, edgecolor="white", lw=.5, zorder=8)
        inset.set_xlim(nx + NEAR_NATURAL["travel"][0], nx + NEAR_NATURAL["travel"][1])
        inset.set_ylim(ny + NEAR_NATURAL["state"][0], ny + NEAR_NATURAL["state"][1])
        inset.set_title("Near the lineage", fontsize=6.5, pad=3)
        inset.tick_params(labelsize=5.5, pad=1.5, length=2)
        inset.locator_params(axis="both", nbins=3)
        for spine in ("top", "right"):
            inset.spines[spine].set_visible(True)
        for spine in inset.spines.values():
            spine.set_color(palette.NEUTRAL["reference"])
            spine.set_linewidth(.6)
        ax.indicate_inset_zoom(inset, edgecolor=palette.NEUTRAL["reference"], alpha=.7, lw=.6)
        title = "Travel measured in pooled geometry" if target == "pooled" else f"Travel measured in embryo {target[-1]}"
        ax.set_title(title, loc="left", fontsize=8.5, pad=6)
        _panel_letter(ax, letter, x=-.2, y=1.04)
        ax.set_xlim(min(-.015, x_all.min() - .015), max(1.025, x_all.max() + .015))
        ax.set_ylim(min(-.015, y_all.min() - .015), max(1.025, y_all.max() + .015))
        ax.set_xlabel(AXIS[0], fontsize=8)
        ax.set_ylabel(AXIS[1], fontsize=8)
        ax.tick_params(labelsize=7)
    return fig
