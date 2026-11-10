"""Panels for the revised layout: plain fronts with direct labels, lineage similarity, rewiring, metric rows."""

from __future__ import annotations

from collections import Counter

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
from publication.figures import terminal_canonical

FRONT_COLOR = palette.MEASURE["front"]
LINEAGE_RED = palette.MEASURE["natural_lineage"]
SIMILARITY_COLOR = palette.MEASURE["lineage_similarity"]
TRAVEL_COLOR = palette.MEASURE["travel"]
STATE_COLOR = palette.MEASURE["cell_state"]
U_COLOR = palette.MEASURE["canonical_position"]
D_COLOR = palette.MEASURE["front_distance"]
R_COLOR = palette.MEASURE["cousin_relative"]
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
    from publication.figures.terminal_subtree_map import terminal_inorder
    from publication.v2.data import edge_changes
    from terminal_pareto import data_loader as dl
    order = terminal_inorder(t.lineage, dl.map_names)
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


# One metric column of the subtree rows: (field, heading, color, axis maximum, reference value or None).
def subtree_metrics(canonical):
    """Figure 5 columns: position u, distance d_LP and the cousin-relative distance r of former S1."""
    dmax = max(.2, np.ceil(float(canonical.d_lp.max()) * 10) / 10)
    return [("u_lineage_lp", f"{nt.label('canonical_position')}\n{nt.math('canonical_position')}", U_COLOR, 1., None),
            ("d_lp", f"Distance\n{nt.math('natural_front_distance')}", D_COLOR, dmax, None),
            ("relative_distance", f"Cousin-relative\ndistance {nt.math('cousin_relative')}", R_COLOR, 1., 1.)]


def metric_columns(fig, left, bottom, height, sections, metrics, *, row_count, pies=None, marks=None, pairs=(),
                   label_width=.1, pie_width=.035, column_width=.085, gap=.028):
    """Subtree rows: [labels | pie | one column per metric] for grouped ``sections``.

    ``sections``: [(heading, [(label, {field: value}, key), ...]), ...]; ``metrics``: see ``subtree_metrics``;
    ``pies``: {key: Counter} or None; ``marks``: {key: (note mark, field)} as ``terminal_canonical.NOTES``;
    ``pairs``: (key, key) left-right counterparts, joined by a bracket left of their labels when adjacent.
    """
    x = left + label_width
    pie_ax = None
    if pies is not None:
        pie_ax = fig.add_axes([x, bottom, pie_width, height])
        x += pie_width + .008
    axes = {}
    for i, (field, title, color, xmax, reference) in enumerate(metrics):
        ax = fig.add_axes([x + i * (column_width + gap), bottom, column_width, height])
        axes[field] = ax
        ax.set_ylim(row_count - .3, -1.2)
        ax.set_xlim(-.07 * xmax, 1.07 * xmax)
        ax.set_xticks([0, xmax / 2, xmax])
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, labelsize=7.5,
                       length=2, pad=2)
        ax.set_yticks([])
        ax.set_title(title, fontsize=8.3, color=color, pad=17, linespacing=1.1)
        ax.spines[["left", "right", "bottom"]].set_visible(False)
        ax.spines["top"].set_color(palette.NEUTRAL["rule"])
        if reference is not None:
            ax.axvline(reference, color=palette.NEUTRAL["reference"], lw=.6, ls=(0, (2, 2)), zorder=0)
    first = axes[metrics[0][0]]
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
            for field, _, color, xmax, _ in metrics:
                ax = axes[field]
                ax.hlines(row, 0, values[field], color=color, lw=.8, alpha=.45)
                ax.scatter(values[field], row, color=color, s=18, edgecolor="white", linewidth=.35, zorder=3)
            if pie_ax is not None and key in pies:
                _pie(pie_ax, (pie_ax.get_xlim()[1] / 2, row), .42, pies[key], PIE_COLORS)
            if key in (marks or {}):
                mark, field = marks[key]
                xmax = next(m[3] for m in metrics if m[0] == field)
                terminal_canonical.mark_value(axes[field], row, values[field], mark, xmax)
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


def fate_legend(fig, types, *, y=.02, x=.5, ncol=6, loc="lower center"):
    handles = [Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor=PIE_COLORS[t], markeredgecolor="white",
                      label=t.replace("_", " ")) for t in types]
    fig.legend(handles=handles, loc=loc, bbox_to_anchor=(x, y), ncol=min(len(types), ncol), fontsize=7.5,
               frameon=False, handletextpad=.2, columnspacing=1.1, title="Terminal fates", title_fontsize=7.5,
               alignment="left")


def subtree_rows(canonical, nodes, *, notes=None, width=7.4, row_inch=.241):
    """Figure 5: every displayed subtree (at least 24 cells, Track D6) with fate pie, position u, distance d_LP
    and cousin-relative distance r, left-right counterparts adjacent and bracketed.

    Two columns of groups: Root & AB and ABa on the left (the fate legend fills the space below), ABp and P1
    on the right. With ``notes`` (the Track D evidence, checked by ``terminal_canonical.check_notes``), ABplp
    carries the dagger of production Figure 5B.
    """
    table = canonical.set_index("subtree")
    order = [name for names in SUBTREE_ORDER.values() for name in names]
    if sorted(order) != sorted(table.index):
        raise ValueError(f"Subtree order and data differ: {sorted(set(order) ^ set(table.index))}")
    metrics = subtree_metrics(canonical)

    def section(heading):
        return heading, [(n, {m[0]: float(table.loc[n, m[0]]) for m in metrics}, n) for n in SUBTREE_ORDER[heading]]
    columns = [[section("Root & AB branch"), section("ABa branch")],
               [section("ABp branch"), section("P1 branch")]]
    row_count = max(sum(len(es) + 2 for _, es in column) - 1 for column in columns)
    pies = {name: nodes[name]["composition"] for name in order if name in nodes}
    marks = None
    if notes:
        terminal_canonical.check_notes(notes)
        marks = terminal_canonical.NOTES
    top, bottom = .62, .2  # inches: column heads above, footnote below
    rows_height = row_inch * (row_count + .9)
    height = rows_height + top + bottom
    fig = plt.figure(figsize=(width, height))
    for left, column in zip((.0, .5), columns):
        metric_columns(fig, left, bottom / height, rows_height / height, column, metrics, row_count=row_count,
                       pies=pies, marks=marks, pairs=SUBTREE_PAIRS)
    types = sorted({t for c in pies.values() for t in c}, key=lambda t: -sum(c[t] for c in pies.values()))
    left_rows = sum(len(es) + 2 for _, es in columns[0]) - 1
    # Rows map to the axes as (row + 1.2) / (row_count + .9) from the top; the legend starts below the left column.
    legend_top = (bottom + rows_height * (row_count - .3 - (left_rows + .2)) / (row_count + .9)) / height
    fate_legend(fig, types, x=.04, y=legend_top, ncol=3, loc="upper left")
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

def layer_position(table, *, dmax=None, size=None):
    """Distance from the lineage to its front and the lowest % of edges differing from it, one point per layer.

    ``table``: layer, born, front_distance, fewest_rewired. Figure 7B (C. elegans) and 8C (C. briggsae) pass
    the same ``dmax`` so their distance axes compare directly.
    """
    from publication.figures import full_tree_layers as ftl
    t = table
    fig, axes = plt.subplots(2, 1, figsize=size or ftl.SUMMARY_SIZE, sharex=True)
    fig.subplots_adjust(left=.2, right=.96, top=.9, bottom=.14, hspace=.22)
    dmax = dmax or max(.3, np.ceil(t.front_distance.max() * 20) / 20 + .05)
    rows = ((axes[0], "front_distance", 1, D_COLOR,
             f"Distance from lineage\nto front, {nt.math('natural_front_distance')}", dmax),
            (axes[1], "fewest_rewired", 100, SIMILARITY_COLOR, "Lowest % of edges\ndiffering, front\nvs lineage", 100))
    for ax, column, scale, color, label, top in rows:
        ftl._size_axis(ax, label=ax is axes[-1])
        ax.plot(t.born, t[column] * scale, color=color, lw=1.2, zorder=3)
        ax.scatter(t.born, t[column] * scale, s=18, color=color, zorder=4, clip_on=False)
        ax.set_ylim(0, top)
        ax.set_ylabel(label)
    ftl._layer_axis(axes[0], t.born, t.layer)
    return fig


def layer_cost(table):
    """Figure 7C: travel and cell-state cost relative to the first-cousin shuffle, with 95% bootstrap intervals."""
    from publication.figures import full_tree_layers as ftl
    t = table
    fig, ax = plt.subplots(figsize=ftl.SUMMARY_SIZE)
    fig.subplots_adjust(left=.2, right=.96, top=.9, bottom=.14)
    ftl._size_axis(ax)
    ax.axhline(1, color=palette.NEUTRAL["reference"], lw=.6, zorder=1)
    ax.text(.97, 1.005, "no better than the shuffle", transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=6.5, color="#777777")
    for column, color, label, shift in (("travel", TRAVEL_COLOR, "Travel", 1 / ftl.OFFSET),
                                        ("state", STATE_COLOR, "Cell state", ftl.OFFSET)):
        x = t.born * shift
        ax.vlines(x, t[f"{column}_low"], t[f"{column}_high"], color=color, lw=.9, zorder=3)
        ax.plot(x, t[column], color=color, lw=1.2, zorder=3, label=label)
        ax.scatter(x, t[column], s=18, color=color, zorder=4)
    ax.set_ylim(*ftl._relative_limits(t))
    ax.set_ylabel(ftl.RELATIVE_LABEL)
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    ftl._layer_axis(ax, t.born, t.layer)
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
