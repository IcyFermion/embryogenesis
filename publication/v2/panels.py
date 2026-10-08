"""Panels for the revised layout: plain fronts with direct labels, lineage similarity, rewiring, metric rows."""

from __future__ import annotations

from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MPath
import numpy as np

from publication import notation as nt
from publication import style

FRONT_COLOR = style.COLORS["black"]
LINEAGE_RED = "#D7191C"
SIMILARITY_COLOR = style.COLORS["blue"]
TRAVEL_COLOR = style.SEMANTIC_COLORS["travel"]
STATE_COLOR = style.SEMANTIC_COLORS["cell_state"]
CLOUD_LABELS = {"first_cousin": "First-cousin shuffle", "second_cousin": "Second-cousin shuffle",
                "third_cousin": "Third-cousin shuffle"}
DISTANCE_CMAP = LinearSegmentedColormap.from_list("tree_distance", ["#F2D49B", "#D9822B", "#A33B3B", "#4B1D4F"])
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
        color = style.NULL_MODEL_COLORS[family]
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
        color = style.NULL_MODEL_COLORS["full_random"]
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
                    color=LINEAGE_RED if lineage else style.NULL_MODEL_COLORS[key],
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
            _, box = _place_label(ax, CLOUD_LABELS[family], anchor, everything, style.NULL_MODEL_COLORS[family],
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


def draw_dendrogram(ax, x, depth, children, *, baseline, color="#555555", lw=.45):
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


def metric_columns(fig, left, bottom, height, sections, *, row_count, dmax, pies=None, label_width=.105,
                   pie_width=.035, column_width=.08, gap=.03):
    """Figure 4B-style rows: [labels | pie | position u | distance d_LP] for grouped ``sections``.

    ``sections``: [(heading, [(label, u, d_lp, key), ...]), ...]; ``pies``: {key: Counter} or None.
    """
    x = left + label_width
    pie_ax = None
    if pies is not None:
        pie_ax = fig.add_axes([x, bottom, pie_width, height])
        x += pie_width + .008
    u_ax = fig.add_axes([x, bottom, column_width, height])
    d_ax = fig.add_axes([x + column_width + gap, bottom, column_width, height])
    for ax, xmax, title, color in ((u_ax, 1., f"{nt.label('canonical_position')}\n{nt.math('canonical_position')}",
                                    style.SEMANTIC_COLORS["canonical_u"]),
                                   (d_ax, dmax, f"Distance\n{nt.math('natural_front_distance')}",
                                    style.SEMANTIC_COLORS["lineage_front_distance"])):
        ax.set_ylim(row_count - .3, -1.2)
        ax.set_xlim(-.07 * xmax, 1.07 * xmax)
        ax.set_xticks([0, xmax / 2, xmax])
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, labelsize=7.5,
                       length=2, pad=2)
        ax.set_yticks([])
        ax.set_title(title, fontsize=8.3, color=color, pad=17, linespacing=1.1)
        ax.spines[["left", "right", "bottom"]].set_visible(False)
        ax.spines["top"].set_color("#BBBBBB")
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
    row = 0
    for heading, entries in sections:
        if heading:
            u_ax.text(heading_x, row, heading, transform=u_ax.get_yaxis_transform(), fontsize=8.5, weight="bold",
                      va="center", clip_on=False)
        for label, u, d, key in entries:
            row += 1
            u_ax.text(label_x, row, label, transform=u_ax.get_yaxis_transform(), ha="right", va="center",
                      fontsize=8.3, clip_on=False)
            for ax, value, color in ((u_ax, u, style.SEMANTIC_COLORS["canonical_u"]),
                                     (d_ax, d, style.SEMANTIC_COLORS["lineage_front_distance"])):
                ax.hlines(row, 0, value, color=color, lw=.8, alpha=.45)
                ax.scatter(value, row, color=color, s=18, edgecolor="white", linewidth=.35, zorder=3)
            if pie_ax is not None and key in pies:
                _pie(pie_ax, (pie_ax.get_xlim()[1] / 2, row), .42, pies[key], PIE_COLORS)
        row += 2
    return u_ax, d_ax


def _pie_colors():
    from publication.figures.terminal_subtree_map import TYPE_COLORS
    return TYPE_COLORS


PIE_COLORS = _pie_colors()


def fate_legend(fig, types, *, y=.02):
    handles = [Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor=PIE_COLORS[t], markeredgecolor="white",
                      label=t.replace("_", " ")) for t in types]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, y), ncol=min(len(types), 6), fontsize=7.5,
               frameon=False, handletextpad=.2, columnspacing=1.1)


def subtree_rows(canonical, nodes, *, size=(7.4, 5.4)):
    """Figure 4B: every subtree with fate pie, position u and distance d_LP (no null)."""
    metrics = canonical
    root = metrics[metrics.subtree.isin(["P0", "AB"])].set_index("subtree")
    root = root.loc[[n for n in ("P0", "AB") if n in root.index]].reset_index()
    groups = {name: metrics[metrics.region == name].sort_values("subtree") for name in ("ABa", "ABp", "P1")}
    p1_order = ["P1", "EMS", "MS", "MSa", "MSaa", "MSap", "MSp", "MSpa", "MSpp", "P2", "C", "P3", "D"]
    rank = {n: i for i, n in enumerate(p1_order)}
    groups["P1"] = groups["P1"].sort_values("subtree", key=lambda s: s.map(rank).fillna(len(rank)))

    def entries(frame):
        return [(r.subtree, r.u_lineage_lp, r.d_lp, r.subtree) for r in frame.itertuples()]
    columns = [[("Root & AB branch", entries(root)), ("ABa branch", entries(groups["ABa"]))],
               [("ABp branch", entries(groups["ABp"]))],
               [("P1 branch", entries(groups["P1"]))]]
    shown = [e[3] for column in columns for _, es in column for e in es]
    if len(shown) != len(set(shown)) or set(shown) != set(metrics.subtree):
        raise ValueError("Every subtree must appear exactly once")
    row_count = max(sum(len(es) + 2 for _, es in column) - 1 for column in columns)
    dmax = max(.2, np.ceil(float(metrics.d_lp.max()) * 10) / 10)
    pies = {name: nodes[name]["composition"] for name in shown if name in nodes}
    fig = plt.figure(figsize=size)
    for left, column in zip((.0, .335, .668), columns):
        metric_columns(fig, left, .1, .8, column, row_count=row_count, dmax=dmax, pies=pies)
    types = sorted({t for c in pies.values() for t in c}, key=lambda t: -sum(c[t] for c in pies.values()))
    fate_legend(fig, types)
    return fig


def layer_rows(metrics, *, size=(3.6, 3.0)):
    """Figure 6C: position u and distance d_LP for each layer (layer 1 = top)."""
    metrics = metrics.sort_values("layer")
    entries = [(f"Layer {int(r.layer)} \u00b7 {int(r.edges)}", r.canonical_position, r.natural_front_distance, int(r.layer))
               for r in metrics.itertuples()]
    fig = plt.figure(figsize=size)
    dmax = max(.2, np.ceil(float(metrics.natural_front_distance.max()) * 10) / 10)
    metric_columns(fig, .02, .06, .78, [("", entries)], row_count=len(entries) + 1, dmax=dmax,
                   label_width=.2, column_width=.3, gap=.08)
    return fig


def layer_similarity(by_layer, *, size=(3.6, 3.0)):
    """Figure 6B: lineage similarity along each layer's front, layer 1 (top) to layer 8 (bottom)."""
    fig, ax = plt.subplots(figsize=size)
    keys = sorted(by_layer)
    colors = style.color_ramp(SIMILARITY_COLOR, len(keys), lightness=(.75, 0.))
    for color, layer in zip(colors, keys):
        curve, metric = by_layer[layer].select("raw3d", "ce_protein")
        ax.plot(curve.D1, curve.edge_retention, color=color, lw=1.3, label=f"Layer {layer}")
        top = curve.iloc[int(metric.maximum_index)]
        ax.scatter(top.D1, top.edge_retention, s=16, color=color, edgecolor="#222222", lw=.5, zorder=4)
    ax.set_xlabel(AXIS[0])
    ax.set_ylabel(similarity_label())
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=6.5, ncol=4, frameon=False, loc="upper center", bbox_to_anchor=(.5, -.24),
              handlelength=1.2, columnspacing=.8)
    fig.tight_layout()
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
    handles = [Line2D([], [], marker="o", ls="", color=style.NULL_MODEL_COLORS[f], markersize=5, label=CLOUD_LABELS[f])
               for f in CLOUD_LABELS]
    handles += [Line2D([], [], marker="o", ls="", color=style.NULL_MODEL_COLORS["full_random"], markersize=5,
                       label="Full-random shuffle (insets)"),
                Line2D([], [], color=FRONT_COLOR, lw=1.4, label="Pareto front"),
                Line2D([], [], marker="X", ls="", color=LINEAGE_RED, markersize=7, label="C. elegans lineage")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7, frameon=False, bbox_to_anchor=(.5, 0),
               handletextpad=.3, columnspacing=1.)
    return fig
