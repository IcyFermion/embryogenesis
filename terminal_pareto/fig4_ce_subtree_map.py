"""Figure 4: C. elegans subtree map of the investigated lineages.

Renders the investigated subtrees as pies in a tidy lineage-order layout
(root P0 at the top): every subtree with at least ``--min-cells``
usable biological terminal cells is a node whose pie shows the terminal fate
composition of its USABLE terminal cells (the same cells used in the Pareto
analysis) and whose circumference style marks whether the natural assignment
is exactly Pareto-optimal (solid) or not (dashed). Subtree label color encodes
the maximum natural-edge retention attained on that subtree's Pareto front.
Edges connect each
subtree to its nearest qualifying descendant. Non-qualifying branches are
omitted entirely: they are not part of this discussion. Renders
``output/publication/fig4_ce_subtree_map_panel.{pdf,png}``.

Visible endpoints are spaced according to pie and label width, parents are
centered over their displayed descendants, and crowded depth levels use
small vertical node offsets. Thus labels stay attached to their pies rather
than being independently displaced.
"""

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle, Rectangle, Wedge

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import plot_style as ps
from terminal_pareto.subtree_analysis import build_type_map

OUT = Path(__file__).resolve().parent / "output" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"

# Terminal fate palette (sentence-case type -> hex color).
TYPE_COLORS = {
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
LABEL_COLOR = "#222222"
ER_NORM = Normalize(vmin=0.6, vmax=1.0)
# Text-safe companion to the shared navy-to-cyan retention scale.  Its pale
# endpoint is compressed so five-point labels retain contrast on white.
ER_TEXT_CMAP = LinearSegmentedColormap.from_list(
    "edge_retention_text_blue", ["#17365D", "#006A9F", "#32A4D2"]
)

# X span of the displayed-tree map. A wider internal coordinate grid gives
# the many qualifying terminal branches distinct positions; the reciprocal
# canvas scale below keeps the final publication width essentially unchanged.
X_SPAN = 20.0
Y_LEVEL = 1.0   # vertical distance per depth level
INCH_PER_UNIT = 0.38   # figure scale: inches of canvas per data unit

# Labels shown on the map. The layout uses this same set when reserving
# horizontal space, so readability is achieved by positioning tree nodes
# rather than detaching text from them after drawing.
LABEL_SET = {
    "P0", "AB", "ABa", "ABp", "P1", "EMS", "MS", "C", "D", "ABal",
    "ABpr", "ABpl", "ABplpp", "ABprp", "ABala", "ABprppp", "ABplppp",
    "MSa", "MSp", "ABalp", "ABara", "ABpra", "ABprpa", "ABprpap",
    "ABplpa", "ABalpa", "ABalpp", "ABarap", "ABpla", "ABplaa",
    "ABpraa", "ABarp", "MSaa", "MSap", "MSpa", "MSpp", "ABaraa",
    "ABprpp", "P3", "P2", "ABalap", "ABar", "ABplp",
}


def terminal_inorder(root):
    """Assign each terminal cell an inorder index 0..n-1 (left-to-right)."""
    order = []

    def dfs(node):
        children = node.get("children", [])
        if not children:
            order.append(dl.map_names(node["did"]))
            return
        for c in children:
            dfs(c)

    dfs(root)
    return order


def collect_nodes(root, tree_index, qualifying, terminal_order, type_map,
                  usable):
    """Positions and compositions for the subtree map.

    Only USABLE biological terminal cells (``usable``, the terminal set used
    by the Pareto analysis) contribute to the pies and counts; all other
    leaves are ignored. Returns (nodes, edges, positions, depth, max_depth):
    nodes is a dict name -> dict(x, y, n, on_front, composition); edges is a
    list of (parent_name, child_name).
    """
    idx = {name: i for i, name in enumerate(terminal_order)}
    n_term = len(terminal_order)
    positions = {}
    comp = {}

    def walk(node):
        name = dl.map_names(node.get("did", ""))
        children = node.get("children", [])
        if not children:
            # Terminal cell: x = its inorder index. Composition counts only
            # biological (usable) terminal cells.
            positions[name] = (idx[name], 0)
            comp[name] = (Counter({type_map.get(name, "other"): 1})
                          if name in usable else Counter())
            return
        xs = []
        for c in children:
            walk(c)
            cname = dl.map_names(c.get("did", ""))
            xs.append(positions[cname][0])
        # x of an internal node = midpoint of its children's xs
        positions[name] = (float(np.mean(xs)) if xs else 0.0, 0.0)
        comp[name] = Counter()
        for c in children:
            cname = dl.map_names(c.get("did", ""))
            comp[name].update(comp.get(cname, Counter()))
        return

    walk(root)

    depth = {}
    for name, info in tree_index.items():
        depth[name] = info["depth"]
    max_depth = max(depth.values()) if depth else 0

    qualifying = {q for q in qualifying}
    nodes = {}
    edges = []
    for name, (x, _y) in positions.items():
        if name not in qualifying:
            continue
        c = comp.get(name, Counter())
        n = sum(c.values())
        nodes[name] = dict(
            x=X_SPAN * (x + 0.5) / (n_term + 1),
            # Root (depth 0) at the TOP: y grows with the distance from the
            # root, so the deepest nodes sit at y = 0 and the root at the top.
            y=(max_depth - depth[name]) * Y_LEVEL,
            n=n,
            composition=c,
        )
    # Edges: parent subtree -> nearest qualifying descendant per branch.
    def parent_of(name):
        return tree_index[name]["parent"]

    for name in nodes:
        p = parent_of(name)
        while p is not None and p not in nodes:
            p = parent_of(p)
        if p is not None:
            edges.append((p, name))
    return nodes, edges, positions, depth, max_depth


def tidy_display_tree(nodes, edges, radii):
    """Lay out the displayed subtree while preserving lineage order.

    Terminal-inorder averages compress branches in proportion to biological
    cell count, which is useful analytically but causes the *displayed*
    qualifying subtrees to pile up. Here visible leaves receive consecutive
    slots with gaps large enough for their pies and labels; internal nodes
    are then centered over the span of their visible descendants. Original
    x coordinates determine sibling order, so no lineage branches cross.
    """
    children = {name: [] for name in nodes}
    parents = {}
    for parent, child in edges:
        children[parent].append(child)
        parents[child] = parent
    for name in children:
        children[name].sort(key=lambda child: nodes[child]["x"])

    roots = sorted((name for name in nodes if name not in parents),
                   key=lambda name: nodes[name]["x"])
    leaves = []

    def collect_leaves(name):
        if not children[name]:
            leaves.append(name)
            return
        for child in children[name]:
            collect_leaves(child)

    for root in roots:
        collect_leaves(root)

    # Approximate the rendered 5-point label width in data units. Reserving
    # this at the node-layout stage makes the later bbox pass a verification
    # and small safety net rather than the primary layout mechanism.
    def half_extent(name):
        pie_extent = radii[name] + 0.08
        if name not in LABEL_SET:
            return pie_extent
        label = f"{name} {nodes[name]['n']}"
        label_extent = 0.043 * len(label) + 0.08
        return max(pie_extent, label_extent)

    leaf_x = {}
    cursor = 0.0
    previous = None
    for leaf in leaves:
        if previous is None:
            cursor = half_extent(leaf)
        else:
            cursor += half_extent(previous) + half_extent(leaf) + 0.16
        leaf_x[leaf] = cursor
        previous = leaf

    laid_out = {}

    def place(name):
        if not children[name]:
            x = leaf_x[name]
        else:
            child_x = [place(child) for child in children[name]]
            # Midrange centering keeps the parent visually above the entire
            # descendant span rather than weighting it by branch density.
            x = 0.5 * (min(child_x) + max(child_x))
        laid_out[name] = x
        return x

    for root in roots:
        place(root)

    xmin = min(laid_out.values())
    xmax = max(laid_out.values())
    if xmax > xmin:
        for name, x in laid_out.items():
            nodes[name]["x"] = 0.35 + (X_SPAN - 0.70) * (x - xmin) / (xmax - xmin)
    else:
        for name in laid_out:
            nodes[name]["x"] = X_SPAN / 2.0


def stagger_crowded_levels(nodes, radii):
    """Move crowded same-depth nodes into two nearby vertical lanes.

    The shift is applied to the node and therefore to its incident edges,
    not merely to its label. A 0.38-unit offset is well below one lineage
    depth but exceeds a five-point label's rendered height at this scale.
    """
    levels = {}
    for name, nd in nodes.items():
        if name in LABEL_SET:
            levels.setdefault(nd["y"], []).append(name)

    def half_label_width(name):
        label = f"{name} {nodes[name]['n']}"
        # Slightly conservative relative to the average glyph width: labels
        # such as ABplpp/ABpraa and MSpa/MSpp contain many wide capitals.
        return max(radii[name], 0.060 * len(label) + 0.12)

    for names in levels.values():
        names.sort(key=lambda name: nodes[name]["x"])
        lane_right = [-float("inf"), -float("inf")]
        for name in names:
            x = nodes[name]["x"]
            half = half_label_width(name)
            left = x - half
            lane = next((i for i, right in enumerate(lane_right)
                         if left >= right + 0.10), None)
            if lane is None:
                lane = int(np.argmin(lane_right))
            nodes[name]["y"] += 0.38 * lane
            lane_right[lane] = x + half


def draw_retention_label_key(ax, x0, x1, y, height=0.18, n_steps=80):
    """Draw the text-safe maximum-retention ramp as one legend row."""
    for i in range(n_steps):
        lo = i / n_steps
        hi = (i + 1) / n_steps
        value = 0.6 + 0.4 * (i + 0.5) / n_steps
        ax.add_patch(Rectangle(
            (x0 + (x1 - x0) * lo, y - height / 2),
            (x1 - x0) * (hi - lo), height,
            facecolor=ER_TEXT_CMAP(ER_NORM(value)), edgecolor="none", zorder=6,
        ))
    ax.add_patch(Rectangle(
        (x0, y - height / 2), x1 - x0, height, fill=False,
        edgecolor="#666666", lw=0.45, zorder=7,
    ))
    ax.text(x0, y - height / 2 - 0.10, "0.6", ha="center", va="top",
            fontsize=5.7, color=LABEL_COLOR)
    ax.text(x1, y - height / 2 - 0.10, "1.0", ha="center", va="top",
            fontsize=5.7, color=LABEL_COLOR)
    ax.text(x1 + 0.36, y,
            "subtree-label color: maximum edge retention",
            ha="left", va="center", fontsize=6.1, color=LABEL_COLOR)


def emit_figure(nodes, edges, min_cells, out_path):
    """Render the subtree map to PDF + PNG with matplotlib.

    Only the investigated subtrees are drawn (non-qualifying nodes are
    omitted), so the view is cropped to the pie region plus the legend.
    """
    rs = {name: 0.10 + 0.085 * math.sqrt(nd["n"] / 12.0)
          for name, nd in nodes.items()}
    tidy_display_tree(nodes, edges, rs)
    stagger_crowded_levels(nodes, rs)
    y_bottom = min(nd["y"] - rs[name] - 0.10 for name, nd in nodes.items())
    y_top = max(nd["y"] + rs[name] for name, nd in nodes.items())
    xlim = (-0.6, X_SPAN + 0.6)
    ylim = (y_bottom - 2.95, y_top + 0.4)
    fig, ax = plt.subplots(
        figsize=((xlim[1] - xlim[0]) * INCH_PER_UNIT,
                 (ylim[1] - ylim[0]) * INCH_PER_UNIT))

    # Pie-to-pie backbone: parent subtree -> nearest qualifying descendant.
    for p, c in edges:
        ax.plot([nodes[p]["x"], nodes[c]["x"]],
                [nodes[p]["y"], nodes[c]["y"]],
                color="#9A9A9A", lw=0.6, zorder=1)

    # Pies + circumference style (solid = on front, dashed = off front).
    for name, nd in nodes.items():
        x, y, r = nd["x"], nd["y"], rs[name]
        comp = nd["composition"]
        total = sum(comp.values())
        # Wedges, starting at 90 deg, sweeping clockwise.
        angle = 90.0
        for tname, cnt in comp.most_common():
            if cnt <= 0:
                continue
            sweep = 360.0 * cnt / total
            a1 = angle - sweep
            ax.add_patch(Wedge((x, y), r, a1, angle,
                               facecolor=TYPE_COLORS[tname],
                               edgecolor="#FFFFFF", lw=0.3, zorder=3))
            angle = a1
        ax.add_patch(Circle((x, y), r, fill=False,
                            edgecolor=LABEL_COLOR, lw=1.1,
                            linestyle="solid" if nd["on_front"] else "dashed",
                            zorder=4))

    # Labels: on-front subtrees, major roots, and the off-front outliers.
    label_artists = []
    for name, nd in nodes.items():
        if name not in LABEL_SET:
            continue
        label_color = ER_TEXT_CMAP(ER_NORM(nd["max_er"]))
        t = ax.text(nd["x"], nd["y"] - rs[name] - 0.10, f"{name} {nd['n']}",
                    ha="center", va="top", fontsize=5.0,
                    fontweight="medium", color=label_color, zorder=6)
        label_artists.append((name, t, nd["x"], nd["y"] - rs[name] - 0.10))

    # De-overlap pass: nudge horizontally until text boxes stop colliding
    # (measured, not heuristic-only). Keeps labels inside the x range.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    xlo, xhi = xlim
    n_collide0 = None
    for it in range(10):
        boxes = [t.get_window_extent(renderer) for _, t, _, _ in label_artists]
        if n_collide0 is None:
            n_collide0 = sum(1 for i in range(len(boxes))
                             for j in range(i + 1, len(boxes))
                             if boxes[i].overlaps(boxes[j]))
        moved = False
        for i in range(len(label_artists)):
            for j in range(i + 1, len(label_artists)):
                if boxes[i].overlaps(boxes[j]):
                    # Push both labels apart horizontally, alternating sides.
                    for k, (sign) in ((i, -1), (j, 1)):
                        if (it + k) % 2:
                            sign = -sign
                        _name, _t, x0, y0 = label_artists[k]
                        xn = float(np.clip(x0 + sign * 0.38,
                                           xlo + 0.5, xhi - 0.5))
                        if abs(xn - x0) > 1e-9:
                            _t.set_position((xn, y0))
                            label_artists[k] = (_name, _t, xn, y0)
                            moved = True
        if not moved:
            break
        fig.canvas.draw()
    boxes = [t.get_window_extent(renderer) for _, t, _, _ in label_artists]
    collision_pairs = [
        (label_artists[i][0], label_artists[j][0])
        for i in range(len(boxes)) for j in range(i + 1, len(boxes))
        if boxes[i].overlaps(boxes[j])
    ]
    n_collide = len(collision_pairs)
    print(f"Label collisions after de-overlap: {n_collide0} -> {n_collide}")
    if collision_pairs:
        print(f"Remaining label pairs: {collision_pairs}")

    # Legend, in data coordinates below the map. The three structural keys
    # share one row across the full width; fate colors use two balanced rows.
    # This keeps the legend visually attached to the map without leaving an
    # empty right-hand half or spending excessive vertical space.
    structural_y = y_bottom - 0.62
    ax.add_patch(Circle((0.45, structural_y), 0.09, fill=False,
                        edgecolor=LABEL_COLOR, lw=1.1, zorder=5))
    ax.text(0.7, structural_y, "natural assignment exactly on front",
            va="center", fontsize=6.2, color=LABEL_COLOR)
    ax.add_patch(Circle((6.30, structural_y), 0.09, fill=False,
                        edgecolor=LABEL_COLOR, lw=1.1, linestyle="dashed",
                        zorder=5))
    ax.text(6.55, structural_y, "natural assignment off front",
            va="center", fontsize=6.2, color=LABEL_COLOR)
    ax.add_patch(Wedge((11.25, structural_y), 0.09, -110, 90,
                       facecolor=TYPE_COLORS["neuron"],
                       edgecolor="#FFFFFF", lw=0.3, zorder=5))
    ax.add_patch(Wedge((11.25, structural_y), 0.09, 90, 250,
                       facecolor=TYPE_COLORS["muscle"],
                       edgecolor="#FFFFFF", lw=0.3, zorder=5))
    ax.text(11.5, structural_y,
            "pie: terminal-cell fate composition; "
                       "radius ∝ √n",
            va="center", fontsize=6.2, color=LABEL_COLOR)
    types = sorted(TYPE_COLORS)
    for i, tname in enumerate(types):
        col, row = i % 6, i // 6
        row_count = min(6, len(types) - row * 6)
        row_offset = 0.1 + 0.5 * (6 - row_count) * 3.35
        x0 = row_offset + col * 3.35
        y0 = y_bottom - 1.35 - 0.48 * row
        ax.add_patch(Rectangle((x0, y0 - 0.07), 0.14, 0.14,
                               facecolor=TYPE_COLORS[tname],
                               edgecolor=LABEL_COLOR, lw=0.3, zorder=5))
        ax.text(x0 + 0.21, y0, tname.replace("_", " "), va="center",
                fontsize=6.0, color=LABEL_COLOR)

    draw_retention_label_key(ax, 5.0, 10.0, y_bottom - 2.43)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ps.save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote {out_path.with_suffix('.pdf')}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-cells", type=int, default=12,
                        help="Minimum usable terminal cells (default 12).")
    args = parser.parse_args(argv)
    ps.configure()

    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    protein_exp = dl.load_protein_expression()
    v_prot = [n for n in valid_ce if n in protein_exp.index]
    tn, tp = dl.collect_terminals(lineage, v_prot)
    type_map = build_type_map(tn)
    usable = set(tn)

    df = pd.read_csv(ANALYSIS_OUT / f"subtree_summary_min{args.min_cells}.csv")
    qualifying = set(df["subtree"])
    on_front = set(df[df["natural_on_front"].astype(bool)]["subtree"])
    summary_by_name = df.set_index("subtree")

    tree_index = lm.build_lineage_tree_index(lineage)
    terminal_order = terminal_inorder(lineage)
    nodes, edges, _positions, _depth, _max_depth = collect_nodes(
        lineage, tree_index, qualifying, terminal_order, type_map, usable)
    # Attach presentation metrics already validated in the subtree summary.
    for name, nd in nodes.items():
        nd["on_front"] = name in on_front
        nd["n"] = sum(nd["composition"].values())
        nd["max_er"] = float(summary_by_name.loc[name, "max_er"])

    OUT.mkdir(parents=True, exist_ok=True)
    emit_figure(nodes, edges, args.min_cells,
                OUT / "fig4_ce_subtree_map_panel")
    print(f"Subtrees: {len(nodes)}, edges: {len(edges)}, "
          f"on-front: {len(on_front)}")


if __name__ == "__main__":
    main()
