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
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from publication import style as ps
from terminal_pareto.subtree_analysis import (
    build_type_map,
    load_validated_subtree_summary,
)
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
)

OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"

from publication.figures.terminal_subtree_map import (  # noqa: E402,F401  (re-exported drawing API)
    ER_NORM, ER_TEXT_CMAP, INCH_PER_UNIT, LABEL_COLOR, LABEL_SET, TYPE_COLORS, X_SPAN, Y_LEVEL,
    draw_retention_label_key, emit_figure, stagger_crowded_levels, tidy_display_tree,
)


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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-cells", type=int, default=12,
                        help="Minimum usable terminal cells (default 12).")
    parser.add_argument(
        "--profile",
        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"),
        help="Opt into an isolated profile-aware run (omission keeps legacy paths).",
    )
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    ps.configure()

    analysis_out = ANALYSIS_OUT
    out = args.out or OUT
    context = None
    if args.profile is None:
        lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
        _xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
        protein_exp = dl.load_protein_expression()
        v_prot = [n for n in valid_ce if n in protein_exp.index]
        tn, _tp = dl.collect_terminals(lineage, v_prot)
    else:
        context = build_analysis_context(
            args.profile, run_id=args.run_id, output_root=args.output_root)
        lineage = context.lineage
        tn = context.terminal_nodes
        analysis_out = context.run_paths.analysis
        if args.out is None:
            out = context.run_paths.display("endpoint")
    type_map = build_type_map(tn)
    usable = set(tn)

    # Validate the existing run and cache-specific identity before loading the
    # table or allowing context.write() to refresh provenance.
    df = load_validated_subtree_summary(
        analysis_out, args.min_cells, context)
    if context is not None:
        context.write()
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

    out.mkdir(parents=True, exist_ok=True)
    emit_figure(nodes, edges, args.min_cells,
                out / "fig4_ce_subtree_map_panel")
    print(f"Subtrees: {len(nodes)}, edges: {len(edges)}, "
          f"on-front: {len(on_front)}")


if __name__ == "__main__":
    main()
