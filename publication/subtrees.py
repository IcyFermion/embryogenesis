"""Subtree helpers shared by the adapter, the manuscript set and the figure families.

Moved on 2026-11-13 from ``publication.figures.terminal_subtree_map`` (lineage order and node layout) and
``publication.figures.terminal_canonical`` (the caption note on ABplp and its build-time check), so that the
manuscript set (``publication.v2``) imports nothing from the figure families. Those modules re-export the names.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

X_SPAN = 20.0
Y_LEVEL = 1.0   # vertical distance per depth level


def terminal_inorder(root, name_of):
    """Assign each terminal cell an inorder index 0..n-1 (left-to-right)."""
    order = []

    def dfs(node):
        children = node.get("children", [])
        if not children:
            order.append(name_of(node["did"]))
            return
        for c in children:
            dfs(c)

    dfs(root)
    return order


def collect_nodes(root, tree_index, qualifying, terminal_order, type_map,
                  usable, name_of):
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
        name = name_of(node.get("did", ""))
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
            cname = name_of(c.get("did", ""))
            xs.append(positions[cname][0])
        # x of an internal node = midpoint of its children's xs
        positions[name] = (float(np.mean(xs)) if xs else 0.0, 0.0)
        comp[name] = Counter()
        for c in children:
            cname = name_of(c.get("did", ""))
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


# Caption note on ABplp (Track D4c): the stored results that support it, checked when a figure is built.
NOTES = {"ABplp": ("†", "d_lp")}


def check_notes(evidence):
    """Raise if the stored Track D results no longer support the caption note on ABplp."""
    claims = {
        "ABplp lies farthest from its front": evidence["farthest"] == "ABplp",
        "ABplp has three apoptotic cells": len(evidence["deaths"]) == 3,
        "only ABplpppapp's birth endpoint puts ABplp on its front":
            {d for d, on in evidence["alone_on_front"].items() if on} == {"ABplpppapp"},
        "matched living cells: 0 of 1000 draws": (evidence["control_draws"], evidence["control_on_front"]) == (1000, 0),
    }
    failed = [claim for claim, ok in claims.items() if not ok]
    if failed:
        raise ValueError(f"Figure 4 caption note (released Figure 5B) no longer holds: {failed}")
