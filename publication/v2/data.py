"""Read-only inputs for the revised manuscript layout (no back-end data changes).

Every loader goes through an existing validated adapter or add-on cache:
the primary terminal analysis (plus its cached per-weight assignments), the
per-round references of the pooled full tree, the Track E per-layer runs, and
the 187-edge cross-species terminal comparison.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Terminal:
    """Primary 275-edge terminal cohort with the cached global assignments."""

    data: object  # publication.contracts.TerminalPrimaryInput
    nodes: list  # terminal cells, cost-matrix order
    parents: np.ndarray  # natural parent of each terminal
    assignments: np.ndarray  # (weights, terminals): slot assigned to terminal j; parent = parents[slot]
    tree_index: dict
    lineage: dict
    retention: np.ndarray  # replayed from assignments; equals the cached front record

    @property
    def travel_optimum(self):
        return len(self.assignments) - 1

    @property
    def state_optimum(self):
        return 0

    @property
    def maximum(self):
        return int(self.data.global_twr["kp"]["max_er_idx"])


def terminal() -> Terminal:
    from publication.adapters import terminal as adapter
    from terminal_pareto.analysis_context import DEFAULT_OUTPUT_ROOT, build_analysis_context
    from terminal_pareto.global_analysis import load_global_analysis
    data = adapter.primary()
    context = build_analysis_context(adapter.PRIMARY_PROFILE, run_id=adapter.PRIMARY_RUN,
                                     output_root=DEFAULT_OUTPUT_ROOT, sweep_intervals=adapter.PRIMARY_INTERVALS)
    result = load_global_analysis(context)
    parents = np.asarray(context.terminal_parents)
    assignments = np.asarray(result.assignments)
    retention = (parents[assignments] == parents[None, :]).mean(axis=1)
    np.testing.assert_allclose(retention, np.asarray(data.global_twr["traditional_er"]), atol=1e-12)
    return Terminal(data=data, nodes=list(context.terminal_nodes), parents=parents, assignments=assignments,
                    tree_index=context.tree_index, lineage=context.lineage, retention=retention)


def tree_distance(tree_index, a, b):
    """Edges on the lineage path between cells ``a`` and ``b``."""
    if a == b:
        return 0
    up_a = [a] + tree_index[a]["ancestors"]
    up_b = [b] + tree_index[b]["ancestors"]
    common = set(up_a) & set(up_b)
    return min(up_a.index(c) for c in common) + min(up_b.index(c) for c in common)


def edge_changes(t: Terminal, weight: int) -> pd.DataFrame:
    """Terminals whose assigned parent differs from their natural parent at one sweep weight."""
    assigned = t.parents[t.assignments[weight]]
    rows = []
    for j, cell in enumerate(t.nodes):
        if assigned[j] != t.parents[j]:
            rows.append(dict(child=cell, natural_parent=t.parents[j], assigned_parent=assigned[j],
                             tree_distance=tree_distance(t.tree_index, t.parents[j], assigned[j])))
    return pd.DataFrame(rows, columns=["child", "natural_parent", "assigned_parent", "tree_distance"])


def layers():
    """Pooled full tree by layer (layer 1 = top-most round): per-round inputs for Figure 6A."""
    from publication.sketches.data import pooled_rounds
    rounds, _, _ = pooled_rounds()
    last = max(rounds)
    return {last + 1 - r: d for r, d in rounds.items()}


def layer_summary():
    """Track E per-layer summaries for Figure 6B-C (layers 3-8), through the validated adapter."""
    from publication.adapters import full_tree as adapter
    return adapter.layers()


def briggsae():
    """187-edge cross-species terminal comparison (C. briggsae RNA is read from it)."""
    from publication.adapters import terminal as adapter
    return adapter.cross_species()
