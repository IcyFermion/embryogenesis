"""Read-only inputs for the revised manuscript layout (no back-end data changes).

Every loader goes through an existing validated adapter or add-on cache:
the primary terminal analysis (plus its cached per-weight assignments), the
per-round references of the pooled full tree, the Track E per-layer runs, and
the 187-edge cross-species terminal comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


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


def briggsae_layers():
    """C. briggsae partial forest by layer (Figure 8C): distance to the front and fewest edges rewired.

    Read from the validated cross-species round references (``layer_hourglass.load_partial`` validates the run)
    and the hash-checked Track E stage table; birth sizes are those of the matched C. elegans cells (Large et al.
    cell table). Degenerate layers (layer 3: 2 cells) are left out.
    """
    from full_tree_pareto import layer_hourglass as lh
    from full_tree_pareto import layer_stages as ls
    from publication.adapters import full_tree as adapter
    stages, _ = adapter._checked_layer_run(ls, adapter.LAYER_RUNS["stages"])
    stage_table = pd.read_csv(stages / "layers.csv")
    refs = lh.load_partial()["references"]
    refs = refs[(refs.config == "cb_rna") & ~refs.degenerate.astype(bool)].assign(layer=lambda f: 9 - f["round"])
    table = refs.merge(stage_table[stage_table.tree == "partial"][["layer", "alive_median"]], on="layer")
    table = table[table.layer >= adapter.FIRST_LAYER].sort_values("layer")
    return pd.DataFrame(dict(layer=table.layer, edges=table.edges, born=table.alive_median,
                             front_distance=table.d_LP, fewest_rewired=1 - table.max_edge_retention,
                             front_points=table.distinct_front_points)).reset_index(drop=True)


FEATURE_RUN = "criteria_metric_20261014"
FEATURE_CANDIDATES = ("top-20", "HV pool", "whole panel")


def feature_selection(run=FEATURE_RUN, n_random=200, seed=20261108):
    """Inputs of the cell-state feature figure (Figure 2): C. elegans protein, cosine, 275-edge cohort.

    * ``ranking``: the published classifier ranking of all 210 reporters (L2 multinomial logistic regression on
      Schema D classes; importance = largest coefficient over classes, and the class it belongs to).
    * ``sisters``: cosine distance of every typed terminal sister pair, same versus different fine type, for
      each candidate; recomputed and checked against the recorded Cohen's d of ``run`` (1e-12).
    * ``separation``: recorded Cohen's d with 95% bootstrap intervals; ``random``: recorded d of the 200
      random 20-reporter sets from the HV pool.
    * ``by_path``: mean cosine distance of all cohort pairs by lineage path length (lengths pooled in twos,
      2-3, 4-5, ...: odd lengths join cells of different depth), for the top-20, the whole panel and
      ``n_random`` random 20-reporter sets from the HV pool, as the recorded random sets (drawn here, seeded).
    * ``memory``: recorded same-type pairs by path length against the within-type permutation null (Methods only;
      the figure panel was dropped on 2026-11-09).
    """
    from scipy.spatial.distance import pdist, squareform
    from terminal_pareto import representation_criteria as rc
    analysis = rc.run_path(run) / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    for name, expected in record["files"].items():
        if rc.file_hash(analysis / name) != expected:
            raise ValueError(f"Track A cache hash mismatch in {run}: {name}")
    ctx = rc.load_context()
    design = ctx.designs["CE protein"]
    criteria = pd.read_csv(analysis / "criteria.csv")
    criteria = criteria[(criteria.dataset == "CE protein") & (criteria.metric == "cosine")].set_index("candidate")
    same, diff = design.sisters["fine"]
    sisters, separation, matrices = [], [], {}
    for candidate in FEATURE_CANDIDATES:
        X = rc.matrix(ctx, "CE protein", ctx.candidates["protein"][candidate])
        D = squareform(pdist(X, "cosine"))
        matrices[candidate] = D
        d = rc.cohen_d(D, same, diff)
        recorded = criteria.loc[candidate]
        if abs(d - recorded.sib_fine) > 1e-12:
            raise ValueError(f"Sibling separation of {candidate} does not reproduce {run}: {d} vs {recorded.sib_fine}")
        for kind, index in (("same", same), ("different", diff)):
            sisters += [dict(candidate=candidate, kind=kind, distance=v) for v in D[index[:, 0], index[:, 1]]]
        separation.append(dict(candidate=candidate, d=d, low=recorded.sib_fine_lo, high=recorded.sib_fine_hi,
                               n_same=int(recorded.n_sisters_fine_same), n_different=int(recorded.n_sisters_fine_diff),
                               memory_rho=recorded.memory_rho, memory_null_mean=recorded.memory_null_mean,
                               memory_null_lo=recorded.memory_null_lo, memory_null_hi=recorded.memory_null_hi))
    random = pd.read_csv(analysis / "random_metrics.csv")
    random = random[(random.dataset == "CE protein") & (random.metric == "cosine")]
    memory = pd.read_csv(analysis / "memory_curves.csv")
    memory = memory[(memory.dataset == "CE protein") & (memory.metric == "cosine")]

    # All cohort pairs by lineage path length, pooled in twos (descriptive; random sets from the HV pool).
    cells = design.cells
    i, j = np.triu_indices(len(cells), 1)
    path = np.array([ctx.lineage.path_length(cells[a], cells[b]) for a, b in zip(i, j)])
    path = path - path % 2
    rng = np.random.default_rng(seed)
    pool = list(ctx.candidates["protein"]["HV pool"])
    draws = {}
    for k in range(n_random):
        X = rc.matrix(ctx, "CE protein", rng.choice(pool, size=20, replace=False))
        draws[k] = squareform(pdist(X, "cosine"))[i, j]
    rows = []
    for length in np.unique(path):
        mask = path == length
        row = dict(path_low=int(length), path_high=int(length) + 1, n_pairs=int(mask.sum()))
        for candidate in ("top-20", "whole panel"):
            row[candidate] = float(matrices[candidate][i[mask], j[mask]].mean())
        values = np.array([np.nanmean(v[mask]) for v in draws.values()])
        row.update(random_median=float(np.median(values)), random_low=float(np.quantile(values, .05)),
                   random_high=float(np.quantile(values, .95)))
        rows.append(row)
    ranking = pd.read_csv(ROOT / "expression_embedding/results/elegans_protein_linear_baseline/protein_feature_ranking.csv")
    if ranking.protein.head(20).tolist() != list(ctx.candidates["protein"]["top-20"]):
        raise ValueError("Ranking file and published top-20 differ")
    return dict(ranking=ranking, sisters=pd.DataFrame(sisters), separation=pd.DataFrame(separation),
                random=random.sib_fine.to_numpy(), by_path=pd.DataFrame(rows), memory=memory, run=run,
                n_cells=len(cells), analysis_id=record.get("analysis_id"))


def tracking_front_distances(t: Terminal) -> dict:
    """d_LP of the natural lineage to its own front in each tracking geometry of released S3 (main-text note).

    Each geometry's front is its own 301-weight sweep (source == target in the S3 projections), in that geometry's
    endpoint coordinates; the pooled value equals the primary global d_LP.
    """
    projections = t.data.s3_projections
    distances = {}
    for target, reference in t.data.s3_references.items():
        L = np.asarray(reference["natural_display"], float)
        front = projections[(projections.target == target) & (projections.source == target)]
        distances[target] = float(np.hypot(front.x - L[0], front.y - L[1]).min())
    primary = float(t.data.canonical_display.set_index("subtree").loc["P0", "d_lp"])
    if abs(distances["pooled"] - primary) > 1e-9:
        raise ValueError(f"Pooled S3 d_LP {distances['pooled']} differs from the primary {primary}")
    return distances
