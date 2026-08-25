"""Cell-type decomposition helpers used by Figures 6 and S2.

The functions in this module attribute natural-edge retention and objective
changes to terminal-cell types along the global Pareto front. They contain no
publication layout code; rendering is handled by the dedicated Figure 6 and
Figure S2 scripts.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.subtree_analysis import first_cousin_null_summary


def decompose_global_front_by_type(lineage, v_prot, tn, tp, xyz_ce, protein_exp,
                                   prot_sel, type_map, tree_index, gp_map,
                                   iteration=300):
    """Per-cell-type edge retention and cost deltas along the global front.

    The preferred cell-type analysis: assignments are those of the single
    global (full-tree) Pareto problem, then decomposed by terminal cell type.
    Every quantity is attributed to the terminal child's type: for child ``k``
    with matched parent slot ``assigned_row[k]`` (``assigned_row[ci] = ri``,
    inverting the Hungarian row/column assignment), the retained-edge test is
    ``parent[assigned_row[k]] == parent[k]``. This matches the global
    biological-edge metric: exchanging two indistinguishable slots belonging
    to the same parent does not remove either parent--child edge. The
    changed-edge tree distance is
    ``d(parent[k], parent[assigned_row[k]])``, and the costs are
    ``cost[assigned_row[k], k]``. Naively testing ``ci[k]`` against child ``k``
    misattributes sibling swaps (114 of 299 parents have two terminal
    children). Returns a tidy DataFrame, one row per (keypoint, type), with
    both total and per-cell cost contributions. Small categories (n <= 4) are
    flagged via ``small``.
    """
    xm, em, _ = pe.build_cost_matrices(tn, tp, xyz_ce, protein_exp, prot_sel)
    # Closed-form first-cousin null moments for the sigma scaling (identical
    # to the summary table); the cousin groups here are the ancestor-based
    # groups used for the first-cousin definition in the publication pipeline.
    rs = first_cousin_null_summary(xm, em, tn, gp_map, seed=42)["_raw"]
    xs = xm / rs["xyz_std"]
    es = em / rs["exp_std"]
    types = [type_map[c] for c in tn]
    idx_by_type = {t: [i for i, ty in enumerate(types) if ty == t]
                   for t in sorted(set(types))}

    # Key points: the two single-objective endpoints and the maximum-edge-
    # retention compromise (the structural compromise used in Figures 2-3).
    # An equal-weight "balanced" point would imply the two objectives are
    # linearly comparable, which the analysis deliberately avoids.
    _xa, _ea, edge_arr, kp = pe.compute_std_scaled_pareto(
        xm, em, tp, rs, iteration=iteration)
    keypoints = [(0, "expr_opt"),
                 (int(kp["max_er_idx"]), "max_er"),
                 (iteration, "spatial_opt")]

    rows = []
    for i, label in keypoints:
        alpha = i / iteration
        ri, ci = pe.linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        assigned_row = np.empty_like(ci)
        assigned_row[ci] = ri
        retained_total = 0
        for t, ix in idx_by_type.items():
            n = len(ix)
            kept = 0
            dsum = 0.0
            tx = te = lx_t = le_t = 0.0
            for k in ix:
                r = int(assigned_row[k])
                tx += xs[r, k]
                te += es[r, k]
                lx_t += xs[k, k]
                le_t += es[k, k]
                if tp[r] == tp[k]:
                    kept += 1
                else:
                    d = lm.lineage_tree_distance(tp[k], tp[r], tree_index)
                    dsum += (d or 0)
            travel = tx - lx_t
            state = te - le_t
            rows.append(dict(
                keypoint=label, type=t, n=n, small=n <= 4, er=kept / n,
                mean_td_changed=dsum / (n - kept) if n - kept else 0.0,
                travel_sigma=travel, state_sigma=state,
                travel_sigma_per_cell=travel / n, state_sigma_per_cell=state / n,
            ))
            retained_total += kept
        if not np.isclose(retained_total / len(tp), edge_arr[i], atol=1e-12):
            raise AssertionError(
                f"Cell-type retention does not recompose global ER at {label}")
    return pd.DataFrame(rows)


def decompose_global_front_retention_by_type(
        lineage, v_prot, tn, tp, xyz_ce, protein_exp, prot_sel, type_map,
        tree_index, gp_map, iteration=300):
    """Cell-type natural-edge retention at every distinct front assignment.

    This is the full-front companion to :func:`decompose_global_front_by_type`.
    It uses the identical cost matrices, first-cousin-null scaling, Hungarian
    assignment direction, and terminal-child type attribution. Consecutive
    sweep solutions with identical objective coordinates are represented once;
    the established maximum-edge-retention keypoint is retained if it belongs
    to such a run. Arc length is measured in the two null-SD objective axes.

    Returns a tidy table with one row per (distinct front assignment, type).
    No interpolation between assignments is performed.
    """
    xm, em, _ = pe.build_cost_matrices(tn, tp, xyz_ce, protein_exp, prot_sel)
    rs = first_cousin_null_summary(xm, em, tn, gp_map, seed=42)["_raw"]
    xs = xm / rs["xyz_std"]
    es = em / rs["exp_std"]
    types = [type_map[c] for c in tn]
    idx_by_type = {t: np.asarray([i for i, ty in enumerate(types) if ty == t])
                   for t in sorted(set(types))}

    xa, ea, edge_arr, kp = pe.compute_std_scaled_pareto(
        xm, em, tp, rs, iteration=iteration)
    max_er_idx = int(kp["max_er_idx"])

    # Collapse only consecutive identical cost points. The assignment at the
    # published max-ER keypoint wins within its run; otherwise use the first
    # deterministic Hungarian solution in the run.
    runs = []
    start = 0
    for i in range(1, iteration + 1):
        if not (np.isclose(xa[i], xa[i - 1], rtol=1e-10, atol=1e-12)
                and np.isclose(ea[i], ea[i - 1], rtol=1e-10, atol=1e-12)):
            runs.append((start, i))
            start = i
    runs.append((start, iteration + 1))
    keep = []
    for lo, hi in runs:
        keep.append(max_er_idx if lo <= max_er_idx < hi else lo)
    keep = np.asarray(keep, dtype=int)

    xk, ek = xa[keep], ea[keep]
    # Match Figure 5's canonical coordinate exactly: normalize each
    # objective by its attainable endpoint range before measuring arc length.
    # The sweep here is state -> travel, so X falls 1 -> 0 and Y rises 0 -> 1.
    denx = float(xk[0] - xk[-1])
    deny = float(ek[-1] - ek[0])
    if not (denx > 0 and deny > 0):
        raise ValueError("Degenerate global-front endpoint range")
    x_norm = (xk - xk[-1]) / denx
    e_norm = (ek - ek[0]) / deny
    steps = np.hypot(np.diff(x_norm), np.diff(e_norm))
    cumulative = np.r_[0.0, np.cumsum(steps)]
    u = cumulative / cumulative[-1] if cumulative[-1] > 0 else cumulative

    parent_array = np.asarray(tp)
    rows = []
    for front_index, (sweep_index, ui) in enumerate(zip(keep, u)):
        alpha = sweep_index / iteration
        ri, ci = pe.linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        assigned_row = np.empty_like(ci)
        assigned_row[ci] = ri
        retained_total = 0
        for cell_type, indices in idx_by_type.items():
            kept_edges = int(np.count_nonzero(
                parent_array[assigned_row[indices]] == parent_array[indices]))
            n = int(len(indices))
            rows.append(dict(
                front_index=front_index,
                sweep_index=int(sweep_index),
                alpha=float(alpha),
                # The raw sweep begins at the cell-state optimum. Figure 5
                # uses travel optimum = 0, so reverse u for visual consistency.
                u=float(1.0 - ui),
                travel_sigma=float(xa[sweep_index]),
                state_sigma=float(ea[sweep_index]),
                is_max_er=bool(sweep_index == max_er_idx),
                type=cell_type,
                n=n,
                small=bool(n <= 4),
                er=float(kept_edges / n),
            ))
            retained_total += kept_edges
        if not np.isclose(retained_total / len(tp), edge_arr[sweep_index],
                          atol=1e-12):
            raise AssertionError(
                "Cell-type retention does not recompose global ER at "
                f"sweep index {sweep_index}")
    return pd.DataFrame(rows).sort_values(["u", "type"]).reset_index(drop=True)


SMALL_TYPES = {"mesoderm", "reproduction", "excretory", "other"}


def _merge_small(ct_df):
    """Merge n <= 4 categories into a single "other (n<=4)" group."""
    out = ct_df.copy()
    small = out["type"].isin(SMALL_TYPES) | out["small"]
    small_rows = out[small]
    merged = []
    for kp in sorted(set(small_rows["keypoint"])):
        g = small_rows[small_rows["keypoint"] == kp]
        n = int(g["n"].sum())
        merged.append({
            "keypoint": kp, "type": "other (n<=4)", "n": n, "small": True,
            "er": float(np.average(g["er"], weights=g["n"])),
            "mean_td_changed": np.nan,
            "travel_sigma": float(g["travel_sigma"].sum()),
            "state_sigma": float(g["state_sigma"].sum()),
            "travel_sigma_per_cell": float(g["travel_sigma"].sum() / n),
            "state_sigma_per_cell": float(g["state_sigma"].sum() / n),
        })
    return pd.concat([out[~small], pd.DataFrame(merged)], ignore_index=True)
