"""Exploratory subtree-level Pareto analysis for publication Figure 4.

Computes a tidy per-subtree summary table for the primary publication
configuration (C. elegans embryo 1, 3D tracking, T <= 255, top-20 protein
features, terminal-only optimization) and writes it under
``terminal_pareto/output/``.

This is an exploratory analysis script, not a publication renderer. It reuses
the exact preprocessing, null construction, std-scaling, and assignment logic
of the completed terminal pipeline (``fig2_fig3_ce_terminal_pareto.py``).

Usage::

    python terminal_pareto/subtree_analysis.py [--min-cells 12]

The ``--min-cells`` threshold controls which subtrees are included; the handoff
directs starting around 10-15 usable terminal cells and testing sensitivity.
"""

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe


OUT = Path(__file__).resolve().parent / "output"

# Analysis parameters (match the completed terminal pipeline).
ITERATION = 300          # Pareto sweep resolution
N_RANDOM_NULL = 1000     # first-cousin null draws (used only for the null cloud)
N_RANDOM_NLAD = 100      # NLAD null draws per subtree
EXACT_NULL_CAP = 200_000   # cousin-shuffle space size below which optimality is enumerated
OPTIMALITY_N_DRAWS = 1000  # Monte Carlo draws for non-enumerable optimality nulls

# Cell-type taxonomy: wormweb lineage -> merged developmental category.
# Join on ``wormweb.lineage``, the convention used across the repo.
_MERGE = pe.MERGE_MAP


def build_type_map(terminal_nodes):
    """Map terminal lineage names to merged cell types via the entropy key.

    Joins on ``wormweb.lineage`` (the convention used in lineage_complexity.ipynb).
    Rows without a type are dropped; duplicated terminals carry identical types
    (verified for the primary terminal set). Terminal cells missing from the key
    are the programmed-death cells: the 32 untyped terminals of the primary
    configuration are exactly the entries of ``data/apoptotic_cells.txt``, and
    ``pareto_engine.decompose_by_cell_type`` already defaults missing types to
    ``programmed_death``. We make that explicit and cross-check it.

    The ``tail`` annotation is kept as its own category: the two tail-typed
    terminals (ABplppppppa, ABprppppppa) are NOT in the apoptotic list, so
    merging them into ``programmed_death`` would mix annotation sources. They
    are n=2 and are merged into the displayed "other (n<=4)" group in figures.
    Returns {lineage: merged_type}.
    """
    ct = dl.load_cell_type_map()
    rows = ct[ct["wormweb.lineage"].astype(str).isin(set(terminal_nodes))]
    rows = rows.dropna(subset=["wormweb.type"])
    out = {}
    for _, row in rows.iterrows():
        name = str(row["wormweb.lineage"])
        raw = str(row["wormweb.type"]).strip()
        out[name] = "tail" if raw == "tail" else _MERGE.get(raw, "other")
    # Programmed-death cells absent from the key: verify against the apoptotic
    # list and record them explicitly.
    apoptotic = set(dl.load_apoptotic())
    missing = [c for c in terminal_nodes if c not in out]
    unverified = [c for c in missing if c not in apoptotic]
    if unverified:
        print(f"Warning: {len(unverified)} terminal cells lack a type and are "
              f"not in the apoptotic list: {unverified[:8]}")
    for c in missing:
        out[c] = "programmed_death"
    return out


def exact_cousin_stats(xm, em, groups):
    """Exact mean/var/cov of the cousin-shuffle total costs.

    For a group ``g`` of size ``k`` under a uniform random permutation ``pi``,
    ``E[x[i, pi(i)]] = rowsum_i / k`` and, for ``i != j``,
    ``E[x[i, pi(i)] x[j, pi(j)]] = (rowsum_i rowsum_j - sum_a x[i,a] x[j,a]) /
    (k (k-1))``. Groups are independent, so totals add per group.
    Returns (mean_x, mean_e, var_x, var_e, cov_xe).
    """
    mean_x = mean_e = var_x = var_e = cov_xe = 0.0
    grouped = set()
    for g in groups:
        idx = np.asarray(g)
        grouped.update(idx.tolist())
        k = len(idx)
        xg = xm[np.ix_(idx, idx)]
        eg = em[np.ix_(idx, idx)]
        rx, re = xg.sum(axis=1), eg.sum(axis=1)
        mx, me = rx.sum() / k, re.sum() / k
        # Sums over all ordered pairs (i, j), including i == j.
        sxx = (xg ** 2).sum()          # sum_{i,a} x[i,a]^2
        see = (eg ** 2).sum()
        sxe = (xg * eg).sum()
        cx = xg.sum(axis=0)
        ce = eg.sum(axis=0)
        # Sum over i != j of sum_a x[i,a] x[j,a] (columns of the same matrix).
        xx_off = (cx ** 2).sum() - sxx
        ee_off = (ce ** 2).sum() - see
        xe_off = (cx * ce).sum() - sxe
        rx2 = (rx ** 2).sum()
        re2 = (re ** 2).sum()
        rxe = (rx * re).sum()
        rx_sum, re_sum = rx.sum(), re.sum()
        denom = k * (k - 1)
        # E[S_x^2] = sum_i E[x^2] + sum_{i!=j} E[x_i x_j].
        ex2 = sxx / k + (rx_sum ** 2 - rx2 - xx_off) / denom
        ee2 = see / k + (re_sum ** 2 - re2 - ee_off) / denom
        exe = sxe / k + (rx_sum * re_sum - rxe - xe_off) / denom
        var_x += ex2 - mx ** 2
        var_e += ee2 - me ** 2
        cov_xe += exe - mx * me
        mean_x += mx
        mean_e += me
    # Ungrouped cells are fixed at their diagonal (identity) cost.
    n = xm.shape[0]
    ungrouped = [i for i in range(n) if i not in grouped]
    if ungrouped:
        mean_x += xm[ungrouped, ungrouped].sum()
        mean_e += em[ungrouped, ungrouped].sum()
    return mean_x, mean_e, var_x, var_e, cov_xe


def first_cousin_null_summary(xm, em, tn, gp_map, n_random=N_RANDOM_NULL, seed=42):
    """First-cousin shuffle null in raw units plus sigma (std) scales.

    The null mean/std/covariance are computed in closed form (``exact_cousin_stats``)
    for every subtree — no Monte Carlo, so no sampling error in the sigma
    scaling. The random-sample arrays are still drawn for the null cloud only.
    Raises if a null standard deviation is not positive, before any division.
    """
    groups = pe.build_cousin_groups(tn, gp_map)
    n_cells = len(tn)
    covered = sum(len(g) for g in groups)
    rs = pe.compute_cousin_random_stats(xm, em, groups, n_random=n_random, seed=seed)
    mx, me_, vx, ve, _cov = exact_cousin_stats(xm, em, groups)
    sx, se = float(np.sqrt(vx)), float(np.sqrt(ve))
    if not (sx > 0 and se > 0):
        raise ValueError(
            f"Null standard deviation is not positive (xyz_std={sx}, "
            f"exp_std={se}); the cousin groups do not randomize these costs. "
            f"Refusing to standardize.")
    rs = dict(rs, xyz_mean=mx, exp_mean=me_, xyz_std=sx, exp_std=se)
    return dict(
        n_cousin_groups=len(groups),
        cousin_covered_frac=covered / n_cells if n_cells else 0.0,
        null_xyz_std=sx,
        null_exp_std=se,
        # Null mean position relative to the natural lineage, in sigma units.
        null_xyz_mean_sigma=(mx - rs["lineage_xyz"]) / sx,
        null_exp_mean_sigma=(me_ - rs["lineage_exp"]) / se,
        _raw=rs,
    )


def analyze_subtree(name, terms, xyz_ce, protein_exp, prot_sel, tree_index,
                    gp_map, type_map, seed=42):
    """Compute the full metric set for one subtree.

    Returns a flat dict of scalars ready for the summary table.
    """
    tn = [c for c, _p in terms]
    tp = [_p for c, _p in terms]
    n = len(tn)

    # ── Structural descriptors ──
    root_depth = tree_index[name]["depth"]
    rel_depths = [tree_index[c]["depth"] - root_depth for c in tn]
    types = [type_map.get(c) for c in tn]
    typed = [t for t in types if t is not None]
    n_types = len(set(typed)) if typed else 0
    entropy = max(pe.type_shannon_entropy(tn, lambda c: type_map.get(c)), 0.0) if typed else np.nan

    # ── Cost matrices and null ──
    xm, em, _ = pe.build_cost_matrices(tn, tp, xyz_ce, protein_exp, prot_sel)
    lineage_xyz_raw = float(xm.diagonal().sum())
    lineage_exp_raw = float(em.diagonal().sum())
    groups = pe.build_cousin_groups(tn, gp_map)
    null = first_cousin_null_summary(xm, em, tn, gp_map, seed=seed)

    # ── Tree-aware Pareto metrics (one call, as in the publication pipeline) ──
    twr = lm.combined_lineage_proximity(
        xm, em, tp, tn, tree_index, null["_raw"],
        iteration=ITERATION, n_random=N_RANDOM_NLAD,
    )
    kp = twr["kp"]
    edge = np.asarray(twr["traditional_er"])
    ct = twr["cost_tree_tradeoff"]
    lmd = np.asarray(twr["lineage_mean_dist"])
    nlad = twr["nlad"]["nlad"]

    # Size-matched optimality baseline: how often are random assignments
    # exactly on the front at this subtree size?
    opt_null = exact_optimality_null(
        xm, em, tp, null["_raw"], groups, twr["xyz_arr"], twr["exp_arr"])

    # Single-objective attainable reductions (positive = saving vs lineage, sigma).
    attain_xyz = float(ct["delta_xyz"][-1])     # toward travel optimum (alpha=1)
    attain_exp = float(ct["delta_exp"][0])      # toward cell-state optimum (alpha=0)

    # Structural change toward each single-objective optimum.
    tree_dist_expr_opt = float(lmd[0])
    tree_dist_spatial_opt = float(lmd[-1])
    nlad_expr_opt = float(nlad[0])
    nlad_spatial_opt = float(nlad[-1])

    # Relative distance from the natural lineage to the Pareto front, with the
    # null mean as the reference point (notebook convention). Also record the
    # two raw components: LP (front proximity, sigma) and NP (distance from the
    # null mean to the same closest front point, sigma).
    dists = np.sqrt(twr["xyz_arr"] ** 2 + twr["exp_arr"] ** 2)
    p_idx = int(np.argmin(dists))
    lp = float(dists[p_idx])
    np_ = float(np.sqrt((null["null_xyz_mean_sigma"] - twr["xyz_arr"][p_idx]) ** 2
                        + (null["null_exp_mean_sigma"] - twr["exp_arr"][p_idx]) ** 2))
    rel_pareto_dist = lp / np_ if np_ > 1e-10 else np.inf

    return dict(
        subtree=name,
        n=n,
        n_typed=len(typed),
        typed_frac=len(typed) / n if n else 0.0,
        n_types=n_types,
        entropy_nats=entropy,
        depth_min=min(rel_depths),
        depth_med=float(np.median(rel_depths)),
        depth_max=max(rel_depths),
        lineage_xyz_raw=lineage_xyz_raw,
        lineage_exp_raw=lineage_exp_raw,
        expr_opt_er=float(edge[0]),
        spatial_opt_er=float(edge[-1]),
        max_er=float(edge[kp["max_er_idx"]]),
        attain_xyz_sigma=attain_xyz,
        attain_exp_sigma=attain_exp,
        tree_dist_expr_opt=tree_dist_expr_opt,
        tree_dist_spatial_opt=tree_dist_spatial_opt,
        nlad_expr_opt=nlad_expr_opt,
        nlad_spatial_opt=nlad_spatial_opt,
        front_proximity_sigma=lp,
        null_to_front_sigma=np_,
        rel_pareto_dist=rel_pareto_dist,
        **opt_null,
        **{k: v for k, v in null.items() if not k.startswith("_")},
    )


def analyze_all(min_cells, lineage, xyz_ce, protein_exp, prot_sel, v_prot,
                tree_index, gp_map, type_map):
    """Run the metric set over every subtree with >= min_cells usable terminals."""
    subtrees = dl.collect_all_subtrees(lineage, v_prot, min_cells=min_cells)
    rows = []
    for name, terms in subtrees:
        row = analyze_subtree(name, terms, xyz_ce, protein_exp, prot_sel,
                              tree_index, gp_map, type_map)
        rows.append(row)
        print(f"  {name:10s} n={row['n']:3d}  types={row['n_types']:2d}  "
              f"entropy={row['entropy_nats']:.2f}  max_er={row['max_er']:.3f}  "
              f"rel_dist={row['rel_pareto_dist']:.3f}  "
              f"depth_med={row['depth_med']:.1f}")
    return pd.DataFrame(rows)


def _cousin_perm(n, groups, rng):
    """Draw one cousin-shuffled permutation (identity on ungrouped cells)."""
    perm = np.arange(n)
    for g in groups:
        if len(g) < 2:
            continue
        sg = np.asarray(g).copy()
        rng.shuffle(sg)
        perm[g] = sg
    return perm


def _wilson_ci(k, n, z=1.96):
    """Wilson 95% score interval for k successes in n Bernoulli draws."""
    if n <= 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / denom
    return (float(centre - half), float(centre + half))


def _on_front_flags(x_costs, e_costs, alphas, m_cost):
    """Boolean per assignment: optimal for some alpha on the grid."""
    # x_costs, e_costs: (n_assignments,) arrays of total costs (sigma units).
    c = (x_costs[:, None] * alphas[None, :]
         + e_costs[:, None] * (1 - alphas)[None, :])
    return np.any(np.isclose(c, m_cost[None, :], rtol=1e-7, atol=1e-9), axis=1)


def exact_optimality_null(xm, em, tp, rs, groups, xyz_arr, exp_arr,
                          iteration=ITERATION, n_draws=OPTIMALITY_N_DRAWS,
                          seed=7):
    """Fraction of null assignments that are exactly Pareto-optimal.

    The natural lineage's "exactly on the front" status (relative distance 0)
    means its cost point is optimal for some alpha on the sweep grid. This
    computes the same property for cousin-shuffled and fully random
    assignments, calibrating how often exact optimality occurs by chance at
    each subtree size.

    The cousin-shuffle space is enumerated exactly when small
    (``prod(|g|!) <= EXACT_NULL_CAP``); otherwise (and always for the full
    random space) ``OPTIMALITY_N_DRAWS`` permutations are sampled and a
    Wilson 95% interval is reported.

    Returns a dict with natural_on_front, cousin (frac, ci, method) and
    full_random (frac, ci, method).
    """
    xs = xm / rs["xyz_std"]
    es = em / rs["exp_std"]
    n = len(tp)
    lin_x = np.diag(xs).sum()
    lin_e = np.diag(es).sum()
    alphas = np.arange(iteration + 1) / iteration
    # Minimum achievable weighted cost at each alpha (reconstructed from the
    # sweep: xyz_arr = m_xyz(alpha) - lin_x, exp_arr = m_exp(alpha) - lin_e).
    m_cost = alphas * (np.asarray(xyz_arr) + lin_x) \
        + (1 - alphas) * (np.asarray(exp_arr) + lin_e)
    rows = np.arange(n)
    rng = np.random.default_rng(seed)

    def _mc_frac(perm_matrix):
        x_costs = xs[rows[None, :], perm_matrix].sum(axis=1)
        e_costs = es[rows[None, :], perm_matrix].sum(axis=1)
        flags = _on_front_flags(x_costs, e_costs, alphas, m_cost)
        k = int(flags.sum())
        lo, hi = _wilson_ci(k, len(perm_matrix))
        return k / len(perm_matrix), lo, hi, "mc"

    def _exact_frac():
        # Per-group permutation cost tables (sigma units, matching m_cost),
        # then product over groups.
        tables = []
        for g in groups:
            idx = np.asarray(g)
            perms = np.array(list(itertools.permutations(idx)))
            x_costs = xs[idx[None, :], perms].sum(axis=1)
            e_costs = es[idx[None, :], perms].sum(axis=1)
            tables.append((x_costs, e_costs))
        # Ungrouped cells are fixed at their diagonal cost.
        grouped = set().union(*[set(g) for g in groups]) if groups else set()
        ungrouped = [i for i in range(n) if i not in grouped]
        base_x = xs[ungrouped, ungrouped].sum()
        base_e = es[ungrouped, ungrouped].sum()
        space = math.prod(len(t[0]) for t in tables)
        x_arr = np.empty(space)
        e_arr = np.empty(space)
        combos = itertools.product(*[range(len(t[0])) for t in tables])
        for cidx, combo in enumerate(combos):
            x_arr[cidx] = base_x + sum(t[0][j] for t, j in zip(tables, combo))
            e_arr[cidx] = base_e + sum(t[1][j] for t, j in zip(tables, combo))
        flags = _on_front_flags(x_arr, e_arr, alphas, m_cost)
        return flags.mean(), np.nan, np.nan, "exact"

    # Cousin-shuffled assignments (same groups as the null scaling).
    space = math.prod(math.factorial(len(g)) for g in groups)
    if 0 < space <= EXACT_NULL_CAP:
        cousin = _exact_frac()
    else:
        cousin_perms = np.stack([
            _cousin_perm(n, groups, rng) for _ in range(n_draws)])
        cousin = _mc_frac(cousin_perms)
    full_perms = np.stack([rng.permutation(n) for _ in range(n_draws)])
    full = _mc_frac(full_perms)

    # Natural lineage = the identity assignment.
    cx = np.diag(xs).sum()
    ce = np.diag(es).sum()
    natural_on_front = bool(np.any(
        np.isclose(alphas * cx + (1 - alphas) * ce, m_cost,
                   rtol=1e-7, atol=1e-9)))
    return dict(
        natural_on_front=natural_on_front,
        cousin_on_front_frac=cousin[0],
        cousin_on_front_ci_low=cousin[1],
        cousin_on_front_ci_high=cousin[2],
        cousin_on_front_method=cousin[3],
        full_random_on_front_frac=full[0],
        full_random_on_front_ci_low=full[1],
        full_random_on_front_ci_high=full[2],
        full_random_on_front_method=full[3],
    )


def nested_null_sharing(subtrees, gp_map):
    """Count subtrees whose cousin-group cell sets coincide.

    When a nested subtree's extra terminals are ungrouped singletons, its
    first-cousin groups are identical to the parent subtree's, so the null
    sigma and centered null mean are literally shared. This is a concrete
    instance of the statistical non-independence of nested subtrees.
    """
    by_groups = {}
    for name, terms in subtrees:
        tn = [c for c, _p in terms]
        groups = pe.build_cousin_groups(tn, gp_map)
        key = frozenset(frozenset(tn[i] for i in g) for g in groups)
        by_groups.setdefault(key, []).append(name)
    return {k: v for k, v in by_groups.items() if len(v) > 1}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-cells", type=int, default=12,
                        help="Minimum usable terminal cells per subtree (default 12).")
    args = parser.parse_args(argv)

    # Primary publication configuration: C. elegans protein, T <= 255, top-20.
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    protein_exp = dl.load_protein_expression()
    prot_sel = dl.load_prot_sel()
    v_prot = [n for n in valid_ce if n in protein_exp.index]
    gp_map = pe.build_grandparent_map(lineage)
    tn_full, _ = dl.collect_terminals(lineage, v_prot)
    type_map = build_type_map(tn_full)
    print(f"Type map: {sum(1 for t in type_map.values() if t is not None)}/"
          f"{len(tn_full)} terminals typed, merged categories "
          f"{sorted(set(type_map.values()))}")

    df = analyze_all(args.min_cells, lineage, xyz_ce, protein_exp, prot_sel,
                     v_prot, tree_index, gp_map, type_map)

    # Nested non-independence: how many subtrees share an identical cousin-group
    # structure (and hence an identical null scale)?
    subtrees_all = dl.collect_all_subtrees(lineage, v_prot,
                                           min_cells=args.min_cells)
    sharing = nested_null_sharing(subtrees_all, gp_map)
    n_shared = sum(len(v) for v in sharing.values())
    print(f"\nNested-null sharing: {len(sharing)} group structures shared by "
          f"{n_shared} subtrees "
          f"(e.g., {[v for v in sharing.values() if len(v) > 1][:3]}); "
          f"their null sigma is identical, so sigma units are not independent "
          f"across nested subtrees.")

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"subtree_summary_min{args.min_cells}.csv"
    df.to_csv(path, index=False)
    print(f"\nWrote {len(df)} rows to {path}")
    return df


if __name__ == "__main__":
    main()
