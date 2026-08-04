"""Pareto optimization engine for terminal-only lineage analysis.

All functions are pure: input arrays/dicts, output arrays/dicts.
No file I/O, no matplotlib — computation only.
"""

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cosine
from collections import defaultdict
from joblib import Parallel, delayed
from tqdm import tqdm

# Import shared utilities from data_loader (no circular dependency:
# data_loader only has I/O functions, pareto_engine has computation)
from terminal_pareto.data_loader import map_names, find_node, collect_terminals

# ── Shared constants ──
SUBTREES = ['AB', 'ABa', 'ABp', 'P1']
MARKERS = {'expr_opt_idx': 's', 'spatial_opt_idx': 'o', 'max_er_idx': 'D'}
MERGE_MAP = {
    'neuron': 'neuron', 'muscle': 'muscle', 'repro': 'reproduction',
    'hypoderm': 'epithelium', 'epithelium': 'epithelium',
    'sheath': 'glial', 'socket': 'glial', 'coelomocyte': 'coelomocyte',
    'excretory': 'excretory', 'mesoderm': 'mesoderm', 'other': 'other',
    'intestine': 'alimentary', 'valve': 'alimentary', 'marginal': 'alimentary',
    'gland': 'alimentary', 'rectal': 'alimentary',
    'tail': 'programmed_death', 'programmed_death': 'programmed_death', 'blast': 'programmed_death',
}


# ═══════════════════════════════════════════════════════════════
# Pareto optimisation
# ═══════════════════════════════════════════════════════════════

def lineage_edge_ratio(ri, ci, terminal_parents):
    """Fraction of real parent→child edges preserved by assignment (ri, ci)."""
    return sum(1 for i, j in zip(ri, ci) if terminal_parents[i] == terminal_parents[j]) / len(terminal_parents)


def pareto_sweep(xyz_mat, exp_mat, terminal_parents, iteration=1000):
    """Run Pareto sweep (min-max normalised) and return key scale-invariant metrics."""
    xyz_n = (xyz_mat - xyz_mat.min()) / (xyz_mat.max() - xyz_mat.min())
    exp_n = (exp_mat - exp_mat.min()) / (exp_mat.max() - exp_mat.min())
    lx, le = xyz_n.diagonal().sum(), exp_n.diagonal().sum()

    best_er, best_dist, best_alpha = 0.0, np.inf, 0.0
    edge_ratios = []

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xyz_n + (1 - alpha) * exp_n)
        er = lineage_edge_ratio(ri, ci, terminal_parents)
        edge_ratios.append(er)
        cx, ce = xyz_n[ri, ci].sum(), exp_n[ri, ci].sum()
        d = np.sqrt((cx - lx)**2 + (ce - le)**2)
        if er > best_er:
            best_er = er
        if d < best_dist:
            best_dist, best_alpha = d, alpha

    # Edge ratios at the two extremes
    ri_e, ci_e = linear_sum_assignment(exp_n)
    er_expr = lineage_edge_ratio(ri_e, ci_e, terminal_parents)
    ri_s, ci_s = linear_sum_assignment(xyz_n)
    er_spatial = lineage_edge_ratio(ri_s, ci_s, terminal_parents)

    return dict(expr_opt_er=er_expr, spatial_opt_er=er_spatial, max_er=best_er,
                closest_dist=best_dist, edge_ratios=np.array(edge_ratios))


# ═══════════════════════════════════════════════════════════════
# Cost matrix construction
# ═══════════════════════════════════════════════════════════════

def build_cost_matrices(terminal_nodes, terminal_parents, xyz_map, exp_df, sel_features=None):
    """Build L2 spatial and cosine expression cost matrices.
    Returns (xyz_cost_mat, exp_cost_mat, full_exp_dict).
    full_exp_dict maps cell names to full (unselected) expression vectors."""
    sel = exp_df[sel_features] if sel_features is not None else exp_df
    n = len(terminal_parents)
    xm, em = np.zeros((n, n)), np.zeros((n, n))
    full_exp_dict = {}
    for i, p in enumerate(terminal_parents):
        for j, c in enumerate(terminal_nodes):
            xm[i, j] = np.linalg.norm(xyz_map[p] - xyz_map[c])
            em[i, j] = cosine(sel.loc[p].values, sel.loc[c].values)
    for node in terminal_nodes:
        full_exp_dict[node] = exp_df.loc[node].values
    for parent in terminal_parents:
        full_exp_dict[parent] = exp_df.loc[parent].values
    return xm, em, full_exp_dict


# ═══════════════════════════════════════════════════════════════
# Grandparent / cousin-group analysis
# ═══════════════════════════════════════════════════════════════

def build_grandparent_map(root):
    """Map each cell name to its grandparent (parent's parent).
    First cousins share the same grandparent but have different parents."""
    gp_map = {}

    def dfs(node, parent_name=None, grandparent_name=None):
        name = map_names(node.get("did", ""))
        gp_map[name] = grandparent_name
        for child in node.get("children", []):
            dfs(child, name, parent_name)

    dfs(root)
    return gp_map


def build_cousin_groups(terminal_nodes, gp_mapping):
    """Group terminal cell indices by grandparent. Only groups with >=2 members."""
    groups = {}
    for idx, node in enumerate(terminal_nodes):
        gp = gp_mapping.get(node)
        if gp is not None:
            groups.setdefault(gp, []).append(idx)
    cousin_groups = [indices for gp, indices in groups.items() if len(indices) >= 2]
    cousin_groups.sort(key=len, reverse=True)
    return cousin_groups


def build_ancestor_groups(terminal_nodes, tree_index, ancestor_steps):
    """Group terminal indices by a shared ancestor at a fixed depth.

    ``ancestor_steps=2`` gives first-cousin groups (shared grandparent),
    ``3`` gives second-cousin groups, and ``4`` gives third-cousin groups.
    """
    groups = {}
    for idx, node in enumerate(terminal_nodes):
        info = tree_index.get(node)
        if info is None or ancestor_steps < 1:
            continue
        ancestors = info.get('ancestors', [])
        if ancestor_steps - 1 >= len(ancestors):
            continue
        ancestor = ancestors[ancestor_steps - 1]
        groups.setdefault(ancestor, []).append(idx)
    result = [indices for indices in groups.values() if len(indices) >= 2]
    result.sort(key=len, reverse=True)
    return result


def compute_group_shuffle_costs(xyz_mat, exp_mat, groups,
                                n_random=1000, seed=42):
    """Sample total costs after permutations restricted to ``groups``."""
    rng = np.random.default_rng(seed)
    n_cells = xyz_mat.shape[0]
    rows = np.arange(n_cells)
    random_xyz = np.empty(n_random)
    random_exp = np.empty(n_random)
    for sample in range(n_random):
        perm = np.arange(n_cells)
        for group in groups:
            shuffled = np.asarray(group).copy()
            rng.shuffle(shuffled)
            perm[group] = shuffled
        random_xyz[sample] = xyz_mat[rows, perm].sum()
        random_exp[sample] = exp_mat[rows, perm].sum()
    return random_xyz, random_exp


def compute_cousin_random_stats(xyz_mat, exp_mat, cousin_groups, n_random=500, seed=42):
    """Compute mean/std of total costs under first-cousin-group permutation."""
    rng = np.random.default_rng(seed)
    n_cells = xyz_mat.shape[0]
    random_xyz, random_exp = [], []
    for _ in range(n_random):
        perm = list(range(n_cells))
        for group in cousin_groups:
            if len(group) < 2:
                continue
            shuffled = group.copy()
            rng.shuffle(shuffled)
            for orig_idx, new_idx in zip(group, shuffled):
                perm[orig_idx] = new_idx
        total_xyz = sum(xyz_mat[i, perm[i]] for i in range(n_cells))
        total_exp = sum(exp_mat[i, perm[i]] for i in range(n_cells))
        random_xyz.append(total_xyz)
        random_exp.append(total_exp)
    random_xyz = np.array(random_xyz)
    random_exp = np.array(random_exp)
    return {
        'xyz_mean': float(random_xyz.mean()), 'xyz_std': float(random_xyz.std()),
        'exp_mean': float(random_exp.mean()), 'exp_std': float(random_exp.std()),
        'lineage_xyz': float(np.diag(xyz_mat).sum()),
        'lineage_exp': float(np.diag(exp_mat).sum()),
        'random_xyz': random_xyz, 'random_exp': random_exp,
    }


def compute_std_scaled_pareto(xyz_mat, exp_mat, terminal_parents, random_stats, iteration=1000):
    """Pareto front in null-SD units, translated so the lineage is at (0, 0).

    Translation does not affect the assignment optimization. One unit remains
    one standard deviation of the first-cousin null distribution.
    """
    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']
    lineage_x = random_stats.get('lineage_xyz', np.diag(xyz_mat).sum()) / random_stats['xyz_std']
    lineage_e = random_stats.get('lineage_exp', np.diag(exp_mat).sum()) / random_stats['exp_std']
    xl, el, edge_list = [], [], []
    for i in range(iteration + 1):
        a = i / iteration
        ri, ci = linear_sum_assignment(a * xs + (1 - a) * es)
        xl.append(xs[ri, ci].sum() - lineage_x)
        el.append(es[ri, ci].sum() - lineage_e)
        edge_list.append(lineage_edge_ratio(ri, ci, terminal_parents))
    edge_list = np.array(edge_list)
    kp = {'expr_opt_idx': 0, 'spatial_opt_idx': iteration, 'xyz_opt_idx': iteration,
          'max_er_idx': int(np.argmax(edge_list))}
    return np.array(xl), np.array(el), edge_list, kp


def get_null_cloud(random_stats, max_points=200):
    """Return null samples in SD units relative to the natural lineage."""
    x = (random_stats['random_xyz'] - random_stats['lineage_xyz']) / random_stats['xyz_std']
    y = (random_stats['random_exp'] - random_stats['lineage_exp']) / random_stats['exp_std']
    if len(x) > max_points:
        rng = np.random.default_rng(42)
        sel = rng.choice(len(x), max_points, replace=False)
        x, y = x[sel], y[sel]
    return x, y


def lineage_std_position(xyz_mat, exp_mat, random_stats):
    """Lineage position in the lineage-centred display coordinates."""
    return 0.0, 0.0


def relative_pareto_distance(xyz_arr, exp_arr, lx, le, nx=0.0, ne=0.0):
    """||LP|| / ||NP|| for lineage L, closest front point P, and null mean N."""
    dists = np.sqrt((xyz_arr - lx)**2 + (exp_arr - le)**2)
    p_idx = np.argmin(dists)
    px, py = xyz_arr[p_idx], exp_arr[p_idx]
    lp = np.sqrt((lx - px)**2 + (le - py)**2)
    np_dist = np.sqrt((px - nx)**2 + (py - ne)**2)
    return lp / np_dist if np_dist > 1e-10 else np.inf


# ═══════════════════════════════════════════════════════════════
# Edge perturbation analysis
# ═══════════════════════════════════════════════════════════════

def edge_perturbation_choose2(xyz_mat, exp_mat):
    """k=2 edge swaps: enumerate all C(n,2) pairs, count cost-improving swaps."""
    n = xyz_mat.shape[0]
    diag_x = np.diag(xyz_mat)
    diag_e = np.diag(exp_mat)
    swap_x = xyz_mat + xyz_mat.T
    swap_e = exp_mat + exp_mat.T
    dx = swap_x - diag_x[:, None] - diag_x[None, :]
    de = swap_e - diag_e[:, None] - diag_e[None, :]
    triu = np.triu_indices(n, k=1)
    dx_u, de_u = dx[triu], de[triu]
    total = len(dx_u)
    saves_xyz = int(np.sum(dx_u < 0))
    saves_exp = int(np.sum(de_u < 0))
    saves_both = int(np.sum((dx_u < 0) & (de_u < 0)))
    return dict(total=total, saves_xyz=saves_xyz, saves_exp=saves_exp, saves_both=saves_both,
                pct_xyz=round(100 * saves_xyz / total, 3),
                pct_exp=round(100 * saves_exp / total, 3),
                pct_both=round(100 * saves_both / total, 3))


def edge_perturbation_choose3(xyz_mat, exp_mat):
    """k=3 cyclic edge swaps: for each triple (i,j,k), test both 3-cycles.
    Uses vectorised inner loop — fast even for n~300."""
    n = xyz_mat.shape[0]
    diag_x = np.diag(xyz_mat)
    diag_e = np.diag(exp_mat)
    delta_x = xyz_mat - diag_x[:, None]
    delta_e = exp_mat - diag_e[:, None]
    saves_xyz = saves_exp = saves_both = total = 0
    for i in range(n - 2):
        dix, die = delta_x[i], delta_e[i]
        for j in range(i + 1, n - 1):
            ks = slice(j + 1, n)
            kc = n - j - 1
            # Forward 3-cycle: i→j, j→k, k→i
            df_x = dix[j] + delta_x[j, ks] + delta_x[ks, i]
            df_e = die[j] + delta_e[j, ks] + delta_e[ks, i]
            # Backward 3-cycle: i→k, k→j, j→i
            db_x = dix[ks] + delta_x[ks, j] + delta_x[j, i]
            db_e = die[ks] + delta_e[ks, j] + delta_e[j, i]
            saves_xyz += int(np.sum(df_x < 0)) + int(np.sum(db_x < 0))
            saves_exp += int(np.sum(df_e < 0)) + int(np.sum(db_e < 0))
            saves_both += int(np.sum((df_x < 0) & (df_e < 0))) + int(np.sum((db_x < 0) & (db_e < 0)))
            total += 2 * kc
    p = lambda v: round(100 * v / total, 3) if total else 0.0
    return dict(total=total, saves_xyz=saves_xyz, saves_exp=saves_exp, saves_both=saves_both,
                pct_xyz=p(saves_xyz), pct_exp=p(saves_exp), pct_both=p(saves_both))


def run_perturbation_test(xyz_mat, exp_mat, label):
    """Run both choose-2 and choose-3 perturbation tests and print summary."""
    r2 = edge_perturbation_choose2(xyz_mat, exp_mat)
    r3 = edge_perturbation_choose3(xyz_mat, exp_mat)
    n = xyz_mat.shape[0]
    print(f"  {label:30s} n={n:4d} | C2={r2['total']:8d}  "
          f"save_xyz={r2['pct_xyz']:6.2f}%  save_exp={r2['pct_exp']:6.2f}%  "
          f"save_both={r2['pct_both']:5.2f}% | "
          f"C3={r3['total']:9d}  save_xyz={r3['pct_xyz']:6.2f}%  "
          f"save_exp={r3['pct_exp']:6.2f}%  save_both={r3['pct_both']:5.2f}%")
    return dict(label=label, n=n, c2=r2, c3=r3)


# ═══════════════════════════════════════════════════════════════
# Z-axis noise null model
# ═══════════════════════════════════════════════════════════════

def make_z_noise_map(xyz_map_3d, seed=42):
    """Replace z-coordinate with Gaussian noise calibrated to xy spatial scale.
    Noise sigma = pooled std of (x,y) coordinates, matching embryo size.
    elegans: [z,y,x], briggsae: [z,x,y] — z is always index 0."""
    rng = np.random.default_rng(seed)
    xy_vals = np.array([v[1:] for v in xyz_map_3d.values()])
    sigma = np.sqrt(np.mean(np.var(xy_vals, axis=0)))
    mu = np.mean([v[1:].mean() for v in xyz_map_3d.values()])
    noisy = {}
    for k, v in xyz_map_3d.items():
        new_z = rng.normal(mu, sigma)
        noisy[k] = np.array([new_z, v[1], v[2]])
    return noisy


def run_z_noise_pareto(xyz_map_3d, tn, tp, em_mat, lineage, gp_map, n_draws=10, base_seed=1234):
    """Run Pareto on 2D+z_noise for n_draws independent noise realisations.
    Returns aggregate dict with mean±std of metrics plus median-draw arrays."""
    results = []
    for d in range(n_draws):
        xyz_noise = make_z_noise_map(xyz_map_3d, seed=base_seed + d)
        n_cells = len(tp)
        xm_nz = np.zeros((n_cells, n_cells))
        for i, p in enumerate(tp):
            for j, c in enumerate(tn):
                xm_nz[i, j] = np.linalg.norm(xyz_noise[p] - xyz_noise[c])
        nz_groups = build_cousin_groups(tn, gp_map)
        nz_rs = compute_cousin_random_stats(xm_nz, em_mat, nz_groups, n_random=300, seed=42)
        xyz_arr, exp_arr, edge_arr, kp = compute_std_scaled_pareto(xm_nz, em_mat, tp, nz_rs)
        lx, le = lineage_std_position(xm_nz, em_mat, nz_rs)
        er_max = edge_arr[kp['max_er_idx']]
        er_expr = edge_arr[0]
        er_spatial = edge_arr[-1]
        dists = np.sqrt((xyz_arr - lx)**2 + (exp_arr - le)**2)
        p_idx = np.argmin(dists)
        rel_dist = (np.sqrt((lx - xyz_arr[p_idx])**2 + (le - exp_arr[p_idx])**2)
                    / np.sqrt(xyz_arr[p_idx]**2 + exp_arr[p_idx]**2))
        results.append(dict(max_er=er_max, expr_er=er_expr, spatial_er=er_spatial, rel_dist=rel_dist,
                            xyz_arr=xyz_arr, exp_arr=exp_arr, edge_arr=edge_arr, kp=kp))
    keys = ['max_er', 'expr_er', 'spatial_er', 'rel_dist']
    agg = {k: (np.mean([r[k] for r in results]), np.std([r[k] for r in results])) for k in keys}
    med_idx = np.argsort([r['max_er'] for r in results])[n_draws // 2]
    agg['median_xyz_arr'] = results[med_idx]['xyz_arr']
    agg['median_exp_arr'] = results[med_idx]['exp_arr']
    agg['median_edge_arr'] = results[med_idx]['edge_arr']
    agg['median_kp'] = results[med_idx]['kp']
    agg['n_draws'] = n_draws
    return agg


# ═══════════════════════════════════════════════════════════════
# Random-feature null model (from terminal_expression.ipynb)
# ═══════════════════════════════════════════════════════════════

def build_random_cost_matrices(terminal_parents, terminal_nodes, full_exp_dict,
                               n_features=20, n_random=1000, seed_base=0):
    """For each of n_random draws, randomly select n_features expression features,
    build a cosine-distance expression cost matrix.  Returns list of n_random matrices."""
    n_total_features = len(full_exp_dict[terminal_nodes[0]])

    def _one_random(seed):
        rng = np.random.default_rng(seed_base + seed)
        rand_indices = rng.choice(n_total_features, size=n_features, replace=False)
        mat = np.zeros((len(terminal_parents), len(terminal_nodes)))
        for row_idx, parent in enumerate(terminal_parents):
            for col_idx, node in enumerate(terminal_nodes):
                pvec = full_exp_dict[parent][rand_indices]
                nvec = full_exp_dict[node][rand_indices]
                mat[row_idx, col_idx] = cosine(pvec, nvec)
        return mat

    random_mats = Parallel(n_jobs=-1, prefer="processes")(
        delayed(_one_random)(seed) for seed in tqdm(range(n_random), desc="Random feature matrices")
    )
    return random_mats


# ═══════════════════════════════════════════════════════════════
# Replicate analysis helper
# ═══════════════════════════════════════════════════════════════

def run_replicate_pareto(xyz_map, valid_names, exp_df, sel_features,
                         lineage_data, gp_map, n_random=300, seed=42):
    """Full pipeline for one tracking replicate against one expression config.
    Returns dict with n, xm, em, xyz_arr, exp_arr, edge_arr, kp, rs."""
    v = [n for n in valid_names if n in exp_df.index]
    tn, tp = collect_terminals(lineage_data, v)
    xm, em, _ = build_cost_matrices(tn, tp, xyz_map, exp_df, sel_features)
    groups = build_cousin_groups(tn, gp_map)
    rs = compute_cousin_random_stats(xm, em, groups, n_random=n_random, seed=seed)
    xyz_arr, exp_arr, edge_arr, kp = compute_std_scaled_pareto(xm, em, tp, rs)
    return dict(n=len(tn), xm=xm, em=em, xyz_arr=xyz_arr, exp_arr=exp_arr,
                edge_arr=edge_arr, kp=kp, rs=rs)


# ═══════════════════════════════════════════════════════════════
# Cell type utilities
# ═══════════════════════════════════════════════════════════════

def type_shannon_entropy(tn, get_ct_func):
    """Shannon entropy of merged cell-type distribution (in nats). 0 = pure single type."""
    counts = defaultdict(int)
    for node in tn:
        counts[get_ct_func(node)] += 1
    total = sum(counts.values())
    probs = np.array([c / total for c in counts.values()])
    return float(-np.sum(probs * np.log(probs + 1e-12)))


def decompose_by_cell_type(tn, tp, type_map, major_types=None):
    """Group terminal cell indices by cell type.
    major_types: set of types to keep separate; others bundled as 'rest (other types)'.
    Returns dict: type_name -> list of indices."""
    if major_types is None:
        major_types = {'neuron', 'muscle', 'epithelium', 'programmed_death'}
    groups = defaultdict(list)
    for i, (node, parent) in enumerate(zip(tn, tp)):
        ct = type_map.get(node, 'programmed_death')
        group = ct if ct in major_types else 'rest (other types)'
        groups[group].append(i)
    return dict(groups)
