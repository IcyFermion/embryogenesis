"""Tree-aware distance metrics for comparing real lineage assignments to Pareto-optimal ones.

Extends the pareto_engine.py framework with metrics that account for the hierarchical
structure of the cell lineage tree. Key insight: not all edge changes are equal —
swapping two first-cousins (tree distance 2-4) is less "difficult" than reassigning
a cell across major branches (tree distance 10+).

All functions are pure: input arrays/dicts, output arrays/dicts. No file I/O, no plotting.
"""

import numpy as np
from scipy.optimize import linear_sum_assignment
from collections import defaultdict
from terminal_pareto.data_loader import map_names


# ═══════════════════════════════════════════════════════════════
# Tree Index
# ═══════════════════════════════════════════════════════════════

def build_lineage_tree_index(root):
    """Build a comprehensive index of the lineage tree for fast distance/LCA queries.

    Returns dict: node_name -> {depth, parent, ancestors (list from parent to root), name}

    Complexity: O(N * max_depth) time, O(N * max_depth) space
    """
    tree_index = {}

    def dfs(node, parent_name=None, depth=0):
        name = node.get('did', node.get('name', ''))
        mapped = map_names(name)
        if mapped not in tree_index:
            tree_index[mapped] = {'depth': depth, 'parent': parent_name, 'name': name}
        for child in node.get('children', []):
            dfs(child, mapped, depth + 1)

    dfs(root)

    # Post-process: fill ancestor lists by walking up through parents
    for node_name, info in tree_index.items():
        ancestors = []
        curr = info['parent']
        while curr is not None and curr in tree_index:
            ancestors.append(curr)
            curr = tree_index[curr]['parent']
        info['ancestors'] = ancestors

    return tree_index


# ═══════════════════════════════════════════════════════════════
# Tree Queries
# ═══════════════════════════════════════════════════════════════

def lineage_tree_distance(n1, n2, tree_index):
    """Path length (#edges) between two nodes going through their LCA.

    Complexity: O(max_depth) — depth is ≤ 18 for C. elegans lineage.
    """
    if n1 == n2:
        return 0
    e1, e2 = tree_index.get(n1), tree_index.get(n2)
    if e1 is None or e2 is None:
        return None

    chain1 = [n1] + e1['ancestors']
    chain2 = [n2] + e2['ancestors']

    # Find LCA: walk up from n1 until we hit n2's chain
    chain2_idx = {n: i for i, n in enumerate(chain2)}
    for i, a1 in enumerate(chain1):
        if a1 in chain2_idx:
            return i + chain2_idx[a1]

    return None  # no common ancestor (should not happen if root is universal)


def lineage_lca(n1, n2, tree_index):
    """Lowest common ancestor of two nodes."""
    if n1 == n2:
        return n1
    e1, e2 = tree_index.get(n1), tree_index.get(n2)
    if e1 is None or e2 is None:
        return None

    chain1 = [n1] + e1['ancestors']
    chain2_set = set([n2] + e2['ancestors'])
    for a in chain1:
        if a in chain2_set:
            return a
    return None


def ancestor_at_level(node_name, level, tree_index):
    """Return the ancestor `level` steps above node (0=self, 1=parent, 2=grandparent, ...).

    Returns None if level exceeds available ancestors (past root).
    """
    if level == 0:
        return node_name
    e = tree_index.get(node_name)
    if e is None:
        return None
    ancestors = e['ancestors']
    if level - 1 < len(ancestors):
        return ancestors[level - 1]
    return None


def cousin_level(n1, n2, tree_index):
    """Determine the cousin relationship level between two terminal cells.

    Returns:
        0: same cell
        1: siblings (same parent) — but for terminal cells each parent has ONE child,
           so same parent never occurs between distinct cells in practice
        2: first cousins (same grandparent)
        3: second cousins (same great-grandparent)
        k: k-th cousins (share ancestor at level k above)
        None: no common ancestor found

    Note: For the terminal cell assignment problem, cells are distinct leaves.
    'cousin_level' between parent p_i and assigned parent q_i tells us how
    far "up" the tree the reassignment goes.
    """
    if n1 == n2:
        return 0
    lca = lineage_lca(n1, n2, tree_index)
    if lca is None:
        return None
    d1, d2 = tree_index[n1]['depth'], tree_index[n2]['depth']
    lca_depth = tree_index[lca]['depth']
    # The max of the two distances to LCA tells us the cousin level
    return max(d1 - lca_depth, d2 - lca_depth)


# ═══════════════════════════════════════════════════════════════
# Approach 1: Tree-Weighted Edge Retention (TWER)
# ═══════════════════════════════════════════════════════════════

def tree_weighted_edge_scores(col_indices, terminal_parents, tree_index, normalize='max_depth'):
    """Compute tree-distance-weighted edge retention score for a single assignment.

    For each cell i, the real parent is p_i = terminal_parents[i].
    The Pareto assignment gives: parent p_i is assigned to child at column col_indices[i],
    which is terminal_nodes[col_indices[i]]. The "lineage edge" for this pair is the
    tree distance between p_i and the actual parent of terminal_nodes[col_indices[i]].

    Wait — in the assignment framework:
      - Row i = parent terminal_parents[i]
      - Column j = child terminal_nodes[j]
      - Assignment: parent i → child col_indices[i]

    The true developmental edge is: parent p_i → child terminal_nodes[i] (the diagonal).

    When the Pareto assignment assigns parent i → child col_indices[i] ≠ i,
    we want to know: how far in tree-space is this reassignment?

    The "edge change distance" for cell i:
      d_i = tree_distance(terminal_parents[i], terminal_parents[col_indices[i]])

    In the true lineage, terminal_nodes[i] and terminal_parents[i] are parent-child
    (tree distance 1). But since terminal_parents[i] is the parent of terminal_nodes[i],
    and we're assigning parent i to a DIFFERENT child, the tree distance between
    terminal_parents[i] and the true parent of the assigned child measures how
    "wrong" the assignment is.

    Parameters:
        col_indices: permutation array from linear_sum_assignment (length n)
        terminal_parents: list of parent cell names (length n)
        tree_index: lineage tree index from build_lineage_tree_index()
        normalize: 'max_depth' (divide by 2*max_depth), 'none', or 'log_depth'

    Returns:
        scores: array of per-cell scores (0 = preserved, 1 = maximally distant)
        mean_score: average across all cells
    """
    n = len(terminal_parents)
    max_dist = max(v['depth'] for v in tree_index.values()) * 2
    scores = np.zeros(n)

    for i in range(n):
        # p_i: the real parent of cell i in the biological lineage
        p_real = terminal_parents[i]
        # The Pareto assignment: parent i → child col_indices[i]
        # p_assigned: the PARENT of the child that was assigned to parent i
        # (terminal_parents[col_indices[i]] is the real parent of the assigned child)
        p_assigned = terminal_parents[col_indices[i]]

        if p_real == p_assigned:
            scores[i] = 0.0  # edge effectively preserved
        else:
            d = lineage_tree_distance(p_real, p_assigned, tree_index)
            if d is not None:
                if normalize == 'max_depth':
                    scores[i] = d / max_dist
                else:
                    scores[i] = d
            else:
                scores[i] = 1.0  # worst case for missing nodes

    return scores, float(np.mean(scores))


def tree_weighted_er_along_pareto(xyz_mat, exp_mat, terminal_parents, terminal_nodes,
                                   tree_index, random_stats, iteration=1000):
    """Compute tree-weighted edge retention along the entire Pareto front.

    Combines compute_std_scaled_pareto with tree_weighted_edge_scores to produce
    a TWER curve that parallels the existing ER curve.

    Returns:
        twer_scores: array of mean TWER at each alpha step
        per_cell_scores: array of shape (iteration+1, n) with per-cell scores
    """
    from terminal_pareto.pareto_engine import compute_std_scaled_pareto

    # Get the std-scaled cost matrices
    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']

    n = len(terminal_parents)
    twer_scores = np.zeros(iteration + 1)
    per_cell_scores = np.zeros((iteration + 1, n))

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        per_cell_scores[i, :], twer_scores[i] = tree_weighted_edge_scores(
            ci, terminal_parents, tree_index)

    return twer_scores, per_cell_scores


# ═══════════════════════════════════════════════════════════════
# Approach 2: Cousin-Stratified Edge Retention (CSER)
# ═══════════════════════════════════════════════════════════════

def stratified_edge_retention(col_indices, terminal_parents, tree_index,
                               max_level=4):
    """Decompose edge retention by cousin distance level.

    For each cell i:
      - Level 0: edge preserved (same parent assigned)
      - Level 1: different parent, tree distance 1 (same parent? unlikely in terminal case)
      - Level 2: tree distance 2 (first cousins: swapped within same grandparent)
      - Level 3: tree distance 3-4
      - Level 4+: beyond (deep cross-branch reassignments)

    Returns:
        level_counts: dict mapping level -> count of cells at that level
        level_fractions: dict mapping level -> fraction of total cells
    """
    n = len(terminal_parents)
    level_counts = defaultdict(int)

    # Map tree distances to levels
    # Level 0 = distance 0 (same parent — edge preserved)
    # Level 1 = distance 1
    # Level 2 = distance 2 (same grandparent — siblings have distance 2)
    # Level 3 = distance 3-4 (close cousins)
    # Level 4+ = distance 5+ (deep changes)

    distance_to_level = {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 3,
    }

    for i in range(n):
        p_real = terminal_parents[i]
        p_assigned = terminal_parents[col_indices[i]]
        if p_real == p_assigned:
            level_counts[0] += 1
        else:
            d = lineage_tree_distance(p_real, p_assigned, tree_index)
            if d is not None:
                level = distance_to_level.get(d, 4)
                level_counts[level] += 1
            else:
                level_counts[4] += 1

    level_fractions = {k: v / n for k, v in level_counts.items()}
    return dict(level_counts), dict(level_fractions)


def cousin_stratified_er_along_pareto(xyz_mat, exp_mat, terminal_parents,
                                       tree_index, random_stats, iteration=1000):
    """Compute cousin-stratified edge retention at key points along the Pareto front.

    Returns dict with stratified fractions at:
      - expression-optimal (alpha=0)
      - balanced (alpha=0.5)
      - spatial-optimal (alpha=1)
      - max-edge-retention point
    """
    from terminal_pareto.pareto_engine import compute_std_scaled_pareto

    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']

    n = len(terminal_parents)
    edge_arr = np.zeros(iteration + 1)
    all_levels = []

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        edge_arr[i] = sum(1 for k, j in zip(ri, ci)
                          if terminal_parents[k] == terminal_parents[j]) / n
        _, fractions = stratified_edge_retention(ci, terminal_parents, tree_index)
        all_levels.append(fractions)

    max_er_idx = int(np.argmax(edge_arr))

    return {
        'expr_opt': all_levels[0],
        'balanced': all_levels[iteration // 2],
        'spatial_opt': all_levels[iteration],
        'max_er': all_levels[max_er_idx],
        'all_levels': all_levels,
        'edge_arr': edge_arr,
        'max_er_idx': max_er_idx,
    }


# ═══════════════════════════════════════════════════════════════
# Approach 3: Edge Change Profile along Pareto Front
# ═══════════════════════════════════════════════════════════════

def edge_change_profile(xyz_mat, exp_mat, terminal_parents, tree_index,
                         random_stats, iteration=1000):
    """Track the tree-distance magnitude of edge changes at each step along the Pareto front.

    Uses a symmetric neighbour comparison: at each α step, computes the mean tree
    distance of edges that differ between the current assignment and the average of
    its two neighbours (previous and next). At endpoints, only one neighbour exists.

    This reveals whether "easy" changes (close cousins) happen at different α values
    than "hard" changes (distant branch reassignments).

    Returns:
        steps: alpha values where changes occurred
        changed_tree_dists: tree distance of edges that changed at each step
        cumulative_changes: cumulative count of edge changes
        mean_tree_dist: mean tree distance of changes at each step (symmetric)
    """
    from terminal_pareto.pareto_engine import compute_std_scaled_pareto

    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']

    n = len(terminal_parents)
    all_ci = []

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        all_ci.append(ci.copy())

    steps = []
    changed_tree_dists = []
    cumulative_changes = []
    mean_tree_dists = []

    def _changed_dists(ci_a, ci_b):
        """Return list of tree distances for edges that differ between two assignments."""
        mask = ci_a != ci_b
        indices = np.where(mask)[0]
        dists = []
        for idx in indices:
            p_a = terminal_parents[ci_a[idx]]
            p_b = terminal_parents[ci_b[idx]]
            if p_a == p_b:
                dists.append(0)
            else:
                d = lineage_tree_distance(p_a, p_b, tree_index)
                dists.append(d if d is not None else 0)
        return indices, dists

    for i in range(iteration + 1):
        alpha = i / iteration
        ci = all_ci[i]

        if i == 0:
            # Only compare with next
            indices, dists = _changed_dists(ci, all_ci[1])
        elif i == iteration:
            # Only compare with previous
            indices, dists = _changed_dists(all_ci[i - 1], ci)
        else:
            # Symmetric: average of (prev→current) and (current→next)
            _, dists_prev = _changed_dists(all_ci[i - 1], ci)
            _, dists_next = _changed_dists(ci, all_ci[i + 1])
            # Pool all changed edges from both comparisons
            # Use union of changed indices for a symmetric view
            mask_prev = all_ci[i - 1] != ci
            mask_next = ci != all_ci[i + 1]
            union_mask = mask_prev | mask_next
            indices = np.where(union_mask)[0]
            dists = []
            for idx in indices:
                ds = []
                if mask_prev[idx]:
                    p_p = terminal_parents[all_ci[i - 1][idx]]
                    p_c = terminal_parents[ci[idx]]
                    ds.append(lineage_tree_distance(p_p, p_c, tree_index) if p_p != p_c else 0)
                if mask_next[idx]:
                    p_c = terminal_parents[ci[idx]]
                    p_n = terminal_parents[all_ci[i + 1][idx]]
                    ds.append(lineage_tree_distance(p_c, p_n, tree_index) if p_c != p_n else 0)
                dists.append(np.mean(ds) if ds else 0)

        if len(indices) > 0:
            steps.append(alpha)
            changed_tree_dists.append(dists)
            mean_tree_dists.append(np.mean(dists) if dists else 0)

    return {
        'steps': np.array(steps),
        'changed_tree_dists': changed_tree_dists,
        'cumulative_changes': np.cumsum(np.array([len(d) for d in changed_tree_dists])),
        'mean_tree_dists': np.array(mean_tree_dists),
    }


# ═══════════════════════════════════════════════════════════════
# Approach 3b: Mean Lineage Distance along Pareto Front
# ═══════════════════════════════════════════════════════════════

def lineage_distance_along_pareto(xyz_mat, exp_mat, terminal_parents, tree_index,
                                   random_stats, iteration=1000):
    """Mean tree distance of ALL assigned edges from the real lineage, at each α.

    For each Pareto point, computes the tree distance between the real parent and
    the Pareto-assigned parent for every cell, then takes the mean.  Unlike TWER
    (which normalizes to a 0-1 retention score), this returns raw mean tree distance.
    Unlike the edge change profile (which only looks at edges that flip between
    neighbouring α steps), this compares every edge against the real lineage.

    Returns:
        mean_dists: array of mean tree distances (length iteration+1)
        per_cell: array of shape (iteration+1, n) with per-cell tree distances
    """
    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']

    n = len(terminal_parents)
    mean_dists = np.zeros(iteration + 1)
    per_cell = np.zeros((iteration + 1, n))

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        for j in range(n):
            p_real = terminal_parents[j]
            p_assigned = terminal_parents[ci[j]]
            if p_real == p_assigned:
                per_cell[i, j] = 0
            else:
                d = lineage_tree_distance(p_real, p_assigned, tree_index)
                per_cell[i, j] = d if d is not None else 0
        mean_dists[i] = per_cell[i, :].mean()

    return mean_dists, per_cell


def lineage_distance_null(terminal_parents, tree_index, n_samples=100000, seed=42):
    """Monte Carlo null distribution for mean tree distance from real lineage.

    Two null models:
      - full random: uniformly random permutation of all cells
      - cousin shuffle: permutation restricted to first-cousin groups

    Returns dict with keys:
        full_mean, full_std: mean/std of full random permutation
        cousin_mean, cousin_std: mean/std of cousin shuffle
    """
    from terminal_pareto.pareto_engine import build_cousin_groups

    n = len(terminal_parents)
    rng = np.random.default_rng(seed)
    cousin_groups = build_cousin_groups_from_tree(terminal_parents, tree_index)

    full_dists = np.zeros(n_samples)
    cousin_dists = np.zeros(n_samples)

    # Precompute pairwise tree distances between all parents
    dist_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                d = lineage_tree_distance(terminal_parents[i], terminal_parents[j], tree_index)
                dist_mat[i, j] = d if d is not None else 0

    for s in range(n_samples):
        # Full random
        perm_full = rng.permutation(n)
        full_dists[s] = dist_mat[np.arange(n), perm_full].mean()

        # Cousin shuffle
        perm_cousin = np.arange(n)
        for group in cousin_groups:
            if len(group) < 2:
                continue
            shuffled = group.copy()
            rng.shuffle(shuffled)
            for orig, new in zip(group, shuffled):
                perm_cousin[orig] = new
        cousin_dists[s] = dist_mat[np.arange(n), perm_cousin].mean()

    return {
        'full_mean': float(full_dists.mean()),
        'full_std':  float(full_dists.std()),
        'cousin_mean': float(cousin_dists.mean()),
        'cousin_std':  float(cousin_dists.std()),
        'full_dists': full_dists,
        'cousin_dists': cousin_dists,
    }


# ═══════════════════════════════════════════════════════════════
# Approach 3c: Cost-vs-Tree-Distance Trade-off
# ═══════════════════════════════════════════════════════════════

def cost_vs_tree_tradeoff(xyz_mat, exp_mat, terminal_parents, tree_index,
                           random_stats, iteration=1000):
    """Compute the trade-off between tree disruption and cost savings along the Pareto front.

    For each α point, relative to the real lineage:
      - Δxyz: xyz cost change (negative = worse spatial, positive = better)
      - Δexp: expression cost change
      - tree_dist: mean tree distance from real lineage
      - marginal_xyz: d(Δxyz)/d(tree_dist) — xyz efficiency
      - marginal_exp: d(Δexp)/d(tree_dist) — expression efficiency

    All costs in σ from null (z-score units).

    Returns dict with all arrays of length iteration+1.
    """
    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']

    # Real lineage costs (diagonal sums, in σ)
    real_xyz = xs.diagonal().sum()
    real_exp = es.diagonal().sum()

    n = len(terminal_parents)
    pareto_xyz = np.zeros(iteration + 1)
    pareto_exp = np.zeros(iteration + 1)
    tree_dists = np.zeros(iteration + 1)

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        pareto_xyz[i] = xs[ri, ci].sum()
        pareto_exp[i] = es[ri, ci].sum()

        # Mean tree distance from real lineage
        dists = np.zeros(n)
        for j in range(n):
            p_real = terminal_parents[j]
            p_assigned = terminal_parents[ci[j]]
            if p_real != p_assigned:
                d = lineage_tree_distance(p_real, p_assigned, tree_index)
                dists[j] = d if d is not None else 0
        tree_dists[i] = dists.mean()

    # Cost changes relative to real lineage (positive = saving / better than real)
    delta_xyz = real_xyz - pareto_xyz
    delta_exp = real_exp - pareto_exp

    # Marginal efficiency: d(cost)/d(tree) — how much cost change per unit tree dist
    marginal_xyz = np.zeros(iteration + 1)
    marginal_exp = np.zeros(iteration + 1)
    for i in range(1, iteration + 1):
        dt = tree_dists[i] - tree_dists[i - 1]
        if abs(dt) > 1e-10:
            marginal_xyz[i] = (delta_xyz[i] - delta_xyz[i - 1]) / dt
            marginal_exp[i] = (delta_exp[i] - delta_exp[i - 1]) / dt

    # Combined marginal (sum of both savings per tree unit)
    marginal_combined = marginal_xyz + marginal_exp

    return {
        'delta_xyz': delta_xyz,
        'delta_exp': delta_exp,
        'tree_dists': tree_dists,
        'marginal_xyz': marginal_xyz,
        'marginal_exp': marginal_exp,
        'marginal_combined': marginal_combined,
        'xyz_arr': np.array([pareto_xyz[i] - real_xyz / xs.shape[0] for i in range(iteration + 1)]),  # placeholder, not used
    }


# ═══════════════════════════════════════════════════════════════
# Approach 4: Normalized Lineage Assignment Distance (NLAD)
# ═══════════════════════════════════════════════════════════════

def lineage_assignment_distance(col_indices, terminal_parents, tree_index):
    """Total tree-distance cost of a permutation relative to the identity (real lineage).

    The real lineage corresponds to the identity permutation (col_indices[i] = i).
    The cost of a permutation π is sum_i tree_distance(terminal_parents[i], terminal_parents[π(i)]).

    This is the "total structural change" — how far in tree-space the cells have been reassigned.
    """
    n = len(terminal_parents)
    total_dist = 0.0
    per_cell = np.zeros(n)
    for i in range(n):
        p_real = terminal_parents[i]
        p_assigned = terminal_parents[col_indices[i]]
        if p_real != p_assigned:
            d = lineage_tree_distance(p_real, p_assigned, tree_index)
            per_cell[i] = d if d is not None else 0
        total_dist += per_cell[i]
    return total_dist, per_cell


def nlad_along_pareto(xyz_mat, exp_mat, terminal_parents, tree_index,
                       random_stats, iteration=1000, n_random=500, seed=42):
    """Normalized Lineage Assignment Distance along the Pareto front.

    The raw total tree distance is normalized by the distribution of total tree
    distances under cousin-random permutations, producing a z-score that's
    comparable across conditions (consistent with the existing std-based framework).

    Returns:
        nlad_scores: z-scored total tree distance at each alpha
        raw_distances: raw total tree distances
        null_mean, null_std: parameters of the null distribution
    """
    from terminal_pareto.pareto_engine import build_cousin_groups

    xs, es = xyz_mat.copy(), exp_mat.copy()
    xs /= random_stats['xyz_std']
    es /= random_stats['exp_std']

    n = len(terminal_parents)

    # Compute null distribution of total tree distances under cousin randomization
    # Build cousin groups from tree_index
    cousin_groups = build_cousin_groups_from_tree(terminal_parents, tree_index)

    rng = np.random.default_rng(seed)
    null_dists = np.zeros(n_random)
    for r in range(n_random):
        perm = np.arange(n)
        for group in cousin_groups:
            if len(group) < 2:
                continue
            shuffled = group.copy()
            rng.shuffle(shuffled)
            for orig, new in zip(group, shuffled):
                perm[orig] = new
        null_dists[r], _ = lineage_assignment_distance(perm, terminal_parents, tree_index)

    null_mean = float(np.mean(null_dists))
    null_std = float(np.std(null_dists))

    # Compute along Pareto
    raw_dists = np.zeros(iteration + 1)
    nlad_scores = np.zeros(iteration + 1)
    per_cell_mat = np.zeros((iteration + 1, n))

    for i in range(iteration + 1):
        alpha = i / iteration
        ri, ci = linear_sum_assignment(alpha * xs + (1 - alpha) * es)
        raw_dists[i], per_cell = lineage_assignment_distance(ci, terminal_parents, tree_index)
        per_cell_mat[i, :] = per_cell
        nlad_scores[i] = (raw_dists[i] - null_mean) / null_std if null_std > 0 else 0

    return {
        'nlad': nlad_scores,
        'raw': raw_dists,
        'null_mean': null_mean,
        'null_std': null_std,
        'null_dists': null_dists,
        'per_cell_mat': per_cell_mat,
    }


def build_cousin_groups_from_tree(terminal_parents, tree_index):
    """Build cousin groups using tree_index directly (no need for grandparent map)."""
    groups = defaultdict(list)
    for idx, parent in enumerate(terminal_parents):
        gp = ancestor_at_level(parent, 2, tree_index)  # grandparent
        if gp is not None:
            groups[gp].append(idx)
    result = [indices for indices in groups.values() if len(indices) >= 2]
    result.sort(key=len, reverse=True)
    return result


# ═══════════════════════════════════════════════════════════════
# Approach 5 (Bonus): Summary — Combined Distance Score
# ═══════════════════════════════════════════════════════════════

def combined_lineage_proximity(xyz_mat, exp_mat, terminal_parents, terminal_nodes,
                                tree_index, random_stats, iteration=1000, n_random=300):
    """Run all tree-aware metrics for a single configuration and return combined results.

    This is the "one-call" entry point for adding tree-aware metrics to the existing
    analysis pipeline.

    Returns a dict with keys:
        - n: number of terminal cells
        - twer: tree-weighted edge retention at each alpha
        - cser: cousin-stratified ER at key points
        - edge_profile: edge change profile
        - nlad: normalized lineage assignment distance at each alpha
        - traditional_er: traditional (flat) edge retention at each alpha
        - xyz_arr, exp_arr: Pareto front coordinates
        - kp: key point indices
    """
    from terminal_pareto.pareto_engine import compute_std_scaled_pareto

    # Traditional Pareto sweep (for ER comparison)
    xyz_arr, exp_arr, edge_arr, kp = compute_std_scaled_pareto(
        xyz_mat, exp_mat, terminal_parents, random_stats, iteration=iteration
    )

    # 1. Tree-Weighted ER
    twer, per_cell_scores = tree_weighted_er_along_pareto(
        xyz_mat, exp_mat, terminal_parents, terminal_nodes,
        tree_index, random_stats, iteration=iteration
    )

    # 2. Cousin-Stratified ER
    cser = cousin_stratified_er_along_pareto(
        xyz_mat, exp_mat, terminal_parents, tree_index, random_stats, iteration=iteration
    )

    # 3. Edge Change Profile
    edge_profile = edge_change_profile(
        xyz_mat, exp_mat, terminal_parents, tree_index, random_stats, iteration=iteration
    )

    # 3b. Mean lineage distance from real lineage (all edges)
    lineage_mean_dist, per_cell_dist = lineage_distance_along_pareto(
        xyz_mat, exp_mat, terminal_parents, tree_index, random_stats, iteration=iteration
    )

    # 3c. Cost-vs-tree-distance trade-off
    cost_tree_tradeoff = cost_vs_tree_tradeoff(
        xyz_mat, exp_mat, terminal_parents, tree_index, random_stats, iteration=iteration
    )

    # Null baselines for lineage distance (full random + cousin shuffle)
    lineage_null = lineage_distance_null(terminal_parents, tree_index)

    # 4. Normalized Lineage Assignment Distance
    nlad = nlad_along_pareto(
        xyz_mat, exp_mat, terminal_parents, tree_index, random_stats,
        iteration=iteration, n_random=n_random
    )

    return {
        'n': len(terminal_parents),
        'twer': twer,
        'per_cell_scores': per_cell_scores,
        'cser': cser,
        'edge_profile': edge_profile,
        'lineage_mean_dist': lineage_mean_dist,
        'per_cell_dist': per_cell_dist,
        'lineage_null': lineage_null,
        'cost_tree_tradeoff': cost_tree_tradeoff,
        'nlad': nlad,
        'traditional_er': edge_arr,
        'xyz_arr': xyz_arr,
        'exp_arr': exp_arr,
        'kp': kp,
    }


# ═══════════════════════════════════════════════════════════════
# Quick Diagnostic
# ═══════════════════════════════════════════════════════════════

def print_comparison(results, label="", config_name=""):
    """Print a summary comparing traditional ER with tree-aware metrics."""
    if label:
        config_name = f"{config_name} ({label})" if config_name else label

    trad_er = results['traditional_er']
    twer = results['twer']
    kp = results['kp']
    nlad = results['nlad']
    cser = results['cser']

    print(f"\n{'='*70}")
    print(f"  {config_name}")
    print(f"{'='*70}")

    # Key points
    for pt_name, pt_idx in [('Expr Opt (α=0)', 0),
                              ('Max ER', kp['max_er_idx']),
                              ('Spatial Opt (α=1)', len(trad_er) - 1)]:
        print(f"  {pt_name:20s}: ER={trad_er[pt_idx]:.3f}  TWER={1-twer[pt_idx]:.3f}  "
              f"NLAD={nlad['nlad'][pt_idx]:+.1f}σ")

    # Cousin-stratified at max ER
    max_er_cser = cser['all_levels'][kp['max_er_idx']]
    print(f"\n  Cousin-stratified ER at Max ER:")
    level_labels = {0: 'Same parent', 1: 'Dist=1', 2: 'Same GP (1st cousins)',
                    3: 'Dist=3-4', 4: 'Deep (≥5)'}
    for level in sorted(max_er_cser.keys()):
        frac = max_er_cser.get(level, 0)
        bar = '█' * int(frac * 50)
        print(f"    L{level} {level_labels.get(level, '?'):25s}: {frac:.3f}  {bar}")

    # Edge change profile summary
    ep = results['edge_profile']
    if len(ep['steps']) > 0:
        print(f"\n  Edge change profile: {len(ep['steps'])} change events")
        print(f"    Mean tree dist of changes: {np.mean(ep['mean_tree_dists']):.1f}")
        print(f"    Min tree dist: {np.min(ep['mean_tree_dists']):.1f}  "
              f"Max tree dist: {np.max(ep['mean_tree_dists']):.1f}")
        # Did early changes tend to be easier?
        half = len(ep['steps']) // 2
        if half > 0:
            early_mean = np.mean(ep['mean_tree_dists'][:half])
            late_mean = np.mean(ep['mean_tree_dists'][half:])
            print(f"    Early (α<0.5) mean dist: {early_mean:.1f}  "
                  f"Late (α>0.5) mean dist: {late_mean:.1f}")
    print()
