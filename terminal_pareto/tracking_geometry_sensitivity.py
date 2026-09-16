"""Exploratory P0 tracking transfer; writes only output/tracking_geometry_sensitivity.

Run from repo root in dev. Assignments are parent-slot -> child permutations;
biological agreement is computed after inversion to child -> parent identity.
No claim of enumeration of the full discrete Pareto set is made.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from terminal_pareto import data_loader as dl, pareto_engine as pe
from terminal_pareto.figS3_ce_tracking_robustness import audit_replicates, _ordered_matched_terms
from terminal_pareto.subtree_analysis import exact_cousin_stats


def solve(cost):
    rows, cols = linear_sum_assignment(cost)
    assert np.array_equal(rows, np.arange(len(rows)))
    assert np.array_equal(np.sort(cols), rows)
    return cols


def evaluate(matrix, assignments):
    return matrix[np.arange(matrix.shape[0])[None, :], assignments].sum(axis=1)


def biological_parents(assignments, parents):
    return parents[np.argsort(assignments, axis=1)]


def sweep(x, e, sx, se, parents, intervals):
    alpha = np.linspace(0, 1, intervals + 1)
    assignments = np.array([solve(a*x/sx + (1-a)*e/se) for a in alpha])
    t, s = evaluate(x, assignments), evaluate(e, assignments)
    span = np.array([t[0]-t[-1], s[-1]-s[0]])
    assert np.all(span > 0)
    normalized = np.column_stack((t-t[-1], s-s[0])) / span
    steps = np.linalg.norm(np.diff(normalized, axis=0), axis=1)
    u = 1 - np.r_[0, np.cumsum(steps)] / steps.sum()
    bp = biological_parents(assignments, parents)
    return dict(alpha=alpha, assignments=assignments, t=t, s=s, u=u,
                bp=bp, retention=(bp == parents).mean(axis=1), sx=sx, se=se,
                j=alpha*t/sx+(1-alpha)*s/se)


def travel_gap_bounds(t, s, target, candidate_t, candidate_s):
    """Bounds on t - min travel subject to expression <= s (raw travel units).

    Feasible candidate solutions bound the constrained optimum from above.
    Exact weighted-sum optima supply supporting-line lower bounds. Unlike
    connecting sampled points, these remain valid with missed discrete points.
    The transferred assignment itself is also a feasible upper bound.
    """
    a = target['alpha'][1:]
    lower_optimum = np.max((target['j'][1:][None, :]
                          - (1-a)[None, :]*s[:, None]/target['se'])
                          / a[None, :], axis=1) * target['sx']
    feasible = candidate_s[None, :] <= s[:, None] + 1e-10
    upper_optimum = np.minimum(t, np.min(np.where(feasible, candidate_t, np.inf), axis=1))
    lower_gap, upper_gap = t-upper_optimum, t-lower_optimum
    assert np.min(upper_gap-lower_gap) > -1e-6
    return np.maximum(0, lower_gap), np.maximum(0, upper_gap)


def selftest():
    """Exhaustive small problem checks for bounds and parent-slot symmetry."""
    rng = np.random.default_rng(91)
    for _ in range(5):
        x, e = rng.uniform(.1, 2, (2, 5, 5))
        parents = np.arange(5)
        f = sweep(x, e, 1, 1, parents, 20)
        assignments = np.array(list(itertools.permutations(range(5))))
        t, s = evaluate(x, assignments), evaluate(e, assignments)
        lo, hi = travel_gap_bounds(t, s, f, f['t'], f['s'])
        truth = np.array([tt-t[s <= ss+1e-10].min() for tt, ss in zip(t, s)])
        assert np.all(lo <= truth+1e-8) and np.all(truth <= hi+1e-8)
    bp = biological_parents(np.array([[0, 1, 2], [1, 0, 2]]), np.array(['p', 'p', 'q']))
    assert np.array_equal(bp[0], bp[1])
    print('Self-test passed: exact enumeration bounds and slot symmetry.', flush=True)


def run(intervals, out):
    out.mkdir(parents=True, exist_ok=True)
    context, audit = audit_replicates()
    terms = _ordered_matched_terms(context)
    assert len(terms) == 275
    children, parents = map(np.array, zip(*terms))
    feature = context['protein'][context['selected']]
    e = cdist(feature.loc[parents].values, feature.loc[children].values, 'cosine')
    assert np.isfinite(e).all()
    groups = pe.build_cousin_groups(children, pe.build_grandparent_map(context['lineage']))
    matrices, fronts, archive, summary = {}, {}, {}, []
    for label, data in context['replicates'].items():
        xyz = data['xyz_map']
        x = cdist(np.array([xyz[p] for p in parents], dtype=float),
                  np.array([xyz[c] for c in children], dtype=float))
        # Verify the vectorized builder against the production builder.
        xx, ee, _ = pe.build_cost_matrices(children[:8], parents[:8], xyz,
                                          context['protein'], context['selected'])
        np.testing.assert_allclose(x[:8, :8], xx, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(e[:8, :8], ee, rtol=1e-12, atol=1e-12)
        mx, me, vx, ve, _ = exact_cousin_stats(x, e, groups)
        assert vx > 0 and ve > 0
        f = sweep(x, e, np.sqrt(vx), np.sqrt(ve), parents, intervals)
        matrices[label], fronts[label] = x, f
        archive[label+'_assignments'] = f['assignments']
        archive[label+'_travel_matrix'] = x
        d = np.hypot((f['t']-np.trace(x))/(f['t'][0]-f['t'][-1]),
                     (f['s']-np.trace(e))/(f['s'][-1]-f['s'][0]))
        i = int(np.argmin(d))
        summary.append(dict(replicate=label, n=len(terms), intervals=intervals,
                            unique_biological_assignments=len(np.unique(f['bp'], axis=0)),
                            u_L=f['u'][i], d_LP=d[i], max_retention=f['retention'].max(),
                            travel_natural=np.trace(x), travel_optimum=f['t'][-1],
                            travel_null_sd=f['sx'], expression_null_sd=f['se']))
        print(label, summary[-1], flush=True)

    archive.update(expression_matrix=e, children=children, natural_parents=parents,
                   alpha=next(iter(fronts.values()))['alpha'])
    np.savez_compressed(out/'assignments_and_costs.npz', **archive)
    pd.DataFrame(summary).to_csv(out/'front_summary.csv', index=False)
    audit.to_csv(out/'matching_audit.csv', index=False)
    pd.DataFrame(dict(child=children, natural_parent=parents)).to_csv(out/'matched_edges.csv', index=False)

    # Per-child neighborhoods of unique biological parents, never duplicated slots.
    unique_parents, parent_rows = np.unique(parents, return_index=True)
    rankings = {r: np.argsort(x[parent_rows], axis=0, kind='stable') for r, x in matrices.items()}
    pair_rows, neighborhood_rows = [], []
    for r, q in itertools.combinations(fronts, 2):
        xr, xq = matrices[r][parent_rows], matrices[q][parent_rows]
        pair_rows.append(dict(source=r, target=q, pearson=np.corrcoef(xr.ravel(), xq.ravel())[0, 1],
                              spearman=spearmanr(xr.ravel(), xq.ravel()).statistic))
        for k in (3, 5, 10, 20):
            for j, child in enumerate(children):
                nr, nq = set(rankings[r][:k, j]), set(rankings[q][:k, j])
                neighborhood_rows.append(dict(source=r, target=q, child=child, k=k,
                                               overlap=len(nr & nq)/k,
                                               chance_overlap=k/len(unique_parents)))
    pd.DataFrame(pair_rows).to_csv(out/'matrix_concordance.csv', index=False)
    pd.DataFrame(neighborhood_rows).to_csv(out/'neighborhoods.csv', index=False)

    transfer_rows, canonical_rows, local_rows = [], [], []
    pool = np.unique(np.vstack([f['assignments'] for f in fronts.values()] + [np.arange(len(terms))[None, :]]), axis=0)
    pool_s = evaluate(e, pool)
    for q, target in fronts.items():
        pool_t = evaluate(matrices[q], pool)
        parent_xyz = np.array([context['replicates'][q]['xyz_map'][p] for p in unique_parents], dtype=float)
        parent_dist = cdist(parent_xyz, parent_xyz)
        local_scale = np.sort(parent_dist, axis=1)[:, 1:6].mean(axis=1)
        unique_index = {p: i for i, p in enumerate(unique_parents)}
        for r, source in fronts.items():
            t = evaluate(matrices[q], source['assignments'])
            s = evaluate(e, source['assignments'])
            np.testing.assert_allclose(s, source['s'], atol=1e-12)
            j = target['alpha']*t/target['sx']+(1-target['alpha'])*s/target['se']
            regret = j-target['j']
            assert regret.min() > -1e-8
            lo, hi = travel_gap_bounds(t, s, target, pool_t, pool_s)
            gain = np.trace(matrices[q])-target['t'][-1]
            assert gain > 0
            for i, alpha in enumerate(source['alpha']):
                transfer_rows.append(dict(source=r, target=q, alpha=alpha, u_source=source['u'][i],
                    travel=t[i], expression=s[i], weighted_regret=max(0, regret[i]),
                    weighted_regret_pct_optimum=100*max(0, regret[i])/target['j'][i],
                    travel_gap_lower=lo[i], travel_gap_upper=hi[i],
                    travel_gap_lower_fraction_available=lo[i]/gain,
                    travel_gap_upper_fraction_available=hi[i]/gain,
                    parent_agreement=np.mean(source['bp'][i] == target['bp'][i]),
                    slot_agreement=np.mean(source['assignments'][i] == target['assignments'][i]),
                    natural_retention=source['retention'][i]))
            if r == q:
                np.testing.assert_allclose(regret, 0, atol=1e-8)
                continue
            for u in np.linspace(0, 1, 101):
                i, k = np.argmin(abs(source['u']-u)), np.argmin(abs(target['u']-u))
                canonical_rows.append(dict(source=r, target=q, u_requested=u,
                    u_source=source['u'][i], u_target=target['u'][k],
                    parent_agreement=np.mean(source['bp'][i] == target['bp'][k])))
            for alpha in (.1, .5, .9, 1.):
                i = int(round(alpha*intervals))
                a = np.array([unique_index[p] for p in source['bp'][i]])
                b = np.array([unique_index[p] for p in target['bp'][i]])
                for child_i, child in enumerate(children):
                    local_rows.append(dict(source=r, target=q, alpha=source['alpha'][i], child=child,
                        same_parent=a[child_i] == b[child_i],
                        replacement_distance_local_scale=parent_dist[a[child_i], b[child_i]]/local_scale[b[child_i]],
                        source_parent_target_spatial_rank=int(np.flatnonzero(rankings[q][:, child_i] == a[child_i])[0])+1))
    transfer = pd.DataFrame(transfer_rows)
    transfer.to_csv(out/'transfer.csv', index=False)
    pd.DataFrame(canonical_rows).to_csv(out/'canonical_agreement.csv', index=False)
    pd.DataFrame(local_rows).to_csv(out/'replacement_locality.csv', index=False)

    # Exact per-edge exclusion margins: all slots of the selected parent must
    # be forbidden. Each solve permits compensating changes anywhere else.
    margins = []
    for r, f in fronts.items():
        for alpha in (.1, .5, .9, 1.):
            i = int(round(alpha*intervals))
            a = f['alpha'][i]
            cost = a*matrices[r]/f['sx']+(1-a)*e/f['se']
            baseline = f['j'][i]
            natural_gain = np.trace(cost)-baseline
            for child_i, child in enumerate(children):
                modified = cost.copy()
                modified[parents == f['bp'][i, child_i], child_i] = np.inf
                alternative = solve(modified)
                margin = cost[np.arange(len(terms)), alternative].sum()-baseline
                assert margin > -1e-8
                margins.append(dict(replicate=r, alpha=a, child=child,
                    selected_parent=f['bp'][i, child_i], margin=max(0, margin),
                    margin_fraction_natural_gain=max(0, margin)/natural_gain))
        print('Exclusion margins complete:', r, flush=True)
    pd.DataFrame(margins).to_csv(out/'edge_exclusion_margins.csv', index=False)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    colors = dict(embryo1='#4477AA', embryo2='#EE6677', embryo3='#228833')
    for ax, (q, f) in zip(axes[0], fronts.items()):
        for r in fronts:
            df = transfer[(transfer.source == r) & (transfer.target == q)]
            ax.plot((df.travel-np.trace(matrices[q]))/f['sx'],
                    (df.expression-np.trace(e))/f['se'], color=colors[r],
                    lw=2 if r == q else 1, alpha=1 if r == q else .6, label=r)
        ax.scatter([0], [0], marker='+', color='black', s=80, zorder=5)
        ax.set(title=f'Target geometry: {q}', xlabel='Travel change (target null SD)',
               ylabel='Cell-state change (null SD)')
    axes[0, 0].legend(fontsize=8)
    for r, q in itertools.combinations(fronts, 2):
        df = transfer[(transfer.source == r) & (transfer.target == q)]
        axes[1, 0].plot(df.alpha, df.parent_agreement, label=f'{r[-1]} vs {q[-1]}')
    axes[1, 0].set(xlabel='Travel weight', ylabel='Biological parent agreement', ylim=(0, 1.03))
    axes[1, 0].legend(fontsize=8)
    for (r, q), df in transfer[transfer.source != transfer.target].groupby(['source', 'target']):
        axes[1, 1].plot(df.alpha, df.weighted_regret_pct_optimum, lw=1, label=f'{r[-1]} to {q[-1]}')
    axes[1, 1].set(xlabel='Travel weight', ylabel='Transfer regret (% target optimum)')
    axes[1, 1].legend(fontsize=7, ncol=2)
    neighborhoods = pd.DataFrame(neighborhood_rows)
    for (r, q), df in neighborhoods.groupby(['source', 'target']):
        g = df.groupby('k').overlap.mean()
        axes[1, 2].plot(g.index, g.values, marker='o', label=f'{r[-1]} vs {q[-1]}')
    axes[1, 2].set(xlabel='Number of nearest distinct parents', ylabel='Mean neighborhood overlap', ylim=(0, 1))
    axes[1, 2].legend(fontsize=8)
    fig.suptitle('Tracking geometry sensitivity: matched 275 terminal edges (exploratory)')
    fig.savefig(out/'tracking_geometry_review.pdf')
    fig.savefig(out/'tracking_geometry_review.png', dpi=150)
    plt.close(fig)
    sources = [Path(p) for _, p, _ in dl.CE_REPLICATES]
    sources += [Path(__file__), ROOT/'terminal_pareto/data_loader.py',
                ROOT/'terminal_pareto/figS3_ce_tracking_robustness.py',
                ROOT/'terminal_pareto/subtree_analysis.py', ROOT/'terminal_pareto/pareto_engine.py',
                ROOT/'data/cell_lineage.json']
    pd.util.hash_pandas_object(feature, index=True).values.tofile(out/'expression_table_hash_input.bin')
    sources += [out/'expression_table_hash_input.bin']
    provenance = dict(intervals=intervals, n=275, numpy=np.__version__, scipy=scipy.__version__,
        selected_features=list(context['selected']),
        sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Wrote', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--intervals', type=int, default=300)
    parser.add_argument('--out', type=Path, default=ROOT/'terminal_pareto/output/tracking_geometry_sensitivity')
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        selftest()
    else:
        if args.intervals < 10 or args.intervals % 10:
            parser.error('--intervals must be a positive multiple of 10, at least 10')
        run(args.intervals, args.out)
