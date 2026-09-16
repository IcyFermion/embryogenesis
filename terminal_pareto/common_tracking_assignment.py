"""Minimax common parent assignment on the existing matched tracking cache.

Minimize the worst travel/natural-travel ratio across three embryos, with
an explicit expression budget. Collapse duplicate parent slots to capacities.
Outputs are exploratory and separate from all publication caches.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix, vstack

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = ROOT/'terminal_pareto/output/tracking_geometry_sensitivity/assignments_and_costs.npz'
DEFAULT_OUT = ROOT/'terminal_pareto/output/common_tracking_assignment'


def optimize_common(travel, expression, capacity, natural_travel, expression_budget,
                    expression_scale, time_limit=45, integral=True):
    """Return the solver result; z is the worst fractional natural travel cost."""
    nr, m, n = travel.shape
    assert expression.shape == (m, n) and sum(capacity) == n
    nv = m*n
    idx = np.arange(nv)
    # Each parent keeps its available child count; every child gets one parent.
    incidence = coo_matrix((np.ones(2*nv),
        (np.r_[np.repeat(np.arange(m), n), m+np.tile(np.arange(n), m)],
         np.r_[idx, idx])), shape=(m+n, nv+1)).tocsr()
    travel_rows = np.column_stack((travel.reshape(nr, nv)/natural_travel[:, None], -np.ones(nr)))
    expression_row = np.r_[expression.ravel()/expression_scale, 0][None, :]
    matrix = vstack((incidence, csr_matrix(travel_rows), csr_matrix(expression_row)), format='csc')
    equal = np.r_[capacity, np.ones(n)]
    lower = np.r_[equal, np.full(nr+1, -np.inf)]
    upper = np.r_[equal, np.zeros(nr), expression_budget/expression_scale]
    objective = np.r_[np.zeros(nv), 1.]
    integrality = np.r_[np.ones(nv, dtype=int), 0] if integral else None
    return milp(objective, integrality=integrality,
                bounds=Bounds(np.zeros(nv+1), np.r_[np.ones(nv), np.inf]),
                constraints=LinearConstraint(matrix, lower, upper),
                options=dict(time_limit=time_limit, mip_rel_gap=1e-7))


def validate_assignment(result, travel, expression, capacity, budget):
    m, n = expression.shape
    assignment = result.x[:-1].reshape(m, n)
    np.testing.assert_allclose(assignment, np.round(assignment), atol=1e-5)
    assignment = np.round(assignment).astype(int)
    assert np.array_equal(assignment.sum(axis=1), capacity)
    assert np.array_equal(assignment.sum(axis=0), np.ones(n))
    t = np.sum(travel*assignment[None, :, :], axis=(1, 2))
    s = np.sum(expression*assignment)
    assert s <= budget+1e-5
    return assignment.argmax(axis=0), t, s


def selftest():
    rng = np.random.default_rng(55)
    # Duplicate biological parent: capacities 2, 1, 1.
    capacity = np.array([2, 1, 1])
    natural = np.array([0, 0, 1, 2])
    for _ in range(3):
        t = rng.uniform(.1, 1, (3, 3, 4))
        e = rng.uniform(.1, 1, (3, 4))
        nt = t[:, natural, np.arange(4)].sum(axis=1)
        budget = e[natural, np.arange(4)].sum()
        candidates = np.unique(list(itertools.permutations(natural)), axis=0)
        values = []
        for p in candidates:
            if e[p, np.arange(4)].sum() <= budget+1e-10:
                values.append(np.max(t[:, p, np.arange(4)].sum(axis=1)/nt))
        result = optimize_common(t, e, capacity, nt, budget, budget)
        assert result.success
        np.testing.assert_allclose(result.fun, min(values), atol=1e-7)
        validate_assignment(result, t, e, capacity, budget)
    print('Exhaustive minimax/capacity self-test passed.', flush=True)


def run(cache, out, reductions, seconds):
    out.mkdir(parents=True, exist_ok=True)
    data = np.load(cache)
    parents = data['natural_parents']
    children = data['children']
    unique, first, inv, capacity = np.unique(parents, return_index=True,
                                            return_inverse=True, return_counts=True)
    labels = ['embryo1', 'embryo2', 'embryo3']
    all_travel = np.array([data[r+'_travel_matrix'] for r in labels])
    all_expression = data['expression_matrix']
    t, e = all_travel[:, first], all_expression[first]
    np.testing.assert_array_equal(all_travel, t[:, inv])
    np.testing.assert_array_equal(all_expression, e[inv])
    nt = np.trace(all_travel, axis1=1, axis2=2)
    ns = np.trace(all_expression)
    results, edge_rows = [], []
    for reduction in reductions:
        budget = ns*(1-reduction/100)
        result = optimize_common(t, e, capacity, nt, budget, ns, seconds)
        row = dict(expression_reduction_required_pct=reduction,
                   expression_budget=budget, solver_status=int(result.status),
                   solver_message=result.message, seconds_limit=seconds,
                   mip_gap=getattr(result, 'mip_gap', None),
                   objective_lower_bound=getattr(result, 'mip_dual_bound', None))
        if result.x is not None:
            assignment, travel, state = validate_assignment(result, t, e, capacity, budget)
            worst = np.max(travel/nt)
            np.testing.assert_allclose(worst, result.fun, atol=1e-6)
            row.update(worst_travel_ratio=worst,
                       guaranteed_travel_reduction_pct=100*(1-worst),
                       expression=state, expression_natural=ns,
                       expression_reduction_pct=100*(1-state/ns),
                       natural_edge_retention=np.mean(unique[assignment] == parents))
            for r, cost, nat in zip(labels, travel, nt):
                row[r+'_travel'] = cost
                row[r+'_natural_travel'] = nat
                row[r+'_travel_reduction_pct'] = 100*(1-cost/nat)
            for child, parent, new_parent in zip(children, parents, unique[assignment]):
                edge_rows.append(dict(expression_reduction_required_pct=reduction,
                                      child=child, natural_parent=parent,
                                      common_parent=new_parent, retained=parent == new_parent))
        results.append(row)
        pd.DataFrame(results).to_csv(out/'common_assignment_summary.csv', index=False)
        pd.DataFrame(edge_rows).to_csv(out/'common_assignments.csv', index=False)
        print(json.dumps(row, default=float), flush=True)
    provenance = dict(cache=str(cache), cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      numpy=np.__version__, scipy=scipy.__version__, reductions=reductions,
                      definition='min max_r T_r(A)/T_r(natural), subject to S(A)<=budget')
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--reductions', nargs='+', type=float, default=[0, 1, 2, 5])
    parser.add_argument('--seconds', type=float, default=45)
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        selftest()
    else:
        run(args.cache, args.out, args.reductions, args.seconds)
