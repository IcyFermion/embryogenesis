"""One-parameter protein branch-variability sensitivity; never updates Figure 8.

Run from the repository root with conda run -n dev python -m
full_tree_pareto.branch_variability_sensitivity. Outputs are isolated under
output/branch_variability/. The covariance is a training-set second-moment
plug-in estimate, NOT jointly optimized with the lognormal parameter.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp, roots_hermitenorm

from full_tree_pareto import publication_analysis as analysis
from full_tree_pareto.clock_sensitivity import edge_manifest, _summarize, _top_share


DEFAULT_OUTPUT = analysis.CACHE_ROOT / "branch_variability"
GAUSSIAN = "Separate-clock Gaussian"
MIXTURE = "Lognormal branch variability"


@lru_cache(maxsize=None)
def quadrature(order):
    nodes, weights = roots_hermitenorm(order)
    keep = weights > 0
    return nodes[keep], np.log(weights[keep]) - .5 * np.log(2 * np.pi)


def second_moment(delta):
    delta = np.asarray(delta, dtype=float)
    if delta.ndim != 2 or not len(delta) or not np.isfinite(delta).all():
        raise ValueError("Expected nonempty, finite edge vectors")
    covariance = delta.T @ delta / len(delta)
    np.linalg.cholesky(covariance)  # Reject singular fits; no silent regularization.
    return covariance


def radial_coordinates(delta, covariance):
    factor = np.linalg.cholesky(covariance)
    whitened = np.linalg.solve(factor, delta.T).T
    return (whitened**2).sum(axis=1), 2 * np.log(np.diag(factor)).sum()


def radial_logpdf(q, dimension, logdet, tau, order=512):
    """Density of full edge vectors, integrated over log(s), not norm density."""
    if tau < 0 or not np.isfinite(tau):
        raise ValueError("tau must be finite and nonnegative")
    normalizer = -.5 * (dimension * np.log(2 * np.pi) + logdet)
    q = np.asarray(q)
    if tau == 0:
        return normalizer - .5 * q
    nodes, logweights = quadrature(order)
    logscale = tau * nodes - .5 * tau**2
    kernel = (-.5 * dimension * logscale[:, None]
              - .5 * np.exp(-logscale[:, None]) * q[None, :])
    return normalizer + logsumexp(kernel + logweights[:, None], axis=0)


def fit_model(delta):
    """Fit tau by conditional vector likelihood, never by matching total cost."""
    covariance = second_moment(delta)
    q, logdet = radial_coordinates(delta, covariance)
    dim = delta.shape[1]
    objective = lambda tau: -radial_logpdf(q, dim, logdet, tau).sum()
    # Bracket the best grid point to avoid assuming global unimodality.
    grid = np.linspace(0, 3, 61)
    scores = np.array([objective(tau) for tau in grid])
    best = int(scores.argmin())
    fitted = minimize_scalar(objective, bounds=(grid[max(0, best-1)],
                                               grid[min(len(grid)-1, best+1)]),
                             method="bounded", options={"xatol": 1e-7})
    if not fitted.success:
        raise RuntimeError("Branch-variability fit did not converge")
    tau = min([0., float(fitted.x), 3.], key=objective)
    if tau >= 2.999:
        raise RuntimeError("tau hit search boundary; inspect before expanding it")
    ll = radial_logpdf(q, dim, logdet, tau)
    check = radial_logpdf(q, dim, logdet, tau, order=1024)
    error = float(np.max(np.abs(ll-check)))
    if error > 1e-5:
        raise RuntimeError(f"Quadrature failed convergence check: {error}")
    return covariance, tau, {
        "n_train": len(delta), "tau": tau, "scale_cv": np.sqrt(np.expm1(tau**2)),
        "expected_norm_multiplier": np.exp(-tau**2 / 8),
        "gaussian_loglik": -objective(0), "mixture_loglik": float(ll.sum()),
        "quadrature_max_abs_logpdf_error": error,
        "fit": "second-moment covariance; conditional vector likelihood for tau",
    }


def gaussian_mean_norm(covariance):
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0)
    def integrand(u):
        if u == 0:
            return 2 * eigenvalues.sum()
        return -np.expm1(-.5 * np.log1p(2 * u*u * eigenvalues).sum()) / (u*u)
    return quad(integrand, 0, np.inf, epsabs=1e-9)[0] / np.sqrt(np.pi)


def paired_norms(covariance, tau, n_edges, n_draws, seed):
    """Exact edge-cost simulation with paired Gaussian draws.

    Independent edge vectors induce unique states on a rooted forest once its
    root states are fixed. Costs depend only on increments, so materializing
    all ancestral states is unnecessary. Scalar mixing leaves directions
    unchanged. Separate RNG streams make the Gaussian baseline tau-invariant.
    """
    if n_edges < 1 or n_draws < 1 or tau < 0:
        raise ValueError("Positive edge/draw counts and nonnegative tau required")
    streams = np.random.SeedSequence(seed).spawn(2)
    normal_rng, scale_rng = [np.random.default_rng(s) for s in streams]
    factor = np.linalg.cholesky(covariance)
    gaussian = np.empty((n_draws, n_edges))
    mixed = np.empty_like(gaussian)
    for start in range(0, n_draws, 100):
        stop = min(start+100, n_draws)
        delta = normal_rng.normal(size=(stop-start, n_edges, len(covariance))) @ factor.T
        gaussian[start:stop] = np.linalg.norm(delta, axis=-1)
        logscale = tau * scale_rng.normal(size=(stop-start, n_edges)) - .5*tau**2
        mixed[start:stop] = gaussian[start:stop] * np.exp(.5 * logscale)
    return {GAUSSIAN: gaussian, MIXTURE: mixed}


def diagnostic_rows(norms_by_model, observed, groups=None, **metadata):
    rows = []
    for model, norms in norms_by_model.items():
        measures = {
            "total": (norms.sum(axis=1), observed.sum()),
            "squared_total": ((norms**2).sum(axis=1), (observed**2).sum()),
            "edge_norm_cv": (norms.std(axis=1)/norms.mean(axis=1), observed.std()/observed.mean()),
            "max_edge_norm": (norms.max(axis=1), observed.max()),
            "top_5pct_raw_sq_share": (_top_share(norms**2), _top_share(observed**2)),
        }
        for group, mask in (groups or {}).items():
            measures[f"mean_edge_norm_{group}"] = (norms[:, mask].mean(axis=1), observed[mask].mean())
        for metric, (simulated, natural) in measures.items():
            rows.append({**metadata, "model": model, "metric": metric,
                         "observed": natural, **_summarize(simulated)})
    return rows


def subtree_folds(context, edges, depth=3):
    """Hold out internal edges of depth-3 subtrees, purging shared-node training.

    Incoming stem edges are not tested. Purging all edges touching any held-out
    node prevents their expression states from entering the covariance fit.
    Feature selection and atlas standardization remain fixed, so this is
    conditional within-atlas validation, not independent-embryo validation.
    """
    nodes = analysis.node_manifest(context).set_index("tree_id")
    parent = np.asarray(context.tree.parent_list)
    ids = edges.tree_id.to_numpy()
    parents = edges.parent_tree_id.to_numpy()
    covered = np.zeros(len(edges), dtype=int)
    for root in nodes.index[nodes.canonical_depth.eq(depth)]:
        descendants = {int(root)}
        for node in nodes.index:
            if parent[node] in descendants:
                descendants.add(int(node))
        child_in = np.isin(ids, list(descendants))
        parent_in = np.isin(parents, list(descendants))
        test = child_in & parent_in
        train = ~(child_in | parent_in)
        if not test.any():
            continue
        train_nodes = set(ids[train]) | set(parents[train])
        test_nodes = set(ids[test]) | set(parents[test])
        if train_nodes & test_nodes:
            raise AssertionError("Held-out and training edges share measured nodes")
        covered += test
        yield str(nodes.loc[root, "cell"]), train, test
    if np.any(covered > 1):
        raise AssertionError("Subtree test sets overlap")


def raw_audit(context, edges, count=20):
    """Audit largest raw-norm edges, not historical time-rate outliers."""
    source = pd.read_csv(analysis.REPO_ROOT / "data/protein/aggregated_all/s3.csv", index_col=0)
    rows, raw_cache = [], {}
    for edge in edges.nlargest(count, "protein_norm").itertuples():
        protein = edge.dominant_protein
        if protein not in raw_cache:
            path = analysis.REPO_ROOT / "data/protein" / f"{protein}.csv"
            raw_cache[protein] = pd.read_csv(path) if path.exists() else None
        raw = raw_cache[protein]
        for role, cell in [("parent", edge.represented_parent), ("child", edge.cell)]:
            subset = raw[raw["Cell-name"].eq(cell)] if raw is not None else None
            present = subset is not None and len(subset) > 0
            rows.append({"tree_id": edge.tree_id, "protein_norm": edge.protein_norm,
                         "role": role, "cell": cell, "protein": protein,
                         "dominant_protein_sq_share": edge.dominant_protein_sq_share,
                         "atlas_value": source.loc[cell, protein],
                         "raw_rows": len(subset) if present else 0,
                         "raw_time_min": subset.Time.min() if present else np.nan,
                         "raw_time_max": subset.Time.max() if present else np.nan,
                         "raw_adjusted_mean": subset["Adjustment-expression"].mean() if present else np.nan,
                         "note": "Reporter timeline is not aligned to tracking; no preprocessing changed"})
    return pd.DataFrame(rows)


def run(context, n_draws=10_000, cv_draws=2_000, seed=242):
    edges = edge_manifest(context)
    ids, parents = edges.tree_id.to_numpy(), edges.parent_tree_id.to_numpy()
    mapping = np.asarray(context.tree.lineage_id_mapping)
    delta = context.optimization.exp_mat[mapping[ids]] - context.optimization.exp_mat[mapping[parents]]
    observed = np.linalg.norm(delta, axis=1)
    # Reuse the official sampler for unchanged travel and its topology audit.
    reference = analysis.separate_clock_reference(context, n_samples=n_draws, seed=seed)
    np.testing.assert_array_equal(ids, reference.evaluated_tree_ids)
    np.testing.assert_allclose(observed.sum(), context.optimization.lineage_exp_cost)
    covariance, tau, fit = fit_model(delta)
    np.testing.assert_allclose(covariance, reference.covariance[3:, 3:])
    norms = paired_norms(covariance, tau, len(ids), n_draws, seed)
    terminal = edges.terminal_in_measured_tree.to_numpy(bool)
    diagnostics = diagnostic_rows(norms, observed, {"internal": ~terminal, "terminal": terminal},
                                  draws=n_draws, seed=seed)
    fits = [{"scope": "all", **fit}]
    draws, quantiles = [], []
    probabilities = np.linspace(0, 1, 101)
    expectation = gaussian_mean_norm(covariance)
    for model, values in norms.items():
        draws.append(pd.DataFrame({"model": model, "replicate": np.arange(n_draws), "seed": seed,
                                   "travel_cost": reference.travel_cost,
                                   "cell_state_cost": values.sum(axis=1)}))
        predicted = len(ids) * expectation * (np.exp(-tau*tau/8) if model == MIXTURE else 1)
        fits[0]["analytic_total_" + ("mixture" if model == MIXTURE else "gaussian")] = predicted
        for q, simulated in zip(probabilities, np.quantile(values, probabilities, axis=1)):
            quantiles.append({"model": model, "quantile": q, "observed": np.quantile(observed, q),
                              **_summarize(simulated)})
    cv, predictive, fold_manifest = [], [], []
    edges["test_fold"] = "stem_not_tested"
    for index, (label, train, test) in enumerate(subtree_folds(context, edges)):
        print(f"Held-out subtree {label}: {train.sum()} train / {test.sum()} test", flush=True)
        cov, fold_tau, info = fit_model(delta[train])
        fits.append({"scope": f"holdout_{label}", **info})
        fold_manifest.append({"fold": label, "n_train": int(train.sum()),
                              "n_test": int(test.sum()), "n_purged": int((~(train | test)).sum()),
                              "seed": seed+100+index, "draws": cv_draws})
        edges.loc[test, "test_fold"] = label
        q, logdet = radial_coordinates(delta[test], cov)
        for model, value in [(GAUSSIAN, 0), (MIXTURE, fold_tau)]:
            ll = radial_logpdf(q, delta.shape[1], logdet, value)
            check = radial_logpdf(q, delta.shape[1], logdet, value, order=1024)
            if np.max(np.abs(ll-check)) > 1e-5:
                raise RuntimeError("Held-out quadrature did not converge")
            for tree_id, logpdf in zip(ids[test], ll):
                predictive.append({"fold": label, "tree_id": tree_id, "model": model, "logpdf": logpdf})
        sim = paired_norms(cov, fold_tau, int(test.sum()), cv_draws, seed+100+index)
        cv.extend(diagnostic_rows(sim, observed[test], fold=label, n_test=int(test.sum()),
                                  draws=cv_draws, seed=seed+100+index))
    # Missing-source sensitivity changes fitting only, not the evaluation scope.
    complete = edges.source_missing_feature_count.eq(0).to_numpy()
    cov_complete, tau_complete, complete_fit = fit_model(delta[complete])
    fits.append({"scope": "fit_excluding_missing_source_edges", **complete_fit,
                 "analytic_total_mixture": len(ids)*gaussian_mean_norm(cov_complete)*np.exp(-tau_complete**2/8),
                 "analytic_total_gaussian": len(ids)*gaussian_mean_norm(cov_complete)})
    audit = pd.DataFrame([{
        "scope": name, "n_edges": len(group),
        "missing_source_edges": int(group.source_missing_feature_count.gt(0).sum()),
        "squared_change_share": (group.protein_norm**2).sum()/(observed**2).sum(),
        "median_dominant_feature_sq_share": group.dominant_protein_sq_share.median(),
    } for name, group in [("all", edges), ("largest_raw_50", edges.nlargest(50, "protein_norm")),
                          ("missing_source", edges[~complete])]])
    summary = []
    for model, values in norms.items():
        for objective, totals, natural in [
            ("protein", values.sum(axis=1), observed.sum()),
            ("travel", reference.travel_cost, context.optimization.lineage_xyz_cost),
        ]:
            summary.append({"model": model, "objective": objective, "draws": n_draws, "seed": seed,
                            "observed": natural, **_summarize(totals),
                            "relative_excess": totals.mean()/natural-1,
                            "interpretation": "conditional plug-in simulation, not a calibrated significance test"})
    return {"fit.csv": pd.DataFrame(fits), "draws.csv": pd.concat(draws, ignore_index=True),
            "summary.csv": pd.DataFrame(summary), "fold_manifest.csv": pd.DataFrame(fold_manifest),
            "diagnostics.csv": pd.DataFrame(diagnostics), "edge_quantiles.csv": pd.DataFrame(quantiles),
            "cross_validation.csv": pd.DataFrame(cv), "heldout_logpdf.csv": pd.DataFrame(predictive),
            "edges.csv": edges, "raw_outlier_windows.csv": raw_audit(context, edges),
            "source_audit.csv": audit,
            "protein_covariance.csv": pd.DataFrame(covariance, columns=context.feature_names[3:]).assign(feature=context.feature_names[3:])}


def plot(outputs, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    colors = {GAUSSIAN: "#2b8c8a", MIXTURE: "#b46531"}
    for model, values in outputs["draws.csv"].groupby("model", sort=False):
        axes[0, 0].hist(values.cell_state_cost, bins=45, density=True, alpha=.55,
                        color=colors[model], label=model)
        q = outputs["edge_quantiles.csv"].query("model == @model")
        axes[0, 1].plot(q["mean"], q["quantile"], color=colors[model])
        axes[0, 1].fill_betweenx(q["quantile"], q.q025, q.q975, color=colors[model], alpha=.15)
        cv = outputs["cross_validation.csv"].query("model == @model and metric == 'total'")
        offset = -.12 if model == GAUSSIAN else .12
        axes[1, 0].errorbar(np.arange(len(cv))+offset, cv["mean"]/cv.observed,
                            yerr=[(cv["mean"]-cv.q025)/cv.observed,
                                  (cv.q975-cv["mean"])/cv.observed],
                            fmt="o", capsize=3, color=colors[model])
    observed = outputs["edges.csv"].protein_norm.sum()
    axes[0, 0].axvline(observed, color="black", label="Observed")
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].set(xlabel="Total protein-state distance", ylabel="Density", title="Full-data plug-in predictions")
    axes[0, 1].plot(q.observed, q["quantile"], color="black")
    axes[0, 1].set(xlabel="Protein edge distance", ylabel="Cumulative fraction",
                    title="Edge quantiles: central 95% simulation bands")
    axes[1, 0].axhline(1, color="black", linestyle="--")
    axes[1, 0].set_xticks(np.arange(len(cv)), cv.fold)
    axes[1, 0].set(ylabel="Predicted / observed total", xlabel="Held-out subtree",
                    title="Held-out totals: central 95% simulation intervals")
    score = outputs["heldout_logpdf.csv"].pivot(index=["fold", "tree_id"], columns="model", values="logpdf")
    gain = (score[MIXTURE]-score[GAUSSIAN]).groupby("fold", sort=False).mean()
    gain = gain.reindex(cv.fold)
    axes[1, 1].bar(gain.index, gain.values, color=colors[MIXTURE])
    axes[1, 1].axhline(0, color="black", linewidth=.8)
    axes[1, 1].set(xlabel="Held-out subtree", ylabel="Mean log-density gain (nats / edge)",
                    title="Positive values favor branch variability")
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Branch-variability sensitivity — Figure 8 unchanged")
    fig.savefig(output_dir / "branch_variability_comparison.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=10_000)
    parser.add_argument("--cv-draws", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=242)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()
    if min(args.draws, args.cv_draws) < 2:
        parser.error("draw counts must be at least two")
    outputs = run(analysis.load_ce_protein_context(max_workers=1), args.draws, args.cv_draws, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(args.output_dir / name, index=False)
    if not args.analysis_only:
        plot(outputs, args.output_dir)
    print(outputs["fit.csv"].to_string(index=False))
    print(outputs["diagnostics.csv"].to_string(index=False))
    print(f"Sensitivity outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
