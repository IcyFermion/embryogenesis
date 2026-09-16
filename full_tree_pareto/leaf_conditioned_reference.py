"""Side comparison: root-and-leaf-conditioned separate-clock Gaussian.

Run from the repository root: conda run -n dev python -m
full_tree_pareto.leaf_conditioned_reference. No publication inputs are changed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import splu

from full_tree_pareto import publication_analysis as analysis
from full_tree_pareto.clock_sensitivity import edge_manifest, _summarize, _top_share


GAUSSIAN = "Root-fixed Gaussian"
CONDITIONED = "Root-and-leaf-fixed Gaussian"
DEFAULT_OUTPUT = analysis.CACHE_ROOT / "leaf_conditioned"


class ConditionalBlock:
    """Exact Gaussian conditioning using the weighted tree Laplacian.

    Nodes use compact local indices. For incidence B and W=diag(1/clock),
    Q=B_U.T W B_U, M=-Q^-1 B_U.T W B_F X_F. The free-node covariance is
    Q^-1 (node factor) times Sigma (feature factor). Sampling solves
    Q H=B_U.T sqrt(W) Z; this gives Cov(H)=Q^-1 without a dense square root.
    """

    def __init__(self, parents, children, clock, values, fixed, covariance):
        self.parents = np.asarray(parents, dtype=int)
        self.children = np.asarray(children, dtype=int)
        self.clock = np.asarray(clock, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.fixed = np.asarray(fixed, dtype=bool)
        self.covariance = np.asarray(covariance, dtype=float)
        n, d = self.values.shape
        e = len(self.children)
        if (self.fixed.shape != (n,) or self.clock.shape != (e,)
                or self.parents.shape != (e,) or self.covariance.shape != (d, d)
                or not np.isfinite(self.values).all()
                or not np.isfinite(self.clock).all() or np.any(self.clock <= 0)):
            raise ValueError("Invalid conditional block dimensions, values or clocks")
        if not self.fixed.any() or self.fixed.all():
            raise ValueError("At least one fixed and one free node required")
        self.factor = np.linalg.cholesky(self.covariance)
        rows = np.repeat(np.arange(e), 2)
        columns = np.column_stack((self.parents, self.children)).ravel()
        incidence = sparse.csr_matrix((np.tile([-1., 1.], e), (rows, columns)), shape=(e, n))
        self.free = ~self.fixed
        self.bu = incidence[:, self.free]
        bf = incidence[:, self.fixed]
        weighted_bu = self.bu.multiply((1/self.clock)[:, None])
        self.precision = (self.bu.T @ weighted_bu).tocsc()
        self.solver = splu(self.precision)
        rhs = -self.bu.T @ ((bf @ self.values[self.fixed])/self.clock[:, None])
        self.mean = self.values.copy()
        self.mean[self.free] = self.solver.solve(np.asarray(rhs))
        self.node_covariance = self.solver.solve(np.eye(self.free.sum()))
        self.edge_variance = np.asarray(
            self.bu.multiply(self.bu @ self.node_covariance).sum(axis=1)
        ).ravel()
        if np.any(self.edge_variance < -1e-10) or np.any(self.edge_variance > self.clock+1e-9):
            raise AssertionError("Conditioning must not increase marginal edge variance")
        np.testing.assert_allclose(self.precision @ self.mean[self.free], rhs, atol=1e-9)

    def sample(self, n_draws, rng):
        if n_draws < 1:
            raise ValueError("Positive draw count required")
        d = self.values.shape[1]
        z = rng.normal(size=(len(self.children), n_draws*d))
        rhs = self.bu.T @ (z / np.sqrt(self.clock)[:, None])
        noise = self.solver.solve(np.asarray(rhs)).reshape(self.free.sum(), n_draws, d)
        noise = noise.transpose(1, 0, 2) @ self.factor.T
        states = np.broadcast_to(self.mean, (n_draws, *self.mean.shape)).copy()
        states[:, self.free] += noise
        # Fixed values are copied verbatim, never snapped after an unconditioned draw.
        np.testing.assert_array_equal(states[:, self.fixed],
                                      np.broadcast_to(self.values[self.fixed], states[:, self.fixed].shape))
        return states

    def analytic_squared_total(self):
        increments = self.mean[self.children] - self.mean[self.parents]
        return float((increments**2).sum() + self.edge_variance.sum()*np.trace(self.covariance))


def prepare(context, covariance):
    edges = edge_manifest(context)
    roots = np.asarray(sorted(tree_id for tree_id, _ in context.optimization.first_internal_layer))
    ids = np.sort(np.concatenate((roots, edges.tree_id.to_numpy())))
    local = {int(tree_id): i for i, tree_id in enumerate(ids)}
    parents = np.array([local[i] for i in edges.parent_tree_id])
    children = np.array([local[i] for i in edges.tree_id])
    leaves = np.array(sorted(set(children)-set(parents)))
    fixed = np.isin(ids, roots)
    fixed[leaves] = True
    nodes = analysis.node_manifest(context).set_index("tree_id").loc[ids].reset_index()
    np.testing.assert_array_equal(ids[leaves],
                                  np.sort(nodes.loc[nodes.terminal_in_measured_tree, "tree_id"]))
    if (len(ids), len(edges), len(roots), len(leaves), int((~fixed).sum())) != (1004, 1000, 4, 504, 496):
        raise AssertionError("Unexpected publication conditioning scope")
    if not edges.canonical_generation_gap.eq(1).all():
        raise AssertionError("Protein clock requires one generation per edge")
    mapping = np.asarray(context.tree.lineage_id_mapping)[ids]
    values = np.hstack((context.optimization.xyz_mat[mapping], context.optimization.exp_mat[mapping]))
    blocks = {
        "travel": ConditionalBlock(parents, children, edges.branch_time.to_numpy(),
                                   values[:, :3], fixed, covariance[:3, :3]),
        "protein": ConditionalBlock(parents, children, np.ones(len(edges)),
                                    values[:, 3:], fixed, covariance[3:, 3:]),
    }
    nodes["fixed_boundary"] = fixed
    return edges, nodes, blocks


def compare(context, n_draws=10_000, seed=242):
    reference = analysis.separate_clock_reference(context, n_samples=n_draws, seed=seed,
                                                 capture_edge_norms=True)
    edges, nodes, blocks = prepare(context, reference.covariance)
    np.testing.assert_array_equal(edges.tree_id, reference.evaluated_tree_ids)
    streams = np.random.SeedSequence(seed).spawn(2)
    generators = {name: np.random.default_rng(stream) for name, stream in zip(blocks, streams)}
    norms = {name: np.empty((n_draws, len(edges))) for name in blocks}
    for start in range(0, n_draws, 100):
        stop = min(start+100, n_draws)
        for name, block in blocks.items():
            states = block.sample(stop-start, generators[name])
            delta = states[:, block.children] - states[:, block.parents]
            norms[name][start:stop] = np.linalg.norm(delta, axis=-1)
        if stop % 2000 == 0:
            print(f"Conditioned draws: {stop:,}/{n_draws:,}", flush=True)
    baseline = {"travel": reference.edge_travel, "protein": reference.edge_cell_state}
    summary, diagnostics, quantiles, analytic = [], [], [], []
    terminal = edges.terminal_in_measured_tree.to_numpy(bool)
    for objective, block in blocks.items():
        observed = edges[f"{objective}_norm"].to_numpy()
        mean_delta = block.mean[block.children]-block.mean[block.parents]
        mean_norm = np.linalg.norm(mean_delta, axis=1)
        analytic.append({"objective": objective, "conditional_mean_state_cost": mean_norm.sum(),
                         "expected_squared_total": block.analytic_squared_total(),
                         "simulated_squared_total": (norms[objective]**2).sum(axis=1).mean(),
                         "squared_total_mc_se": (norms[objective]**2).sum(axis=1).std()/np.sqrt(n_draws)})
        edges[f"{objective}_conditional_mean_norm"] = mean_norm
        edges[f"{objective}_conditional_variance_factor"] = block.edge_variance
        for model, values in [(GAUSSIAN, baseline[objective]), (CONDITIONED, norms[objective])]:
            totals = values.sum(axis=1)
            summary.append({"model": model, "objective": objective, "draws": n_draws, "seed": seed,
                            "observed": observed.sum(), **_summarize(totals),
                            "relative_excess": totals.mean()/observed.sum()-1,
                            "interpretation": "conditional plug-in reference, not a calibrated significance test"})
            measures = {
                "edge_norm_cv": (values.std(axis=1)/values.mean(axis=1), observed.std()/observed.mean()),
                "max_edge_norm": (values.max(axis=1), observed.max()),
                "top_5pct_raw_sq_share": (_top_share(values**2), _top_share(observed**2)),
                "mean_internal": (values[:, ~terminal].mean(axis=1), observed[~terminal].mean()),
                "mean_terminal": (values[:, terminal].mean(axis=1), observed[terminal].mean()),
            }
            for metric, (sim, natural) in measures.items():
                diagnostics.append({"model": model, "objective": objective, "metric": metric,
                                    "observed": natural, **_summarize(sim)})
            probabilities = np.linspace(0, 1, 101)
            for q, sim in zip(probabilities, np.quantile(values, probabilities, axis=1)):
                quantiles.append({"model": model, "objective": objective, "quantile": q,
                                  "observed": np.quantile(observed, q), **_summarize(sim)})
    mean_states = nodes.copy()
    for name, block in blocks.items():
        features = context.feature_names[:3] if name == "travel" else context.feature_names[3:]
        for j, feature in enumerate(features):
            mean_states[f"conditional_mean_{feature}"] = block.mean[:, j]
    draws = pd.concat([pd.DataFrame({"model": model, "replicate": np.arange(n_draws), "seed": seed,
                                     "travel_cost": data["travel"].sum(axis=1),
                                     "cell_state_cost": data["protein"].sum(axis=1)})
                       for model, data in [(GAUSSIAN, baseline), (CONDITIONED, norms)]], ignore_index=True)
    covariance = pd.DataFrame(reference.covariance, columns=context.feature_names)
    covariance.insert(0, "feature", context.feature_names)
    return {"summary.csv": pd.DataFrame(summary), "diagnostics.csv": pd.DataFrame(diagnostics),
            "edge_quantiles.csv": pd.DataFrame(quantiles), "analytic_checks.csv": pd.DataFrame(analytic),
            "draws.csv": draws, "edges.csv": edges, "conditional_mean_states.csv": mean_states,
            "covariance.csv": covariance}


def plot(outputs, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    colors = {GAUSSIAN: "#2b8c8a", CONDITIONED: "#8061a5"}
    for model, draws in outputs["draws.csv"].groupby("model", sort=False):
        shown = draws.iloc[np.linspace(0, len(draws)-1, min(1000, len(draws)), dtype=int)]
        axes[0, 0].scatter(shown.travel_cost, shown.cell_state_cost, s=6, alpha=.22,
                           color=colors[model], label=model)
        axes[0, 1].hist(draws.travel_cost, bins=45, density=True, alpha=.5, color=colors[model])
        axes[1, 0].hist(draws.cell_state_cost, bins=45, density=True, alpha=.5, color=colors[model])
        q = outputs["edge_quantiles.csv"].query("model == @model and objective == 'protein'")
        axes[1, 1].plot(q["mean"], q["quantile"], color=colors[model])
        axes[1, 1].fill_betweenx(q["quantile"], q.q025, q.q975, color=colors[model], alpha=.15)
    observed = outputs["summary.csv"].groupby("objective")["observed"].first()
    mean_cost = outputs["analytic_checks.csv"].set_index("objective").conditional_mean_state_cost
    axes[0, 0].scatter(observed.travel, observed.protein, c="black", marker="x", s=75, label="Observed")
    axes[0, 0].scatter(mean_cost.travel, mean_cost.protein, c="#b46531", marker="D", s=35,
                       label="Conditional mean states (not minimum cost)")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].axvline(observed.travel, color="black")
    axes[1, 0].axvline(observed.protein, color="black")
    axes[1, 1].plot(q.observed, q["quantile"], color="black")
    specs = [("Travel cost", "Protein-state cost", "Same topology; different endpoint constraints"),
             ("Travel cost", "Density", "Travel totals"),
             ("Protein-state cost", "Density", "Protein totals"),
             ("Protein edge distance", "Cumulative fraction", "Edge quantiles: central 95% simulation bands")]
    for ax, (xlabel, ylabel, title) in zip(axes.flat, specs):
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Leaf-conditioned separate-clock Gaussian — Figure 8 unchanged")
    fig.savefig(output_dir / "leaf_conditioned_comparison.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=242)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()
    if args.draws < 2:
        parser.error("--draws must be at least two")
    outputs = compare(analysis.load_ce_protein_context(max_workers=1), args.draws, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(args.output_dir / name, index=False)
    if not args.analysis_only:
        plot(outputs, args.output_dir)
    print(outputs["summary.csv"].to_string(index=False))
    print(outputs["analytic_checks.csv"].to_string(index=False))
    print(f"Side-analysis outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
