"""Compare the historical Brownian model with the publication separate-clock reference.

Run from the repository root with conda run -n dev python
full_tree_pareto/clock_sensitivity.py. This writes separate sensitivity outputs;
it does not replace Figure 8 or its publication caches.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from . import publication_analysis as analysis
except ImportError:
    import publication_analysis as analysis


BASELINE = "Tracking-time Brownian (full covariance)"
SEPARATE = "Tracking-time spatial / per-transition protein"
DEFAULT_OUTPUT = analysis.CACHE_ROOT / "clock_sensitivity"


def edge_manifest(context: analysis.FullTreeContext) -> pd.DataFrame:
    """Observed increments, clock metadata and upstream missing-value flags."""
    nodes = analysis.node_manifest(context).set_index("tree_id")
    ids = np.sort(np.concatenate([layer.child_tree_ids for layer in context.layers]))
    parents = np.asarray(context.tree.parent_list)[ids]
    mapping = np.asarray(context.tree.lineage_id_mapping)
    values = np.hstack((context.optimization.xyz_mat, context.optimization.exp_mat))[mapping]
    delta = values[ids] - values[parents]
    frame = nodes.loc[ids].reset_index()
    frame["parent_tree_id"] = parents
    frame["branch_time"] = np.asarray(context.tree.branch_time_length)[ids]
    frame["canonical_generation_gap"] = (
        nodes.loc[ids, "canonical_depth"].to_numpy()
        - nodes.loc[parents, "canonical_depth"].to_numpy()
    )
    frame["at_tracking_cutoff"] = frame["tracking_time"] == analysis.TRACKING_CUTOFF
    frame["travel_norm"] = np.linalg.norm(delta[:, :3], axis=1)
    frame["protein_norm"] = np.linalg.norm(delta[:, 3:], axis=1)
    frame["protein_time_rate_sq"] = frame["protein_norm"] ** 2 / frame["branch_time"]
    frame["protein_transition_rate_sq"] = frame["protein_norm"] ** 2
    features = context.feature_names[3:]
    source = pd.read_csv(
        analysis.REPO_ROOT / "data/protein/aggregated_all/s3.csv", index_col=0
    )[features]
    child_missing = source.loc[frame["cell"]].isna().to_numpy()
    parent_missing = source.loc[frame["represented_parent"]].isna().to_numpy()
    frame["source_missing_feature_count"] = (child_missing | parent_missing).sum(axis=1)
    contribution = delta[:, 3:] ** 2
    top_feature = contribution.argmax(axis=1)
    frame["dominant_protein"] = np.asarray(features)[top_feature]
    frame["dominant_protein_sq_share"] = np.divide(
        contribution[np.arange(len(ids)), top_feature], contribution.sum(axis=1),
        out=np.zeros(len(ids)), where=contribution.sum(axis=1) > 0,
    )
    return frame


def observation_window_audit(context, edges, count=12):
    """Read dominant-protein raw windows for influential edges, without aligning clocks.

    Raw reporter Time and tracking t belong to different observation timelines;
    their endpoints are recorded as provenance, never substituted as durations.
    """
    source = pd.read_csv(
        analysis.REPO_ROOT / "data/protein/aggregated_all/s3.csv", index_col=0
    )
    rows = []
    raw_by_protein = {}
    for edge in edges.nlargest(count, "protein_time_rate_sq").itertuples():
        protein = edge.dominant_protein
        if protein not in raw_by_protein:
            path = analysis.REPO_ROOT / "data/protein" / f"{protein}.csv"
            raw_by_protein[protein] = pd.read_csv(path) if path.exists() else None
        raw = raw_by_protein[protein]
        for role, cell in [("parent", edge.represented_parent), ("child", edge.cell)]:
            subset = raw[raw["Cell-name"] == cell] if raw is not None else None
            present = subset is not None and len(subset) > 0
            rows.append({
                "child_tree_id": edge.tree_id, "role": role, "cell": cell,
                "protein": protein, "tracking_last_time": context.tracking_times[cell],
                "edge_tracking_duration": edge.branch_time,
                "atlas_value": source.loc[cell, protein],
                "atlas_value_missing": pd.isna(source.loc[cell, protein]),
                "raw_rows": len(subset) if present else 0,
                "raw_time_min": subset["Time"].min() if present else np.nan,
                "raw_time_max": subset["Time"].max() if present else np.nan,
                "raw_adjusted_mean": subset["Adjustment-expression"].mean() if present else np.nan,
                "clock_note": "Reporter Time is not aligned to tracking t",
            })
    return pd.DataFrame(rows)


def _summarize(values):
    return {
        "mean": float(np.mean(values)), "std": float(np.std(values)),
        "q025": float(np.quantile(values, .025)),
        "q975": float(np.quantile(values, .975)),
    }


def _top_share(squared):
    k = max(1, int(np.ceil(squared.shape[-1] * .05)))
    return np.sort(squared, axis=-1)[..., -k:].sum(axis=-1) / squared.sum(axis=-1)


def compare(context, n_samples=10_000, seed=242):
    edges = edge_manifest(context)
    frames, summaries, diagnostics, distributions, covariances = [], [], [], [], []
    groups = {
        "all": np.ones(len(edges), dtype=bool),
        "internal": ~edges["terminal_in_measured_tree"].to_numpy(bool),
        "terminal": edges["terminal_in_measured_tree"].to_numpy(bool),
        "cutoff": edges["at_tracking_cutoff"].to_numpy(bool),
        "before_cutoff": ~edges["at_tracking_cutoff"].to_numpy(bool),
    }
    for model, sampler in [
        (BASELINE, analysis.parametric_brownian_reference),
        (SEPARATE, analysis.separate_clock_reference),
    ]:
        print(f"Simulating {n_samples:,} draws: {model}", flush=True)
        ref = sampler(context, n_samples=n_samples, seed=seed, capture_edge_norms=True)
        np.testing.assert_array_equal(ref.evaluated_tree_ids, edges["tree_id"])
        frame = pd.DataFrame({
            "model": model, "replicate": np.arange(n_samples), "seed": seed,
            "travel_cost": ref.travel_cost, "cell_state_cost": ref.cell_state_cost,
        })
        frames.append(frame)
        covariance = pd.DataFrame(ref.covariance, index=context.feature_names,
                                  columns=context.feature_names).rename_axis("feature").reset_index()
        covariance.insert(0, "model", model)
        covariance.insert(1, "spatial_clock", "tracking time")
        covariance.insert(2, "protein_clock", "tracking time" if model == BASELINE else "one transition")
        covariances.append(covariance)
        for objective, norms, observed, total in [
            ("travel", ref.edge_travel, edges["travel_norm"].to_numpy(), ref.travel_cost),
            ("protein", ref.edge_cell_state, edges["protein_norm"].to_numpy(), ref.cell_state_cost),
        ]:
            np.testing.assert_allclose(norms.sum(axis=1), total, rtol=1e-12)
            summaries.append({
                "model": model, "objective": objective, "draws": n_samples, "seed": seed,
                "observed": observed.sum(), **_summarize(total),
                "relative_excess": total.mean() / observed.sum() - 1,
                "observed_lower_tail_count": int((total <= observed.sum()).sum()),
                "interpretation": "plug-in sensitivity; tail counts are not calibrated p-values",
            })
            divisor = (edges["branch_time"].to_numpy()
                       if objective == "travel" or model == BASELINE else np.ones(len(edges)))
            measures = {
                "edge_norm_cv": (norms.std(axis=1) / norms.mean(axis=1), observed.std() / observed.mean()),
                "max_edge_norm": (norms.max(axis=1), observed.max()),
                "top_5pct_raw_sq_share": (_top_share(norms**2), _top_share(observed**2)),
                "top_5pct_model_rate_sq_share": (
                    _top_share(norms**2 / divisor), _top_share(observed**2 / divisor)),
            }
            for group, mask in groups.items():
                measures[f"mean_edge_norm_{group}"] = (norms[:, mask].mean(axis=1), observed[mask].mean())
            for metric, (simulated, natural) in measures.items():
                diagnostics.append({"model": model, "objective": objective, "metric": metric,
                                    "observed": natural, **_summarize(simulated)})
            # Distribution of within-replicate edge quantiles, not independent
            # pseudo-replicates pooled across every simulated edge.
            for q, simulated in zip(np.linspace(0, 1, 101), np.quantile(norms, np.linspace(0, 1, 101), axis=1)):
                distributions.append({"model": model, "objective": objective,
                                      "quantile": q, "observed": np.quantile(observed, q),
                                      **_summarize(simulated)})
    return {
        "draws.csv": pd.concat(frames, ignore_index=True),
        "summary.csv": pd.DataFrame(summaries),
        "diagnostics.csv": pd.DataFrame(diagnostics),
        "edge_quantiles.csv": pd.DataFrame(distributions),
        "covariances.csv": pd.concat(covariances, ignore_index=True),
        "edges.csv": edges,
        "observation_windows.csv": observation_window_audit(context, edges),
    }


def plot_comparison(outputs, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {BASELINE: "#917ca8", SEPARATE: "#2b8c8a"}
    labels = {BASELINE: "Historical Brownian", SEPARATE: "Separate clocks"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    draws, summary = outputs["draws.csv"], outputs["summary.csv"]
    for model, group in draws.groupby("model", sort=False):
        shown = group.iloc[np.linspace(0, len(group)-1, min(1000, len(group)), dtype=int)]
        axes[0, 0].scatter(shown.travel_cost, shown.cell_state_cost, s=5, alpha=.2,
                           color=colors[model], label=labels[model])
        axes[0, 1].hist(group.cell_state_cost, bins=40, density=True, alpha=.5,
                        color=colors[model], label=labels[model])
        q = outputs["edge_quantiles.csv"].query("model == @model and objective == 'protein'")
        axes[1, 0].plot(q["mean"], q["quantile"], color=colors[model], label=labels[model])
        axes[1, 0].fill_betweenx(q["quantile"], q.q025, q.q975, color=colors[model], alpha=.15)
        diag = outputs["diagnostics.csv"].query("model == @model and objective == 'protein'")
        groups = ["internal", "terminal", "cutoff", "before_cutoff"]
        d = diag.set_index("metric").loc[[f"mean_edge_norm_{g}" for g in groups]]
        offset = -.12 if model == BASELINE else .12
        axes[1, 1].errorbar(np.arange(4)+offset, d["mean"],
                            yerr=[d["mean"]-d.q025, d.q975-d["mean"]], fmt="o",
                            color=colors[model], capsize=3, label=labels[model])
    observed = summary.groupby("objective")["observed"].first()
    axes[0, 0].scatter(observed["travel"], observed["protein"], marker="x", s=90,
                       color="black", linewidths=2, label="Observed")
    axes[0, 1].axvline(observed["protein"], color="black", label="Observed")
    axes[1, 0].plot(q.observed, q["quantile"], color="black", label="Observed")
    axes[1, 1].scatter(np.arange(4), d.observed, color="black", marker="x", s=45, label="Observed")
    axes[1, 1].set_xticks(np.arange(4), ["Internal", "Terminal", "At cutoff", "Before cutoff"])
    specs = [
        ("Total travel cost", "Total protein-state cost", "Fixed-topology reference clouds"),
        ("Total protein-state cost", "Density", "Protein totals"),
        ("Protein edge length", "Cumulative fraction", "Protein edge distributions (95% replicate bands)"),
        ("Edge group (groups overlap)", "Mean protein edge length", "Grouped means (95% simulation intervals)"),
    ]
    for ax, (xlabel, ylabel, title) in zip(axes.flat, specs):
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(markerscale=1.5, fontsize=8)
    fig.suptitle("Clock sensitivity: measured topology, four roots, 1,000 edges")
    for ext in ["pdf", "png"]:
        fig.savefig(output_dir / f"clock_comparison.{ext}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=242)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()
    if args.draws < 2:
        parser.error("--draws must be at least 2 for comparison diagnostics")
    outputs = compare(analysis.load_ce_protein_context(), args.draws, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(args.output_dir / name, index=False)
    if not args.analysis_only:
        plot_comparison(outputs, args.output_dir)
    print(outputs["summary.csv"].to_string(index=False))
    print(f"Sensitivity outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
