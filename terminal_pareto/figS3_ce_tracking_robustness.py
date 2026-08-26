"""Three-replicate C. elegans tracking robustness analysis and figure.

This pipeline is intentionally separate from the accepted Figures 2--6
renderers.  It compares the three stage-matched tracking replicates using a
strict intersection of natural terminal parent--child edges, while holding
the lineage tree, top-20 protein representation, null construction, and
Pareto sweep fixed.

Outputs
-------
Analysis caches in ``terminal_pareto/output/``:

* ``ce_tracking_replicate_audit.csv``
* ``ce_tracking_replicate_cell_manifest.csv``
* ``ce_tracking_replicate_subtree_manifest.csv``
* ``ce_tracking_replicate_fronts.csv``
* ``ce_tracking_replicate_null_clouds.csv``
* ``ce_tracking_replicate_major_metrics.csv``

Publication assets in ``terminal_pareto/output/publication/``:

* ``figS3A_ce_tracking_replicate_fronts.{pdf,png}``
* ``figS3B_ce_tracking_replicate_metrics.{pdf,png}``

Usage::

    python terminal_pareto/figS3_ce_tracking_robustness.py
    python terminal_pareto/figS3_ce_tracking_robustness.py --audit-only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto import plot_style as ps
from terminal_pareto.fig5_table1_ce_canonical_metrics import (
    compute_front_landmarks,
)
from terminal_pareto.subtree_analysis import first_cousin_null_summary


ANALYSIS_OUT = Path(__file__).resolve().parent / "output"
PUBLICATION_OUT = ANALYSIS_OUT / "publication"
MAJOR_SUBTREES = ("P0", "AB", "ABa", "ABp", "P1")
MIN_CELLS = 12
# ``compute_std_scaled_pareto`` includes both endpoints, so 300 intervals
# produce 301 sampled weights.
ITERATION = 300
NULL_DRAWS = 1000
EXPECTED_MATCHED_EDGES = 279
EXPECTED_COMMON_SUBTREES = 42

# Replicate is a local categorical encoding. Recurrent biological quantities
# retain the manuscript-wide semantic colors below.
REPLICATE_STYLES = {
    "embryo1": dict(color="#4477AA", ls="-"),
    "embryo2": dict(color="#EE6677", ls=(0, (5, 2))),
    "embryo3": dict(color="#228833", ls=(0, (1.5, 1.5))),
}
NULL_COLOR = ps.SEMANTIC_COLORS["first_cousin_null"]
RETENTION_COLOR = ps.COLORS["blue"]
INK = ps.COLORS["black"]


def _tree_nodes(tree_index):
    return set(tree_index)


def _major_region(name, tree_index):
    """Major lineage containing ``name``; used only as manifest metadata."""
    if name == "P0":
        return "P0"
    current = name
    while current not in {"AB", "ABa", "ABp", "P1"}:
        info = tree_index.get(current)
        if info is None or info["parent"] is None:
            return "P0"
        current = info["parent"]
    return current


def audit_replicates():
    """Load predefined stage-matched replicates and enumerate terminal edges."""
    lineage = dl.load_json(str(REPO_ROOT / "data" / "cell_lineage.json"))
    tree_index = lm.build_lineage_tree_index(lineage)
    protein = dl.load_protein_expression()
    selected = dl.load_prot_sel()
    protein_names = set(protein.index)
    replicate_data = {}
    audit_rows = []

    for label, path_string, cutoff in dl.CE_REPLICATES:
        path = Path(path_string)
        raw = pd.read_csv(path, sep="\t")
        xyz_map, valid_names = dl.load_elegans_tracking(cutoff, path=path_string)
        expression_valid = [name for name in valid_names if name in protein_names]
        terminal_nodes, terminal_parents = dl.collect_terminals(
            lineage, expression_valid
        )
        terms = list(zip(terminal_nodes, terminal_parents))
        if len(terms) != len(set(terms)):
            raise AssertionError(f"Duplicate terminal edges in {label}")

        def terminal_count_at(candidate_cutoff):
            """Loader-equivalent count from the already loaded track table."""
            through_cutoff = raw[raw["t"] <= candidate_cutoff]
            time_summary = through_cutoff.groupby("name")["t"].agg(
                ["count", "first"]
            )
            present = set(time_summary.index[
                (time_summary["count"] != 1)
                | (time_summary["first"] != candidate_cutoff)
            ])
            expression_present = present & protein_names
            candidate_nodes, _ = dl.collect_terminals(
                lineage, expression_present
            )
            return len(candidate_nodes)

        previous_count = terminal_count_at(cutoff - 1)
        selected_count = terminal_count_at(cutoff)
        next_count = terminal_count_at(cutoff + 1)
        if selected_count != len(terms):
            raise AssertionError(
                f"Raw timing audit disagrees with loader for {label}"
            )
        replicate_data[label] = dict(
            cutoff=cutoff,
            xyz_map=xyz_map,
            term_set=set(terms),
        )
        audit_rows.append(dict(
            replicate=label,
            source_path=str(path.relative_to(REPO_ROOT)),
            cutoff=cutoff,
            source_time_min=int(raw["t"].min()),
            source_time_max=int(raw["t"].max()),
            source_unique_cells=int(raw["name"].nunique()),
            valid_cells_at_cutoff=len(valid_names),
            expression_valid_cells=len(expression_valid),
            terminal_edges_available=len(terms),
            terminal_edges_previous_frame=previous_count,
            terminal_edges_next_frame=next_count,
        ))

    term_sets = [d["term_set"] for d in replicate_data.values()]
    common_terms = set.intersection(*term_sets)
    union_terms = set.union(*term_sets)
    if not common_terms:
        raise ValueError("The three tracking replicates have no shared terminal edges")

    audit = pd.DataFrame(audit_rows)
    audit["matched_terminal_edges"] = len(common_terms)
    audit["excluded_for_strict_match"] = (
        audit["terminal_edges_available"] - len(common_terms)
    )
    audit["matched_fraction"] = (
        len(common_terms) / audit["terminal_edges_available"]
    )
    audit["cutoff_basis"] = "predefined stage-matched cutoff"

    context = dict(
        lineage=lineage,
        tree_index=tree_index,
        protein=protein,
        selected=selected,
        replicates=replicate_data,
        common_terms=common_terms,
        union_terms=union_terms,
    )
    return context, audit


def build_manifests(context, audit):
    """Write edge-level and subtree-level availability/matching manifests."""
    replicate_data = context["replicates"]
    tree_index = context["tree_index"]
    common_terms = context["common_terms"]
    union_terms = context["union_terms"]

    cell_rows = []
    for cell, parent in sorted(union_terms):
        row = dict(
            terminal_cell=cell,
            natural_parent=parent,
            lineage_region=_major_region(cell, tree_index),
        )
        for label, data in replicate_data.items():
            row[f"present_{label}"] = (cell, parent) in data["term_set"]
        row["present_all_replicates"] = (cell, parent) in common_terms
        cell_rows.append(row)
    cells = pd.DataFrame(cell_rows)

    # Traverse the complete reference tree once; exact edge-set intersections
    # then avoid manufacturing an edge merely because its endpoints occur in
    # different replicate-specific subsets.
    all_subtrees = dl.collect_all_subtrees(
        context["lineage"], _tree_nodes(tree_index), min_cells=1
    )
    subtree_rows = []
    for name, terms in all_subtrees:
        descendant_edges = set(terms)
        counts = {
            label: len(descendant_edges & data["term_set"])
            for label, data in replicate_data.items()
        }
        common_n = len(descendant_edges & common_terms)
        if not (common_n >= MIN_CELLS
                or max(counts.values(), default=0) >= MIN_CELLS
                or name in MAJOR_SUBTREES):
            continue
        row = dict(
            subtree=name,
            lineage_region=_major_region(name, tree_index),
            reference_terminal_edges=len(terms),
            **{f"n_{label}": count for label, count in counts.items()},
            n_common=common_n,
            eligible_common_min12=common_n >= MIN_CELLS,
            major_subtree=name in MAJOR_SUBTREES,
        )
        subtree_rows.append(row)
    subtrees = pd.DataFrame(subtree_rows).sort_values(
        ["n_common", "subtree"], ascending=[False, True]
    )

    ANALYSIS_OUT.mkdir(parents=True, exist_ok=True)
    audit.to_csv(ANALYSIS_OUT / "ce_tracking_replicate_audit.csv", index=False)
    cells.to_csv(
        ANALYSIS_OUT / "ce_tracking_replicate_cell_manifest.csv", index=False
    )
    subtrees.to_csv(
        ANALYSIS_OUT / "ce_tracking_replicate_subtree_manifest.csv", index=False
    )
    return cells, subtrees


def _ordered_matched_terms(context):
    """Strictly matched edges in reference-lineage traversal order."""
    terminal_nodes, terminal_parents = dl.collect_terminals(
        context["lineage"], _tree_nodes(context["tree_index"])
    )
    ordered = [term for term in zip(terminal_nodes, terminal_parents)
               if term in context["common_terms"]]
    if set(ordered) != context["common_terms"]:
        raise AssertionError("Reference traversal did not recover the common edge set")
    return ordered


def _subtree_term_map(context, ordered_common):
    """Map major subtree names to matched descendant edges in common order."""
    reference_subtrees = dict(dl.collect_all_subtrees(
        context["lineage"], _tree_nodes(context["tree_index"]), min_cells=1
    ))
    out = {}
    for name in MAJOR_SUBTREES:
        descendants = set(reference_subtrees[name])
        out[name] = [term for term in ordered_common if term in descendants]
    return out


def run_analysis(context, iteration=ITERATION, null_draws=NULL_DRAWS):
    """Compute matched global fronts and major-subtree canonical metrics."""
    ordered_common = _ordered_matched_terms(context)
    subtree_terms = _subtree_term_map(context, ordered_common)
    gp_map = pe.build_grandparent_map(context["lineage"])
    front_rows = []
    null_rows = []
    metric_rows = []
    expression_reference = {}

    for label, data in context["replicates"].items():
        for subtree in MAJOR_SUBTREES:
            terms = subtree_terms[subtree]
            tn = [cell for cell, _parent in terms]
            tp = [parent for _cell, parent in terms]
            xm, em, _ = pe.build_cost_matrices(
                tn, tp, data["xyz_map"], context["protein"], context["selected"]
            )
            if subtree in expression_reference:
                if not np.array_equal(em, expression_reference[subtree]):
                    raise AssertionError(
                        f"Expression geometry changed across replicates for {subtree}"
                    )
            else:
                expression_reference[subtree] = em.copy()

            null = first_cousin_null_summary(
                xm, em, tn, gp_map, n_random=null_draws, seed=42
            )
            xa, ea, edge, kp = pe.compute_std_scaled_pareto(
                xm, em, tp, null["_raw"], iteration=iteration
            )
            landmarks = compute_front_landmarks(
                xa, ea, edge,
                nx=null["null_xyz_mean_sigma"],
                ny=null["null_exp_mean_sigma"],
            )
            if not landmarks["endpoint_ok"]:
                raise ValueError(f"Degenerate endpoint geometry for {label}/{subtree}")
            metric_rows.append(dict(
                replicate=label,
                cutoff=data["cutoff"],
                subtree=subtree,
                n=len(terms),
                u_L=landmarks["u_lineage_lp"],
                d_LP=landmarks["d_lp"],
                d_NP=landmarks["d_np"],
                max_edge_retention=landmarks["max_er"],
                natural_on_front=landmarks["exactly_on_front"],
                sampled_front_points=landmarks["n_front_points"],
                null_travel_mean_sigma=null["null_xyz_mean_sigma"],
                null_cell_state_mean_sigma=null["null_exp_mean_sigma"],
                travel_endpoint_span_sigma=landmarks["x_scale"],
                cell_state_endpoint_span_sigma=landmarks["y_scale"],
            ))

            if subtree == "P0":
                max_idx = int(kp["max_er_idx"])
                for i, (x, y, er) in enumerate(zip(xa, ea, edge)):
                    front_rows.append(dict(
                        replicate=label,
                        cutoff=data["cutoff"],
                        n=len(terms),
                        sweep_index=i,
                        alpha_travel=i / iteration,
                        travel_sigma=x,
                        cell_state_sigma=y,
                        edge_retention=er,
                        null_travel_mean_sigma=null["null_xyz_mean_sigma"],
                        null_cell_state_mean_sigma=null["null_exp_mean_sigma"],
                        is_max_retention=i == max_idx,
                        is_cell_state_optimum=i == 0,
                        is_travel_optimum=i == iteration,
                    ))
                raw_null = null["_raw"]
                nx = ((np.asarray(raw_null["random_xyz"])
                       - raw_null["lineage_xyz"]) / raw_null["xyz_std"])
                ny = ((np.asarray(raw_null["random_exp"])
                       - raw_null["lineage_exp"]) / raw_null["exp_std"])
                for i, (x, y) in enumerate(zip(nx, ny)):
                    null_rows.append(dict(
                        replicate=label,
                        draw=i,
                        travel_sigma=x,
                        cell_state_sigma=y,
                    ))

    fronts = pd.DataFrame(front_rows)
    nulls = pd.DataFrame(null_rows)
    metrics = pd.DataFrame(metric_rows)
    fronts.to_csv(
        ANALYSIS_OUT / "ce_tracking_replicate_fronts.csv", index=False
    )
    nulls.to_csv(
        ANALYSIS_OUT / "ce_tracking_replicate_null_clouds.csv", index=False
    )
    metrics.to_csv(
        ANALYSIS_OUT / "ce_tracking_replicate_major_metrics.csv", index=False
    )
    return fronts, nulls, metrics


def _panel_letter(ax, letter, x=-0.13, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", ha="left", va="top", clip_on=False)


def plot_fronts(fronts, nulls, audit):
    """Panel A: matched all-terminal fronts in common null-SD axes."""
    fig, ax = plt.subplots(figsize=(7.15, 4.45))
    fig.subplots_adjust(left=0.105, right=0.79, bottom=0.13, top=0.92)
    ax.axhline(0, color="#777777", lw=0.7, ls=":", zorder=0)
    ax.axvline(0, color="#777777", lw=0.7, ls=":", zorder=0)

    # The same first-cousin construction is applied separately to each
    # replicate's travel geometry. Gold remains the shared null encoding.
    for label in REPLICATE_STYLES:
        cloud = nulls[nulls["replicate"] == label]
        ax.scatter(cloud["travel_sigma"], cloud["cell_state_sigma"],
                   s=7, color=NULL_COLOR, alpha=0.045, edgecolors="none",
                   rasterized=True, zorder=1)
        front = fronts[fronts["replicate"] == label]
        ax.scatter(front["null_travel_mean_sigma"].iloc[0],
                   front["null_cell_state_mean_sigma"].iloc[0],
                   marker="+", s=42,
                   color=NULL_COLOR, lw=1.2, zorder=5)

    for label, style in REPLICATE_STYLES.items():
        front = fronts[fronts["replicate"] == label].sort_values("sweep_index")
        cutoff = int(front["cutoff"].iloc[0])
        ax.plot(front["travel_sigma"], front["cell_state_sigma"],
                color=style["color"], ls=style["ls"], lw=1.8,
                label=f"{label} (T={cutoff})", zorder=3)
        key = front[front["is_max_retention"]]
        ax.scatter(key["travel_sigma"], key["cell_state_sigma"],
                   marker="D", s=42, facecolor=RETENTION_COLOR,
                   edgecolor=INK, lw=0.6, zorder=6)

    ax.scatter([0], [0], marker="X", s=82, facecolor=INK,
               edgecolor="white", lw=0.7, zorder=8)
    n_matched = int(audit["matched_terminal_edges"].iloc[0])
    ax.set_xlabel("Travel distance (first-cousin-null σ; natural lineage = 0)")
    ax.set_ylabel("Cell-state distance (first-cousin-null σ; natural lineage = 0)")
    ax.set_title(f"Matched terminal-cell Pareto fronts (n={n_matched} edges)",
                 loc="left", pad=5)
    ax.grid(True)
    _panel_letter(ax, "A", x=-0.11, y=1.07)

    cutoffs = audit.set_index("replicate")["cutoff"].astype(int)
    replicate_handles = [
        Line2D([0], [0], color=style["color"], ls=style["ls"], lw=1.8,
               label=f"{label} (T={cutoffs[label]})")
        for label, style in REPLICATE_STYLES.items()
    ]
    semantic_handles = [
        Line2D([0], [0], marker="X", ls="", markersize=7,
               markerfacecolor=INK, markeredgecolor="white",
               label="Natural lineage"),
        Line2D([0], [0], marker="D", ls="", markersize=5.5,
               markerfacecolor=RETENTION_COLOR, markeredgecolor=INK,
               label="Maximum edge retention"),
        Line2D([0], [0], marker="o", ls="", markersize=5,
               markerfacecolor=NULL_COLOR, markeredgecolor="none", alpha=0.55,
               label="First-cousin shuffle"),
        Line2D([0], [0], marker="+", ls="", markersize=7,
               color=NULL_COLOR, label="Null mean"),
    ]
    legend_reps = ax.legend(
        handles=replicate_handles, title="Tracking replicate",
        loc="upper left", bbox_to_anchor=(1.015, 1.0), borderaxespad=0,
        title_fontsize=7.3,
    )
    ax.add_artist(legend_reps)
    ax.legend(handles=semantic_handles, loc="upper left",
              bbox_to_anchor=(1.015, 0.70), borderaxespad=0)
    return fig


def plot_metrics(metrics):
    """Panel B: short horizontal marks in categorical subtree coordinates."""
    metric_specs = [
        ("u_L", r"Canonical position, $u_L$", (0, 1)),
        ("d_LP", r"Front proximity, $d_{LP}$", None),
        ("max_edge_retention", "Maximum edge retention", (0, 1)),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.95), sharex=True)
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.25, top=0.82,
                        wspace=0.33)
    x = np.arange(len(MAJOR_SUBTREES), dtype=float)
    half_width = 0.13

    for ax, (column, title, ylim) in zip(axes, metric_specs):
        for label, style in REPLICATE_STYLES.items():
            rows = metrics[metrics["replicate"] == label].set_index("subtree")
            values = rows.loc[list(MAJOR_SUBTREES), column].to_numpy()
            ax.hlines(values, x - half_width, x + half_width,
                      color=style["color"], lw=1.35, alpha=0.60,
                      zorder=3, label=label)
        ax.set_title(title, pad=5)
        ax.set_xticks(x, MAJOR_SUBTREES, rotation=35, ha="right")
        values = metrics[column]
        span = float(values.max() - values.min())
        minimum_pad = {"u_L": 0.035, "d_LP": 0.004,
                       "max_edge_retention": 0.035}[column]
        pad = max(0.18 * span, minimum_pad)
        lower = float(values.min() - pad)
        upper = float(values.max() + pad)
        if ylim is not None:
            lower = max(ylim[0], lower)
            upper = min(ylim[1] + (0.015 if ylim[1] == 1 else 0), upper)
        ax.set_ylim(lower, upper)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    _panel_letter(axes[0], "B", x=-0.31, y=1.22)

    handles = [
        Line2D([0], [0], marker="_", ls="", markersize=11,
               markeredgewidth=2.2, color=style["color"], alpha=0.60,
               label=label)
        for label, style in REPLICATE_STYLES.items()
    ]
    fig.legend(handles=handles, title="Tracking replicate", ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 0.995),
               title_fontsize=7.3, columnspacing=1.5, handletextpad=0.5)
    return fig


def save_panels(fronts, nulls, metrics, audit):
    PUBLICATION_OUT.mkdir(parents=True, exist_ok=True)
    fig_a = plot_fronts(fronts, nulls, audit)
    ps.save_figure(
        fig_a, PUBLICATION_OUT / "figS3A_ce_tracking_replicate_fronts.png"
    )
    plt.close(fig_a)
    fig_b = plot_metrics(metrics)
    ps.save_figure(
        fig_b, PUBLICATION_OUT / "figS3B_ce_tracking_replicate_metrics.png"
    )
    plt.close(fig_b)


def validate(context, audit, cells, subtrees, fronts=None, nulls=None,
             metrics=None):
    checks = []

    def check(name, condition, detail=""):
        checks.append(bool(condition))
        print(f"  [{'PASS' if condition else 'FAIL'}] {name} {detail}")

    expected_labels = [row[0] for row in dl.CE_REPLICATES]
    check("replicate order", list(context["replicates"]) == expected_labels)
    common_n = len(context["common_terms"])
    check("cell manifest common count",
          int(cells["present_all_replicates"].sum()) == common_n,
          f"n={common_n}")
    check("release matched-edge count",
          common_n == EXPECTED_MATCHED_EDGES,
          f"expected={EXPECTED_MATCHED_EDGES}, observed={common_n}")
    check("strict match retains >=90% per replicate",
          bool((audit["matched_fraction"] >= 0.90).all()),
          ", ".join(f"{row.replicate}={row.matched_fraction:.3f}"
                    for row in audit.itertuples()))
    check("adjacent frames change terminal coverage",
          bool(((audit["terminal_edges_previous_frame"]
                 != audit["terminal_edges_available"])
                & (audit["terminal_edges_next_frame"]
                   != audit["terminal_edges_available"])).all()),
          ", ".join(
              f"{row.replicate}={row.terminal_edges_previous_frame}/"
              f"{row.terminal_edges_available}/"
              f"{row.terminal_edges_next_frame}"
              for row in audit.itertuples()
          ))
    major = subtrees[subtrees["major_subtree"]].set_index("subtree")
    check("all major subtrees present", set(MAJOR_SUBTREES) <= set(major.index))
    check("major subtree counts are matched and usable",
          bool((major.loc[list(MAJOR_SUBTREES), "n_common"] >= MIN_CELLS).all()))
    eligible_n = int(subtrees["eligible_common_min12"].sum())
    check("release qualifying-subtree count",
          eligible_n == EXPECTED_COMMON_SUBTREES,
          f"expected={EXPECTED_COMMON_SUBTREES}, observed={eligible_n}")

    if fronts is not None and nulls is not None and metrics is not None:
        check("endpoint-inclusive 300-interval publication sweep",
              bool((fronts.groupby("replicate").size() == ITERATION + 1).all()))
        check("displayed null-cloud draw count",
              bool((nulls.groupby("replicate").size() == NULL_DRAWS).all()))
        check("one maximum-retention marker per front",
              bool((fronts.groupby("replicate")["is_max_retention"].sum()
                    == 1).all()))
        check("major metrics complete",
              len(metrics) == len(expected_labels) * len(MAJOR_SUBTREES))
        expected_pairs = {
            (label, subtree)
            for label in expected_labels for subtree in MAJOR_SUBTREES
        }
        observed_pairs = set(zip(metrics["replicate"], metrics["subtree"]))
        check("replicate/subtree metric grid complete",
              observed_pairs == expected_pairs)
        expected_n = major["n_common"].to_dict()
        check("metric sample sizes match subtree manifest",
              bool((metrics["n"]
                    == metrics["subtree"].map(expected_n)).all()))
        numeric = metrics[["u_L", "d_LP", "d_NP", "max_edge_retention"]]
        check("metrics finite", bool(np.isfinite(numeric.to_numpy()).all()))
        check("distance metrics non-negative",
              bool(((metrics["d_LP"] >= 0) & (metrics["d_NP"] >= 0)).all()))
        check("bounded u and retention",
              bool(((metrics["u_L"] >= 0) & (metrics["u_L"] <= 1)
                    & (metrics["max_edge_retention"] >= 0)
                    & (metrics["max_edge_retention"] <= 1)).all()))

    if not all(checks):
        raise AssertionError(f"{sum(not c for c in checks)} validation checks failed")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-only", action="store_true",
                        help="Write availability/matching manifests and stop")
    return parser.parse_args()


def main():
    args = parse_args()
    ps.configure()
    context, audit = audit_replicates()
    cells, subtrees = build_manifests(context, audit)
    print("Replicate audit:")
    for row in audit.itertuples():
        print(f"  {row.replicate}: T={row.cutoff}, "
              f"available={row.terminal_edges_available}, "
              f"matched={row.matched_terminal_edges} "
              f"({row.matched_fraction:.1%}); adjacent frames="
              f"{row.terminal_edges_previous_frame}/"
              f"{row.terminal_edges_next_frame}")
    print(f"  Common qualifying subtrees (n>={MIN_CELLS}): "
          f"{int(subtrees['eligible_common_min12'].sum())}")
    if args.audit_only:
        validate(context, audit, cells, subtrees)
        return

    fronts, nulls, metrics = run_analysis(context)
    validate(context, audit, cells, subtrees, fronts, nulls, metrics)
    save_panels(fronts, nulls, metrics, audit)
    print("Major-subtree replicate ranges:")
    for column in ["u_L", "d_LP", "d_NP", "max_edge_retention"]:
        spreads = metrics.groupby("subtree")[column].agg(lambda x: x.max() - x.min())
        print(f"  {column}: maximum spread={spreads.max():.3f} "
              f"({spreads.idxmax()})")
    print("Wrote Figure S3 analysis caches and publication panels")


if __name__ == "__main__":
    main()
