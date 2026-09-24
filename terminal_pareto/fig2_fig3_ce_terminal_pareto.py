"""Generate C. elegans terminal-cell Pareto Figures 2 and 3."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto import plot_style as ps
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
)
from terminal_pareto.front_coordinates import (
    build_endpoint_transform,
    null_sd_coordinates,
    percent_natural_coordinates,
)
from terminal_pareto.global_analysis import get_or_compute_global_analysis


OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "publication"
DIAGNOSTIC_OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "ce_protein"
EDGE_RETENTION_CMAP = LinearSegmentedColormap.from_list(
    "edge_retention_blue", ["#17365D", "#0072B2", "#72C7EC"]
)
TRAVEL_COLOR = ps.SEMANTIC_COLORS["travel"]
CELL_STATE_COLOR = ps.SEMANTIC_COLORS["cell_state"]
TREE_DISTANCE_COLOR = ps.SEMANTIC_COLORS["tree_distance"]
FIRST_COUSIN_COLOR = ps.SEMANTIC_COLORS["first_cousin_null"]
SECOND_COUSIN_COLOR = ps.NULL_MODEL_COLORS["second_cousin"]
THIRD_COUSIN_COLOR = ps.NULL_MODEL_COLORS["third_cousin"]
FULL_RANDOM_COLOR = ps.NULL_MODEL_COLORS["full_random"]
# Panels B and C are assembled side by side in Figure 3. Their top edges are
# already identical; fixing the axes baseline keeps their horizontal axes on
# the same level despite B having a two-line x label and C a one-line label.
SUPPORT_AXES_BOTTOM = 0.217


def _standardize(x, y, reference):
    return (
        (x - reference["lineage_xyz"]) / reference["xyz_std"],
        (y - reference["lineage_exp"]) / reference["exp_std"],
    )


def _load_analysis():
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_map, valid = dl.load_elegans_tracking(dl.T_CE)
    protein = dl.load_protein_expression()
    selected = dl.load_prot_sel()
    valid = [name for name in valid if name in protein.index]
    terminal_nodes, terminal_parents = dl.collect_terminals(lineage, valid)
    xyz_mat, exp_mat, _ = pe.build_cost_matrices(
        terminal_nodes, terminal_parents, xyz_map, protein, selected
    )
    first_groups = pe.build_ancestor_groups(terminal_nodes, tree_index, 2)
    first_stats = pe.compute_cousin_random_stats(
        xyz_mat, exp_mat, first_groups, n_random=1000, seed=42
    )
    twr = lm.combined_lineage_proximity(
        xyz_mat, exp_mat, terminal_parents, terminal_nodes, tree_index,
        first_stats, iteration=300, n_random=100,
    )

    nulls = {}
    for degree, steps, seed in [(1, 2, 42), (2, 3, 43), (3, 4, 44)]:
        groups = pe.build_ancestor_groups(terminal_nodes, tree_index, steps)
        raw_x, raw_y = pe.compute_group_shuffle_costs(
            xyz_mat, exp_mat, groups, n_random=1000, seed=seed
        )
        nulls[degree] = _standardize(raw_x, raw_y, first_stats)
    rows = np.arange(len(terminal_nodes))
    rng = np.random.default_rng(45)
    full_x = np.empty(2000)
    full_y = np.empty(2000)
    for i in range(2000):
        perm = rng.permutation(len(terminal_nodes))
        full_x[i] = xyz_mat[rows, perm].sum()
        full_y[i] = exp_mat[rows, perm].sum()
    nulls["full"] = _standardize(full_x, full_y, first_stats)
    return twr, nulls


def _profile_analysis(context, *, force=False):
    """Load one cache-keyed numerical result and prepare its null-SD view."""
    result = get_or_compute_global_analysis(context, force=force)
    natural = result.natural_costs
    stds = result.null_stds
    key_map = {
        1: "first_cousin",
        2: "second_cousin",
        3: "third_cousin",
        "full": "full_random",
    }
    nulls = {
        key: null_sd_coordinates(
            *result.null_raw[cache_key], natural_costs=natural,
            null_stds=stds)
        for key, cache_key in key_map.items()
    }
    return result.twr, nulls, result


def _display_data(twr, nulls, *, mode, result=None, context=None):
    """Return plot coordinates without modifying saved numerical results."""
    reference = context.cache_key if context is not None else "legacy"
    if mode == "null_sd":
        return dict(
            x=np.asarray(twr["xyz_arr"], dtype=float),
            y=np.asarray(twr["exp_arr"], dtype=float),
            nulls=nulls,
            natural=(0.0, 0.0),
            xlabel=("Travel distance\n(null standard deviations; "
                    "natural lineage = 0)"),
            ylabel=("Cell-state distance\n(null standard deviations; "
                    "natural lineage = 0)"),
            metadata={"display_mode": "null_sd",
                      "reference_analysis_id": reference},
        )
    if result is None or context is None:
        raise ValueError(f"Display mode {mode!r} requires a saved profile result")
    key_map = {
        1: "first_cousin", 2: "second_cousin",
        3: "third_cousin", "full": "full_random",
    }
    if mode == "endpoint":
        transform = build_endpoint_transform(
            result.raw_front_travel, result.raw_front_state,
            travel_optimum_index=context.spec.sweep_intervals,
            state_optimum_index=0,
            reference_analysis_id=result.analysis_cache_key,
            assignment_ids=[f"sweep:{index}" for index in range(
                context.spec.sweep_intervals + 1)],
        )
        x, y = transform.transform(
            result.raw_front_travel, result.raw_front_state)
        natural_x, natural_y = transform.transform(
            np.asarray([result.natural_costs[0]]),
            np.asarray([result.natural_costs[1]]))
        display_nulls = {
            key: transform.transform(*result.null_raw[raw_key])
            for key, raw_key in key_map.items()
        }
        return dict(
            x=x, y=y, nulls=display_nulls,
            natural=(float(natural_x[0]), float(natural_y[0])),
            xlabel="Travel distance\n(fraction of endpoint cost span)",
            ylabel="Cell-state distance\n(fraction of endpoint cost span)",
            metadata=transform.metadata(clipping=False),
        )
    if mode == "percent_natural":
        x, y = percent_natural_coordinates(
            result.raw_front_travel, result.raw_front_state,
            natural_costs=result.natural_costs)
        display_nulls = {
            key: percent_natural_coordinates(
                *result.null_raw[raw_key], natural_costs=result.natural_costs)
            for key, raw_key in key_map.items()
        }
        return dict(
            x=x, y=y, nulls=display_nulls, natural=(0.0, 0.0),
            xlabel="Travel-distance change from natural lineage (%)",
            ylabel="Cell-state-distance change from natural lineage (%)",
            metadata={"display_mode": "percent_natural",
                      "reference_analysis_id": result.analysis_cache_key,
                      "natural_costs": list(result.natural_costs)},
        )
    raise ValueError(f"Unsupported display mode: {mode}")


def _save(fig, stem, dpi=400, out_dir=OUT, fixed_canvas=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Side-by-side panels use the full, identical figure canvas so independent
    # tight cropping cannot undo their shared axes baseline.
    bbox = fig.bbox_inches if fixed_canvas else "tight"
    fig.savefig(out_dir / f"{stem}.pdf", facecolor="white",
                bbox_inches=bbox)
    fig.savefig(out_dir / f"{stem}.png", dpi=dpi, facecolor="white",
                bbox_inches=bbox)


def proportional_limits(values, pad_fraction=0.12):
    """Return finite data limits with padding proportional to their range."""
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("Cannot determine limits from non-finite data")
    lo = float(finite.min())
    hi = float(finite.max())
    span = hi - lo
    if span <= np.finfo(float).eps * max(abs(lo), abs(hi), 1.0):
        span = max(abs(lo), abs(hi), 1.0)
    pad = float(pad_fraction) * span
    return lo - pad, hi + pad


def plot_main(twr, nulls, *, display=None, out_dir=OUT,
              panel_letter=None):
    """Full-width main figure: Pareto front and progressively relaxed nulls."""
    x = np.asarray(display["x"] if display else twr["xyz_arr"])
    y = np.asarray(display["y"] if display else twr["exp_arr"])
    shown_nulls = display["nulls"] if display else nulls
    ct = twr["cost_tree_tradeoff"]
    kp = twr["kp"]
    if display:
        lineage_x, lineage_y = display["natural"]
    else:
        lineage_x = x[0] + ct["delta_xyz"][0]
        lineage_y = y[0] + ct["delta_exp"][0]

    fig, ax = plt.subplots(figsize=(7.15, 4.85))
    fig.subplots_adjust(left=0.12, right=0.88, bottom=0.15, top=0.95)
    ax.axhline(0, color="#777777", lw=0.7, ls=":", zorder=0)
    ax.axvline(0, color="#777777", lw=0.7, ls=":", zorder=0)

    null_styles = [
        (1, "First-cousin shuffle", FIRST_COUSIN_COLOR),
        (2, "Second-cousin shuffle", SECOND_COUSIN_COLOR),
        (3, "Third-cousin shuffle", THIRD_COUSIN_COLOR),
    ]
    for degree, label, color in null_styles:
        nx, ny = shown_nulls[degree]
        ax.scatter(nx, ny, s=10, color=color, alpha=0.19,
                   edgecolors="none", label=label, zorder=1)
        mean = (float(np.mean(nx)), float(np.mean(ny)))
        ax.scatter(*mean, marker="+", s=38, color=color, lw=1.2, zorder=4)

    edge_retention = np.asarray(twr["traditional_er"])
    norm = Normalize(vmin=float(np.nanmin(edge_retention)), vmax=1.0)
    edge_cmap = EDGE_RETENTION_CMAP
    ax.plot(x, y, color="#0072B2", lw=1.3, alpha=0.5, zorder=2)
    front = ax.scatter(x, y, c=edge_retention, cmap=edge_cmap, norm=norm, s=20,
                       edgecolors="none", zorder=3)
    ax.scatter([lineage_x], [lineage_y], marker="X", s=78, color="#222222",
               edgecolor="white", lw=0.6, label="Natural lineage", zorder=7)
    k = kp["max_er_idx"]
    max_er_color = edge_cmap(norm(edge_retention[k]))
    ax.scatter([x[k]], [y[k]], marker="o", s=32,
               facecolor=max_er_color, edgecolor="#222222", lw=1.0,
               label="Maximum edge retention", zorder=7)
    ax.plot([lineage_x, x[k]], [lineage_y, y[k]], color="#222222",
            lw=0.8, ls=(0, (2, 2)), zorder=5)

    full_x, full_y = shown_nulls["full"]
    inset = ax.inset_axes([0.66, 0.71, 0.31, 0.24])
    inset.scatter(full_x, full_y, s=6, color=FULL_RANDOM_COLOR, alpha=0.18,
                  edgecolors="none", rasterized=True)
    full_mean = (float(np.mean(full_x)), float(np.mean(full_y)))
    inset.scatter(*full_mean, marker="+", s=30, color=FULL_RANDOM_COLOR,
                  lw=1.1)
    inset.set_xlim(*proportional_limits(full_x))
    inset.set_ylim(*proportional_limits(full_y))
    inset.set_title("Full-random shuffle", fontsize=7, pad=2)
    inset.set_xlabel("Travel distance", fontsize=6, labelpad=1)
    inset.set_ylabel("Cell-state distance", fontsize=6, labelpad=1)
    inset.tick_params(labelsize=5.5, length=2, pad=1)
    inset.grid(True, alpha=0.25)

    full_handle = Line2D([0], [0], marker="o", ls="", markersize=5,
                         markerfacecolor=FULL_RANDOM_COLOR,
                         markeredgecolor="none",
                         alpha=0.7, label="Full-random shuffle (inset)")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend_labels = [
        "First-cousin shuffle", "Second-cousin shuffle",
        "Third-cousin shuffle", "Full-random shuffle (inset)",
        "Natural lineage", "Maximum edge retention",
    ]
    by_label[full_handle.get_label()] = full_handle
    ax.legend([by_label[label] for label in legend_labels], legend_labels,
              loc="upper left", bbox_to_anchor=(0.66, 0.59), fontsize=6.2,
              frameon=True, facecolor="white", edgecolor="none",
              framealpha=0.92, borderaxespad=0, borderpad=0.25,
              labelspacing=0.22, handletextpad=0.4)
    ax.set_xlabel(display["xlabel"] if display else
                  "Travel distance\n(σ; natural lineage = 0)")
    ax.set_ylabel(display["ylabel"] if display else
                  "Cell-state distance\n(σ; natural lineage = 0)")
    if display:
        visible_x = np.concatenate(
            [x, np.asarray([lineage_x])] +
            [np.asarray(shown_nulls[key][0]) for key in (1, 2, 3)])
        visible_y = np.concatenate(
            [y, np.asarray([lineage_y])] +
            [np.asarray(shown_nulls[key][1]) for key in (1, 2, 3)])
        xpad = max(0.04 * np.ptp(visible_x), 0.04)
        ypad = max(0.04 * np.ptp(visible_y), 0.04)
        ax.set_xlim(visible_x.min() - xpad, visible_x.max() + xpad)
        ax.set_ylim(visible_y.min() - ypad, visible_y.max() + ypad)
    else:
        ax.set_xlim(min(-8, x.min() - 2), max(90, x.max() + 3))
        ax.set_ylim(min(-45, y.min() - 3), max(50, nulls[3][1].max() + 3))
    ax.grid(True, alpha=0.32)
    if panel_letter:
        ax.text(-0.16, 1.04, panel_letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", ha="left", va="top",
                clip_on=False)
    cbar = fig.colorbar(front, ax=ax, fraction=0.047, pad=0.025)
    cbar.set_label("Edge retention")
    _save(fig, "fig2_ce_terminal_pareto_front", out_dir=out_dir)
    return fig


def plot_support_b(twr, distance_mode="changed_edges", *, display=None,
                   out_dir=OUT, diagnostic_out=DIAGNOSTIC_OUT):
    x = np.asarray(display["x"] if display else twr["xyz_arr"])
    er = np.asarray(twr["traditional_er"])
    null = twr["lineage_null"]
    if distance_mode == "changed_edges":
        per_cell = np.asarray(twr["per_cell_dist"])
        td = np.array([
            row[row > 0].mean() if np.any(row > 0) else 0
            for row in per_cell
        ])
        full_mean = null["full_changed_mean"]
        cousin_mean = null["cousin_changed_mean"]
        distance_label = "Mean tree distance\nof changed edges"
    elif distance_mode == "all_edges":
        td = np.asarray(twr["lineage_mean_dist"])
        full_mean = null["full_mean"]
        cousin_mean = null["cousin_mean"]
        distance_label = "Mean tree distance\nacross all edges"
    else:
        raise ValueError(f"Unknown tree-distance mode: {distance_mode}")
    ct = twr["cost_tree_tradeoff"]
    k = twr["kp"]["max_er_idx"]
    lineage_x = (display["natural"][0] if display
                 else x[0] + ct["delta_xyz"][0])
    fig, ax = plt.subplots(figsize=(3.45, 2.85))
    ax_td = ax.twinx()
    er_line, = ax.plot(x, er, color="#0072B2", lw=2, label="Edge retention")
    tree_color = TREE_DISTANCE_COLOR
    td_line, = ax_td.plot(x, td, color=tree_color, lw=1.6,
                          label="Tree distance")
    td_max = max(td.max(), full_mean) * 1.06
    # Leave matching headroom on both axes: ER=1 and tree distance=0 occupy
    # exactly the same vertical level, with the lineage cross fully visible.
    ax_td.set_ylim(td_max, -0.04 * td_max)
    for value, label, color, ls in [
        (full_mean, "Full random", FULL_RANDOM_COLOR, "--"),
        (cousin_mean, "First cousin", FIRST_COUSIN_COLOR, ":"),
    ]:
        ax_td.plot([0.965, 1], [value, value],
                   transform=ax_td.get_yaxis_transform(), color=color,
                   lw=1.1, ls=ls, clip_on=False)
        ax_td.text(0.955, value, f"{label} mean",
                   transform=ax_td.get_yaxis_transform(), ha="right",
                   va="center", fontsize=6, color=color)
    natural = ax.scatter([lineage_x], [1], marker="X", s=48,
                         color="#222222", edgecolor="white", lw=0.5,
                         zorder=7, label="Natural lineage")
    er_norm = Normalize(vmin=float(np.nanmin(er)), vmax=1.0)
    max_er_color = EDGE_RETENTION_CMAP(er_norm(er[k]))
    ax.scatter([x[k]], [er[k]], marker="o", s=30,
               facecolor=max_er_color, edgecolor="#222222", lw=0.9, zorder=7,
               label="Maximum edge retention")
    min_td_idx = int(np.nanargmin(td))
    ax_td.scatter([x[min_td_idx]], [td[min_td_idx]], marker="o", s=35,
                  facecolor=tree_color, edgecolor="white", lw=0.6,
                  zorder=6, label="Minimum tree distance")
    ax.set_xlabel(display["xlabel"] if display else
                  "Travel distance\n(σ; natural lineage = 0)")
    ax.set_ylabel("Edge retention", color="#0072B2")
    ax_td.set_ylabel(distance_label, color=tree_color)
    ax.tick_params(axis="y", colors="#0072B2")
    ax_td.tick_params(axis="y", colors=tree_color)
    ax_td.spines["right"].set_color(tree_color)
    ax.set_ylim(0, 1.04)
    marker_handles = [
        natural,
        Line2D([0], [0], marker="o", ls="", markersize=5,
               markerfacecolor=max_er_color, markeredgecolor="#222222",
               label="Maximum edge retention"),
        Line2D([0], [0], marker="o", ls="", markersize=5,
               markerfacecolor=tree_color, markeredgecolor="white",
               label="Minimum tree distance"),
    ]
    ax.legend(handles=marker_handles, loc="lower left", fontsize=5.7,
              frameon=False, labelspacing=0.25, handletextpad=0.4)
    ax.text(-0.16, 1.08, "B", transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")
    ax.grid(True, alpha=0.32)
    fig.tight_layout()
    for support_ax in (ax, ax_td):
        pos = support_ax.get_position()
        support_ax.set_position([
            pos.x0, SUPPORT_AXES_BOTTOM, pos.width,
            pos.y1 - SUPPORT_AXES_BOTTOM,
        ])
    if distance_mode == "changed_edges":
        _save(fig, "fig3B_ce_edge_retention_tree_distance",
              fixed_canvas=True, out_dir=out_dir)
    else:
        # Retain the selectable all-edge definition as a diagnostic without
        # adding non-manuscript variants to output/publication.
        _save(fig, "fig3B_ce_edge_retention_tree_distance_all_edges",
              out_dir=diagnostic_out)
    return fig


def plot_support_c(twr, *, display=None, out_dir=OUT):
    x = np.asarray(display["x"] if display else twr["xyz_arr"])
    y = np.asarray(display["y"] if display else twr["exp_arr"])
    tree_distance = np.asarray(twr["lineage_mean_dist"])
    edge_retention = np.asarray(twr["traditional_er"])
    k = twr["kp"]["max_er_idx"]
    fig, ax = plt.subplots(figsize=(3.45, 2.85))

    trajectories = [
        ("Toward travel optimum", np.arange(k, len(x)), x, TRAVEL_COLOR),
        ("Toward cell-state optimum", np.arange(k, -1, -1), y,
         CELL_STATE_COLOR),
    ]
    final_prices = []
    for label, indices, cost, color in trajectories:
        distance_reduction = cost[k] - cost[indices]
        reduction_fraction = distance_reduction / distance_reduction[-1]
        retention = 100 * edge_retention[indices]
        ax.plot(100 * reduction_fraction, retention, color=color, lw=2.0,
                label=label, zorder=3)
        threshold_idx = np.flatnonzero(reduction_fraction >= 0.95)[0]
        final_prices.append((
            retention[threshold_idx] - retention[-1],
            tree_distance[indices[-1]] - tree_distance[indices[threshold_idx]],
            retention[threshold_idx], retention[-1], color,
        ))
        ax.scatter([100 * reduction_fraction[threshold_idx], 100],
                   [retention[threshold_idx], retention[-1]], s=[18, 24],
                   color=color, edgecolor="white", lw=0.4, zorder=5)

    ax.axvspan(95, 100, color="#BDBDBD", alpha=0.20, lw=0,
               label="Final 5% of attainable reduction")
    ax.scatter([0], [100 * edge_retention[k]], marker="D", s=30,
               facecolor="white",
               edgecolor="#222222", lw=0.8, zorder=5)
    for ypos, (edge_loss, tree_gain, _, _, color) in zip((63, 28), final_prices):
        ax.text(102.0, ypos,
                f"Final 5%:\n−{edge_loss:.1f} pp edges\n+{tree_gain:.1f} tree distance",
                color=color, fontsize=5.3, ha="left", va="center")

    # Minimal map back to the Pareto front: the diamond is the common starting
    # point, and colored arrows show the two directions represented above.
    inset = ax.inset_axes([0.10, 0.10, 0.38, 0.29])
    inset.plot(x, y, color="#777777", lw=1.0)
    inset.scatter([x[k]], [y[k]], marker="D", s=20, facecolor="white",
                  edgecolor="#222222", lw=0.7, zorder=4)
    travel_target = len(x) - 1
    state_target = 0
    inset.annotate("", xy=(x[travel_target], y[travel_target]),
                   xytext=(x[k], y[k]),
                   arrowprops=dict(arrowstyle="->", color=TRAVEL_COLOR, lw=1.5,
                                   connectionstyle="arc3,rad=-0.12"))
    inset.annotate("", xy=(x[state_target], y[state_target]),
                   xytext=(x[k], y[k]),
                   arrowprops=dict(arrowstyle="->", color=CELL_STATE_COLOR, lw=1.5,
                                   connectionstyle="arc3,rad=0.12"))
    inset.text(0.05, 0.93, "Travel", transform=inset.transAxes,
               color=TRAVEL_COLOR, fontsize=4.8, va="top")
    inset.text(0.98, 0.14, "Cell state", transform=inset.transAxes,
               color=CELL_STATE_COLOR, fontsize=4.8, ha="right")
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_xlabel("Travel distance", fontsize=4.8, labelpad=1)
    inset.set_ylabel("Cell-state distance", fontsize=4.8, labelpad=1)
    inset.grid(False)

    ax.set_title("Moving along the Pareto front", fontsize=8, pad=5)
    ax.set_xlabel("Attainable distance reduction achieved (%)")
    ax.set_ylabel("Natural edges retained (%)")
    ax.set_xlim(-2, 124)
    ax.set_ylim(15, 86)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.legend(loc="upper right", fontsize=5.2, frameon=False,
              labelspacing=0.25, handlelength=1.7)
    ax.text(-0.16, 1.08, "C", transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")
    ax.grid(True, alpha=0.32)
    fig.tight_layout()
    pos = ax.get_position()
    ax.set_position([
        pos.x0, SUPPORT_AXES_BOTTOM, pos.width,
        pos.y1 - SUPPORT_AXES_BOTTOM,
    ])
    _save(fig, "fig3C_ce_structural_retention", fixed_canvas=True,
          out_dir=out_dir)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the terminal-cell publication figure panels."
    )
    parser.add_argument(
        "--tree-distance-mode",
        choices=("changed_edges", "all_edges", "both"),
        default="changed_edges",
        help=("Tree-distance definition(s) for supporting panel B "
              "(default: changed_edges; all_edges is diagnostic-only)."),
    )
    parser.add_argument(
        "--profile",
        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"),
        help="Opt into an isolated profile-aware run (omission keeps legacy paths).",
    )
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--display", choices=("null_sd", "endpoint",
                                               "percent_natural"),
                        default="null_sd")
    parser.add_argument(
        "--main-panel-letter",
        help="Optional panel letter for the Figure 2 front component.",
    )
    parser.add_argument("--force-analysis", action="store_true")
    args = parser.parse_args(argv)
    ps.configure()
    display = None
    out_dir = OUT
    diagnostic_out = DIAGNOSTIC_OUT
    if args.profile is None:
        twr, nulls = _load_analysis()
    else:
        context = build_analysis_context(
            args.profile, run_id=args.run_id, output_root=args.output_root)
        context.write()
        twr, nulls, result = _profile_analysis(
            context, force=args.force_analysis)
        display = _display_data(
            twr, nulls, mode=args.display, result=result, context=context)
        out_dir = context.run_paths.display(args.display)
        diagnostic_out = context.run_paths.analysis / "diagnostics" / args.display
        out_dir.mkdir(parents=True, exist_ok=True)
        diagnostic_out.mkdir(parents=True, exist_ok=True)
        metadata = dict(display["metadata"])
        metadata.update({
            "profile": context.profile,
            "context_cache_key": context.cache_key,
            "analysis_cache_key": result.analysis_cache_key,
            "assignment_ids_hash": hashlib.sha256(
                result.assignments.tobytes()).hexdigest(),
        })
        (out_dir / "fig2_fig3_display_manifest.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    plot_main(twr, nulls, display=display, out_dir=out_dir,
              panel_letter=args.main_panel_letter)
    modes = ("changed_edges", "all_edges") if args.tree_distance_mode == "both" else (args.tree_distance_mode,)
    for mode in modes:
        plot_support_b(twr, distance_mode=mode, display=display,
                       out_dir=out_dir, diagnostic_out=diagnostic_out)
    plot_support_c(twr, display=display, out_dir=out_dir)
    print("Publication figure panels written to", out_dir)


if __name__ == "__main__":
    main()
