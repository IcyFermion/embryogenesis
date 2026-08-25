"""Figure 6 panels B--C: optimization within terminal cell types.

This analysis asks a different question from ``fig6a_figs2_ce_cell_types.py``. Instead
of decomposing assignments from the unrestricted embryo-wide Pareto front,
it forbids assignments across terminal cell-type groups and solves the
assignment problem independently inside each group.

Two outputs are produced:

1. Individual within-type fronts in shared, embryo-level null-SD units per
   cell. This avoids pretending that degenerate or very sparse within-type
   cousin nulls define comparable scales.
2. The aggregate type-preserving front versus the unrestricted global front.
   Because the restricted feasible set is a subset of the global one, the
   restricted single-objective endpoints must never outperform the global
   endpoints; this is asserted as an implementation check.

The established first-cousin-null proximity is also computed per type when
both objective variances are positive and stored in the summary table. A
full-random-within-type sensitivity value is stored separately, never used as
a silent replacement for a degenerate cousin null.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
from scipy.optimize import linear_sum_assignment

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto import plot_style as ps
from terminal_pareto.subtree_analysis import (
    build_type_map,
    exact_cousin_stats,
    first_cousin_null_summary,
)


OUT = Path(__file__).resolve().parent / "output" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"
ITERATION = 300
MIN_DISPLAY_N = 12
TYPE_ORDER = [
    "programmed_death", "epithelium", "neuron", "muscle",
    "alimentary", "glial",
]
TYPE_COLORS = {
    "programmed_death": ps.COLORS["purple"],
    "epithelium": ps.COLORS["orange"],
    "neuron": ps.COLORS["blue"],
    "muscle": ps.COLORS["vermillion"],
    "alimentary": ps.COLORS["green"],
    "glial": ps.COLORS["sky"],
}
TRAVEL_COLOR = ps.SEMANTIC_COLORS["travel"]
CELL_STATE_COLOR = ps.SEMANTIC_COLORS["cell_state"]
UNRESTRICTED_COLOR = "#444444"
TYPE_RESTRICTED_COLOR = ps.COLORS["orange"]


def load_primary_data():
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    protein_exp = dl.load_protein_expression()
    prot_sel = dl.load_prot_sel()
    v_prot = [name for name in valid_ce if name in protein_exp.index]
    tn, tp = dl.collect_terminals(lineage, v_prot)
    return dict(
        lineage=lineage, tree_index=tree_index, xyz=xyz_ce,
        expression=protein_exp, features=prot_sel, tn=tn, tp=tp,
        type_map=build_type_map(tn), gp_map=pe.build_grandparent_map(lineage),
    )


def exact_null_stats(xm, em, groups):
    """Exact permutation-null moments in the random-stats API shape."""
    mx, me, vx, ve, _ = exact_cousin_stats(xm, em, groups)
    return dict(
        xyz_mean=float(mx), exp_mean=float(me),
        xyz_std=float(np.sqrt(max(vx, 0.0))),
        exp_std=float(np.sqrt(max(ve, 0.0))),
        lineage_xyz=float(np.diag(xm).sum()),
        lineage_exp=float(np.diag(em).sum()),
    )


def relative_front_record(xm, em, tp, stats, iteration):
    """Natural-to-front ratio using one explicitly supplied null model."""
    if not (stats["xyz_std"] > 1e-10 and stats["exp_std"] > 1e-10):
        return dict(relative_distance=np.nan, closest_sigma=np.nan,
                    null_to_front_sigma=np.nan, natural_on_front=False)
    xa, ea, _edge, _kp = pe.compute_std_scaled_pareto(
        xm, em, tp, stats, iteration=iteration)
    distance = np.hypot(xa, ea)
    closest = int(np.argmin(distance))
    nx = (stats["xyz_mean"] - stats["lineage_xyz"]) / stats["xyz_std"]
    ny = (stats["exp_mean"] - stats["lineage_exp"]) / stats["exp_std"]
    lp = float(distance[closest])
    np_dist = float(np.hypot(nx - xa[closest], ny - ea[closest]))
    return dict(
        relative_distance=lp / np_dist if np_dist > 1e-10 else np.nan,
        closest_sigma=lp,
        null_to_front_sigma=np_dist,
        natural_on_front=bool(np.any((np.abs(xa) < 1e-9)
                                     & (np.abs(ea) < 1e-9))),
    )


def solve_block_sweep(xs, es, groups, lineage_x, lineage_e, iteration):
    """Pareto sweep with assignments constrained to cell-type blocks."""
    travel = np.empty(iteration + 1)
    state = np.empty(iteration + 1)
    for step in range(iteration + 1):
        alpha = step / iteration
        total_x = total_e = 0.0
        for indices in groups.values():
            ix = np.asarray(indices, dtype=int)
            xb = xs[np.ix_(ix, ix)]
            eb = es[np.ix_(ix, ix)]
            ri, ci = linear_sum_assignment(alpha * xb + (1 - alpha) * eb)
            total_x += float(xb[ri, ci].sum())
            total_e += float(eb[ri, ci].sum())
        travel[step] = total_x - lineage_x
        state[step] = total_e - lineage_e
    return travel, state


def analyze(iteration=ITERATION):
    data = load_primary_data()
    tn, tp = data["tn"], data["tp"]
    xm, em, _ = pe.build_cost_matrices(
        tn, tp, data["xyz"], data["expression"], data["features"])
    global_null = first_cousin_null_summary(
        xm, em, tn, data["gp_map"], seed=42)["_raw"]
    xstd, estd = global_null["xyz_std"], global_null["exp_std"]
    xs, es = xm / xstd, em / estd
    lineage_x = float(np.diag(xs).sum())
    lineage_e = float(np.diag(es).sum())

    global_x, global_e, _edge, _kp = pe.compute_std_scaled_pareto(
        xm, em, tp, global_null, iteration=iteration)
    types = [data["type_map"][cell] for cell in tn]
    indices_by_type = {
        cell_type: [i for i, value in enumerate(types) if value == cell_type]
        for cell_type in sorted(set(types))
    }
    restricted_x, restricted_e = solve_block_sweep(
        xs, es, indices_by_type, lineage_x, lineage_e, iteration)

    # A type-preserving feasible set is a subset of the unrestricted set.
    tol = 1e-8
    if restricted_x[-1] < global_x[-1] - tol:
        raise AssertionError("Type-restricted travel optimum beats global optimum")
    if restricted_e[0] < global_e[0] - tol:
        raise AssertionError("Type-restricted state optimum beats global optimum")

    front_rows = []
    summary_rows = []
    for cell_type, indices in indices_by_type.items():
        ix = np.asarray(indices, dtype=int)
        n = len(ix)
        xb = xm[np.ix_(ix, ix)]
        eb = em[np.ix_(ix, ix)]
        xsb = xs[np.ix_(ix, ix)]
        esb = es[np.ix_(ix, ix)]
        lx = float(np.diag(xsb).sum())
        le = float(np.diag(esb).sum())
        local_x, local_e = solve_block_sweep(
            xsb, esb, {cell_type: list(range(n))}, lx, le, iteration)
        local_x /= n
        local_e /= n
        for step, (dx, de) in enumerate(zip(local_x, local_e)):
            front_rows.append(dict(
                type=cell_type, n=n, step=step, alpha=step / iteration,
                travel_sigma_per_cell=dx, state_sigma_per_cell=de,
            ))

        tnodes = [tn[i] for i in ix]
        tparents = [tp[i] for i in ix]
        cousin_groups = pe.build_cousin_groups(tnodes, data["gp_map"])
        cousin_stats = exact_null_stats(xb, eb, cousin_groups)
        cousin_rec = relative_front_record(
            xb, eb, tparents, cousin_stats, iteration)
        random_stats = exact_null_stats(
            xb, eb, [list(range(n))] if n >= 2 else [])
        random_rec = relative_front_record(
            xb, eb, tparents, random_stats, iteration)
        closest = float(np.min(np.hypot(local_x, local_e)))
        summary_rows.append(dict(
            type=cell_type,
            n=n,
            cousin_groups=len(cousin_groups),
            cousin_covered=sum(map(len, cousin_groups)),
            cousin_covered_frac=sum(map(len, cousin_groups)) / n,
            cousin_xyz_std=cousin_stats["xyz_std"],
            cousin_state_std=cousin_stats["exp_std"],
            cousin_null_valid=bool(cousin_stats["xyz_std"] > 1e-10
                                   and cousin_stats["exp_std"] > 1e-10),
            cousin_relative_distance=cousin_rec["relative_distance"],
            cousin_natural_on_front=cousin_rec["natural_on_front"],
            full_random_relative_distance=random_rec["relative_distance"],
            shared_scale_closest_per_cell=closest,
            travel_reduction_per_cell=-float(local_x[-1]),
            state_reduction_per_cell=-float(local_e[0]),
        ))

    aggregate = pd.DataFrame({
        "step": np.arange(iteration + 1),
        "alpha": np.arange(iteration + 1) / iteration,
        "global_travel_sigma": global_x,
        "global_state_sigma": global_e,
        "restricted_travel_sigma": restricted_x,
        "restricted_state_sigma": restricted_e,
    })
    endpoints = dict(
        travel_penalty_sigma=float(restricted_x[-1] - global_x[-1]),
        state_penalty_sigma=float(restricted_e[0] - global_e[0]),
        global_travel_reduction_sigma=float(-global_x[-1]),
        restricted_travel_reduction_sigma=float(-restricted_x[-1]),
        global_state_reduction_sigma=float(-global_e[0]),
        restricted_state_reduction_sigma=float(-restricted_e[0]),
    )
    return pd.DataFrame(front_rows), pd.DataFrame(summary_rows), aggregate, endpoints


def unique_front(data):
    """Remove consecutive duplicate display coordinates."""
    xy = data[["travel_sigma_per_cell", "state_sigma_per_cell"]].to_numpy()
    keep = np.r_[True, np.any(np.abs(np.diff(xy, axis=0)) > 1e-12, axis=1)]
    return data.loc[keep]


def save_panel_crops(fig, upper_axes, lower_ax, b_artists, c_artists,
                     out_dir):
    """Export Figure 6B and 6C independently from their shared canvas."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_inches = fig.dpi_scale_trans.inverted()
    width, height = fig.get_size_inches()

    upper_box = Bbox.union(
        [axis.get_tightbbox(renderer) for axis in upper_axes]
        + [artist.get_tightbbox(renderer) for artist in b_artists]
    ).transformed(to_inches)
    lower_box = Bbox.union(
        [lower_ax.get_tightbbox(renderer)]
        + [artist.get_tightbbox(renderer) for artist in c_artists]
    ).transformed(to_inches)
    crops = [
        ("fig6B_ce_within_type_fronts", upper_axes, b_artists,
         Bbox.from_extents(0.0, max(0.0, upper_box.y0 - 0.04),
                           width, min(height, upper_box.y1 + 0.04))),
        ("fig6C_ce_type_restricted_aggregate", [lower_ax], c_artists,
         Bbox.from_extents(0.0, max(0.0, lower_box.y0 - 0.04),
                           width, min(height, lower_box.y1 + 0.04))),
    ]

    original_axis_visibility = {axis: axis.get_visible() for axis in fig.axes}
    figure_artists = list(b_artists) + list(c_artists)
    original_artist_visibility = {
        artist: artist.get_visible() for artist in figure_artists
    }
    for stem, visible_axes, visible_artists, crop in crops:
        visible_axes = set(visible_axes)
        visible_artists = set(visible_artists)
        for axis in fig.axes:
            axis.set_visible(axis in visible_axes)
        for artist in figure_artists:
            artist.set_visible(artist in visible_artists)
        fig.savefig(out_dir / f"{stem}.pdf", facecolor="white",
                    bbox_inches=crop, pad_inches=0)
        fig.savefig(out_dir / f"{stem}.png", dpi=300, facecolor="white",
                    bbox_inches=crop, pad_inches=0)

    for axis, was_visible in original_axis_visibility.items():
        axis.set_visible(was_visible)
    for artist, was_visible in original_artist_visibility.items():
        artist.set_visible(was_visible)


def plot_within_type(front_df, summary_df, aggregate, endpoints, out_dir=OUT):
    """Render the within-type panels used below the retention heatmap."""
    fig = plt.figure(figsize=(7.15, 5.75))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.56)
    upper = outer[0].subgridspec(2, 3, hspace=0.58, wspace=0.34)
    axes = [fig.add_subplot(upper[i, j]) for i in range(2) for j in range(3)]

    for ax, cell_type in zip(axes, TYPE_ORDER):
        data = unique_front(front_df[front_df["type"] == cell_type])
        row = summary_df.set_index("type").loc[cell_type]
        color = TYPE_COLORS[cell_type]
        ax.plot(data["travel_sigma_per_cell"], data["state_sigma_per_cell"],
                color=color, lw=1.45, zorder=2)
        ax.scatter(data["travel_sigma_per_cell"], data["state_sigma_per_cell"],
                   s=5, color=color, alpha=0.45, edgecolors="none", zorder=3)
        ax.scatter([0], [0], marker="X", s=34, color="#222222",
                   edgecolor="white", lw=0.5, zorder=5)
        ax.scatter([data.iloc[-1]["travel_sigma_per_cell"]],
                   [data.iloc[-1]["state_sigma_per_cell"]],
                   marker="o", s=22, facecolor=TRAVEL_COLOR,
                   edgecolor="#222222", lw=0.45, zorder=5)
        ax.scatter([data.iloc[0]["travel_sigma_per_cell"]],
                   [data.iloc[0]["state_sigma_per_cell"]],
                   marker="s", s=22, facecolor=CELL_STATE_COLOR,
                   edgecolor="#222222", lw=0.45, zorder=5)
        ax.axhline(0, color="#999999", lw=0.55, ls=":")
        ax.axvline(0, color="#999999", lw=0.55, ls=":")
        title = cell_type.replace("_", " ")
        ax.set_title(f"{title} (n={int(row['n'])})", fontsize=7.4,
                     color="#222222")
        ax.grid(True, alpha=0.22)
        ax.tick_params(labelsize=6.1)

    for ax in axes[3:]:
        ax.set_xlabel("Travel change\n(global σ per cell)", fontsize=6.7)
    for ax in (axes[0], axes[3]):
        ax.set_ylabel("Cell-state change\n(global σ per cell)", fontsize=6.7)
    b_letter = fig.text(0.012, 0.975, "B", fontsize=10, fontweight="bold",
                        ha="left", va="top")
    b_heading = fig.text(
        0.055, 0.972, "Pareto fronts within terminal cell types",
        fontsize=9.2, fontweight="semibold", ha="left", va="top")

    ax = fig.add_subplot(outer[1])
    ax.plot(aggregate["global_travel_sigma"], aggregate["global_state_sigma"],
            color=UNRESTRICTED_COLOR, lw=1.8,
            label="Unrestricted assignment",
            zorder=3)
    ax.plot(aggregate["restricted_travel_sigma"],
            aggregate["restricted_state_sigma"],
            color=TYPE_RESTRICTED_COLOR, lw=1.8,
            label="Assignments restricted within cell type", zorder=3)
    ax.scatter([0], [0], marker="X", s=48, color="#222222",
               edgecolor="white", lw=0.55, zorder=5, label="Natural lineage")
    for xcol, ycol, color in [
        ("global_travel_sigma", "global_state_sigma", UNRESTRICTED_COLOR),
        ("restricted_travel_sigma", "restricted_state_sigma",
         TYPE_RESTRICTED_COLOR),
    ]:
        ax.scatter([aggregate.iloc[-1][xcol]], [aggregate.iloc[-1][ycol]],
                   marker="o", s=27, facecolor=color, edgecolor="#222222",
                   lw=0.5, zorder=5)
        ax.scatter([aggregate.iloc[0][xcol]], [aggregate.iloc[0][ycol]],
                   marker="s", s=27, facecolor=color, edgecolor="#222222",
                   lw=0.5, zorder=5)
    ax.axhline(0, color="#888888", lw=0.6, ls=":")
    ax.axvline(0, color="#888888", lw=0.6, ls=":")
    ax.set_xlabel("Travel-distance change (global null σ)")
    ax.set_ylabel("Cell-state-distance change\n(global null σ)")
    ax.set_title("Aggregate consequence of forbidding assignments across cell types",
                 loc="left", pad=7)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=6.4)
    note = (f"Endpoint penalty of type restriction:  "
            f"travel +{endpoints['travel_penalty_sigma']:.2f}σ   ·   "
            f"cell state +{endpoints['state_penalty_sigma']:.2f}σ")
    ax.text(0.015, 0.035, note, transform=ax.transAxes, fontsize=6.5,
            color="#555555", ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="#CCCCCC", lw=0.5, alpha=0.9))
    c_letter = fig.text(0.012, 0.385, "C", fontsize=10, fontweight="bold",
                        ha="left", va="top")

    endpoint_handles = [
        Line2D([0], [0], marker="o", ls="", ms=4.7,
               markerfacecolor=TRAVEL_COLOR, markeredgecolor="#222222",
               label="Travel optimum"),
        Line2D([0], [0], marker="s", ls="", ms=4.7,
               markerfacecolor=CELL_STATE_COLOR, markeredgecolor="#222222",
               label="Cell-state optimum"),
    ]
    endpoint_legend = fig.legend(
        handles=endpoint_handles, loc="upper right", ncol=2,
        bbox_to_anchor=(0.975, 0.975), fontsize=6.2, frameon=False)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.09)
    save_panel_crops(
        fig, axes, ax,
        b_artists=[b_letter, b_heading, endpoint_legend],
        c_artists=[c_letter], out_dir=out_dir,
    )
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iteration", type=int, default=ITERATION)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args(argv)
    ps.configure()
    fronts, summary, aggregate, endpoints = analyze(iteration=args.iteration)
    fronts.to_csv(ANALYSIS_OUT / "within_type_fronts.csv", index=False)
    summary.to_csv(ANALYSIS_OUT / "within_type_summary.csv", index=False)
    aggregate.to_csv(ANALYSIS_OUT / "type_preserving_aggregate_front.csv",
                     index=False)
    plot_within_type(fronts, summary, aggregate, endpoints, out_dir=args.out)
    print(summary[["type", "n", "cousin_null_valid",
                   "cousin_relative_distance", "full_random_relative_distance",
                   "travel_reduction_per_cell", "state_reduction_per_cell"]]
          .to_string(index=False))
    print("Endpoint comparison:", endpoints)
    print("Wrote within-type analysis and Figure 6 panels B--C to", args.out)


if __name__ == "__main__":
    main()
