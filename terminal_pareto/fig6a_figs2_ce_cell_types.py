"""Figure 6A and Figure S2: fate-specific terminal-front changes.

Panel A is a display-smoothed heatmap of natural-edge retention for each
terminal cell type across the observed Pareto assignments. Horizontal
position is normalized arc length along the same travel-to-cell-state
direction used by Figure 5. For presentation only, the irregularly spaced
assignments are resampled to a uniform grid and locally averaged; all
analysis and marked keypoints retain their exact unsmoothed values.

The endpoint cost-gain ledger is retained as a supplementary figure. It shows
the per-cell change relative to the natural lineage at the travel optimum,
maximum-retention compromise, and cell-state optimum. Travel and cell-state
axes remain separate and are never combined or compared linearly.

The established optimization and terminal-child attribution are unchanged.
The full-front retention cache simply evaluates the existing Pareto sweep at
all distinct assignments rather than keeping only three keypoints.
"""

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
from terminal_pareto.fig2_fig3_ce_terminal_pareto import EDGE_RETENTION_CMAP
from terminal_pareto.subtree_analysis import build_type_map
from terminal_pareto.subtree_explore import (
    SMALL_TYPES,
    _merge_small,
    decompose_global_front_by_type,
    decompose_global_front_retention_by_type,
)

OUT = Path(__file__).resolve().parent / "output" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"
FULL_RETENTION_CACHE = ANALYSIS_OUT / "global_front_retention_by_cell_type.csv"
KEYPOINT_CACHE = ANALYSIS_OUT / "global_front_by_cell_type.csv"

# Shared row order: strongest fate-specific contrasts first, abundance-heavy
# types next, then smaller responses. This order is fixed across both panels.
ORDER = ["programmed_death", "epithelium", "neuron", "muscle",
         "alimentary", "glial", "other (n<=4)"]

KEYPOINT_COLORS = {
    "spatial_opt": ps.SEMANTIC_COLORS["travel"],
    "max_er": ps.COLORS["blue"],
    "expr_opt": ps.SEMANTIC_COLORS["cell_state"],
}
KEYPOINT_LABELS = {
    "spatial_opt": "Travel-minimizing assignment",
    "max_er": "Maximum-retention compromise",
    "expr_opt": "Cell-state-minimizing assignment",
}
KEYPOINT_MARKERS = {"spatial_opt": "o", "max_er": "D", "expr_opt": "o"}

# Heatmap rendering parameters only. These do not affect the assignments,
# keypoints, statistics, or cached retention values used by the analysis.
HEATMAP_BINS = 120
HEATMAP_SMOOTH_BANDWIDTH = 0.018


def compute_cell_type_caches(iteration=300):
    """Evaluate full-front retention and the three-keypoint decomposition."""
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    protein_exp = dl.load_protein_expression()
    prot_sel = dl.load_prot_sel()
    v_prot = [name for name in valid_ce if name in protein_exp.index]
    gp_map = pe.build_grandparent_map(lineage)
    tn, tp = dl.collect_terminals(lineage, v_prot)
    type_map = build_type_map(tn)
    full = decompose_global_front_retention_by_type(
        lineage, v_prot, tn, tp, xyz_ce, protein_exp, prot_sel, type_map,
        tree_index, gp_map, iteration=iteration)
    keypoints = decompose_global_front_by_type(
        lineage, v_prot, tn, tp, xyz_ce, protein_exp, prot_sel, type_map,
        tree_index, gp_map, iteration=iteration)
    full.to_csv(FULL_RETENTION_CACHE, index=False)
    keypoints.to_csv(KEYPOINT_CACHE, index=False)
    return full, keypoints


def merge_small_retention(full_df):
    """Merge n<=4 fate categories within every distinct front assignment."""
    small = full_df["type"].isin(SMALL_TYPES) | full_df["small"].astype(bool)
    retained = full_df.loc[~small].copy()
    merged_rows = []
    for front_index, group in full_df.loc[small].groupby("front_index"):
        first = group.iloc[0]
        n = int(group["n"].sum())
        merged_rows.append(dict(
            front_index=int(front_index),
            sweep_index=int(first["sweep_index"]),
            alpha=float(first["alpha"]),
            u=float(first["u"]),
            travel_sigma=float(first["travel_sigma"]),
            state_sigma=float(first["state_sigma"]),
            is_max_er=bool(first["is_max_er"]),
            type="other (n<=4)",
            n=n,
            small=True,
            er=float(np.average(group["er"], weights=group["n"])),
        ))
    return pd.concat([retained, pd.DataFrame(merged_rows)], ignore_index=True)


def ordered_types(df):
    present = set(df["type"])
    return [cell_type for cell_type in ORDER if cell_type in present]


def display_type_name(cell_type):
    """Publication label without repeating the small-group threshold."""
    if cell_type == "other (n<=4)":
        return "other"
    return cell_type.replace("_", " ")


def smooth_retention_profiles(u, values, n_bins=HEATMAP_BINS,
                              bandwidth=HEATMAP_SMOOTH_BANDWIDTH):
    """Resample irregular profiles and apply narrow display-only averaging.

    Linear resampling prevents densely sampled sweep regions from receiving
    extra visual weight. Gaussian averaging is then performed on the uniform
    canonical-position grid. Edge padding avoids artificial endpoint decay.
    """
    u = np.asarray(u, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(u)
    u = u[order]
    values = values[:, order]

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    resampled = np.vstack([
        np.interp(centers, u, profile)
        for profile in values
    ])

    sigma_bins = bandwidth * n_bins
    radius = max(1, int(np.ceil(3.0 * sigma_bins)))
    offsets = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (offsets / sigma_bins) ** 2)
    kernel /= kernel.sum()
    smoothed = np.vstack([
        np.convolve(np.pad(profile, radius, mode="edge"), kernel,
                    mode="valid")
        for profile in resampled
    ])
    return edges, np.clip(smoothed, 0.0, 1.0)


def plot_retention_heatmap(full_df):
    """Panel A: natural-edge retention across every distinct assignment."""
    data = merge_small_retention(full_df)
    order = ordered_types(data)
    fronts = (data.drop_duplicates("front_index")
              .sort_values("u").reset_index(drop=True))
    front_ids = fronts["front_index"].to_numpy()
    u = fronts["u"].to_numpy(dtype=float)
    observed_values = np.vstack([
        (data[data["type"] == cell_type]
         .set_index("front_index").reindex(front_ids)["er"].to_numpy())
        for cell_type in order
    ])
    display_edges, display_values = smooth_retention_profiles(
        u, observed_values)
    # Anchor the color scale at zero, but use the full palette over the
    # observed display range rather than the unused 0--1 interval. Rounding
    # upward to 0.05 keeps the colorbar stable and guarantees no clipping.
    color_max = min(1.0, max(0.10,
                    np.ceil(float(np.nanmax(display_values)) / 0.05) * 0.05))

    fig, ax = plt.subplots(figsize=(7.15, 2.65))
    mesh = ax.imshow(
        display_values, extent=(display_edges[0], display_edges[-1],
                                len(order), 0),
        aspect="auto", interpolation="none",
        cmap=EDGE_RETENTION_CMAP, vmin=0.0, vmax=color_max,
        rasterized=True,
    )
    ax.set_ylim(len(order), 0)
    ax.set_xlim(0, 1)
    ax.set_yticks(np.arange(len(order)) + 0.5)
    n_by_type = (data.drop_duplicates("type").set_index("type")["n"].to_dict())
    ax.set_yticklabels([
        f"{display_type_name(cell_type)}  (n={n_by_type[cell_type]})"
        for cell_type in order
    ])
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xlabel("Position along the Pareto front (u)\n"
                  "Travel optimum (u = 0) → Cell-state optimum (u = 1)")
    ax.set_title("Natural-edge retention across the Pareto front",
                 loc="left", pad=9)

    max_rows = fronts[fronts["is_max_er"].astype(bool)]
    if not max_rows.empty:
        max_u = float(max_rows.iloc[0]["u"])
        ax.axvline(max_u, color=KEYPOINT_COLORS["max_er"], lw=1.2,
                   ls=(0, (3, 2)), zorder=4)
        ax.scatter([max_u], [-0.12], marker="D", s=19,
                   facecolor=KEYPOINT_COLORS["max_er"], edgecolor="#222222",
                   lw=0.45, clip_on=False, zorder=5)
        max_handle = Line2D(
            [0], [0], marker="D", ls=(0, (3, 2)), lw=1.0,
            color=KEYPOINT_COLORS["max_er"],
            markerfacecolor=KEYPOINT_COLORS["max_er"],
            markeredgecolor="#222222", markeredgewidth=0.45,
            label="Maximum retention",
        )
        ax.legend(handles=[max_handle], loc="lower right",
                  bbox_to_anchor=(1.0, 1.015), fontsize=6.2,
                  borderaxespad=0, handlelength=2.2)

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.022, pad=0.018)
    cbar.set_label("Natural edges retained", fontsize=7, labelpad=2)
    cbar.set_ticks([0.0, color_max / 2.0, color_max])
    cbar.ax.set_yticklabels([
        "0", f"{color_max / 2.0:.2f}", f"{color_max:.2f}"
    ])
    cbar.ax.tick_params(labelsize=6)
    ax.tick_params(axis="y", length=0)
    ax.grid(False)
    fig.text(0.012, 0.965, "A", fontsize=10, fontweight="bold",
             ha="left", va="top")
    fig.subplots_adjust(left=0.22, right=0.94, top=0.84, bottom=0.25)
    ps.save_figure(fig, OUT / "fig6A_ce_retention_heatmap")
    plt.close(fig)


def marker_area(n, max_n):
    """Compressed cell-count encoding for endpoint-ledger markers."""
    return 18.0 + 42.0 * np.sqrt(float(n) / float(max_n))


def plot_endpoint_ledger(keypoint_df):
    """Supplement: paired signed changes for each objective and cell type."""
    data = _merge_small(keypoint_df)
    order = ordered_types(data)
    y = np.arange(len(order))
    n_by_type = (data[data["keypoint"] == "expr_opt"]
                 .set_index("type")["n"].to_dict())
    max_n = max(n_by_type.values())

    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.15), sharey=True,
                             gridspec_kw={"wspace": 0.12})
    panels = [
        (axes[0], "travel_sigma_per_cell",
         "Travel-distance change per cell (σ)", False,
         "Travel reduction", "Travel increase",
         ps.SEMANTIC_COLORS["travel"]),
        (axes[1], "state_sigma_per_cell",
         "Cell-state-distance change per cell (σ)", True,
         "Cell-state increase", "Cell-state reduction",
         ps.SEMANTIC_COLORS["cell_state"]),
    ]
    indexed = {kp: data[data["keypoint"] == kp].set_index("type")
               for kp in KEYPOINT_COLORS}

    for ax, column, xlabel, invert_axis, left_header, right_header, cue_color in panels:
        ax.axvline(0, color="#777777", lw=0.75, ls=":", zorder=1)
        all_values = []
        for yi, cell_type in zip(y, order):
            endpoint_values = [indexed[kp].loc[cell_type, column]
                               for kp in ("spatial_opt", "expr_opt")]
            ax.plot(endpoint_values, [yi, yi], color="#A0A0A0",
                    lw=0.8, zorder=2)
            for kp in ("spatial_opt", "max_er", "expr_opt"):
                value = float(indexed[kp].loc[cell_type, column])
                all_values.append(value)
                ax.scatter(
                    [value], [yi], marker=KEYPOINT_MARKERS[kp],
                    s=marker_area(n_by_type[cell_type], max_n),
                    facecolor=KEYPOINT_COLORS[kp], edgecolor="#222222",
                    lw=0.55, zorder=4 if kp == "max_er" else 3,
                )
        lo, hi = min(all_values + [0.0]), max(all_values + [0.0])
        span = max(hi - lo, 0.05)
        pad = 0.08 * span
        # A restrained background cue makes the sign convention visible
        # before the reader inspects individual markers. Negative change is
        # always a reduction, even though the cell-state axis is inverted to
        # preserve travel-to-cell-state reading order.
        ax.axvspan(lo - pad, 0, color=cue_color, alpha=0.055, zorder=0)
        ax.axvspan(0, hi + pad, color="#B8B8B8", alpha=0.075, zorder=0)
        if invert_axis:
            # Preserve the established Pareto-front reading direction:
            # travel optimum at left, cell-state optimum at right. Numeric
            # signs consequently run in the opposite direction in this panel.
            ax.set_xlim(hi + pad, lo - pad)
        else:
            ax.set_xlim(lo - pad, hi + pad)
        ax.set_xlabel(xlabel)
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.grid(True, axis="x", alpha=0.3)
        left_color = cue_color if "reduction" in left_header.lower() else "#666666"
        right_color = cue_color if "reduction" in right_header.lower() else "#666666"
        ax.text(0.01, 1.015, left_header, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=6.2, color=left_color,
                fontweight="semibold" if left_color == cue_color else "normal")
        ax.text(0.99, 1.015, right_header, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=6.2, color=right_color,
                fontweight="semibold" if right_color == cue_color else "normal")
        ax.annotate("", xy=(0.66, 1.035), xytext=(0.34, 1.035),
                    xycoords="axes fraction", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", color="#999999",
                                    lw=0.6, mutation_scale=6))

    axes[0].set_yticks(y)
    axes[0].set_yticklabels([
        f"{display_type_name(cell_type)}  (n={n_by_type[cell_type]})"
        for cell_type in order
    ])
    axes[0].tick_params(axis="y", length=0)
    axes[1].tick_params(axis="y", left=False, labelleft=False)

    handles = [
        Line2D([0], [0], marker=KEYPOINT_MARKERS[kp], ls="", ms=5.2,
               markerfacecolor=KEYPOINT_COLORS[kp], markeredgecolor="#222222",
               markeredgewidth=0.5, label=KEYPOINT_LABELS[kp])
        for kp in ("spatial_opt", "max_er", "expr_opt")
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.57, 0.965), fontsize=6.4,
               handletextpad=0.35, columnspacing=1.1)
    fig.suptitle("Fate-specific consequences of moving along the Pareto front",
                 x=0.57, y=0.995, fontsize=9.2, fontweight="semibold")
    fig.text(0.57, 0.045,
             "Both panels run from the travel optimum (left) to the "
             "cell-state optimum (right). Dotted zero = natural lineage; "
             "colored shading = reduction; grey shading = increase.",
             ha="center", va="bottom", fontsize=6.2, color="#555555")
    fig.subplots_adjust(left=0.22, right=0.98, top=0.82, bottom=0.18)
    ps.save_figure(fig, OUT / "figS2_ce_cell_type_cost_gain_panel")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iteration", type=int, default=300,
                        help="Existing Pareto sweep resolution (default 300).")
    parser.add_argument("--reuse-cache", action="store_true",
                        help="Reuse the full-front retention CSV if present.")
    # Kept for command-line compatibility with the former renderer. The cell-
    # type analysis is global and does not use a subtree-size threshold.
    parser.add_argument("--min-cells", type=int, default=12,
                        help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    ps.configure()

    if (args.reuse_cache and FULL_RETENTION_CACHE.exists()
            and KEYPOINT_CACHE.exists()):
        full_df = pd.read_csv(FULL_RETENTION_CACHE)
        keypoint_df = pd.read_csv(KEYPOINT_CACHE)
        print(f"Reused {len(full_df)} full-front and "
              f"{len(keypoint_df)} keypoint cell-type rows.")
    else:
        full_df, keypoint_df = compute_cell_type_caches(
            iteration=args.iteration)
        print(f"Wrote {len(full_df)} full-front rows to "
              f"{FULL_RETENTION_CACHE} and {len(keypoint_df)} keypoint rows "
              f"to {KEYPOINT_CACHE}.")

    plot_retention_heatmap(full_df)
    plot_endpoint_ledger(keypoint_df)
    print("Figure 6A and Figure S2 panels written to", OUT)


if __name__ == "__main__":
    main()
