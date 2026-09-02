"""Generate the C. elegans protein full-tree Pareto publication figures."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from full_tree_pareto import publication_analysis as analysis
from terminal_pareto import plot_style as ps


OUT = analysis.PUBLICATION_ROOT
EDGE_RETENTION_CMAP = LinearSegmentedColormap.from_list(
    "full_tree_edge_retention",
    ["#17365D", ps.COLORS["blue"], "#72C7EC"],
)
LAYER_LINE = "#326A86"
SLOT_NULL_COLOR = ps.SEMANTIC_COLORS["first_cousin_null"]
HEURISTIC_COLORS = {
    "Bottom-up by layer": "#56B4E9",
    "Layerwise assignment": "#0072B2",
    "Paired bottom-up": "#CC79A7",
    "Degree-constrained spanning forest": "#009E73",
    "Top-down rebuild": "#D55E00",
    "Terminal-only rebuild": "#E69F00",
}
NULL_COLORS = {
    "First-cousin shuffle": ps.SEMANTIC_COLORS["first_cousin_null"],
    "Internal-layer shuffle": "#A6611A",
    "Full assignment shuffle": "#7B4F2C",
    "Random rebuild": "#6F3B5C",
    "Parametric Brownian reference": "#6B6ECF",
}


def _save(fig, stem: str, dpi: int = 400) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(
        OUT / f"{stem}.png",
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
    )


def _padded_limits(values: np.ndarray, include_zero: bool = True) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    low = float(np.nanmin(values))
    high = float(np.nanmax(values))
    if include_zero:
        low = min(low, 0.0)
        high = max(high, 0.0)
    span = max(high - low, 1.0)
    return low - 0.08 * span, high + 0.08 * span


def _plot_retention_front(ax, frame: pd.DataFrame, norm: Normalize):
    frame = frame.sort_values("alpha_travel")
    nondominated = frame[frame["sampled_nondominated"]]
    x = frame["travel_standardized"].to_numpy(float)
    y = frame["cell_state_standardized"].to_numpy(float)
    ax.plot(x, y, color=LAYER_LINE, lw=1.15, alpha=0.55, zorder=2)
    scatter = ax.scatter(
        nondominated["travel_standardized"],
        nondominated["cell_state_standardized"],
        c=nondominated["edge_retention"],
        cmap=EDGE_RETENTION_CMAP,
        norm=norm,
        s=13,
        edgecolors="none",
        zorder=3,
    )
    natural = ax.scatter(
        [0], [0], marker="X", s=43, color=ps.COLORS["black"],
        edgecolor="white", lw=0.45, zorder=7,
    )
    max_retention = frame[frame["edge_retention"] == frame["edge_retention"].max()]
    marker_row = max_retention.iloc[len(max_retention) // 2]
    max_marker = ax.scatter(
        [marker_row["travel_standardized"]],
        [marker_row["cell_state_standardized"]],
        marker="o",
        s=33,
        facecolor=EDGE_RETENTION_CMAP(norm(marker_row["edge_retention"])),
        edgecolor=ps.COLORS["black"],
        lw=0.9,
        zorder=6,
    )
    ax.axhline(0, color="#888888", lw=0.55, ls=":", zorder=0)
    ax.axvline(0, color="#888888", lw=0.55, ls=":", zorder=0)
    ax.set_xlim(_padded_limits(x))
    ax.set_ylim(_padded_limits(y))
    ax.grid(True, alpha=0.27)
    return scatter, natural, max_marker


def plot_layerwise_rounds(fronts: pd.DataFrame, manifest: pd.DataFrame):
    """Figure 7A: all bottom-up contraction-round Pareto fronts."""
    rounds = fronts[fronts["scope"] != "aggregate"]
    norm = Normalize(vmin=0.0, vmax=1.0)
    fig, axes = plt.subplots(
        2,
        4,
        figsize=(7.15, 4.25),
        constrained_layout=True,
        sharex=False,
        sharey=False,
    )
    axes_flat = axes.ravel()
    last_scatter = None
    for round_index, ax in enumerate(axes_flat, start=1):
        frame = rounds[rounds["round"] == round_index]
        row = manifest[manifest["round"] == round_index].iloc[0]
        last_scatter, _, _ = _plot_retention_front(ax, frame, norm)
        location = "bottom" if round_index == 1 else (
            "top" if round_index == 8 else ""
        )
        location_text = f" · {location}" if location else ""
        ax.set_title(
            f"Round {round_index}{location_text}\n{int(row['edges'])} edges",
            fontsize=8.1,
            pad=3,
        )
    cbar = fig.colorbar(
        last_scatter,
        ax=axes_flat.tolist(),
        fraction=0.025,
        pad=0.015,
        aspect=30,
    )
    cbar.set_label("Natural edges retained")
    fig.supxlabel(
        "Travel distance (parent-slot-null σ; natural lineage = 0)",
        fontsize=8.2,
    )
    fig.supylabel(
        "Cell-state distance (parent-slot-null σ; natural lineage = 0)",
        fontsize=8.2,
    )
    fig.text(
        0.5,
        1.012,
        "Bottom of measured tree  →  top of measured tree",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color=ps.COLORS["grey"],
    )
    fig.text(0.003, 1.015, "A", fontsize=11, fontweight="bold", va="bottom")
    _save(fig, "fig7A_ce_full_tree_layerwise_rounds")
    return fig


def plot_layerwise_aggregate(fronts: pd.DataFrame, nulls: pd.DataFrame):
    """Figure 7B: aggregate of the eight independently optimized rounds."""
    frame = fronts[fronts["scope"] == "aggregate"]
    null = nulls[nulls["scope"] == "aggregate"]
    norm = Normalize(vmin=0.0, vmax=1.0)
    fig, ax = plt.subplots(figsize=(7.15, 3.85))
    fig.subplots_adjust(left=0.115, right=0.87, bottom=0.18, top=0.94)
    scatter, natural, max_marker = _plot_retention_front(ax, frame, norm)

    inset = ax.inset_axes([0.67, 0.67, 0.30, 0.27])
    inset.scatter(
        null["travel_standardized"],
        null["cell_state_standardized"],
        s=5,
        color=SLOT_NULL_COLOR,
        alpha=0.16,
        edgecolors="none",
        rasterized=True,
    )
    inset.scatter(
        [null["travel_standardized"].mean()],
        [null["cell_state_standardized"].mean()],
        marker="+",
        s=30,
        color=SLOT_NULL_COLOR,
        lw=1.0,
        zorder=4,
    )
    inset.set_title("Independent parent-slot shuffle", fontsize=6.5, pad=2)
    inset.set_xlabel("Travel (σ)", fontsize=5.7, labelpad=1)
    inset.set_ylabel("Cell state (σ)", fontsize=5.7, labelpad=1)
    inset.tick_params(labelsize=5.2, length=2, pad=1)
    inset.grid(True, alpha=0.22)

    null_handle = Line2D(
        [0], [0], marker="o", ls="", markersize=5,
        markerfacecolor=SLOT_NULL_COLOR, markeredgecolor="none", alpha=0.7,
        label="Parent-slot shuffle (inset)",
    )
    ax.legend(
        [natural, max_marker, null_handle],
        ["Natural lineage", "Maximum edge retention", "Parent-slot shuffle (inset)"],
        loc="upper left",
        bbox_to_anchor=(0.67, 0.55),
        fontsize=6.4,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.9,
    )
    ax.set_xlabel("Travel distance (aggregate parent-slot-null σ; natural lineage = 0)")
    ax.set_ylabel(
        "Cell-state distance\n(aggregate parent-slot-null σ; natural lineage = 0)"
    )
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.042, pad=0.025)
    cbar.set_label("Natural edges retained")
    ax.text(-0.14, 1.04, "B", transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")
    _save(fig, "fig7B_ce_full_tree_layerwise_aggregate")
    return fig


def _null_inset(ax, frame: pd.DataFrame, models: list[str], title: str):
    for model in models:
        data = frame[frame["null_model"] == model]
        ax.scatter(
            data["travel_standardized"],
            data["cell_state_standardized"],
            s=4,
            color=NULL_COLORS[model],
            alpha=0.15,
            edgecolors="none",
            rasterized=True,
        )
        ax.scatter(
            [data["travel_standardized"].mean()],
            [data["cell_state_standardized"].mean()],
            marker="+",
            s=25,
            color=NULL_COLORS[model],
            lw=1.0,
            zorder=4,
        )
    ax.set_title(title, fontsize=6.3, pad=2)
    ax.set_xlabel("Travel (σ)", fontsize=5.5, labelpad=1)
    ax.set_ylabel("Cell state (σ)", fontsize=5.5, labelpad=1)
    ax.tick_params(labelsize=5.0, length=2, pad=1)
    ax.grid(True, alpha=0.22)


def plot_collective(
    heuristics: pd.DataFrame,
    nulls: pd.DataFrame,
):
    """Figure 8: collective sampled front across validated heuristics."""
    fig, ax = plt.subplots(figsize=(7.15, 5.05))
    fig.subplots_adjust(left=0.115, right=0.98, bottom=0.14, top=0.96)

    cousin = nulls[nulls["null_model"] == "First-cousin shuffle"]
    ax.scatter(
        cousin["travel_standardized"],
        cousin["cell_state_standardized"],
        s=8,
        color=NULL_COLORS["First-cousin shuffle"],
        alpha=0.17,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )
    ax.scatter(
        [cousin["travel_standardized"].mean()],
        [cousin["cell_state_standardized"].mean()],
        marker="+",
        s=35,
        color=NULL_COLORS["First-cousin shuffle"],
        lw=1.1,
        zorder=4,
    )
    brownian = nulls[nulls["null_model"] == "Parametric Brownian reference"]
    brownian_display = brownian[brownian["displayed"]]
    ax.scatter(
        brownian_display["travel_standardized"],
        brownian_display["cell_state_standardized"],
        s=8,
        color=NULL_COLORS["Parametric Brownian reference"],
        alpha=0.14,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )
    ax.scatter(
        [brownian["travel_standardized"].mean()],
        [brownian["cell_state_standardized"].mean()],
        marker="+",
        s=35,
        color=NULL_COLORS["Parametric Brownian reference"],
        lw=1.1,
        zorder=4,
    )

    line_handles = []
    for name in HEURISTIC_COLORS:
        frame = heuristics[
            (heuristics["heuristic"] == name)
            & heuristics["within_heuristic_nondominated"]
        ].sort_values("travel_standardized")
        line, = ax.plot(
            frame["travel_standardized"],
            frame["cell_state_standardized"],
            color=HEURISTIC_COLORS[name],
            lw=1.25,
            alpha=0.72,
            label=name,
            zorder=2,
        )
        line_handles.append(line)

    natural = ax.scatter(
        [0], [0], marker="X", s=72, color=ps.COLORS["black"],
        edgecolor="white", lw=0.55, label="Natural lineage", zorder=8,
    )
    ax.axhline(0, color="#888888", lw=0.6, ls=":", zorder=0)
    ax.axvline(0, color="#888888", lw=0.6, ls=":", zorder=0)

    broad_models = [
        "Internal-layer shuffle",
        "Full assignment shuffle",
        "Random rebuild",
    ]
    broad = ax.inset_axes([0.665, 0.70, 0.31, 0.25])
    _null_inset(broad, nulls, broad_models, "Broad random nulls")
    broad.text(
        0.03, 0.08, "Internal layer", transform=broad.transAxes,
        color=NULL_COLORS["Internal-layer shuffle"], fontsize=5.0,
    )
    broad.text(
        0.97, 0.91, "Full assignment", transform=broad.transAxes,
        color=NULL_COLORS["Full assignment shuffle"], fontsize=5.0,
        ha="right", va="top",
    )
    broad.text(
        0.97, 0.81, "Random rebuild", transform=broad.transAxes,
        color=NULL_COLORS["Random rebuild"], fontsize=5.0,
        ha="right", va="top",
    )
    cousin_handle = Line2D(
        [0], [0], marker="o", ls="", markersize=4.5,
        markerfacecolor=NULL_COLORS["First-cousin shuffle"],
        markeredgecolor="none", alpha=0.75,
        label="First-cousin shuffle",
    )
    brownian_handle = Line2D(
        [0], [0], marker="o", ls="", markersize=4.5,
        markerfacecolor=NULL_COLORS["Parametric Brownian reference"],
        markeredgecolor="none", alpha=0.75,
        label="Brownian reference (fixed topology)",
    )
    first_legend = ax.legend(
        [natural, cousin_handle, brownian_handle] + line_handles,
        [
            "Natural lineage",
            "First-cousin shuffle",
            "Brownian reference (fixed topology)",
        ] + list(HEURISTIC_COLORS),
        loc="upper left",
        fontsize=6.0,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.9,
        borderpad=0.3,
        labelspacing=0.22,
        handletextpad=0.45,
    )
    ax.add_artist(first_legend)

    main_x = np.concatenate(
        [
            heuristics.loc[
                heuristics["within_heuristic_nondominated"],
                "travel_standardized",
            ].to_numpy(float),
            cousin["travel_standardized"].to_numpy(float),
            brownian_display["travel_standardized"].to_numpy(float),
            np.array([0.0]),
        ]
    )
    main_y = np.concatenate(
        [
            heuristics.loc[
                heuristics["within_heuristic_nondominated"],
                "cell_state_standardized",
            ].to_numpy(float),
            cousin["cell_state_standardized"].to_numpy(float),
            brownian_display["cell_state_standardized"].to_numpy(float),
            np.array([0.0]),
        ]
    )
    ax.set_xlim(_padded_limits(main_x))
    ax.set_ylim(_padded_limits(main_y))
    ax.set_xlabel("Travel distance (first-cousin-null σ; natural lineage = 0)")
    ax.set_ylabel("Cell-state distance\n(first-cousin-null σ; natural lineage = 0)")
    ax.grid(True, alpha=0.28)
    _save(fig, "fig8_ce_full_tree_collective_panel")
    return fig


def report_findings(outputs: dict[str, pd.DataFrame]) -> None:
    manifest = outputs["ce_full_tree_layer_manifest.csv"]
    heuristics = outputs["ce_full_tree_collective_heuristics.csv"]
    aggregate = manifest[manifest["scope"] == "aggregate"].iloc[0]
    print(
        "Layerwise aggregate: "
        f"{int(aggregate['sampled_unique_solutions'])} unique sampled solutions; "
        f"maximum edge retention {aggregate['maximum_edge_retention']:.3f}; "
        f"natural dominated={bool(aggregate['natural_sampled_dominated'])}"
    )
    print("Collective sampled-front contributions:")
    for name in HEURISTIC_COLORS:
        count = int(
            (
                (heuristics["heuristic"] == name)
                & heuristics["collective_nondominated"]
            ).sum()
        )
        print(f"  {name}: {count}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Write and validate caches without rendering panels.",
    )
    args = parser.parse_args()
    ps.configure()
    outputs = analysis.build_publication_caches()
    report_findings(outputs)
    if args.analysis_only:
        return
    plot_layerwise_rounds(
        outputs["ce_full_tree_layerwise_fronts.csv"],
        outputs["ce_full_tree_layer_manifest.csv"],
    )
    plot_layerwise_aggregate(
        outputs["ce_full_tree_layerwise_fronts.csv"],
        outputs["ce_full_tree_layerwise_nulls.csv"],
    )
    plot_collective(
        outputs["ce_full_tree_collective_heuristics.csv"],
        outputs["ce_full_tree_collective_nulls.csv"],
    )
    print(f"Publication panels written to {OUT}")


if __name__ == "__main__":
    main()
