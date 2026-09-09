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
from full_tree_pareto.heuristic_inventory import write_inventory
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
    analysis.PUBLICATION_REFERENCE_LABEL: "#6B6ECF",
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


MAIN_HEURISTICS = ("Layerwise assignment", "Degree-constrained spanning forest")
MAIN_NULLS = ("First-cousin shuffle", analysis.PUBLICATION_REFERENCE_LABEL, "Random rebuild")
# The terminal-only cache is retained for provenance, not plotted: its current
# generator changes internal positions and can omit the final scored pair.
SUPPLEMENT_HEURISTICS = tuple(name for name in HEURISTIC_COLORS if name != "Terminal-only rebuild")


def _plot_selected_comparison(heuristics, nulls, models, null_models, stem, summary=None):
    reference = nulls[nulls["null_model"] == analysis.PUBLICATION_REFERENCE_LABEL]
    if reference.empty or not reference["source"].eq(analysis.PUBLICATION_REFERENCE_SOURCE).all():
        raise AssertionError("Figure 8 requires the separate-clock Gaussian reference")
    if not set(models).issubset(set(heuristics["heuristic"])):
        raise AssertionError("Selected heuristic is absent from the cache")
    if not set(null_models).issubset(set(nulls["null_model"])):
        raise AssertionError("Selected null is absent from the cache")
    fig, ax = plt.subplots(figsize=(7.15, 5.05))
    fig.subplots_adjust(left=0.115, right=0.98, bottom=0.14, top=0.96)
    coordinates = [np.array([[0., 0.]])]
    handles = []
    for name in models:
        frame = heuristics[
            (heuristics["heuristic"] == name) & heuristics["within_heuristic_nondominated"]
        ].sort_values("travel_standardized")
        xy = frame[["travel_standardized", "cell_state_standardized"]].to_numpy(float)
        coordinates.append(xy)
        line, = ax.plot(xy[:, 0], xy[:, 1], color=HEURISTIC_COLORS[name],
                        lw=1.65 if name in MAIN_HEURISTICS else 1.1,
                        ls="-" if name in MAIN_HEURISTICS else "--",
                        alpha=0.95 if name in MAIN_HEURISTICS else 0.75,
                        label=name, zorder=3)
        handles.append(line)
    for name in null_models:
        group = nulls[nulls["null_model"] == name]
        shown = group[group["displayed"]]
        xy = shown[["travel_standardized", "cell_state_standardized"]].to_numpy(float)
        coordinates.append(xy)
        ax.scatter(xy[:, 0], xy[:, 1], s=8, alpha=.14, color=NULL_COLORS[name],
                   edgecolors="none", rasterized=True, zorder=1)
        if summary is None:
            mean = group[["travel_standardized", "cell_state_standardized"]].mean().to_numpy()
        else:
            row = summary.set_index("null_model").loc[name]
            mean = row[["travel_mean_standardized", "cell_state_mean_standardized"]].to_numpy(float)
        coordinates.append(mean[None, :])
        marker = ax.scatter(*mean, marker="+", s=40, color=NULL_COLORS[name],
                            lw=1.15, label=name, zorder=4)
        handles.append(marker)
    natural = ax.scatter(0, 0, marker="X", s=72, color=ps.COLORS["black"],
                         edgecolor="white", lw=.55, label="Natural lineage", zorder=8)
    handles.append(natural)
    ax.axhline(0, color="#888888", lw=.6, ls=":", zorder=0)
    ax.axvline(0, color="#888888", lw=.6, ls=":", zorder=0)
    xy = np.concatenate(coordinates)
    ax.set_xlim(_padded_limits(xy[:, 0]))
    ax.set_ylim(_padded_limits(xy[:, 1]))
    ax.set_xlabel("Travel distance (first-cousin-null σ; natural lineage = 0)")
    ax.set_ylabel("Cell-state distance\n(first-cousin-null σ; natural lineage = 0)")
    ax.grid(True, alpha=.22)
    ax.legend(handles=handles, loc="upper left", fontsize=6.0, frameon=True,
              facecolor="white", edgecolor="none", framealpha=.9,
              borderpad=.3, labelspacing=.25, handletextpad=.45)
    _save(fig, stem)
    return fig


def plot_collective(heuristics, nulls, summary=None):
    """Figure 8: two complementary reconstructions, three clouds, continuous axes."""
    return _plot_selected_comparison(
        heuristics, nulls, MAIN_HEURISTICS, MAIN_NULLS,
        "fig8_ce_full_tree_collective_panel", summary,
    )


def plot_supplementary(heuristics, nulls, summary=None):
    """Five cached comparisons and all nulls; terminal-only is inventory-only."""
    return _plot_selected_comparison(
        heuristics, nulls, SUPPLEMENT_HEURISTICS, tuple(NULL_COLORS),
        "figs_ce_full_tree_heuristics_panel", summary,
    )


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
    print("Heuristic cache inventory (not a common-feasible-set ranking):")
    for name in HEURISTIC_COLORS:
        count = int((heuristics["heuristic"] == name).sum())
        placement = "main" if name in MAIN_HEURISTICS else "supplement" if name in SUPPLEMENT_HEURISTICS else "held; table only"
        print(f"  {name}: {count} cached weights; {placement}")


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
    write_inventory(outputs["ce_full_tree_collective_heuristics.csv"], OUT)
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
        outputs["ce_full_tree_null_summary.csv"],
    )
    plot_supplementary(
        outputs["ce_full_tree_collective_heuristics.csv"],
        outputs["ce_full_tree_collective_nulls.csv"],
        outputs["ce_full_tree_null_summary.csv"],
    )
    print(f"Publication panels written to {OUT}")


if __name__ == "__main__":
    main()
