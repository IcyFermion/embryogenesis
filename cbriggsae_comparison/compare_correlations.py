"""
Core comparison: match timepoints between CSV and SBD datasets,
interpolate SBD cell positions, then compute pairwise distance
correlations between datasets.

Key insight: cell overlap is maximized at mid-late timepoints
(CSV t≈140-150, SBD frames≈300-350), not at the very end.
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from parse_data import load_all_data

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def interpolate_cell_positions(cells, target_frame):
    """Linearly interpolate each cell's (x, y) position at target_frame.

    Only interpolates if target_frame is within the cell's observed frame range.
    """
    positions = {}
    for c in cells:
        frames = c["frame"]
        if len(frames) < 2:
            continue
        f_min, f_max = frames.min(), frames.max()
        if not (f_min <= target_frame <= f_max):
            continue
        x_interp = np.interp(target_frame, frames, c["x"])
        y_interp = np.interp(target_frame, frames, c["y"])
        positions[c["cell"]] = (x_interp, y_interp)
    return positions


def get_csv_positions(csv_df, timepoint):
    """Get (x, y) positions from CSV at a given timepoint."""
    subset = csv_df[csv_df.time == timepoint]
    return {row["cell"]: (row["x"], row["y"]) for _, row in subset.iterrows()}


def find_best_matching_frame(cells, csv_cells, max_f, step=10):
    """Find the SBD frame that maximizes cell overlap with a CSV timepoint."""
    best_frame, best_count = 0, 0
    for frame in range(50, int(max_f), step):
        cells_at_frame = set()
        for c in cells:
            if c["frame"].min() <= frame <= c["frame"].max():
                cells_at_frame.add(c["cell"])
        common = len(cells_at_frame & csv_cells)
        if common > best_count:
            best_count = common
            best_frame = frame
    return best_frame, best_count


def compute_pairwise_distance_matrix(positions, cell_list):
    """Compute pairwise Euclidean distance matrix for ordered cell_list."""
    coords = np.array([positions[c] for c in cell_list])
    n = len(cell_list)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sqrt(np.sum((coords[i] - coords[j]) ** 2))
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d
    return dist_matrix


def compare_at_timepoint(csv_df, sbd_cells, csv_t, sbd_frame, label):
    """Compare pairwise distances at a matched timepoint pair."""
    csv_pos = get_csv_positions(csv_df, csv_t)
    sbd_pos = interpolate_cell_positions(sbd_cells, sbd_frame)

    common_cells = sorted(set(csv_pos.keys()) & set(sbd_pos.keys()))

    if len(common_cells) < 10:
        return {"common_cells": common_cells, "n": len(common_cells),
                "csv_t": csv_t, "sbd_frame": sbd_frame, "label": label}

    csv_dist = compute_pairwise_distance_matrix(csv_pos, common_cells)
    sbd_dist = compute_pairwise_distance_matrix(sbd_pos, common_cells)

    # Extract upper triangle (excluding diagonal)
    triu_idx = np.triu_indices_from(csv_dist, k=1)
    csv_flat = csv_dist[triu_idx]
    sbd_flat = sbd_dist[triu_idx]

    # Correlations
    spear_r, spear_p = spearmanr(csv_flat, sbd_flat)
    pear_r, pear_p = pearsonr(csv_flat, sbd_flat)
    log_pear_r, log_pear_p = pearsonr(np.log1p(csv_flat), np.log1p(sbd_flat))

    return {
        "common_cells": common_cells,
        "n": len(common_cells),
        "csv_t": csv_t,
        "sbd_frame": sbd_frame,
        "label": label,
        "csv_dist": csv_dist,
        "sbd_dist": sbd_dist,
        "csv_flat": csv_flat,
        "sbd_flat": sbd_flat,
        "spearman_r": spear_r,
        "spearman_p": spear_p,
        "pearson_r": pear_r,
        "pearson_p": pear_p,
        "log_pearson_r": log_pear_r,
        "log_pearson_p": log_pear_p,
    }


def run_full_comparison():
    """Run the full comparison pipeline with optimal timepoint matching."""
    csv_df, sbd2_cells, sbd3_cells = load_all_data()
    active2 = [c for c in sbd2_cells if c["active"]]
    active3 = [c for c in sbd3_cells if c["active"]]

    max_f2 = int(max(c["frame"].max() for c in active2))
    max_f3 = int(max(c["frame"].max() for c in active3))

    # ================================================================
    # Part A: Correlation across multiple timepoints (stability analysis)
    # ================================================================
    print("=" * 60)
    print("Part A: Correlation stability across time")
    print("=" * 60)

    sweep_results = []
    for csv_t in range(10, 201, 20):
        csv_t_cells = set(csv_df[csv_df.time == csv_t].cell.unique())

        for cells, sbd_label, max_f in [(active2, "SBD2", max_f2),
                                         (active3, "SBD3", max_f3)]:
            best_f, best_c = find_best_matching_frame(cells, csv_t_cells, max_f)
            result = compare_at_timepoint(csv_df, cells, csv_t, best_f, sbd_label)
            if result["n"] >= 20:
                sweep_results.append(result)
                if csv_t in [50, 100, 140, 150, 160, 200]:
                    print(f"  CSV t={csv_t} ↔ {sbd_label} f={best_f}: "
                          f"n={result['n']}, Spearman r={result.get('spearman_r', float('nan')):.3f}")

    # ================================================================
    # Part B: Primary comparison at optimal endpoint
    # ================================================================
    print("\n" + "=" * 60)
    print("Part B: Primary comparison (late timepoint, good coverage)")
    print("=" * 60)

    # Use CSV t=150 as sweet spot (75% development, ~243 common cells)
    primary_t = 150
    csv_150 = set(csv_df[csv_df.time == primary_t].cell.unique())

    # Auto-find best frame, or use the 80%-coverage frame
    # The 80%-of-peak-coverage approach from task instructions:
    def find_endpoint_80pct(cells):
        frame_counts = {}
        for c in cells:
            for f in c["frame"]:
                frame_counts[f] = frame_counts.get(f, 0) + 1
        max_count = max(frame_counts.values())
        threshold = max_count * 0.8
        eligible = sorted([f for f, cnt in frame_counts.items() if cnt >= threshold])
        return eligible[-1]

    sbd2_endpoint = find_endpoint_80pct(active2)
    sbd3_endpoint = find_endpoint_80pct(active3)

    # Also find best match for t=150
    best_f2, best_c2 = find_best_matching_frame(active2, csv_150, max_f2)
    best_f3, best_c3 = find_best_matching_frame(active3, csv_150, max_f3)

    print(f"SBD2 80%-coverage frame: {sbd2_endpoint}")
    print(f"SBD2 best match for CSV t=150: frame {best_f2} ({best_c2} common)")
    print(f"SBD3 80%-coverage frame: {sbd3_endpoint}")
    print(f"SBD3 best match for CSV t=150: frame {best_f3} ({best_c3} common)")

    # Run primary comparison with both methods
    primary_results = {}

    # Method 1: Best matching frame
    res1 = compare_at_timepoint(csv_df, active2, primary_t, best_f2, "SBD2")
    res2 = compare_at_timepoint(csv_df, active3, primary_t, best_f3, "SBD3")
    primary_results["best_match"] = {"SBD2": res1, "SBD3": res2}

    print(f"\n=== Method 1: Best-matched frame ===")
    for res in [res1, res2]:
        print(f"  {res['label']}: n={res['n']}, "
              f"Spearman r={res['spearman_r']:.4f} (p={res['spearman_p']:.2e}), "
              f"Pearson r={res['pearson_r']:.4f} (p={res['pearson_p']:.2e}), "
              f"log-Pearson r={res['log_pearson_r']:.4f}")

    # Method 2: 80%-coverage frame
    res3 = compare_at_timepoint(csv_df, active2, primary_t, sbd2_endpoint, "SBD2")
    res4 = compare_at_timepoint(csv_df, active3, primary_t, sbd3_endpoint, "SBD3")
    primary_results["80pct"] = {"SBD2": res3, "SBD3": res4}

    print(f"\n=== Method 2: 80%-coverage frame ===")
    for res in [res3, res4]:
        print(f"  {res['label']}: n={res['n']}, "
              f"Spearman r={res['spearman_r']:.4f} (p={res['spearman_p']:.2e}), "
              f"Pearson r={res['pearson_r']:.4f} (p={res['pearson_p']:.2e})")

    # ================================================================
    # Part C: SBD2 vs SBD3 (same source, different embryos)
    # ================================================================
    print("\n" + "=" * 60)
    print("Part C: SBD2 vs SBD3 comparison (same source, different embryos)")
    print("=" * 60)

    sbd2_pos_best = interpolate_cell_positions(active2, best_f2)
    sbd3_pos_best = interpolate_cell_positions(active3, best_f3)
    common_23 = sorted(set(sbd2_pos_best.keys()) & set(sbd3_pos_best.keys()))

    sbd2_dist_23 = compute_pairwise_distance_matrix(sbd2_pos_best, common_23)
    sbd3_dist_23 = compute_pairwise_distance_matrix(sbd3_pos_best, common_23)

    triu_idx = np.triu_indices_from(sbd2_dist_23, k=1)
    sbd2_flat_23 = sbd2_dist_23[triu_idx]
    sbd3_flat_23 = sbd3_dist_23[triu_idx]

    spear_23, spear_23_p = spearmanr(sbd2_flat_23, sbd3_flat_23)
    pear_23, pear_23_p = pearsonr(sbd2_flat_23, sbd3_flat_23)
    log_pear_23, log_pear_23_p = pearsonr(
        np.log1p(sbd2_flat_23), np.log1p(sbd3_flat_23))

    print(f"  n={len(common_23)} cells")
    print(f"  Spearman r={spear_23:.4f} (p={spear_23_p:.2e})")
    print(f"  Pearson r={pear_23:.4f} (p={pear_23_p:.2e})")
    print(f"  log-Pearson r={log_pear_23:.4f} (p={log_pear_23_p:.2e})")

    sbd23_result = {
        "common_cells": common_23,
        "n": len(common_23),
        "sbd2_flat": sbd2_flat_23,
        "sbd3_flat": sbd3_flat_23,
        "spearman_r": spear_23,
        "pearson_r": pear_23,
        "log_pearson_r": log_pear_23,
    }

    return sweep_results, primary_results, sbd23_result


def plot_results(sweep_results, primary_results, sbd23_result):
    """Generate comprehensive visualization."""
    fig = plt.figure(figsize=(20, 14))

    # ----------------------------------------------------------
    # Panel 1: Correlation vs time (sweep)
    # ----------------------------------------------------------
    ax1 = fig.add_subplot(2, 3, 1)
    sbd2_sweep = [r for r in sweep_results if r["label"] == "SBD2"]
    sbd3_sweep = [r for r in sweep_results if r["label"] == "SBD3"]

    for sweep, color, marker, lbl in [
        (sbd2_sweep, "blue", "o", "CSV↔SBD2"),
        (sbd3_sweep, "red", "s", "CSV↔SBD3"),
    ]:
        ts = [r["csv_t"] for r in sweep]
        rs = [r["spearman_r"] for r in sweep]
        ns = [r["n"] for r in sweep]
        ax1.scatter(ts, rs, c=color, marker=marker, s=[n / 5 for n in ns],
                    alpha=0.6, label=lbl)
    ax1.set_xlabel("CSV timepoint")
    ax1.set_ylabel("Spearman rank correlation")
    ax1.set_title("Pairwise distance correlation vs time\n(marker size ∝ n common cells)")
    ax1.legend()
    ax1.set_ylim(0, 1)
    ax1.grid(True, alpha=0.3)

    # ----------------------------------------------------------
    # Panel 2: SBD2 vs CSV scatter (primary)
    # ----------------------------------------------------------
    ax2 = fig.add_subplot(2, 3, 2)
    res = primary_results["best_match"]["SBD2"]
    hb = ax2.hexbin(res["csv_flat"], res["sbd_flat"], gridsize=50,
                     cmap="viridis", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax2, label="log10(count)")
    max_val = max(res["csv_flat"].max(), res["sbd_flat"].max())
    ax2.plot([0, max_val], [0, max_val], "r--", alpha=0.5, label="y=x")
    ax2.set_xlabel("CSV pairwise distance (px)")
    ax2.set_ylabel("SBD2 pairwise distance (px)")
    ax2.set_title(f"CSV t=150 ↔ SBD2 f={res['sbd_frame']}\n"
                  f"n={res['n']} cells, "
                  f"Spearman r={res['spearman_r']:.3f}, "
                  f"Pearson r={res['pearson_r']:.3f}")
    ax2.legend()

    # ----------------------------------------------------------
    # Panel 3: SBD3 vs CSV scatter (primary)
    # ----------------------------------------------------------
    ax3 = fig.add_subplot(2, 3, 3)
    res = primary_results["best_match"]["SBD3"]
    hb = ax3.hexbin(res["csv_flat"], res["sbd_flat"], gridsize=50,
                     cmap="viridis", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax3, label="log10(count)")
    max_val = max(res["csv_flat"].max(), res["sbd_flat"].max())
    ax3.plot([0, max_val], [0, max_val], "r--", alpha=0.5, label="y=x")
    ax3.set_xlabel("CSV pairwise distance (px)")
    ax3.set_ylabel("SBD3 pairwise distance (px)")
    ax3.set_title(f"CSV t=150 ↔ SBD3 f={res['sbd_frame']}\n"
                  f"n={res['n']} cells, "
                  f"Spearman r={res['spearman_r']:.3f}, "
                  f"Pearson r={res['pearson_r']:.3f}")
    ax3.legend()

    # ----------------------------------------------------------
    # Panel 4: SBD2 vs SBD3 scatter
    # ----------------------------------------------------------
    ax4 = fig.add_subplot(2, 3, 4)
    hb = ax4.hexbin(sbd23_result["sbd2_flat"], sbd23_result["sbd3_flat"],
                     gridsize=50, cmap="viridis", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax4, label="log10(count)")
    max_val = max(sbd23_result["sbd2_flat"].max(), sbd23_result["sbd3_flat"].max())
    ax4.plot([0, max_val], [0, max_val], "r--", alpha=0.5, label="y=x")
    ax4.set_xlabel("SBD2 pairwise distance (px)")
    ax4.set_ylabel("SBD3 pairwise distance (px)")
    ax4.set_title(f"SBD2 vs SBD3 (same source)\n"
                  f"n={sbd23_result['n']} cells, "
                  f"Spearman r={sbd23_result['spearman_r']:.3f}, "
                  f"Pearson r={sbd23_result['pearson_r']:.3f}")
    ax4.legend()

    # ----------------------------------------------------------
    # Panel 5: Distance distributions
    # ----------------------------------------------------------
    ax5 = fig.add_subplot(2, 3, 5)
    res2 = primary_results["best_match"]["SBD2"]
    res3 = primary_results["best_match"]["SBD3"]
    # Normalize to [0,1] for shape comparison
    for data, label, ls in [
        (res2["csv_flat"], "CSV (ref)", "-"),
        (res2["sbd_flat"], "SBD2", "--"),
        (res3["sbd_flat"], "SBD3", "--"),
        (sbd23_result["sbd2_flat"], "SBD2 (alt frame)", ":"),
        (sbd23_result["sbd3_flat"], "SBD3 (alt frame)", ":"),
    ]:
        vals = data / data.max() if len(data) > 0 else data
        ax5.hist(vals, bins=80, density=True, alpha=0.3, label=label,
                 histtype="stepfilled" if "--" in ls else "step",
                 linewidth=1.5, linestyle=ls)
    ax5.set_xlabel("Normalized pairwise distance")
    ax5.set_ylabel("Density")
    ax5.set_title("Distance distribution comparison")
    ax5.legend(fontsize=7)

    # ----------------------------------------------------------
    # Panel 6: Spatial layout with matched positions
    # ----------------------------------------------------------
    ax6 = fig.add_subplot(2, 3, 6)
    res = primary_results["best_match"]["SBD2"]
    csv_pos = get_csv_positions(
        pd.read_csv(
            Path(__file__).parent.parent / "data/c_briggsae/CD140715HLH1cbp1.csv"),
        150)
    sbd2_cells = load_all_data()[1]
    active2 = [c for c in sbd2_cells if c["active"]]
    sbd2_pos = interpolate_cell_positions(active2, res["sbd_frame"])

    common = sorted(set(csv_pos.keys()) & set(sbd2_pos.keys()))
    csv_coords = np.array([csv_pos[c] for c in common])
    sbd_coords = np.array([sbd2_pos[c] for c in common])

    # Normalize both to zero mean, unit variance
    csv_norm = (csv_coords - csv_coords.mean(axis=0)) / csv_coords.std(axis=0)
    sbd_norm = (sbd_coords - sbd_coords.mean(axis=0)) / sbd_coords.std(axis=0)

    ax6.scatter(csv_norm[:, 0], csv_norm[:, 1], c="blue", alpha=0.4, s=8,
                label="CSV (normalized)")
    ax6.scatter(sbd_norm[:, 0], sbd_norm[:, 1], c="red", alpha=0.4, s=8,
                label=f"SBD2 (normalized)")
    ax6.set_xlabel("Normalized X")
    ax6.set_ylabel("Normalized Y")
    ax6.set_title(f"Aligned spatial positions\n"
                  f"({len(common)} cells, centered & scaled)")
    ax6.legend()
    ax6.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "pairwise_comparison.png", dpi=150)
    plt.close()
    print("Saved pairwise_comparison.png")

    # ----------------------------------------------------------
    # Extra: comprehensive sweep plot
    # ----------------------------------------------------------
    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))

    # Number of common cells vs time
    ax = axes2[0]
    for sweep, color, lbl in [
        (sbd2_sweep, "blue", "CSV↔SBD2"),
        (sbd3_sweep, "red", "CSV↔SBD3"),
    ]:
        ts = [r["csv_t"] for r in sweep]
        ns = [r["n"] for r in sweep]
        ax.plot(ts, ns, "-o", color=color, label=lbl, markersize=5)
    ax.set_xlabel("CSV timepoint")
    ax.set_ylabel("Common cells")
    ax.set_title("Cell overlap across time")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Correlation vs n_cells
    ax = axes2[1]
    for sweep, color, marker, lbl in [
        (sbd2_sweep, "blue", "o", "CSV↔SBD2"),
        (sbd3_sweep, "red", "s", "CSV↔SBD3"),
    ]:
        ns = [r["n"] for r in sweep]
        rs = [r["spearman_r"] for r in sweep]
        ax.scatter(ns, rs, c=color, marker=marker, alpha=0.6, label=lbl)
    ax.set_xlabel("Number of common cells")
    ax.set_ylabel("Spearman r")
    ax.set_title("Correlation vs. sample size")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Correlation bar chart summary
    ax = axes2[2]
    labels = []
    values = []
    colors = []
    for method_name, method_results in primary_results.items():
        for sbd_label, res in method_results.items():
            labels.append(f"{sbd_label}\n({method_name})")
            values.append(res["spearman_r"])
            colors.append("steelblue" if sbd_label == "SBD2" else "coral")
    # Add SBD2 vs SBD3
    labels.append("SBD2↔SBD3")
    values.append(sbd23_result["spearman_r"])
    colors.append("purple")

    bars = ax.bar(range(len(labels)), values, color=colors, alpha=0.7)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Spearman rank correlation")
    ax.set_title("Pairwise distance correlation summary")
    ax.set_ylim(0, 1)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:.3f}", ha="center", fontsize=9, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "correlation_summary.png", dpi=150)
    plt.close()
    print("Saved correlation_summary.png")


if __name__ == "__main__":
    sweep_results, primary_results, sbd23_result = run_full_comparison()
    plot_results(sweep_results, primary_results, sbd23_result)
