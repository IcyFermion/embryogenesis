#!/usr/bin/env python3
"""Visualize the A-P axis determination results."""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d, Axes3D
import seaborn as sns

from ap_axis import (
    compute_ap_axis, pca_axis,
    convex_hull_inertia_axis, lineage_separation_axis,
    is_ab_lineage, is_p_lineage, quick_check,
)


def load_tracking_data(path, time_cutoff=None):
    """Load tracking data, optionally cutoff at time_cutoff."""
    df = pd.read_csv(path, sep="\t")
    if time_cutoff is not None:
        df = df[df["t"] <= time_cutoff]
    return df

TRACKS_PATH = "../data/embryo1/tracks.txt"
TIME_CUTOFF = 255
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Style
plt.rcParams.update({
    "figure.dpi": 150,
    "figure.figsize": (8, 6),
    "font.size": 10,
})

# ==============================================================================
# Figure 1: AB vs P1 at t=0 — The Biological Ground Truth
# ==============================================================================

def plot_initial_state(df):
    """Show AB and P1 at t=0 with the A-P vector."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    t0 = df[df["t"] == 0]
    ab = t0[t0["name"] == "AB"].iloc[0]
    p1 = t0[t0["name"] == "P1"].iloc[0]

    projections = [("x", "y"), ("x", "z"), ("y", "z")]
    for ax, (h, v) in zip(axes, projections):
        ax.scatter(ab[h], ab[v], c="C0", s=300, marker="o",
                   label=f"AB (anterior, R={ab['radius']:.0f}μm)", zorder=5,
                   edgecolors="black", linewidth=1)
        ax.scatter(p1[h], p1[v], c="C3", s=150, marker="o",
                   label=f"P1 (posterior, R={p1['radius']:.0f}μm)", zorder=5,
                   edgecolors="black", linewidth=1)

        # Arrow from P1 to AB (posterior → anterior)
        ax.annotate("", xy=(ab[h], ab[v]), xytext=(p1[h], p1[v]),
                    arrowprops=dict(arrowstyle="->", color="C0", lw=2,
                                   connectionstyle="arc3,rad=0"))

        # Mark which direction is anterior
        mid_x = (ab[h] + p1[h]) / 2
        mid_y = (ab[v] + p1[v]) / 2
        ax.annotate("Anterior", xy=(ab[h], ab[v]),
                    xytext=(ab[h] - 30, ab[v] + 20 if v != "z" else ab[v] - 20),
                    fontsize=10, color="C0", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="C0", lw=0.5))
        ax.annotate("Posterior", xy=(p1[h], p1[v]),
                    xytext=(p1[h] + 20, p1[v] - 15),
                    fontsize=10, color="C3", fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="C3", lw=0.5))

        ax.set_xlabel(h.upper())
        ax.set_ylabel(v.upper())
        ax.legend(fontsize=8, loc="best")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)

    fig.suptitle("C. elegans Embryo at t=0: AB (Anterior) vs P1 (Posterior)\n"
                f"AB→P1 vector: dx={ab['x']-p1['x']:.0f}, dy={ab['y']-p1['y']:.0f}, dz={ab['z']-p1['z']:.0f}",
                fontweight="bold", fontsize=13)
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/fig1_initial_state.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1_initial_state.png")


# ==============================================================================
# Figure 2: PCA at multiple timepoints
# ==============================================================================

def plot_pca_timepoints(df, timepoints):
    """Show PCA of nuclear positions at several timepoints."""
    n = len(timepoints)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = axes.flatten()

    axis_names = ["X", "Y", "Z"]

    for idx, t in enumerate(timepoints):
        ax = axes[idx]
        cells = df[df["t"] == t]
        positions = cells[["x", "y", "z"]].values

        # PCA
        centroid, pc1, var_ratio, all_evecs = pca_axis(positions)
        abs_comp = np.abs(pc1)

        # Project onto PC1 and PC2
        proj1 = (positions - centroid) @ pc1
        proj2 = (positions - centroid) @ all_evecs[:, 1]

        # Color by lineage
        colors = []
        for name in cells["name"]:
            if is_ab_lineage(name):
                colors.append("C0")  # AB = anterior
            elif is_p_lineage(name):
                colors.append("C3")  # P = posterior
            else:
                colors.append("C2")  # intermediate

        ax.scatter(proj1, proj2, c=colors, s=5, alpha=0.7)

        # Draw PC1 arrow showing anterior direction
        ab_mask = np.array([is_ab_lineage(n) for n in cells["name"]])
        p_mask = np.array([is_p_lineage(n) for n in cells["name"]])
        if ab_mask.sum() > 0 and p_mask.sum() > 0:
            ab_proj = proj1[ab_mask].mean()
            p_proj = proj1[p_mask].mean()
            # Arrow from P to AB along PC1
            ax.annotate("", xy=(ab_proj, 0), xytext=(p_proj, 0),
                       arrowprops=dict(arrowstyle="->", color="black", lw=2))

        ax.set_title(f"t={t} ({len(cells)} cells)\n"
                    f"PC1 = {axis_names[np.argmax(abs_comp)]} "
                    f"({max(abs_comp):.3f})",
                    fontsize=9)
        ax.set_xlabel("PC1 (A-P axis)")
        ax.set_ylabel("PC2")
        ax.axhline(0, color="gray", lw=0.5, alpha=0.3)
        ax.axvline(0, color="gray", lw=0.5, alpha=0.3)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="C0",
               markersize=8, label="AB lineage (anterior)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="C3",
               markersize=8, label="P lineage (posterior)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="C2",
               markersize=8, label="Intermediate"),
    ]
    axes[-1].legend(handles=legend_elements, loc="center", fontsize=9)
    axes[-1].axis("off")

    fig.suptitle("PCA of Nuclear Positions at Key Timepoints\n"
                "(PC1 = A-P axis, colored by lineage)",
                fontweight="bold", fontsize=14)
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/fig2_pca_timepoints.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig2_pca_timepoints.png")


# ==============================================================================
# Figure 3: A-P Axis Temporal Stability
# ==============================================================================

def plot_temporal_stability(df):
    """Show how the A-P axis direction is stable across time."""
    valid_times = [t for t in sorted(df["t"].unique())
                   if len(df[df["t"] == t]) >= 10]

    n = len(valid_times)
    pc1_components = np.zeros((n, 3))
    hull_components = np.zeros((n, 3))
    cell_counts = np.zeros(n)
    ab_p_separations = np.zeros(n)

    for i, t in enumerate(valid_times):
        cells = df[df["t"] == t]
        positions = cells[["x", "y", "z"]].values
        cell_counts[i] = len(cells)

        # PCA
        _, pc1, _, _ = pca_axis(positions)
        pc1_components[i] = np.abs(pc1)

        # Hull+inertia
        try:
            _, ap_axis, _, _ = convex_hull_inertia_axis(positions)
            if ap_axis is not None:
                hull_components[i] = np.abs(ap_axis)
        except Exception:
            hull_components[i] = np.nan

        # AB-P separation along PCA axis
        ab_mask = np.array([is_ab_lineage(n) for n in cells["name"]])
        p_mask = np.array([is_p_lineage(n) for n in cells["name"]])
        if ab_mask.sum() > 0 and p_mask.sum() > 0:
            proj = (positions - positions.mean(axis=0)) @ pc1
            ab_p_separations[i] = proj[ab_mask].mean() - proj[p_mask].mean()

    # Normalize so first timepoint is positive
    reference = pc1_components[0]
    for i in range(n):
        if np.dot(pc1_components[i], reference) < 0:
            pc1_components[i] *= -1

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top left: PCA components over time
    ax = axes[0, 0]
    ax.plot(valid_times, np.abs(pc1_components[:, 0]), label="|PC1·X|", lw=1.5)
    ax.plot(valid_times, np.abs(pc1_components[:, 1]), label="|PC1·Y|", lw=1.5)
    ax.plot(valid_times, np.abs(pc1_components[:, 2]), label="|PC1·Z|", lw=1.5)
    ax.set_xlabel("Time (t)")
    ax.set_ylabel("Absolute component")
    ax.set_title("PCA: PC1 Alignment with Microscopy Axes")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)

    # Top right: Hull components over time
    ax = axes[0, 1]
    ax.plot(valid_times, np.abs(hull_components[:, 0]), label="|AP·X|", lw=1.5)
    ax.plot(valid_times, np.abs(hull_components[:, 1]), label="|AP·Y|", lw=1.5)
    ax.plot(valid_times, np.abs(hull_components[:, 2]), label="|AP·Z|", lw=1.5)
    ax.set_xlabel("Time (t)")
    ax.set_ylabel("Absolute component")
    ax.set_title("Hull+Inertia: Long Axis Alignment with Microscopy Axes")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)

    # Bottom left: AB-P separation over time
    ax = axes[1, 0]
    mask = np.abs(ab_p_separations) > 1e-10
    ax.plot(np.array(valid_times)[mask], np.abs(ab_p_separations[mask]),
            color="C4", lw=1.5)
    ax.set_xlabel("Time (t)")
    ax.set_ylabel("|AB projection - P projection| (μm)")
    ax.set_title("AB vs P Lineage Separation along PC1")
    ax.grid(True, alpha=0.3)

    # Bottom right: Cell count over time
    ax = axes[1, 1]
    ax.plot(valid_times, cell_counts, color="C5", lw=1.5)
    ax.axvline(x=TIME_CUTOFF, color="red", ls="--", alpha=0.5,
              label=f"Time cutoff (t={TIME_CUTOFF})")
    ax.set_xlabel("Time (t)")
    ax.set_ylabel("Number of cells")
    ax.set_title("Cell Count Progression")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("Temporal Stability of A-P Axis Determination\n"
                f"Data: embryo1, t=0–{TIME_CUTOFF} (pre-twitching)",
                fontweight="bold", fontsize=13)
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/fig3_temporal_stability.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig3_temporal_stability.png")


# ==============================================================================
# Figure 4: 3D View of Embryo with A-P Axis
# ==============================================================================

def plot_3d_embryo(df, t=200):
    """3D scatter plot of nuclear positions at a key timepoint with A-P axis."""
    from mpl_toolkits.mplot3d import Axes3D

    cells = df[df["t"] == t]
    positions = cells[["x", "y", "z"]].values
    centroid, pc1, var_ratio, all_evecs = pca_axis(positions)

    # Determine anterior direction from lineage
    projections = (positions - centroid) @ pc1
    ab_mask = np.array([is_ab_lineage(n) for n in cells["name"]])
    p_mask = np.array([is_p_lineage(n) for n in cells["name"]])

    fig = plt.figure(figsize=(14, 6))

    # 3D View
    ax = fig.add_subplot(121, projection="3d")
    colors = []
    for name in cells["name"]:
        if is_ab_lineage(name):
            colors.append("C0")
        elif is_p_lineage(name):
            colors.append("C3")
        else:
            colors.append("C2")
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
              c=colors, s=3, alpha=0.6)

    # Draw A-P axis line
    axis_len = np.linalg.norm(positions.std(axis=0)) * 3
    line = np.outer(np.linspace(-1, 1, 100), pc1) * axis_len + centroid
    ax.plot(line[:, 0], line[:, 1], line[:, 2], color="black", lw=2, label="A-P axis (PC1)")

    # Mark AB and P centroids
    ab_cells = positions[ab_mask]
    p_cells = positions[p_mask]
    if len(ab_cells) > 0:
        ab_cent = ab_cells.mean(axis=0)
        ax.scatter(*ab_cent, c="C0", s=100, marker="*", edgecolors="black",
                  label="AB centroid (anterior)")
    if len(p_cells) > 0:
        p_cent = p_cells.mean(axis=0)
        ax.scatter(*p_cent, c="C3", s=100, marker="*", edgecolors="black",
                  label="P centroid (posterior)")

    ax.set_xlabel("X (μm)")
    ax.set_ylabel("Y (μm)")
    ax.set_zlabel("Z (μm)")
    ax.set_title(f"3D View at t={t} ({len(cells)} cells)\n"
                f"PC1 X-comp={abs(pc1[0]):.4f}, Y={abs(pc1[1]):.4f}, Z={abs(pc1[2]):.4f}")
    ax.legend(fontsize=7, loc="upper left")

    # Histogram of cells along A-P axis (projected)
    ax2 = fig.add_subplot(122)

    # Use -X as the A-P axis for the histogram (anterior is negative X)
    ap_positions = -positions[:, 0]  # flip X so anterior is positive
    ax2.hist(ap_positions[ab_mask], bins=30, alpha=0.5, color="C0",
            label=f"AB lineage ({ab_mask.sum()} cells)", density=True)
    ax2.hist(ap_positions[p_mask], bins=20, alpha=0.5, color="C3",
            label=f"P lineage ({p_mask.sum()} cells)", density=True)

    ax2.axvline(x=ap_positions[ab_mask].mean(), color="C0", ls="--", lw=1.5,
               label=f"AB mean ({ap_positions[ab_mask].mean():.0f})")
    ax2.axvline(x=ap_positions[p_mask].mean(), color="C3", ls="--", lw=1.5,
               label=f"P mean ({ap_positions[p_mask].mean():.0f})")

    ax2.set_xlabel("A-P Position (-X, μm) → anterior")
    ax2.set_ylabel("Density")
    ax2.set_title(f"Cell Distribution Along A-P Axis (t={t})")
    ax2.legend(fontsize=8)

    fig.suptitle("3D Embryo with Computed A-P Axis",
                fontweight="bold", fontsize=13)
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/fig4_3d_embryo.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig4_3d_embryo.png")


# ==============================================================================
# Main
# ==============================================================================

def main():
    print("Loading data...")
    df = load_tracking_data(TRACKS_PATH, time_cutoff=TIME_CUTOFF)
    print(f"Loaded {len(df)} records (t <= {TIME_CUTOFF})")

    # Run quick check
    qc = quick_check(TRACKS_PATH, TIME_CUTOFF)
    print(f"\nQuick check: A-P axis = {qc['ap_axis']}, "
          f"anterior = {qc['anterior_direction']}, "
          f"consistent = {qc['consistent']}")

    # Generate figures
    print("\nGenerating figures...")

    key_timepoints = [0, 7, 14, 28, 50, 100, 150, 200]

    print("  Figure 1: Initial state...")
    plot_initial_state(df)

    print("  Figure 2: PCA at timepoints...")
    plot_pca_timepoints(df, key_timepoints)

    print("  Figure 3: Temporal stability...")
    plot_temporal_stability(df)

    print("  Figure 4: 3D embryo view...")
    plot_3d_embryo(df, t=200)

    print(f"\nAll figures saved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
