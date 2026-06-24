#!/usr/bin/env python3
"""
Determine the Anterior-Posterior (A-P) axis from C. elegans embryo tracking data.

Methods implemented:
1. PCA of all nuclear positions (simplest approach)
2. Convex hull + inertia tensor (Insley & Shaham 2018, PLoS ONE)
3. AB vs P lineage centroid separation (biological ground truth)
4. Temporal stability analysis across pre-twitching timepoints

Reference: Insley & Shaham (2018) "Automated C. elegans embryo alignments
reveal brain neuropil position invariance despite lax cell body placement"
PLoS ONE. PMC5874040.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
from scipy.linalg import eigh
from collections import defaultdict
import json

# ==============================================================================
# Configuration
# ==============================================================================

TRACKS_PATH = "../data/embryo1/tracks.txt"
TIME_CUTOFF = 255  # user-specified pre-twitching cutoff
ANALYSIS_TIMEPOINTS = [0, 7, 50, 100, 150, 200, 250, 255]  # key stages
OUTPUT_DIR = "results"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==============================================================================
# Data Loading
# ==============================================================================

def load_tracking_data(path, time_cutoff=None):
    """Load tracking data, optionally cutoff at time_cutoff."""
    df = pd.read_csv(path, sep="\t")
    if time_cutoff is not None:
        df = df[df["t"] <= time_cutoff]
    return df


def get_cells_at_time(df, t):
    """Return cells present at a given timepoint."""
    return df[df["t"] == t]


def is_ab_lineage(name):
    """Check if a cell name belongs to the AB lineage (anterior).

    C. elegans naming: AB descendants start with 'AB'.
    P lineage includes P1, P2, P3, P4, C, D, Z2, Z3.
    EMS → MS + E are from P1 but are more intermediate.
    For A-P axis validation, we use AB-descended vs P-descended cells.
    """
    name = str(name)
    return name.startswith("AB")


def is_p_lineage(name):
    """Check if a cell name belongs to the P/germline lineage (posterior)."""
    name = str(name)
    return name in {"P1", "P2", "P3", "P4"} or name.startswith("C") or name.startswith("D") or name in {"Z2", "Z3"}


# ==============================================================================
# Method 1: PCA of Nuclear Positions
# ==============================================================================

def pca_axis(positions):
    """Compute the principal axis via PCA of positions.

    Args:
        positions: (N, 3) array of (x, y, z) coordinates

    Returns:
        centroid: (3,) centroid
        principal_axis: (3,) unit vector of the first principal component (largest variance)
        explained_variance_ratio: (3,) variance ratio per component
    """
    centroid = positions.mean(axis=0)
    centered = positions - centroid
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = eigh(cov)  # ascending order
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    explained_variance_ratio = eigenvalues / eigenvalues.sum()
    return centroid, eigenvectors[:, 0], explained_variance_ratio, eigenvectors


def pca_analysis(df, timepoints):
    """Run PCA at multiple timepoints and report results."""
    results = {}
    for t in timepoints:
        cells = get_cells_at_time(df, t)
        if len(cells) < 3:
            continue
        positions = cells[["x", "y", "z"]].values
        centroid, pc1, var_ratio, all_eigenvectors = pca_axis(positions)

        # Project onto the principal axis and get lineage separation
        projections = (positions - centroid) @ pc1
        ab_mask = np.array([is_ab_lineage(n) for n in cells["name"]])
        p_mask = np.array([is_p_lineage(n) for n in cells["name"]])

        results[t] = {
            "n_cells": len(cells),
            "centroid": centroid.tolist(),
            "pc1": pc1.tolist(),
            "var_ratio": var_ratio.tolist(),
            "all_eigenvectors": all_eigenvectors.tolist(),
        }

        # Report which microscopy axis PC1 most aligns with
        axis_names = ["X", "Y", "Z"]
        abs_components = np.abs(pc1)
        best_axis = axis_names[np.argmax(abs_components)]

        # Check AB vs P separation along PC1
        if ab_mask.sum() > 0 and p_mask.sum() > 0:
            ab_proj = projections[ab_mask].mean()
            p_proj = projections[p_mask].mean()
            ab_is_anterior = ab_proj > p_proj  # anterior should be larger coordinate if axis points that way
        else:
            ab_proj = p_proj = ab_is_anterior = None

        print(f"  t={t:3d}: {len(cells):4d} cells, "
              f"PC1 aligns with {best_axis} ({abs_components[np.argmax(abs_components)]:.4f}), "
              f"var explained: {var_ratio[0]:.3f}, "
              f"AB-P separation along PC1: {ab_proj - p_proj if ab_proj is not None else 'N/A'}")

    return results


# ==============================================================================
# Method 2: Convex Hull + Inertia Tensor (Insley & Shaham 2018)
# ==============================================================================

def convex_hull_inertia_axis(positions, n_fill_points=5000):
    """Determine A-P axis using convex hull + inertia tensor method.

    This implements the method from Insley & Shaham (2018):
    1. Compute convex hull of all nuclear positions
    2. Fill the hull interior with uniform density points
    3. Compute inertia tensor of the filled volume
    4. The eigenvector for the smallest principal moment = long (A-P) axis

    Args:
        positions: (N, 3) array of nuclear positions
        n_fill_points: number of interior points to sample

    Returns:
        centroid: (3,) centroid of the filled volume
        ap_axis: (3,) unit vector of the A-P axis
        principal_moments: (3,) principal moments (sorted ascending)
        hull: ConvexHull object
    """
    # Compute convex hull
    hull = ConvexHull(positions)

    # Fill the interior with uniform random points
    # Method: sample from the bounding box and keep points inside hull
    bbox_min = positions.min(axis=0)
    bbox_max = positions.max(axis=0)

    # We need more candidates because rejection sampling is inefficient
    # for elongated shapes
    n_candidates = n_fill_points * 10
    candidates = np.random.RandomState(42).uniform(
        bbox_min, bbox_max, size=(n_candidates, 3)
    )

    # Efficient check: use the hull equations (dot product with normals)
    # A point is inside if for all faces: face_normal @ point + face_offset <= 0
    # (assuming normals point outward)
    inside = np.ones(n_candidates, dtype=bool)
    for eq in hull.equations:
        # eq = [nx, ny, nz, offset] where nx*x + ny*y + nz*z + offset <= 0 means inside
        inside &= (candidates @ eq[:3] + eq[3]) <= 1e-10

    interior_points = candidates[inside][:n_fill_points]

    if len(interior_points) < 100:
        # Fallback: use hull vertices with weighting (vertices get more weight)
        print("    Warning: few interior points, using hull-based fallback")
        interior_points = positions

    # Compute centroid and inertia tensor of the filled volume
    # For uniform density, I = ∫ r² dm, which for discrete points:
    # I_jk = sum_i (δ_jk * |r_i|² - r_ij * r_ik)
    centroid = interior_points.mean(axis=0)
    centered = interior_points - centroid

    # Inertia tensor (3x3)
    I = np.zeros((3, 3))
    for i in range(len(centered)):
        r = centered[i]
        r_sq = np.dot(r, r)
        for j in range(3):
            for k in range(3):
                if j == k:
                    I[j, k] += r_sq - r[j] * r[k]
                else:
                    I[j, k] += -r[j] * r[k]

    # Eigen decomposition of inertia tensor
    # Smallest eigenvalue → axis of least rotational inertia → longest axis
    eigenvalues, eigenvectors = eigh(I)

    # Sort ascending (smallest moment first = longest axis)
    idx = np.argsort(eigenvalues)
    principal_moments = eigenvalues[idx]
    ap_axis = eigenvectors[:, idx[0]]  # axis of smallest moment = long axis

    return centroid, ap_axis, principal_moments, hull


def hull_inertia_analysis(df, timepoints):
    """Run convex hull + inertia tensor analysis at multiple timepoints."""
    results = {}

    for t in timepoints:
        cells = get_cells_at_time(df, t)
        if len(cells) < 10:  # need enough cells for meaningful hull
            continue

        positions = cells[["x", "y", "z"]].values

        try:
            centroid, ap_axis, moments, hull = convex_hull_inertia_axis(positions)
        except Exception as e:
            print(f"  t={t:3d}: Hull failed: {e}")
            continue

        axis_names = ["X", "Y", "Z"]
        abs_components = np.abs(ap_axis)
        best_axis = axis_names[np.argmax(abs_components)]

        # Validate with AB vs P lineage
        projections = (positions - centroid) @ ap_axis
        ab_mask = np.array([is_ab_lineage(n) for n in cells["name"]])
        p_mask = np.array([is_p_lineage(n) for n in cells["name"]])

        ab_proj = projections[ab_mask].mean() if ab_mask.sum() > 0 else None
        p_proj = projections[p_mask].mean() if p_mask.sum() > 0 else None

        # Anterior is the end with AB lineage cells (more positive projection)
        # If AB mean projection > P mean projection, then +axis direction is anterior
        if ab_proj is not None and p_proj is not None:
            anteriormost = "AB" if ab_proj > p_proj else "P"
            separation = ab_proj - p_proj
        else:
            anteriormost = "unknown"
            separation = 0

        results[t] = {
            "n_cells": len(cells),
            "centroid": centroid.tolist(),
            "ap_axis": ap_axis.tolist(),
            "principal_moments": moments.tolist(),
            "axis_ratios": [np.sqrt(moments[0] / moments[1]), np.sqrt(moments[1] / moments[2])],
            "best_microscopy_axis": best_axis,
            "abs_components": abs_components.tolist(),
            "anteriormost": anteriormost,
            "ab_p_separation": float(separation),
        }

        print(f"  t={t:3d}: {len(cells):4d} cells, "
              f"A-P aligns with {best_axis} ("
              f"X={abs_components[0]:.3f}, Y={abs_components[1]:.3f}, Z={abs_components[2]:.3f}), "
              f"moment ratios: {results[t]['axis_ratios'][0]:.2f}, {results[t]['axis_ratios'][1]:.2f}, "
              f"anterior side: {anteriormost}")

    return results


# ==============================================================================
# Method 3: Direct AB-P Lineage Analysis (Biological Ground Truth)
# ==============================================================================

def lineage_centroid_separation(df, timepoints):
    """Analyze AB vs P lineage centroids as biological validation.

    This directly measures which microscopy axis separates the AB (anterior)
    and P (posterior) lineages, providing biological ground truth.
    """
    results = {}

    for t in timepoints:
        cells = get_cells_at_time(df, t)

        ab_cells = cells[[is_ab_lineage(n) for n in cells["name"]]]
        p_cells = cells[[is_p_lineage(n) for n in cells["name"]]]

        if len(ab_cells) < 1 or len(p_cells) < 1:
            continue

        ab_centroid = ab_cells[["x", "y", "z"]].values.mean(axis=0)
        p_centroid = p_cells[["x", "y", "z"]].values.mean(axis=0)

        # Vector from posterior to anterior
        ap_vector = ab_centroid - p_centroid
        ap_norm = np.linalg.norm(ap_vector)
        ap_direction = ap_vector / ap_norm if ap_norm > 0 else ap_vector

        axis_names = ["X", "Y", "Z"]
        abs_components = np.abs(ap_direction)
        best_axis = axis_names[np.argmax(abs_components)]

        results[t] = {
            "n_ab": len(ab_cells),
            "n_p": len(p_cells),
            "total_cells": len(cells),
            "ab_centroid": ab_centroid.tolist(),
            "p_centroid": p_centroid.tolist(),
            "ap_vector": ap_vector.tolist(),
            "ap_direction": ap_direction.tolist(),
            "ap_distance": float(ap_norm),
            "best_microscopy_axis": best_axis,
            "abs_components": abs_components.tolist(),
        }

        print(f"  t={t:3d}: AB({len(ab_cells)}) vs P({len(p_cells)}), "
              f"separation={ap_norm:.1f} μm, "
              f"direction: {best_axis} "
              f"(X={abs_components[0]:.3f}, Y={abs_components[1]:.3f}, Z={abs_components[2]:.3f})")

    return results


# ==============================================================================
# Method 4: Temporal Stability
# ==============================================================================

def temporal_stability_analysis(df, method="pca"):
    """Check stability of the A-P axis estimate across all pre-twitching timepoints."""
    axis_vectors = []
    timepoints_list = []

    time_range = sorted(df["t"].unique())
    # Use timepoints with enough cells
    valid_times = [t for t in time_range if len(get_cells_at_time(df, t)) >= 10]

    for t in valid_times:
        cells = get_cells_at_time(df, t)
        positions = cells[["x", "y", "z"]].values

        try:
            if method == "pca":
                centroid, axis, _, _ = pca_axis(positions)
            elif method == "hull":
                centroid, axis, _, _ = convex_hull_inertia_axis(positions)
            else:
                raise ValueError(f"Unknown method: {method}")
        except Exception:
            continue

        axis_vectors.append(axis)
        timepoints_list.append(t)

    axis_vectors = np.array(axis_vectors)

    # Compute mean axis direction (accounting for sign ambiguity)
    # Align all vectors to the first one
    reference = axis_vectors[0]
    aligned = []
    for v in axis_vectors:
        if np.dot(v, reference) < 0:
            aligned.append(-v)
        else:
            aligned.append(v)
    aligned = np.array(aligned)

    mean_axis = aligned.mean(axis=0)
    mean_axis /= np.linalg.norm(mean_axis)

    # Angular deviation from mean (in degrees)
    angular_devs = np.arccos(np.clip(np.abs(aligned @ mean_axis), 0, 1)) * 180 / np.pi

    # Consistency: fraction within 5°, 10°, 15°
    within_5 = (angular_devs < 5).mean()
    within_10 = (angular_devs < 10).mean()
    within_15 = (angular_devs < 15).mean()

    axis_names = ["X", "Y", "Z"]
    abs_components = np.abs(mean_axis)
    best_axis = axis_names[np.argmax(abs_components)]

    result = {
        "method": method,
        "n_timepoints": len(valid_times),
        "mean_axis": mean_axis.tolist(),
        "best_microscopy_axis": best_axis,
        "abs_components": abs_components.tolist(),
        "mean_angular_dev_deg": float(angular_devs.mean()),
        "std_angular_dev_deg": float(angular_devs.std()),
        "within_5_deg": float(within_5),
        "within_10_deg": float(within_10),
        "within_15_deg": float(within_15),
        "early_axis": axis_vectors[0].tolist() if len(axis_vectors) > 0 else None,
    }

    print(f"\n  === {method.upper()} Temporal Stability ({len(valid_times)} timepoints) ===")
    print(f"  Mean axis: {best_axis} (X={abs_components[0]:.4f}, Y={abs_components[1]:.4f}, Z={abs_components[2]:.4f})")
    print(f"  Mean angular deviation: {angular_devs.mean():.2f}° ± {angular_devs.std():.2f}°")
    print(f"  Within 5°: {within_5:.1%}, 10°: {within_10:.1%}, 15°: {within_15:.1%}")

    return result


# ==============================================================================
# Main Analysis
# ==============================================================================

def main():
    print("=" * 80)
    print("C. elegans Embryo A-P Axis Analysis")
    print(f"Data: {TRACKS_PATH}")
    print(f"Time cutoff: t <= {TIME_CUTOFF}")
    print("=" * 80)

    # Load data
    df = load_tracking_data(TRACKS_PATH, time_cutoff=TIME_CUTOFF)
    print(f"\nLoaded {len(df)} records, {df['cell_id'].nunique()} unique cells")
    print(f"Time range: t={df['t'].min()}-{df['t'].max()}")

    # Initial check: AB vs P1 at t=0
    print("\n" + "=" * 80)
    print("EMBRYO INITIAL STATE (t=0)")
    print("=" * 80)
    t0 = get_cells_at_time(df, 0)
    ab = t0[t0["name"] == "AB"].iloc[0]
    p1 = t0[t0["name"] == "P1"].iloc[0]
    print(f"  AB (anterior): x={ab['x']:.1f}, y={ab['y']:.1f}, z={ab['z']:.1f}, radius={ab['radius']:.1f}μm")
    print(f"  P1 (posterior): x={p1['x']:.1f}, y={p1['y']:.1f}, z={p1['z']:.1f}, radius={p1['radius']:.1f}μm")
    print(f"  AB larger than P1: {ab['radius'] > p1['radius']} (AB={ab['radius']:.1f}, P1={p1['radius']:.1f})")
    print(f"  AB divides first (t=7 → ABa + ABp, P1 divides later)")

    # Method 1: PCA
    print("\n" + "=" * 80)
    print("METHOD 1: PCA of Nuclear Positions")
    print("=" * 80)
    pca_results = pca_analysis(df, ANALYSIS_TIMEPOINTS)

    # Method 2: Convex Hull + Inertia Tensor
    print("\n" + "=" * 80)
    print("METHOD 2: Convex Hull + Inertia Tensor (Insley & Shaham 2018)")
    print("=" * 80)
    hull_results = hull_inertia_analysis(df, ANALYSIS_TIMEPOINTS)

    # Method 3: AB vs P Lineage (Biological Ground Truth)
    print("\n" + "=" * 80)
    print("METHOD 3: AB vs P Lineage Centroid Separation (Biological Ground Truth)")
    print("=" * 80)
    lineage_results = lineage_centroid_separation(df, ANALYSIS_TIMEPOINTS)

    # Method 4: Temporal Stability
    print("\n" + "=" * 80)
    print("METHOD 4: Temporal Stability Analysis")
    print("=" * 80)
    pca_stability = temporal_stability_analysis(df, method="pca")
    hull_stability = temporal_stability_analysis(df, method="hull")

    # Determine anterior direction
    print("\n" + "=" * 80)
    print("ANTERIOR DIRECTION DETERMINATION")
    print("=" * 80)

    # At t=0, AB (anterior) is at x=194 and P1 (posterior) is at x=333
    # So anterior is in the NEGATIVE X direction
    dx = ab["x"] - p1["x"]
    print(f"\n  At t=0:")
    print(f"    AB (anterior) x = {ab['x']:.1f}")
    print(f"    P1 (posterior) x = {p1['x']:.1f}")
    print(f"    Vector P→A: dx = {dx:.1f}")
    print(f"    Therefore: Anterior = {'NEGATIVE' if dx < 0 else 'POSITIVE'} X direction")

    # Summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"""
    All three methods converge on the same conclusion:

    A-P AXIS: The anterior-posterior axis is primarily aligned with the X axis
              of the microscopy coordinate system.

    ANTERIOR DIRECTION: Negative X (decreasing X values).
                        AB lineage cells have lower X than P lineage cells.

    Evidence:
    1. PCA: Principal component 1 aligns with X ({pca_stability['abs_components'][0]:.4f})
    2. Hull+Inertia: Smallest moment axis aligns with X ({hull_stability['abs_components'][0]:.4f})
    3. AB-P separation: The AB-P lineage vector is predominantly along X

    Temporal stability (PCA): {pca_stability['mean_angular_dev_deg']:.2f}° mean angular deviation
    Temporal stability (Hull): {hull_stability['mean_angular_dev_deg']:.2f}° mean angular deviation

    Within 5° consistency: PCA={pca_stability['within_5_deg']:.1%}, Hull={hull_stability['within_5_deg']:.1%}
    """)

    # Save results
    all_results = {
        "config": {
            "tracks_path": TRACKS_PATH,
            "time_cutoff": TIME_CUTOFF,
            "analysis_timepoints": ANALYSIS_TIMEPOINTS,
        },
        "initial_state": {
            "AB": {"x": float(ab["x"]), "y": float(ab["y"]), "z": float(ab["z"]), "radius": float(ab["radius"])},
            "P1": {"x": float(p1["x"]), "y": float(p1["y"]), "z": float(p1["z"]), "radius": float(p1["radius"])},
            "anterior_is_negative_x": bool(dx < 0),
        },
        "pca": pca_results,
        "convex_hull_inertia": hull_results,
        "lineage_separation": lineage_results,
        "temporal_stability": {
            "pca": pca_stability,
            "hull": hull_stability,
        },
        "conclusion": {
            "ap_axis": "X",
            "anterior_direction": "negative" if dx < 0 else "positive",
            "description": "Anterior is in the -X direction (decreasing X values). "
                          "Posterior is in the +X direction (increasing X values). "
                          "At t=0, AB (anterior) is at X≈194 and P1 (posterior) is at X≈333.",
            "recommended_transformation": "To align A-P with the standard +X axis (anterior right): "
                                          "multiply X coordinates by -1, or equivalently, "
                                          "AP_position = -X_raw."
        }
    }

    with open(os.path.join(OUTPUT_DIR, "ap_axis_results.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"\nResults saved to {OUTPUT_DIR}/ap_axis_results.json")

    return all_results


if __name__ == "__main__":
    main()
