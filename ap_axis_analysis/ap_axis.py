"""
A-P Axis Utility for C. elegans Embryo Tracking Data.

Determines the anterior-posterior (A-P) axis from 3D nuclear tracking data.
The A-P axis is established at first cleavage: AB (anterior, larger cell) and
P1 (posterior, smaller cell). This module implements multiple methods:

1. **PCA** — Principal component of nuclear positions (simplest, fastest)
2. **Convex Hull + Inertia Tensor** — The Insley & Shaham (2018) method
3. **Lineage Centroid** — Direct AB vs P lineage separation (biological ground truth)

Based on analysis of the embryo1/2/3 tracking data, the A-P axis is aligned
with the **microscopy X axis**, with **anterior in the -X direction**.

Reference:
    Insley & Shaham (2018) "Automated C. elegans embryo alignments reveal brain
    neuropil position invariance despite lax cell body placement." PLoS ONE.
    https://pmc.ncbi.nlm.nih.gov/articles/PMC5874040/

Usage:
    >>> from ap_axis import compute_ap_axis, get_ap_position, transform_to_ap
    >>> import pandas as pd
    >>> df = pd.read_csv("data/embryo1/tracks.txt", sep="\\t")
    >>>
    >>> # Get the A-P axis from tracking data
    >>> result = compute_ap_axis(df, method="pca", time_cutoff=255)
    >>> print(f"A-P axis: {result.axis_name}, "
    ...       f"anterior direction: {result.anterior_direction}")
    >>>
    >>> # Transform positions so A-P is aligned with coordinate axes
    >>> df_ap = transform_to_ap(df)
    >>> # Now X_ap increases from posterior to anterior
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, Literal
import json


# ==============================================================================
# Data Structures
# ==============================================================================

@dataclass
class APAaxisResult:
    """Result of A-P axis computation.

    Attributes:
        axis_name: Which microscopy axis the A-P axis aligns with ("X", "Y", or "Z").
        axis_index: Index of that axis in (x, y, z) ordering (0=X, 1=Y, 2=Z).
        axis_vector: Unit vector of the A-P axis in (x, y, z) coordinates.
        anterior_direction: +1 if anterior is in the +direction of axis_name,
                            -1 if anterior is in the -direction.
        method: Which method was used ("pca", "hull_inertia", "lineage").
        centroid: The centroid of all cells used in computation (x, y, z).
        explained_variance_ratio: For PCA, variance explained per component.
        principal_moments: For hull method, principal moments of inertia.
        time_cutoff: The time cutoff used.
        n_cells: Number of cells used.
        n_timepoints: Number of timepoints averaged over.
        metadata: Additional method-specific information.
    """
    axis_name: str  # "X", "Y", or "Z"
    axis_index: int  # 0, 1, 2
    axis_vector: np.ndarray  # (3,) unit vector in (x, y, z) order
    anterior_direction: int  # +1 or -1
    method: str
    centroid: np.ndarray  # (3,)
    explained_variance_ratio: Optional[np.ndarray] = None
    principal_moments: Optional[np.ndarray] = None
    time_cutoff: Optional[int] = None
    n_cells: int = 0
    n_timepoints: int = 1
    metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Serialize to a JSON-compatible dict."""
        d = {
            "axis_name": self.axis_name,
            "axis_index": self.axis_index,
            "axis_vector": self.axis_vector.tolist() if self.axis_vector is not None else None,
            "anterior_direction": self.anterior_direction,
            "method": self.method,
            "centroid": self.centroid.tolist() if self.centroid is not None else None,
            "time_cutoff": self.time_cutoff,
            "n_cells": self.n_cells,
            "n_timepoints": self.n_timepoints,
            "metadata": self.metadata,
        }
        if self.explained_variance_ratio is not None:
            d["explained_variance_ratio"] = self.explained_variance_ratio.tolist()
        if self.principal_moments is not None:
            d["principal_moments"] = self.principal_moments.tolist()
        return d

    def save(self, path: str):
        """Save result to JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str) -> "APAaxisResult":
        """Load result from JSON file."""
        with open(path) as f:
            d = json.load(f)
        return cls(
            axis_name=d["axis_name"],
            axis_index=d["axis_index"],
            axis_vector=np.array(d["axis_vector"]),
            anterior_direction=d["anterior_direction"],
            method=d["method"],
            centroid=np.array(d["centroid"]),
            explained_variance_ratio=(
                np.array(d["explained_variance_ratio"])
                if d.get("explained_variance_ratio") else None
            ),
            principal_moments=(
                np.array(d["principal_moments"])
                if d.get("principal_moments") else None
            ),
            time_cutoff=d.get("time_cutoff"),
            n_cells=d.get("n_cells", 0),
            n_timepoints=d.get("n_timepoints", 1),
            metadata=d.get("metadata", {}),
        )

    def __repr__(self) -> str:
        dir_str = f"{'-' if self.anterior_direction < 0 else '+'}{self.axis_name}"
        return (f"APAaxisResult(axis={self.axis_name}, "
                f"anterior={dir_str}, "
                f"method={self.method}, "
                f"n_cells={self.n_cells})")


# ==============================================================================
# Data Loading
# ==============================================================================

def load_tracking_data(path: str, time_cutoff: Optional[int] = None) -> "pd.DataFrame":
    """Load tracking data from a tab-separated file.

    Args:
        path: Path to tracks.txt file.
        time_cutoff: If provided, filter to t <= time_cutoff.

    Returns:
        DataFrame with columns: t, z, y, x, cell_id, parent_id, track_id, radius, name, div_state.
    """
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    if time_cutoff is not None:
        df = df[df["t"] <= time_cutoff]
    return df


# ==============================================================================
# Lineage Utilities
# ==============================================================================

def is_ab_lineage(name: str) -> bool:
    """Check if a cell belongs to the AB lineage (anterior).

    In C. elegans, AB is the anterior blastomere from the first division.
    All AB descendants have names starting with 'AB'.
    """
    return str(name).startswith("AB")


def is_p_lineage(name: str) -> bool:
    """Check if a cell belongs to the P/germline lineage (posterior).

    Includes P1, P2, P3, P4 and their somatic descendants (C, D)
    plus the primordial germ cells (Z2, Z3).
    Note: EMS, MS, and E are intermediate (from P1 but in the middle of the embryo)
    and are NOT classified as posterior here.
    """
    name = str(name)
    return (name in {"P1", "P2", "P3", "P4"}
            or name.startswith("C")
            or name.startswith("D")
            or name in {"Z2", "Z3"})


# ==============================================================================
# Core Computation Methods
# ==============================================================================

def pca_axis(positions: np.ndarray) -> tuple:
    """Compute the principal axis via PCA.

    Args:
        positions: (N, 3) array of (x, y, z) coordinates.

    Returns:
        centroid: (3,) centroid.
        pc1: (3,) unit vector of PC1 (direction of maximum variance).
        explained_variance_ratio: (3,) variance fraction per component.
        all_eigenvectors: (3, 3) eigenvectors as columns, sorted by eigenvalue descending.
    """
    centroid = positions.mean(axis=0)
    centered = positions - centroid
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    explained_variance_ratio = eigenvalues / eigenvalues.sum()
    return centroid, eigenvectors[:, 0], explained_variance_ratio, eigenvectors


def convex_hull_inertia_axis(
    positions: np.ndarray,
    n_fill_points: int = 5000,
    random_seed: int = 42,
) -> tuple:
    """Determine the long axis via convex hull + inertia tensor.

    Implements the method from Insley & Shaham (2018):
    1. Compute convex hull of nuclear positions.
    2. Fill the interior with uniformly distributed points.
    3. Compute the inertia tensor of the filled volume.
    4. The eigenvector for the *smallest* principal moment = long (A-P) axis.

    Requires at least 10 cells for a meaningful hull.

    Args:
        positions: (N, 3) array of (x, y, z) coordinates.
        n_fill_points: Number of interior points to sample.
        random_seed: Seed for reproducible interior sampling.

    Returns:
        centroid: (3,) centroid of the filled volume.
        ap_axis: (3,) unit vector of the long axis.
        principal_moments: (3,) principal moments (sorted ascending).
        hull: The scipy ConvexHull object (or None if hull failed).
    """
    from scipy.spatial import ConvexHull

    if len(positions) < 10:
        return None, None, None, None

    try:
        hull = ConvexHull(positions)
    except Exception:
        return None, None, None, None

    # Fill interior with uniform random points
    bbox_min = positions.min(axis=0)
    bbox_max = positions.max(axis=0)

    n_candidates = n_fill_points * 10
    rng = np.random.RandomState(random_seed)
    candidates = rng.uniform(bbox_min, bbox_max, size=(n_candidates, 3))

    # Efficient inside-hull check using half-space equations
    inside = np.ones(n_candidates, dtype=bool)
    for eq in hull.equations:
        inside &= (candidates @ eq[:3] + eq[3]) <= 1e-10

    interior_points = candidates[inside][:n_fill_points]

    if len(interior_points) < 100:
        # Fallback: use original positions directly
        interior_points = positions

    # Compute centroid and inertia tensor
    centroid = interior_points.mean(axis=0)
    centered = interior_points - centroid

    I = np.zeros((3, 3))
    for i in range(len(centered)):
        r = centered[i]
        r_sq = np.dot(r, r)
        I += r_sq * np.eye(3) - np.outer(r, r)

    # Eigen decomposition — smallest eigenvalue = least rotational inertia = longest axis
    eigenvalues, eigenvectors = np.linalg.eigh(I)
    idx = np.argsort(eigenvalues)
    principal_moments = eigenvalues[idx]
    ap_axis = eigenvectors[:, idx[0]]

    return centroid, ap_axis, principal_moments, hull


def lineage_separation_axis(
    cells_df,
) -> Optional[tuple]:
    """Determine the A-P axis from AB vs P lineage centroid separation.

    This is the biological ground truth: AB is the anterior lineage,
    P is the posterior lineage. The vector from P centroid to AB centroid
    defines the A-P axis direction.

    Args:
        cells_df: DataFrame with columns ['x', 'y', 'z', 'name'].

    Returns:
        (centroid, ap_direction, ap_distance, ab_centroid, p_centroid) or None.
    """
    ab_mask = np.array([is_ab_lineage(n) for n in cells_df["name"]])
    p_mask = np.array([is_p_lineage(n) for n in cells_df["name"]])

    if ab_mask.sum() == 0 or p_mask.sum() == 0:
        return None

    ab_cells = cells_df[ab_mask][["x", "y", "z"]].values
    p_cells = cells_df[p_mask][["x", "y", "z"]].values

    ab_centroid = ab_cells.mean(axis=0)
    p_centroid = p_cells.mean(axis=0)
    global_centroid = np.concatenate([ab_cells, p_cells]).mean(axis=0)

    ap_vector = ab_centroid - p_centroid
    ap_distance = np.linalg.norm(ap_vector)

    if ap_distance < 1e-10:
        return None

    ap_direction = ap_vector / ap_distance

    return global_centroid, ap_direction, ap_distance, ab_centroid, p_centroid


# ==============================================================================
# High-Level API
# ==============================================================================

def compute_ap_axis(
    df,
    method: Literal["pca", "hull_inertia", "lineage"] = "pca",
    time_cutoff: Optional[int] = 255,
    timepoints: Optional[list] = None,
    validate_with_lineage: bool = True,
) -> APAaxisResult:
    """Compute the A-P axis from tracking data.

    Args:
        df: Tracking DataFrame with columns ['t', 'x', 'y', 'z', 'name'].
        method: Computation method:
            - "pca": PCA of nuclear positions (fast, ~1° accuracy).
            - "hull_inertia": Convex hull + inertia tensor (literature method).
            - "lineage": Direct AB vs P centroid separation (biological ground truth).
        time_cutoff: Only use t <= time_cutoff (pre-twitching).
        timepoints: Specific timepoints to average over. If None, uses a range
                    of timepoints with sufficient cell counts.
        validate_with_lineage: If True, uses AB/P lineage to determine which
                               end is anterior (resolving axis sign ambiguity).

    Returns:
        APAaxisResult with the computed axis information.

    Example:
        >>> result = compute_ap_axis(df, method="pca")
        >>> print(result.axis_name)       # "X"
        >>> print(result.anterior_direction)  # -1 (anterior = -X)
    """
    # Apply time cutoff
    work_df = df.copy()
    if time_cutoff is not None:
        work_df = work_df[work_df["t"] <= time_cutoff]

    axis_names = ["X", "Y", "Z"]

    # Select timepoints
    if timepoints is None:
        all_times = sorted(work_df["t"].unique())
        # Use timepoints with at least 10 cells, spread across development
        valid_times = [t for t in all_times
                       if len(work_df[work_df["t"] == t]) >= 10]
        # Evenly sample ~20 timepoints
        if len(valid_times) > 20:
            step = len(valid_times) // 20
            timepoints = valid_times[::step]
        else:
            timepoints = valid_times

    # Collect axis estimates across timepoints
    axis_estimates = []
    centroids = []

    for t in timepoints:
        cells = work_df[work_df["t"] == t]
        if len(cells) < 3:
            continue
        positions = cells[["x", "y", "z"]].values.copy()

        if method == "pca":
            centroid, axis_vec, var_ratio, _ = pca_axis(positions)
            axis_estimates.append(axis_vec)
            centroids.append(centroid)

        elif method == "hull_inertia":
            result = convex_hull_inertia_axis(positions)
            if result[0] is not None:
                centroid, axis_vec, moments, hull = result
                axis_estimates.append(axis_vec)
                centroids.append(centroid)

        elif method == "lineage":
            result = lineage_separation_axis(cells)
            if result is not None:
                centroid, axis_vec, _, _, _ = result
                axis_estimates.append(axis_vec)
                centroids.append(centroid)

        else:
            raise ValueError(f"Unknown method: {method}")

    if not axis_estimates:
        raise ValueError("No valid timepoints found for axis computation")

    # Average axis estimates (handling sign ambiguity)
    axis_estimates = np.array(axis_estimates)
    centroids = np.array(centroids)

    # Align all estimates to the first one
    reference = axis_estimates[0]
    aligned = []
    for v in axis_estimates:
        if np.dot(v, reference) < 0:
            aligned.append(-v)
        else:
            aligned.append(v)
    aligned = np.array(aligned)

    mean_axis = aligned.mean(axis=0)
    mean_axis /= np.linalg.norm(mean_axis)
    mean_centroid = centroids.mean(axis=0)

    # Determine which microscopy axis this corresponds to
    abs_components = np.abs(mean_axis)
    best_idx = int(np.argmax(abs_components))
    axis_name = axis_names[best_idx]

    # Determine anterior direction using lineage
    anterior_direction = 1  # default: anterior = +direction
    explained_variance_ratio = None
    principal_moments = None

    if validate_with_lineage and method != "lineage":
        # Use the lineage method at the timepoint with most AB and P cells
        best_ab_p = None
        best_t = None
        for t in timepoints:
            cells = work_df[work_df["t"] == t]
            n_ab = sum(1 for n in cells["name"] if is_ab_lineage(n))
            n_p = sum(1 for n in cells["name"] if is_p_lineage(n))
            if n_ab > 0 and n_p > 0:
                if best_ab_p is None or (n_ab + n_p) > best_ab_p[0]:
                    best_ab_p = (n_ab + n_p, t, n_ab, n_p)

        if best_ab_p is not None:
            _, t_val, _, _ = best_ab_p
            cells = work_df[work_df["t"] == t_val]

            # Use direct coordinate comparison along the dominant axis
            # to avoid eigenvector sign ambiguity
            ab_mask = np.array([is_ab_lineage(n) for n in cells["name"]])
            p_mask = np.array([is_p_lineage(n) for n in cells["name"]])

            coord_names = ["x", "y", "z"]
            dominant_coord = coord_names[best_idx]
            ab_mean = cells[ab_mask][dominant_coord].mean()
            p_mean = cells[p_mask][dominant_coord].mean()

            # If AB (anterior) is at lower coordinate than P (posterior),
            # then anterior is in the -direction of that axis
            anterior_direction = 1 if ab_mean > p_mean else -1

    # Also compute full PCA at the largest timepoint for var ratio
    if method == "pca":
        largest_t = timepoints[-1]
        cells = work_df[work_df["t"] == largest_t]
        if len(cells) >= 3:
            positions = cells[["x", "y", "z"]].values
            _, _, var_ratio, _ = pca_axis(positions)
            explained_variance_ratio = var_ratio
    elif method == "hull_inertia":
        largest_t = timepoints[-1]
        cells = work_df[work_df["t"] == largest_t]
        if len(cells) >= 10:
            positions = cells[["x", "y", "z"]].values
            result = convex_hull_inertia_axis(positions)
            if result[0] is not None:
                _, _, moments, _ = result
                principal_moments = moments

    return APAaxisResult(
        axis_name=axis_name,
        axis_index=best_idx,
        axis_vector=mean_axis,
        anterior_direction=anterior_direction,
        method=method,
        centroid=mean_centroid,
        explained_variance_ratio=explained_variance_ratio,
        principal_moments=principal_moments,
        time_cutoff=time_cutoff,
        n_cells=int(sum(len(work_df[work_df["t"] == t]) for t in timepoints)),
        n_timepoints=len(timepoints),
        metadata={
            "timepoints_used": timepoints[:10] + ["..."] + timepoints[-3:]
            if len(timepoints) > 15 else timepoints,
        },
    )


# ==============================================================================
# Coordinate Transformations
# ==============================================================================

def get_ap_position(
    df_or_row,
    ap_result: Optional[APAaxisResult] = None,
    centered: bool = False,
) -> np.ndarray:
    """Extract A-P position(s) from tracking data.

    If ap_result is provided, uses the computed axis. Otherwise uses the
    canonical result: A-P axis = X, anterior = -X, so ap_position = -x.

    Args:
        df_or_row: DataFrame with ['x', 'y', 'z'] columns, or a dict/series.
        ap_result: Optional APAaxisResult from compute_ap_axis().
        centered: If True, subtract the centroid.

    Returns:
        A-P position(s). Scalar for a single row, (N,) array for a DataFrame.
    """
    import pandas as pd

    if isinstance(df_or_row, pd.DataFrame):
        xyz = df_or_row[["x", "y", "z"]].values
    elif isinstance(df_or_row, dict):
        xyz = np.array([df_or_row["x"], df_or_row["y"], df_or_row["z"]])
    elif isinstance(df_or_row, pd.Series):
        xyz = np.array([df_or_row["x"], df_or_row["y"], df_or_row["z"]])
    else:
        xyz = np.asarray(df_or_row)

    if ap_result is None:
        # Canonical: A-P = -X (anterior is negative X)
        ap_values = -xyz[..., 0]
    else:
        axis = ap_result.axis_vector
        direction = ap_result.anterior_direction
        centroid = ap_result.centroid if centered else np.zeros(3)

        # Project onto axis, oriented so anterior is positive
        ap_values = (xyz - centroid) @ axis * direction

    return ap_values


def transform_to_ap(
    df,
    ap_result: Optional[APAaxisResult] = None,
    centered: bool = False,
) -> "pd.DataFrame":
    """Transform tracking data so the A-P axis is aligned with coordinate axes.

    Creates new columns:
        - ap: A-P position (anterior = positive)
        - dv: dorsal-ventral position (second principal component)
        - lr: left-right position (third principal component)

    If ap_result is None, uses the canonical result and renames:
        - x → ap (flipped: -x, so anterior is positive)
        - Other axes are not redefined (they may not align with D-V or L-R).

    Args:
        df: Tracking DataFrame.
        ap_result: Optional APAaxisResult.
        centered: If True, center positions on the embryo centroid.

    Returns:
        DataFrame with added 'ap', 'dv', 'lr' columns.
    """
    import pandas as pd

    result = df.copy()

    if ap_result is None:
        # Simple canonical transform
        result["ap"] = -result["x"]  # anterior = -X, flip so anterior is positive
        if centered:
            result["ap"] -= result["ap"].mean()
        result["dv"] = result["y"]  # placeholder
        result["lr"] = result["z"]  # placeholder
        return result

    # Full transform using the computed PCA axes
    xyz = result[["x", "y", "z"]].values
    centroid = ap_result.centroid if centered else np.zeros(3)

    # We have the A-P axis. For D-V and L-R, we can use the other PCA components
    # stored in explained_variance_ratio or metadata, or leave them as-is.
    ap_direction = ap_result.axis_vector * ap_result.anterior_direction

    # Project
    result["ap"] = (xyz - centroid) @ ap_direction
    result["dv"] = result["y"]  # kept as-is unless full basis is available
    result["lr"] = result["z"]

    return result


# ==============================================================================
# Quick Diagnostics
# ==============================================================================

def quick_check(tracks_path: str, time_cutoff: int = 255) -> dict:
    """Quick diagnostic check of a tracking file.

    Returns a dict with the primary axis, anterior direction, and confidence.

    Example:
        >>> quick_check("data/embryo1/tracks.txt")
        {'ap_axis': 'X', 'anterior_direction': 'negative',
         'pca_x_component': 0.9989, 'consistent': True}
    """
    import pandas as pd

    df = pd.read_csv(tracks_path, sep="\t")
    df = df[df["t"] <= time_cutoff]

    # Check t=0 AB/P1
    ab_p1_info = {}
    t0 = df[df["t"] == 0]
    if len(t0) >= 2:
        ab = t0[t0["name"] == "AB"]
        p1 = t0[t0["name"] == "P1"]
        if len(ab) > 0 and len(p1) > 0:
            ab = ab.iloc[0]
            p1 = p1.iloc[0]
            dx, dy, dz = ab["x"] - p1["x"], ab["y"] - p1["y"], ab["z"] - p1["z"]
            ab_p1_info["ab_p1_dx"] = float(dx)
            ab_p1_info["ab_p1_dy"] = float(dy)
            ab_p1_info["ab_p1_dz"] = float(dz)

    # PCA at max timepoint
    max_t = df["t"].max()
    cells = df[df["t"] == max_t]
    positions = cells[["x", "y", "z"]].values
    _, pc1, var_ratio, _ = pca_axis(positions)

    abs_x, abs_y, abs_z = np.abs(pc1)
    axis_names = ["X", "Y", "Z"]
    best_idx = int(np.argmax([abs_x, abs_y, abs_z]))
    best_axis = axis_names[best_idx]

    # Determine anterior direction using raw coordinates (independent of PCA sign)
    # AB lineage cells have lower X than P lineage cells → anterior is -X
    ab_x = cells[[is_ab_lineage(n) for n in cells["name"]]]["x"].mean()
    p_x = cells[[is_p_lineage(n) for n in cells["name"]]]["x"].mean()
    posterior_to_anterior_dx = ab_x - p_x
    # If AB is at lower X than P, then anterior is in the -X direction
    anterior_sign = -1 if posterior_to_anterior_dx < 0 else 1
    anterior_dir = f"{'-' if anterior_sign < 0 else '+'}{best_axis}"

    return {
        "ap_axis": best_axis,
        "anterior_direction": anterior_dir,
        "anterior_sign": anterior_sign,
        "pca_component": {axis_names[i]: float(np.abs(pc1)[i]) for i in range(3)},
        "variance_explained": float(var_ratio[0]),
        "mean_ab_x": float(ab_x),
        "mean_p_x": float(p_x),
        "ab_p1_initial": ab_p1_info,
        "consistent": abs_x > 0.95,  # X is strongly dominant
    }
