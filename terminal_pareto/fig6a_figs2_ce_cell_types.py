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
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from publication import style as ps
from publication.figures import terminal_cell_types as cell_type_figures
from publication.figures.terminal_cell_types import (  # noqa: F401  (re-exported)
    HEATMAP_BINS, HEATMAP_SMOOTH_BANDWIDTH, KEYPOINT_COLORS, KEYPOINT_LABELS, KEYPOINT_MARKERS, ORDER,
    display_type_name, marker_area, ordered_types, smooth_retention_profiles,
)
from terminal_pareto.subtree_analysis import build_type_map
from terminal_pareto.subtree_explore import (
    SMALL_TYPES,
    _merge_small,
    decompose_global_front_by_type,
    decompose_global_front_retention_by_type,
)
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    AnalysisContext,
    build_analysis_context,
)

OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"
FULL_RETENTION_CACHE = ANALYSIS_OUT / "global_front_retention_by_cell_type.csv"
KEYPOINT_CACHE = ANALYSIS_OUT / "global_front_by_cell_type.csv"
CELL_TYPE_CACHE_MANIFEST = "fig6a_cell_type_cache_manifest.json"
CELL_TYPE_CACHE_VERSION = "fig6a-cell-type-cache-v1"

def _file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cell_type_cache_identity(iteration, context):
    """Configuration fields that must match before cached tables are read."""
    return {
        "version": CELL_TYPE_CACHE_VERSION,
        "profile": context.profile if context is not None else "embryo1_legacy",
        "context_cache_key": context.cache_key if context is not None else None,
        "iteration": int(iteration),
        "terminal_edges": (len(context.terminal_nodes)
                           if context is not None else 299),
    }


def write_cell_type_cache_manifest(analysis_out, iteration, context):
    """Record cache identity and hashes only after both tables are complete."""
    analysis_out = Path(analysis_out)
    full_cache = analysis_out / FULL_RETENTION_CACHE.name
    keypoint_cache = analysis_out / KEYPOINT_CACHE.name
    payload = cell_type_cache_identity(iteration, context)
    payload["files"] = {
        full_cache.name: _file_sha256(full_cache),
        keypoint_cache.name: _file_sha256(keypoint_cache),
    }
    path = analysis_out / CELL_TYPE_CACHE_MANIFEST
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def _validate_cell_type_tables(full_df, keypoint_df, iteration,
                               expected_edges):
    required_full = {
        "front_index", "sweep_index", "alpha", "type", "n", "er",
    }
    required_keypoints = {"keypoint", "type", "n", "er"}
    if not required_full.issubset(full_df.columns):
        raise ValueError("Full-front cell-type cache has incomplete columns")
    if not required_keypoints.issubset(keypoint_df.columns):
        raise ValueError("Keypoint cell-type cache has incomplete columns")
    expected_alpha = full_df["sweep_index"].to_numpy(dtype=float) / iteration
    if not np.allclose(full_df["alpha"], expected_alpha,
                       atol=1e-12, rtol=0):
        raise ValueError("Cell-type cache sweep resolution does not match request")
    if (full_df["sweep_index"].min() < 0
            or full_df["sweep_index"].max() > iteration):
        raise ValueError("Cell-type cache contains out-of-range sweep indices")
    front_counts = full_df.groupby("front_index")["n"].sum().to_numpy()
    keypoint_counts = keypoint_df.groupby("keypoint")["n"].sum().to_numpy()
    if (not np.all(front_counts == expected_edges)
            or not np.all(keypoint_counts == expected_edges)):
        raise ValueError("Cell-type cache cohort does not match requested profile")


def load_validated_cell_type_caches(analysis_out, iteration, context):
    """Validate identity and file hashes before loading reusable CSV tables."""
    analysis_out = Path(analysis_out)
    full_cache = analysis_out / FULL_RETENTION_CACHE.name
    keypoint_cache = analysis_out / KEYPOINT_CACHE.name
    missing = [path.name for path in (full_cache, keypoint_cache)
               if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing cell-type cache files: {missing}")

    expected = cell_type_cache_identity(iteration, context)
    manifest_path = analysis_out / CELL_TYPE_CACHE_MANIFEST
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        mismatched = {
            key: (manifest.get(key), value)
            for key, value in expected.items()
            if manifest.get(key) != value
        }
        if mismatched:
            raise ValueError(
                f"Cell-type cache configuration mismatch: {mismatched}")
        recorded_hashes = manifest.get("files", {})
        for path in (full_cache, keypoint_cache):
            if recorded_hashes.get(path.name) != _file_sha256(path):
                raise ValueError(f"Cell-type cache hash mismatch: {path.name}")
    elif context is not None:
        raise ValueError(
            "Profile-specific cell-type caches lack an identity manifest; "
            "rerun without --reuse-cache")

    # Loading happens only after configuration and byte-level provenance pass.
    full_df = pd.read_csv(full_cache)
    keypoint_df = pd.read_csv(keypoint_cache)
    _validate_cell_type_tables(
        full_df, keypoint_df, iteration, expected["terminal_edges"])
    return full_df, keypoint_df


def compute_cell_type_caches(iteration=300, *,
                             context: AnalysisContext | None = None,
                             analysis_out=ANALYSIS_OUT):
    """Evaluate full-front retention and the three-keypoint decomposition."""
    if context is None:
        lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
        tree_index = lm.build_lineage_tree_index(lineage)
        xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
        protein_exp = dl.load_protein_expression()
        prot_sel = dl.load_prot_sel()
        v_prot = [name for name in valid_ce if name in protein_exp.index]
        gp_map = pe.build_grandparent_map(lineage)
        tn, tp = dl.collect_terminals(lineage, v_prot)
        prepared = None
    else:
        lineage = context.lineage
        tree_index = context.tree_index
        xyz_ce = None
        protein_exp = context.protein_exp
        prot_sel = context.prot_sel
        v_prot = context.v_prot
        gp_map = context.gp_map
        tn, tp = context.terminal_nodes, context.terminal_parents
        prepared = context.cost_matrices()
    type_map = build_type_map(tn)
    full = decompose_global_front_retention_by_type(
        lineage, v_prot, tn, tp, xyz_ce, protein_exp, prot_sel, type_map,
        tree_index, gp_map, iteration=iteration,
        prepared_matrices=prepared)
    keypoints = decompose_global_front_by_type(
        lineage, v_prot, tn, tp, xyz_ce, protein_exp, prot_sel, type_map,
        tree_index, gp_map, iteration=iteration,
        prepared_matrices=prepared)
    analysis_out = Path(analysis_out)
    analysis_out.mkdir(parents=True, exist_ok=True)
    full.to_csv(analysis_out / FULL_RETENTION_CACHE.name, index=False)
    keypoints.to_csv(analysis_out / KEYPOINT_CACHE.name, index=False)
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


def plot_retention_heatmap(full_df, out_dir=OUT):
    """Figure 6A from the raw cached table (merges n<=4 groups first)."""
    return cell_type_figures.plot_retention_heatmap(merge_small_retention(full_df), out_dir=out_dir)


def plot_endpoint_ledger(keypoint_df, out_dir=OUT):
    """Figure S2 from the raw cached keypoint table (merges n<=4 groups first)."""
    return cell_type_figures.plot_endpoint_ledger(_merge_small(keypoint_df), out_dir=out_dir)


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
    parser.add_argument(
        "--profile",
        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"),
        help="Opt into an isolated profile-aware run (omission keeps legacy paths).",
    )
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    ps.configure()

    context = None
    analysis_out = ANALYSIS_OUT
    out = args.out or OUT
    if args.profile is not None:
        context = build_analysis_context(
            args.profile, run_id=args.run_id, output_root=args.output_root,
            sweep_intervals=args.iteration)
        analysis_out = context.run_paths.analysis
        if args.out is None:
            out = context.run_paths.display("endpoint")
    out.mkdir(parents=True, exist_ok=True)
    full_cache = analysis_out / FULL_RETENTION_CACHE.name
    keypoint_cache = analysis_out / KEYPOINT_CACHE.name
    if args.reuse_cache:
        # Validate the existing run before context.write() can replace its
        # provenance. A mismatch raises without touching manifests or tables.
        full_df, keypoint_df = load_validated_cell_type_caches(
            analysis_out, args.iteration, context)
        if context is not None:
            context.write()
        print(f"Reused {len(full_df)} full-front and "
              f"{len(keypoint_df)} keypoint cell-type rows.")
    else:
        if context is not None:
            context.write()
        full_df, keypoint_df = compute_cell_type_caches(
            iteration=args.iteration, context=context,
            analysis_out=analysis_out)
        manifest = write_cell_type_cache_manifest(
            analysis_out, args.iteration, context)
        print(f"Wrote {len(full_df)} full-front rows to "
              f"{full_cache} and {len(keypoint_df)} keypoint rows "
              f"to {keypoint_cache}; manifest {manifest}.")

    # These panels retain their canonical-position and interpretable ledger
    # axes; endpoint display normalization is intentionally not applied.
    plot_retention_heatmap(full_df, out_dir=out)
    plot_endpoint_ledger(keypoint_df, out_dir=out)
    print("Figure 6A and Figure S2 panels written to", out)


if __name__ == "__main__":
    main()
