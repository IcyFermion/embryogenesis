"""Primary and supplementary canonical subtree summaries.

This presentation-only renderer uses the validated per-subtree canonical
metrics produced by ``fig5_table1_ce_canonical_metrics.py``.

Profile-aware Figure 5 uses named rows: a paired natural/null comparison for
five major subtrees, then canonical position and natural distance for every
surveyed subtree, grouped by branch. Shapes do not encode lineage identity.

Historical and supplementary panels:
  A. Symbolic definition of canonical position u and front distances.
  B. Lineage- and first-cousin-null distances for the five major subtrees.
  C. All eligible lineage-to-front distances, with the five major-subtree
     first-cousin-null distances overlaid as references.

Drawing code lives in ``publication/figures/terminal_canonical.py``; this
module keeps the validated loader and command-line entrypoint.

The primary figure uses null-independent endpoint-normalized distance d_LP.
The supplementary version preserves first-cousin-null-relative r and adds a
correlation inset. No analysis is recomputed and no maximum-retention
position is assigned a second canonical coordinate.
"""

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from publication import style as ps
from publication.figures.terminal_canonical import (  # noqa: F401  (re-exported drawing API)
    LP_COLOR, MAJOR_MARKERS, MAJOR_SUBTREES, NP_COLOR, REGION_LABELS, REGION_MARKERS, REGION_ORDER,
    STATE_COLOR, TRAVEL_COLOR, U_COLOR, add_colored_fraction, add_metric_correlation_inset,
    add_panel_letter, annotate_abplpp_pair, plot_all_dlp_with_major_dnp, plot_definition_dlp,
    plot_definition_r, plot_major_dlp_dnp, plot_primary_rows, plot_summary, save_component_crops,
    scatter_records,
)
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
)
from terminal_pareto.fig5_table1_ce_canonical_metrics import (
    ITERATION,
    MIN_CELLS,
    load_validated_canonical_metrics,
)


OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "publication"
DEFAULT_METRICS = (Path(__file__).resolve().parent / "output"
                   / "ce_subtree_canonical_metrics.csv")


def load_metrics(path, min_cells=MIN_CELLS, iteration=ITERATION, context=None):
    """Load the established canonical metrics without recomputing analysis."""
    data = load_validated_canonical_metrics(
        path, min_cells=min_cells, iteration=iteration, context=context)
    return data[data["endpoint_ok"].astype(bool)
                & data["u_monotone"].astype(bool)].copy()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--min-cells", type=int, default=MIN_CELLS)
    parser.add_argument("--iteration", type=int, default=ITERATION)
    parser.add_argument("--primary-only", action="store_true",
                        help="Render Figure 5 without rebuilding Figure S1.")
    parser.add_argument(
        "--profile",
        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"),
        help="Read validated metrics from an isolated profile run.",
    )
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args(argv)
    ps.configure()
    context = None
    if args.profile is not None:
        context = build_analysis_context(
            args.profile, run_id=args.run_id, output_root=args.output_root,
            sweep_intervals=args.iteration)
        if args.metrics == DEFAULT_METRICS:
            args.metrics = (context.run_paths.analysis
                            / "ce_subtree_canonical_metrics.csv")
        if args.out == OUT:
            args.out = context.run_paths.display("endpoint")
    # Rendering validates the existing cache and never rewrites its provenance.
    metrics = load_metrics(
        args.metrics, min_cells=args.min_cells,
        iteration=args.iteration, context=context)
    if args.profile is None:
        primary_names = (
            "fig5A_ce_canonical_definition",
            "fig5B_ce_canonical_major_subtrees",
            "fig5C_ce_canonical_all_subtrees",
        )
        primary_letters = ("A", "B", "C")
        plot_summary(metrics, metric="dlp", component_names=primary_names,
                     out_dir=args.out, panel_letters=primary_letters)
    else:
        plot_primary_rows(metrics, args.out)
    if not args.primary_only:
        plot_summary(metrics, metric="r",
                     component_names=(
                         "figS1A_ce_cousin_r_definition",
                         "figS1B_ce_cousin_r_major_subtrees",
                         "figS1C_ce_cousin_r_all_subtrees",
                     ),
                     out_dir=args.out)
    print(f"Wrote canonical summaries for {len(metrics)} subtrees to {args.out}")


if __name__ == "__main__":
    main()
