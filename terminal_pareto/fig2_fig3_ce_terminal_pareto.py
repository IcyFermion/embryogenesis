"""Generate C. elegans terminal-cell Pareto Figures 2 and 3."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from publication import style as ps
from publication.artists import proportional_limits  # noqa: F401  (re-exported)
from publication.figures.terminal_global import (  # noqa: F401  (re-exported drawing API)
    CELL_STATE_COLOR, FIRST_COUSIN_COLOR, FULL_RANDOM_COLOR, SECOND_COUSIN_COLOR, SUPPORT_AXES_BOTTOM,
    THIRD_COUSIN_COLOR, TRAVEL_COLOR, TREE_DISTANCE_COLOR, plot_main, plot_support_b, plot_support_c,
)
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
)
from terminal_pareto.front_coordinates import (
    build_endpoint_transform,
    null_sd_coordinates,
    percent_natural_coordinates,
)
from terminal_pareto.global_analysis import get_or_compute_global_analysis


OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "publication"
DIAGNOSTIC_OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "ce_protein"
EDGE_RETENTION_CMAP = ps.EDGE_RETENTION_CMAP


def _standardize(x, y, reference):
    return (
        (x - reference["lineage_xyz"]) / reference["xyz_std"],
        (y - reference["lineage_exp"]) / reference["exp_std"],
    )


def _load_analysis():
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_map, valid = dl.load_elegans_tracking(dl.T_CE)
    protein = dl.load_protein_expression()
    selected = dl.load_prot_sel()
    valid = [name for name in valid if name in protein.index]
    terminal_nodes, terminal_parents = dl.collect_terminals(lineage, valid)
    xyz_mat, exp_mat, _ = pe.build_cost_matrices(
        terminal_nodes, terminal_parents, xyz_map, protein, selected
    )
    first_groups = pe.build_ancestor_groups(terminal_nodes, tree_index, 2)
    first_stats = pe.compute_cousin_random_stats(
        xyz_mat, exp_mat, first_groups, n_random=1000, seed=42
    )
    twr = lm.combined_lineage_proximity(
        xyz_mat, exp_mat, terminal_parents, terminal_nodes, tree_index,
        first_stats, iteration=300, n_random=100,
    )

    nulls = {}
    for degree, steps, seed in [(1, 2, 42), (2, 3, 43), (3, 4, 44)]:
        groups = pe.build_ancestor_groups(terminal_nodes, tree_index, steps)
        raw_x, raw_y = pe.compute_group_shuffle_costs(
            xyz_mat, exp_mat, groups, n_random=1000, seed=seed
        )
        nulls[degree] = _standardize(raw_x, raw_y, first_stats)
    rows = np.arange(len(terminal_nodes))
    rng = np.random.default_rng(45)
    full_x = np.empty(2000)
    full_y = np.empty(2000)
    for i in range(2000):
        perm = rng.permutation(len(terminal_nodes))
        full_x[i] = xyz_mat[rows, perm].sum()
        full_y[i] = exp_mat[rows, perm].sum()
    nulls["full"] = _standardize(full_x, full_y, first_stats)
    return twr, nulls


def _profile_analysis(context, *, force=False):
    """Load one cache-keyed numerical result and prepare its null-SD view."""
    result = get_or_compute_global_analysis(context, force=force)
    natural = result.natural_costs
    stds = result.null_stds
    key_map = {
        1: "first_cousin",
        2: "second_cousin",
        3: "third_cousin",
        "full": "full_random",
    }
    nulls = {
        key: null_sd_coordinates(
            *result.null_raw[cache_key], natural_costs=natural,
            null_stds=stds)
        for key, cache_key in key_map.items()
    }
    return result.twr, nulls, result


def _display_data(twr, nulls, *, mode, result=None, context=None):
    """Return plot coordinates without modifying saved numerical results."""
    reference = context.cache_key if context is not None else "legacy"
    if mode == "null_sd":
        return dict(
            x=np.asarray(twr["xyz_arr"], dtype=float),
            y=np.asarray(twr["exp_arr"], dtype=float),
            nulls=nulls,
            natural=(0.0, 0.0),
            xlabel=("Travel distance\n(null standard deviations; "
                    "natural lineage = 0)"),
            ylabel=("Cell-state distance\n(null standard deviations; "
                    "natural lineage = 0)"),
            metadata={"display_mode": "null_sd",
                      "reference_analysis_id": reference},
        )
    if result is None or context is None:
        raise ValueError(f"Display mode {mode!r} requires a saved profile result")
    key_map = {
        1: "first_cousin", 2: "second_cousin",
        3: "third_cousin", "full": "full_random",
    }
    if mode == "endpoint":
        transform = build_endpoint_transform(
            result.raw_front_travel, result.raw_front_state,
            travel_optimum_index=context.spec.sweep_intervals,
            state_optimum_index=0,
            reference_analysis_id=result.analysis_cache_key,
            assignment_ids=[f"sweep:{index}" for index in range(
                context.spec.sweep_intervals + 1)],
        )
        x, y = transform.transform(
            result.raw_front_travel, result.raw_front_state)
        natural_x, natural_y = transform.transform(
            np.asarray([result.natural_costs[0]]),
            np.asarray([result.natural_costs[1]]))
        display_nulls = {
            key: transform.transform(*result.null_raw[raw_key])
            for key, raw_key in key_map.items()
        }
        return dict(
            x=x, y=y, nulls=display_nulls,
            natural=(float(natural_x[0]), float(natural_y[0])),
            xlabel="Travel distance\n(fraction of endpoint cost span)",
            ylabel="Cell-state distance\n(fraction of endpoint cost span)",
            metadata=transform.metadata(clipping=False),
        )
    if mode == "percent_natural":
        x, y = percent_natural_coordinates(
            result.raw_front_travel, result.raw_front_state,
            natural_costs=result.natural_costs)
        display_nulls = {
            key: percent_natural_coordinates(
                *result.null_raw[raw_key], natural_costs=result.natural_costs)
            for key, raw_key in key_map.items()
        }
        return dict(
            x=x, y=y, nulls=display_nulls, natural=(0.0, 0.0),
            xlabel="Travel-distance change from natural lineage (%)",
            ylabel="Cell-state-distance change from natural lineage (%)",
            metadata={"display_mode": "percent_natural",
                      "reference_analysis_id": result.analysis_cache_key,
                      "natural_costs": list(result.natural_costs)},
        )
    raise ValueError(f"Unsupported display mode: {mode}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the terminal-cell publication figure panels."
    )
    parser.add_argument(
        "--tree-distance-mode",
        choices=("changed_edges", "all_edges", "both"),
        default="changed_edges",
        help=("Tree-distance definition(s) for supporting panel B "
              "(default: changed_edges; all_edges is diagnostic-only)."),
    )
    parser.add_argument(
        "--profile",
        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"),
        help="Opt into an isolated profile-aware run (omission keeps legacy paths).",
    )
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--display", choices=("null_sd", "endpoint",
                                               "percent_natural"),
                        default="null_sd")
    parser.add_argument(
        "--main-panel-letter",
        help="Optional panel letter for the Figure 2 front component.",
    )
    parser.add_argument("--force-analysis", action="store_true")
    args = parser.parse_args(argv)
    ps.configure()
    display = None
    out_dir = OUT
    diagnostic_out = DIAGNOSTIC_OUT
    if args.profile is None:
        twr, nulls = _load_analysis()
    else:
        context = build_analysis_context(
            args.profile, run_id=args.run_id, output_root=args.output_root)
        context.write()
        twr, nulls, result = _profile_analysis(
            context, force=args.force_analysis)
        display = _display_data(
            twr, nulls, mode=args.display, result=result, context=context)
        out_dir = context.run_paths.display(args.display)
        diagnostic_out = context.run_paths.analysis / "diagnostics" / args.display
        out_dir.mkdir(parents=True, exist_ok=True)
        diagnostic_out.mkdir(parents=True, exist_ok=True)
        metadata = dict(display["metadata"])
        metadata.update({
            "profile": context.profile,
            "context_cache_key": context.cache_key,
            "analysis_cache_key": result.analysis_cache_key,
            "assignment_ids_hash": hashlib.sha256(
                result.assignments.tobytes()).hexdigest(),
        })
        (out_dir / "fig2_fig3_display_manifest.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    plot_main(twr, nulls, display=display, out_dir=out_dir,
              panel_letter=args.main_panel_letter)
    modes = ("changed_edges", "all_edges") if args.tree_distance_mode == "both" else (args.tree_distance_mode,)
    for mode in modes:
        plot_support_b(twr, distance_mode=mode, display=display,
                       out_dir=out_dir, diagnostic_out=diagnostic_out)
    plot_support_c(twr, display=display, out_dir=out_dir)
    print("Publication figure panels written to", out_dir)


if __name__ == "__main__":
    main()
