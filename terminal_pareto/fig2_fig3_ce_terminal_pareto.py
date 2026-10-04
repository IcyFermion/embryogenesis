"""Global terminal-cell Pareto analysis behind Figures 2 and 3 (numerical only).

Figures are drawn by ``python -m publication build --family terminal-primary``.
``_load_analysis`` keeps the historical no-profile reconstruction used by the
migration validator.
"""

import argparse
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
)
from terminal_pareto.front_coordinates import null_sd_coordinates
from terminal_pareto.global_analysis import get_or_compute_global_analysis




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


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute or validate the cached global terminal front.")
    parser.add_argument("--profile", required=True,
                        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"))
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--force-analysis", action="store_true")
    args = parser.parse_args(argv)
    context = build_analysis_context(args.profile, run_id=args.run_id, output_root=args.output_root)
    context.write()
    _, _, result = _profile_analysis(context, force=args.force_analysis)
    print(f"Global front cache {result.analysis_cache_key} in {context.run_paths.analysis}")


if __name__ == "__main__":
    main()
