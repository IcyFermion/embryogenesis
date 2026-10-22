"""Shared test helpers."""

from __future__ import annotations

from contextlib import ExitStack, contextmanager
import sys
from unittest import mock

import numpy as np
import scipy.optimize


@contextmanager
def no_experiments():
    """Make any solver or null/reference generation call fail loudly."""
    def forbidden(*args, **kwargs):
        raise AssertionError("Presentation build reached scientific computation")

    from full_tree_pareto import cousin_references, cross_species_analysis as fcs, pooled_analysis
    from terminal_pareto import cross_species_analysis as tcs
    from terminal_pareto import global_analysis, fig6bc_ce_within_type as within, fig5_table1_ce_canonical_metrics as canonical
    from terminal_pareto import fig6a_figs2_ce_cell_types as cell_types
    from terminal_pareto import feature_sensitivity, representation_criteria
    from terminal_pareto import apoptosis_sensitivity, apoptosis_single_deaths, apoptosis_subtree_control
    from terminal_pareto import subtree_symmetry
    from terminal_pareto import subtree_symmetry_robustness
    targets = [(tcs, "analyze"), (tcs, "write_analysis"), (tcs, "endpoint_assignment"),
               (feature_sensitivity, "analyze"), (feature_sensitivity, "draw_feature_sets"),
               (feature_sensitivity, "run_survey"),
               (representation_criteria, "compute"), (representation_criteria, "criteria"),
               (representation_criteria, "front"), (representation_criteria, "draw_random_sets"),
               (subtree_symmetry, "front_metrics"), (subtree_symmetry, "compute_compare"),
               (subtree_symmetry_robustness, "compute"),
               (apoptosis_sensitivity, "compute"), (apoptosis_sensitivity, "controls"),
               (apoptosis_subtree_control, "compute"), (apoptosis_subtree_control, "controls"),
               (apoptosis_single_deaths, "compute"), (apoptosis_single_deaths, "single"),
               (fcs, "build"), (fcs, "solve_sweep"), (fcs, "random_rebuilds"),
               (cousin_references, "sample_permutations"), (pooled_analysis, "_solve"),
               (pooled_analysis, "random_nulls"), (pooled_analysis, "gaussian"),
               (global_analysis, "compute_global_analysis"), (global_analysis, "get_or_compute_global_analysis"),
               (within, "analyze"), (canonical, "analyze_all"), (cell_types, "compute_cell_type_caches"),
               (scipy.optimize, "linear_sum_assignment"), (np.random, "default_rng")]
    targets += [(module, "linear_sum_assignment") for module in list(sys.modules.values())
                if getattr(module, "linear_sum_assignment", None) is scipy.optimize.linear_sum_assignment]
    with ExitStack() as stack:
        for module, name in targets:
            stack.enter_context(mock.patch.object(module, name, forbidden))
        yield
