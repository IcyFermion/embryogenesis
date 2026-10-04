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
    targets = [(tcs, "analyze"), (tcs, "write_analysis"), (tcs, "endpoint_assignment"),
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
