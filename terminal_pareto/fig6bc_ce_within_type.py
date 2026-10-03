"""Figure 6 panels B--C: optimization within terminal cell types.

This analysis asks a different question from ``fig6a_figs2_ce_cell_types.py``. Instead
of decomposing assignments from the unrestricted embryo-wide Pareto front,
it forbids assignments across terminal cell-type groups and solves the
assignment problem independently inside each group.

Two outputs are produced:

1. Individual within-type fronts in shared, embryo-level null-SD units per
   cell. This avoids pretending that degenerate or very sparse within-type
   cousin nulls define comparable scales.
2. The aggregate type-preserving front versus the unrestricted global front.
   Because the restricted feasible set is a subset of the global one, the
   restricted single-objective endpoints must never outperform the global
   endpoints; this is asserted as an implementation check.

The established first-cousin-null proximity is also computed per type when
both objective variances are positive and stored in the summary table. A
full-random-within-type sensitivity value is stored separately, never used as
a silent replacement for a degenerate cousin null.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from publication import style as ps
from terminal_pareto.subtree_analysis import (
    build_type_map,
    exact_cousin_stats,
    first_cousin_null_summary,
)
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    AnalysisContext,
    build_analysis_context,
)
from terminal_pareto.front_coordinates import (
    DegenerateEndpointSpan,
    EndpointTransform,
)


OUT = Path(__file__).resolve().parent / "output" / "legacy" / "rebuild" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"
ITERATION = 300
MIN_DISPLAY_N = 12
from publication.figures.terminal_within_type import (  # noqa: E402,F401  (re-exported drawing API)
    CELL_STATE_COLOR, TRAVEL_COLOR, TYPE_COLORS, TYPE_ORDER, TYPE_RESTRICTED_COLOR, UNRESTRICTED_COLOR,
    plot_within_type, save_panel_crops, unique_front,
)


def load_primary_data(context: AnalysisContext | None = None):
    if context is not None:
        tn, tp = context.terminal_nodes, context.terminal_parents
        return dict(
            lineage=context.lineage, tree_index=context.tree_index,
            xyz=None, expression=context.protein_exp,
            features=context.prot_sel, tn=tn, tp=tp,
            type_map=build_type_map(tn), gp_map=context.gp_map,
            prepared_matrices=context.cost_matrices(),
        )
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    protein_exp = dl.load_protein_expression()
    prot_sel = dl.load_prot_sel()
    v_prot = [name for name in valid_ce if name in protein_exp.index]
    tn, tp = dl.collect_terminals(lineage, v_prot)
    return dict(
        lineage=lineage, tree_index=tree_index, xyz=xyz_ce,
        expression=protein_exp, features=prot_sel, tn=tn, tp=tp,
        type_map=build_type_map(tn), gp_map=pe.build_grandparent_map(lineage),
        prepared_matrices=None,
    )


def exact_null_stats(xm, em, groups):
    """Exact permutation-null moments in the random-stats API shape."""
    mx, me, vx, ve, _ = exact_cousin_stats(xm, em, groups)
    return dict(
        xyz_mean=float(mx), exp_mean=float(me),
        xyz_std=float(np.sqrt(max(vx, 0.0))),
        exp_std=float(np.sqrt(max(ve, 0.0))),
        lineage_xyz=float(np.diag(xm).sum()),
        lineage_exp=float(np.diag(em).sum()),
    )


def relative_front_record(xm, em, tp, stats, iteration):
    """Natural-to-front ratio using one explicitly supplied null model."""
    if not (stats["xyz_std"] > 1e-10 and stats["exp_std"] > 1e-10):
        return dict(relative_distance=np.nan, closest_sigma=np.nan,
                    null_to_front_sigma=np.nan, natural_on_front=False)
    xa, ea, _edge, _kp = pe.compute_std_scaled_pareto(
        xm, em, tp, stats, iteration=iteration)
    distance = np.hypot(xa, ea)
    closest = int(np.argmin(distance))
    nx = (stats["xyz_mean"] - stats["lineage_xyz"]) / stats["xyz_std"]
    ny = (stats["exp_mean"] - stats["lineage_exp"]) / stats["exp_std"]
    lp = float(distance[closest])
    np_dist = float(np.hypot(nx - xa[closest], ny - ea[closest]))
    return dict(
        relative_distance=lp / np_dist if np_dist > 1e-10 else np.nan,
        closest_sigma=lp,
        null_to_front_sigma=np_dist,
        natural_on_front=bool(np.any((np.abs(xa) < 1e-9)
                                     & (np.abs(ea) < 1e-9))),
    )


def solve_block_sweep(xs, es, groups, lineage_x, lineage_e, iteration):
    """Pareto sweep with assignments constrained to cell-type blocks."""
    travel = np.empty(iteration + 1)
    state = np.empty(iteration + 1)
    for step in range(iteration + 1):
        alpha = step / iteration
        total_x = total_e = 0.0
        for indices in groups.values():
            ix = np.asarray(indices, dtype=int)
            xb = xs[np.ix_(ix, ix)]
            eb = es[np.ix_(ix, ix)]
            ri, ci = linear_sum_assignment(alpha * xb + (1 - alpha) * eb)
            total_x += float(xb[ri, ci].sum())
            total_e += float(eb[ri, ci].sum())
        travel[step] = total_x - lineage_x
        state[step] = total_e - lineage_e
    return travel, state


def analyze(iteration=ITERATION, *, context: AnalysisContext | None = None):
    data = load_primary_data(context)
    tn, tp = data["tn"], data["tp"]
    if data["prepared_matrices"] is None:
        xm, em, _ = pe.build_cost_matrices(
            tn, tp, data["xyz"], data["expression"], data["features"])
    else:
        xm, em = data["prepared_matrices"]
    global_null = first_cousin_null_summary(
        xm, em, tn, data["gp_map"], seed=42)["_raw"]
    xstd, estd = global_null["xyz_std"], global_null["exp_std"]
    xs, es = xm / xstd, em / estd
    lineage_x = float(np.diag(xs).sum())
    lineage_e = float(np.diag(es).sum())

    global_x, global_e, _edge, _kp = pe.compute_std_scaled_pareto(
        xm, em, tp, global_null, iteration=iteration)
    types = [data["type_map"][cell] for cell in tn]
    indices_by_type = {
        cell_type: [i for i, value in enumerate(types) if value == cell_type]
        for cell_type in sorted(set(types))
    }
    restricted_x, restricted_e = solve_block_sweep(
        xs, es, indices_by_type, lineage_x, lineage_e, iteration)

    # A type-preserving feasible set is a subset of the unrestricted set.
    tol = 1e-8
    if restricted_x[-1] < global_x[-1] - tol:
        raise AssertionError("Type-restricted travel optimum beats global optimum")
    if restricted_e[0] < global_e[0] - tol:
        raise AssertionError("Type-restricted state optimum beats global optimum")

    front_rows = []
    summary_rows = []
    for cell_type, indices in indices_by_type.items():
        ix = np.asarray(indices, dtype=int)
        n = len(ix)
        xb = xm[np.ix_(ix, ix)]
        eb = em[np.ix_(ix, ix)]
        xsb = xs[np.ix_(ix, ix)]
        esb = es[np.ix_(ix, ix)]
        lx = float(np.diag(xsb).sum())
        le = float(np.diag(esb).sum())
        local_x, local_e = solve_block_sweep(
            xsb, esb, {cell_type: list(range(n))}, lx, le, iteration)
        local_x /= n
        local_e /= n
        reference_id = (f"{context.cache_key}:{cell_type}"
                        if context is not None else f"legacy:{cell_type}")
        try:
            local_transform = EndpointTransform.from_endpoints(
                reference_analysis_id=reference_id,
                travel_optimum_assignment_id=f"{cell_type}:sweep:{iteration}",
                state_optimum_assignment_id=f"{cell_type}:sweep:0",
                travel_optimum_costs=(local_x[-1], local_e[-1]),
                state_optimum_costs=(local_x[0], local_e[0]),
            )
            local_endpoint_x, local_endpoint_y = local_transform.transform(
                local_x, local_e)
            local_natural_x, local_natural_y = local_transform.transform(
                np.asarray([0.0]), np.asarray([0.0]))
            endpoint_valid = True
            display_reference = reference_id
        except DegenerateEndpointSpan:
            # Preserve the analysis and report an explicit null-SD fallback;
            # never divide by zero or silently invent another endpoint span.
            local_endpoint_x, local_endpoint_y = local_x.copy(), local_e.copy()
            local_natural_x = local_natural_y = np.asarray([0.0])
            endpoint_valid = False
            display_reference = f"{reference_id}:null_sd_fallback"
        for step, (dx, de, endpoint_x, endpoint_y) in enumerate(zip(
                local_x, local_e, local_endpoint_x, local_endpoint_y)):
            front_rows.append(dict(
                type=cell_type, n=n, step=step, alpha=step / iteration,
                travel_sigma_per_cell=dx, state_sigma_per_cell=de,
                endpoint_travel=float(endpoint_x),
                endpoint_state=float(endpoint_y),
                endpoint_natural_travel=float(local_natural_x[0]),
                endpoint_natural_state=float(local_natural_y[0]),
                endpoint_reference=display_reference,
                endpoint_valid=endpoint_valid,
            ))

        tnodes = [tn[i] for i in ix]
        tparents = [tp[i] for i in ix]
        cousin_groups = pe.build_cousin_groups(tnodes, data["gp_map"])
        cousin_stats = exact_null_stats(xb, eb, cousin_groups)
        cousin_rec = relative_front_record(
            xb, eb, tparents, cousin_stats, iteration)
        random_stats = exact_null_stats(
            xb, eb, [list(range(n))] if n >= 2 else [])
        random_rec = relative_front_record(
            xb, eb, tparents, random_stats, iteration)
        closest = float(np.min(np.hypot(local_x, local_e)))
        summary_rows.append(dict(
            type=cell_type,
            n=n,
            cousin_groups=len(cousin_groups),
            cousin_covered=sum(map(len, cousin_groups)),
            cousin_covered_frac=sum(map(len, cousin_groups)) / n,
            cousin_xyz_std=cousin_stats["xyz_std"],
            cousin_state_std=cousin_stats["exp_std"],
            cousin_null_valid=bool(cousin_stats["xyz_std"] > 1e-10
                                   and cousin_stats["exp_std"] > 1e-10),
            cousin_relative_distance=cousin_rec["relative_distance"],
            cousin_natural_on_front=cousin_rec["natural_on_front"],
            full_random_relative_distance=random_rec["relative_distance"],
            endpoint_display_valid=endpoint_valid,
            endpoint_display_reference=display_reference,
            shared_scale_closest_per_cell=closest,
            travel_reduction_per_cell=-float(local_x[-1]),
            state_reduction_per_cell=-float(local_e[0]),
        ))

    aggregate = pd.DataFrame({
        "step": np.arange(iteration + 1),
        "alpha": np.arange(iteration + 1) / iteration,
        "global_travel_sigma": global_x,
        "global_state_sigma": global_e,
        "restricted_travel_sigma": restricted_x,
        "restricted_state_sigma": restricted_e,
    })
    shared_transform = EndpointTransform.from_endpoints(
        reference_analysis_id=(context.cache_key if context is not None
                               else "legacy_global_unrestricted"),
        travel_optimum_assignment_id=f"unrestricted:sweep:{iteration}",
        state_optimum_assignment_id="unrestricted:sweep:0",
        travel_optimum_costs=(global_x[-1], global_e[-1]),
        state_optimum_costs=(global_x[0], global_e[0]),
    )
    aggregate["global_endpoint_travel"], aggregate["global_endpoint_state"] = (
        shared_transform.transform(global_x, global_e))
    (aggregate["restricted_endpoint_travel"],
     aggregate["restricted_endpoint_state"]) = shared_transform.transform(
         restricted_x, restricted_e)
    natural_endpoint = shared_transform.transform(
        np.asarray([0.0]), np.asarray([0.0]))
    endpoints = dict(
        travel_penalty_sigma=float(restricted_x[-1] - global_x[-1]),
        state_penalty_sigma=float(restricted_e[0] - global_e[0]),
        global_travel_reduction_sigma=float(-global_x[-1]),
        restricted_travel_reduction_sigma=float(-restricted_x[-1]),
        global_state_reduction_sigma=float(-global_e[0]),
        restricted_state_reduction_sigma=float(-restricted_e[0]),
        endpoint_reference=shared_transform.reference_analysis_id,
        endpoint_natural_travel=float(natural_endpoint[0][0]),
        endpoint_natural_state=float(natural_endpoint[1][0]),
        endpoint_metadata=shared_transform.metadata(clipping=False),
    )
    return pd.DataFrame(front_rows), pd.DataFrame(summary_rows), aggregate, endpoints


def build_parser():
    """Build the CLI parser so legacy display defaults are regression-tested."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iteration", type=int, default=ITERATION)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument(
        "--profile",
        choices=("embryo1_legacy", "embryo1_matched", "pooled_tracking_v1"),
        help="Opt into an isolated profile-aware run (omission keeps legacy paths).",
    )
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--display", choices=("null_sd", "endpoint"), default="null_sd",
        help=("Figure coordinate system. Legacy-compatible null_sd is the "
              "default; endpoint must be requested explicitly."),
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    ps.configure()
    context = None
    analysis_out = ANALYSIS_OUT
    if args.profile is not None:
        context = build_analysis_context(
            args.profile, run_id=args.run_id, output_root=args.output_root,
            sweep_intervals=args.iteration)
        context.write()
        analysis_out = context.run_paths.analysis
        if args.out == OUT:
            args.out = context.run_paths.display(args.display)
    analysis_out.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    fronts, summary, aggregate, endpoints = analyze(
        iteration=args.iteration, context=context)
    fronts.to_csv(analysis_out / "within_type_fronts.csv", index=False)
    summary.to_csv(analysis_out / "within_type_summary.csv", index=False)
    aggregate.to_csv(analysis_out / "type_preserving_aggregate_front.csv",
                     index=False)
    display_manifest = {
        "profile": context.profile if context is not None else "embryo1_legacy",
        "display_mode": args.display,
        "aggregate_reference": (
            endpoints["endpoint_metadata"] if args.display == "endpoint"
            else {
                "display_mode": "null_sd",
                "reference_analysis_id": (
                    context.cache_key if context is not None
                    else "legacy_global_first_cousin_null"),
                "natural_lineage": [0.0, 0.0],
            }
        ),
        "type_references": summary[[
            "type", "endpoint_display_valid", "endpoint_display_reference"
        ]].to_dict(orient="records"),
    }
    (args.out / "fig6bc_display_manifest.json").write_text(
        json.dumps(display_manifest, indent=2, sort_keys=True) + "\n")
    plot_within_type(
        fronts, summary, aggregate, endpoints, out_dir=args.out,
        display_mode=args.display)
    print(summary[["type", "n", "cousin_null_valid",
                   "cousin_relative_distance", "full_random_relative_distance",
                   "travel_reduction_per_cell", "state_reduction_per_cell"]]
          .to_string(index=False))
    print("Endpoint comparison:", endpoints)
    print("Wrote within-type analysis and Figure 6 panels B--C to", args.out)


if __name__ == "__main__":
    main()
