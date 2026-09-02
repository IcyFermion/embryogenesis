"""Validated analysis helpers for the full-tree publication figures.

The exploratory notebook remains useful for development history.  This module
turns its C. elegans protein configuration into a deterministic, testable
publication pipeline without importing or executing notebook state.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path
import sys

import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pareto_core import LineageOptimization, LineageTree
from utils import lineage_name_mapping as map_names, load_json


DATA_ROOT = MODULE_DIR / "output" / "internal_opt" / "celegans_prot" / "embryo1"
CACHE_ROOT = MODULE_DIR / "output"
PUBLICATION_ROOT = CACHE_ROOT / "publication"
TRACKING_CUTOFF = 255
TRACKING_SCALE = 0.1625
EARLY_ANCHORS = ("P0", "AB", "P1")
ITERATION = 300
N_WEIGHTS = ITERATION + 1
N_LAYER_NULL = 1_000
N_BROWNIAN_REFERENCE = 10_000
N_DISPLAY_NULL = 1_000
SEED = 42
EXPECTED_TREE_NODES = 1_007
EXPECTED_OPTIMIZED_NODES = 1_004
EXPECTED_INTERNAL = 500
EXPECTED_TERMINAL = 504
EXPECTED_EDGES = 1_000
EXPECTED_LAYER_SIZES = (504, 230, 126, 68, 38, 19, 10, 5)


HEURISTIC_SPECS = (
    ("Bottom-up by layer", "bottom_up_by_layer.npz"),
    ("Layerwise assignment", None),
    ("Paired bottom-up", "paired_bottom_up.npz"),
    ("Degree-constrained spanning forest", None),
    ("Top-down rebuild", "top_down_rebuild.npz"),
    ("Terminal-only rebuild", "terminal_only_no_replacement.npz"),
)

NULL_SPECS = (
    ("First-cousin shuffle", "random_cousin.npz"),
    ("Internal-layer shuffle", "random_by_layer.npz"),
    ("Full assignment shuffle", "random_full.npz"),
    ("Random rebuild", "random_rebuild.npz"),
    # The historical phylo_bm.npz cache was generated for the older 10-PC
    # experiment and is not compatible with the top-20 z-scored analysis.
    ("Parametric Brownian reference", None),
)


@dataclass(frozen=True)
class ContractionLayer:
    """One bottom-up contraction round in the measured lineage forest."""

    round_index: int
    child_tree_ids: np.ndarray
    parent_tree_ids: np.ndarray


@dataclass
class FullTreeContext:
    """Measured C. elegans protein lineage and optimization state."""

    lineage_data: dict
    names: list[str]
    feature_names: list[str]
    tracking_times: dict[str, int]
    measured_tracking: set[str]
    measured_expression: set[str]
    tree: LineageTree
    optimization: LineageOptimization
    layers: list[ContractionLayer]


@dataclass(frozen=True)
class BrownianReference:
    """All-edge plug-in Brownian simulations and their fitted covariance."""

    travel_cost: np.ndarray
    cell_state_cost: np.ndarray
    covariance: np.ndarray
    evaluated_tree_ids: np.ndarray
    root_tree_ids: np.ndarray
    seed: int


def _canonical_maps(lineage_data: dict) -> tuple[dict[str, str | None], dict[str, int]]:
    parents: dict[str, str | None] = {}
    depths: dict[str, int] = {}

    def visit(node: dict, parent: str | None, depth: int) -> None:
        name = map_names(node["did"])
        parents[name] = parent
        depths[name] = depth
        for child in node.get("children", []):
            visit(child, name, depth + 1)

    visit(lineage_data, None, 0)
    return parents, depths


def _contraction_layers(opt: LineageOptimization) -> list[ContractionLayer]:
    """Reproduce the notebook's asynchronous bottom-up contraction rounds."""
    current = list(opt.terminal_tree_ids)
    first_roots = {tree_id for tree_id, _ in opt.first_internal_layer}
    seen = set(current)
    layers: list[ContractionLayer] = []

    while True:
        parents = [opt.lineage_tree.parent_list[tree_id] for tree_id in current]
        layers.append(
            ContractionLayer(
                round_index=len(layers) + 1,
                child_tree_ids=np.asarray(current, dtype=int),
                parent_tree_ids=np.asarray(parents, dtype=int),
            )
        )
        next_layer: list[int] = []
        for parent_id in parents:
            if parent_id in first_roots:
                continue
            if parent_id in seen:
                next_layer.append(parent_id)
            else:
                seen.add(parent_id)
        if not next_layer:
            break
        if len(next_layer) != len(set(next_layer)):
            raise AssertionError("A contraction round contains duplicate child nodes")
        current = next_layer

    evaluated = np.concatenate([layer.child_tree_ids for layer in layers])
    expected = {
        tree_id
        for tree_id in range(opt.lineage_tree.size)
        if tree_id not in first_roots
        and tree_id not in {
            opt.lineage_tree.reverse_lineage_id_mapping[i]
            for i in opt.untracked_lineage_ids
        }
    }
    if len(evaluated) != EXPECTED_EDGES or len(set(evaluated)) != EXPECTED_EDGES:
        raise AssertionError("Contraction rounds do not partition 1,000 edges")
    if tuple(len(layer.child_tree_ids) for layer in layers) != EXPECTED_LAYER_SIZES:
        raise AssertionError("Unexpected contraction-round sizes")
    if set(evaluated) != expected:
        raise AssertionError("Contraction rounds do not match the evaluated forest")
    return layers


def load_ce_protein_context(max_workers: int = 1) -> FullTreeContext:
    """Load the exact notebook configuration and construct its measured tree."""
    lineage_data = load_json(REPO_ROOT / "data" / "cell_lineage.json")

    expression = pd.read_csv(
        REPO_ROOT / "data" / "protein" / "aggregated_all" / "s3_zscore.csv",
        index_col=0,
    ).T
    selected = pd.read_csv(
        REPO_ROOT
        / "expression_embedding"
        / "results"
        / "elegans_protein_linear_baseline"
        / "top20_protein_names.csv"
    )["protein"].tolist()
    available_features = [feature for feature in selected if feature in expression.columns]
    if len(available_features) != 20:
        raise AssertionError("The full-tree protein configuration requires 20 features")
    expression = expression[available_features].fillna(0)

    tracks = pd.read_csv(REPO_ROOT / "data" / "embryo1" / "tracks.txt", sep="\t")
    tracks = tracks[tracks["t"] <= TRACKING_CUTOFF]
    xyz_rows: dict[str, np.ndarray] = {}
    tracking_times: dict[str, int] = {}
    for name in tracks["name"].unique():
        if name not in expression.index:
            continue
        cell_tracks = tracks[tracks["name"] == name]
        time_points = cell_tracks["t"].to_numpy()
        if len(time_points) == 1 and time_points[0] == TRACKING_CUTOFF:
            continue
        xyz_rows[name] = (
            cell_tracks[["x", "y", "z"]].to_numpy(dtype=float)[-1]
            * TRACKING_SCALE
        )
        tracking_times[name] = int(time_points[-1])

    xyz = pd.DataFrame.from_dict(xyz_rows, orient="index", columns=["x", "y", "z"])
    measured_tracking = set(xyz.index)
    measured_expression = set(expression.index)
    for name, time in zip(EARLY_ANCHORS, (-20, -3, -2)):
        tracking_times[name] = time
        xyz.loc[name] = np.repeat(np.nan, 3)
        expression.loc[name] = np.repeat(np.nan, expression.shape[1])
    names = list(EARLY_ANCHORS) + [name for name in xyz_rows]
    name_to_lid = {name: idx for idx, name in enumerate(names)}

    tree = LineageTree()
    tree.add_node(0, -1)
    tree.root = 0
    tree.add_node(1, 0)
    tree.add_node(2, 0)
    queue = deque([(lineage_data, -1)])
    while queue:
        node, represented_parent_lid = queue.popleft()
        name = map_names(node["did"])
        if name in name_to_lid:
            lineage_id = name_to_lid[name]
            if lineage_id >= len(EARLY_ANCHORS):
                parent_tree_id = tree.reverse_lineage_id_mapping[represented_parent_lid]
                tree.add_node(lineage_id, parent_tree_id)
            for child in node.get("children", []):
                queue.append((child, lineage_id))
        else:
            for child in node.get("children", []):
                queue.append((child, represented_parent_lid))

    last_times = [tracking_times.get(name, -1) for name in names]
    tree.record_last_tracked_time(last_times)
    xyz_matrix = np.nan_to_num(xyz.loc[names].to_numpy(dtype=float), nan=0.0)
    exp_matrix = np.nan_to_num(expression.loc[names].to_numpy(dtype=float), nan=0.0)

    first_layer = [
        (child_id, 2)
        for parent_id in (1, 2)
        for child_id in tree.children_list[parent_id]
    ]
    terminal_lids = [
        tree.lineage_id_mapping[tree_id]
        for tree_id in range(tree.size)
        if not tree.children_list[tree_id]
    ]
    type_codes = {names[lineage_id]: 0 for lineage_id in terminal_lids}
    opt = LineageOptimization(
        xyz_matrix,
        exp_matrix,
        tree,
        first_internal_layer=first_layer,
        lineage_names=names,
        lineage_type_code_dict=type_codes,
        exp_norm=2,
        max_workers=max_workers,
    )
    layers = _contraction_layers(opt)
    context = FullTreeContext(
        lineage_data=lineage_data,
        names=names,
        feature_names=["x", "y", "z"] + available_features,
        tracking_times=tracking_times,
        measured_tracking=measured_tracking,
        measured_expression=measured_expression,
        tree=tree,
        optimization=opt,
        layers=layers,
    )
    validate_context(context)
    return context


def validate_context(context: FullTreeContext) -> None:
    opt = context.optimization
    if context.tree.size != EXPECTED_TREE_NODES:
        raise AssertionError(f"Expected {EXPECTED_TREE_NODES} represented nodes")
    if len(opt.internal_tree_ids) != EXPECTED_INTERNAL:
        raise AssertionError(f"Expected {EXPECTED_INTERNAL} evaluated internal cells")
    if len(opt.terminal_tree_ids) != EXPECTED_TERMINAL:
        raise AssertionError(f"Expected {EXPECTED_TERMINAL} evaluated terminal cells")
    if len(opt.internal_tree_ids) + len(opt.terminal_tree_ids) != EXPECTED_OPTIMIZED_NODES:
        raise AssertionError("Unexpected evaluated full-tree node count")
    if not np.isclose(opt.lineage_xyz_cost, 4465.5476125941705):
        raise AssertionError("Natural travel cost differs from the audited notebook value")
    if not np.isclose(opt.lineage_exp_cost, 3102.414696174934):
        raise AssertionError("Natural cell-state cost differs from the audited notebook value")


def node_manifest(context: FullTreeContext) -> pd.DataFrame:
    """Return one row per represented cell, including measurement provenance."""
    canonical_parent, canonical_depth = _canonical_maps(context.lineage_data)
    tree = context.tree
    round_by_tree_id = {
        int(tree_id): layer.round_index
        for layer in context.layers
        for tree_id in layer.child_tree_ids
    }
    represented_depth = {tree.root: 0}
    queue = deque([tree.root])
    while queue:
        parent_id = queue.popleft()
        for child_id in tree.children_list[parent_id]:
            represented_depth[child_id] = represented_depth[parent_id] + 1
            queue.append(child_id)
    first_roots = {tree_id for tree_id, _ in context.optimization.first_internal_layer}
    rows = []
    for tree_id in range(tree.size):
        lineage_id = tree.lineage_id_mapping[tree_id]
        name = context.names[lineage_id]
        parent_tree_id = tree.parent_list[tree_id]
        represented_parent = None
        if parent_tree_id != -1:
            represented_parent = context.names[tree.lineage_id_mapping[parent_tree_id]]
        rows.append(
            {
                "tree_id": tree_id,
                "lineage_id": lineage_id,
                "cell": name,
                "canonical_parent": canonical_parent.get(name),
                "represented_parent": represented_parent,
                "canonical_depth": canonical_depth.get(name),
                "represented_depth": represented_depth[tree_id],
                "tracking_time": context.tracking_times.get(name),
                "measured_tracking": name in context.measured_tracking,
                "measured_expression": name in context.measured_expression,
                "early_anchor": name in EARLY_ANCHORS,
                "terminal_in_measured_tree": not tree.children_list[tree_id],
                "optimization_root": tree_id in first_roots,
                "edge_evaluated": tree_id in round_by_tree_id,
                "contraction_round": round_by_tree_id.get(tree_id),
            }
        )
    result = pd.DataFrame(rows)
    if int(result["edge_evaluated"].sum()) != EXPECTED_EDGES:
        raise AssertionError("Node manifest does not contain exactly 1,000 evaluated edges")
    return result


def _nondominated_mask(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Mark one representative of each sampled non-dominated cost pair."""
    points = np.column_stack((np.asarray(x, float), np.asarray(y, float)))
    unique, first = np.unique(points, axis=0, return_index=True)
    order = np.lexsort((unique[:, 1], unique[:, 0]))
    keep_unique = np.zeros(len(unique), dtype=bool)
    best_y = np.inf
    for unique_idx in order:
        candidate_y = unique[unique_idx, 1]
        if candidate_y < best_y - 1e-10:
            keep_unique[unique_idx] = True
            best_y = candidate_y
    mask = np.zeros(len(points), dtype=bool)
    mask[first[keep_unique]] = True
    return mask


def compute_layerwise_analysis(
    context: FullTreeContext,
    iteration: int = ITERATION,
    n_null: int = N_LAYER_NULL,
    seed: int = SEED,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute exact per-round assignments and matched slot-shuffle nulls."""
    opt = context.optimization
    weights = np.linspace(0.0, 1.0, iteration + 1)
    layer_fronts: list[pd.DataFrame] = []
    layer_nulls: list[pd.DataFrame] = []
    manifest_rows = []
    aggregate_xyz = np.zeros(len(weights))
    aggregate_exp = np.zeros(len(weights))
    aggregate_retained = np.zeros(len(weights))
    aggregate_null_xyz = np.zeros(n_null)
    aggregate_null_exp = np.zeros(n_null)

    for layer in context.layers:
        child_lids = np.asarray(
            [context.tree.lineage_id_mapping[t] for t in layer.child_tree_ids],
            dtype=int,
        )
        parent_lids = np.asarray(
            [context.tree.lineage_id_mapping[t] for t in layer.parent_tree_ids],
            dtype=int,
        )
        n_edges = len(child_lids)
        rows = np.arange(n_edges)
        xyz_cost = opt.xyz_cost_mat[np.ix_(child_lids, parent_lids)]
        exp_cost = opt.exp_cost_mat[np.ix_(child_lids, parent_lids)]
        natural_xyz = float(xyz_cost[rows, rows].sum())
        natural_exp = float(exp_cost[rows, rows].sum())

        front_xyz = np.empty(len(weights))
        front_exp = np.empty(len(weights))
        retention = np.empty(len(weights))
        for weight_idx, alpha in enumerate(weights):
            assignment_cost = alpha * xyz_cost + (1.0 - alpha) * exp_cost
            row_idx, col_idx = linear_sum_assignment(assignment_cost)
            front_xyz[weight_idx] = xyz_cost[row_idx, col_idx].sum()
            front_exp[weight_idx] = exp_cost[row_idx, col_idx].sum()
            assigned_parents = parent_lids[col_idx]
            retention[weight_idx] = np.mean(assigned_parents == parent_lids[row_idx])

        rng = np.random.default_rng(seed + layer.round_index)
        null_xyz = np.empty(n_null)
        null_exp = np.empty(n_null)
        for draw in range(n_null):
            permutation = rng.permutation(n_edges)
            null_xyz[draw] = xyz_cost[rows, permutation].sum()
            null_exp[draw] = exp_cost[rows, permutation].sum()
        xyz_std = float(null_xyz.std())
        exp_std = float(null_exp.std())
        if xyz_std <= 0 or exp_std <= 0:
            raise AssertionError("Layer shuffle has non-positive objective variance")

        scope = f"round_{layer.round_index}"
        layer_df = pd.DataFrame(
            {
                "scope": scope,
                "round": layer.round_index,
                "weight_index": np.arange(len(weights)),
                "alpha_travel": weights,
                "travel_cost": front_xyz,
                "cell_state_cost": front_exp,
                "edge_retention": retention,
                "travel_standardized": (front_xyz - natural_xyz) / xyz_std,
                "cell_state_standardized": (front_exp - natural_exp) / exp_std,
            }
        )
        layer_df["sampled_nondominated"] = _nondominated_mask(front_xyz, front_exp)
        layer_fronts.append(layer_df)
        layer_nulls.append(
            pd.DataFrame(
                {
                    "scope": scope,
                    "round": layer.round_index,
                    "draw": np.arange(n_null),
                    "travel_cost": null_xyz,
                    "cell_state_cost": null_exp,
                    "travel_standardized": (null_xyz - natural_xyz) / xyz_std,
                    "cell_state_standardized": (null_exp - natural_exp) / exp_std,
                }
            )
        )

        parent_counts = Counter(layer.parent_tree_ids.tolist())
        dominated_natural = bool(
            np.any(
                (front_xyz <= natural_xyz + 1e-10)
                & (front_exp <= natural_exp + 1e-10)
                & (
                    (front_xyz < natural_xyz - 1e-10)
                    | (front_exp < natural_exp - 1e-10)
                )
            )
        )
        manifest_rows.append(
            {
                "scope": scope,
                "round": layer.round_index,
                "position": "bottom" if layer.round_index == 1 else (
                    "top" if layer.round_index == len(context.layers) else "intermediate"
                ),
                "edges": n_edges,
                "unique_parents": len(parent_counts),
                "parents_with_one_slot": sum(count == 1 for count in parent_counts.values()),
                "parents_with_two_slots": sum(count == 2 for count in parent_counts.values()),
                "natural_travel_cost": natural_xyz,
                "natural_cell_state_cost": natural_exp,
                "null_travel_mean": float(null_xyz.mean()),
                "null_travel_std": xyz_std,
                "null_cell_state_mean": float(null_exp.mean()),
                "null_cell_state_std": exp_std,
                "natural_sampled_dominated": dominated_natural,
                "sampled_unique_solutions": int(
                    len(np.unique(np.column_stack((front_xyz, front_exp)), axis=0))
                ),
                "maximum_edge_retention": float(retention.max()),
            }
        )
        aggregate_xyz += front_xyz
        aggregate_exp += front_exp
        aggregate_retained += retention * n_edges
        aggregate_null_xyz += null_xyz
        aggregate_null_exp += null_exp

    natural_xyz = float(sum(row["natural_travel_cost"] for row in manifest_rows))
    natural_exp = float(sum(row["natural_cell_state_cost"] for row in manifest_rows))
    if not np.isclose(natural_xyz, opt.lineage_xyz_cost):
        raise AssertionError("Layer natural travel costs do not recompose the tree")
    if not np.isclose(natural_exp, opt.lineage_exp_cost):
        raise AssertionError("Layer natural cell-state costs do not recompose the tree")
    aggregate_xyz_std = float(aggregate_null_xyz.std())
    aggregate_exp_std = float(aggregate_null_exp.std())
    aggregate_df = pd.DataFrame(
        {
            "scope": "aggregate",
            "round": 0,
            "weight_index": np.arange(len(weights)),
            "alpha_travel": weights,
            "travel_cost": aggregate_xyz,
            "cell_state_cost": aggregate_exp,
            "edge_retention": aggregate_retained / EXPECTED_EDGES,
            "travel_standardized": (aggregate_xyz - natural_xyz) / aggregate_xyz_std,
            "cell_state_standardized": (aggregate_exp - natural_exp) / aggregate_exp_std,
        }
    )
    aggregate_df["sampled_nondominated"] = _nondominated_mask(
        aggregate_xyz, aggregate_exp
    )
    layer_fronts.append(aggregate_df)
    layer_nulls.append(
        pd.DataFrame(
            {
                "scope": "aggregate",
                "round": 0,
                "draw": np.arange(n_null),
                "travel_cost": aggregate_null_xyz,
                "cell_state_cost": aggregate_null_exp,
                "travel_standardized": (
                    aggregate_null_xyz - natural_xyz
                ) / aggregate_xyz_std,
                "cell_state_standardized": (
                    aggregate_null_exp - natural_exp
                ) / aggregate_exp_std,
            }
        )
    )
    aggregate_dominated = bool(
        np.any(
            (aggregate_xyz <= natural_xyz + 1e-10)
            & (aggregate_exp <= natural_exp + 1e-10)
            & (
                (aggregate_xyz < natural_xyz - 1e-10)
                | (aggregate_exp < natural_exp - 1e-10)
            )
        )
    )
    manifest_rows.append(
        {
            "scope": "aggregate",
            "round": 0,
            "position": "full tree",
            "edges": EXPECTED_EDGES,
            "unique_parents": EXPECTED_INTERNAL,
            "parents_with_one_slot": np.nan,
            "parents_with_two_slots": np.nan,
            "natural_travel_cost": natural_xyz,
            "natural_cell_state_cost": natural_exp,
            "null_travel_mean": float(aggregate_null_xyz.mean()),
            "null_travel_std": aggregate_xyz_std,
            "null_cell_state_mean": float(aggregate_null_exp.mean()),
            "null_cell_state_std": aggregate_exp_std,
            "natural_sampled_dominated": aggregate_dominated,
            "sampled_unique_solutions": int(
                len(np.unique(np.column_stack((aggregate_xyz, aggregate_exp)), axis=0))
            ),
            "maximum_edge_retention": float(
                (aggregate_retained / EXPECTED_EDGES).max()
            ),
        }
    )

    fronts = pd.concat(layer_fronts, ignore_index=True)
    nulls = pd.concat(layer_nulls, ignore_index=True)
    manifest = pd.DataFrame(manifest_rows)
    validate_layerwise(context, fronts, nulls, manifest, iteration, n_null)
    return fronts, nulls, manifest


def validate_layerwise(
    context: FullTreeContext,
    fronts: pd.DataFrame,
    nulls: pd.DataFrame,
    manifest: pd.DataFrame,
    iteration: int,
    n_null: int,
) -> None:
    scopes = [f"round_{i}" for i in range(1, 9)] + ["aggregate"]
    if set(fronts["scope"]) != set(scopes) or set(nulls["scope"]) != set(scopes):
        raise AssertionError("Layerwise outputs are missing a required scope")
    if not all(len(fronts[fronts["scope"] == scope]) == iteration + 1 for scope in scopes):
        raise AssertionError("Every layerwise scope must contain every weight")
    if not all(len(nulls[nulls["scope"] == scope]) == n_null for scope in scopes):
        raise AssertionError("Every layerwise scope must contain every null draw")
    if len(manifest) != len(scopes):
        raise AssertionError("Unexpected layer manifest length")
    numeric = fronts[
        [
            "travel_cost",
            "cell_state_cost",
            "edge_retention",
            "travel_standardized",
            "cell_state_standardized",
        ]
    ].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise AssertionError("Layerwise front contains non-finite values")
    if not fronts["edge_retention"].between(0, 1).all():
        raise AssertionError("Layerwise edge retention is outside [0, 1]")
    aggregate = fronts[fronts["scope"] == "aggregate"].sort_values("weight_index")
    per_round = [
        fronts[fronts["scope"] == f"round_{i}"].sort_values("weight_index")
        for i in range(1, 9)
    ]
    if not np.allclose(
        aggregate["travel_cost"],
        np.sum([frame["travel_cost"].to_numpy() for frame in per_round], axis=0),
    ):
        raise AssertionError("Aggregate travel front does not equal layer sum")
    if not np.allclose(
        aggregate["cell_state_cost"],
        np.sum([frame["cell_state_cost"].to_numpy() for frame in per_round], axis=0),
    ):
        raise AssertionError("Aggregate cell-state front does not equal layer sum")

    legacy = np.load(
        DATA_ROOT / "heuristics" / "layerwise_pareto_assignment.npz",
        allow_pickle=True,
    )
    for weight_idx, legacy_idx in ((0, 0), (iteration // 2, 500), (iteration, 1000)):
        row = aggregate.iloc[weight_idx]
        if not np.isclose(row["travel_cost"], legacy["xyz"][legacy_idx]):
            raise AssertionError("Layerwise travel landmarks differ from legacy cache")
        if not np.isclose(row["cell_state_cost"], legacy["exp"][legacy_idx]):
            raise AssertionError("Layerwise cell-state landmarks differ from legacy cache")


def degree_constrained_spanning_forest(
    context: FullTreeContext,
    iteration: int = ITERATION,
) -> pd.DataFrame:
    """Construct a valid four-root binary forest with a Kruskal heuristic.

    The notebook's saved unconstrained-MST cache includes both internal and
    terminal costs, but its reconstructed graph is a single unrooted component
    and can give one internal cell more than two children.  The current legacy
    ``mst_rebuild`` implementation is additionally inconsistent with that cache:
    an early return reports terminal-assignment costs only.  Here, greedy
    degree-constrained Kruskal preserves four fixed roots, prevents components
    from joining two roots, caps every rooted parent at two children, and then
    assigns all terminals to the remaining openings.  Every sampled solution
    therefore contains exactly the same 1,000 evaluated edges as the natural
    lineage.
    """
    opt = context.optimization
    internal = np.asarray(
        [context.tree.lineage_id_mapping[t] for t in opt.internal_tree_ids],
        dtype=int,
    )
    terminals = np.asarray(
        [context.tree.lineage_id_mapping[t] for t in opt.terminal_tree_ids],
        dtype=int,
    )
    roots = {
        context.tree.lineage_id_mapping[t]
        for t, _ in opt.first_internal_layer
    }
    local_root = np.asarray([lineage_id in roots for lineage_id in internal])
    edge_u, edge_v = np.triu_indices(len(internal), 1)
    weights = np.linspace(0.0, 1.0, iteration + 1)
    records = []

    natural_parent = {}
    for tree_id in opt.internal_tree_ids + opt.terminal_tree_ids:
        if tree_id in {t for t, _ in opt.first_internal_layer}:
            continue
        child_lid = context.tree.lineage_id_mapping[tree_id]
        parent_lid = context.tree.lineage_id_mapping[context.tree.parent_list[tree_id]]
        natural_parent[child_lid] = parent_lid

    for weight_idx, alpha in enumerate(weights):
        combined = alpha * opt.xyz_cost_mat + (1.0 - alpha) * opt.exp_cost_mat
        edge_weights = combined[internal[edge_u], internal[edge_v]]
        order = np.argsort(edge_weights, kind="stable")

        union_parent = np.arange(len(internal))
        union_rank = np.zeros(len(internal), dtype=np.int8)
        degree = np.zeros(len(internal), dtype=np.int8)
        roots_per_component = local_root.astype(np.int8)

        def find(local_id: int) -> int:
            while union_parent[local_id] != local_id:
                union_parent[local_id] = union_parent[union_parent[local_id]]
                local_id = int(union_parent[local_id])
            return local_id

        selected_local: list[tuple[int, int]] = []
        for order_idx in order:
            u = int(edge_u[order_idx])
            v = int(edge_v[order_idx])
            u_cap = 2 if local_root[u] else 3
            v_cap = 2 if local_root[v] else 3
            if degree[u] >= u_cap or degree[v] >= v_cap:
                continue
            root_u = find(u)
            root_v = find(v)
            if root_u == root_v:
                continue
            if roots_per_component[root_u] + roots_per_component[root_v] > 1:
                continue
            if union_rank[root_u] < union_rank[root_v]:
                root_u, root_v = root_v, root_u
            union_parent[root_v] = root_u
            roots_per_component[root_u] += roots_per_component[root_v]
            if union_rank[root_u] == union_rank[root_v]:
                union_rank[root_u] += 1
            degree[u] += 1
            degree[v] += 1
            selected_local.append((u, v))
            if len(selected_local) == EXPECTED_INTERNAL - len(roots):
                break

        graph = nx.Graph()
        graph.add_nodes_from(internal.tolist())
        graph.add_edges_from(
            (int(internal[u]), int(internal[v])) for u, v in selected_local
        )
        components = list(nx.connected_components(graph))
        if len(selected_local) != EXPECTED_INTERNAL - len(roots):
            raise AssertionError("Spanning-forest Kruskal did not reach four components")
        if len(components) != len(roots):
            raise AssertionError("Spanning forest does not have four components")
        if any(sum(node in roots for node in component) != 1 for component in components):
            raise AssertionError("Every spanning-forest component must contain one root")

        child_count: dict[int, int] = {}
        directed_internal_edges: list[tuple[int, int]] = []
        for root in roots:
            component = nx.node_connected_component(graph, root)
            subgraph = graph.subgraph(component)
            for parent, children in nx.bfs_successors(subgraph, root):
                child_count[int(parent)] = len(children)
                for child in children:
                    child_count.setdefault(int(child), 0)
                    directed_internal_edges.append((int(parent), int(child)))
        if max(child_count.values()) > 2:
            raise AssertionError("Spanning forest violates the binary child cap")
        openings = np.asarray(
            [
                lineage_id
                for lineage_id in internal
                for _ in range(2 - child_count.get(int(lineage_id), 0))
            ],
            dtype=int,
        )
        if len(openings) != EXPECTED_TERMINAL:
            raise AssertionError("Spanning forest must expose 504 terminal openings")

        terminal_cost = combined[np.ix_(terminals, openings)]
        terminal_rows, terminal_cols = linear_sum_assignment(terminal_cost)
        assigned_openings = openings[terminal_cols]
        internal_u = np.asarray([edge[0] for edge in directed_internal_edges], dtype=int)
        internal_v = np.asarray([edge[1] for edge in directed_internal_edges], dtype=int)
        travel = float(
            opt.xyz_cost_mat[internal_u, internal_v].sum()
            + opt.xyz_cost_mat[terminals[terminal_rows], assigned_openings].sum()
        )
        cell_state = float(
            opt.exp_cost_mat[internal_u, internal_v].sum()
            + opt.exp_cost_mat[terminals[terminal_rows], assigned_openings].sum()
        )
        retained = sum(
            natural_parent.get(child) == parent
            for parent, child in directed_internal_edges
        )
        retained += sum(
            natural_parent.get(int(child)) == int(parent)
            for child, parent in zip(terminals[terminal_rows], assigned_openings)
        )
        records.append(
            {
                "weight_index": weight_idx,
                "alpha_travel": alpha,
                "travel_cost": travel,
                "cell_state_cost": cell_state,
                "edge_retention": retained / EXPECTED_EDGES,
                "internal_edges": len(directed_internal_edges),
                "terminal_edges": len(terminals),
                "components": len(components),
                "maximum_children": max(child_count.values()),
            }
        )

    result = pd.DataFrame(records)
    result["sampled_nondominated"] = _nondominated_mask(
        result["travel_cost"].to_numpy(), result["cell_state_cost"].to_numpy()
    )
    if not (result["internal_edges"] + result["terminal_edges"] == EXPECTED_EDGES).all():
        raise AssertionError("Spanning-forest edge totals are not comparable")
    if not (result["components"] == 4).all() or not (result["maximum_children"] <= 2).all():
        raise AssertionError("Spanning-forest structural validation failed")
    return result


def _mark_collective_nondominated(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["collective_nondominated"] = False
    mask = _nondominated_mask(
        result["travel_cost"].to_numpy(), result["cell_state_cost"].to_numpy()
    )
    result.loc[mask, "collective_nondominated"] = True
    return result


def parametric_brownian_reference(
    context: FullTreeContext,
    n_samples: int = N_BROWNIAN_REFERENCE,
    seed: int = SEED + 200,
) -> BrownianReference:
    """Generate an all-edge plug-in Brownian reference on the fixed topology.

    The full 23-dimensional rate covariance is estimated from the 1,000
    observed, branch-standardized increments.  Each replicate fixes the four
    optimization-root states and draws independent Gaussian increments down
    every scored edge.  This is a parametric reference for Pareto geometry,
    not an ancestral-state reconstruction or a realistic embryogenesis model.
    """
    tree = context.tree
    opt = context.optimization
    evaluated_tree_ids = np.sort(
        np.concatenate([layer.child_tree_ids for layer in context.layers])
    ).astype(int)
    root_tree_ids = np.asarray(
        sorted({tree_id for tree_id, _ in opt.first_internal_layer}), dtype=int
    )
    if len(evaluated_tree_ids) != EXPECTED_EDGES or len(root_tree_ids) != 4:
        raise AssertionError("Brownian reference scope does not match Figure 8")
    if set(evaluated_tree_ids) & set(root_tree_ids):
        raise AssertionError("Brownian roots must not also be scored descendants")

    tree_to_lineage = np.asarray(tree.lineage_id_mapping, dtype=int)
    values = np.hstack(
        (opt.xyz_mat[tree_to_lineage], opt.exp_mat[tree_to_lineage])
    )
    n_features = values.shape[1]
    if n_features != len(context.feature_names) or n_features != 23:
        raise AssertionError("Brownian reference requires 3 spatial + 20 protein features")

    parent_tree_ids = np.asarray(
        [tree.parent_list[int(tree_id)] for tree_id in evaluated_tree_ids],
        dtype=int,
    )
    branch_length = np.asarray(tree.branch_time_length, dtype=float)
    evaluated_time = branch_length[evaluated_tree_ids]
    if np.any(evaluated_time <= 0) or not np.isfinite(evaluated_time).all():
        raise AssertionError("Brownian sampling requires positive branch lengths")

    measured_increment = values[evaluated_tree_ids] - values[parent_tree_ids]
    direct_travel = np.linalg.norm(measured_increment[:, :3], axis=1).sum()
    direct_cell_state = np.linalg.norm(measured_increment[:, 3:], axis=1).sum()
    if not np.isclose(direct_travel, opt.lineage_xyz_cost):
        raise AssertionError("Brownian edge manifest does not reproduce travel cost")
    if not np.isclose(direct_cell_state, opt.lineage_exp_cost):
        raise AssertionError("Brownian edge manifest does not reproduce cell-state cost")

    standardized = measured_increment / np.sqrt(evaluated_time)[:, None]
    covariance = standardized.T @ standardized / len(standardized)
    covariance = (covariance + covariance.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = np.finfo(float).eps * n_features * eigenvalues[-1]
    if not np.isfinite(covariance).all() or eigenvalues[0] <= tolerance:
        raise AssertionError("Fitted Brownian covariance is not positive definite")
    covariance_factor = eigenvectors * np.sqrt(eigenvalues)
    if not np.allclose(
        covariance_factor @ covariance_factor.T,
        covariance,
        rtol=1e-10,
        atol=1e-12,
    ):
        raise AssertionError("Brownian covariance square root failed validation")

    rng = np.random.default_rng(seed)
    travel_cost = np.empty(n_samples, dtype=float)
    cell_state_cost = np.empty(n_samples, dtype=float)
    batch_size = min(100, n_samples)
    for start in range(0, n_samples, batch_size):
        stop = min(start + batch_size, n_samples)
        batch = stop - start
        states = np.full((tree.size, batch, n_features), np.nan, dtype=float)
        states[root_tree_ids] = values[root_tree_ids, None, :]
        batch_travel = np.zeros(batch)
        batch_cell_state = np.zeros(batch)
        for tree_id in evaluated_tree_ids:
            parent_id = tree.parent_list[tree_id]
            if not np.isfinite(states[parent_id]).all():
                raise AssertionError("Brownian simulation order is not root-to-tip")
            increment = (
                rng.normal(size=(batch, n_features)) @ covariance_factor.T
            ) * np.sqrt(branch_length[tree_id])
            states[tree_id] = states[parent_id] + increment
            batch_travel += np.linalg.norm(increment[:, :3], axis=1)
            batch_cell_state += np.linalg.norm(increment[:, 3:], axis=1)
        travel_cost[start:stop] = batch_travel
        cell_state_cost[start:stop] = batch_cell_state

    if not np.isfinite(travel_cost).all() or not np.isfinite(cell_state_cost).all():
        raise AssertionError("Parametric Brownian reference contains non-finite costs")
    return BrownianReference(
        travel_cost=travel_cost,
        cell_state_cost=cell_state_cost,
        covariance=covariance,
        evaluated_tree_ids=evaluated_tree_ids,
        root_tree_ids=root_tree_ids,
        seed=seed,
    )


def collective_analysis(
    context: FullTreeContext,
    layerwise_fronts: pd.DataFrame,
    spanning_forest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, BrownianReference]:
    """Combine validated heuristic fronts and all notebook null models."""
    natural_xyz = context.optimization.lineage_xyz_cost
    natural_exp = context.optimization.lineage_exp_cost
    cousin = np.load(DATA_ROOT / "null" / "random_cousin.npz")
    xyz_std = float(cousin["xyz"].std())
    exp_std = float(cousin["exp"].std())
    if xyz_std <= 0 or exp_std <= 0:
        raise AssertionError("First-cousin null has non-positive variance")

    heuristic_frames = []
    for name, filename in HEURISTIC_SPECS:
        if name == "Layerwise assignment":
            source = layerwise_fronts[layerwise_fronts["scope"] == "aggregate"].copy()
            x = source["travel_cost"].to_numpy(float)
            y = source["cell_state_cost"].to_numpy(float)
            alpha = source["alpha_travel"].to_numpy(float)
            source_cache = "publication layerwise analysis"
        elif name == "Degree-constrained spanning forest":
            x = spanning_forest["travel_cost"].to_numpy(float)
            y = spanning_forest["cell_state_cost"].to_numpy(float)
            alpha = spanning_forest["alpha_travel"].to_numpy(float)
            source_cache = "publication spanning-forest analysis"
        else:
            cache = np.load(DATA_ROOT / "heuristics" / filename, allow_pickle=True)
            x = np.asarray(cache["xyz"], dtype=float)
            y = np.asarray(cache["exp"], dtype=float)
            alpha = np.linspace(0.0, 1.0, len(x))
            source_cache = filename
        frame = pd.DataFrame(
            {
                "heuristic": name,
                "weight_index": np.arange(len(x)),
                "alpha_travel": alpha,
                "travel_cost": x,
                "cell_state_cost": y,
                "travel_standardized": (x - natural_xyz) / xyz_std,
                "cell_state_standardized": (y - natural_exp) / exp_std,
                "source_cache": source_cache,
            }
        )
        frame["within_heuristic_nondominated"] = _nondominated_mask(x, y)
        heuristic_frames.append(frame)
    heuristics = _mark_collective_nondominated(
        pd.concat(heuristic_frames, ignore_index=True)
    )

    brownian = parametric_brownian_reference(context)
    null_frames = []
    null_summary = []
    for null_index, (name, filename) in enumerate(NULL_SPECS):
        if name == "Parametric Brownian reference":
            x = brownian.travel_cost
            y = brownian.cell_state_cost
            sample_idx = np.arange(len(x))
            displayed = sample_idx < min(N_DISPLAY_NULL, len(x))
            source = "all-edge fixed-topology parametric bootstrap"
        else:
            cache = np.load(DATA_ROOT / "null" / filename)
            x = np.asarray(cache["xyz"], dtype=float)
            y = np.asarray(cache["exp"], dtype=float)
            rng = np.random.default_rng(SEED + 100 + null_index)
            sample_idx = np.sort(
                rng.choice(len(x), size=min(N_DISPLAY_NULL, len(x)), replace=False)
            )
            displayed = np.ones(len(sample_idx), dtype=bool)
            source = filename
        null_frames.append(
            pd.DataFrame(
                {
                    "null_model": name,
                    "draw": sample_idx,
                    "travel_cost": x[sample_idx],
                    "cell_state_cost": y[sample_idx],
                    "travel_standardized": (x[sample_idx] - natural_xyz) / xyz_std,
                    "cell_state_standardized": (y[sample_idx] - natural_exp) / exp_std,
                    "displayed": displayed,
                    "source": source,
                }
            )
        )
        null_summary.append(
            {
                "null_model": name,
                "draws": len(x),
                "travel_mean": float(x.mean()),
                "travel_std": float(x.std()),
                "cell_state_mean": float(y.mean()),
                "cell_state_std": float(y.std()),
                "travel_mean_standardized": float((x.mean() - natural_xyz) / xyz_std),
                "cell_state_mean_standardized": float((y.mean() - natural_exp) / exp_std),
                "correlation": float(np.corrcoef(x, y)[0, 1]),
            }
        )
    nulls = pd.concat(null_frames, ignore_index=True)
    summary = pd.DataFrame(null_summary)
    validate_collective(heuristics, nulls, summary)
    return heuristics, nulls, summary, brownian


def validate_collective(
    heuristics: pd.DataFrame,
    nulls: pd.DataFrame,
    null_summary: pd.DataFrame,
) -> None:
    expected_heuristics = {name for name, _ in HEURISTIC_SPECS}
    expected_nulls = {name for name, _ in NULL_SPECS}
    if set(heuristics["heuristic"]) != expected_heuristics:
        raise AssertionError("Collective cache is missing a heuristic")
    if set(nulls["null_model"]) != expected_nulls:
        raise AssertionError("Collective cache is missing a null model")
    if set(null_summary["null_model"]) != expected_nulls:
        raise AssertionError("Null summary is incomplete")
    numeric = heuristics[
        [
            "travel_cost",
            "cell_state_cost",
            "travel_standardized",
            "cell_state_standardized",
        ]
    ].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise AssertionError("Collective heuristic cache contains non-finite values")
    if not heuristics["collective_nondominated"].any():
        raise AssertionError("Collective sampled front is empty")
    brownian = nulls[nulls["null_model"] == "Parametric Brownian reference"]
    if not (
        brownian["source"] == "all-edge fixed-topology parametric bootstrap"
    ).all():
        raise AssertionError("Collective analysis reused the legacy Brownian cache")
    if len(brownian) != N_BROWNIAN_REFERENCE:
        raise AssertionError("Parametric Brownian reference has the wrong draw count")
    if int(brownian["displayed"].sum()) != N_DISPLAY_NULL:
        raise AssertionError("Parametric Brownian display subsample has the wrong size")
    # The audited replacement must not silently fall back to the historical
    # partial-cost MST cache.
    if heuristics["source_cache"].astype(str).str.contains("mst_rebuild").any():
        raise AssertionError("Invalid legacy MST costs entered publication analysis")


def build_publication_caches(
    iteration: int = ITERATION,
    n_layer_null: int = N_LAYER_NULL,
) -> dict[str, pd.DataFrame]:
    """Run the validated analysis and write publication-oriented CSV caches."""
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    context = load_ce_protein_context()
    nodes = node_manifest(context)
    layer_fronts, layer_nulls, layer_manifest = compute_layerwise_analysis(
        context, iteration=iteration, n_null=n_layer_null
    )
    spanning = degree_constrained_spanning_forest(context, iteration=iteration)
    heuristics, collective_nulls, null_summary, brownian = collective_analysis(
        context, layer_fronts, spanning
    )
    covariance = pd.DataFrame(
        brownian.covariance,
        index=context.feature_names,
        columns=context.feature_names,
    ).rename_axis("feature").reset_index()
    eigenvalues = np.linalg.eigvalsh(brownian.covariance)
    evaluated_time = np.asarray(context.tree.branch_time_length, dtype=float)[
        brownian.evaluated_tree_ids
    ]
    tree_to_lineage = np.asarray(context.tree.lineage_id_mapping, dtype=int)
    measured_state = np.hstack(
        (
            context.optimization.xyz_mat[tree_to_lineage],
            context.optimization.exp_mat[tree_to_lineage],
        )
    )
    evaluated_parent = np.asarray(
        [context.tree.parent_list[int(i)] for i in brownian.evaluated_tree_ids]
    )
    standardized_increment = (
        measured_state[brownian.evaluated_tree_ids]
        - measured_state[evaluated_parent]
    ) / np.sqrt(evaluated_time)[:, None]
    travel_rate_sq = np.square(standardized_increment[:, :3]).sum(axis=1)
    cell_state_rate_sq = np.square(standardized_increment[:, 3:]).sum(axis=1)

    def top_share(values: np.ndarray, fraction: float) -> float:
        count = max(1, int(np.ceil(fraction * len(values))))
        return float(np.sort(values)[-count:].sum() / values.sum())

    brownian_diagnostics = pd.DataFrame(
        [
            ("generator", "all-edge fixed-topology parametric bootstrap"),
            ("interpretation", "geometric reference; not a realistic embryogenesis model"),
            ("draws", len(brownian.travel_cost)),
            ("displayed_draws", N_DISPLAY_NULL),
            ("seed", brownian.seed),
            ("evaluated_edges", len(brownian.evaluated_tree_ids)),
            ("fixed_roots", len(brownian.root_tree_ids)),
            ("features", brownian.covariance.shape[0]),
            ("covariance_rank", np.linalg.matrix_rank(brownian.covariance)),
            ("covariance_min_eigenvalue", eigenvalues[0]),
            ("covariance_max_eigenvalue", eigenvalues[-1]),
            ("covariance_condition_number", np.linalg.cond(brownian.covariance)),
            ("minimum_branch_time", evaluated_time.min()),
            ("maximum_branch_time", evaluated_time.max()),
            (
                "standardized_mean_norm_travel",
                np.linalg.norm(standardized_increment[:, :3].mean(axis=0)),
            ),
            (
                "standardized_mean_norm_cell_state",
                np.linalg.norm(standardized_increment[:, 3:].mean(axis=0)),
            ),
            (
                "branch_time_vs_travel_rate_sq_correlation",
                np.corrcoef(evaluated_time, travel_rate_sq)[0, 1],
            ),
            (
                "branch_time_vs_cell_state_rate_sq_correlation",
                np.corrcoef(evaluated_time, cell_state_rate_sq)[0, 1],
            ),
            ("top_1pct_cell_state_rate_sq_share", top_share(cell_state_rate_sq, 0.01)),
            ("top_5pct_cell_state_rate_sq_share", top_share(cell_state_rate_sq, 0.05)),
            ("top_10pct_cell_state_rate_sq_share", top_share(cell_state_rate_sq, 0.10)),
            ("natural_travel_cost", context.optimization.lineage_xyz_cost),
            ("natural_cell_state_cost", context.optimization.lineage_exp_cost),
            ("reference_travel_mean", brownian.travel_cost.mean()),
            ("reference_travel_std", brownian.travel_cost.std()),
            ("reference_cell_state_mean", brownian.cell_state_cost.mean()),
            ("reference_cell_state_std", brownian.cell_state_cost.std()),
            (
                "reference_cost_correlation",
                np.corrcoef(brownian.travel_cost, brownian.cell_state_cost)[0, 1],
            ),
        ],
        columns=["metric", "value"],
    )
    outputs = {
        "ce_full_tree_node_manifest.csv": nodes,
        "ce_full_tree_layer_manifest.csv": layer_manifest,
        "ce_full_tree_layerwise_fronts.csv": layer_fronts,
        "ce_full_tree_layerwise_nulls.csv": layer_nulls,
        "ce_full_tree_spanning_forest.csv": spanning,
        "ce_full_tree_collective_heuristics.csv": heuristics,
        "ce_full_tree_collective_nulls.csv": collective_nulls,
        "ce_full_tree_null_summary.csv": null_summary,
        "ce_full_tree_brownian_covariance.csv": covariance,
        "ce_full_tree_brownian_diagnostics.csv": brownian_diagnostics,
    }
    for filename, frame in outputs.items():
        frame.to_csv(CACHE_ROOT / filename, index=False)
    return outputs
