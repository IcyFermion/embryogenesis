"""Cacheable global terminal analysis shared by numerical and display layers."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto.analysis_context import AnalysisContext
from terminal_pareto.subtree_analysis import first_cousin_null_summary


GLOBAL_ANALYSIS_VERSION = "global-terminal-v1"


@dataclass
class GlobalAnalysisResult:
    context_cache_key: str
    analysis_cache_key: str
    twr: dict
    first_stats: dict
    raw_front_travel: np.ndarray
    raw_front_state: np.ndarray
    assignments: np.ndarray
    null_raw: dict[str, tuple[np.ndarray, np.ndarray]]

    @property
    def natural_costs(self) -> tuple[float, float]:
        return (float(self.first_stats["lineage_xyz"]),
                float(self.first_stats["lineage_exp"]))

    @property
    def null_stds(self) -> tuple[float, float]:
        return (float(self.first_stats["xyz_std"]),
                float(self.first_stats["exp_std"]))


def _analysis_cache_key(context: AnalysisContext, *, n_lineage_null: int) -> str:
    payload = {
        "version": GLOBAL_ANALYSIS_VERSION,
        "context_cache_key": context.cache_key,
        "sweep_intervals": context.spec.sweep_intervals,
        "null_draws": context.spec.null_draws,
        "lineage_null_draws": int(n_lineage_null),
        "seed": context.spec.seed,
        "scaling_policy": context.spec.scaling_policy,
    }
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def compute_global_analysis(context: AnalysisContext, *,
                            n_lineage_null: int = 100) -> GlobalAnalysisResult:
    """Compute global fronts/nulls without consulting any legacy cache."""
    travel, state = context.cost_matrices()
    nodes = context.terminal_nodes
    parents = context.terminal_parents
    first_groups = pe.build_ancestor_groups(nodes, context.tree_index, 2)
    if context.spec.scaling_policy == "legacy_sampled_1000":
        first_stats = pe.compute_cousin_random_stats(
            travel, state, first_groups,
            n_random=context.spec.null_draws, seed=context.spec.seed)
    else:
        first_stats = first_cousin_null_summary(
            travel, state, nodes, context.gp_map,
            n_random=context.spec.null_draws,
            seed=context.spec.seed)["_raw"]

    twr = lm.combined_lineage_proximity(
        travel, state, parents, nodes, context.tree_index,
        first_stats, iteration=context.spec.sweep_intervals,
        n_random=n_lineage_null,
    )

    travel_scaled = travel / first_stats["xyz_std"]
    state_scaled = state / first_stats["exp_std"]
    raw_front_travel = []
    raw_front_state = []
    assignments = []
    intervals = context.spec.sweep_intervals
    for step in range(intervals + 1):
        alpha = step / intervals
        rows, columns = linear_sum_assignment(
            alpha * travel_scaled + (1.0 - alpha) * state_scaled)
        assigned_row = np.empty_like(columns)
        assigned_row[columns] = rows
        assignments.append(assigned_row)
        raw_front_travel.append(float(travel[rows, columns].sum()))
        raw_front_state.append(float(state[rows, columns].sum()))

    null_raw: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    null_raw["first_cousin"] = (
        np.asarray(first_stats["random_xyz"], dtype=float),
        np.asarray(first_stats["random_exp"], dtype=float),
    )
    for label, steps, seed in (
        ("second_cousin", 3, context.spec.seed + 1),
        ("third_cousin", 4, context.spec.seed + 2),
    ):
        groups = pe.build_ancestor_groups(nodes, context.tree_index, steps)
        null_raw[label] = pe.compute_group_shuffle_costs(
            travel, state, groups,
            n_random=context.spec.null_draws, seed=seed)
    rng = np.random.default_rng(context.spec.seed + 3)
    rows = np.arange(len(nodes))
    full_x = np.empty(2000)
    full_y = np.empty(2000)
    for draw in range(2000):
        permutation = rng.permutation(len(nodes))
        full_x[draw] = travel[rows, permutation].sum()
        full_y[draw] = state[rows, permutation].sum()
    null_raw["full_random"] = (full_x, full_y)

    return GlobalAnalysisResult(
        context_cache_key=context.cache_key,
        analysis_cache_key=_analysis_cache_key(
            context, n_lineage_null=n_lineage_null),
        twr=twr,
        first_stats=first_stats,
        raw_front_travel=np.asarray(raw_front_travel),
        raw_front_state=np.asarray(raw_front_state),
        assignments=np.asarray(assignments, dtype=int),
        null_raw=null_raw,
    )


def _encode(value: Any, arrays: dict[str, np.ndarray], prefix: str) -> Any:
    if isinstance(value, np.ndarray):
        key = prefix.replace("/", "__")
        arrays[key] = value
        return {"__array__": key}
    if isinstance(value, dict):
        return {str(key): _encode(item, arrays, f"{prefix}/{key}")
                for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_encode(item, arrays, f"{prefix}/{index}")
                for index, item in enumerate(value)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"Cannot encode cache value {prefix}: {type(value).__name__}")


def _decode(value: Any, arrays: Any) -> Any:
    if isinstance(value, dict) and set(value) == {"__array__"}:
        return np.asarray(arrays[value["__array__"]])
    if isinstance(value, dict):
        return {key: _decode(item, arrays) for key, item in value.items()}
    if isinstance(value, list):
        return [_decode(item, arrays) for item in value]
    return value


def save_global_analysis(result: GlobalAnalysisResult,
                         context: AnalysisContext) -> tuple[Path, Path]:
    """Write arrays and JSON metadata to the profile's analysis directory."""
    if result.context_cache_key != context.cache_key:
        raise ValueError("Result belongs to a different analysis context")
    out = context.run_paths.create().analysis
    arrays = {
        "raw_front_travel": result.raw_front_travel,
        "raw_front_state": result.raw_front_state,
        "assignments": result.assignments,
    }
    tree = {
        "twr": _encode(result.twr, arrays, "twr"),
        "first_stats": _encode(result.first_stats, arrays, "first_stats"),
        "null_raw": _encode(result.null_raw, arrays, "null_raw"),
    }
    npz_path = out / "global_terminal_analysis.npz"
    metadata_path = out / "global_terminal_analysis.json"
    np.savez_compressed(npz_path, **arrays)
    metadata = {
        "version": GLOBAL_ANALYSIS_VERSION,
        "profile": context.profile,
        "context_cache_key": context.cache_key,
        "analysis_cache_key": result.analysis_cache_key,
        "tree": tree,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2,
                                        sort_keys=True) + "\n")
    return npz_path, metadata_path


def load_global_analysis(context: AnalysisContext) -> GlobalAnalysisResult:
    """Load only a complete cache matching this context and configuration."""
    out = context.run_paths.analysis
    npz_path = out / "global_terminal_analysis.npz"
    metadata_path = out / "global_terminal_analysis.json"
    if not npz_path.exists() or not metadata_path.exists():
        raise FileNotFoundError(
            f"Missing run-specific global cache under {out}; refusing legacy fallback")
    metadata = json.loads(metadata_path.read_text())
    expected = _analysis_cache_key(context, n_lineage_null=100)
    if metadata.get("version") != GLOBAL_ANALYSIS_VERSION:
        raise ValueError("Unsupported or incomplete global-analysis cache")
    if metadata.get("context_cache_key") != context.cache_key:
        raise ValueError("Global cache context does not match requested profile")
    if metadata.get("analysis_cache_key") != expected:
        raise ValueError("Global cache settings do not match requested analysis")
    with np.load(npz_path, allow_pickle=False) as arrays:
        tree = _decode(metadata["tree"], arrays)
        return GlobalAnalysisResult(
            context_cache_key=context.cache_key,
            analysis_cache_key=metadata["analysis_cache_key"],
            twr=tree["twr"],
            first_stats=tree["first_stats"],
            raw_front_travel=np.asarray(arrays["raw_front_travel"]),
            raw_front_state=np.asarray(arrays["raw_front_state"]),
            assignments=np.asarray(arrays["assignments"]),
            null_raw={
                key: (np.asarray(value[0]), np.asarray(value[1]))
                for key, value in tree["null_raw"].items()
            },
        )


def get_or_compute_global_analysis(context: AnalysisContext, *,
                                   force: bool = False) -> GlobalAnalysisResult:
    if not force:
        try:
            return load_global_analysis(context)
        except FileNotFoundError:
            pass
    result = compute_global_analysis(context)
    save_global_analysis(result, context)
    return result
