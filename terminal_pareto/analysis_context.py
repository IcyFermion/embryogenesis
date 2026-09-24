"""Profile-aware inputs and isolated output paths for terminal Pareto analyses.

This module is deliberately an adapter around the existing matrix-based
analysis engine.  It owns cohort selection, replicate aggregation, provenance,
and run paths; it does not implement a second optimizer.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe


MODULE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = MODULE_ROOT.parent
DEFAULT_OUTPUT_ROOT = MODULE_ROOT / "output"
CONTEXT_VERSION = "pooled-travel-context-v1"


@dataclass(frozen=True)
class TrackingInput:
    label: str
    path: str
    cutoff: int


@dataclass(frozen=True)
class AnalysisSpec:
    """Everything that gives a numerical terminal analysis its identity."""

    profile: str
    tracking: tuple[TrackingInput, ...]
    cohort_rule: str
    expression_representation: str = "top20_protein_cosine"
    travel_aggregation: str = "single_raw"
    travel_normalization: str = "none"
    null_family: str = "first_cousin_shuffle"
    scaling_policy: str = "legacy_sampled_1000"
    sweep_intervals: int = 300
    null_draws: int = 1000
    seed: int = 42


def _tracking_inputs() -> tuple[TrackingInput, ...]:
    return tuple(
        TrackingInput(label, str(Path(path).resolve()), int(cutoff))
        for label, path, cutoff in dl.CE_REPLICATES
    )


def get_analysis_spec(profile: str, *, sweep_intervals: int = 300,
                      null_draws: int = 1000, seed: int = 42) -> AnalysisSpec:
    """Return one of the three explicit migration profiles."""
    tracking = _tracking_inputs()
    common = dict(sweep_intervals=sweep_intervals, null_draws=null_draws,
                  seed=seed)
    if profile == "embryo1_legacy":
        return AnalysisSpec(
            profile=profile,
            tracking=(tracking[0],),
            cohort_rule="single_replicate_natural_edges",
            travel_aggregation="single_raw",
            travel_normalization="none",
            scaling_policy="legacy_sampled_1000",
            **common,
        )
    if profile == "embryo1_matched":
        return AnalysisSpec(
            profile=profile,
            tracking=tracking,
            cohort_rule="strict_natural_edge_intersection",
            travel_aggregation="embryo1_raw",
            travel_normalization="none",
            scaling_policy="exact_first_cousin_moments",
            **common,
        )
    if profile == "pooled_tracking_v1":
        return AnalysisSpec(
            profile=profile,
            tracking=tracking,
            cohort_rule="strict_natural_edge_intersection",
            travel_aggregation="equal_mean_normalized_components",
            travel_normalization="fixed_global_natural_total_per_replicate",
            scaling_policy="exact_first_cousin_moments",
            **common,
        )
    raise ValueError(
        f"Unknown analysis profile {profile!r}; expected embryo1_legacy, "
        "embryo1_matched, or pooled_tracking_v1"
    )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_array(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode())
    digest.update(str(array.shape).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def _stable_hash(payload: Mapping) -> str:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


@dataclass(frozen=True)
class RunPaths:
    """All numerical, figure, and validation paths for one isolated run."""

    output_root: Path
    profile: str
    run_id: str

    @property
    def root(self) -> Path:
        return self.output_root / "runs" / self.profile / self.run_id

    @property
    def analysis(self) -> Path:
        return self.root / "analysis"

    @property
    def publication(self) -> Path:
        return self.root / "publication"

    @property
    def validation(self) -> Path:
        return self.root / "validation"

    def display(self, mode: str) -> Path:
        if mode not in {"endpoint", "null_sd", "percent_natural"}:
            raise ValueError(f"Unsupported display mode: {mode}")
        return self.publication / mode

    def create(self) -> "RunPaths":
        for path in (self.analysis, self.publication, self.validation):
            path.mkdir(parents=True, exist_ok=True)
        return self


@dataclass
class AnalysisContext:
    """Ordered matrices plus biological/provenance metadata for one profile."""

    spec: AnalysisSpec
    lineage: dict
    tree_index: dict
    protein_exp: pd.DataFrame
    prot_sel: list[str]
    terminal_nodes: list[str]
    terminal_parents: list[str]
    travel_matrix: np.ndarray
    expression_matrix: np.ndarray
    component_travel_matrices: dict[str, np.ndarray]
    component_natural_totals: dict[str, float]
    component_xyz_maps: dict[str, dict[str, np.ndarray]]
    replicate_edge_sets: dict[str, set[tuple[str, str]]]
    reference_subtrees: dict[str, list[tuple[str, str]]]
    gp_map: dict
    source_hashes: dict[str, str]
    cache_key: str
    run_paths: RunPaths

    @property
    def profile(self) -> str:
        return self.spec.profile

    @property
    def ordered_edges(self) -> list[tuple[str, str]]:
        return list(zip(self.terminal_nodes, self.terminal_parents))

    @property
    def v_prot(self) -> list[str]:
        """Compatibility list for reference-tree traversal helpers."""
        names = []
        seen = set()
        for child, parent in self.ordered_edges:
            for name in (child, parent):
                if name not in seen:
                    names.append(name)
                    seen.add(name)
        return names

    @property
    def xyz_ce(self) -> dict[str, np.ndarray]:
        """Legacy coordinate-map adapter; invalid for pooled travel."""
        if self.profile == "pooled_tracking_v1":
            raise AttributeError(
                "Pooled travel has no single coordinate map; use cost_matrices()"
            )
        return self.component_xyz_maps["embryo1"]

    def edge_indices(self, terms: Sequence[tuple[str, str]]) -> np.ndarray:
        lookup = {edge: i for i, edge in enumerate(self.ordered_edges)}
        missing = [edge for edge in terms if edge not in lookup]
        if missing:
            raise KeyError(f"Edges are outside profile {self.profile}: {missing[:5]}")
        return np.asarray([lookup[edge] for edge in terms], dtype=int)

    def cost_matrices(self, terms: Sequence[tuple[str, str]] | None = None
                      ) -> tuple[np.ndarray, np.ndarray]:
        """Return global costs or an indexed submatrix with fixed normalizers."""
        if terms is None:
            return self.travel_matrix, self.expression_matrix
        indices = self.edge_indices(terms)
        selector = np.ix_(indices, indices)
        return self.travel_matrix[selector], self.expression_matrix[selector]

    def component_cost_matrices(
        self, terms: Sequence[tuple[str, str]] | None = None
    ) -> dict[str, np.ndarray]:
        if terms is None:
            return {key: value.copy()
                    for key, value in self.component_travel_matrices.items()}
        indices = self.edge_indices(terms)
        selector = np.ix_(indices, indices)
        return {key: value[selector]
                for key, value in self.component_travel_matrices.items()}

    def subtrees(self, min_cells: int = 12
                 ) -> list[tuple[str, list[tuple[str, str]]]]:
        cohort = set(self.ordered_edges)
        rows = []
        for name, reference_terms in self.reference_subtrees.items():
            descendant = set(reference_terms)
            terms = [edge for edge in self.ordered_edges
                     if edge in descendant and edge in cohort]
            if len(terms) >= min_cells:
                rows.append((name, terms))
        rows.sort(key=lambda item: (-len(item[1]), item[0]))
        return rows

    def manifest(self) -> dict:
        return {
            "context_version": CONTEXT_VERSION,
            "profile": self.profile,
            "run_id": self.run_paths.run_id,
            "cache_key": self.cache_key,
            "spec": {
                **asdict(self.spec),
                "tracking": [asdict(item) for item in self.spec.tracking],
            },
            "n_edges": len(self.terminal_nodes),
            "n_subtrees_min12": len(self.subtrees(12)),
            "ordered_edge_hash": _stable_hash({
                "edges": self.ordered_edges,
            }),
            "matrix_hashes": {
                "travel": _sha256_array(self.travel_matrix),
                "expression": _sha256_array(self.expression_matrix),
                **{
                    f"travel_component_{key}": _sha256_array(value)
                    for key, value in self.component_travel_matrices.items()
                },
            },
            "component_natural_totals": self.component_natural_totals,
            "source_hashes": self.source_hashes,
        }

    def write(self) -> None:
        """Write a complete, cache-keyed context fixture to the run directory."""
        paths = self.run_paths.create()
        arrays = {
            "travel_matrix": self.travel_matrix,
            "expression_matrix": self.expression_matrix,
            "children": np.asarray(self.terminal_nodes),
            "natural_parents": np.asarray(self.terminal_parents),
            "component_labels": np.asarray(
                list(self.component_travel_matrices), dtype=str),
            "component_natural_totals": np.asarray([
                self.component_natural_totals[label]
                for label in self.component_travel_matrices
            ]),
        }
        for label, matrix in self.component_travel_matrices.items():
            arrays[f"travel_component_{label}"] = matrix
        np.savez_compressed(paths.analysis / "analysis_context.npz", **arrays)
        (paths.analysis / "analysis_manifest.json").write_text(
            json.dumps(self.manifest(), indent=2, sort_keys=True) + "\n"
        )
        self.cell_manifest().to_csv(
            paths.analysis / "cell_manifest.csv", index=False)
        self.subtree_manifest().to_csv(
            paths.analysis / "subtree_manifest.csv", index=False)
    def cell_manifest(self) -> pd.DataFrame:
        union = set().union(*self.replicate_edge_sets.values())
        cohort = set(self.ordered_edges)
        rows = []
        for child, parent in sorted(union):
            row = {
                "terminal_cell": child,
                "natural_parent": parent,
                "included_profile": (child, parent) in cohort,
                "exclusion_reason": (
                    "" if (child, parent) in cohort
                    else "not_in_required_natural_edge_cohort"
                ),
            }
            for label, edge_set in self.replicate_edge_sets.items():
                row[f"present_{label}"] = (child, parent) in edge_set
            rows.append(row)
        return pd.DataFrame(rows)

    def subtree_manifest(self) -> pd.DataFrame:
        cohort = set(self.ordered_edges)
        rows = []
        for name, reference_terms in self.reference_subtrees.items():
            descendant = set(reference_terms)
            n = sum(edge in cohort for edge in descendant)
            rows.append({
                "subtree": name,
                "reference_terminal_edges": len(reference_terms),
                "n_profile": n,
                "eligible_min12": n >= 12,
            })
        return pd.DataFrame(rows).sort_values(
            ["n_profile", "subtree"], ascending=[False, True])


def validate_existing_context_manifest(
    context: AnalysisContext,
    analysis_out: Path | str | None = None,
) -> dict:
    """Validate an existing run identity before any provenance is rewritten.

    Downstream renderers call this before loading reusable tables and before
    ``context.write()``.  The cache key already covers the complete analysis
    specification, source hashes, ordered cohort, and cost matrices; the
    additional fields make failures easier to diagnose.
    """
    analysis_out = (context.run_paths.analysis if analysis_out is None
                    else Path(analysis_out))
    path = analysis_out / "analysis_manifest.json"
    if not path.exists():
        raise ValueError(
            f"Profile-specific caches lack the context manifest {path}; "
            "rerun the upstream analysis before rendering"
        )
    try:
        actual = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as error:
        raise ValueError(f"Unreadable analysis context manifest: {path}") from error

    expected = context.manifest()
    identity_fields = (
        "context_version",
        "profile",
        "cache_key",
        "n_edges",
        "ordered_edge_hash",
        "matrix_hashes",
        "source_hashes",
    )
    mismatched = {
        key: (actual.get(key), expected.get(key))
        for key in identity_fields
        if actual.get(key) != expected.get(key)
    }
    if mismatched:
        raise ValueError(f"Analysis context manifest mismatch: {mismatched}")
    return actual


def _load_replicates(spec: AnalysisSpec, lineage: dict,
                     protein_names: set[str]) -> tuple[
                         dict[str, dict[str, np.ndarray]],
                         dict[str, set[tuple[str, str]]],
                     ]:
    xyz_maps = {}
    edge_sets = {}
    for tracking in spec.tracking:
        xyz_map, valid = dl.load_elegans_tracking(
            tracking.cutoff, path=tracking.path)
        valid_expression = [name for name in valid if name in protein_names]
        children, parents = dl.collect_terminals(lineage, valid_expression)
        edges = set(zip(children, parents))
        if len(edges) != len(children):
            raise AssertionError(f"Duplicate natural edges in {tracking.label}")
        xyz_maps[tracking.label] = {
            name: np.asarray(value, dtype=float) for name, value in xyz_map.items()
        }
        edge_sets[tracking.label] = edges
    return xyz_maps, edge_sets


def build_analysis_context(
    profile: str = "embryo1_legacy", *, run_id: str | None = None,
    output_root: Path | str = DEFAULT_OUTPUT_ROOT, sweep_intervals: int = 300,
    null_draws: int = 1000, seed: int = 42,
) -> AnalysisContext:
    """Build production matrices directly from source data."""
    spec = get_analysis_spec(
        profile, sweep_intervals=sweep_intervals,
        null_draws=null_draws, seed=seed)
    lineage_path = REPO_ROOT / "data" / "cell_lineage.json"
    protein_path = (REPO_ROOT / "data" / "protein" / "aggregated_all"
                    / "s3_zscore.csv")
    selected_path = (REPO_ROOT / "expression_embedding" / "results"
                     / "elegans_protein_linear_baseline"
                     / "top20_protein_names.csv")
    lineage = dl.load_json(str(lineage_path))
    tree_index = lm.build_lineage_tree_index(lineage)
    protein = dl.load_protein_expression(str(protein_path))
    selected = dl.load_prot_sel(str(selected_path))
    xyz_maps, edge_sets = _load_replicates(spec, lineage, set(protein.index))

    if spec.cohort_rule == "single_replicate_natural_edges":
        cohort = set(edge_sets[spec.tracking[0].label])
    else:
        cohort = set.intersection(*(edge_sets[item.label]
                                    for item in spec.tracking))

    # Determine terminal eligibility by edge intersection, while preserving
    # full reference-tree traversal order and duplicate biological parent slots.
    ref_children, ref_parents = dl.collect_terminals(
        lineage, set(tree_index))
    ordered_edges = [edge for edge in zip(ref_children, ref_parents)
                     if edge in cohort]
    if set(ordered_edges) != cohort:
        missing = sorted(cohort - set(ordered_edges))
        raise AssertionError(
            f"Reference traversal did not recover profile edges: {missing[:5]}")
    terminal_nodes = [child for child, _parent in ordered_edges]
    terminal_parents = [parent for _child, parent in ordered_edges]

    component_matrices = {}
    expression_matrix = None
    for item in spec.tracking:
        matrix, expression, _ = pe.build_cost_matrices(
            terminal_nodes, terminal_parents, xyz_maps[item.label],
            protein, selected)
        component_matrices[item.label] = np.asarray(matrix, dtype=float)
        if expression_matrix is None:
            expression_matrix = np.asarray(expression, dtype=float)
        elif not np.array_equal(expression_matrix, expression):
            raise AssertionError("Expression matrix changed between tracking inputs")
    assert expression_matrix is not None

    natural_totals = {
        label: float(np.trace(matrix))
        for label, matrix in component_matrices.items()
    }
    if not all(value > 0 for value in natural_totals.values()):
        raise ValueError("All global natural travel totals must be positive")
    if spec.travel_aggregation in {"single_raw", "embryo1_raw"}:
        travel_matrix = component_matrices["embryo1"].copy()
    elif spec.travel_aggregation == "equal_mean_normalized_components":
        travel_matrix = np.mean([
            component_matrices[item.label] / natural_totals[item.label]
            for item in spec.tracking
        ], axis=0)
    else:
        raise ValueError(f"Unsupported aggregation: {spec.travel_aggregation}")

    all_subtrees = dict(dl.collect_all_subtrees(
        lineage, set(tree_index), min_cells=1))
    source_paths = [lineage_path, protein_path, selected_path] + [
        Path(item.path) for item in spec.tracking
    ]
    source_hashes = {
        str(path.relative_to(REPO_ROOT)): _sha256_file(path)
        for path in source_paths
    }
    identity = {
        "context_version": CONTEXT_VERSION,
        "spec": {
            **asdict(spec),
            "tracking": [asdict(item) for item in spec.tracking],
        },
        "source_hashes": source_hashes,
        "ordered_edges": ordered_edges,
        "travel_hash": _sha256_array(travel_matrix),
        "expression_hash": _sha256_array(expression_matrix),
        "component_hashes": {
            label: _sha256_array(matrix)
            for label, matrix in component_matrices.items()
        },
        "component_natural_totals": natural_totals,
    }
    cache_key = _stable_hash(identity)
    resolved_run_id = run_id or cache_key[:12]
    paths = RunPaths(Path(output_root).resolve(), profile, resolved_run_id)
    return AnalysisContext(
        spec=spec,
        lineage=lineage,
        tree_index=tree_index,
        protein_exp=protein,
        prot_sel=selected,
        terminal_nodes=terminal_nodes,
        terminal_parents=terminal_parents,
        travel_matrix=travel_matrix,
        expression_matrix=expression_matrix,
        component_travel_matrices=component_matrices,
        component_natural_totals=natural_totals,
        component_xyz_maps=xyz_maps,
        replicate_edge_sets=edge_sets,
        reference_subtrees=all_subtrees,
        gp_map=pe.build_grandparent_map(lineage),
        source_hashes=source_hashes,
        cache_key=cache_key,
        run_paths=paths,
    )


def write_hash_manifest(root: Path | str,
                        destination: Path | str | None = None) -> Path:
    """Hash an artifact tree without modifying any archived file."""
    root = Path(root).resolve()
    destination = (Path(destination).resolve() if destination is not None
                   else root.parent / f"{root.name}_sha256.json")
    files = sorted(path for path in root.rglob("*") if path.is_file()
                   and path.resolve() != destination)
    payload = {
        "root": str(root),
        "files": {
            str(path.relative_to(root)): _sha256_file(path) for path in files
        },
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return destination


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="pooled_tracking_v1",
                        choices=("embryo1_legacy", "embryo1_matched",
                                 "pooled_tracking_v1"))
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--hash-tree", type=Path)
    parser.add_argument("--hash-destination", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> AnalysisContext | Path:
    args = _parse_args(argv)
    if args.hash_tree is not None:
        path = write_hash_manifest(args.hash_tree, args.hash_destination)
        print(f"Wrote hash manifest to {path}")
        return path
    context = build_analysis_context(
        args.profile, run_id=args.run_id, output_root=args.output_root)
    context.write()
    print(json.dumps(context.manifest(), indent=2, sort_keys=True))
    return context


if __name__ == "__main__":
    main()
