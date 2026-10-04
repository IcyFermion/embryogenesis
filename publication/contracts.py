"""Minimal validated inputs consumed by publication layouts.

Adapters build these records from back-end caches after the back end's own
scientific validator has passed. Layouts read only these fields; they never
infer scope from a filename, figure number or row count.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


# Persisted metric columns -> notation concept keys used by layouts.
CANONICAL_FIELDS = {
    "u_L": "canonical_position",
    "d_LP": "natural_front_distance",
    "d_NP": "null_front_distance",
}
FRONT_COLUMNS = ("geometry", "config", "sweep_index", "D1", "D2", "edge_retention")
CLOUD_COLUMNS = ("geometry", "config", "family", "D1", "D2")
METRIC_COLUMNS = ("geometry", "config", "natural_D1", "natural_D2", "null_D1", "null_D2",
                  "nearest_index", "maximum_index", *CANONICAL_FIELDS.values())


@dataclass(frozen=True)
class FrontComparisonInput:
    """Base-cohort fronts and canonical metrics for one comparison scope.

    ``fronts`` rows are ordered by sweep index within each geometry/config;
    ``metrics`` has one row per geometry/config with canonical metrics under
    concept keys (see ``CANONICAL_FIELDS``).
    """

    family: str
    run: Path
    analysis_id: str
    edges: int
    configs: tuple[str, ...]
    geometries: tuple[str, ...]
    fronts: pd.DataFrame
    metrics: pd.DataFrame
    null_mean: str  # "analytic first-cousin" etc.; stated in captions, not inferred
    caption: dict = field(default_factory=dict)
    validation: dict = field(default_factory=dict)
    input_files: dict = field(default_factory=dict)
    # Reference clouds (geometry, config, family, D1, D2) and how the scope displays them:
    # main (families in main axes), inset, inset_title, inset_legend, display_draws (None = all).
    clouds: pd.DataFrame | None = None
    references: dict = field(default_factory=dict)
    # Scope-specific companions, e.g. terminal tracking transfers {config: [(label, D1, D2), ...]}.
    extras: dict = field(default_factory=dict)

    def __post_init__(self):
        for name, frame, columns in (("fronts", self.fronts, FRONT_COLUMNS),
                                     ("metrics", self.metrics, METRIC_COLUMNS)):
            missing = set(columns) - set(frame.columns)
            if missing:
                raise ValueError(f"{self.family} {name} missing {sorted(missing)}")
        if self.clouds is not None:
            missing = set(CLOUD_COLUMNS) - set(self.clouds.columns)
            if missing:
                raise ValueError(f"{self.family} clouds missing {sorted(missing)}")
            declared = set(self.references.get("main", ())) | {self.references.get("inset")} - {None}
            if declared - set(self.clouds.family):
                raise ValueError(f"{self.family} clouds lack declared reference families")
        expected = {(g, c) for g in self.geometries for c in self.configs}
        found = set(zip(self.metrics.geometry, self.metrics.config))
        if found != expected or len(self.metrics) != len(expected):
            raise ValueError(f"{self.family} metrics must have one row per geometry/config")

    def select(self, geometry: str, config: str):
        mask = (self.fronts.geometry == geometry) & (self.fronts.config == config)
        curve = self.fronts.loc[mask].sort_values("sweep_index")
        metric = self.metrics[(self.metrics.geometry == geometry) & (self.metrics.config == config)].iloc[0]
        return curve, metric

    def cloud(self, geometry: str, config: str):
        return self.clouds[(self.clouds.geometry == geometry) & (self.clouds.config == config)]


@dataclass(frozen=True)
class TerminalPrimaryInput:
    """Validated inputs for the primary 275-edge terminal figures (1-6, S1-S3) and Table 1.

    Each field is the cached result the corresponding layout consumes; the
    adapter has already run the back end's validators and prepared display
    coordinates (endpoint transforms) without recomputing any analysis.
    """

    run: Path
    analysis_id: str  # pooled analysis-context cache key
    edges: int
    subtrees: int
    min_cells: int
    global_twr: dict  # cached global front record (Figures 2-3)
    global_nulls: dict  # null-SD clouds keyed 1/2/3/"full"
    global_display: dict  # endpoint display coordinates for the global front
    canonical: pd.DataFrame  # validated canonical metrics (all endpoint-valid rows; Table 1)
    canonical_display: pd.DataFrame  # rows drawn in Figures 5/S1 (monotone front parameterization)
    subtree_summary: pd.DataFrame
    subtree_nodes: dict  # Figure 4 display nodes with validated attributes
    subtree_edges: list
    cell_type_retention: pd.DataFrame  # n<=4 groups merged (Figure 6A)
    cell_type_keypoints: pd.DataFrame  # n<=4 groups merged (Figure S2)
    within_type: dict  # fronts, summary, aggregate, endpoints (Figure 6B-C)
    s3_projections: pd.DataFrame
    s3_references: dict
    validation: dict = field(default_factory=dict)
    input_files: dict = field(default_factory=dict)


@dataclass(frozen=True)
class FullTreePooledInput:
    """Validated inputs for pooled full-tree Figure 7 and provisional supplement S4.

    978 measured nodes, 974 scored edges, Euclidean molecular distance. Fronts
    and references are replay-validated by the back end; endpoint transforms
    are per round (Figure 7A) and aggregate layerwise (Figure 7B and S4).
    """

    run: Path
    analysis_id: str
    edges: int
    leaves: int
    internal: int
    round_edges: tuple[int, ...]
    natural: tuple[float, float]  # aggregate natural (travel, state)
    fronts: pd.DataFrame  # method, scope, weight_index, travel, state, retention, nondominated
    layers: pd.DataFrame  # scope, edges, natural_travel, natural_state
    references: dict  # display name -> (draws, 2) raw costs
    methods: tuple[str, ...]  # first is the main layerwise method
    main_references: tuple[str, ...]
    supplement_references: tuple[str, ...]
    inset_reference: str
    round_transforms: dict  # scope -> endpoint transform (``.transform``, ``.metadata``)
    aggregate_transform: object
    settings: dict  # weights, draws, display_draws
    validation: dict = field(default_factory=dict)
    input_files: dict = field(default_factory=dict)
