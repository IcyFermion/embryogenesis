"""Single source of displayed notation for publication figures, captions and tables.

Concepts are keyed by semantic identity, never by their printed symbol. The
persisted scientific field names that carry each concept (for example
``d_NP``/``d_np`` in cached CSVs) are recorded for adapters and documentation
only; changing a displayed symbol here never renames a cache field.

Each symbol is stored once as a math body (no delimiters). Backends add their
own delimiters: ``math`` for Matplotlib mathtext, ``tex`` for LaTeX documents
and ``plain`` for log/CSV-adjacent text. Context-dependent scientific
qualifications (analytic vs sampled null means, cohort scope) belong in the
caption templates, not here.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class Concept:
    key: str
    symbol: str | None  # math body shared by mathtext and LaTeX, e.g. r"d_{LP}"
    plain: str | None  # readable ASCII fallback, e.g. "d_LP"
    label: str  # short readable name used before the symbol in headings
    definition: str
    precision: int | None = None
    persisted_fields: tuple[str, ...] = field(default_factory=tuple)


_DEFAULTS = (
    Concept("travel_axis", "D_1", "D1", "Travel distance",
            "Travel cost normalized by the endpoint cost span; travel optimum at 0",
            persisted_fields=("D1",)),
    Concept("cell_state_axis", "D_2", "D2", "Cell-state distance",
            "Cell-state cost normalized by the endpoint cost span; cell-state optimum at 0",
            persisted_fields=("D2",)),
    Concept("travel_optimum", "T", "T", "Travel optimum",
            "Attained assignment minimizing travel; endpoint coordinates (0, 1)"),
    Concept("cell_state_optimum", "S", "S", "Cell-state optimum",
            "Attained assignment minimizing cell-state cost; endpoint coordinates (1, 0)"),
    Concept("travel_cost", "T", "T", "Travel cost",
            "Total travel cost of an assignment before endpoint normalization"),
    Concept("cell_state_cost", "S", "S", "Cell-state cost",
            "Total cell-state cost of an assignment before endpoint normalization"),
    Concept("natural_lineage", "L", "L", "Natural lineage",
            "The observed natural assignment in endpoint coordinates",
            persisted_fields=("natural_D1", "natural_D2")),
    Concept("null_mean", "C", "C", "First-cousin null mean",
            "Mean endpoint coordinates of the first-cousin shuffle null",
            persisted_fields=("null_D1", "null_D2")),
    Concept("closest_point", "P^*", "P*", "Closest attained point",
            "Attained sampled front assignment closest to L in endpoint coordinates",
            persisted_fields=("nearest_index",)),
    Concept("canonical_position", "u", "u", "Position",
            "Front arc length from the travel optimum to P*, divided by total arc length",
            precision=3, persisted_fields=("u_L", "u_lineage_lp")),
    Concept("natural_front_distance", "d_{LP}", "d_LP", "Natural distance",
            "Euclidean endpoint-coordinate distance from L to P*",
            precision=4, persisted_fields=("d_LP", "d_lp")),
    Concept("null_front_distance", "d_{CP}", "d_CP", "Null-mean distance",
            "Euclidean endpoint-coordinate distance from the first-cousin null mean to the same P*",
            precision=4, persisted_fields=("d_NP", "d_np")),
    Concept("retention", None, None, "Edge retention",
            "Fraction of edges whose biological parent matches natural lineage",
            precision=3, persisted_fields=("edge_retention",)),
    Concept("cousin_relative", "r", "r", "Cousin-relative statistic",
            "Separately null-standardized sensitivity statistic; not d_LP divided by the null distance",
            precision=3),
)

_REGISTRY: dict[str, Concept] = {concept.key: concept for concept in _DEFAULTS}


def concept(key: str) -> Concept:
    try:
        return _REGISTRY[key]
    except KeyError:
        raise KeyError(f"Unknown notation concept: {key}") from None


def symbol(key: str) -> str:
    value = concept(key).symbol
    if value is None:
        raise ValueError(f"{key} is displayed by label only")
    return value


def math(key: str) -> str:
    """Matplotlib mathtext, e.g. ``$d_{LP}$``."""
    return f"${symbol(key)}$"


def tex(key: str) -> str:
    r"""Inline LaTeX, e.g. ``\(d_{LP}\)``."""
    return rf"\({symbol(key)}\)"


def plain(key: str) -> str:
    return concept(key).plain


def label(key: str) -> str:
    return concept(key).label


def heading(key: str) -> str:
    """Readable label followed by its mathtext symbol, e.g. ``Position, $u$``."""
    return f"{label(key)}, {math(key)}"


def snapshot() -> dict[str, dict]:
    """Serializable registry state for presentation provenance."""
    return {key: dict(symbol=c.symbol, plain=c.plain, label=c.label) for key, c in _REGISTRY.items()}


@contextmanager
def override(**changes: dict):
    """Temporarily replace concept fields, e.g. ``override(null_front_distance=dict(symbol="x"))``.

    Intended for propagation tests; the committed registry above is the
    publication notation.
    """
    saved = dict(_REGISTRY)
    try:
        for key, fields in changes.items():
            _REGISTRY[key] = replace(concept(key), **fields)
        yield
    finally:
        _REGISTRY.clear()
        _REGISTRY.update(saved)


# Display names for configuration/geometry identifiers shared by both back ends.
CONFIG_LABELS = {"ce_protein": "C. elegans protein", "ce_rna": "C. elegans RNA", "cb_rna": "C. briggsae RNA"}
GEOMETRY_LABELS = {"raw3d": "3D", "xy": "2D (XY)"}
