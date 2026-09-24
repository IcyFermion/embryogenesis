"""Pure display-coordinate transforms for saved terminal Pareto results."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping, Sequence

import numpy as np


TRANSFORM_VERSION = "endpoint-affine-v1"


class DegenerateEndpointSpan(ValueError):
    """Raised when endpoint coordinates do not define two positive spans."""


@dataclass(frozen=True)
class EndpointTransform:
    """Affine map sending travel optimum A to (0,1), state optimum B to (1,0)."""

    reference_analysis_id: str
    travel_optimum_assignment_id: str
    state_optimum_assignment_id: str
    travel_at_travel_optimum: float
    state_at_travel_optimum: float
    travel_at_state_optimum: float
    state_at_state_optimum: float
    travel_span: float
    state_span: float
    tolerance: float
    version: str = TRANSFORM_VERSION

    @classmethod
    def from_endpoints(
        cls, *, reference_analysis_id: str,
        travel_optimum_assignment_id: str,
        state_optimum_assignment_id: str,
        travel_optimum_costs: Sequence[float],
        state_optimum_costs: Sequence[float],
        rtol: float = 1e-12, atol: float = 1e-15,
    ) -> "EndpointTransform":
        ta, sa = map(float, travel_optimum_costs)
        tb, sb = map(float, state_optimum_costs)
        travel_span = tb - ta
        state_span = sa - sb
        scale = max(abs(ta), abs(sa), abs(tb), abs(sb), 1.0)
        tolerance = max(float(atol), float(rtol) * scale)
        if not (travel_span > tolerance and state_span > tolerance):
            raise DegenerateEndpointSpan(
                "Endpoint normalization requires positive, non-degenerate "
                f"spans; travel_span={travel_span:.12g}, "
                f"state_span={state_span:.12g}, tolerance={tolerance:.12g}. "
                "Use a clearly labeled null-SD or raw-cost view instead."
            )
        return cls(
            reference_analysis_id=reference_analysis_id,
            travel_optimum_assignment_id=str(travel_optimum_assignment_id),
            state_optimum_assignment_id=str(state_optimum_assignment_id),
            travel_at_travel_optimum=ta,
            state_at_travel_optimum=sa,
            travel_at_state_optimum=tb,
            state_at_state_optimum=sb,
            travel_span=travel_span,
            state_span=state_span,
            tolerance=tolerance,
        )

    def transform(self, travel, state) -> tuple[np.ndarray, np.ndarray]:
        travel = np.asarray(travel, dtype=float)
        state = np.asarray(state, dtype=float)
        if travel.shape != state.shape:
            raise ValueError("Travel and state arrays must have identical shapes")
        return (
            (travel - self.travel_at_travel_optimum) / self.travel_span,
            (state - self.state_at_state_optimum) / self.state_span,
        )

    def inverse(self, x, y) -> tuple[np.ndarray, np.ndarray]:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.shape != y.shape:
            raise ValueError("Endpoint-coordinate arrays must have identical shapes")
        return (
            self.travel_at_travel_optimum + x * self.travel_span,
            self.state_at_state_optimum + y * self.state_span,
        )

    def metadata(self, *, clipping: bool = False,
                 axis_limits: Mapping[str, Sequence[float]] | None = None) -> dict:
        payload = asdict(self)
        payload.update({
            "display_mode": "endpoint",
            "clipping": bool(clipping),
            "axis_limits": dict(axis_limits or {}),
        })
        return payload


def null_sd_coordinates(travel, state, *, natural_costs: Sequence[float],
                        null_stds: Sequence[float]
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Natural-centred coordinates in first-cousin-null SD units."""
    travel = np.asarray(travel, dtype=float)
    state = np.asarray(state, dtype=float)
    natural_travel, natural_state = map(float, natural_costs)
    travel_sd, state_sd = map(float, null_stds)
    if not (travel_sd > 0 and state_sd > 0):
        raise ValueError("Null standard deviations must be positive")
    return ((travel - natural_travel) / travel_sd,
            (state - natural_state) / state_sd)


def percent_natural_coordinates(travel, state, *, natural_costs: Sequence[float]
                                ) -> tuple[np.ndarray, np.ndarray]:
    """Percentage cost change relative to the natural lineage (nature = 0)."""
    travel = np.asarray(travel, dtype=float)
    state = np.asarray(state, dtype=float)
    natural_travel, natural_state = map(float, natural_costs)
    if not (natural_travel > 0 and natural_state > 0):
        raise ValueError("Natural costs must be positive")
    return (100.0 * (travel / natural_travel - 1.0),
            100.0 * (state / natural_state - 1.0))


def build_endpoint_transform(
    travel, state, *, travel_optimum_index: int, state_optimum_index: int,
    reference_analysis_id: str, assignment_ids: Sequence[str] | None = None,
) -> EndpointTransform:
    """Fit the declared endpoint map from one saved analysis result."""
    travel = np.asarray(travel, dtype=float)
    state = np.asarray(state, dtype=float)
    if travel.shape != state.shape or travel.ndim != 1:
        raise ValueError("Front costs must be same-length one-dimensional arrays")
    if assignment_ids is None:
        assignment_ids = [f"sweep:{index}" for index in range(len(travel))]
    if len(assignment_ids) != len(travel):
        raise ValueError("assignment_ids length does not match front arrays")
    return EndpointTransform.from_endpoints(
        reference_analysis_id=reference_analysis_id,
        travel_optimum_assignment_id=assignment_ids[travel_optimum_index],
        state_optimum_assignment_id=assignment_ids[state_optimum_index],
        travel_optimum_costs=(travel[travel_optimum_index],
                              state[travel_optimum_index]),
        state_optimum_costs=(travel[state_optimum_index],
                             state[state_optimum_index]),
    )
