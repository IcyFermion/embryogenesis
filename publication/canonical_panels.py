"""Canonical-metric replay checks and named-row metric panels.

The closest attained point P* is replayed from the cached front and checked
against the cached metrics before anything is drawn. Both distances use this
same P*, never an interpolated or maximum-retention point.
"""

from __future__ import annotations

import numpy as np

from publication import notation as nt
from publication import style

METRIC_KEYS = ("canonical_position", "natural_front_distance", "null_front_distance")


def canonical_point(curve, metric):
    """Return P* and its replayed canonical metrics; fail if the cache disagrees."""
    distances = np.hypot(curve.D1 - metric.natural_D1, curve.D2 - metric.natural_D2)
    nearest = int(np.argmin(distances.to_numpy()))
    if nearest != int(metric.nearest_index):
        raise ValueError("Cached nearest point does not match the attained front")
    closest = curve.iloc[nearest]
    arc = np.r_[0., np.cumsum(np.hypot(np.diff(curve.D1.to_numpy()[::-1]),
                                     np.diff(curve.D2.to_numpy()[::-1])))]
    if arc[-1] <= 0:
        raise ValueError("Canonical position requires a nondegenerate front")
    values = dict(canonical_position=float((arc / arc[-1])[::-1][nearest]),
                  natural_front_distance=float(distances.iloc[nearest]),
                  null_front_distance=float(np.hypot(closest.D1 - metric.null_D1, closest.D2 - metric.null_D2)))
    np.testing.assert_allclose([values[key] for key in METRIC_KEYS], [metric[key] for key in METRIC_KEYS],
                               rtol=1e-11, atol=1e-12)
    return closest, values


def canonical_records(data):
    records = []
    for geometry in data.geometries:
        for config in data.configs:
            curve, metric = data.select(geometry, config)
            closest, values = canonical_point(curve, metric)
            records.append(dict(geometry=geometry, config=config, **values,
                                nearest_index=int(metric.nearest_index),
                                maximum_index=int(metric.maximum_index),
                                closest_D1=float(closest.D1), closest_D2=float(closest.D2)))
    return records


def _column_max(values, step):
    return max(step, np.ceil(max(values) / step) * step)


def canonical_summary(fig, data, *, label="C"):
    """Figure 5B-style labelled rows: configurations grouped by tracking geometry."""
    records = canonical_records(data)
    per_group = len(data.configs)
    rows = np.array([1. + 5 * g + i for g in range(len(data.geometries)) for i in range(per_group)])
    fig.text(.07, .405, f"{label}   Canonical metrics at the closest attained assignment, {nt.math('closest_point')}",
             fontsize=10.5, weight="bold", va="top")
    columns = (
        (.30, "canonical_position", style.SEMANTIC_COLORS["canonical_u"], 1.),
        (.545, "natural_front_distance", style.SEMANTIC_COLORS["lineage_front_distance"],
         _column_max([r["natural_front_distance"] for r in records], .02)),
        (.79, "null_front_distance", style.SEMANTIC_COLORS["first_cousin_null"],
         _column_max([r["null_front_distance"] for r in records], .05)),
    )
    axes = []
    for left, key, color, xmax in columns:
        ax = fig.add_axes([left, .075, .18, .24])
        axes.append(ax)
        ax.set_ylim(rows[-1] + .65, -.55)
        ax.set_xlim(-.07 * xmax, 1.07 * xmax)
        ax.set_xticks([0., xmax / 2, xmax])
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False,
                       labelbottom=False, labelsize=8, length=2, pad=3)
        ax.set_yticks([])
        ax.set_title(nt.heading(key), fontsize=9.5, color=color, pad=21)
        ax.spines[["left", "right", "bottom"]].set_visible(False)
        ax.spines["top"].set_color("#CCCCCC")
        ax.grid(axis="x", color="#E8E8E8", lw=.55)
        values = [record[key] for record in records]
        ax.hlines(rows, 0, values, color=color, lw=.8, alpha=.45)
        points = ax.scatter(values, rows, color=color, s=30, edgecolor="white", linewidth=.4, zorder=3)
        points.set_gid(f"canonical:{key}")
    for g, geometry in enumerate(data.geometries):
        axes[0].text(-.08, 5. * g, f"{nt.GEOMETRY_LABELS[geometry]} tracking coordinates",
                     transform=axes[0].get_yaxis_transform(),
                     ha="right", va="center", fontsize=9.5, weight="bold", clip_on=False)
    for row, record in zip(rows, records):
        axes[0].text(-.08, row, nt.CONFIG_LABELS[record["config"]], transform=axes[0].get_yaxis_transform(),
                     ha="right", va="center", fontsize=9.5, clip_on=False)
    u = nt.symbol("canonical_position")
    fig.text(.07, .025, f"${u}=0$: travel optimum; ${u}=1$: cell-state optimum.  "
             f"Both distances use the same {nt.math('closest_point')}; "
             f"{nt.math('null_mean')} is the {data.null_mean} null mean.",
             fontsize=8.5, color="#444444")
    return axes
