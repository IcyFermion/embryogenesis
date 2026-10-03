"""Figure 2 style for the pooled terminal molecular/species comparison."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from terminal_pareto import plot_style as ps
from terminal_pareto.cross_species_analysis import LABELS, PRIMARY, file_hash
from terminal_pareto.fig2_fig3_ce_terminal_pareto import EDGE_RETENTION_CMAP, proportional_limits
from terminal_pareto.front_coordinates import EndpointTransform


TRACKING_CAVEAT = (
    "All C. briggsae 3D panels are subject to limitations of traditional embryo-tracking techniques, "
    "particularly for z-axis measurements."
)
TRACKING_CAPTION_STATUS = (
    "This description is provisional and will be refined with experimental collaborators."
)


def save(fig, out, stem):
    fig.savefig(out / f"{stem}.png", dpi=300, facecolor="white")
    fig.savefig(out / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def select(fronts, metrics, geometry, config):
    mask = (fronts.geometry == geometry) & (fronts.config == config) & (fronts.cohort == "base")
    curve = fronts.loc[mask].sort_values("sweep_index")
    metric = metrics[(metrics.geometry == geometry) & (metrics.config == config) & (metrics.cohort == "base")].iloc[0]
    return curve, metric


def base_axes(ax):
    ax.axhline(0, color="#777777", lw=.7, ls=":", zorder=0)
    ax.axvline(0, color="#777777", lw=.7, ls=":", zorder=0)
    ax.grid(True, alpha=.32)
    ax.spines[["top", "right"]].set_visible(True)


def comparison_heading(ax, config, geometry, panel):
    """Keep configuration first and make tracking geometry a clear second line."""
    ax.set_title(LABELS[config], fontsize=11, pad=25)
    geometry_label = "3D tracking" if geometry == "raw3d" else "2D (XY) tracking"
    subtitle = ax.text(.5, 1.025, geometry_label, transform=ax.transAxes,
                       fontsize=11, fontweight="bold", color=ps.COLORS["black"],
                       ha="center", va="bottom")
    subtitle.set_gid(f"comparison_geometry:{geometry}")
    ax.text(-.025, 1.025, panel, transform=ax.transAxes,
            fontsize=10, fontweight="bold", ha="left", va="bottom")


def comparison(fronts, metrics, clouds, geometry, out, n):
    """Compose one row or the primary two-row 3D/XY comparison."""
    geometries = ("raw3d", "xy") if geometry is None else (geometry,)
    combined = len(geometries) == 2
    fig, axes = plt.subplots(len(geometries), 3, figsize=(14.4, 9.8 if combined else 5.1),
                             sharex=True, sharey=True, squeeze=False)
    if combined:
        fig.subplots_adjust(left=.06, right=.91, bottom=.16, top=.87, hspace=.32, wspace=.13)
    else:
        fig.subplots_adjust(left=.06, right=.91, bottom=.24, top=.80, wspace=.13)
    fig.suptitle(f"Terminal-only comparison | {n} shared terminal edges", y=.97, fontsize=12)
    visible_x, visible_y = [], []
    panels = [(row, col, ax, g, config) for row, g in enumerate(geometries)
              for col, (ax, config) in enumerate(zip(axes[row], PRIMARY))]
    for row, col, ax, g, config in panels:
        base_axes(ax)
        curve, m = select(fronts, metrics, g, config)
        null = clouds[(clouds.geometry == g) & (clouds.config == config)]
        for family in ("first_cousin", "second_cousin", "third_cousin"):
            sample = null[null.family == family]
            color = ps.NULL_MODEL_COLORS[family]
            ax.scatter(sample.D1, sample.D2, s=7, color=color, alpha=.19,
                       edgecolors="none", rasterized=True, zorder=1)
            # First-cousin mean is analytic, matching the canonical metric.
            mean = (m.null_D1, m.null_D2) if family == "first_cousin" else (sample.D1.mean(), sample.D2.mean())
            ax.scatter(*mean, marker="+", s=34, color=color, lw=1.1, zorder=4)
            visible_x.extend(sample.D1.to_list())
            visible_y.extend(sample.D2.to_list())
        ax.plot(curve.D1, curve.D2, color=ps.COLORS["blue"], lw=1.3, alpha=.5, zorder=2)
        front = ax.scatter(curve.D1, curve.D2, c=curve.edge_retention, cmap=EDGE_RETENTION_CMAP,
                           norm=Normalize(0, 1), s=18, edgecolors="none", zorder=3)
        ax.scatter(m.natural_D1, m.natural_D2, marker="X", s=78, color="#222222",
                   edgecolor="white", lw=.6, zorder=7)
        maximum = curve.iloc[int(m.maximum_index)]
        ax.scatter(maximum.D1, maximum.D2, s=32, facecolor=EDGE_RETENTION_CMAP(maximum.edge_retention),
                   edgecolor="#222222", lw=1.0, zorder=7)
        visible_x.extend(curve.D1.to_list() + [m.natural_D1])
        visible_y.extend(curve.D2.to_list() + [m.natural_D2])
        full = null[null.family == "full_random"]
        inset = ax.inset_axes([.74, .69, .23, .26])
        inset.scatter(full.D1, full.D2, s=4, color=ps.NULL_MODEL_COLORS["full_random"],
                      alpha=.18, edgecolors="none", rasterized=True)
        inset.scatter(full.D1.mean(), full.D2.mean(), marker="+", s=26,
                      color=ps.NULL_MODEL_COLORS["full_random"], lw=1.0)
        inset.set_xlim(*proportional_limits(full.D1))
        inset.set_ylim(*proportional_limits(full.D2))
        inset.set_title("Full-random shuffle", fontsize=7, pad=2)
        inset.set_xlabel("Travel distance", fontsize=6, labelpad=1)
        inset.set_ylabel("Cell-state distance", fontsize=6, labelpad=1)
        inset.tick_params(labelsize=5.5, length=2, pad=1)
        inset.locator_params(axis="both", nbins=3)
        inset.grid(True, alpha=.25)
        inset.spines[["top", "right"]].set_visible(True)
        comparison_heading(ax, config, g, "ABCDEF"[3 * row + col] if combined else "ABC"[col])
        if row == len(geometries) - 1:
            ax.set_xlabel("Travel distance\n(fraction of endpoint cost span)")
        if col == 0:
            ax.set_ylabel("Cell-state distance\n(fraction of endpoint cost span)")
    for values, setter in ((visible_x, axes[0, 0].set_xlim), (visible_y, axes[0, 0].set_ylim)):
        lo, hi = min(values), max(values)
        pad = max(.04 * (hi - lo), .04)
        setter(lo - pad, hi + pad)
    color_axis = fig.add_axes([.934, .16 if combined else .24, .013, .71 if combined else .56])
    bar = fig.colorbar(front, cax=color_axis)
    bar.set_label("Edge retention")
    handles = [Line2D([], [], marker="o", ls="", color=ps.NULL_MODEL_COLORS[family], alpha=.7,
                      markersize=5, markeredgecolor="none", label=label) for family, label in (
        ("first_cousin", "First-cousin shuffle"), ("second_cousin", "Second-cousin shuffle"),
        ("third_cousin", "Third-cousin shuffle"), ("full_random", "Full-random shuffle (insets)"))]
    handles.extend([Line2D([], [], marker="X", ls="", color="#222222", markersize=7, label="Natural lineage"),
                    Line2D([], [], marker="o", ls="", color=ps.COLORS["blue"], markeredgecolor="#222222",
                           markersize=5, label="Maximum edge retention")])
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.49, .025), ncol=3, fontsize=8)
    stem = ("terminal_cross_species_comparison" if combined else
            "terminal_cross_species_existing_3d" if geometry == "raw3d" else "terminal_cross_species_xy")
    save(fig, out, stem)


def overlay_panel(ax, fronts, metrics, geometry, *, zoom_limits=None, legend=True):
    """Draw one geometry; optional common zoom limits support 3D/XY parity."""
    base_axes(ax)
    handles = []
    styles = ("-", "--", "-.")
    visible_x, visible_y = [], []
    for config, style in zip(PRIMARY, styles):
        curve, m = select(fronts, metrics, geometry, config)
        color = ps.SPECIES_COLORS[config]
        ax.plot(curve.D1, curve.D2, color=color, lw=1.6, ls=style)
        ax.scatter(m.natural_D1, m.natural_D2, marker="X", s=70, color=color, edgecolor="white", lw=.6, zorder=7)
        maximum = curve.iloc[int(m.maximum_index)]
        ax.scatter(maximum.D1, maximum.D2, s=34, color=color, edgecolor="#222222", lw=1, zorder=6)
        closest, _ = canonical_point(curve, m)
        connector, = ax.plot([m.natural_D1, closest.D1], [m.natural_D2, closest.D2],
                             color=color, lw=.9, ls=(0, (2, 2)), zorder=5)
        connector.set_gid(f"natural_to_closest:{geometry}:{config}")
        ax.scatter(closest.D1, closest.D2, marker="D", s=30, facecolor="white",
                   edgecolor=color, lw=1.1, zorder=8)
        handles.append(Line2D([], [], color=color, ls=style, label=LABELS[config]))
        visible_x.extend(curve.D1.to_list() + [m.natural_D1])
        visible_y.extend(curve.D2.to_list() + [m.natural_D2])
    ax.set(xlabel="Travel distance\n(fraction of endpoint cost span)",
           ylabel="Cell-state distance\n(fraction of endpoint cost span)",
           xlim=(min(visible_x) - .04, max(visible_x) + .04),
           ylim=(min(visible_y) - .04, max(visible_y) + .04))
    # Zoom shares the exact coordinates and anchors of the main axes.
    inset = ax.inset_axes([.48, .40, .47, .47])
    natural_x, natural_y = [], []
    for config, style in zip(PRIMARY, styles):
        curve, m = select(fronts, metrics, geometry, config)
        color = ps.SPECIES_COLORS[config]
        inset.plot(curve.D1, curve.D2, color=color, ls=style, lw=1.4)
        inset.scatter(m.natural_D1, m.natural_D2, marker="X", s=45, color=color, edgecolor="white", lw=.5, zorder=4)
        closest, _ = canonical_point(curve, m)
        connector, = inset.plot([m.natural_D1, closest.D1], [m.natural_D2, closest.D2],
                                color=color, lw=1, ls=(0, (2, 2)), zorder=3)
        connector.set_gid(f"natural_to_closest:{geometry}:{config}")
        inset.scatter(closest.D1, closest.D2, marker="D", s=25, facecolor="white",
                      edgecolor=color, lw=1, zorder=5)
        natural_x.append(m.natural_D1)
        natural_y.append(m.natural_D2)
    inset.set_xlim(*(zoom_limits[0] if zoom_limits else (min(natural_x) - .035, max(natural_x) + .035)))
    inset.set_ylim(*(zoom_limits[1] if zoom_limits else (min(natural_y) - .08, max(natural_y) + .08)))
    inset.set_title("Near natural lineages", fontsize=7.5)
    inset.tick_params(labelsize=6)
    inset.grid(True, alpha=.3)
    inset.spines[["top", "right"]].set_visible(True)
    handles.extend([Line2D([], [], marker="X", ls="", color="#222222", label="Natural lineage"),
                    Line2D([], [], marker="o", ls="", color="white", markeredgecolor="#222222", label="Maximum edge retention"),
                    Line2D([], [], marker="D", ls=(0, (2, 2)), color="#666666", markerfacecolor="white",
                           label=r"Closest attained point, $P^*$")])
    if legend:
        ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(.98, .32), fontsize=7)
    return handles


def overlay(fronts, metrics, out, n, geometry):
    fig, ax = plt.subplots(figsize=(7.15, 4.85))
    fig.subplots_adjust(left=.12, right=.96, bottom=.17, top=.86)
    overlay_panel(ax, fronts, metrics, geometry)
    fig.suptitle(f"Pooled terminal fronts | {n} shared edges", y=.97, fontsize=10.5)
    subtitle = ("3D: pooled distances" if geometry == "raw3d" else
                "2D (XY): axial coordinates omitted")
    fig.text(.5, .91, subtitle, ha="center", fontsize=7.5, color="#555555")
    suffix = "existing_3d" if geometry == "raw3d" else "xy"
    save(fig, out, f"terminal_cross_species_overlay_{suffix}")


def canonical_point(curve, metric):
    """Replay the Figure 1 metrics at an attained point, never at max ER."""
    distances = np.hypot(curve.D1 - metric.natural_D1, curve.D2 - metric.natural_D2)
    nearest = int(np.argmin(distances.to_numpy()))
    if nearest != int(metric.nearest_index):
        raise ValueError("Cached nearest point does not match the attained front")
    closest = curve.iloc[nearest]
    arc = np.r_[0., np.cumsum(np.hypot(np.diff(curve.D1.to_numpy()[::-1]),
                                     np.diff(curve.D2.to_numpy()[::-1])))]
    if arc[-1] <= 0:
        raise ValueError("Canonical position requires a nondegenerate front")
    values = dict(u_L=float((arc / arc[-1])[::-1][nearest]),
                  d_LP=float(distances.iloc[nearest]),
                  d_NP=float(np.hypot(closest.D1 - metric.null_D1, closest.D2 - metric.null_D2)))
    np.testing.assert_allclose([values[key] for key in values],
                               [metric[key] for key in values], rtol=1e-11, atol=1e-12)
    return closest, values


def canonical_records(fronts, metrics):
    records = []
    for geometry in ("raw3d", "xy"):
        for config in PRIMARY:
            curve, metric = select(fronts, metrics, geometry, config)
            closest, values = canonical_point(curve, metric)
            records.append(dict(geometry=geometry, config=config, **values,
                                nearest_index=int(metric.nearest_index),
                                maximum_index=int(metric.maximum_index),
                                closest_D1=float(closest.D1), closest_D2=float(closest.D2)))
    return records


def canonical_summary(fig, fronts, metrics):
    """Figure 5B-style labelled rows, grouped by 3D and XY geometry."""
    records = canonical_records(fronts, metrics)
    rows = np.array([1., 2., 3., 6., 7., 8.])
    fig.text(.07, .405, r"C   Canonical metrics at the closest attained assignment, $P^*$",
             fontsize=10.5, weight="bold", va="top")
    axes = []
    for left, field, title, color, xmax in (
        (.30, "u_L", r"Position, $u$", ps.SEMANTIC_COLORS["canonical_u"], 1.),
        (.545, "d_LP", r"Natural distance, $d_{LP}$", ps.SEMANTIC_COLORS["lineage_front_distance"],
         max(.02, np.ceil(max(row["d_LP"] for row in records) / .02) * .02)),
        (.79, "d_NP", r"Null-mean distance, $d_{NP}$", ps.SEMANTIC_COLORS["first_cousin_null"],
         max(.05, np.ceil(max(row["d_NP"] for row in records) / .05) * .05)),
    ):
        ax = fig.add_axes([left, .075, .18, .24])
        axes.append(ax)
        ax.set_ylim(8.65, -.55)
        ax.set_xlim(-.07 * xmax, 1.07 * xmax)
        ax.set_xticks([0., xmax / 2, xmax])
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False,
                       labelbottom=False, labelsize=8, length=2, pad=3)
        ax.set_yticks([])
        ax.set_title(title, fontsize=9.5, color=color, pad=21)
        ax.spines[["left", "right", "bottom"]].set_visible(False)
        ax.spines["top"].set_color("#CCCCCC")
        ax.grid(axis="x", color="#E8E8E8", lw=.55)
        values = [record[field] for record in records]
        ax.hlines(rows, 0, values, color=color, lw=.8, alpha=.45)
        points = ax.scatter(values, rows, color=color, s=30, edgecolor="white", linewidth=.4, zorder=3)
        points.set_gid(f"canonical:{field}")
    for y, label in ((0., "3D tracking coordinates"), (5., "2D (XY) tracking coordinates")):
        axes[0].text(-.08, y, label, transform=axes[0].get_yaxis_transform(),
                     ha="right", va="center", fontsize=9.5, weight="bold", clip_on=False)
    for row, record in zip(rows, records):
        axes[0].text(-.08, row, LABELS[record["config"]], transform=axes[0].get_yaxis_transform(),
                     ha="right", va="center", fontsize=9.5, clip_on=False)
    fig.text(.07, .025, r"$u=0$: travel optimum; $u=1$: cell-state optimum.  "
             r"Both distances use the same $P^*$; $N$ is the analytic first-cousin null mean.",
             fontsize=8.5, color="#444444")
    return axes


def overlay_comparison(fronts, metrics, out, n):
    fig = plt.figure(figsize=(12.8, 8.8))
    first = fig.add_axes([.07, .59, .42, .30])
    axes = [first, fig.add_axes([.56, .59, .42, .30], sharex=first, sharey=first)]
    primary = metrics[(metrics.cohort == "base") & metrics.config.isin(PRIMARY)]
    zoom = ((primary.natural_D1.min() - .035, primary.natural_D1.max() + .035),
            (primary.natural_D2.min() - .08, primary.natural_D2.max() + .08))
    for ax, geometry, title in zip(axes, ("raw3d", "xy"), ("A   3D", "B   2D (XY)")):
        handles = overlay_panel(ax, fronts, metrics, geometry, zoom_limits=zoom, legend=False)
        ax.set_title(title, fontsize=10, pad=12)
    axes[1].set_ylabel("")
    fig.suptitle(f"Pooled terminal fronts and canonical metrics | {n} shared edges", y=.985, fontsize=12)
    fig.text(.5, .947, "Overlay main axes and near-natural zoom limits shared",
             ha="center", fontsize=8, color="#555555")
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.52, .465), ncol=3, fontsize=8)
    canonical_summary(fig, fronts, metrics)
    save(fig, out, "terminal_cross_species_overlay")


def tracking_comparison(fronts, metrics, record, arrays, out):
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.3), sharex=True, sharey=True)
    fig.subplots_adjust(left=.06, right=.98, top=.77, bottom=.28, wspace=.16)
    children = np.arange(record["cohort_size"])
    colors = [ps.COLORS["blue"], ps.COLORS["orange"], ps.COLORS["green"]]
    all_x, all_y = [], []
    for i, (ax, config) in enumerate(zip(axes, PRIMARY)):
        base_axes(ax)
        curve, m = select(fronts, metrics, "raw3d", config)
        ax.plot(curve.D1, curve.D2, color="#222222", lw=1.8, label="Pooled optimum")
        ax.scatter(m.natural_D1, m.natural_D2, color="#222222", marker="X", s=55, zorder=5)
        all_x.extend(curve.D1.to_list() + [m.natural_D1])
        all_y.extend(curve.D2.to_list() + [m.natural_D2])
        key = f"raw3d__{config}"
        meta = record["endpoint_transforms"][key]
        transform = EndpointTransform(**{k: v for k, v in meta.items() if k not in ("display_mode", "clipping", "axis_limits")})
        species = "cb" if config == "cb_rna" else "ce"
        for j, replica in enumerate(record["tracking"][species]):
            perms = arrays[f"{key}__{replica['label']}__permutations"]
            t = arrays[f"{key}__travel"][children, perms].sum(axis=1)
            s = arrays[f"{key}__state"][children, perms].sum(axis=1)
            x, y = transform.transform(t, s)
            all_x.extend(x)
            all_y.extend(y)
            label = f"AF16 {'p1' if j == 0 else 'p5'}" if species == "cb" else f"Embryo {j + 1}"
            ax.plot(x, y, color=colors[j], lw=1.2, ls="--", alpha=.85, label=label)
        ax.set_title(f"{'ABC'[i]}   {LABELS[config]}", fontsize=9, pad=10)
        ax.set_xlabel("Travel distance\n(pooled endpoint coordinates)")
        ax.legend(loc="upper right", bbox_to_anchor=(.98, .65), fontsize=7)
    axes[0].set_ylabel("Cell-state distance\n(pooled endpoint coordinates)")
    axes[0].set_xlim(min(all_x) - .04, max(all_x) + .04)
    axes[0].set_ylim(min(all_y) - .04, max(all_y) + .04)
    fig.suptitle("Tracking sensitivity: all assignments evaluated using pooled travel", fontsize=11, y=.95)
    fig.text(.5, .87, f"{record['cohort_size']} shared edges | 3D pooled distances", ha="center", fontsize=8, color="#555555")
    fig.text(.5, .06, "Solid: pooled optimization. Dashed: individual-embryo optimization transferred to the pooled objective.\n"
             "Every panel uses its pooled front's anchors; molecular matrices and biological-parent capacities stay fixed.",
             ha="center", fontsize=8)
    save(fig, out, "terminal_cross_species_tracking_sensitivity")


def render(run):
    ps.configure()
    analysis, out = run / "analysis", run / "figures"
    out.mkdir(parents=True, exist_ok=True)
    record = json.loads((analysis / "provenance.json").read_text())
    fronts, metrics, clouds = (pd.read_csv(analysis / name) for name in ("fronts.csv", "metrics.csv", "null_clouds.csv"))
    for geometry in ("raw3d", "xy"):
        comparison(fronts, metrics, clouds, geometry, out, record["cohort_size"])
        overlay(fronts, metrics, out, record["cohort_size"], geometry)
    comparison(fronts, metrics, clouds, None, out, record["cohort_size"])
    overlay_comparison(fronts, metrics, out, record["cohort_size"])
    with np.load(analysis / "assignments_and_costs.npz", allow_pickle=False) as arrays:
        tracking_comparison(fronts, metrics, record, arrays, out)
    caption = (
        "Terminal parentage across CE protein, CE RNA and CB AF16 RNA. The same 187 natural terminal edges and parent capacities "
        "are used in all panels. In the main comparison, A-C are existing 3D and D-F are 2D (XY); columns are CE protein, CE RNA "
        "and CB RNA. All six panels share axis limits and one retention scale. The separate overlay figure places 3D and XY side "
        "by side in A-B, sharing main-axis and near-natural zoom limits, species colors and line styles. "
        "Hollow diamonds mark the closest attained assignment P*, with dashed connectors from natural lineage to P*, "
        "not to maximum retention. Panel C groups the three configurations under 3D and XY and shows canonical "
        "position u, natural distance d_LP and first-cousin-null-mean distance d_NP in separate metric columns. "
        "Both distances use the same P*, and u is front arc length from the travel optimum normalized by total arc length. "
        "Travel is the equal mean of "
        "per-embryo pairwise distances normalized by fixed natural totals: "
        "three CE embryos at cutoffs 255/247/225 on the published 275-edge reference, and two CB AF16 embryos at cutoffs 148/156 "
        f"on their {record['cb_reference_size']}-edge tracking/RNA intersection. Protein uses the frozen top-20 z-scored reporters; "
        "RNA uses the frozen shared 20 TFs, cosine distance on stored values. Each front uses its own endpoint spans; within a "
        "panel all nulls and landmarks share those anchors. Colors encode biological-parent edge retention, crosses natural lineage, "
        "and outlined circles observed maximum retention. The first-cousin mean is analytic; other mean marks summarize display draws. "
        "Full-random nulls use insets with the same coordinates and separate limits. Optimization uses exact first-cousin null SDs "
        "and 301 weights with checked endpoint tie handling; 1,201-weight checks are cached separately. Lines connect attained "
        f"assignments. {TRACKING_CAVEAT} {TRACKING_CAPTION_STATUS} "
        "Developmental alignment and RNA measurement provenance remain unresolved. The XY row "
        "drops z and recomputes travel, nulls, endpoints and assignments; tracking spread does not estimate molecular uncertainty.\n"
    )
    (out / "caption.txt").write_text(caption)
    manifest = dict(analysis_id=record["analysis_id"], shared_axes=True, retention_color_scale=[0, 1],
        primary_comparison="terminal_cross_species_comparison", comparison_rows=["raw3d", "xy"],
        comparison_columns=list(PRIMARY), primary_overlay="terminal_cross_species_overlay", overlay_panels=["raw3d", "xy"],
        comparison_geometry_labels="bold second line beneath each configuration title",
        comparison_title_tracking_counts=False,
        shared_overlay_zoom_limits=True,
        comparison_natural_to_maximum_connectors=False,
        overlay_closest_point_connectors=True, canonical_panel="C",
        canonical_metrics=canonical_records(fronts, metrics),
        endpoint_normalization="separate configuration anchors; shared anchors for transfers within each panel",
        plotting_sources={str(path.relative_to(Path(__file__).resolve().parents[1])): file_hash(path) for path in (
            Path(__file__), Path(ps.__file__), Path(__file__).with_name("fig2_fig3_ce_terminal_pareto.py"))},
        files={path.name: file_hash(path) for path in sorted(out.iterdir()) if path.suffix in (".png", ".pdf", ".txt")})
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
