"""Round-by-round sketches: Figure 7 as Figure-2-style round panels, and per-round cross-species comparisons."""

from __future__ import annotations

from pathlib import Path
import shutil
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np

from publication import notation as nt
from publication import style
from publication.artists import AXIS_LABELS, proportional_limits
from publication.figures.cross_species import FAMILY_LABELS, _display

CONFIG_MARKERS = {"ce_protein": "o", "ce_rna": "s", "cb_rna": "^"}


def round_title(r, edges, last):
    suffix = " (bottom)" if r == 1 else " (top)" if r == last else ""
    return f"Round {r}{suffix} · {edges} edges"


def round_panel(ax, data, geometry, config, *, inset=True, marker_scale=1.):
    """One Figure-2-style panel: cousin clouds, retention-colored front, L, maximum retention, full-random inset."""
    refs = data.references
    curve, metric = data.select(geometry, config)
    null = data.cloud(geometry, config)
    ax.axhline(0, color="#777777", lw=.6, ls=":", zorder=0)
    ax.axvline(0, color="#777777", lw=.6, ls=":", zorder=0)
    xs, ys = list(curve.D1) + [metric.natural_D1], list(curve.D2) + [metric.natural_D2]
    for family in refs["main"]:
        full = null[null.family == family]
        sample = _display(full, refs["display_draws"])
        color = style.NULL_MODEL_COLORS[family]
        ax.scatter(sample.D1, sample.D2, s=6*marker_scale, color=color, alpha=.19, edgecolors="none", rasterized=True, zorder=1)
        mean = (metric.null_D1, metric.null_D2) if family == "first_cousin" else (full.D1.mean(), full.D2.mean())
        ax.scatter(*mean, marker="+", s=30*marker_scale, color=color, lw=1.1, zorder=4)
        xs += list(sample.D1)
        ys += list(sample.D2)
    ax.plot(curve.D1, curve.D2, color=style.COLORS["blue"], lw=1.1, alpha=.5, zorder=2)
    front = ax.scatter(curve.D1, curve.D2, c=curve.edge_retention, cmap=style.EDGE_RETENTION_CMAP, norm=Normalize(0, 1),
                       s=12*marker_scale, edgecolors="none", zorder=3)
    maximum = curve.iloc[int(metric.maximum_index)]
    ax.plot([metric.natural_D1, maximum.D1], [metric.natural_D2, maximum.D2], color="#222222", lw=.7,
            ls=(0, (2, 2)), zorder=5)
    ax.scatter(maximum.D1, maximum.D2, s=26*marker_scale, facecolor=style.EDGE_RETENTION_CMAP(maximum.edge_retention),
               edgecolor="#222222", lw=.9, zorder=7)
    ax.scatter(metric.natural_D1, metric.natural_D2, marker="X", s=60*marker_scale, color="#222222", edgecolor="white",
               lw=.5, zorder=8)
    for values, setter, headroom in ((xs, ax.set_xlim, 1.), (ys, ax.set_ylim, 1.55 if inset else 1.)):
        low, high = min(values), max(values)
        pad = max(.05 * (high - low), .04)
        # Raise the top so the data stay below the inset in the upper-right corner.
        setter(low - pad, low + headroom * (high - low) + pad)
    ax.grid(True, alpha=.3)
    if inset:
        full = null[null.family == refs["inset"]]
        color = style.NULL_MODEL_COLORS["full_random"]
        box = ax.inset_axes([.66, .64, .31, .31])
        box.scatter(full.D1, full.D2, s=3, color=color, alpha=.18, edgecolors="none", rasterized=True)
        box.scatter(full.D1.mean(), full.D2.mean(), marker="+", s=22, color=color, lw=1.)
        box.set_xlim(*proportional_limits(full.D1))
        box.set_ylim(*proportional_limits(full.D2))
        box.tick_params(labelsize=5, length=1.5, pad=1)
        box.locator_params(axis="both", nbins=3)
        box.grid(True, alpha=.25)
        box.spines[["top", "right"]].set_visible(True)
    return front


def legend_handles(refs):
    handles = [Line2D([], [], marker="o", ls="", color=style.NULL_MODEL_COLORS[f], alpha=.7, markersize=5,
                      markeredgecolor="none", label=FAMILY_LABELS[f]) for f in refs["main"]]
    handles.append(Line2D([], [], marker="o", ls="", color=style.NULL_MODEL_COLORS["full_random"], alpha=.7,
                          markersize=5, markeredgecolor="none", label=refs["inset_legend"]))
    handles.append(Line2D([], [], marker="+", ls="", color=style.NULL_MODEL_COLORS["first_cousin"], markersize=7,
                          label=f"{nt.label('null_mean')} {nt.math('null_mean')} (exact)"))
    handles += [Line2D([], [], marker="X", ls="", color="#222222", markersize=7, label=nt.label("natural_lineage")),
                Line2D([], [], marker="o", ls="", color=style.COLORS["blue"], markeredgecolor="#222222", markersize=5,
                       label="Maximum edge retention")]
    return handles


def round_grid(rounds, *, geometry="raw3d", config="ce_protein", title=None, columns=4):
    """All rounds of one configuration as Figure-2-style panels."""
    keys = sorted(rounds)
    nrows = int(np.ceil(len(keys) / columns))
    fig, axes = plt.subplots(nrows, columns, figsize=(3.55 * columns, 3.3 * nrows + .9), squeeze=False)
    fig.subplots_adjust(left=.06, right=.92, bottom=.15 if nrows > 1 else .25, top=.9, hspace=.42, wspace=.28)
    for ax, r in zip(axes.flat, keys):
        front = round_panel(ax, rounds[r], geometry, config, marker_scale=.8)
        ax.set_title(round_title(r, rounds[r].edges, max(keys)), fontsize=9)
    for ax in axes.flat[len(keys):]:
        ax.set_visible(False)
    fig.supxlabel(AXIS_LABELS[0].replace("\n", " "), y=.085 if nrows > 1 else .13, fontsize=9)
    fig.supylabel(AXIS_LABELS[1].replace("\n", " "), x=.015, fontsize=9)
    if title:
        fig.suptitle(title, fontsize=11.5, y=.975)
    bar = fig.colorbar(front, cax=fig.add_axes([.935, .2, .012, .65]))
    bar.set_label(nt.label("retention"))
    fig.legend(handles=legend_handles(rounds[keys[0]].references), loc="lower center", ncol=6, fontsize=7.5,
               bbox_to_anchor=(.49, .005))
    return fig


def figure2_single(data, out, stem, *, geometry="raw3d", config="ce_protein"):
    """Exact Figure 2 layout (``terminal_global.plot_main``) applied to one round."""
    from publication.figures import terminal_global
    curve, metric = data.select(geometry, config)
    null = data.cloud(geometry, config)
    pick = lambda fam, n: _display(null[null.family == fam], n)
    shown = {k: (lambda s: (s.D1.to_numpy(), s.D2.to_numpy()))(pick(fam, 1000 if k != "full" else None))
             for k, fam in ((1, "first_cousin"), (2, "second_cousin"), (3, "third_cousin"), ("full", "full_random"))}
    display = dict(x=curve.D1.to_numpy(), y=curve.D2.to_numpy(), natural=(metric.natural_D1, metric.natural_D2),
                   nulls=shown, xlabel=AXIS_LABELS[0], ylabel=AXIS_LABELS[1])
    twr = dict(traditional_er=curve.edge_retention.to_numpy(), cost_tree_tradeoff=None,
               kp=dict(max_er_idx=int(metric.maximum_index)))
    with tempfile.TemporaryDirectory() as tmp:
        fig = terminal_global.plot_main(twr, None, out_dir=tmp, display=display)
        plt.close(fig)
        for ext in ("pdf", "png"):
            shutil.move(Path(tmp) / f"fig2_ce_terminal_pareto_front.{ext}", Path(out) / f"{stem}.{ext}")


def round_metrics(metrics, *, configs, geometries, title, aggregate=None, degenerate_note=None):
    """Canonical metrics by round (rows), one marker per configuration; geometries stacked."""
    keys = (("canonical_position", style.SEMANTIC_COLORS["canonical_u"]),
            ("natural_front_distance", style.SEMANTIC_COLORS["lineage_front_distance"]),
            ("null_front_distance", style.SEMANTIC_COLORS["first_cousin_null"]),
            ("max_edge_retention", style.COLORS["blue"]))
    rounds = sorted(metrics["round"].unique())
    labels = [f"Round {r} · {int(metrics.loc[metrics['round'] == r, 'edges'].iloc[0])}" for r in rounds]
    if aggregate is not None:
        labels.append("Aggregate")
    n = len(labels)
    fig, axes = plt.subplots(len(geometries), len(keys), figsize=(12.5, (.32 * n + 1.1) * len(geometries) + .9),
                             sharey=True, squeeze=False)
    fig.subplots_adjust(left=.12, right=.98, top=1 - .9 / fig.get_figheight(), bottom=.9 / fig.get_figheight(),
                        hspace=.45, wspace=.12)
    offsets = np.linspace(-.22, .22, len(configs)) if len(configs) > 1 else [0.]
    for g, geometry in enumerate(geometries):
        for k, (key, color) in enumerate(keys):
            ax = axes[g, k]
            for c, (config, offset) in enumerate(zip(configs, offsets)):
                rows = metrics[(metrics.geometry == geometry) & (metrics.config == config)].set_index("round")
                values = [rows[key].get(r, np.nan) for r in rounds]
                if aggregate is not None:
                    a = aggregate[(aggregate.geometry == geometry) & (aggregate.config == config)]
                    values.append(float(a[key].iloc[0]) if len(a) else np.nan)
                y = np.arange(n) + offset
                single = len(configs) == 1
                ax.hlines(y, 0, values, color=color if single else style.SPECIES_COLORS[config], lw=.8, alpha=.4)
                ax.scatter(values, y, s=26, marker=CONFIG_MARKERS[config], zorder=3, edgecolor="white", lw=.4,
                           color=color if single else style.SPECIES_COLORS[config],
                           label=nt.CONFIG_LABELS[config] if (g, k) == (0, 0) else None)
            if aggregate is not None:
                ax.axhline(n - 1.5, color="#BBBBBB", lw=.6, ls="--")
            ax.set_ylim(n - .4, -.6)
            ax.set_yticks(range(n), labels)
            ax.tick_params(axis="y", labelleft=k == 0)
            ax.grid(axis="x", color="#E8E8E8", lw=.55)
            ax.set_xlim(left=-.02 * max(1e-9, np.nanmax(metrics[key])))
            heading = nt.heading(key) if key != "max_edge_retention" else "Maximum edge retention"
            ax.set_title(heading if g == 0 else "", fontsize=9, color=color)
            if k == 0:
                ax.set_ylabel(f"{nt.GEOMETRY_LABELS[geometry]} tracking", fontsize=9, weight="bold")
    fig.suptitle(title, fontsize=11, y=.995)
    u = nt.symbol("canonical_position")
    note = (f"${u}=0$: travel optimum; ${u}=1$: cell-state optimum. Both distances use each round's own "
            f"{nt.math('closest_point')} and exact first-cousin null mean {nt.math('null_mean')}.")
    if degenerate_note:
        note += " " + degenerate_note
    fig.text(.12, .012, note, fontsize=7.8, color="#444444")
    if len(configs) > 1:
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center", ncol=len(configs), fontsize=8,
                   bbox_to_anchor=(.55, 1 - .35 / fig.get_figheight()))
    return fig
