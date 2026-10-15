"""Shared sketch panels: named-row canonical metrics with one marker per configuration."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from publication import notation as nt
from publication import style

CONFIG_MARKERS = {"ce_protein": "o", "ce_rna": "s", "cb_rna": "^"}
METRIC_COLUMNS = (("canonical_position", style.SEMANTIC_COLORS["canonical_u"]),
                  ("natural_front_distance", style.SEMANTIC_COLORS["lineage_front_distance"]),
                  ("null_front_distance", style.SEMANTIC_COLORS["first_cousin_null"]),
                  ("max_edge_retention", style.COLORS["blue"]))


def metric_rows(table, *, row_key, rows, configs, geometries, title, note="", separators=(), headings=None,
                bold=(), row_height=.32, columns=METRIC_COLUMNS):
    """Figure 5B-style named rows: ``rows`` is [(row id, label)]; geometries stacked, metrics side by side.

    ``table`` has geometry, config, ``row_key`` and metric columns under notation concept keys
    (missing or NaN values are left blank). ``headings`` maps a row index to a group title drawn above it.
    """
    headings = headings or {}
    # Each group heading takes its own gap row so it never overlaps a row label.
    position = np.arange(len(rows)) + np.array([sum(1 for h in headings if h <= i) for i in range(len(rows))])
    n = int(position[-1]) + 1
    height = (row_height * n + 1.1) * len(geometries) + .9
    fig, axes = plt.subplots(len(geometries), len(columns), figsize=(3.1 * len(columns) + .1, height), sharey=True,
                             squeeze=False)
    fig.subplots_adjust(left=.13, right=.98, top=1 - 1. / height, bottom=.75 / height, hspace=.3 if n > 15 else .45,
                        wspace=.12)
    offsets = np.linspace(-.22, .22, len(configs)) if len(configs) > 1 else [0.]
    single = len(configs) == 1
    for g, geometry in enumerate(geometries):
        for k, (key, color) in enumerate(columns):
            ax = axes[g, k]
            for config, offset in zip(configs, offsets):
                frame = table[(table.geometry == geometry) & (table.config == config)].set_index(row_key)
                values = np.array([frame[key].get(r, np.nan) if key in frame else np.nan for r, _ in rows], float)
                y = position + offset
                tint = color if single else style.SPECIES_COLORS[config]
                ax.hlines(y, 0, values, color=tint, lw=.8, alpha=.4)
                ax.scatter(values, y, s=24, marker=CONFIG_MARKERS[config], zorder=3, edgecolor="white", lw=.4,
                           color=tint, label=nt.CONFIG_LABELS[config] if (g, k) == (0, 0) else None)
            for index in separators:
                ax.axhline(position[index] - .5, color="#BBBBBB", lw=.6, ls="--")
            ax.set_ylim(n - .4, -.6)
            ax.set_yticks(position, [label for _, label in rows])
            ax.tick_params(axis="y", labelleft=k == 0, labelsize=7.5 if n > 15 else 8)
            for tick, (r, _) in zip(ax.get_yticklabels(), rows):
                if r in bold:
                    tick.set_fontweight("bold")
            ax.grid(axis="x", color="#E8E8E8", lw=.55)
            finite = table[key].to_numpy(float) if key in table else np.array([1.])
            ax.set_xlim(-.03 * max(1e-9, np.nanmax(finite)), 1.05 * max(1e-9, np.nanmax(finite)))
            heading = nt.heading(key) if key != "max_edge_retention" else "Maximum edge retention"
            ax.set_title(heading if g == 0 else "", fontsize=9, color=color)
            if k == 0:
                ax.set_ylabel(f"{nt.GEOMETRY_LABELS[geometry]} tracking", fontsize=9, weight="bold")
            if k == 0:
                for index, text in headings.items():
                    ax.text(-.02, position[index] - 1, text, transform=ax.get_yaxis_transform(), ha="right",
                            va="center", fontsize=8, weight="bold", color="#555555", clip_on=False)
    fig.suptitle(title, fontsize=11, y=1 - .12 / height)
    u = nt.symbol("canonical_position")
    fig.text(.13, .1 / height, f"${u}=0$: travel optimum; ${u}=1$: cell-state optimum. " + note, fontsize=7.8,
             color="#444444", wrap=True)
    if not single:
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center", ncol=len(configs), fontsize=8,
                   bbox_to_anchor=(.55, 1 - .45 / height))
    return fig
