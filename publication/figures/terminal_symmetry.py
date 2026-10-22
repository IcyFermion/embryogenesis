"""Mirrored subtree pairs: pair numbers for the Figure 5B rows and the Figure 5C pair rows.

Input is ``TerminalPrimaryInput.symmetry`` (Track D, hash-checked by the
adapter): one row per mirrored pair with both partners' readouts and reference
percentiles, plus the evidence the caption notes cite. ``check_notes`` fails
the build when the stored results no longer support those notes. The helpers
take plain axes, so the revised subtree rows (``publication/v2``) can reuse
them.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

from publication import notation as nt
from publication import style as ps

SIDE_COLORS = dict(left=ps.COLORS["blue"], right=ps.COLORS["orange"])
SHADE = "#E6E6E6"
HIGHLIGHT = 0.9  # reference percentile at or above which a within-pair difference is shaded
ZERO = 1e-9  # d_SD at or below this is on the front (Track D4b tie tolerance)
GROUPS = ("Posterior AB", "MS", "Anterior AB (exploratory)")
# Caption notes: mark and the readouts it explains, keyed by the left partner.
NOTES = {"ABplpa": ("†", ("u",)), "MSaa": ("‡", ("u", "z_state")),
         "ABplp": ("§", ("d_lp",)), "ABplpp": ("§", ("d_lp",)), "ABplppp": ("§", ("d_lp",))}
# Figure 5C columns: readout key, notation concept, title label, axis limits and ticks.
COLUMNS = (("u", "canonical_position", "Position", (0, 1), (0, .5, 1)),
           ("d_lp", "natural_front_distance", "Distance", None, None),
           ("d_sd", "natural_front_distance_sd", "Distance", (0, 2), (0, 1, 2)),
           ("z_travel", "natural_travel_position", "Travel", (-4, 0), (-4, -2, 0)),
           ("z_state", "natural_state_position", "Cell state", (-3, 1), (-3, -1, 1)))
NULL_SD_COLUMNS = 3  # the last three columns share the null-SD unit


def _group(pair):
    if pair.tier == "exploratory":
        return GROUPS[2]
    return GROUPS[1] if pair.left.startswith("MS") else GROUPS[0]


def ordered(pairs):
    """Pairs in display order (group, then left name) with their pair number."""
    out = pairs.assign(group=[_group(p) for p in pairs.itertuples()])
    out["rank"] = out.group.map(GROUPS.index)
    out = out.sort_values(["rank", "left"]).drop(columns="rank").reset_index(drop=True)
    out["number"] = range(1, len(out) + 1)
    return out


def row_tags(pairs):
    """Per paired subtree: pair number, side, partner and the note on the readout it explains."""
    tags = {}
    for p in ordered(pairs).itertuples():
        mark, metrics = NOTES.get(p.left, (None, ()))
        for side, name, mate in (("left", p.left, p.right), ("right", p.right, p.left)):
            tags[name] = dict(number=p.number, side=side, partner=mate, mark=mark, metrics=metrics)
    return tags


def check_notes(symmetry):
    """Raise if the stored Track D results no longer support the caption notes."""
    pairs, ev = symmetry["pairs"].set_index("left"), symmetry["evidence"]
    keys = [c[0] for c in COLUMNS]
    shaded = {left for left, p in pairs.iterrows() if any(p[f"{k}_percentile"] >= HIGHLIGHT for k in keys)}
    strict = set(pairs.index[pairs.tier != "exploratory"])
    claims = {
        "every shaded pair has a note": shaded <= set(NOTES),
        "two strict pairs agree throughout": len(strict - set(NOTES)) == 2 and not (strict - set(NOTES)) & shaded,
        "seven strict pairs": len(strict) == 7,
        "ABplpa: numerical, travel range 0.04 null SD": pairs.at["ABplpa", "label"] == "numerical"
            and round(ev["abplpa_travel_span_sd"], 2) == .04,
        "ABplpa: one apoptotic cell; without it the front is a point": ev["abplpa_deaths"] == 1
            and set(ev["abplpa_without_death"].values()) == {"degenerate_span"},
        "ABplpa and ABprpa on their fronts": max(pairs.at["ABplpa", "d_sd_left"], pairs.at["ABplpa", "d_sd_right"]) <= ZERO,
        "MSaa/MSpa: data, 16 against 12 edges": pairs.at["MSaa", "label"] == "data"
            and (pairs.at["MSaa", "edges_left"], pairs.at["MSaa", "edges_right"]) == (16, 12),
        "MSpa: late-born cells past the cutoff, null on one group": {"late_born_one_side", "low_null_support"}
            <= {r.strip() for r in ev["msaa_reasons"].split(";")} and ev["msaa_null_support_share"] == .25,
        "ABplppp: one apoptotic cell, ABplpppapp": ev["abplppp_deaths"] == ["ABplpppapp"],
        "birth endpoints put ABplp, ABplpp and ABplppp on their fronts": all(ev["birth_endpoint_on_front"].values()),
        "ABplp control: 0 of 1000 matched draws": (ev["abplp_control_draws"], ev["abplp_control_on_front"]) == (1000, 0),
    }
    failed = [claim for claim, ok in claims.items() if not ok]
    if failed:
        raise ValueError(f"Figure 5 symmetry notes no longer hold: {failed}")


def tag_row(ax, row, tag, *, x=-.30, fontsize=6.5):
    """Side dot and pair number left of the first metric column (``x`` in axes fraction)."""
    transform = ax.get_yaxis_transform()
    ax.scatter([x], [row], s=11, color=SIDE_COLORS[tag["side"]], transform=transform, clip_on=False,
               linewidth=0, zorder=4)
    ax.text(x + .055, row, str(tag["number"]), transform=transform, fontsize=fontsize, color="#444444",
            va="center", ha="left", clip_on=False)


def mark_value(ax, row, value, mark, xmax, *, fontsize=8):
    """Note mark just right of a plotted readout."""
    ax.text(value + .06 * xmax, row, mark, fontsize=fontsize, color="#222222", va="center", ha="left",
            clip_on=False, zorder=4)


def _title(concept, label):
    return f"{label}\n{nt.math(concept)}"


def plot_mirrored_pairs(symmetry, out_dir, *, dmax=.2, stem="fig5C_ce_mirrored_pairs"):
    """Figure 5C: one row per mirrored pair, left and right partner on each readout."""
    check_notes(symmetry)
    pairs = ordered(symmetry["pairs"])
    slots, row = [], 0
    for group in GROUPS:
        members = pairs[pairs.group == group]
        if len(members):
            slots.append((row, group, None))
            slots += [(row + 1 + i, None, p) for i, p in enumerate(members.itertuples())]
            row += len(members) + 1
    fig = plt.figure(figsize=(7.4, 3.4))
    fig.text(.015, .985, "C", fontsize=12, weight="bold", va="top")
    fig.text(.065, .985, "Mirrored subtree pairs: left and right partner on one row",
             fontsize=10.5, weight="bold", va="top")
    left, width, gap = .275, .117, .0275
    axes = []
    for j, (key, concept, label, limits, ticks) in enumerate(COLUMNS):
        limits = limits or (0, dmax)
        ticks = ticks or (0, dmax / 2, dmax)
        ax = fig.add_axes([left + j * (width + gap), .135, width, .6])
        span = limits[1] - limits[0]
        ax.set_xlim(limits[0] - .07 * span, limits[1] + .07 * span)
        ax.set_ylim(row - .4, -.6)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}".replace("-", "−")))
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, labelsize=7.5,
                       length=2, pad=3)
        ax.set_yticks([])
        ax.set_title(_title(concept, label), fontsize=8, pad=19, linespacing=1.15)
        ax.spines[["left", "right", "bottom"]].set_visible(False)
        ax.spines["top"].set_color("#CCCCCC")
        ax.grid(axis="x", color="#E8E8E8", lw=.55)
        if limits[0] < 0 < limits[1]:
            ax.axvline(0, color="#AAAAAA", lw=.7, ls=(0, (3, 2)), zorder=1)
        axes.append((ax, key, span))
    start = left + (len(COLUMNS) - NULL_SD_COLUMNS) * (width + gap)
    end = left + len(COLUMNS) * (width + gap) - gap
    fig.text((start + end) / 2, .905, "In units of the first-cousin null SD", fontsize=8, ha="center",
             va="bottom", color="#444444")
    fig.add_artist(Line2D([start, end], [.9, .9], transform=fig.transFigure, color="#AAAAAA", lw=.6))
    label_x = (.02 - left) / width  # first-column axes fraction of the label block
    for slot, group, p in slots:
        first = axes[0][0]
        if group:
            first.text(label_x, slot, group, transform=first.get_yaxis_transform(), fontsize=8.3,
                       weight="bold", va="center", clip_on=False)
            continue
        mark = NOTES.get(p.left, ("", ()))[0]
        first.text(label_x, slot, f"{p.number}", transform=first.get_yaxis_transform(), fontsize=7,
                   color="#666666", va="center", clip_on=False)
        first.text(label_x + .03 / width, slot, f"{p.left} · {p.right}",
                   transform=first.get_yaxis_transform(), fontsize=8.3, color="#222222", va="center",
                   clip_on=False)
        if mark:
            first.text(-.12, slot, mark, transform=first.get_yaxis_transform(), fontsize=8.5,
                       color="#222222", va="center", ha="right", clip_on=False)
        hollow = p.tier == "exploratory"
        for ax, key, span in axes:
            a, b = getattr(p, f"{key}_left"), getattr(p, f"{key}_right")
            if getattr(p, f"{key}_percentile") >= HIGHLIGHT:
                ax.axhspan(slot - .42, slot + .42, color=SHADE, lw=0, zorder=0)
            ax.hlines(slot, min(a, b), max(a, b), color="#B5B5B5", lw=1.1, zorder=2)
            offset = .17 if abs(a - b) < .035 * span else 0
            for value, side, dy in ((b, "right", offset), (a, "left", -offset)):
                color = SIDE_COLORS[side]
                ax.scatter(value, slot + dy, s=20, zorder=3, linewidth=.9 if hollow else .35,
                           facecolor="white" if hollow else color, edgecolor=color if hollow else "white")
    handles = [Line2D([], [], marker="o", ls="", ms=5, color=SIDE_COLORS["left"], label="Left partner"),
               Line2D([], [], marker="o", ls="", ms=5, color=SIDE_COLORS["right"], label="Right partner"),
               Line2D([], [], marker="o", ls="", ms=5, markerfacecolor="white", markeredgecolor="#555555",
                      label="Exploratory pair"),
               Patch(facecolor=SHADE, edgecolor="none", label="Difference ≥ 90th percentile of matched pairs")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.015, .0), ncol=4, frameon=False,
               fontsize=7.8, handletextpad=.3, columnspacing=1.6, borderaxespad=.2)
    out_dir = Path(out_dir)
    for extension in ("pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{extension}", dpi=300, facecolor="white",
                    bbox_inches=Bbox.from_bounds(0, 0, *fig.get_size_inches()))
    return fig
