"""Reserve note: one apoptotic cell's track endpoint and the ABplp/ABprp difference (Track D4-D4c).

Answers the question a reader of Figure 5C (note §) may ask: is the posterior
ABpl/ABpr difference a break in left-right symmetry? Every number comes from
hash-checked Track D tables, and ``check`` fails the build if a sentence no
longer holds.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from publication import assembly, notation as nt, style
from publication.figures.terminal_symmetry import SIDE_COLORS
from publication.provenance import sha256

STEM = "noteA_apoptosis_abplp"
LABEL = "A"
RUNS = dict(compare="compare_20261015", deaths="apoptosis_20261019", control="apoptosis_subtrees_20261021",
            single="apoptosis_single_20261022")
TABLES = dict(compare=("subtree_metrics", "pair_comparison"), deaths=("decision",), control=("decision",),
              single=("single", "endpoints", "readouts"))
CELL = "ABplpppapp"
LEFT = ("ABplp", "ABplpp", "ABplppp", "ABpl")
PAIRS = ("ABplp", "ABplpp", "ABplppp")  # the kept pairs marked in Figure 5C
DYING = (("ABplpppapp", "ABprpppapp"), ("ABplppaaap", "ABprppaaap"))
ZERO = 1e-9  # Track D4b tie tolerance


def mirror(name):
    return name[:3] + {"l": "r", "r": "l"}[name[3]] + name[4:]


def load():
    """Hash-checked Track D tables, keyed (stage, name), and the input files read."""
    from terminal_pareto import subtree_symmetry as ss
    tables, files = {}, {}
    for stage, run_id in RUNS.items():
        folder = ss.run_path(run_id) / "analysis"
        record = json.loads((folder / "provenance.json").read_text())
        for key, expected in (("compare_run", RUNS["compare"]), ("apoptosis_run", RUNS["deaths"])):
            if key in record and record[key] != expected:
                raise ValueError(f"{run_id} was built from {record[key]}, not {expected}")
        files[str(folder.relative_to(ss.RUN_ROOT) / "provenance.json")] = sha256(folder / "provenance.json")
        for name in TABLES[stage]:
            path = folder / f"{name}.csv"
            if sha256(path) != record["table_hashes"][name]:
                raise ValueError(f"{path} differs from its provenance record")
            tables[stage, name] = pd.read_csv(path)
            files[str(path.relative_to(ss.RUN_ROOT))] = sha256(path)
    return tables, files


def facts(tables):
    """The numbers the note prints, gathered in one place."""
    single = tables["single", "single"]
    t = single.set_index(["subtree", "treatment", "death", "variant"]).sort_index()
    as_is = single[single.treatment == "none"].set_index("subtree").d_SD
    rows = []
    for s in LEFT:
        others = single[(single.subtree == s) & (single.treatment == "all_but") & (single.death == CELL)]
        rows.append(dict(
            subtree=s, mirror=mirror(s), left=as_is[s], right=as_is[mirror(s)],
            alone={v: (t.loc[(s, "only", CELL, v), "status"], t.loc[(s, "only", CELL, v), "d_SD"])
                   for v in ("drop_deaths", "birth_endpoint")},
            others={r.variant: (r.status, r.d_SD) for r in others.itertuples()},
            right_max=float(single[(single.subtree == mirror(s)) & (single.status == "ok")].d_SD.max())))
    ends = tables["single", "endpoints"]
    cell = ends[ends.cell == CELL].sort_values("embryo")
    pc = tables["compare", "pair_comparison"]
    z = pc[pc.subtree_x.isin(PAIRS) & pc.metric.isin(["natural_travel_z", "natural_state_z"])]
    control = tables["control", "decision"].set_index(["subtree", "variant"]).loc[("ABplp", "drop_deaths")]
    whole = tables["deaths", "decision"].set_index("variant").loc["drop_deaths"]
    return dict(
        rows=rows, readouts=tables["single", "readouts"].set_index("readout").holds.to_dict(),
        travel=list(zip(cell.travel, cell.mirror_travel)),
        frames=((int(cell.frames_to_endpoint.min()), int(cell.frames_to_endpoint.max())),
                (int(cell.mirror_frames_to_endpoint.min()), int(cell.mirror_frames_to_endpoint.max()))),
        second=[(float(a), float(b)) for a, b in ends[ends.cell == DYING[1][0]].sort_values("embryo")[
            ["travel", "mirror_travel"]].itertuples(index=False)],
        z_percentiles=(z.reference_percentile.min(), z.reference_percentile.max()),
        control=dict(draws=int(control.n_ok), on_front=int(round(control.d_SD_p * (control.n_ok + 1))) - 1,
                     median=float(control.d_SD_control_median)),
        whole=dict(n_deaths=30, d_SD_p=float(whole.d_SD_p), retention_p=float(whole.retention_at_nearest_p)),
        endpoints=ends)


def check(f):
    """Raise if a sentence of the note no longer holds."""
    rows = {r["subtree"]: r for r in f["rows"]}
    on = lambda status_value: status_value[0] == "ok" and status_value[1] <= ZERO
    claims = {
        "stored D4c readouts hold": all(f["readouts"].values()),
        "left subtrees off, right partners on, as is": all(r["left"] > ZERO and r["right"] <= ZERO for r in rows.values()),
        "ABplpppapp alone puts each left subtree on its front":
            all(on(r["alone"]["birth_endpoint"]) for r in rows.values())
            and all(on(r["alone"]["drop_deaths"]) or r["alone"]["drop_deaths"][0] == "below_min_edges"
                    for r in rows.values()),
        "the other deaths leave each left subtree off its front":
            all(s == "ok" and v > ZERO for r in rows.values() for s, v in r["others"].values()),
        "no treatment moves a right partner off its front": all(r["right_max"] <= ZERO for r in rows.values()),
        "left endpoint farther in every embryo": all(a > b for a, b in f["travel"]),
        "ABplppaaap farther in two of three embryos": sum(a > b for a, b in f["second"]) == 2,
        "ABplp control: 0 of 1000": (f["control"]["draws"], f["control"]["on_front"]) == (1000, 0),
        "whole cohort: deaths not special": f["whole"]["d_SD_p"] > .05 and f["whole"]["retention_p"] > .05,
    }
    failed = [claim for claim, ok in claims.items() if not ok]
    if failed:
        raise ValueError(f"Note {LABEL} no longer holds: {failed}")


def endpoint_figure(ends, out: Path):
    """Parent division to track end, left cell against its mirror, per embryo."""
    style.configure()
    fig, ax = plt.subplots(figsize=(3.3, 1.95))
    rows, labels = [], []
    for k, (left, right) in enumerate(DYING):
        for e, embryo in enumerate(("embryo1", "embryo2", "embryo3")):
            r = ends[(ends.cell == left) & (ends.embryo == embryo)].iloc[0]
            y = k * 4.2 + e
            ax.hlines(y, min(r.travel, r.mirror_travel), max(r.travel, r.mirror_travel), color="#B5B5B5", lw=1.1)
            ax.scatter([r.mirror_travel], [y], s=18, color=SIDE_COLORS["right"], edgecolor="white", lw=.35, zorder=3)
            ax.scatter([r.travel], [y], s=18, color=SIDE_COLORS["left"], edgecolor="white", lw=.35, zorder=3)
            rows.append(y)
            labels.append(f"embryo {e + 1}")
        ax.text(-.02, k * 4.2 - 1.05, f"{left} / {right}", transform=ax.get_yaxis_transform(),
                fontsize=6.8, weight="bold", ha="left", va="center")
    ax.set_yticks(rows, labels, fontsize=6.5)
    ax.set_ylim(max(rows) + .7, -1.6)
    ax.set_xlim(0, 9.6)
    ax.set_xlabel("Parent division to track end (µm)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5, length=2)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", lw=.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color=SIDE_COLORS["left"], label="left cell"),
                       Line2D([], [], marker="o", ls="", ms=4, color=SIDE_COLORS["right"], label="mirror cell")],
              loc="lower right", bbox_to_anchor=(1, 1.02), ncol=2, fontsize=6.3, frameon=False,
              handletextpad=.2, columnspacing=.8, borderaxespad=0)
    fig.tight_layout(pad=.2)
    fig.savefig(out / f"{STEM}_endpoints.pdf")
    plt.close(fig)
    return f"{STEM}_endpoints.pdf"


def _tt(name):
    return name.replace("_", r"\_")


def _d(status_value):
    status, value = status_value
    if status == "below_min_edges":
        return r"{\footnotesize $<$12 edges}"
    return "0" if value <= ZERO else f"{value:.2f}"


def table(f):
    body = []
    for r in f["rows"]:
        others = (" / ".join(_d(r["others"][v]) for v in ("drop_deaths", "birth_endpoint"))
                  if r["others"] else "---")
        body.append(f"{r['subtree']} & {_d(('ok', r['left']))} & {_d(r['alone']['drop_deaths'])} / "
                    f"{_d(r['alone']['birth_endpoint'])} & {others} & {r['mirror']} & {_d(('ok', r['right']))} & "
                    f"{_d(('ok', r['right_max']))} \\\\")
    d = nt.tex('natural_front_distance_sd')
    return (r"\begin{tabular}{@{}l c c c l c c@{}}" "\n" r"\toprule" "\n"
            rf" & \multicolumn{{3}}{{c}}{{Left subtree, {d}}} & \multicolumn{{3}}{{c}}{{Mirror subtree, {d}}} \\" "\n"
            r"\cmidrule(lr){2-4}\cmidrule(l){5-7}" "\n"
            rf"Subtree & as is & {CELL} alone & other deaths & Subtree & as is & any treatment \\" "\n"
            r"\midrule" "\n" + "\n".join(body) + "\n" r"\bottomrule" "\n" r"\end{tabular}")


def tex(f, figure):
    t = nt.tex
    (l1, l2, l3), (r1, r2, r3) = zip(*[(a, b) for a, b in f["travel"]])
    (lf0, lf1), (rf0, rf1) = f["frames"]
    zp0, zp1 = (round(100 * p) for p in f["z_percentiles"])
    rows = {r["subtree"]: r for r in f["rows"]}
    c, w = f["control"], f["whole"]
    return rf"""\documentclass[10pt]{{article}}
\usepackage[margin=0.75in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{graphicx}}
\usepackage{{amsmath}}
\usepackage[T1]{{fontenc}}
\usepackage{{textcomp}}
\pagestyle{{empty}}
\setlength{{\parskip}}{{4pt}}
\setlength{{\parindent}}{{0pt}}
\begin{{document}}
{{\large\textbf{{Supplementary Note {LABEL}: One apoptotic cell's track endpoint accounts for\\
the \mbox{{ABplp--ABprp}} difference}}}}

\textbf{{Question.}} In Figure~5C (note~\S), ABplp, ABplpp and ABplppp lie off their terminal
fronts ({t('natural_front_distance_sd')} = {rows['ABplp']['left']:.2f}, {rows['ABplpp']['left']:.2f} and
{rows['ABplppp']['left']:.2f}) while their mirror subtrees ABprp, ABprpp and ABprppp lie on theirs; ABpl
against ABpr shows the same ({rows['ABpl']['left']:.2f} against 0). Is this a break in left--right symmetry?
The readouts that do not depend on the front, the travel and cell-state positions {t('natural_travel_position')}
and {t('natural_state_position')}, differ between partners no more than between random subtree pairs
matched for size and depth ({zp0}th--{zp1}th percentile).

\textbf{{One cell.}} Each left subtree was recomputed with each of its apoptotic cells treated alone, and
with all its other apoptotic cells treated together. A cell is treated either by removing it from the
cohort or by holding its track endpoint at its birth position (table, removal / birth endpoint).
Treating {CELL} alone puts every left subtree on its front; treating all the other apoptotic cells
leaves every left subtree off its front. The same treatments never move a mirror subtree off its front.
The effect is specific to this cell: in ABplp, applying the same change to three living cells of the
subtree, matched on depth and birth time, never places the lineage on its front ({c['on_front']} of
{c['draws']} draws; median {t('natural_front_distance_sd')} {c['median']:.2f}).

\begin{{center}}
{table(f)}
\end{{center}}

\begin{{minipage}}[t]{{.5\textwidth}}
\vspace{{0pt}}
\textbf{{The endpoint.}} {CELL} and its mirror cell ABprpppapp both undergo apoptosis, and both
tracks continue for a similar time after birth ({lf0}--{lf1} and {rf0}--{rf1} frames). In every
embryo the left track ends farther from the parent's division: {l1:.1f}, {l2:.1f} and {l3:.1f}~\textmu m
against {r1:.1f}, {r2:.1f} and {r3:.1f}~\textmu m (right). Because holding the endpoint at birth removes the
difference, the travel cost carries it, not the cell state. The other apoptotic cell of ABplpp,
ABplppaaap (lower rows), also ends farther than its mirror in {sum(a > b for a, b in f['second'])} of 3
embryos, but treating it moves no subtree onto its front.
\end{{minipage}}\hfill
\begin{{minipage}}[t]{{.46\textwidth}}
\vspace{{0pt}}
\centering
\includegraphics[width=\linewidth]{{{figure}}}
\end{{minipage}}

\textbf{{Reading.}} The left and right lineages are symmetric; the difference lies in where one dying
cell's track ends. It recurs in all three embryos, so it may reflect a real left--right difference in
how the corpse is displaced or engulfed rather than tracking noise; these data cannot separate the two,
and no conclusion of the paper depends on it. Apoptotic cells remain in the analysis: across the whole
cohort, removing all {w['n_deaths']} of them changes {t('natural_front_distance_sd')} and the edge
retention at the nearest front point no more than removing {w['n_deaths']} living cells matched on depth
and birth time ($p = {w['d_SD_p']:.2f}$ and ${w['retention_p']:.2f}$).

{{\footnotesize\textbf{{Data.}} Track D runs \texttt{{{_tt(RUNS['single'])}}} (single-cell treatments,
endpoints), \texttt{{{_tt(RUNS['control'])}}} (matched control within ABplp), \texttt{{{_tt(RUNS['deaths'])}}}
(whole cohort) and \texttt{{{_tt(RUNS['compare'])}}} (pair comparison), each replay-verified; pooled travel
over three embryos, top-20 protein cosine cell-state cost, first-cousin null.}}
\end{{document}}
"""


def build(out: Path) -> dict:
    out = Path(out)
    tables, inputs = load()
    f = facts(tables)
    check(f)
    figure = endpoint_figure(f["endpoints"], out)
    path = out / f"{STEM}.tex"
    path.write_text(tex(f, figure))
    compiled = assembly.compile_wrapper(path, label=LABEL, kind="Supplementary Note")
    return dict(file=f"{STEM}.pdf", compiled=compiled, inputs=inputs,
                source=sha256(Path(__file__)), notation=nt.snapshot())
