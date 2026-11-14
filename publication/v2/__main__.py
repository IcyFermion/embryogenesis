"""Build the revised manuscript layout: ``python -m publication.v2 --output-dir DIR``.

Draft pages carry panels and panel letters only (captions: ``publication.v2.captions``).
Figure 1 is the schematic figure: panels A-C and D are the authors' OmniGraffle
files from ``paper/FIGS/``, included unchanged. Colours come from
``publication.palette``; the build also draws the palette sheet (``palette.pdf``).
Numbering of 2026-11-10: the cell-state feature figure is Supplementary Figure S2
(next to S1, travel pooled across embryos), so the main figures are 1-7 again.
Figure 7 (C. briggsae) is its front alone; its comparisons with C. elegans (the
former 8B-C) are built as a backlog page, ``backlog_fig7BC_briggsae_counterparts``.
Since 2026-11-13 the build makes every figure of the manuscript set, one file per
figure, and Table 1; nothing is taken from a released figure set. Figure 5 and
S1 keep the released panels, drawn by the released figure code
(``publication.figures``) in its own style.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from publication import notation as nt
from publication import palette, provenance, style, tables
from publication.figures import terminal_cell_types, terminal_s3, terminal_within_type
from publication.provenance import sha256
from publication.v2 import data, panels

ROOT = Path(__file__).resolve().parents[2]
ASSETS = Path(__file__).resolve().parent / "assets"
EXTERNAL = {"Fig1_overview.pdf": ROOT / "paper/FIGS/Fig1_overview.pdf",
            "FigX_distances_schematic.pdf": ROOT / "paper/FIGS/FigX_distances_schematic.pdf"}
# Collaborator suggestion adopted for this candidate (2026-10-08); one-line switch in notation.py on adoption.
NOTATION = dict(retention=dict(label="Lineage similarity"))
# Figure 7 clouds sit side by side, so its labels are placed by hand (data coordinates).
BRIGGSAE_LABELS = {"lineage": (.075, .385, "left", "top"), "first_cousin": (.135, .60, "center", "bottom"),
                   "second_cousin": (.265, .44, "left", "top"), "third_cousin": (.575, .60, "left", "center")}

PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{graphicx}
\usepackage{tikz}
\pagestyle{empty}
\newcommand{\panel}[2]{\par\noindent{\textbf{\large\sffamily #1}}\par\vspace{-1.5ex}\noindent\hspace{5mm}#2\par\medskip}
\begin{document}
"""
# Manuscript mode: the figure alone, cropped to its content, at the draft pages' text width (8.5in - 2 x 0.55in).
MANUSCRIPT_PREAMBLE = r"""\documentclass[varwidth=7.4in,border=4pt]{standalone}
\usepackage{graphicx}
\usepackage{tikz}
\newcommand{\panel}[2]{\par\noindent{\textbf{\large\sffamily #1}}\par\vspace{-1.5ex}\noindent\hspace{5mm}#2\par\medskip}
\begin{document}
"""


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return f"{stem}.pdf"


def page(out, stem, title, body, note="Draft layout; caption to follow.", manuscript=False):
    """A draft page with its title and note or, in ``manuscript`` mode, the figure alone, cropped to its content."""
    tex = out / f"{stem}.tex"
    head = "" if manuscript else (f"\\noindent{{\\sffamily\\bfseries {title}}}"
                                  f"\\hfill{{\\sffamily\\small\\color{{gray}}{note}}}\\par\\bigskip\n")
    tex.write_text((MANUSCRIPT_PREAMBLE if manuscript else PREAMBLE) + head + body + "\n\\end{document}\n")
    result = subprocess.run(["tectonic", tex.name], cwd=out, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{tex.name} failed to compile:\n{result.stdout[-3000:]}")
    return f"{stem}.pdf"


def layer_ceiling(*tables):
    """Shared d_LP axis of Figure 6B and backlog 7C, so the two species compare directly if 7C is ever used."""
    return max(.3, np.ceil(max(float(t.front_distance.max()) for t in tables) * 20) / 20 + .05)


def figure6(out, manuscript=False):
    """Figure 6 page (full tree by layer); callable alone (the other figures need the terminal analysis context)."""
    # A: every layer's front (layers 1-2 show why the summaries start at layer 3).
    # B-C: Track E's per-layer summaries on the birth-size axis (author decision 2026-11-05), drawn with the palette.
    summary = data.layer_summary().layers
    a = save(panels.layer_grid(data.layers()), out, "fig6A_layer_fronts")
    b = save(panels.layer_position(summary, dmax=layer_ceiling(summary, data.briggsae_layers())), out,
             "fig6B_layer_position")
    c = save(panels.layer_cost(summary), out, "fig6C_layer_cost")
    return page(out, "fig6_full_tree_layers", "Figure 6",
        f"\\panel{{A}}{{\\includegraphics[width=.95\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{B}}{{\\includegraphics[width=.9\\textwidth]{{{b}}}}}\\end{{minipage}}"
        f"\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{C}}{{\\includegraphics[width=.9\\textwidth]{{{c}}}}}\\end{{minipage}}",
        manuscript=manuscript)


def released_panels(out, render, names):
    """Draw panels with the released figure code (its own style) in a scratch folder and keep ``names``
    ({released stem: stem here}, PDF and PNG) in ``out``. Returns the kept PDFs in order."""
    scratch = out / "_released_panels"
    scratch.mkdir()
    try:
        with plt.rc_context():
            style.configure()
            render(scratch)
            plt.close("all")
        for old, new in names.items():
            for ext in ("pdf", "png"):
                (scratch / f"{old}.{ext}").replace(out / f"{new}.{ext}")
    finally:
        shutil.rmtree(scratch)
    return [f"{new}.pdf" for new in names.values()]


def figure5(out, t, manuscript=False):
    """Figure 5 page (cell types): the released panels, drawn by the released figure code and composed here like every
    other figure (author, 2026-11-13: one build, one figure file each). A restyle is pending."""
    def render(scratch):
        terminal_cell_types.plot_retention_heatmap(t.data.cell_type_retention, out_dir=scratch)
        within = t.data.within_type
        terminal_within_type.plot_within_type(within["fronts"], within["summary"], within["aggregate"],
                                              within["endpoints"], out_dir=scratch, display_mode="endpoint")
    a, b, c = released_panels(out, render, {"fig6A_ce_retention_heatmap": "fig5A_cell_type_retention",
                                            "fig6B_ce_within_type_fronts": "fig5B_within_type_fronts",
                                            "fig6C_ce_type_restricted_aggregate": "fig5C_type_restricted_aggregate"})
    return page(out, "fig5_cell_types", "Figure 5",
        "{\\centering\n"
        f"\\includegraphics[width=.95\\textwidth]{{{a}}}\\\\[1mm]\n"
        f"\\includegraphics[width=.95\\textwidth]{{{b}}}\\\\[1mm]\n"
        f"\\includegraphics[width=.95\\textwidth]{{{c}}}\\par}}", manuscript=manuscript)


def figure_s1(out, t, manuscript=False):
    """Supplementary Figure S1 page (terminal assignments in four tracking geometries): the released panel (released
    S3), drawn by the released figure code."""
    (a,) = released_panels(out, lambda scratch: terminal_s3.render_publication(
        t.data.s3_projections, t.data.s3_references, out_dir=scratch),
        {"figS3_four_geometries_with_insets": "figS1_four_geometries"})
    return page(out, "figS1_tracking_robustness", "Figure S1", f"\\noindent\\includegraphics[width=\\textwidth]{{{a}}}",
                manuscript=manuscript)


def table1(out, t):
    """Table 1 (subtree statistics) as a standalone captioned TeX table; its caption cites the subtree figure (4)."""
    return tables.write_subtree_statistics(out, t.data.canonical, t.data.subtree_summary, t.data.min_cells,
                                           stem="table1_subtree_statistics", subtree_figure="4").name


def figure_s2(out, manuscript=False):
    """Supplementary Figure S2 page (cell-state features, Track A); callable alone."""
    a = save(panels.feature_figure(data.feature_selection()), out, "figS2_cell_state_features_panels")
    return page(out, "figS2_cell_state_features", "Figure S2", f"\\noindent\\includegraphics[width=\\textwidth]{{{a}}}",
                manuscript=manuscript)


def figure7(out, manuscript=False):
    """Figure 7 page (C. briggsae): the terminal front alone (author, 2026-11-10). Its C. elegans counterparts,
    lineage similarity along the front (B) and the partial forest by layer (C), are the backlog page that
    ``briggsae_counterparts`` builds. Returns both stems."""
    cs = data.briggsae()
    curve, metric = cs.select("raw3d", "cb_rna")
    curve = curve.reset_index(drop=True)
    cloud = cs.cloud("raw3d", "cb_rna")
    pick = lambda f: (cloud[cloud.family == f].D1.to_numpy(), cloud[cloud.family == f].D2.to_numpy())
    clouds = {f: pick(f) for f in panels.CLOUD_LABELS}
    a = save(panels.front_figure((curve.D1.to_numpy(), curve.D2.to_numpy()), (metric.natural_D1, metric.natural_D2),
                                 clouds, pick("full_random"), lineage_label="C. briggsae lineage",
                                 label_positions=BRIGGSAE_LABELS),
             out, "fig7_briggsae_front")
    main = page(out, "fig7_briggsae", "Figure 7", f"\\noindent\\includegraphics[width=.78\\textwidth]{{{a}}}",
                manuscript=manuscript)
    return main, briggsae_counterparts(out, curve, metric, manuscript=manuscript)


def briggsae_counterparts(out, curve, metric, manuscript=False):
    """Backlog page: C. briggsae counterparts of Figures 3A and 6B (the former 8B-C), kept for reviewers."""
    maximum = np.flatnonzero(curve.sweep_index.to_numpy() == int(metric.maximum_index))
    if len(maximum) != 1 or not np.isclose(curve.edge_retention[maximum[0]], metric.max_edge_retention):
        raise ValueError("C. briggsae maximum-similarity point does not match its metrics")
    b = save(panels.similarity_figure(curve.D1.to_numpy(), curve.edge_retention.to_numpy(), metric.natural_D1,
                                      travel_index=int(curve.alpha_travel.idxmax()),
                                      state_index=int(curve.alpha_travel.idxmin()), maximum_index=int(maximum[0])),
             out, "backlog_fig7B_briggsae_similarity")
    layers = data.briggsae_layers()
    c = save(panels.layer_position(layers, dmax=layer_ceiling(data.layer_summary().layers, layers)), out,
             "backlog_fig7C_briggsae_layers")
    return page(out, "backlog_fig7BC_briggsae_counterparts", "Backlog: Figure 7B-C",
        f"\\noindent\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{B}}{{\\includegraphics[width=.9\\textwidth]{{{b}}}}}\\end{{minipage}}"
        f"\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{C}}{{\\includegraphics[width=.9\\textwidth]{{{c}}}}}\\end{{minipage}}",
        note="Kept in reserve: C. briggsae against C. elegans.", manuscript=manuscript)


def presentation_sources() -> dict[str, str]:
    """Hashes of the publication sources the revised layout renders from: the v2 build, its captions and the manuscript
    set, with every publication module they import. The build
    records them; ``publication.v2.manuscript`` refuses a build, and reports a set as stale, when they changed."""
    files = provenance.module_closure(["publication.v2.__main__", "publication.v2.captions", "publication.v2.manuscript"])
    return {str(p.relative_to(ROOT)): sha256(p) for p in sorted({*files, *ASSETS.rglob("*")}) if p.is_file()}


def figure_order(name):
    """Main figures by number, then supplementary figures (``fig7`` before ``figS2``)."""
    number = name[3:]
    return (number.startswith("S"), int(number.lstrip("S")))


def build(out: Path, manuscript=False):
    """Revised layout, numbering of 2026-11-10: 1 schematics, 2 terminal front, 3 lineage similarity, 4 subtrees,
    5 cell types, 6 full tree by layer, 7 C. briggsae front; S1 terminal assignments in four tracking geometries,
    S2 the cell-state features; Table 1; backlog: C. briggsae 7B-C."""
    out.mkdir(parents=True, exist_ok=True)
    panels.configure()
    for name, source in EXTERNAL.items():
        shutil.copy2(source, out / name)
    shutil.copy2(ASSETS / "null_models.tex", out / "fig2B_null_models.tex")
    t = data.terminal()
    display = t.data.global_display
    clouds = {f: display["nulls"][k] for k, f in ((1, "first_cousin"), (2, "second_cousin"), (3, "third_cousin"))}
    figures = {}

    # Figure 2: terminal Pareto front (A) and null models (B).
    a = save(panels.front_figure((display["x"], display["y"]), display["natural"], clouds, display["nulls"]["full"],
                                 lineage_label="C. elegans lineage"), out, "fig2A_terminal_front")
    figures["fig2"] = page(out, "fig2_terminal_pareto", "Figure 2",
        f"\\panel{{A}}{{\\includegraphics[width=.78\\textwidth]{{{a}}}}}\n"
        "\\panel{B}{\\resizebox{.9\\textwidth}{!}{\\input{fig2B_null_models.tex}}}", manuscript=manuscript)

    # Figure 3: lineage similarity along the front (A) and lineage rewiring at T, M, S (B-D).
    a = save(panels.similarity_figure(display["x"], t.retention, display["natural"][0], travel_index=t.travel_optimum,
                                      state_index=t.state_optimum, maximum_index=t.maximum),
             out, "fig3A_lineage_similarity")
    points = [("B", f"Travel optimum {nt.symbol('travel_optimum')}", t.travel_optimum),
              ("C", "Maximum similarity M", t.maximum),
              ("D", f"Cell-state optimum {nt.symbol('cell_state_optimum')}", t.state_optimum)]
    fig, _ = panels.rewiring_figure(t, points)
    b = save(fig, out, "fig3BD_lineage_rewiring")
    for _, name, weight in points:
        data.edge_changes(t, weight).to_csv(out / f"fig3_edge_changes_{name.split()[0].lower()}_{weight}.csv", index=False)
    figures["fig3"] = page(out, "fig3_lineage_similarity", "Figure 3",
        f"\\panel{{A}}{{\\includegraphics[width=.5\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\includegraphics[width=\\textwidth]{{{b}}}", manuscript=manuscript)

    # Figure 4: subtree rows with fate pies, u, and d_LP with d_CP on one axis (2026-11-10), counterparts adjacent.
    a = save(panels.subtree_rows(t.data.canonical_display, t.data.subtree_nodes, notes=t.data.notes), out,
             "fig4_subtree_rows")
    figures["fig4"] = page(out, "fig4_subtrees", "Figure 4",
        f"\\noindent\\includegraphics[width=\\textwidth]{{{a}}}", manuscript=manuscript)

    figures["fig5"] = figure5(out, t, manuscript=manuscript)
    figures["figS1"] = figure_s1(out, t, manuscript=manuscript)
    figures["figS2"] = figure_s2(out, manuscript=manuscript)
    figures["fig6"] = figure6(out, manuscript=manuscript)
    figures["fig7"], backlog = figure7(out, manuscript=manuscript)

    # Figure 1: schematics. Placeholder assembly of the two authors' files.
    figures["fig1"] = page(out, "fig1_overview", "Figure 1",
        "\\noindent\\includegraphics[width=\\textwidth]{Fig1_overview.pdf}\\par\\bigskip\n"
        "\\panel{D}{\\includegraphics[width=.62\\textwidth]{FigX_distances_schematic.pdf}}",
        note="Authors' files, unchanged: B-C still read \\emph{expression}.", manuscript=manuscript)
    palette_sheet = save(palette.sheet(), out, "palette")
    record = dict(figures=dict(sorted(figures.items(), key=lambda kv: figure_order(kv[0]))),
                  table1=table1(out, t), backlog=dict(fig7BC=backlog), notation=NOTATION,
                  manuscript=manuscript, palette=palette_sheet, palette_groups=palette.GROUPS,
                  main_text_notes=dict(tracking_front_distance=data.tracking_front_distances(t)),
                  external={name: sha256(out / name) for name in EXTERNAL}, sources=presentation_sources(),
                  files={p.name: sha256(p) for p in sorted(out.iterdir())
                         if p.suffix in (".pdf", ".png", ".csv") or p.name.startswith("table1")})
    (out / "build.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manuscript", action="store_true",
                        help="figures alone, cropped, without the draft title and note (publication.v2.manuscript)")
    args = parser.parse_args(argv)
    with nt.override(**NOTATION):
        record = build(args.output_dir, manuscript=args.manuscript)
    print(json.dumps(record["figures"], indent=2))


if __name__ == "__main__":
    main()
