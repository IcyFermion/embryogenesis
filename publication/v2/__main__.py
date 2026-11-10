"""Build the revised manuscript layout: ``python -m publication.v2 --output-dir DIR``.

Draft pages carry panels and panel letters only (captions: ``publication.v2.captions``).
Figure 1 is the schematic figure: panels A-C and D are the authors' OmniGraffle
files from ``paper/FIGS/``, included unchanged. Colours come from
``publication.palette``; the build also draws the palette sheet (``palette.pdf``).
Numbering of 2026-11-08: Figure 2 (cell-state features) is new, and the former
Figures 2-7 are now 3-8.
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
from publication import palette
from publication.provenance import sha256
from publication.v2 import data, panels

ROOT = Path(__file__).resolve().parents[2]
ASSETS = Path(__file__).resolve().parent / "assets"
EXTERNAL = {"Fig1_overview.pdf": ROOT / "paper/FIGS/Fig1_overview.pdf",
            "FigX_distances_schematic.pdf": ROOT / "paper/FIGS/FigX_distances_schematic.pdf"}
# Collaborator suggestion adopted for this candidate (2026-10-08); one-line switch in notation.py on adoption.
NOTATION = dict(retention=dict(label="Lineage similarity"))
# Figure 8A clouds sit side by side, so its labels are placed by hand (data coordinates).
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
    """Shared d_LP axis of Figures 7B and 8C, so the two species compare directly."""
    return max(.3, np.ceil(max(float(t.front_distance.max()) for t in tables) * 20) / 20 + .05)


def figure7(out, manuscript=False):
    """Figure 7 page (full tree by layer); callable alone (the other figures need the terminal analysis context)."""
    # A: every layer's front (layers 1-2 show why the summaries start at layer 3).
    # B-C: Track E's per-layer summaries on the birth-size axis (author decision 2026-11-05), drawn with the palette.
    summary = data.layer_summary().layers
    a = save(panels.layer_grid(data.layers()), out, "fig7A_layer_fronts")
    b = save(panels.layer_position(summary, dmax=layer_ceiling(summary, data.briggsae_layers())), out,
             "fig7B_layer_position")
    c = save(panels.layer_cost(summary), out, "fig7C_layer_cost")
    return page(out, "fig7_full_tree_layers", "Figure 7",
        f"\\panel{{A}}{{\\includegraphics[width=.95\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{B}}{{\\includegraphics[width=.9\\textwidth]{{{b}}}}}\\end{{minipage}}"
        f"\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{C}}{{\\includegraphics[width=.9\\textwidth]{{{c}}}}}\\end{{minipage}}",
        manuscript=manuscript)


def figure2(out, manuscript=False):
    """Figure 2 page (cell-state features, Track A); callable alone."""
    a = save(panels.feature_figure(data.feature_selection()), out, "fig2_cell_state_features_panels")
    return page(out, "fig2_cell_state_features", "Figure 2", f"\\noindent\\includegraphics[width=\\textwidth]{{{a}}}",
                manuscript=manuscript)


def figure8(out, manuscript=False):
    """Figure 8 page (C. briggsae): A terminal front; B lineage similarity along it; C the partial forest by layer."""
    cs = data.briggsae()
    curve, metric = cs.select("raw3d", "cb_rna")
    curve = curve.reset_index(drop=True)
    cloud = cs.cloud("raw3d", "cb_rna")
    pick = lambda f: (cloud[cloud.family == f].D1.to_numpy(), cloud[cloud.family == f].D2.to_numpy())
    clouds = {f: pick(f) for f in panels.CLOUD_LABELS}
    a = save(panels.front_figure((curve.D1.to_numpy(), curve.D2.to_numpy()), (metric.natural_D1, metric.natural_D2),
                                 clouds, pick("full_random"), lineage_label="C. briggsae lineage",
                                 label_positions=BRIGGSAE_LABELS),
             out, "fig8A_briggsae_front")
    maximum = np.flatnonzero(curve.sweep_index.to_numpy() == int(metric.maximum_index))
    if len(maximum) != 1 or not np.isclose(curve.edge_retention[maximum[0]], metric.max_edge_retention):
        raise ValueError("C. briggsae maximum-similarity point does not match its metrics")
    b = save(panels.similarity_figure(curve.D1.to_numpy(), curve.edge_retention.to_numpy(), metric.natural_D1,
                                      travel_index=int(curve.alpha_travel.idxmax()),
                                      state_index=int(curve.alpha_travel.idxmin()), maximum_index=int(maximum[0])),
             out, "fig8B_briggsae_similarity")
    layers = data.briggsae_layers()
    c = save(panels.layer_position(layers, dmax=layer_ceiling(data.layer_summary().layers, layers)), out,
             "fig8C_briggsae_layers")
    return page(out, "fig8_briggsae", "Figure 8",
        f"\\panel{{A}}{{\\includegraphics[width=.78\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{B}}{{\\includegraphics[width=.9\\textwidth]{{{b}}}}}\\end{{minipage}}"
        f"\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{C}}{{\\includegraphics[width=.9\\textwidth]{{{c}}}}}\\end{{minipage}}",
        manuscript=manuscript)


def build(out: Path, manuscript=False):
    """Revised layout, numbering of 2026-11-08: 1 schematics, 2 cell-state features (new), 3 terminal front,
    4 lineage similarity, 5 subtrees, 6 cell types (released panels, placed by ``manuscript``), 7 full tree by
    layer, 8 C. briggsae."""
    out.mkdir(parents=True, exist_ok=True)
    panels.configure()
    for name, source in EXTERNAL.items():
        shutil.copy2(source, out / name)
    shutil.copy2(ASSETS / "null_models.tex", out / "fig3B_null_models.tex")
    t = data.terminal()
    display = t.data.global_display
    clouds = {f: display["nulls"][k] for k, f in ((1, "first_cousin"), (2, "second_cousin"), (3, "third_cousin"))}
    figures = {}

    # Figure 3: terminal Pareto front (A) and null models (B).
    a = save(panels.front_figure((display["x"], display["y"]), display["natural"], clouds, display["nulls"]["full"],
                                 lineage_label="C. elegans lineage"), out, "fig3A_terminal_front")
    figures["fig3"] = page(out, "fig3_terminal_pareto", "Figure 3",
        f"\\panel{{A}}{{\\includegraphics[width=.78\\textwidth]{{{a}}}}}\n"
        "\\panel{B}{\\resizebox{.9\\textwidth}{!}{\\input{fig3B_null_models.tex}}}", manuscript=manuscript)

    # Figure 4: lineage similarity along the front (A) and lineage rewiring at T, M, S (B-D).
    a = save(panels.similarity_figure(display["x"], t.retention, display["natural"][0], travel_index=t.travel_optimum,
                                      state_index=t.state_optimum, maximum_index=t.maximum),
             out, "fig4A_lineage_similarity")
    points = [("B", f"Travel optimum {nt.symbol('travel_optimum')}", t.travel_optimum),
              ("C", "Maximum similarity M", t.maximum),
              ("D", f"Cell-state optimum {nt.symbol('cell_state_optimum')}", t.state_optimum)]
    fig, _ = panels.rewiring_figure(t, points)
    b = save(fig, out, "fig4BD_lineage_rewiring")
    for _, name, weight in points:
        data.edge_changes(t, weight).to_csv(out / f"fig4_edge_changes_{name.split()[0].lower()}_{weight}.csv", index=False)
    figures["fig4"] = page(out, "fig4_lineage_similarity", "Figure 4",
        f"\\panel{{A}}{{\\includegraphics[width=.5\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\includegraphics[width=\\textwidth]{{{b}}}", manuscript=manuscript)

    # Figure 5: subtree rows with fate pies, u, d_LP and r (former S1), counterparts adjacent (2026-11-08).
    a = save(panels.subtree_rows(t.data.canonical_display, t.data.subtree_nodes, notes=t.data.notes), out,
             "fig5_subtree_rows")
    figures["fig5"] = page(out, "fig5_subtrees", "Figure 5",
        f"\\noindent\\includegraphics[width=\\textwidth]{{{a}}}", manuscript=manuscript)

    figures["fig2"] = figure2(out, manuscript=manuscript)
    figures["fig7"] = figure7(out, manuscript=manuscript)
    figures["fig8"] = figure8(out, manuscript=manuscript)

    # Figure 1: schematics. Placeholder assembly of the two authors' files.
    figures["fig1"] = page(out, "fig1_overview", "Figure 1",
        "\\noindent\\includegraphics[width=\\textwidth]{Fig1_overview.pdf}\\par\\bigskip\n"
        "\\panel{D}{\\includegraphics[width=.62\\textwidth]{FigX_distances_schematic.pdf}}",
        note="Authors' files, unchanged: B-C still read \\emph{expression}.", manuscript=manuscript)
    palette_sheet = save(palette.sheet(), out, "palette")
    record = dict(figures=dict(sorted(figures.items(), key=lambda kv: int(kv[0][3:]))), notation=NOTATION,
                  manuscript=manuscript, palette=palette_sheet, palette_groups=palette.GROUPS,
                  main_text_notes=dict(tracking_front_distance=data.tracking_front_distances(t)),
                  external={name: sha256(out / name) for name in EXTERNAL},
                  files={p.name: sha256(p) for p in sorted(out.iterdir()) if p.suffix in (".pdf", ".png", ".csv")})
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
