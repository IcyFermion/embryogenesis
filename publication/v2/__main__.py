"""Build the revised manuscript layout: ``python -m publication.v2 --output-dir DIR``.

Draft pages carry panels and panel letters only; captions follow once the
layout is approved. Figure 1 and Figure 4A are the authors' OmniGraffle files
from ``paper/FIGS/``, included unchanged.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import matplotlib.pyplot as plt

from publication import notation as nt
from publication.provenance import sha256
from publication.v2 import data, panels

ROOT = Path(__file__).resolve().parents[2]
ASSETS = Path(__file__).resolve().parent / "assets"
EXTERNAL = {"Fig1_overview.pdf": ROOT / "paper/FIGS/Fig1_overview.pdf",
            "FigX_distances_schematic.pdf": ROOT / "paper/FIGS/FigX_distances_schematic.pdf"}
# Collaborator suggestion adopted for this candidate (2026-10-08); one-line switch in notation.py on adoption.
NOTATION = dict(retention=dict(label="Lineage similarity"))
# Figure 7A clouds sit side by side, so its labels are placed by hand (data coordinates).
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


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return f"{stem}.pdf"


def page(out, stem, title, body, note="Draft layout; caption to follow."):
    tex = out / f"{stem}.tex"
    tex.write_text(PREAMBLE + f"\\noindent{{\\sffamily\\bfseries {title}}}\\hfill{{\\sffamily\\small\\color{{gray}}{note}}}"
                   "\\par\\bigskip\n" + body + "\n\\end{document}\n")
    result = subprocess.run(["tectonic", tex.name], cwd=out, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{tex.name} failed to compile:\n{result.stdout[-3000:]}")
    return f"{stem}.pdf"


def build(out: Path):
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
        "\\panel{B}{\\resizebox{.9\\textwidth}{!}{\\input{fig2B_null_models.tex}}}")

    # Figure 3: lineage similarity along the front (A) and lineage rewiring at T, M, S (B-D).
    fig, ax = plt.subplots(figsize=(3.6, 2.7))
    panels.similarity_curve(ax, display["x"], t.retention, display["natural"][0], travel_index=t.travel_optimum,
                            state_index=t.state_optimum, maximum_index=t.maximum)
    fig.tight_layout()
    a = save(fig, out, "fig3A_lineage_similarity")
    points = [("B", f"Travel optimum {nt.symbol('travel_optimum')}", t.travel_optimum),
              ("C", "Maximum similarity M", t.maximum),
              ("D", f"Cell-state optimum {nt.symbol('cell_state_optimum')}", t.state_optimum)]
    fig, _ = panels.rewiring_figure(t, points)
    b = save(fig, out, "fig3BD_lineage_rewiring")
    for _, name, weight in points:
        data.edge_changes(t, weight).to_csv(out / f"fig3_edge_changes_{name.split()[0].lower()}_{weight}.csv", index=False)
    figures["fig3"] = page(out, "fig3_lineage_similarity", "Figure 3",
        f"\\panel{{A}}{{\\includegraphics[width=.5\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\includegraphics[width=\\textwidth]{{{b}}}")

    # Figure 4: distance schematic (A) and subtree rows with fate pies (B).
    b = save(panels.subtree_rows(t.data.canonical_display, t.data.subtree_nodes), out, "fig4B_subtree_rows")
    figures["fig4"] = page(out, "fig4_subtrees", "Figure 4",
        "\\panel{A}{\\includegraphics[width=.62\\textwidth]{FigX_distances_schematic.pdf}}\n"
        f"\\panel{{B}}{{\\includegraphics[width=.95\\textwidth]{{{b}}}}}")

    # Figure 6: full tree by layer.
    by_layer, layer_metrics = data.layers()
    a = save(panels.layer_grid(by_layer), out, "fig6A_layer_fronts")
    b = save(panels.layer_similarity(by_layer), out, "fig6B_layer_similarity")
    c = save(panels.layer_rows(layer_metrics), out, "fig6C_layer_metrics")
    figures["fig6"] = page(out, "fig6_full_tree_layers", "Figure 6",
        f"\\panel{{A}}{{\\includegraphics[width=.95\\textwidth]{{{a}}}}}\n"
        f"\\noindent\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{B}}{{\\includegraphics[width=.9\\textwidth]{{{b}}}}}\\end{{minipage}}"
        f"\\begin{{minipage}}[t]{{.5\\textwidth}}\\panel{{C}}{{\\includegraphics[width=.9\\textwidth]{{{c}}}}}\\end{{minipage}}")

    # Figure 7: C. briggsae terminal front (A); B to follow.
    cs = data.briggsae()
    curve, metric = cs.select("raw3d", "cb_rna")
    cloud = cs.cloud("raw3d", "cb_rna")
    pick = lambda f: (cloud[cloud.family == f].D1.to_numpy(), cloud[cloud.family == f].D2.to_numpy())
    a = save(panels.front_figure((curve.D1.to_numpy(), curve.D2.to_numpy()), (metric.natural_D1, metric.natural_D2),
                                 {f: pick(f) for f in clouds}, pick("full_random"), lineage_label="C. briggsae lineage",
                                 label_positions=BRIGGSAE_LABELS),
             out, "fig7A_briggsae_front")
    figures["fig7"] = page(out, "fig7_briggsae", "Figure 7",
        f"\\panel{{A}}{{\\includegraphics[width=.78\\textwidth]{{{a}}}}}\n\\panel{{B}}{{\\sffamily\\small To follow.}}")

    figures["fig1"] = page(out, "fig1_overview", "Figure 1",
        "\\noindent\\includegraphics[width=\\textwidth]{Fig1_overview.pdf}",
        note="Authors' file, unchanged: panels B-C still read \\emph{expression}.")
    record = dict(figures=figures, notation=NOTATION,
                  external={name: sha256(out / name) for name in EXTERNAL},
                  files={p.name: sha256(p) for p in sorted(out.iterdir()) if p.suffix in (".pdf", ".png", ".csv")})
    (out / "build.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    with nt.override(**NOTATION):
        record = build(args.output_dir)
    print(json.dumps(record["figures"], indent=2))


if __name__ == "__main__":
    main()
