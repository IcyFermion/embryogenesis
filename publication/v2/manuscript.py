"""Assemble the manuscript figure set: ``python -m publication.v2.manuscript --v2-build DIR --output-dir DIR``.

Numbering of 2026-11-08. Main Figures 1-5, 7 and 8 come from a v2 build made
with ``--manuscript`` (figures alone, cropped). Figure 6 holds the released
cell-type panels (production Figure 6) until they are restyled. The supplement
is released S3, renumbered S1; released S1 is retired (its r is in Figure 5) and
released S2, S4 and S5 are kept in ``backlog/`` with their captions, outside
``review.pdf``. Table 1 is the released page. Production must verify cleanly first.

Every figure file is a captioned page (``Figure1.pdf`` ... ``FigureS1.pdf``,
``Table1.pdf``, plus ``Table1.tex``). Main-figure captions are agreed or, until
the author approves them, drafts in the v1 style (``publication.v2.captions``).
``graphics/`` holds the main figures' bare graphics for LaTeX ``includegraphics``.
``review.tex`` is the single source of ``review.pdf`` (``tectonic review.tex`` in
the folder recompiles it); the main figures' pages are split out of it.
``captions.tex``, ``README.md`` and ``manifest.json`` complete the folder.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from publication import notation as nt
from publication.provenance import ROOT, sha256
from publication.release import PRODUCTION, verify_production
from publication.v2 import captions
from publication.v2.__main__ import MANUSCRIPT_PREAMBLE, NOTATION

MAIN = {"1": "fig1_overview.pdf", "2": "fig2_cell_state_features.pdf", "3": "fig3_terminal_pareto.pdf",
        "4": "fig4_lineage_similarity.pdf", "5": "fig5_subtrees.pdf", "7": "fig7_full_tree_layers.pdf",
        "8": "fig8_briggsae.pdf"}
# Cell types until restyled: the released cell-type figure (production Figure 6), panels and caption.
CELL_TYPES = dict(number="6", wrapper="fig6_ce_cell_types.tex", panels=("fig6A_ce_retention_heatmap.pdf",
                                                                       "fig6B_ce_within_type_fronts.pdf",
                                                                       "fig6C_ce_type_restricted_aggregate.pdf"))
# Supplement of the manuscript set (new number -> released stem); author decisions of 2026-11-08.
SUPPLEMENT = {"S1": "figS3_ce_tracking_robustness"}
BACKLOG = {"S2": "figS2_ce_cell_type_cost_gain", "S4": "figs_ce_full_tree_heuristics",
           "S5": "figS5_full_tree_layer_travel"}
RETIRED = {"S1": ("figS1_ce_canonical_summary_cousin_r", "r is a column of Figure 5; panel A repeated Figure 1D")}
TABLE1 = "table1_ce_subtree_statistics"
# Released figures that the revised layout does not place (listed in the README, not copied).
NOT_PLACED = {"7": "full-tree layerwise rounds and aggregate front (7A is redrawn as Figure 7A; 7B is not placed)",
              "8": "terminal cross-species comparison (its C. briggsae front is redrawn as Figure 8A)",
              "9": "terminal cross-species overlays",
              "10": "partial-forest cross-species comparison (Figure 8C draws its C. briggsae layers)",
              "11": "partial-forest cross-species overlays"}
# Released captions that cite the released numbering.
OLD_REFERENCES = {}
REVIEW_PREAMBLE = r"""% Review copy of the manuscript figure set. Recompile in this folder: tectonic review.tex
% (or pdflatex review.tex). Main-figure captions can be edited here; the supplement and
% Table 1 are the released pages, edited in their production sources.
\documentclass[10pt]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{amsmath}
\usepackage{pdfpages}
\captionsetup{font=small,labelfont=bf}
\pagestyle{empty}
\begin{document}
"""


def caption_body(tex: str) -> str:
    """The argument of the first ``\\caption{...}`` in a wrapper (brace matched)."""
    start = tex.index(r"\caption{") + len(r"\caption{")
    depth = 1
    for i in range(start, len(tex)):
        depth += {"{": 1, "}": -1}.get(tex[i], 0) if tex[i - 1] != "\\" else 0
        if depth == 0:
            return tex[start:i].strip()
    raise ValueError("Unbalanced caption")


def tectonic(tex: Path) -> Path:
    result = subprocess.run(["tectonic", tex.name], cwd=tex.parent, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{tex.name} failed to compile:\n{result.stdout[-3000:]}")
    return tex.with_suffix(".pdf")


def cell_types(work: Path) -> Path:
    for name in CELL_TYPES["panels"]:
        shutil.copy2(PRODUCTION / name, work / name)
    body = "\\\\[1mm]\n".join(f"\\includegraphics[width=.95\\textwidth]{{{name}}}" for name in CELL_TYPES["panels"])
    tex = work / "cell_types.tex"
    tex.write_text(MANUSCRIPT_PREAMBLE + "\\noindent\\centering\n" + body + "\n\\end{document}\n")
    return tectonic(tex)


def git_state() -> dict:
    run = lambda *args: subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True).stdout.strip()
    return dict(commit=run("rev-parse", "HEAD"), branch=run("branch", "--show-current"),
                dirty=bool(run("status", "--porcelain", "--untracked-files=no")))


def notes(record: dict) -> list[str]:
    """Numbers the main text should carry, recorded by the v2 build (``main_text_notes``)."""
    d = record["main_text_notes"]["tracking_front_distance"]
    embryos = [d[f"embryo{i}"] for i in (1, 2, 3)]
    return [
        "- **Supplementary Figure S1 (released S3), d_LP from pooled tracking.** The natural lineage's distance to "
        "its own front, each in that geometry's endpoint coordinates: embryo 1 "
        + ", embryo 2 ".join(f"{v:.3f}" for v in embryos[:2]) + f", embryo 3 {embryos[2]:.3f} "
        f"(mean {np.mean(embryos):.3f}); pooled travel {d['pooled']:.3f} (the global d_LP of Figure 3). Pooling "
        "gives a smaller distance than two of the three single embryos and than their mean, but not the smallest "
        "(embryo 3); the text should say the pooled estimate averages out single-embryo tracking variation, not "
        "that it is closest to the front.",
    ]


def assemble(v2_build: Path, out: Path) -> dict:
    v2_build, out = Path(v2_build), Path(out)
    record = json.loads((v2_build / "build.json").read_text())
    if not record.get("manuscript"):
        raise ValueError(f"{v2_build} is a draft build; rebuild it with `python -m publication.v2 --manuscript`")
    production = verify_production()
    stale = {k: v for k, v in production.get("stale_presentation_sources", {}).items() if v}
    if stale:
        raise ValueError(f"Production is not current for: {sorted(stale)}; release a fresh build first")
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory must be new or empty: {out}")
    work, graphics = out / "_tex", out / "graphics"
    work.mkdir(parents=True)
    graphics.mkdir()

    cell = CELL_TYPES["number"]
    bare = {f"Figure{n}.pdf": v2_build / name for n, name in MAIN.items()} | {f"Figure{cell}.pdf": cell_types(work)}
    for name, source in bare.items():
        shutil.copy2(source, graphics / name)
    released_pages = {f"Figure{n}.pdf": PRODUCTION / f"{stem}.pdf" for n, stem in SUPPLEMENT.items()}
    released_pages |= {"Table1.pdf": PRODUCTION / f"{TABLE1}.pdf", "Table1.tex": PRODUCTION / f"{TABLE1}.tex"}
    for name, source in released_pages.items():
        shutil.copy2(source, out / name)
    (out / "backlog").mkdir()
    for n, stem in BACKLOG.items():  # kept in reserve, under their released numbers
        shutil.copy2(PRODUCTION / f"{stem}.pdf", out / "backlog" / f"Figure{n}_released.pdf")

    with nt.override(**NOTATION):
        agreed = {n: make() for n, make in captions.AGREED.items()}
        drafts = {n: make() for n, make in captions.DRAFT.items()}
    main_captions = {n: agreed.get(n) or drafts.get(n) for n in sorted([*MAIN, cell], key=int)}
    main_captions[cell] = caption_body((PRODUCTION / CELL_TYPES["wrapper"]).read_text())
    status = {n: "agreed" if n in agreed else "released" if n == cell else "draft" if n in drafts else "missing"
              for n in main_captions}
    if "missing" in status.values():
        raise ValueError(f"No caption for main figures {[n for n, v in status.items() if v == 'missing']}")
    released = {n: caption_body((PRODUCTION / f"{stem}.tex").read_text()) for n, stem in SUPPLEMENT.items()}
    released["Table 1"] = caption_body((PRODUCTION / f"{TABLE1}.tex").read_text())
    backlog = {n: caption_body((PRODUCTION / f"{stem}.tex").read_text()) for n, stem in BACKLOG.items()}

    blocks = [f"% Figure {n} ({status[n]})\n\\setcounter{{figure}}{{{int(n) - 1}}}\n\\begin{{figure}}[p]\n\\centering\n"
              f"\\includegraphics[width=\\textwidth,height=.7\\textheight,keepaspectratio]{{graphics/Figure{n}.pdf}}\n"
              f"\\caption{{{text}}}\n\\end{{figure}}\n\\clearpage\n" for n, text in main_captions.items()]
    pages = [f"Figure{n}.pdf" for n in SUPPLEMENT] + ["Table1.pdf"]
    blocks.append("% Supplement and Table 1: released pages\n"
                  + "".join(f"\\includepdf[pages=-]{{{name}}}\n" for name in pages))
    (out / "review.tex").write_text(REVIEW_PREAMBLE + "\n".join(blocks) + "\\end{document}\n")
    tectonic(out / "review.tex")
    for page, n in enumerate(main_captions, start=1):  # each main figure's captioned page, as in review.pdf
        subprocess.run(["qpdf", str(out / "review.pdf"), "--pages", ".", str(page), "--", str(out / f"Figure{n}.pdf")],
                       check=True)

    blocks = ["% Captions for the manuscript figure set (TeX caption bodies).",
              "% Main figures: agreed captions, or drafts in the v1 style awaiting the author's review.",
              f"% Figure {cell} and the supplement: released captions, unchanged (S1 is released S3).",
              "% Backlog (released S2, S4, S5) follows at the end, under the released numbers."]
    for n, text in main_captions.items():
        label = "released (production Figure 6)" if n == cell else status[n]
        blocks.append(f"\n% Figure {n} ({label})\n{text}")
    for n, text in released.items():
        source = f", released {SUPPLEMENT[n].split('_')[0][3:]}" if n in SUPPLEMENT else ""
        blocks.append(f"\n% {'Figure ' + n if n != 'Table 1' else n} (released{source})\n{text}")
    for n, text in backlog.items():
        blocks.append(f"\n% Backlog: released Figure {n}\n{text}")
    (out / "captions.tex").write_text("\n".join(blocks) + "\n")

    shutil.rmtree(work)
    families = json.loads((PRODUCTION / "release_manifest.json").read_text())["families"]
    readme = [
        "# Manuscript figure set (revised layout)", "",
        f"Built by `python -m publication.v2.manuscript` from `{v2_build.name}` (v2 manuscript-mode build) and "
        f"production released on {max(f['released_at'] for f in families.values())[:8]}. "
        "Every figure file is a page with its caption; `review.pdf` joins them all and is compiled from "
        "`review.tex` (`tectonic review.tex` here recompiles it). `graphics/` holds the main figures without "
        "captions, for `\\includegraphics`, and `captions.tex` holds every caption as TeX.", "",
        "| File | Figure | Source | Caption |", "|---|---|---|---|",
        *(f"| Figure{n}.pdf | {n} | revised layout (v2) | {status[n]} |"
          if n != cell else
          f"| Figure{cell}.pdf | {cell} | released Figure 6 panels (cell types), not yet restyled | released, renumbered |"
          for n in main_captions),
        *(f"| Figure{n}.pdf | {n} | released `{stem}` | released, renumbered |" for n, stem in SUPPLEMENT.items()),
        f"| Table1.pdf, Table1.tex | Table 1 | released `{TABLE1}` | released |", "",
        "## Notes for the main text", "",
        *notes(record), "",
        "## Open items", "",
        "- Figure 1 is the authors' schematic file, unchanged: panels B-C still read \"expression\", and 1D writes "
        "d_LP*/d_CP* where the other figures write d_LP/d_CP.",
        f"- Figure {cell} is a placeholder until the cell-type figure is restyled; its caption is the released one.",
        "- Captions of Figures " + ", ".join(n for n, v in status.items() if v == "draft")
        + " are drafts in the v1 style, for the author's review (`publication/v2/captions.py`).",
        *(f"- Released caption of {n} {why}." for n, why in OLD_REFERENCES.items()), "",
        "## Backlog (kept in reserve, not in review.pdf)", "",
        *(f"- `backlog/Figure{n}_released.pdf`: released `{stem}`; caption at the end of `captions.tex`."
          for n, stem in BACKLOG.items()),
        "- Terminal cohort sensitivity (Track D, D7; author kept the 275-edge cohort on 2026-11-09): 46 cells end "
        "embryogenesis undivided but are internal in the reference lineage (10 E nuclei, which is why Figure 4 "
        "shows two E cells; 11 each in ABar, ABpl, ABpr; 3 in MS). With them (321 edges) the global d_LP is 0.031 "
        "instead of 0.015, u 0.204, maximum lineage similarity 0.75; the lineage stays far better than its "
        "first-cousin null (travel z -9.9, cell-state z -5.1). Run "
        "`terminal_pareto/output/runs/cohort_sensitivity_v1/embryonic_terminals_20261108`.", "",
        "## Retired", "",
        *(f"- Released {n} (`{stem}`): {why}." for n, (stem, why) in RETIRED.items()), "",
        "## Released figures not placed in the revised layout", "",
        *(f"- Released Figure {n}: {what}." for n, what in NOT_PLACED.items()), "",
    ]
    (out / "README.md").write_text("\n".join(readme))
    manifest = dict(version="manuscript-set-2", git=git_state(), v2_build=str(v2_build),
                    v2_build_record=sha256(v2_build / "build.json"),
                    production={k: dict(released_at=v["released_at"], presentation_id=v["presentation_id"])
                                for k, v in families.items()},
                    captions={f"Figure {n}": v for n, v in status.items()},
                    files={str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--v2-build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = assemble(args.v2_build, args.output_dir)
    print(json.dumps(dict(captions=manifest["captions"], files=len(manifest["files"])), indent=2))


if __name__ == "__main__":
    main()
