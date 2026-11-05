"""Assemble the manuscript figure set: ``python -m publication.v2.manuscript --v2-build DIR --output-dir DIR``.

Main Figures 1-4, 6 and 7 come from a v2 build made with ``--manuscript``
(figures alone, cropped). Figure 5 holds the released cell-type panels
(production Figure 6) under their new number until they are restyled.
Supplementary Figures S1-S5 and Table 1 are the released production pages,
with their released captions. Production must verify cleanly first.

Every figure file is a captioned page (``Figure1.pdf`` ... ``FigureS5.pdf``,
``Table1.pdf``, plus ``Table1.tex``); main figures without an agreed caption
say "Caption to follow". ``graphics/`` holds the main figures' bare graphics
for LaTeX ``includegraphics``. ``review.pdf`` joins every page; ``captions.tex``,
``README.md`` and ``manifest.json`` complete the folder.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from publication import notation as nt
from publication.provenance import ROOT, sha256
from publication.release import PRODUCTION, verify_production
from publication.v2 import captions
from publication.v2.__main__ import MANUSCRIPT_PREAMBLE, NOTATION

MAIN = {"1": "fig1_overview.pdf", "2": "fig2_terminal_pareto.pdf", "3": "fig3_lineage_similarity.pdf",
        "4": "fig4_subtrees.pdf", "6": "fig6_full_tree_layers.pdf", "7": "fig7_briggsae.pdf"}
# Figure 5 until it is restyled: the released cell-type figure (production Figure 6), panels and caption.
FIGURE5 = dict(wrapper="fig6_ce_cell_types.tex", panels=("fig6A_ce_retention_heatmap.pdf",
                                                         "fig6B_ce_within_type_fronts.pdf",
                                                         "fig6C_ce_type_restricted_aggregate.pdf"))
SUPPLEMENT = {"S1": "figS1_ce_canonical_summary_cousin_r", "S2": "figS2_ce_cell_type_cost_gain",
              "S3": "figS3_ce_tracking_robustness", "S4": "figs_ce_full_tree_heuristics",
              "S5": "figS5_full_tree_layer_travel"}
TABLE1 = "table1_ce_subtree_statistics"
# Released figures that the revised layout does not place (listed in the README, not copied).
NOT_PLACED = {"7": "full-tree layerwise rounds and aggregate front (7A is redrawn as Figure 6A; 7B is not placed)",
              "8": "terminal cross-species comparison", "9": "terminal cross-species overlays",
              "10": "partial-forest cross-species comparison", "11": "partial-forest cross-species overlays"}
# Released captions that cite the released numbering.
OLD_REFERENCES = {"S4": "cites Figure 7 and 7B (released full-tree figure; 7B is not placed)",
                  "Table 1": "cites Figure 5 (released canonical summary; its subtree rows are now Figure 4)"}
REVIEW_PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{amsmath}
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


def figure5(work: Path) -> Path:
    for name in FIGURE5["panels"]:
        shutil.copy2(PRODUCTION / name, work / name)
    body = "\\\\[1mm]\n".join(f"\\includegraphics[width=.95\\textwidth]{{{name}}}" for name in FIGURE5["panels"])
    tex = work / "figure5.tex"
    tex.write_text(MANUSCRIPT_PREAMBLE + "\\noindent\\centering\n" + body + "\n\\end{document}\n")
    return tectonic(tex)


def git_state() -> dict:
    run = lambda *args: subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True).stdout.strip()
    return dict(commit=run("rev-parse", "HEAD"), branch=run("branch", "--show-current"),
                dirty=bool(run("status", "--porcelain", "--untracked-files=no")))


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

    bare = {f"Figure{n}.pdf": v2_build / name for n, name in MAIN.items()} | {"Figure5.pdf": figure5(work)}
    for name, source in bare.items():
        shutil.copy2(source, graphics / name)
    released_pages = {f"Figure{n}.pdf": PRODUCTION / f"{stem}.pdf" for n, stem in SUPPLEMENT.items()}
    released_pages |= {"Table1.pdf": PRODUCTION / f"{TABLE1}.pdf", "Table1.tex": PRODUCTION / f"{TABLE1}.tex"}
    for name, source in released_pages.items():
        shutil.copy2(source, out / name)

    with nt.override(**NOTATION):
        agreed = {n: make() for n, make in captions.AGREED.items()}
    main_captions = {n: agreed.get(n) for n in sorted([*MAIN, "5"], key=int)}
    main_captions["5"] = caption_body((PRODUCTION / FIGURE5["wrapper"]).read_text())
    released = {n: caption_body((PRODUCTION / f"{stem}.tex").read_text()) for n, stem in SUPPLEMENT.items()}
    released["Table 1"] = caption_body((PRODUCTION / f"{TABLE1}.tex").read_text())

    for n, text in main_captions.items():  # one captioned page per main figure, like the released pages
        shutil.copy2(graphics / f"Figure{n}.pdf", work / f"graphic{n}.pdf")
        text = text or r"{\itshape Caption to follow.}"
        tex = work / f"Figure{n}.tex"
        tex.write_text(REVIEW_PREAMBLE + f"\\setcounter{{figure}}{{{int(n) - 1}}}\n"
                       "\\begin{figure}[p]\n\\centering\n"
                       f"\\includegraphics[width=\\textwidth,height=.7\\textheight,keepaspectratio]{{graphic{n}.pdf}}\n"
                       f"\\caption{{{text}}}\n\\end{{figure}}\n\\end{{document}}\n")
        shutil.copy2(tectonic(tex), out / f"Figure{n}.pdf")
    order = [f"Figure{n}.pdf" for n in main_captions] + [f"Figure{n}.pdf" for n in SUPPLEMENT] + ["Table1.pdf"]
    subprocess.run(["pdfunite", *(str(out / name) for name in order), str(out / "review.pdf")], check=True)

    blocks = ["% Captions for the manuscript figure set (TeX caption bodies).",
              "% Main figures: agreed captions; 'to follow' where none is agreed yet.",
              "% Figure 5 and the supplement: released captions, unchanged."]
    for n, text in main_captions.items():
        status = "agreed" if n in agreed else ("released (production Figure 6)" if n == "5" else "to follow")
        blocks.append(f"\n% Figure {n} ({status})\n{text or '% to follow'}")
    for n, text in released.items():
        blocks.append(f"\n% {'Figure ' + n if n != 'Table 1' else n} (released)\n{text}")
    (out / "captions.tex").write_text("\n".join(blocks) + "\n")

    shutil.rmtree(work)
    families = json.loads((PRODUCTION / "release_manifest.json").read_text())["families"]
    readme = [
        "# Manuscript figure set (revised layout)", "",
        f"Built by `python -m publication.v2.manuscript` from `{v2_build.name}` (v2 manuscript-mode build) and "
        f"production released on {max(f['released_at'] for f in families.values())[:8]}. "
        "Every figure file is a page with its caption; `review.pdf` joins them all. `graphics/` holds the "
        "main figures without captions, for `\\includegraphics`, and `captions.tex` holds every caption as TeX.", "",
        "| File | Figure | Source | Caption |", "|---|---|---|---|",
        *(f"| Figure{n}.pdf | {n} | revised layout (v2) | {'agreed' if n in agreed else 'to follow'} |"
          if n != "5" else
          "| Figure5.pdf | 5 | released Figure 6 panels (cell types), not yet restyled | released, renumbered |"
          for n in main_captions),
        *(f"| Figure{n}.pdf | {n} | released `{stem}` | released |" for n, stem in SUPPLEMENT.items()),
        f"| Table1.pdf, Table1.tex | Table 1 | released `{TABLE1}` | released |", "",
        "## Open items", "",
        "- Figure 1 is the authors' schematic file, unchanged: panels B-C still read \"expression\", and 1D writes "
        "d_LP*/d_CP* where the other figures write d_LP/d_CP.",
        "- Figure 5 is a placeholder until the cell-type figure is restyled; its caption is the released one.",
        "- Figure 7B is to follow.",
        *(f"- Released caption of {n} {why}." for n, why in OLD_REFERENCES.items()), "",
        "## Released figures not placed in the revised layout", "",
        *(f"- Released Figure {n}: {what}." for n, what in NOT_PLACED.items()), "",
    ]
    (out / "README.md").write_text("\n".join(readme))
    manifest = dict(version="manuscript-set-1", git=git_state(), v2_build=str(v2_build),
                    v2_build_record=sha256(v2_build / "build.json"),
                    production={k: dict(released_at=v["released_at"], presentation_id=v["presentation_id"])
                                for k, v in families.items()},
                    captions={("Figure " + n): ("agreed" if n in agreed else "released" if n == "5" else "to follow")
                              for n in main_captions},
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
