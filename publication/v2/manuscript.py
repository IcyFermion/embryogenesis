"""Assemble the manuscript release set: ``python -m publication.v2.manuscript --v2-build DIR --output-dir DIR``.

The release set is laid out like ``publication/output/production`` (author request, 2026-11-09).
Each figure and table has a standalone LaTeX wrapper, ``figN_<stem>.tex``, that holds its caption
and label. Next to it sit the compiled captioned page (``figN_<stem>.pdf``) and the graphics it
includes.

Numbering of 2026-11-08:

- Figures 1-5, 7 and 8 come from a v2 build made with ``--manuscript`` (figures alone, cropped).
  Each wrapper includes that figure as ``figN_<stem>_panel.pdf``, as released Figures 4 and 8-11
  include theirs. Captions are agreed or, until the author approves them, drafts in the v1 style
  (``publication.v2.captions``).
- Figure 6 is the released cell-type figure (wrapper, panels and caption) until it is restyled.
- The supplement is released S3, renumbered S1 (wrapper and graphic renamed ``figS1_*``).
- Table 1 is the released table.
- Released S1 is retired (its r is in Figure 5). Released S2, S4 and S5 are copied unchanged, with
  their wrappers and graphics, into ``backlog/``.

Every wrapper is compiled here. Each must give one page, print its own number (read from the
``.aux``) and raise no TeX warning other than the tolerated float overflow (``publication.assembly``).
``review.tex`` joins the compiled pages into ``review.pdf``. To change a caption, edit the figure's
``.tex``, recompile it, then recompile ``review.tex``. ``README.md`` and ``manifest.json`` (with
file hashes) complete the folder. Production must verify cleanly first.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from publication import assembly
from publication import notation as nt
from publication.provenance import ROOT, sha256
from publication.release import PRODUCTION, verify_production
from publication.v2 import captions
from publication.v2.__main__ import NOTATION

VERSION = "manuscript-release-1"
# Main figures of the revised layout: number -> stem of the v2 manuscript-mode build.
MAIN = {"1": "fig1_overview", "2": "fig2_cell_state_features", "3": "fig3_terminal_pareto",
        "4": "fig4_lineage_similarity", "5": "fig5_subtrees", "7": "fig7_full_tree_layers", "8": "fig8_briggsae"}
# Released figures placed in the set: number -> (released stem, stem in the set). Figure 6 holds the
# released cell-type figure until it is restyled; S1 is released S3 (author decisions of 2026-11-08).
ADOPTED = {"6": ("fig6_ce_cell_types", "fig6_ce_cell_types"),
           "S1": ("figS3_ce_tracking_robustness", "figS1_ce_tracking_robustness")}
SUPPLEMENT = ("S1",)
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
# The main figure fills the text width and at most 70% of the text height, so its caption fits on the page.
GRAPHIC_SIZE = r"width=\textwidth,height=.7\textheight,keepaspectratio"
INCLUDED = re.compile(r"\\(?:includegraphics(?:\[[^\]]*\])?|input)\{([^}]+)\}")
NUMBERING = re.compile(r"\\setcounter\{figure\}\{\d+\}\n|\\renewcommand\{\\thefigure\}\{[^}]*\}\n")
LABEL = re.compile(r"\\label\{([^}]+)\}")
REVIEW = r"""% Review copy: every page of the release set, in order. Recompile here: tectonic review.tex
% Captions live in each figure's .tex; recompile that file first, then this one.
\documentclass{article}
\usepackage{pdfpages}
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


def page_count(pdf: Path) -> int:
    if shutil.which("qpdf"):
        return int(subprocess.check_output(["qpdf", "--show-npages", str(pdf)], text=True))
    if shutil.which("pdfinfo"):
        info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
        return int(next(line.split(":")[1] for line in info.splitlines() if line.startswith("Pages:")))
    raise RuntimeError("qpdf or pdfinfo is required to count pages")


def preview(pdf: Path) -> Path | None:
    """A 300 dpi PNG of a one-page graphic, when poppler is installed (the released sets carry them)."""
    if not shutil.which("pdftoppm"):
        return None
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile", str(pdf), str(pdf.with_suffix(""))], check=True)
    return pdf.with_suffix(".png")


def compile_page(tex: Path, *, number: str, label: str) -> dict:
    """Compile a wrapper with tectonic and check it as a release does: one page, its own number, no TeX warning
    other than the tolerated float overflow. The printed number is read from the label in the ``.aux``."""
    result = subprocess.run(["tectonic", "--keep-logs", "--keep-intermediates", tex.name], cwd=tex.parent, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{tex.name} failed to compile:\n{result.stdout[-3000:]}")
    aux = tex.with_suffix(".aux").read_text(errors="replace")
    printed = re.search(r"\\newlabel\{" + re.escape(label) + r"\}\{\{([^}]*)\}", aux)
    if not printed or printed.group(1) != number:
        raise ValueError(f"{tex.name}: prints {printed.group(1) if printed else 'no label'}, expected {number}")
    pages = page_count(tex.with_suffix(".pdf"))
    if pages != 1:
        raise ValueError(f"{tex.name}: {pages} pages, expected 1")
    log = tex.with_suffix(".log").read_text(errors="replace") + result.stdout
    warnings = sorted({line.strip() for line in log.splitlines() if assembly.TEX_WARNINGS.search(line)})
    unexpected = [w for w in warnings if not any(ok in w for ok in assembly.TOLERATED_WARNINGS)]
    if unexpected:
        raise ValueError(f"{tex.name}: TeX warnings: {unexpected}")
    for side in tex.parent.glob(f"{tex.stem}.*"):  # intermediates; the set keeps .tex and .pdf
        if side.suffix not in (".tex", ".pdf"):
            side.unlink()
    return dict(number=number, label=label, pages=pages, warnings=warnings)


def copy_released(name: str, target: Path) -> list[Path]:
    """A released graphic and, if production has one, its PNG preview."""
    copied = [shutil.copy2(PRODUCTION / name, target)]
    png = (PRODUCTION / name).with_suffix(".png")
    if png.suffix != Path(name).suffix and png.exists():
        copied.append(shutil.copy2(png, Path(target).with_suffix(".png")))
    return [Path(p) for p in copied]


def main_figure(v2_build: Path, out: Path, number: str, stem: str, caption: str, status: str) -> dict:
    graphic = f"{stem}_panel.pdf"
    shutil.copy2(v2_build / f"{stem}.pdf", out / graphic)
    preview(out / graphic)
    label = "fig:" + stem.split("_", 1)[1]
    note = ("agreed with the author" if status == "agreed" else
            "draft in the v1 style, awaiting the author's review")
    body = (f"% Caption: {note}; source publication/v2/captions.py.\n"
            f"% Graphic: {stem}.pdf of the v2 manuscript-mode build {v2_build.name}.\n"
            "\\begin{figure}[p]\n\\centering\n"
            f"\\includegraphics[{GRAPHIC_SIZE}]{{{graphic}}}\n"
            f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{figure}}")
    tex = assembly.write_wrapper(out, stem, number=number, body=body)
    return dict(stem=stem, source=f"v2 build ({stem}.pdf)", caption=status, **compile_page(tex, number=number,
                                                                                            label=label))


def renumber(tex: str, released: str, stem: str, number: str) -> tuple[str, dict]:
    """A released wrapper under a new number, its included files renamed to the new prefix
    (``figS3_*``, ``figS3A_*`` -> ``figS1_*``, ``figS1A_*``). Returns the TeX and the renames."""
    old, new = released.split("_", 1)[0], stem.split("_", 1)[0]
    prefix = re.compile(re.escape(old) + r"(?=[A-Z]?_)")
    names = {name: prefix.sub(new, name, count=1) for name in INCLUDED.findall(tex)}
    tex = INCLUDED.sub(lambda m: m.group(0).replace(m.group(1), names[m.group(1)]), tex)
    note = f"% Released caption of production {released}, unchanged; numbered {number} here.\n"
    tex, count = NUMBERING.subn(lambda _: assembly.numbering(number) + note, tex)
    if count != 1:
        raise ValueError(f"{released}: expected one numbering line, found {count}")
    return tex, names


def adopt(out: Path, number: str, released: str, stem: str) -> dict:
    """A released figure under its new number: wrapper, caption and graphics, renamed with the new prefix."""
    tex, names = renumber((PRODUCTION / f"{released}.tex").read_text(), released, stem, number)
    for name, new in names.items():
        copy_released(name, out / new)
    label = LABEL.search(tex).group(1)
    (out / f"{stem}.tex").write_text(tex)
    return dict(stem=stem, source=f"released {released}", caption="released",
                **compile_page(out / f"{stem}.tex", number=number, label=label))


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
    out.mkdir(parents=True, exist_ok=True)

    with nt.override(**NOTATION):
        agreed = {n: make() for n, make in captions.AGREED.items()}
        drafts = {n: make() for n, make in captions.DRAFT.items()}
    missing = [n for n in MAIN if n not in agreed and n not in drafts]
    if missing:
        raise ValueError(f"No caption for main figures {missing}")
    figures = {}
    for n in sorted([*MAIN, *ADOPTED], key=lambda k: (k.startswith("S"), int(k.lstrip("S")))):
        if n in MAIN:
            figures[n] = main_figure(v2_build, out, n, MAIN[n], agreed.get(n) or drafts[n],
                                     "agreed" if n in agreed else "draft")
        else:
            figures[n] = adopt(out, n, *ADOPTED[n])
    for suffix in (".tex", ".pdf"):  # Table 1: the released page, unchanged
        shutil.copy2(PRODUCTION / f"{TABLE1}{suffix}", out / f"{TABLE1}{suffix}")
    backlog = out / "backlog"
    backlog.mkdir()
    for stem in BACKLOG.values():  # kept in reserve, under their released numbers
        tex = (PRODUCTION / f"{stem}.tex").read_text()
        for suffix in (".tex", ".pdf"):
            shutil.copy2(PRODUCTION / f"{stem}{suffix}", backlog / f"{stem}{suffix}")
        for name in INCLUDED.findall(tex):
            copy_released(name, backlog / name)

    pages = [f"{figures[n]['stem']}.pdf" for n in figures] + [f"{TABLE1}.pdf"]
    (out / "review.tex").write_text(REVIEW + "".join(f"\\includepdf[pages=-]{{{p}}}\n" for p in pages)
                                    + "\\end{document}\n")
    result = subprocess.run(["tectonic", "review.tex"], cwd=out, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode or page_count(out / "review.pdf") != len(pages):
        raise RuntimeError(f"review.tex failed:\n{result.stdout[-3000:]}")

    families = json.loads((PRODUCTION / "release_manifest.json").read_text())["families"]
    overflow = [n for n, f in figures.items() if f["warnings"]]
    readme = [
        "# Manuscript release set (revised layout)", "",
        f"Built by `python -m publication.v2.manuscript` from `{v2_build.name}` (v2 manuscript-mode build) and "
        f"production released on {max(f['released_at'] for f in families.values())[:8]}. Laid out like "
        "`publication/output/production`: every figure and table has a standalone LaTeX wrapper (`.tex`) with its "
        "caption and label, the compiled captioned page (`.pdf`) and the graphics it includes. The main figures of "
        "the revised layout include their graphic as `<stem>_panel.pdf`. To change a caption, edit the `.tex` and "
        "run `tectonic <stem>.tex` here. `review.pdf` joins every page (`tectonic review.tex` rebuilds it).", "",
        "| Wrapper | Figure | Label | Source | Caption |", "|---|---|---|---|---|",
        *(f"| `{f['stem']}.tex` | {n} | `{f['label']}` | {f['source']} | {f['caption']} |" for n, f in figures.items()),
        f"| `{TABLE1}.tex` | Table 1 | `tab:ce_subtree_statistics` | released `{TABLE1}` | released |", "",
        "Every wrapper compiled here to one page that prints its own number"
        + (f"; figures {', '.join(overflow)} report a tolerated float overflow." if overflow else
           ", with no TeX warning."),
        "PNG previews of the main figures are written only where poppler (`pdftoppm`) is installed.", "",
        "## Notes for the main text", "",
        *notes(record), "",
        "## Open items", "",
        "- Figure 1 is the authors' schematic file, unchanged: panels B-C still read \"expression\", and 1D writes "
        "d_LP*/d_CP* where the other figures write d_LP/d_CP.",
        "- Figure 6 is a placeholder until the cell-type figure is restyled; its caption is the released one.",
        "- Captions of Figures " + ", ".join(n for n, f in figures.items() if f["caption"] == "draft")
        + " are drafts in the v1 style, for the author's review (`publication/v2/captions.py`).",
        "- Table 1 is the released table: its column header and caption still say maximum edge retention, which "
        "the figures now call lineage similarity.",
        *(f"- Released caption of {n} {why}." for n, why in OLD_REFERENCES.items()), "",
        "## Backlog (kept in reserve, not in review.pdf)", "",
        *(f"- `backlog/{stem}.tex`: released Figure {n}, under its released number, with its graphics."
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
    manifest = dict(version=VERSION, git=git_state(), v2_build=str(v2_build),
                    v2_build_record=sha256(v2_build / "build.json"),
                    production={k: dict(released_at=v["released_at"], presentation_id=v["presentation_id"])
                                for k, v in families.items()},
                    figures=figures, table1=dict(stem=TABLE1, source=f"released {TABLE1}"),
                    backlog=BACKLOG, retired={n: stem for n, (stem, _) in RETIRED.items()},
                    files={str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--v2-build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = assemble(args.v2_build, args.output_dir)
    summary = {n: dict(stem=f["stem"], caption=f["caption"], warnings=len(f["warnings"]))
               for n, f in manifest["figures"].items()}
    print(json.dumps(dict(figures=summary, files=len(manifest["files"])), indent=2))


if __name__ == "__main__":
    main()
