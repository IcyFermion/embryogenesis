"""Assemble the manuscript set: ``python -m publication.v2.manuscript --v2-build DIR --output-dir DIR``.

The set is free-standing (author, 2026-11-13): every figure and Table 1 comes from one v2 build made
with ``--manuscript`` (figures alone, cropped), and nothing is taken from another figure set. Each
figure and the table has a standalone LaTeX wrapper, ``figN_<stem>.tex``, that holds its caption and
label; next to it sit the compiled captioned page (``figN_<stem>.pdf``) and the one graphic it
includes, ``figN_<stem>_panel.pdf``.

Numbering of 2026-11-10: main Figures 1-7 (1 schematics, 2 terminal front, 3 lineage similarity,
4 subtrees, 5 cell types, 6 full tree by layer, 7 C. briggsae front), S1 (terminal assignments in four
tracking geometries), S2 (cell-state features) and Table 1. Captions come from
``publication.v2.captions``: agreed, drafts in the v1 style awaiting review, or the released captions
of Figure 5 and S1 carried over unchanged. ``backlog/`` holds the v2 page 7B-C (C. briggsae against
C. elegans); every other figure kept as a backup is recorded, with its command, in
``publication/BACKLOG.md``.

Every wrapper is compiled here. Each must give one page, print its own number (read from the
``.aux``) and raise no TeX warning other than the tolerated float overflow (``publication.assembly``).
``review.tex`` joins the compiled pages into ``review.pdf``. ``README.md`` and ``manifest.json`` (file
hashes, the presentation sources and external files the v2 build recorded) complete the folder. The
v2 build must be current first.

The official set (author, 2026-11-12) is the one ``publication/output/manuscript/current`` points to::

    python -m publication.v2.manuscript --verify [SET]     # files against the manifest, then currency
    python -m publication.v2.manuscript --promote SET      # verify, then point ``current`` at SET

``--verify`` fails on a missing, extra or modified file. It reports a set as stale when presentation
sources or Figure 1's external files changed since assembly. ``--promote`` accepts only a set that is
current and was assembled from committed code; it appends to ``history.json``. To change a caption,
edit ``publication.v2.captions`` and assemble a new set (a hand-edited wrapper fails verification).
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from publication import assembly
from publication import notation as nt
from publication.provenance import ROOT, changed_sources, sha256
from publication.v2 import captions
from publication.v2.__main__ import EXTERNAL, NOTATION, presentation_sources

VERSION = "manuscript-release-3"  # 2: sources and external files recorded (2026-11-12); 3: free-standing (2026-11-13)
SETS = ROOT / "publication/output/manuscript"
CURRENT = "current"  # symlink in SETS to the official set
HISTORY = "history.json"  # promotions, oldest first
# Figures of the set, all drawn by the v2 manuscript-mode build: number -> stem.
MAIN = {"1": "fig1_overview", "2": "fig2_terminal_pareto", "3": "fig3_lineage_similarity", "4": "fig4_subtrees",
        "5": "fig5_cell_types", "6": "fig6_full_tree_layers", "7": "fig7_briggsae",
        "S1": "figS1_tracking_robustness", "S2": "figS2_cell_state_features"}
# v2 pages kept in reserve (author, 2026-11-10): number printed -> (v2 build stem, stem in backlog/).
BACKLOG_V2 = {"7B--C": ("backlog_fig7BC_briggsae_counterparts", "fig7BC_briggsae_counterparts")}
TABLE1 = "table1_subtree_statistics"
TABLE1_LABEL = "tab:ce_subtree_statistics"
# The main figure fills the text width and at most 70% of the text height, so its caption fits on the page.
GRAPHIC_SIZE = r"width=\textwidth,height=.7\textheight,keepaspectratio"
REVIEW = r"""% Review copy: every page of the set, in order. Recompile here: tectonic review.tex
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


def main_figure(v2_build: Path, out: Path, number: str, stem: str, caption: str, status: str,
                source: str | None = None) -> dict:
    """A v2 figure in a captioned wrapper; ``source`` is its stem in the v2 build when that differs."""
    source = source or stem
    graphic = f"{stem}_panel.pdf"
    shutil.copy2(v2_build / f"{source}.pdf", out / graphic)
    preview(out / graphic)
    label = "fig:" + stem.split("_", 1)[1]
    note = {"agreed": "agreed with the author", "draft": "draft in the v1 style, awaiting the author's review",
            "released": "the released caption, unchanged", "backlog": "kept in reserve, not in review.pdf"}[status]
    body = (f"% Caption: {note}; source publication/v2/captions.py.\n"
            f"% Graphic: {source}.pdf of the v2 manuscript-mode build {v2_build.name}.\n"
            "\\begin{figure}[p]\n\\centering\n"
            f"\\includegraphics[{GRAPHIC_SIZE}]{{{graphic}}}\n"
            f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{figure}}")
    tex = assembly.write_wrapper(out, stem, number=number, body=body)
    return dict(stem=stem, source=f"v2 build ({source}.pdf)", caption=status, **compile_page(tex, number=number,
                                                                                              label=label))


def git_state() -> dict:
    run = lambda *args: subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True).stdout.strip()
    return dict(commit=run("rev-parse", "HEAD"), branch=run("branch", "--show-current"),
                dirty=bool(run("status", "--porcelain", "--untracked-files=no")))


def notes(record: dict) -> list[str]:
    """Numbers the main text should carry, recorded by the v2 build (``main_text_notes``)."""
    d = record["main_text_notes"]["tracking_front_distance"]
    embryos = [d[f"embryo{i}"] for i in (1, 2, 3)]
    return [
        "- **Supplementary Figure S1, d_LP from pooled tracking.** The natural lineage's distance to "
        "its own front, each in that geometry's endpoint coordinates: embryo 1 "
        + ", embryo 2 ".join(f"{v:.3f}" for v in embryos[:2]) + f", embryo 3 {embryos[2]:.3f} "
        f"(mean {np.mean(embryos):.3f}); pooled travel {d['pooled']:.3f} (the global d_LP of Figure 2). Pooling "
        "gives a smaller distance than two of the three single embryos and than their mean, but not the smallest "
        "(embryo 3); the text should say the pooled estimate averages out single-embryo tracking variation, not "
        "that it is closest to the front.",
    ]


def assemble(v2_build: Path, out: Path) -> dict:
    v2_build, out = Path(v2_build), Path(out)
    record = json.loads((v2_build / "build.json").read_text())
    if not record.get("manuscript"):
        raise ValueError(f"{v2_build} is a draft build; rebuild it with `python -m publication.v2 --manuscript`")
    built = {Path(p).stem for p in record["figures"].values()} | {Path(record.get("table1", "")).stem}
    if missing := sorted({*MAIN.values(), TABLE1} - built):
        raise ValueError(f"{v2_build} does not hold {missing}; rebuild it")
    if sorted(stem for stem, _ in BACKLOG_V2.values()) != sorted(Path(p).stem for p in record.get("backlog", {}).values()):
        raise ValueError(f"{v2_build} does not hold the backlog pages {sorted(BACKLOG_V2)}; rebuild it")
    if "sources" not in record:
        raise ValueError(f"{v2_build} records no presentation sources; rebuild it")
    if changed := changed_sources(record["sources"], presentation_sources()):
        raise ValueError(f"{v2_build} is stale, sources changed since it was built: {changed}; rebuild it")
    if changed := external_changes(record["external"]):
        raise ValueError(f"{v2_build} is stale, external files changed since it was built: {changed}; rebuild it")
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory must be new or empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    with nt.override(**NOTATION):
        status = {**{n: "agreed" for n in captions.AGREED}, **{n: "draft" for n in captions.DRAFT},
                  **{n: "released" for n in captions.RELEASED}}
        text = {n: make() for n, make in {**captions.AGREED, **captions.DRAFT, **captions.RELEASED}.items()}
        reserve = {n: make() for n, make in captions.BACKLOG.items()}
    if missing := [n for n in MAIN if n not in text]:
        raise ValueError(f"No caption for figures {missing}")
    figures = {n: main_figure(v2_build, out, n, MAIN[n], text[n], status[n])
               for n in sorted(MAIN, key=lambda k: (k.startswith("S"), int(k.lstrip("S"))))}
    shutil.copy2(v2_build / f"{TABLE1}.tex", out / f"{TABLE1}.tex")
    table = dict(stem=TABLE1, source=f"v2 build ({TABLE1}.tex)", caption="released, citing Figure 4",
                 **compile_page(out / f"{TABLE1}.tex", number="1", label=TABLE1_LABEL))
    backlog = out / "backlog"
    backlog.mkdir()
    reserved = {number: main_figure(v2_build, backlog, number, stem, reserve[number.replace("--", "-")], "backlog",
                                    source=source)
                for number, (source, stem) in BACKLOG_V2.items()}

    pages = [f"{figures[n]['stem']}.pdf" for n in figures] + [f"{TABLE1}.pdf"]
    (out / "review.tex").write_text(REVIEW + "".join(f"\\includepdf[pages=-]{{{p}}}\n" for p in pages)
                                    + "\\end{document}\n")
    result = subprocess.run(["tectonic", "review.tex"], cwd=out, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode or page_count(out / "review.pdf") != len(pages):
        raise RuntimeError(f"review.tex failed:\n{result.stdout[-3000:]}")

    overflow = [n for n, f in figures.items() if f["warnings"]]
    readme = [
        "# Manuscript set (revised layout)", "",
        f"Built by `python -m publication.v2.manuscript` from the v2 manuscript-mode build `{v2_build.name}`, which "
        "draws every figure and Table 1; nothing is taken from another figure set. Every figure and the table has a "
        "standalone LaTeX wrapper (`.tex`) with its caption and label, the compiled captioned page (`.pdf`) and the "
        "one graphic it includes, `<stem>_panel.pdf`. `review.pdf` joins every page.", "",
        "The official set is the one `publication/output/manuscript/current` points to. `python -m "
        "publication.v2.manuscript --verify <set>` checks these files against `manifest.json` and reports whether "
        "the set is still current. Captions come from `publication/v2/captions.py`: change them there and assemble "
        "a new set; a wrapper edited here fails verification.", "",
        "| Wrapper | Figure | Label | Source | Caption |", "|---|---|---|---|---|",
        *(f"| `{f['stem']}.tex` | {n} | `{f['label']}` | {f['source']} | {f['caption']} |" for n, f in figures.items()),
        f"| `{TABLE1}.tex` | Table 1 | `{TABLE1_LABEL}` | {table['source']} | {table['caption']} |", "",
        "Every wrapper compiled here to one page that prints its own number"
        + (f"; figures {', '.join(overflow)} report a tolerated float overflow." if overflow else
           ", with no TeX warning."),
        "PNG previews of the figures are written only where poppler (`pdftoppm`) is installed.", "",
        "## Notes for the main text", "",
        *notes(record), "",
        "## Open items", "",
        "- Figure 1 is the authors' schematic file, unchanged: panels B-C still read \"expression\", and 1D writes "
        "d_LP*/d_CP* where the other figures write d_LP/d_CP.",
        "- Figures 5 and S1 were redrawn in the v2 style on 2026-11-13 and keep their released captions (Figure 5's "
        "says natural-edge retention).",
        "- Captions of Figures " + ", ".join(n for n, f in figures.items() if f["caption"] == "draft")
        + " are drafts in the v1 style, for the author's review (`publication/v2/captions.py`).",
        "- Table 1 keeps the released caption, now citing Figure 4; its column header and caption still say maximum "
        "edge retention, which the figures call lineage similarity.", "",
        "## Backlog (kept in reserve, not in review.pdf)", "",
        *(f"- `backlog/{f['stem']}.tex`: Figure {n.replace('--', '-')}, C. briggsae counterparts of Figures 3A and 6B "
          "(the former 8B-C), with its caption. Moved out of Figure 7 on 2026-11-10: the gaps in data resolution "
          "between the two species make these comparisons hard to justify; kept for reviewers."
          for n, f in reserved.items()),
        "- Every other figure kept as a backup (released figures not placed here, sketches, reserve notes) is "
        "recorded with its regeneration command in `publication/BACKLOG.md`.", "",
    ]
    (out / "README.md").write_text("\n".join(readme))
    manifest = dict(version=VERSION, git=git_state(), v2_build=str(v2_build),
                    v2_build_record=sha256(v2_build / "build.json"), sources=record["sources"],
                    external=record["external"], figures=figures, table1=table, backlog=dict(v2=reserved),
                    files={p.relative_to(out).as_posix(): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def external_changes(recorded: dict) -> list[str]:
    """External files (Figure 1's sources in ``paper/FIGS``) missing or changed since ``recorded``."""
    return sorted(n for n, digest in recorded.items()
                  if n not in EXTERNAL or not EXTERNAL[n].is_file() or sha256(EXTERNAL[n]) != digest)


def check_files(set_dir: Path, files: dict) -> None:
    """The set holds exactly the recorded files (``manifest.json`` aside), each with its recorded hash."""
    present = {p.relative_to(set_dir).as_posix() for p in set_dir.rglob("*") if p.is_file()} - {"manifest.json"}
    if present != set(files):
        raise ValueError(f"{set_dir}: missing {sorted(set(files) - present)}, unexpected {sorted(present - set(files))}")
    if modified := sorted(n for n, digest in files.items() if sha256(set_dir / n) != digest):
        raise ValueError(f"{set_dir}: modified since assembly: {modified}")


def verify(set_dir: Path = SETS / CURRENT) -> dict:
    """Integrity (raises) and currency (reported in ``stale``): presentation sources and Figure 1's external files,
    each against what the set recorded."""
    set_dir = Path(set_dir)
    manifest = json.loads((set_dir / "manifest.json").read_text())
    check_files(set_dir, manifest["files"])
    stale = dict(
        sources=(changed_sources(manifest["sources"], presentation_sources()) if "sources" in manifest
                 else [f"not recorded ({manifest['version']})"]),
        external=external_changes(manifest["external"]) if "external" in manifest else ["not recorded"],
        assembly=[] if "production" not in manifest else ["reuses a released figure set (before 2026-11-13)"])
    return dict(set=str(set_dir.resolve()), version=manifest["version"], files=len(manifest["files"]),
                git=manifest["git"], stale={k: v for k, v in stale.items() if v})


def point_current(sets: Path, set_dir: Path, record: dict) -> dict:
    """Point ``sets/current`` at ``set_dir`` (a folder directly in ``sets``) and append the promotion to the history."""
    sets, set_dir = Path(sets), Path(set_dir)
    if set_dir.resolve().parent != sets.resolve() or not set_dir.is_dir():
        raise ValueError(f"{set_dir} must be a set folder directly in {sets}")
    link, history = sets / CURRENT, sets / HISTORY
    previous = os.readlink(link) if link.is_symlink() else None
    if link.exists() and not link.is_symlink():
        raise ValueError(f"{link} exists and is not a symlink")
    entry = dict(set=set_dir.name, previous=previous, promoted_at=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
                 **record)
    entries = json.loads(history.read_text()) if history.exists() else []
    staged = sets / f".{CURRENT}.tmp"
    if staged.is_symlink():
        staged.unlink()
    staged.symlink_to(set_dir.name)
    os.replace(staged, link)
    history.write_text(json.dumps([*entries, entry], indent=2) + "\n")
    return entry


def promote(set_dir: Path, sets: Path = SETS) -> dict:
    """Make ``set_dir`` the official set: it must verify, be current and come from committed code."""
    result = verify(set_dir)
    if result["stale"]:
        raise ValueError(f"{set_dir} is not current: {result['stale']}; assemble a new set")
    if result["git"]["dirty"]:
        raise ValueError(f"{set_dir} was assembled from uncommitted changes; commit, rebuild and assemble again")
    return point_current(sets, set_dir, dict(commit=result["git"]["commit"], branch=result["git"]["branch"],
                                             version=result["version"], files=result["files"]))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--v2-build", type=Path, help="assemble a set from this v2 manuscript-mode build")
    parser.add_argument("--output-dir", type=Path, help="new or empty folder for the assembled set")
    parser.add_argument("--verify", type=Path, nargs="?", const=SETS / CURRENT, metavar="SET",
                        help="check a set (default: the official set)")
    parser.add_argument("--promote", type=Path, metavar="SET", help="verify SET and make it the official set")
    args = parser.parse_args(argv)
    if args.verify or args.promote:
        if args.v2_build or args.output_dir or (args.verify and args.promote):
            parser.error("--verify and --promote take no other option")
        print(json.dumps(verify(args.verify) if args.verify else promote(args.promote), indent=2))
        return
    if not (args.v2_build and args.output_dir):
        parser.error("assembling needs --v2-build and --output-dir")
    manifest = assemble(args.v2_build, args.output_dir)
    summary = {n: dict(stem=f["stem"], caption=f["caption"], warnings=len(f["warnings"]))
               for n, f in manifest["figures"].items()}
    print(json.dumps(dict(figures=summary, files=len(manifest["files"])), indent=2))


if __name__ == "__main__":
    main()
