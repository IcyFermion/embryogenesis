"""Standalone captioned figure wrappers: preamble, numbering, compilation and checks.

Generated TeX is self-contained and does not depend on the repository at
compile time.
"""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess

PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{amsmath}
\usepackage{tikz}
\captionsetup{font=small,labelfont=bf}
\pagestyle{empty}
\begin{document}
"""
TEX_WARNINGS = re.compile(r"Overfull|Underfull|LaTeX Warning|Missing character|undefined", re.IGNORECASE)


def numbering(number: str) -> str:
    """Make a standalone wrapper print its manuscript label ("9", "S3", "1 amendment")."""
    if number.isdigit():
        return rf"\setcounter{{figure}}{{{int(number) - 1}}}" + "\n"
    return rf"\renewcommand{{\thefigure}}{{{number}}}" + "\n"


def write_figure_wrapper(out: Path, stem: str, *, number: str, graphic: str, caption: str, label: str,
                         width: str = r"\textwidth") -> Path:
    body = ("\\begin{figure}[p]\n\\centering\n"
            f"\\includegraphics[width={width}]{{{graphic}}}\n"
            f"\\caption{{{caption}}}\n"
            f"\\label{{{label}}}\n\\end{{figure}}")
    return write_wrapper(out, stem, number=number, body=body)


def compile_wrapper(wrapper: Path, *, label: str, pages: int = 1, known_warnings: tuple[str, ...] = (),
                    kind: str = "Figure") -> dict:
    """Compile with tectonic (or pdflatex); check page count, printed label and TeX warnings.

    ``known_warnings`` lists exact substrings of pre-existing, reviewed warnings
    that are recorded rather than fatal; any other warning fails the build.
    """
    compiler = shutil.which("tectonic") or shutil.which("pdflatex")
    if compiler is None:
        raise RuntimeError("Activate dev: tectonic or pdflatex is required")
    wrapper = Path(wrapper)
    args = ([compiler, "--keep-logs", wrapper.name] if Path(compiler).name == "tectonic" else
            [compiler, "-interaction=nonstopmode", "-halt-on-error", wrapper.name])
    result = subprocess.run(args, cwd=wrapper.parent, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    build_log = wrapper.with_suffix(".build.log")
    build_log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"Compilation failed for {wrapper.name}; see {build_log}")
    pdf = wrapper.with_suffix(".pdf")
    info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    found = int(next(line.split(":")[1] for line in info.splitlines() if line.startswith("Pages:")))
    if found != pages:
        raise ValueError(f"{wrapper.stem}: expected {pages} page(s), found {found}")
    text = subprocess.check_output(["pdftotext", str(pdf), "-"], text=True)
    if f"{kind} {label}:" not in text:
        raise ValueError(f"{wrapper.stem}: compiled label '{kind} {label}:' missing")
    tex_log = wrapper.with_suffix(".log")
    log_text = (tex_log.read_text(errors="replace") if tex_log.exists() else "") + result.stdout
    warnings = sorted({line.strip() for line in log_text.splitlines() if TEX_WARNINGS.search(line)})
    unexpected = [w for w in warnings if not any(known in w for known in known_warnings)]
    if unexpected:
        raise ValueError(f"{wrapper.stem}: TeX warnings: {unexpected}")
    return dict(pages=found, label=f"{kind} {label}", warnings=warnings)


def write_wrapper(out: Path, stem: str, *, number: str, body: str, preamble: str = PREAMBLE,
                  numbering_tex: str | None = None) -> Path:
    """Standalone wrapper around a complete ``figure`` environment."""
    path = Path(out) / f"{stem}.tex"
    counter = numbering(number) if numbering_tex is None else numbering_tex
    path.write_text(preamble + counter + body.strip() + "\n\\end{document}\n")
    return path
