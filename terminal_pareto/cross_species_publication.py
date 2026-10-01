"""Assemble numbered Figures 8 and 9 from validated terminal comparison panels.

This presentation-only build writes to the comparison run's publication/
directory. It never modifies numerical caches or the accepted release.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto.cross_species_analysis import DEFAULT_RUN, file_hash, run_path, validate_run
from terminal_pareto.fig_terminal_cross_species import TRACKING_CAPTION_STATUS, TRACKING_CAVEAT
from terminal_pareto.publication_wrappers import _write

FIGURES = {
    "fig8_terminal_cross_species_comparison": (8, "terminal_cross_species_comparison"),
    "fig9_terminal_cross_species_overlays": (9, "terminal_cross_species_overlay"),
}


def captions(record):
    """Publication prose without species rankings or absolute-cost claims."""
    caveat = TRACKING_CAVEAT.replace("C. briggsae", r"\emph{C.\@ briggsae}")
    n, cb_reference = record["cohort_size"], record["cb_reference_size"]
    return {
        8: rf"""\textbf{{Terminal-only trade-offs across molecular representations and species.}}
The same {n} natural terminal parent--child edges and biological-parent
capacities are used throughout. Columns show \emph{{C.\@ elegans}} protein,
\emph{{C.\@ elegans}} RNA and \emph{{C.\@ briggsae}} AF16 RNA.
\textbf{{(A--C)}} Three-dimensional tracking distances.
\textbf{{(D--F)}} Two-dimensional XY distances, with z omitted and travel,
nulls, optima and endpoint coordinates recomputed.
Travel is the equal mean of per-embryo pairwise displacements divided by
fixed natural-lineage totals: three \emph{{C.\@ elegans}} embryos on the
275-edge reference and two AF16 embryos on their {cb_reference}-edge reference,
before restriction to the matched cohort. Displacement is not integrated
trajectory length. Cell-state cost is cosine dissimilarity across the frozen
20 z-scored protein reporters or the shared 20 RNA transcription factors
on stored values. Each sampled front connects 301 weighted assignment
solutions; connecting segments need not be attainable assignments.
Each panel uses its own endpoint cost spans, mapping the travel optimum to
\((0,1)\) and the cell-state optimum to \((1,0)\); natural lineage and nulls
share those anchors. All panels share axis limits and the retention scale.
Blue encodes biological-parent edge retention, black crosses mark natural
lineage, and outlined circles mark sampled maximum retention. Duplicate slots
of one parent do not count as different edges. First-, second- and third-cousin
shuffles are shown in the main axes; full-random shuffles use insets with the
same coordinate transform and separate limits. Plus signs mark the analytic
first-cousin mean or display-sample means for the other nulls.
These are descriptive comparisons, not calibrated significance tests or
absolute cost comparisons between datasets. {caveat} {TRACKING_CAPTION_STATUS}""",
        9: rf"""\textbf{{Endpoint-normalized terminal fronts and canonical metrics.}}
The matched {n}-edge cohort, pooled travel, molecular representations and
sampled assignment sweeps are the same as in Figure~8.
\textbf{{(A)}} Three-dimensional distances. \textbf{{(B)}} XY distances.
Blue solid lines show \emph{{C.\@ elegans}} protein, orange dashed lines
\emph{{C.\@ elegans}} RNA, and green dash-dotted lines \emph{{C.\@ briggsae}}
AF16 RNA. Colored crosses mark natural lineage; outlined circles mark
sampled maximum biological-parent retention. Hollow diamonds mark the closest
attained front assignment \(P^*\); dashed connectors run from natural lineage
to \(P^*\), not to maximum retention. Insets enlarge the region near
the natural lineages using the same coordinates as the main axes.
Both panels share main-axis and inset limits, colors and line styles.
Each dataset/geometry retains its own endpoint spans, with travel and
cell-state optima at \((0,1)\) and \((1,0)\), respectively.
\textbf{{(C)}} Named rows group the three configurations under 3D and XY,
following the separate-metric-column layout of Figure~5B. Gray shows canonical
position \(u\), purple natural-lineage distance \(d_{{LP}}\), and gold distance
\(d_{{NP}}\) from the analytic first-cousin null mean \(N\) to the same \(P^*\).
Here \(P^*\) minimizes Euclidean distance to natural lineage \(L\) among the
attained sampled assignments in endpoint coordinates;
\(d_{{LP}}=\lVert D(L)-D(P^*)\rVert_2\) and
\(d_{{NP}}=\lVert D(N)-D(P^*)\rVert_2\).
Position \(u\) is front arc length from the travel optimum to \(P^*\), divided
by total front arc length. Scales are shared across all six rows within each
metric; neither the closest point nor the distance is defined by maximum
retention. Overlays and metrics compare relative trade-off shapes and
natural-lineage proximity, not absolute travel or molecular costs.
Curves connect 301 attained weighted
solutions rather than enumerate the complete discrete Pareto set; connecting
segments need not be attainable assignments. Null clouds are omitted here
and shown in Figure~8. Protein/RNA differences are not a controlled modality
effect, and tracking replicates share molecular measurements.
{caveat} {TRACKING_CAPTION_STATUS}""",
    }


def write_wrappers(out, record):
    result = []
    for stem, (number, _) in FIGURES.items():
        body = (
            "\\begin{figure}[p]\n\\centering\n"
            f"\\includegraphics[width=\\textwidth]{{{stem}_panel.pdf}}\n"
            f"\\caption{{{captions(record)[number]}}}\n"
            f"\\label{{fig:terminal_cross_species_{number}}}\n\\end{{figure}}"
        )
        result.append(_write(Path(out), stem, body))
    return result


def archive_previous(run):
    source = run / "publication"
    if not source.exists() or not any(source.iterdir()):
        return None
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = run / "layout_history" / stamp / "publication"
    hashes = {str(p.relative_to(source)): file_hash(p) for p in source.rglob("*") if p.is_file()}
    shutil.copytree(source, target)
    for name, digest in hashes.items():
        if file_hash(target / name) != digest:
            raise ValueError(f"Layout archive verification failed: {name}")
    (target.parent / "manifest.json").write_text(json.dumps(dict(files=hashes), indent=2) + "\n")
    return target.parent


def compile_wrappers(wrappers):
    compiler = shutil.which("tectonic") or shutil.which("pdflatex")
    if compiler is None:
        raise RuntimeError("Activate dev: tectonic or pdflatex is required")
    for wrapper in wrappers:
        args = ([compiler, "--keep-logs", wrapper.name] if Path(compiler).name == "tectonic" else
                [compiler, "-interaction=nonstopmode", "-halt-on-error", wrapper.name])
        result = subprocess.run(args, cwd=wrapper.parent, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = wrapper.with_suffix(".build.log")
        log.write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(f"Compilation failed for {wrapper.name}; see {log}")
        info = subprocess.check_output(["pdfinfo", str(wrapper.with_suffix(".pdf"))], text=True)
        pages = int(next(line.split(":")[1] for line in info.splitlines() if line.startswith("Pages:")))
        if pages != 1:
            raise ValueError(f"Expected one-page Figure: {wrapper.stem} has {pages} pages")
        content = subprocess.check_output(["pdftotext", str(wrapper.with_suffix(".pdf")), "-"], text=True)
        number = FIGURES[wrapper.stem][0]
        if f"Figure {number}:" not in content:
            raise ValueError(f"Compiled label missing: Figure {number}")


def build(run_id=DEFAULT_RUN):
    run = run_path(run_id)
    report = validate_run(run)
    record = json.loads((run / "analysis/provenance.json").read_text())
    panels = run / "figures"
    manifest_path = panels / "figure_manifest.json"
    rendered = json.loads(manifest_path.read_text())
    if rendered["analysis_id"] != report["analysis_id"]:
        raise ValueError("Panel/numerical analysis identity mismatch")
    for name, digest in rendered["plotting_sources"].items():
        if file_hash(ROOT / name) != digest:
            raise ValueError(f"Changed plotting source: {name}; rerun --render-only first")
    for name, digest in rendered["files"].items():
        if file_hash(panels / name) != digest:
            raise ValueError(f"Changed rendered asset: {name}")
    archive = archive_previous(run)
    out = run / "publication"
    out.mkdir(exist_ok=True)
    for stem, (_, source) in FIGURES.items():
        for suffix in (".pdf", ".png"):
            shutil.copy2(panels / f"{source}{suffix}", out / f"{stem}_panel{suffix}")
    wrappers = write_wrappers(out, record)
    compile_wrappers(wrappers)
    manifest = dict(
        analysis_id=report["analysis_id"], validation=report,
        figure_numbers={stem: number for stem, (number, _) in FIGURES.items()},
        caption_tracking_wording_provisional=True, release_promoted=False,
        rendered_manifest_hash=file_hash(manifest_path),
        previous_layout=str(archive.relative_to(run)) if archive else None,
        sources={str(p.relative_to(ROOT)): file_hash(p) for p in (
            Path(__file__), ROOT / "terminal_pareto/publication_wrappers.py")},
        analysis_files={p.name: file_hash(p) for p in sorted((run / "analysis").iterdir()) if p.is_file()},
        files={p.name: file_hash(p) for p in sorted(out.iterdir()) if p.suffix in (".pdf", ".png", ".tex")},
    )
    (out / "publication_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    args = parser.parse_args()
    print(f"Numbered Figures 8 and 9 ready: {build(args.run_id)}")


if __name__ == "__main__":
    main()
