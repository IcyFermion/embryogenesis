"""Numbered Figures 10/11 and explicit, additive production promotion."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from full_tree_pareto.cross_species_analysis import DEFAULT_RUN, ROOT, digest, run_path, validate_run, check_hashes
from full_tree_pareto.fig_cross_species import COMPARISON, OVERLAY, archive_previous
from terminal_pareto.fig_terminal_cross_species import TRACKING_CAVEAT, TRACKING_CAPTION_STATUS
from terminal_pareto.publication_wrappers import _write

FIGURES = {
    "fig10_full_tree_cross_species_comparison": (10, COMPARISON),
    "fig11_full_tree_cross_species_overlays": (11, OVERLAY),
}
PRODUCTION = ROOT / "full_tree_pareto/output/publication"
RELEASE_MANIFEST = "cross_species_release_manifest.json"
ASSEMBLY_COPY = "cross_species_publication_manifest.json"
LEGACY_FIG8 = (
    "fig8_ce_full_tree_collective.pdf", "fig8_ce_full_tree_collective.tex",
    "fig8_ce_full_tree_collective_panel.pdf", "fig8_ce_full_tree_collective_panel.png",
)


def publication_files():
    return {f"{stem}{suffix}" for stem in FIGURES for suffix in (".pdf", ".tex", "_panel.pdf", "_panel.png")}


def inventory(directory):
    directory = Path(directory)
    return {str(path.relative_to(directory)): digest(path) for path in directory.rglob("*") if path.is_file()}


def captions(record):
    n, nodes, terminals, roots = (record[key] for key in ("edges", "cohort_size", "terminal_count", "roots"))
    rounds = ", ".join(map(str, record["round_edges"]))
    weights, draws = record["settings"]["intervals"]+1, record["settings"]["draws"]
    caveat = TRACKING_CAVEAT.replace("C. briggsae", r"\emph{C.\@ briggsae}")
    return {
        COMPARISON: rf"""\textbf{{Layerwise trade-offs on a matched terminal-anchored partial forest.}}
Columns show \emph{{C.\@ elegans}} protein, \emph{{C.\@ elegans}} RNA and
\emph{{C.\@ briggsae}} AF16 RNA; \textbf{{(A--C)}} use 3D tracking and
\textbf{{(D--F)}} use XY, recomputing distances, nulls and assignments after omitting z.
Starting at the same {terminals} matched biological terminal edges as Figures~8--9,
each branch ascends through observed canonical ancestors and stops at its first
measurement gap. The resulting {nodes} measured cells and {n} edges form a
partial forest with {roots} fixed boundary roots, not a complete embryonic tree.
No missing state is imputed and no ancestor is skipped. Six asynchronous
contraction rounds ({rounds} edges) preserve observed parent-slot multiplicities,
including one-child boundaries. Exact linear assignments at {weights} shared
weights use global analytic first-cousin-null SD scaling; aggregate costs sum
all rounds. This optimizes a product of round-wise assignment spaces, not
unrestricted tree reconstruction. Travel averages pairwise displacements from
three CE embryos (cutoffs 255/247/225) or two AF16 embryos (148/156), each divided
by its fixed natural total on these {n} edges. Coordinates are not averaged.
Cell-state distance is Euclidean across 20 z-scored protein reporters or the
shared 20 RNA TFs on stored values, unlike terminal cosine distance.
Each panel maps its own travel optimum to \((0,1)\) and state optimum to \((1,0)\);
lineage and nulls share those anchors, with common main-axis limits.
Blue indicates biological-parent retention; black crosses mark natural lineage
and outlined circles sampled maximum retention. Cousin shuffles permute only
the {terminals} terminal identities within canonical ancestor groups two, three
or four transitions back, keeping internal states fixed and scoring all {n}
edges. Random rebuild fixes roots and observed capacities and appears in insets
with the same transform but separate limits. Dots show 1,000 of {draws:,} draws;
plus signs show the analytic first-cousin mean or full-draw means for other
references. These descriptive comparisons are not calibrated significance
tests or absolute-cost species rankings. {caveat} {TRACKING_CAPTION_STATUS}""",
        OVERLAY: rf"""\textbf{{Partial-forest front overlays and canonical metrics.}}
The same {n}-edge, {nodes}-cell terminal-anchored cohort, pooled travel and
layerwise assignment spaces are used as in Figure~10.
\textbf{{(A)}} 3D tracking. \textbf{{(B)}} XY tracking. Blue solid lines show
\emph{{C.\@ elegans}} protein, orange dashed lines \emph{{C.\@ elegans}} RNA,
and green dash-dotted lines \emph{{C.\@ briggsae}} AF16 RNA. Crosses mark natural
lineage and outlined circles sampled maximum biological-parent retention.
Hollow diamonds mark the closest attained assignment \(P^*\), and dashed
connectors join natural lineage to \(P^*\), not maximum retention. Insets
enlarge the natural-lineage region without changing coordinates. Both panels
share main-axis and zoom limits; each configuration retains its own aggregate
layerwise endpoint spans. \textbf{{(C)}} Named rows group configurations under
3D and XY, using separate metric columns as in Figures~5B and 9C.
Gray shows position \(u\), purple distance \(d_{{LP}}\) from natural lineage
\(L\) to \(P^*\), and gold distance \(d_{{NP}}\) from the analytic partial-forest
first-cousin-null mean \(N\) to the same \(P^*\).
\(P^*\) minimizes Euclidean distance to \(L\) among sampled attained
assignments in endpoint coordinates; \(u\) is front arc length from the travel
optimum to \(P^*\), divided by total arc length. Metric scales are shared
across all six rows within each column. Curves connect {weights} attained
weighted solutions, not the complete discrete Pareto set; connecting segments
need not be attainable. Separate {record['settings']['dense_intervals']+1:,}-weight
checks are retained in the numerical cache. These panels compare relative
trade-off shapes and lineage proximity, not absolute travel or molecular costs.
Coverage determines the partial-tree boundary; developmental alignment and
RNA measurement provenance remain unresolved. Protein/RNA differences are not
a controlled modality effect, and tracking replicates share molecular matrices.
{caveat} {TRACKING_CAPTION_STATUS}""",
    }


def write_wrappers(out, record):
    result = []
    for stem, (number, source) in FIGURES.items():
        body = ("\\begin{figure}[p]\n\\centering\n"
                + rf"\includegraphics[width=\textwidth]{{{stem}_panel.pdf}}" + "\n"
                + rf"\caption{{{captions(record)[source]}}}" + "\n"
                + rf"\label{{fig:full_tree_cross_species_{number}}}" + "\n\\end{figure}\n")
        result.append(_write(Path(out), stem, body))
    return result


def compile_wrappers(wrappers):
    compiler = shutil.which("tectonic") or shutil.which("pdflatex")
    if compiler is None:
        raise RuntimeError("Activate dev: tectonic or pdflatex is required")
    for wrapper in wrappers:
        args = ([compiler, "--keep-logs", wrapper.name] if Path(compiler).name == "tectonic" else
                [compiler, "-interaction=nonstopmode", "-halt-on-error", wrapper.name])
        result = subprocess.run(args, cwd=wrapper.parent, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        wrapper.with_suffix(".build.log").write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(f"Failed compilation: {wrapper}")
        info = subprocess.check_output(["pdfinfo", str(wrapper.with_suffix(".pdf"))], text=True)
        pages = int(next(line.split(":")[1] for line in info.splitlines() if line.startswith("Pages:")))
        if pages != 1:
            raise ValueError(f"Figure must fit one page: {wrapper.stem} has {pages}")
        text = subprocess.check_output(["pdftotext", str(wrapper.with_suffix(".pdf")), "-"], text=True)
        if f"Figure {FIGURES[wrapper.stem][0]}:" not in text:
            raise ValueError("Missing numbered figure label")


def verify(run, *, publication=None):
    run = Path(run)
    report = validate_run(run)
    figures = json.loads((run / "figures/figure_manifest.json").read_text())
    publication = Path(publication) if publication is not None else run / "publication"
    published = json.loads((publication / "publication_manifest.json").read_text())
    if report["analysis_id"] != figures["analysis_id"] or report["analysis_id"] != published["analysis_id"]:
        raise ValueError("Numerical/rendered/assembled identity mismatch")
    check_hashes(ROOT, figures["plotting_sources"])
    check_hashes(run / "figures", figures["files"])
    check_hashes(ROOT, published["sources"])
    if set(published["files"]) != publication_files() or published["figure_numbers"] != {stem: number for stem, (number, _) in FIGURES.items()}:
        raise ValueError("Invalid Figures 10/11 assembly inventory or numbering")
    check_hashes(publication, published["files"])
    check_hashes(run / "analysis", published["analysis_files"])
    if digest(run / "figures/figure_manifest.json") != published["rendered_manifest_hash"]:
        raise ValueError("Changed rendered manifest")
    return report


def build(run):
    run = Path(run)
    report = validate_run(run)
    record = json.loads((run / "analysis/provenance.json").read_text())
    panels = run / "figures"
    figure_manifest = json.loads((panels / "figure_manifest.json").read_text())
    if figure_manifest["analysis_id"] != report["analysis_id"]:
        raise ValueError("Panel/numerical identity mismatch")
    check_hashes(ROOT, figure_manifest["plotting_sources"])
    check_hashes(panels, figure_manifest["files"])
    out = run / "publication"
    temporary = Path(tempfile.mkdtemp(prefix=".numbered-assembly-", dir=run))
    stage, backup = temporary / "publication", temporary / "previous"
    stage.mkdir()
    installed = False
    try:
        for stem, (_, source) in FIGURES.items():
            for extension in ("pdf", "png"):
                shutil.copy2(panels / f"{source}.{extension}", stage / f"{stem}_panel.{extension}")
        compile_wrappers(write_wrappers(stage, record))
        previous = archive_previous(run, "publication")
        manifest = dict(analysis_id=report["analysis_id"], validation=report, figure_numbers_assigned=True,
            figure_numbers={stem: number for stem, (number, _) in FIGURES.items()},
            publication_target="isolated numbered build; production release managed separately",
            caption_tracking_wording_provisional=True,
            previous_layout=str(previous.relative_to(run)) if previous else None,
            rendered_manifest_hash=digest(panels / "figure_manifest.json"),
            sources={str(p.relative_to(ROOT)): digest(p) for p in (
                Path(__file__), ROOT / "terminal_pareto/publication_wrappers.py",
                ROOT / "terminal_pareto/fig_terminal_cross_species.py", ROOT / "full_tree_pareto/fig_cross_species.py")},
            analysis_files={p.name: digest(p) for p in (run / "analysis").iterdir() if p.is_file()},
            files={name: digest(stage / name) for name in sorted(publication_files())})
        (stage / "publication_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        verify(run, publication=stage)
        if out.exists():
            out.rename(backup)
        stage.rename(out)
        installed = True
        verify(run)
    except Exception:
        if installed:
            out.rename(temporary / "failed")
        if backup.exists():
            backup.rename(out)
        raise
    finally:
        shutil.rmtree(temporary)
    return out


def verify_release(publication=PRODUCTION, *, check_archive=True):
    publication = Path(publication)
    record = json.loads((publication / RELEASE_MANIFEST).read_text())
    if set(record["files"]) != publication_files() | {ASSEMBLY_COPY}:
        raise ValueError("Invalid cross-species production inventory")
    if record["figure_numbers"] != {stem: number for stem, (number, _) in FIGURES.items()}:
        raise ValueError("Invalid production figure numbers")
    check_hashes(publication, record["files"])
    check_hashes(publication, record["preserved_files"])
    if set(inventory(publication)) != set(record["files"]) | set(record["preserved_files"]) | {RELEASE_MANIFEST}:
        raise ValueError("Unexpected or missing production files")
    if any((publication / name).exists() for name in LEGACY_FIG8):
        raise ValueError("Retired full-tree Figure 8 remains in production")
    run = Path(record["run"])
    report = verify(run)
    if report["analysis_id"] != record["analysis_id"]:
        raise ValueError("Production/numerical analysis mismatch")
    if digest(run / "publication/publication_manifest.json") != record["files"][ASSEMBLY_COPY]:
        raise ValueError("Production differs from its numbered source assembly")
    if check_archive and record["previous_publication"]:
        archive = Path(record["previous_publication"])
        check_hashes(archive / "publication", record["previous_files"])
        if json.loads((archive / "manifest.json").read_text())["files"] != record["previous_files"]:
            raise ValueError("Previous production archive manifest mismatch")
    return record


def promote(run, *, production=PRODUCTION):
    """Add Figures 10/11, retire only known Figure 8 assets, preserve other files."""
    if Path(production).is_symlink():
        raise ValueError("Production must be an existing real directory")
    run, target = Path(run).resolve(), Path(production).resolve()
    report = verify(run)
    source = run / "publication"
    if not target.is_dir() or target.is_symlink():
        raise ValueError("Production must be an existing real directory")
    if any(path.is_symlink() for path in target.rglob("*")):
        raise ValueError("Do not copy/retire a production symlink")
    unexpected = [path.name for path in target.iterdir() if path.name.startswith("fig8") and path.name not in LEGACY_FIG8]
    if unexpected:
        raise ValueError(f"Unrecognized legacy Figure 8 assets; inspect before retirement: {unexpected}")
    before = inventory(target)
    retired = {name: before[name] for name in LEGACY_FIG8 if name in before}
    # Existing unrelated release manifests must not be silently invalidated.
    generic_manifest = target / "release_manifest.json"
    if generic_manifest.exists() and set(json.loads(generic_manifest.read_text())["files"]) & set(retired):
        raise ValueError("An existing release pins retired Figure 8; coordinate its manifest first")
    managed = publication_files() | {ASSEMBLY_COPY, RELEASE_MANIFEST}
    preserved = {name: value for name, value in before.items() if name not in managed and name not in retired}
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = target.parent / "legacy/cross_species_releases" / stamp
    stage = Path(tempfile.mkdtemp(prefix=".cross-species-stage-", dir=target.parent))
    previous_moved = installed = False
    try:
        shutil.copytree(target, stage, dirs_exist_ok=True)
        check_hashes(stage, before)
        for name in retired:
            (stage / name).unlink()  # Only staging copies; originals are archived below.
        for name in publication_files():
            shutil.copy2(source / name, stage / name)
        shutil.copy2(source / "publication_manifest.json", stage / ASSEMBLY_COPY)
        record = dict(version="full-tree-cross-species-release-1", run=str(run),
            analysis_id=report["analysis_id"], published_at=stamp, release_promoted=True,
            figure_numbers={stem: number for stem, (number, _) in FIGURES.items()},
            files={name: digest(stage / name) for name in sorted(publication_files() | {ASSEMBLY_COPY})},
            preserved_files=preserved, retired_figure8_files=retired,
            previous_publication=str(archive), previous_files=before,
            scope="add Figures 10/11 and retire four obsolete full-tree Figure 8 assets; no Figure 7 promotion")
        (stage / RELEASE_MANIFEST).write_text(json.dumps(record, indent=2) + "\n")
        verify_release(stage, check_archive=False)
        if inventory(target) != before:
            raise ValueError("Production changed during staging; leave it untouched")
        archive.mkdir(parents=True)
        target.rename(archive / "publication")
        previous_moved = True
        (archive / "manifest.json").write_text(json.dumps(dict(files=before, retired_figure8_files=retired), indent=2) + "\n")
        check_hashes(archive / "publication", before)
        stage.rename(target)
        installed = True
        verify_release(target)
    except Exception:
        if installed:
            shutil.rmtree(target)
        if previous_moved:
            (archive / "publication").rename(target)
            shutil.rmtree(archive)
        raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    print(f"Promoted Figures 10/11; retired {len(retired)} Figure 8 files; preserved {len(preserved)} other files.\nPrevious production: {archive}", flush=True)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--publish", action="store_true", help="Explicitly add Figures 10/11 to production and archive/retire legacy Figure 8")
    parser.add_argument("--verify-production", action="store_true")
    args = parser.parse_args()
    run = run_path(args.run_id)
    if args.verify_production:
        record = verify_release()
        print(f"Verified production Figures 10/11, preserved assets and archive: {record['published_at']}", flush=True)
    elif args.verify:
        print(json.dumps(verify(run), indent=2), flush=True)
    else:
        out = build(run)
        if args.publish:
            out = promote(run)
        print(f"Numbered Figures 10 and 11 ready: {out}", flush=True)


if __name__ == "__main__":
    main()
