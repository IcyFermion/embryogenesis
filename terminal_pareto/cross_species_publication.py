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
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto.cross_species_analysis import DEFAULT_RUN, file_hash, run_path, validate_run
from publication import assembly, provenance
from publication.captions import cross_species as cs_captions

FIGURES = {
    "fig8_terminal_cross_species_comparison": (8, "terminal_cross_species_comparison"),
    "fig9_terminal_cross_species_overlays": (9, "terminal_cross_species_overlay"),
}


def captions(record):
    """Publication prose without species rankings or absolute-cost claims."""
    meta = dict(edges=record["cohort_size"], cb_reference_size=record["cb_reference_size"])
    return {8: cs_captions.terminal_comparison(meta), 9: cs_captions.terminal_overlay(meta)}


def write_wrappers(out, record):
    text = captions(record)
    return [assembly.write_figure_wrapper(Path(out), stem, number=str(number), graphic=f"{stem}_panel.pdf",
                                          caption=text[number], label=f"fig:terminal_cross_species_{number}")
            for stem, (number, _) in FIGURES.items()]


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
    for wrapper in wrappers:
        assembly.compile_wrapper(wrapper, label=str(FIGURES[wrapper.stem][0]))


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
        sources={str(Path(__file__).relative_to(ROOT)): file_hash(Path(__file__)),
                 **provenance.presentation_sources()},
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
