"""Validate and publish a frozen pooled figure bundle; preserve earlier releases."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from terminal_pareto.analysis_context import DEFAULT_OUTPUT_ROOT
from publication.provenance import stale_files

WRAPPERS = (
    "fig1_ce_endpoint_normalization_amendment", "fig2_ce_terminal_pareto_main",
    "fig3_ce_terminal_pareto_supporting", "fig4_ce_subtree_map",
    "fig5_ce_canonical_summary", "fig6_ce_cell_types",
    "figS1_ce_canonical_summary_cousin_r", "figS2_ce_cell_type_cost_gain",
    "figS3_ce_tracking_robustness",
)
CROSS_SPECIES_STEMS = ("fig8_terminal_cross_species_comparison", "fig9_terminal_cross_species_overlays")
CROSS_SPECIES_MANIFEST = "cross_species_publication_manifest.json"


def cross_species_files():
    return {f"{stem}{suffix}" for stem in CROSS_SPECIES_STEMS
            for suffix in (".pdf", ".tex", "_panel.pdf", "_panel.png")} | {CROSS_SPECIES_MANIFEST}


def file_inventory(directory):
    directory = Path(directory)
    return {str(p.relative_to(directory)): digest(p) for p in directory.rglob("*") if p.is_file()}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(output_root=DEFAULT_OUTPUT_ROOT, *, publication=None):
    publication = Path(publication) if publication is not None else Path(output_root) / "publication"
    record = json.loads((publication / "release_manifest.json").read_text())
    if record["profile"] != "pooled_tracking_v1" or record["display"] != "endpoint":
        raise ValueError("Unexpected publication profile/display")
    mismatches = [name for name, value in record["files"].items()
                  if not (publication / name).is_file()
                  or digest(publication / name) != value]
    for name, value in record["analysis_files"].items():
        path = Path(output_root) / name
        if not path.is_file() or digest(path) != value:
            mismatches.append(name)
    if mismatches:
        raise ValueError(f"Published release mismatch: {mismatches}")
    if "cross_species" in record:
        comparison = record["cross_species"]
        if set(comparison["files"]) != cross_species_files():
            raise ValueError("Invalid cross-species release inventory")
        if comparison["figure_numbers"] != dict(zip(CROSS_SPECIES_STEMS, (8, 9))):
            raise ValueError("Invalid cross-species figure numbers")
        if any(record["files"].get(name) != sha for name, sha in comparison["files"].items()):
            raise ValueError("Cross-species artifact hashes disagree with release")
        if any(record["analysis_files"].get(name) != sha for name, sha in comparison["analysis_files"].items()):
            raise ValueError("Cross-species analysis hashes disagree with release")
        assembly = json.loads((publication / CROSS_SPECIES_MANIFEST).read_text())
        if assembly["analysis_id"] != comparison["analysis_id"] or assembly["figure_numbers"] != comparison["figure_numbers"]:
            raise ValueError("Cross-species assembly identity mismatch")
        if assembly["files"] != {name: sha for name, sha in comparison["files"].items() if name != CROSS_SPECIES_MANIFEST}:
            raise ValueError("Cross-species assembly asset hashes mismatch")
    return record


def validate_cross_species_source(run):
    """Replay the comparison, then verify its frozen panels and numbered assembly."""
    from terminal_pareto.cross_species_analysis import validate_run
    run = Path(run)
    report = validate_run(run)
    source = run / "publication"
    record = json.loads((source / "publication_manifest.json").read_text())
    rendered = json.loads((run / "figures/figure_manifest.json").read_text())
    if record["analysis_id"] != report["analysis_id"] or rendered["analysis_id"] != report["analysis_id"]:
        raise ValueError("Comparison numerical/rendered/assembly identity mismatch")
    if set(record["files"]) != cross_species_files() - {CROSS_SPECIES_MANIFEST}:
        raise ValueError("Invalid numbered comparison assembly inventory")
    if record["figure_numbers"] != dict(zip(CROSS_SPECIES_STEMS, (8, 9))):
        raise ValueError("Invalid comparison numbering")
    # Frozen panels, assembly and numerical inputs must be unchanged. Live
    # presentation source may have moved on since assembly; report it only.
    for base, hashes in ((run / "figures", rendered["files"]), (source, record["files"]),
                         (run / "analysis", record["analysis_files"])):
        for name, sha in hashes.items():
            if digest(base / name) != sha:
                raise ValueError(f"Changed comparison source/asset: {base / name}")
    stale = sorted(set(stale_files(record["sources"])) | set(stale_files(rendered["plotting_sources"])))
    if stale:
        print(f"Presentation source changed since the comparison assembly (artifacts unchanged): {stale}", flush=True)
    if digest(run / "figures/figure_manifest.json") != record["rendered_manifest_hash"]:
        raise ValueError("Changed comparison render manifest")
    if file_inventory(run / "analysis") != record["analysis_files"]:
        raise ValueError("Changed comparison analysis inventory")
    for stem, number in zip(CROSS_SPECIES_STEMS, (8, 9)):
        pdf = source / f"{stem}.pdf"
        info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
        pages = next(line for line in info.splitlines() if line.startswith("Pages:"))
        content = subprocess.check_output(["pdftotext", str(pdf), "-"], text=True)
        if int(pages.split(":", 1)[1]) != 1 or f"Figure {number}:" not in content:
            raise ValueError(f"Invalid one-page numbered comparison: {pdf}")
    return record


def promote_cross_species(run_id=None, output_root=DEFAULT_OUTPUT_ROOT):
    """Add only Figures 8/9 to the existing release; keep the source run intact."""
    from terminal_pareto.cross_species_analysis import DEFAULT_RUN
    output_root = Path(output_root).resolve()
    run = output_root / "runs/cross_species_terminal_v1" / (run_id or DEFAULT_RUN)
    assembled = validate_cross_species_source(run)
    publication = output_root / "publication"
    if publication.is_symlink() or any(p.is_symlink() for p in publication.rglob("*")):
        raise ValueError("Do not replace a publication symlink")
    record = verify(output_root)
    before = file_inventory(publication)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = output_root / "legacy/releases" / stamp
    stage = Path(tempfile.mkdtemp(prefix=".comparison-stage-", dir=output_root))
    previous_moved = installed = False
    try:
        shutil.copytree(publication, stage, dirs_exist_ok=True)
        if file_inventory(stage) != before:
            raise ValueError("Staged publication copy mismatch")
        for name in assembled["files"]:
            shutil.copy2(run / "publication" / name, stage / name)
        shutil.copy2(run / "publication/publication_manifest.json", stage / CROSS_SPECIES_MANIFEST)
        comparison = dict(run=str(run.relative_to(output_root)), analysis_id=assembled["analysis_id"],
            published_at_utc=stamp, figure_numbers=assembled["figure_numbers"],
            validation=assembled["validation"], previous_publication=str((archive / "publication").relative_to(output_root)),
            files={name: digest(stage / name) for name in sorted(cross_species_files())},
            analysis_files={str((run / "analysis" / name).relative_to(output_root)): sha
                            for name, sha in assembled["analysis_files"].items()})
        record["cross_species"] = comparison
        record["files"].update(comparison["files"])
        record["analysis_files"].update(comparison["analysis_files"])
        (stage / "release_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
        verify(output_root, publication=stage)
        if file_inventory(publication) != before:
            raise ValueError("Publication changed during staging; leave it untouched")
        archive.mkdir(parents=True)
        publication.rename(archive / "publication")
        previous_moved = True
        (archive / "manifest.json").write_text(json.dumps(dict(files=before), indent=2) + "\n")
        if file_inventory(archive / "publication") != before:
            raise ValueError("Previous release archive mismatch")
        stage.rename(publication)
        installed = True
        verify(output_root)
    except Exception:
        if installed:
            publication.rename(stage)
        if previous_moved:
            (archive / "publication").rename(publication)
            shutil.rmtree(archive)
        raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    print(f"Published Figures 8/9 alongside existing terminal figures: {publication}\nPrevious release: {archive}")
    return publication


def promote(run_id, output_root=DEFAULT_OUTPUT_ROOT):
    from terminal_pareto.validate_pooled_migration import validate
    output_root = Path(output_root).resolve()
    report = validate(run_id, output_root)
    run = output_root / "runs" / "pooled_tracking_v1" / run_id
    source = run / "publication" / "endpoint"
    for stem in WRAPPERS:
        pdf = source / f"{stem}.pdf"
        info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
        pages = next(line for line in info.splitlines() if line.startswith("Pages:"))
        if int(pages.split(":", 1)[1]) != 1:
            raise ValueError(f"Expected one-page wrapper: {pdf}")
        if not (source / f"{stem}.tex").is_file():
            raise FileNotFoundError(stem)
    publication = output_root / "publication"
    # A later pooled-layout release must not silently discard the comparison.
    existing = (verify(output_root) if (publication / "release_manifest.json").exists()
                and "cross_species" in json.loads((publication / "release_manifest.json").read_text()) else None)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    previous = None
    stage = Path(tempfile.mkdtemp(prefix=".publication-stage-", dir=output_root))
    try:
        for path in source.iterdir():
            if path.is_file() and path.suffix in {".pdf", ".png", ".svg", ".tex", ".csv", ".json"}:
                shutil.copy2(path, stage / path.name)
        record = {
            "published_at_utc": stamp, "profile": "pooled_tracking_v1",
            "display": "endpoint", "run_id": run_id,
            "source": str(source.relative_to(output_root)),
            "context_cache_key": report["context_cache_key"],
            "validation_checks": len(report["checks"]),
            "files": {p.name: digest(p) for p in sorted(stage.iterdir())},
            "analysis_files": {str(p.relative_to(output_root)): digest(p)
                               for p in sorted((run / "analysis").rglob("*")) if p.is_file()},
        }
        if existing is not None:
            comparison = existing["cross_species"]
            for name in comparison["files"]:
                shutil.copy2(publication / name, stage / name)
            record["cross_species"] = comparison
            record["files"].update(comparison["files"])
            record["analysis_files"].update(comparison["analysis_files"])
        (stage / "release_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
        if publication.exists():
            if (publication / "release_manifest.json").exists():
                previous = output_root / "legacy" / "releases" / stamp / "publication"
            else:
                previous = output_root / "legacy" / "embryo1" / "publication"
            if previous.exists():
                raise FileExistsError(f"Will not overwrite legacy material: {previous}")
            previous.parent.mkdir(parents=True, exist_ok=True)
            publication.rename(previous)
        try:
            stage.rename(publication)
            verify(output_root)
        except Exception:
            if publication.exists():
                publication.rename(stage)
            if previous is not None:
                previous.rename(publication)
            raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    print(f"Published pooled endpoint figures: {publication}")
    return publication


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", help="Run in the selected pooled or cross-species profile")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--verify", action="store_true", help="Check current release hashes without rebuilding")
    parser.add_argument("--cross-species", action="store_true", help="Add validated Figures 8/9 without rebuilding or replacing other figures")
    args = parser.parse_args()
    if args.verify:
        record = verify(args.output_root)
        print(f"Verified {len(record['files'])} published files and all recorded analysis files")
    elif args.cross_species:
        promote_cross_species(args.run_id, args.output_root)
    else:
        promote(args.run_id or "migration_candidate_20260920", args.output_root)


if __name__ == "__main__":
    main()
