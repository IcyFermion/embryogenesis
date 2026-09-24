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

WRAPPERS = (
    "fig1_ce_endpoint_normalization_amendment", "fig2_ce_terminal_pareto_main",
    "fig3_ce_terminal_pareto_supporting", "fig4_ce_subtree_map",
    "fig5_ce_canonical_summary", "fig6_ce_cell_types",
    "figS1_ce_canonical_summary_cousin_r", "figS2_ce_cell_type_cost_gain",
    "figS3_ce_tracking_robustness",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(output_root=DEFAULT_OUTPUT_ROOT):
    publication = Path(output_root) / "publication"
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
    return record


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
    parser.add_argument("--run-id", default="migration_candidate_20260920")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--verify", action="store_true", help="Check current release hashes without rebuilding")
    args = parser.parse_args()
    if args.verify:
        record = verify(args.output_root)
        print(f"Verified {len(record['files'])} published files and all recorded analysis files")
    else:
        promote(args.run_id, args.output_root)


if __name__ == "__main__":
    main()
