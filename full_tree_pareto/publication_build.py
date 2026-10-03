"""Build Figure 7 and its supplement from isolated, validated pooled caches.

python -m full_tree_pareto.publication_build [--layout-only] [--publish]

Drawing and captions live in ``publication/`` (``figures/full_tree.py``,
``captions/full_tree.py``); this module keeps the sweep settings, cache
building, working-layout archive and release mechanics.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import cousin_references as cr
from publication import assembly
from publication.adapters import full_tree as full_tree_adapter
from publication.captions import full_tree as ft_captions
from publication.figures import full_tree as full_tree_figures
from publication.registry import FULL_TREE_KNOWN_WARNINGS
from publication.figures.full_tree import plot_comparison, plot_rounds, retention_front  # noqa: F401  (re-exported)
from terminal_pareto.front_coordinates import EndpointTransform

# Single source for the sweep/draw settings used by builds, checkpoints and
# captions. pooled_analysis.py is hash-pinned in cache identities, so these live
# here rather than there; changing them requires a new run directory.
INTERVALS = 300
DRAWS = 10000
SEED = 42
DEFAULT_WORKERS = 24
DISPLAY_DRAWS = 1000
WRAPPERS = ("fig7_ce_full_tree_layerwise", "figs_ce_full_tree_heuristics")
MAIN_REFERENCES = (pa.NULLS[0], *cr.REFERENCES)
SUPPLEMENT_REFERENCES = (*MAIN_REFERENCES, *pa.NULLS[1:])


def endpoint(frame, scope):
    a = frame.loc[frame.weight_index.idxmax()]
    b = frame.loc[frame.weight_index.idxmin()]
    return EndpointTransform.from_endpoints(
        reference_analysis_id=f"{pa.PROFILE}:layerwise:{scope}",
        travel_optimum_assignment_id=f"layerwise:{scope}:{int(a.weight_index)}",
        state_optimum_assignment_id=f"layerwise:{scope}:{int(b.weight_index)}",
        travel_optimum_costs=(a.travel, a.state),
        state_optimum_costs=(b.travel, b.state))


def write_wrappers(ctx, out):
    """Figure 7 and S4 wrappers; captions live in ``publication/captions/full_tree.py``."""
    meta = dict(edges=len(ctx.edges), leaves=len(ctx.leaves), internal=len(ctx.internal),
                round_edges=[len(r) for r in ctx.layers], weights=INTERVALS + 1, draws=DRAWS,
                display_draws=DISPLAY_DRAWS)
    for stem, numbering, body in ((WRAPPERS[0], ft_captions.FIGURE7_NUMBERING, ft_captions.figure7(meta)),
                                  (WRAPPERS[1], ft_captions.SUPPLEMENT_NUMBERING, ft_captions.supplement(meta))):
        assembly.write_wrapper(out, stem, number="", body=body, preamble=ft_captions.PREAMBLE, numbering_tex=numbering)


def compile_wrappers(out):
    for stem, label in zip(WRAPPERS, ("7", "S4")):
        assembly.compile_wrapper(Path(out) / f"{stem}.tex", label=label,
                                 known_warnings=FULL_TREE_KNOWN_WARNINGS.get(stem, ()))


def archive_working_publication(run):
    """Retain a hash-verified layout snapshot before overwriting review assets."""
    run = Path(run)
    source = run / "publication"
    if not source.exists() or not any(source.iterdir()):
        return None
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = run / "layout_history" / stamp
    hashes = {str(p.relative_to(source)): pa.digest(p)
              for p in source.rglob("*") if p.is_file()}
    shutil.copytree(source, archive / "publication")
    for name, sha in hashes.items():
        if pa.digest(archive / "publication" / name) != sha:
            raise ValueError("Working layout archive verification failed")
    pa.write_json(archive / "manifest.json", dict(files=hashes,
                  analysis_files={str(p.relative_to(run / "analysis")): pa.digest(p)
                                  for p in (run / "analysis").rglob("*") if p.is_file()}))
    print(f"Preserved preceding layout: {archive}", flush=True)
    return archive


def render(ctx, fronts, layers, nulls, run):
    data = full_tree_adapter.pooled_from_results(run, ctx, fronts, layers, nulls)
    archive_working_publication(run)
    out = Path(run) / "publication"
    out.mkdir(parents=True, exist_ok=True)
    transforms = full_tree_figures.render_panels(data, out)
    pa.write_json(out / "display_transforms.json", transforms)
    write_wrappers(ctx, out)
    compile_wrappers(out)
    return out


def verify_release(publication):
    publication = Path(publication)
    record = json.loads((publication / "release_manifest.json").read_text())
    for filename, sha in record["files"].items():
        if pa.digest(publication / filename) != sha:
            raise ValueError(f"Changed publication file: {filename}")
    run = Path(record["run"])
    for filename, sha in record["analysis_files"].items():
        if pa.digest(run / "analysis" / filename) != sha:
            raise ValueError(f"Changed analysis file: {filename}")
    return record


def sweep_settings():
    return dict(intervals=INTERVALS, draws=DRAWS, seed=SEED)


def promote(run):
    """Explicit publication replacement; archive/hash all former assets first."""
    run = Path(run).resolve()
    source = run / "publication"
    root = pa.ROOT / "full_tree_pareto/output"
    target = root / "publication"
    if (target / "cross_species_release_manifest.json").exists():
        raise ValueError("Production includes Figures 10/11. Coordinate an additive Figure 7 release "
                         "that preserves them and updates their preserved-file hashes.")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    # Revalidate source/config and replay parent arrays before promotion.
    ctx, _, _, _ = pa.build(run=run, layout_only=True, **sweep_settings())
    cr.build(ctx, run, draws=DRAWS, seed=SEED, layout_only=True)
    for stem in WRAPPERS:
        info = subprocess.check_output(["pdfinfo", str(source / f"{stem}.pdf")], text=True)
        if int(next(x.split(":")[1] for x in info.splitlines() if x.startswith("Pages:"))) != 1:
            raise ValueError("Invalid wrapper page count")
    # Stage and verify the complete replacement before moving the live bundle.
    # Every failure path removes the stage and any partial archive.
    stage = Path(tempfile.mkdtemp(prefix=".pooled-stage-", dir=root))
    archive = root / "legacy/releases" / stamp
    previous = None
    created_archive = swapped = False
    try:
        for p in source.iterdir():
            if p.suffix in (".pdf", ".png", ".tex", ".json"):
                shutil.copy2(p, stage / p.name)
        record = dict(profile=pa.PROFILE, run=str(run), published_at=stamp,
                      files={p.name: pa.digest(p) for p in stage.iterdir() if p.is_file()},
                      analysis_files={str(p.relative_to(run / "analysis")): pa.digest(p)
                                      for p in (run / "analysis").rglob("*") if p.is_file()})
        pa.write_json(stage / "release_manifest.json", record)
        verify_release(stage)
        if target.exists():
            archive.mkdir(parents=True)
            created_archive = True
            hashes = {str(p.relative_to(target)): pa.digest(p) for p in sorted(target.rglob("*")) if p.is_file()}
            previous = archive / "publication"
            target.rename(previous)
            pa.write_json(archive / "manifest.json", hashes)
            for name, sha in hashes.items():
                if pa.digest(previous / name) != sha:
                    raise ValueError("Legacy archive verification failed")
        stage.rename(target)
        swapped = True
        verify_release(target)
    except Exception:
        if swapped:
            shutil.rmtree(target)
        if previous is not None and previous.exists():
            previous.rename(target)
        if created_archive:
            shutil.rmtree(archive, ignore_errors=True)
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return target


def positive_workers(value):
    """Validate the shared CLI worker option without affecting cache identity."""
    try:
        workers = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("workers must be a positive integer") from exc
    if workers < 1:
        raise argparse.ArgumentTypeError("workers must be a positive integer")
    return workers


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=pa.DEFAULT_RUN)
    parser.add_argument("--workers", type=positive_workers, default=DEFAULT_WORKERS,
                        help=f"Parallel solver processes (default: {DEFAULT_WORKERS}; use 1 for serial execution).")
    parser.add_argument("--layout-only", action="store_true")
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--verify", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.verify:
        verify_release(pa.ROOT / "full_tree_pareto/output/publication")
        print("Verified pooled full-tree publication and analysis hashes")
        return
    result = pa.build(run=args.run, workers=args.workers, layout_only=args.layout_only,
                      **sweep_settings())
    result[3].update(cr.build(result[0], args.run, draws=DRAWS, seed=SEED,
                              layout_only=args.layout_only))
    out = render(*result, args.run)
    if args.publish:
        out = promote(args.run)
    print(f"Figures written to {out}")


if __name__ == "__main__":
    main()
