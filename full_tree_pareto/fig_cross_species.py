"""Compatibility entrypoint for the partial-forest cross-species figures (Figures 10/11).

Presentation lives in ``publication`` (``figures/cross_species.py``). ``render``
archives the previous working layout and writes the historical panels and
manifest into the run's ``figures/`` directory for
``cross_species_analysis --render-only`` and the numbered assembly.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

from full_tree_pareto.cross_species_analysis import PRIMARY, GEOMETRIES, ROOT, digest
from publication import notation as nt
from publication import provenance, style
from publication.adapters import full_tree as adapter
from publication.canonical_panels import canonical_records
from publication.figures import cross_species as layouts
from publication.registry import FAMILIES

LABELS = nt.CONFIG_LABELS
COMPARISON = "partial_forest_cross_species_comparison"
OVERLAY = "partial_forest_cross_species_overlay"
COMPARISON_OPTIONS = FAMILIES["full-tree-cross-species"].figures[0].render.keywords
OVERLAY_OPTIONS = FAMILIES["full-tree-cross-species"].figures[1].render.keywords


def archive_previous(run, directory):
    source = Path(run) / directory
    if not source.exists() or not any(source.iterdir()):
        return None
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = Path(run) / "layout_history" / stamp / directory
    hashes = {str(p.relative_to(source)): digest(p) for p in source.rglob("*") if p.is_file()}
    shutil.copytree(source, target)
    for name, expected in hashes.items():
        if digest(target / name) != expected:
            raise ValueError(f"Failed layout archive: {name}")
    (target.parent / "manifest.json").write_text(json.dumps(dict(directory=directory, files=hashes), indent=2) + "\n")
    return target.parent


def render(run):
    style.configure()
    run = Path(run)
    data = adapter.cross_species(run)
    archive = archive_previous(run, "figures")
    out = run / "figures"
    out.mkdir(exist_ok=True)
    comparison_limits = layouts.comparison_figure(data, out, COMPARISON, **COMPARISON_OPTIONS)
    comparison_limits["random_rebuild_insets"] = comparison_limits.pop("reference_insets")
    overlay_limits = layouts.overlay_figure(data, out, OVERLAY, **OVERLAY_OPTIONS)
    manifest = dict(analysis_id=data.analysis_id, comparison=COMPARISON, overlay=OVERLAY,
        comparison_rows=list(GEOMETRIES), comparison_columns=list(PRIMARY), retention_scale=[0, 1],
        comparison_geometry_labels="bold second line beneath each configuration title",
        comparison_title_tracking_counts=False,
        comparison_limits=comparison_limits, overlay_limits=overlay_limits,
        comparison_natural_to_maximum_connectors=False, overlay_closest_point_connectors=True,
        canonical_panel="C", canonical_metrics=canonical_records(data), notation=nt.snapshot(),
        previous_layout=str(archive.relative_to(run)) if archive else None,
        plotting_sources={str(Path(__file__).relative_to(ROOT)): digest(Path(__file__)),
                          **provenance.presentation_sources()},
        files={p.name: digest(p) for p in out.iterdir() if p.suffix in (".png", ".pdf")})
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return out
