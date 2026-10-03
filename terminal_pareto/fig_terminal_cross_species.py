"""Compatibility entrypoint for the terminal cross-species figures (Figures 8/9).

Presentation lives in ``publication`` (``figures/cross_species.py``,
``captions/cross_species.py``). ``render(run)`` keeps writing the historical
panel set, caption and manifest into the run's ``figures/`` directory for
``cross_species_analysis.py --render-only`` and the numbered assembly.
"""

from __future__ import annotations

import json
from pathlib import Path

from publication import notation as nt
from publication import provenance, style
from publication.adapters import terminal as adapter
from publication.canonical_panels import canonical_records
from publication.captions.cross_species import TRACKING_CAPTION_STATUS, TRACKING_CAVEAT  # noqa: F401
from publication.figures import cross_species as layouts
from publication.registry import FAMILIES
from terminal_pareto.cross_species_analysis import LABELS, PRIMARY, file_hash  # noqa: F401

COMPARISON_OPTIONS = FAMILIES["terminal-cross-species"].figures[0].render.keywords
OVERLAY_OPTIONS = FAMILIES["terminal-cross-species"].figures[1].render.keywords


def caption(record) -> str:
    d_lp, d_cp = nt.plain("natural_front_distance"), nt.plain("null_front_distance")
    return (
        "Terminal parentage across CE protein, CE RNA and CB AF16 RNA. The same 187 natural terminal edges and parent capacities "
        "are used in all panels. In the main comparison, A-C are existing 3D and D-F are 2D (XY); columns are CE protein, CE RNA "
        "and CB RNA. All six panels share axis limits and one retention scale. The separate overlay figure places 3D and XY side "
        "by side in A-B, sharing main-axis and near-natural zoom limits, species colors and line styles. "
        f"Hollow diamonds mark the closest attained assignment {nt.plain('closest_point')}, with dashed connectors from natural lineage "
        f"to {nt.plain('closest_point')}, not to maximum retention. Panel C groups the three configurations under 3D and XY and shows canonical "
        f"position u, natural distance {d_lp} and first-cousin-null-mean distance {d_cp} in separate metric columns. "
        f"Both distances use the same {nt.plain('closest_point')}, and u is front arc length from the travel optimum normalized by total arc length. "
        "Travel is the equal mean of "
        "per-embryo pairwise distances normalized by fixed natural totals: "
        "three CE embryos at cutoffs 255/247/225 on the published 275-edge reference, and two CB AF16 embryos at cutoffs 148/156 "
        f"on their {record['cb_reference_size']}-edge tracking/RNA intersection. Protein uses the frozen top-20 z-scored reporters; "
        "RNA uses the frozen shared 20 TFs, cosine distance on stored values. Each front uses its own endpoint spans; within a "
        "panel all nulls and landmarks share those anchors. Colors encode biological-parent edge retention, crosses natural lineage, "
        "and outlined circles observed maximum retention. The first-cousin mean is analytic; other mean marks summarize display draws. "
        "Full-random nulls use insets with the same coordinates and separate limits. Optimization uses exact first-cousin null SDs "
        "and 301 weights with checked endpoint tie handling; 1,201-weight checks are cached separately. Lines connect attained "
        f"assignments. {TRACKING_CAVEAT} {TRACKING_CAPTION_STATUS} "
        "Developmental alignment and RNA measurement provenance remain unresolved. The XY row "
        "drops z and recomputes travel, nulls, endpoints and assignments; tracking spread does not estimate molecular uncertainty.\n"
    )


def render(run):
    style.configure()
    run = Path(run)
    out = run / "figures"
    out.mkdir(parents=True, exist_ok=True)
    data = adapter.cross_species(run)
    record = json.loads((run / "analysis/provenance.json").read_text())
    for geometry, suffix in (("raw3d", "existing_3d"), ("xy", "xy")):
        layouts.comparison_figure(data, out, f"terminal_cross_species_{suffix}", geometries=(geometry,),
                                  **COMPARISON_OPTIONS)
        layouts.save(layouts.compose_single_overlay(data, geometry), out, f"terminal_cross_species_overlay_{suffix}")
    layouts.comparison_figure(data, out, "terminal_cross_species_comparison", **COMPARISON_OPTIONS)
    layouts.overlay_figure(data, out, "terminal_cross_species_overlay", **OVERLAY_OPTIONS)
    layouts.save(layouts.compose_tracking(data), out, "terminal_cross_species_tracking_sensitivity")
    (out / "caption.txt").write_text(caption(record))
    manifest = dict(analysis_id=record["analysis_id"], shared_axes=True, retention_color_scale=[0, 1],
        primary_comparison="terminal_cross_species_comparison", comparison_rows=["raw3d", "xy"],
        comparison_columns=list(PRIMARY), primary_overlay="terminal_cross_species_overlay", overlay_panels=["raw3d", "xy"],
        comparison_geometry_labels="bold second line beneath each configuration title",
        comparison_title_tracking_counts=False,
        shared_overlay_zoom_limits=True,
        comparison_natural_to_maximum_connectors=False,
        overlay_closest_point_connectors=True, canonical_panel="C",
        canonical_metrics=canonical_records(data),
        endpoint_normalization="separate configuration anchors; shared anchors for transfers within each panel",
        notation=nt.snapshot(),
        plotting_sources={"terminal_pareto/fig_terminal_cross_species.py": file_hash(Path(__file__)),
                          **provenance.presentation_sources()},
        files={path.name: file_hash(path) for path in sorted(out.iterdir()) if path.suffix in (".png", ".pdf", ".txt")})
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
