"""Publication families: numbering, stems, input loaders and owned assets.

A family groups figures built from one validated back-end result. Builds are
cache-only and always write to an explicit, new candidate directory.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import partial
from importlib import import_module
import json
from pathlib import Path
from typing import Callable

from publication import assembly, provenance, style, tables
from publication.canonical_panels import canonical_records
from publication.captions import cross_species as cs_captions
from publication.figures import cross_species as cs_figures


@dataclass(frozen=True)
class FigureSpec:
    number: str
    stem: str
    tex_label: str
    render: Callable  # (data, out, panel_stem) -> layout metadata
    caption: Callable  # (caption metadata) -> TeX caption body

    @property
    def files(self):
        return tuple(f"{self.stem}{suffix}" for suffix in (".pdf", ".tex", "_panel.pdf", "_panel.png"))


@dataclass(frozen=True)
class Family:
    key: str
    description: str
    adapter: str  # "module:function" loaded lazily so listing never imports back ends
    builder: Callable  # (family, data, out) -> dict(artifacts=[paths], figures={stem: metadata}, extra={})
    status: str
    figures: tuple[FigureSpec, ...] = ()  # single-graphic figures (captioned from a spec)
    wrappers: tuple[tuple[str, str], ...] = ()  # (stem, printed number) for multi-panel figures

    def load(self, run_id: str | None = None):
        module, function = self.adapter.split(":")
        loader = getattr(import_module(module), function)
        return loader() if run_id is None else loader(run_id)

    @property
    def numbers(self) -> dict[str, str]:
        return {spec.stem: spec.number for spec in self.figures} | dict(self.wrappers)


def _build_specs(family, data, out):
    """Single-graphic figures: render the panel, then caption, wrap and compile it."""
    artifacts, figures = [], {}
    for spec in family.figures:
        layout = spec.render(data, out, f"{spec.stem}_panel")
        wrapper = assembly.write_figure_wrapper(out, spec.stem, number=spec.number, graphic=f"{spec.stem}_panel.pdf",
                                                caption=spec.caption(data.caption), label=spec.tex_label)
        compiled = assembly.compile_wrapper(wrapper, label=spec.number)
        artifacts += [out / name for name in spec.files]
        figures[spec.stem] = dict(number=spec.number, layout=layout, compiled=compiled)
    return dict(artifacts=artifacts, figures=figures, extra=dict(canonical_metrics=canonical_records(data)))


def _build_terminal_primary(family, data, out):
    """Render every panel, then write and compile each multi-panel wrapper."""
    from publication.captions.terminal import CAPTIONED_FIGURES
    from publication.figures import terminal
    terminal.render_panels(data, out)
    figures = {}
    for stem, number in family.wrappers:
        wrapper = assembly.write_wrapper(out, stem, number=number, body=CAPTIONED_FIGURES[stem][1]())
        figures[stem] = dict(number=number, compiled=assembly.compile_wrapper(wrapper, label=number))
    table = out / f"{tables.TABLE1_STEM}.tex"
    figures[table.stem] = dict(number="1", compiled=assembly.compile_wrapper(
        table, label="1", kind="Table", known_warnings=TABLE1_KNOWN_WARNINGS))
    artifacts = sorted(p for p in out.iterdir() if p.suffix in (".pdf", ".png", ".svg", ".tex"))
    return dict(artifacts=artifacts, figures=figures, extra={})


# Reviewed, pre-existing warnings (one-page outputs): the approved 2026-09-30 Figure 7 build,
# and published Table 1, whose standalone TeX was never compiled in production.
TABLE1_KNOWN_WARNINGS = ("Float too large for page by 3.75319pt",)
FULL_TREE_KNOWN_WARNINGS = {"fig7_ce_full_tree_layerwise": ("Float too large for page by 8.2687pt",)}


def _build_full_tree_pooled(family, data, out):
    """Figure 7 and provisional S4: panels, wrappers (own preamble/numbering) and compilation."""
    from publication.adapters.full_tree import caption_meta
    from publication.captions import full_tree as captions
    from publication.figures import full_tree
    transforms = full_tree.render_panels(data, out)
    meta, figures = caption_meta(data), {}
    for (stem, number), numbering, body in zip(family.wrappers, (captions.FIGURE7_NUMBERING, captions.SUPPLEMENT_NUMBERING),
                                               (captions.figure7, captions.supplement)):
        wrapper = assembly.write_wrapper(out, stem, number=number, body=body(meta), preamble=captions.PREAMBLE,
                                         numbering_tex=numbering)
        figures[stem] = dict(number=number, compiled=assembly.compile_wrapper(
            wrapper, label=number, known_warnings=FULL_TREE_KNOWN_WARNINGS.get(stem, ())))
    artifacts = sorted(p for p in out.iterdir() if p.suffix in (".pdf", ".png", ".tex"))
    return dict(artifacts=artifacts, figures=figures, extra=dict(display_transforms=transforms))


TERMINAL_WRAPPERS = (
    ("fig1_ce_endpoint_normalization_amendment", "1 amendment"), ("fig2_ce_terminal_pareto_main", "2"),
    ("fig3_ce_terminal_pareto_supporting", "3"), ("fig4_ce_subtree_map", "4"), ("fig5_ce_canonical_summary", "5"),
    ("fig6_ce_cell_types", "6"), ("figS1_ce_canonical_summary_cousin_r", "S1"),
    ("figS2_ce_cell_type_cost_gain", "S2"), ("figS3_ce_tracking_robustness", "S3"),
)

FAMILIES = {family.key: family for family in (
    Family("terminal-primary",
           "Pooled 275-edge terminal analysis, 42 subtrees (cosine molecular distance); Figures 1-6, S1-S3, Table 1",
           "publication.adapters.terminal:primary", _build_terminal_primary,
           status="active: published 2026-09-23 (Figure 1 amendment pending co-author integration)",
           wrappers=TERMINAL_WRAPPERS),
    Family("full-tree-pooled",
           "Pooled full tree, 978 nodes/974 edges, layerwise and four other reconstructions (Euclidean molecular distance)",
           "publication.adapters.full_tree:pooled", _build_full_tree_pooled,
           status="approved working build 2026-09-30; NOT promoted (S4 numbering provisional)",
           wrappers=(("fig7_ce_full_tree_layerwise", "7"), ("figs_ce_full_tree_heuristics", "S4"))),
    Family("terminal-cross-species",
           "Pooled 187-edge terminal CE protein/CE RNA/CB RNA comparison (cosine molecular distance)",
           "publication.adapters.terminal:cross_species", _build_specs,
           status="active: Figures 8-9 published 2026-10-01",
           figures=(FigureSpec("8", "fig8_terminal_cross_species_comparison", "fig:terminal_cross_species_8",
                               partial(cs_figures.comparison_figure,
                                       title="Terminal-only comparison | {edges} shared terminal edges"),
                               cs_captions.terminal_comparison),
                    FigureSpec("9", "fig9_terminal_cross_species_overlays", "fig:terminal_cross_species_9",
                               partial(cs_figures.overlay_figure, title="Pooled terminal fronts and canonical metrics",
                                       zoom_includes_closest=False),
                               cs_captions.terminal_overlay))),
    Family("full-tree-cross-species",
           "Terminal-anchored partial forest, 485 cells/454 edges (Euclidean molecular distance)",
           "publication.adapters.full_tree:cross_species", _build_specs,
           status="active: Figures 10-11 published 2026-10-01",
           figures=(FigureSpec("10", "fig10_full_tree_cross_species_comparison", "fig:full_tree_cross_species_10",
                               partial(cs_figures.comparison_figure,
                                       title="Partial-forest layerwise comparison | {edges} shared edges",
                                       front_alpha=.65, front_size=16,
                                       colorbar_label="Biological-parent edge retention"),
                               cs_captions.full_tree_comparison),
                    FigureSpec("11", "fig11_full_tree_cross_species_overlays", "fig:full_tree_cross_species_11",
                               partial(cs_figures.overlay_figure, title="Partial-forest fronts and canonical metrics",
                                       zoom_includes_closest=True),
                               cs_captions.full_tree_overlay))),
)}

PROTECTED = tuple(provenance.ROOT / path for path in (
    "terminal_pareto/output", "full_tree_pareto/output"))


def _check_destination(out: Path):
    out = Path(out).resolve()
    if any(out == p.resolve() or p.resolve() in out.parents for p in PROTECTED):
        raise ValueError(f"Refusing to write a presentation build inside back-end output: {out}")
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory must be new or empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    return out


def build(family_key: str, output_dir: Path, *, run_id: str | None = None) -> dict:
    """Render, caption and compile one family from validated caches only."""
    family = FAMILIES[family_key]
    data = family.load(run_id)
    out = _check_destination(output_dir)
    style.configure()
    result = family.builder(family, data, out)
    return provenance.write_manifest(
        out, family=family.key, artifacts=result["artifacts"],
        inputs=dict(run=str(data.run), analysis_id=data.analysis_id, input_files=data.input_files),
        extra=dict(validation=data.validation, figures=result["figures"], figure_numbers=family.numbers,
                   release_promoted=False, **result["extra"]))


def describe() -> list[dict]:
    return [dict(family=f.key, description=f.description, status=f.status,
                 figures={number: stem for stem, number in f.numbers.items()}) for f in FAMILIES.values()]


def dumps(value) -> str:
    return json.dumps(value, indent=2, default=str)
