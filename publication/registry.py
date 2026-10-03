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

from publication import assembly, provenance, style
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
    figures: tuple[FigureSpec, ...]
    status: str

    def load(self, run_id: str | None = None):
        module, function = self.adapter.split(":")
        loader = getattr(import_module(module), function)
        return loader() if run_id is None else loader(run_id)


FAMILIES = {family.key: family for family in (
    Family("terminal-cross-species",
           "Pooled 187-edge terminal CE protein/CE RNA/CB RNA comparison (cosine molecular distance)",
           "publication.adapters.terminal:cross_species",
           (FigureSpec("8", "fig8_terminal_cross_species_comparison", "fig:terminal_cross_species_8",
                       partial(cs_figures.comparison_figure,
                               title="Terminal-only comparison | {edges} shared terminal edges"),
                       cs_captions.terminal_comparison),
            FigureSpec("9", "fig9_terminal_cross_species_overlays", "fig:terminal_cross_species_9",
                       partial(cs_figures.overlay_figure, title="Pooled terminal fronts and canonical metrics",
                               zoom_includes_closest=False),
                       cs_captions.terminal_overlay),),
           status="active: Figures 8-9 published 2026-10-01"),
    Family("full-tree-cross-species",
           "Terminal-anchored partial forest, 485 cells/454 edges (Euclidean molecular distance)",
           "publication.adapters.full_tree:cross_species",
           (FigureSpec("10", "fig10_full_tree_cross_species_comparison", "fig:full_tree_cross_species_10",
                       partial(cs_figures.comparison_figure,
                               title="Partial-forest layerwise comparison | {edges} shared edges",
                               front_alpha=.65, front_size=16,
                               colorbar_label="Biological-parent edge retention"),
                       cs_captions.full_tree_comparison),
            FigureSpec("11", "fig11_full_tree_cross_species_overlays", "fig:full_tree_cross_species_11",
                       partial(cs_figures.overlay_figure, title="Partial-forest fronts and canonical metrics",
                               zoom_includes_closest=True),
                       cs_captions.full_tree_overlay),),
           status="active: Figures 10-11 published 2026-10-01"),
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
    artifacts, figures = [], {}
    for spec in family.figures:
        layout = spec.render(data, out, f"{spec.stem}_panel")
        wrapper = assembly.write_figure_wrapper(out, spec.stem, number=spec.number, graphic=f"{spec.stem}_panel.pdf",
                                                caption=spec.caption(data.caption), label=spec.tex_label)
        compiled = assembly.compile_wrapper(wrapper, label=spec.number)
        artifacts += [out / name for name in spec.files]
        figures[spec.stem] = dict(number=spec.number, layout=layout, compiled=compiled)
    return provenance.write_manifest(
        out, family=family.key, artifacts=artifacts,
        inputs=dict(run=str(data.run), analysis_id=data.analysis_id, input_files=data.input_files),
        extra=dict(validation=data.validation, figures=figures, canonical_metrics=canonical_records(data),
                   release_promoted=False))


def describe() -> list[dict]:
    return [dict(family=f.key, description=f.description, status=f.status,
                 figures={s.number: s.stem for s in f.figures}) for f in FAMILIES.values()]


def dumps(value) -> str:
    return json.dumps(value, indent=2, default=str)
