"""Write pooled-profile LaTeX wrappers beside generated publication panels.

Compatibility module: captions live in ``publication/captions/terminal.py`` and
the preamble/numbering in ``publication/assembly.py``.
"""

from __future__ import annotations

from pathlib import Path
import shutil

from publication import assembly
from publication.assembly import PREAMBLE  # noqa: F401  (re-exported)
from publication.captions.terminal import CAPTIONED_FIGURES, NULL_SCHEMATIC


def _write(out_dir: Path, stem: str, body: str) -> Path:
    """Legacy helper: derive the printed number from the stem (fig1 -> "1 amendment")."""
    number = stem.split("_")[0].removeprefix("fig")
    return assembly.write_wrapper(out_dir, stem, number="1 amendment" if number == "1" else number, body=body)


def write_pooled_wrappers(out_dir: Path | str) -> list[Path]:
    """Generate candidate wrappers without touching accepted legacy files."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    null_schematic = NULL_SCHEMATIC
    null_schematic_copy = out_dir / null_schematic.name
    if null_schematic.resolve() != null_schematic_copy.resolve():
        shutil.copy2(null_schematic, null_schematic_copy)
    return [assembly.write_wrapper(out_dir, stem, number=number, body=caption())
            for stem, (number, caption) in CAPTIONED_FIGURES.items()]
