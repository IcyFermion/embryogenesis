"""Retired: Figures 8/9 are drawn by the publication front end.

``cross_species_analysis.py --render-only`` (whose source is hash-pinned by the
numerical cache) still imports ``render``; it now explains where figures live.
"""
from __future__ import annotations

COMMAND = "python -m publication build --family terminal-cross-species --output-dir PATH"


def render(run):
    raise SystemExit(f"Figures are no longer drawn into {run}; build them with `{COMMAND}`.")
