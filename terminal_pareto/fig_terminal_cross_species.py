"""Retired: Figures 8/9 are drawn by the publication front end.

``cross_species_analysis.py --render-only`` (whose source is hash-pinned by the
numerical cache) still imports ``render``; it now explains where figures live.
"""
from __future__ import annotations

COMMAND = "python -m publication build --family terminal-cross-species --output-dir PATH"


def render(run):
    # The hash-pinned caller reaches this hook after a successful computation
    # as well as for --render-only. Stop before its obsolete figure-path message
    # without turning a completed numerical run into a failed shell command.
    print(f"Numerical results are ready in {run}; build figures with `{COMMAND}`.")
    raise SystemExit(0)
