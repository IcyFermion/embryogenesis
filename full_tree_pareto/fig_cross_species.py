"""Retired: Figures 10/11 are drawn by the publication front end.

``cross_species_analysis --render-only`` (whose source is hash-pinned by the
numerical cache) still imports ``render``; it now explains where figures live.
"""
from __future__ import annotations

COMMAND = "python -m publication build --family full-tree-cross-species --output-dir PATH"


def render(run):
    # The hash-pinned caller also reaches this hook after a successful build.
    # Exit before its obsolete figure-path message, with a successful status.
    print(f"Numerical results are ready in {run}; build figures with `{COMMAND}`.")
    raise SystemExit(0)
