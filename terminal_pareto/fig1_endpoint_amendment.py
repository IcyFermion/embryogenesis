"""Export the pooled coordinate schematic for co-author Figure 1 integration.

This is an independent component with no assumed final panel letter. It uses
the same canonical definition drawing as the subtree summary figures; the
drawing lives in ``publication/figures/terminal_canonical.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from publication.figures.terminal_canonical import FIGURE1_STEM as STEM, render_figure1_amendment
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
    validate_existing_context_manifest,
)


def render(out_dir: Path, *, context=None) -> None:
    render_figure1_amendment(out_dir)
    if context is not None:
        (out_dir / "fig1_endpoint_amendment_manifest.json").write_text(
            json.dumps({
                "status": "standalone candidate component for co-author Figure 1",
                "profile": context.profile,
                "run_id": context.run_paths.run_id,
                "context_cache_key": context.cache_key,
                "final_panel_letter": None,
                "component_stem": STEM,
            }, indent=2, sort_keys=True) + "\n")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("pooled_tracking_v1",), required=True)
    parser.add_argument("--run-id")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    context = build_analysis_context(
        args.profile, run_id=args.run_id, output_root=args.output_root)
    validate_existing_context_manifest(context)
    out_dir = (args.out if args.out is not None else
               context.run_paths.display("endpoint"))
    render(out_dir, context=context)
    print(f"Wrote {out_dir / (STEM + '.pdf')}")


if __name__ == "__main__":
    main()
