"""Export the pooled coordinate schematic for co-author Figure 1 integration.

This is an independent component with no assumed final panel letter. It uses
the same canonical definition drawing as the subtree summary figures.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto import plot_style as ps
from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
    validate_existing_context_manifest,
)
from terminal_pareto.fig5_figs1_ce_canonical_summary import plot_definition_dlp

STEM = "fig1_endpoint_normalization_amendment"


def render(out_dir: Path, *, context=None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    ps.configure()
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    fig.subplots_adjust(left=.16, right=.84, bottom=.18, top=.84)
    plot_definition_dlp(ax, panel_letter="")
    for label in ax.texts:
        label.set_fontsize(label.get_fontsize() * 1.35)
    ax.set_title("Endpoint coordinates and canonical metrics",
                 fontsize=13, loc="left", pad=12)
    with plt.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none"}):
        for extension in ("pdf", "png", "svg"):
            fig.savefig(out_dir / f"{STEM}.{extension}", dpi=300,
                        bbox_inches="tight", facecolor="white")
    plt.close(fig)
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
