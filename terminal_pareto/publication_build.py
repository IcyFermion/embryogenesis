"""Build the pooled terminal publication figures in a versioned work area.

Build stages write only to the selected run directory.
Use --publish to validate and promote the build to the current publication set.
Legacy no-argument entrypoints write only to the legacy rebuild area.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from terminal_pareto.analysis_context import (
    DEFAULT_OUTPUT_ROOT,
    build_analysis_context,
    validate_existing_context_manifest,
)
from terminal_pareto.publication_wrappers import write_pooled_wrappers


def _run(command: list[str], *, log: Path, env: dict[str, str],
         cwd: Path = ROOT) -> None:
    result = subprocess.run(
        command, cwd=cwd, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(result.stdout)
    if result.returncode:
        raise subprocess.CalledProcessError(
            result.returncode, command, output=result.stdout)


def _script(name: str, *args: str) -> list[str]:
    return [sys.executable, str(ROOT / "terminal_pareto" / name), *args]


def _render_single_panel_fig2(context, out_dir: Path) -> None:
    """Reuse the production global cache to export Figure 2 without a letter."""
    import matplotlib.pyplot as plt
    from terminal_pareto.fig2_fig3_ce_terminal_pareto import (
        _display_data,
        plot_main,
    )
    from terminal_pareto.front_coordinates import null_sd_coordinates
    from terminal_pareto.global_analysis import load_global_analysis

    result = load_global_analysis(context)
    names = {
        1: "first_cousin", 2: "second_cousin",
        3: "third_cousin", "full": "full_random",
    }
    nulls = {
        key: null_sd_coordinates(
            *result.null_raw[name], natural_costs=result.natural_costs,
            null_stds=result.null_stds)
        for key, name in names.items()
    }
    display = _display_data(
        result.twr, nulls, mode="endpoint", result=result, context=context)
    fig = plot_main(
        result.twr, nulls, display=display, out_dir=out_dir,
        panel_letter=None)
    plt.close(fig)


def build(profile: str, run_id: str, output_root: Path, *,
          compile_tex: bool = True, layout_only: bool = False) -> Path:
    context = build_analysis_context(
        profile, run_id=run_id, output_root=output_root)
    if (context.run_paths.analysis / "analysis_manifest.json").exists():
        validate_existing_context_manifest(context)
    elif layout_only:
        raise FileNotFoundError(
            "Layout-only build requires the existing pooled analysis cache")
    else:
        context.write()
    paths = context.run_paths.create()
    primary = paths.display("endpoint")
    primary.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    mplconfig = paths.validation / "matplotlib"
    mplconfig.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mplconfig)
    common = ["--profile", profile, "--run-id", run_id,
              "--output-root", str(output_root)]
    commands = [
        ("fig2_fig3_endpoint", _script(
            "fig2_fig3_ce_terminal_pareto.py", *common,
            "--display", "endpoint")),
        ("fig2_fig3_null_sd", _script(
            "fig2_fig3_ce_terminal_pareto.py", *common,
            "--display", "null_sd")),
        ("fig2_fig3_percent_natural", _script(
            "fig2_fig3_ce_terminal_pareto.py", *common,
            "--display", "percent_natural")),
        ("subtree_analysis", _script(
            "subtree_analysis.py", *common, "--min-cells", "12")),
        ("canonical_metrics", _script(
            "fig5_table1_ce_canonical_metrics.py", *common,
            "--min-cells", "12", "--iteration", "300")),
        ("fig4", _script(
            "fig4_ce_subtree_map.py", *common, "--min-cells", "12")),
        ("fig5_figs1", _script(
            "fig5_figs1_ce_canonical_summary.py", *common)),
        ("fig6a_figs2", _script(
            "fig6a_figs2_ce_cell_types.py", *common,
            "--iteration", "300")),
        ("fig6bc", _script(
            "fig6bc_ce_within_type.py", *common,
            "--iteration", "300", "--display", "endpoint")),
    ]
    if layout_only:
        _render_single_panel_fig2(context, primary)
        _run(_script("fig5_figs1_ce_canonical_summary.py", *common, "--primary-only"),
             log=paths.validation / "fig5_layout.log", env=env)
    else:
        for label, command in commands:
            _run(command, log=paths.validation / f"{label}.log", env=env)

    organization_commands = [
        ("fig1_endpoint_amendment", _script(
            "fig1_endpoint_amendment.py", *common, "--out", str(primary))),
        ("figs3_cross_geometry", _script(
            "figS3_cross_geometry.py",
            "--run-id", run_id, "--output-root", str(output_root),
            "--tracking-cache", str(ROOT / "terminal_pareto" / "output"
                                      / "tracking_geometry_sensitivity"
                                      / "assignments_and_costs.npz"),
            "--out", str(primary), "--publication")),
    ]
    for label, command in organization_commands:
        _run(command, log=paths.validation / f"{label}.log", env=env)

    wrappers = write_pooled_wrappers(primary)
    if compile_tex:
        pdflatex = shutil.which("pdflatex")
        tectonic = (shutil.which("tectonic")
                    or str(Path("/opt/codex-desktop/resources/plugins"
                                "/openai-bundled/plugins/latex/bin/tectonic")))
        if pdflatex is None and not Path(tectonic).is_file():
            raise RuntimeError(
                "Publication wrapper compilation needs pdflatex or tectonic")
        for wrapper in wrappers:
            if pdflatex is not None:
                command = [pdflatex, "-interaction=nonstopmode",
                           "-halt-on-error", wrapper.name]
            else:
                command = [tectonic, "--keep-logs", wrapper.name]
            _run(
                command,
                log=paths.validation / f"latex_{wrapper.stem}.log",
                env=env | {"TEXINPUTS": f"{primary}:"},
                cwd=primary,
            )
    return paths.root


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="pooled_tracking_v1",
                        choices=("pooled_tracking_v1",))
    parser.add_argument("--run-id", default="migration_candidate_20260920")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--no-compile", action="store_true")
    parser.add_argument("--layout-only", action="store_true",
                        help="Use validated numerical caches and update only layout assets.")
    parser.add_argument("--publish", action="store_true",
                        help="Validate and promote the build to output/publication.")
    args = parser.parse_args(argv)
    if args.publish and args.no_compile:
        parser.error("--publish requires compiled wrappers; omit --no-compile")
    root = build(args.profile, args.run_id, args.output_root,
                 compile_tex=not args.no_compile,
                 layout_only=args.layout_only)
    if args.publish:
        from terminal_pareto.publication_release import promote
        promote(args.run_id, args.output_root)
    print(f"Complete pooled publication build: {root}")


if __name__ == "__main__":
    main()
