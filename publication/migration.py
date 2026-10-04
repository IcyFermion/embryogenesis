"""One-off move of figure material out of the back ends (author decision 2026-10-04).

    python -m publication.migration --archive        # copy + hash-verify into publication/output/archive/
    python -m publication.migration --verify         # archive still matches its manifest
    python -m publication.migration --purge --confirm  # after the PR: re-verify, then delete the originals

Only figure-production material is listed: rendered panels, wrappers, layout
snapshots and render logs written into run folders, plus the two old
production folders. Numerical results (``analysis/``, ``checkpoints/``, replay
and migration-validation reports) stay in ``terminal_pareto/`` and
``full_tree_pareto/``. The S3 projection table and provenance were copied into
the pooled run's ``analysis/s3_cross_geometry/`` before archiving.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from publication.provenance import ROOT, sha256

ARCHIVE = ROOT / "publication/output/archive/migrated_20261004"
MANIFEST = "manifest.json"
TERMINAL_RUNS = "terminal_pareto/output/runs"
FULL_TREE_RUNS = "full_tree_pareto/output/runs"
POOLED_VALIDATION = f"{TERMINAL_RUNS}/pooled_tracking_v1/migration_candidate_20260920/validation"

FIGURE_PATHS = (
    # Old production folders, superseded by publication/output/production/.
    "terminal_pareto/output/publication",
    "full_tree_pareto/output/publication",
    # Working figures written beside numerical runs.
    f"{TERMINAL_RUNS}/pooled_tracking_v1/migration_candidate_20260920/publication",
    f"{TERMINAL_RUNS}/embryo1_legacy/migration_comparison_20260920/publication",
    f"{TERMINAL_RUNS}/embryo1_matched/migration_comparison_20260920/publication",
    f"{TERMINAL_RUNS}/cross_species_terminal_v1/pooled_comparison_20260930/figures",
    f"{TERMINAL_RUNS}/cross_species_terminal_v1/pooled_comparison_20260930/publication",
    f"{TERMINAL_RUNS}/cross_species_terminal_v1/pooled_comparison_20260930/layout_history",
    f"{FULL_TREE_RUNS}/pooled_full_tree_v1/terminal_clamped_20260927/publication",
    f"{FULL_TREE_RUNS}/pooled_full_tree_v1/terminal_clamped_20260927/layout_history",
    f"{FULL_TREE_RUNS}/cross_species_layerwise_v1/terminal_anchored_20260930/figures",
    f"{FULL_TREE_RUNS}/cross_species_layerwise_v1/terminal_anchored_20260930/publication",
    f"{FULL_TREE_RUNS}/cross_species_layerwise_v1/terminal_anchored_20260930/layout_history",
    # Figure-only logs and layout snapshots in the pooled run's validation folder.
    f"{POOLED_VALIDATION}/fig5_before_redesign_20260923",
    f"{POOLED_VALIDATION}/s3_legend_before_20260923",
    f"{POOLED_VALIDATION}/fig5_redesign_preview.png",
    f"{POOLED_VALIDATION}/matplotlib",
    f"{POOLED_VALIDATION}/fig1_endpoint_amendment.log",
    f"{POOLED_VALIDATION}/fig4.log",
    f"{POOLED_VALIDATION}/fig5_figs1.log",
    f"{POOLED_VALIDATION}/latex_fig1_ce_endpoint_normalization_amendment.log",
    f"{POOLED_VALIDATION}/latex_fig1_endpoint_normalization_amendment.log",
    f"{POOLED_VALIDATION}/latex_fig2_ce_terminal_pareto_main.log",
    f"{POOLED_VALIDATION}/latex_fig3_ce_terminal_pareto_supporting.log",
    f"{POOLED_VALIDATION}/latex_fig4_ce_subtree_map.log",
    f"{POOLED_VALIDATION}/latex_fig5_ce_canonical_summary.log",
    f"{POOLED_VALIDATION}/latex_fig6_ce_cell_types.log",
    f"{POOLED_VALIDATION}/latex_figS1_ce_canonical_summary_cousin_r.log",
    f"{POOLED_VALIDATION}/latex_figS2_ce_cell_type_cost_gain.log",
    f"{POOLED_VALIDATION}/latex_figS3_ce_tracking_robustness.log",
)


def _files(path: Path) -> list[Path]:
    if path.is_symlink():
        raise ValueError(f"Refusing to migrate a symlink: {path}")
    return [path] if path.is_file() else sorted(p for p in path.rglob("*") if p.is_file())


def present() -> list[str]:
    return [name for name in FIGURE_PATHS if (ROOT / name).exists()]


def archive() -> dict:
    """Copy every listed path into the archive and verify each file's hash."""
    if (ARCHIVE / MANIFEST).exists():
        raise FileExistsError(f"{ARCHIVE} already holds a migration archive")
    hashes = {}
    for name in present():
        for source in _files(ROOT / name):
            relative = str(source.relative_to(ROOT))
            target = ARCHIVE / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            hashes[relative] = sha256(source)
            if sha256(target) != hashes[relative]:
                raise ValueError(f"Archive copy mismatch: {relative}")
    record = dict(version="figure-migration-1", paths=present(), files=hashes)
    (ARCHIVE / MANIFEST).write_text(json.dumps(record, indent=2) + "\n")
    return dict(paths=len(record["paths"]), files=len(hashes))


def verify() -> dict:
    record = json.loads((ARCHIVE / MANIFEST).read_text())
    bad = [name for name, digest in record["files"].items()
           if not (ARCHIVE / name).is_file() or sha256(ARCHIVE / name) != digest]
    if bad:
        raise ValueError(f"Migration archive mismatch: {bad[:10]}")
    return dict(paths=len(record["paths"]), files=len(record["files"]))


def purge() -> dict:
    """Delete migrated originals once the archive verifies and still covers them exactly."""
    record = json.loads((ARCHIVE / MANIFEST).read_text())
    verify()
    removed = []
    for name in record["paths"]:
        path = ROOT / name
        if not path.exists():
            continue
        current = {str(p.relative_to(ROOT)): sha256(p) for p in _files(path)}
        unarchived = [f for f, digest in current.items() if record["files"].get(f) != digest]
        if unarchived:
            raise ValueError(f"{name} changed since archiving; re-archive first: {unarchived[:5]}")
        shutil.rmtree(path) if path.is_dir() else path.unlink()
        removed.append(name)
    return dict(removed=removed)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--archive", action="store_true")
    action.add_argument("--verify", action="store_true")
    action.add_argument("--purge", action="store_true")
    parser.add_argument("--confirm", action="store_true", help="Required with --purge")
    args = parser.parse_args(argv)
    if args.purge and not args.confirm:
        parser.error("--purge deletes the originals; add --confirm")
    result = archive() if args.archive else verify() if args.verify else purge()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
