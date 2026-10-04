"""Shared release mechanics: inventories, hash-checked archives and staged transactions.

A release stages a complete copy of a production directory, applies one
family's changes to the copy, verifies it, archives the live directory with
hashes, installs the stage and verifies again; any failure restores the
original directory byte-for-byte. Ownership checks guarantee that a family
release changes only that family's assets (plus explicitly coordinated
manifests). Family-specific scientific validation and manifest semantics stay
with the family release functions below.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile
from typing import Callable

from publication import provenance
from publication.provenance import sha256


def stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def inventory(directory: Path) -> dict[str, str]:
    directory = Path(directory)
    if any(path.is_symlink() for path in directory.rglob("*")):
        raise ValueError(f"Release directories must not contain symlinks: {directory}")
    return {str(p.relative_to(directory)): sha256(p) for p in sorted(directory.rglob("*")) if p.is_file()}


def check(directory: Path, hashes: dict[str, str]) -> None:
    changed = [name for name, digest in hashes.items()
               if not (Path(directory) / name).is_file() or sha256(Path(directory) / name) != digest]
    if changed:
        raise ValueError(f"Hash mismatch in {directory}: {changed}")


def check_ownership(before: dict, after: dict, *, owned: set[str], retired: set[str]) -> dict:
    """Every added, changed or removed file must be owned by the release or explicitly retired."""
    added = set(after) - set(before)
    removed = set(before) - set(after)
    changed = {name for name in set(before) & set(after) if before[name] != after[name]}
    foreign = sorted((added | changed) - owned) + sorted(removed - retired - owned)
    if foreign:
        raise ValueError(f"Release would alter assets it does not own: {foreign}")
    return dict(added=sorted(added), changed=sorted(changed), removed=sorted(removed))


def transaction(target: Path, mutate: Callable[[Path, dict], dict], verify: Callable[[Path], None], *,
                archive_root: Path, when: str | None = None) -> dict:
    """Stage, verify, archive and install; restore the original directory on any failure."""
    target = Path(target)
    if target.is_symlink() or not target.is_dir():
        raise ValueError("Production must be an existing real directory")
    before = inventory(target)
    when = when or stamp()
    archive = Path(archive_root) / when
    if archive.exists():
        raise FileExistsError(archive)
    stage = Path(tempfile.mkdtemp(prefix=".release-stage-", dir=target.parent))
    moved = installed = False
    try:
        shutil.copytree(target, stage, dirs_exist_ok=True)
        check(stage, before)
        summary = mutate(stage, before)
        verify(stage)
        if inventory(target) != before:
            raise ValueError("Production changed during staging; leaving it untouched")
        archive.mkdir(parents=True)
        target.rename(archive / "publication")
        moved = True
        (archive / "manifest.json").write_text(json.dumps(dict(files=before), indent=2) + "\n")
        check(archive / "publication", before)
        stage.rename(target)
        installed = True
        verify(target)
    except Exception:
        if installed:
            shutil.rmtree(target)
        if moved:
            (archive / "publication").rename(target)
        if archive.exists():
            shutil.rmtree(archive)
        raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return dict(summary or {}, archive=str(archive), published_at=when)


# ---------------------------------------------------------------- full tree

FULL_TREE_RELEASE_MANIFEST = "release_manifest.json"
FULL_TREE_PRESENTATION_COPY = "full_tree_pooled_presentation_manifest.json"
CROSS_SPECIES_RELEASE_MANIFEST = "cross_species_release_manifest.json"
# Historical single-embryo Figure 7/supplement assets superseded by the pooled
# family without a same-named replacement, plus the phylogenetic-reference
# methods notes the author chose to drop from production for now (2026-10-03);
# every retired file stays recoverable in the hash-checked release archive.
FULL_TREE_POOLED_RETIRES = (
    "fig7B_ce_full_tree_layerwise_aggregate.pdf", "fig7B_ce_full_tree_layerwise_aggregate.png",
    "table_ce_full_tree_heuristics.pdf", "ce_full_tree_heuristic_inventory.csv",
    "ce_full_tree_heuristic_inventory_rows.tex",
    "brownian_covariance_mle_derivation.pdf", "parametric_brownian_bootstrap.pdf",
    "separate_clock_reference.pdf",
)


def _default_full_tree_verifiers():
    from full_tree_pareto import cross_species_publication as cross_species
    from full_tree_pareto import publication_build as pooled
    return pooled.verify_release, cross_species.verify_release


def release_full_tree_pooled(build: Path, production: Path, *, archive_root: Path | None = None,
                             verifiers=None, when: str | None = None) -> dict:
    """Mixed full-tree release: pooled Figure 7/S4 in, Figures 10/11 preserved with refreshed records.

    ``build`` is a verified ``full-tree-pooled`` presentation build. The pooled
    ``release_manifest.json`` keeps the legacy format read by
    ``publication_build --verify``; the cross-species release manifest's
    preserved-file record is rewritten to the new inventory (previous record
    kept under ``preserved_files_history``) so both verifiers agree.
    """
    build, production = Path(build), Path(production)
    status = provenance.verify_build(build)
    record = json.loads((build / provenance.MANIFEST).read_text())
    if record["family"] != "full-tree-pooled":
        raise ValueError(f"Not a full-tree-pooled build: {record['family']}")
    verify_pooled, verify_cross_species = verifiers or _default_full_tree_verifiers()
    when = when or stamp()
    archive_root = Path(archive_root) if archive_root else production.parent / "legacy/releases"
    from publication.registry import FAMILIES
    family_assets = set(FAMILIES["full-tree-pooled"].owned_assets)
    build_files = sorted(record["files"])
    if set(build_files) - family_assets:
        raise ValueError(f"Build lists assets the full-tree-pooled family does not own: "
                         f"{sorted(set(build_files) - family_assets)}")
    owned = family_assets | {FULL_TREE_RELEASE_MANIFEST, FULL_TREE_PRESENTATION_COPY, CROSS_SPECIES_RELEASE_MANIFEST}
    retired = set(FULL_TREE_POOLED_RETIRES)
    run = Path(record["run"])

    def mutate(stage: Path, before: dict) -> dict:
        removed = {name: before[name] for name in FULL_TREE_POOLED_RETIRES if name in before}
        for name in removed:
            (stage / name).unlink()
        for name in build_files:
            shutil.copy2(build / name, stage / name)
        shutil.copy2(build / provenance.MANIFEST, stage / FULL_TREE_PRESENTATION_COPY)
        files = {name: sha256(stage / name) for name in [*build_files, FULL_TREE_PRESENTATION_COPY]}
        pooled_record = dict(
            version="full-tree-pooled-release-2", profile="pooled_full_tree_v1", run=str(run), published_at=when,
            presentation_id=record["presentation_id"], figure_numbers=record["figure_numbers"], files=files,
            analysis_files={name.removeprefix("analysis/"): digest for name, digest in record["input_files"].items()
                            if name.startswith("analysis/")},
            retired_files=removed, previous_publication=str(archive_root / when))
        (stage / FULL_TREE_RELEASE_MANIFEST).write_text(json.dumps(pooled_record, indent=2) + "\n")
        cross = stage / CROSS_SPECIES_RELEASE_MANIFEST
        if cross.exists():
            cross_record = json.loads(cross.read_text())
            current = inventory(stage)
            preserved = {name: digest for name, digest in current.items()
                         if name not in cross_record["files"] and name != CROSS_SPECIES_RELEASE_MANIFEST}
            history = cross_record.get("preserved_files_history", [])
            history.append(dict(updated_at=when, by="full-tree-pooled release", previous=cross_record["preserved_files"]))
            cross_record.update(preserved_files=preserved, preserved_files_history=history)
            cross.write_text(json.dumps(cross_record, indent=2) + "\n")
        return check_ownership(before, inventory(stage), owned=owned, retired=retired)

    def verify(directory: Path) -> None:
        verify_pooled(directory)
        check(directory, {name: digest for name, digest in record["files"].items()})
        if (directory / CROSS_SPECIES_RELEASE_MANIFEST).exists():
            verify_cross_species(directory, check_archive=False)

    result = transaction(production, mutate, verify, archive_root=archive_root, when=when)
    return dict(result, build=str(build), stale_presentation_sources=status["stale_sources"])


def rehearse_full_tree_pooled(build: Path, production: Path, *, keep: Path | None = None) -> dict:
    """Run the mixed release against a scratch copy of ``production`` with the real verifiers.

    Production itself is only read. With ``keep``, the rehearsed directory is
    left there for inspection; otherwise it is removed.
    """
    production = Path(production)
    before = inventory(production)
    scratch = Path(keep) if keep else Path(tempfile.mkdtemp(prefix="publication-rehearsal-"))
    copy = scratch / "output/publication"
    if copy.exists():
        raise FileExistsError(copy)
    shutil.copytree(production, copy)
    try:
        result = release_full_tree_pooled(build, copy, archive_root=scratch / "output/legacy/releases")
        if (copy / CROSS_SPECIES_RELEASE_MANIFEST).exists():
            from full_tree_pareto import cross_species_publication
            cross_species_publication.verify_release(copy, check_archive=True)
        if inventory(production) != before:
            raise ValueError("Production changed during the rehearsal")
        return dict(result, rehearsal=str(copy) if keep else None, production_unchanged=True)
    finally:
        if not keep:
            shutil.rmtree(scratch)
