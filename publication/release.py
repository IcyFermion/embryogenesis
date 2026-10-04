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


# ---------------------------------------------------------- production

PRODUCTION = provenance.PACKAGE / "output" / "production"
RELEASE_ARCHIVE = provenance.PACKAGE / "output" / "archive" / "releases"
PRODUCTION_MANIFEST = "release_manifest.json"


def family_manifest(key: str) -> str:
    return f"{key}.presentation_manifest.json"


def builds_in(path: Path) -> dict[str, Path]:
    """Family key -> build directory, for one family build or a build set."""
    from publication.registry import BUILD_SET
    path = Path(path)
    if (path / BUILD_SET).exists():
        families = json.loads((path / BUILD_SET).read_text())["families"]
        return {key: path / entry["directory"] for key, entry in families.items()}
    return {json.loads((path / provenance.MANIFEST).read_text())["family"]: path}


def release(build: Path, *, families=None, production: Path = PRODUCTION, archive_root: Path = RELEASE_ARCHIVE,
            when: str | None = None) -> dict:
    """Release family builds into the single flat production folder.

    Each family may add, replace or retire only its registry-owned assets
    (plus its manifest copy and the shared release manifest); every other
    production file must stay byte-identical. Builds must be intact and built
    from the current presentation source.
    """
    from publication.registry import FAMILIES
    available = builds_in(build)
    selected = list(families or available)
    if unknown := set(selected) - set(available):
        raise ValueError(f"Families not in this build: {sorted(unknown)}")
    records = {}
    for key in selected:
        status = provenance.verify_build(available[key])
        if not status["ready"]:
            raise ValueError(f"{key} build is stale ({status['stale_sources']}); rebuild before releasing")
        record = json.loads((available[key] / provenance.MANIFEST).read_text())
        if extra := set(record["files"]) - set(FAMILIES[key].owned_assets):
            raise ValueError(f"{key} build lists assets the family does not own: {sorted(extra)}")
        records[key] = record
    owned = {PRODUCTION_MANIFEST} | {family_manifest(k) for k in selected}
    owned |= {name for k in selected for name in FAMILIES[k].owned_assets}
    production = Path(production)
    production.mkdir(parents=True, exist_ok=True)
    when = when or stamp()

    def mutate(stage: Path, before: dict) -> dict:
        path = stage / PRODUCTION_MANIFEST
        manifest = (json.loads(path.read_text()) if path.exists()
                    else dict(version="publication-release-1", families={}, history=[]))
        retired = set()
        for key in selected:
            record, source = records[key], available[key]
            previous = set(manifest["families"].get(key, {}).get("files", {}))
            for name in previous - set(record["files"]) - {family_manifest(key)}:
                (stage / name).unlink()
                retired.add(name)
            for name in record["files"]:
                shutil.copy2(source / name, stage / name)
            shutil.copy2(source / provenance.MANIFEST, stage / family_manifest(key))
            manifest["families"][key] = dict(
                released_at=when, status=FAMILIES[key].status, build=str(source.resolve()),
                presentation_id=record["presentation_id"], run=record["run"], analysis_id=record["analysis_id"],
                figure_numbers=record["figure_numbers"], input_files=record["input_files"],
                files={name: sha256(stage / name) for name in [*sorted(record["files"]), family_manifest(key)]})
        manifest["history"].append(dict(released_at=when, families=selected,
                                        previous=str(Path(archive_root) / when) if before else None))
        path.write_text(json.dumps(manifest, indent=2) + "\n")
        return dict(check_ownership(before, inventory(stage), owned=owned, retired=retired), families=selected)

    return transaction(production, mutate, lambda d: verify_production(d), archive_root=archive_root, when=when)


def verify_production(production: Path = PRODUCTION, *, check_inputs: bool = True) -> dict:
    """Integrity of released files (and their numerical inputs); stale presentation source is reported."""
    production = Path(production)
    manifest = json.loads((production / PRODUCTION_MANIFEST).read_text())
    expected = {PRODUCTION_MANIFEST}
    stale, summary = set(), {}
    for key, entry in manifest["families"].items():
        check(production, entry["files"])
        expected |= set(entry["files"])
        if check_inputs:
            check(Path(entry["run"]), entry["input_files"])
        recorded = json.loads((production / family_manifest(key)).read_text())["presentation_sources"]
        stale |= set(provenance.stale_files(recorded))
        summary[key] = dict(files=len(entry["files"]), released_at=entry["released_at"])
    if set(inventory(production)) != expected:
        raise ValueError(f"Unexpected or missing production files: "
                         f"{sorted(set(inventory(production)) ^ expected)}")
    return dict(families=summary, stale_presentation_sources=sorted(stale))


def rehearse(build: Path, *, families=None, production: Path = PRODUCTION, keep: Path | None = None) -> dict:
    """Release into a scratch copy of production with the real checks; production is only read."""
    production = Path(production)
    before = inventory(production) if production.exists() else {}
    scratch = Path(keep) if keep else Path(tempfile.mkdtemp(prefix="publication-rehearsal-"))
    copy = scratch / "production"
    if copy.exists():
        raise FileExistsError(copy)
    if production.exists():
        shutil.copytree(production, copy)
    try:
        result = release(build, families=families, production=copy, archive_root=scratch / "archive")
        verify_production(copy)
        if (inventory(production) if production.exists() else {}) != before:
            raise ValueError("Production changed during the rehearsal")
        return dict(result, rehearsal=str(copy) if keep else None, production_unchanged=True)
    finally:
        if not keep:
            shutil.rmtree(scratch)
