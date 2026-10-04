"""Presentation identity: dependency hashes, build manifests and verification.

Three identities are kept separate. *Scientific* identity belongs to back-end
caches and is checked by their validators. *Presentation* identity is the hash
set of every publication source plus the notation snapshot used for a build.
*Release* identity (production manifests) is handled elsewhere.

Verification distinguishes artifact integrity (files still match the build
manifest) from readiness (live presentation source still matches what built
them). A frozen build can be integrity-verified after source changes without
treating that source as having produced it.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from publication import notation

PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE.parent
MANIFEST = "presentation_manifest.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stale_files(hashes: dict[str, str], base: Path = ROOT) -> list[str]:
    """Recorded files whose live bytes differ from (or no longer exist at) their hash.

    Use for *readiness* questions about source code. Frozen artifacts must be
    checked strictly instead; stale source never invalidates an existing build.
    """
    return sorted(name for name, digest in hashes.items()
                  if not (Path(base) / name).is_file() or sha256(Path(base) / name) != digest)


# Maintenance tools that never influence rendered assets.
NON_RENDERING = {"migration.py", "parity.py"}


def presentation_sources() -> dict[str, str]:
    """Every publication source and asset that can affect a build; any change invalidates rendered assets."""
    files = [*PACKAGE.rglob("*.py"), *(PACKAGE / "assets").rglob("*")]
    return {str(path.relative_to(ROOT)): sha256(path) for path in sorted(files)
            if path.is_file() and path.relative_to(PACKAGE).as_posix() not in NON_RENDERING
            and not {"tests", "output", "__pycache__"} & set(path.relative_to(PACKAGE).parts)}


def presentation_id(sources: dict[str, str], notation_snapshot: dict) -> str:
    payload = json.dumps(dict(sources=sources, notation=notation_snapshot), sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def write_manifest(out: Path, *, family: str, inputs: dict, artifacts: list[Path], extra: dict | None = None) -> dict:
    """Record a family build; ``inputs`` maps run-relative cache files to hashes."""
    sources, snapshot = presentation_sources(), notation.snapshot()
    record = dict(
        version="publication-build-1", family=family,
        presentation_id=presentation_id(sources, snapshot),
        presentation_sources=sources, notation=snapshot,
        **inputs,
        files={path.name: sha256(path) for path in sorted(artifacts)},
        **(extra or {}))
    (Path(out) / MANIFEST).write_text(json.dumps(record, indent=2) + "\n")
    return record


def verify_build(out: Path, *, check_inputs: bool = True) -> dict:
    """Return integrity/readiness status for one family build directory.

    Raises on artifact or numerical-input mismatches (integrity). Changed live
    presentation source is reported as ``stale_sources`` (not ready to rebuild
    identically) rather than silently blessed or treated as corruption.
    """
    out = Path(out)
    record = json.loads((out / MANIFEST).read_text())
    changed = [name for name, digest in record["files"].items()
               if not (out / name).is_file() or sha256(out / name) != digest]
    if changed:
        raise ValueError(f"Changed or missing build artifacts in {out}: {changed}")
    if check_inputs:
        run = Path(record["run"])
        changed = [name for name, digest in record["input_files"].items()
                   if not (run / name).is_file() or sha256(run / name) != digest]
        if changed:
            raise ValueError(f"Numerical inputs changed since build: {changed}")
    live = presentation_sources()
    stale = sorted(name for name in set(live) | set(record["presentation_sources"])
                   if live.get(name) != record["presentation_sources"].get(name))
    current = presentation_id(live, notation.snapshot()) == record["presentation_id"]
    return dict(family=record["family"], files=len(record["files"]),
                integrity=True, ready=current and not stale, stale_sources=stale)
