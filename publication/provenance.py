"""Presentation identity: dependency hashes, build manifests and verification.

Three identities are kept separate. *Scientific* identity belongs to back-end
caches and is checked by their validators. *Presentation* identity is the hash
set of the publication sources and assets a family's build depends on, plus
the notation snapshot. Dependencies are derived per family (registry entry
points followed through ``publication`` imports), so a change that cannot
affect one family's figures does not mark that family stale.
*Release* identity (production manifests) is handled elsewhere.

Verification distinguishes artifact integrity (files still match the build
manifest) from readiness (live presentation source still matches what built
them). A frozen build can be integrity-verified after source changes without
treating that source as having produced it.
"""

from __future__ import annotations

import ast
import hashlib
import inspect
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
# Shared by every build through ``registry.build`` (style, assembly, manifests).
BUILD_MODULES = ("publication.assembly", "publication.provenance", "publication.style", "publication.notation")
# Holds layout options (titles, colorbar labels) for every family: hashed, not traversed.
REGISTRY = "publication.registry"


def _module_file(name: str) -> Path | None:
    path = ROOT / Path(*name.split("."))
    for candidate in (path.with_suffix(".py"), path / "__init__.py"):
        if candidate.is_file():
            return candidate
    return None


def _publication_imports(tree: ast.AST) -> set[str]:
    """Every ``publication`` module imported anywhere in ``tree`` (including inside functions)."""
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found |= {alias.name for alias in node.names if alias.name.startswith("publication")}
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and (node.module or "").startswith("publication"):
            found.add(node.module)
            found |= {f"{node.module}.{alias.name}" for alias in node.names
                      if _module_file(f"{node.module}.{alias.name}") is not None}
    return found


def module_closure(entries, *, files_only=()) -> set[Path]:
    """Source files of ``entries`` and every ``publication`` module they import, transitively.

    ``files_only`` modules are hashed but not traversed (their relevant
    imports were already resolved by the caller).
    """
    pending, seen, files = list(entries), set(), {_module_file(name) for name in files_only}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        path = _module_file(name)
        if path is None:
            continue
        files.add(path)
        if name != REGISTRY and name not in files_only:
            pending.extend(_publication_imports(ast.parse(path.read_text())))
    return files


def _function_imports(module: str, function: str) -> set[str]:
    """Top-level imports of ``module`` plus those inside ``function`` and same-module functions it calls.

    Adapter modules serve several families; this keeps one family's adapter
    from inheriting the lazy imports of another family's loader.
    """
    tree = ast.parse(_module_file(module).read_text())
    functions = {node.name: node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    found = _publication_imports(ast.Module(body=[node for node in tree.body if node.__class__.__name__ not in
                                                  ("FunctionDef", "AsyncFunctionDef", "ClassDef")], type_ignores=[]))
    pending, seen = [function], set()
    while pending:
        name = pending.pop()
        if name in seen or name not in functions:
            continue
        seen.add(name)
        found |= _publication_imports(functions[name])
        pending.extend(node.func.id for node in ast.walk(functions[name])
                       if isinstance(node, ast.Call) and isinstance(node.func, ast.Name))
    return found


def family_entry_modules(family) -> set[str]:
    """Registry entry points of one family: adapter, figure/caption functions, builder imports."""
    module, function = family.adapter.split(":")
    entries = {*_function_imports(module, function), *BUILD_MODULES, REGISTRY}
    for spec in family.figures:
        entries.add(getattr(spec.render, "func", spec.render).__module__)
        entries.add(spec.caption.__module__)
    entries |= _publication_imports(ast.parse(inspect.getsource(family.builder)))
    return entries


def _hashes(files) -> dict[str, str]:
    return {str(path.relative_to(ROOT)): sha256(path) for path in sorted(files)
            if path.is_file() and path.relative_to(PACKAGE).as_posix() not in NON_RENDERING
            and not {"tests", "output", "__pycache__"} & set(path.relative_to(PACKAGE).parts)}


def presentation_sources(family: str | None = None) -> dict[str, str]:
    """Sources and assets a family's build depends on (every candidate when ``family`` is None)."""
    if family is None:
        return _hashes([*PACKAGE.rglob("*.py"), *(PACKAGE / "assets").rglob("*")])
    from publication.registry import FAMILIES
    spec = FAMILIES[family]
    adapter = spec.adapter.split(":")[0]
    entries = family_entry_modules(spec) - {adapter}
    return _hashes([*module_closure(entries, files_only=(adapter,)), *(PACKAGE / asset for asset in spec.assets)])


def changed_sources(recorded: dict[str, str], live: dict[str, str]) -> list[str]:
    """Dependencies added, removed or changed since ``recorded``."""
    return sorted(name for name in set(live) | set(recorded) if live.get(name) != recorded.get(name))


def presentation_id(sources: dict[str, str], notation_snapshot: dict) -> str:
    payload = json.dumps(dict(sources=sources, notation=notation_snapshot), sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def write_manifest(out: Path, *, family: str, inputs: dict, artifacts: list[Path], extra: dict | None = None) -> dict:
    """Record a family build; ``inputs`` maps run-relative cache files to hashes."""
    sources, snapshot = presentation_sources(family), notation.snapshot()
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
    live = presentation_sources(record["family"])
    stale = changed_sources(record["presentation_sources"], live)
    current = presentation_id(live, notation.snapshot()) == record["presentation_id"]
    return dict(family=record["family"], files=len(record["files"]),
                integrity=True, ready=current and not stale, stale_sources=stale)
