"""Cache-only presentation CLI.

    python -m publication list
    python -m publication build --family terminal-cross-species --output-dir PATH
    python -m publication build --all --output-dir PATH
    python -m publication verify --build PATH | --production
    python -m publication release --build PATH [--family KEY ...] --rehearse | --apply

Builds read validated caches only and write to a new directory under
publication/output/candidates/. Releases copy family builds into the single
production folder publication/output/production/ (archiving the previous state
under publication/output/archive/releases/); --rehearse uses a scratch copy.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from publication import registry  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m publication", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="Show publication families and their figures")
    build = commands.add_parser("build", help="Render one family (or --all) from validated caches into a new directory")
    target = build.add_mutually_exclusive_group(required=True)
    target.add_argument("--family", choices=sorted(registry.FAMILIES))
    target.add_argument("--all", action="store_true")
    build.add_argument("--output-dir", required=True, type=Path)
    build.add_argument("--run-id", help="Back-end run identity (defaults to the family's published run)")
    verify = commands.add_parser("verify", help="Check a build's artifacts, inputs and source readiness")
    verify.add_argument("--build", type=Path)
    verify.add_argument("--production", action="store_true", help="Verify publication/output/production instead")
    rel = commands.add_parser("release", help="Rehearse, or apply, a release of family builds into production")
    rel.add_argument("--build", required=True, type=Path, help="A family build or a build set (build --all)")
    rel.add_argument("--family", action="append", choices=sorted(registry.FAMILIES),
                     help="Release only these families from a build set (repeatable)")
    mode = rel.add_mutually_exclusive_group(required=True)
    mode.add_argument("--rehearse", action="store_true", help="Release into a scratch copy of production only")
    mode.add_argument("--apply", action="store_true", help="Archive and update production (staged, rolled back on failure)")
    rel.add_argument("--keep", type=Path, help="With --rehearse: keep the rehearsed production here")
    args = parser.parse_args(argv)
    if args.command == "list":
        print(registry.dumps(registry.describe()))
    elif args.command == "build" and args.all:
        manifest = registry.build_all(args.output_dir)
        print(f"Built {len(manifest['families'])} families in {args.output_dir}")
    elif args.command == "build":
        record = registry.build(args.family, args.output_dir, run_id=args.run_id)
        print(f"Built {record['family']} ({len(record['files'])} files) in {args.output_dir}")
    elif args.command == "verify":
        if args.production == bool(args.build):
            parser.error("verify needs exactly one of --build or --production")
        from publication import release
        print(registry.dumps(release.verify_production() if args.production else registry.verify_any(args.build)))
    else:
        from publication import release
        if args.apply:
            print(registry.dumps(release.release(args.build, families=args.family)))
        else:
            print(registry.dumps(release.rehearse(args.build, families=args.family, keep=args.keep)))


if __name__ == "__main__":
    main()
