"""Cache-only presentation CLI.

    python -m publication list
    python -m publication build --family terminal-cross-species --output-dir PATH
    python -m publication build --all --output-dir PATH
    python -m publication verify --build PATH

Builds read validated caches only and write to a new directory under
publication/output/candidates/. Since 2026-11-13 the families only regenerate
backlog figures (publication/BACKLOG.md); the manuscript's figures are the set
built by publication.v2, and family releases into a production folder are retired.
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
    verify.add_argument("--build", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "list":
        print(registry.dumps(registry.describe()))
    elif args.command == "build" and args.all:
        manifest = registry.build_all(args.output_dir)
        print(f"Built {len(manifest['families'])} families in {args.output_dir}")
    elif args.command == "build":
        record = registry.build(args.family, args.output_dir, run_id=args.run_id)
        print(f"Built {record['family']} ({len(record['files'])} files) in {args.output_dir}")
    else:
        print(registry.dumps(registry.verify_any(args.build)))


if __name__ == "__main__":
    main()
