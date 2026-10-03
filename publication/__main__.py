"""Cache-only presentation CLI.

    python -m publication list
    python -m publication build --family terminal-cross-species --output-dir PATH
    python -m publication verify --build PATH
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from publication import provenance, registry  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m publication", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="Show publication families and their figures")
    build = commands.add_parser("build", help="Render one family from validated caches into a new directory")
    build.add_argument("--family", required=True, choices=sorted(registry.FAMILIES))
    build.add_argument("--output-dir", required=True, type=Path)
    build.add_argument("--run-id", help="Back-end run identity (defaults to the family's published run)")
    verify = commands.add_parser("verify", help="Check a build's artifacts, inputs and source readiness")
    verify.add_argument("--build", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.command == "list":
        print(registry.dumps(registry.describe()))
    elif args.command == "build":
        record = registry.build(args.family, args.output_dir, run_id=args.run_id)
        print(f"Built {record['family']} ({len(record['files'])} files) in {args.output_dir}")
    else:
        print(registry.dumps(provenance.verify_build(args.build)))


if __name__ == "__main__":
    main()
