"""Build the reserve notes: ``python -m publication.notes --output-dir DIR``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from publication.notes import apoptosis

NOTES = {"apoptosis_abplp": apoptosis.build}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"{args.output_dir} is not empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    record = {name: build(args.output_dir) for name, build in NOTES.items()}
    (args.output_dir / "build.json").write_text(json.dumps(record, indent=2, default=str) + "\n")
    print(json.dumps(record, indent=2, default=str))


if __name__ == "__main__":
    main()
