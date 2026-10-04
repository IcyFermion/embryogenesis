"""Compare a candidate build with reference (e.g. production) assets.

    python -m publication.parity --candidate BUILD_DIR --reference PRODUCTION_DIR [--names a.pdf b.png ...]

TeX is compared byte-for-byte; PNG panels and rasterized PDF pages pixel by
pixel. Reports differing-pixel counts so near-identical renders can be
inspected rather than judged by PDF bytes (whose metadata changes).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image


def _pixels(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("RGB")).astype(np.int16)


def _raster(pdf: Path, folder: Path, dpi: int) -> list[Path]:
    prefix = folder / pdf.stem
    subprocess.run(["pdftoppm", "-r", str(dpi), "-png", str(pdf), str(prefix)], check=True)
    return sorted(folder.glob(f"{pdf.stem}*.png"))


def _image_diff(a: np.ndarray, b: np.ndarray) -> dict:
    if a.shape != b.shape:
        return dict(identical=False, shape=[list(a.shape), list(b.shape)])
    changed = np.any(a != b, axis=2)
    return dict(identical=not changed.any(), differing_pixels=int(changed.sum()),
                differing_fraction=float(changed.mean()), max_channel_diff=int(np.abs(a - b).max()))


def compare(candidate: Path, reference: Path, names=None, *, dpi: int = 100) -> dict:
    candidate, reference = Path(candidate), Path(reference)
    names = names or sorted(p.name for p in candidate.iterdir() if p.suffix in (".tex", ".png", ".pdf"))
    report = {}
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "a").mkdir()
        (Path(tmp) / "b").mkdir()
        for name in names:
            a, b = candidate / name, reference / name
            if not b.exists():
                report[name] = dict(missing_reference=True)
            elif a.suffix == ".tex":
                report[name] = dict(identical=a.read_bytes() == b.read_bytes())
            elif a.suffix == ".png":
                report[name] = _image_diff(_pixels(a), _pixels(b))
            else:
                pages_a, pages_b = _raster(a, Path(tmp) / "a", dpi), _raster(b, Path(tmp) / "b", dpi)
                if len(pages_a) != len(pages_b):
                    report[name] = dict(identical=False, pages=[len(pages_a), len(pages_b)])
                else:
                    diffs = [_image_diff(_pixels(x), _pixels(y)) for x, y in zip(pages_a, pages_b)]
                    report[name] = diffs[0] if len(diffs) == 1 else dict(pages=diffs)
                for page in pages_a + pages_b:
                    page.unlink()
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--names", nargs="*")
    parser.add_argument("--dpi", type=int, default=100)
    args = parser.parse_args(argv)
    print(json.dumps(compare(args.candidate, args.reference, args.names, dpi=args.dpi), indent=2))


if __name__ == "__main__":
    main()
