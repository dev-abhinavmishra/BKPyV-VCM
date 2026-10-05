#!/usr/bin/env python3
"""Generate the reviewable BKPyV figures and evidence manifest."""
from pathlib import Path
import sys
import argparse

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import matplotlib

matplotlib.use("Agg")  # headless-friendly for Windows/CI environments

from vcm.viz.bkpyv_review import build_review_bundle


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="outputs/bkpyv_review")
    parser.add_argument("--days", type=float, default=60.0)
    args = parser.parse_args()
    paths = build_review_bundle(output_dir=args.output_dir, days=args.days)
    print("BKPyV review bundle generated:")
    for name, path in paths.items():
        print(f"  {name}: {path}")
