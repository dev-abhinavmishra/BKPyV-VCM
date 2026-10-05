#!/usr/bin/env python3
"""One-command download of the GSE317012 raw 10x archive and series matrix.

Fetches the public GEO supplementary archive ``GSE317012_RAW.tar`` from NCBI
FTP over HTTPS with resume support, verifies its sha256 against the value
recorded in docs/GSE317012_ANALYSIS.md provenance, and optionally runs the
existing safe extractor/organizer so ``data/raw/GSE317012_RAW/`` is ready for
``scripts/analyze_gse317012.py``. No manual steps.

Usage:
    python scripts/download_gse317012.py                # download + verify + extract
    python scripts/download_gse317012.py --no-extract   # download + verify only
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
RAW_TAR = DATA_DIR / "GSE317012_RAW.tar"
RAW_ROOT = DATA_DIR / "raw" / "GSE317012_RAW"
RESEARCH_DIR = DATA_DIR / "research"

RAW_TAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/suppl/"
    "GSE317012_RAW.tar"
)
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/matrix/"
    "GSE317012_series_matrix.txt.gz"
)

# sha256 recorded in docs/GSE317012_ANALYSIS.md provenance (2026-09-19).
EXPECTED_RAW_TAR_SHA256 = (
    "1EA67E31A190AA05F39A7664721BC0CAF3D219997DFE67C2FC7920C4854327A8"
)
EXPECTED_SERIES_MATRIX_SHA256 = (
    "221b8892bfb3e4277e467ecc8aad13acf4f4a4348d2b3e119ffcc6db1c58b83e"
)

CHUNK = 1 << 20  # 1 MiB


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def download_with_resume(url: str, dest: Path, expected_sha256: str | None) -> dict:
    """Stream url to dest, resuming a partial file if present."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    fetched = datetime.now(timezone.utc).isoformat(timespec="seconds")

    if (
        dest.exists()
        and expected_sha256
        and sha256_file(dest) == expected_sha256.lower()
    ):
        print(f"already complete and verified: {dest}")
        return {
            "url": url,
            "path": str(dest.relative_to(REPO_ROOT)),
            "sha256": expected_sha256,
            "bytes": dest.stat().st_size,
            "fetched_utc": "pre-existing (verified)",
        }

    # Download to a temp sibling so a failed run never leaves a "complete" file.
    tmp = dest.with_suffix(dest.suffix + ".part")
    mode = "ab" if tmp.exists() else "wb"
    have = tmp.stat().st_size if tmp.exists() else 0

    headers = {"User-Agent": "vcm-analysis/1.0"}
    if have:
        headers["Range"] = f"bytes={have}-"
        print(f"resuming {dest.name} at byte {have}")

    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=120) as resp:
        # If the server ignored Range and sent 200, restart from scratch.
        if have and resp.status == 200:
            print("server ignored Range; restarting download")
            mode = "wb"
            have = 0
        total = resp.headers.get("Content-Length")
        total = int(total) + have if total else None
        done = have
        with open(tmp, mode) as fh:
            while True:
                chunk = resp.read(CHUNK)
                if not chunk:
                    break
                fh.write(chunk)
                done += len(chunk)
                if total:
                    pct = 100.0 * done / total
                    print(
                        f"\r{dest.name}: {done / 1e6:,.0f} / {total / 1e6:,.0f} MB"
                        f" ({pct:.1f}%)",
                        end="",
                        flush=True,
                    )
    print()

    tmp.rename(dest)
    digest = sha256_file(dest)
    if expected_sha256 and digest != expected_sha256.lower():
        raise ValueError(
            f"sha256 mismatch for {dest}: got {digest}, expected {expected_sha256}"
        )
    print(f"verified sha256: {digest}")
    return {
        "url": url,
        "path": str(dest.relative_to(REPO_ROOT)),
        "sha256": digest,
        "bytes": dest.stat().st_size,
        "fetched_utc": fetched,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="download + verify only; skip extraction into data/raw/",
    )
    parser.add_argument(
        "--skip-series-matrix",
        action="store_true",
        help="skip the (small) series-matrix download",
    )
    args = parser.parse_args(argv)

    manifest = []

    if not args.skip_series_matrix:
        manifest.append(
            download_with_resume(
                SERIES_MATRIX_URL,
                RESEARCH_DIR / "gse317012_series_matrix.txt.gz",
                EXPECTED_SERIES_MATRIX_SHA256,
            )
        )

    manifest.append(download_with_resume(RAW_TAR_URL, RAW_TAR, EXPECTED_RAW_TAR_SHA256))

    if not args.no_extract:
        if not any(RAW_ROOT.glob("GSM*")):
            print("extracting archive into data/raw/GSE317012_RAW/ ...")
            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "reorganize_gse317012.py"),
                    str(RAW_ROOT),
                    "--extract",
                    str(RAW_TAR),
                    "--move",
                    "--verify",
                ],
                check=True,
            )
        else:
            print(f"extracted tree already present: {RAW_ROOT}")

    print("\nDownload manifest:")
    for entry in manifest:
        print(f"  {entry['path']}  sha256={entry['sha256']}  bytes={entry['bytes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
