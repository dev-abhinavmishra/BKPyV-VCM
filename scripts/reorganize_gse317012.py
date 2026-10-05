#!/usr/bin/env python3
"""Group flat GSE317012 10x files into GSM-named directories.

The utility is deliberately file-format agnostic: it only copies the three
standard 10x artifacts and never imports scanpy or rewrites their contents.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import tarfile
from collections.abc import Iterable, Mapping
from pathlib import Path, PurePosixPath

TENX_FILE = re.compile(
    r"^(?P<sample>GSM\d+)(?:_[^.]+?)?_"
    r"(?P<prefix>raw_)?"
    r"(?P<kind>matrix|barcodes|features)\.(?P<extension>mtx|tsv)"
    r"(?P<compressed>\.gz)?$"
)
REQUIRED_KINDS = ("matrix", "barcodes", "features")

MISSING_INPUT_NOTE = """\
ERROR: required input is missing: {path}

The raw archive (data/GSE317012_RAW.tar) and the extracted tree
(data/raw/GSE317012_RAW/) were deleted on 2026-09-19 under user-authorized
storage reclamation. To restore: re-download GSE317012_RAW.tar from NCBI
GEO FTP
(https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/suppl/GSE317012_RAW.tar),
verify sha256 == 1EA67E31A190AA05F39A7664721BC0CAF3D219997DFE67C2FC7920C4854327A8,
then re-run with --extract <tar> --move --verify.
"""


def _kind_key(match: re.Match) -> str:
    """Full kind key including the raw_ prefix so filtered and raw coexist."""
    return f"{match.group('prefix') or ''}{match.group('kind')}"


def extract_tar(tar_path: Path, dest_dir: Path) -> list:
    """Extract a tar archive safely.

    Rejects traversal/absolute member paths, skips members whose destination
    already exists with identical content (idempotent), and refuses to
    overwrite destination files holding different content.
    Returns the list of paths written this run.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_root = dest_dir.resolve()
    written = []
    with tarfile.open(tar_path) as tf:
        for member in tf.getmembers():
            name = member.name
            rel = PurePosixPath(name)
            target = (dest_root / rel).resolve()
            if (
                rel.is_absolute()
                or ".." in rel.parts
                or re.match(r"^[A-Za-z]:", rel.parts[0] if rel.parts else "")
                or not target.is_relative_to(dest_root)
            ):
                raise ValueError(f"unsafe tar member path (traversal/absolute): {name}")
            if not member.isfile():
                continue
            payload = tf.extractfile(member).read()
            if target.exists():
                existing = target.read_bytes()
                if existing == payload:
                    continue
                raise ValueError(
                    f"refusing to overwrite existing file with different content: {target}"
                )
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
            written.append(target)
    return written


def discover_samples(source_dir: Path) -> dict[str, dict[str, Path]]:
    """Return ``{GSM id: {10x kind: source path}}`` for a flat directory."""
    if not source_dir.is_dir():
        raise FileNotFoundError(f"Source directory does not exist: {source_dir}")

    samples: dict[str, dict[str, Path]] = {}
    for path in sorted(source_dir.iterdir()):
        if not path.is_file():
            continue
        match = TENX_FILE.match(path.name)
        if not match:
            continue
        sample = match.group("sample")
        kind = _kind_key(match)
        previous = samples.setdefault(sample, {}).get(kind)
        if previous is not None:
            raise ValueError(f"Duplicate {kind} file for {sample}: {previous.name} and {path.name}")
        samples[sample][kind] = path
    return samples


def _destination_path(output_dir: Path, sample: str, source: Path) -> Path:
    match = TENX_FILE.match(source.name)
    if match is None:
        raise ValueError(f"Not a recognized 10x file: {source.name}")
    suffix = f".{match.group('extension')}{match.group('compressed') or ''}"
    return output_dir / sample / f"{_kind_key(match)}{suffix}"


def verify_samples(
    samples: Mapping[str, Mapping[str, Path]],
    output_dir: Path,
    dry_run: bool = False,
) -> None:
    """Verify that every sample has all required 10x artifacts."""
    incomplete = {
        sample: sorted(set(REQUIRED_KINDS) - set(files))
        for sample, files in samples.items()
        if set(REQUIRED_KINDS) - set(files)
    }
    if incomplete:
        details = ", ".join(f"{sample}: {', '.join(kinds)}" for sample, kinds in incomplete.items())
        raise ValueError(f"Incomplete 10x sample groups: {details}")

    if dry_run:
        return

    missing = []
    for sample, files in samples.items():
        for source in files.values():
            destination = _destination_path(output_dir, sample, source)
            if not destination.is_file():
                missing.append(str(destination))
    if missing:
        raise FileNotFoundError("Missing reorganized files: " + ", ".join(missing))


def reorganize(
    source_dir: Path,
    output_dir: Path | None = None,
    *,
    dry_run: bool = False,
    verify: bool = False,
    move: bool = False,
) -> dict[str, dict[str, Path]]:
    """Copy (or move) flat 10x files into ``output_dir/GSM.../`` directories.

    Every discovered kind (filtered and ``raw_`` prefixed) is organized.
    Existing destinations are skipped when identical, and the run refuses to
    overwrite a destination that holds different content.
    """
    source_dir = Path(source_dir)
    output_dir = Path(output_dir) if output_dir is not None else source_dir
    samples = discover_samples(source_dir)

    for sample, files in samples.items():
        for kind, source in sorted(files.items()):
            destination = _destination_path(output_dir, sample, source)
            verb = "move" if move else "copy"
            print(
                f"{'Would ' + verb if dry_run else verb.capitalize() + 'ing'} "
                f"{source.name} -> {destination}"
            )
            if dry_run:
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                if destination.read_bytes() == source.read_bytes():
                    if move and source != destination:
                        source.unlink()
                    continue
                raise ValueError(
                    f"refusing to overwrite existing file with different content: {destination}"
                )
            if move:
                shutil.move(str(source), destination)
            else:
                shutil.copy2(source, destination)

    if verify:
        verify_samples(samples, output_dir, dry_run=dry_run)
        print("Verification passed" + (" (dry-run)" if dry_run else ""))
    return samples


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", nargs="?", type=Path, default=Path("data/raw/GSE317012_RAW"))
    parser.add_argument("--output-dir", type=Path, help="Destination root; defaults to source_dir")
    parser.add_argument(
        "--extract",
        type=Path,
        metavar="TAR",
        help="Safely extract TAR into source_dir before organizing",
    )
    parser.add_argument(
        "--move",
        action="store_true",
        help="Move files instead of copying (no flat duplicates retained)",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print planned copies without changing files"
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Require complete sample groups and verify destinations",
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.extract is not None:
            written = extract_tar(args.extract, args.source_dir)
            print(f"Extracted {len(written)} file(s) from {args.extract} into {args.source_dir}")
        samples = reorganize(
            args.source_dir,
            args.output_dir,
            dry_run=args.dry_run,
            verify=args.verify,
            move=args.move,
        )
    except FileNotFoundError as exc:
        print(MISSING_INPUT_NOTE.format(path=exc.filename or args.source_dir), file=sys.stderr)
        return 1
    print(f"Discovered {len(samples)} GSM sample group(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
