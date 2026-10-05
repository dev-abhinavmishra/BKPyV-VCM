"""Tests for the scanpy-free GSE317012 file organizer."""

from pathlib import Path

import pytest

from scripts.reorganize_gse317012 import reorganize


def _write_fixture(source: Path, sample: str = "GSM1234567") -> None:
    for suffix in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
        (source / f"{sample}_T_11_F2_{suffix}").write_bytes(b"fixture")


def test_reorganize_groups_files_and_preserves_contents(tmp_path):
    source = tmp_path / "flat"
    output = tmp_path / "grouped"
    source.mkdir()
    _write_fixture(source)

    reorganize(source, output, verify=True)

    sample_dir = output / "GSM1234567"
    assert sorted(path.name for path in sample_dir.iterdir()) == [
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
    ]
    assert (sample_dir / "matrix.mtx.gz").read_bytes() == b"fixture"


def test_dry_run_does_not_create_destination(tmp_path):
    source = tmp_path / "flat"
    output = tmp_path / "grouped"
    source.mkdir()
    _write_fixture(source)

    reorganize(source, output, dry_run=True, verify=True)

    assert not output.exists()


def test_verify_rejects_incomplete_sample(tmp_path):
    source = tmp_path / "flat"
    source.mkdir()
    (source / "GSM1234567_matrix.mtx.gz").write_bytes(b"fixture")

    with pytest.raises(ValueError, match="Incomplete"):
        reorganize(source, tmp_path / "grouped", dry_run=True, verify=True)
