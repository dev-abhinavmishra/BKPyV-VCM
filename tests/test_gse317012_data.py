"""Tests for GSE317012 data integration: organizer, extraction, metadata, analyzer."""

import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

from scripts import reorganize_gse317012 as reorg


def _write_flat_sample(directory: Path, gsm: str = "GSM100", suffix: str = "T_1_F1"):
    """Create one GSM's six files: filtered + raw_ triplets."""
    for kind, ext in (("matrix", "mtx"), ("barcodes", "tsv"), ("features", "tsv")):
        (directory / f"{gsm}_{suffix}_{kind}.{ext}.gz").write_bytes(f"{kind}".encode())
        (directory / f"{gsm}_{suffix}_raw_{kind}.{ext}.gz").write_bytes(f"raw_{kind}".encode())


def _make_tar(path: Path, members: dict):
    with tarfile.open(path, "w") as tf:
        for name, data in members.items():
            payload = data if isinstance(data, bytes) else data.encode()
            ti = tarfile.TarInfo(name)
            ti.size = len(payload)
            tf.addfile(ti, io.BytesIO(payload))


# ---------- organizer: raw vs filtered regression ----------


def test_discover_samples_raw_filtered_coexist(tmp_path):
    """Regression: raw_ triplet must not collapse into filtered kinds."""
    _write_flat_sample(tmp_path)
    samples = reorg.discover_samples(tmp_path)
    kinds = set(samples["GSM100"])
    assert kinds == {
        "matrix",
        "barcodes",
        "features",
        "raw_matrix",
        "raw_barcodes",
        "raw_features",
    }


def test_organize_move_six_artifacts_no_flat_dupes(tmp_path):
    src = tmp_path / "flat"
    src.mkdir()
    _write_flat_sample(src)
    reorg.reorganize(src, move=True)
    dest = src / "GSM100"
    names = sorted(p.name for p in dest.iterdir())
    assert names == [
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
        "raw_barcodes.tsv.gz",
        "raw_features.tsv.gz",
        "raw_matrix.mtx.gz",
    ]
    assert not list(src.glob("*.gz")), "flat duplicate copies must not remain after move"


# ---------- extraction safety ----------


def test_extract_rejects_dotdot_traversal(tmp_path):
    tar = tmp_path / "evil.tar"
    _make_tar(tar, {"../escape.txt": b"x_mat"})
    with pytest.raises(ValueError, match="traversal|unsafe|escape"):
        reorg.extract_tar(tar, tmp_path / "out")


def test_extract_rejects_absolute_path(tmp_path):
    tar = tmp_path / "abs.tar"
    _make_tar(tar, {"/abs/evil.txt": b"x_mat"})
    with pytest.raises(ValueError, match="traversal|unsafe|absolute"):
        reorg.extract_tar(tar, tmp_path / "out")


def test_extract_refuses_foreign_overwrite(tmp_path):
    tar = tmp_path / "data.tar"
    _make_tar(tar, {"a.txt": b"original"})
    dest = tmp_path / "out"
    reorg.extract_tar(tar, dest)
    (dest / "a.txt").write_bytes(b"foreign-change")
    with pytest.raises(ValueError, match="refuse|overwrite|differ"):
        reorg.extract_tar(tar, dest)


def test_extract_idempotent_identical(tmp_path):
    tar = tmp_path / "data.tar"
    _make_tar(tar, {"a.txt": b"original"})
    dest = tmp_path / "out"
    reorg.extract_tar(tar, dest)
    before = (dest / "a.txt").read_bytes()
    written = reorg.extract_tar(tar, dest)
    assert written == [] or all(p == dest / "a.txt" for p in written)
    assert (dest / "a.txt").read_bytes() == before


def test_organize_repeat_run_is_noop(tmp_path):
    src = tmp_path / "flat"
    src.mkdir()
    _write_flat_sample(src)
    reorg.reorganize(src, move=True)
    second = reorg.reorganize(src, move=True)
    assert second == {}
    names = sorted(p.name for p in (src / "GSM100").iterdir())
    assert len(names) == 6


# ---------- metadata parsing ----------


SERIES_MATRIX_FIXTURE = (
    '!Series_title\t"fake"\n'
    '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n'
    '!Sample_title\t"sample one"\t"sample two"\t"sample three"\n'
    '!Sample_characteristics_ch1\t"disease phase: control"\t'
    '"disease phase: peaking"\t"disease phase: resolving"\n'
    '!Sample_source_name_ch1\t"kidney biopsy"\t"kidney biopsy"\t"kidney biopsy"\n'
)


def test_parse_series_matrix_phases(tmp_path):
    f = tmp_path / "series.txt"
    f.write_text(SERIES_MATRIX_FIXTURE)
    from scripts import analyze_gse317012 as anz

    parsed = anz.parse_series_matrix(f)
    assert parsed["GSM1"]["phase"] == "Control"
    assert parsed["GSM2"]["phase"] == "Peaking"
    assert parsed["GSM3"]["phase"] == "Resolving"
    assert parsed["GSM1"]["source_name"] == "kidney biopsy"


def test_metadata_crosscheck_mismatch_raises(tmp_path):
    f = tmp_path / "series.txt"
    f.write_text(SERIES_MATRIX_FIXTURE)
    from scripts import analyze_gse317012 as anz

    parsed = anz.parse_series_matrix(f)
    with pytest.raises(ValueError, match="GSM"):
        anz.crosscheck_gsm_sets(set(parsed), {"GSM1", "GSM2", "GSM9"})


# ---------- analyzer math on tiny sparse fixtures ----------


def _tiny_matrix():
    """3 genes x_mat 4 cells; col sums chosen so normalization is hand-checkable."""
    # genes: MKI67 (index0), LRP2 (1), PTPRC (2)
    x_mat = sp.csr_matrix(
        np.array(
            [
                [100, 0, 50, 0],
                [10, 10, 0, 40],
                [0, 90, 0, 10],
            ],
            dtype=float,
        )
    )
    features = ["MKI67", "LRP2", "PTPRC"]
    barcodes = ["c1", "c2", "c3", "c4"]
    return x_mat, features, barcodes


def test_sparse_score_and_fraction_correctness():
    from scripts import analyze_gse317012 as anz

    x_mat, features, _ = _tiny_matrix()
    norm = anz.normalize_log1p(x_mat)
    scores = anz.gene_set_scores(norm, features, {"sg2m": ["MKI67"]})
    col_sums = np.asarray(x_mat.sum(axis=0)).ravel()
    expected = np.log1p(np.array([100, 0, 50, 0]) / col_sums * 1e4)
    assert np.allclose(scores["sg2m"], expected)
    fracs = anz.marker_fractions(x_mat, features, {"epi": ["LRP2"], "imm": ["PTPRC"]})
    assert fracs["epi"] == pytest.approx(3 / 4)
    assert fracs["imm"] == pytest.approx(2 / 4)


def test_qc_exclusion_accounting():
    from scripts import analyze_gse317012 as anz

    x_mat = sp.csr_matrix(
        np.array(
            [
                [500, 10, 300, 200],
                [500, 5, 300, 200],
                [0, 10, 0, 0],
            ],
            dtype=float,
        )
    )
    features = ["A", "B", "MT-CO1"]
    keep, stats = anz.qc_filter(x_mat, features, min_genes=2, max_mt_frac=0.2)
    assert stats["cells_pre"] == 4
    assert stats["cells_post"] == 3
    assert stats["excluded"] == 1
    assert keep.shape[1] == 3


def test_output_schema_and_determinism(tmp_path):
    from scripts import analyze_gse317012 as anz

    x_mat, features, barcodes = _tiny_matrix()
    sample_dir = tmp_path / "GSMX"
    sample_dir.mkdir()
    (sample_dir / "features.tsv").write_text("MKI67\nLRP2\nPTPRC\n")
    (sample_dir / "barcodes.tsv").write_text("c1\nc2\nc3\nc4\n")
    scipy.io.mmwrite(str(sample_dir / "matrix.mtx"), x_mat.tocoo())

    meta = {"GSMX": {"phase": "Control", "title": "t", "source_name": "s"}}
    row1 = anz.analyze_sample(sample_dir, meta["GSMX"])
    row2 = anz.analyze_sample(sample_dir, meta["GSMX"])
    assert pd.Series(row1).equals(pd.Series(row2))
    for col in (
        "gsm",
        "phase",
        "cells_pre_qc",
        "cells_post_qc",
        "frac_epithelial_proxy",
        "frac_immune_proxy",
        "score_s_g2m",
        "score_ddr",
        "score_mito",
        "score_translation",
        "score_antigen_presentation",
        "missing_markers",
    ):
        assert col in row1
