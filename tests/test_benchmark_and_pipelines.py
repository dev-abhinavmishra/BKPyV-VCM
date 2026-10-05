"""Tests for the Phase-2 viral-load benchmark and pipeline helpers.

Covers the units added by scripts/benchmark_viral_load.py,
scripts/download_gse317012.py, and scripts/cluster_gse317012.py — the
pure helpers (metrics, hash verification, deterministic holdout split)
plus a fast end-to-end benchmark smoke run. No network access required.
"""

import hashlib
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))


def _load_script(name):
    """Import a scripts/*.py file as a module."""
    spec = importlib.util.spec_from_file_location(
        name.replace(".py", ""), REPO / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


benchmark = _load_script("benchmark_viral_load.py")
download = _load_script("download_gse317012.py")


# ---------------------------------------------------------------------------
# benchmark_viral_load.py — pure metric helpers
# ---------------------------------------------------------------------------

class TestBandExcessRmse:
    def test_inside_band_is_zero(self):
        model = np.array([100.0, 50.0, 10.0])
        lo = np.array([1.0, 1.0, 1.0])
        hi = np.array([1000.0, 1000.0, 1000.0])
        rmse, frac = benchmark.band_excess_rmse(model, lo, hi)
        assert rmse == 0.0
        assert frac == 1.0

    def test_above_band_penalized(self):
        model = np.array([1e5])          # log10 = 5
        lo = np.array([1e3])             # log10 = 3
        hi = np.array([1e4])             # log10 = 4
        rmse, frac = benchmark.band_excess_rmse(model, lo, hi)
        assert rmse == pytest.approx(1.0)  # 1 decade above the band
        assert frac == 0.0

    def test_below_band_penalized(self):
        model = np.array([1.0])          # log10 = 0
        lo = np.array([100.0])
        hi = np.array([1000.0])
        rmse, _ = benchmark.band_excess_rmse(model, lo, hi)
        assert rmse == pytest.approx(2.0)


class TestHalfLifeAndThreshold:
    def test_half_life_exact(self):
        t = np.linspace(0, 10, 101)
        copies = 1e4 * 0.5 ** (t / 2.0)      # true t1/2 = 2 days
        assert benchmark.clearance_half_life(t, copies) == pytest.approx(2.0, abs=0.11)

    def test_half_life_nan_when_never_halves(self):
        t = np.linspace(0, 10, 101)
        copies = np.full_like(t, 5000.0)
        assert np.isnan(benchmark.clearance_half_life(t, copies))

    def test_time_below_threshold(self):
        t = np.linspace(0, 10, 101)
        copies = np.where(t < 5, 2000.0, 500.0)
        assert benchmark.time_below(t, copies, 1000.0) == pytest.approx(5.0, abs=0.11)

    def test_time_below_never(self):
        t = np.linspace(0, 10, 11)
        copies = np.full_like(t, 5000.0)
        assert benchmark.time_below(t, copies, 1000.0) == np.inf


class TestPublishedTable:
    def test_table_schema_and_citations(self):
        df = pd.DataFrame(benchmark.PUBLISHED)
        for col in ("source", "citation", "quantity", "value", "unit", "context"):
            assert col in df.columns
        assert (df["citation"].str.len() > 20).all()      # every row cites a source
        assert {"funk2006", "funk2008"}.issubset(set(df["source"]))

    def test_writes_csv(self, tmp_path):
        out = benchmark.write_published_table(tmp_path)
        df = pd.read_csv(out)
        assert len(df) == len(benchmark.PUBLISHED)


# ---------------------------------------------------------------------------
# Benchmark end-to-end smoke (fast; the scientific checks live in the script)
# ---------------------------------------------------------------------------

class TestBenchmarkIntegration:
    @pytest.mark.slow
    def test_quick_benchmark_produces_report(self, tmp_path):
        sys.argv = ["benchmark_viral_load.py", "--quick",
                    "--outdir", str(tmp_path), "--research-dir", str(tmp_path)]
        benchmark.main()
        import json
        report = json.loads((tmp_path / "viral_load_benchmark.json").read_text())
        for key in ("calibration", "validation", "parameter_uncertainty",
                    "data_limitation", "overall_pass"):
            assert key in report
        # qualitative check that survives --quick horizons:
        # 50% curtailment must not achieve sustained clearance
        res = report["validation"]["funk2008_curtailment"]["results"]
        assert res["curtail_50"]["sustained_clearance"] is False


# ---------------------------------------------------------------------------
# download_gse317012.py — sha256 verification helper
# ---------------------------------------------------------------------------

class TestDownloadHelpers:
    def test_sha256_of_file(self, tmp_path):
        f = tmp_path / "x.bin"
        f.write_bytes(b"hello world")
        expected = hashlib.sha256(b"hello world").hexdigest()
        assert download.sha256_file(f) == expected

    def test_verify_case_insensitive(self, tmp_path):
        f = tmp_path / "x.bin"
        f.write_bytes(b"hello world")
        expected = hashlib.sha256(b"hello world").hexdigest()
        assert download.sha256_matches(f, expected.upper()) is True
        assert download.sha256_matches(f, expected.lower()) is True
        assert download.sha256_matches(f, "deadbeef") is False


# ---------------------------------------------------------------------------
# cluster_gse317012.py — deterministic holdout split (no data needed)
# ---------------------------------------------------------------------------

class TestClusteringHelpers:
    @staticmethod
    def _meta():
        """Real signature: dict of gsm -> {'phase': ...}."""
        phases = (["control"] * 12 + ["peaking"] * 5 + ["resolving"] * 9)
        return {f"GSM{i}": {"phase": p} for i, p in enumerate(phases)}

    def test_split_deterministic(self):
        cluster = _load_script("cluster_gse317012.py")
        meta = self._meta()
        a = cluster.split_samples(meta, seed=7, frac=1.0 / 3.0)
        b = cluster.split_samples(meta, seed=7, frac=1.0 / 3.0)
        assert a.equals(b)
        assert set(a["split"]) <= {"calibration", "validation"}
        # every phase contributes to both splits
        for phase in ("control", "peaking", "resolving"):
            sub = a[a["phase"] == phase]
            assert {"calibration", "validation"}.issubset(set(sub["split"]))

    def test_split_changes_with_seed(self):
        cluster = _load_script("cluster_gse317012.py")
        meta = self._meta()
        a = cluster.split_samples(meta, seed=7, frac=1.0 / 3.0)
        b = cluster.split_samples(meta, seed=11, frac=1.0 / 3.0)
        assert not a["split"].equals(b["split"])
