"""Clinical bridge (copies/mL) integrity tests.

The bridge is piecewise log-linear through assumption-labelled anchors; the
contract is that anchors are reproduced EXACTLY (previous Hill-based version
missed its own anchors by 3-4 orders of magnitude).
"""

import numpy as np
import pytest

from vcm.clinical.viral_load_mapper import ViralLoadMapper


@pytest.fixture(scope="module")
def mapper(tmp_path_factory):
    d = tmp_path_factory.mktemp("bridge")
    return ViralLoadMapper(params_path=str(d / "params.json"))


def test_anchors_reproduced_exactly(mapper):
    for anchor in mapper.params["anchors"]:
        v = anchor["viral_load_v"]
        expected = anchor["copies_per_ml"]
        got = mapper.normalized_to_copies(v)
        assert abs(got - expected) <= max(1e-6, 1e-6 * expected)


def test_bridge_monotonic(mapper):
    vs = np.linspace(0.0, 4.0, 81)
    copies = [mapper.normalized_to_copies(v) for v in vs]
    diffs = np.diff(copies)
    assert (diffs >= -1e-9).all(), "bridge is not monotone non-decreasing"


def test_bridge_below_detection(mapper):
    assert mapper.normalized_to_copies(0.0) == 0.0
    assert mapper.normalized_to_copies(0.0) < 100.0


def test_drug_ordering_through_bridge(mapper):
    """infection < sirolimus < ... < tacrolimus for peak copies (regime test)."""
    peaks = {}
    for scenario in ("infection", "tacrolimus", "sirolimus"):
        df = mapper.simulate_clinical_trajectory(scenario, weeks=30)
        peaks[scenario] = float(df["copies_per_ml"].max())
    assert peaks["sirolimus"] < peaks["infection"] < peaks["tacrolimus"]


def test_weekly_grid_shape(mapper):
    df = mapper.simulate_clinical_trajectory("infection", weeks=12)
    assert len(df) == 13
    assert df["week"].tolist() == list(range(13))
    assert (df["viral_load_norm"].values >= 0).all()


def test_no_silent_hill_collapse(mapper):
    """A mid-scale value should sit between its neighbouring anchors."""
    # V=0.6 lies between anchors 0.2 (1k) and 1.0 (10k): log-linear ~3.16k
    v = mapper.normalized_to_copies(0.6)
    assert 1000.0 < v < 10000.0
