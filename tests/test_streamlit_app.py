"""Streamlit AppTest coverage for the BKPyV UI pages."""

from pathlib import Path
from types import SimpleNamespace

from streamlit.testing.v1 import AppTest

from vcm.ui import streamlit_app

REPO_ROOT = Path(__file__).resolve().parents[1]
APP_PATH = REPO_ROOT / "src" / "vcm" / "ui" / "streamlit_app.py"




def _app():
    return AppTest.from_file(str(APP_PATH), default_timeout=120)


def _open_page(at, name):
    """Navigate the sidebar radio to the named page (order-insensitive)."""
    options = list(at.sidebar.radio[0].options)
    target = next((o for o in options if name in o), None)
    assert target is not None, f"page {name!r} not in sidebar options {options}"
    at.sidebar.radio[0].set_value(target)
    at.run()
    return at


def _fake_result(viral_loads, timestep=1.0):
    steps = [
        SimpleNamespace(
            timestamp=i * timestep,
            cell_state=SimpleNamespace(metadata={"viral_load": v}),
        )
        for i, v in enumerate(viral_loads)
    ]
    return SimpleNamespace(steps=steps)


def _entry(viral_loads, scenario="infection_no_drug", nccr="archetype", timestep=1.0, horizon=None):
    horizon = horizon if horizon is not None else timestep * (len(viral_loads) - 1)
    config = {
        "experiment_id": f"bkpyv_{scenario}",
        "plugin": "transplant.bk_polyomavirus",
        "simulator": "bkpyv_ode",
        "simulation_length": float(horizon),
        "timestep": float(timestep),
        "output_path": "outputs/bkpyv/ui/",
        "scenario": scenario,
        "nccr_variant": nccr,
        "simulator_parameters": {},
    }
    return {"result": _fake_result(viral_loads, timestep), "config": config, "label": "seeded"}


# ---------- all-page coverage ----------


def test_all_ten_pages_render_with_expected_content():
    """Every sidebar page loads with no exception and page-specific content."""
    at = _app()
    at.run()
    assert len(at.exception) == 0
    options = list(at.sidebar.radio[0].options)
    assert len(options) == 10

    expectations = {
        "Home": lambda a: any("Welcome to the BKPyV" in m.value for m in a.markdown),
        "Simulation": lambda a: any("Choose simulation scenario" == s.label for s in a.selectbox),
        "Visualization": lambda a: any(
            "No simulation results available" in w.value for w in a.warning
        ),
        "Single-Cell": lambda a: any("Single-Cell Analysis" in m.value for m in a.markdown),
        "Viral-Load Validation": lambda a: any(
            "Viral-Load Validation" in m.value for m in a.markdown
        ),
        "Risk Prediction": lambda a: any("Patient Clinical Data" == s.value for s in a.subheader),
        "Comparison": lambda a: any("populate this page" in i.value for i in a.info),
        "Parameters & Assumptions": lambda a: any(
            "Parameters & Assumptions" in m.value for m in a.markdown
        ),
        "Review Bundle": lambda a: any("Generate review bundle" in b.label for b in a.button),
        "Documentation": lambda a: any("Clinical Background" in s.value for s in a.subheader),
    }
    for name, check in expectations.items():
        at2 = _app()
        at2.run()
        _open_page(at2, name)
        assert len(at2.exception) == 0, f"{name} raised: {at2.exception}"
        assert check(at2), f"{name} missing expected content"


# ---------- comparison page: empty / insufficient / populated ----------


def test_comparison_empty_state():
    at = _app()
    at.run()
    _open_page(at, "Comparison")
    assert len(at.exception) == 0
    assert any("populate this page" in i.value for i in at.info)
    assert any("Stored runs available" in m.value for m in at.markdown)
    assert len(at.dataframe) == 0


def test_comparison_insufficient_selection_warns():
    at = _app()
    at.run()
    at.session_state["simulation_history"] = {
        "run 1: a": _entry([0.0, 0.5, 0.2]),
        "run 2: b": _entry([0.0, 0.3, 0.7]),
    }
    _open_page(at, "Comparison")
    ms = at.multiselect[0]
    ms.set_value(list(ms.options)[:1])
    at.run()
    assert len(at.exception) == 0
    assert any("at least two stored runs" in w.value for w in at.warning)


def test_comparison_populated_metrics_equal_stored_values():
    loads_a = [0.0, 0.4, 0.9, 0.3]
    loads_b = [0.0, 0.2, 0.1]
    at = _app()
    at.run()
    at.session_state["simulation_history"] = {
        "run 1: infection_no_drug | NCCR archetype": _entry(loads_a, scenario="infection_no_drug"),
        "run 2: tacrolimus_exposure | NCCR rearranged": _entry(
            loads_b, scenario="tacrolimus_exposure", nccr="rearranged"
        ),
    }
    _open_page(at, "Comparison")
    assert len(at.exception) == 0
    assert len(at.dataframe) >= 1
    df = at.dataframe[0].value
    assert set(df["run"]) == {
        "run 1: infection_no_drug | NCCR archetype",
        "run 2: tacrolimus_exposure | NCCR rearranged",
    }
    row_a = df[df["run"].str.contains("infection_no_drug")].iloc[0]
    row_b = df[df["run"].str.contains("tacrolimus")].iloc[0]
    assert row_a["peak_virtual_load"] == max(loads_a)
    assert row_b["peak_virtual_load"] == max(loads_b)
    assert row_a["endpoint_load"] == loads_a[-1]
    assert row_b["endpoint_load"] == loads_b[-1]
    assert row_a["peak_day"] == 2.0
    assert row_b["peak_day"] == 1.0
    assert row_a["scenario"] == "infection_no_drug"
    assert row_b["nccr_variant"] == "rearranged"
    assert any("shared day window" in c.value for c in at.caption), "alignment caption missing"


def test_comparison_clinical_bridge_carries_assumption_label():
    at = _app()
    at.run()
    at.session_state["simulation_history"] = {
        "run 1: a": _entry([0.0, 0.5, 0.2]),
        "run 2: b": _entry([0.0, 0.3, 0.7]),
    }
    _open_page(at, "Comparison")
    df = at.dataframe[0].value
    assert "peak_plasma_copies_per_ml" not in df.columns
    assert "peak_risk_category" not in df.columns
    at.checkbox[0].set_value(True)
    at.run()
    assert len(at.exception) == 0
    df = at.dataframe[0].value
    assert "peak_plasma_copies_per_ml" in df.columns
    assert "peak_risk_category" in df.columns
    assert any(
        "Assumption-labelled bridge" in c.value and "not a calibration" in c.value
        for c in at.caption
    )


# ---------- pure-helper unit tests ----------


def test_shared_time_window_overlap_and_none():
    e1 = _entry([0.0, 0.1, 0.2], timestep=1.0)  # 0..2
    e2 = _entry([0.0, 0.1, 0.2, 0.3, 0.4], timestep=1.0)  # 0..4
    assert streamlit_app._shared_time_window([e1, e2]) == (0.0, 2.0)
    e3 = _entry([0.0, 0.1], timestep=1.0)
    e3["result"].steps = [
        SimpleNamespace(timestamp=10.0 + i, cell_state=SimpleNamespace(metadata={"viral_load": v}))
        for i, v in enumerate([0.5, 0.6])
    ]
    assert streamlit_app._shared_time_window([e1, e3]) is None


def test_windowed_series_restricts_to_window():
    traj = {"timepoints": [0.0, 1.0, 2.0, 3.0, 4.0], "viral_loads": [0, 1, 2, 3, 4]}
    xs, ys = streamlit_app._windowed_series(traj, (1.0, 3.0))
    assert xs == [1.0, 2.0, 3.0]
    assert ys == [1, 2, 3]


def test_comparison_metrics_truthful_to_stored():
    entry = _entry([0.0, 0.7, 0.2], timestep=0.5)
    metrics = streamlit_app._comparison_metrics("r", entry)
    assert metrics["peak_virtual_load"] == 0.7
    assert metrics["peak_day"] == 0.5
    assert metrics["endpoint_load"] == 0.2
    assert metrics["timestep_days"] == 0.5
    assert metrics["engine"] == "bkpyv_ode"


# ---------- real simulation stores a comparable run ----------


def test_simulation_run_stores_history_entry():
    """A real UI run must land in simulation_history for the Comparison page."""
    at = _app()
    at.run()
    _open_page(at, "Simulation")
    run_buttons = [b for b in at.button if "Run Simulation" in b.label]
    assert run_buttons, "run button not found"
    run_buttons[0].click()
    at.run()
    assert len(at.exception) == 0
    history = at.session_state["simulation_history"]
    assert len(history) == 1
    entry = next(iter(history.values()))
    assert entry["result"] is not None
    assert entry["config"]["simulator"] == "bkpyv_ode"
    assert "infection" in entry["label"] or "run 1" in entry["label"]
