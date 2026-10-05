"""Runner plumbing + reproducibility regression tests.

Covered invariants (from the 2026-09 audit):
- 'bkpyv_ode' is registered in the ExperimentRunner registry.
- simulator_parameters from configs actually reach the simulator constructor
  (previous versions silently dropped them).
- The legacy 'simulator_params' YAML key is accepted as an alias.
- The full perturbation list (not only the first entry) reaches simulators
  that support it.
- Perturbation.name defaults to the perturbation id when omitted.
"""

import numpy as np
import pytest

from vcm.core.models import ExperimentConfig, Perturbation, PerturbationType
from vcm.experiments.runner import ExperimentRunner
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.base import BaseSimulator


class _ProbeSimulator(BaseSimulator):
    """Captures the constructor config and perturbations actually passed."""

    captured = None

    def __init__(self, config=None):
        super().__init__(config)
        type(self).captured = dict(config or {})

    def step(self, current_state, perturbation=None, environment=None, timestep=1.0):
        return current_state

    def simulate(self, initial_state, perturbation=None, perturbations=None,
                 environment=None, n_steps=100, timestep=1.0):
        type(self).perturbation_ids = [
            p.id for p in (perturbations or ([perturbation] if perturbation else []))
        ]
        # Minimal fake result object
        from vcm.core.models import SimulationResult

        return SimulationResult(
            experiment_id="probe",
            simulator_type="probe",
            plugin=initial_state.cell_type,
            config_id="probe",
            steps=[],
            final_state=initial_state,
        )


@pytest.fixture()
def runner():
    # Ensure the plugin used by these tests is registered regardless of
    # test execution order (PluginRegistry is process-global).
    BKPolyomavirusPlugin().register()
    return ExperimentRunner()


@pytest.fixture()
def bkpyv_state():
    plugin = BKPolyomavirusPlugin()
    return plugin.create_initial_state()


def test_bkpyv_ode_registered(runner):
    assert "bkpyv_ode" in runner.simulator_registry
    assert "bkpyv_specific" in runner.simulator_registry


def test_simulator_parameters_reach_simulator(runner, bkpyv_state):
    runner.simulator_registry["probe"] = _ProbeSimulator
    config = ExperimentConfig(
        experiment_id="probe",
        plugin="transplant.bk_polyomavirus",
        simulator="probe",
        simulation_length=1.0,
        timestep=1.0,
        simulator_parameters={"p": 12.0, "mtor_inhibition": 0.25},
        metadata={"human_label": "not-a-simulator-param"},
    )
    runner.run_experiment(config, initial_state=bkpyv_state)

    assert _ProbeSimulator.captured["p"] == 12.0
    assert _ProbeSimulator.captured["mtor_inhibition"] == 0.25
    # Free-form metadata must NOT be passed to the simulator as config
    assert "human_label" not in _ProbeSimulator.captured


def test_legacy_simulator_params_alias_merges():
    config = ExperimentConfig(
        experiment_id="alias",
        plugin="transplant.bk_polyomavirus",
        simulator="bkpyv_ode",
        simulation_length=1.0,
        timestep=1.0,
        simulator_params={"magnitude_test": 3.3},  # legacy key
        simulator_parameters={"other": 1.0},
    )
    assert config.simulator_parameters["magnitude_test"] == 3.3
    assert config.simulator_parameters["other"] == 1.0


def test_all_perturbations_passed(runner, bkpyv_state):
    runner.simulator_registry["probe"] = _ProbeSimulator
    perts = [
        Perturbation(id="infection", perturbation_type=PerturbationType.VIRAL_INFECTION, timing=21.0),
        Perturbation(id="drug", perturbation_type=PerturbationType.DRUG_TREATMENT,
                     target_id="FKBP1A", magnitude=1.0, timing=0.0),
    ]
    config = ExperimentConfig(
        experiment_id="probe",
        plugin="transplant.bk_polyomavirus",
        simulator="probe",
        simulation_length=1.0,
        timestep=1.0,
        perturbations=perts,
    )
    runner.run_experiment(config, initial_state=bkpyv_state)
    assert _ProbeSimulator.perturbation_ids == ["infection", "drug"]


def test_perturbation_name_defaults_to_id():
    p = Perturbation(id="my_perturbation", perturbation_type=PerturbationType.DRUG_TREATMENT)
    assert p.name == "my_perturbation"


def test_time_unit_defaults_to_days():
    config = ExperimentConfig(
        experiment_id="x", plugin="p", simulator="mechanistic"
    )
    assert config.time_unit == "days"


def test_nccr_variant_flows_from_config_to_state_and_ode(runner):
    """simulator_parameters.nccr_variant must reach the plugin state and ODE params."""
    from vcm.cli.main import register_default_plugins

    register_default_plugins()
    config = ExperimentConfig(
        experiment_id="nccr_test",
        plugin="transplant.bk_polyomavirus",
        simulator="bkpyv_ode",
        simulation_length=20.0,
        timestep=1.0,
        simulator_parameters={"nccr_variant": "rearranged"},
        perturbations=[
            Perturbation(
                id="infection",
                perturbation_type=PerturbationType.VIRAL_INFECTION,
                magnitude=1.0,
                timing=1.0,
            )
        ],
    )
    result = runner.run_experiment(config)
    assert result.final_state is not None
    assert result.final_state.metadata["nccr_variant"] == "rearranged"


# ---------------------------------------------------------------------------
# ODE-level behaviour invariants


def test_ode_tolerances_accept_string_floats():
    """YAML 1.1 parses 1e-6 as a string; the simulator must coerce, not crash."""
    from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator

    sim = BKPyVODESimulator({"rtol": "1e-6", "atol": "1e-8", "ode_solver": "LSODA"})
    assert sim.rtol == pytest.approx(1e-6)
    assert sim.atol == pytest.approx(1e-8)


def test_cell_cycle_permissiveness_bounded_and_non_ratcheting():
    """CC is a bounded permissiveness variable; it must not monotonically lock at 1."""
    from scipy.integrate import solve_ivp

    from vcm.simulators.ode_system import BKPyVODESystem

    system = BKPyVODESystem({})
    y0 = system.get_infection_conditions(0.5)
    sol = solve_ivp(
        lambda t, y: system.ode_system(t, y, None),
        (0, 120), y0, method="LSODA", t_eval=np.arange(0.0, 121.0, 1.0),
    )
    cc = sol.y[6]
    assert np.all(cc <= 1.0 + 1e-6), "CC exceeded its [0,1] bound"
    assert np.all(cc >= -1e-6), "CC went negative"
    # No one-way ratchet: over the last 30 days CC must show downward movement
    tail = cc[-31:]
    assert (np.diff(tail) < -1e-6).any() or tail.max() < 0.9, (
        "CC appears locked at its ceiling (regression of the reset bug)"
    )
