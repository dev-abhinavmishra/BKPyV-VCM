"""Tests for simulator-backed BKPyV sensitivity analysis."""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from vcm.analysis.sensitivity import SensitivityAnalyzer


class _Mapper:
    @staticmethod
    def normalized_to_copies(value):
        return float(value) * 10.0


def test_sensitivity_runs_one_fresh_simulator_per_sample():
    simulator_configs = []
    perturbation_batches = []

    class FakeSimulator:
        def __init__(self, config):
            simulator_configs.append(dict(config))

        def simulate(self, initial_state, perturbations, environment, n_steps, timestep):
            perturbation_batches.append(perturbations)
            value = simulator_configs[-1]["t_antigen_replication_threshold"]
            step = SimpleNamespace(
                cell_state=SimpleNamespace(metadata={"viral_load": value})
            )
            return SimpleNamespace(steps=[step])

    with patch("vcm.analysis.sensitivity.ViralLoadMapper", return_value=_Mapper()):
        with patch("vcm.analysis.sensitivity.BKPyVODESimulator", FakeSimulator):
            analyzer = SensitivityAnalyzer(n_samples=3, random_state=7)
            result = analyzer.run_sensitivity_analysis(
                "t_antigen_replication_threshold", scenario="infection"
            )

    assert len(simulator_configs) == 3
    assert len(perturbation_batches) == 3
    # Infection at day 21 (week 3), matching the ODE/clinical day convention
    assert [batch[0].timing for batch in perturbation_batches] == [21.0] * 3
    assert result.outcome_values == [value * 10.0 for value in result.parameter_values]
    assert [config["t_antigen_replication_threshold"] for config in simulator_configs] == result.parameter_values


def test_sensitivity_rejects_unknown_parameter_and_scenario():
    with patch("vcm.analysis.sensitivity.ViralLoadMapper", return_value=_Mapper()):
        analyzer = SensitivityAnalyzer(n_samples=1)

    with pytest.raises(ValueError, match="not found"):
        analyzer.run_sensitivity_analysis("not_a_parameter")
    with pytest.raises(ValueError, match="Invalid scenario"):
        analyzer.run_sensitivity_analysis("t_antigen_replication_threshold", "unknown")
