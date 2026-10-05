"""Tests for simulator implementations."""

import numpy as np
import pytest

from vcm.core.models import CellState, Environment, Gene, Protein
from vcm.simulators.base import BaseSimulator
from vcm.simulators.hybrid import HybridSimulator
from vcm.simulators.mechanistic import MechanisticSimulator
from vcm.simulators.ml_based import MLBasedSimulator


def test_base_simulator_interface():
    assert hasattr(BaseSimulator, "simulate")
    assert hasattr(BaseSimulator, "step")
    assert hasattr(BaseSimulator, "get_simulator_info")


def test_mechanistic_simulator_creation():
    sim = MechanisticSimulator()
    assert sim.decay_rate == 0.05
    assert sim.growth_rate == 0.1
    sim_with_config = MechanisticSimulator({"decay_rate": 0.1, "growth_rate": 0.2})
    assert sim_with_config.decay_rate == 0.1
    assert sim_with_config.growth_rate == 0.2


def test_mechanistic_simulator_step():
    sim = MechanisticSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )
    new_state = sim.step(state, None, Environment(), 1.0)
    assert new_state.timestamp == 1.0
    assert new_state.cell_id == "cell1"
    assert new_state.genes["gene1"].expression_level != 1.0


def test_mechanistic_simulator_simulate():
    sim = MechanisticSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )
    result = sim.simulate(state, None, Environment(), n_steps=10, timestep=1.0)
    assert result.simulator_type == "mechanistic"
    assert len(result.steps) == 10
    assert result.final_state.timestamp == 10.0


def test_ml_based_simulator_creation():
    sim = MLBasedSimulator()
    assert sim.model is None


def test_ml_based_simulator_training():
    teacher = MechanisticSimulator()
    sim = MLBasedSimulator()
    sim.train(teacher, n_steps=50)
    assert sim.model is not None
    assert sim._input_size is not None
    assert sim._output_size is not None
    assert hasattr(sim.model, "predict")


def test_ml_based_simulator_step_after_training():
    teacher = MechanisticSimulator()
    sim = MLBasedSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )
    sim.train(teacher, n_steps=50, initial_state=state)
    new_state = sim.step(state, None, Environment(), 1.0)
    assert new_state.timestamp == 1.0
    assert new_state.cell_id == "cell1"


def test_ml_based_simulator_auto_train_on_simulate():
    sim = MLBasedSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )
    result = sim.simulate(state, None, Environment(), n_steps=10, timestep=1.0)
    assert sim.model is not None
    assert result.simulator_type == "ml_based"
    assert len(result.steps) == 10
    assert result.final_state is not None


def test_ml_based_simulator_save_load(tmp_path):
    teacher = MechanisticSimulator()
    sim = MLBasedSimulator()
    sim.train(teacher, n_steps=50)
    model_path = str(tmp_path / "model.pkl")
    sim.save_model(model_path)

    loaded = MLBasedSimulator()
    loaded.load_model(model_path)
    assert loaded.model is not None
    assert loaded._input_size == sim._input_size
    assert loaded._output_size == sim._output_size


def test_ml_based_simulator_prediction_error_bounded():
    sim = MLBasedSimulator()
    teacher = MechanisticSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )

    sim.train(teacher, n_steps=300, initial_state=state)

    teacher_result = teacher.simulate(state, n_steps=10)
    for i in range(len(teacher_result.steps) - 1):
        cur = teacher_result.steps[i].cell_state
        nxt = teacher_result.steps[i + 1].cell_state

        pred = sim.step(cur)
        pred_vec = pred.get_state_vector()
        target_vec = nxt.get_state_vector()

        mse = sum((a - b) ** 2 for a, b in zip(pred_vec, target_vec)) / max(len(pred_vec), 1)
        assert mse < 5.0, f"MSE too high at step {i}: {mse}"


def test_ml_based_simulator_get_simulator_info():
    sim = MLBasedSimulator()
    info = sim.get_simulator_info()
    assert info["simulator_type"] == "MLBasedSimulator"
    assert info["model_trained"] is False

    teacher = MechanisticSimulator()
    sim.train(teacher, n_steps=30)
    info = sim.get_simulator_info()
    assert info["model_trained"] is True
    assert info["input_size"] is not None
    assert info["output_size"] is not None


def test_hybrid_simulator_creation():
    sim = HybridSimulator()
    assert sim.mechanistic_weight == 0.7
    assert sim.ml_weight == 0.3


def test_hybrid_simulator_step():
    sim = HybridSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )
    new_state = sim.step(state, None, Environment(), 1.0)
    assert new_state.timestamp == 1.0
    assert new_state.cell_id == "cell1"


def test_hybrid_simulator_simulate():
    sim = HybridSimulator()
    state = CellState(
        cell_id="cell1", cell_type="test_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={}, pathways={},
    )
    result = sim.simulate(state, None, Environment(), n_steps=10, timestep=1.0)
    assert result.simulator_type == "hybrid"
    assert len(result.steps) == 10
    assert result.final_state is not None


def test_simulator_info():
    sim = MechanisticSimulator()
    info = sim.get_simulator_info()
    assert "simulator_type" in info
    assert info["simulator_type"] == "MechanisticSimulator"
    assert "config" in info
