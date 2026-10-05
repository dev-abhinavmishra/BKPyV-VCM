"""ML-based simulator using a trained neural network for state transitions."""

import copy
import pickle
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from sklearn.neural_network import MLPRegressor

from vcm.core.models import (
    CellState,
    Environment,
    Perturbation,
    SimulationResult,
    SimulationStep,
)
from vcm.simulators.base import BaseSimulator


class MLBasedSimulator(BaseSimulator):
    """ML-based simulator using a trained neural network.

    Learns state transitions from a teacher simulator via supervised
    learning. An MLPRegressor predicts the delta (change) between
    consecutive states given the current state vector.

    Usage:
        sim = MLBasedSimulator()
        sim.train(teacher_simulator, n_steps=500)
        result = sim.simulate(initial_state, n_steps=100)
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self.model: Optional[MLPRegressor] = None
        self._input_size: Optional[int] = None
        self._output_size: Optional[int] = None

    def train(
        self,
        teacher: BaseSimulator,
        n_steps: int = 500,
        initial_state: Optional[CellState] = None,
    ) -> "MLBasedSimulator":
        """Train the ML model on state transitions from a teacher simulator.

        Runs the teacher simulator to generate trajectory data, then
        trains an MLPRegressor to predict state deltas.

        Args:
            teacher: A trained simulator to learn from.
            n_steps: Number of simulation steps for training data.
            initial_state: Starting state. Uses a default if not given.

        Returns:
            Self, with trained model.
        """
        if initial_state is None:
            initial_state = _default_cell_state()

        result = teacher.simulate(initial_state, n_steps=n_steps)
        if len(result.steps) < 2:
            raise ValueError("Need at least 2 steps for training")

        X, y = _build_training_data(result)
        self._input_size = X.shape[1]
        self._output_size = y.shape[1]

        self.model = MLPRegressor(
            hidden_layer_sizes=(max(64, self._input_size), max(32, self._input_size // 2)),
            activation="relu",
            solver="adam",
            max_iter=1000,
            random_state=42,
            early_stopping=True,
            validation_fraction=0.1,
            n_iter_no_change=20,
            verbose=False,
        )
        self.model.fit(X, y)
        return self

    def simulate(
        self,
        initial_state: CellState,
        perturbation: Optional[Perturbation] = None,
        environment: Optional[Environment] = None,
        n_steps: int = 100,
        timestep: float = 1.0,
    ) -> SimulationResult:
        if environment is None:
            environment = Environment()

        if self.model is None:
            self._auto_train(initial_state)

        result = SimulationResult(
            experiment_id=f"ml_based_{initial_state.cell_id}",
            simulator_type="ml_based",
            plugin=initial_state.cell_type,
            config_id="default",
        )

        current_state = copy.deepcopy(initial_state)
        current_state.timestamp = 0.0

        for step_num in range(n_steps):
            active_perturbations = []
            if perturbation and perturbation.timing is not None:
                if (
                    perturbation.timing <= current_state.timestamp
                    and (
                        perturbation.duration is None
                        or current_state.timestamp < perturbation.timing + perturbation.duration
                    )
                ):
                    active_perturbations.append(perturbation)

            step = SimulationStep(
                step_number=step_num,
                timestamp=current_state.timestamp,
                cell_state=copy.deepcopy(current_state),
                applied_perturbations=active_perturbations,
                environment=environment,
            )
            result.steps.append(step)

            current_state = self.step(current_state, perturbation, environment, timestep)

        result.final_state = current_state
        return result

    def step(
        self,
        current_state: CellState,
        perturbation: Optional[Perturbation] = None,
        environment: Optional[Environment] = None,
        timestep: float = 1.0,
    ) -> CellState:
        if environment is None:
            environment = Environment()

        if self.model is None:
            self._auto_train(current_state)

        new_state = copy.deepcopy(current_state)
        new_state.timestamp += timestep

        state_vec = np.array(current_state.get_state_vector(), dtype=np.float64).reshape(1, -1)
        n_features = state_vec.shape[1]

        if self._input_size is not None and n_features != self._input_size:
            state_vec = _resize_vector(state_vec, self._input_size)

        delta = self.model.predict(state_vec)[0]

        if self._output_size is not None and len(delta) != self._output_size:
            delta = _resize_1d(delta, self._output_size)

        index = 0
        for gene in new_state.genes.values():
            if index < len(delta):
                gene.expression_level = max(0.0, gene.expression_level + delta[index] * timestep * 0.1)
                index += 1
        for protein in new_state.proteins.values():
            if index < len(delta):
                protein.concentration = max(0.0, protein.concentration + delta[index] * timestep * 0.1)
                index += 1
        for metabolite in new_state.metabolites.values():
            if index < len(delta):
                metabolite.concentration = max(0.0, metabolite.concentration + delta[index] * timestep * 0.05)
                index += 1
        for pathway in new_state.pathways.values():
            if index < len(delta):
                pathway.flux = max(0.0, pathway.flux + delta[index] * timestep * 0.1)
                index += 1

        if perturbation and perturbation.timing is not None:
            if perturbation.timing <= current_state.timestamp:
                _apply_perturbation_ml(new_state, perturbation, timestep)

        env_factor = (environment.temperature - 37.0) / 37.0
        for gene in new_state.genes.values():
            gene.expression_level = max(0.0, gene.expression_level * (1.0 + env_factor * 0.05))
        for protein in new_state.proteins.values():
            protein.concentration = max(0.0, protein.concentration * (1.0 + env_factor * 0.05))

        new_state.update_state_vector()
        return new_state

    def save_model(self, path: str) -> None:
        """Save trained model to disk."""
        if self.model is None:
            raise ValueError("No trained model to save")
        payload = {
            "model": self.model,
            "input_size": self._input_size,
            "output_size": self._output_size,
        }
        with open(path, "wb") as f:
            pickle.dump(payload, f)

    def load_model(self, path: str) -> "MLBasedSimulator":
        """Load a trained model from disk."""
        with open(path, "rb") as f:
            payload = pickle.load(f)
        self.model = payload["model"]
        self._input_size = payload["input_size"]
        self._output_size = payload["output_size"]
        return self

    def _auto_train(self, state: CellState) -> None:
        """Auto-train on a MechanisticSimulator when no model is provided."""
        from vcm.simulators.mechanistic import MechanisticSimulator
        self.train(MechanisticSimulator(), n_steps=200, initial_state=state)

    def get_simulator_info(self) -> Dict[str, Any]:
        info = super().get_simulator_info()
        info["model_trained"] = self.model is not None
        if self.model is not None:
            info["input_size"] = self._input_size
            info["output_size"] = self._output_size
            info["n_layers"] = len(self.model.coefs_)
        return info


def _build_training_data(result: SimulationResult):
    """Build feature/label arrays from simulation trajectory."""
    rows = []
    for i in range(len(result.steps) - 1):
        cur = result.steps[i].cell_state.get_state_vector()
        nxt = result.steps[i + 1].cell_state.get_state_vector()
        delta = [nxt[j] - cur[j] for j in range(len(cur))]
        rows.append((cur, delta))
    X = np.array([r[0] for r in rows], dtype=np.float64)
    y = np.array([r[1] for r in rows], dtype=np.float64)
    return X, y


def _resize_vector(vec: np.ndarray, target: int) -> np.ndarray:
    """Pad or truncate a 2-D feature vector to target size."""
    n = vec.shape[1]
    if n < target:
        out = np.zeros((vec.shape[0], target))
        out[:, :n] = vec
        return out
    return vec[:, :target]


def _resize_1d(arr: np.ndarray, target: int) -> np.ndarray:
    """Pad or truncate a 1-D array to target size."""
    n = len(arr)
    if n < target:
        out = np.zeros(target)
        out[:n] = arr
        return out
    return arr[:target]


def _apply_perturbation_ml(state: CellState, perturbation: Perturbation, timestep: float) -> None:
    """Apply perturbation effects to an ML-updated state."""
    mag = perturbation.magnitude
    tid = perturbation.target_id
    ptype = perturbation.perturbation_type.value if perturbation.perturbation_type else ""

    if "knockout" in ptype or "inhibition" in ptype:
        if tid and tid in state.genes:
            state.genes[tid].expression_level *= max(0.0, 1.0 - mag * timestep * 0.5)
        if tid and tid in state.proteins:
            state.proteins[tid].concentration *= max(0.0, 1.0 - mag * timestep * 0.5)
            state.proteins[tid].active = False
    elif "overexpression" in ptype and tid and tid in state.genes:
        state.genes[tid].expression_level *= (1.0 + mag * timestep * 0.5)
    elif tid and tid in state.metabolites:
        if "depletion" in ptype:
            state.metabolites[tid].concentration *= max(0.0, 1.0 - mag * timestep * 0.5)
        else:
            state.metabolites[tid].concentration += mag * timestep * 0.1


def _default_cell_state() -> CellState:
    """Create a minimal default cell state for training."""
    from vcm.core.models import Gene, Metabolite, Pathway, Protein
    return CellState(
        cell_id="default_train",
        cell_type="training_cell",
        genes={"gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0)},
        proteins={"prot1": Protein(id="prot1", name="Protein 1", concentration=1.0, gene_id="gene1")},
        metabolites={"atp": Metabolite(id="atp", name="ATP", concentration=5.0)},
        pathways={"metabolism": Pathway(id="metabolism", name="Metabolism", flux=1.0)},
    )
