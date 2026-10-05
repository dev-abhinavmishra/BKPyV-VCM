"""Experiment runner for managing simulation experiments."""

import inspect
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from vcm.core.models import CellState, ExperimentConfig, SimulationResult
from vcm.data.loader import DataLoader
from vcm.plugins.base import PluginRegistry
from vcm.simulators.hybrid import HybridSimulator
from vcm.simulators.mechanistic import MechanisticSimulator
from vcm.simulators.ml_based import MLBasedSimulator
from vcm.simulators.bkpyv_simulator import BKPyVSimulator
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator


class ExperimentRunner:
    """Run simulation experiments with plugins and configs."""

    def __init__(self, data_dir: Optional[Path] = None):
        """Initialize the experiment runner.

        Args:
            data_dir: Base directory for data
        """
        self.data_loader = DataLoader(data_dir)
        self.simulator_registry = {
            "mechanistic": MechanisticSimulator,
            "ml_based": MLBasedSimulator,
            "hybrid": HybridSimulator,
            "bkpyv_specific": BKPyVSimulator,
            "bkpyv_ode": BKPyVODESimulator,
        }

    def run_experiment(
        self,
        config: ExperimentConfig,
        initial_state: Optional[CellState] = None,
    ) -> SimulationResult:
        """Run a single experiment.

        Args:
            config: Experiment configuration
            initial_state: Optional initial state (if None, loads from plugin)

        Returns:
            SimulationResult
        """
        # Load plugin
        plugin = PluginRegistry.get(config.plugin)
        if plugin is None:
            raise ValueError(f"Plugin '{config.plugin}' not found")

        # Get initial state; pass simulator_parameters so state-level knobs
        # (e.g., BKPyV nccr_variant) flow from the config into the initial state.
        if initial_state is None:
            initial_state = plugin.create_initial_state(config.simulator_parameters or None)

        # Get simulator
        simulator_class = self.simulator_registry.get(config.simulator)
        if simulator_class is None:
            raise ValueError(f"Simulator '{config.simulator}' not found")

        # Simulator receives the declared simulator_parameters (plus plugin
        # kwargs such as nccr_variant), never the human-facing metadata.
        simulator_config = dict(config.simulator_parameters)
        if initial_state is not None:
            nccr_variant = initial_state.metadata.get("nccr_variant")
            if nccr_variant and "nccr_variant" not in simulator_config:
                simulator_config["nccr_variant"] = nccr_variant
        simulator = simulator_class(simulator_config)

        # Calculate number of steps
        n_steps = int(config.simulation_length / config.timestep)

        # Pass the full perturbation list to simulators that support it;
        # older simulators only accept a single perturbation.
        simulate_params = inspect.signature(simulator.simulate).parameters
        if "perturbations" in simulate_params:
            result = simulator.simulate(
                initial_state=initial_state,
                perturbations=list(config.perturbations),
                environment=config.environment,
                n_steps=n_steps,
                timestep=config.timestep,
            )
        else:
            perturbation_sequence = config.perturbations or [None]
            for perturbation in perturbation_sequence:
                result = simulator.simulate(
                    initial_state=initial_state,
                    perturbation=perturbation,
                    environment=config.environment,
                    n_steps=n_steps,
                    timestep=config.timestep,
                )
                initial_state = result.final_state or initial_state

        result.end_time = datetime.now(timezone.utc)
        result.config_id = config.experiment_id

        return result

    def save_result(self, result: SimulationResult, output_path: Path) -> None:
        """Save simulation result to files.

        Args:
            result: SimulationResult to save
            output_path: Directory to save results
        """
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save as JSON with custom serialization
        json_path = output_path / f"{result.experiment_id}.json"
        result_dict = {
            "experiment_id": result.experiment_id,
            "simulator_type": result.simulator_type,
            "plugin": result.plugin,
            "config_id": result.config_id,
            "start_time": result.start_time.isoformat() if result.start_time else None,
            "end_time": result.end_time.isoformat() if result.end_time else None,
            "success": result.success,
            "error_message": result.error_message,
            "metadata": result.metadata,
            "n_steps": len(result.steps),
        }
        self.data_loader.save_json(result_dict, json_path)

        # Save trajectory as CSV
        if result.steps:
            trajectory_data = []
            for step in result.steps:
                row = {
                    "step": step.step_number,
                    "timestamp": step.timestamp,
                }
                # Add gene expressions
                for gene_id, gene in step.cell_state.genes.items():
                    row[f"gene_{gene_id}"] = gene.expression_level
                # Add protein concentrations
                for protein_id, protein in step.cell_state.proteins.items():
                    row[f"protein_{protein_id}"] = protein.concentration
                # Add metabolite concentrations
                for metabolite_id, metabolite in step.cell_state.metabolites.items():
                    row[f"metabolite_{metabolite_id}"] = metabolite.concentration
                # Add pathway fluxes
                for pathway_id, pathway in step.cell_state.pathways.items():
                    row[f"pathway_{pathway_id}"] = pathway.flux

                trajectory_data.append(row)

            import pandas as pd

            df = pd.DataFrame(trajectory_data)
            csv_path = output_path / f"{result.experiment_id}_trajectory.csv"
            df.to_csv(csv_path, index=False)

    def load_config(self, config_path: Path) -> ExperimentConfig:
        """Load experiment configuration from YAML file.

        Args:
            config_path: Path to YAML config file

        Returns:
            ExperimentConfig
        """
        with open(config_path, "r") as f:
            config_data = yaml.safe_load(f)

        # Convert perturbation dicts to Perturbation objects
        perturbations = []
        if "perturbations" in config_data:
            for p_data in config_data["perturbations"]:
                from vcm.core.models import Perturbation, PerturbationType

                p_data["perturbation_type"] = PerturbationType(p_data["perturbation_type"])
                perturbations.append(Perturbation(**p_data))

        config_data["perturbations"] = perturbations

        # Convert environment dict to Environment object
        if "environment" in config_data:
            from vcm.core.models import Environment

            config_data["environment"] = Environment(**config_data["environment"])

        return ExperimentConfig(**config_data)

    def save_config(self, config: ExperimentConfig, config_path: Path) -> None:
        """Save experiment configuration to YAML file.

        Args:
            config: ExperimentConfig to save
            config_path: Path to save config file
        """
        config_path = Path(config_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)

        config_data = config.model_dump(mode="json")

        # Convert enum to string
        if config_data.get("perturbations"):
            for p in config_data["perturbations"]:
                if "perturbation_type" in p:
                    p["perturbation_type"] = p["perturbation_type"].value

        with open(config_path, "w") as f:
            yaml.dump(config_data, f, default_flow_style=False)


class ExperimentComparator:
    """Compare results from multiple experiments."""

    def __init__(self):
        """Initialize the experiment comparator."""
        self.results: Dict[str, SimulationResult] = {}

    def add_result(self, name: str, result: SimulationResult) -> None:
        """Add a result for comparison.

        Args:
            name: Name/identifier for this result
            result: SimulationResult to add
        """
        self.results[name] = result

    def compare_states(
        self,
        experiment1: str,
        experiment2: str,
        timepoint: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Compare cell states between two experiments.

        Args:
            experiment1: Name of first experiment
            experiment2: Name of second experiment
            timepoint: Optional specific timepoint to compare (default: final state)

        Returns:
            Dictionary with comparison data
        """
        if experiment1 not in self.results or experiment2 not in self.results:
            raise ValueError("Both experiments must be added before comparison")

        result1 = self.results[experiment1]
        result2 = self.results[experiment2]

        # Get states to compare
        if timepoint is None:
            state1 = result1.final_state
            state2 = result2.final_state
        else:
            # Find closest timepoint
            state1 = None
            state2 = None
            for step in result1.steps:
                if abs(step.timestamp - timepoint) < 0.01:
                    state1 = step.cell_state
                    break
            for step in result2.steps:
                if abs(step.timestamp - timepoint) < 0.01:
                    state2 = step.cell_state
                    break

            if state1 is None or state2 is None:
                raise ValueError(f"Timepoint {timepoint} not found in results")

        # Compare genes
        gene_diffs = {}
        for gene_id in state1.genes:
            if gene_id in state2.genes:
                diff = state2.genes[gene_id].expression_level - state1.genes[gene_id].expression_level
                gene_diffs[gene_id] = diff

        # Compare proteins
        protein_diffs = {}
        for protein_id in state1.proteins:
            if protein_id in state2.proteins:
                diff = (
                    state2.proteins[protein_id].concentration - state1.proteins[protein_id].concentration
                )
                protein_diffs[protein_id] = diff

        # Compare metabolites
        metabolite_diffs = {}
        for metabolite_id in state1.metabolites:
            if metabolite_id in state2.metabolites:
                diff = (
                    state2.metabolites[metabolite_id].concentration
                    - state1.metabolites[metabolite_id].concentration
                )
                metabolite_diffs[metabolite_id] = diff

        return {
            "experiment1": experiment1,
            "experiment2": experiment2,
            "timepoint": timepoint if timepoint else "final",
            "gene_differences": gene_diffs,
            "protein_differences": protein_diffs,
            "metabolite_differences": metabolite_diffs,
        }

    def compare_trajectories(
        self,
        experiment1: str,
        experiment2: str,
        entity_id: str,
        entity_type: str,
    ) -> Dict[str, Any]:
        """Compare trajectories for a specific entity between experiments.

        Args:
            experiment1: Name of first experiment
            experiment2: Name of second experiment
            entity_id: ID of entity to compare
            entity_type: Type of entity (gene, protein, metabolite, pathway)

        Returns:
            Dictionary with trajectory comparison data
        """
        if experiment1 not in self.results or experiment2 not in self.results:
            raise ValueError("Both experiments must be added before comparison")

        result1 = self.results[experiment1]
        result2 = self.results[experiment2]

        trajectory1 = result1.get_trajectory(entity_id, entity_type)
        trajectory2 = result2.get_trajectory(entity_id, entity_type)
        timepoints = result1.get_timepoints()

        # Calculate differences
        if len(trajectory1) == len(trajectory2):
            differences = [t2 - t1 for t1, t2 in zip(trajectory1, trajectory2)]
        else:
            differences = []

        return {
            "experiment1": experiment1,
            "experiment2": experiment2,
            "entity_id": entity_id,
            "entity_type": entity_type,
            "timepoints": timepoints,
            "trajectory1": trajectory1,
            "trajectory2": trajectory2,
            "differences": differences,
        }
