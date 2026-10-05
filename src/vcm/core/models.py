"""Core domain models for the Virtual Cell Model platform."""

from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, model_validator


class PerturbationType(str, Enum):
    """Types of perturbations that can be applied to a cell."""

    GENE_KNOCKOUT = "gene_knockout"
    GENE_OVEREXPRESSION = "gene_overexpression"
    PROTEIN_INHIBITION = "protein_inhibition"
    METABOLITE_ADDITION = "metabolite_addition"
    METABOLITE_DEPLETION = "metabolite_depletion"
    ENVIRONMENTAL_CHANGE = "environmental_change"
    VIRAL_INFECTION = "viral_infection"
    DRUG_TREATMENT = "drug_treatment"
    RADIATION = "radiation"
    CUSTOM = "custom"


class Gene(BaseModel):
    """Represents a gene in the cell model."""

    id: str
    name: str
    expression_level: float = Field(default=1.0, ge=0.0, description="Expression level (arbitrary units)")
    sequence: Optional[str] = Field(default=None, description="DNA/RNA sequence")
    function: Optional[str] = Field(default=None, description="Biological function")
    is_essential: bool = Field(default=False, description="Whether gene is essential for survival")


class Protein(BaseModel):
    """Represents a protein in the cell model."""

    id: str
    name: str
    concentration: float = Field(default=1.0, ge=0.0, description="Protein concentration")
    gene_id: Optional[str] = Field(default=None, description="Associated gene ID")
    function: Optional[str] = Field(default=None, description="Biological function")
    is_enzyme: bool = Field(default=False, description="Whether protein is an enzyme")
    active: bool = Field(default=True, description="Whether protein is active")


class Metabolite(BaseModel):
    """Represents a metabolite in the cell model."""

    id: str
    name: str
    concentration: float = Field(default=1.0, ge=0.0, description="Metabolite concentration")
    compartment: str = Field(default="cytoplasm", description="Cellular compartment")
    is_essential: bool = Field(default=False, description="Whether metabolite is essential")


class Pathway(BaseModel):
    """Represents a biological pathway."""

    id: str
    name: str
    genes: List[str] = Field(default_factory=list, description="Gene IDs in pathway")
    proteins: List[str] = Field(default_factory=list, description="Protein IDs in pathway")
    metabolites: List[str] = Field(default_factory=list, description="Metabolite IDs in pathway")
    flux: float = Field(default=1.0, description="Pathway flux/activity level")


class Environment(BaseModel):
    """Represents the external environment of the cell."""

    temperature: float = Field(default=37.0, description="Temperature in Celsius")
    ph: float = Field(default=7.4, description="pH level")
    oxygen_level: float = Field(default=21.0, description="Oxygen percentage")
    glucose_concentration: float = Field(default=5.0, description="Glucose concentration (mM)")
    nutrient_availability: Dict[str, float] = Field(
        default_factory=dict, description="Availability of various nutrients"
    )
    stress_factors: Dict[str, float] = Field(
        default_factory=dict, description="Various stress factors (e.g., oxidative stress)"
    )


class Perturbation(BaseModel):
    """Represents a perturbation to be applied to the cell."""

    id: str
    name: Optional[str] = Field(default=None, description="Human-readable name (defaults to id)")
    perturbation_type: PerturbationType
    target_id: Optional[str] = Field(default=None, description="Target gene/protein/metabolite ID")
    magnitude: float = Field(default=1.0, description="Magnitude of perturbation")
    duration: Optional[float] = Field(default=None, description="Duration in time units")
    timing: Optional[float] = Field(default=0.0, description="When to apply perturbation")
    parameters: Dict[str, Any] = Field(default_factory=dict, description="Additional parameters")

    @model_validator(mode="after")
    def _default_name(self) -> "Perturbation":
        if self.name is None or self.name == "":
            self.name = self.id
        return self


class CellState(BaseModel):
    """Represents the state of a cell at a given time point.

    This is the central data structure that captures the complete cellular state
    including genes, proteins, metabolites, and pathways.
    """

    cell_id: str
    cell_type: str
    timestamp: float = Field(default=0.0, description="Simulation time")
    genes: Dict[str, Gene] = Field(default_factory=dict, description="Gene registry")
    proteins: Dict[str, Protein] = Field(default_factory=dict, description="Protein registry")
    metabolites: Dict[str, Metabolite] = Field(default_factory=dict, description="Metabolite registry")
    pathways: Dict[str, Pathway] = Field(default_factory=dict, description="Pathway registry")
    state_vector: Optional[List[float]] = Field(
        default=None, description="Vector representation of state for ML models"
    )
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Additional metadata")

    # BKPyV-specific metadata properties for ergonomic access
    @property
    def viral_load(self) -> float:
        """Get viral load from metadata (0.0-1.0 scale)."""
        return self.metadata.get("viral_load", 0.0)

    @property
    def infection_status(self) -> str:
        """Get infection status: 'uninfected', 'latent', or 'active_lytic'."""
        return self.metadata.get("infection_status", "uninfected")

    @property
    def t_antigen_level(self) -> str:
        """Get T antigen level: 'none', 'low', 'medium', or 'high'."""
        return self.metadata.get("t_antigen_level", "none")

    @property
    def viral_replication_phase(self) -> str:
        """Get viral replication phase: 'none', 'early', or 'late'."""
        return self.metadata.get("viral_replication_phase", "none")

    @property
    def cell_cycle_phase(self) -> str:
        """Get cell cycle phase: 'G0/G1', 'S', or 'G2/M'."""
        return self.metadata.get("cell_cycle_phase", "G0/G1")

    @property
    def host_dna_synthesis(self) -> str:
        """Get host DNA synthesis status: 'active', 'inactive', or 'suppressed'."""
        return self.metadata.get("host_dna_synthesis", "inactive")

    @property
    def mitochondrial_stress(self) -> float:
        """Get mitochondrial stress level (0.0-1.0)."""
        return self.metadata.get("mitochondrial_stress", 0.0)

    @property
    def antigen_presentation(self) -> float:
        """Get antigen presentation level (0.0-1.0, 1.0 = normal)."""
        return self.metadata.get("antigen_presentation", 1.0)

    @property
    def innate_immune_suppression(self) -> float:
        """Get innate immune suppression level (0.0-1.0)."""
        return self.metadata.get("innate_immune_suppression", 0.0)

    @property
    def time_since_infection(self) -> float:
        """Get time since infection in hours."""
        return self.metadata.get("time_since_infection", 0.0)

    @property
    def drug_effects(self) -> Dict[str, float]:
        """Get current drug effect levels."""
        return self.metadata.get("drug_effects", {"tacrolimus": 0.0, "sirolimus": 0.0, "everolimus": 0.0})

    @property
    def drug_exposure_duration(self) -> Dict[str, float]:
        """Get drug exposure duration in hours."""
        return self.metadata.get("drug_exposure_duration", {"tacrolimus": 0.0, "sirolimus": 0.0})

    @property
    def pathway_activities(self) -> Dict[str, float]:
        """Get pathway activity levels."""
        return self.metadata.get("pathway_activities", {})

    def get_state_vector(self) -> List[float]:
        """Convert state to a vector representation."""
        if self.state_vector is not None:
            return self.state_vector

        # Simple concatenation of all quantities
        vector: List[float] = []
        vector.extend([g.expression_level for g in self.genes.values()])
        vector.extend([p.concentration for p in self.proteins.values()])
        vector.extend([m.concentration for m in self.metabolites.values()])
        vector.extend([p.flux for p in self.pathways.values()])
        return vector

    def update_state_vector(self) -> None:
        """Update the state vector from current entity values."""
        self.state_vector = self.get_state_vector()


class SimulationStep(BaseModel):
    """Represents a single step in a simulation."""

    step_number: int
    timestamp: float
    cell_state: CellState
    applied_perturbations: List[Perturbation] = Field(default_factory=list)
    environment: Environment = Field(default_factory=Environment)


class SimulationResult(BaseModel):
    """Represents the results of a simulation run."""

    experiment_id: str
    simulator_type: str
    plugin: str
    config_id: str
    steps: List[SimulationStep] = Field(default_factory=list)
    final_state: Optional[CellState] = Field(default=None)
    metadata: Dict[str, Any] = Field(default_factory=dict)
    start_time: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    end_time: Optional[datetime] = Field(default=None)
    success: bool = Field(default=True)
    error_message: Optional[str] = Field(default=None)

    def get_trajectory(self, entity_id: str, entity_type: str) -> List[float]:
        """Extract time series trajectory for a specific entity."""
        trajectory: List[float] = []
        for step in self.steps:
            state = step.cell_state
            if entity_type == "gene" and entity_id in state.genes:
                trajectory.append(state.genes[entity_id].expression_level)
            elif entity_type == "protein" and entity_id in state.proteins:
                trajectory.append(state.proteins[entity_id].concentration)
            elif entity_type == "metabolite" and entity_id in state.metabolites:
                trajectory.append(state.metabolites[entity_id].concentration)
            elif entity_type == "pathway" and entity_id in state.pathways:
                trajectory.append(state.pathways[entity_id].flux)
        return trajectory

    def get_timepoints(self) -> List[float]:
        """Get all timepoints from the simulation."""
        return [step.timestamp for step in self.steps]


class ExperimentConfig(BaseModel):
    """Configuration for an experiment.

    Time unit discipline: ``time_unit`` declares the unit of
    ``simulation_length``, ``timestep``, and perturbation ``timing``/``duration``.
    The canonical unit for BKPyV experiments is "days" (the ODE simulator's
    rate constants are per day). The legacy discrete ``bkpyv_specific``
    simulator historically interpreted timesteps as hours; new work should use
    the ODE simulator with days.
    """

    experiment_id: str
    plugin: str
    simulator: str = Field(default="mechanistic", description="Simulator type to use")
    simulation_length: float = Field(default=100.0, description="Total simulation time (in time_unit)")
    timestep: float = Field(default=1.0, description="Time step size (in time_unit)")
    time_unit: str = Field(default="days", description="Unit of all times in this config: 'days', 'hours', or 'weeks'")
    simulator_parameters: Dict[str, Any] = Field(
        default_factory=dict,
        description="Parameters passed to the simulator constructor (e.g., rate constants)",
    )
    perturbations: List[Perturbation] = Field(default_factory=list)
    environment: Environment = Field(default_factory=Environment)
    output_path: str = Field(default="outputs/", description="Where to save results")
    save_interval: int = Field(default=1, description="Save every N steps")
    seed: Optional[int] = Field(default=None, description="Random seed for reproducibility")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Additional metadata (not passed to the simulator)")

    @model_validator(mode="before")
    @classmethod
    def _merge_legacy_simulator_params(cls, data: Any) -> Any:
        # Older configs used the key ``simulator_params``; accept it as an
        # alias rather than silently dropping parameters.
        if isinstance(data, dict) and isinstance(data.get("simulator_params"), dict):
            legacy = data.pop("simulator_params")
            merged = {**legacy, **(data.get("simulator_parameters") or {})}
            data["simulator_parameters"] = merged
        return data
