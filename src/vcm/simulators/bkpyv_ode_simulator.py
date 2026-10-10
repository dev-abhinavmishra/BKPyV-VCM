"""ODE-based BKPyV simulator with scipy integration.

This simulator provides a true ODE-based implementation of BKPyV dynamics
using scipy's numerical ODE solvers, replacing the discrete-time multiplicative
update rules with continuous-time differential equations.

Research Grounding:
- Drug mechanisms: AJT-16-821.pdf (tacrolimus activates via FKBP-12, sirolimus inhibits via mTOR)
- Viral half-lives: abstract-16323135.txt (1-2h fast, 20-38h moderate clearance)
- Single-cell host biology: Weissbach et al., J Virol 2024;98(12):e01382-24
"""

import copy
import numpy as np
from typing import Any, Dict, Optional
from scipy.integrate import solve_ivp

from vcm.core.models import CellState, Environment, Perturbation, SimulationResult, SimulationStep
from vcm.simulators.base import BaseSimulator
from vcm.simulators.ode_system import BKPyVODESystem


class BKPyVODESimulator(BaseSimulator):
    """ODE-based BK polyomavirus simulator using scipy integration.
    
    This simulator uses a system of 20 coupled ODEs to model:
    - Viral dynamics (viral load, T antigen, viral gene expression)
    - Host cell states (healthy, infected, dead cells)
    - Cell cycle progression and DNA synthesis
    - Immune response (effector cells, interferon, antiviral state)
    - BKPyV-specific adaptive T cells (naive → effector; calcineurin-sensitive)
    - Urinary/urothelial compartment with kidney↔bladder cross-feeding (Funk 2008)
    - Drug pharmacodynamics (tacrolimus, sirolimus)
    - Pathway activities (DNA replication, innate immune)
    
    The ODE parameters are mapped from the existing research-validated discrete
    parameters to ensure continuity with existing validation results.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize the ODE-based BKPyV simulator.

        Args:
            config: Configuration dictionary with parameters:
                - ode_solver: Solver method ('LSODA', 'RK45', 'BDF', etc.)
                - rtol: Relative tolerance for adaptive solver (default: 1e-6)
                - atol: Absolute tolerance for adaptive solver (default: 1e-8)
                - max_step: Maximum step size, days (default: 1.0)
                - research_params: Dictionary of research-validated parameters

        Note: YAML floats like ``1e-6`` are parsed as strings by PyYAML
        (YAML 1.1); numeric-looking strings are coerced defensively here.
        """
        super().__init__(config)

        # ODE solver configuration (coerce YAML-1.1 string forms like "1e-6")
        def _num(key: str, default: float) -> float:
            raw = config.get(key, default) if config else default
            try:
                return float(raw)
            except (TypeError, ValueError):
                return default

        self.ode_solver = str(config.get("ode_solver", "LSODA")) if config else "LSODA"
        self.rtol = _num("rtol", 1e-6)
        self.atol = _num("atol", 1e-8)
        self.max_step = _num("max_step", 1.0)
        
        # Initialize ODE system with research parameters
        research_params = config.get("research_params", {}) if config else {}
        self.ode_system = BKPyVODESystem(research_params)
        
        # Map existing discrete parameters to ODE parameters
        self._map_discrete_to_ode_parameters()
    
    def _map_discrete_to_ode_parameters(self):
        """Map legacy discrete-simulator parameter names to ODE parameters.

        Note on ``tacrolimus_enhancement_factor``: in the legacy discrete
        simulator this directly multiplied viral production. In the ODE model
        tacrolimus instead weakens immune control, so values > 1.0 are mapped
        onto ``innate_immune_suppression`` as ``min(1, factor - 1)`` (the
        default 1.5 maps to the default 0.5). A directly supplied
        ``innate_immune_suppression`` (or ``innate_immune_suppression_factor``)
        value always wins.
        """
        # Get the discrete parameters that were used in the original simulator
        discrete_params = {
            'tacrolimus_enhancement_factor': 1.5,    # Retained for compatibility; now affects immune control
            'mtor_inhibition_factor': 0.5,          # From AJT-16-821.pdf
            't_antigen_replication_threshold': 0.5,  # Research-validated
            'cell_cycle_s_phase_bonus': 2.0,        # Research-validated
            'dna_replication_coupling': 0.8,       # Research-validated
            'innate_immune_suppression_factor': 0.5,  # Clinical-validated
            'dna_damage_response_enhancement': 1.3,   # Clinical-validated
            'translation_enhancement_factor': 2.0,     # Single-cell validated
            'mitochondrial_function_importance': 0.8,   # Single-cell validated
            'protein_degradation_inhibition': 0.3,      # Pathway-validated
            'nccr_early_expression_multiplier': 1.0,
            'nccr_capsid_expression_multiplier': 1.0,
        }
        
        # Override config if provided
        if self.config:
            for key, value in discrete_params.items():
                if key in self.config:
                    discrete_params[key] = self.config[key]
        
        # Map to ODE parameters
        ode_param_mapping = {
            'tac_enhancement': discrete_params['tacrolimus_enhancement_factor'],
            'mtor_inhibition': discrete_params['mtor_inhibition_factor'],
            't_threshold': discrete_params['t_antigen_replication_threshold'],
            's_phase_bonus': discrete_params['cell_cycle_s_phase_bonus'],
            'dna_coupling': discrete_params['dna_replication_coupling'],
            'innate_immune_suppression': discrete_params['innate_immune_suppression_factor'],
            'ddr_enhancement': discrete_params['dna_damage_response_enhancement'],
            'translation_enhancement': discrete_params['translation_enhancement_factor'],
            'mitochondrial_importance': discrete_params['mitochondrial_function_importance'],
            'protein_degradation_inhibition': discrete_params['protein_degradation_inhibition'],
            'nccr_early_expression_multiplier': discrete_params['nccr_early_expression_multiplier'],
            'nccr_capsid_expression_multiplier': discrete_params['nccr_capsid_expression_multiplier'],
        }

        # Route the legacy direct-replication multiplier onto the immune-control
        # strength UNLESS an explicit immune-suppression value was configured.
        if self.config:
            has_explicit_immune = (
                "innate_immune_suppression" in self.config
                or "innate_immune_suppression_factor" in self.config
            )
            if "tacrolimus_enhancement_factor" in self.config and not has_explicit_immune:
                factor = float(self.config["tacrolimus_enhancement_factor"])
                ode_param_mapping["innate_immune_suppression"] = min(1.0, max(0.0, factor - 1.0))

        # Update ODE system parameters (legacy mapped names first)
        self.ode_system.params.update(ode_param_mapping)

        # Then allow direct ODE-parameter overrides by name (e.g. {'p': 12.0}).
        # Solver settings and registry-only names are not ODE params; skip them.
        solver_keys = {"ode_solver", "rtol", "atol", "max_step", "research_params"}
        if self.config:
            legacy_names = set(self._legacy_discrete_param_names())
            for key, value in self.config.items():
                if key in solver_keys or key in legacy_names:
                    continue
                if key in self.ode_system.params:
                    try:
                        self.ode_system.params[key] = float(value)
                    except (TypeError, ValueError):
                        # Keep YAML-string values (e.g. a mis-written "1e-6")
                        # from poisoning the numeric parameter dictionary.
                        pass

    @staticmethod
    def _legacy_discrete_param_names():
        return (
            'tacrolimus_enhancement_factor', 'mtor_inhibition_factor',
            't_antigen_replication_threshold', 'cell_cycle_s_phase_bonus',
            'dna_replication_coupling', 'innate_immune_suppression_factor',
            'dna_damage_response_enhancement', 'translation_enhancement_factor',
            'mitochondrial_function_importance', 'protein_degradation_inhibition',
        )
    
    def simulate(
        self,
        initial_state: CellState,
        perturbation: Optional[Perturbation] = None,
        perturbations: Optional[list[Perturbation]] = None,
        environment: Optional[Environment] = None,
        n_steps: int = 100,
        timestep: float = 1.0,
    ) -> SimulationResult:
        """Run an ODE-based BKPyV simulation.
        
        Args:
            initial_state: Starting cell state
            perturbation: Optional single perturbation (for backward compatibility)
            perturbations: Optional list of perturbations (new preferred interface)
            environment: Environmental conditions
            n_steps: Number of simulation steps
            timestep: Time step size
        
        Returns:
            SimulationResult with full trajectory
        """
        if environment is None:
            environment = Environment()
        
        # Combine perturbations
        all_perturbations = []
        if perturbation:
            all_perturbations.append(perturbation)
        if perturbations:
            all_perturbations.extend(perturbations)
        
        result = SimulationResult(
            experiment_id=f"bkpyv_ode_{initial_state.cell_id}",
            simulator_type="bkpyv_ode",
            plugin=initial_state.cell_type,
            config_id="default",
        )
        
        # Convert CellState to ODE initial conditions
        y0 = self._cellstate_to_ode(initial_state)

        # Continuous drug administration is part of the ODE right-hand side
        # (dosing context); only impulsive events (infection inoculum) are
        # applied discretely at segment boundaries so the adaptive solver never
        # sees a discontinuous state mutation mid-step.
        dosing_context = self._build_dosing_context(all_perturbations)
        end_time = n_steps * timestep
        t_eval = np.linspace(0.0, end_time, n_steps + 1)
        event_times = sorted({
            float(pert.timing)
            for pert in all_perturbations
            if pert.perturbation_type.value == "viral_infection"
            and pert.timing is not None
            and 0.0 < pert.timing < end_time
        })
        boundaries = [0.0, *event_times, end_time]

        def apply_events(state_vector: np.ndarray, event_time: float) -> np.ndarray:
            updated = state_vector.copy()
            for pert in all_perturbations:
                if pert.perturbation_type.value != "viral_infection":
                    continue
                if pert.timing is None or abs(float(pert.timing) - event_time) > 1e-9:
                    continue
                updated[0] += pert.magnitude
                updated[2] += pert.magnitude * 0.5
                updated[1] = max(updated[1], 0.3)
                updated[4] += pert.magnitude * 0.1
                updated[3] = max(0.0, updated[3] - pert.magnitude * 0.1)
            return updated

        y_current = apply_events(y0, 0.0)
        solved_times: list[float] = []
        solved_states: list[np.ndarray] = []
        for start, stop in zip(boundaries[:-1], boundaries[1:]):
            segment_times = t_eval[(t_eval >= start) & (t_eval <= stop)]
            if len(segment_times) == 0 or segment_times[-1] < stop:
                segment_times = np.append(segment_times, stop)
            sol = self._solve_ode_system(y_current, (start, stop), segment_times, dosing_context)
            for segment_time, segment_state in zip(sol.t, sol.y.T):
                if solved_times and abs(float(segment_time) - solved_times[-1]) < 1e-9:
                    solved_states[-1] = segment_state.copy()
                else:
                    solved_times.append(float(segment_time))
                    solved_states.append(segment_state.copy())
            y_current = apply_events(sol.y[:, -1], stop)
            if solved_times and abs(solved_times[-1] - stop) < 1e-9:
                solved_states[-1] = y_current.copy()
        
        # Convert ODE solution back to CellState steps
        for i, (t, y) in enumerate(zip(solved_times, solved_states)):
            cell_state = self._ode_to_cellstate(y, initial_state, t)
            
            # Determine active perturbations at this time
            active_perturbations = []
            for pert in all_perturbations:
                if pert.timing is not None:
                    if (pert.timing <= t <= (pert.timing + pert.duration) if pert.duration else pert.timing <= t):
                        active_perturbations.append(pert)
            
            step = SimulationStep(
                step_number=i,
                timestamp=t,
                cell_state=cell_state,
                applied_perturbations=active_perturbations,
                environment=environment,
            )
            result.steps.append(step)
        
        result.final_state = result.steps[-1].cell_state
        return result
    
    def step(
        self,
        current_state: CellState,
        perturbation: Optional[Perturbation] = None,
        environment: Optional[Environment] = None,
        timestep: float = 1.0,
    ) -> CellState:
        """Perform a single ODE simulation step.
        
        Args:
            current_state: Current cell state
            perturbation: Optional perturbation to apply
            environment: Environmental conditions
            timestep: Time step size
        
        Returns:
            Updated cell state after one step
        """
        if environment is None:
            environment = Environment()
        
        # Convert to ODE state
        y0 = self._cellstate_to_ode(current_state)

        # Setup time span for single step
        t_span = (current_state.timestamp, current_state.timestamp + timestep)
        t_eval = [current_state.timestamp + timestep]

        # Dosing context for this step (continuous infusion semantics).
        perts = [p for p in (perturbation,) if p is not None]
        dosing_context = self._build_dosing_context(perts)

        # Solve ODE for single step
        sol = self._solve_ode_system(y0, t_span, t_eval, dosing_context)

        # Convert back to CellState
        new_state = self._ode_to_cellstate(sol.y.T[0], current_state, current_state.timestamp + timestep)
        return new_state

    @staticmethod
    def _build_dosing_context(perturbations: list[Perturbation]) -> Dict[str, Any]:
        """Build the continuous-dosing schedule used by the ODE right-hand side.

        Each drug entry is ``{'start': day, 'stop': day|None, 'target': magnitude}``.
        Multiple treatments of the same drug merge to the earliest start,
        latest stop, and strongest target. ``stop=None`` means ongoing dosing.
        """
        context: Dict[str, Any] = {}
        for pert in perturbations:
            if pert.perturbation_type.value != "drug_treatment":
                continue
            drug = {"FKBP1A": "tacrolimus", "MTOR": "sirolimus"}.get(pert.target_id)
            if drug is None:
                continue
            start = float(pert.timing) if pert.timing is not None else 0.0
            stop = None if pert.duration is None else start + float(pert.duration)
            magnitude = float(pert.magnitude)
            if drug not in context:
                context[drug] = {"start": start, "stop": stop, "target": magnitude, "continuous": pert.duration is None}
            else:
                entry = context[drug]
                entry["start"] = min(entry["start"], start)
                if entry["continuous"] or pert.duration is None:
                    entry["stop"] = None
                    entry["continuous"] = True
                elif entry["stop"] is not None and stop is not None:
                    entry["stop"] = max(entry["stop"], stop)
                entry["target"] = max(entry["target"], magnitude)
        for entry in context.values():
            if entry.pop("continuous", False):
                entry["stop"] = None
        return context

    def _solve_ode_system(self, y0: np.ndarray, t_span: tuple,
                         t_eval: np.ndarray, dosing_context: Dict[str, Any]):
        """Solve the ODE system using scipy.

        Args:
            y0: Initial conditions
            t_span: Time span (t_start, t_end) in days
            t_eval: Time points for evaluation
            dosing_context: Continuous drug administration schedule; see
                :meth:`BKPyVODESystem.ode_system`.

        Returns:
            scipy ODE solution object
        """
        def ode_func(t, y):
            return self.ode_system.ode_system(t, y, dosing_context)

        return solve_ivp(
            ode_func,
            t_span,
            y0,
            t_eval=t_eval,
            method=self.ode_solver,
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step,
        )
    
    def _cellstate_to_ode(self, cell_state: CellState) -> np.ndarray:
        """Convert CellState to ODE state vector.
        
        Args:
            cell_state: CellState object
        
        Returns:
            ODE state vector
        """
        # Extract viral load from metadata
        viral_load = cell_state.metadata.get("viral_load", 0.0)
        nccr_variant = cell_state.metadata.get("nccr_variant", "archetype")
        
        # Extract T antigen level
        t_antigen = cell_state.genes.get("viral_LT", Gene(id="viral_LT", name="viral_LT")).expression_level
        
        # Extract viral gene expression
        viral_genes = ["viral_LT", "viral_ST", "viral_VP1", "viral_VP2", "viral_VP3"]
        viral_gene_expr = sum(cell_state.genes.get(g, Gene(id=g, name=g)).expression_level for g in viral_genes)
        
        # Extract cell counts (simplified)
        infection_status = cell_state.metadata.get("infection_status", "uninfected")
        if infection_status == "uninfected":
            healthy_cells = 1.0
            infected_cells = 0.0
        else:
            healthy_cells = 0.7
            infected_cells = 0.3
        
        # Extract cell cycle phase
        cell_cycle_phase = cell_state.metadata.get("cell_cycle_phase", "G0/G1")
        cc_map = {"G0/G1": 0.1, "S": 0.5, "G2/M": 0.9}
        cc = cc_map.get(cell_cycle_phase, 0.1)
        
        # Extract DNA synthesis status
        dna_status = cell_state.metadata.get("dna_synthesis_status", "active")
        dna_map = {"active": 0.8, "inactive": 0.2, "suppressed": 0.1}
        dna = dna_map.get(dna_status, 0.5)
        
        # Extract pathway activities
        pathway_activities = cell_state.metadata.get("pathway_activities", {})
        p_rep = pathway_activities.get("dna_replication", 0.7)
        p_immune = pathway_activities.get("innate_immune", 0.5)
        
        # Build state vector
        y0 = self.ode_system.get_initial_conditions()
        y0[0] = viral_load  # V
        y0[1] = t_antigen   # T
        y0[2] = viral_gene_expr * 0.5  # G_v
        y0[3] = healthy_cells  # C
        y0[4] = infected_cells  # I
        y0[6] = cc  # CC
        y0[7] = dna  # DNA
        y0[13] = p_rep  # P_rep
        y0[14] = p_immune  # P_immune
        # Variant presets are overridable for calibration/sensitivity work.
        # This makes it possible to distinguish a biological scenario from a
        # fitted coefficient instead of silently overwriting the coefficient.
        if "nccr_early_expression_multiplier" not in self.config:
            self.ode_system.params["nccr_early_expression_multiplier"] = 2.0 if nccr_variant == "rearranged" else 1.0
        if "nccr_capsid_expression_multiplier" not in self.config:
            self.ode_system.params["nccr_capsid_expression_multiplier"] = 0.5 if nccr_variant == "rearranged" else 1.0
        
        return y0
    
    def _ode_to_cellstate(self, y: np.ndarray, template_state: CellState, timestamp: float) -> CellState:
        """Convert ODE state vector to CellState.
        
        Args:
            y: ODE state vector
            template_state: Template CellState for structure
            timestamp: Current timestamp
        
        Returns:
            CellState object
        """
        # Unpack state vector (15-dim intracellular block + appended
        # T-cell and urothelial compartments; extras defensively defaulted
        # for any legacy 15-vector caller)
        V, T, G_v, C, I, D, CC, DNA, E, IFN, AK, D_tac, D_sir, P_rep, P_immune = y[:15]
        T_naive = float(y[15]) if len(y) > 15 else 0.0
        T_eff = float(y[16]) if len(y) > 16 else 0.0
        C_u = float(y[17]) if len(y) > 17 else 1.0
        I_u = float(y[18]) if len(y) > 18 else 0.0
        V_u = float(y[19]) if len(y) > 19 else 0.0
        
        # Create new state based on template
        new_state = copy.deepcopy(template_state)
        new_state.timestamp = timestamp
        
        # Update viral genes
        new_state.genes["viral_LT"].expression_level = T
        new_state.genes["viral_ST"].expression_level = G_v * 0.6
        capsid_multiplier = self.ode_system.params.get("nccr_capsid_expression_multiplier", 1.0)
        new_state.genes["viral_VP1"].expression_level = G_v * 0.3 * capsid_multiplier
        new_state.genes["viral_VP2"].expression_level = G_v * 0.05
        new_state.genes["viral_VP3"].expression_level = G_v * 0.05
        
        # Update viral proteins
        new_state.proteins["LT"].concentration = T
        new_state.proteins["LT"].active = T > 0.1
        new_state.proteins["ST"].concentration = G_v * 0.6
        new_state.proteins["VP1"].concentration = G_v * 0.3 * capsid_multiplier
        
        # Update metadata (clamp solver overshoot: a state can dip slightly
        # below zero without biological meaning)
        V = max(0.0, V)
        new_state.metadata["viral_load"] = V
        new_state.metadata["t_antigen_level"] = T
        new_state.metadata["healthy_cell_count"] = C
        new_state.metadata["infected_cell_count"] = I
        new_state.metadata["dead_cell_count"] = D
        
        # Update cell cycle phase
        if CC < 0.3:
            new_state.metadata["cell_cycle_phase"] = "G0/G1"
        elif CC < 0.7:
            new_state.metadata["cell_cycle_phase"] = "S"
        else:
            new_state.metadata["cell_cycle_phase"] = "G2/M"
        
        # Update DNA synthesis status
        if DNA > 0.6:
            new_state.metadata["dna_synthesis_status"] = "active"
        elif DNA > 0.3:
            new_state.metadata["dna_synthesis_status"] = "inactive"
        else:
            new_state.metadata["dna_synthesis_status"] = "suppressed"
        
        # Update infection status
        if V > 0.01:
            if T > self.ode_system.params['t_threshold']:
                new_state.metadata["infection_status"] = "active_lytic"
            else:
                new_state.metadata["infection_status"] = "latent"
        else:
            new_state.metadata["infection_status"] = "uninfected"
        
        # Update pathway activities
        pathway_activities = new_state.metadata.get("pathway_activities", {})
        pathway_activities["dna_replication"] = P_rep
        pathway_activities["innate_immune"] = P_immune
        pathway_activities["interferon_response"] = min(1.0, IFN * 2.0)
        pathway_activities["cell_cycle"] = CC
        pathway_activities["mTOR_signaling"] = max(0.0, 1.0 - D_sir * 0.7)
        new_state.metadata["pathway_activities"] = pathway_activities
        
        # Update drug effects metadata
        new_state.metadata["tacrolimus_effect"] = D_tac
        new_state.metadata["sirolimus_effect"] = D_sir
        new_state.metadata["nccr_variant"] = template_state.metadata.get("nccr_variant", "archetype")

        # Adaptive T-cell arm and urinary compartment (Funk 2008)
        new_state.metadata["bkpyv_tcell_naive"] = T_naive
        new_state.metadata["bkpyv_tcell_effector"] = T_eff
        new_state.metadata["urothelial_healthy_cells"] = C_u
        new_state.metadata["urothelial_infected_cells"] = I_u
        new_state.metadata["urine_viral_load"] = max(0.0, V_u)
        new_state.metadata["intracellular_replication_flux"] = float(max(0.0, V * (T / (T + 0.5))))
        new_state.metadata["viral_production_rate"] = float(max(0.0, P_rep * T))
        new_state.metadata["immune_control_index"] = float(max(0.0, min(1.0, P_immune)))
        
        # Update viral replication pathway flux
        new_state.pathways["viral_replication"].flux = V * (T / (T + 0.5))
        
        # Update state vector
        new_state.update_state_vector()
        
        return new_state
    
    def _create_drug_events(self, perturbations: list[Perturbation]) -> Dict[str, Any]:
        """Create drug event dictionary from perturbations.
        
        Args:
            perturbations: List of perturbations
        
        Returns:
            Dictionary mapping time to drug effects
        """
        drug_events = {}
        for pert in perturbations:
            if pert.perturbation_type.value == "drug_treatment" and pert.timing is not None:
                if pert.target_id == "FKBP1A":
                    drug_events[pert.timing] = {'tacrolimus': pert.magnitude}
                elif pert.target_id == "MTOR":
                    drug_events[pert.timing] = {'sirolimus': pert.magnitude}
        return drug_events


# Need to import Gene for the _cellstate_to_ode method
from vcm.core.models import Gene
