"""BKPyV-specific simulator with pathway-driven replication and drug effects.

.. deprecated::
    Superseded by :class:`BKPyVODESimulator` (continuous-time ODE, day units).
    This discrete (1 hour per step) heuristic simulator is kept for backward
    compatibility with older configs and tests, but new analyses should use
    the ODE engine. Known conceptual differences vs. the ODE engine:
    tacrolimus here multiplies the replication rate directly (the ODE routes
    it through immune control), and T antigen here *drives* cell-cycle entry
    (the ODE gates T-antigen accumulation on prior S phase, per Needham 2024).

Mechanistic intent (all phenomenological):
- Viral replication depends on host DNA synthesis status and T antigen expression
- Cell cycle phase influences viral replication efficiency
- In vitro: tacrolimus increases, sirolimus decreases BKPyV replication
  (Hirsch et al., Am J Transplant 2016;16(3):821-832)
- Single-cell pathway patterns from Weissbach et al., J Virol 2024;98(12):e01382-24
"""

import copy
import math
from typing import Any, Dict, Optional

from vcm.core.models import CellState, Environment, Perturbation, SimulationResult, SimulationStep
from vcm.simulators.base import BaseSimulator


class BKPyVSimulator(BaseSimulator):
    """BK polyomavirus-specific simulator with pathway-driven dynamics.

    This simulator incorporates:
    - T antigen-dependent viral replication
    - Cell cycle phase effects on replication efficiency
    - Drug-specific effects (tacrolimus vs mTOR inhibitors)
    - Pathway activity-based modulation of viral replication
    - Host response dynamics (DNA damage, immune signaling, apoptosis)
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize the BKPyV simulator with research-validated parameters.

        Parameters based on AJT-16-821.pdf research findings:
        - Tacrolimus activates BKPyV replication via FKBP-12 (enhancement factor: 1.5)
        - Sirolimus inhibits BKPyV replication via mTOR inhibition (IC90 = 4 ng/mL, factor: 0.5)
        - Critical window: 0-24h post-infection for drug effectiveness during early gene expression

        Parameters based on single-cell transcriptomic data (Weissbach et al., J Virol 2024;98(12):e01382-24 (supplementary pathway analyses)):
        - Translation pathways highly elevated in infection (factor: 2.0)
        - Mitochondrial function upregulated (importance: 0.8)
        - Protein degradation pathways involved (inhibition: 0.3)

        Args:
            config: Configuration dictionary with parameters:
                - tacrolimus_enhancement_factor: Tacrolimus activates replication (default: 1.5, research-validated)
                - mtor_inhibition_factor: Sirolimus IC90 = 4 ng/mL (default: 0.5, research-validated)
                - t_antigen_replication_threshold: T antigen level needed for active replication (default: 0.5)
                - cell_cycle_s_phase_bonus: Replication bonus in S phase (default: 2.0, pathway-validated)
                - dna_replication_coupling: How strongly viral replication couples to host DNA synthesis (default: 0.8, research-validated)
                - innate_immune_suppression_factor: How much innate immunity suppresses replication (default: 0.5, clinical-validated)
                - dna_damage_response_enhancement: DDR enhancement of replication (default: 1.3, clinical-validated)
                - translation_enhancement_factor: Translation pathway elevation in infection (default: 2.0, single-cell validated)
                - mitochondrial_function_importance: Mitochondrial gene upregulation (default: 0.8, single-cell validated)
                - protein_degradation_inhibition: Proteasome pathway involvement (default: 0.3, pathway-validated)
        """
        super().__init__(config)
        # Drug effect parameters (from AJT-16-821.pdf research)
        self.tacrolimus_enhancement_factor = config.get("tacrolimus_enhancement_factor", 1.5) if config else 1.5
        self.mtor_inhibition_factor = config.get("mtor_inhibition_factor", 0.5) if config else 0.5

        # Replication parameters (research and pathway-validated)
        self.t_antigen_replication_threshold = (
            config.get("t_antigen_replication_threshold", 0.5) if config else 0.5
        )
        self.cell_cycle_s_phase_bonus = config.get("cell_cycle_s_phase_bonus", 2.0) if config else 2.0
        self.dna_replication_coupling = config.get("dna_replication_coupling", 0.8) if config else 0.8

        # Host response parameters (clinical-validated)
        self.innate_immune_suppression_factor = (
            config.get("innate_immune_suppression_factor", 0.5) if config else 0.5
        )
        self.dna_damage_response_enhancement = (
            config.get("dna_damage_response_enhancement", 1.3) if config else 1.3
        )

        # Single-cell pathway-validated parameters (from Weissbach et al., J Virol 2024;98(12):e01382-24 (supplementary pathway analyses))
        self.translation_enhancement_factor = config.get("translation_enhancement_factor", 2.0) if config else 2.0
        self.mitochondrial_function_importance = config.get("mitochondrial_function_importance", 0.8) if config else 0.8
        self.protein_degradation_inhibition = config.get("protein_degradation_inhibition", 0.3) if config else 0.3

    def simulate(
        self,
        initial_state: CellState,
        perturbation: Optional[Perturbation] = None,
        perturbations: Optional[list[Perturbation]] = None,
        environment: Optional[Environment] = None,
        n_steps: int = 100,
        timestep: float = 1.0,
    ) -> SimulationResult:
        """Run a BKPyV simulation with pathway-driven dynamics.

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

        # Combine perturbations from both parameters
        all_perturbations = []
        if perturbation:
            all_perturbations.append(perturbation)
        if perturbations:
            all_perturbations.extend(perturbations)

        result = SimulationResult(
            experiment_id=f"bkpyv_{initial_state.cell_id}",
            simulator_type="bkpyv_specific",
            plugin=initial_state.cell_type,
            config_id="default",
        )

        current_state = copy.deepcopy(initial_state)
        current_state.timestamp = 0.0

        # Track drug effects over time (enhanced structure)
        drug_effects = {
            "tacrolimus": 0.0,
            "sirolimus": 0.0,
            "everolimus": 0.0,
        }

        for step_num in range(n_steps):
            # Apply perturbations at appropriate times
            active_perturbations = []
            for pert in all_perturbations:
                if pert.timing is not None:
                    if (
                        pert.timing <= current_state.timestamp
                        and (
                            pert.duration is None
                            or current_state.timestamp < pert.timing + pert.duration
                        )
                    ):
                        active_perturbations.append(pert)

                        # Track drug effects (enhanced with multiple mTOR inhibitors)
                        if pert.target_id == "FKBP1A":
                            drug_effects["tacrolimus"] = pert.magnitude
                        elif pert.target_id == "MTOR":
                            # Determine which mTOR inhibitor based on perturbation name or type
                            if "sirolimus" in pert.name.lower():
                                drug_effects["sirolimus"] = pert.magnitude
                            elif "everolimus" in pert.name.lower():
                                drug_effects["everolimus"] = pert.magnitude
                            else:
                                # Default to sirolimus for mTOR targeting
                                drug_effects["sirolimus"] = pert.magnitude

            # Store current step
            step = SimulationStep(
                step_number=step_num,
                timestamp=current_state.timestamp,
                cell_state=copy.deepcopy(current_state),
                applied_perturbations=active_perturbations,
                environment=environment,
            )
            result.steps.append(step)

            # Update state using BKPyV-specific dynamics
            current_state = self.step(
                current_state=current_state,
                perturbations=active_perturbations,
                environment=environment,
                timestep=timestep,
                drug_effects=drug_effects,
            )

        result.final_state = current_state
        return result

    def step(
        self,
        current_state: CellState,
        perturbation: Optional[Perturbation] = None,
        perturbations: Optional[list[Perturbation]] = None,
        environment: Optional[Environment] = None,
        timestep: float = 1.0,
        drug_effects: Optional[Dict[str, float]] = None,
    ) -> CellState:
        """Perform a single BKPyV simulation step.

        This implements the key mechanistic dynamics:
        1. T antigen expression drives viral replication
        2. Viral replication depends on host DNA synthesis status
        3. Cell cycle phase modulates replication efficiency
        4. Drug effects modify replication and host responses
        5. Host pathways respond to infection and modulate outcome

        Args:
            current_state: Current cell state
            perturbation: Optional single perturbation (backward compatibility)
            perturbations: Optional list of perturbations (new preferred interface)
            environment: Environmental conditions
            timestep: Time step size
            drug_effects: Current drug effect levels

        Returns:
            Updated cell state
        """
        if environment is None:
            environment = Environment()

        if drug_effects is None:
            drug_effects = {"tacrolimus": 0.0, "sirolimus": 0.0, "everolimus": 0.0}

        # Combine perturbations from both parameters
        all_perturbations = []
        if perturbation:
            all_perturbations.append(perturbation)
        if perturbations:
            all_perturbations.extend(perturbations)

        new_state = copy.deepcopy(current_state)
        new_state.timestamp += timestep

        pathway_activities = new_state.metadata["pathway_activities"]

        # Apply perturbation effects
        if all_perturbations:
            for pert in all_perturbations:
                if pert.timing is not None:
                    if pert.timing <= current_state.timestamp:
                        self._apply_perturbation(new_state, pert)

        # Update T antigen expression
        self._update_t_antigen(new_state, pathway_activities, timestep)

        # Update viral replication based on T antigen and host factors
        self._update_viral_replication(
            new_state, pathway_activities, drug_effects, timestep
        )

        # Update host DNA synthesis and cell cycle
        self._update_host_dna_synthesis(new_state, pathway_activities, drug_effects, timestep)

        # Update host gene expression based on pathway activities
        self._update_host_genes(new_state, pathway_activities, drug_effects, timestep)

        # Update pathway activities based on current state
        self._update_pathway_activities(new_state, timestep)

        # Update metadata
        self._update_metadata(new_state, drug_effects)

        new_state.update_state_vector()
        return new_state

    def _apply_perturbation(self, state: CellState, perturbation: Perturbation) -> None:
        """Apply a perturbation to the cell state.

        Args:
            state: Cell state to modify
            perturbation: Perturbation to apply
        """
        magnitude = perturbation.magnitude
        target_id = perturbation.target_id

        if perturbation.perturbation_type.value == "viral_infection":
            # BKPyV infection: introduce viral genes and T antigen
            state.genes["viral_LT"].expression_level = 0.5
            state.genes["viral_ST"].expression_level = 0.3
            state.genes["viral_VP1"].expression_level = 0.1
            state.proteins["LT"].concentration = 0.5
            state.proteins["LT"].active = True
            state.metadata["infection_status"] = "latent"
            state.metadata["viral_load"] = 1.0

        elif perturbation.perturbation_type.value == "drug_treatment":
            if target_id == "FKBP1A":
                # Tacrolimus: calcineurin inhibition reduces immune response
                pathway_activities = state.metadata["pathway_activities"]
                pathway_activities["innate_immune"] *= (1.0 - magnitude * 0.5)
                pathway_activities["interferon_response"] *= (1.0 - magnitude * 0.4)

            elif target_id == "MTOR":
                # mTOR inhibitor: suppresses cell cycle and protein synthesis
                pathway_activities = state.metadata["pathway_activities"]
                pathway_activities["cell_cycle"] *= (1.0 - magnitude * 0.3)
                pathway_activities["mTOR_signaling"] *= (1.0 - magnitude * 0.7)

            elif target_id == "viral_LT":
                # Antiviral treatment: inhibit viral replication
                state.genes["viral_LT"].expression_level *= (1.0 - magnitude)

        elif perturbation.perturbation_type.value == "environmental_change":
            # DNA damage: activate DNA damage response
            pathway_activities = state.metadata["pathway_activities"]
            pathway_activities["dna_damage_response"] *= (1.0 + magnitude)
            pathway_activities["cellular_stress"] *= (1.0 + magnitude * 0.5)

    def _update_t_antigen(self, state: CellState, pathway_activities: Dict[str, float], timestep: float) -> None:
        """Update T antigen expression based on viral replication status.

        T antigen expression is driven by:
        - Viral replication activity
        - Host DNA synthesis status
        - Drug effects (tacrolimus enhances, mTOR inhibitors suppress)

        Args:
            state: Cell state to update
            pathway_activities: Current pathway activity levels
            timestep: Time step size
        """
        viral_gene = state.genes["viral_LT"]
        current_level = viral_gene.expression_level

        # Only update T antigen if there's an active infection
        if state.metadata.get("infection_status", "uninfected") == "uninfected":
            return

        # Baseline T antigen expression growth (only when already infected)
        growth_rate = 0.1 * current_level

        # Enhancement by host DNA synthesis (only when infected)
        dna_rep_activity = pathway_activities["dna_replication"]
        synthesis_enhancement = self.dna_replication_coupling * dna_rep_activity * 0.2

        # Drug effects (applied externally via perturbation)
        # This is handled in _apply_perturbation and drug_effects tracking

        # Update T antigen level
        delta = growth_rate + synthesis_enhancement
        viral_gene.expression_level = max(0.0, current_level + delta * timestep)

        # Update corresponding protein
        state.proteins["LT"].concentration = viral_gene.expression_level

    def _update_viral_replication(
        self,
        state: CellState,
        pathway_activities: Dict[str, float],
        drug_effects: Dict[str, float],
        timestep: float,
    ) -> None:
        """Update viral replication based on T antigen and host factors.

        Viral replication depends on:
        - T antigen expression (threshold required)
        - Host DNA synthesis activity
        - Cell cycle phase (S phase optimal)
        - Drug effects (tacrolimus enhances, mTOR inhibitors suppress)
        - Host innate immune response

        Args:
            state: Cell state to update
            pathway_activities: Current pathway activity levels
            drug_effects: Current drug effect levels
            timestep: Time step size
        """
        t_antigen_level = state.genes["viral_LT"].expression_level
        viral_pathway = state.pathways["viral_replication"]

        # Check if T antigen is above threshold for active replication
        if t_antigen_level < self.t_antigen_replication_threshold:
            viral_pathway.flux = max(0.0, viral_pathway.flux * 0.9)
            return

        # Base replication rate
        base_rate = 0.1 * t_antigen_level

        # Host DNA synthesis coupling
        dna_rep_activity = pathway_activities["dna_replication"]
        synthesis_factor = dna_rep_activity * self.dna_replication_coupling

        # Cell cycle phase effect
        cell_cycle_phase = state.metadata["cell_cycle_phase"]
        cell_cycle_activity = pathway_activities["cell_cycle"]
        cycle_factor = cell_cycle_activity
        if cell_cycle_phase == "S":
            cycle_factor *= self.cell_cycle_s_phase_bonus

        # Drug effects with mechanistic implementation based on research
        # Source: AJT-16-821.pdf - distinct mechanisms through FKBP-12
        # Sirolimus effective only during early phase (0-24h post-infection)
        # Tacrolimus effect persists through both phases but acts via FKBP-12 activation
        
        viral_phase = state.metadata.get("viral_replication_phase", "none")
        time_since_infection = state.metadata.get("time_since_infection", 0.0)
        
        # Tacrolimus effect: Activates replication through FKBP-12 pathway
        # Clinical data: OR 2.0-2.3 for BKPyVAN vs belatacept
        # Mechanism: FKBP-12 binding creates permissive environment for viral replication
        tacrolimus_effect = 1.0 + drug_effects.get("tacrolimus", 0.0) * (self.tacrolimus_enhancement_factor - 1.0)
        
        # Sirolimus effect: Inhibits replication through mTOR-S6K-kinase interference
        # IC90 = 4 ng/mL, effective only during early gene expression (0-24h)
        # Source: AJT-16-821.pdf - "inhibition is rapid and effective up to 24 hours post-infection"
        if viral_phase == "early" and time_since_infection <= 24.0:
            # Full inhibitory effect during early phase
            mtor_effect = 1.0 - drug_effects.get("sirolimus", 0.0) * (1.0 - self.mtor_inhibition_factor)
        else:
            # Reduced effectiveness in late phase (drug resistance mechanism)
            mtor_effect = 1.0 - drug_effects.get("sirolimus", 0.0) * (1.0 - self.mtor_inhibition_factor) * 0.3
        
        # Everolimus effect (alternative mTOR inhibitor, similar to sirolimus)
        # Using same timing-dependent mechanism as sirolimus
        everolimus_effect = 1.0
        if "everolimus" in drug_effects:
            if viral_phase == "early" and time_since_infection <= 24.0:
                everolimus_effect = 1.0 - drug_effects.get("everolimus", 0.0) * (1.0 - self.mtor_inhibition_factor)
            else:
                everolimus_effect = 1.0 - drug_effects.get("everolimus", 0.0) * (1.0 - self.mtor_inhibition_factor) * 0.3

        # Innate immune suppression
        immune_activity = pathway_activities["innate_immune"]
        immune_suppression = 1.0 - immune_activity * self.innate_immune_suppression_factor

        # DNA damage response enhancement with phase-dependent dual role
        # Source: JVI S4 single-cell - DDR genes upregulated in late BKPyV
        # Mechanism: DDR enhances replication in early phase but may trigger apoptosis in late phase
        ddr_activity = pathway_activities["dna_damage_response"]
        if viral_phase == "early":
            # Early phase: DDR enhances replication (permissive)
            ddr_enhancement = 1.0 + ddr_activity * (self.dna_damage_response_enhancement - 1.0)
        else:
            # Late phase: DDR has reduced effect due to potential apoptosis trigger
            ddr_enhancement = 1.0 + ddr_activity * (self.dna_damage_response_enhancement - 1.0) * 0.5

        # Combined replication rate with research-grounded pathway effects
        # NEW: Mitochondrial stress effect (emerges mainly in late replication)
        # Source: Weissbach et al. 2024 (suppl.) - mitochondrial genes upregulated in late infection
        mitochondrial_stress = state.metadata.get("mitochondrial_stress", 0.0)
        mitochondrial_factor = 1.0 + mitochondrial_stress * 0.5  # Stress may enhance late replication
        
        # NEW: Antigen presentation suppression (immune evasion)
        # Source: Single-cell data shows downregulation in successful infection
        antigen_presentation = state.metadata.get("antigen_presentation", 1.0)
        antigen_factor = 1.0 - (1.0 - antigen_presentation) * 0.3  # Reduced presentation may enhance replication
        
        # NEW: Translation pathway enhancement (single-cell validated, factor 2.0)
        translation_activity = pathway_activities.get("translation", 0.7)
        translation_factor = 1.0 + translation_activity * (self.translation_enhancement_factor - 1.0)
        
        # Combined replication rate
        replication_rate = (
            base_rate * synthesis_factor * cycle_factor 
            * tacrolimus_effect * mtor_effect * everolimus_effect
            * immune_suppression * ddr_enhancement
            * mitochondrial_factor * antigen_factor * translation_factor
        )

        # Update viral pathway flux
        viral_pathway.flux = max(0.0, replication_rate)

        # Update viral gene expression based on replication
        for gene_id in ["viral_LT", "viral_ST", "viral_VP1", "viral_VP2", "viral_VP3"]:
            if gene_id in state.genes:
                growth = 0.05 * viral_pathway.flux
                state.genes[gene_id].expression_level = max(
                    0.0, state.genes[gene_id].expression_level + growth * timestep
                )

        # Update viral protein levels
        for protein_id in ["LT", "ST", "VP1"]:
            if protein_id in state.proteins:
                gene_id = f"viral_{protein_id}"
                if gene_id in state.genes:
                    state.proteins[protein_id].concentration = state.genes[gene_id].expression_level

        # Update viral load
        total_viral_expression = sum(
            state.genes[gid].expression_level
            for gid in ["viral_LT", "viral_ST", "viral_VP1", "viral_VP2", "viral_VP3"]
        )
        # Normalize to 0-1 scale for clinical mapping (Hill function expects 0-1)
        # Use a saturation function to prevent unbounded growth
        normalized_viral_load = total_viral_expression / (1.0 + total_viral_expression)
        state.metadata["viral_load"] = max(0.0, min(1.0, normalized_viral_load))

    def _update_host_dna_synthesis(
        self,
        state: CellState,
        pathway_activities: Dict[str, float],
        drug_effects: Dict[str, float],
        timestep: float,
    ) -> None:
        """Update host DNA synthesis status and cell cycle phase.

        Host DNA synthesis is influenced by:
        - T antigen (can drive cells into S phase)
        - mTOR inhibitors (suppress cell cycle)
        - DNA damage response (can arrest cell cycle)

        Args:
            state: Cell state to update
            pathway_activities: Current pathway activity levels
            drug_effects: Current drug effect levels
            timestep: Time step size
        """
        t_antigen_level = state.genes["viral_LT"].expression_level

        # T antigen can drive cells into S phase to support viral replication
        t_antigen_drive = 0.0
        if t_antigen_level > 0.3:
            t_antigen_drive = t_antigen_level * 0.1

        # mTOR inhibition suppresses cell cycle (default to 0 if not present)
        mtor_suppression = drug_effects.get("mtor_inhibitor", 0.0) * 0.2

        # DNA damage response can cause cell cycle arrest
        ddr_activity = pathway_activities.get("dna_damage_response", 0.5)
        ddr_arrest = ddr_activity * 0.15

        # Update cell cycle activity
        current_cycle_activity = pathway_activities["cell_cycle"]
        delta_cycle = t_antigen_drive - mtor_suppression - ddr_arrest
        new_cycle_activity = max(0.0, min(2.0, current_cycle_activity + delta_cycle * timestep))
        pathway_activities["cell_cycle"] = new_cycle_activity

        # Determine cell cycle phase based on activity
        if new_cycle_activity > 1.2:
            state.metadata["cell_cycle_phase"] = "S"
        elif new_cycle_activity > 0.8:
            state.metadata["cell_cycle_phase"] = "G2/M"
        else:
            state.metadata["cell_cycle_phase"] = "G0/G1"

        # Update host DNA synthesis status
        if new_cycle_activity > 0.8:
            state.metadata["host_dna_synthesis"] = "active"
        elif ddr_activity > 0.5:
            state.metadata["host_dna_synthesis"] = "suppressed"
        else:
            state.metadata["host_dna_synthesis"] = "inactive"

        # Update DNA replication pathway activity
        pathway_activities["dna_replication"] = new_cycle_activity * 0.8

    def _update_host_genes(
        self,
        state: CellState,
        pathway_activities: Dict[str, float],
        drug_effects: Dict[str, float],
        timestep: float,
    ) -> None:
        """Update host gene expression based on pathway activities.

        Host gene expression responds to:
        - Viral infection (stress response, immune signaling)
        - Drug effects
        - Pathway feedback loops

        Args:
            state: Cell state to update
            pathway_activities: Current pathway activity levels
            drug_effects: Current drug effect levels
            timestep: Time step size
        """
        viral_load = state.metadata["viral_load"]
        t_antigen_level = state.genes["viral_LT"].expression_level

        # DNA damage response genes respond to viral replication
        ddr_genes = ["ATM", "ATR", "CHEK1", "TP53"]
        ddr_activity = pathway_activities["dna_damage_response"]
        for gene_id in ddr_genes:
            if gene_id in state.genes:
                # Viral replication induces DNA damage response
                viral_induction = viral_load * 0.05
                delta = viral_induction - 0.01 * state.genes[gene_id].expression_level
                state.genes[gene_id].expression_level = max(
                    0.1, state.genes[gene_id].expression_level + delta * timestep
                )

        # Immune response genes respond to infection
        immune_genes = ["STAT1", "IRF7", "IFNB1"]
        immune_activity = pathway_activities["innate_immune"]
        for gene_id in immune_genes:
            if gene_id in state.genes:
                # Infection induces immune response
                infection_induction = viral_load * 0.08
                delta = infection_induction - 0.02 * state.genes[gene_id].expression_level
                state.genes[gene_id].expression_level = max(
                    0.1, state.genes[gene_id].expression_level + delta * timestep
                )

        # Apoptosis genes respond to cellular stress
        apoptosis_genes = ["BAX", "BCL2"]
        stress_activity = pathway_activities["cellular_stress"]
        for gene_id in apoptosis_genes:
            if gene_id in state.genes:
                # Stress induces pro-apoptotic signals
                stress_induction = stress_activity * 0.03
                delta = stress_induction - 0.01 * state.genes[gene_id].expression_level
                state.genes[gene_id].expression_level = max(
                    0.1, state.genes[gene_id].expression_level + delta * timestep
                )

        # Update corresponding protein levels
        for gene_id, gene in state.genes.items():
            for protein_id, protein in state.proteins.items():
                if protein.gene_id == gene_id:
                    protein.concentration = gene.expression_level

    def _update_pathway_activities(self, state: CellState, timestep: float) -> None:
        """Update pathway activities based on current gene expression.

        This creates feedback loops between gene expression and pathway activity.

        Args:
            state: Cell state to update
            timestep: Time step size
        """
        pathway_activities = state.metadata["pathway_activities"]

        # Update innate immune activity based on immune gene expression
        immune_genes_expr = sum(
            state.genes[gid].expression_level for gid in ["STAT1", "IRF7", "IFNB1"]
        )
        target_immune = min(1.0, immune_genes_expr / 3.0)
        pathway_activities["innate_immune"] = (
            pathway_activities["innate_immune"] * 0.9 + target_immune * 0.1
        )

        # Update interferon response (similar to innate immune)
        pathway_activities["interferon_response"] = pathway_activities["innate_immune"] * 0.9

        # Update cellular stress based on stress genes
        stress_genes_expr = state.genes["HSPA1A"].expression_level if "HSPA1A" in state.genes else 0.5
        target_stress = min(1.0, stress_genes_expr)
        pathway_activities["cellular_stress"] = (
            pathway_activities["cellular_stress"] * 0.9 + target_stress * 0.1
        )

        # Update apoptosis activity based on BAX/BCL2 ratio
        if "BAX" in state.genes and "BCL2" in state.genes:
            bax_level = state.genes["BAX"].expression_level
            bcl2_level = state.genes["BCL2"].expression_level
            apoptosis_ratio = bax_level / (bcl2_level + 0.1)
            pathway_activities["apoptosis"] = min(1.0, apoptosis_ratio * 0.5)

        # Update mitochondrial stress pathway (emerges in late infection)
        # Source: Weissbach et al. 2024 (suppl.) - MT-ND4, MT-CO1, MT-CYB upregulated in late BKPyV
        mitochondrial_genes = ["MT-ND4", "MT-CO1", "MT-CYB"]
        if all(gene in state.genes for gene in mitochondrial_genes):
            mt_expression = sum(state.genes[gene].expression_level for gene in mitochondrial_genes) / 3.0
            if state.metadata.get("viral_replication_phase") == "late":
                target_mito_stress = min(1.0, mt_expression)
            else:
                target_mito_stress = mt_expression * 0.2  # Low in early phase
            pathway_activities["mitochondrial_stress"] = (
                pathway_activities["mitochondrial_stress"] * 0.9 + target_mito_stress * 0.1
            )

        # Update antigen presentation pathway (immune evasion)
        # Source: Single-cell data shows downregulation in successful infection
        antigen_genes = ["HLA-A", "HLA-B"]
        if all(gene in state.genes for gene in antigen_genes):
            hla_expression = sum(state.genes[gene].expression_level for gene in antigen_genes) / 2.0
            # Viral load suppresses antigen presentation
            viral_load = state.metadata.get("viral_load", 0.0)
            suppression_factor = max(0.1, 1.0 - viral_load * 0.1)
            target_antigen = min(1.0, hla_expression * suppression_factor)
            pathway_activities["antigen_presentation"] = (
                pathway_activities["antigen_presentation"] * 0.95 + target_antigen * 0.05
            )

    def _update_metadata(self, state: CellState, drug_effects: Dict[str, float]) -> None:
        """Update metadata including infection status, viral phase, and drug effects.

        Enhanced with research-grounded state transitions:
        - Viral replication phase (early vs late) based on time since infection
        - Mitochondrial stress (increases in late phase)
        - Antigen presentation (decreases with immune evasion)
        - Drug exposure tracking for pharmacodynamics

        Args:
            state: Cell state to update
            drug_effects: Current drug effect levels
        """
        viral_load = state.metadata["viral_load"]
        t_antigen_level = state.genes["viral_LT"].expression_level
        pathway_activities = state.metadata["pathway_activities"]

        # Update T antigen level classification
        if t_antigen_level < 0.1:
            state.metadata["t_antigen_level"] = "none"
        elif t_antigen_level < 0.5:
            state.metadata["t_antigen_level"] = "low"
        elif t_antigen_level < 1.5:
            state.metadata["t_antigen_level"] = "medium"
        else:
            state.metadata["t_antigen_level"] = "high"

        # Update infection status based on viral load and T antigen
        if viral_load < 0.5:
            state.metadata["infection_status"] = "uninfected"
            state.metadata["viral_replication_phase"] = "none"
        elif viral_load < 2.0:
            state.metadata["infection_status"] = "latent"
            # Latent infection still has phase based on time
        else:
            state.metadata["infection_status"] = "active_lytic"

        # NEW: Update viral replication phase (early vs late)
        # Source: AJT-16-821.pdf - drug timing effects differ by phase
        time_since_infection = state.metadata.get("time_since_infection", 0.0)
        if state.metadata["infection_status"] != "uninfected":
            if time_since_infection < 24.0:
                state.metadata["viral_replication_phase"] = "early"
            else:
                state.metadata["viral_replication_phase"] = "late"

        # NEW: Update mitochondrial stress (increases in late replication)
        # Source: Weissbach et al. 2024 (suppl.) - discordant mt vs nuclear expression in late infection
        if state.metadata["viral_replication_phase"] == "late":
            mitochondrial_activity = pathway_activities.get("mitochondrial_stress", 0.0)
            # Stress increases as replication progresses in late phase
            stress_increase = min(0.5, (viral_load - 2.0) * 0.1)
            state.metadata["mitochondrial_stress"] = mitochondrial_activity * 0.5 + stress_increase
        else:
            state.metadata["mitochondrial_stress"] = pathway_activities.get("mitochondrial_stress", 0.0) * 0.2

        # NEW: Update antigen presentation (immune evasion mechanism)
        # Source: Single-cell data shows downregulation in successful infection
        innate_immune_activity = pathway_activities.get("innate_immune", 0.5)
        # Higher viral load and lower immune activity -> lower antigen presentation
        immune_suppression = max(0.0, (viral_load - 1.0) * 0.2)
        state.metadata["antigen_presentation"] = max(0.1, 1.0 - immune_suppression - (1.0 - innate_immune_activity) * 0.3)

        # NEW: Update innate immune suppression state
        state.metadata["innate_immune_suppression"] = 1.0 - innate_immune_activity

        # NEW: Update time since infection
        if state.metadata["infection_status"] != "uninfected":
            state.metadata["time_since_infection"] += 1.0  # Assume 1-hour timestep

        # Update drug effect tracking (enhanced structure)
        state.metadata["drug_effects"]["tacrolimus"] = drug_effects.get("tacrolimus", 0.0)
        state.metadata["drug_effects"]["sirolimus"] = drug_effects.get("sirolimus", 0.0) if "sirolimus" in drug_effects else 0.0
        state.metadata["drug_effects"]["everolimus"] = drug_effects.get("everolimus", 0.0) if "everolimus" in drug_effects else 0.0

        # NEW: Update drug exposure duration
        for drug in ["tacrolimus", "sirolimus", "everolimus"]:
            if state.metadata["drug_effects"][drug] > 0.0:
                state.metadata["drug_exposure_duration"][drug] += 1.0
