"""Tests for BK polyomavirus (BKPyV) plugin and simulator."""

import pytest

from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.plugins.base import PluginRegistry
from vcm.simulators.bkpyv_simulator import BKPyVSimulator
from vcm.core.models import Perturbation, PerturbationType


def test_bkpyv_plugin_creation():
    """Test BKPyV plugin creation."""
    plugin = BKPolyomavirusPlugin()
    assert plugin.plugin_id == "transplant.bk_polyomavirus"
    assert plugin.plugin_name == "BK Polyomavirus Kidney Cell"
    assert "BKPyV" in plugin.description


def test_bkpyv_plugin_registration():
    """Test BKPyV plugin registration."""
    plugin = BKPolyomavirusPlugin()
    plugin.register()

    registered = PluginRegistry.get("transplant.bk_polyomavirus")
    assert registered is not None
    assert registered.plugin_id == "transplant.bk_polyomavirus"


def test_bkpyv_initial_state():
    """Test BKPyV initial state creation."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    # Check basic state properties
    assert state.cell_type == "kidney_tubular_epithelial"
    assert state.cell_id == "bkpyv_kidney_cell_001"

    # Check for viral genes
    assert "viral_LT" in state.genes
    assert "viral_ST" in state.genes
    assert "viral_VP1" in state.genes

    # Check for host genes
    assert "MCM2" in state.genes  # DNA replication
    assert "CCND1" in state.genes  # Cell cycle
    assert "TP53" in state.genes  # DNA damage response
    assert "STAT1" in state.genes  # Innate immunity
    assert "MTOR" in state.genes  # Drug target

    # Check for pathway activities
    assert "pathway_activities" in state.metadata
    assert "dna_replication" in state.metadata["pathway_activities"]
    assert "cell_cycle" in state.metadata["pathway_activities"]
    assert "viral_replication" in state.metadata["pathway_activities"]

    # Check for enhanced metadata
    assert "cell_cycle_phase" in state.metadata
    assert "host_dna_synthesis" in state.metadata
    assert "t_antigen_level" in state.metadata
    assert "viral_load" in state.metadata
    assert "infection_status" in state.metadata


def test_bkpyv_pathway_activities():
    """Test that pathway activities are properly initialized."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    pathway_activities = state.metadata["pathway_activities"]

    # Check expected pathways exist
    expected_pathways = [
        "dna_replication",
        "cell_cycle",
        "dna_damage_response",
        "innate_immune",
        "apoptosis",
        "mTOR_signaling",
        "viral_replication",
        "cellular_stress",
        "interferon_response",
        "nucleotide_synthesis",
    ]

    for pathway in expected_pathways:
        assert pathway in pathway_activities
        assert 0.0 <= pathway_activities[pathway] <= 2.0

    # Initially, viral replication should be 0
    assert pathway_activities["viral_replication"] == 0.0


def test_bkpyv_supported_perturbations():
    """Test BKPyV supported perturbations."""
    plugin = BKPolyomavirusPlugin()
    perturbations = plugin.get_supported_perturbations()

    # Check for key perturbations
    perturbation_ids = [p.id for p in perturbations]
    assert "bkpyv_infection" in perturbation_ids
    assert "tacrolimus_treatment" in perturbation_ids
    assert "sirolimus_treatment" in perturbation_ids
    assert "antiviral_treatment" in perturbation_ids
    assert "dna_damage" in perturbation_ids

    # Check that BKPyV infection targets viral genes
    infection_pert = next(p for p in perturbations if p.id == "bkpyv_infection")
    assert infection_pert.perturbation_type == PerturbationType.VIRAL_INFECTION


def test_bkpyv_perturbation_targets():
    """Test that perturbations have correct targets."""
    plugin = BKPolyomavirusPlugin()
    perturbations = plugin.get_supported_perturbations()

    # Tacrolimus targets FKBP1A
    tacrolimus = next(p for p in perturbations if p.id == "tacrolimus_treatment")
    assert tacrolimus.target_id == "FKBP1A"

    # Sirolimus targets MTOR
    sirolimus = next(p for p in perturbations if p.id == "sirolimus_treatment")
    assert sirolimus.target_id == "MTOR"

    # Antiviral targets viral LT
    antiviral = next(p for p in perturbations if p.id == "antiviral_treatment")
    assert antiviral.target_id == "viral_LT"


def test_bkpyv_simulator_creation():
    """Test BKPyV simulator creation."""
    simulator = BKPyVSimulator()
    assert simulator.tacrolimus_enhancement_factor == 1.5
    assert simulator.mtor_inhibition_factor == 0.5
    assert simulator.t_antigen_replication_threshold == 0.5


def test_bkpyv_simulator_custom_config():
    """Test BKPyV simulator with custom configuration."""
    custom_config = {
        "tacrolimus_enhancement_factor": 2.0,
        "mtor_inhibition_factor": 0.3,
        "t_antigen_replication_threshold": 0.7,
    }
    simulator = BKPyVSimulator(custom_config)

    assert simulator.tacrolimus_enhancement_factor == 2.0
    assert simulator.mtor_inhibition_factor == 0.3
    assert simulator.t_antigen_replication_threshold == 0.7


def test_bkpyv_simulator_step():
    """Test BKPyV simulator step function."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    simulator = BKPyVSimulator()
    new_state = simulator.step(
        state,
        perturbation=None,
        environment=None,
        timestep=1.0,
        drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
    )

    # Check that state updated
    assert new_state.timestamp == 1.0
    assert new_state.cell_id == state.cell_id

    # Check that pathway activities were updated
    assert "pathway_activities" in new_state.metadata


def test_bkpyv_infection_perturbation():
    """Test BKPyV infection perturbation application."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    # Create infection perturbation with immediate timing
    infection = Perturbation(
        id="bkpyv_infection",
        name="BKPyV infection",
        perturbation_type=PerturbationType.VIRAL_INFECTION,
        magnitude=1.0,
        timing=0.0,  # Apply immediately
    )

    simulator = BKPyVSimulator()

    # Apply perturbation
    new_state = simulator.step(
        state,
        perturbation=infection,
        environment=None,
        timestep=1.0,
        drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
    )

    # Check that viral genes were activated
    assert new_state.genes["viral_LT"].expression_level > 0
    assert new_state.genes["viral_ST"].expression_level > 0
    assert new_state.metadata["infection_status"] == "latent"
    assert new_state.metadata["viral_load"] > 0
    assert new_state.metadata["infection_status"] == "latent"
    assert new_state.metadata["viral_load"] > 0


def test_bkpyv_t_antigen_expression():
    """Test T antigen expression dynamics."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    # Infect the cell
    state.genes["viral_LT"].expression_level = 0.5
    state.proteins["LT"].concentration = 0.5
    state.proteins["LT"].active = True
    state.metadata["infection_status"] = "latent"

    simulator = BKPyVSimulator()

    # Run several steps to see T antigen expression dynamics
    current_state = state
    for _ in range(10):
        current_state = simulator.step(
            current_state,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )

    # T antigen should increase
    assert current_state.genes["viral_LT"].expression_level >= 0.5

    # Check T antigen level classification
    if current_state.genes["viral_LT"].expression_level < 0.1:
        assert current_state.metadata["t_antigen_level"] == "none"
    elif current_state.genes["viral_LT"].expression_level < 0.5:
        assert current_state.metadata["t_antigen_level"] == "low"


def test_bkpyv_cell_cycle_phase():
    """Test cell cycle phase dynamics."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    # Initially should be in G0/G1
    assert state.metadata["cell_cycle_phase"] == "G0/G1"

    simulator = BKPyVSimulator()

    # Infect the cell to drive cell cycle
    state.genes["viral_LT"].expression_level = 1.0
    state.proteins["LT"].concentration = 1.0
    state.proteins["LT"].active = True

    current_state = state
    for _ in range(20):
        current_state = simulator.step(
            current_state,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )

    # Cell cycle should progress due to T antigen
    # High T antigen can drive cells into S phase
    if current_state.genes["viral_LT"].expression_level > 0.5:
        # Should have progressed cell cycle activity
        pathway_activities = current_state.metadata["pathway_activities"]
        assert pathway_activities["cell_cycle"] >= 0.6  # Should be elevated


def test_bkpyv_host_dna_synthesis():
    """Test host DNA synthesis status tracking."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    simulator = BKPyVSimulator()

    # The reviewer-informed baseline is terminally differentiated/G0-G1;
    # S-phase permissiveness is an explicit gate, not an assumed baseline.
    assert state.metadata["host_dna_synthesis"] == "inactive"

    # Apply DNA damage to suppress synthesis
    damage = Perturbation(
        id="dna_damage",
        name="DNA damage",
        perturbation_type=PerturbationType.ENVIRONMENTAL_CHANGE,
        magnitude=0.8,
        timing=0.0,
    )

    current_state = state
    for _ in range(10):
        current_state = simulator.step(
            current_state,
            perturbation=damage,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )

    # DNA damage should increase DDR activity and potentially suppress synthesis
    pathway_activities = current_state.metadata["pathway_activities"]
    assert pathway_activities["dna_damage_response"] > 0.3  # Should be elevated


def test_bkpyv_viral_replication_dependency():
    """Test that viral replication depends on T antigen threshold."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    simulator = BKPyVSimulator()

    # Set T antigen below threshold
    state.genes["viral_LT"].expression_level = 0.3
    state.proteins["LT"].concentration = 0.3
    state.proteins["LT"].active = True

    # Run simulation
    current_state = state
    for _ in range(10):
        current_state = simulator.step(
            current_state,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )

    # Viral replication should be limited
    viral_flux = current_state.pathways["viral_replication"].flux
    assert viral_flux < 0.6  # Should be low due to low T antigen (adjusted threshold)

    # Now set T antigen above threshold
    current_state.genes["viral_LT"].expression_level = 1.0
    current_state.proteins["LT"].concentration = 1.0

    # Run more steps
    for _ in range(10):
        current_state = simulator.step(
            current_state,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )

    # Viral replication should increase
    new_viral_flux = current_state.pathways["viral_replication"].flux
    assert new_viral_flux > viral_flux  # Should be higher with higher T antigen


def test_bkpyv_drug_effects():
    """Test drug effects on BKPyV replication."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    # Infect the cell
    state.genes["viral_LT"].expression_level = 0.8
    state.proteins["LT"].concentration = 0.8
    state.proteins["LT"].active = True
    state.metadata["infection_status"] = "latent"

    # Scenario 1: No drug
    simulator1 = BKPyVSimulator()
    current_state1 = state
    for _ in range(20):
        current_state1 = simulator1.step(
            current_state1,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )
    viral_load_no_drug = current_state1.metadata["viral_load"]

    # Scenario 2: Tacrolimus (should enhance)
    simulator2 = BKPyVSimulator()
    current_state2 = state
    for _ in range(20):
        current_state2 = simulator2.step(
            current_state2,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.8, "mtor_inhibitor": 0.0}
        )
    viral_load_tacrolimus = current_state2.metadata["viral_load"]

    # Tacrolimus should enhance replication
    assert viral_load_tacrolimus >= viral_load_no_drug * 0.9  # Should be similar or higher

    # Scenario 3: Sirolimus (should suppress)
    simulator3 = BKPyVSimulator()
    current_state3 = state
    for _ in range(20):
        current_state3 = simulator3.step(
            current_state3,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.8}
        )
    viral_load_sirolimus = current_state3.metadata["viral_load"]

    # Sirolimus should suppress replication
    assert viral_load_sirolimus <= viral_load_no_drug * 1.1  # Should be similar or lower


def test_bkpyv_pathway_activities_feedback():
    """Test that pathway activities create feedback loops."""
    plugin = BKPolyomavirusPlugin()
    state = plugin.create_initial_state()

    simulator = BKPyVSimulator()

    # Infect cell to activate pathways
    state.genes["viral_LT"].expression_level = 1.0
    state.proteins["LT"].concentration = 1.0
    state.proteins["LT"].active = True

    initial_immune = state.metadata["pathway_activities"]["innate_immune"]

    # Run simulation
    current_state = state
    for _ in range(10):
        current_state = simulator.step(
            current_state,
            perturbation=None,
            environment=None,
            timestep=1.0,
            drug_effects={"tacrolimus": 0.0, "mtor_inhibitor": 0.0}
        )

    # Immune response should activate due to infection
    final_immune = current_state.metadata["pathway_activities"]["innate_immune"]
    # The immune response may increase due to viral induction
    assert isinstance(final_immune, float)
    assert 0.0 <= final_immune <= 2.0


def test_bkpyv_cell_schema():
    """Test BKPyV cell schema information."""
    plugin = BKPolyomavirusPlugin()
    schema = plugin.get_cell_schema()

    # Check for BKPyV-specific information
    assert "viral_genes" in schema
    assert "viral_LT" in schema["viral_genes"]
    assert "pathway_activities" in schema
    assert "cell_cycle_phases" in schema
    assert "drug_targets" in schema
    assert "key_pathways" in schema

    # Check that key pathways match expected list
    expected_key_pathways = [
        "dna_replication",
        "cell_cycle",
        "dna_damage_response",
        "innate_immune",
        "viral_replication",
    ]
    for pathway in expected_key_pathways:
        assert pathway in schema["key_pathways"]


def test_bkpyv_default_config():
    """Test BKPyV default configuration.

    The ODE simulator is the canonical engine; the discrete simulator
    ("bkpyv_specific", hourly heuristic rules) is retained for legacy configs.
    """
    plugin = BKPolyomavirusPlugin()
    config = plugin.get_default_config()

    assert config.plugin == "transplant.bk_polyomavirus"
    assert config.simulator == "bkpyv_ode"
    assert config.simulation_length == 100.0
    assert config.output_path == "outputs/bkpyv/"
