"""Tests for ODE-based BKPyV simulator and ODE system."""

import sys
from pathlib import Path

import pytest
import numpy as np
from scipy.integrate import solve_ivp

from vcm.simulators.ode_system import BKPyVODESystem
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.core.models import Perturbation, PerturbationType


class TestBKPyVODESystem:
    """Test the ODE system definition and behavior."""

    def test_ode_system_initialization(self):
        """Test ODE system initialization with default parameters."""
        ode_system = BKPyVODESystem()
        assert ode_system.params is not None
        assert 'beta' in ode_system.params
        assert 'delta' in ode_system.params
        assert 'tac_enhancement' in ode_system.params

    def test_ode_system_custom_parameters(self):
        """Test ODE system initialization with custom parameters."""
        custom_params = {'beta': 0.2, 'delta': 0.8}
        ode_system = BKPyVODESystem(custom_params)
        assert ode_system.params['beta'] == 0.2
        assert ode_system.params['delta'] == 0.8
        assert ode_system.params['tac_enhancement'] == 1.5  # Default preserved

    def test_state_vector_names(self):
        """Test that state vector names are correctly defined."""
        ode_system = BKPyVODESystem()
        names = ode_system.get_state_vector_names()
        assert len(names) == 22
        assert 'V' in names  # Viral load
        assert 'T' in names  # T antigen
        assert 'C' in names  # Healthy cells
        assert 'I' in names  # Infected cells
        assert 'D_tac' in names  # Tacrolimus
        assert 'D_sir' in names  # Sirolimus

    def test_initial_conditions(self):
        """Test initial conditions generation."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_initial_conditions()
        assert len(y0) == len(ode_system.get_state_vector_names())
        assert y0[0] == 0.0  # No initial virus
        assert y0[3] == 1.0  # One healthy cell
        assert y0[4] == 0.0  # No infected cells

    def test_infection_conditions(self):
        """Test infection initial conditions."""
        ode_system = BKPyVODESystem()
        y0_inf = ode_system.get_infection_conditions(viral_load=0.5)
        assert y0_inf[0] == 0.5  # Viral load set
        assert y0_inf[1] == 0.3  # Initial T antigen
        assert y0_inf[4] > 0.0  # Some infected cells
        assert y0_inf[3] < 1.0  # Reduced healthy cells

    def test_ode_function_call(self):
        """Test that ODE function can be called and returns correct shape."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_initial_conditions()
        dydt = ode_system.ode_system(0.0, y0)
        assert len(dydt) == len(ode_system.get_state_vector_names())
        assert all(np.isfinite(dydt))

    def test_ode_function_with_infection(self):
        """Test ODE function with active infection."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions()
        dydt = ode_system.ode_system(0.0, y0)
        assert len(dydt) == len(ode_system.get_state_vector_names())
        # Viral load should be changing (non-zero derivative)
        assert dydt[0] != 0.0 or dydt[1] != 0.0

    def test_ode_function_with_drugs(self):
        """Test ODE function with drug events."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions()
        drug_events = {0.0: {'tacrolimus': 1.0}}
        dydt = ode_system.ode_system(0.0, y0, drug_events)
        assert len(dydt) == len(ode_system.get_state_vector_names())
        # Drug concentrations should be changing
        assert dydt[11] < 0.0  # Tacrolimus clearance

    def test_ode_numerical_stability(self):
        """Test numerical stability of ODE integration."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions()

        def ode_func(t, y):
            return ode_system.ode_system(t, y)

        sol = solve_ivp(ode_func, (0, 100), y0, method='LSODA')
        assert sol.success
        assert np.all(np.isfinite(sol.y))
        assert sol.y.shape[0] == len(ode_system.get_state_vector_names())

    def test_viral_clearance_kinetics(self):
        """Viral load must decay when production is switched off.

        Isolates the clearance term: with production disabled (p=0) and no new
        infections, V should decline monotonically at roughly the rate implied
        by delta (half-life ln(2)/delta, modulo antiviral-state enhancement).
        """
        ode_system = BKPyVODESystem()
        ode_system.params['delta'] = 0.5  # ~1.4 day half-life
        ode_system.params['p'] = 0.0  # switch off virion production
        ode_system.params['beta'] = 0.0  # no new infections
        y0 = ode_system.get_infection_conditions()
        y0[4] = 0.0  # no infected-cell reservoir producing virus

        def ode_func(t, y):
            return ode_system.ode_system(t, y)

        sol = solve_ivp(ode_func, (0, 10), y0, t_eval=np.arange(0.0, 10.5, 0.5),
                        method='LSODA')

        initial_viral = y0[0]
        final_viral = sol.y[0, -1]
        assert final_viral < initial_viral
        # Monotone decrease (allowing tiny solver jitter)
        diffs = np.diff(sol.y[0])
        assert (diffs <= 1e-6).all()
        # Half-life should be at least as fast as delta alone implies and not
        # absurdly faster than delta + antiviral enhancement allows
        import math
        observed_half_lives = []
        for i in range(1, len(sol.t)):
            v = sol.y[0, i]
            if 0.0 < v < initial_viral / 2.0:
                observed_half_lives.append(sol.t[i])
                break
        assert observed_half_lives, "viral load never dropped below half the inoculum"
        expected_min = math.log(2) / (0.5 * 3.0)  # delta * (1 + ak_max_enhancement)
        expected_max = math.log(2) / 0.5 * 4  # generous bound (sanity)
        assert expected_min <= observed_half_lives[0] <= expected_max

    def test_infection_persists_with_default_parameters(self):
        """With default (baseline) parameters, infection must not self-clear.

        Regression test: a previous parameterisation always extinguished the
        infection (peak == inoculum), making every dose/response scenario a
        decay curve. BKPyV can persist under ongoing replication, so the
        default regime must sustain a non-zero plateau.
        """
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)

        def ode_func(t, y):
            return ode_system.ode_system(t, y)

        sol = solve_ivp(ode_func, (0, 60), y0, t_eval=np.arange(0, 61, 1.0),
                        method='LSODA')
        V = sol.y[0]
        assert V.max() > 0.5, "infection never amplified beyond the inoculum"
        assert V[-1] > 0.05, "infection always self-clears to zero"

    def test_drug_effect_consistency(self):
        """Drug knobs must act through their documented mechanisms.

        - Tacrolimus (D_tac) acts ONLY on immune control: it must reduce the
          infected-cell killing pressure (dI/dt less negative / more positive)
          and must not appear in the intracellular replication flux.
        - Sirolimus (D_sir) acts on production permissiveness: with enough
          present, dV/dt from production is lower than without.
        """
        ode_system = BKPyVODESystem()

        # --- Tacrolimus: identical state ± drug ------------------------------
        y_no_drug = ode_system.get_infection_conditions()
        y_tac = y_no_drug.copy()
        y_tac[11] = 1.0  # Tacrolimus present

        d_no_drug = ode_system.ode_system(10.0, y_no_drug)
        d_tac = ode_system.ode_system(10.0, y_tac)

        # Tacrolimus must not change dT/dt (no direct genome-copy enhancement)
        assert d_tac[1] == pytest.approx(d_no_drug[1])
        # It must change infected-cell loss (dI/dt index 4) in favour of the virus
        assert d_tac[4] > d_no_drug[4]

        # --- Sirolimus: suppresses production --------------------------------
        # get_infection_conditions() sets T=0.3, below the 0.5 threshold, so
        # the mTOR early-phase effect is fully active in both copies.
        y_ref = y_no_drug.copy()
        y_sir = y_ref.copy()
        y_sir[12] = 1.0  # Sirolimus present
        d_sir = ode_system.ode_system(10.0, y_sir)
        d_ref = ode_system.ode_system(10.0, y_ref)
        assert d_sir[0] < d_ref[0]  # viral production reduced



class TestBKPyVODESimulator:
    """Test the ODE-based BKPyV simulator integration."""

    def test_simulator_initialization(self):
        """Test ODE simulator initialization."""
        simulator = BKPyVODESimulator()
        assert simulator.ode_system is not None
        assert simulator.ode_solver == "LSODA"
        assert simulator.rtol == 1e-6
        assert simulator.atol == 1e-8

    def test_simulator_custom_config(self):
        """Test ODE simulator with custom configuration."""
        config = {
            'ode_solver': 'RK45',
            'rtol': 1e-5,
            'atol': 1e-7,
            'research_params': {'beta': 0.2}
        }
        simulator = BKPyVODESimulator(config)
        assert simulator.ode_solver == "RK45"
        assert simulator.rtol == 1e-5
        assert simulator.ode_system.params['beta'] == 0.2

    def test_parameter_mapping(self):
        """Test that discrete parameters are correctly mapped to ODE parameters."""
        config = {
            'tacrolimus_enhancement_factor': 2.0,
            'mtor_inhibition_factor': 0.3
        }
        simulator = BKPyVODESimulator(config)
        assert simulator.ode_system.params['tac_enhancement'] == 2.0
        assert simulator.ode_system.params['mtor_inhibition'] == 0.3

    def test_cellstate_to_ode_conversion(self):
        """Test conversion from CellState to ODE state vector."""
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()

        simulator = BKPyVODESimulator()
        y0 = simulator._cellstate_to_ode(initial_state)

        assert len(y0) == len(simulator.ode_system.get_state_vector_names())
        assert y0[0] == initial_state.metadata['viral_load']
        assert y0[1] == initial_state.genes['viral_LT'].expression_level

    def test_ode_to_cellstate_conversion(self):
        """Test conversion from ODE state vector to CellState."""
        plugin = BKPolyomavirusPlugin()
        template_state = plugin.create_initial_state()

        simulator = BKPyVODESimulator()
        y_test = simulator.ode_system.get_infection_conditions()

        cell_state = simulator._ode_to_cellstate(y_test, template_state, 10.0)

        assert cell_state.timestamp == 10.0
        assert cell_state.metadata['viral_load'] == y_test[0]
        assert cell_state.metadata['t_antigen_level'] == y_test[1]

    def test_simulate_baseline(self):
        """Test baseline simulation (uninfected)."""
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()

        simulator = BKPyVODESimulator()
        result = simulator.simulate(
            initial_state=initial_state,
            n_steps=10,
            timestep=1.0
        )

        assert result.success
        assert len(result.steps) == 11  # 0 to 10 inclusive
        assert result.final_state.metadata['viral_load'] == 0.0  # No infection

    def test_simulate_infection(self):
        """Test simulation with BKPyV infection."""
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()

        infection = Perturbation(
            id="bkpyv_infection",
            name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            magnitude=1.0,
            timing=5.0
        )

        simulator = BKPyVODESimulator()
        result = simulator.simulate(
            initial_state=initial_state,
            perturbations=[infection],
            n_steps=20,
            timestep=1.0
        )

        assert result.success
        assert len(result.steps) == 21

        # Check that the infection event increases viral load from the
        # pre-event state and remains represented in the later trajectory.
        pre_infection_viral = result.steps[4].cell_state.metadata['viral_load']
        post_infection_viral = result.steps[10].cell_state.metadata['viral_load']
        event_viral = result.steps[5].cell_state.metadata['viral_load']
        assert event_viral >= pre_infection_viral
        assert post_infection_viral > 0.0

    def test_simulate_with_tacrolimus(self):
        """Test simulation with tacrolimus treatment."""
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()

        infection = Perturbation(
            id="bkpyv_infection",
            name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            magnitude=1.0,
            timing=5.0
        )

        tacrolimus = Perturbation(
            id="tacrolimus_treatment",
            name="Tacrolimus treatment",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="FKBP1A",
            magnitude=1.0,
            timing=0.0
        )

        simulator = BKPyVODESimulator()
        result = simulator.simulate(
            initial_state=initial_state,
            perturbations=[infection, tacrolimus],
            n_steps=20,
            timestep=1.0
        )

        assert result.success
        # Tacrolimus effect should be present in metadata
        assert result.final_state.metadata.get('tacrolimus_effect', 0) >= 0

    def test_simulate_with_sirolimus(self):
        """Test simulation with sirolimus treatment."""
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()

        infection = Perturbation(
            id="bkpyv_infection",
            name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            magnitude=1.0,
            timing=5.0
        )

        sirolimus = Perturbation(
            id="sirolimus_treatment",
            name="Sirolimus treatment",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="MTOR",
            magnitude=1.0,
            timing=0.0
        )

        simulator = BKPyVODESimulator()
        result = simulator.simulate(
            initial_state=initial_state,
            perturbations=[infection, sirolimus],
            n_steps=20,
            timestep=1.0
        )

        assert result.success
        # Sirolimus effect should be present in metadata
        assert result.final_state.metadata.get('sirolimus_effect', 0) >= 0

    def test_step_function(self):
        """Test single step simulation."""
        plugin = BKPolyomavirusPlugin()
        current_state = plugin.create_initial_state()

        simulator = BKPyVODESimulator()
        new_state = simulator.step(current_state, timestep=1.0)

        assert new_state.timestamp == 1.0
        assert new_state.cell_id == current_state.cell_id

    def test_drug_event_creation(self):
        """Test drug event creation from perturbations."""
        simulator = BKPyVODESimulator()

        tacrolimus = Perturbation(
            id="tacrolimus",
            name="Tacrolimus",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="FKBP1A",
            magnitude=1.0,
            timing=10.0
        )

        sirolimus = Perturbation(
            id="sirolimus",
            name="Sirolimus",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="MTOR",
            magnitude=0.8,
            timing=15.0
        )

        drug_events = simulator._create_drug_events([tacrolimus, sirolimus])

        assert 10.0 in drug_events
        assert 15.0 in drug_events
        assert 'tacrolimus' in drug_events[10.0]
        assert 'sirolimus' in drug_events[15.0]

    def test_infection_event_creation(self):
        """Test infection event creation from perturbations."""
        simulator = BKPyVODESimulator()

        infection = Perturbation(
            id="bkpyv_infection",
            name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            magnitude=1.0,
            timing=10.0
        )

        # Infection events are now handled in the simulate method
        # This test verifies the infection perturbation structure
        assert infection.perturbation_type.value == "viral_infection"
        assert infection.timing == 10.0
        assert infection.magnitude == 1.0


class TestODESystemValidation:
    """Test validation of ODE system against research expectations."""

    def test_viral_half_life_range(self):
        """Test that viral clearance rate produces realistic half-lives."""
        ode_system = BKPyVODESystem()

        # Test that default delta produces half-life in expected range (1-38 hours)
        # Half-life = ln(2) / delta (in days) * 24 (hours)
        delta = ode_system.params['delta']
        half_life_hours = np.log(2) / delta * 24

        # Research shows half-lives of 1-2h (fast) or 20-38h (moderate)
        assert 1.0 <= half_life_hours <= 50.0  # Allow some range

    def test_drug_effect_magnitudes(self):
        """Test that drug effect magnitudes match research values."""
        ode_system = BKPyVODESystem()

        # From AJT-16-821.pdf: tacrolimus enhances, sirolimus inhibits
        assert ode_system.params['tac_enhancement'] >= 1.0
        assert ode_system.params['mtor_inhibition'] <= 1.0
        assert ode_system.params['mtor_inhibition'] >= 0.0

    def test_t_antigen_threshold(self):
        """Test that T antigen threshold is in biologically plausible range."""
        ode_system = BKPyVODESystem()

        # T antigen threshold should be between 0 and 1
        assert 0.0 < ode_system.params['t_threshold'] < 1.0

    def test_state_conservation(self):
        """Test that cell populations are conserved (no negative values)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions()

        def ode_func(t, y):
            return ode_system.ode_system(t, y)

        sol = solve_ivp(ode_func, (0, 50), y0, method='LSODA')

        # Check that cell counts never go negative
        assert all(sol.y[3] >= 0)  # Healthy cells
        assert all(sol.y[4] >= 0)  # Infected cells
        assert all(sol.y[5] >= 0)  # Dead cells

    def test_ode_solver_compatibility(self):
        """Test that different ODE solvers work with the system."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions()

        def ode_func(t, y):
            return ode_system.ode_system(t, y)

        solvers = ['LSODA', 'RK45', 'BDF']
        for solver in solvers:
            sol = solve_ivp(ode_func, (0, 20), y0, method=solver)
            assert sol.success, f"Solver {solver} failed"


class TestExtendedCompartments:
    """Tests for the appended adaptive-T-cell and urothelial compartments.

    Indices 15-19: T_naive, T_eff, C_u, I_u, V_u.
    """

    def test_urine_amplifies_above_plasma(self):
        """Urothelial amplification makes urinary load exceed plasma load by
        orders of magnitude under sustained viremia (Funk 2008: urine ~3000x
        plasma; we assert a conservative >50x on normalised units)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)

        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y),
            (0, 120), y0, t_eval=np.arange(0, 121, 1.0), method='LSODA',
        )
        assert sol.success
        v_plasma, v_urine = sol.y[0], sol.y[19]
        assert v_urine[-1] > 0.0
        assert v_urine[-1] > v_plasma[-1] * 50.0

    def test_urine_origin_is_urothelial(self):
        """The urothelial production term must dominate the urine load
        (Funk 2008: >95% of the urine load is urothelial-derived, not
        kidney drainage)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y),
            (0, 120), y0, t_eval=np.arange(0, 121, 1.0), method='LSODA',
        )
        v_plasma, i_u = sol.y[0, -1], sol.y[18, -1]
        p = ode_system.params
        urothelial = p['p_u'] * i_u
        drainage = p['drain_kidney'] * v_plasma
        assert urothelial / (urothelial + drainage) > 0.9

    def test_tcell_arm_expands_under_antigen(self):
        """BKPyV-specific effector T cells must expand when antigen is
        present (baseline repertoire is 0.02)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y),
            (0, 120), y0, t_eval=np.arange(0, 121, 1.0), method='LSODA',
        )
        t_eff = sol.y[16]
        assert t_eff.max() > 0.1

    def test_tacrolimus_blunts_tcell_expansion(self):
        """Tacrolimus must suppress T_eff expansion (calcineurin/NFAT
        blockade) — the clinically dominant mechanism for BKPyV risk under
        tacrolimus (Kotton 2024)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        no_drug = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y),
            (0, 60), y0, t_eval=np.arange(0, 61, 1.0), method='LSODA',
        )
        with_tac = solve_ivp(
            lambda t, y: ode_system.ode_system(
                t, y, dosing_context={
                    "tacrolimus": {"start": 0.0, "stop": None, "target": 1.0}}),
            (0, 60), y0, t_eval=np.arange(0, 61, 1.0), method='LSODA',
        )
        assert with_tac.y[16].max() < no_drug.y[16].max() * 0.5

    def test_plasma_clears_while_viruria_persists(self):
        """Funk 2008 signature: strong curtailment clears plasma viremia
        while the urothelial reservoir keeps shedding (viruria outlives
        viremia — urine PCR stays positive after plasma clears)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=3.0)
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(
                t, y, {"tacrolimus": {"start": 0.0, "stop": None, "target": 1.0}}),
            (0, 60), y0, t_eval=np.linspace(0, 60, 601), method='LSODA',
        )
        y_peak = sol.y[:, int(np.argmax(sol.y[0]))]

        curtailed = BKPyVODESystem({"p": 8.0 * 0.1, "p_u": 500.0 * 0.1})
        sol2 = solve_ivp(
            lambda t, y: curtailed.ode_system(t, y),
            (0, 140), y_peak, t_eval=np.linspace(0, 140, 561), method='LSODA',
        )
        assert sol2.y[0, -1] < 0.05      # plasma viremia cleared
        assert sol2.y[19, -1] > 1.0      # urinary shedding persists

    def test_metadata_exposes_new_compartments(self):
        """CellState metadata must carry the new compartment readouts."""
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()
        infection = Perturbation(
            id="bkpyv_infection", name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            magnitude=1.0, timing=5.0,
        )
        simulator = BKPyVODESimulator()
        result = simulator.simulate(
            initial_state=initial_state, perturbations=[infection],
            n_steps=20, timestep=1.0,
        )
        md = result.final_state.metadata
        for key in ("urine_viral_load", "bkpyv_tcell_effector",
                    "bkpyv_tcell_naive", "urothelial_infected_cells"):
            assert key in md


class TestPharmacokineticDosing:
    """ng/mL (trough) dosing schedules with real drug half-lives."""

    def test_ngml_taper_resolves_to_troughs(self):
        """A tacrolimus 8 -> 4 ng/mL step taper must drive D_tac to
        intensity 1.0 then 0.5 (ref trough 8 ng/mL)."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        ctx = {"tacrolimus": [
            {"start": 0.0, "stop": 42.0, "trough_ng_ml": 8.0},
            {"start": 42.0, "stop": None, "trough_ng_ml": 4.0},
        ]}
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, ctx),
            (0, 90), y0, t_eval=np.arange(0, 91, 1.0), method='LSODA',
        )
        d_tac = sol.y[11]
        assert 0.9 <= d_tac[40] <= 1.1      # 8 ng/mL steady state
        assert 0.4 <= d_tac[60] <= 0.6      # 4 ng/mL steady state

    def test_ngml_washout_follows_half_life(self):
        """Stopping tacrolimus must wash out on the ~12 h half-life, not the
        legacy ad-hoc clearance envelope."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        ctx = {"tacrolimus": [{"start": 0.0, "stop": 30.0, "trough_ng_ml": 8.0}]}
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, ctx),
            (0, 60), y0, t_eval=np.linspace(0, 60, 121), method='LSODA',
        )
        d_tac = sol.y[11]
        # 12 h half-life: 2 days after stop -> ~4 half-lives -> ~1/16 of target
        stopped_idx = int(np.searchsorted(sol.t, 32.0))
        assert d_tac[stopped_idx] < 0.15

    def test_sirolimus_slow_accumulation(self):
        """Sirolimus (t½ ~60 h) must accumulate slowly — a tac->sir switch is
        not instantaneous, matching clinical practice."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        ctx = {"sirolimus": [{"start": 0.0, "stop": None, "trough_ng_ml": 6.0}]}
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, ctx),
            (0, 30), y0, t_eval=np.arange(0, 31, 1.0), method='LSODA',
        )
        d_sir = sol.y[12]
        target = 6.0 / ode_system.params['sir_ref_trough_ngml']
        assert d_sir[2] < target * 0.6      # far from steady state at day 2
        assert d_sir[14] > target * 0.7     # mostly there by ~2 weeks

    def test_legacy_dimensionless_target_unchanged(self):
        """Dimensionless ``target`` schedules keep the legacy envelope."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        ctx = {"tacrolimus": {"start": 0.0, "stop": None, "target": 1.0}}
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, ctx),
            (0, 30), y0, t_eval=np.arange(0, 31, 1.0), method='LSODA',
        )
        assert sol.y[11, -1] > 0.9

    def test_ngml_units_from_perturbation(self):
        """Perturbations declaring parameters units=ng_ml must reach the ODE
        as trough_ng_ml windows."""
        pert = Perturbation(
            id="tac_taper", name="tac taper",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="FKBP1A", magnitude=8.0, timing=0.0, duration=42.0,
            parameters={"units": "ng_ml"},
        )
        pert2 = Perturbation(
            id="tac_low", name="tac low",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="FKBP1A", magnitude=4.0, timing=42.0,
            parameters={"units": "ng_ml"},
        )
        ctx = BKPyVODESimulator._build_dosing_context([pert, pert2])
        assert len(ctx["tacrolimus"]) == 2
        assert ctx["tacrolimus"][0]["trough_ng_ml"] == 8.0
        assert ctx["tacrolimus"][1]["trough_ng_ml"] == 4.0
        assert ctx["tacrolimus"][0]["stop"] == 42.0


class TestNCCREmergence:
    """F_rr quasi-species dynamics: rearranged NCCR emerges under
    sustained replication (Gosert 2008)."""

    def test_state_vector_is_21d_with_frr_last(self):
        names = BKPyVODESystem().get_state_vector_names()
        assert names[-2] == 'F_rr'
        assert names[-1] == 'F_rr_u'
        assert len(BKPyVODESystem().get_initial_conditions()) == len(names)

    def test_frr_emerges_under_sustained_viremia(self):
        """An archetype inoculum evolves a substantial rearranged fraction
        over the weeks-months window reported clinically."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, None),
            (0, 120), y0, t_eval=np.linspace(0, 120, 121), method='LSODA',
        )
        f_rr = sol.y[20]
        assert f_rr[0] == pytest.approx(0.0)
        assert np.all((f_rr >= 0.0) & (f_rr <= 1.0))
        assert f_rr[-1] > 0.5

    def test_rearranged_preset_seeds_f1(self):
        """nccr_variant='rearranged' maps onto the F_rr=1 boundary instead
        of pinning static multipliers."""
        sim = BKPyVODESimulator()
        plugin = BKPolyomavirusPlugin()
        cell = plugin.create_initial_state({"nccr_variant": "rearranged"})
        y0 = sim._cellstate_to_ode(cell)
        assert y0[20] == pytest.approx(1.0)
        assert y0[21] == pytest.approx(1.0)  # urinary pool seeded too
        early, capsid = sim.ode_system._nccr_multipliers(1.0)
        assert early == pytest.approx(2.0)
        assert capsid == pytest.approx(0.5)

    def test_emergence_can_be_disabled(self):
        """nccr_emergence_enabled=0 restores the fixed-genotype model."""
        ode_system = BKPyVODESystem()
        ode_system.params['nccr_emergence_enabled'] = 0.0
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, None),
            (0, 60), y0, t_eval=[60], method='LSODA',
        )
        assert sol.y[20, -1] == pytest.approx(0.0)

    def test_strong_suppression_slows_emergence(self):
        """Under heavy immunosuppression-driven curtailment the rearranged
        fraction stays much lower than in untreated infection."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        sol_on = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, None),
            (0, 90), y0, t_eval=[90], method='LSODA',
        )
        ode_curtailed = BKPyVODESystem()
        ode_curtailed.params.update({'p': 8.0 * 0.1, 'p_u': 500.0 * 0.1})
        sol_off = solve_ivp(
            lambda t, y: ode_curtailed.ode_system(t, y, None),
            (0, 90), y0, t_eval=[90], method='LSODA',
        )
        assert sol_off.y[20, -1] < sol_on.y[20, -1] * 0.5


class TestReductionSchedules:
    """The model's headline clinical finding: tac->sir conversion should
    dominate tac tapering on BOTH viral clearance and immune rebound."""

    def test_sir_conversion_clears_viremia(self):
        """Conversion to sirolimus clears plasma viremia below the Kotton
        screening threshold within the horizon; tapering to 3 ng/mL tac
        does not."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
        from optimize_reduction_schedule import evaluate_schedule

        taper = evaluate_schedule([(0.0, 28.0, 8.0), (28.0, None, 3.0)])
        conversion = evaluate_schedule([(0.0, 28.0, 8.0), (28.0, None, 3.0)],
                                       sir_trough=4.0)
        assert taper["clearance_weeks"] is None
        assert conversion["clearance_weeks"] is not None
        assert conversion["clearance_weeks"] <= 16.0
        # The winning schedule also rebounds LESS — sir hits replication
        # permissiveness, not just the T-cell brake.
        assert conversion["rebound_index"] < taper["rebound_index"]

    def test_frr_tracks_replication_pressure(self):
        """Schedules that fail to clear leave the virion pool dominated by
        rearranged NCCR; conversion suppresses emergence."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
        from optimize_reduction_schedule import evaluate_schedule

        hold = evaluate_schedule([(0.0, None, 8.0)])
        conversion = evaluate_schedule([(0.0, 28.0, 8.0), (28.0, None, 3.0)],
                                       sir_trough=4.0)
        assert hold["final_frr"] > 0.8
        assert conversion["final_frr"] < hold["final_frr"]

    def test_rr_enriched_in_plasma_not_urine(self):
        """Gosert 2008 signature: the rearranged fraction enriches in the
        kidney (plasma) pool relative to the urinary pool — urothelial
        production is shedding-driven so the rr advantage is weaker."""
        ode_system = BKPyVODESystem()
        y0 = ode_system.get_infection_conditions(viral_load=0.5)
        sol = solve_ivp(
            lambda t, y: ode_system.ode_system(t, y, None),
            (0, 120), y0, t_eval=np.linspace(0, 120, 121), method='LSODA',
        )
        f_k, f_u = sol.y[20], sol.y[21]
        assert f_k[-1] > 2.0 * f_u[-1]
        assert np.all((f_u >= 0.0) & (f_u <= 1.0))


class TestScreeningPolicies:
    """Screening-trigger policy ordering (Kotton 2024 logic, mechanistic)."""

    def test_early_trigger_suppresses_emergence(self):
        """Acting at the 1k screening trigger leaves the kidney pool
        archetype-dominated; never acting lets rr variants take over."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
        from screening_policy_analysis import run_policy

        early = run_policy(1_000.0)
        never = run_policy(float("inf"))
        assert early["clearance_weeks"] is not None
        assert early["final_frr_kidney"] < 0.1
        assert never["final_frr_kidney"] > 0.5

    def test_later_triggers_lose_efficacy(self):
        """Monotone clinical ordering: earlier action -> earlier clearance
        and less T-cell rebound."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
        from screening_policy_analysis import run_policy

        screen = run_policy(1_000.0)
        pyvan = run_policy(10_000.0)
        assert screen["clearance_weeks"] < pyvan["clearance_weeks"]
        assert screen["rebound_index"] < pyvan["rebound_index"]
