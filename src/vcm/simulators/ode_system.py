"""ODE system for BKPyV viral-host dynamics.

This module defines the system of ordinary differential equations (ODEs) that model
BK polyomavirus infection in kidney tubular epithelial cells, including drug effects
and host response dynamics.

The ODE system is based on standard virus-host modeling approaches extended with:
- Drug-specific effects (tacrolimus vs sirolimus)
- Pathway activity dynamics
- Cell cycle progression
- Immune response kinetics
- Research-validated parameters from clinical studies

State Variables (15-dimensional vector):
1. V: Viral load (normalised scale; see clinical.viral_load_mapper for the bridge to copies/mL)
2. T: Large T-antigen level (arbitrary units)
3. G_v: Viral gene expression level (arbitrary units)
4. C: Healthy target cells (fraction of carrying capacity)
5. I: Infected cells (fraction)
6. D: Dead/damaged cells (fraction)
7. CC: Cell-cycle permissiveness (bounded 0-1; relaxation toward replication-pathway drive,
        NOT an oscillating clock — previous versions used an unbounded ramp that never reset)
8. DNA: Host DNA synthesis activity
9. E: Effector immune cells (BKPyV-specific cellular response proxy)
10. IFN: Interferon concentration
11. AK: Antiviral state
12. D_tac: Tacrolimus dosing intensity (dimensionless, 0 = none, ~1 = full immunosuppression;
          NOT a plasma ng/mL concentration)
13. D_sir: Sirolimus dosing intensity (same convention)
14. P_rep: DNA replication pathway activity
15. P_immune: Innate immune pathway activity

Time unit: all rate constants are PER DAY. Clinical interpretation should compare against
weeks-scale plasma DNAemia (Funk 2006, PMID 16323135) rather than in-vitro hours.

Research Grounding:
- Drug mechanisms: Hirsch et al., Am J Transplant 2016;16(3):821-832 (in-vitro:
  sirolimus inhibits BKPyV replication via mTOR, IC90 ~ 4 ng/mL; tacrolimus
  increases replication in primary RPTECs). Modelled here as immune-control
  weakening (tacrolimus) and S-phase/production permissiveness (sirolimus).
- Viral clearance kinetics: Funk et al., J Infect Dis 2006;193:80-87
  (PMID 16323135; decay half-lives in patients).
- Host-cell pathway biology: Weissbach et al., J Virol 2024;98(12):e01382-24
  (single-cell transcriptomics of BKPyV in primary RPTECs); Needham et al.,
  PLoS Pathog 2024;20(12):e1012663 (host S phase precedes large T antigen).
"""

import numpy as np
from typing import Dict, Any, Optional


class BKPyVODESystem:
    """ODE system for BKPyV viral-host dynamics.
    
    This class defines the system of differential equations that govern
    BKPyV infection dynamics, drug effects, and host response.
    """
    
    def __init__(self, params: Optional[Dict[str, Any]] = None):
        """Initialize the ODE system with research-validated parameters.
        
        Args:
            params: Dictionary of ODE parameters (uses defaults if not provided)
        
        Default parameters are based on:
        - AJT-16-821.pdf: Drug effect mechanisms
        - abstract-16323135.txt: Viral clearance kinetics
        - Single-cell studies: Pathway activity dynamics
        """
        self.params = self._get_default_params()
        if params:
            self.params.update(params)
        self.applied_infections = set()  # Track applied infections
    
    def _get_default_params(self) -> Dict[str, float]:
        """Get default ODE parameters based on research findings.
        
        Returns:
            Dictionary of default parameter values
        """
        return {
            # Viral dynamics parameters
            'beta': 0.3,           # Infection rate (1/day), target-cell limited
            'delta': 0.4,          # Viral clearance rate (1/day) ≈ t1/2 1.7 d; within
                                   # the range reported after intervention change by
                                   # Funk 2006 (t1/2 6h-17d). NOT the 1-2h
                                   # post-nephrectomy phase.
            'p': 8.0,              # Viral production per infected cell (1/day);
                                   # tuned (with beta and immune_kill) so the default
                                   # regime supports persistent viremia (final V ~0.4,
                                   # peak ~3x inoculum) while clearance alone
                                   # (p small) self-resolves
            'c': 0.3,              # (legacy alias of infected-cell cytopathic death)

            # T antigen dynamics
            't_prod': 0.2,        # T antigen production rate (1/day)
            't_decay': 0.1,       # T antigen decay rate (1/day)
            't_threshold': 0.5,   # T antigen threshold separating early/late behaviour
            'nccr_early_expression_multiplier': 1.0,  # Archetype baseline; rearranged is a scenario
            'nccr_capsid_expression_multiplier': 1.0,  # Kept separate to expose the NCCR trade-off

            # Viral gene expression
            'g_prod': 0.15,        # Viral gene production rate (1/day)
            'g_decay': 0.08,      # Viral gene decay rate (1/day)

            # Cell dynamics
            'lambda_cell': 0.1,   # Constant target-cell regeneration source (cells/day,
                                  # Nowak-May form: dC/dt = lambda - d*C - beta*V*C)
            'd_cell': 0.05,        # Natural cell death rate (1/day)
            'd_infected': 0.3,     # Infected cell cytopathic death rate (1/day)
            'immune_kill': 0.8,    # Additional infected-cell killing per unit effector E

            # Cell cycle parameters (CC is a bounded permissiveness variable, 0-1)
            'cc_rate': 0.2,        # Relaxation rate of CC toward its set point (1/day)
            'cc_s_phase': 0.5,     # Hill half-point of the S-phase permissiveness gate

            # DNA synthesis dynamics
            'dna_prod': 0.3,      # DNA synthesis production rate (1/day), saturating
            'dna_decay': 0.15,     # DNA synthesis decay rate (1/day)
            'dna_coupling': 0.8,   # Host-viral replication coupling

            # Immune response parameters
            'e_prod': 0.05,       # Effector cell production rate (1/day)
            'e_decay': 0.02,      # Effector cell decay rate (1/day)
            'ifn_prod': 0.1,      # Interferon production rate (1/day)
            'ifn_decay': 0.08,     # Interferon decay rate (1/day)
            'ak_prod': 0.15,      # Antiviral state production rate (1/day)
            'ak_decay': 0.1,       # Antiviral state decay rate (1/day)
            'ak_max_enhancement': 2.0,  # Max extra viral clearance from antiviral state
            'ak_half': 1.0,            # Half-saturation of the antiviral clearance effect

            # Drug dosing (dimensionless intensity; see class docstring)
            'tac_absorption': 0.5,    # Tacrolimus onset rate while treatment active (1/day)
            'tac_clearance': 0.3,     # Tacrolimus washout rate after treatment stops (1/day)
            'sir_absorption': 0.4,    # Sirolimus onset rate while treatment active (1/day)
            'sir_clearance': 0.25,    # Sirolimus washout rate after treatment stops (1/day)

            # Research-grounded drug effects (Hirsch 2016, in vitro)
            'tac_enhancement': 1.5,    # Compatibility parameter; immune-control effect is used below
            'mtor_inhibition': 0.5,   # Residual permissiveness under full sirolimus (~IC90 4 ng/mL)
            'innate_immune_suppression': 0.5,  # Immune-control weakening per unit tacrolimus
            'sir_half': 0.5,          # Half-saturation of sirolimus effect
            'sir_late_weight': 0.3,   # Residual sirolimus effect after productive phase onset

            # Pathway dynamics
            'p_rep_prod': 0.2,     # DNA replication pathway production rate (1/day), saturating
            'p_rep_decay': 0.1,    # DNA replication pathway decay rate (1/day)
            'p_immune_prod': 0.15, # Immune pathway production rate (1/day)
            'p_immune_decay': 0.08, # Immune pathway decay rate (1/day)

            # Cell cycle and replication enhancement
            's_phase_bonus': 2.0,  # Extra DNA synthesis rate while the S-phase gate is open
            'ddr_enhancement': 1.3,  # DNA damage response enhancement (currently a scaling placeholder)

            # Single-cell pathway scaling (Weissbach 2024; phenomenological magnitudes)
            'translation_enhancement': 2.0,     # Translation pathway elevation
            'mitochondrial_importance': 0.8,   # Mitochondrial function importance
            'protein_degradation_inhibition': 0.3,  # Proteasome pathway involvement

            # Hill function parameters
            'hill_coeff': 4.0,     # Hill coefficient for smooth gates (CC/DNA)
            'half_saturation': 0.5,  # Half-saturation constant (T-antigen→production map)
        }
    
    def get_state_vector_names(self) -> list[str]:
        """Get names of state variables in order.
        
        Returns:
            List of state variable names
        """
        return [
            'V',      # Viral load
            'T',      # T antigen concentration
            'G_v',    # Viral gene expression
            'C',      # Healthy target cells
            'I',      # Infected cells
            'D',      # Dead/damaged cells
            'CC',     # Cell cycle phase
            'DNA',    # DNA synthesis activity
            'E',      # Effector immune cells
            'IFN',    # Interferon concentration
            'AK',     # Antiviral state
            'D_tac',  # Tacrolimus concentration
            'D_sir',  # Sirolimus concentration
            'P_rep',  # DNA replication pathway activity
            'P_immune' # Innate immune pathway activity
        ]
    
    def ode_system(self, t: float, y: np.ndarray,
                   dosing_context: Optional[Dict[str, Any]] = None) -> np.ndarray:
        """Define the ODE system.

        Args:
            t: Current time (days)
            y: State vector [V, T, G_v, C, I, D, CC, DNA, E, IFN, AK, D_tac, D_sir, P_rep, P_immune]
            dosing_context: Optional dict describing drug administration:
                {'tacrolimus': {'start': float, 'stop': float|None, 'target': float},
                 'sirolimus': {...}}
                ``target`` is the dimensionless dosing intensity (≈1.0 = full
                immunosuppression). While active, D approaches ``target`` with
                the drug's absorption rate; afterwards it decays with the
                drug's clearance rate. (Legacy ``{time: {drug: magnitude}}``
                bolus dicts are also accepted and applied exactly once when
                ``t`` crosses the key time.)

        Returns:
            Derivatives dy/dt (all per day)
        """
        # Unpack state variables
        V, T, G_v, C, I, D, CC, DNA, E, IFN, AK, D_tac, D_sir, P_rep, P_immune = y

        p = self.params

        # Safety floor for numerical excursions (states are formulated to stay
        # non-negative; the max() guards against adaptive-solver overshoot).
        V = max(0.0, V)
        T = max(0.0, T)
        G_v = max(0.0, G_v)
        C = max(0.0, C)
        I = max(0.0, I)
        D = max(0.0, D)
        CC = min(1.0, max(0.0, CC))
        DNA = max(0.0, DNA)
        E = max(0.0, E)
        IFN = max(0.0, IFN)
        AK = max(0.0, AK)
        D_tac = max(0.0, D_tac)
        D_sir = max(0.0, D_sir)
        P_rep = max(0.0, P_rep)
        P_immune = max(0.0, P_immune)

        # --- Drug administration -------------------------------------------
        # Continuous dosing: while within [start, stop) the intensity relaxes
        # toward `target`; afterwards it washes out. This replaces the old
        # single-bolus model, which let tacrolimus vanish within days.
        tac_active = sir_active = False
        tac_target = sir_target = 0.0
        tac_bolus = sir_bolus = 0.0
        if dosing_context:
            for drug, sched in dosing_context.items():
                if isinstance(sched, dict) and "target" in sched:
                    start = float(sched.get("start", 0.0))
                    stop = sched.get("stop", None)
                    active = start <= t and (stop is None or t < float(stop))
                    if drug == "tacrolimus":
                        tac_active, tac_target = active, float(sched["target"])
                    elif drug == "sirolimus":
                        sir_active, sir_target = active, float(sched["target"])
            # Legacy bolus form {time: {drug: magnitude}} (kept for backward
            # compatibility with old callers/tests)
            if all(isinstance(k, (int, float)) for k in dosing_context):
                for event_time, drugs in dosing_context.items():
                    if isinstance(drugs, dict) and abs(t - float(event_time)) < 1e-9:
                        tac_bolus += float(drugs.get("tacrolimus", 0.0))
                        sir_bolus += float(drugs.get("sirolimus", 0.0))

        D_tac_eff = D_tac + tac_bolus
        D_sir_eff = D_sir + sir_bolus

        # --- Drug mechanisms ------------------------------------------------
        # Tacrolimus: weakens immune control (saturating; never negative).
        tac_immune_effect = 1.0 / (1.0 + 2.0 * p['innate_immune_suppression'] * D_tac_eff)

        # Sirolimus: mTOR inhibition. Effect saturates with dose; strongest
        # before productive (late) phase, modelled by a smooth weight on T.
        sir_sat = D_sir_eff / (D_sir_eff + p['sir_half'])
        phase_weight = p['sir_late_weight'] + (1.0 - p['sir_late_weight']) / (
            1.0 + np.exp((T - p['t_threshold']) / 0.05)
        )
        sir_effect = 1.0 - sir_sat * (1.0 - p['mtor_inhibition']) * phase_weight
        sir_effect = max(0.0, sir_effect)
        mtor_rep_effect = 1.0 - 0.7 * sir_sat  # mTOR also drives host DNA-replication pathway

        # --- Viral dynamics -------------------------------------------------
        # dV/dt: production from infected cells, gated by T antigen supply and
        # sirolimus; clearance is increased (saturating) by the antiviral state.
        viral_production = p['p'] * I * (T / (T + p['half_saturation'])) * sir_effect
        ak_boost = p['ak_max_enhancement'] * AK / (AK + p['ak_half'])
        viral_clearance = p['delta'] * V * (1.0 + ak_boost)
        dVdt = viral_production - viral_clearance

        # dT/dt: T antigen requires host S phase and DNA synthesis first
        # (Needham 2024: S phase precedes T-antigen accumulation). Smooth
        # Hill gates keep the RHS Lipschitz for the adaptive solver.
        n = p['hill_coeff']
        gate_cc = CC**n / (p['cc_s_phase']**n + CC**n)
        gate_dna = DNA**n / (0.5**n + DNA**n)
        t_production = p['t_prod'] * G_v * P_rep * p['dna_coupling'] * gate_cc * gate_dna
        t_decay = p['t_decay'] * T
        dTdt = t_production - t_decay

        # dG_v/dt: Viral gene expression (saturating to keep G_v finite)
        g_production = (
            p['g_prod'] * V * p['translation_enhancement']
            * p['nccr_early_expression_multiplier'] / (1.0 + G_v)
        )
        g_decay = p['g_decay'] * G_v
        dG_vdt = g_production - g_decay

        # --- Cell populations ------------------------------------------------
        # dC/dt: constant regeneration source (Nowak-May target-cell supply),
        # infection loss, baseline turnover. A logistic C*(1-C) form collapses
        # permanently once C is depleted, which is why the previous model could
        # never sustain infection.
        infection_rate = p['beta'] * V * C
        cell_death = p['d_cell'] * C
        dCdt = p['lambda_cell'] - infection_rate - cell_death

        # dI/dt: infection gain minus cytopathic death and effector killing
        # (effector killing is weakened by tacrolimus via tac_immune_effect)
        infected_gain = infection_rate
        infected_death = p['d_infected'] * I + p['immune_kill'] * E * I * tac_immune_effect
        dIdt = infected_gain - infected_death

        # dD/dt: dead/damaged cells accumulate and clear slowly
        dDdt = cell_death + infected_death - 0.1 * D

        # --- Cell cycle permissiveness --------------------------------------
        # CC relaxes toward a set point driven by the DNA-replication pathway
        # and reduced by mTOR inhibition. Bounded in [0, 1) by construction.
        cc_setpoint = (P_rep * mtor_rep_effect) / (1.0 + P_rep * mtor_rep_effect)
        dCCdt = p['cc_rate'] * (cc_setpoint - CC)

        # --- Host DNA synthesis ---------------------------------------------
        # Saturating production driven by the replication pathway, with a bonus
        # while the S-phase gate is open.
        dna_production = p['dna_prod'] * P_rep * (1.0 + p['s_phase_bonus'] * gate_cc) / (1.0 + DNA)
        dna_decay = p['dna_decay'] * DNA
        dDNAdt = dna_production - dna_decay

        # --- Immune response -------------------------------------------------
        # dE/dt: effector cells expand with the innate immune pathway, with
        # self-limiting growth (without saturation the E <-> IFN <-> P_immune
        # loop is a runaway positive feedback that always clears the virus)
        dEdt = p['e_prod'] * P_immune / (1.0 + E) - p['e_decay'] * E

        # dIFN/dt: interferon induced by virus, sensed through effector tone
        dIFNdt = p['ifn_prod'] * E * V - p['ifn_decay'] * IFN

        # dAK/dt: antiviral state follows interferon
        dAKdt = p['ak_prod'] * IFN - p['ak_decay'] * AK

        # --- Drug pharmacokinetics ------------------------------------------
        if tac_active:
            dD_tacdt = p['tac_absorption'] * (tac_target - D_tac)
        else:
            # Decay applies to any present amount (incl. legacy bolus doses).
            dD_tacdt = -p['tac_clearance'] * D_tac_eff
        if sir_active:
            dD_sirdt = p['sir_absorption'] * (sir_target - D_sir)
        else:
            dD_sirdt = -p['sir_clearance'] * D_sir_eff

        # --- Pathway dynamics -----------------------------------------------
        # dP_rep/dt: saturating production from DNA synthesis, scaled by mTOR tone
        dP_repdt = p['p_rep_prod'] * (DNA / (1.0 + DNA)) * mtor_rep_effect - p['p_rep_decay'] * P_rep

        # dP_immune/dt: interferon-driven, suppressed by tacrolimus
        dP_immunedt = p['p_immune_prod'] * IFN * tac_immune_effect - p['p_immune_decay'] * P_immune

        return np.array([dVdt, dTdt, dG_vdt, dCdt, dIdt, dDdt, dCCdt, dDNAdt,
                        dEdt, dIFNdt, dAKdt, dD_tacdt, dD_sirdt, dP_repdt, dP_immunedt])
    
    def get_initial_conditions(self, cell_count: float = 1.0) -> np.ndarray:
        """Get initial conditions for the ODE system.
        
        Args:
            cell_count: Initial healthy cell count
        
        Returns:
            Initial state vector
        """
        return np.array([
            0.0,          # V: No initial virus
            0.0,          # T: No T antigen
            0.0,          # G_v: No viral gene expression
            cell_count,   # C: Healthy cells
            0.0,          # I: No infected cells
            0.0,          # D: No dead cells
            0.1,          # CC: Start in early G0/G1
            0.5,          # DNA: Baseline DNA synthesis
            0.1,          # E: Baseline immune cells
            0.0,          # IFN: No interferon
            0.0,          # AK: No antiviral state
            0.0,          # D_tac: No tacrolimus
            0.0,          # D_sir: No sirolimus
            0.7,          # P_rep: Baseline DNA replication pathway
            0.5           # P_immune: Baseline immune pathway
        ])
    
    def get_infection_conditions(self, viral_load: float = 0.5) -> np.ndarray:
        """Get initial conditions with active BKPyV infection.
        
        Args:
            viral_load: Initial viral load (0-1 scale)
        
        Returns:
            Initial state vector with infection
        """
        y0 = self.get_initial_conditions()
        y0[0] = viral_load  # Set viral load
        y0[2] = viral_load * 0.5  # Initial viral gene expression
        y0[4] = viral_load * 0.1  # Some initially infected cells
        y0[1] = 0.3  # Initial T antigen
        y0[3] -= y0[4]  # Reduce healthy cells
        return y0
