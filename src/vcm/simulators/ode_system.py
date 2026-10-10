"""ODE system for BKPyV viral-host dynamics.

This module defines the system of ordinary differential equations (ODEs) that model
BK polyomavirus infection in kidney tubular epithelial cells, including drug effects
and host response dynamics.

The ODE system is based on standard virus-host modeling approaches extended with:
- Drug-specific effects (tacrolimus vs sirolimus)
- Pathway activity dynamics
- Cell cycle progression
- Immune response kinetics
- BKPyV-specific adaptive T-cell arm (activation → clonal expansion → CTL killing)
- A urinary/urothelial compartment with kidney↔bladder cross-feeding
- Research-validated parameters from clinical studies

State Variables (20-dimensional vector):
1. V: Viral load (normalised scale; see clinical.viral_load_mapper for the bridge to copies/mL)
2. T: Large T-antigen level (arbitrary units)
3. G_v: Viral gene expression level (arbitrary units)
4. C: Healthy target cells (fraction of carrying capacity)
5. I: Infected cells (fraction)
6. D: Dead/damaged cells (fraction)
7. CC: Cell-cycle permissiveness (bounded 0-1; relaxation toward replication-pathway drive,
        NOT an oscillating clock — previous versions used an unbounded ramp that never reset)
8. DNA: Host DNA synthesis activity
9. E: Effector immune cells (innate/nonspecific cellular response proxy)
10. IFN: Interferon concentration
11. AK: Antiviral state
12. D_tac: Tacrolimus dosing intensity (dimensionless, 0 = none, ~1 = full immunosuppression;
          NOT a plasma ng/mL concentration)
13. D_sir: Sirolimus dosing intensity (same convention)
14. P_rep: DNA replication pathway activity
15. P_immune: Innate immune pathway activity
16. T_naive: BKPyV-specific naive/precursor T cells (repertoire fraction)
17. T_eff: BKPyV-specific effector T cells (the antigen-driven CTL arm; distinct
          from the nonspecific E proxy — clinically the arm measured by T-cell
          monitoring assays and boosted by virus-specific T-cell therapy)
18. C_u: Healthy urothelial cells (urinary/bladder compartment, fraction)
19. I_u: Infected urothelial cells (fraction)
20. V_u: Urinary viral load (virions in the bladder/urine compartment)

Time unit: all rate constants are PER DAY. Clinical interpretation should compare against
weeks-scale plasma DNAemia (Funk 2006, PMID 16323135) rather than in-vitro hours.

Research Grounding:
- Drug mechanisms: Hirsch et al., Am J Transplant 2016;16(3):821-832 (in-vitro:
  sirolimus inhibits BKPyV replication via mTOR, IC90 ~ 4 ng/mL; tacrolimus
  increases replication in primary RPTECs). Modelled here as immune-control
  weakening (tacrolimus) and S-phase/production permissiveness (sirolimus).
  Tacrolimus additionally suppresses the virus-specific T-cell arm more
  strongly than the innate arm — calcineurin/NFAT blockade's primary clinical
  target is T-cell activation (Kotton et al., Transplantation 2024 consensus).
- Viral clearance kinetics: Funk et al., J Infect Dis 2006;193:80-87
  (PMID 16323135; decay half-lives in patients).
- Host-cell pathway biology: Weissbach et al., J Virol 2024;98(12):e01382-24
  (single-cell transcriptomics of BKPyV in primary RPTECs); Needham et al.,
  PLoS Pathog 2024;20(12):e1012663 (host S phase precedes large T antigen).
- Urinary tract dynamics: Funk et al., Am J Transplant 2008;8 (within-host
  two-compartment model: kidney replication seeds the urothelium, which
  amplifies and dominates the urinary load — >95% of urine load is urothelial,
  urine ~3000× plasma — with cross-feeding back to the kidney). The urothelial
  compartment is deliberately population-level (no per-cell T-antigen gates):
  production is an aggregate amplification term, matching Funk's formulation.
- Virus-specific T-cell control is the clinical determinant of clearance
  (Kotton 2024 consensus; the basis for virus-specific T-cell therapy), so the
  T-cell arm is modelled explicitly rather than folded into the generic E pool.
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
            'nccr_early_expression_multiplier': 1.0,  # Archetype baseline (F_rr=0 endpoint)
            'nccr_capsid_expression_multiplier': 1.0,  # Archetype baseline (F_rr=0 endpoint)

            # NCCR quasi-species dynamics (Gosert 2008; Broekema 2021):
            # archetype is the transmitted/persistent form; rearranged NCCR
            # variants emerge in vivo under sustained replication and
            # outcompete (early-gene overexpression -> faster genome
            # copying). F_rr is the rearranged fraction of the virion pool;
            # the discrete archetype/rearranged presets are the F=0/F=1
            # boundary conditions of this dynamics.
            'rr_early_gain': 2.0,          # Early-gene expression at F_rr=1 (2x archetype)
            'rr_capsid_fraction': 0.5,     # Capsid expression at F_rr=1 (0.5x archetype)
            'nccr_emergence_rate': 0.01,   # Rearrangement supply per unit replication (1/day)
            'nccr_selection_rate': 0.8,    # Competitive advantage of rr under replication (1/day)
            'nccr_reversion_rate': 0.0005, # Slow drift back toward archetype (1/day)
            'nccr_emergence_enabled': 1.0, # 0 restores the fixed-genotype model

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

            # Pharmacokinetics in real units (used when a dosing schedule
            # specifies ``trough_ng_ml`` instead of dimensionless ``target``).
            # One-compartment approximation at maintenance dosing: the state
            # relaxes toward the target trough with the drug's elimination
            # half-life; peak-trough oscillation is averaged out.
            'tac_pk_half_life_days': 0.5,   # Tacrolimus t½ ~12 h
            'sir_pk_half_life_days': 2.5,   # Sirolimus t½ ~60 h (slow washout —
                                            # why tac→sir switches take weeks)
            'tac_ref_trough_ngml': 8.0,     # Trough producing unit dosing
                                            # intensity (typical early-maintenance
                                            # trough 5-10 ng/mL)
            'sir_ref_trough_ngml': 4.0,     # Trough producing unit dosing
                                            # intensity; anchored to the in-vitro
                                            # antiviral IC90 ~4 ng/mL (Hirsch 2016),
                                            # so real-world 5-10 ng/mL troughs
                                            # saturate the mTOR effect

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

            # --- BKPyV-specific T-cell arm (adaptive immunity) ------------
            # Naive cells are replenished slowly and primed by antigen drive
            # (free virions + infected-cell display); effectors expand
            # antigen-dependently and kill infected cells. Tacrolimus blocks
            # activation/expansion (calcineurin→NFAT), not the cytolytic hit.
            'lambda_tcell': 0.005,        # Naive T-cell replenishment (1/day)
            'a_tcell': 0.3,               # Antigen-driven priming rate (1/day)
            'd_tcell': 0.005,             # Naive T-cell turnover (1/day; long-lived)
            'prolif_tcell': 0.6,          # Effector clonal expansion rate (1/day)
            'tcell_carry': 1.0,           # Effector carrying capacity (fraction scale)
            'd_teff': 0.05,               # Effector T-cell decay (1/day; t½ ~14 d)
            'tcell_kill': 0.25,           # Infected-cell killing per unit T_eff (1/day)
            'tcell_antigen_i_weight': 0.5, # Infected-cell contribution to antigen drive
            'tac_tcell_suppression': 12.0, # Tacrolimus blockade of T-cell priming/
                                          # expansion — much stronger than the
                                          # innate-arm effect (calcineurin's
                                          # primary target; clinical troughs are
                                          # several× the T-cell-proliferation IC50)

            # --- Urinary/urothelial compartment (Funk 2008) ----------------
            # Population-level bladder compartment: kidney virions drain into
            # urine and seed urothelial infection; the urothelium amplifies
            # (>95% of the urine load is urothelial-derived); a small fraction
            # feeds back onto the kidney (cross-feeding). No intracellular
            # gates here — production is aggregate, per Funk's formulation.
            'lambda_u': 0.05,             # Urothelial regeneration source (cells/day)
            'd_cell_u': 0.03,             # Urothelial natural turnover (1/day)
            'd_infected_u': 0.2,          # Infected urothelial death (1/day)
            'beta_u': 0.8,                # Urothelial infection rate (1/day)
            'seed_kidney_u': 0.05,        # Kidney→urothelium seeding pressure on V
            'p_u': 500.0,                 # Urothelial virion production per infected
                                          # cell (1/day) — amplification compartment
            'drain_kidney': 0.5,          # Kidney virions appearing in urine (1/day)
            'delta_u': 0.5,               # Urinary clearance/washout (1/day)
            'immune_kill_u': 0.1,         # Sparse bladder immune surveillance per E
            'cross_feed': 0.0005,         # Bladder→kidney reinfection pressure on V_u;
                                          # deliberately weak — viruria without
                                          # viremia is clinically common, so the
                                          # bladder reservoir can sustain only
                                          # low-level reseeding
        }
    
    def get_state_vector_names(self) -> list[str]:
        """Get names of state variables in order.
        
        Returns:
            List of state variable names
        """
        return [
            'V',       # Viral load
            'T',       # T antigen concentration
            'G_v',     # Viral gene expression
            'C',       # Healthy target cells
            'I',       # Infected cells
            'D',       # Dead/damaged cells
            'CC',      # Cell cycle phase
            'DNA',     # DNA synthesis activity
            'E',       # Effector immune cells
            'IFN',     # Interferon concentration
            'AK',      # Antiviral state
            'D_tac',   # Tacrolimus concentration
            'D_sir',   # Sirolimus concentration
            'P_rep',   # DNA replication pathway activity
            'P_immune',# Innate immune pathway activity
            'T_naive', # BKPyV-specific naive T cells
            'T_eff',   # BKPyV-specific effector T cells
            'C_u',     # Healthy urothelial cells
            'I_u',     # Infected urothelial cells
            'V_u',     # Urinary viral load
            'F_rr',    # Rearranged-NCCR fraction of the virion pool
        ]
    
    def ode_system(self, t: float, y: np.ndarray,
                   dosing_context: Optional[Dict[str, Any]] = None) -> np.ndarray:
        """Define the ODE system.

        Args:
            t: Current time (days)
            y: State vector [V, T, G_v, C, I, D, CC, DNA, E, IFN, AK, D_tac,
               D_sir, P_rep, P_immune, T_naive, T_eff, C_u, I_u, V_u, F_rr]
            dosing_context: Optional dict describing drug administration:
                {'tacrolimus': {'start': float, 'stop': float|None,
                                'target': float | 'trough_ng_ml': float} | [ ... ],
                 'sirolimus': {...}}
                Each drug maps to one schedule dict or a LIST of schedule
                dicts (stepwise regimens: e.g. a tacrolimus taper is
                [{start 0, stop 42, trough_ng_ml 8}, {start 42, trough 4}]).
                When windows overlap, the strongest intensity wins.
                ``target`` is the dimensionless dosing intensity (≈1.0 = full
                immunosuppression) driven with the ad-hoc absorption/clearance
                envelope. ``trough_ng_ml`` is a real steady-state trough
                concentration in ng/mL — converted to intensity through the
                drug's reference trough and driven with first-order kinetics
                at the drug's elimination half-life (tac ~12 h, sir ~60 h).
                If both are present, ``trough_ng_ml`` wins. (Legacy
                ``{time: {drug: magnitude}}`` bolus dicts are also accepted
                and applied exactly once when ``t`` crosses the key time.)

        Returns:
            Derivatives dy/dt (all per day)
        """
        # Unpack state variables (kidney/intracellular block first, then the
        # appended adaptive-T-cell and urothelial compartments)
        V, T, G_v, C, I, D, CC, DNA, E, IFN, AK, D_tac, D_sir, P_rep, P_immune = y[:15]
        T_naive, T_eff, C_u, I_u, V_u = y[15:20]
        F_rr = float(y[20]) if len(y) > 20 else 0.0

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
        T_naive = max(0.0, T_naive)
        T_eff = max(0.0, T_eff)
        C_u = max(0.0, C_u)
        I_u = max(0.0, I_u)
        V_u = max(0.0, V_u)
        F_rr = min(1.0, max(0.0, F_rr))

        # --- Drug administration -------------------------------------------
        # Continuous dosing: while within [start, stop) the intensity relaxes
        # toward `target`; afterwards it washes out. This replaces the old
        # single-bolus model, which let tacrolimus vanish within days.
        tac_active = sir_active = False
        tac_target = sir_target = 0.0
        tac_bolus = sir_bolus = 0.0
        tac_pk_mode = sir_pk_mode = False
        if dosing_context:
            for drug, scheds in dosing_context.items():
                # Each drug maps to a schedule LIST (a single dict is treated
                # as a one-window list) so stepwise regimens — e.g. a
                # tacrolimus taper 8 -> 4 ng/mL — are expressible. When
                # several windows are active the strongest intensity wins.
                if isinstance(scheds, dict):
                    scheds = [scheds]
                if not isinstance(scheds, list):
                    continue
                for sched in scheds:
                    if not isinstance(sched, dict):
                        continue
                    start = float(sched.get("start", 0.0))
                    stop = sched.get("stop", None)
                    active = start <= t and (stop is None or t < float(stop))
                    # ng/mL schedules: ``trough_ng_ml`` is the steady-state
                    # target trough; converted to dosing intensity via the
                    # drug's reference trough and driven with first-order PK
                    # kinetics. A dimensionless ``target`` keeps the legacy
                    # intensity path.
                    trough = sched.get("trough_ng_ml", None)
                    if drug == "tacrolimus":
                        if trough is not None:
                            tac_pk_mode = True
                            if active:
                                tac_active = True
                                tac_target = max(tac_target,
                                                 float(trough) / p['tac_ref_trough_ngml'])
                        elif "target" in sched and active:
                            tac_active = True
                            tac_target = max(tac_target, float(sched["target"]))
                    elif drug == "sirolimus":
                        if trough is not None:
                            sir_pk_mode = True
                            if active:
                                sir_active = True
                                sir_target = max(sir_target,
                                                 float(trough) / p['sir_ref_trough_ngml'])
                        elif "target" in sched and active:
                            sir_active = True
                            sir_target = max(sir_target, float(sched["target"]))
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

        # dG_v/dt: Viral gene expression (saturating to keep G_v finite).
        # The early-gene multiplier is F_rr-weighted: as the rearranged
        # fraction grows, early-gene expression approaches rr_early_gain x
        # the archetype baseline (Gosert 2008: rr-NCCR overexpresses early
        # genes).
        nccr_early_eff, nccr_capsid_eff = self._nccr_multipliers(F_rr)
        g_production = (
            p['g_prod'] * V * p['translation_enhancement']
            * nccr_early_eff / (1.0 + G_v)
        )
        g_decay = p['g_decay'] * G_v
        dG_vdt = g_production - g_decay

        # --- Cell populations ------------------------------------------------
        # dC/dt: constant regeneration source (Nowak-May target-cell supply),
        # infection loss, baseline turnover. A logistic C*(1-C) form collapses
        # permanently once C is depleted, which is why the previous model could
        # never sustain infection.
        # Kidney infection pressure includes a small bladder→kidney
        # cross-feeding term (Funk 2008: urothelial amplification reseeds the
        # graft) — kept deliberately small so it is a persistence mechanism,
        # not a primary driver.
        infection_rate = p['beta'] * (V + p['cross_feed'] * V_u) * C
        cell_death = p['d_cell'] * C
        dCdt = p['lambda_cell'] - infection_rate - cell_death

        # dI/dt: infection gain minus cytopathic death, nonspecific effector
        # killing (weakened by tacrolimus via tac_immune_effect), and
        # BKPyV-specific CTL killing by T_eff. The T-cell hit itself is NOT
        # attenuated by tacrolimus — calcineurin blockade acts on priming and
        # expansion (below), not on the cytolytic synapse.
        infected_gain = infection_rate
        infected_death = (p['d_infected'] * I
                          + p['immune_kill'] * E * I * tac_immune_effect
                          + p['tcell_kill'] * T_eff * I)
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
        # ng/mL schedules use the drug's real elimination half-life for both
        # approach-to-steady-state and washout (one-compartment approximation;
        # dosing-frequency oscillation averaged out). Dimensionless schedules
        # keep the legacy absorption/clearance envelope.
        tac_pk_rate = np.log(2.0) / p['tac_pk_half_life_days']
        sir_pk_rate = np.log(2.0) / p['sir_pk_half_life_days']
        if tac_active:
            rate = tac_pk_rate if tac_pk_mode else p['tac_absorption']
            dD_tacdt = rate * (tac_target - D_tac)
        else:
            # Decay applies to any present amount (incl. legacy bolus doses).
            rate = tac_pk_rate if tac_pk_mode else p['tac_clearance']
            dD_tacdt = -rate * D_tac_eff
        if sir_active:
            rate = sir_pk_rate if sir_pk_mode else p['sir_absorption']
            dD_sirdt = rate * (sir_target - D_sir)
        else:
            rate = sir_pk_rate if sir_pk_mode else p['sir_clearance']
            dD_sirdt = -rate * D_sir_eff

        # --- Pathway dynamics -----------------------------------------------
        # dP_rep/dt: saturating production from DNA synthesis, scaled by mTOR tone
        dP_repdt = p['p_rep_prod'] * (DNA / (1.0 + DNA)) * mtor_rep_effect - p['p_rep_decay'] * P_rep

        # dP_immune/dt: interferon-driven, suppressed by tacrolimus
        dP_immunedt = p['p_immune_prod'] * IFN * tac_immune_effect - p['p_immune_decay'] * P_immune

        # --- BKPyV-specific T-cell arm --------------------------------------
        # Antigen drive: free virions plus infected-cell display. Tacrolimus
        # blunts priming and clonal expansion (calcineurin/NFAT) with a
        # saturating coefficient stronger than its innate-arm effect — this is
        # the clinically dominant reason tacrolimus regimens carry the highest
        # BKPyV risk (Demey 2018 meta-analysis; Kotton 2024 consensus).
        antigen_drive = V + p['tcell_antigen_i_weight'] * I
        tac_tcell_effect = 1.0 / (1.0 + p['tac_tcell_suppression'] * D_tac_eff)

        # dT_naive/dt: slow replenishment, antigen-driven recruitment into the
        # effector pool (consumes naive cells), baseline turnover.
        tcell_priming = p['a_tcell'] * antigen_drive * T_naive * tac_tcell_effect
        dT_naivedt = p['lambda_tcell'] - tcell_priming - p['d_tcell'] * T_naive

        # dT_eff/dt: priming input plus antigen-dependent expansion with
        # logistic ceiling, minus effector decay.
        tcell_expansion = (p['prolif_tcell'] * antigen_drive * T_eff
                           * (1.0 - T_eff / p['tcell_carry']) * tac_tcell_effect)
        dT_effdt = tcell_priming + tcell_expansion - p['d_teff'] * T_eff

        # --- Urinary/urothelial compartment (Funk 2008) ---------------------
        # Urothelial infection pressure: local spread (beta_u * V_u) plus
        # seeding by kidney-derived virions draining into the bladder.
        uro_infection_rate = (p['beta_u'] * V_u + p['seed_kidney_u'] * V) * C_u
        dC_udt = p['lambda_u'] - uro_infection_rate - p['d_cell_u'] * C_u
        dI_udt = (uro_infection_rate
                  - p['d_infected_u'] * I_u
                  - p['immune_kill_u'] * E * I_u * tac_immune_effect)
        # Urinary viral load = aggregate urothelial production + kidney
        # drainage − urinary washout. Amplification here is what makes urine
        # loads exceed plasma by orders of magnitude (>95% urothelial origin).
        dV_udt = p['p_u'] * I_u + p['drain_kidney'] * V - p['delta_u'] * V_u

        # --- NCCR quasi-species dynamics ------------------------------------
        # Rearranged variants arise in proportion to replication activity
        # (mutation supply) and are competitively favoured by it (their
        # early-gene overexpression copies genomes faster); slow reversion
        # toward archetype. F_rr therefore stays ~0 in quiescent infection
        # and rises toward 1 under sustained high viremia — the observed
        # clinical pattern (rr-NCCR marks high-load plasma).
        replication_pressure = I * (T / (T + p['half_saturation']))
        if p['nccr_emergence_enabled'] >= 0.5:
            dF_rrdt = (p['nccr_emergence_rate'] * replication_pressure * (1.0 - F_rr)
                       + p['nccr_selection_rate'] * replication_pressure * F_rr * (1.0 - F_rr)
                       - p['nccr_reversion_rate'] * F_rr)
        else:
            dF_rrdt = 0.0

        return np.array([dVdt, dTdt, dG_vdt, dCdt, dIdt, dDdt, dCCdt, dDNAdt,
                        dEdt, dIFNdt, dAKdt, dD_tacdt, dD_sirdt, dP_repdt, dP_immunedt,
                        dT_naivedt, dT_effdt, dC_udt, dI_udt, dV_udt, dF_rrdt])
    
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
            0.5,          # P_immune: Baseline immune pathway
            0.02,         # T_naive: small primed repertoire (seropositive recipient)
            0.02,         # T_eff: low pre-existing virus-specific memory
            1.0,          # C_u: full urothelial compartment
            0.0,          # I_u: no urothelial infection
            0.0,          # V_u: no urinary viral load
            0.0,          # F_rr: archetype NCCR at transmission
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

    def _nccr_multipliers(self, F_rr: float) -> tuple:
        """Effective early-gene and capsid expression multipliers at a given
        rearranged fraction. The configured ``nccr_*_multiplier`` params are
        the F_rr=0 (archetype) endpoints; F_rr=1 approaches ``rr_early_gain``x
        early expression and ``rr_capsid_fraction``x capsid expression
        (Gosert 2008: rr-NCCR overexpresses early genes ~2x, capsid ~0.5x).
        """
        p = self.params
        early = p['nccr_early_expression_multiplier'] * (
            1.0 + F_rr * (p['rr_early_gain'] - 1.0))
        capsid = p['nccr_capsid_expression_multiplier'] * (
            1.0 - F_rr * (1.0 - p['rr_capsid_fraction']))
        return early, capsid
