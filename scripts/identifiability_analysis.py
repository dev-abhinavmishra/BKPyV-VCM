#!/usr/bin/env python3
"""Practical identifiability analysis for the BKPyV ODE model.

Question that matters for the project (and for the patient-data
collaboration): given ONLY serial plasma qPCR — the clinically available
measurement — which model parameters are actually learnable, and which
require the urine series or a different observable?

Method: normalized sensitivity analysis (the standard first-pass
identifiability screen for ODE models, e.g. Raue 2009; Wu 2008). For
each parameter theta we compute

    S_theta(t) = d log V_obs(t) / d log theta

by central finite differences (+-5%), then take the L2 norm over the
simulation horizon as the parameter's information content for that
observable. Observables examined separately: plasma V, urine V_u, and
the intracellular T-antigen proxy (unmeasurable clinically — included to
show what a tissue readout would add).

A parameter with near-zero sensitivity under EVERY observable is
practically unidentifiable — its value cannot be learned from data and
must come from literature priors. Parameters with collinear sensitivity
profiles are flagged as a confounded pair (only their combination is
learnable). Results are assumption-labelled: ranks depend on the
baseline regimen simulated.
"""

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402

HORIZON = 90.0
STEP = 0.05  # +-5% log-perturbation

# Parameters probed — grouped by which biological block they live in.
PROBE_PARAMS = [
    # viral kinetics (Nowak-May block)
    "beta", "p", "delta", "c", "p_u", "beta_u", "delta_u",
    "lambda_u", "d_infected_u", "seed_kidney_u", "drain_kidney", "cross_feed",
    # intracellular gates
    "t_prod", "t_decay", "g_prod", "g_decay", "cc_rate", "dna_coupling",
    "dna_decay", "dna_prod",
    # innate immunity
    "immune_kill", "e_prod", "ak_prod", "ak_max_enhancement",
    # adaptive (BKPyV T-cell arm)
    "a_tcell", "prolif_tcell", "d_teff", "tcell_kill", "tac_tcell_suppression",
    "tcell_antigen_i_weight", "lambda_tcell", "d_tcell", "tcell_carry",
    # NCCR quasi-species
    "nccr_emergence_rate", "nccr_selection_rate", "nccr_reversion_rate",
    "rr_early_gain", "rr_capsid_fraction",
    # drug PD
    "tac_immune_effect", "mtor_inhibition", "tac_enhancement",
]
OBSERVABLES = {"plasma_V": 0, "urine_V": 19, "t_antigen": 1, "F_rr": 20}
COLLINEARITY_THRESHOLD = 0.98  # |corr| of sensitivity profiles


def simulate_traj(params_overrides=None, y0=None):
    ode = BKPyVODESystem()
    if params_overrides:
        ode.params.update(params_overrides)
    if y0 is None:
        y0 = ode.get_infection_conditions(0.5)
    t_eval = np.linspace(0, HORIZON, int(HORIZON * 4) + 1)
    sol = solve_ivp(lambda t, y: ode.ode_system(t, y, None),
                    (0, HORIZON), y0, t_eval=t_eval, method="LSODA")
    return sol


def sensitivity(param, obs_idx, base_params, y0):
    """d log obs / d log param via +-5% central differences."""
    theta = base_params[param]
    if theta <= 0:
        return None
    d = theta * STEP
    sol_hi = simulate_traj({param: theta + d}, y0)
    sol_lo = simulate_traj({param: theta - d}, y0)
    obs_hi = np.maximum(sol_hi.y[obs_idx], 1e-12)
    obs_lo = np.maximum(sol_lo.y[obs_idx], 1e-12)
    return (np.log(obs_hi) - np.log(obs_lo)) / (2.0 * np.log(1.0 + STEP))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true",
                    help="write sensitivity table to outputs/identifiability.json")
    ap.add_argument("--fast", action="store_true",
                    help="probe a reduced parameter set for smoke tests")
    args = ap.parse_args()

    ode = BKPyVODESystem()
    base = ode.params
    y0 = ode.get_infection_conditions(0.5)
    params = PROBE_PARAMS[:12] if args.fast else PROBE_PARAMS

    # Sensitivity tensor: param -> observable -> |S|_L2
    table = {}
    profiles = {}
    for name in params:
        if name not in base:
            continue
        row = {}
        for obs_name, obs_idx in OBSERVABLES.items():
            s = sensitivity(name, obs_idx, base, y0)
            if s is None:
                continue
            row[obs_name] = float(np.linalg.norm(s) / np.sqrt(len(s)))
            profiles[(name, obs_name)] = s
        table[name] = row

    # Collinearity on the plasma-V observable: pairs whose normalized
    # sensitivity shapes are |corr| > threshold are confounded pairwise.
    confounded = []
    plasma_names = [n for n in table if "plasma_V" in table[n]]
    for a, b in combinations(plasma_names, 2):
        sa, sb = profiles[(a, "plasma_V")], profiles[(b, "plasma_V")]
        na, nb = np.linalg.norm(sa), np.linalg.norm(sb)
        if na < 1e-9 or nb < 1e-9:
            continue
        corr = abs(float(np.dot(sa, sb) / (na * nb)))
        if corr > COLLINEARITY_THRESHOLD:
            confounded.append((a, b, round(corr, 3)))

    ranked = sorted(table.items(),
                    key=lambda kv: -kv[1].get("plasma_V", 0.0))
    print("Identifiability screen — RMS sensitivity of each observable to a")
    print("+-5% log parameter perturbation (baseline: untreated infection, 90 d)\n")
    hdr = f"{'parameter':<24} {'plasma_V':>9} {'urine_V':>9} {'T_ag':>9} {'F_rr':>9}"
    print(hdr)
    print("-" * len(hdr))
    for name, row in ranked:
        print(f"{name:<24} {row.get('plasma_V', 0):>9.3f} {row.get('urine_V', 0):>9.3f} "
              f"{row.get('t_antigen', 0):>9.3f} {row.get('F_rr', 0):>9.3f}")

    print("\nConfounded parameter pairs on plasma V alone (|corr| > "
          f"{COLLINEARITY_THRESHOLD}):")
    for a, b, c in confounded[:12]:
        print(f"  {a:<22} <-> {b:<22} corr={c}")

    if args.json:
        out = Path("outputs")
        out.mkdir(exist_ok=True)
        payload = {
            "method": "central-difference normalized sensitivity, +-5%",
            "horizon_days": HORIZON,
            "rms_sensitivity": table,
            "confounded_pairs_plasma": confounded,
            "note": ("Ranks are baseline-regimen dependent; near-zero rows are "
                     "practically unidentifiable and must come from literature "
                     "priors, not data."),
        }
        (out / "identifiability.json").write_text(json.dumps(payload, indent=2))
        print(f"\nWrote {out/'identifiability.json'}")


if __name__ == "__main__":
    main()
