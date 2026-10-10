#!/usr/bin/env python3
"""Screening-trigger policy analysis: WHEN to act on rising DNAemia.

Kotton 2024 consensus screens plasma BKPyV DNAemia monthly (then q3mo)
and recommends acting at >=1,000 cp/mL persisting or >=10,000 cp/mL
(presumptive PyVAN). This script makes the model answer the question a
transplant program actually asks: does intervening at the screening
threshold beat waiting for presumptive PyVAN — and by how much?

Design: untreated infection rises from inoculum; each policy watches the
plasma trajectory at weekly resolution and fires the SAME intervention
(tacrolimus 8 -> 3 ng/mL + sirolimus 4 ng/mL conversion — the regimen
the optimizer found dominant) the week its trigger is crossed. Policies
differ only in trigger level, matching clinical practice. Reported:
peak DNAemia, weeks to sustained clearance, rearranged-NCCR fraction at
day 180, and T-cell rebound AUC.

Emergent result (mechanistic hypothesis, assumption-labelled): earlier
triggers cap both peak viremia and quasi-species emergence — the model's
mechanistic justification for intensive early screening.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402
from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402

HORIZON = 180.0
SCREENING_V = 0.2034  # ~1,000 cp/mL anchor
POLICIES = [
    ("act at 1k cp/mL (screening)", 1_000.0),
    ("act at 10k cp/mL (presumptive PyVAN)", 10_000.0),
    ("act at 100k cp/mL (late)", 100_000.0),
    ("never act", np.inf),
]


def run_policy(trigger_cpml, intervention="conversion", horizon=HORIZON):
    """Phase 1: untreated rise, sampled weekly. When plasma crosses the
    trigger, phase 2 applies the intervention schedule from that day."""
    mapper = ViralLoadMapper()
    ode = BKPyVODESystem()
    y0 = ode.get_infection_conditions(0.5)
    t_eval = np.linspace(0, horizon, int(horizon * 4) + 1)

    sol1 = solve_ivp(lambda t, y: ode.ode_system(t, y, None),
                     (0, horizon), y0, t_eval=t_eval, method="LSODA")
    cp1 = np.array([mapper.normalized_to_copies(float(v))
                    for v in np.maximum(sol1.y[0], 1e-9)])
    weekly = np.arange(0, horizon + 1e-9, 7.0)
    cp_weekly = np.interp(weekly, sol1.t, cp1)
    trigger_days = weekly[cp_weekly >= trigger_cpml]
    t_trigger = float(trigger_days[0]) if len(trigger_days) else np.inf

    if not np.isfinite(t_trigger):
        return _metrics(sol1, cp1, t_trigger, mapper)

    # continue from the pre-trigger state at the trigger day
    i_trig = int(np.searchsorted(sol1.t, t_trigger))
    y_trig = sol1.y[:, i_trig].copy()
    t2 = np.linspace(t_trigger, horizon,
                     max(2, int((horizon - t_trigger) * 4) + 1))
    if intervention == "conversion":
        dosing = {
            "tacrolimus": [{"start": 0.0, "stop": t_trigger, "trough_ng_ml": 8.0},
                           {"start": t_trigger, "stop": None, "trough_ng_ml": 3.0}],
            "sirolimus": [{"start": t_trigger, "stop": None, "trough_ng_ml": 4.0}],
        }
    else:  # taper only
        dosing = {"tacrolimus": [{"start": 0.0, "stop": t_trigger,
                                  "trough_ng_ml": 8.0},
                                 {"start": t_trigger, "stop": None,
                                  "trough_ng_ml": 3.0}]}
    sol2 = solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                     (t_trigger, horizon), y_trig, t_eval=t2,
                     method="LSODA")
    cp2 = np.array([mapper.normalized_to_copies(float(v))
                    for v in np.maximum(sol2.y[0], 1e-9)])
    return _metrics(sol2, np.concatenate([cp1[:i_trig], cp2]),
                    t_trigger, mapper)


def _metrics(sol, cp_full, t_trigger, mapper):
    V, V_u, T_eff = sol.y[0], sol.y[19], sol.y[16]
    F_k = sol.y[20] if sol.y.shape[0] > 20 else np.zeros_like(V)
    F_u = sol.y[21] if sol.y.shape[0] > 21 else np.zeros_like(V)
    below = cp_full < 1_000.0
    clear_weeks = np.inf
    tt = np.linspace(0, HORIZON, len(cp_full))
    for i in range(len(cp_full)):
        if below[i:].all():
            clear_weeks = tt[i] / 7.0
            break
    return {
        "trigger_day": None if np.isinf(t_trigger) else round(t_trigger, 1),
        "peak_log10_cpml": round(float(np.log10(max(cp_full.max(), 1.0))), 2),
        "clearance_weeks": None if np.isinf(clear_weeks) else round(float(clear_weeks), 1),
        "final_frr_kidney": round(float(F_k[-1]), 3),
        "final_frr_urine": round(float(F_u[-1]), 3),
        "rebound_index": round(float(np.trapezoid(T_eff, sol.t) / HORIZON), 3),
        "final_u_p_ratio": round(float(V_u[-1] / max(V[-1], 1e-12)), 1),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--taper", action="store_true",
                    help="intervene with tac taper alone instead of conversion")
    args = ap.parse_args()
    intervention = "taper" if args.taper else "conversion"

    print(f"Screening-trigger policy analysis ({intervention} at trigger)")
    print("(mechanistic hypothesis generator — assumption-labelled)\n")
    hdr = (f"{'policy':<40} {'trig d':>7} {'peak log10':>11} {'clear@wk':>9} "
           f"{'F_rr k':>7} {'F_rr u':>7} {'rebound':>8}")
    print(hdr)
    print("-" * len(hdr))
    results = []
    for label, trig in POLICIES:
        m = run_policy(trig, intervention=intervention)
        m["policy"] = label
        results.append(m)
        trig_d = "-" if m["trigger_day"] is None else f"{m['trigger_day']:.0f}"
        clr = "never" if m["clearance_weeks"] is None else f"{m['clearance_weeks']:.1f}"
        print(f"{label:<40} {trig_d:>7} {m['peak_log10_cpml']:>11} {clr:>9} "
              f"{m['final_frr_kidney']:>7} {m['final_frr_urine']:>7} {m['rebound_index']:>8}")

    if args.json:
        out = Path("outputs")
        out.mkdir(exist_ok=True)
        (out / "screening_policies.json").write_text(json.dumps(results, indent=2))
        print(f"\nWrote {out/'screening_policies.json'}")


if __name__ == "__main__":
    main()
