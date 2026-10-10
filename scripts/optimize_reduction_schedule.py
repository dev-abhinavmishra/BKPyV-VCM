#!/usr/bin/env python3
"""Model-guided immunosuppression-reduction schedules for BKPyV.

The clinical dilemma (Funk 2008; Kotton 2024): reducing immunosuppression
is the only effective lever against BKPyV viremia, but every increment of
reduction buys rejection risk. This script turns the 21-dimensional ODE
into a decision-support tool:

  * Schedules are specified in REAL units — tacrolimus/sirolimus troughs
    in ng/mL with stepwise tapers — using the PK layer in
    ``BKPyVODESystem.ode_system``.
  * Viral efficacy is reported in clinical terms: weeks until sustained
    plasma viremia < 1,000 copies/mL (Kotton screening threshold), and
    the re-equilibrated load at the horizon.
  * Rejection counterweight: the integrated BKPyV-relevant T-cell
    rebound after reduction. The model cannot distinguish BKPyV-specific
    from alloreactive T cells, so T_eff expansion is the honest proxy —
    the same adaptive recovery that clears virus is what threatens the
    graft. Reported as "rebound index" = AUC(T_eff) over the horizon,
    normalized to the untreated reference.

Outputs a ranked table of candidate schedules (and optional JSON under
outputs/). Assumption-labelled throughout: these are mechanistic
hypotheses for protocol design, NOT fitted patient recommendations.
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

HORIZON_DAYS = 180.0
SCREENING_V = 0.2034  # V corresponding to 1,000 cp/mL (ViralLoadMapper)


def evaluate_schedule(schedule, sir_trough=None, horizon=HORIZON_DAYS,
                      n_points=721):
    """Simulate one taper schedule; return clinical metrics dict.

    ``schedule``: list of (start_day, stop_day|None, tac_trough_ng_ml).
    ``sir_trough``: optional sirolimus trough in ng/mL started at the
    first taper step (a tac->sir conversion scenario).
    """
    ode = BKPyVODESystem()
    y0 = ode.get_infection_conditions(0.5)
    t_eval = np.linspace(0, horizon, n_points)
    dosing = {
        "tacrolimus": [
            {"start": s, "stop": e, "trough_ng_ml": trough}
            for s, e, trough in schedule
        ],
    }
    if sir_trough is not None:
        dosing["sirolimus"] = [
            {"start": schedule[-1][0], "stop": None,
             "trough_ng_ml": sir_trough},
        ]
    sol = solve_ivp(
        lambda t, y: ode.ode_system(t, y, dosing),
        (0, horizon), y0, t_eval=t_eval, method="LSODA",
    )
    V, V_u, T_eff, F_rr = sol.y[0], sol.y[19], sol.y[16], sol.y[20]

    below = V < SCREENING_V
    # sustained clearance: never re-rises above threshold afterwards
    clear_weeks = np.inf
    for i, t in enumerate(sol.t):
        if below[i:].all():
            clear_weeks = t / 7.0
            break

    rebound_index = float(np.trapezoid(T_eff, sol.t) / horizon)
    return {
        "clearance_weeks": None if np.isinf(clear_weeks) else round(float(clear_weeks), 1),
        "final_log10_cpml": round(float(np.log10(max(ViralLoadMapper()
                               .normalized_to_copies(V[-1]), 1.0))), 2),
        "peak_teff": round(float(T_eff.max()), 3),
        "rebound_index": round(rebound_index, 3),
        "final_frr": round(float(F_rr[-1]), 3),
        # Model-internal urine:plasma ratio — a well-defined quantity
        # with no unit-conversion assumptions (Funk 2008: urine >> plasma).
        "urine_plasma_ratio": round(float(V_u[-1] / max(V[-1], 1e-12)), 1),
        "schedule": [(s, e, tr) for s, e, tr in schedule],
        "sir_trough_ngml": sir_trough,
    }


def candidate_schedules():
    """Baseline 8 ng/mL tac maintenance, step-down to a lower trough at a
    chosen week — the shape of real reduction protocols — plus tac->sir
    conversions."""
    scenarios = []
    for step_week in (2, 4, 6):
        for trough in (3.0, 4.0, 5.0, 6.0):
            scenarios.append((
                [(0.0, step_week * 7.0, 8.0), (step_week * 7.0, None, trough)],
                None,
                f"tac 8 -> {trough:g} ng/mL at wk{step_week}",
            ))
    for trough in (4.0, 6.0):
        scenarios.append((
            [(0.0, 28.0, 8.0), (28.0, None, 3.0)],
            trough,
            f"tac 8 -> 3 ng/mL + sir {trough:g} ng/mL at wk4",
        ))
    scenarios.append(([(0.0, None, 8.0)], None, "reference: hold tac 8 ng/mL"))
    return scenarios


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true",
                        help="also write results to outputs/reduction_schedules.json")
    args = parser.parse_args()

    print("BKPyV immunosuppression-reduction schedule sweep")
    print("(mechanistic hypothesis generator — assumption-labelled, not patient-fitted)\n")
    header = (f"{'schedule':<42} {'clear@wk':>9} {'final log10':>11} "
              f"{'peakTeff':>9} {'rebound':>8} {'F_rr':>6} {'u:p':>6}")
    print(header)
    print("-" * len(header))
    results = []
    for schedule, sir, label in candidate_schedules():
        m = evaluate_schedule(schedule, sir_trough=sir)
        m["label"] = label
        results.append(m)
        clear = "never" if m["clearance_weeks"] is None else f"{m['clearance_weeks']:.1f}"
        # plasma -> ~0 under clearing schedules makes u:p diverge; cap display
        up = ">1e4" if m["urine_plasma_ratio"] > 1e4 else f"{m['urine_plasma_ratio']:.0f}"
        print(f"{label:<42} {clear:>9} {m['final_log10_cpml']:>11} "
              f"{m['peak_teff']:>9} {m['rebound_index']:>8} {m['final_frr']:>6} "
              f"{up:>6}")

    # Pareto view: schedules that clear, ranked by lowest rebound first
    clearing = [r for r in results if r["clearance_weeks"] is not None]
    clearing.sort(key=lambda r: (r["rebound_index"], r["clearance_weeks"]))
    print("\nPareto front (clearing schedules, lowest rebound first):")
    for r in clearing[:5]:
        print(f"  {r['label']:<42} clear wk{r['clearance_weeks']:<6} "
              f"rebound {r['rebound_index']}")

    if args.json:
        out = Path("outputs")
        out.mkdir(exist_ok=True)
        (out / "reduction_schedules.json").write_text(json.dumps(results, indent=2))
        print(f"\nWrote {out/'reduction_schedules.json'}")


if __name__ == "__main__":
    main()
