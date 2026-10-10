#!/usr/bin/env python3
"""Stochastic reactivation model: WHEN does post-transplant viremia start?

Clinically, BKPyV DNAemia doesn't begin at transplant day 0 — it emerges
from latently infected cells that reactivate under immunosuppression,
and onset clusters ~weeks 4-16 (median ~8; consistent with the cohort
patterns motivating Kotton's monthly screening start). No published
BKPyV kinetic model generates an onset-time DISTRIBUTION; this script
produces one mechanistically.

Model: each reactivation draw samples a reactivation day from an
exponential hazard that scales with tacrolimus trough (deeper
immunosuppression -> earlier reactivation). A reactivated cell injects
a small inoculum; the canonical 23-dim ODE then carries it to the 1k
cp/mL screening threshold — onset day = threshold crossing. Monte-Carlo
over draws gives the onset distribution per drug level.

Emergent predictions: higher trough -> earlier AND tighter onset;
the distribution's right tail explains late-onset cases under
immunosuppression reduction.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402
from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402

ONSET_CP_ML = 1_000.0
H_MAX = 180.0          # observation window (days)
INOCULUM = 0.05        # a single reactivated focus: well below threshold
# Hazard model: base rate on minimal immunosuppression + linear boost per
# ng/mL above a 3 ng/mL maintenance floor. Chosen so median onset lands
# inside the clinical 4-16 wk cluster at tac ~8 ng/mL.
H_BASE = 0.005         # latent reactivation hazard at trough ~3 ng/mL (1/day)
H_PER_NGML = 0.0025    # hazard gain per ng/mL of tacrolimus trough


def hazard(tac_ngml):
    return H_BASE + H_PER_NGML * max(0.0, tac_ngml - 3.0)


def onset_days(tac_ngml, n=400, seed=7, horizon=H_MAX):
    """Monte-Carlo onset distribution at a constant tacrolimus trough."""
    rng = np.random.default_rng(seed)
    mapper = ViralLoadMapper()
    h = hazard(tac_ngml)
    dosing = {"tacrolimus": [{"start": 0.0, "stop": None,
                              "trough_ng_ml": tac_ngml}]}
    onsets = []
    for _ in range(n):
        t_react = rng.exponential(1.0 / h)
        if t_react > horizon:
            onsets.append(np.inf)
            continue
        ode = BKPyVODESystem()
        y0 = ode.get_infection_conditions(INOCULUM)
        t_eval = np.linspace(t_react, horizon, max(2, int((horizon - t_react) * 2) + 1))
        sol = solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                        (t_react, horizon), y0, t_eval=t_eval, method="LSODA")
        cp = np.array([mapper.normalized_to_copies(float(v))
                       for v in np.maximum(sol.y[0], 1e-9)])
        crossed = t_eval[cp >= ONSET_CP_ML]
        onsets.append(float(crossed[0]) if len(crossed) else np.inf)
    onsets = np.asarray(onsets)
    finite = onsets[np.isfinite(onsets)]
    return {
        "tac_ngml": tac_ngml,
        "hazard_per_day": round(h, 4),
        "pct_reactivated": round(100.0 * len(finite) / n, 1),
        "median_onset_weeks": (None if not len(finite)
                               else round(float(np.median(finite)) / 7.0, 1)),
        "p10_weeks": (None if not len(finite)
                      else round(float(np.percentile(finite, 10)) / 7.0, 1)),
        "p90_weeks": (None if not len(finite)
                      else round(float(np.percentile(finite, 90)) / 7.0, 1)),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=200)
    args = ap.parse_args()
    print("Reactivation-onset distribution vs tacrolimus trough")
    print("(mechanistic hypothesis generator — assumption-labelled)\n")
    hdr = (f"{'tac ng/mL':>10} {'hazard/d':>9} {'% onset':>8} "
           f"{'p10 wk':>7} {'med wk':>7} {'p90 wk':>7}")
    print(hdr)
    print("-" * len(hdr))
    for tac in (3.0, 5.0, 8.0, 12.0):
        m = onset_days(tac, n=args.n)
        med = "-" if m["median_onset_weeks"] is None else f"{m['median_onset_weeks']}"
        p10 = "-" if m["p10_weeks"] is None else f"{m['p10_weeks']}"
        p90 = "-" if m["p90_weeks"] is None else f"{m['p90_weeks']}"
        print(f"{tac:>10} {m['hazard_per_day']:>9} {m['pct_reactivated']:>8} "
              f"{p10:>7} {med:>7} {p90:>7}")


if __name__ == "__main__":
    main()
