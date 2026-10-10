#!/usr/bin/env python3
"""Early-window forecast: a digital-twin proof of concept.

Clinical question: at first detection of DNAemia, can we predict WHEN a
patient clears under a given intervention — from only the first three
weekly qPCR points? Method: a virtual cohort (each patient is a
parameter-jittered individual: beta, delta, production, inoculum) gives
the true trajectory; the forecaster sees only the first three weekly
measurements above baseline, perturbed by realistic qPCR noise
(0.25 log10), infers (beta, delta) on a coarse joint grid, then
forward-simulates under the intervention and reports predicted vs true
clearance week.

Honest by construction: the forecaster uses a DIFFERENT parameterization
than the generator's full jitter set (it infers only beta+delta while
production and inoculum also vary), and reports error, not just numbers.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402
from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402

HORIZON = 180.0
MAPPER = ViralLoadMapper()
QPCR_NOISE_LOG10 = 0.25   # realistic assay variability
OBS_DAYS = (4, 11, 18)    # first three weekly checks post-detection
ACT_DAY = 21.0            # intervention decision point (wk 3)


def _simulate(ode, y0, dosing, horizon=HORIZON, dt=1.0):
    t_eval = np.arange(0, horizon + 1e-9, dt)
    sol = solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                    (0, horizon), y0, t_eval=t_eval, method="LSODA")
    cp = np.array([MAPPER.normalized_to_copies(float(v))
                   for v in np.maximum(sol.y[0], 1e-9)])
    return t_eval, cp, sol


def _clearance_week(t, cp):
    below = cp < 1_000.0
    for i in range(len(cp)):
        if below[i:].all():
            return t[i] / 7.0
    return np.inf


def _conversion(day):
    return {"tacrolimus": [{"start": 0.0, "stop": day, "trough_ng_ml": 8.0},
                           {"start": day, "stop": None, "trough_ng_ml": 3.0}],
            "sirolimus": [{"start": day, "stop": None, "trough_ng_ml": 4.0}]}


def forecast_patient(obs_days, obs_cp, act_day=ACT_DAY):
    """Infer (beta, delta) from the noisy early points, predict clearance."""
    best = None
    for b in np.linspace(0.18, 0.5, 9):
        for d in np.linspace(0.28, 0.55, 7):
            ode = BKPyVODESystem(params={"beta": float(b),
                                         "delta": float(d)})
            y0 = ode.get_infection_conditions(0.1)
            tt, cp, _ = _simulate(
                ode, y0,
                {"tacrolimus": [{"start": 0, "stop": None,
                                 "trough_ng_ml": 8.0}]},
                horizon=float(obs_days[-1]) + 1)
            pred = np.interp(obs_days, tt, cp)
            err = np.mean((np.log10(np.maximum(pred, 1)) -
                           np.log10(np.maximum(obs_cp, 1))) ** 2)
            if best is None or err < best[0]:
                best = (err, b, d)
    ode = BKPyVODESystem(params={"beta": best[1], "delta": best[2]})
    y0 = ode.get_infection_conditions(0.1)
    tt, cp, _ = _simulate(ode, y0, _conversion(act_day))
    return (best[1], best[2]), _clearance_week(tt, cp)


def virtual_cohort(n=24, seed=11, noise=QPCR_NOISE_LOG10):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        truth = {"beta": float(rng.uniform(0.2, 0.5)),
                 "delta": float(rng.uniform(0.3, 0.55)),
                 "p": float(rng.uniform(6.0, 11.0))}
        inoc = float(rng.uniform(0.05, 0.2))
        ode = BKPyVODESystem(params=truth)
        y0 = ode.get_infection_conditions(inoc)
        tt, cp_true, _ = _simulate(ode, y0, _conversion(ACT_DAY))
        # natural-history obs: same patient trajectory before ACT_DAY but
        # under steady tac 8 (the conversion only kicks at ACT_DAY anyway)
        obs = cp_true[np.searchsorted(tt, np.array(OBS_DAYS))]
        obs = 10 ** (np.log10(np.maximum(obs, 1.0))
                     + rng.normal(0, noise, size=len(obs)))
        (b_hat, d_hat), pred_wk = forecast_patient(OBS_DAYS, obs)
        true_wk = _clearance_week(tt, cp_true)
        rows.append({"beta_true": truth["beta"], "beta_hat": round(b_hat, 3),
                     "true_clearance_weeks": (None if np.isinf(true_wk)
                                              else round(true_wk, 1)),
                     "pred_clearance_weeks": (None if np.isinf(pred_wk)
                                              else round(pred_wk, 1))})
    return rows


def main():
    print("Early-window forecast — virtual cohort (n=24, 3 noisy weekly points)")
    print("(mechanistic hypothesis generator — assumption-labelled)\n")
    rows = virtual_cohort()
    paired = [r for r in rows
              if r["pred_clearance_weeks"] is not None
              and r["true_clearance_weeks"] is not None]
    errs = [abs(r["pred_clearance_weeks"] - r["true_clearance_weeks"])
            for r in paired]
    discordant = sum(1 for r in rows
                     if (r["pred_clearance_weeks"] is None)
                     != (r["true_clearance_weeks"] is None))
    print(f"{'beta true':>9} {'beta hat':>9} {'true wk':>8} {'pred wk':>8}")
    print("-" * 38)
    for r in rows[:12]:
        tw = "never" if r["true_clearance_weeks"] is None else r["true_clearance_weeks"]
        pw = "never" if r["pred_clearance_weeks"] is None else r["pred_clearance_weeks"]
        print(f"{r['beta_true']:>9.3f} {r['beta_hat']:>9.3f} {tw:>8} {pw:>8}")
    print(f"... ({len(rows)} patients)")
    if errs:
        print(f"\nClearance-week forecast MAE: {np.mean(errs):.2f} wk "
              f"(max {np.max(errs):.1f}) on {len(errs)} resolved pairs")
    print(f"Clear/not-clear discordance: {discordant}/{len(rows)}")


if __name__ == "__main__":
    main()
