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
from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402

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


def forecast_posterior(obs_days, obs_cp, act_day=ACT_DAY,
                       n_mcmc=1500, n_forward=40, seed=3):
    """Full uncertainty quantification: Metropolis posterior over
    (beta, delta) given the noisy early points, propagated through the
    conversion simulation -> clearance-week DISTRIBUTION.

    Returns dict(median_weeks, lo90, hi90, clear_prob, n_forward). Unlike
    `forecast_patient` (point estimate on a grid), this admits how much
    the data constrain the prediction — the honest answer a clinician
    needs ("clear by week X with 90% probability").

    Likelihood: Gaussian in log10 cp/mL with the assay-noise sd; uniform
    prior over the plausible (beta, delta) box. Random-walk proposals.
    Each retained posterior draw is forward-simulated under conversion.
    `clear_prob` = fraction of draws clearing within HORIZON; never-
    clearing draws contribute HORIZON+ to the interval (conservative).
    """
    rng = np.random.default_rng(seed)
    lo, hi = np.array([0.18, 0.28]), np.array([0.50, 0.55])

    def neglogpost(theta):
        if np.any(theta < lo) or np.any(theta > hi):
            return np.inf
        ode = BKPyVODESystem(params={"beta": float(theta[0]),
                                     "delta": float(theta[1])})
        y0 = ode.get_infection_conditions(0.1)
        tt, cp, _ = _simulate(
            ode, y0, {"tacrolimus": [{"start": 0, "stop": None,
                                      "trough_ng_ml": 8.0}]},
            horizon=float(obs_days[-1]) + 1)
        pred = np.interp(obs_days, tt, cp)
        sse = np.mean((np.log10(np.maximum(pred, 1))
                       - np.log10(np.maximum(obs_cp, 1))) ** 2)
        return sse / (2 * QPCR_NOISE_LOG10 ** 2)

    theta = np.array([0.35, 0.42])
    lp = neglogpost(theta)
    chain = []
    step = np.array([0.02, 0.02])
    for _ in range(n_mcmc):
        prop = theta + rng.normal(0, step)
        lp_prop = neglogpost(prop)
        if np.log(rng.uniform()) < lp - lp_prop:
            theta, lp = prop, lp_prop
        chain.append(theta.copy())
    chain = np.array(chain[n_mcmc // 3:])          # burn-in
    draws = chain[rng.choice(len(chain), size=min(n_forward, len(chain)),
                             replace=False)]
    weeks = []
    for b, d in draws:
        ode = BKPyVODESystem(params={"beta": float(b), "delta": float(d)})
        y0 = ode.get_infection_conditions(0.1)
        tt, cp, _ = _simulate(ode, y0, _conversion(act_day))
        weeks.append(_clearance_week(tt, cp))
    weeks = np.array(weeks)
    finite = weeks[np.isfinite(weeks)]
    clear_prob = len(finite) / len(weeks)
    fill = np.where(np.isfinite(weeks), weeks, HORIZON / 7.0 + 1)
    return {"median_weeks": float(np.median(fill)),
            "lo90": float(np.quantile(fill, 0.05)),
            "hi90": float(np.quantile(fill, 0.95)),
            "clear_prob": clear_prob,
            "n_forward": len(weeks)}


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

    print("\nPosterior forecast (Metropolis UQ) for first 2 patients:")
    rng = np.random.default_rng(11)
    for i in range(2):
        truth = {"beta": float(rng.uniform(0.2, 0.5)),
                 "delta": float(rng.uniform(0.3, 0.55)),
                 "p": float(rng.uniform(6.0, 11.0))}
        rng.uniform(0.05, 0.2)  # keep the rng stream aligned with cohort
        ode = BKPyVODESystem(params=truth)
        y0 = ode.get_infection_conditions(0.1)
        tt, cp_true, _ = _simulate(ode, y0, _conversion(ACT_DAY))
        obs = cp_true[np.searchsorted(tt, np.array(OBS_DAYS))]
        obs = 10 ** (np.log10(np.maximum(obs, 1.0))
                     + rng.normal(0, QPCR_NOISE_LOG10, size=len(obs)))
        r = forecast_posterior(OBS_DAYS, obs)
        tw = _clearance_week(tt, cp_true)
        print(f"  patient {i}: truth {tw:.1f} wk — posterior median "
              f"{r['median_weeks']:.1f} wk, 90% CI [{r['lo90']:.1f}, "
              f"{r['hi90']:.1f}], P(clear) {r['clear_prob']:.0%}")


if __name__ == "__main__":
    main()
