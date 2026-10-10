"""Fit the 23-dim BKPyV ODE to REAL longitudinal patient series.

Data: digitized plasma BKV trajectories from Funk et al. 2008
(validation_data/funk2008/, CC BY-NC-ND; see that folder's README).

Method (honest construction): the patient's week-0 is NOT the infection
onset — it is a mid-course observation. We therefore pre-roll the model
under high immunosuppression so the joint state sits on the model's own
quasi-steady manifold when observation begins, then let five parameters
vary: (beta, d_infected, pre_days, act_week, tac_lo) where the
intervention is an immunosuppression reduction (tac trough T_HI -> tac_lo)
at act_week — the clinically documented response for these patients.

Two scores per patient:
  * full fit: Nelder-Mead on all points -> RMSE / R2 in log10 space
  * holdout: fit on the first 60% of points, predict the remaining 40%
    -> MAE vs the naive baseline (per-patient mean of train points),
    the falsifiable "does early data predict later course" test.

This is the model's first contact with real patient data — residuals are
reported, not tuned away.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402

MAPPER = ViralLoadMapper()
DATA_DIR = (Path(__file__).resolve().parent.parent
            / "validation_data" / "funk2008")
T_HI = 9.0          # high-IS tac trough before intervention (ng/mL)
INOC = 0.5          # fixed inoculum; pre-roll absorbs its scale


def _log_cp(v_arr):
    cp = np.array([MAPPER.normalized_to_copies(float(v))
                   for v in np.maximum(v_arr, 1e-9)])
    return np.log10(np.maximum(cp, 1.0))


def load_series(path):
    """Return (weeks post-transplant, log10 cp/mL) as stored in the CSV.
    Callers convert weeks -> days exactly once (x7) where needed."""
    df = pd.read_csv(path)
    return df["week"].to_numpy(), df["plasma_log10_cp_ml"].to_numpy()


def simulate_patient(theta, obs_days):
    """theta = (beta, delta, pre_days, act_week, tac_lo).
    Returns log10 cp/mL at obs_days (post-transplant timeline)."""
    theta = np.asarray(theta, dtype=float)
    if np.any(~np.isfinite(theta)):
        return None
    beta, delta, pre_days, act_wk, tac_lo = theta
    act_day = pre_days + act_wk * 7.0
    end_day = pre_days + obs_days[-1] + 14.0
    if not (np.isfinite(end_day) and end_day > pre_days >= 0):
        return None
    ode = BKPyVODESystem(params={"beta": float(beta),
                                 "delta": float(delta)})
    y0 = ode.get_infection_conditions(INOC)
    dosing = {"tacrolimus": [{"start": 0.0, "stop": act_day,
                              "trough_ng_ml": T_HI},
                             {"start": act_day, "stop": None,
                              "trough_ng_ml": tac_lo}]}
    eval_days = np.unique(pre_days + obs_days)
    try:
        sol = solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                        (0.0, end_day), y0, method="LSODA",
                        t_eval=eval_days)
    except (ValueError, RuntimeError):
        return None
    if not sol.success or sol.y.shape[1] != len(eval_days):
        return None
    pred = np.interp(pre_days + obs_days, eval_days, sol.y[0])
    return _log_cp(pred)


BOUNDS_LO = np.array([0.10, 0.10, 20.0, 0.0, 0.0])


def _sse(theta, obs_days, obs_log, hi=None):
    theta = np.clip(np.asarray(theta, dtype=float), BOUNDS_LO, hi)
    pred = simulate_patient(theta, obs_days)
    if pred is None or np.any(~np.isfinite(pred)):
        return 1e6
    return float(np.sum((pred - obs_log) ** 2))


def fit_series(obs_weeks, obs_log, n_starts=3, seed=7):
    """Nelder-Mead fit of (beta, delta, pre_days, act_week, tac_lo)."""
    obs_days = obs_weeks * 7.0
    bounds_lo = BOUNDS_LO
    bounds_hi = np.array([1.50, 1.50, 400.0, obs_weeks[-1] + 4.0, 9.0])
    rng = np.random.default_rng(seed)
    starts = [np.array([0.6, 0.40, 120.0, obs_weeks[-1] * 0.5, 3.0])]
    for _ in range(n_starts - 1):
        starts.append(bounds_lo + rng.random(5) * (bounds_hi - bounds_lo))
    best = None
    for s0 in starts:
        res = minimize(_sse, s0, args=(obs_days, obs_log, bounds_hi),
                       method="Nelder-Mead",
                       options={"maxiter": 250, "xatol": 1e-3,
                                "fatol": 1e-4})
        x = np.clip(res.x, bounds_lo, bounds_hi)
        val = _sse(x, obs_days, obs_log, bounds_hi)
        if best is None or val < best[0]:
            best = (val, x)
    sse, theta = best
    pred = simulate_patient(theta, obs_days)
    rmse = float(np.sqrt(sse / len(obs_log)))
    ss_res = sse
    ss_tot = float(np.sum((obs_log - obs_log.mean()) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-9)
    return {"theta": theta, "rmse": rmse, "r2": r2, "pred": pred}


def holdout_score(obs_weeks, obs_log, frac=0.6, seed=7):
    """Fit on first frac of points, MAE on the rest vs mean-baseline."""
    n_train = max(3, int(np.ceil(len(obs_weeks) * frac)))
    tr_wk, te_wk = obs_weeks[:n_train], obs_weeks[n_train:]
    tr_y, te_y = obs_log[:n_train], obs_log[n_train:]
    fit = fit_series(tr_wk, tr_y, n_starts=2, seed=seed)
    # extrapolate: pre_days/act_week fitted on train stay absolute in
    # patient time, so simulate the fitted patient through test days
    theta = fit["theta"].copy()
    pred_te = simulate_patient(theta, te_wk * 7.0)
    mae = float(np.mean(np.abs(pred_te - te_y)))
    base_mae = float(np.mean(np.abs(te_y - tr_y.mean())))
    return {"n_train": n_train, "n_test": len(te_wk), "mae": mae,
            "baseline_mae": base_mae, "theta": theta}


def main():
    print("Fitting 23-dim ODE to digitized Funk-2008 patient plasma series")
    print("(real-data holdout validation — residuals reported honestly)\n")
    print(f"{'patient':>10} {'fullRMSE':>9} {'fullR2':>7} "
          f"{'holdMAE':>8} {'baseMAE':>8}  fitted (beta,delta,pre,actwk,taclo)")
    print("-" * 92)
    results = {}
    for csv in sorted(DATA_DIR.glob("pat_*.csv")):
        name = csv.stem.replace("pat_", "")
        wk, y = load_series(csv)
        full = fit_series(wk, y)
        ho = holdout_score(wk, y)
        results[name] = (full, ho)
        th = full["theta"]
        print(f"{name:>10} {full['rmse']:>9.3f} {full['r2']:>7.3f} "
              f"{ho['mae']:>8.3f} {ho['baseline_mae']:>8.3f}  "
              f"({th[0]:.2f},{th[1]:.3f},{th[2]:.0f},{th[3]:.1f},{th[4]:.1f})")
    mae_all = np.mean([r[1]["mae"] for r in results.values()])
    base_all = np.mean([r[1]["baseline_mae"] for r in results.values()])
    print("-" * 92)
    print(f"cohort mean holdout MAE {mae_all:.3f} log10 vs baseline "
          f"{base_all:.3f} log10")
    _plot_fits(results)


def _plot_fits(results):
    """6-panel figure: observed points + full-fit curve; holdout split
    marked. Saved to outputs/figures/patient_fits.png."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable — skipping patient_fits figure")
        return
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=True)
    for ax, (csv, (name, (full, ho))) in zip(
            axes.ravel(),
            sorted(zip(sorted(DATA_DIR.glob("pat_*.csv")),
                       results.items()))):
        wk, y = load_series(csv)
        theta = full["theta"]
        grid_wk = np.linspace(0, wk[-1] + 4, 300)
        pred = simulate_patient(theta, grid_wk * 7.0)
        ax.plot(grid_wk, pred, lw=1.5, label="model fit")
        ntr = ho["n_train"]
        ax.plot(wk[:ntr], y[:ntr], "o", ms=4, label="observed (train)")
        ax.plot(wk[ntr:], y[ntr:], "s", ms=4, mfc="none",
                label="observed (holdout)")
        ax.axvspan(wk[ntr - 1], wk[-1] + 4, alpha=0.08)
        ax.set_title(f"Pat {name}  R2={full['r2']:.2f}  "
                     f"holdout {ho['mae']:.2f} vs {ho['baseline_mae']:.2f}",
                     fontsize=9)
        ax.set_ylim(0, 11)
    for ax in axes[-1]:
        ax.set_xlabel("weeks post-transplant")
    for ax in axes[:, 0]:
        ax.set_ylabel("log10 cp/mL")
    axes[0, 0].legend(fontsize=7, loc="upper right")
    fig.suptitle("ODE fits to digitized patient plasma series "
                 "(Funk 2008); shaded = holdout region", fontsize=11)
    fig.tight_layout()
    out = (Path(__file__).resolve().parent.parent
           / "outputs" / "figures" / "patient_fits.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170)
    print(f"figure -> {out}")


if __name__ == "__main__":
    main()
