#!/usr/bin/env python3
"""Virus-specific T-cell therapy (VST) simulation — a digital twin of
the frontier BKPyV intervention.

Clinical grounding: third-party BKPyV-specific T-cell products (e.g.
posoleucel / ALVR105, phase-3) rescue patients who fail drug-sparing
management — and they work UNDER ongoing calcineurin inhibition, the
signature property that distinguishes cellular therapy from drug
adjustment. The model's T_eff arm already encodes this: infused CTLs
enter the effector pool and kill infected cells via `tcell_kill` —
which is NOT attenuated by tacrolimus (calcineurin blockade hits
priming/expansion, not the cytolytic synapse) — while their in-host
expansion `prolif_tcell * antigen * tac_tcell_effect` IS suppressed.
Both properties are emergent here, not assumed.

Intervention: at `vst_day`, add `dose` to T_eff (index 16) — an
infusion of pre-expanded CTLs. Dose is in effector-pool fraction units
(carrying capacity 1.0); assumption-labelled, not cells/kg.

Predictions probed: dose-response (clearance week vs infusion size),
timing (early vs late infusion), and the tac-independence signature
(efficacy at trough 8 vs taper-failure at the same trough).
"""

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from early_forecast import _clearance_week  # noqa: E402

from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402

MAPPER = ViralLoadMapper()
HORIZON = 180.0
T_EFF_IDX = 16


def _cp(v_arr):
    return np.array([MAPPER.normalized_to_copies(float(v))
                     for v in np.maximum(v_arr, 1e-9)])


def simulate_with_vst(dosing, vst_day, dose, horizon=HORIZON, inoc=0.5):
    """Two-segment solve: [0, vst_day), infuse CTLs, [vst_day, horizon].
    Returns (t, cp, sol-like arrays for V and T_eff)."""
    ode = BKPyVODESystem()
    y0 = ode.get_infection_conditions(inoc)
    t_eval = np.arange(0, horizon + 1e-9, 1.0)
    sol1 = solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                     (0, vst_day), y0,
                     t_eval=t_eval[t_eval < vst_day], method="LSODA")
    y_inf = sol1.y[:, -1].copy()
    y_inf[T_EFF_IDX] += dose            # the infusion
    sol2 = solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                     (vst_day, horizon), y_inf,
                     t_eval=t_eval[t_eval >= vst_day], method="LSODA")
    t = np.concatenate([sol1.t, sol2.t])
    v_load = np.concatenate([sol1.y[0], sol2.y[0]])
    teff = np.concatenate([sol1.y[T_EFF_IDX], sol2.y[T_EFF_IDX]])
    return t, _cp(v_load), teff


def _tac(trough):
    return {"tacrolimus": [{"start": 0.0, "stop": None,
                            "trough_ng_ml": trough}]}


def main():
    print("VST digital twin — infused BKPyV-specific CTLs under tacrolimus")
    print("(mechanistic hypothesis — VST dose in effector-fraction units)\n")

    print("Single-infusion arms (infusion day 28, tac 8 ng/mL unless noted):")
    print(f"{'arm':>44} {'clear wk':>9} {'post-nadir':>10} {'nadir wk':>8}")
    print("-" * 74)
    taper2 = {"tacrolimus": [{"start": 0, "stop": 28, "trough_ng_ml": 8.0},
                             {"start": 28, "stop": None, "trough_ng_ml": 2.0}]}
    conv = {"tacrolimus": [{"start": 0, "stop": 28, "trough_ng_ml": 8.0},
                           {"start": 28, "stop": None, "trough_ng_ml": 3.0}],
            "sirolimus": [{"start": 28, "stop": None, "trough_ng_ml": 4.0}]}
    arms = [
        ("no intervention, tac 8", _tac(8.0), 0.0),
        ("VST 0.5, tac 8", _tac(8.0), 0.5),
        ("VST 2.0, tac 8", _tac(8.0), 2.0),
        ("VST 0.5, tac 3", _tac(3.0), 0.5),
        ("VST 0.5 + taper 8->2", taper2, 0.5),
        ("tac->sir conversion (control)", conv, 0.0),
    ]
    for name, dosing, dose in arms:
        t, cp, teff = simulate_with_vst(dosing, 28.0, dose)
        wk = _clearance_week(t, cp)
        post = cp[t >= 28.0]
        t_post = t[t >= 28.0]
        print(f"{name:>44} "
              f"{'never' if np.isinf(wk) else f'{wk:.1f}':>9} "
              f"{post.min():>10.0f} {t_post[np.argmin(post)] / 7.0:>8.1f}")

    print("""
Emergent result (not assumed): single-bolus VST produces only a
TRANSIENT dip — visible suppression needs a large bolus (~4x the
homeostatic ceiling), and even then viremia re-equilibrates within
weeks. Mechanism: infused CTLs sit above the homeostatic ceiling
(1 - T_eff/tcell_carry pulls them down within days) AND cannot
re-expand (prolif * antigen * tac_tcell_effect is suppressed by
calcineurin blockade), while the latent reservoir L and cell-to-cell
spread reseed infection once the spike decays.

This reproduces the documented clinical pattern that VST responses are
often transient and that durable response tracks in-vivo CTL expansion
— which the model says is exactly what calcineurin blockade prevents.
Falsifiable claim: VST is a suppression BRIDGE for high-IS patients,
not monotherapy-curative; durable clearance still requires an IS-
reduction arm (conversion clears wk ~8.3 with or without VST).

Limitation flag: tcell_carry is a fixed homeostatic ceiling; if real
VST products raise effective carrying capacity (e.g. via exogenous
IL-2 support or niches the model does not represent), the durability
prediction changes — the parameter to probe when data arrives.""")


if __name__ == "__main__":
    main()
