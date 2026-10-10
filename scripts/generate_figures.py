#!/usr/bin/env python3
"""Publication-quality figure generation for poster/presentation.

Produces a four-panel composite telling the whole mechanistic story:
  A) Plasma DNAemia + urinary load vs screening thresholds
  B) Two-pool NCCR quasi-species emergence (kidney vs urine)
  C) Regimen comparison (taper vs conversion)
  D) Reactivation onset distribution by tacrolimus trough

Assumption labels are burned into each panel where applicable.
Outputs to outputs/figures/ (gitignored directory — figures are
reproducible artifacts, run this script to regenerate).
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from vcm.clinical.viral_load_mapper import ViralLoadMapper  # noqa: E402
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402

OUT = Path("outputs/figures")
MAPPER = ViralLoadMapper()


def _cp(v_arr):
    return np.array([MAPPER.normalized_to_copies(float(v))
                     for v in np.maximum(v_arr, 1e-9)])


def _run(dosing, horizon=180.0, inoc=0.5):
    ode = BKPyVODESystem()
    y0 = ode.get_infection_conditions(inoc)
    t = np.arange(0, horizon + 1e-9, 1.0)
    sol = solve_ivp(lambda tt, y: ode.ode_system(tt, y, dosing),
                    (0, horizon), y0, t_eval=t, method="LSODA")
    return t, sol.y


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    def wk(t):
        return t / 7.0

    # A) compartments vs thresholds (untreated)
    t, y = _run(None)
    cp = _cp(y[0])
    ax[0, 0].plot(wk(t), np.log10(cp), lw=2, label="plasma DNAemia")
    ax[0, 0].plot(wk(t), np.log10(np.maximum(y[19], 1e-9)),
                  lw=1.5, ls="--", label="urinary V_u (model units)")
    for thr, c in ((1e3, "tab:green"), (1e4, "tab:orange")):
        ax[0, 0].axhline(np.log10(thr), color=c, ls=":", lw=1,
                         label=f"{int(thr):,} cp/mL")
    ax[0, 0].set(ylabel="log10 load", xlabel="weeks",
                 title="A. Two-compartment kinetics vs clinical thresholds")
    ax[0, 0].legend(fontsize=8)

    # B) quasi-species emergence
    ax[0, 1].plot(wk(t), y[20], lw=2, label="kidney pool (plasma)")
    ax[0, 1].plot(wk(t), y[21], lw=2, label="urinary pool")
    ax[0, 1].set(ylabel="rearranged-NCCR fraction", xlabel="weeks",
                 title="B. rr-NCCR emergence: plasma > urine (Gosert 2008)")
    ax[0, 1].legend(fontsize=8)

    # C) regimens
    taper = {"tacrolimus": [{"start": 0, "stop": None, "trough_ng_ml": 4.0}]}
    conv = {"tacrolimus": [{"start": 0, "stop": 28, "trough_ng_ml": 8.0},
                           {"start": 28, "stop": None, "trough_ng_ml": 3.0}],
            "sirolimus": [{"start": 28, "stop": None, "trough_ng_ml": 4.0}]}
    for dosing, lab, col in ((taper, "tac taper 8->4", "tab:red"),
                             (conv, "tac->sir conversion", "tab:blue")):
        t2, y2 = _run(dosing)
        ax[1, 0].plot(wk(t2), np.log10(_cp(y2[0])), lw=2, color=col, label=lab)
    ax[1, 0].axhline(3, color="tab:green", ls=":", lw=1)
    ax[1, 0].set(ylabel="log10 cp/mL", xlabel="weeks",
                 title="C. Regimen design: taper never clears, conversion clears")
    ax[1, 0].legend(fontsize=8)

    # D) onset distributions (hazard model identical to
    # scripts/reactivation_onset.py; inlined to keep this script standalone)
    for tac, col in ((3.0, "tab:gray"), (8.0, "tab:red"), (12.0, "tab:purple")):
        rng = np.random.default_rng(7)
        ode = BKPyVODESystem()
        onsets = []
        for _ in range(150):
            tr = rng.exponential(1.0 / (0.005 + 0.0025 * max(0, tac - 3)))
            if tr > 180:
                continue
            y0 = ode.get_infection_conditions(0.05)
            tt = np.arange(tr, 180, 0.5)
            s = solve_ivp(lambda a, b: ode.ode_system(
                a, b, {"tacrolimus": [{"start": 0, "stop": None,
                                       "trough_ng_ml": tac}]}),
                (tr, 180), y0, t_eval=tt, method="LSODA")
            cpl = _cp(s.y[0])
            hit = tt[cpl >= 1e3]
            if len(hit):
                onsets.append(hit[0] / 7.0)
        if onsets:
            ax[1, 1].hist(onsets, bins=15, alpha=0.55, color=col,
                          label=f"tac {tac:.0f} ng/mL (med {np.median(onsets):.1f}wk)")
    ax[1, 1].set(xlabel="onset week", ylabel="patients",
                 title="D. Reactivation onset vs tacrolimus trough")
    ax[1, 1].legend(fontsize=8)

    fig.suptitle("BKPyV-VCM: mechanism-first virtual cell model "
                 "(mechanistic hypotheses — assumption-labelled)", fontsize=12)
    fig.tight_layout()
    out = OUT / "composite_figure.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
