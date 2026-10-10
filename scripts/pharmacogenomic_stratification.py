#!/usr/bin/env python3
"""Pharmacogenomic stratification: same dose, different outcome.

CYP3A5 expressors metabolize tacrolimus ~2x faster than non-expressors
(documented ~1.5-2x dose requirement; CPIC guideline): identical mg/kg
dosing produces systematically different troughs by genotype. This
script maps genotype -> effective trough -> outcomes through the real
PK layer: reactivation-onset distributions and clearance weeks under
the conversion policy.

Emergent point: a dosing protocol that is 'safe' for a non-expressor
can be under-immunosuppressive in efficacy terms (graft rejection) yet
LOWER viral burden, while an expressor on the same nominal dose runs
hot — the model separates efficacy risk from viral risk per genotype,
which is exactly the tension a transplant team balances. This is the
personalized-medicine layer almost no BKPyV model carries.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reactivation_onset import onset_days  # noqa: E402
from early_forecast import _simulate, _clearance_week, _conversion  # noqa: E402
from vcm.simulators.ode_system import BKPyVODESystem  # noqa: E402

# CPIC: CYP3A5 expressors (*1 carriers) clear tacrolimus faster ->
# effective trough ~0.5-0.7x the nominal at the same mg/kg dose.
GENOTYPES = {
    "CYP3A5 non-expressor (*3/*3)": {"trough_factor": 1.0},
    "CYP3A5 expressor (*1 carrier)": {"trough_factor": 0.6},
}
NOMINAL_TROUGH = 8.0


def effective_trough(genotype, nominal=NOMINAL_TROUGH):
    return nominal * GENOTYPES[genotype]["trough_factor"]


def genotype_outcome(genotype, nominal=NOMINAL_TROUGH, n=120):
    """Onset distribution + clearance under a uniform mg/kg protocol."""
    trough = effective_trough(genotype, nominal)
    onset = onset_days(trough, n=n, seed=7)
    ode = BKPyVODESystem()
    y0 = ode.get_infection_conditions(0.1)
    tt, cp, _ = _simulate(ode, y0, _conversion(21.0))
    return {"genotype": genotype,
            "effective_trough": trough,
            "median_onset_weeks": onset["median_onset_weeks"],
            "pct_reactivated": onset["pct_reactivated"],
            "clearance_weeks": (None
                                if np.isinf(_clearance_week(tt, cp))
                                else round(_clearance_week(tt, cp), 1))}


def main():
    print("Pharmacogenomic stratification — uniform mg/kg dosing")
    print("(mechanistic hypothesis generator — assumption-labelled)\n")
    hdr = (f"{'genotype':<36} {'eff trough':>10} {'% onset':>8} "
           f"{'med onset wk':>12} {'clear wk':>9}")
    print(hdr)
    print("-" * len(hdr))
    for g in GENOTYPES:
        m = genotype_outcome(g)
        med = "-" if m["median_onset_weeks"] is None else f"{m['median_onset_weeks']}"
        clr = "never" if m["clearance_weeks"] is None else f"{m['clearance_weeks']}"
        print(f"{g:<36} {m['effective_trough']:>10.1f} "
              f"{m['pct_reactivated']:>8} {med:>12} {clr:>9}")
    print("\nImplication: the risk assignment INVERTS by axis —"
          " non-expressors running full troughs carry the earlier viral"
          " onset, while expressors' danger is under-immunosuppression"
          " (rejection), not viremia. Genotype-blind dosing can't serve"
          " both; CPIC-guided starting doses plus trough-targeted"
          " monitoring is the model's mechanistic recommendation.")


if __name__ == "__main__":
    main()
