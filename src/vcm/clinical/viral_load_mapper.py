"""Clinical Viral Load Mapper for BKPyV Plugin.

Bridges the model's dimensionless viral-load variable V (ODE units, unbounded)
to indicative plasma BKPyV-DNAemia in copies/mL.

IMPORTANT INTERPRETATION BOUNDARY
---------------------------------
There is no measured patient dataset in this repository to fit against, so this
bridge is an explicitly-labelled ASSUMPTION SET, not a calibration: the model's
V is mapped through piecewise log-linear interpolation between anchor points
that correspond to consensus clinical thresholds (1,000 and 10,000 copies/mL;
AST Infectious Diseases Community of Practice 2019 guideline, and Kotton et
al., Second International Consensus Guidelines, Transplantation 2024). Treat
the resulting copies/mL as an order-of-magnitude visualisation aid.

Previous versions used a Hill function "fitted" to anchors (0.3 -> 1,000,
0.6 -> 10,000, 1.0 -> 1e7). Those anchors were impossible for a saturating
Hill with Vmax = 1e7, so the shipped fit missed them by one to four orders of
magnitude while reporting the anchors as if they were reproduced. The current
interpolation scheme reproduces every anchor exactly, by construction.

Reference Sources:
- AST Infectious Diseases Community of Practice: Hirsch HH, Randhawa PS.
  BK polyomavirus in solid organ transplantation. Clin Transplant
  2019;33(9):e13528.
- Kotton CN, et al. The Second International Consensus Guidelines on BK
  polyomavirus in kidney transplantation. Transplantation 2024;108(9):1834-1866.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict


def hill_function(vl, Vmax, K, n):
    """Hill equation: copies_per_mL = Vmax * (vl^n) / (K^n + vl^n).

    DEPRECATED for the default mapping (kept for API compatibility and the
    existing test-suite): the historical Hill anchors were not jointly
    satisfiable, so the mapper now uses piecewise log-linear interpolation
    (see ``ViralLoadMapper``).
    """
    return Vmax * (vl**n) / (K**n + vl**n)


# Anchor points of the clinical bridge, in ODE viral-load units V -> copies/mL.
# These are ASSUMPTIONS chosen so the bridge passes through consensus clinical
# thresholds at biologically plausible points of the simulated trajectory;
# they are not fitted to patient data.
DEFAULT_ANCHORS = [
    (0.0, 0.0),       # no virus
    (0.02, 100.0),    # ~assay detection limit
    (0.2, 1000.0),    # screening threshold (AST IDCOP 2019 / Kotton 2024)
    (1.0, 10000.0),   # presumptive-PyVAN threshold (same guidelines)
    (3.0, 1e6),       # established viremia, order-of-magnitude
    (5.0, 1e7),       # high-level viremia, order-of-magnitude
]


def build_bridge_parameters():
    """Build the (deterministic) bridge-parameter record.

    The mapping is piecewise log-linear between anchors, so every anchor is
    reproduced exactly by construction. The record advertises that fact in its
    ``fit_quality`` block instead of pretending to be a statistical fit.
    """
    anchors_vl = np.array([a[0] for a in DEFAULT_ANCHORS], dtype=float)
    anchors_cp = np.array([a[1] for a in DEFAULT_ANCHORS], dtype=float)
    predicted = np.array([piecewise_loglinear(v, anchors_vl, anchors_cp) for v in anchors_vl])
    percent_error = np.where(
        anchors_cp > 0, np.abs(predicted - anchors_cp) / np.maximum(anchors_cp, 1e-12) * 100, 0.0
    )
    return {
        "model": "piecewise_log_linear",
        "anchors": [
            {"viral_load_v": float(v), "copies_per_ml": float(c)}
            for v, c in DEFAULT_ANCHORS
        ],
        "source": "Consensus thresholds: AST IDCOP 2019 (Hirsch & Randhawa, "
                  "Clin Transplant 33(9):e13528) and Kotton et al. 2024 "
                  "(Transplantation 108(9):1834-1866). Anchor placement on the "
                  "model V axis is an assumption, not a fit.",
        "fit_method": "piecewise log-linear interpolation through anchors (exact by construction)",
        "fit_quality": {
            "predicted_values": predicted.tolist(),
            "actual_values": anchors_cp.tolist(),
            "percent_error": [float(x) for x in percent_error],
        },
        "caveat": "Bridging function for qualitative visualisation; not a "
                  "patient-calibrated copy-number prediction.",
    }


def piecewise_loglinear(vl: float, anchors_vl: np.ndarray, anchors_cp: np.ndarray) -> float:
    """Interpolate log10(copies) linearly between anchors; flat below/above."""
    vl = max(0.0, float(vl))
    if vl <= anchors_vl[0]:
        return float(anchors_cp[0])
    if vl >= anchors_vl[-1]:
        # Extrapolate with the last segment's log-slope
        x0, x1 = anchors_vl[-2], anchors_vl[-1]
        y0 = np.log10(max(anchors_cp[-2], 1e-12))
        y1 = np.log10(max(anchors_cp[-1], 1e-12))
        y = y1 + (vl - x1) * (y1 - y0) / (x1 - x0)
        return float(10.0 ** y)
    for i in range(1, len(anchors_vl)):
        if vl <= anchors_vl[i]:
            x0, x1 = anchors_vl[i - 1], anchors_vl[i]
            c0, c1 = anchors_cp[i - 1], anchors_cp[i]
            if c0 <= 0.0 or c1 <= 0.0:
                # Linear (not log) interpolation when a segment touches zero
                frac = (vl - x0) / (x1 - x0)
                return float(c0 + frac * (c1 - c0))
            y0, y1 = np.log10(c0), np.log10(c1)
            frac = (vl - x0) / (x1 - x0)
            return float(10.0 ** (y0 + frac * (y1 - y0)))
    return float(anchors_cp[-1])


def save_hill_parameters(params: Dict, output_path: str):
    """Save bridge parameters to JSON (function name kept for compatibility)."""
    with open(output_path, 'w') as f:
        json.dump(params, f, indent=2)


def fit_hill_parameters():
    """DEPRECATED: delegating to the exact bridge.

    Kept so legacy callers still work; returns the piecewise bridge record
    (which contains ``model='piecewise_log_linear'``).
    """
    return build_bridge_parameters()


class ViralLoadMapper:
    """Mapper for converting simulator viral_load (0-1) to clinical plasma copies/mL.
    
    Uses a Hill function fitted to clinical anchor points to map normalized
    viral load from the BKPyV simulator to clinically relevant plasma viral load
    values in copies/mL.
    """
    
    CLINICAL_THRESHOLDS = {
        'screening': 1000.0,    # copies/mL - screening threshold
        'treatment': 10000.0,    # copies/mL - treatment threshold
        'severe': 1e7            # copies/mL - severe nephropathy
    }
    
    def __init__(self, params_path: Optional[str] = None):
        """Initialize ViralLoadMapper.

        Args:
            params_path: Path to a JSON bridge-parameter file. If the file is
                missing (or is an old record without ``model`` ==
                "piecewise_log_linear"), the default assumption anchors are
                used and saved.
        """
        if params_path is None:
            params_path = "data/processed/viral_load_mapper_params.json"

        params = None
        if Path(params_path).exists():
            with open(params_path, 'r') as f:
                params = json.load(f)
            if params.get("model") != "piecewise_log_linear":
                # Stale record from the old (broken) Hill fit — rebuild rather
                # than silently use parameters that do not reproduce anchors.
                params = None
        if params is None:
            params = build_bridge_parameters()
            try:
                Path(params_path).parent.mkdir(parents=True, exist_ok=True)
                save_hill_parameters(params, params_path)
            except OSError:
                pass  # read-only location: still usable with in-memory params
        self.params = params

        anchors = [(a["viral_load_v"], a["copies_per_ml"]) for a in self.params["anchors"]]
        self._anchors_vl = np.array([a[0] for a in anchors], dtype=float)
        self._anchors_cp = np.array([a[1] for a in anchors], dtype=float)

        # Backwards-compatible attributes (old Hill record exposed these)
        self.Vmax = float(self._anchors_cp[-1])
        self.K = 0.5
        self.n = 1.0

    def normalized_to_copies(self, viral_load: float) -> float:
        """Convert the ODE model's dimensionless viral load to an indicative
        copies/mL via the anchored piecewise log-linear bridge.

        Args:
            viral_load: V from the ODE model (>= 0; may exceed 1.0)

        Returns:
            Indicative plasma viral load in copies/mL (not patient-calibrated).
        """
        return piecewise_loglinear(viral_load, self._anchors_vl, self._anchors_cp)
    
    def copies_to_risk_category(self, copies_per_ml: float) -> str:
        """Return risk category based on clinical viral load.
        
        Args:
            copies_per_ml: Plasma viral load in copies/mL
            
        Returns:
            risk_category: One of 'undetectable', 'low_risk', 'screening', 
                          'treatment', 'severe'
        """
        if copies_per_ml < 100:  # Detection limit
            return 'undetectable'
        elif copies_per_ml < self.CLINICAL_THRESHOLDS['screening']:
            return 'low_risk'
        elif copies_per_ml < self.CLINICAL_THRESHOLDS['treatment']:
            return 'screening'
        elif copies_per_ml < self.CLINICAL_THRESHOLDS['severe']:
            return 'treatment'
        else:
            return 'severe'
    
    # Infection is introduced at day 21 (~3 weeks post-transplant), matching
    # the clinical schedule in which DNAemia is screened from week 2-4 onward.
    INFECTION_DAY = 21.0

    def simulate_clinical_trajectory(
        self,
        scenario: str,
        weeks: int = 52,
        simulator_config: Optional[Dict] = None,
    ) -> pd.DataFrame:
        """Run the BKPyV ODE model and return a weekly copies/mL trajectory.

        Uses the ODE simulator (days as the time unit) with continuous drug
        dosing, then converts V to indicative copies/mL via the anchored
        bridge. This replaces the legacy implementation that simulated in
        1-hour steps for 52 weeks with the heuristic discrete simulator.

        Args:
            scenario: One of 'baseline', 'infection', 'tacrolimus', 'sirolimus'
            weeks: Number of weeks to simulate (default: 52)
            simulator_config: Optional dict of ODE parameter overrides
                (e.g. {'p': 12.0}) forwarded to ``BKPyVODESimulator``.

        Returns:
            DataFrame with columns: week, viral_load_norm, copies_per_ml,
                                 risk_category, drug_scenario
        """
        from vcm.plugins.transplant.bk_polyomavirus.bk_polyomavirus import (
            BKPolyomavirusPlugin,
        )
        from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator
        from vcm.core.models import Environment, Perturbation, PerturbationType

        valid_scenarios = ['baseline', 'infection', 'tacrolimus', 'sirolimus']
        if scenario not in valid_scenarios:
            raise ValueError(f"Invalid scenario: {scenario}. Must be one of {valid_scenarios}")

        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()
        environment = Environment()
        simulator = BKPyVODESimulator(simulator_config or {})

        perturbations = []
        if scenario in ('infection', 'tacrolimus', 'sirolimus'):
            perturbations.append(Perturbation(
                id='bkpyv_infection',
                name='BKPyV infection',
                perturbation_type=PerturbationType.VIRAL_INFECTION,
                target_id='viral_entry',
                magnitude=1.0,
                timing=self.INFECTION_DAY,  # days
                duration=None,
            ))
        if scenario == 'tacrolimus':
            perturbations.append(Perturbation(
                id='tacrolimus_treatment',
                name='Tacrolimus treatment',
                perturbation_type=PerturbationType.DRUG_TREATMENT,
                target_id='FKBP1A',
                magnitude=1.0,
                timing=0.0,      # from transplant day
                duration=None,   # continuous dosing
            ))
        elif scenario == 'sirolimus':
            perturbations.append(Perturbation(
                id='sirolimus_treatment',
                name='Sirolimus treatment',
                perturbation_type=PerturbationType.DRUG_TREATMENT,
                target_id='MTOR',
                magnitude=1.0,
                timing=0.0,
                duration=None,
            ))

        total_days = max(int(weeks * 7), int(self.INFECTION_DAY) + 7)
        result = simulator.simulate(
            initial_state=initial_state,
            perturbations=perturbations,
            environment=environment,
            n_steps=total_days,   # 1 day per step
            timestep=1.0,
        )

        daily_loads = [
            max(0.0, step.cell_state.metadata.get("viral_load", 0.0))
            for step in result.steps
        ]
        daily_days = [step.timestamp for step in result.steps]

        # Weekly sampling (week w -> day 7w)
        weekly_viral_loads = []
        for week in range(weeks + 1):
            target_day = min(7 * week, daily_days[-1])
            idx = min(range(len(daily_days)), key=lambda i: abs(daily_days[i] - target_day))
            weekly_viral_loads.append(daily_loads[idx])

        results = []
        for week, vl_norm in enumerate(weekly_viral_loads):
            copies = self.normalized_to_copies(vl_norm)
            results.append({
                'week': week,
                'viral_load_norm': vl_norm,
                'copies_per_ml': copies,
                'risk_category': self.copies_to_risk_category(copies),
                'drug_scenario': scenario,
            })

        return pd.DataFrame(results)

    def simulate_viral_load_trajectory(self, virtual_loads, timepoints, **kwargs):
        """Convert an existing normalized trajectory for legacy UI callers.

        This compatibility method intentionally preserves the historical
        dataframe shape while using the mapper's canonical Hill conversion and
        risk taxonomy.
        """
        frame = pd.DataFrame({"week": list(timepoints), "viral_load_norm": list(virtual_loads)})
        frame["copies_per_ml"] = [
            self.normalized_to_copies(float(value)) for value in virtual_loads
        ]
        frame["plasma_viral_load"] = frame["copies_per_ml"]
        frame["risk_category"] = [
            ClinicalThresholds.get_risk_category(float(value))
            for value in frame["copies_per_ml"]
        ]
        return frame


def generate_clinical_summary(trajectories: Dict[str, pd.DataFrame]) -> Dict:
    """Generate clinical summary statistics for each scenario.
    
    Args:
        trajectories: Dict mapping scenario names to DataFrames
        
    Returns:
        Dict with summary statistics for each scenario
    """
    summary = {}
    
    for scenario, df in trajectories.items():
        # Time to first detectable viremia (weeks post-transplant)
        detectable = df[df['copies_per_ml'] >= 100]
        if len(detectable) > 0:
            time_to_detectable = int(detectable['week'].min())
        else:
            time_to_detectable = None
        
        # Peak viral load (copies/mL)
        peak_viral_load = float(df['copies_per_ml'].max())
        
        # Week of peak
        week_of_peak = int(df.loc[df['copies_per_ml'].idxmax(), 'week'])
        
        # Time above screening threshold (weeks)
        above_screening = df[df['copies_per_ml'] >= 1000]
        time_above_screening = int(len(above_screening))
        
        # Time above treatment threshold (weeks)
        above_treatment = df[df['copies_per_ml'] >= 10000]
        time_above_treatment = int(len(above_treatment))
        
        # Final viral load at week 52
        final_viral_load = float(df[df['week'] == 52]['copies_per_ml'].values[0] if len(df[df['week'] == 52]) > 0 else df['copies_per_ml'].iloc[-1])
        
        summary[scenario] = {
            'time_to_first_detectable_viremia_weeks': time_to_detectable,
            'peak_viral_load_copies_per_ml': peak_viral_load,
            'week_of_peak': week_of_peak,
            'time_above_screening_threshold_weeks': time_above_screening,
            'time_above_treatment_threshold_weeks': time_above_treatment,
            'final_viral_load_week_52_copies_per_ml': final_viral_load
        }
    
    return summary

class ClinicalThresholds:
    DETECTION_LIMIT = 100
    SCREENING_POSITIVE = 1_000
    TREATMENT_POSITIVE = 10_000
    SEVERE = 1e7
    # Historical threshold name retained for callers that use the UI plot.
    HIGH_RISK = TREATMENT_POSITIVE

    @staticmethod
    def get_risk_category(copies_per_ml: float) -> str:
        if copies_per_ml < ClinicalThresholds.DETECTION_LIMIT:
            return "undetectable"
        if copies_per_ml < ClinicalThresholds.SCREENING_POSITIVE:
            return "low_risk"
        if copies_per_ml < ClinicalThresholds.TREATMENT_POSITIVE:
            return "screening"
        if copies_per_ml < ClinicalThresholds.SEVERE:
            return "treatment"
        return "severe"


# Backward-compatible name used by the original Streamlit dashboard.
ClinicalViralLoadMapper = ViralLoadMapper
