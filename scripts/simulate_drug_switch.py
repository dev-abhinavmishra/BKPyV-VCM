#!/usr/bin/env python3
"""
Drug switching simulation for BKPyV VCM project.

Simulates the clinically motivated scenario:
- Patient starts on tacrolimus -> develops viremia -> switches to sirolimus
at the screening threshold (1,000 copies/mL) or treatment threshold
(10,000 copies/mL) weeks.

A switch is modelled as a SINGLE continuous ODE trajectory: tacrolimus dosing
from day 0 to the switch day, sirolimus dosing from the switch day onward.
(The previous implementation spliced two independent simulations, which
restarted the infection clock and is replaced here.)
"""

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcm.clinical.viral_load_mapper import ViralLoadMapper
from vcm.core.models import Perturbation, PerturbationType
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator


def simulate_drug_switch(switch_week=None, weeks: int = 52):
    """Run one continuous ODE trajectory with a tac->sir switch at `switch_week`.

    Args:
        switch_week: Week at which tacrolimus stops and sirolimus starts
            (None = tacrolimus throughout).
        weeks: Horizon in weeks.

    Returns:
        DataFrame with weekly columns: week, viral_load_norm, copies_per_ml.
    """
    mapper = ViralLoadMapper()
    total_days = weeks * 7
    infection_day = ViralLoadMapper.INFECTION_DAY

    perturbations = [
        Perturbation(
            id="bkpyv_infection",
            name="BKPyV infection",
            perturbation_type=PerturbationType.VIRAL_INFECTION,
            target_id="viral_entry",
            magnitude=1.0,
            timing=infection_day,
        ),
        Perturbation(
            id="tacrolimus_baseline",
            name="Tacrolimus (baseline immunosuppression)",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="FKBP1A",
            magnitude=1.0,
            timing=0.0,
            duration=(None if switch_week is None else switch_week * 7.0),
        ),
    ]
    if switch_week is not None:
        perturbations.append(Perturbation(
            id="sirolimus_switch",
            name="Sirolimus after switch",
            perturbation_type=PerturbationType.DRUG_TREATMENT,
            target_id="MTOR",
            magnitude=1.0,
            timing=float(switch_week * 7),
            duration=None,
        ))

    plugin = BKPolyomavirusPlugin()
    sim = BKPyVODESimulator({})
    result = sim.simulate(
        initial_state=plugin.create_initial_state(),
        perturbations=perturbations,
        n_steps=total_days,
        timestep=1.0,
    )

    daily = pd.DataFrame(
        {
            "day": [s.timestamp for s in result.steps],
            "viral_load_norm": [
                max(0.0, s.cell_state.metadata.get("viral_load", 0.0)) for s in result.steps
            ],
        }
    )
    weekly = daily.groupby((daily["day"] // 7).astype(int)).last().reset_index(drop=True)
    weekly.insert(0, "week", weekly.index)
    weekly = weekly[weekly["week"] <= weeks].reset_index(drop=True)
    weekly["copies_per_ml"] = weekly["viral_load_norm"].apply(mapper.normalized_to_copies)
    weekly["risk_category"] = weekly["copies_per_ml"].apply(mapper.copies_to_risk_category)
    weekly["drug_scenario"] = "no_switch" if switch_week is None else f"switch_week_{switch_week}"
    weekly["intervention"] = (
        "No switch" if switch_week is None else f"Switch at week {switch_week}"
    )
    return weekly


def main():
    """Main function to run drug switching simulation."""
    print("=" * 80)
    print("DRUG SWITCHING SIMULATION")
    print("=" * 80)
    print()
    
    # Run three scenarios
    print("Running simulation scenarios...")
    
    # Scenario A: No intervention (tacrolimus only)
    print("  Scenario A: No intervention (tacrolimus only)")
    df_no_switch = simulate_drug_switch(switch_week=None)
    print(f"    Peak viral load: {df_no_switch['copies_per_ml'].max():,.0f} copies/mL")
    print(f"    Viral load at week 52: {df_no_switch['copies_per_ml'].iloc[-1]:,.0f} copies/mL")
    
    # Scenario B: Switch at week 8 (treatment threshold)
    print("  Scenario B: Switch to sirolimus at week 8 (treatment threshold)")
    df_switch_week8 = simulate_drug_switch(switch_week=8)
    print(f"    Peak viral load: {df_switch_week8['copies_per_ml'].max():,.0f} copies/mL")
    print(f"    Viral load at week 52: {df_switch_week8['copies_per_ml'].iloc[-1]:,.0f} copies/mL")
    
    # Scenario C: Switch at week 4 (early intervention)
    print("  Scenario C: Switch to sirolimus at week 4 (screening threshold)")
    df_switch_week4 = simulate_drug_switch(switch_week=4)
    print(f"    Peak viral load: {df_switch_week4['copies_per_ml'].max():,.0f} copies/mL")
    print(f"    Viral load at week 52: {df_switch_week4['copies_per_ml'].iloc[-1]:,.0f} copies/mL")
    
    print()
    
    # Generate comparison figure
    print("Generating drug switch comparison figure...")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot trajectories
    ax.plot(df_no_switch['week'], df_no_switch['copies_per_ml'], 
            'r-', linewidth=2, label='No switch (tacrolimus only)', alpha=0.8)
    ax.plot(df_switch_week8['week'], df_switch_week8['copies_per_ml'], 
            'g-', linewidth=2, label='Switch at week 8 (treatment threshold)', alpha=0.8)
    ax.plot(df_switch_week4['week'], df_switch_week4['copies_per_ml'], 
            'b-', linewidth=2, label='Switch at week 4 (early intervention)', alpha=0.8)
    
    # Add intervention lines
    ax.axvline(x=8, color='g', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axvline(x=4, color='b', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Add intervention annotations
    ax.text(8.2, ax.get_ylim()[1] * 0.95, 'Week 8:\nTreatment\nthreshold', 
            fontsize=9, color='green', va='top')
    ax.text(4.2, ax.get_ylim()[1] * 0.95, 'Week 4:\nScreening\nthreshold', 
            fontsize=9, color='blue', va='top')
    
    # Add clinical thresholds
    ax.axhline(y=10000, color='orange', linestyle=':', linewidth=1, alpha=0.5)
    ax.axhline(y=1000, color='yellow', linestyle=':', linewidth=1, alpha=0.5)
    
    # Log scale
    ax.set_yscale('log')
    ax.set_ylim(100, 1e7)
    
    # Labels and title
    ax.set_xlabel('Weeks Post-Transplant', fontsize=12)
    ax.set_ylabel('Viral Load (copies/mL)', fontsize=12)
    ax.set_title('Drug Switching Simulation: BKPyV Viral Load Trajectories', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Add clinical threshold labels
    ax.text(1, 10000, 'Treatment threshold (10,000)', fontsize=8, color='orange', va='bottom')
    ax.text(1, 1000, 'Screening threshold (1,000)', fontsize=8, color='orange', va='bottom')
    
    plt.tight_layout()
    
    # Save figure
    output_dir = Path("outputs/clinical")
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / "drug_switch_comparison.png", dpi=300, bbox_inches='tight')
    print(f"Saved to {output_dir / 'drug_switch_comparison.png'}")
    plt.close()
    
    # Summary statistics
    print()
    print("=" * 80)
    print("DRUG SWITCHING SIMULATION SUMMARY")
    print("=" * 80)
    print()
    print("Scenario Comparison:")
    print(f"  {'Scenario':<30} {'Peak (copies/mL)':<20} {'Week 52 (copies/mL)':<20}")
    print(f"  {'-'*30} {'-'*20} {'-'*20}")
    print(f"  {'No switch (tacrolimus)':<30} {df_no_switch['copies_per_ml'].max():>15,.0f} {df_no_switch['copies_per_ml'].iloc[-1]:>15,.0f}")
    print(f"  {'Switch at week 8':<30} {df_switch_week8['copies_per_ml'].max():>15,.0f} {df_switch_week8['copies_per_ml'].iloc[-1]:>15,.0f}")
    print(f"  {'Switch at week 4':<30} {df_switch_week4['copies_per_ml'].max():>15,.0f} {df_switch_week4['copies_per_ml'].iloc[-1]:>15,.0f}")
    print()
    
    # Calculate reduction vs no switch (guard against a zero denominator)
    denom = max(float(df_no_switch['copies_per_ml'].iloc[-1]), 1e-9)
    reduction_week8 = (df_no_switch['copies_per_ml'].iloc[-1] - df_switch_week8['copies_per_ml'].iloc[-1]) / denom * 100
    reduction_week4 = (df_no_switch['copies_per_ml'].iloc[-1] - df_switch_week4['copies_per_ml'].iloc[-1]) / denom * 100
    
    print("Viral load reduction at week 52 vs no switch:")
    print(f"  Switch at week 8: {reduction_week8:.1f}% reduction")
    print(f"  Switch at week 4: {reduction_week4:.1f}% reduction")
    print()
    
    # Save results to JSON
    results = {
        'no_switch': {
            'peak_viral_load': float(df_no_switch['copies_per_ml'].max()),
            'week_52_viral_load': float(df_no_switch['copies_per_ml'].iloc[-1])
        },
        'switch_week_8': {
            'peak_viral_load': float(df_switch_week8['copies_per_ml'].max()),
            'week_52_viral_load': float(df_switch_week8['copies_per_ml'].iloc[-1]),
            'reduction_vs_no_switch': reduction_week8
        },
        'switch_week_4': {
            'peak_viral_load': float(df_switch_week4['copies_per_ml'].max()),
            'week_52_viral_load': float(df_switch_week4['copies_per_ml'].iloc[-1]),
            'reduction_vs_no_switch': reduction_week4
        },
        'clinical_thresholds': {
            'screening_threshold': 1000,
            'treatment_threshold': 10000
        },
        'interpretation': (
            "Model-internal prediction only: switching earlier reduces simulated "
            "viral burden under this model's assumptions (assumption-labelled "
            "bridge, immune-control formulation). Not a clinical recommendation."
        )
    }
    
    import json
    results_path = output_dir / "drug_switch_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved results to {results_path}")
    print()
    
    print("=" * 80)
    print("Drug switching simulation complete!")
    print("=" * 80)

if __name__ == "__main__":
    main()
