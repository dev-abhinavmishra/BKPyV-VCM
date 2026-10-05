#!/usr/bin/env python3
"""
Sensitivity analysis for BKPyV simulator parameters.

Performs one-at-a-time (OAT) sensitivity analysis on all model parameters
to assess how parameter uncertainty affects model outputs (peak viral load).
"""

import matplotlib

matplotlib.use("Agg")  # headless backend: tkinter is unavailable/crashy in some environments

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcm.clinical.viral_load_mapper import ViralLoadMapper

def main():
    """Main function to perform sensitivity analysis."""
    print("=" * 80)
    print("BKPyV SIMULATOR SENSITIVITY ANALYSIS")
    print("=" * 80)
    print()
    
    # The mapper runs the ODE simulator; simulator_config overrides are passed
    # through on every call so perturbations actually reach the model
    # (the previous version constructed a perturbed simulator that was never
    # used, producing a flat tornado chart).
    mapper = ViralLoadMapper()
    
    # Parameters wired into the ODE right-hand side. Unwired placeholders
    # (ddr_enhancement, mitochondrial_importance, protein_degradation_inhibition)
    # are excluded: they are registered but not yet dynamic, so their
    # sensitivity is identically zero by construction.
    default_params = {
        'tacrolimus_enhancement_factor': 1.5,
        'mtor_inhibition_factor': 0.5,
        't_antigen_replication_threshold': 0.5,
        'cell_cycle_s_phase_bonus': 2.0,
        'dna_replication_coupling': 0.8,
        'innate_immune_suppression_factor': 0.5,
        'translation_enhancement_factor': 2.0,
    }
    unwired = ['dna_damage_response_enhancement',
               'mitochondrial_function_importance',
               'protein_degradation_inhibition']
    
    print("Parameters for sensitivity analysis:")
    for param, value in default_params.items():
        print(f"  {param}: {value}")
    print("Skipped (registered but not wired into the ODE dynamics):")
    for param in unwired:
        print(f"  {param}: excluded from sweep")
    print()
    
    # Define perturbation levels: ±20%, ±40%
    perturbations = [-0.4, -0.2, 0.0, 0.2, 0.4]
    perturbation_labels = ['-40%', '-20%', '0%', '+20%', '+40%']
    
    # Baseline runs (no perturbations)
    print("Running baseline simulations...")
    baseline_results = {}
    scenarios = ['infection', 'tacrolimus', 'sirolimus']
    
    for scenario in scenarios:
        df = mapper.simulate_clinical_trajectory(scenario, weeks=52)
        baseline_results[scenario] = df['copies_per_ml'].max()
        print(f"  {scenario}: {baseline_results[scenario]:,.0f} copies/mL")
    
    print()
    print("Running sensitivity analysis (OAT)...")
    
    # Store results
    sensitivity_results = []
    
    for param_name, default_value in default_params.items():
        param_results = []
        
        for perturbation, label in zip(perturbations, perturbation_labels):
            # Calculate perturbed value
            if perturbation == 0.0:
                perturbed_value = default_value
            else:
                perturbed_value = default_value * (1 + perturbation)
            
            # Skip invalid values (negative, > reasonable bounds)
            if perturbed_value <= 0:
                param_results.append({
                    'parameter': param_name,
                    'perturbation': label,
                    'perturbation_value': perturbed_value,
                    'infection_peak': np.nan,
                    'tacrolimus_peak': np.nan,
                    'sirolimus_peak': np.nan
                })
                continue
            
            # Perturbed parameter value is passed THROUGH to the ODE simulator
            config = {param_name: perturbed_value}

            try:
                # Run simulations for each scenario
                scenario_peaks = {}
                for scenario in scenarios:
                    # Map scenario to drug settings
                    if scenario == 'tacrolimus':
                        mapper_scenario = 'tacrolimus'
                    elif scenario == 'sirolimus':
                        mapper_scenario = 'sirolimus'
                    else:
                        mapper_scenario = 'infection'

                    df = mapper.simulate_clinical_trajectory(
                        mapper_scenario, weeks=52, simulator_config=config
                    )
                    scenario_peaks[scenario] = df['copies_per_ml'].max()
                
                param_results.append({
                    'parameter': param_name,
                    'perturbation': label,
                    'perturbation_value': perturbed_value,
                    'infection_peak': scenario_peaks['infection'],
                    'tacrolimus_peak': scenario_peaks['tacrolimus'],
                    'sirolimus_peak': scenario_peaks['sirolimus']
                })
                
            except Exception as e:
                print(f"    Error with {param_name} at {label}: {e}")
                param_results.append({
                    'parameter': param_name,
                    'perturbation': label,
                    'perturbation_value': perturbed_value,
                    'infection_peak': np.nan,
                    'tacrolimus_peak': np.nan,
                    'sirolimus_peak': np.nan
                })
        
        sensitivity_results.extend(param_results)
        print(f"  Completed: {param_name}")
    
    # Convert to DataFrame
    df = pd.DataFrame(sensitivity_results)
    
    # Calculate percentage change from baseline
    for scenario in scenarios:
        baseline = baseline_results[scenario]
        df[f'{scenario}_pct_change'] = (df[f'{scenario}_peak'] - baseline) / baseline * 100
    
    # Save results
    output_dir = Path("outputs/sensitivity")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(output_dir / "sensitivity_results.csv", index=False)
    print(f"\nSaved results to {output_dir / 'sensitivity_results.csv'}")
    
    # Generate tornado plot
    print("\nGenerating tornado plot...")
    
    # Calculate sensitivity (max absolute % change across all perturbations)
    sensitivity_summary = []
    for param_name in default_params.keys():
        param_data = df[df['parameter'] == param_name]
        
        for scenario in scenarios:
            # Get max absolute % change for this parameter in this scenario
            max_change = param_data[f'{scenario}_pct_change'].abs().max()
            
            sensitivity_summary.append({
                'parameter': param_name,
                'scenario': scenario,
                'max_pct_change': max_change,
                'direction': 'increase' if param_data.loc[param_data[f'{scenario}_pct_change'].abs().idxmax(), f'{scenario}_pct_change'] > 0 else 'decrease'
            })
    
    sens_df = pd.DataFrame(sensitivity_summary)
    
    # Sort by sensitivity (absolute change)
    sens_df = sens_df.sort_values('max_pct_change', ascending=False)
    
    # Create tornado plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create horizontal bar chart (tacrolimus vs sirolimus per parameter)
    parameters = sens_df['parameter'].unique()
    y_pos = np.arange(len(parameters))
    
    # Get values for each scenario
    tac_changes = sens_df[sens_df['scenario'] == 'tacrolimus'].set_index('parameter').loc[parameters, 'max_pct_change']
    sir_changes = sens_df[sens_df['scenario'] == 'sirolimus'].set_index('parameter').loc[parameters, 'max_pct_change']
    
    # Plot bars
    bar_height = 0.35
    ax.barh(y_pos + bar_height/2, tac_changes, bar_height, label='Tacrolimus', color='red', alpha=0.7)
    ax.barh(y_pos - bar_height/2, sir_changes, bar_height, label='Sirolimus', color='blue', alpha=0.7)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(parameters, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('Maximum % Change in Peak Viral Load', fontsize=12)
    ax.set_title('Parameter Sensitivity Analysis: BKPyV Simulator', fontsize=14, fontweight='bold')
    ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
    ax.legend(loc='best', fontsize=11)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "tornado_plot.png", dpi=300, bbox_inches='tight')
    print(f"Saved tornado plot to {output_dir / 'tornado_plot.png'}")
    plt.close()
    
    # Print summary
    print("\n" + "=" * 80)
    print("SENSITIVITY ANALYSIS SUMMARY")
    print("=" * 80)
    print("\nTop 5 most sensitive parameters (tacrolimus scenario):")
    top_tac = sens_df[sens_df['scenario'] == 'tacrolimus'].head(5)
    for _, row in top_tac.iterrows():
        print(f"  {row['parameter']}: {row['max_pct_change']:.2f}% ({row['direction']})")
    
    print("\nTop 5 most sensitive parameters (sirolimus scenario):")
    top_sir = sens_df[sens_df['scenario'] == 'sirolimus'].head(5)
    for _, row in top_sir.iterrows():
        print(f"  {row['parameter']}: {row['max_pct_change']:.2f}% ({row['direction']})")
    
    print("\n" + "=" * 80)
    print("Sensitivity analysis complete!")
    print("=" * 80)

if __name__ == "__main__":
    main()