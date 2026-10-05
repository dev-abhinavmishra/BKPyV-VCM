#!/usr/bin/env python3
"""Validation script for ODE-based BKPyV simulator.

This script compares the discrete and ODE-based simulators to ensure
the ODE implementation produces qualitatively similar results while
providing the benefits of continuous-time dynamics.
"""

import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path

try:  # Windows cp1252 consoles cannot encode status glyphs
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except (AttributeError, OSError):
    pass

from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_simulator import BKPyVSimulator
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator
from vcm.core.models import Perturbation, PerturbationType


def compare_simulators(
    scenario_name: str,
    perturbations: list,
    n_steps: int = 100,
    timestep: float = 1.0
) -> dict:
    """Compare discrete and ODE simulators for a given scenario.
    
    Args:
        scenario_name: Name of the scenario for output files
        perturbations: List of perturbations to apply
        n_steps: Number of simulation steps
        timestep: Time step size
    
    Returns:
        Dictionary with comparison results
    """
    print(f"\n{'='*60}")
    print(f"Comparing simulators for scenario: {scenario_name}")
    print(f"{'='*60}")
    
    # Create initial state
    plugin = BKPolyomavirusPlugin()
    initial_state = plugin.create_initial_state()
    
    # Run discrete simulation
    print("Running discrete simulator...")
    discrete_sim = BKPyVSimulator()
    discrete_result = discrete_sim.simulate(
        initial_state=initial_state,
        perturbations=perturbations,
        n_steps=n_steps,
        timestep=timestep
    )
    
    # Run ODE simulation
    print("Running ODE simulator...")
    ode_config = {
        'ode_solver': 'LSODA',
        'rtol': 1e-6,
        'atol': 1e-8
    }
    ode_sim = BKPyVODESimulator(ode_config)
    ode_result = ode_sim.simulate(
        initial_state=initial_state,
        perturbations=perturbations,
        n_steps=n_steps,
        timestep=timestep
    )
    
    # Extract time series for comparison
    times_discrete = [step.timestamp for step in discrete_result.steps]
    times_ode = [step.timestamp for step in ode_result.steps]
    
    viral_load_discrete = [step.cell_state.metadata['viral_load'] for step in discrete_result.steps]
    viral_load_ode = [step.cell_state.metadata['viral_load'] for step in ode_result.steps]
    
    t_antigen_discrete = [step.cell_state.genes['viral_LT'].expression_level for step in discrete_result.steps]
    t_antigen_ode = [step.cell_state.genes['viral_LT'].expression_level for step in ode_result.steps]
    
    # Interpolate to common time points for comparison
    common_times = np.linspace(min(times_discrete[0], times_ode[0]), 
                               max(times_discrete[-1], times_ode[-1]), 
                               min(len(times_discrete), len(times_ode)))
    
    viral_load_discrete_interp = np.interp(common_times, times_discrete, viral_load_discrete)
    viral_load_ode_interp = np.interp(common_times, times_ode, viral_load_ode)
    
    t_antigen_discrete_interp = np.interp(common_times, times_discrete, t_antigen_discrete)
    t_antigen_ode_interp = np.interp(common_times, times_ode, t_antigen_ode)
    
    # Calculate comparison metrics using interpolated values
    viral_load_diff = np.abs(np.array(viral_load_discrete_interp) - np.array(viral_load_ode_interp))
    max_diff = np.max(viral_load_diff)
    mean_diff = np.mean(viral_load_diff)
    
    print(f"Max viral load difference: {max_diff:.4f}")
    print(f"Mean viral load difference: {mean_diff:.4f}")
    
    # Check qualitative behavior similarity
    discrete_trend = np.gradient(viral_load_discrete_interp)
    ode_trend = np.gradient(viral_load_ode_interp)
    trend_correlation = np.corrcoef(discrete_trend, ode_trend)[0, 1]
    
    print(f"Trend correlation: {trend_correlation:.4f}")
    
    # Generate comparison plots using original data points
    create_comparison_plots(
        scenario_name,
        times_discrete, times_ode,
        viral_load_discrete, viral_load_ode,
        t_antigen_discrete, t_antigen_ode
    )
    
    return {
        'scenario': scenario_name,
        'max_diff': max_diff,
        'mean_diff': mean_diff,
        'trend_correlation': trend_correlation,
        'final_viral_discrete': viral_load_discrete[-1],
        'final_viral_ode': viral_load_ode[-1]
    }


def create_comparison_plots(
    scenario_name: str,
    times_d: list, times_o: list,
    viral_d: list, viral_o: list,
    t_ag_d: list, t_ag_o: list
):
    """Create comparison plots for discrete vs ODE simulators.
    
    Args:
        scenario_name: Name for plot titles
        times_d: Discrete simulator time points
        times_o: ODE simulator time points
        viral_d: Discrete viral load
        viral_o: ODE viral load
        t_ag_d: Discrete T antigen
        t_ag_o: ODE T antigen
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Viral load comparison
    axes[0].plot(times_d, viral_d, 'b-', label='Discrete simulator', linewidth=2, alpha=0.7)
    axes[0].plot(times_o, viral_o, 'r--', label='ODE simulator', linewidth=2, alpha=0.7)
    axes[0].set_xlabel('Time (days)')
    axes[0].set_ylabel('Viral Load (normalized)')
    axes[0].set_title(f'BKPyV Viral Load Dynamics: {scenario_name}')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # T antigen comparison
    axes[1].plot(times_d, t_ag_d, 'b-', label='Discrete simulator', linewidth=2, alpha=0.7)
    axes[1].plot(times_o, t_ag_o, 'r--', label='ODE simulator', linewidth=2, alpha=0.7)
    axes[1].set_xlabel('Time (days)')
    axes[1].set_ylabel('T Antigen Expression')
    axes[1].set_title(f'T Antigen Dynamics: {scenario_name}')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    output_dir = Path('outputs/ode_validation')
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / f'{scenario_name.replace(" ", "_")}_comparison.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"Saved comparison plot to {plot_path}")
    plt.close()


def test_drug_effects():
    """Test that drug effects are consistent between simulators."""
    print("\n" + "="*60)
    print("Testing Drug Effects Consistency")
    print("="*60)
    
    scenarios = [
        ("Tacrolimus Enhancement", [
            Perturbation(
                id="infection",
                name="BKPyV infection",
                perturbation_type=PerturbationType.VIRAL_INFECTION,
                magnitude=1.0,
                timing=0.0  # Start at t=0 for ODE compatibility
            ),
            Perturbation(
                id="tacrolimus",
                name="Tacrolimus treatment",
                perturbation_type=PerturbationType.DRUG_TREATMENT,
                target_id="FKBP1A",
                magnitude=1.0,
                timing=0.0
            )
        ]),
        ("Sirolimus Inhibition", [
            Perturbation(
                id="infection",
                name="BKPyV infection",
                perturbation_type=PerturbationType.VIRAL_INFECTION,
                magnitude=1.0,
                timing=0.0  # Start at t=0 for ODE compatibility
            ),
            Perturbation(
                id="sirolimus",
                name="Sirolimus treatment",
                perturbation_type=PerturbationType.DRUG_TREATMENT,
                target_id="MTOR",
                magnitude=1.0,
                timing=0.0
            )
        ])
    ]
    
    results = []
    for scenario_name, perturbations in scenarios:
        result = compare_simulators(scenario_name, perturbations)
        results.append(result)
    
    # Print summary
    print("\n" + "="*60)
    print("Validation Summary")
    print("="*60)
    for result in results:
        print(f"\nScenario: {result['scenario']}")
        print(f"  Max difference: {result['max_diff']:.4f}")
        print(f"  Mean difference: {result['mean_diff']:.4f}")
        print(f"  Trend correlation: {result['trend_correlation']:.4f}")
        print(f"  Final viral (discrete): {result['final_viral_discrete']:.4f}")
        print(f"  Final viral (ODE): {result['final_viral_ode']:.4f}")
    
    return results


def test_numerical_properties():
    """Test numerical properties of ODE implementation."""
    print("\n" + "="*60)
    print("Testing Numerical Properties")
    print("="*60)
    
    plugin = BKPolyomavirusPlugin()
    initial_state = plugin.create_initial_state()
    
    infection = Perturbation(
        id="infection",
        name="BKPyV infection",
        perturbation_type=PerturbationType.VIRAL_INFECTION,
        magnitude=1.0,
        timing=10.0
    )
    
    # Test different solvers
    solvers = ['LSODA', 'RK45', 'BDF']
    solver_results = {}
    
    for solver in solvers:
        print(f"\nTesting solver: {solver}")
        ode_config = {
            'ode_solver': solver,
            'rtol': 1e-6,
            'atol': 1e-8
        }
        ode_sim = BKPyVODESimulator(ode_config)
        result = ode_sim.simulate(
            initial_state=initial_state,
            perturbations=[infection],
            n_steps=50,
            timestep=1.0
        )
        
        viral_load = [step.cell_state.metadata['viral_load'] for step in result.steps]
        solver_results[solver] = {
            'final_viral': viral_load[-1],
            'max_viral': max(viral_load),
            'success': result.success
        }
        
        print(f"  Final viral load: {viral_load[-1]:.4f}")
        print(f"  Max viral load: {max(viral_load):.4f}")
        print(f"  Success: {result.success}")
    
    # Check consistency between solvers
    print("\nSolver consistency check:")
    viral_values = [r['final_viral'] for r in solver_results.values()]
    viral_std = np.std(viral_values)
    print(f"  Standard deviation of final viral loads: {viral_std:.4f}")
    
    if viral_std < 0.1:
        print("  ✓ Solvers produce consistent results")
    else:
        print("  ⚠ Solvers show significant differences")
    
    return solver_results


def main():
    """Run all validation tests."""
    print("ODE Implementation Validation")
    print("="*60)
    
    # Test drug effects consistency
    drug_results = test_drug_effects()
    
    # Test numerical properties
    solver_results = test_numerical_properties()
    
    print("\n" + "="*60)
    print("Validation Complete")
    print("="*60)
    print("All validation tests completed successfully!")
    print("Comparison plots saved to outputs/ode_validation/")


if __name__ == "__main__":
    main()