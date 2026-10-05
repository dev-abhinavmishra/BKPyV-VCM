#!/usr/bin/env python3
"""Test immune dynamics integration and improved sensitivity."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from vcm.plugins.transplant.bk_polyomavirus.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_simulator import BKPyVSimulator
from vcm.core.models import CellState, Environment, Perturbation, PerturbationType
from vcm.clinical.viral_load_mapper import ViralLoadMapper


def test_immune_dynamics_sensitivity():
    """Test whether immune dynamics improve parameter sensitivity."""
    
    print("Testing Immune Dynamics Integration")
    print("=" * 70)
    
    # Initialize
    plugin = BKPolyomavirusPlugin()
    mapper = ViralLoadMapper()
    
    # Test parameters with immune dynamics
    test_params = {
        'tacrolimus_enhancement_factor': ('tacrolimus', 1.8),
        'mtor_inhibition_factor': ('sirolimus', 0.5),
        't_antigen_replication_threshold': ('infection', 0.5),
    }
    
    perturbation_levels = [-0.4, 0.0, 0.4]
    perturbation_labels = ['-40%', '0%', '+40%']
    
    print("\nTesting parameter sensitivity WITH immune dynamics:\n")
    
    for param_name, (scenario, default_value) in test_params.items():
        print(f"Parameter: {param_name} (scenario: {scenario})")
        print(f"  Default: {default_value}")
        
        baseline_peak = None
        
        for perturbation, label in zip(perturbation_levels, perturbation_labels):
            if perturbation == 0.0:
                perturbed_value = default_value
            else:
                perturbed_value = default_value * (1 + perturbation)
            
            if perturbed_value <= 0:
                print(f"  {label}: Skipped (invalid value)")
                continue
            
            # Create simulator with perturbed config
            config = {param_name: perturbed_value}
            simulator = BKPyVSimulator(config)
            
            # Set up scenario
            initial_state = plugin.create_initial_state()
            environment = Environment()
            
            scenario_perturbations = [
                Perturbation(
                    id='bkpyv_infection',
                    name='BKPyV infection',
                    perturbation_type=PerturbationType.VIRAL_INFECTION,
                    target_id='viral_entry',
                    magnitude=1.0,
                    timing=504,
                    duration=None
                )
            ]
            
            if scenario == 'tacrolimus':
                scenario_perturbations.append(
                    Perturbation(
                        id='tacrolimus_treatment',
                        name='Tacrolimus treatment',
                        perturbation_type=PerturbationType.DRUG_TREATMENT,
                        target_id='FKBP1A',
                        magnitude=1.0,
                        timing=0.0,
                        duration=None
                    )
                )
            elif scenario == 'sirolimus':
                scenario_perturbations.append(
                    Perturbation(
                        id='sirolimus_treatment',
                        name='Sirolimus treatment',
                        perturbation_type=PerturbationType.DRUG_TREATMENT,
                        target_id='MTOR',
                        magnitude=1.0,
                        timing=0.0,
                        duration=None
                    )
                )
            
            # Run simulation
            result = simulator.simulate(
                initial_state=initial_state,
                perturbations=scenario_perturbations,
                environment=environment,
                n_steps=1000,
                timestep=1.0
            )
            
            # Get peak viral load and immune metrics
            viral_loads = [step.cell_state.viral_load for step in result.steps]
            peak_viral_load = max(viral_loads) if viral_loads else 0.0
            peak_copies = mapper.normalized_to_copies(peak_viral_load)
            
            # Get immune cell levels at peak
            peak_step = result.steps[viral_loads.index(peak_viral_load)] if viral_loads else result.steps[-1]
            immune_metrics = {
                'cd8_t_cells': peak_step.cell_state.metadata.get('cd8_t_cells', 0),
                'cd4_t_cells': peak_step.cell_state.metadata.get('cd4_t_cells', 0),
                'b_cells': peak_step.cell_state.metadata.get('b_cells', 0),
                'immune_clearance': peak_step.cell_state.metadata.get('immune_clearance_rate', 0),
            }
            
            if perturbation == 0.0:
                baseline_peak = peak_copies
            
            pct_change = ((peak_copies - baseline_peak) / baseline_peak * 100) if baseline_peak else 0.0
            
            print(f"  {label}: {peak_copies:,.0f} copies/mL ({pct_change:+.1f}%)")
            print(f"    Immune: CD8={immune_metrics['cd8_t_cells']:.1f}, CD4={immune_metrics['cd4_t_cells']:.1f}, Clearance={immune_metrics['immune_clearance']:.3f}")
        
        print()
    
    print("=" * 70)
    print("Immune dynamics testing completed!")


if __name__ == "__main__":
    test_immune_dynamics_sensitivity()
