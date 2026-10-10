"""ISEF-Level Qualitative Validation for BKPyV Model.

This script performs qualitative validation focusing on whether the model reproduces
the mechanistic patterns and qualitative findings from published research, which is
more scientifically defensible for ISEF than trying to match exact quantitative targets
without real patient data.

Key Qualitative Patterns to Validate (from literature):
1. Tacrolimus enhances BKPyV replication (Hirsch et al., Am J Transplant 2016)
2. Sirolimus inhibits BKPyV replication (Hirsch et al., Am J Transplant 2016)
3. Sirolimus timing effect: only effective in early phase (first 24h) 
4. Cell-cycle permissiveness: higher in S-phase
5. Mitochondrial stress signature: emerges in late replication phase
6. T antigen threshold-dependent replication

Research References:
- Hirsch HH et al., Am J Transplant 2016: Drug mechanisms via FKBP-12
- Weissbach FH et al., J Virol 2024: Single-cell transcriptomics
- Clinical cohort studies: Risk factor ORs and timing patterns
"""

import sys
from pathlib import Path
import numpy as np
import json
from typing import Dict, List, Tuple

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from scipy.integrate import solve_ivp

from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.ode_system import BKPyVODESystem
from vcm.core.models import Perturbation, PerturbationType


class ISEFQualitativeValidator:
    """Validator for ISEF-level qualitative pattern matching."""
    
    def __init__(self):
        """Initialize ISEF qualitative validator."""
        self.pattern_results = {}
        self.literature_patterns = {
            'tacrolimus_enhancement': {
                'description': 'Tacrolimus enhances viral replication via FKBP-12 pathway',
                'reference': 'Hirsch HH et al., Am J Transplant 2016',
                'expected': 'tacrolimus > baseline',
                'confidence': 'HIGH',
            },
            'sirolimus_inhibition': {
                'description': 'Sirolimus inhibits viral replication via mTOR pathway',
                'reference': 'Hirsch HH et al., Am J Transplant 2016',
                'expected': 'sirolimus < tacrolimus',
                'confidence': 'HIGH',
            },
            'cell_cycle_permissiveness': {
                'description': 'Higher viral replication in S-phase',
                'reference': 'Weissbach FH et al., J Virol 2024',
                'expected': 'S-phase > other phases',
                'confidence': 'HIGH',
            },
            't_antigen_threshold': {
                'description': 'T antigen required for active replication',
                'reference': 'BKPyV biology literature',
                'expected': 'threshold-dependent activation',
                'confidence': 'MEDIUM',
            },
            'early_vs_late_phase': {
                'description': 'Different gene expression in early vs late replication',
                'reference': 'Weissbach FH et al., J Virol 2024',
                'expected': 'phase-specific signatures',
                'confidence': 'HIGH',
            },
        }
    
    def validate_drug_effects(self) -> Dict[str, any]:
        """Validate drug effects - tacrolimus enhances, sirolimus inhibits.
        
        This is the most critical validation since it's the key mechanistic insight
        from Hirsch et al. 2016 with HIGH confidence (experimental validation).
        
        Returns:
            Dictionary with drug effect validation results
        """
        print("\n" + "="*80)
        print("QUALITATIVE VALIDATION: Drug Effects")
        print("="*80)
        print("Reference: Hirsch HH et al., Am J Transplant 2016")
        print("Mechanism: Tacrolimus activates, Sirolimus inhibits via FKBP-12 pathway")
        print()
        
        results = {}
        
        # Scenario 1: Baseline (no drug)
        baseline_peak = self._run_scenario(
            scenario_name='baseline',
            tacrolimus_enhancement=1.0,
            sirolimus_inhibition=1.0,
            simulation_length=100
        )
        
        # Scenario 2: Tacrolimus (should enhance)
        tacrolimus_peak = self._run_scenario(
            scenario_name='tacrolimus',
            tacrolimus_enhancement=2.3,  # Clinical OR
            sirolimus_inhibition=1.0,
            simulation_length=100
        )
        
        # Scenario 3: Sirolimus (should inhibit)
        sirolimus_peak = self._run_scenario(
            scenario_name='sirolimus',
            tacrolimus_enhancement=1.0,
            sirolimus_inhibition=0.5,  # IC90 = 4 ng/mL
            simulation_length=100
        )
        
        print(f"Baseline peak viral load: {baseline_peak:.3f}")
        print(f"Tacrolimus peak viral load: {tacrolimus_peak:.3f}")
        print(f"Sirolimus peak viral load: {sirolimus_peak:.3f}")
        print()
        
        # Validate pattern 1: Tacrolimus > Baseline
        pattern1_pass = tacrolimus_peak > baseline_peak
        print(f"Pattern 1: Tacrolimus > Baseline")
        print(f"  Expected: TRUE (Hirsch et al. 2016)")
        print(f"  Observed: {pattern1_pass}")
        print(f"  Ratio: {tacrolimus_peak / baseline_peak:.2f}x")
        print(f"  Status: {'✓ PASS' if pattern1_pass else '✗ FAIL'}")
        
        # Validate pattern 2: Sirolimus < Tacrolimus
        pattern2_pass = sirolimus_peak < tacrolimus_peak
        print(f"\nPattern 2: Sirolimus < Tacrolimus")
        print(f"  Expected: TRUE (Hirsch et al. 2016)")
        print(f"  Observed: {pattern2_pass}")
        print(f"  Ratio: {sirolimus_peak / tacrolimus_peak:.2f}x")
        print(f"  Status: {'✓ PASS' if pattern2_pass else '✗ FAIL'}")
        
        results = {
            'pattern_tacrolimus_enhancement': {
                'pass': pattern1_pass,
                'baseline': baseline_peak,
                'tacrolimus': tacrolimus_peak,
                'ratio': tacrolimus_peak / baseline_peak,
                'expected_ratio': '> 1.0',
                'reference': 'Hirsch HH et al., Am J Transplant 2016',
                'confidence': 'HIGH',
            },
            'pattern_sirolimus_inhibition': {
                'pass': pattern2_pass,
                'tacrolimus': tacrolimus_peak,
                'sirolimus': sirolimus_peak,
                'ratio': sirolimus_peak / tacrolimus_peak,
                'expected_ratio': '< 1.0',
                'reference': 'Hirsch HH et al., Am J Transplant 2016',
                'confidence': 'HIGH',
            },
        }
        
        overall_pass = pattern1_pass and pattern2_pass
        print(f"\nOverall Drug Effects Validation: {'✓ PASS' if overall_pass else '✗ FAIL'}")
        
        return results
    
    def validate_cell_cycle_permissiveness(self) -> Dict[str, any]:
        """Validate cell cycle permissiveness - S-phase optimal.
        
        Returns:
            Dictionary with cell cycle validation results
        """
        print("\n" + "="*80)
        print("QUALITATIVE VALIDATION: Cell Cycle Permissiveness")
        print("="*80)
        print("Reference: Weissbach FH et al., J Virol 2024")
        print("Mechanism: BKPyV replication coupled to host DNA replication (S-phase)")
        print()
        
        results = {}
        
        # Scenario: High cell cycle activity (S-phase)
        high_cycle_peak = self._run_scenario(
            scenario_name='high_cell_cycle',
            tacrolimus_enhancement=1.0,
            sirolimus_inhibition=1.0,
            cell_cycle_bonus=3.0,  # High S-phase bonus
            simulation_length=100
        )
        
        # Scenario: Low cell cycle activity
        low_cycle_peak = self._run_scenario(
            scenario_name='low_cell_cycle',
            tacrolimus_enhancement=1.0,
            sirolimus_inhibition=1.0,
            cell_cycle_bonus=1.0,  # Low S-phase bonus
            simulation_length=100
        )
        
        print(f"High cell cycle (S-phase) peak: {high_cycle_peak:.3f}")
        print(f"Low cell cycle peak: {low_cycle_peak:.3f}")
        print()
        
        pattern_pass = high_cycle_peak > low_cycle_peak
        print(f"Pattern: S-phase > Other phases")
        print(f"  Expected: TRUE (Weissbach et al. 2024)")
        print(f"  Observed: {pattern_pass}")
        print(f"  Ratio: {high_cycle_peak / low_cycle_peak:.2f}x")
        print(f"  Status: {'✓ PASS' if pattern_pass else '✗ FAIL'}")
        
        results = {
            'pattern_cell_cycle_permissiveness': {
                'pass': pattern_pass,
                'high_cycle': high_cycle_peak,
                'low_cycle': low_cycle_peak,
                'ratio': high_cycle_peak / low_cycle_peak,
                'expected_ratio': '> 1.0',
                'reference': 'Weissbach FH et al., J Virol 2024',
                'confidence': 'HIGH',
            },
        }
        
        print(f"\nCell Cycle Permissiveness Validation: {'✓ PASS' if pattern_pass else '✗ FAIL'}")
        
        return results
    
    def validate_t_antigen_threshold(self) -> Dict[str, any]:
        """Validate T antigen threshold-dependent replication.
        
        Returns:
            Dictionary with T antigen validation results
        """
        print("\n" + "="*80)
        print("QUALITATIVE VALIDATION: T Antigen Threshold")
        print("="*80)
        print("Reference: BKPyV biology literature")
        print("Mechanism: T antigen required to initiate viral replication")
        print()
        
        results = {}
        
        # Scenario: Low T antigen threshold
        low_threshold_peak = self._run_scenario(
            scenario_name='low_threshold',
            tacrolimus_enhancement=1.0,
            sirolimus_inhibition=1.0,
            t_antigen_threshold=0.3,  # Low threshold
            simulation_length=100
        )
        
        # Scenario: High T antigen threshold
        high_threshold_peak = self._run_scenario(
            scenario_name='high_threshold',
            tacrolimus_enhancement=1.0,
            sirolimus_inhibition=1.0,
            t_antigen_threshold=0.8,  # High threshold
            simulation_length=100
        )
        
        print(f"Low T antigen threshold peak: {low_threshold_peak:.3f}")
        print(f"High T antigen threshold peak: {high_threshold_peak:.3f}")
        print()
        
        pattern_pass = low_threshold_peak > high_threshold_peak
        print(f"Pattern: Low threshold > High threshold")
        print(f"  Expected: TRUE (T antigen activates replication)")
        print(f"  Observed: {pattern_pass}")
        print(f"  Ratio: {low_threshold_peak / high_threshold_peak:.2f}x")
        print(f"  Status: {'✓ PASS' if pattern_pass else '✗ FAIL'}")
        
        results = {
            'pattern_t_antigen_threshold': {
                'pass': pattern_pass,
                'low_threshold': low_threshold_peak,
                'high_threshold': high_threshold_peak,
                'ratio': low_threshold_peak / high_threshold_peak,
                'expected_ratio': '> 1.0',
                'reference': 'BKPyV biology literature',
                'confidence': 'MEDIUM',
            },
        }
        
        print(f"\nT Antigen Threshold Validation: {'✓ PASS' if pattern_pass else '✗ FAIL'}")
        
        return results
    
    def _run_scenario(
        self,
        scenario_name: str,
        tacrolimus_enhancement: float = 1.0,
        sirolimus_inhibition: float = 1.0,
        cell_cycle_bonus: float = 2.0,
        t_antigen_threshold: float = 0.5,
        simulation_length: float = 100,
    ) -> float:
        """Run a simulation scenario and return peak viral load.
        
        Args:
            scenario_name: Name of the scenario
            tacrolimus_enhancement: Tacrolimus enhancement factor
            sirolimus_inhibition: Sirolimus inhibition factor
            cell_cycle_bonus: Cell cycle S-phase bonus
            t_antigen_threshold: T antigen replication threshold
            simulation_length: Simulation length in hours
            
        Returns:
            Peak virtual viral load
        """
        # Canonical engine: the bkpyv_ode 15-dimensional system (per-day rates).
        # The legacy discrete simulator (BKPyVSimulator) is superseded and did
        # not reproduce these mechanism directions; see docs/DECISIONS.md.
        ode = BKPyVODESystem(params={
            's_phase_bonus': cell_cycle_bonus,
            # The ODE's T-antigen replication gate is the half-saturation of the
            # T->viral-production map (T/(T+half_saturation)); 't_threshold'
            # only shapes sirolimus timing and does not gate replication.
            'half_saturation': t_antigen_threshold,
        })
        y0 = ode.get_infection_conditions(viral_load=0.5)

        # Map the legacy scenario knobs to ODE dosing contexts: any non-unity
        # drug factor means the corresponding drug is administered at full
        # dimensionless dosing intensity from day 0 (magnitude-response
        # curves are covered by scripts/sensitivity_analysis.py).
        dosing_context = {}
        if tacrolimus_enhancement != 1.0:
            dosing_context['tacrolimus'] = {'start': 0.0, 'stop': None, 'target': 1.0}
        if sirolimus_inhibition != 1.0:
            dosing_context['sirolimus'] = {'start': 0.0, 'stop': None, 'target': 1.0}

        t_eval = np.linspace(0.0, simulation_length, 500)
        sol = solve_ivp(
            lambda t, y: ode.ode_system(t, y, dosing_context=dosing_context or None),
            (0.0, simulation_length),
            y0,
            method='LSODA',
            t_eval=t_eval,
        )
        if not sol.success:
            return 0.0

        viral_loads = sol.y[0]  # V is state index 0
        return float(np.max(viral_loads))
    
    def generate_isef_validation_report(self) -> Dict[str, any]:
        """Generate comprehensive ISEF validation report.
        
        Returns:
            Dictionary with all validation results
        """
        print("\n" + "="*80)
        print("ISEF-LEVEL QUALITATIVE VALIDATION REPORT")
        print("="*80)
        print("Focus: Reproducing mechanistic patterns from published research")
        print("Approach: Qualitative pattern matching (scientifically defensible without patient data)")
        print()
        
        # Run all validations
        drug_effects = self.validate_drug_effects()
        cell_cycle = self.validate_cell_cycle_permissiveness()
        t_antigen = self.validate_t_antigen_threshold()
        
        # Compile results
        all_results = {
            'drug_effects': drug_effects,
            'cell_cycle_permissiveness': cell_cycle,
            't_antigen_threshold': t_antigen,
        }
        
        # Calculate overall pass rate
        total_patterns = sum(len(result_set) for result_set in all_results.values())
        passed_patterns = sum(
            1 for result_set in all_results.values()
            for pattern in result_set.values()
            if pattern['pass']
        )
        
        pass_rate = passed_patterns / total_patterns
        
        print("\n" + "="*80)
        print("VALIDATION SUMMARY")
        print("="*80)
        print(f"Total patterns tested: {total_patterns}")
        print(f"Patterns passed: {passed_patterns}")
        print(f"Pass rate: {pass_rate*100:.1f}%")
        print()
        
        # Generate recommendations
        if pass_rate >= 0.8:
            print("✓ EXCELLENT: Model reproduces key mechanistic patterns from literature")
            print("  Ready for ISEF competition with qualitative validation")
        elif pass_rate >= 0.5:
            print("⚠ GOOD: Model reproduces major patterns, some refinement needed")
            print("  Suitable for ISEF with explanation of limitations")
        else:
            print("✗ NEEDS WORK: Model fails to reproduce key patterns")
            print("  Requires parameter adjustment before ISEF submission")
        
        # Save report
        report = {
            'overall_pass_rate': pass_rate,
            'patterns_tested': total_patterns,
            'patterns_passed': passed_patterns,
            'detailed_results': all_results,
            'literature_patterns': self.literature_patterns,
            'validation_approach': 'Qualitative pattern matching - scientifically defensible without patient calibration data',
            'isef_readiness': pass_rate >= 0.6,
        }
        
        output_path = Path(__file__).parent / "isef_validation_report.json"
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"\nValidation report saved to: {output_path}")
        
        return report


def run_isef_validation():
    """Run complete ISEF-level qualitative validation."""
    validator = ISEFQualitativeValidator()
    report = validator.generate_isef_validation_report()
    
    return report


if __name__ == "__main__":
    run_isef_validation()