"""Sensitivity analysis for BKPyV simulator parameters.

This module provides lightweight uncertainty quantification through Monte Carlo sampling
over key parameters defined in the BKPyV parameter registry. This adds scientific rigor by
quantifying how parameter uncertainty affects simulation outcomes.
"""

import copy
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from vcm.plugins.transplant.bk_polyomavirus.parameters import BKPyVParameterRegistry
from vcm.clinical.viral_load_mapper import ViralLoadMapper
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator
from vcm.core.models import Environment, Perturbation, PerturbationType

INFECTION_DAY = 21.0  # infection at ~week 3 post-transplant (days; ODE unit)


@dataclass
class SensitivityResult:
    """Results from sensitivity analysis."""
    parameter_name: str
    parameter_values: List[float]
    outcome_values: List[float]
    sensitivity_index: float
    correlation_coefficient: float
    confidence_interval_low: float
    confidence_interval_high: float


class SensitivityAnalyzer:
    """Lightweight sensitivity analysis for BKPyV parameters.
    
    This class performs Monte Carlo sampling over key BKPyV parameters to quantify
    how parameter uncertainty affects simulation outcomes, particularly peak viral load.
    
    Approach:
    1. Select 3-5 key parameters based on clinical importance and uncertainty
    2. Sample within their documented ranges (uniform distribution for simplicity)
    3. Run simulations for each sampled parameter value
    4. Calculate sensitivity metrics (correlation, elasticity, confidence intervals)
    5. Provide visualization-ready outputs
    
    Parameters are grounded in the BKPyV parameter registry with evidence sources.
    """
    
    def __init__(self, n_samples: int = 50, random_state: Optional[int] = None):
        """Initialize the sensitivity analyzer.
        
        Args:
            n_samples: Number of Monte Carlo samples per parameter
            random_state: Random seed for reproducibility
        """
        self.n_samples = n_samples
        self.random_state = random_state
        self.parameter_registry = BKPyVParameterRegistry()
        self.mapper = ViralLoadMapper()
        self.base_simulator_config: Dict[str, Any] = {}
        
        if random_state is not None:
            np.random.seed(random_state)
    
    def get_key_parameters(self) -> List[str]:
        """Get the 5 most clinically important BKPyV parameters.
        
        Selection based on:
        - Clinical impact (drug effects)
        - Literature support (high/medium confidence)
        - Model sensitivity (large effect on outcomes)
        
        Returns:
            List of parameter names
        """
        # Select parameters based on clinical importance and uncertainty
        key_parameters = [
            'tacrolimus_enhancement_factor',  # HIGH confidence, large clinical impact
            'mtor_inhibition_factor',           # HIGH confidence, drug effect
            'cell_cycle_s_phase_bonus',          # MEDIUM confidence, pathway coupling
            'dna_replication_coupling',           # MEDIUM confidence, core mechanism
            'innate_immune_suppression_factor'    # MEDIUM confidence, immune evasion
        ]
        return key_parameters
    
    def sample_parameter_values(self, parameter_name: str) -> List[float]:
        """Sample parameter values within documented range.
        
        Args:
            parameter_name: Name of parameter to sample
            
        Returns:
            List of sampled values
        """
        param = self.parameter_registry.parameters.get(parameter_name)
        if not param:
            raise ValueError(f"Parameter {parameter_name} not found in registry")
        
        # Uniform sampling within documented range
        min_val, max_val = param.value_range
        samples = np.random.uniform(min_val, max_val, self.n_samples)
        return samples.tolist()
    
    def run_sensitivity_analysis(self, 
                              parameter_name: str,
                              scenario: str = 'tacrolimus') -> SensitivityResult:
        """Run sensitivity analysis for a single parameter.
        
        Args:
            parameter_name: Parameter to analyze
            scenario: Scenario to run (baseline, tacrolimus, sirolimus)
            
        Returns:
            SensitivityResult with analysis metrics
        """
        valid_scenarios = {"baseline", "infection", "tacrolimus", "sirolimus"}
        if scenario not in valid_scenarios:
            raise ValueError(f"Invalid scenario: {scenario}. Must be one of {sorted(valid_scenarios)}")

        if parameter_name not in self.parameter_registry.parameters:
            raise ValueError(f"Parameter {parameter_name} not found in registry")

        # Sample parameter values
        param_values = self.sample_parameter_values(parameter_name)
        
        # Store outcomes
        peak_viral_loads = []
        
        for value in param_values:
            # Every sampled point receives its own simulator instance and
            # configuration. This is intentionally a deep copy so future
            # simulator defaults or stateful configuration cannot leak between
            # samples.
            simulator_config = copy.deepcopy(getattr(self, "base_simulator_config", {}))
            simulator_config[parameter_name] = value
            result = self._run_scenario(simulator_config, scenario)
            normalized_loads = [
                step.cell_state.metadata.get("viral_load", 0.0)
                for step in result.steps
            ]
            if not normalized_loads:
                peak_viral_loads.append(0.0)
            else:
                peak_viral_loads.append(
                    max(self.mapper.normalized_to_copies(load) for load in normalized_loads)
                )
        
        # Calculate sensitivity metrics
        sensitivity_index = self._calculate_sensitivity_index(param_values, peak_viral_loads)
        correlation = self._calculate_correlation(param_values, peak_viral_loads)
        ci_low, ci_high = self._calculate_confidence_interval(peak_viral_loads)
        
        return SensitivityResult(
            parameter_name=parameter_name,
            parameter_values=param_values,
            outcome_values=peak_viral_loads,
            sensitivity_index=sensitivity_index,
            correlation_coefficient=correlation,
            confidence_interval_low=ci_low,
            confidence_interval_high=ci_high
        )

    def _run_scenario(
        self,
        simulator_config: Dict[str, Any],
        scenario: str,
        days: int = 364,
    ):
        """Run one complete ODE trajectory for a sensitivity sample.

        Time unit is DAYS (the ODE simulator's native unit); 364 days = 52
        weeks. Infection is introduced at day 21 (week 3) and drug dosing is
        continuous from day 0.
        """
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state()
        perturbations = []

        if scenario != "baseline":
            perturbations.append(
                Perturbation(
                    id="bkpyv_infection",
                    name="BKPyV infection",
                    perturbation_type=PerturbationType.VIRAL_INFECTION,
                    target_id="viral_entry",
                    magnitude=1.0,
                    timing=INFECTION_DAY,
                )
            )
        if scenario == "tacrolimus":
            perturbations.append(
                Perturbation(
                    id="tacrolimus_treatment",
                    name="Tacrolimus treatment",
                    perturbation_type=PerturbationType.DRUG_TREATMENT,
                    target_id="FKBP1A",
                    magnitude=1.0,
                    timing=0.0,
                )
            )
        elif scenario == "sirolimus":
            perturbations.append(
                Perturbation(
                    id="sirolimus_treatment",
                    name="Sirolimus treatment",
                    perturbation_type=PerturbationType.DRUG_TREATMENT,
                    target_id="MTOR",
                    magnitude=1.0,
                    timing=0.0,
                )
            )

        return BKPyVODESimulator(simulator_config).simulate(
            initial_state=initial_state,
            perturbations=perturbations,
            environment=Environment(),
            n_steps=days,
            timestep=1.0,
        )
    
    def _calculate_sensitivity_index(self, param_values: List[float], 
                                outcome_values: List[float]) -> float:
        """Calculate elasticity sensitivity index.
        
        Elasticity: (% change in outcome) / (% change in parameter)
        
        Args:
            param_values: List of parameter values
            outcome_values: List of outcome values (peak viral loads)
            
        Returns:
            Sensitivity index
        """
        param_array = np.array(param_values)
        outcome_array = np.array(outcome_values)

        # Avoid division by zero
        if (
            param_array.size == 0
            or outcome_array.size == 0
            or param_array.std() < 1e-6
            or abs(param_array.mean()) < 1e-12
        ):
            return 0.0

        # Calculate percentage changes
        param_pct_changes = (param_array - param_array.mean()) / param_array.mean()
        if np.abs(outcome_array.mean()) < 1e-12:
            return 0.0
        outcome_pct_changes = (outcome_array - outcome_array.mean()) / outcome_array.mean()

        # Elasticity
        valid = np.abs(param_pct_changes) >= 1e-6
        if not valid.any():
            return 0.0

        elasticity = np.mean(outcome_pct_changes[valid] / param_pct_changes[valid])
        
        return float(elasticity)
    
    def _calculate_correlation(self, param_values: List[float],
                           outcome_values: List[float]) -> float:
        """Calculate Pearson correlation coefficient.
        
        Args:
            param_values: List of parameter values
            outcome_values: List of outcome values
            
        Returns:
            Correlation coefficient
        """
        if len(param_values) < 2 or np.std(param_values) < 1e-12 or np.std(outcome_values) < 1e-12:
            return 0.0
        correlation = np.corrcoef(param_values, outcome_values)[0, 1]
        return float(correlation)
    
    def _calculate_confidence_interval(self, values: List[float],
                                    confidence: float = 0.95) -> Tuple[float, float]:
        """Calculate confidence interval for outcome values.
        
        Args:
            values: List of outcome values
            confidence: Confidence level (0-1)
            
        Returns:
            Tuple of (lower_bound, upper_bound)
        """
        array = np.array(values)
        mean = np.mean(array)
        std = np.std(array)
        
        # Use t-distribution for small samples
        n = len(array)
        if n < 2:
            return (float(mean), float(mean))
        
        from scipy import stats
        t_critical = stats.t.ppf((1 + confidence) / 2, n - 1)
        
        margin = t_critical * std / np.sqrt(n)
        
        return (float(mean - margin), float(mean + margin))
    
    def run_comprehensive_analysis(self, 
                                scenario: str = 'tacrolimus',
                                parameters: Optional[List[str]] = None) -> List[SensitivityResult]:
        """Run comprehensive sensitivity analysis on multiple parameters.
        
        Args:
            scenario: Scenario to analyze
            parameters: List of parameters to analyze (if None, uses key parameters)
            
        Returns:
            List of SensitivityResult objects
        """
        if parameters is None:
            parameters = self.get_key_parameters()
        
        results = []
        for param in parameters:
            try:
                result = self.run_sensitivity_analysis(param, scenario)
                results.append(result)
            except Exception as e:
                print(f"Warning: Could not analyze parameter {param}: {e}")
                continue
        
        return results
    
    def get_summary_statistics(self, results: List[SensitivityResult]) -> pd.DataFrame:
        """Get summary statistics from sensitivity analysis results.
        
        Args:
            results: List of SensitivityResult objects
            
        Returns:
            DataFrame with summary statistics
        """
        summary_data = []
        
        for result in results:
            param_info = self.parameter_registry.parameters.get(result.parameter_name)
            summary_data.append({
                'Parameter': result.parameter_name,
                'Sensitivity Index': result.sensitivity_index,
                'Correlation': result.correlation_coefficient,
                'Confidence Interval (95%)': f"[{result.confidence_interval_low:.0f}, {result.confidence_interval_high:.0f}]",
                'Confidence': param_info.confidence if param_info else 'UNKNOWN',
                'Evidence-Based': param_info.is_evidence_based if param_info else False
            })
        
        return pd.DataFrame(summary_data)
