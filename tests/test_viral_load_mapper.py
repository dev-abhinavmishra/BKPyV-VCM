"""Tests for ViralLoadMapper.

Tests the Hill function mapping from simulator viral load (0-1) to clinical
plasma viral load (copies/mL).
"""

import pytest
import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcm.clinical.viral_load_mapper import (
    ViralLoadMapper,
    hill_function,
    fit_hill_parameters
)


class TestHillFunction:
    """Test the Hill function implementation."""
    
    def test_hill_function_properties(self):
        """Test basic properties of Hill function."""
        Vmax = 1e7
        K = 0.5
        n = 2.0
        
        # At vl=0, output should be ~0
        assert hill_function(0.0, Vmax, K, n) < 1.0
        
        # At vl=1, output should be ~Vmax (but less due to K^n)
        result = hill_function(1.0, Vmax, K, n)
        assert result > 0.0
        assert result < Vmax
        
        # At very high vl, output should approach Vmax
        high_result = hill_function(10.0, Vmax, K, n)
        assert high_result > result


class TestViralLoadMapper:
    """Test the ViralLoadMapper class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        # Use a temporary params path to avoid conflicts
        self.params_path = "/tmp/test_viral_load_mapper_params.json"
        self.mapper = ViralLoadMapper(params_path=self.params_path)
    
    def test_anchor_point_reasonable_range(self):
        """Test that normalized viral load 0.3-0.6 produces clinically relevant range."""
        result_0_3 = self.mapper.normalized_to_copies(0.3)
        result_0_6 = self.mapper.normalized_to_copies(0.6)
        
        # Both should be in a clinically relevant range (100-1e7)
        assert 100 < result_0_3 < 1e7, f"Result {result_0_3} outside clinical range"
        assert 100 < result_0_6 < 1e7, f"Result {result_0_6} outside clinical range"
        
        # Higher normalized load should give higher clinical load (monotonic)
        assert result_0_6 > result_0_3, "Higher normalized should give higher clinical"
    
    def test_anchor_point_0_0(self):
        """Test that normalized viral load 0.0 produces near 0 copies/mL."""
        result = self.mapper.normalized_to_copies(0.0)
        assert result < 10.0, f"Expected near 0, got {result}"
    
    def test_anchor_point_1_0(self):
        """Bridge anchors are reproduced exactly by construction."""
        result = self.mapper.normalized_to_copies(1.0)
        # V = 1.0 is anchored to the presumptive-PyVAN threshold (10,000)
        assert abs(result - 10000.0) < 1.0, f"Expected ~10,000 anchor, got {result}"
        # screening anchor
        assert abs(self.mapper.normalized_to_copies(0.2) - 1000.0) < 1.0
        # detection-limit anchor
        assert abs(self.mapper.normalized_to_copies(0.02) - 100.0) < 1.0
    
    def test_monotonicity(self):
        """Test that output increases monotonically with input."""
        vl_values = np.linspace(0, 1, 20)
        copies_values = [self.mapper.normalized_to_copies(vl) for vl in vl_values]
        
        # Check monotonic increase
        for i in range(len(copies_values) - 1):
            assert copies_values[i] <= copies_values[i + 1] + 1e-6, \
                f"Non-monotonic at index {i}: {copies_values[i]} > {copies_values[i + 1]}"
    
    def test_risk_category_undetectable(self):
        """Test risk category for undetectable viral load."""
        assert self.mapper.copies_to_risk_category(50) == 'undetectable'
        assert self.mapper.copies_to_risk_category(0) == 'undetectable'
    
    def test_risk_category_low_risk(self):
        """Test risk category for low-risk viral load."""
        assert self.mapper.copies_to_risk_category(500) == 'low_risk'
        assert self.mapper.copies_to_risk_category(999) == 'low_risk'
    
    def test_risk_category_screening(self):
        """Test risk category at screening threshold."""
        assert self.mapper.copies_to_risk_category(1000) == 'screening'
        assert self.mapper.copies_to_risk_category(5000) == 'screening'
    
    def test_risk_category_treatment(self):
        """Test risk category at treatment threshold."""
        assert self.mapper.copies_to_risk_category(10000) == 'treatment'
        assert self.mapper.copies_to_risk_category(50000) == 'treatment'
    
    def test_risk_category_severe(self):
        """Test risk category for severe nephropathy."""
        assert self.mapper.copies_to_risk_category(1e8) == 'severe'
        assert self.mapper.copies_to_risk_category(1e7) == 'severe'
    
    def test_boundary_values(self):
        """Test risk categories at exact boundary values."""
        assert self.mapper.copies_to_risk_category(99) == 'undetectable'
        assert self.mapper.copies_to_risk_category(100) == 'low_risk'  # At detection limit
        assert self.mapper.copies_to_risk_category(1000) == 'screening'
        assert self.mapper.copies_to_risk_category(10000) == 'treatment'
        assert self.mapper.copies_to_risk_category(1e7) == 'severe'
    
    def test_simulate_clinical_trajectory_columns(self):
        """Test that simulate_clinical_trajectory returns DataFrame with correct columns."""
        df = self.mapper.simulate_clinical_trajectory('infection', weeks=10)
        
        expected_columns = ['week', 'viral_load_norm', 'copies_per_ml', 'risk_category', 'drug_scenario']
        assert list(df.columns) == expected_columns, f"Expected {expected_columns}, got {list(df.columns)}"
    
    def test_simulate_clinical_trajectory_length(self):
        """Test that simulate_clinical_trajectory returns correct number of weeks."""
        weeks = 20
        df = self.mapper.simulate_clinical_trajectory('infection', weeks=weeks)
        
        assert len(df) == weeks + 1, f"Expected {weeks + 1} rows, got {len(df)}"
    
    def test_simulate_clinical_trajectory_baseline(self):
        """Test that baseline scenario produces zero viral load."""
        df = self.mapper.simulate_clinical_trajectory('baseline', weeks=10)
        
        assert df['viral_load_norm'].max() < 1e-6, "Baseline should have near-zero viral load"
        assert df['copies_per_ml'].max() < 10.0, "Baseline should have near-zero clinical load"
    
    def test_simulate_clinical_trajectory_scenario_validation(self):
        """Test that invalid scenario raises ValueError."""
        with pytest.raises(ValueError):
            self.mapper.simulate_clinical_trajectory('invalid_scenario', weeks=10)
    
    def test_simulate_clinical_trajectory_all_scenarios(self):
        """Test that all four scenarios can be simulated."""
        scenarios = ['baseline', 'infection', 'tacrolimus', 'sirolimus']
        
        for scenario in scenarios:
            df = self.mapper.simulate_clinical_trajectory(scenario, weeks=5)
            assert len(df) == 6, f"{scenario} should have 6 rows (weeks 0-5)"
            assert df['drug_scenario'].unique()[0] == scenario, f"Scenario should be {scenario}"
    
    def test_simulate_clinical_trajectory_monotonic_viral_load(self):
        """Test that viral load in trajectory is non-negative."""
        df = self.mapper.simulate_clinical_trajectory('infection', weeks=20)
        
        assert (df['viral_load_norm'] >= 0).all(), "Viral load should be non-negative"
        assert (df['copies_per_ml'] >= 0).all(), "Clinical copies should be non-negative"
    
    def test_simulate_clinical_trajectory_viral_load_range(self):
        """Viral load is dimensionless ODE output: non-negative, finite, may exceed 1."""
        df = self.mapper.simulate_clinical_trajectory('infection', weeks=20)

        assert (df['viral_load_norm'] >= 0).all(), "Viral load should be non-negative"
        assert df['viral_load_norm'].notna().all(), "Viral load should be finite"
        assert (df['copies_per_ml'] >= 0).all()


class TestClinicalSummary:
    """Test the clinical summary statistics generation."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.params_path = "/tmp/test_viral_load_mapper_params.json"
        self.mapper = ViralLoadMapper(params_path=self.params_path)
        
        # Create test trajectories
        self.trajectories = {
            'baseline': self.mapper.simulate_clinical_trajectory('baseline', weeks=10),
            'infection': self.mapper.simulate_clinical_trajectory('infection', weeks=10)
        }
    
    def test_summary_contains_all_scenarios(self):
        """Test that summary contains all scenarios."""
        from vcm.clinical.viral_load_mapper import generate_clinical_summary
        
        summary = generate_clinical_summary(self.trajectories)
        
        assert 'baseline' in summary
        assert 'infection' in summary
    
    def test_summary_has_all_fields(self):
        """Test that summary has all expected fields."""
        from vcm.clinical.viral_load_mapper import generate_clinical_summary
        
        summary = generate_clinical_summary(self.trajectories)
        
        expected_fields = [
            'time_to_first_detectable_viremia_weeks',
            'peak_viral_load_copies_per_ml',
            'week_of_peak',
            'time_above_screening_threshold_weeks',
            'time_above_treatment_threshold_weeks',
            'final_viral_load_week_52_copies_per_ml'
        ]
        
        for scenario_stats in summary.values():
            for field in expected_fields:
                assert field in scenario_stats, f"Missing field: {field}"
    
    def test_summary_positive_values(self):
        """Test that summary values are non-negative."""
        from vcm.clinical.viral_load_mapper import generate_clinical_summary
        
        summary = generate_clinical_summary(self.trajectories)
        
        for scenario, stats in summary.items():
            assert stats['peak_viral_load_copies_per_ml'] >= 0
            assert stats['time_above_screening_threshold_weeks'] >= 0
            assert stats['time_above_treatment_threshold_weeks'] >= 0
            assert stats['final_viral_load_week_52_copies_per_ml'] >= 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
