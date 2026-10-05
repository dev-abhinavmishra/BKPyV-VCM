"""Clinical Thresholds Management Module for BKPyV.

This module provides comprehensive management of clinical thresholds for BK viremia
monitoring, risk stratification, and clinical decision-making. It incorporates
guidelines from KDIGO, AST consensus statements, and published clinical studies.

Research Grounding:
- KDIGO Clinical Practice Guidelines for Kidney Transplant Recipients
- American Society of Transplantation (AST) consensus statements
- Published cohort studies on BK viremia thresholds
- Clinical screening protocols and monitoring intervals

Threshold Levels:
- Screening Positive: ≥1,000 copies/mL
- High Risk: ≥10,000 copies/mL
- Clinical Nephropathy Risk: ≥10,000 copies/mL (persistent)
- Intervention Threshold: ≥1,000 copies/mL (rising)
"""

import math
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum
import json
from pathlib import Path


class RiskCategory(Enum):
    """BKPyV risk categories based on viral load."""
    NEGATIVE = "negative"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class InterventionLevel(Enum):
    """Clinical intervention levels based on risk."""
    MONITORING = "monitoring"
    CONSIDER_INTERVENTION = "consider_intervention"
    ACTIVE_INTERVENTION = "active_intervention"
    URGENT_INTERVENTION = "urgent_intervention"


@dataclass
class ClinicalThreshold:
    """Clinical threshold for BK viremia monitoring."""
    
    name: str
    threshold_copies_ml: float
    risk_category: RiskCategory
    intervention_level: InterventionLevel
    clinical_significance: str
    guideline_reference: str
    monitoring_frequency: str  # e.g., "Monthly", "Biweekly", "Weekly"
    action_recommendations: List[str]
    evidence_level: str  # "A", "B", "C" based on strength of evidence


@dataclass
class ViralLoadResult:
    """Result of viral load interpretation."""
    
    viral_load_copies_ml: float
    viral_load_log10: float
    risk_category: RiskCategory
    intervention_level: InterventionLevel
    threshold_crossings: List[str]
    clinical_recommendations: List[str]
    monitoring_frequency: str
    clinical_context: str


class ClinicalThresholdsManager:
    """Manager for clinical thresholds and risk stratification.
    
    This class provides:
    1. Standard clinical thresholds from guidelines
    2. Risk stratification based on viral load
    3. Clinical intervention recommendations
    4. Monitoring frequency guidelines
    5. Trend analysis and interpretation
    """
    
    def __init__(self):
        """Initialize clinical thresholds manager with guideline-based thresholds."""
        self.thresholds = self._initialize_thresholds()
        self.guideline_references = self._initialize_guidelines()
    
    def _initialize_thresholds(self) -> Dict[str, ClinicalThreshold]:
        """Initialize clinical thresholds from published guidelines.
        
        Returns:
            Dictionary of clinical thresholds
        """
        return {
            'screening_positive': ClinicalThreshold(
                name="Screening Positive",
                threshold_copies_ml=1000.0,
                risk_category=RiskCategory.LOW,
                intervention_level=InterventionLevel.MONITORING,
                clinical_significance="First detectable viremia, requires close monitoring",
                guideline_reference="KDIGO 2023, AST 2019",
                monitoring_frequency="Biweekly",
                action_recommendations=[
                    "Increase monitoring frequency to biweekly",
                    "Review immunosuppression regimen",
                    "Assess for other causes of renal dysfunction",
                    "Consider baseline renal function trends",
                ],
                evidence_level="A",
            ),
            'high_risk': ClinicalThreshold(
                name="High Risk",
                threshold_copies_ml=10000.0,
                risk_category=RiskCategory.HIGH,
                intervention_level=InterventionLevel.ACTIVE_INTERVENTION,
                clinical_significance="Significant viremia, requires intervention to prevent nephropathy",
                guideline_reference="KDIGO 2023, AST 2019, Multiple cohort studies",
                monitoring_frequency="Weekly",
                action_recommendations=[
                    "Reduce immunosuppression (calcineurin inhibitor reduction)",
                    "Consider switching from tacrolimus to sirolimus",
                    "Weekly monitoring of viral load",
                    "Close renal function monitoring",
                    "Early nephrology consultation",
                    "Consider renal biopsy if renal function declines",
                ],
                evidence_level="A",
            ),
            'critical': ClinicalThreshold(
                name="Critical",
                threshold_copies_ml=100000.0,
                risk_category=RiskCategory.CRITICAL,
                intervention_level=InterventionLevel.URGENT_INTERVENTION,
                clinical_significance="Very high viremia, urgent intervention required",
                guideline_reference="Clinical experience, case series",
                monitoring_frequency="Twice weekly",
                action_recommendations=[
                    "Significant immunosuppression reduction",
                    "Switch to mTOR inhibitor (sirolimus)",
                    "Twice weekly viral load monitoring",
                    "Immediate nephrology consultation",
                    "Consider renal biopsy for definitive diagnosis",
                    "Evaluate for adjunctive antiviral therapy",
                ],
                evidence_level="B",
            ),
            'declining': ClinicalThreshold(
                name="Declining Trend",
                threshold_copies_ml=0.0,  # Dynamic threshold
                risk_category=RiskCategory.LOW,
                intervention_level=InterventionLevel.MONITORING,
                clinical_significance="Viral load declining, intervention effective",
                guideline_reference="Clinical practice",
                monitoring_frequency="Monthly",
                action_recommendations=[
                    "Continue current immunosuppression reduction",
                    "Return to monthly monitoring",
                    "Monitor for viral rebound",
                    "Assess renal function recovery",
                ],
                evidence_level="B",
            ),
        }
    
    def _initialize_guidelines(self) -> Dict[str, str]:
        """Initialize guideline references.
        
        Returns:
            Dictionary of guideline references
        """
        return {
            'KDIGO_2023': "KDIGO Clinical Practice Guidelines for Kidney Transplant Recipients (2023)",
            'AST_2019': "American Society of Transplantation Consensus Statement on BK Virus (2019)",
            'AST_2021': "American Society of Transplantation Guidelines for BK Virus Nephropathy (2021)",
            'KDIGO_2024': "KDIGO Clinical Practice Guidelines for the Care of Kidney Transplant Recipients (2024)",
        }
    
    def interpret_viral_load(
        self,
        viral_load_copies_ml: float,
        trend: Optional[str] = None,
        prior_values: Optional[List[float]] = None,
        immunosuppression_status: Optional[str] = None,
    ) -> ViralLoadResult:
        """Interpret viral load result with clinical context.
        
        Args:
            viral_load_copies_ml: Viral load in copies/mL
            trend: Viral load trend ('rising', 'falling', 'stable', None)
            prior_values: Prior viral load measurements for trend analysis
            immunosuppression_status: Current immunosuppression status
            
        Returns:
            ViralLoadResult with interpretation
        """
        # Calculate log10 viral load
        viral_load_log10 = math.log10(viral_load_copies_ml) if viral_load_copies_ml > 0 else 0.0
        
        # Determine risk category and intervention level
        risk_category, intervention_level = self._determine_risk_category(
            viral_load_copies_ml, trend
        )
        
        # Identify threshold crossings
        threshold_crossings = self._identify_threshold_crossings(viral_load_copies_ml)
        
        # Generate clinical recommendations
        clinical_recommendations = self._generate_recommendations(
            risk_category, intervention_level, trend, immunosuppression_status
        )
        
        # Determine monitoring frequency
        monitoring_frequency = self._get_monitoring_frequency(intervention_level)
        
        # Clinical context
        clinical_context = self._generate_clinical_context(
            viral_load_copies_ml, trend, prior_values
        )
        
        return ViralLoadResult(
            viral_load_copies_ml=viral_load_copies_ml,
            viral_load_log10=viral_load_log10,
            risk_category=risk_category,
            intervention_level=intervention_level,
            threshold_crossings=threshold_crossings,
            clinical_recommendations=clinical_recommendations,
            monitoring_frequency=monitoring_frequency,
            clinical_context=clinical_context,
        )
    
    def _determine_risk_category(
        self,
        viral_load: float,
        trend: Optional[str]
    ) -> Tuple[RiskCategory, InterventionLevel]:
        """Determine risk category and intervention level.
        
        Args:
            viral_load: Viral load in copies/mL
            trend: Viral load trend
            
        Returns:
            Tuple of (risk_category, intervention_level)
        """
        if viral_load >= self.thresholds['critical'].threshold_copies_ml:
            return RiskCategory.CRITICAL, InterventionLevel.URGENT_INTERVENTION
        
        if viral_load >= self.thresholds['high_risk'].threshold_copies_ml:
            return RiskCategory.HIGH, InterventionLevel.ACTIVE_INTERVENTION
        
        if viral_load >= self.thresholds['screening_positive'].threshold_copies_ml:
            if trend == 'rising':
                return RiskCategory.MEDIUM, InterventionLevel.CONSIDER_INTERVENTION
            else:
                return RiskCategory.LOW, InterventionLevel.MONITORING
        
        return RiskCategory.NEGATIVE, InterventionLevel.MONITORING
    
    def _identify_threshold_crossings(self, viral_load: float) -> List[str]:
        """Identify which thresholds have been crossed.
        
        Args:
            viral_load: Viral load in copies/mL
            
        Returns:
            List of threshold crossing descriptions
        """
        crossings = []
        
        if viral_load >= self.thresholds['critical'].threshold_copies_ml:
            crossings.append(f"Critical threshold ({self.thresholds['critical'].threshold_copies_ml:.0f} copies/mL) crossed")
        
        if viral_load >= self.thresholds['high_risk'].threshold_copies_ml:
            crossings.append(f"High-risk threshold ({self.thresholds['high_risk'].threshold_copies_ml:.0f} copies/mL) crossed")
        
        if viral_load >= self.thresholds['screening_positive'].threshold_copies_ml:
            crossings.append(f"Screening positive threshold ({self.thresholds['screening_positive'].threshold_copies_ml:.0f} copies/mL) crossed")
        
        return crossings
    
    def _generate_recommendations(
        self,
        risk_category: RiskCategory,
        intervention_level: InterventionLevel,
        trend: Optional[str],
        immunosuppression_status: Optional[str]
    ) -> List[str]:
        """Generate clinical recommendations based on risk.
        
        Args:
            risk_category: Current risk category
            intervention_level: Required intervention level
            trend: Viral load trend
            immunosuppression_status: Current immunosuppression status
            
        Returns:
            List of clinical recommendations
        """
        recommendations = []
        
        # Get base recommendations from threshold
        if risk_category == RiskCategory.HIGH:
            recommendations.extend(self.thresholds['high_risk'].action_recommendations)
        elif risk_category == RiskCategory.CRITICAL:
            recommendations.extend(self.thresholds['critical'].action_recommendations)
        elif risk_category == RiskCategory.LOW:
            recommendations.extend(self.thresholds['screening_positive'].action_recommendations)
        
        # Add trend-specific recommendations
        if trend == 'rising' and risk_category != RiskCategory.CRITICAL:
            recommendations.append("Consider more aggressive immunosuppression reduction")
            recommendations.append("Monitor more frequently than standard recommendations")
        elif trend == 'falling':
            recommendations.append("Continue current management strategy")
            recommendations.append("Monitor for viral rebound")
        
        # Add immunosuppression-specific recommendations
        if immunosuppression_status == 'tacrolimus' and risk_category in [RiskCategory.HIGH, RiskCategory.CRITICAL]:
            recommendations.append("Strongly consider switching from tacrolimus to sirolimus")
            recommendations.append("Tacrolimus may be enhancing viral replication (OR 2.0-2.3)")
        
        return recommendations
    
    def _get_monitoring_frequency(self, intervention_level: InterventionLevel) -> str:
        """Get recommended monitoring frequency.
        
        Args:
            intervention_level: Required intervention level
            
        Returns:
            Monitoring frequency string
        """
        frequency_map = {
            InterventionLevel.MONITORING: "Monthly",
            InterventionLevel.CONSIDER_INTERVENTION: "Biweekly",
            InterventionLevel.ACTIVE_INTERVENTION: "Weekly",
            InterventionLevel.URGENT_INTERVENTION: "Twice weekly",
        }
        
        return frequency_map.get(intervention_level, "Monthly")
    
    def _generate_clinical_context(
        self,
        viral_load: float,
        trend: Optional[str],
        prior_values: Optional[List[float]]
    ) -> str:
        """Generate clinical context description.
        
        Args:
            viral_load: Current viral load
            trend: Viral load trend
            prior_values: Prior viral load measurements
            
        Returns:
            Clinical context description
        """
        context = f"Viral load: {viral_load:.0f} copies/mL"
        
        if trend:
            context += f" ({trend} trend)"
        
        if prior_values and len(prior_values) > 1:
            # Calculate rate of change
            prior = prior_values[-2]
            if prior > 0:
                change_rate = (viral_load - prior) / prior
                if change_rate > 0:
                    context += f", increasing by {change_rate*100:.1f}% from prior measurement"
                else:
                    context += f", decreasing by {abs(change_rate)*100:.1f}% from prior measurement"
        
        return context
    
    def analyze_trend(
        self,
        viral_loads: List[float],
        timepoints: List[float],
        min_measurements: int = 3
    ) -> Dict[str, Any]:
        """Analyze viral load trend over time.
        
        Args:
            viral_loads: List of viral load measurements
            timepoints: Corresponding time points
            min_measurements: Minimum measurements required for trend analysis
            
        Returns:
            Dictionary with trend analysis results
        """
        if len(viral_loads) < min_measurements:
            return {
                'trend': 'insufficient_data',
                'slope': None,
                'rate_of_change': None,
                'doubling_time': None,
                'recommendation': 'Insufficient data for trend analysis',
            }
        
        # Calculate simple linear trend
        import numpy as np
        
        timepoints_array = np.array(timepoints)
        viral_loads_array = np.array(viral_loads)
        
        # Log-transform viral loads for linear analysis
        log_loads = np.log10(viral_loads_array + 1)  # Add 1 to avoid log(0)
        
        # Linear regression
        slope, intercept = np.polyfit(timepoints_array, log_loads, 1)
        
        # Determine trend direction
        if slope > 0.01:  # Positive slope indicates rising
            trend = 'rising'
        elif slope < -0.01:  # Negative slope indicates falling
            trend = 'falling'
        else:
            trend = 'stable'
        
        # Calculate rate of change (guard against degenerate time axis)
        if len(viral_loads) >= 2 and timepoints[-1] != timepoints[0]:
            rate_of_change = (viral_loads[-1] - viral_loads[0]) / (timepoints[-1] - timepoints[0])
        else:
            rate_of_change = None
        
        # Calculate doubling time if rising. The slope was fitted to
        # log10-transformed loads, so doubling time uses log10(2), not ln(2).
        doubling_time = None
        if trend == 'rising' and slope > 0:
            doubling_time = math.log10(2) / slope
        
        # Generate recommendation
        if trend == 'rising' and slope > 0.05:  # Steep rise
            recommendation = "Rapidly rising viral load - urgent intervention recommended"
        elif trend == 'rising':
            recommendation = "Rising viral load - consider intervention"
        elif trend == 'falling':
            recommendation = "Falling viral load - continue current strategy"
        else:
            recommendation = "Stable viral load - continue monitoring"
        
        return {
            'trend': trend,
            'slope': slope,
            'rate_of_change': rate_of_change,
            'doubling_time': doubling_time,
            'recommendation': recommendation,
        }
    
    def generate_monitoring_schedule(
        self,
        risk_category: RiskCategory,
        intervention_level: InterventionLevel,
        post_transplant_days: int
    ) -> Dict[str, str]:
        """Generate recommended monitoring schedule.
        
        Args:
            risk_category: Current risk category
            intervention_level: Required intervention level
            post_transplant_days: Days post-transplant
            
        Returns:
            Dictionary with monitoring schedule recommendations
        """
        # Base monitoring frequency
        base_frequency = self._get_monitoring_frequency(intervention_level)
        
        # Adjust for post-transplant timeline
        if post_transplant_days < 90:  # First 3 months
            schedule_note = "Enhanced monitoring in early post-transplant period"
            if intervention_level == InterventionLevel.MONITORING:
                base_frequency = "Biweekly"  # More frequent early on
        elif post_transplant_days < 365:  # First year
            schedule_note = "Standard first-year monitoring"
        else:
            schedule_note = "Long-term monitoring protocol"
        
        # Add specific recommendations
        if risk_category in [RiskCategory.HIGH, RiskCategory.CRITICAL]:
            additional_tests = [
                "Serum creatinine",
                "eGFR monitoring",
                "Urine protein",
                "Consider renal biopsy if renal function declines",
            ]
        else:
            additional_tests = [
                "Serum creatinine",
                "eGFR monitoring",
            ]
        
        return {
            'monitoring_frequency': base_frequency,
            'schedule_note': schedule_note,
            'additional_tests': additional_tests,
            'next_monitoring': self._calculate_next_monitoring_date(base_frequency),
        }
    
    def _calculate_next_monitoring_date(self, frequency: str) -> str:
        """Calculate next monitoring date based on frequency.
        
        Args:
            frequency: Monitoring frequency string
            
        Returns:
            Next monitoring date description
        """
        frequency_days = {
            'Twice weekly': 3.5,
            'Weekly': 7,
            'Biweekly': 14,
            'Monthly': 30,
        }
        
        days = frequency_days.get(frequency, 30)
        return f"In {days} days (or {frequency.lower()})"
    
    def export_thresholds_config(self, output_path: str) -> None:
        """Export thresholds configuration to file.
        
        Args:
            output_path: Path to save configuration
        """
        config = {
            'thresholds': {},
            'guideline_references': self.guideline_references,
        }
        
        for threshold_name, threshold in self.thresholds.items():
            config['thresholds'][threshold_name] = {
                'name': threshold.name,
                'threshold_copies_ml': threshold.threshold_copies_ml,
                'risk_category': threshold.risk_category.value,
                'intervention_level': threshold.intervention_level.value,
                'clinical_significance': threshold.clinical_significance,
                'guideline_reference': threshold.guideline_reference,
                'monitoring_frequency': threshold.monitoring_frequency,
                'action_recommendations': threshold.action_recommendations,
                'evidence_level': threshold.evidence_level,
            }
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"Thresholds configuration exported to {output_file}")
    
    def import_thresholds_config(self, input_path: str) -> None:
        """Import thresholds configuration from file.
        
        Args:
            input_path: Path to load configuration from
        """
        input_file = Path(input_path)
        
        if not input_file.exists():
            raise FileNotFoundError(f"Configuration file not found: {input_path}")
        
        with open(input_file, 'r') as f:
            config = json.load(f)
        
        # Reconstruct thresholds
        for threshold_name, threshold_data in config['thresholds'].items():
            self.thresholds[threshold_name] = ClinicalThreshold(
                name=threshold_data['name'],
                threshold_copies_ml=threshold_data['threshold_copies_ml'],
                risk_category=RiskCategory(threshold_data['risk_category']),
                intervention_level=InterventionLevel(threshold_data['intervention_level']),
                clinical_significance=threshold_data['clinical_significance'],
                guideline_reference=threshold_data['guideline_reference'],
                monitoring_frequency=threshold_data['monitoring_frequency'],
                action_recommendations=threshold_data['action_recommendations'],
                evidence_level=threshold_data['evidence_level'],
            )
        
        print(f"Thresholds configuration imported from {input_path}")


# Convenience functions
def get_default_thresholds_manager() -> ClinicalThresholdsManager:
    """Get default clinical thresholds manager.
    
    Returns:
        ClinicalThresholdsManager with guideline-based thresholds
    """
    return ClinicalThresholdsManager()


def interpret_bk_viral_load(
    viral_load_copies_ml: float,
    trend: Optional[str] = None,
) -> Dict[str, Any]:
    """Convenience function to interpret BK viral load.
    
    Args:
        viral_load_copies_ml: Viral load in copies/mL
        trend: Viral load trend (optional)
        
    Returns:
        Dictionary with interpretation results
    """
    manager = get_default_thresholds_manager()
    result = manager.interpret_viral_load(viral_load_copies_ml, trend)
    
    return {
        'viral_load': result.viral_load_copies_ml,
        'viral_load_log': result.viral_load_log10,
        'risk_category': result.risk_category.value,
        'intervention_level': result.intervention_level.value,
        'threshold_crossings': result.threshold_crossings,
        'recommendations': result.clinical_recommendations,
        'monitoring_frequency': result.monitoring_frequency,
        'clinical_context': result.clinical_context,
    }


if __name__ == "__main__":
    # Example usage
    manager = ClinicalThresholdsManager()
    
    print("Clinical Thresholds Manager")
    print("=" * 80)
    
    # Test interpretation
    test_cases = [
        (500.0, "stable"),
        (1500.0, "rising"),
        (15000.0, "rising"),
        (50000.0, "rising"),
    ]
    
    print("\nViral Load Interpretation Examples:")
    for viral_load, trend in test_cases:
        result = manager.interpret_viral_load(viral_load, trend)
        print(f"\nViral Load: {viral_load:.0f} copies/mL ({trend})")
        print(f"Risk Category: {result.risk_category.value.upper()}")
        print(f"Intervention Level: {result.intervention_level.value}")
        print(f"Monitoring: {result.monitoring_frequency}")
        print(f"Context: {result.clinical_context}")
    
    # Export configuration
    manager.export_thresholds_config("clinical_thresholds_config.json")
    
    print("\nClinical thresholds configuration exported to clinical_thresholds_config.json")