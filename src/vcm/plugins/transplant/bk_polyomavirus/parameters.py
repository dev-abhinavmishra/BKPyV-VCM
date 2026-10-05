"""Central, auditable BKPyV parameter registry.

Values are intentionally labelled as evidence-based or phenomenological so the
model does not imply that every coefficient is clinically calibrated.
"""
from dataclasses import dataclass
from typing import Dict, Optional

@dataclass(frozen=True)
class ParameterDefinition:
    name: str
    description: str
    default_value: float
    value_range: tuple[float, float]
    evidence_source: str
    confidence: str
    is_evidence_based: bool
    calibration_notes: str
    requires_refinement: bool = False

class BKPyVParameterRegistry:
    def __init__(self):
        self.parameters = {
            "tacrolimus_enhancement_factor": ParameterDefinition("tacrolimus_enhancement_factor", "Immune escape multiplier on infected-cell survival/production", 1.5, (1.0, 2.5), "Clinical association; not a direct intracellular replication rate", "MEDIUM", False, "Applied to immune clearance, not genome-copying flux.", True),
            "mtor_inhibition_factor": ParameterDefinition("mtor_inhibition_factor", "Residual permissiveness under sirolimus", 0.5, (0.0, 1.0), "Hirsch et al. 2016 (primary RPTECs, IC90 ~4 ng/mL); mTOR magnitudes are strongly cell-type dependent (HEK293 vs RPTEC titers differ by >1 log), so treat as direction-only", "MEDIUM", True, "Phenomenological until calibrated in primary RPTEC data.", True),
            "t_antigen_replication_threshold": ParameterDefinition("t_antigen_replication_threshold", "LT threshold for genome replication", 0.5, (0.0, 1.0), "Single-cell biology", "MEDIUM", False, "Threshold is a modeling abstraction, not a measured concentration.", True),
            "cell_cycle_s_phase_bonus": ParameterDefinition("cell_cycle_s_phase_bonus", "Relative S-phase permissiveness", 2.0, (1.0, 3.0), "Single-cell S-phase coupling", "MEDIUM", True, "Relative effect pending quantitative calibration.", True),
            "dna_replication_coupling": ParameterDefinition("dna_replication_coupling", "Coupling of host DNA synthesis to LT/replication", 0.8, (0.0, 1.0), "Single-cell and mechanistic evidence", "MEDIUM", True, "Supports the S-phase gate.", True),
            "innate_immune_suppression_factor": ParameterDefinition("innate_immune_suppression_factor", "Immune clearance sensitivity", 0.5, (0.0, 1.0), "Host-response hypothesis", "LOW", False, "Requires patient and primary-cell calibration.", True),
        }
    def get_parameter(self, name: str) -> Optional[ParameterDefinition]: return self.parameters.get(name)
    def get_all_parameters(self) -> Dict[str, ParameterDefinition]: return self.parameters
    def get_evidence_based_parameters(self): return {k:v for k,v in self.parameters.items() if v.is_evidence_based}
    def get_heuristic_parameters(self): return {k:v for k,v in self.parameters.items() if not v.is_evidence_based}
    def get_parameters_requiring_refinement(self): return {k:v for k,v in self.parameters.items() if v.requires_refinement}
    def get_parameter_dict(self): return {k:v.default_value for k,v in self.parameters.items()}

_registry = BKPyVParameterRegistry()
def get_registry(): return _registry
def get_parameter(name: str):
    item = _registry.get_parameter(name)
    return item.default_value if item else None
def get_all_parameters(): return _registry.get_parameter_dict()

def get_bkpyv_defaults():
    """Return defaults for the interactive dashboard."""
    return _registry.get_parameter_dict()
