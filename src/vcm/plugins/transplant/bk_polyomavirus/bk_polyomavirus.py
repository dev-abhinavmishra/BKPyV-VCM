"""BK polyomavirus renal-cell plugin.

The plugin is intentionally a state/schema layer. Numerical dynamics live in
``BKPyVODESimulator`` and the discrete simulator; keeping the state definition
here makes each experiment inspectable and reproducible.
"""

from typing import Any, Dict, List, Optional

from vcm.core.models import (
    CellState,
    ExperimentConfig,
    Gene,
    Metabolite,
    Pathway,
    Perturbation,
    PerturbationType,
    Protein,
)
from vcm.plugins.base import BasePlugin


class BKPolyomavirusPlugin(BasePlugin):
    """State schema for BKPyV infection in renal tubular epithelial cells."""

    NCCR_VARIANTS = ("archetype", "rearranged")

    def __init__(self) -> None:
        super().__init__()
        self.plugin_id = "transplant.bk_polyomavirus"
        self.plugin_name = "BK Polyomavirus Kidney Cell"
        self.description = (
            "BKPyV infection in renal tubular epithelial cells with S-phase gating, "
            "immune-control drug effects, and NCCR scenarios"
        )

    def get_default_config(self) -> ExperimentConfig:
        return ExperimentConfig(
            experiment_id="bkpyv_default",
            plugin=self.plugin_id,
            simulator="bkpyv_ode",  # canonical engine; discrete legacy: "bkpyv_specific"
            simulation_length=100.0,
            timestep=1.0,
            time_unit="days",
            output_path="outputs/bkpyv/",
        )

    @staticmethod
    def _host_genes() -> Dict[str, Gene]:
        definitions = {
            "MCM2": ("MCM2 helicase", 1.0),
            "MCM5": ("MCM5 helicase", 1.0),
            "PCNA": ("PCNA", 1.2),
            "DNA_POL_ALPHA": ("DNA polymerase alpha", 1.0),
            "CDC45": ("CDC45", 1.0),
            "CCND1": ("Cyclin D1", 0.8),
            "CDK2": ("CDK2", 1.0),
            "RB1": ("RB1", 1.2),
            "TP53": ("p53", 1.0),
            "ATM": ("ATM", 0.9),
            "ATR": ("ATR", 0.9),
            "CHEK1": ("CHK1", 0.8),
            "BRCA1": ("BRCA1", 0.7),
            "BRCA2": ("BRCA2", 0.7),
            "PRKDC": ("PRKDC", 0.9),
            "FANCI": ("FANCI", 0.6),
            "MMS22L": ("MMS22L", 0.6),
            "STAT1": ("STAT1", 0.7),
            "IRF7": ("IRF7", 0.5),
            "IFNB1": ("IFNB1", 0.3),
            "BAX": ("BAX", 0.5),
            "BCL2": ("BCL2", 0.6),
            "HSPA1A": ("HSP70", 0.7),
            "FKBP1A": ("FKBP12", 1.0),
            "MTOR": ("mTOR", 1.0),
            "MT-ND4": ("MT-ND4", 0.8),
            "MT-CO1": ("MT-CO1", 0.8),
            "MT-CYB": ("MT-CYB", 0.8),
            "HLA-A": ("HLA-A", 1.0),
            "HLA-B": ("HLA-B", 1.0),
        }
        return {gene_id: Gene(id=gene_id, name=name, expression_level=level) for gene_id, (name, level) in definitions.items()}

    @staticmethod
    def _viral_genes(nccr_variant: str) -> Dict[str, Gene]:
        # Rearranged NCCR is represented as an early-gene bias scenario. The
        # ODE simulator applies the corresponding capsid trade-off explicitly.
        lt_baseline = 0.12 if nccr_variant == "rearranged" else 0.0
        definitions = {
            "viral_LT": ("Large T antigen", lt_baseline),
            "viral_ST": ("Small T antigen", 0.0),
            "viral_VP1": ("VP1", 0.0),
            "viral_VP2": ("VP2", 0.0),
            "viral_VP3": ("VP3", 0.0),
        }
        return {gene_id: Gene(id=gene_id, name=name, expression_level=level) for gene_id, (name, level) in definitions.items()}

    @staticmethod
    def _pathways() -> Dict[str, Pathway]:
        names = [
            "dna_replication", "cell_cycle", "dna_damage_response",
            "innate_immune", "apoptosis", "mTOR_signaling",
            "viral_replication", "cellular_stress", "interferon_response",
            "nucleotide_synthesis", "mitochondrial_stress",
            "antigen_presentation", "translation",
        ]
        return {
            name: Pathway(
                id=name,
                name=name.replace("_", " ").title(),
                flux=0.0 if name == "viral_replication" else 0.5,
            )
            for name in names
        }

    def create_initial_state(self, config: Optional[Dict[str, Any]] = None) -> CellState:
        config = config or {}
        nccr_variant = config.get("nccr_variant", "archetype")
        if nccr_variant not in self.NCCR_VARIANTS:
            raise ValueError(f"nccr_variant must be one of {self.NCCR_VARIANTS}")

        genes = {**self._host_genes(), **self._viral_genes(nccr_variant)}
        lt_level = genes["viral_LT"].expression_level
        proteins = {
            "LT": Protein(id="LT", name="Large T antigen", concentration=lt_level, gene_id="viral_LT", is_enzyme=True, active=lt_level > 0),
            "ST": Protein(id="ST", name="Small T antigen", concentration=0.0, gene_id="viral_ST"),
            "VP1": Protein(id="VP1", name="VP1", concentration=0.0, gene_id="viral_VP1"),
        }
        metabolites = {
            nucleotide: Metabolite(id=nucleotide, name=nucleotide, concentration=0.5, compartment="nucleus", is_essential=True)
            for nucleotide in ("dATP", "dCTP", "dGTP", "dTTP")
        }
        metabolites["atp"] = Metabolite(id="atp", name="ATP", concentration=3.0, is_essential=True)

        pathways = self._pathways()
        pathway_activities = {name: pathway.flux for name, pathway in pathways.items()}
        metadata = {
            "nccr_variant": nccr_variant,
            "nccr_description": (
                "Rearranged NCCR: hypothesized early-gene bias"
                if nccr_variant == "rearranged"
                else "Archetype NCCR: presumed persistent/transmissible baseline"
            ),
            "cell_cycle_phase": "G0/G1",
            "host_dna_synthesis": "inactive",
            "dna_synthesis_status": "inactive",
            "t_antigen_level": "none",
            "viral_load": 0.0,
            "infection_status": "uninfected",
            "viral_replication_phase": "none",
            "time_since_infection": 0.0,
            "pathway_activities": pathway_activities,
            "drug_effects": {"tacrolimus": 0.0, "sirolimus": 0.0, "everolimus": 0.0},
            "drug_exposure_duration": {"tacrolimus": 0.0, "sirolimus": 0.0, "everolimus": 0.0},
            "mitochondrial_stress": 0.0,
            "antigen_presentation": 1.0,
            "innate_immune_suppression": 0.0,
            "healthy_cell_count": 1.0,
            "infected_cell_count": 0.0,
            "dead_cell_count": 0.0,
        }
        return CellState(
            cell_id="bkpyv_kidney_cell_001",
            cell_type="kidney_tubular_epithelial",
            genes=genes,
            proteins=proteins,
            metabolites=metabolites,
            pathways=pathways,
            metadata=metadata,
        )

    def get_supported_perturbations(self) -> List[Perturbation]:
        return [
            Perturbation(id="bkpyv_infection", name="BKPyV infection", perturbation_type=PerturbationType.VIRAL_INFECTION, magnitude=1.0),
            Perturbation(id="tacrolimus_treatment", name="Tacrolimus treatment (immune escape)", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="FKBP1A", magnitude=1.0),
            Perturbation(id="sirolimus_treatment", name="Sirolimus treatment (S-phase gate)", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="MTOR", magnitude=1.0),
            Perturbation(id="antiviral_treatment", name="Antiviral treatment", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="viral_LT", magnitude=1.0),
            Perturbation(id="dna_damage", name="DNA damage", perturbation_type=PerturbationType.ENVIRONMENTAL_CHANGE, magnitude=0.8),
        ]

    def get_cell_schema(self) -> Dict[str, Any]:
        schema = super().get_cell_schema()
        schema.update(
            {
                "viral_genes": ["viral_LT", "viral_ST", "viral_VP1", "viral_VP2", "viral_VP3"],
                "pathway_activities": ["dna_replication", "cell_cycle", "viral_replication"],
                "cell_cycle_phases": ["G0/G1", "S", "G2/M"],
                "drug_targets": ["FKBP1A", "MTOR"],
                "key_pathways": ["dna_replication", "cell_cycle", "dna_damage_response", "innate_immune", "viral_replication"],
                "nccr_variants": list(self.NCCR_VARIANTS),
            }
        )
        return schema
