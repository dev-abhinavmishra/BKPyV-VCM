#!/usr/bin/env python3
"""
ARC-style evaluation metrics for virtual cell model validation.

Implements metrics inspired by the ARC Virtual Cell Challenge:
- Differential Expression Score (DES)
- Perturbation Discrimination Score (PDS) 
- Mean Absolute Error (MAE)
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


@dataclass
class EvaluationResults:
    """Results from ARC-style evaluation metrics."""
    des: float  # Differential Expression Score
    pds: float  # Perturbation Discrimination Score
    mae: float  # Mean Absolute Error
    details: Dict[str, Any]


class ARCEvaluator:
    """ARC-style evaluator for virtual cell model validation."""
    
    def __init__(self):
        self.gene_weights = {
            'viral_genes': 2.0,  # Higher weight for viral response genes
            'immune_genes': 1.5,  # Higher weight for immune response genes
            'housekeeping': 1.0,  # Standard weight for housekeeping genes
        }
    
    def calculate_des(
        self,
        predicted_de: Dict[str, float],
        actual_de: Dict[str, float],
        top_n: int = 100
    ) -> float:
        """
        Calculate Differential Expression Score (DES).
        
        DES captures whether the model recovers the correct set of 
        differentially expressed genes.
        
        Args:
            predicted_de: Predicted differentially expressed genes {gene: fold_change}
            actual_de: Actual differentially expressed genes {gene: fold_change}
            top_n: Number of top genes to consider
            
        Returns:
            DES score (0-1, higher is better)
        """
        if not predicted_de or not actual_de or top_n <= 0:
            return 0.0

        # Get top N predicted and actual genes by absolute fold change. The
        # gene-name tie-break makes results reproducible across runs.
        predicted_top = sorted(
            predicted_de.items(), 
            key=lambda x: (-abs(x[1]), str(x[0])),
        )[:top_n]
        actual_top = sorted(
            actual_de.items(), 
            key=lambda x: (-abs(x[1]), str(x[0])),
        )[:top_n]
        
        predicted_genes = {gene for gene, _ in predicted_top}
        actual_genes = {gene for gene, _ in actual_top}
        
        # Calculate overlap
        overlap = predicted_genes & actual_genes
        
        denominator = min(top_n, len(predicted_genes), len(actual_genes))
        des = len(overlap) / denominator if denominator else 0.0
        
        return des
    
    def calculate_pds(
        self,
        perturbation_effects: Dict[str, Dict[str, float]],
        true_perturbations: Dict[str, str]
    ) -> float:
        """
        Calculate Perturbation Discrimination Score (PDS).
        
        PDS measures whether the model assigns the correct effect 
        to the correct perturbation.
        
        Args:
            perturbation_effects: Predicted effects {perturbation_id: {metric: value}}
            true_perturbations: True perturbation assignments {sample_id: perturbation_id}
            
        Returns:
            PDS score (0-1, higher is better)
        """
        if not perturbation_effects or not true_perturbations:
            return 0.0

        correct_assignments = 0
        total_assignments = 0
        
        for sample_id, true_pert in true_perturbations.items():
            if sample_id not in perturbation_effects:
                continue
                
            predicted_effects = perturbation_effects[sample_id]
            
            # Find perturbation with highest predicted effect
            if predicted_effects:
                max_effect_pert = sorted(
                    predicted_effects.items(),
                    key=lambda x: (-x[1], str(x[0])),
                )[0][0]
                
                if max_effect_pert == true_pert:
                    correct_assignments += 1
                
                total_assignments += 1
        
        pds = correct_assignments / total_assignments if total_assignments > 0 else 0.0
        return pds
    
    def calculate_mae(
        self,
        predicted_expression: Dict[str, float],
        actual_expression: Dict[str, float],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> float:
        """Calculate weighted MAE over aligned genes.

        Viral/BK genes use ``viral_genes``, immune-marker names use
        ``immune_genes``, and all other genes use ``housekeeping``. Explicit
        per-gene entries in ``gene_weights`` override category defaults.
        ``inf`` is returned for invalid/non-finite inputs, no shared genes, or
        weights that cannot produce a positive denominator.

        Args:
            predicted_expression: Predicted gene expression by gene.
            actual_expression: Actual gene expression by gene.
            gene_weights: Optional category or per-gene weight overrides.

        Returns:
            MAE score (lower is better).
        """
        if not isinstance(predicted_expression, dict) or not isinstance(actual_expression, dict):
            return float('inf')

        weights = dict(self.gene_weights)
        if gene_weights:
            weights.update(gene_weights)
        
        # Align genes
        common_genes = set(predicted_expression.keys()) & set(actual_expression.keys())
        
        if not common_genes:
            return float('inf')
        
        # Calculate weighted absolute errors
        weighted_errors = []
        total_weight = 0.0
        
        for gene in sorted(common_genes):
            try:
                pred = float(predicted_expression[gene])
                actual = float(actual_expression[gene])
            except (TypeError, ValueError):
                return float('inf')
            if not np.isfinite(pred) or not np.isfinite(actual):
                return float('inf')

            # Determine weight based on gene category
            weight = gene_weights.get(gene) if gene_weights else None
            if weight is None:
                weight = weights.get('housekeeping', 1.0)
                if 'viral' in gene.lower() or 'bk' in gene.lower():
                    weight = weights.get('viral_genes', weight)
                elif any(immune_term in gene.lower() for immune_term in ['cd', 'interferon', 'cytokine', 'immune']):
                    weight = weights.get('immune_genes', weight)
            
            error = abs(pred - actual) * weight
            weighted_errors.append(error)
            total_weight += weight
        
        mae = sum(weighted_errors) / total_weight if total_weight > 0 else float('inf')
        return mae
    
    def evaluate_full(
        self,
        predicted_de: Dict[str, float],
        actual_de: Dict[str, float],
        perturbation_effects: Dict[str, Dict[str, float]],
        true_perturbations: Dict[str, str],
        predicted_expression: Dict[str, float],
        actual_expression: Dict[str, float]
    ) -> EvaluationResults:
        """
        Calculate all ARC evaluation metrics.
        
        Args:
            predicted_de: Predicted differentially expressed genes
            actual_de: Actual differentially expressed genes
            perturbation_effects: Predicted perturbation effects
            true_perturbations: True perturbation assignments
            predicted_expression: Predicted gene expression
            actual_expression: Actual gene expression
            
        Returns:
            EvaluationResults containing all metrics
        """
        des = self.calculate_des(predicted_de, actual_de)
        pds = self.calculate_pds(perturbation_effects, true_perturbations)
        mae = self.calculate_mae(predicted_expression, actual_expression)
        
        details = {
            'des_details': {
                'top_n': len(predicted_de),
                'predicted_genes': len(predicted_de),
                'actual_genes': len(actual_de)
            },
            'pds_details': {
                'total_samples': len(true_perturbations),
                'evaluated_samples': len(perturbation_effects)
            },
            'mae_details': {
                'common_genes': len(set(predicted_expression.keys()) & set(actual_expression.keys())),
                'predicted_genes': len(predicted_expression),
                'actual_genes': len(actual_expression)
            }
        }
        
        return EvaluationResults(
            des=des,
            pds=pds, 
            mae=mae,
            details=details
        )
    
    def calculate_cell_cycle_scores(
        self,
        gene_expression: Dict[str, float],
        s_phase_genes: List[str],
        g2m_phase_genes: List[str]
    ) -> Tuple[float, float, str]:
        """
        Calculate ARC-style cell cycle scores.
        
        Args:
            gene_expression: Gene expression levels {gene: expression}
            s_phase_genes: List of S-phase marker genes
            g2m_phase_genes: List of G2M-phase marker genes
            
        Returns:
            Tuple of (s_score, g2m_score, predicted_phase)
        """
        # Calculate S phase score
        s_score = 0.0
        if s_phase_genes:
            s_gene_expr = [gene_expression.get(gene, 0.0) for gene in s_phase_genes]
            s_score = np.mean(s_gene_expr) if s_gene_expr else 0.0
        
        # Calculate G2M phase score
        g2m_score = 0.0
        if g2m_phase_genes:
            g2m_gene_expr = [gene_expression.get(gene, 0.0) for gene in g2m_phase_genes]
            g2m_score = np.mean(g2m_gene_expr) if g2m_gene_expr else 0.0
        
        # Predict cell cycle phase
        if s_score > g2m_score and s_score > 0.5:
            phase = "S"
        elif g2m_score > s_score and g2m_score > 0.5:
            phase = "G2M"
        else:
            phase = "G1"
        
        return s_score, g2m_score, phase
