#!/usr/bin/env python3
"""
Generate risk prediction figures for BKPyVAN risk prediction.

Creates:
- Figure 1: ROC curves comparison (clinical vs clinical + VCM)
- Figure 2: Feature importance plot for VCM-enhanced model
- Figure 3: Calibration curves for both models
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcm.clinical.risk_prediction import RiskPredictor
from sklearn.metrics import roc_curve, auc
from sklearn.calibration import calibration_curve


def plot_roc_comparison(baseline_results, vcm_results, df, output_path):
    """Generate ROC curve comparison plot.
    
    Args:
        baseline_results: Results from train_clinical_baseline
        vcm_results: Results from train_vcm_enhanced
        df: DataFrame with data
        output_path: Path to save figure
    """
    # Get predictions
    clinical_features = ['age', 'sex', 'prior_transplant', 'diabetes', 'tacrolimus_use',
                        'hla_mismatch', 'donor_age']
    vcm_features = ['peak_viral_load_copies', 'weeks_above_1k', 'weeks_above_10k',
                   'area_under_curve_log', 'time_to_peak_weeks']
    
    X_clinical = df[clinical_features].values
    X_vcm = df[vcm_features].copy()
    X_vcm['peak_viral_load_copies'] = np.log10(df['peak_viral_load_copies'] + 1)
    X_vcm = X_vcm.values
    X_combined = np.concatenate([X_clinical, X_vcm], axis=1)
    
    y = df['bkypan_outcome'].values
    
    # Get probability predictions (in-sample, i.e. predictions on the same
    # data the models were trained on; cross-validated AUCs live in the
    # results JSON, not in this figure)
    baseline_probs = baseline_results['model'].predict_proba(X_clinical)[:, 1]
    vcm_probs = vcm_results['model'].predict_proba(X_combined)[:, 1]

    # Calculate ROC curves
    baseline_fpr, baseline_tpr, _ = roc_curve(y, baseline_probs)
    vcm_fpr, vcm_tpr, _ = roc_curve(y, vcm_probs)

    baseline_auc = auc(baseline_fpr, baseline_tpr)
    vcm_auc = auc(vcm_fpr, vcm_tpr)
    
    # Create plot
    plt.figure(figsize=(8, 8))
    plt.plot(baseline_fpr, baseline_tpr, 'b-', linewidth=2,
             label=f'Clinical only (AUC={baseline_auc:.3f})')
    plt.plot(vcm_fpr, vcm_tpr, 'r-', linewidth=2,
             label=f'Clinical + VCM features (AUC={vcm_auc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random chance')
    
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curves: BKPyVAN Risk Prediction', fontsize=14, fontweight='bold')
    plt.legend(loc='lower right', fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    
    # Add note - be explicit that these curves are in-sample
    plt.text(0.05, 0.95, 'n=500 synthetic patients\nROC shown in-sample (training data);\nsee results JSON for 5-fold CV AUCs',
             transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved ROC comparison plot to {output_path}")


def plot_feature_importance(vcm_results, output_path):
    """Generate feature importance plot for VCM-enhanced model.
    
    Args:
        vcm_results: Results from train_vcm_enhanced
        output_path: Path to save figure
    """
    # Extract feature importances
    importances = vcm_results['feature_importances']
    feature_names = vcm_results['feature_names']
    
    # Sort by absolute coefficient value
    sorted_idx = np.argsort([abs(importances[f]) for f in feature_names])
    sorted_features = [feature_names[i] for i in sorted_idx]
    sorted_values = [importances[f] for f in sorted_features]
    
    # Determine colors (VCM features vs clinical features)
    vcm_feature_names = ['peak_viral_load_copies_log', 'weeks_above_1k', 'weeks_above_10k',
                         'area_under_curve_log', 'time_to_peak_weeks']
    colors = ['red' if f in vcm_feature_names else 'steelblue' for f in sorted_features]
    
    # Create horizontal bar plot
    plt.figure(figsize=(10, 8))
    y_pos = np.arange(len(sorted_features))
    
    bars = plt.barh(y_pos, sorted_values, color=colors, alpha=0.7)
    
    # Add vertical line at x=0
    plt.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
    
    plt.yticks(y_pos, sorted_features, fontsize=11)
    plt.xlabel('Logistic Regression Coefficient', fontsize=12)
    plt.title('Feature Importance: BKPyVAN Risk Prediction (Clinical + VCM)',
              fontsize=14, fontweight='bold')
    plt.grid(axis='x', alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='steelblue', alpha=0.7, label='Clinical features'),
        Patch(facecolor='red', alpha=0.7, label='VCM-derived features')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved feature importance plot to {output_path}")


def plot_calibration_curves(baseline_results, vcm_results, df, output_path):
    """Generate calibration curves for both models.
    
    Args:
        baseline_results: Results from train_clinical_baseline
        vcm_results: Results from train_vcm_enhanced
        df: DataFrame with data
        output_path: Path to save figure
    """
    # Get predictions
    clinical_features = ['age', 'sex', 'prior_transplant', 'diabetes', 'tacrolimus_use',
                        'hla_mismatch', 'donor_age']
    vcm_features = ['peak_viral_load_copies', 'weeks_above_1k', 'weeks_above_10k',
                   'area_under_curve_log', 'time_to_peak_weeks']
    
    X_clinical = df[clinical_features].values
    X_vcm = df[vcm_features].copy()
    X_vcm['peak_viral_load_copies'] = np.log10(df['peak_viral_load_copies'] + 1)
    X_vcm = X_vcm.values
    X_combined = np.concatenate([X_clinical, X_vcm], axis=1)
    
    y = df['bkypan_outcome'].values
    
    # Get probability predictions
    baseline_probs = baseline_results['model'].predict_proba(X_clinical)[:, 1]
    vcm_probs = vcm_results['model'].predict_proba(X_combined)[:, 1]
    
    # Calculate calibration curves (10 bins)
    baseline_prob_true, baseline_prob_pred = calibration_curve(y, baseline_probs, n_bins=10)
    vcm_prob_true, vcm_prob_pred = calibration_curve(y, vcm_probs, n_bins=10)
    
    # Create plot
    plt.figure(figsize=(8, 8))
    
    # Plot calibration curves
    plt.plot(baseline_prob_pred, baseline_prob_true, 'b-o', linewidth=2, markersize=6,
             label=f'Clinical only (Brier={baseline_results["brier_score"]:.3f})')
    plt.plot(vcm_prob_pred, vcm_prob_true, 'r-s', linewidth=2, markersize=6,
             label=f'Clinical + VCM (Brier={vcm_results["brier_score"]:.3f})')
    
    # Plot perfect calibration line
    plt.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Perfect calibration')
    
    plt.xlabel('Mean Predicted Probability', fontsize=12)
    plt.ylabel('Fraction of Positives', fontsize=12)
    plt.title('Calibration Curves: BKPyVAN Risk Prediction',
              fontsize=14, fontweight='bold')
    plt.legend(loc='upper left', fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    
    # Add note
    plt.text(0.05, 0.05, 'n=500 simulated patients, 10 bins',
             transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved calibration plot to {output_path}")


def main():
    """Main function to generate all risk prediction figures."""
    print("Generating risk prediction figures...")
    
    # Create output directory
    output_dir = Path("outputs/clinical")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load synthetic cohort data
    cohort_path = Path("data/processed/synthetic_patient_cohort.csv")
    if not cohort_path.exists():
        print(f"Error: Cohort data not found at {cohort_path}")
        print("Please run scripts/generate_risk_dataset.py first")
        return
    
    df = pd.read_csv(cohort_path)
    print(f"Loaded cohort data: {df.shape}")
    
    # Initialize risk predictor
    predictor = RiskPredictor(random_state=42)
    
    # Train models
    print("Training baseline clinical model...")
    baseline_results = predictor.train_clinical_baseline(df)
    print(f"Baseline AUC: {baseline_results['auc_mean']:.3f} ± {baseline_results['auc_std']:.3f}")
    
    print("Training VCM-enhanced model...")
    vcm_results = predictor.train_vcm_enhanced(df)
    print(f"VCM-enhanced AUC: {vcm_results['auc_mean']:.3f} ± {vcm_results['auc_std']:.3f}")
    
    # Generate figures
    print("\nGenerating figures...")
    
    # Figure 1: ROC comparison
    roc_path = output_dir / "roc_comparison.png"
    plot_roc_comparison(baseline_results, vcm_results, df, roc_path)
    
    # Figure 2: Feature importance
    importance_path = output_dir / "feature_importance.png"
    plot_feature_importance(vcm_results, importance_path)
    
    # Figure 3: Calibration curves
    calibration_path = output_dir / "calibration_plot.png"
    plot_calibration_curves(baseline_results, vcm_results, df, calibration_path)
    
    print("\nAll figures generated successfully!")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()