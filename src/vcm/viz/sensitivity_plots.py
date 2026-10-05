"""Visualization functions for sensitivity analysis results."""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import List, Optional
from pathlib import Path

from vcm.analysis.sensitivity import SensitivityResult


def plot_sensitivity_scatter(result: SensitivityResult,
                           output_path: Optional[Path] = None,
                           show: bool = False) -> None:
    """Plot parameter vs outcome scatter plot for sensitivity analysis.
    
    Args:
        result: SensitivityResult to visualize
        output_path: Optional path to save the plot
        show: Whether to display the plot
    """
    plt.figure(figsize=(10, 6))
    
    # Scatter plot
    plt.scatter(result.parameter_values, result.outcome_values, 
                alpha=0.6, s=50, color='steelblue', edgecolors='darkblue')
    
    # Add trend line
    z = np.polyfit(result.parameter_values, result.outcome_values, 1)
    p = np.poly1d(z)
    x_line = np.linspace(min(result.parameter_values), max(result.parameter_values), 100)
    plt.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2, label='Trend line')
    
    # Add correlation annotation
    plt.text(0.05, 0.95, f'Correlation: {result.correlation_coefficient:.3f}\nSensitivity Index: {result.sensitivity_index:.3f}',
             transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.xlabel(f'{result.parameter_name} (Parameter Value)', fontsize=12)
    plt.ylabel('Peak Viral Load (copies/mL)', fontsize=12)
    plt.title(f'Sensitivity Analysis: {result.parameter_name}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_tornado_chart(results: List[SensitivityResult],
                      output_path: Optional[Path] = None,
                      show: bool = False) -> None:
    """Plot tornado chart for comparing parameter sensitivities.
    
    Args:
        results: List of SensitivityResult objects
        output_path: Optional path to save the plot
        show: Whether to display the plot
    """
    # Sort results by absolute sensitivity index
    sorted_results = sorted(results, key=lambda x: abs(x.sensitivity_index))
    
    # Extract data
    param_names = [r.parameter_name.replace('_', ' ').title() for r in sorted_results]
    sensitivity_indices = [r.sensitivity_index for r in sorted_results]
    
    # Create colors based on sign
    colors = ['red' if idx < 0 else 'green' for idx in sensitivity_indices]
    
    # Horizontal bar chart
    plt.figure(figsize=(12, 8))
    y_pos = np.arange(len(param_names))
    
    plt.barh(y_pos, sensitivity_indices, color=colors, alpha=0.7, edgecolor='black')
    
    # Add vertical line at zero
    plt.axvline(x=0, color='black', linestyle='-', linewidth=1)
    
    # Add value labels
    for i, (idx, name) in enumerate(zip(sensitivity_indices, param_names)):
        plt.text(idx, i, f'{idx:.3f}', va='center', ha='left' if idx > 0 else 'right',
                fontsize=10, fontweight='bold')
    
    plt.yticks(y_pos, param_names)
    plt.xlabel('Sensitivity Index (Elasticity)', fontsize=12)
    plt.ylabel('Parameters', fontsize=12)
    plt.title('Tornado Chart: Parameter Sensitivity Comparison\n(Higher absolute value = more sensitive)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_confidence_intervals(results: List[SensitivityResult],
                             output_path: Optional[Path] = None,
                             show: bool = False) -> None:
    """Plot confidence intervals for parameter sensitivity results.
    
    Args:
        results: List of SensitivityResult objects
        output_path: Optional path to save the plot
        show: Whether to display the plot
    """
    # Extract data
    param_names = [r.parameter_name.replace('_', ' ').title() for r in results]
    mean_outcomes = [np.mean(r.outcome_values) for r in results]
    ci_lows = [r.confidence_interval_low for r in results]
    ci_highs = [r.confidence_interval_high for r in results]
    
    # Calculate error bars
    yerr_low = np.array(mean_outcomes) - np.array(ci_lows)
    yerr_high = np.array(ci_highs) - np.array(mean_outcomes)
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    plt.errorbar(param_names, mean_outcomes, yerr=[yerr_low, yerr_high],
                 fmt='o', color='steelblue', ecolor='red', elinewidth=2, capsize=5, capthick=2,
                 markersize=8, markeredgecolor='darkblue', markeredgewidth=2)
    
    plt.xlabel('Parameters', fontsize=12)
    plt.ylabel('Peak Viral Load (copies/mL)', fontsize=12)
    plt.title('Parameter Sensitivity: 95% Confidence Intervals\n(Uncertainty in outcomes due to parameter uncertainty)', 
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_correlation_heatmap(results: List[SensitivityResult],
                            output_path: Optional[Path] = None,
                            show: bool = False) -> None:
    """Plot correlation heatmap for parameter sensitivities.
    
    Args:
        results: List of SensitivityResult objects
        output_path: Optional path to save the plot
        show: Whether to display the plot
    """
    # Extract correlation coefficients
    param_names = [r.parameter_name.replace('_', ' ').title() for r in results]
    correlations = [r.correlation_coefficient for r in results]
    
    # Create single-row heatmap
    corr_matrix = np.array(correlations).reshape(1, -1)
    
    plt.figure(figsize=(14, 4))
    
    # Create heatmap
    sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdBu_r',
                center=0, vmin=-1, vmax=1,
                xticklabels=param_names, yticklabels=['Correlation'],
                cbar_kws={'label': 'Correlation Coefficient'},
                linewidths=1, linecolor='black')
    
    plt.title('Parameter-Outcome Correlation Heatmap\n(Correlation between parameter value and peak viral load)', 
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()


def plot_parameter_distribution(result: SensitivityResult,
                               output_path: Optional[Path] = None,
                               show: bool = False) -> None:
    """Plot distribution of parameter values and corresponding outcomes.
    
    Args:
        result: SensitivityResult to visualize
        output_path: Optional path to save the plot
        show: Whether to display the plot
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Parameter distribution
    ax1.hist(result.parameter_values, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.set_xlabel(f'{result.parameter_name} (Parameter Value)', fontsize=11)
    ax1.set_ylabel('Frequency', fontsize=11)
    ax1.set_title(f'Parameter Value Distribution', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Outcome distribution
    ax2.hist(result.outcome_values, bins=20, color='coral', alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Peak Viral Load (copies/mL)', fontsize=11)
    ax2.set_ylabel('Frequency', fontsize=11)
    ax2.set_title('Outcome Distribution', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle(f'Sensitivity Analysis: {result.parameter_name}', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()