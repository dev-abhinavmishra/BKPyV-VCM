#!/usr/bin/env python3
"""Generate clinical viral load trajectory figures.

This script runs BKPyV simulations under different immunosuppression scenarios
and generates publication-quality figures of clinical viral load trajectories.
"""

import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import json


def create_trajectories_figure(trajectories):
    """Create the main viral load trajectories figure.
    
    Args:
        trajectories: Dict mapping scenario names to DataFrames
    """
    plt.figure(figsize=(12, 8), dpi=300)
    
    # Colors for each scenario
    colors = {
        'baseline': '#2ecc71',      # Green
        'infection': '#3498db',     # Blue
        'tacrolimus': '#e74c3c',    # Red
        'sirolimus': '#9b59b6'      # Purple
    }
    
    # Plot trajectories
    for scenario, df in trajectories.items():
        plt.plot(df['week'], df['copies_per_ml'], 
                color=colors[scenario], 
                linewidth=2.5, 
                label=scenario.capitalize(),
                alpha=0.8)
    
    # Add clinical threshold lines
    plt.axhline(y=1000, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='Screening threshold (1,000 copies/mL)')
    plt.axhline(y=10000, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Treatment threshold (10,000 copies/mL)')
    
    # Shade high risk zone
    plt.axhspan(10000, plt.ylim()[1], alpha=0.15, color='red', label='High risk zone')
    
    # Annotate tacrolimus peak and sirolimus suppression
    if 'tacrolimus' in trajectories:
        tac_df = trajectories['tacrolimus']
        peak_idx = tac_df['copies_per_ml'].idxmax()
        peak_week = tac_df.loc[peak_idx, 'week']
        peak_load = tac_df.loc[peak_idx, 'copies_per_ml']
        
        plt.annotate('Tacrolimus peak', 
                    xy=(peak_week, peak_load),
                    xytext=(peak_week + 5, peak_load * 0.5),
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                    fontsize=10, color='red')
    
    if 'sirolimus' in trajectories:
        siro_df = trajectories['sirolimus']
        max_siro = siro_df['copies_per_ml'].max()
        max_week = siro_df.loc[siro_df['copies_per_ml'].idxmax(), 'week']
        
        plt.annotate('Sirolimus suppression',
                    xy=(max_week, max_siro),
                    xytext=(max_week - 10, max_siro * 5),
                    arrowprops=dict(arrowstyle='->', color='purple', lw=1.5),
                    fontsize=10, color='purple')
    
    # Format axes
    plt.xlabel('Weeks Post-Transplant', fontsize=12, fontweight='bold')
    plt.ylabel('Plasma BKPyV DNA (copies/mL)', fontsize=12, fontweight='bold')
    plt.yscale('log')
    plt.title('Simulated BKPyV Plasma Viral Load Under Different Immunosuppression Regimens',
              fontsize=14, fontweight='bold', pad=20)
    
    plt.xlim(0, 52)
    plt.ylim(10, 1e8)
    
    # Legend
    plt.legend(loc='upper left', fontsize=10, framealpha=0.9)
    
    plt.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    
    return plt.gcf()


def create_risk_timeline_figure(trajectories):
    """Create risk category timeline heatmap.
    
    Args:
        trajectories: Dict mapping scenario names to DataFrames
    """
    # Define risk category colors
    risk_colors = {
        'undetectable': '#2ecc71',   # Green
        'low_risk': '#3498db',       # Blue
        'screening': '#f39c12',      # Orange
        'treatment': '#e67e22',      # Dark orange
        'severe': '#e74c3c'          # Red
    }
    
    # Create matrix of risk categories
    scenarios = list(trajectories.keys())
    weeks = 53  # 0-52 weeks
    
    risk_matrix = np.zeros((len(scenarios), weeks))
    for i, scenario in enumerate(scenarios):
        df = trajectories[scenario]
        for week in range(weeks):
            if week < len(df):
                risk_cat = df.iloc[week]['risk_category']
                # Map to numeric for color mapping
                if risk_cat == 'undetectable':
                    risk_matrix[i, week] = 0
                elif risk_cat == 'low_risk':
                    risk_matrix[i, week] = 1
                elif risk_cat == 'screening':
                    risk_matrix[i, week] = 2
                elif risk_cat == 'treatment':
                    risk_matrix[i, week] = 3
                elif risk_cat == 'severe':
                    risk_matrix[i, week] = 4
    
    # Create colormap
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([risk_colors['undetectable'], risk_colors['low_risk'], 
                          risk_colors['screening'], risk_colors['treatment'], 
                          risk_colors['severe']])
    
    fig, ax = plt.subplots(figsize=(14, 6), dpi=300)
    
    im = ax.imshow(risk_matrix, aspect='auto', cmap=cmap, vmin=0, vmax=4)
    
    # Set ticks
    ax.set_xticks(np.arange(0, weeks, 4))
    ax.set_xticklabels(np.arange(0, weeks, 4))
    ax.set_yticks(np.arange(len(scenarios)))
    ax.set_yticklabels([s.capitalize() for s in scenarios])
    
    # Labels
    ax.set_xlabel('Weeks Post-Transplant', fontsize=12, fontweight='bold')
    ax.set_ylabel('Immunosuppression Scenario', fontsize=12, fontweight='bold')
    ax.set_title('BKPyV Risk Category Timeline by Scenario',
                fontsize=14, fontweight='bold', pad=20)
    
    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_ticks([0, 1, 2, 3, 4])
    cbar.set_ticklabels(['Undetectable', 'Low Risk', 'Screening', 'Treatment', 'Severe'])
    
    plt.tight_layout()
    
    return plt.gcf()


def main():
    """Main function to generate clinical figures."""
    print("Generating clinical viral load trajectory figures...")
    
    # Add src to path for imports
    src_path = Path(__file__).parent.parent / "src"
    sys.path.insert(0, str(src_path))
    
    from vcm.clinical.viral_load_mapper import ViralLoadMapper, generate_clinical_summary
    
    # Create output directory
    output_dir = Path("outputs/clinical")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize mapper
    mapper = ViralLoadMapper()
    
    # Scenarios to simulate
    scenarios = ['baseline', 'infection', 'tacrolimus', 'sirolimus']
    
    # Simulate trajectories
    trajectories = {}
    for scenario in scenarios:
        print(f"Simulating {scenario} scenario...")
        trajectories[scenario] = mapper.simulate_clinical_trajectory(scenario, weeks=52)
    
    # Generate main trajectories figure
    print("Generating trajectories figure...")
    fig1 = create_trajectories_figure(trajectories)
    fig1.savefig(output_dir / "viral_load_trajectories.png", dpi=300, bbox_inches='tight')
    fig1.savefig(output_dir / "viral_load_trajectories.svg", format='svg', bbox_inches='tight')
    plt.close(fig1)
    print(f"Saved to {output_dir / 'viral_load_trajectories.png'}")
    print(f"Saved to {output_dir / 'viral_load_trajectories.svg'}")
    
    # Generate risk timeline figure
    print("Generating risk timeline figure...")
    fig2 = create_risk_timeline_figure(trajectories)
    fig2.savefig(output_dir / "risk_timeline.png", dpi=300, bbox_inches='tight')
    plt.close(fig2)
    print(f"Saved to {output_dir / 'risk_timeline.png'}")
    
    # Generate clinical summary statistics
    print("Generating clinical summary statistics...")
    summary = generate_clinical_summary(trajectories)
    
    summary_path = output_dir / "clinical_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved to {summary_path}")
    
    # Print summary
    print("\nClinical Summary Statistics:")
    print("=" * 80)
    for scenario, stats in summary.items():
        print(f"\n{scenario.capitalize()}:")
        print(f"  Time to first detectable viremia: {stats['time_to_first_detectable_viremia_weeks']} weeks")
        print(f"  Peak viral load: {stats['peak_viral_load_copies_per_ml']:.0f} copies/mL")
        print(f"  Week of peak: {stats['week_of_peak']}")
        print(f"  Time above screening threshold: {stats['time_above_screening_threshold_weeks']} weeks")
        print(f"  Time above treatment threshold: {stats['time_above_treatment_threshold_weeks']} weeks")
        print(f"  Final viral load (week 52): {stats['final_viral_load_week_52_copies_per_ml']:.0f} copies/mL")
    
    print("\n" + "=" * 80)
    print("Clinical figure generation complete!")


if __name__ == "__main__":
    main()
