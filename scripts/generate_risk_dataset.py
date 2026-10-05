#!/usr/bin/env python3
"""
Generate synthetic training dataset for BKPyVAN risk prediction.

Creates a cohort of 500 simulated patients with:
- Clinical covariates (from Fang et al. 2022 risk factors, PMC9428263)
- VCM simulation features (run actual simulator for each patient)
- Outcome: BKPyVAN within 1 year using logistic model with literature-based ORs
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vcm.clinical.viral_load_mapper import ViralLoadMapper


def generate_clinical_covariates(n_patients=500, random_state=42):
    """Generate clinical covariates based on Fang et al. 2022.
    
    Args:
        n_patients: Number of patients to generate
        random_state: Random seed for reproducibility
        
    Returns:
        DataFrame with clinical covariates
    """
    np.random.seed(random_state)
    
    data = {
        # Age: random, mean=45, std=12
        'age': np.random.normal(45, 12, n_patients),
        
        # Sex: binary (0=female, 1=male), 60% male
        'sex': np.random.binomial(1, 0.60, n_patients),
        
        # Prior transplant: binary, 20% prevalence
        'prior_transplant': np.random.binomial(1, 0.20, n_patients),
        
        # Diabetes: binary, 25% prevalence
        'diabetes': np.random.binomial(1, 0.25, n_patients),
        
        # Tacrolimus use: binary, 75% prevalence (vs sirolimus/other)
        'tacrolimus_use': np.random.binomial(1, 0.75, n_patients),
        
        # HLA mismatch: integer 0-6, mean=3.2
        'hla_mismatch': np.random.normal(3.2, 1.5, n_patients).clip(0, 6).astype(int),
        
        # Donor age: random, mean=42, std=15
        'donor_age': np.random.normal(42, 15, n_patients)
    }
    
    # Ensure age is positive and reasonable
    data['age'] = data['age'].clip(18, 80)
    data['donor_age'] = data['donor_age'].clip(10, 75)
    
    return pd.DataFrame(data)


def run_vcm_simulation_for_patient(covariates, mapper):
    """Run VCM simulation for a single patient and extract features.
    
    Args:
        covariates: Dict with clinical covariates for one patient
        mapper: ViralLoadMapper instance
        
    Returns:
        Dict with VCM features
    """
    # Determine scenario based on tacrolimus_use
    if covariates['tacrolimus_use'] == 1:
        scenario = 'tacrolimus'
    else:
        # Assume sirolimus or other (use sirolimus as proxy)
        scenario = 'sirolimus'
    
    # Run 52-week simulation
    df = mapper.simulate_clinical_trajectory(scenario, weeks=52)
    
    # Extract features
    peak_viral_load_copies = df['copies_per_ml'].max()
    weeks_above_1k = len(df[df['copies_per_ml'] >= 1000])
    weeks_above_10k = len(df[df['copies_per_ml'] >= 10000])
    
    # Area under curve (log scale)
    # Approximate using trapezoidal rule
    log_viral_loads = np.log10(df['copies_per_ml'] + 1)  # +1 to avoid log(0)
    weeks = df['week'].values
    area_under_curve_log = np.trapezoid(log_viral_loads, weeks)
    
    # Time to peak (weeks)
    peak_idx = df['copies_per_ml'].idxmax()
    time_to_peak_weeks = df.loc[peak_idx, 'week']
    
    return {
        'peak_viral_load_copies': peak_viral_load_copies,
        'weeks_above_1k': weeks_above_1k,
        'weeks_above_10k': weeks_above_10k,
        'area_under_curve_log': area_under_curve_log,
        'time_to_peak_weeks': time_to_peak_weeks
    }


def generate_outcomes(df):
    """Generate BKPyVAN outcome from CLINICAL COVARIATES ONLY (+ noise).

    Integrity note: previous versions included a +0.5 * log(peak_viral_load)
    term, i.e. outcomes depended on the VCM features themselves. That made the
    downstream "does the VCM improve prediction?" comparison circular (the
    synthetic truth contained the feature being evaluated). The synthetic
    outcome now depends only on clinical covariates and noise. Any predictive
    value of VCM features then comes only through their correlation with those
    covariates - the honest null setup.

    Direction of risk factors is broadly consistent with published systematic
    reviews (e.g. Demey et al. 2018): tacrolimus, male sex, prior transplant,
    older age. Exact magnitudes are illustrative, not extracted estimates.

    Args:
        df: DataFrame with clinical and VCM features

    Returns:
        Series with outcomes (1=BKPyVAN, 0=no BKPyVAN)
    """
    np.random.seed(42)  # For reproducible noise

    # Logistic model coefficients (log-odds scale); magnitudes illustrative
    intercept = -6.3
    coef_tacrolimus = np.log(2.3)
    coef_prior_transplant = np.log(2.1)
    coef_male = np.log(1.6)
    coef_age = 0.01
    coef_diabetes = np.log(1.3)
    coef_hla_mismatch = 0.1

    # NO VCM-feature term - see docstring (this is the integrity fix).
    X = intercept + \
        coef_age * (df['age'] - 45) / 12 + \
        coef_male * df['sex'] + \
        coef_prior_transplant * df['prior_transplant'] + \
        coef_diabetes * df['diabetes'] + \
        coef_tacrolimus * df['tacrolimus_use'] + \
        coef_hla_mismatch * (df['hla_mismatch'] - 3.2) / 1.5
    
    # Convert to probability
    prob = 1 / (1 + np.exp(-X))
    
    # Add random noise for realism
    noise = np.random.normal(0, 0.1, len(df))
    prob = prob + noise
    prob = np.clip(prob, 0.01, 0.99)
    
    # Generate binary outcome
    outcomes = np.random.binomial(1, prob)
    
    return pd.Series(outcomes, index=df.index)


def main():
    """Main function to generate synthetic patient cohort."""
    print("Generating synthetic patient cohort for BKPyVAN risk prediction...")
    
    # Create output directory
    output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate clinical covariates
    print("Generating clinical covariates (n=500)...")
    df = generate_clinical_covariates(n_patients=500, random_state=42)
    
    # Initialize ViralLoadMapper
    print("Initializing ViralLoadMapper...")
    mapper = ViralLoadMapper()
    
    # Run VCM simulation for each patient
    print("Running VCM simulations for each patient...")
    vcm_features = []
    
    for idx, row in df.iterrows():
        if idx % 50 == 0:
            print(f"  Progress: {idx}/500")
        
        vcm_features.append(run_vcm_simulation_for_patient(row.to_dict(), mapper))
    
    # Add VCM features to DataFrame
    vcm_df = pd.DataFrame(vcm_features)
    df = pd.concat([df, vcm_df], axis=1)
    
    # Generate outcomes
    print("Generating BKPyVAN outcomes using logistic model...")
    df['bkypan_outcome'] = generate_outcomes(df)
    
    # Calculate prevalence
    prevalence = df['bkypan_outcome'].mean()
    print(f"Prevalence: {prevalence:.2%}")
    
    # Save cohort data
    cohort_path = output_dir / "synthetic_patient_cohort.csv"
    df.to_csv(cohort_path, index=False)
    print(f"Saved cohort data to {cohort_path}")
    
    # Save generation parameters (recorded to match the code exactly)
    params = {
        'n_patients': 500,
        'random_state': 42,
        'integrity_note': (
            'Synthetic cohort. Outcomes depend on clinical covariates plus '
            'noise ONLY - VCM features are NOT in the outcome generator, so '
            'this cohort can honestly test whether VCM features add value. '
            'Covariate distributions and ORs are illustrative, broadly '
            'consistent with Demey et al. 2018 (systematic review). '
            'Do not treat coefficients as extracted clinical estimates.'
        ),
        'clinical_covariates': {
            'age': {'mean': 45, 'std': 12},
            'sex': {'male_prevalence': 0.60},
            'prior_transplant': {'prevalence': 0.20},
            'diabetes': {'prevalence': 0.25},
            'tacrolimus_use': {'prevalence': 0.75},
            'hla_mismatch': {'mean': 3.2, 'std': 1.5},
            'donor_age': {'mean': 42, 'std': 15}
        },
        'vcm_simulation': {
            'mapper': 'ViralLoadMapper with piecewise log-linear anchor bridge (assumption-labelled)',
            'scenario_mapping': {
                'tacrolimus_use=1': 'tacrolimus scenario',
                'tacrolimus_use=0': 'sirolimus scenario (proxy)'
            },
            'features_extracted': [
                'peak_viral_load_copies',
                'weeks_above_1k',
                'weeks_above_10k',
                'area_under_curve_log',
                'time_to_peak_weeks'
            ]
        },
        'outcome_model': {
            'type': 'logistic probability + Gaussian noise + Bernoulli draw',
            'intercept': -6.3,
            'coefficients': {
                'age_standardized': 0.01,
                'sex': float(np.log(1.6)),
                'prior_transplant': float(np.log(2.1)),
                'diabetes': float(np.log(1.3)),
                'tacrolimus_use': float(np.log(2.3)),
                'hla_mismatch_standardized': 0.1,
            },
            'vcm_feature_coefficient': 0.0,
            'achieved_prevalence': float(prevalence)
        }
    }
    
    params_path = output_dir / "cohort_generation_params.json"
    with open(params_path, 'w') as f:
        json.dump(params, f, indent=2)
    print(f"Saved generation parameters to {params_path}")
    
    print("\nDataset summary:")
    print(f"  Total patients: {len(df)}")
    print(f"  BKPyVAN cases: {df['bkypan_outcome'].sum()}")
    print(f"  Prevalence: {prevalence:.2%}")
    print(f"  Clinical features: 7")
    print(f"  VCM features: 5")
    print(f"  Total features: 12")
    
    print("\nVCM feature statistics:")
    print(df[['peak_viral_load_copies', 'weeks_above_1k', 'weeks_above_10k', 
              'area_under_curve_log', 'time_to_peak_weeks']].describe())
    
    print("\nSynthetic cohort generation complete!")


if __name__ == "__main__":
    main()
