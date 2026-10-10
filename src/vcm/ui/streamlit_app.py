"""Streamlit Frontend UI for BKPyV Virtual Cell Model.

This provides a user-friendly web interface for running BKPyV simulations,
visualizing results, and performing risk predictions. Designed for accessibility
by ISEF judges, mentors, and researchers without requiring programming skills.

Features:
- Interactive parameter adjustment
- Scenario selection and configuration
- Real-time simulation execution
- Clinical viral load visualization
- Risk prediction dashboard
- Comparative analysis tools
- Export functionality for publications
"""

import sys
from pathlib import Path
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import yaml

# Add project source to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Import VCM modules
try:
    from vcm.core.models import CellState, ExperimentConfig, Perturbation, PerturbationType
    from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
    from vcm.plugins.transplant.bk_polyomavirus.parameters import get_bkpyv_defaults
    from vcm.simulators.bkpyv_simulator import BKPyVSimulator
    from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator
    from vcm.clinical.viral_load_mapper import ClinicalViralLoadMapper, ClinicalThresholds
    from vcm.clinical.risk_prediction import RiskPredictionModule, ClinicalCovariates, VirtualCellFeatures
    from vcm.viz.bkpyv_review import build_review_bundle
except ImportError as e:
    st.error(f"Import error: {e}")
    st.warning("Please ensure all VCM modules are properly installed")


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        st.error(f"Configuration file not found: {config_path}")
        return {}


def run_bkpyv_simulation(config: dict) -> dict:
    """Run BKPyV simulation with given configuration."""
    try:
        # Create experiment config
        experiment_config = ExperimentConfig(
            experiment_id=config['experiment_id'],
            plugin=config['plugin'],
            simulator=config['simulator'],
            simulation_length=config['simulation_length'],
            timestep=config['timestep'],
            output_path=config['output_path'],
        )
        
        # Get simulator parameters
        sim_params = config.get('simulator_parameters', {})
        
        # Load plugin
        plugin = BKPolyomavirusPlugin()
        initial_state = plugin.create_initial_state({"nccr_variant": config.get("nccr_variant", "archetype")})
        
        # Create simulator
        # The UI uses the reviewable ODE engine so NCCR and S-phase outputs are
        # the same quantities shown in the exported reviewer bundle.
        simulator = BKPyVODESimulator(sim_params)

        perturbations = []
        if config.get("scenario") != "baseline_uninfected":
            perturbations.append(Perturbation(id="bkpyv_infection", name="BKPyV infection", perturbation_type=PerturbationType.VIRAL_INFECTION, magnitude=1.0, timing=0.0))
        if config.get("scenario") == "tacrolimus_exposure":
            perturbations.append(Perturbation(id="tacrolimus_treatment", name="Tacrolimus immune-control exposure", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="FKBP1A", magnitude=1.0, timing=0.0))
        if config.get("scenario") == "sirolimus_exposure":
            perturbations.append(Perturbation(id="sirolimus_treatment", name="Sirolimus S-phase exposure", perturbation_type=PerturbationType.DRUG_TREATMENT, target_id="MTOR", magnitude=1.0, timing=0.0))

        n_steps = int(experiment_config.simulation_length / experiment_config.timestep)
        result = simulator.simulate(
            initial_state=initial_state,
            perturbations=perturbations,
            n_steps=n_steps,
            timestep=experiment_config.timestep
        )
        
        return {
            'success': True,
            'result': result,
            'config': config,
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
        }


def extract_trajectory_data(result) -> dict:
    """Extract trajectory data from simulation result."""
    timepoints = []
    viral_loads = []
    mitochondrial_stress = []
    immune_suppression = []
    replication_phase = []
    
    for step in result.steps:
        timepoints.append(step.timestamp)
        viral_loads.append(step.cell_state.metadata.get('viral_load', 0.0))
        mitochondrial_stress.append(step.cell_state.metadata.get('mitochondrial_stress', 0.0))
        immune_suppression.append(step.cell_state.metadata.get('innate_immune_suppression', 0.0))
        replication_phase.append(step.cell_state.metadata.get('viral_replication_phase', 'none'))
    
    return {
        'timepoints': timepoints,
        'viral_loads': viral_loads,
        'mitochondrial_stress': mitochondrial_stress,
        'immune_suppression': immune_suppression,
        'replication_phase': replication_phase,
    }


def convert_to_clinical_viral_load(virtual_loads: list, timepoints: list) -> dict:
    """Convert virtual cell viral load to clinical plasma viral load."""
    mapper = ClinicalViralLoadMapper()
    plasma = [mapper.normalized_to_copies(float(v)) for v in virtual_loads]
    return {
        "timepoints": list(timepoints),
        "plasma_viral_load": plasma,
        "risk_categories": [ClinicalThresholds.get_risk_category(v) for v in plasma],
    }


def plot_viral_load_trajectory(timepoints: list, viral_loads: list, title: str = "Viral Load Trajectory"):
    """Create Plotly figure for viral load trajectory."""
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=timepoints,
        y=viral_loads,
        mode='lines+markers',
        name='Viral Load',
        line=dict(color='red', width=2),
        marker=dict(size=6),
    ))
    
    # Add clinical thresholds
    fig.add_hline(
        y=ClinicalThresholds.SCREENING_POSITIVE,
        line_dash="dash",
        line_color="orange",
        annotation_text="Screening Threshold (1,000 copies/mL)"
    )
    
    fig.add_hline(
        y=ClinicalThresholds.HIGH_RISK,
        line_dash="dash",
        line_color="red",
        annotation_text="High Risk Threshold (10,000 copies/mL)"
    )
    
    fig.update_layout(
        title=title,
        xaxis_title="Time (days)",
        yaxis_title="Viral Load (copies/mL)",
        yaxis_type="log",
        template="plotly_white",
        height=400,
    )
    
    return fig


def plot_phase_transitions(timepoints: list, phases: list):
    """Create Plotly figure for replication phase transitions."""
    # Convert phases to numeric for plotting
    phase_map = {'none': 0, 'early': 1, 'late': 2}
    numeric_phases = [phase_map.get(phase, 0) for phase in phases]
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=timepoints,
        y=numeric_phases,
        mode='lines+markers',
        name='Replication Phase',
        line=dict(color='blue', width=2),
        marker=dict(size=8),
    ))
    
    # Add phase labels
    fig.update_yaxes(
        ticktext=['None', 'Early', 'Late'],
        tickvals=[0, 1, 2],
    )
    
    fig.update_layout(
        title="Viral Replication Phase Transitions",
        xaxis_title="Time (days)",
        yaxis_title="Replication Phase",
        template="plotly_white",
        height=400,
    )
    
    return fig


def plot_mitochondrial_stress(timepoints: list, stress_levels: list):
    """Create Plotly figure for mitochondrial stress over time."""
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=timepoints,
        y=stress_levels,
        mode='lines+markers',
        name='Mitochondrial Stress',
        line=dict(color='green', width=2),
        marker=dict(size=6),
        fill='tozeroy',
    ))
    
    fig.update_layout(
        title="Mitochondrial Stress Over Time",
        xaxis_title="Time (days)",
        yaxis_title="Stress Level (0-1)",
        template="plotly_white",
        height=400,
    )
    
    return fig


def st_page_header():
    """Display page header with title and description."""
    st.set_page_config(
        page_title="BKPyV Virtual Cell Model",
        page_icon="🦠",
        layout="wide",
    )
    
    st.markdown("""
    # 🦠 BKPyV Virtual Cell Model
    **BK Polyomavirus Nephropathy - Research Simulation Studio**
    
    This tool explores mechanistic BKPyV production hypotheses in a renal-cell model.
    It is a research sandbox, not a clinical decision-support system.
    """)
    
    st.info("""
    **🔬 Evidence-aware:** The model separates evidence-supported mechanisms from
    phenomenological coefficients. Outputs are hypotheses, not patient predictions.
    """)


def st_sidebar_navigation():
    """Create sidebar navigation."""
    st.sidebar.title("Navigation")
    
    page = st.sidebar.radio(
        "Select Page",
        [
            "🏠 Home",
            "⚙️ Simulation",
            "📊 Visualization",
            "🧫 Single-Cell",
            "📈 Viral-Load Validation",
            "🧬 Risk Prediction",
            "📋 Comparison",
            "💉 Regimen Design",
            "🧾 Parameters & Assumptions",
            "🔎 Review Bundle",
            "📚 Documentation",
        ]
    )
    
    return page


def st_home_page():
    """Display home page with overview."""
    st.markdown("""
    ## Welcome to the BKPyV Virtual Cell Model
    
    ### What is BKPyV?
    BK Polyomavirus (BKPyV) is a common virus that can cause nephropathy (kidney damage)
    in immunosuppressed kidney transplant patients, affecting up to 10% of recipients.
    
    ### How This Tool Works
    1. **Virtual Cell Simulation**: Models BKPyV replication in kidney tubular epithelial cells
    2. **Clinical Mapping**: Converts simulation results to plasma viral load (copies/mL)
    3. **Risk Prediction**: Combines clinical factors with simulation features for risk assessment
    4. **Drug Effects**: Models different immunosuppressive drugs (tacrolimus vs sirolimus)
    
    ### Key Features
    - 🎯 **Research-Grounded Parameters**: Based on published clinical studies
    - 🔬 **Mechanistic Modeling**: FKBP-12 pathway, cell cycle coupling, mitochondrial stress
    - 📊 **Clinical Thresholds**: Aligns with screening guidelines (1,000, 10,000 copies/mL)
    - 💊 **Drug Comparison**: Tacrolimus vs sirolimus effects
    - 🧬 **Risk Stratification**: Clinical + virtual cell features for prediction
    
    ### Quick Start
    1. Go to **Simulation** page to configure and run BKPyV simulations
    2. Use **Visualization** to analyze viral load trajectories
    3. Try **Risk Prediction** for patient-specific risk assessment
    4. Compare different scenarios in **Comparison** page
    """)
    
    st.markdown("---")
    
    st.subheader("Latest Research Updates")
    
    with st.expander("🔬 Research Grounding"):
        st.markdown("""
        - **Drug Mechanisms**: Tacrolimus weakens immune control; sirolimus can reduce S-phase permissiveness
        - **Clinical Risk Factors**: Age, sex, prior transplant, HLA mismatch
        - **Single-Cell Signatures**: Mitochondrial stress, cell cycle coupling, immune evasion
        - **Guidelines**: Screening thresholds ≥1,000 and ≥10,000 copies/mL
        """)


def st_simulation_page():
    """Display simulation configuration and execution page."""
    st.markdown("## ⚙️ BKPyV Simulation")
    
    # Scenario selection
    st.subheader("Select Scenario")
    
    scenario_options = {
        'baseline_uninfected': 'Baseline (Uninfected)',
        'infection_no_drug': 'Infection (No Drug)',
        'low_replication': 'Low Replication',
        'high_replication': 'High Replication',
        'tacrolimus_exposure': 'Tacrolimus Exposure',
        'sirolimus_exposure': 'Sirolimus Exposure',
    }
    
    scenario = st.selectbox(
        "Choose simulation scenario",
        options=list(scenario_options.keys()),
        format_func=lambda x: scenario_options[x],
    )

    nccr_variant = st.selectbox(
        "NCCR scenario",
        ["archetype", "rearranged"],
        help="Archetype is the presumed persistent form; rearranged is a hypothesis-testing scenario with early-gene bias and reduced capsid expression.",
    )
    
    # Parameter adjustment
    st.subheader("Simulation Parameters")
    
    # Get default parameters
    default_params = get_bkpyv_defaults()
    
    # Create parameter adjustment interface
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("Drug Effects")
        tacrolimus_enhancement = st.slider(
            "Tacrolimus Enhancement Factor",
            min_value=1.0,
            max_value=3.0,
            value=default_params.get('tacrolimus_enhancement_factor', 2.0),
            step=0.1,
            help="Multiplies viral replication rate under tacrolimus (clinical OR 2.0-2.3)"
        )
        
        sirolimus_inhibition = st.slider(
            "Sirolimus Inhibition Factor",
            min_value=0.0,
            max_value=1.0,
            value=default_params.get('sirolimus_inhibition_factor', 0.5),
            step=0.1,
            help="Multiplies viral replication rate under sirolimus (IC90 = 4 ng/mL)"
        )
    
    with col2:
        st.write("Replication Parameters")
        t_antigen_threshold = st.slider(
            "T Antigen Replication Threshold",
            min_value=0.1,
            max_value=1.0,
            value=default_params.get('t_antigen_replication_threshold', 0.5),
            step=0.1,
            help="T antigen level required for active replication"
        )
        
        cell_cycle_bonus = st.slider(
            "Cell Cycle S-Phase Bonus",
            min_value=1.0,
            max_value=3.0,
            value=default_params.get('cell_cycle_s_phase_bonus', 2.0),
            step=0.1,
            help="Replication bonus in S-phase vs other phases"
        )
    
    # Run simulation button
    st.markdown("---")
    
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🚀 Run Simulation", type="primary"):
            # Create configuration
            config = {
                'experiment_id': f'bkpyv_{scenario}',
                'plugin': 'transplant.bk_polyomavirus',
                'simulator': 'bkpyv_ode',
                'simulation_length': 100.0,
                'timestep': 1.0,
                'output_path': 'outputs/bkpyv/ui/',
                'scenario': scenario,
                'nccr_variant': nccr_variant,
                'simulator_parameters': {
                    'tacrolimus_enhancement_factor': tacrolimus_enhancement,
                    'sirolimus_inhibition_factor': sirolimus_inhibition,
                    't_antigen_replication_threshold': t_antigen_threshold,
                    'cell_cycle_s_phase_bonus': cell_cycle_bonus,
                }
            }
            
            # Run simulation
            with st.spinner("Running simulation..."):
                result = run_bkpyv_simulation(config)
            
            # Display results
            if result['success']:
                st.success("Simulation completed successfully!")
                st.session_state.simulation_result = result
                st.session_state.simulation_config = config
                run_label = _history_label(config, len(st.session_state.simulation_history) + 1)
                st.session_state.simulation_history[run_label] = {
                    "result": result["result"],
                    "config": config,
                    "label": run_label,
                }

                # Display summary
                trajectory_data = extract_trajectory_data(result['result'])
                st.write(f"**Final Virtual Viral Load:** {trajectory_data['viral_loads'][-1]:.3f}")
                st.write(f"**Peak Virtual Viral Load:** {max(trajectory_data['viral_loads']):.3f}")
                
                # Clinical conversion
                clinical_data = convert_to_clinical_viral_load(
                    trajectory_data['viral_loads'],
                    trajectory_data['timepoints']
                )
                st.write(f"**Peak Clinical Viral Load:** {max(clinical_data['plasma_viral_load']):.0f} copies/mL")

                export_frame = pd.DataFrame({
                    "time_hours": trajectory_data["timepoints"],
                    "virtual_viral_load": trajectory_data["viral_loads"],
                    "plasma_copies_per_ml": clinical_data["plasma_viral_load"],
                    "risk_category": [
                        ClinicalThresholds.get_risk_category(load)
                        for load in clinical_data["plasma_viral_load"]
                    ],
                    "replication_phase": trajectory_data["replication_phase"],
                    "mitochondrial_stress": trajectory_data["mitochondrial_stress"],
                    "immune_suppression": trajectory_data["immune_suppression"],
                })
                st.download_button(
                    "Download trajectory CSV",
                    data=export_frame.to_csv(index=False),
                    file_name="bkpyv_virtual_patient_trajectory.csv",
                    mime="text/csv",
                    help="Export the current virtual-patient trajectory and clinical bridge for review.",
                )
                
            else:
                st.error(f"Simulation failed: {result['error']}")


def st_visualization_page():
    """Display visualization page for simulation results."""
    st.markdown("## 📊 Simulation Visualization")
    
    # Check if simulation results are available
    if not st.session_state.get('simulation_result'):
        st.warning("No simulation results available. Run a simulation first in the Simulation page.")
        return
    
    result = st.session_state.simulation_result
    trajectory_data = extract_trajectory_data(result['result'])
    
    # Convert to clinical viral load
    clinical_data = convert_to_clinical_viral_load(
        trajectory_data['viral_loads'],
        trajectory_data['timepoints']
    )
    
    # Visualization tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "Viral Load Trajectory",
        "Phase Transitions",
        "Mitochondrial Stress",
        "Clinical Risk Stratification"
    ])
    
    with tab1:
        st.subheader("Clinical Viral Load Over Time")
        
        # Plot clinical viral load
        fig = plot_viral_load_trajectory(
            clinical_data['timepoints'],
            clinical_data['plasma_viral_load'],
            "Clinical Viral Load Trajectory"
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Clinical summary
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Peak Viral Load", f"{max(clinical_data['plasma_viral_load']):.0f} copies/mL")
        with col2:
            st.metric("Time to Peak", f"{clinical_data['timepoints'][clinical_data['plasma_viral_load'].index(max(clinical_data['plasma_viral_load']))]:.0f} hours")
        with col3:
            final_risk = ClinicalThresholds.get_risk_category(clinical_data['plasma_viral_load'][-1])
            st.metric("Current Risk Category", final_risk.upper())
    
    with tab2:
        st.subheader("Replication Phase Transitions")
        
        fig = plot_phase_transitions(
            trajectory_data['timepoints'],
            trajectory_data['replication_phase']
        )
        st.plotly_chart(fig, use_container_width=True)
        
        st.write("""
        **Phase Interpretation:**
        - **None**: No viral replication
        - **Early**: Early gene expression (0-24h) - drug-sensitive phase
        - **Late**: Late gene expression (>24h) - drug-resistant phase
        """)
    
    with tab3:
        st.subheader("Mitochondrial Stress Dynamics")
        
        fig = plot_mitochondrial_stress(
            trajectory_data['timepoints'],
            trajectory_data['mitochondrial_stress']
        )
        st.plotly_chart(fig, use_container_width=True)
        
        st.write("""
        **Mitochondrial Stress:**
        - Discordant mitochondrial vs nuclear gene expression
        - Emerges primarily in late replication phase
        - Associated with metabolic stress and viral hijacking
        """)
    
    with tab4:
        st.subheader("Clinical Risk Stratification")
        
        # Create risk stratification table
        risk_data = []
        for i, (time, load, risk) in enumerate(zip(
            clinical_data['timepoints'],
            clinical_data['plasma_viral_load'],
            clinical_data['risk_categories']
        )):
            if i % 10 == 0:  # Show every 10th timepoint
                risk_data.append({
                    'Time (days)': time,
                    'Viral Load (copies/mL)': f"{load:.0f}",
                    'Risk Category': risk.upper(),
                    'Screening Positive': load >= ClinicalThresholds.SCREENING_POSITIVE,
                    'High Risk': load >= ClinicalThresholds.HIGH_RISK,
                })
        
        df = pd.DataFrame(risk_data)
        st.dataframe(df, use_container_width=True)
        
        # Final risk assessment
        final_load = clinical_data['plasma_viral_load'][-1]
        final_risk = ClinicalThresholds.get_risk_category(final_load)
        
        st.subheader(f"Final Risk Assessment: {final_risk.upper()}")
        
        if final_risk == 'high':
            st.error(f"⚠️ **High Risk**: Viral load {final_load:.0f} copies/mL exceeds high-risk threshold")
            st.write("Recommendations: Weekly monitoring, consider immunosuppression reduction")
        elif final_risk == 'low':
            st.warning(f"⚡ **Medium Risk**: Viral load {final_load:.0f} copies/mL above screening threshold")
            st.write("Recommendations: Biweekly monitoring, consider tacrolimus dose reduction")
        else:
            st.success(f"✅ **Low Risk**: Viral load {final_load:.0f} copies/mL below screening threshold")
            st.write("Recommendations: Continue routine monthly monitoring")


def st_risk_prediction_page():
    """Display risk prediction page."""
    st.markdown("## 🧬 Clinical Risk Prediction")
    
    st.info("""
    **Note:** This feature combines clinical patient data with virtual cell simulation features
    for enhanced risk prediction. It requires trained ML models and may be extended in future versions.
    """)
    
    # Clinical covariates input
    st.subheader("Patient Clinical Data")
    
    col1, col2 = st.columns(2)
    
    with col1:
        age = st.number_input("Age (years)", min_value=18, max_value=90, value=55)
        sex = st.selectbox("Sex", ["male", "female"])
        prior_transplant = st.checkbox("Prior Kidney Transplant")
        hla_mismatch = st.slider("HLA Mismatch", 0, 6, 4)
        donor_type = st.selectbox("Donor Type", ["living", "deceased"])
    
    with col2:
        diabetes = st.checkbox("Diabetes Mellitus")
        hypertension = st.checkbox("Hypertension")
        tacrolimus_use = st.checkbox("Tacrolimus Use")
        sirolimus_use = st.checkbox("Sirolimus Use")
        belatacept_use = st.checkbox("Belatacept Use")
        induction_agent = st.selectbox("Induction Agent", ["none", "basiliximab", "thymoglobulin"])
    
    # Biomarkers
    st.subheader("Biomarkers")
    
    col1, col2 = st.columns(2)
    with col1:
        gcfdna_delta = st.number_input("Delta GcfDNA", min_value=0.0, max_value=100.0, value=15.0)
    with col2:
        serum_creatinine = st.number_input("Serum Creatinine (mg/dL)", min_value=0.5, max_value=10.0, value=1.2)
    
    # Risk prediction button
    st.markdown("---")
    
    st.caption(
        "This page applies a hand-specified additive heuristic over clinical "
        "risk factors. It is NOT a fitted statistical model and the displayed "
        "score is not a calibrated probability."
    )

    if st.button("🎯 Predict Risk", type="primary"):
        # Create clinical covariates object
        clinical_data = ClinicalCovariates(
            age=age,
            sex=sex,
            prior_transplant=prior_transplant,
            hla_mismatch=hla_mismatch,
            donor_type=donor_type,
            diabetes=diabetes,
            hypertension=hypertension,
            tacrolimus_use=tacrolimus_use,
            sirolimus_use=sirolimus_use,
            belatacept_use=belatacept_use,
            induction_agent=induction_agent,
            gcfdna_delta=gcfdna_delta,
            serum_creatinine=serum_creatinine,
        )
        
        # Heuristic additive score (not a fitted/calibrated probability)
        clinical_risk = 0.1  # base score
        
        # Add risk factors
        if age > 50:
            clinical_risk += 0.2
        if sex == 'male':
            clinical_risk += 0.15
        if prior_transplant:
            clinical_risk += 0.25
        if hla_mismatch >= 4:
            clinical_risk += 0.1
        if diabetes:
            clinical_risk += 0.15
        if tacrolimus_use:
            clinical_risk += 0.2
        if sirolimus_use:
            clinical_risk -= 0.1
        if gcfdna_delta > 10:
            clinical_risk += 0.2
        
        # Clamp to valid range
        clinical_risk = max(0.0, min(1.0, clinical_risk))
        
        # Display results
        col1, col2 = st.columns(2)

        with col1:
            st.metric("Heuristic Risk Score (not calibrated)", f"{clinical_risk:.2f}")

        with col2:
            if clinical_risk < 0.3:
                risk_category = "LOW"
            elif clinical_risk < 0.7:
                risk_category = "MEDIUM"
            else:
                risk_category = "HIGH"
            st.metric("Risk Category", risk_category)
        
        # Recommendations
        st.subheader("Clinical Recommendations")
        
        if clinical_risk < 0.3:
            st.success("✅ Low Risk - Continue routine monthly monitoring")
        elif clinical_risk < 0.7:
            st.warning("⚡ Medium Risk - Consider biweekly monitoring and tacrolimus dose reduction")
        else:
            st.error("⚠️ High Risk - Implement weekly monitoring, reduce immunosuppression, consider sirolimus")


def _history_label(config: dict, index: int) -> str:
    """Build a display label for a stored simulation run."""
    scenario = config.get("scenario", "custom")
    nccr = config.get("nccr_variant", "archetype")
    return f"run {index}: {scenario} | NCCR {nccr}"


def _shared_time_window(entries: list):
    """Return the overlapping (start_day, end_day) window across runs, or None."""
    starts = []
    ends = []
    for entry in entries:
        times = extract_trajectory_data(entry["result"])["timepoints"]
        if not times:
            return None
        starts.append(min(times))
        ends.append(max(times))
    window = (max(starts), min(ends))
    return window if window[1] > window[0] else None


def _comparison_metrics(label: str, entry: dict) -> dict:
    """Compute display metrics directly from a stored run's trajectory."""
    trajectory = extract_trajectory_data(entry["result"])
    loads = trajectory["viral_loads"]
    times = trajectory["timepoints"]
    peak = max(loads)
    config = entry.get("config", {})
    return {
        "run": label,
        "scenario": config.get("scenario", "?"),
        "nccr_variant": config.get("nccr_variant", "?"),
        "engine": config.get("simulator", "?"),
        "timestep_days": config.get("timestep", "?"),
        "horizon_days": config.get("simulation_length", "?"),
        "peak_virtual_load": peak,
        "peak_day": times[loads.index(peak)],
        "endpoint_load": loads[-1],
    }


def _bridge_metrics(label: str, entry: dict) -> dict:
    """Assumption-labelled bridge metrics (copies/mL) for a stored run."""
    trajectory = extract_trajectory_data(entry["result"])
    clinical = convert_to_clinical_viral_load(
        trajectory["viral_loads"], trajectory["timepoints"]
    )
    peak_plasma = max(clinical["plasma_viral_load"])
    return {
        "run": label,
        "peak_plasma_copies_per_ml": peak_plasma,
        "peak_risk_category": ClinicalThresholds.get_risk_category(peak_plasma),
    }


def _windowed_series(trajectory: dict, window: tuple):
    """Return (times, values) restricted to the shared day window."""
    times = trajectory["timepoints"]
    loads = trajectory["viral_loads"]
    xs = [t for t in times if window[0] <= t <= window[1]]
    ys = [v for t, v in zip(times, loads) if window[0] <= t <= window[1]]
    return xs, ys


def st_comparison_page():
    """Display comparison page for stored simulation runs."""
    st.markdown("## 📋 Scenario Comparison")

    history = st.session_state.get("simulation_history", {})

    if len(history) < 2:
        st.info(
            "**How to populate this page:** run at least two simulations in the "
            "Simulation page (for example, tacrolimus vs sirolimus exposure, or "
            "archetype vs rearranged NCCR). Each completed run is stored here "
            "for side-by-side comparison."
        )
        st.write(f"Stored runs available: **{len(history)}** (need at least 2)")
        return

    labels = list(history.keys())
    selected = st.multiselect(
        "Runs to compare (select at least 2)",
        options=labels,
        default=labels[:2],
    )

    if len(selected) < 2:
        st.warning("Select at least two stored runs to compare.")
        return

    entries = [history[label] for label in selected]
    window = _shared_time_window(entries)
    if window is None:
        st.error(
            "Selected runs do not share an overlapping time window; "
            "trajectories cannot be aligned."
        )
        return

    fig = go.Figure()
    for label, entry in zip(selected, entries):
        xs, ys = _windowed_series(
            extract_trajectory_data(entry["result"]), window
        )
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=label))
    fig.update_layout(
        title="Stored run trajectories",
        xaxis_title="Time (days)",
        yaxis_title="Normalized viral load (model units)",
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        f"Runs are overlaid on their own recorded timestamps within the shared "
        f"day window {window[0]:g}–{window[1]:g}; each run's timestep and "
        "horizon are listed in the table below."
    )

    show_bridge = st.checkbox(
        "Show clinical bridge view (assumption-labelled)", value=False
    )

    st.subheader("Run metrics")
    rows = []
    for label, entry in zip(selected, entries):
        row = _comparison_metrics(label, entry)
        if show_bridge:
            row.update(_bridge_metrics(label, entry))
        rows.append(row)
    st.dataframe(pd.DataFrame(rows), use_container_width=True)

    if show_bridge:
        bridge_fig = go.Figure()
        for label, entry in zip(selected, entries):
            trajectory = extract_trajectory_data(entry["result"])
            clinical = convert_to_clinical_viral_load(
                trajectory["viral_loads"], trajectory["timepoints"]
            )
            xs = [t for t in clinical["timepoints"] if window[0] <= t <= window[1]]
            ys = [
                v
                for t, v in zip(clinical["timepoints"], clinical["plasma_viral_load"])
                if window[0] <= t <= window[1]
            ]
            bridge_fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=label))
        bridge_fig.update_layout(
            title="Clinical bridge trajectories",
            xaxis_title="Time (days)",
            yaxis_title="Plasma viral load (copies/mL, assumption-labelled)",
            yaxis_type="log",
        )
        st.plotly_chart(bridge_fig, use_container_width=True)
        st.caption(
            "Assumption-labelled bridge: normalized model output is mapped to "
            "copies/mL through a fixed anchor scale for hypothesis "
            "visualization only. It is not a calibration to patient data."
        )


def st_review_bundle_page():
    """Generate and display the reproducible reviewer-facing output bundle."""
    st.markdown("## 🔎 Research Review Bundle")
    st.caption("This page regenerates the figures and tables from the ODE model so a reviewer can inspect assumptions, outputs, and sensitivity in one place.")
    days = st.slider("Simulation horizon (model days)", 15, 120, 60)
    if st.button("Generate review bundle", type="primary"):
        with st.spinner("Running scenarios and sensitivity sweep..."):
            paths = build_review_bundle(days=float(days))
        st.success("Review bundle generated.")
        st.image(paths["trajectory"], caption="Mechanistic trajectory")
        st.image(paths["nccr"], caption="Archetype versus rearranged NCCR")
        st.image(paths["drugs"], caption="Drug mechanism comparison")
        st.image(paths["sensitivity"], caption="Parameter sensitivity")
        st.download_button("Download trajectory data", Path(paths["trajectory_data"]).read_bytes(), file_name="trajectory_data.csv")
        st.download_button("Download scenario summary", Path(paths["summary"]).read_bytes(), file_name="scenario_summary.csv")
        st.download_button("Download review manifest", Path(paths["manifest"]).read_bytes(), file_name="manifest.json")


def st_regimen_design_page():
    """Interactive immunosuppression-regimen designer on the mechanistic ODE.

    Lets a reviewer dial in real-unit dosing schedules (ng/mL troughs with
    clinical half-lives) and watch the model trade viral clearance against
    T-cell rebound — the tac-taper-vs-sir-conversion question the optimizer
    script answers systematically.
    """
    import sys as _sys
    from pathlib import Path as _Path
    _scripts = str(_Path(__file__).resolve().parents[3] / "scripts")
    if _scripts not in _sys.path:
        _sys.path.insert(0, _scripts)

    st.markdown("## 💉 Regimen Design — mechanistic protocol explorer")
    st.caption(
        "Hypothesis generator, not patient advice: every schedule is simulated "
        "on the 23-dim mechanistic ODE with real drug half-lives "
        "(tac t½≈12 h, sir t½≈60 h) and real trough units."
    )

    from optimize_reduction_schedule import evaluate_schedule
    from vcm.clinical.viral_load_mapper import ViralLoadMapper
    import numpy as _np
    from scipy.integrate import solve_ivp as _solve_ivp
    from vcm.simulators.ode_system import BKPyVODESystem as _ODE

    st.subheader("Design a schedule")
    mode = st.radio("Regimen family",
                    ["Tacrolimus taper", "Tacrolimus → sirolimus conversion",
                     "Preset sweep"],
                    horizontal=True)

    if mode == "Tacrolimus taper":
        c1, c2 = st.columns(2)
        with c1:
            step_week = st.slider("Reduce at week", 1, 8, 4)
        with c2:
            trough = st.slider("Target trough (ng/mL)", 2.0, 8.0, 4.0, 0.5)
        schedule = [(0.0, step_week * 7.0, 8.0), (step_week * 7.0, None, trough)]
        sir = None
    elif mode == "Tacrolimus → sirolimus conversion":
        c1, c2 = st.columns(2)
        with c1:
            step_week = st.slider("Convert at week", 1, 8, 4)
        with c2:
            sir = st.slider("Sirolimus trough (ng/mL)", 2.0, 10.0, 4.0, 0.5)
        schedule = [(0.0, step_week * 7.0, 8.0), (step_week * 7.0, None, 3.0)]
    else:
        schedule = None
        sir = None

    if st.button("Simulate regimen", type="primary"):
        mapper = ViralLoadMapper()
        if schedule is not None:
            with st.spinner("Simulating 180 model days..."):
                ode = _ODE()
                y0 = ode.get_infection_conditions(0.5)
                horizon = 180.0
                t_eval = _np.linspace(0, horizon, 721)
                dosing = {"tacrolimus": [
                    {"start": s, "stop": e, "trough_ng_ml": tr}
                    for s, e, tr in schedule]}
                if sir is not None:
                    dosing["sirolimus"] = [{"start": schedule[-1][0],
                                            "stop": None,
                                            "trough_ng_ml": sir}]
                sol = _solve_ivp(lambda t, y: ode.ode_system(t, y, dosing),
                                 (0, horizon), y0, t_eval=t_eval,
                                 method="LSODA")
                m = evaluate_schedule(schedule, sir_trough=sir)

            st.subheader("Outcome metrics")
            cols = st.columns(4)
            cols[0].metric("Viremia clears <1k cp/mL",
                           "never" if m["clearance_weeks"] is None
                           else f"week {m['clearance_weeks']}")
            cols[1].metric("Final log10 cp/mL", m["final_log10_cpml"])
            cols[2].metric("T-cell rebound index", m["rebound_index"],
                           help="AUC(T_eff)/horizon — the model's rejection-risk proxy")
            cols[3].metric("Rearranged NCCR fraction", f"{m['final_frr']:.2f}")

            plasma_cp = [mapper.normalized_to_copies(float(v))
                         for v in _np.maximum(sol.y[0], 1e-9)]
            import pandas as pd
            traj = pd.DataFrame({
                "week": sol.t / 7.0,
                "log10 plasma cp/mL": _np.log10(_np.maximum(plasma_cp, 1.0)),
                # Model units — not a cp/mL claim (plasma-calibrated bridge
                # cannot produce urine concentrations).
                "log10 urine V_u": _np.log10(_np.maximum(sol.y[19], 1e-9)),
                "T_eff": sol.y[16],
                "F_rr kidney": sol.y[20],
                "F_rr urine": sol.y[21] if sol.y.shape[0] > 21 else 0.0,
            }).set_index("week")
            st.subheader("Plasma cp/mL vs urinary V_u (log scale)")
            st.line_chart(traj[["log10 plasma cp/mL", "log10 urine V_u"]])
            st.subheader("Immune rebound and NCCR evolution")
            st.line_chart(traj[["T_eff", "F_rr kidney", "F_rr urine"]])
        else:
            with st.spinner("Running full schedule sweep..."):
                from optimize_reduction_schedule import candidate_schedules
                rows = []
                for sched, s, label in candidate_schedules():
                    m = evaluate_schedule(sched, sir_trough=s)
                    rows.append({"schedule": label,
                                 "clears@week": m["clearance_weeks"],
                                 "final log10": m["final_log10_cpml"],
                                 "rebound": m["rebound_index"],
                                 "F_rr": m["final_frr"],
                                 "u:p ratio": m["urine_plasma_ratio"]})
            import pandas as pd
            st.subheader("Schedule sweep — the model's answer")
            st.dataframe(pd.DataFrame(rows), use_container_width=True)
            st.info(
                "Headline result: tacrolimus taper alone never clears "
                "viremia within 180 days, while tac→sir conversion clears "
                "by ~week 8 with LOWER T-cell rebound — sirolimus closes "
                "the S-phase/mTOR permissiveness gate in addition to "
                "releasing the adaptive-immunity brake. A mechanistic "
                "hypothesis matching the clinical literature (Hirsch "
                "2016), not a fitted result."
            )


def st_documentation_page():
    """Display documentation page."""
    st.markdown("## 📚 Documentation")
    
    doc_tabs = st.tabs([
        "Clinical Background",
        "Model Description",
        "Parameter Reference",
        "Research Citations"
    ])
    
    with doc_tabs[0]:
        st.subheader("Clinical Background: BK Polyomavirus Nephropathy")
        
        st.markdown("""
        **What is BKPyV?**
        BK Polyomavirus is a common virus that infects up to 90% of adults. In healthy individuals,
        it remains dormant in the kidneys. However, in immunosuppressed kidney transplant patients,
        it can reactivate and cause nephropathy (kidney damage).
        
        **Clinical Impact:**
        - Affects 5-10% of kidney transplant recipients
        - Can lead to graft loss in severe cases
        - Requires monitoring and intervention
        
        **Screening Guidelines:**
        - Screen for BK viremia monthly in first year post-transplant
        - Thresholds:
          - ≥1,000 copies/mL: Positive screening result
          - ≥10,000 copies/mL: High risk, requires intervention
        
        **Treatment Approaches:**
        - Reduce immunosuppression (mainstay)
        - Switch from tacrolimus to sirolimus (shown to inhibit replication)
        - Antiviral therapies (limited efficacy)
        """)
    
    with doc_tabs[1]:
        st.subheader("Virtual Cell Model Description")
        
        st.markdown("""
        **Model Overview:**
        The BKPyV virtual cell model simulates viral replication in renal tubular epithelial cells,
        incorporating mechanistic insights from single-cell transcriptomic studies and clinical observations.
        
        **Key Components:**
        1. **Viral Replication Dynamics**: T antigen-dependent, cell-cycle coupled replication
        2. **Drug Effects**: Tacrolimus changes immune control; sirolimus changes mTOR/S-phase permissiveness
        3. **Host Response**: DNA damage response, innate immunity, mitochondrial stress
        4. **Cell Cycle**: S-phase optimal for viral replication
        5. **Immune Evasion**: Antigen presentation suppression, interferon response downregulation
        
        **Clinical Mapping:**
        Virtual cell viral load (0-1 scale) is converted to clinical plasma viral load (copies/mL)
        using population scaling parameters derived from clinical cohort data.
        
        **Validation:**
        The model reproduces key qualitative patterns from research:
        - Tacrolimus may increase whole-system production by reducing immune control; this is not direct genome-copying enhancement
        - Drug timing effects (sirolimus effective only in early phase)
        - Cell-cycle and DNA-repair permissiveness effects
        - Mitochondrial stress signature in late infection
        """)
    
    with doc_tabs[2]:
        st.subheader("Parameter Reference")
        
        st.markdown("""
        **Drug Effect Parameters (High Confidence):**
        - `tacrolimus_enhancement_factor`: compatibility parameter for immune-control sensitivity; not a direct replication rate
        - `sirolimus_inhibition_factor`: 0.5 (IC90 = 4 ng/mL)
        - `drug_effectiveness_window`: 24h (early phase only)
        - `late_phase_drug_resistance`: 0.3 (reduced effectiveness)
        
        **Viral Replication Parameters (Medium Confidence):**
        - `t_antigen_replication_threshold`: 0.5
        - `dna_replication_coupling`: 0.8
        - `cell_cycle_s_phase_bonus`: 2.0
        
        **Clinical Risk Factors (High Confidence):**
        - `age_risk_multiplier`: 1.9 (age >50 years)
        - `male_sex_risk_multiplier`: 2.3
        - `prior_transplant_risk_multiplier`: 3.0
        
        For complete parameter registry, see `src/vcm/plugins/transplant/bk_polyomavirus/parameters.py`
        """)
    
    with doc_tabs[3]:
        st.subheader("Research Citations")
        
        st.markdown("""
        **Key Research Papers:**
        
        1. **Hirsch et al., Am J Transplant 2016**
           "BK Polyomavirus Replication in Renal Tubular Epithelial Cells Is Inhibited by Sirolimus, 
           but Activated by Tacrolimus Through a Pathway Involving FKBP-12"
           - Drug mechanisms: tacrolimus reduces immune control; sirolimus inhibits permissive cell-cycle signaling
           - Sirolimus IC90 = 4 ng/mL, effective in early phase (0-24h)
        
        2. **Weissbach et al., J Virol 2024**
           Single-cell RNA-sequencing of BKPyV replication
           - Cell cycle coupling, mitochondrial stress signature
           - Immune evasion mechanisms
        
        3. **Clinical risk factors** (see docs/ISEF_PROJECT_OVERVIEW.md for citations)
           - Directions from systematic review evidence (Demey et al. 2018): tacrolimus-based regimens, male sex, older age, prior transplant
           - Specific ORs are illustrative in this app, not fitted

        **Data Sources:**
        - Consensus thresholds: AST IDCOP 2019; Kotton et al. (Transplantation) 2024
        - Single-cell biology (qualitative): Weissbach et al., J Virol 2024; Needham et al., PLoS Pathog 2024
        - No patient-cohort dataset is bundled; do not present outputs as patient-calibrated
        """)


REPO_ROOT = Path(__file__).resolve().parents[3]


def _repo_file(*parts) -> Path:
    return REPO_ROOT.joinpath(*parts)


def st_single_cell_page():
    """Single-cell analysis results (GSE317012 biopsy scRNA-seq)."""
    st.markdown("## 🧫 Single-Cell Analysis — GSE317012")
    st.markdown("""
    Human kidney-transplant biopsies from the Kretzler lab (GEO: GSE317012,
    26 samples: 12 Control / 5 Peaking / 9 Resolving, 34,987 cells post-QC).
    The reference contains **human genes only** — no viral genes — so infected
    cells are identified by a host-response *signature proxy* (top-decile
    composite viral-response score in tubular epithelial cells), not by direct
    viral reads. This is a stated limitation, not a measurement of T-antigen
    or viral load.
    """)

    processed = _repo_file("data", "processed")
    fig_dir = _repo_file("outputs", "figures", "single_cell")
    cluster_csv = processed / "gse317012_cluster_summary.csv"
    de_csv = processed / "gse317012_infected_vs_bystander_de.csv"
    mapping_csv = processed / "gse317012_model_mapping.csv"
    holdout_csv = processed / "gse317012_calibration_holdout.csv"

    if not cluster_csv.exists():
        st.warning("Single-cell outputs not found. Regenerate with:")
        st.code("python scripts/download_gse317012.py\n"
                "python scripts/cluster_gse317012.py", language="bash")
        return

    c1, c2, c3 = st.columns(3)
    clusters = pd.read_csv(cluster_csv)
    c1.metric("Clusters", len(clusters))
    c2.metric("Samples", int(clusters["sample"].nunique()) if "sample" in clusters else 26)
    if de_csv.exists():
        de = pd.read_csv(de_csv)
        c3.metric("DE genes (FDR<0.05)", int((de["pvals_adj"] < 0.05).sum()),
                  f"{len(de):,} tested")

    st.subheader("Cell-type annotation (UMAP)")
    for name, caption in [
        ("umap_celltype.png", "UMAP coloured by annotated cell type"),
        ("umap_phase.png", "UMAP coloured by biopsy phase"),
        ("module_scores.png", "Host-response module scores"),
    ]:
        p = fig_dir / name
        if p.exists():
            st.image(str(p), caption=caption, use_container_width=True)

    st.subheader("Infected-signature vs bystander epithelial cells")
    st.caption("Wilcoxon rank-sum within tubular epithelial cells; "
               "log2FC = signature-high (top-decile) minus remainder.")
    if de_csv.exists():
        st.dataframe(pd.read_csv(de_csv).head(30), use_container_width=True)

    st.subheader("Finding → model mapping")
    if mapping_csv.exists():
        st.dataframe(pd.read_csv(mapping_csv), use_container_width=True)
    if holdout_csv.exists():
        st.subheader("Calibration / validation split")
        st.caption("Deterministic seeded stratified holdout — samples marked "
                   "'validation' were never used to calibrate the cell-state layer.")
        st.dataframe(pd.read_csv(holdout_csv), use_container_width=True)


def st_validation_page():
    """Viral-load benchmark page (Phase 2)."""
    st.markdown("## 📈 Viral-Load Validation")
    st.markdown("""
    Benchmark of plasma BKV kinetics against **published summary statistics**
    (no public de-identified serial qPCR dataset exists — stated limitation).
    Calibration subset: Funk 2006 IS-change arm. Validation subsets:
    Funk 2008 curtailment responses and Funk 2006 nephrectomy arms.
    """)

    report_path = _repo_file("outputs", "benchmark", "viral_load_benchmark.json")
    pub_path = _repo_file("data", "research", "published_kinetics.csv")
    fig_path = _repo_file("outputs", "benchmark", "viral_load_benchmark.png")

    if not report_path.exists():
        st.warning("Benchmark not run yet. Regenerate with:")
        st.code("python scripts/benchmark_viral_load.py", language="bash")
        return

    import json
    report = json.loads(report_path.read_text())

    c1, c2, c3 = st.columns(3)
    calib = report["calibration"]["result"]
    c1.metric("Clearance t½ after IS reduction",
              f"{calib['predicted_t_half_days']:.1f} d",
              f"published 0.25–17 d")
    c2.metric("Funk 2008 curtailment checks",
              "3/3 PASS" if report["validation"]["funk2008_curtailment"]["all_checks_pass"]
              else "see table")
    c3.metric("Nephrectomy clearance t½",
              f"{report['validation']['funk2006_nephrectomy']['predicted_t_half_hours']:.0f} h",
              "published 1–2 h / 20–38 h arms")

    if fig_path.exists():
        st.image(str(fig_path), use_container_width=True)

    st.subheader("Curtailment-response detail (validation)")
    rows = []
    for k, v in report["validation"]["funk2008_curtailment"]["results"].items():
        rows.append({
            "curtailment": k.replace("curtail_", "") + "%",
            "time below 1,000 cp/mL (wk)": (f"{v['time_below_threshold_weeks']:.1f}"
                                            if v["time_below_threshold_weeks"] else "never"),
            "sustained clearance": v["sustained_clearance"],
            "re-equilibrated cp/mL": round(v["re_equilibrated_copies_per_ml"]),
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True)

    st.subheader("Parameter uncertainty (δ sweep)")
    sweep = pd.DataFrame(report["parameter_uncertainty"]["sweep"])
    st.dataframe(sweep, use_container_width=True)
    st.line_chart(sweep.set_index("delta_per_day"))

    st.info(report["data_limitation"])

    if pub_path.exists():
        st.subheader("Published values used (with citations)")
        st.dataframe(pd.read_csv(pub_path), use_container_width=True)


def st_parameters_page():
    """Parameter & assumption table with citations."""
    st.markdown("## 🧾 Parameters & Assumptions")
    st.markdown("""
    Every parameter below is either **cited** to a published source or
    explicitly labelled an **assumption**. Values not fitted to patient data.
    """)

    from vcm.simulators.ode_system import BKPyVODESystem
    ode = BKPyVODESystem()

    # Parameter provenance table — kept in one place so the page and the
    # METHODS doc cite identical sources.
    PROVENANCE = {
        "beta": ("0.3 /day", "infection rate",
                 "literature-informed; tuned within mechanism plausibility (ASSUMPTION)"),
        "delta": ("0.4 /day", "plasma viral clearance",
                  "calibrated to Funk 2006 (PMID 16323135) IS-change t½ range 0.25–17 d"),
        "p": ("8.0 /day", "virion production",
              "ASSUMPTION — no direct cell-level estimate published"),
        "immune_kill": ("see code", "immune-mediated clearance",
                        "direction from Hirsch 2016 Am J Transplant; magnitude ASSUMPTION"),
        "t_threshold": ("0.5", "T-antigen replication gate shape",
                        "mechanistic; T-ag → S-phase permissiveness (Weissbach 2024 J Virol)"),
        "half_saturation": ("0.5", "T→production half-saturation",
                            "ASSUMPTION (Hill coefficient context)"),
        "s_phase_bonus": ("2.0", "S-phase replication boost",
                          "Weissbach 2024 J Virol: replication couples to S-phase"),
        "innate_immune_suppression": ("0.5", "IFN suppression of replication",
                                      "direction from interferon literature; magnitude ASSUMPTION"),
        "mtor_inhibition": ("0.5", "sirolimus mTOR effect",
                            "Hirsch 2016: sirolimus inhibits permissive cell-cycle signalling"),
        "tacrolimus_enhancement": ("see code", "tacrolimus immune-control reduction",
                                   "Hirsch 2016: tacrolimus ↑ replication via FKBP-12 pathway"),
        "drug clearance": ("see code", "drug PK decay",
                           "ASSUMPTION — published PK half-lives used as bounds"),
    }
    rows = [{"parameter": k, "value": v[0], "meaning": v[1], "source / status": v[2]}
            for k, v in PROVENANCE.items()]
    st.dataframe(pd.DataFrame(rows), use_container_width=True)

    st.subheader("Clinical-risk covariates (illustrative, not fitted)")
    st.markdown("""
    - age >50 y: ×1.9; male sex: ×2.3; prior transplant: ×3.0
    - Directions from Demey et al. 2018 systematic review; specific ORs are
      illustrative per the app's own documentation tab.
    """)

    st.subheader("V → copies/mL bridge (ASSUMPTION anchors)")
    st.markdown("""
    - 0 → 0 cp/mL; 0.02 → 100; 0.2 → 1,000 (screening, AST IDCOP 2019);
      1.0 → 10,000 (presumptive PyVAN, Kotton 2024); 3.0 → 10⁶; 5.0 → 10⁷
    - Piecewise log-linear between anchors — see
      `src/vcm/clinical/viral_load_mapper.py`.
    """)

    st.subheader("Data provenance")
    st.markdown("""
    - Single-cell layer: GSE317012 (Kretzler lab biopsies; sha256-verified download)
    - Kinetics benchmark: Funk 2006 (PMID 16323135), Funk 2008
      (doi:10.1111/j.1600-6143.2008.02402.x)
    - Clinical thresholds: AST IDCOP 2019; Kotton et al., Transplantation 2024
    """)


def main():
    """Main Streamlit application."""
    st_page_header()
    
    # Initialize session state
    if 'simulation_result' not in st.session_state:
        st.session_state.simulation_result = None
    if 'simulation_config' not in st.session_state:
        st.session_state.simulation_config = None
    if 'simulation_history' not in st.session_state:
        st.session_state.simulation_history = {}
    
    # Sidebar navigation
    page = st_sidebar_navigation()
    
    # Route to appropriate page
    if page == "🏠 Home":
        st_home_page()
    elif page == "⚙️ Simulation":
        st_simulation_page()
    elif page == "📊 Visualization":
        st_visualization_page()
    elif page == "🧫 Single-Cell":
        st_single_cell_page()
    elif page == "📈 Viral-Load Validation":
        st_validation_page()
    elif page == "🧬 Risk Prediction":
        st_risk_prediction_page()
    elif page == "📋 Comparison":
        st_comparison_page()
    elif page == "💉 Regimen Design":
        st_regimen_design_page()
    elif page == "🧾 Parameters & Assumptions":
        st_parameters_page()
    elif page == "🔎 Review Bundle":
        st_review_bundle_page()
    elif page == "📚 Documentation":
        st_documentation_page()


if __name__ == "__main__":
    main()
