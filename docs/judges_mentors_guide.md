> **NOTE (post-audit):** Sections describing parameter values, calibrations, or
> numbers may be historical. The canonical, current artifacts are:
> README.md, docs/BKPYV_MODEL_CARD.md, docs/ISEF_PROJECT_OVERVIEW.md, and
> the code itself (src/vcm/simulators/ode_system.py is the source of truth
> for the ODE model and its parameters). Time unit: days for the ODE engine.

# BKPyV Clinical Guide for Judges and Mentors

## Executive Summary

This project develops a **virtual cell model** to predict BK polyomavirus (BKPyV) replication in kidney transplant patients. BKPyV is a common virus that can cause kidney damage in transplant recipients, affecting up to 10% of patients. Our model uses **mechanistic simulations** of individual kidney cells, combined with **clinical data**, to predict viral replication kinetics and patient risk. This approach provides insights that could help clinicians make better treatment decisions.

**Impact**: This tool could help prevent kidney damage in transplant patients by predicting which patients are at highest risk and which treatments work best.

---

## What is BKPyV and Why Does It Matter?

### The Clinical Problem

**BK Polyomavirus (BKPyV)** is a common virus that infects most adults without causing symptoms. However, in **kidney transplant patients** who take immunosuppressive drugs to prevent organ rejection, the virus can reactivate and cause:

- **BKPyV Nephropathy (BKPyVAN)**: Kidney damage caused by the virus
- **Graft Loss**: In severe cases, patients may lose their transplanted kidney
- **Need for Dialysis**: Failed grafts require return to dialysis
- **Increased Healthcare Costs**: Treatment of complications is expensive

### Current Clinical Challenge

Doctors currently use:
- **Monthly blood tests** to detect viral load (copies of virus per mL of blood)
- **Guidelines** that define risk thresholds:
  - **≥1,000 copies/mL**: Requires closer monitoring
  - **≥10,000 copies/mL**: High risk, requires treatment
- **Treatment**: Reduce immunosuppressive drugs, but this increases rejection risk

**The Problem**: Current guidelines don't account for:
- Individual patient differences
- How different drugs affect viral replication
- When treatments should start for best outcomes
- How to balance infection risk vs rejection risk

---

## Our Solution: Virtual Cell Model

### What is a Virtual Cell Model?

Think of it as a **digital laboratory** where we can simulate how BKPyV behaves in individual kidney cells, similar to running experiments but in a computer instead of a wet lab.

**Key Components:**

1. **Virtual Kidney Cell**: A computer model of a kidney tubular epithelial cell (the cell type BKPyV infects)

2. **Viral Replication Simulation**: Models how the virus replicates inside the cell, including:
   - When viral genes turn on/off
   - How the virus uses the cell's machinery
   - Drug effects on viral replication

3. **Clinical Translation**: Converts cell-level results to blood-level viral load that doctors actually measure

4. **Risk Prediction**: Combines simulation results with patient data to predict individual risk

### How It Works

```
Patient Data → Virtual Cell Simulation → Clinical Viral Load → Risk Prediction
     ↓                ↓                      ↓                 ↓
Age, sex, drugs   Viral replication    Blood copies/mL    Risk category
HLA matching      kinetics            (1,000-10,000)     (low/med/high)
```

---

## Research Foundation

Our model is built on **published clinical research**, not assumptions:

### Key Research Findings Incorporated:

1. **Drug effects** (Hirsch et al., Am J Transplant 2016;16(3):821-832, *in vitro* in primary renal cells):
   - **Tacrolimus promotes viral *production***: it weakens immune control of infected cells, so more virions accumulate (it does not directly promote genome copying). Wording follows expert review of this project (S. Thompson lab, UAB, Sep 2026).
   - **Sirolimus may prevent viral *production*** by reducing S-phase entry/T-antigen expression (mTOR/S6K pathway; IC90 ≈ 4 ng/mL).
   - Clinical standpoint: tacrolimus-based regimens are associated with higher BKPyV risk than some comparators; whether mTOR-inhibitor regimens *protect* at the patient level remains debated in meta-analyses - the model claims only direction and mechanism framing.
   - Terminology note (per expert review): within this project, *replication* means intracellular genome copying; system-level changes in virions/viremia are *production*.

2. **Cell cycle dependence** (Needham et al., PLoS Pathogens 2024; Weissbach et al., J Virol 2024):
   - Host S phase *precedes* large T-antigen accumulation (the "S-phase gate" in the model)
   - BKPyV replication requires host DNA-replication machinery
   - Infected cells show a mitochondrial-stress signature

3. **Clinical risk factors** (Demey et al., J Clin Virol 2018 systematic review; Yamauchi et al., Renal Failure 2025):
   - Older age, male sex, prior transplant, tacrolimus-based regimen (directionally consistent factors; specific odds-ratio magnitudes are study-dependent and no particular OR set is hard-coded as ground truth)

4. **Drug timing** (Hirsch 2016, in vitro):
   - Sirolimus inhibition was strongest when applied before productive replication; the ODE implements this as a T-antigen-level--gated window rather than a hard 24-hour clock.
   - This timing is critical for treatment decisions

### Data Sources:

- **Clinical tables** extracted from research papers (7 tables from 4 DOCX files)
- **Single-cell RNA sequencing** data from BKPyV-infected cells
- **Published cohort studies** on BKPyV risk factors
- **Clinical guidelines** (KDIGO, AST consensus statements)

---

## Model Features

### 1. Mechanistic Drug Effects

Our model shows **why** different drugs have different effects:

**Tacrolimus**:
- Potent calcineurin inhibitor that weakens T-cell immune targeting of infected cells
- In the model, this weaker immune pressure **promotes viral production** relative to baseline
- Literature odds ratio used in the synthetic cohort model: OR ≈ 2.3 (Fang 2022, PMC9428263)

**Sirolimus**:  
- mTOR inhibitor; in the model it **may prevent production** via the S-phase entry gate (host-cell permissiveness)
- In vitro IC90 ≈ 4 ng/mL reported by Hirsch 2016 (AJT-16-821); the model treats this as a mechanistic prior, not a clinical dosing claim
- Drug timing matters qualitatively (early viral-gene-expression window); no hard cutoff or drug-resistance claim

**Why This Matters**:
- Helps doctors choose the right immunosuppressive regimen
- Predicts which patients will respond to drug changes
- Provides mechanistic explanation for clinical observations

### 2. Clinical Viral Load Mapping

Our model converts **cell-level viral load** (0-1 scale) to **clinical blood viral load** (copies/mL):

**Conversion Process**:
1. Calculate infected cell count from virtual load
2. Estimate viral copies per infected cell
3. Apply plasma distribution factors
4. Account for kidney clearance
5. Match clinical thresholds (1,000 and 10,000 copies/mL)

**Validation**: The model reproduces clinical patterns:
- Time to first detectable viremia (~1 month post-transplant)
- Peak viral load timing (~2 months)
- Clearance rates after treatment

### 3. Risk Prediction

Our model combines **traditional clinical factors** with **virtual cell simulation features**:

**Traditional Clinical Factors**:
- Age, sex, prior transplant
- HLA mismatch, diabetes
- Drug regimen, biomarkers (GcfDNA)

**Virtual Cell Features**:
- Simulated peak viral load
- Time to peak
- Viral clearance rate
- Replication phase duration
- Mitochondrial stress level

**Prediction Models**:
- **Clinical-only model**: Uses just traditional factors (baseline)
- **Enhanced model**: Clinical + virtual cell features (improved accuracy)

### 4. Clinical Decision Support

The system provides **actionable recommendations** based on viral load:

**Low Risk (<1,000 copies/mL)**:
- Continue routine monthly monitoring
- Maintain current immunosuppression

**Medium Risk (1,000-10,000 copies/mL)**:
- Increase monitoring to biweekly
- Consider tacrolimus dose reduction
- Monitor kidney function closely

**High Risk (≥10,000 copies/mL)**:
- Weekly monitoring
- Reduce immunosuppression significantly
- Consider switching from tacrolimus to sirolimus
- Early nephrology consultation

---

## How Judges Can Evaluate This Project

### Technical Rigor

✅ **Evidence-Based Parameters**: All 24 parameters documented with research sources and confidence levels

✅ **Validation**: Model reproduces qualitative patterns from literature (drug effects, timing dependence)

✅ **Code Quality**: Clean architecture, type hints, comprehensive documentation

✅ **Reproducibility**: Configuration files, parameter registry, validation scripts

### Scientific Method

✅ **Literature Review**: Comprehensive analysis of BKPyV research papers

✅ **Data Extraction**: Systematic extraction of clinical data from published tables

✅ **Model Building**: Mechanistic simulation based on biological understanding

✅ **Validation**: Qualitative validation against published findings

✅ **Clinical Relevance**: Aligns with clinical guidelines and thresholds

### Originality

✅ **Novel Approach**: Combines virtual cell simulation with clinical risk prediction

✅ **Mechanistic Insight**: Explains drug effects at molecular level (FKBP-12 pathway)

✅ **Timing Effects**: Models drug timing dependence not captured by clinical models

✅ **Integration**: Bridges cellular-level biology to patient-level outcomes

### Real-World Impact

✅ **Clinical Decision Support**: Could help doctors choose better treatments

✅ **Risk Stratification**: Identifies high-risk patients early

✅ **Personalized Medicine**: Accounts for individual patient differences

✅ **Drug Selection**: Helps choose tacrolimus vs sirolimus based on patient profile

---

## How to Use the System

### For Judges: Quick Demo

**Option 1: Streamlit Web Interface** (Recommended)

```bash
# Navigate to project directory
cd virtual-cell-model

# Install required packages
pip install streamlit plotly pandas numpy pyyaml

# Run the web interface
streamlit run src/vcm/ui/streamlit_app.py
```

**Interface Features**:
- **Home**: Overview and quick start guide
- **Simulation**: Configure and run BKPyV simulations
- **Visualization**: Interactive plots of viral load trajectories
- **Risk Prediction**: Patient-specific risk assessment
- **Comparison**: Compare different drug scenarios
- **Documentation**: Complete reference guide

**Option 2: Validation Script**

```bash
python validate_bkpyv.py
```

This runs automated validation showing the model reproduces research findings:
- Tacrolimus > Control in viral replication
- Sirolimus < Tacrolimus in early replication  
- High cell-cycle permissiveness > Low permissiveness
- Late replication > Early in mitochondrial stress

### For Mentors: Technical Overview

**Project Structure**:
```
virtual-cell-model/
├── src/vcm/
│   ├── clinical/              # Clinical integration modules
│   │   ├── viral_load_mapper.py
│   │   ├── calibration.py
│   │   ├── risk_prediction.py
│   │   └── thresholds.py
│   ├── plugins/
│   │   └── transplant/
│   │       └── bk_polyomavirus/
│   │           └── parameters.py
│   ├── simulators/
│   │   └── bkpyv_simulator.py
│   └── ui/
│       └── streamlit_app.py
├── data/                     # Research data
├── configs/                  # Simulation configurations
├── docs/                     # Documentation
└── validate_bkpyv.py        # Validation script
```

**Key Files to Review**:

1. **`docs/bkpyv_research_to_model_map.md`**: Maps research findings to model implementation

2. **`src/vcm/plugins/transplant/bk_polyomavirus/parameters.py`**: Parameter registry with evidence sources

3. **`src/vcm/clinical/viral_load_mapper.py`**: Converts cell load to clinical viral load

4. **`validate_bkpyv.py`**: Automated validation script

**Running Experiments**:

```python
# Example: Run tacrolimus vs sirolimus comparison
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_simulator import BKPyVSimulator
from vcm.clinical.viral_load_mapper import ClinicalViralLoadMapper

# Load plugin and simulator
plugin = BKPolyomavirusPlugin()
initial_state = plugin.create_initial_state()

# Tacrolimus simulation
sim_tac = BKPyVSimulator({'tacrolimus_enhancement_factor': 2.3})
result_tac = sim_tac.simulate(initial_state)

# Sirolimus simulation  
sim_siro = BKPyVSimulator({'sirolimus_inhibition_factor': 0.5})
result_siro = sim_siro.simulate(initial_state)

# Compare results
mapper = ClinicalViralLoadMapper()
# Convert and compare viral loads...
```

---

## Key Scientific Questions Answered

### 1. Why does tacrolimus increase BKPyV risk?

**Answer**: Tacrolimus is a calcineurin inhibitor that weakens T-cell immune targeting of infected cells. In the model this weaker immune pressure **promotes viral production** relative to baseline. The synthetic-cohort risk model uses a literature odds ratio of OR ≈ 2.3 for tacrolimus use (Fang 2022, PMC9428263).

### 2. Can sirolimus protect against BKPyV?

**Answer**: In the model, sirolimus **may prevent production** via the S-phase entry gate on host-cell permissiveness (mTOR-dependent cell-cycle state). Its effect is time-windowed around early viral gene expression — a graded mechanistic window, not a hard cutoff, and the model does not claim drug resistance emerges.

### 3. What factors increase BKPyV risk?

**Answer**: The synthetic-cohort risk model uses literature-sourced odds ratios (Fang 2022, PMC9428263):
- Tacrolimus use (OR ≈ 2.3)
- Prior kidney transplant (OR ≈ 2.1)
- Male sex (OR ≈ 1.6)
- Diabetes (OR ≈ 1.3)
- Higher simulated peak viral load (derived from VCM features; higher peak → higher risk)

### 4. How does the model translate cell simulation to clinical outcomes?

**Answer**: We use a **population scaling model** that:
- Estimates infected cell count from virtual load
- Calculates viral copies per cell
- Applies plasma distribution and kidney clearance factors
- Converts to blood viral load (copies/mL) that doctors measure

### 5. How does this compare to existing clinical models?

**Answer**: Traditional models use only clinical factors (age, sex, drugs). Our model adds **virtual cell simulation features** (peak viral load, clearance rate, mitochondrial stress) that capture the biology of viral replication, potentially improving prediction accuracy.

---

## Limitations and Future Work

### Current Limitations

1. **Parameter Calibration**: Some parameters are phenomenological (estimated) rather than directly measured
2. **Single-Cell Focus**: Models individual cells, not full tissue architecture
3. **Simplified Pharmacokinetics**: Drug levels are binary (present/absent), not concentration-dependent
4. **No Clinical Validation**: Model matches qualitative patterns but needs clinical patient data for quantitative validation

### Future Improvements

1. **Calibrate with Real Patient Data**: Use actual BKPyV patient viral load trajectories
2. **Add Multi-Cell Modeling**: Include immune cells for tissue-level dynamics
3. **Drug Pharmacokinetics**: Model concentration-dependent drug effects
4. **Clinical Trial**: Validate predictions against actual patient outcomes

---

## How This Project Could Be Used

### Clinical Applications

1. **Pre-Transplant Planning**: Help choose immunosuppression regimen based on patient risk factors
2. **Early Risk Detection**: Identify high-risk patients before viral load becomes detectable
3. **Drug Selection**: Guide choice between tacrolimus and sirolimus
4. **Treatment Timing**: Determine optimal timing for drug changes
5. **Patient Education**: Explain BKPyV risk in individualized terms

### Research Applications

1. **Mechanistic Studies**: Test hypotheses about BKPyV biology
2. **Drug Development**: Simulate new drug effects
3. **Clinical Trial Design**: Identify optimal patient subgroups
4. **Personalized Medicine**: Develop individualized treatment protocols

---

## Questions for the Student

### About the Model

- **How did you choose which biological processes to include?**
  - We prioritized processes with strong experimental evidence and clear clinical relevance
  
- **How do you validate your model against real data?**
  - We use qualitative validation (does it reproduce known patterns?) and are working on quantitative validation with patient data

- **What happens if the model's predictions don't match real patient outcomes?**
  - We can calibrate parameters using the clinical calibration module to improve accuracy

### About the Implementation

- **Why did you choose Python for implementation?**
  - Python has excellent scientific computing libraries and is widely used in computational biology
  
- **How long did it take to build this system?**
  - [Student's actual answer - this shows project timeline]
  
- **What was the most challenging technical problem?**
  - [Student's reflection on main challenges]

### About the Science

- **How does this model compare to other BKPyV prediction tools?**
  - Most existing tools use statistical models; our model adds mechanistic understanding
  
- **Could this approach be applied to other viruses?**
  - Yes, the framework is designed to be extendable to other transplant-related viruses

- **What would you need to make this clinically ready?**
  - Validation with patient data, regulatory approval, integration with hospital systems

---

## Resources for Further Reading

### Key Research Papers

1. **Hirsch et al., Am J Transplant 2016**
   - Drug mechanisms: tacrolimus activates, sirolimus inhibits
   - FKBP-12 pathway explanation

2. **Weissbach et al., J Virol 2024**
   - Single-cell transcriptomics of BKPyV infection
   - Cell cycle coupling and mitochondrial stress

3. **Clinical Cohort Studies**
   - Risk factor analysis (age, sex, prior transplant)
   - Clinical outcomes and treatment responses

### Clinical Guidelines

- **KDIGO 2023**: Kidney transplant recipient guidelines
- **AST 2019**: BK virus consensus statement
- **AST 2021**: BK virus nephropathy guidelines

### Project Documentation

- **`docs/bkpyv_research_to_model_map.md`**: Detailed research-to-model mapping
- **`docs/bkpyv_progress_update.md`**: Technical progress and implementation details
- **`data/README.md`**: Data sources and processing methods

---

## Conclusion

This project represents a **serious, research-informed approach** to BKPyV prediction that bridges **mechanistic biology** with **clinical decision-making**. By using virtual cell simulations grounded in published research, we can:

1. **Explain** drug effects at the molecular level
2. **Predict** individual patient risk
3. **Guide** clinical decision-making
4. **Personalize** treatment strategies

The system is **research-grounded**, **clinically relevant**, and **technically rigorous**, making it suitable for ISEF-level competition while maintaining the potential for real-world impact in transplant patient care.

---

## Contact and Support

For questions about this project:
- Review the comprehensive documentation in the `docs/` folder
- Examine the parameter registry in `src/vcm/plugins/transplant/bk_polyomavirus/parameters.py`
- Run the validation script to see the model reproduce research findings
- Use the Streamlit interface for interactive exploration

**Best of luck with your evaluation!**