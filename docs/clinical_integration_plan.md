> **NOTE (post-audit):** Sections describing parameter values, calibrations, or
> numbers may be historical. The canonical, current artifacts are:
> README.md, docs/BKPYV_MODEL_CARD.md, docs/ISEF_PROJECT_OVERVIEW.md, and
> the code itself (src/vcm/simulators/ode_system.py is the source of truth
> for the ODE model and its parameters). Time unit: days for the ODE engine.

# Comprehensive Clinical Integration & Frontend UI Development Plan

## Overview
Transform the BKPyV plugin into a clinically relevant, user-friendly system that bridges virtual cell simulations to real-world clinical decision-making.

## Phase 1: Clinical Viral Load Integration (HIGH PRIORITY)
### 1.1 Design Viral Load Mapping System
**Goal**: Map virtual cell viral load (0-1 scale) to plasma viral load (copies/mL)

**Approach**:
- Use clinical cohort data to establish conversion factors
- Implement population-level scaling (single cell → plasma copies/mL)
- Incorporate log-scale transformations for clinical relevance
- Account for renal clearance and plasma distribution

**Components**:
- `ClinicalViralLoadMapper` class in BKPyV plugin
- Population scaling parameters (cells/mL, viral copies/cell)
- Log-scale conversion formulas
- Clinical threshold integration (1,000, 10,000 copies/mL)

**Research Data Sources**:
- PubMed studies on BK viremia kinetics
- PMC articles on DNAemia thresholds
- Clinical guidelines for screening protocols

### 1.2 Implement Clinical Parameter Calibration
**Goal**: Calibrate simulation parameters to match clinical cohort summaries

**Approach**:
- Use extracted clinical table data from DOCX files
- Implement optimization algorithm for parameter fitting
- Match timing: time to first viremia, peak timing
- Match magnitude: peak viral load, clinical thresholds
- Match kinetics: rise/fall patterns with immunosuppression changes

**Components**:
- `ClinicalCalibration` class with optimization methods
- Parameter sensitivity analysis
- Clinical trajectory matching functions
- Calibration metrics (MSE, AUC comparison)

**Key Parameters to Calibrate**:
- Infection rate (cells/time)
- Viral replication rate
- Immune clearance rate
- Drug effect kinetics
- Population scaling factors

### 1.3 Integrate Clinical Thresholds
**Goal**: Implement clinical decision thresholds (1,000 and 10,000 copies/mL)

**Approach**:
- Map virtual cell states to clinical risk categories
- Implement automated risk stratification
- Add clinical alerting system
- Generate clinical reports

**Components**:
- `ClinicalRiskStratifier` class
- Risk category logic (low, medium, high)
- Clinical alert system
- Report generation templates

## Phase 2: Risk Prediction Module (HIGH PRIORITY)
### 2.1 Build ML Risk Prediction System
**Goal**: Create risk prediction model combining clinical covariates + virtual cell features

**Approach**:
- Implement logistic regression baseline (clinical covariates only)
- Implement gradient boosting model (clinical + virtual cell features)
- Use published risk models as benchmarks
- Compare performance on held-out data

**Clinical Covariates** (from extracted data):
- Age (OR 1.75-1.99)
- Sex (OR 2.22-2.42 for male)
- Prior transplant (OR 2.79-3.28)
- HLA mismatch (OR 1.30 for 4-6 mismatch)
- Diabetes status
- Tacrolimus use
- GcfDNA levels (delta GcfDNA AUC 0.83)

**Virtual Cell Features**:
- Simulated peak viral load
- Area under viral load curve (AUC)
- Time to partial clearance
- Viral replication phase duration
- Mitochondrial stress trajectory
- Antigen presentation suppression
- Drug effectiveness metrics

**Components**:
- `RiskPredictionModule` class
- ML model implementations (logistic regression, gradient boosting)
- Feature engineering pipeline
- Model comparison framework
- Performance metrics (AUC, calibration, Brier score)

### 2.2 Clinical Model Benchmarking
**Goal**: Compare virtual cell-enhanced models against clinical-only models

**Approach**:
- Recreate baseline BKPyVAN risk score from literature
- Implement dynamic prediction models
- Cross-validation on clinical data
- Document experimental setup clearly

**Components**:
- Baseline clinical risk model
- Enhanced risk model with virtual cell features
- Model comparison pipeline
- Statistical significance testing
- Documentation framework

## Phase 3: Frontend UI Development (HIGH PRIORITY)
### 3.1 Create Web-Based UI
**Goal**: User-friendly interface for running BKPyV simulations and interpreting results

**Technology Stack**:
- Streamlit (Python-based, easy to deploy)
- FastAPI (backend API)
- React/Vue (advanced frontend, optional)

**UI Components**:
- Configuration panel (parameter adjustment)
- Scenario selector (baseline, drug exposures, risk profiles)
- Real-time simulation visualization
- Viral load trajectory plotting
- Clinical risk display
- Parameter sensitivity analysis
- Comparative scenario analysis

### 3.2 Interactive Visualization
**Goal**: Provide clear, publication-ready visualizations

**Visualization Types**:
- Viral load trajectory plots (log scale)
- Drug comparison plots (tacrolimus vs sirolimus)
- Phase transition diagrams (early vs late)
- Mitochondrial stress over time
- Risk probability heatmaps
- Parameter sensitivity tornado plots
- Clinical outcome predictions

**Components**:
- `VisualizationModule` class
- Plotly/Altair integration
- Interactive chart generation
- Export functionality (PNG, PDF, SVG)

### 3.3 User Documentation
**Goal**: Make the system accessible to judges, mentors, and researchers

**Documentation Structure**:
- Quick start guide (step-by-step instructions)
- Clinical background (BKPyV biology, nephropathy)
- Model explanation (virtual cell to clinical mapping)
- Parameter reference (research-grounded values)
- Case studies (typical use cases)
- Troubleshooting guide

## Phase 4: VCM Architecture Refactoring (MEDIUM PRIORITY)
### 4.1 Modular Plugin System
**Goal**: Enable easy addition of new disease modules (e.g., RSV)

**Current Issues**:
- BKPyV-specific code mixed with core VCM
- Hard-coded assumptions about cell types
- No clear plugin interface

**Refactoring Goals**:
- Define plugin interface with CellState schema declarations
- Separate simulator logic from plugin logic
- Make experiment runner plugin-agnostic
- Enable hot-swapping of disease modules

**Components**:
- `BasePlugin` interface with required methods
- `CellStateSchema` definition system
- Plugin registry with dependency resolution
- Simulator factory pattern
- Configuration schema validation

### 4.2 Experiment Runner Enhancement
**Goal**: Support multi-scenario experiments with standardized outputs

**Features**:
- Batch scenario execution
- Parameter sweep capability
- Parallel execution support
- Standardized JSON/CSV outputs
- Result aggregation and comparison
- Automated report generation

**Components**:
- `ExperimentRunner` with scenario management
- `ParameterSweep` class
- `ResultAggregator` for multi-scenario analysis
- Standard output format schemas
- Report generation templates

## Phase 5: Documentation Suite (MEDIUM PRIORITY)
### 5.1 Technical Documentation
**Goals**:
- BKPyV biology and clinical problem explanation
- Virtual cell to clinical load mapping
- Step-by-step experiment execution
- Parameter reference with research citations
- Architecture and API documentation

**Components**:
- BKPyV Clinical Guide (for judges/mentors)
- Technical Reference (for developers)
- API Documentation (for advanced users)
- Parameter Registry (already created)
- Research-to-Model Map (already created)

### 5.2 User Documentation
**Goals**:
- Easy onboarding for new users
- Clear step-by-step instructions
- Example workflows
- Troubleshooting common issues

**Components**:
- Quick Start Guide
- Tutorial Notebooks
- Example Configurations
- FAQ Document
- Video Tutorials (optional)

## Implementation Priority Order

### IMMEDIATE (This Session)
1. ✅ Create comprehensive plan (this document)
2. 🔄 Implement Clinical Viral Load Mapper
3. 🔄 Implement Clinical Calibration Script
4. 🔄 Create Streamlit Frontend UI

### SHORT TERM (Next 1-2 sessions)
5. Build Risk Prediction Module with ML
6. Integrate Clinical Thresholds
7. Create Interactive Visualizations
8. Write Quick Start Guide

### MEDIUM TERM (Future sessions)
9. Refactor VCM Architecture for modularity
10. Enhance Experiment Runner
11. Write comprehensive documentation suite
12. Add second disease module (RSV) as test case

## Success Criteria

### Clinical Integration
- ✅ Viral load mapping reproduces clinical timing and magnitude
- ✅ Calibration script matches cohort summaries within 10-20%
- ✅ Clinical thresholds (1,000, 10,000 copies/mL) correctly stratified
- ✅ Risk prediction model outperforms clinical-only baseline

### Frontend UI
- ✅ Web UI allows easy configuration and execution
- ✅ Visualizations are clear and publication-ready
- ✅ Non-expert users can run simulations
- ✅ Results are interpretable for clinical decision-making

### Architecture
- ✅ Plugin system supports easy addition of new diseases
- ✅ Experiment runner handles multi-scenario batches
- ✅ Standardized outputs enable easy analysis
- ✅ Documentation is comprehensive and accessible

## File Structure Plan

```
virtual-cell-model/
├── src/vcm/
│   ├── clinical/
│   │   ├── viral_load_mapper.py      # NEW: Virtual cell → plasma copies/mL
│   │   ├── calibration.py              # NEW: Clinical parameter calibration
│   │   ├── risk_prediction.py          # NEW: ML risk prediction module
│   │   └── thresholds.py               # NEW: Clinical threshold management
│   ├── plugins/
│   │   ├── base.py                     # REFACTOR: Improved plugin interface
│   │   └── transplant/
│   │       └── bk_polyomavirus/
│   │           ├── viral_load_mapper.py # NEW: BKPyV-specific viral load mapping
│   │           └── risk_features.py    # NEW: Virtual cell feature extraction
│   ├── ui/
│   │   ├── streamlit_app.py            # NEW: Main Streamlit application
│   │   ├── visualization.py            # NEW: Visualization components
│   │   └── api.py                      # NEW: FastAPI backend
│   └── experiments/
│       └── runner.py                   # ENHANCE: Multi-scenario support
├── docs/
│   ├── clinical_guide.md               # NEW: Clinical background for judges
│   ├── quick_start.md                  # NEW: Step-by-step user guide
│   └── api_reference.md                # NEW: API documentation
├── notebooks/
│   ├── calibration_demo.ipynb          # NEW: Calibration demonstration
│   └── risk_prediction_demo.ipynb       # NEW: Risk prediction demonstration
├── ui/                                 # NEW: Frontend assets
│   ├── static/
│   └── templates/
└── validate_clinical.py                # NEW: Clinical validation script
```

## Next Steps

I will now implement the highest priority components in this order:

1. **Clinical Viral Load Mapper** - Bridging virtual cell to clinical relevance
2. **Clinical Calibration Script** - Matching simulation to clinical data
3. **Streamlit Frontend UI** - Making the system usable
4. **Risk Prediction Module** - Demonstrating clinical utility

This approach ensures we create immediate value while building toward the comprehensive vision.