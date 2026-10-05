> **NOTE (post-audit):** Sections describing parameter values, calibrations, or
> numbers may be historical. The canonical, current artifacts are:
> README.md, docs/BKPYV_MODEL_CARD.md, docs/ISEF_PROJECT_OVERVIEW.md, and
> the code itself (src/vcm/simulators/ode_system.py is the source of truth
> for the ODE model and its parameters). Time unit: days for the ODE engine.

# BKPyV Virtual Cell Model - Technical Documentation

## Table of Contents

1. [Installation Guide](#installation-guide)
2. [Quick Start Tutorial](#quick-start-tutorial)
3. [Architecture Overview](#architecture-overview)
4. [API Reference](#api-reference)
5. [Configuration Guide](#configuration-guide)
6. [Development Guide](#development-guide)
7. [Troubleshooting](#troubleshooting)

---

## Installation Guide

### System Requirements

- **Python**: >=3.10 (pyproject `requires-python`)
- **Operating System**: macOS, Linux, Windows
- **Memory**: 4GB RAM minimum, 8GB recommended
- **Storage**: 500MB for project files (plus ~2GB if extracting the GSE317012 raw archive)

### Prerequisites

```bash
# Python >=3.10
python --version

# pip (Python package manager)
pip --version
```

### Installation Steps

#### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/virtual-cell-model.git
cd virtual-cell-model
```

#### 2. Create Virtual Environment

```bash
# Create the isolated environment (the project's .venv is gitignored)
python -m venv .venv

# Activate it
# On macOS/Linux:
source .venv/bin/activate
# On Windows:
.venv\Scripts\activate
```

#### 3. Install Dependencies

```bash
# Editable install with the dev + UI extras (pinned in pyproject.toml)
pip install -e ".[dev,ui]"
```

#### 4. Verify Installation

```bash
# Test Python imports (inside the activated .venv)
python -c "import vcm; print('VCM imported successfully')"

# Run validation script
python validate_bkpyv.py
```

### Development Installation

For contributors who want to modify the code:

```bash
# Install in development mode
pip install -e .

# Install development tools
pip install pytest black flake8 mypy
```

---

## Quick Start Tutorial

### Running Your First Simulation

#### Option 1: Using Configuration Files (Recommended)

```python
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_simulator import BKPyVSimulator
from vcm.core.models import ExperimentConfig
import yaml

# Load configuration
with open('configs/bkpyv_infection.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Create plugin and simulator
plugin = BKPolyomavirusPlugin()
initial_state = plugin.create_initial_state()

# Create simulator with parameters
sim_params = config.get('simulator_parameters', {})
simulator = BKPyVSimulator(sim_params)

# Run simulation
result = simulator.simulate(
    initial_state=initial_state,
    n_steps=100,
    timestep=1.0
)

# Access results
final_viral_load = result.final_state.metadata['viral_load']
print(f"Final viral load: {final_viral_load}")
```

#### Option 2: Using the Streamlit Interface (Easiest)

```bash
streamlit run src/vcm/ui/streamlit_app.py
```

This opens a web interface at `http://localhost:8501` with:
- Interactive parameter adjustment
- Scenario selection
- Real-time visualization
- Risk prediction tools

#### Option 3: Using the Validation Script

```bash
python validate_bkpyv.py
```

This runs automated validation showing the model reproduces research findings.

### Converting Virtual Load to Clinical Viral Load

```python
from vcm.clinical.viral_load_mapper import ClinicalViralLoadMapper

# Create mapper
mapper = ClinicalViralLoadMapper()

# Convert virtual cell load (0-1 scale) to clinical load (copies/mL)
virtual_load = 0.8
clinical_load = mapper.virtual_to_plasma_viral_load(virtual_load)

print(f"Virtual load {virtual_load} → Clinical load {clinical_load:.0f} copies/mL")

# Get risk stratification
risk_assessment = mapper.get_clinical_risk_stratification(clinical_load)
print(f"Risk category: {risk_assessment['risk_category']}")
print(f"Recommendations: {risk_assessment['recommendations']}")
```

### Running Risk Prediction

```python
from vcm.clinical.risk_prediction import (
    RiskPredictionModule,
    ClinicalCovariates,
    VirtualCellFeatures
)

# Create risk prediction module
module = RiskPredictionModule()

# Define patient clinical data
clinical_data = ClinicalCovariates(
    age=55.0,
    sex='male',
    prior_transplant=True,
    hla_mismatch=4,
    donor_type='deceased',
    diabetes=True,
    tacrolimus_use=True,
    sirolimus_use=False,
    induction_agent='basiliximab',
)

# Run risk prediction (requires trained models)
# For demo, this uses a simplified version
# result = module.predict_risk(clinical_data, virtual_cell_data)
```

---

## Architecture Overview

### System Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    User Interface Layer                      │
│  (Streamlit Web Interface | Python API | Jupyter Notebooks) │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                   Clinical Integration Layer                  │
│  Viral Load Mapper | Calibration | Risk Prediction          │
│  Thresholds Management | Clinical Validation                │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                      Core Simulation Layer                   │
│  Experiment Runner | BKPyV Simulator | Configuration      │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                       Plugin Layer                            │
│  BKPyV Plugin | Parameter Registry | State Management      │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                    Data Layer                                │
│  Research Data | Clinical Tables | Configuration Files    │
└─────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

#### 1. Clinical Integration Layer
- **`viral_load_mapper.py`**: Converts virtual cell loads to clinical viral load
- **`calibration.py`**: Calibrates simulation parameters to match clinical data
- **`risk_prediction.py`**: ML-based risk prediction combining clinical + simulation features
- **`thresholds.py`**: Manages clinical thresholds and risk stratification

#### 2. Core Simulation Layer
- **`bkpyv_simulator.py`**: BKPyV-specific simulation engine
- **`runner.py`**: Experiment orchestration and result management

#### 3. Plugin Layer
- **`bk_polyomavirus.py`**: BKPyV-specific state definition and parameters
- **`parameters.py`**: Central parameter registry with evidence documentation

#### 4. Data Layer
- **`data/research/`**: Extracted clinical tables and research papers
- **`configs/`**: YAML configuration files for different scenarios

### Data Flow

```
Patient Data → Clinical Covariates → Risk Prediction Module → Risk Score
Configuration → BKPyV Plugin → Initial State → Simulator → Trajectory
Trajectory → Viral Load Mapper → Clinical Viral Load → Risk Stratification
Clinical Data → Calibration Module → Optimized Parameters → Improved Simulations
```

---

## API Reference

### BKPyV Plugin API

#### `BKPolyomavirusPlugin`

```python
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin

plugin = BKPolyomavirusPlugin()

# Create initial cell state
initial_state = plugin.create_initial_state()

# Get supported perturbations
perturbations = plugin.get_supported_perturbations()

# Get default configuration
config = plugin.get_default_config()
```

#### Key Methods

- **`create_initial_state(config=None)`**: Creates initial cell state with BKPyV-specific genes, proteins, pathways
- **`get_supported_perturbations()`**: Returns list of available perturbations (infection, drugs)
- **`get_default_config()`**: Returns default experiment configuration

### Simulator API

#### `BKPyVSimulator`

```python
from vcm.simulators.bkpyv_simulator import BKPyVSimulator

# Initialize with parameters
simulator = BKPyVSimulator({
    'tacrolimus_enhancement_factor': 2.0,
    'sirolimus_inhibition_factor': 0.5,
    't_antigen_replication_threshold': 0.5,
})

# Run simulation
result = simulator.simulate(
    initial_state=initial_state,
    perturbation=None,
    n_steps=100,
    timestep=1.0
)
```

#### Key Parameters

- **`tacrolimus_enhancement_factor`**: Multiplier for replication under tacrolimus (default: 2.0)
- **`sirolimus_inhibition_factor`**: Multiplier for replication under sirolimus (default: 0.5)
- **`t_antigen_replication_threshold`**: T antigen level needed for replication (default: 0.5)
- **`cell_cycle_s_phase_bonus`**: Replication bonus in S-phase (default: 2.0)
- **`dna_replication_coupling`**: Host-viral coupling strength (default: 0.8)

### Clinical Integration API

#### `ClinicalViralLoadMapper`

```python
from vcm.clinical.viral_load_mapper import ClinicalViralLoadMapper

mapper = ClinicalViralLoadMapper()

# Convert virtual to clinical load
clinical_load = mapper.virtual_to_plasma_viral_load(
    virtual_viral_load=0.8,
    infected_cell_fraction=0.01,
    time_hours=24.0,
    drug_effect=1.0,
    use_log_scale=True
)

# Get risk stratification
risk = mapper.get_clinical_risk_stratification(clinical_load)
```

#### Key Methods

- **`virtual_to_plasma_viral_load()`**: Converts virtual cell load to clinical copies/mL
- **`get_clinical_risk_stratification()`**: Returns risk category and recommendations
- **`simulate_viral_load_trajectory()`**: Generates clinical trajectory from simulation

#### `RiskPredictionModule`

```python
from vcm.clinical.risk_prediction import (
    RiskPredictionModule,
    ClinicalCovariates,
    VirtualCellFeatures
)

module = RiskPredictionModule()

# Create clinical covariates
clinical_data = ClinicalCovariates(
    age=55.0,
    sex='male',
    prior_transplant=True,
    hla_mismatch=4,
    tacrolimus_use=True,
)

# Run risk prediction
result = module.predict_risk(clinical_data, virtual_cell_data)
```

### Configuration API

#### YAML Configuration Format

```yaml
experiment_id: bkpyv_infection
plugin: transplant.bk_polyomavirus
simulator: bkpyv_specific
simulation_length: 100.0
timestep: 1.0
output_path: outputs/bkpyv/

simulator_parameters:
  tacrolimus_enhancement_factor: 2.0
  sirolimus_inhibition_factor: 0.5
  t_antigen_replication_threshold: 0.5
  cell_cycle_s_phase_bonus: 2.0
  dna_replication_coupling: 0.8

perturbations:
  - id: bkpyv_infection
    name: "BKPyV infection"
    perturbation_type: viral_infection
    magnitude: 1.0
    timing: 10.0
```

---

## Configuration Guide

### Parameter Categories

#### Drug Effect Parameters (HIGH Confidence)

These are based on direct experimental measurements or clinical data:

```yaml
# legacy discrete-simulator parameters (bkpyv_specific, hours/step) - the ODE
# engine (bkpyv_ode, days) uses the S-phase gate + graded early-window instead
simulator_parameters:
  tacrolimus_enhancement_factor: 2.0    # literature OR ≈ 2.3 (Fang 2022, PMC9428263)
  sirolimus_inhibition_factor: 0.5      # in-vitro IC90 ≈ 4 ng/mL (Hirsch 2016)
  drug_effectiveness_window: 24.0        # early-phase window (hours)
  late_phase_drug_resistance: 0.3        # reduced late-phase effectiveness factor
```

#### Viral Replication Parameters (MEDIUM Confidence)

These are based on biological consensus but require calibration:

```yaml
simulator_parameters:
  t_antigen_replication_threshold: 0.5   # T antigen requirement
  dna_replication_coupling: 0.8          # Host-viral coupling
  cell_cycle_s_phase_bonus: 2.0          # S-phase optimal
  viral_replication_rate: 0.1            # Base replication rate
  viral_clearance_rate: 0.1              # Base clearance rate
```

#### Clinical Risk Factors (HIGH Confidence)

These are derived from clinical cohort studies:

```yaml
simulator_parameters:
  age_risk_multiplier: 1.9               # Age >50 years
  male_sex_risk_multiplier: 2.3          # Male sex
  prior_transplant_risk_multiplier: 3.0    # Prior transplant
  hla_mismatch_risk_multiplier: 1.3      # HLA 4-6 mismatch
```

### Scenario Configurations

#### Baseline Uninfected

```yaml
configs/bkpyv_baseline.yaml
# Healthy kidney cell with no infection
```

#### Infection No Drug

```yaml
configs/bkpyv_infection.yaml  
# BKPyV infection without immunosuppression
```

#### Drug Exposures

```yaml
configs/bkpyv_tacrolimus.yaml   # Tacrolimus-based regimen
configs/bkpyv_sirolimus.yaml     # Sirolimus-based regimen
```

#### Risk Profiles

```yaml
configs/bkpyv_low_replication.yaml   # Low-risk patient scenario
configs/bkpyv_high_replication.yaml  # High-risk patient scenario
```

---

## Development Guide

### Adding a New Disease Module

The BKPyV plugin serves as a template for adding new disease modules:

#### 1. Create Plugin Directory

```bash
mkdir -p src/vcm/plugins/respiratory/rsv
```

#### 2. Implement Plugin Class

```python
# src/vcm/plugins/respiratory/rsv/rsv.py
from vcm.plugins.base import BasePlugin

class RSVPlugin(BasePlugin):
    def __init__(self):
        super().__init__()
        self.plugin_id = "respiratory.rsv"
        self.plugin_name = "RSV Respiratory Model"
    
    def create_initial_state(self, config=None):
        # Define RSV-specific state
        # Genes, proteins, pathways
        pass
    
    def get_supported_perturbations(self):
        # Define RSV-specific perturbations
        pass
```

#### 3. Create Simulator

```python
# src/vcm/simulators/rsv_simulator.py
from vcm.simulators.base import BaseSimulator

class RSVSimulator(BaseSimulator):
    def simulate(self, initial_state, perturbation, n_steps, timestep):
        # Implement RSV-specific dynamics
        pass
```

#### 4. Create Parameter Registry

```python
# src/vcm/plugins/respiratory/rsv/parameters.py
# Document RSV-specific parameters with evidence sources
```

#### 5. Register Plugin

Add to plugin registry in `src/vcm/plugins/__init__.py`

### Running Tests

```bash
# Run all tests
pytest tests/

# Run specific test
pytest tests/test_bkpyv_plugin.py

# Run with coverage
pytest --cov=vcm tests/
```

### Code Style

```bash
# Format code
black src/

# Check linting
flake8 src/

# Type checking
mypy src/
```

### Contributing Guidelines

1. **Fork the repository** and create a feature branch
2. **Write tests** for new functionality
3. **Update documentation** for any API changes
4. **Submit pull request** with clear description of changes

---

## Troubleshooting

### Common Issues

#### Import Errors

**Problem**: `ModuleNotFoundError: No module named 'vcm'`

**Solution**: Ensure you're in the project root and have activated the virtual environment

```bash
cd /path/to/virtual-cell-model
source venv/bin/activate  # or venv\Scripts\activate on Windows
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

#### Streamlit Won't Start

**Problem**: Streamlit fails to start or displays errors

**Solution**: Ensure all dependencies are installed

```bash
pip install streamlit plotly pandas numpy pyyaml
```

#### Simulation Results Unexpected

**Problem**: Simulation produces unrealistic viral load values

**Solution**: 
- Check parameter values are in valid ranges
- Verify configuration file format is correct
- Review parameter registry for expected value ranges
- Run validation script: `python validate_bkpyv.py`

#### Data Files Not Found

**Problem**: "Configuration file not found" or "Data file not found"

**Solution**: 
- Verify you're running from project root directory
- Check file paths are relative to project root
- Use absolute paths if relative paths aren't working

### Debug Mode

Enable detailed logging for troubleshooting:

```python
import logging

logging.basicConfig(level=logging.DEBUG)
```

### Performance Issues

**Problem**: Simulations run slowly

**Solution**:
- Reduce simulation length or timestep
- Use simpler pathway models
- Disable unnecessary logging
- Consider using faster hardware for large-scale simulations

---

## Additional Resources

### Documentation Files

- **`docs/judges_mentors_guide.md`**: Guide for ISEF judges and mentors
- **`docs/bkpyv_research_to_model_map.md`**: Detailed research-to-model mapping
- **`docs/clinical_integration_plan.md`**: Clinical integration roadmap
- **`data/README.md`**: Data sources and processing methods

### Key Source Files

- **`src/vcm/plugins/transplant/bk_polyomavirus/parameters.py`**: Parameter registry
- **`src/vcm/clinical/viral_load_mapper.py`**: Clinical load conversion
- **`src/vcm/clinical/risk_prediction.py`**: Risk prediction module
- **`src/vcm/ui/streamlit_app.py`**: Web interface

### Community and Support

- **GitHub Issues**: Report bugs and request features
- **Documentation**: See `docs/` folder for detailed guides
- **Examples**: Check `notebooks/` for tutorial notebooks

---

## Citation

If you use this model in your research, please cite:

```
BKPyV Virtual Cell Model - A Mechanistic Simulation of BK Polyomavirus
Replication in Kidney Transplant Recipients

Version 1.0
https://github.com/yourusername/virtual-cell-model

Key References:
- Hirsch HH et al., Am J Transplant 2016
- Weissbach FH et al., J Virol 2024
- KDIGO Clinical Practice Guidelines 2023
- AST Consensus Statement 2019
```