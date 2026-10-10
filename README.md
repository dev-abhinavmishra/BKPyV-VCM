# Virtual Cell Model (VCM)

A modular, extensible platform for cellular simulation designed for ISEF research projects. VCM provides a framework for mechanistic and ML-based cell modeling, with support for various biological systems including minimal bacterial cells, virus-host interactions, immune cells, and cancer cells.

## What is VCM?

The Virtual Cell Model (VCM) is a Python-based computational biology platform that enables researchers to:

- Simulate cellular behavior using mechanistic, ML-based, or hybrid approaches
- Model different biological systems through a plugin architecture
- Apply and compare perturbations (drug treatments, genetic modifications, environmental changes)
- Track and visualize simulation results
- Build reproducible experiments with configurable parameters

VCM is designed as a foundation for computational biology research, particularly suitable for high-school science fair projects (ISEF) that require a technically rigorous, extensible codebase.

### BKPyV research mode

The BKPyV model includes an evidence-aware renal-cell mode shaped by primary-cell single-cell biology:

- host S phase/DNA synthesis gates large T-antigen accumulation (Needham et al., PLoS Pathog 2024);
- intracellular replication flux is kept separate from system-level viral production;
- tacrolimus is modeled as reduced immune control, not direct genome-copying enhancement (mechanism motif from Hirsch et al., AJT 2016, in vitro);
- archetype and rearranged NCCR scenarios expose a testable early-gene/capsid trade-off.

These are mechanistic hypotheses and visualization outputs, not patient-specific clinical predictions. See [the BKPyV model card](docs/BKPYV_MODEL_CARD.md) for interpretation boundaries and proposed experiments.

**Time units**: the ODE simulator (`bkpyv_ode`) uses **days** with per-day rate constants. The legacy discrete simulator (`bkpyv_specific`) uses **hours** per step and is superseded for analysis. `configs/*.yaml` declare `time_unit` explicitly.

## Features

### Core Capabilities
- **Modular Architecture**: Plugin-based system for different biological models
- **Multiple Simulation Engines**: Mechanistic, ML-based, and hybrid simulators
- **Experiment Management**: Run, compare, and track simulation experiments
- **Data Integration**: Support for JSON/CSV data loading and validation
- **Visualization**: Time-series plots, state comparisons, and pathway networks
- **CLI Interface**: Command-line tools for running experiments

### Included Plugins
1. **Minimal Bacterial Cell** (`bacteria.minimal_cell`): Simplified bacterial model (e.g., Mycoplasma-like)
2. **SARS-CoV-2 Virus-Host** (`virus_host.sars_cov2`): Virus-host interaction model
3. **Immune T-Cell** (`mammalian.immune_tcell`): T-cell activation and checkpoint modeling
4. **Cancer Cell** (`mammalian.cancer_cell`): Cancer cell model with oncogenes and tumor suppressors
5. **BK Polyomavirus Kidney Cell** (`transplant.bk_polyomavirus`): BKPyV infection in renal tubular epithelial cells with drug effects

## Installation

### Prerequisites
- Python 3.10 or higher
- pip or uv (recommended)

### Setup

1. Clone the repository:
```bash
git clone <repository-url>
cd virtual-cell-model
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -e .
# optional: dashboard UI
pip install -e ".[ui]"
# optional: single-cell analysis pipeline (scanpy/anndata/umap-learn/pyarrow)
pip install -e ".[analysis]"
# or with uv
uv pip install -e .
```

4. Verify installation:
```bash
python -m vcm list-plugins
```

### Dashboard (one command)

```bash
streamlit run src/vcm/ui/streamlit_app.py
```

Pages: model overview (Home), scenario simulator (choose drug, dose,
start/stop day), single-cell results (GSE317012), viral-load validation
(benchmark vs published kinetics), risk prediction, comparison, a
parameter/assumption table with citations, and a review-bundle export.

### Reproducing the analyses

All pipeline outputs under `data/` and `outputs/` are generated artifacts
(gitignored) — regenerate them with the scripts below. A clean clone plus
these commands reproduces every figure and table:

```bash
# 1. GSE317012 biopsy scRNA-seq (1.9 GB download, sha256-verified, no manual steps)
python scripts/download_gse317012.py          # fetch + extract + verify
python scripts/analyze_gse317012.py           # biopsy-level pathway analysis
python scripts/cluster_gse317012.py           # cell-level QC/clustering/DE/holdout

# 2. Plasma viral-load benchmark (published kinetics, ~1 min)
python scripts/benchmark_viral_load.py

# 3. Tests
python -m pytest tests/ -q
```

Fixed seeds make outputs deterministic. See `docs/METHODS.md` for the
validation design and `docs/STATE.md` for the audit log.

## Quick Start

### Running Your First Simulation

1. List available plugins:
```bash
python -m vcm list-plugins
```

2. Inspect a plugin:
```bash
python -m vcm inspect --plugin bacteria.minimal_cell
```

3. Run a baseline simulation:
```bash
python -m vcm run --plugin bacteria.minimal_cell --config configs/minimal_cell_baseline.yaml
```

4. Run a simulation with perturbation:
```bash
python -m vcm run --plugin bacteria.minimal_cell --config configs/minimal_cell_antibiotic.yaml
```

5. Compare two scenarios:
```bash
python -m vcm compare \
  --plugin bacteria.minimal_cell \
  --scenario1 configs/minimal_cell_baseline.yaml \
  --scenario2 configs/minimal_cell_antibiotic.yaml
```

### Python API Usage

```python
from vcm.core.models import CellState, ExperimentConfig, Perturbation, PerturbationType
from vcm.plugins.bacteria.minimal_cell import MinimalCellPlugin
from vcm.simulators.mechanistic import MechanisticSimulator
from vcm.experiments.runner import ExperimentRunner

# Load plugin
plugin = MinimalCellPlugin()
plugin.register()

# Create initial state
initial_state = plugin.create_initial_state()

# Create perturbation
perturbation = Perturbation(
    id="antibiotic_stress",
    name="Antibiotic stress",
    perturbation_type=PerturbationType.ENVIRONMENTAL_CHANGE,
    magnitude=1.0,
    timing=10.0
)

# Configure experiment
config = ExperimentConfig(
    experiment_id="my_experiment",
    plugin="bacteria.minimal_cell",
    simulator="mechanistic",
    simulation_length=50.0,
    timestep=0.5,
    perturbations=[perturbation],
    output_path="outputs/my_experiment/"
)

# Run experiment
runner = ExperimentRunner()
result = runner.run_experiment(config, initial_state)

# Save results
runner.save_result(result, Path(config.output_path))
```

## Architecture

### Core Components

#### 1. Domain Models (`src/vcm/core/models.py`)
- `CellState`: Complete cellular state including genes, proteins, metabolites, pathways
- `Gene`, `Protein`, `Metabolite`, `Pathway`: Biological entity representations
- `Perturbation`: Interventions (knockouts, drugs, environmental changes)
- `Environment`: External conditions (temperature, pH, nutrients)
- `SimulationResult`: Complete simulation trajectory and metadata
- `ExperimentConfig`: Simulation configuration

#### 2. Simulation Engines (`src/vcm/simulators/`)
- `BaseSimulator`: Abstract interface for all simulators
- `MechanisticSimulator`: Deterministic update rules (placeholder for ODE systems)
- `MLBasedSimulator`: Neural network-like state transitions
- `HybridSimulator`: Combines mechanistic and ML approaches

#### 3. Plugin System (`src/vcm/plugins/`)
- `BasePlugin`: Abstract plugin interface
- `PluginRegistry`: Global plugin management
- Plugin implementations for different biological systems

#### 4. Data Layer (`src/vcm/data/`)
- `DataLoader`: Load/save JSON, CSV, and CellState objects
- `DataValidator`: Schema validation for biological entities
- Support for raw, processed, and mock data

#### 5. Experiment Framework (`src/vcm/experiments/`)
- `ExperimentRunner`: Execute and manage experiments
- `ExperimentComparator`: Compare results between scenarios

#### 6. Visualization (`src/vcm/viz/`)
- `plot_trajectory`: Time-series plots for entities
- `plot_state_comparison`: Compare states between simulations
- `plot_pathway_network`: Visualize pathway interactions

#### 7. CLI (`src/vcm/cli/`)
- `run`: Execute simulation experiments
- `compare`: Compare two scenarios
- `list-plugins`: List available plugins
- `inspect`: Inspect plugin details

## Project Structure

```
virtual-cell-model/
├── README.md
├── pyproject.toml
├── requirements.txt
├── configs/                     # Experiment configurations (BKPyV ODE + legacy scenarios)
├── data/                        # Data directory (gitignored)
│   ├── processed/               # Real processed outputs (GSE317012 pathway scores,
│   │                            #   contrasts, sample metadata, microarray inventory,
│   │                            #   synthetic patient cohort + generation params)
│   ├── research/                # Downloaded GEO series matrices + extracted parameters
│   └── mock/                    # Reserved (currently empty - no mock JSONs shipped)
├── notebooks/                   # Jupyter notebooks (executed workflow)
├── scripts/                     # Analysis/regeneration scripts (review bundle,
│                                #   sensitivity, drug-switch, GSE317012 analyzer)
├── src/
│   └── vcm/                     # Main package
│       ├── core/                # Domain models
│       ├── simulators/          # Simulation engines (bkpyv_ode is canonical)
│       ├── plugins/             # Biological system plugins
│       │   ├── bacteria/
│       │   ├── virus_host/
│       │   ├── mammalian/
│       │   └── transplant/      # bk_polyomavirus
│       ├── data/                # Data loading utilities
│       ├── analysis/            # Analysis helpers
│       ├── clinical/            # Assumption-labelled clinical bridge + risk model
│       ├── evaluation/          # ARC evaluation
│       ├── experiments/         # Experiment runner
│       ├── viz/                 # Visualization tools
│       ├── cli/                 # Command-line interface
│       └── ui/                  # Streamlit app
├── tests/                       # Test suite
└── outputs/                     # Simulation results (gitignored)
```

## Adding a New Plugin

1. Create a new plugin class in `src/vcm/plugins/<category>/<plugin_name>.py`:

```python
from vcm.plugins.base import BasePlugin
from vcm.core.models import CellState, ExperimentConfig, Perturbation, Gene

class MyCustomPlugin(BasePlugin):
    def __init__(self):
        super().__init__()
        self.plugin_id = "category.my_custom"
        self.plugin_name = "My Custom Model"
        self.description = "Description of your model"

    def get_default_config(self) -> ExperimentConfig:
        return ExperimentConfig(
            experiment_id="my_custom_default",
            plugin=self.plugin_id,
            simulator="mechanistic",
            simulation_length=100.0,
            timestep=1.0,
            output_path="outputs/my_custom/"
        )

    def create_initial_state(self, config=None) -> CellState:
        # Define your genes, proteins, metabolites, pathways
        genes = {
            "gene1": Gene(id="gene1", name="Gene 1", expression_level=1.0),
            # Add more genes...
        }
        return CellState(
            cell_id="my_cell_001",
            cell_type="my_custom_cell",
            genes=genes,
            proteins={},
            metabolites={},
            pathways={}
        )

    def get_supported_perturbations(self):
        # Define supported perturbations
        return [
            Perturbation(
                id="pert1",
                name="My Perturbation",
                perturbation_type=PerturbationType.GENE_KNOCKOUT,
                target_id="gene1",
                magnitude=1.0
            )
        ]
```

2. Register the plugin in the CLI (`src/vcm/cli/main.py`):
```python
from vcm.plugins.category.my_custom import MyCustomPlugin

def register_default_plugins():
    plugins = [
        # ... existing plugins ...
        MyCustomPlugin(),
    ]
```

3. Create a config file in `configs/`:
```yaml
experiment_id: my_custom_baseline
plugin: category.my_custom
simulator: mechanistic
simulation_length: 100.0
timestep: 1.0
output_path: outputs/my_custom/baseline/
environment:
  temperature: 37.0
  ph: 7.4
perturbations: []
```

4. Test your plugin:
```bash
python -m vcm inspect --plugin category.my_custom
python -m vcm run --plugin category.my_custom --config configs/my_custom_baseline.yaml
```

## Configuration Reference

Experiment configurations use YAML format with the following structure:

```yaml
# Required fields
experiment_id: unique_experiment_name
plugin: plugin_id
simulator: mechanistic|ml_based|hybrid|bkpyv_specific|bkpyv_ode
simulation_length: 100.0
timestep: 1.0
time_unit: days          # canonical unit; ODE simulator = days,
                          # legacy discrete simulator = hours
output_path: outputs/experiment/

# Passed to the simulator constructor (rate constants, solver settings)
simulator_parameters: {}

# Optional fields
save_interval: 1
seed: 42

# Environment settings
environment:
  temperature: 37.0
  ph: 7.4
  oxygen_level: 21.0
  glucose_concentration: 5.0
  nutrient_availability: {}
  stress_factors: {}

# Perturbations
perturbations:
  - id: perturbation_id
    name: "Perturbation name"
    perturbation_type: gene_knockout|protein_inhibition|etc.
    target_id: entity_id
    magnitude: 1.0
    timing: 10.0
    duration: null
    parameters: {}

# Metadata
metadata:
  description: "Experiment description"
  researcher: "Your name"
```

## Roadmap: Next Steps for Real Biology Integration

This v1 foundation provides the architecture. Here's how to extend it with real biological data:

### 1. Adding Real Omics Data
- Integrate gene expression data from RNA-seq
- Load proteomics data for protein concentrations
- Incorporate metabolomics measurements
- Use `DataLoader` to validate schemas against real datasets
- Add preprocessing pipelines for normalization and batch correction

### 2. Viral Sequence Evolution
- Implement viral genome models with mutation rates
- Add quasispecies dynamics
- Incorporate selective pressure from immune response
- Track viral evolution over simulation time
- Add phylogenetic analysis tools

### 3. Mechanistic ODE Modules
- Replace toy update rules with real ODE systems
- Implement rate constants from literature
- Add parameter estimation from experimental data
- Include stochastic simulation algorithms (Gillespie)
- Add sensitivity analysis tools

### 4. Single-Cell Perturbation Modeling
- Integrate single-cell RNA-seq data
- Model cell-to-cell heterogeneity
- Add population-level simulations
- Incorporate lineage tracing
- Add dimensionality reduction (PCA, UMAP) for visualization

### 5. Validation Pipelines
- Compare simulation results to experimental data
- Add statistical validation metrics
- Implement parameter fitting/optimization
- Add cross-validation for ML components
- Create benchmark datasets

### 6. Additional Plugins
- Add more virus-host systems (influenza, HIV, RSV)
- Implement additional transplant-specific models (CMV)
- Add stem cell differentiation models
- Incorporate microbiome interactions
- Add multi-organ/system models

## Development

### Running Tests
```bash
pytest tests/
```

### Generate the reviewer bundle

The reviewer bundle runs the ODE engine for archetype and rearranged NCCR scenarios and exports readable figures plus the underlying tables:

```bash
python scripts/generate_bkpyv_review_bundle.py
```

Outputs are written to `outputs/bkpyv_review/`:

- `mechanistic_trajectory.png` â€” clinical bridge, S-phase/LT gate, and separated replication/production/immune-control traces
- `nccr_comparison.png` â€” archetype versus rearranged NCCR hypothesis comparison
- `drug_mechanisms.png` â€” tacrolimus versus sirolimus separated into immune control, replication, and production
- `sensitivity_heatmap.png` â€” one-at-a-time sensitivity of peak production
- `trajectory_data.csv`, `scenario_summary.csv`, `sensitivity_results.csv` â€” reviewable source tables
- `manifest.json` â€” model horizon and interpretation caveats

### Consistency checks and analyses

```bash
python validate_bkpyv.py       # config harness: ODE configs run + ordering sanity
python validate_isef.py        # mechanism direction checks (drugs, S-phase gate, persistence)
python validate_clinical.py    # regenerates clinical_validation_report.json from the bridge
python scripts/sensitivity_analysis.py   # OAT sweep that actually perturbs the ODE
python scripts/simulate_drug_switch.py   # continuous tac->siro switch scenario
python scripts/compare_patient_series.py # model-vs-patient comparison (runs when a real CSV is given)
python scripts/compare_patient_series.py --make-template   # writes data/research/patient_series_template.csv
```

The Streamlit **Simulation** page also exports the current virtual-patient
trajectory as `bkpyv_virtual_patient_trajectory.csv`.

To launch the interactive interface:

```bash
python -m streamlit run src/vcm/ui/streamlit_app.py
```

### Code Style
The project uses:
- Black for formatting
- Ruff for linting
- mypy for type checking

```bash
black src/vcm/
ruff check src/vcm/
mypy src/vcm/
```

## Contributing

This is an ISEF research project. Contributions should:
1. Follow the existing code style and architecture
2. Include tests for new features
3. Update documentation
4. Use meaningful commit messages

## License

MIT License - See LICENSE file for details

## Acknowledgments

Developed for ISEF research project on computational biology and virtual cell modeling.

## Contact

For questions or suggestions, please open an issue in the repository.
