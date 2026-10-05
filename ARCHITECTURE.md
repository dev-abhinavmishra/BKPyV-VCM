# VCM Architecture Overview

> This document describes the platform architecture. For the current BKPyV
> implementation status, runnable commands, and reviewer outputs, see
> [README.md](README.md).

## Design Philosophy

The Virtual Cell Model (VCM) is designed with the following principles:

1. **Modularity**: Components are loosely coupled through well-defined interfaces
2. **Extensibility**: New biological systems can be added as plugins without modifying core code
3. **Reproducibility**: Experiments are fully configurable and version-controlled
4. **Simplicity**: The codebase is understandable by strong student researchers
5. **Scientific Rigor**: Type hints, validation, and testing ensure reliability

## Component Architecture

### Layered Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     CLI & API Layer                          │
│  (User interface: commands, Python API, visualization)       │
└────────────────────────┬────────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────────┐
│                  Experiment Layer                             │
│  (ExperimentRunner, ExperimentComparator, Config management)  │
└────────────────────────┬────────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────────┐
│                  Simulation Layer                            │
│  (BaseSimulator, MechanisticSimulator, MLBasedSimulator,    │
│   HybridSimulator)                                           │
└────────────────────────┬────────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────────┐
│                   Plugin Layer                               │
│  (BasePlugin, PluginRegistry, biological system plugins)     │
└────────────────────────┬────────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────────┐
│                    Core Layer                                │
│  (Domain models: CellState, Gene, Protein, etc.)             │
└────────────────────────┬────────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────────┐
│                    Data Layer                                │
│  (DataLoader, DataValidator, file I/O)                       │
└─────────────────────────────────────────────────────────────┘
```

## Key Design Decisions

### 1. Plugin Architecture

**Decision**: Use a plugin system instead of hardcoding biological models

**Rationale**:
- Allows researchers to add new systems without modifying core code
- Enables sharing of models between projects
- Facilitates testing and validation of individual models
- Makes the codebase more maintainable

**Implementation**:
- `BasePlugin` defines the interface
- `PluginRegistry` manages plugin lifecycle
- Each plugin provides: default config, initial state, supported perturbations

### 2. Simulator Interface

**Decision**: Define a common interface for all simulators

**Rationale**:
- Allows swapping simulation engines without changing experiment code
- Enables comparison of different approaches
- Facilitates adding new simulator types (e.g., stochastic, agent-based)

**Implementation**:
- `BaseSimulator` abstract class with `simulate()` and `step()` methods
- Three implementations: Mechanistic, ML-based, Hybrid
- Consistent input/output: CellState → SimulationResult

### 3. Pydantic Models

**Decision**: Use Pydantic for data validation

**Rationale**:
- Automatic validation of biological entities
- Clear schema definitions
- JSON serialization/deserialization built-in
- Type safety with runtime validation

**Implementation**:
- All core models inherit from Pydantic `BaseModel` or use `@dataclass`
- Field constraints ensure biologically plausible values
- Validation errors are clear and actionable

### 4. Configuration Files

**Decision**: Use YAML for experiment configurations

**Rationale**:
- Human-readable and editable
- Supports comments and complex structures
- Version-control friendly
- Easy to share and modify

**Implementation**:
- `ExperimentConfig` model for validation
- YAML loading/saving in `ExperimentRunner`
- Example configs provided for each plugin

### 5. Vector Representation

**Decision**: Include both named entities and vector representation

**Rationale**:
- Named entities are interpretable and biologically meaningful
- Vector representation enables ML approaches
- Hybrid approach supports both mechanistic and ML methods

**Implementation**:
- `CellState.get_state_vector()` concatenates all quantities
- Named entities (genes, proteins, etc.) remain accessible
- Automatic synchronization between representations

## Data Flow

### Experiment Execution Flow

```
1. Load Config (YAML)
   ↓
2. Load Plugin
   ↓
3. Create Initial State
   ↓
4. Initialize Simulator
   ↓
5. For each step:
   - Apply perturbations (if timing matches)
   - Update state using simulator
   - Store step in result
   ↓
6. Return SimulationResult
   ↓
7. Save results (JSON + CSV)
   ↓
8. Generate visualizations (optional)
```

### Comparison Flow

```
1. Load two configs
   ↓
2. Run both experiments
   ↓
3. Add results to comparator
   ↓
4. Compare final states or trajectories
   ↓
5. Generate comparison plots
   ↓
6. Output comparison data
```

## Extension Points

### Adding a New Simulator

1. Inherit from `BaseSimulator`
2. Implement `simulate()` and `step()` methods
3. Register in `ExperimentRunner.simulator_registry`
4. Add tests

### Adding a New Perturbation Type

1. Add enum value to `PerturbationType`
2. Implement effect in simulator(s)
3. Add to plugin's supported perturbations
4. Update documentation

### Adding Visualization Types

1. Add function in `src/vcm/viz/plots.py`
2. Follow existing naming convention
3. Add to `__init__.py`
4. Update CLI or Python API if needed

### Adding Data Formats

1. Add method to `DataLoader`
2. Implement validation in `DataValidator`
3. Add tests
4. Update documentation

## Testing Strategy

### Unit Tests
- Test individual components in isolation
- Mock dependencies where appropriate
- Cover edge cases and error conditions

### Integration Tests
- Test component interactions
- Use mock data for consistency
- Verify end-to-end workflows

### Test Coverage
- Core models: 100% coverage
- Simulators: High coverage (update rules)
- Plugins: Test state creation and perturbations
- CLI: Test command parsing and execution

## Performance Considerations

### Current State
- The BKPyV path includes both the pathway simulator and a reviewable adaptive ODE engine
- The reviewer bundle records solver tolerances and generates source tables beside figures
- Simulations remain in-memory unless an export path is explicitly requested

### Future Optimizations
- Lazy loading of large datasets
- Parallel execution of independent simulations
- Caching of simulation results
- Efficient data structures for large networks

## Security Considerations

### Current State
- No network access required
- No external API calls
- All data is local

### Future Considerations
- Validate user inputs in CLI
- Sanitize file paths
- Secure handling of sensitive biological data
- Rate limiting for external APIs (if added)

## Documentation Strategy

### Code Documentation
- Docstrings for all public functions/classes
- Type hints for all function signatures
- Inline comments for complex logic

### User Documentation
- README with quick start
- Architecture overview (this document)
- Plugin development guide
- API reference (generated from docstrings)

### Example Documentation
- Example configs for each plugin
- Jupyter notebooks demonstrating workflows
- Tutorial scripts for common tasks

## Post-audit invariants (2026-09)

The following invariants are enforced by tests and must hold in any change:

1. **`bkpyv_ode` is registered** in `ExperimentRunner.simulator_registry` and
   is the canonical BKPyV dynamics engine. The discrete `bkpyv_specific`
   simulator is retained for backward compatibility only.
2. **Time units are explicit.** ODE rate constants are per **day**; configs
   declare `time_unit:` (days for ODE, hours for the legacy discrete engine).
   The continuous drug dosing uses `timing`/`duration` in the config's unit.
3. **Config parameters actually reach simulators.** `simulator_parameters:`
   (and the legacy alias `simulator_params:`) are merged and passed to the
   simulator constructor; YAML `1e-6`-style floats must be written as `1.0e-6`
   (PyYAML parses the former as strings; the simulator also coerces).
4. **No free-floating success banners.** Scripts under repo root
   (`validate_*.py`) exit nonzero on failure and never print unconditional
   success text.
5. **The copies/mL bridge is assumption-labelled.** All clinical-facing
   outputs carry explicit caveat text.

## Error Handling

### Validation Errors
- Pydantic validation for model inputs
- Clear error messages with context
- Suggestions for fixing common issues

### Runtime Errors
- Graceful degradation where possible
- Informative error messages
- Logging of failures for debugging

### File I/O Errors
- Check file existence before operations
- Handle permission errors
- Provide clear paths in error messages

## Dependencies

### Core Dependencies
- `pydantic`: Data validation
- `pyyaml`: Config file parsing
- `numpy`: Numerical operations
- `pandas`: Data manipulation

### Optional Dependencies
- `matplotlib`: Static plotting
- `plotly`: Interactive plotting
- `click`: CLI framework
- `rich`: CLI formatting

### Development Dependencies
- `pytest`: Testing
- `black`: Code formatting
- `ruff`: Linting
- `mypy`: Type checking

## Future Architecture Enhancements

### Potential Improvements
1. **Database Integration**: Store results in SQLite/PostgreSQL
2. **Web Interface**: FastAPI/Streamlit frontend
3. **Distributed Computing**: Run simulations on clusters
4. **Model Zoo**: Central repository for shared plugins
5. **Version Control**: Track model versions and experiments
6. **CI/CD**: Automated testing and deployment

### Architectural Debt
- Add logging framework
- Improve error handling consistency
- Add configuration validation before execution
- Standardize file naming conventions
- Add progress bars for long simulations

## Conclusion

This architecture provides a solid foundation for computational biology research while remaining simple enough for student researchers to understand and extend. The modular design enables incremental addition of real biological data and sophisticated simulation methods as the project evolves.
