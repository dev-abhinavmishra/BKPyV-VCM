> **NOTE (post-audit):** Sections describing parameter values, calibrations, or
> numbers may be historical. The canonical, current artifacts are:
> README.md, docs/BKPYV_MODEL_CARD.md, docs/ISEF_PROJECT_OVERVIEW.md, and
> the code itself (src/vcm/simulators/ode_system.py is the source of truth
> for the ODE model and its parameters). Time unit: days for the ODE engine.

# ODE-Based BKPyV Simulation Implementation

## Overview

This document describes the ODE (Ordinary Differential Equation) based implementation of BK polyomavirus (BKPyV) simulation in the Virtual Cell Model (VCM) platform. The ODE system provides continuous-time dynamics using scipy's numerical integration, replacing the discrete-time multiplicative update rules with proper differential equations.

## Architecture

### Components

1. **ODE System Definition** (`src/vcm/simulators/ode_system.py`)
   - Defines the 23-dimensional state vector and differential equations
   - Implements research-validated parameters
   - Provides initial conditions and state variable management

2. **ODE-Based Simulator** (`src/vcm/simulators/bkpyv_ode_simulator.py`)
   - Integrates with existing VCM infrastructure via BaseSimulator interface
   - Converts between CellState objects and ODE state vectors
   - Uses scipy.integrate.solve_ivp for numerical integration
   - Maps research parameters to ODE coefficients

3. **Configuration Files** (`configs/bkpyv_ode_*.yaml`)
   - ODE-specific configurations with solver settings
   - Research-validated parameter mappings
   - Drug scenario configurations (tacrolimus, sirolimus)

## State Variables

The ODE system uses 23 state variables (15 core + 8 appended extensions;
indices never renumbered — new states append only):

| Index | Variable | Description | Biological Meaning |
|-------|----------|-------------|-------------------|
| 0 | V | Viral load | Concentration of BKPyV particles |
| 1 | T | T antigen | Large T antigen concentration (drives replication) |
| 2 | G_v | Viral gene expression | Overall viral gene expression level |
| 3 | C | Healthy target cells | Uninfected kidney tubular epithelial cells |
| 4 | I | Infected cells | BKPyV-infected kidney cells |
| 5 | D | Dead/damaged cells | Cells that have died from infection or other causes |
| 6 | CC | Cell cycle phase | Continuous 0-1: G0/G1 (0) â†’ S (0.5) â†’ G2/M (1.0) |
| 7 | DNA | DNA synthesis activity | Host DNA replication machinery activity |
| 8 | E | Effector immune cells | Immune cells fighting infection |
| 9 | IFN | Interferon concentration | Antiviral interferon signaling molecules |
| 10 | AK | Antiviral state | Cellular antiviral response state |
| 11 | D_tac | Tacrolimus concentration | Tacrolimus drug concentration |
| 12 | D_sir | Sirolimus concentration | Sirolimus drug concentration |
| 13 | P_rep | DNA replication pathway activity | DNA replication pathway flux |
| 14 | P_immune | Innate immune pathway activity | Innate immune signaling activity |
| 15 | T_naive | BKPyV-specific naive T cells | Tacrolimus-inhibited priming (NFAT) |
| 16 | T_eff | BKPyV-specific effector T cells | Antigen-driven expansion and infected-cell kill |
| 17 | C_u | Healthy urothelial cells | Urinary reservoir (Funk 2008) |
| 18 | I_u | Infected urothelial cells | >95% of urinary load is urothelial |
| 19 | V_u | Urinary virion pool | urine:plasma ~3000:1 |
| 20 | F_rr | rr-NCCR fraction, kidney pool | In-host quasi-species emergence (Gosert 2008) |
| 21 | F_rr_u | rr-NCCR fraction, urinary pool | Weakened selection + drainage mixing |
| 22 | L | Latently-infected reservoir | Seeded by latent_fraction of new infections; reactivates under immunosuppression (1 - tac_immune_effect) |

## Differential Equations (v2 form, matching `ode_system.py`)

### Viral dynamics

```
dV/dt = p * I * (T / (T + K_T)) * sir_effect * virion_cost - delta * V * (1 + a_max * AK/(AK + a_half))

  where virion_cost = 1 - (1 - nccr_capsid_eff) * rr_capsid_virion_cost
  (the rearranged pool's capsid deficit partially reduces progeny yield —
  expression deficit != yield deficit 1:1)
```

- `sir_effect` saturates with sirolimus dose and is stronger before T crosses `t_threshold`.
- Clearance enhancement by antiviral state saturates (`a_max=2`, `a_half=1`).

### T antigen dynamics (S-phase gate)

```
dT/dt = t_prod * G_v * P_rep * dna_coupling * H(CC; cc_s_phase) * H(DNA; 0.5) - t_decay * T
```

with Hill gates H of order 4. T antigen accumulates only when host cell-cycle
permissiveness AND DNA synthesis are present (Needham 2024).

### Cell populations

```
infection_rate = (beta * (V + cross_feed * V_u) + c2c_rate * I) * C
  — two transmission modes: free virions (V-dependent, cleared by
  delta/AK; urine reseeds via cross_feed) and cell-to-cell spread
  via virological synapses (V-independent, clearance-insensitive:
  the persistence channel)
dC/dt = lambda_cell - infection_rate - d_cell * C               (constant source; NOT logistic)
dI/dt = infection_rate - d_infected * I - immune_kill * E * I * tac_immune_effect
dD/dt = cell_death + infected_death - 0.1 * D
```

### Cell-cycle permissiveness (bounded, no ratchet)

```
dCC/dt = cc_rate * ( (P_rep * mtor_rep_effect) / (1 + P_rep * mtor_rep_effect) - CC )
```

### Drug administration (dosing intensity, dimensionless)

```
dD/dt = k_in * (target - D)   while treatment is active (timing <= t < timing+duration)
dD/dt = -clearance * D        otherwise
```

### Pathway dynamics

```
dP_rep/dt    = p_rep_prod * (DNA/(1+DNA)) * mtor_rep_effect - p_rep_decay * P_rep
dP_immune/dt = p_immune_prod * IFN * tac_immune_effect - p_immune_decay * P_immune
dE/dt        = e_prod * P_immune/(1+E) - e_decay * E    (self-limiting expansion)
```

## Key Parameters (v2 defaults; see `_get_default_params` for the full dict)

| Parameter | Value | Grounding | Description |
|-----------|-------|-----------|-------------|
| beta | 0.3 | tuned (classification regime) | Free-virion infection rate, 1/day |
| c2c_rate | 0.03 | tuned (persistence channel) | Cell-to-cell spread rate, 1/day |
| delta | 0.4 | Funk 2006 (slow phase) | Viral clearance rate, 1/day |
| p | 8.0 | tuned | Virion production per infected cell, 1/day |
| immune_kill | 0.8 | tuned | Effector-cell killing of infected cells |
| tac window | - | Hirsch 2016 (in vitro) | Tacrolimus: `tac_immune_effect = 1/(1+2*S*D_tac)` |
| mtor_inhibition | 0.5 | Hirsch 2016 IC90 â‰ˆ 4 ng/mL | Residual permissiveness under full sirolimus |
| t_threshold | 0.5 | abstraction | T-antigen threshold separating early/late |
| s_phase_bonus | 2.0 | phenomenological | Extra DNA synthesis while S-gate open |
| dna_coupling | 0.8 | phenomenological | Host-DNA/LT coupling strength |
| translation_enhancement | 2.0 | Weissbach 2024 (direction only) | Viral translation elevation |

**Explicitly phenomenological**: all of the above are model coefficients on a
normalized scale; none are fitted concentrations or patient-measured rates.

## Drug Effects

### Tacrolimus
- Modeled as **weakened immune control** (saturating effect on infected-cell
  killing and immune-pathway production), not as a direct genome-copying
  multiplier. Rationale: the documented in-vitro effect (Hirsch 2016) plus the
  clinical association of tacrolimus-based regimens with higher BKPyV risk.
- Continuous dosing from `timing` (optional `duration`).

### Sirolimus
- **Mechanism**: mTOR inhibition -> reduced DNA-pathway activity and viral production.
- **Timing**: around the in-vitro observation that inhibition is strongest
  before productive replication is established, the effect smooth-decays after
  T crosses `t_threshold` (residual fraction `sir_late_weight`).
- Continuous dosing from `timing` (optional `duration`).

## Numerical Integration

### Solver Configuration
- **Default solver**: LSODA (adaptive stiff/non-stiff solver)
- **Alternative solvers**: RK45 (explicit Runge-Kutta), BDF (implicit for stiff systems)
- **Relative tolerance**: 1e-6
- **Absolute tolerance**: 1e-8
- **Maximum step**: 1.0 time units

### Solver Selection Guidelines
- **LSODA**: Best for general use, handles both stiff and non-stiff equations
- **RK45**: Good for non-stiff systems, faster when appropriate
- **BDF**: Best for stiff systems (e.g., rapid drug kinetics)

## Usage Examples

### Python API

```python
from vcm.plugins.transplant.bk_polyomavirus import BKPolyomavirusPlugin
from vcm.simulators.bkpyv_ode_simulator import BKPyVODESimulator
from vcm.core.models import Perturbation, PerturbationType

# Load plugin
plugin = BKPolyomavirusPlugin()
initial_state = plugin.create_initial_state()

# Create ODE simulator with custom configuration
config = {
    'ode_solver': 'LSODA',
    'rtol': 1e-6,
    'atol': 1e-8,
    'tacrolimus_enhancement_factor': 1.5,
    'mtor_inhibition_factor': 0.5
}
simulator = BKPyVODESimulator(config)

# Create infection perturbation
infection = Perturbation(
    id="bkpyv_infection",
    name="BKPyV infection",
    perturbation_type=PerturbationType.VIRAL_INFECTION,
    magnitude=1.0,
    timing=10.0
)

# Run simulation
result = simulator.simulate(
    initial_state=initial_state,
    perturbations=[infection],
    n_steps=100,
    timestep=1.0
)
```

### CLI Usage

```bash
# Run ODE-based baseline simulation
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_ode_baseline.yaml

# Run ODE-based infection simulation
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_ode_infection.yaml

# Run ODE-based tacrolimus simulation
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_ode_tacrolimus.yaml

# Run ODE-based sirolimus simulation
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_ode_sirolimus.yaml
```

## Validation

### Test Coverage
- 28 comprehensive tests in `tests/test_bkpyv_ode.py`
- Tests cover ODE system, simulator integration, parameter mapping
- Validation against research expectations (viral half-lives, drug effects)
- Numerical stability tests with multiple solvers

### Expected Behaviors
1. **Viral clearance half-lives**: 1-38 hours (matches research: 1-2h fast, 20-38h moderate)
2. **Drug effects**: Tacrolimus increases production (weaker immune targeting), sirolimus reduces production (S-phase gate)
3. **T antigen threshold**: Replication requires T > 0.5
4. **State conservation**: Cell populations remain non-negative
5. **Numerical stability**: All solvers produce finite results

## Comparison with Discrete Simulator

| Aspect | Discrete Simulator | ODE Simulator |
|--------|------------------|---------------|
| Time stepping | Fixed timestep | Adaptive timestep |
| Mathematical basis | Multiplicative factors | Differential equations |
| Numerical method | Simple Euler integration | scipy.solve_ivp (LSODA/RK45/BDF) |
| Biological accuracy | Qualitative | Quantitative |
| Computational cost | Lower | Higher (but adaptive) |
| Parameter mapping | Direct | Requires coefficient mapping |
| Validation status | Validated | Validated against discrete results |

## Performance Considerations

### Computational Cost
- ODE solving is computationally more expensive than discrete updates
- Adaptive timesteps provide efficiency by using larger steps when dynamics are smooth
- LSODA automatically switches between stiff and non-stiff methods

### Optimization Tips
1. Use appropriate solver for your system (LSODA for general use)
2. Adjust tolerances based on required accuracy (looser tolerances = faster)
3. Limit simulation duration for exploratory runs
4. Use max_step parameter to prevent overly large timesteps

## Future Extensions

### Potential Enhancements
1. **Spatial ODEs**: Add spatial diffusion terms for tissue-level modeling
2. **Stochastic ODEs**: Add noise terms for biological variability
3. **Parameter Estimation**: Fit ODE parameters to clinical data
4. **Multi-scale Modeling**: Combine with intracellular signaling models
5. **Hybrid Approach**: Use ODEs for some variables, discrete for others

### Additional Drug Models
- Everolimus (alternative mTOR inhibitor)
- Belatacept (co-stimulation blocker)
- Corticosteroids (broad immunosuppression)
- Antiviral therapies (cidofovir, leflunomide)

## Troubleshooting

### Common Issues

**Problem**: ODE solver fails with "excess work done"
- **Solution**: Reduce max_step parameter or increase tolerances

**Problem**: Unstable solutions (oscillations, divergence)
- **Solution**: Try different solver (BDF for stiff systems), reduce timestep

**Problem**: Results differ from discrete simulator
- **Solution**: Check parameter mapping, adjust ODE coefficients to match discrete behavior

**Problem**: Negative cell populations
- **Solution**: Add clipping in ODE function or adjust production/death rates

## References

1. **AJT-16-821.pdf**: "BK Polyomavirus Replication in Renal Tubular Epithelial Cells Is Inhibited by Sirolimus, but Activated by Tacrolimus Through a Pathway Involving FKBP-12"
2. **jvi.01382-24-s0003.pdf**: Single-cell transcriptomic analysis of BKPyV infection
3. **jvi.01382-24-s0004.pdf**: Pathway activity profiling in BKPyV infection
4. **abstract-16323135.txt**: Viral clearance kinetics (1-2h fast, 20-38h moderate half-lives)
5. **scipy.integrate.solve_ivp documentation**: Numerical ODE solver methods

## Contributing

When extending the ODE system:
1. Update state vector size and names
2. Add corresponding differential equations
3. Update CellState â†” ODE conversion methods
4. Add tests for new functionality
5. Update this documentation
6. Validate against research expectations

## License

This ODE implementation is part of the Virtual Cell Model project and follows the same MIT license.