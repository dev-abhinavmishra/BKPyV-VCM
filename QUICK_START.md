# Quick Start Guide - Virtual Cell Model

> Historical quick-start snapshot. For the current source of truth, use
> [README.md](README.md), especially the reviewer-bundle and Streamlit sections.
> Test counts are intentionally not hard-coded because the suite grows with the
> review tooling.

---

## Current Project Location
Run commands from the repository root (the directory containing `README.md`).

## âœ… Verified Working Commands

### 1. Run Tests
```bash
python -m pytest tests/ -v
```

### 2. List Available Plugins (5 plugins)
```bash
python -m vcm list-plugins
```

### 3. Inspect BKPyV Plugin (Flagship)
```bash
python -m vcm inspect --plugin transplant.bk_polyomavirus
```

### 4. Run Baseline Simulation
```bash
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_baseline.yaml
```

### 5. Run Advanced Drug Taper Scenario
```bash
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_drug_taper.yaml
```

### 6. Compare Drug Scenarios
```bash
python -m vcm compare --plugin transplant.bk_polyomavirus --scenario1 configs/bkpyv_baseline.yaml --scenario2 configs/bkpyv_tacrolimus.yaml
```

### 7. Launch Streamlit Dashboard
```bash
python -m streamlit run src/vcm/ui/streamlit_app.py
```
**Note:** This will launch the dashboard at http://localhost:8501
Press Ctrl+C to stop

---

## âŒ Common Mistakes to Avoid

### DON'T use inline comments in commands
**Wrong:**
```bash
python -m pytest tests/ -v  # Run all tests  â† This breaks the command
```

**Right:**
```bash
python -m pytest tests/ -v
```

### Run the dashboard through Python's module syntax
```bash
python -m streamlit run src/vcm/ui/streamlit_app.py
```
(This requires the optional UI dependencies: `pip install -e ".[ui]"`)

### DON'T forget to navigate to project first
```bash
cd <repo-root>
```

---

## ðŸŽ¯ Quick Test Sequence

```bash
# Navigate to project
cd <repo-root>

# Run tests (should take ~25 seconds)
python -m pytest tests/ -v

# List plugins (should show 5 plugins)
python -m vcm list-plugins

# Run a simulation (should complete in ~2 seconds)
python -m vcm run --plugin transplant.bk_polyomavirus --config configs/bkpyv_baseline.yaml

# Launch dashboard (opens browser at localhost:8501)
python -m streamlit run src/vcm/ui/streamlit_app.py
```

---

## ðŸ“Š Expected Results

### Tests
The command should finish with zero failures; the number of discovered tests is
environment- and revision-dependent.

### Plugin List
```
Available Plugins:
â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”³â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”“
â”ƒ Plugin ID               â”ƒ Plugin Name              â”ƒ
â”¡â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â•‡â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â•‡
â”‚ bacteria.minimal_cell   â”‚ Minimal Bacterial Cell   â”‚
â”‚ virus_host.sars_cov2    â”‚ SARS-CoV-2 Virus-Host    â”‚
â”‚ mammalian.immune_tcell  â”‚ Immune T-Cell            â”‚
â”‚ mammalian.cancer_cell   â”‚ Cancer Cell              â”‚
â”‚ transplant.bk_polyomavâ€¦ â”‚ BK Polyomavirus Kidney   â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”´â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
```

### Dashboard
```
You can now view your Streamlit app in your browser.
Local URL: http://localhost:8501
```

---

## ðŸ”§ Troubleshooting

### If tests fail with "ModuleNotFoundError"
```bash
pip3 install -e .
```

### If streamlit command not found
```bash
pip install -e ".[ui]"  # installs the optional UI dependencies (streamlit)
```

### If CLI command fails
```bash
# Make sure you're in the correct directory
pwd
# Should show the repository root containing README.md
```

---

## ðŸ“š Full Documentation

For detailed instructions, see:
- **`README.md`** - Canonical current status, commands, and review outputs

---

Run the commands from the repository root after installing the project dependencies.
