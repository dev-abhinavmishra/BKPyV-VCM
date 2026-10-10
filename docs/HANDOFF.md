# HANDOFF — where this project stands and what's next

Updated: 2026-10-09. Branch `devin/finish-vcm` (tracks origin; several
commits ahead of main, PR pending).

## What was built (this arc, 4 commits)

1. **21-dimensional ODE** (`11b3161`): added BKPyV-specific T-cell arm
   (T_naive priming + T_eff expansion, tac-attenuated via
   `tac_tcell_suppression=12` — calcineurin hits priming/proliferation,
   NOT the kill itself) and a population-level urothelial compartment
   (C_u/I_u/V_u) cross-feeding the kidney. Reproduces the Funk 2008
   signature: urine:plasma ~384×, plasma clears under 90% curtailment
   while viruria persists.
2. **Real PK dosing** (`800ad67`): `trough_ng_ml` schedule windows with
   actual half-lives (tac 12h, sir 60h) and reference troughs (tac 8,
   sir 4 = Hirsch IC90). Schedule LISTS express tapers; `parameters:
   {"units": "ng_ml"}` on a Perturbation routes through
   `_build_dosing_context`. Legacy dict/target/bolus dosing unchanged.
3. **NCCR emergence** (`8cf9170`): F_rr (index 20) — mutation supply ∝
   replication pressure, logistic rr advantage, slow reversion. The old
   discrete archetype/rearranged presets became F=0/F=1 initial
   conditions; effective multipliers interpolate (early ×1→2, capsid
   ×1→0.5). F_rr≥0.5 by ~wk 14 under sustained viremia (Gosert window).
4. **Schedule optimizer** (`5c5376e`): `scripts/optimize_reduction_schedule.py`
   — taper/conversion sweep scored on clearance week vs T_eff rebound
   AUC (rejection proxy). HEADLINE RESULT: tac taper alone never clears
   within 180d; tac→sir conversion clears ~wk8 with LOWER rebound —
   emergent from mechanism, matching Hirsch 2016 clinical observations.
5. **Identifiability** (`scripts/identifiability_analysis.py`):
   per-observable sensitivity screen — plasma-V sees only the kinetics
   core; urine params need urine data; drug-PD needs drug-perturbation
   data; intracellular gate chain is pairwise-confounded. This is the
   honest answer to "what can patient series fit."

Also: benchmark extended with compartment-signature checks (all PASS,
run `scripts/benchmark_viral_load.py --quick`); 💉 Regimen Design
Streamlit page; 208 tests green.

## Ironclad invariants (do not break)

- Indices APPEND only — never renumber the state vector.
- T-antigen never drives dCC/dt; replication=intracellular copying vs
  production=system virion output; drug wording is production-side.
- copies/mL bridge stays assumption-labelled; "hypothesis generator",
  never "fitted to patients".
- `validate_*.py` exit nonzero on failure; YAML floats `1.0e-6`;
  ODE=days; ~1.6k grandfathered ruff findings — touched files stay clean.
- ISEF **2027**. `/data`,`/outputs` gitignored. No fake/placeholder data.

## Run it

```
.venv/bin/python -m pytest tests/ -x -q          # 208 pass, 8 skip
.venv/bin/python scripts/benchmark_viral_load.py --quick   # OVERALL: PASS
.venv/bin/python scripts/optimize_reduction_schedule.py    # regimen table
.venv/bin/python scripts/identifiability_analysis.py       # sensitivity ranks
streamlit run src/vcm/ui/streamlit_app.py        # 11 pages
```

## What's next (highest value first)

1. **Open the PR** for this whole arc to main — it's the coherent
   multi-compartment story: mechanism → emergence → clinical regimen
   question answered mechanistically.
2. **Patient-series fitting** when Thompson's longitudinal data arrives:
   fit ONLY the plasma-visible parameter subset (per identifiability
   ranks), use urine series for urothelial params. `compare_patient_series.py`
   is the data-gated runner.
3. **Sensitivity figure for the paper**: extend the identifiability
   script into a pub-quality tornado/heatmap figure (F_rr×V phase plot
   is the money visual — emergence under viremia, suppression under
   conversion).
4. **Urine:plasma rr fraction**: single-F simplification can't show
   Gosert's compartment asymmetry; a second F_rr per compartment is the
   honest upgrade if a judge asks.
5. Abstract/poster updates: the tac→sir conversion result is the
   headline; the identifiability table is the rigor signal.
