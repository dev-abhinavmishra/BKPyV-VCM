# BKPyV-VCM — Everything, in one file

*The complete nitty-gritty reference for this project: what it is, how it
works, every mechanism, every result, every number, every file. Written to
be read top-to-bottom by someone who has never seen the repo — and to be
the fastest possible on-ramp for a fresh Devin session or a judge.*
Last synced: Oct 2026, branch `devin/novel-extensions` @ `80df021`
(suite: 219 tests pass / 8 skipped).

---

## 1. What this project is

A **mechanistic virtual cell model of BK polyomavirus (BKPyV) infection in
kidney transplant recipients**, built for ISEF 2027 by Abhinav Mishra
(Health Careers High School, San Antonio). It is a Python platform
(`src/vcm`) whose flagship is a **23-dimensional ODE** that simulates the
virus inside the graft, its shedding into urine, the host immune response,
immunosuppressive drugs in real ng/mL units, viral quasi-species
(archetype → rearranged NCCR), a latent reservoir, and cell-to-cell spread.

**The pitch in one sentence:** prior BKPyV models reproduce population
kinetics; this one emits *falsifiable, individual-level predictions* —
when viremia starts (as a distribution, not a date), when a specific
patient clears, and which immunosuppression change is optimal — and it has
already been tested against real patient trajectories.

**Why it's novel:** no published BKPyV kinetic model produces an onset-time
distribution, an individual clearance forecast with credible intervals, or
a *derived* optimal policy. All three exist here, plus eight emergent
predictions that were not fitted.

**Working relationship:** Sunnie Thompson's lab (UAB Microbiology — the
group behind the S-phase single-cell finding, Needham 2024) reviewed the
model conventions, corrected the drug wording (production, not
replication), suggested the NCCR axis (now a headline result), and offered
longitudinal patient series (urine + plasma) — the next dataset.

---

## 2. Repository layout

```
src/vcm/
  core/models.py            CellState, ExperimentConfig, domain objects
  simulators/
    ode_system.py           THE ENGINE — 23-dim BKPyVODESystem (LSODA)
    bkpyv_ode_simulator.py  bridge ODE <-> CellState domain objects
    bkpyv_specific.py       legacy hourly discrete sim (kept, not canonical)
  clinical/
    viral_load_mapper.py    V (model units) -> copies/mL bridge, ASSUMPTION-LABELLED
    risk_prediction.py      logistic regression on synthetic cohort (honest null)
  ui/streamlit_app.py       12-page Streamlit app
  plugins/                  plugin registry (5 modules; bkpyv flagship)
  cli/main.py               `vcm` CLI (run configs, validate, etc.)
scripts/                  ~20 analysis/validation/generation scripts (see §6)
tests/                    pytest suite — 219 pass / 8 skipped
validation_data/funk2008/ 6 digitized patient CSVs + README (committed)
docs/                     ~25 docs incl. CLINICAL_PAPER.md (living paper),
                          HANDOFF.md (resume path), LIMITATIONS.md,
                          REFERENCES.md (full citation list), this file
outputs/                  gitignored artifacts (figures, run JSONs)
data/                     gitignored; GEO scRNA data restored on demand
```

Environments: `.venv/bin/python` for everything; streamlit 1.65; no CI
on the repo (verification is local pytest + validators).

---

## 3. The engine: 23-dim ODE, state by state

`solve_ivp(LSODA)`, all rates **per day**, states are fractions unless noted.
Index order is append-only (never renumber — tests pin indices).

| # | State | Meaning |
|---|-------|---------|
| 0 | V | Free virions, kidney pool (model units) |
| 1 | T | Intracellular T-antigen |
| 2 | G_v | Viral gene expression |
| 3 | C | Healthy tubular cells (fraction of carrying capacity) |
| 4 | I | Productively infected cells |
| 5 | D | Dead/damaged cells |
| 6 | CC | Cell-cycle permissiveness, 0–1, relaxes toward pathway drive (NOT an oscillating clock) |
| 7 | DNA | Host DNA-synthesis activity |
| 8 | E | Nonspecific effector immune proxy |
| 9 | IFN | Interferon |
| 10 | AK | Antiviral state |
| 11 | D_tac | Tacrolimus effect intensity (dimensionless; PK layer converts ng/mL → this) |
| 12 | D_sir | Sirolimus effect intensity (same) |
| 13 | P_rep | Replication-pathway activity |
| 14 | P_immune | Innate-immune pathway activity |
| 15 | T_naive | BKPyV-specific naïve T cells (repertoire fraction) |
| 16 | T_eff | BKPyV-specific effector T cells — the antigen-driven CTL arm |
| 17 | C_u | Healthy urothelial cells (urinary compartment) |
| 18 | I_u | Infected urothelial cells |
| 19 | V_u | Urinary viral load (model units) |
| 20 | F_rr | Rearranged-NCCR fraction of the **kidney** pool (0–1) |
| 21 | F_rr_u | Rearranged-NCCR fraction of the **urinary** pool (0–1) |
| 22 | L | Latently infected reservoir |

Infection is applied as discrete inoculum events at solver segment
boundaries (`get_infection_conditions(inoc)`). Dosing contexts:
`{"drug": [{"start","stop","target"}]}` dimensionless, or
`{"start","stop","trough_ng_ml"}` for real PK dosing (`trough_ng_ml` wins).

## 4. Mechanisms, one by one

**S-phase gate (invariant #1).** T-antigen production is Hill-gated by
host CC × DNA-synthesis permissiveness (Needham 2024: in primary RPTE
cells, no host S phase → no TAg → no replication). **T-antigen never
appears in dCC/dt** — it does NOT drive cell-cycle entry. Thompson
confirmed this reading and added that something else (not TAg) drives
re-entry; her lab has a manuscript under review on what.

**Production vs replication wording (invariant #2, Thompson's fix).**
Tacrolimus permits viral *production* by blunting immune control of
infected cells; sirolimus can suppress production by removing S-phase
permissiveness. Neither drug touches intracellular genome copying, so
docs say "production" — "replication" is reserved for genome copying.

**Pharmacokinetics.** `trough_ng_ml` schedules relax the drug-effect
state toward `trough/ref_trough` with the drug's elimination half-life
(`tac_pk_half_life_days`, `sir_pk_half_life_days`; reference troughs
`tac_ref_trough_ngml` = 8, `sir_ref_trough_ngml` = 4). Legacy
dimensionless `target` contexts still work.

**Immune arms.** Two layers: nonspecific E/IFN/AK (innate proxy) and the
BKPyV-specific T-cell arm: T_naive → T_eff via priming + antigen-driven
expansion capped by `tcell_carry` = 1.0 homeostatic ceiling (expansion
goes negative above carry → an infused bolus contracts within days);
`d_teff` = 0.05 decay; `tcell_kill` = 0.25 kill rate NOT
tacrolimus-attenuated; `tac_tcell_suppression` = 12.0 gates priming and
expansion. This ceiling is exactly why the VST result is transient-only
(see §7, #8).

**Urinary compartment.** Urothelial C_u/I_u/V_u with kidney→urine drainage
(`drain_kidney`) and bidirectional cross-feeding; urothelial production is
shedding-driven rather than S-phase-coupled (so the rr early-gene
advantage is weakened there — the basis for the plasma>urine split).

**NCCR quasi-species (two pools).** F_rr (kidney) and F_rr_u (urinary)
track the rearranged fraction: spontaneous rearrangement during
replication (`nccr_emergence_rate`), early-gene overexpression gain
(`rr_early_gain` ≈2.0 on production), capsid-suppression cost
(`nccr_capsid_fraction` ≈0.5 on expression, applied to virion yield via
partial coupling `rr_capsid_virion_cost` = 0.4 — full 1.0 coupling breaks
the Funk curtailment signature), reversion (`nccr_reversion_rate`), and a
blunted advantage in the urinary pool (`uro_rr_advantage`). Result:
rr enriches in plasma ~2.7× over urine — Gosert's clinical direction,
emergent rather than coded.

**Latent reservoir L.** `latent_fraction` = 0.02 of new infections abort
into latency (nonproductive, non-cytopathic); `lat_decay` drains slowly;
reactivation flux `lat_reactivation`·L·(1 − tac_immune_effect) feeds I —
the mechanistic reservoir behind the reactivation-onset model and the
reason viremia resurges when IS is interrupted mid-clearance.

**Cell-to-cell spread.** `c2c_rate` = 0.03 (≈10% of β by design):
direct contact transmission independent of extracellular V — structurally
shielded from free-virus clearance; the persistence channel.

**Copies/mL bridge (invariant #3).** `ViralLoadMapper.normalized_to_copies`
is **scalar-only** (numpy arrays crash it — the Regimen Design chart bug
from the first audit). Piecewise log-linear anchor bridge, explicitly
labelled assumption; anchors: 10³ cp/mL screening, 10⁴ presumptive-PVAN.
The bridge is assumption-labelled everywhere — never present readouts as
patient-calibrated truth.

## 5. The scripts, what each does, and the verified numbers

| Script | Does | Verified output |
|--------|------|-----------------|
| `benchmark_viral_load.py` | 6 published signatures | OVERALL PASS: t½ 6.6d (Funk 0.25–17); curtail 50% ineff / ≥80% ≤7wk / ≥90% ≤3wk; u:p 445× (Funk ~3000×, gap documented); rr plasma:urine 2.7× (≥2×; Gosert ~5×); F_rr 0.759@120d; T_eff 0.29× |
| `reactivation_onset.py` | Poisson hazard h = 0.005 + 0.0025·(tac−3), 150 MC patients/trough | Onset DISTRIBUTION: median ≈6.3–6.4 wk @tac 8 (p10–p90 ≈1–17 wk); P(onset) 57%→96% across troughs 3→8; higher trough → earlier AND tighter |
| `early_forecast.py` | jittered virtual cohort + 0.25-log₁₀ noise; (β,δ) from first 3 weekly points; Metropolis posterior → conversion | MAE 0.17 wk (~1.2 d); 0/24 clear/not-clear discordance; posterior median 7.0–7.2 wk vs truth 6.7–7.1; honest 90% CIs + P(clear) |
| `optimal_control_analysis.py` | (conversion day, sir trough) continuous controls, Nelder-Mead | Optimum hits feasibility boundary — convert IMMEDIATELY at max mTOR signal; delay costs ~2.4 objective-weeks/week |
| `optimize_reduction_schedule.py` | trigger×regimen grid with rejection counterweight | Conversion Pareto-dominant (clears ~wk 7.8–8.3, rebound 0.103); TAPER PARADOX: taper at any trigger → MORE rr emergence than no action (0.99 vs 0.93) |
| `screening_policy_analysis.py` | act at 10³ vs 10⁴ vs PVAN | Cadence itself changes clearance timing + emergence — monitoring is an intervention |
| `pharmacogenomic_stratification.py` | CYP3A5 → tac clearance (CPIC) | Expressor trough ~4.8 vs 8.0 → onset ~10 vs ~6.3 wk; genotype-blind dosing assigns non-expressors viral-onset risk, expressors rejection risk |
| `tcell_therapy_simulation.py` | VST CTL bolus under maintained tac | Transient dips only at ~4× homeostatic ceiling; NEVER durable clearance — ceiling + calcineurin-blocked expansion + L/c2c reseeding (matches clinical VST transience) |
| `identifiability_analysis.py` | which params patient data could constrain | Identifiability report — bounds what fits can mean |
| `fit_patient_series.py` | 5-param fit (β, δ, pre_roll, act_week, tac_lo) to real CSVs; 60/40 honest holdout | Cohort holdout MAE **0.82 vs 1.36** log₁₀; wins 5/6; per-patient table in §8; figure → `outputs/figures/patient_fits.png` |
| `generate_figures.py` | 4-panel composite | `outputs/figures/composite_figure.png` |
| `validate_bkpyv.py`, `validate_isef.py` | pipeline gates | exit 0 |
| `download_gse317012.py`, `cluster_gse317012.py`, `analyze_gse317012.py` | scRNA pipeline | 26 biopsies; 34,987 QC cells; MWU p=0.16 (descriptive) |
| `compare_patient_series.py` | Thompson-series comparer | data-gated, waiting for series |

## 6. Streamlit app — 12 sidebar pages

Home · Simulation · Visualization · Single-Cell · Viral-Load Validation ·
Risk Prediction · Comparison · Parameters & Assumptions · Review Bundle ·
Regimen Design · Documentation · **🏆 Findings** (headline signatures,
live posterior-forecast demo, onset explorer). AppTest-verified, 0
exceptions; `test_streamlit_app.py` asserts all 12 + per-page content.

## 7. The eight emergent predictions (the paper's §3.2)

1. Onset is a distribution (median ≈6 wk @tac8, in the 4–16 wk cluster).
2. Clearance forecastable from 3 early points (MAE 0.17 wk, 0 discordance).
3. Optimal policy derived: convert immediately, max tolerated mTOR signal.
4. CYP3A5 genotype orthogonalizes onset-vs-rejection risk.
5. rr-NCCR segregates plasma>urine (~2.7×, Gosert direction).
6. Monitoring cadence is itself an intervention.
7. Persistence has two channels (c2c + L resurge on interrupted IS).
8. VST durability gated by residual IS (honest negative result).

## 8. Real-data validation — the Funk-2008 digitization

Source: Funk et al. 2008, *Am J Transplant* 8(11):2368–2377,
doi:10.1111/j.1600-6143.2008.02402.x, CC BY-NC-ND. No open qPCR dataset
exists anywhere (PMC supplements, Zenodo, Dryad, GEO, edoc all checked);
Fig 3's six patient panels were digitized programmatically —
pixel-calibrated against the figure's own log₁₀ 4/7 dashed reference
lines; plasma vs urine markers separated by blob core intensity (dark
fill = plasma). Precision ≈ ±0.15 log₁₀, ±0.5 wk. CSVs +
`validation_data/funk2008/README.md` are committed (`.gitignore` carries
`!` negations because `/data/` and `*.csv` are ignored).

**Fit design:** patient week-0 is mid-course, not infection onset → the
model pre-rolls under high-IS (`pre_days`) so the state sits on its own
quasi-steady manifold; the intervention is IS reduction (T_HI 9 ng/mL →
fitted `tac_lo`) at fitted `act_week`.

**Results (post-audit, true timescale):**

| Patient | RMSE | R² | Holdout MAE | Baseline | Verdict |
|---------|------|-----|-------------|----------|---------|
| 00107 | 0.47 | 0.53 | 0.29 | 0.62 | wins |
| 00885 | 1.07 | 0.17 | 1.49 | 1.48 | tie (refractory) |
| 06119 | 1.68 | −0.40 | 0.43 | 2.03 | wins big |
| 15207 | 0.06 | **0.98** | 0.42 | 0.68 | wins |
| 19771 | 0.92 | 0.53 | 1.69 | 2.74 | wins |
| 27447 | 0.59 | 0.52 | 0.61 | 0.61 | wins (marginal) |
| **Cohort** | | | **0.82** | **1.36** | **5/6** |

Failure modes, named not hidden: 06119's ~10⁷·⁷ spike exceeds the ~10⁶
plasma ceiling (structural limit, full-R² negative despite a 4.7× holdout
win); 00885 is refractory fluctuation the model declines to chase; the
surges fits can't reach are exactly the L/stochastic-release events —
the mechanism-pointer argument in the paper.

**Audit catch worth remembering:** `load_series` originally returned
weeks×7 then callers multiplied by 7 again → 49×-stretched timeline —
self-consistent (MAE/R² were valid) but the figure axis showed
days-as-weeks and fitted params were distorted. Fixed; the true-timescale
refit *improved* every headline number. Lesson in `~/memory`.

## 9. Invariants — things that must never break

1. T-antigen never drives dCC/dt (S-phase gate is host-driven only).
2. "Production" for drug effects; "replication" only for genome copying.
3. Copies/mL bridge always assumption-labelled; mapper is scalar-only.
4. State indices append-only (23 = L; never renumber).
5. `validate_bkpyv.py`/`validate_isef.py` must exit nonzero on failure.
6. YAML floats written `1.0e-6`.
7. ISEF **2027** everywhere.
8. ~1.6k ruff findings are grandfathered debt — NEVER `ruff --fix`
   whole files; touched files stay clean.
9. `data/` and `outputs/` gitignored; `validation_data/` committed via
   `!` negations. No fake/placeholder data anywhere.
10. PR #2 stays **unmerged** until Abhinav confirms the arc.
11. `rr_capsid_virion_cost` stays 0.4 — full coupling breaks curtailment.

## 10. Tests

219 pass / 8 skipped (32.9s). Structure: `test_bkpyv_ode.py` (engine +
invariants + VST-transience + patient-harness), `test_streamlit_app.py`
(12-page AppTest coverage), `test_gse317012_data.py`,
`test_viral_load_mapper.py`, `test_risk_prediction.py`, plugin/CLI tests.
Run: `.venv/bin/python -m pytest tests/ -x -q`.

## 11. Docs map

- `CLINICAL_PAPER.md` — the living competition paper (v2.0)
- `PROJECT_EVERYTHING.md` — this file
- `ISEF_ABSTRACT.md` — 218-word abstract
- `PRESENTATION_POINTS.md` — judge narrative + Q&A
- `HANDOFF.md` — resume path + research queue for next session
- `METHODS.md`, `ODE_IMPLEMENTATION.md`, `LIMITATIONS.md`, `STATE.md`,
  `DECISIONS.md`, `ROADMAP.md`, `PROGRESS.md`, `BKPYV_MODEL_CARD.md`,
  `GSE317012_ANALYSIS.md`, `bkpyv_research_to_model_map.md`,
  `ISEF_PROJECT_OVERVIEW.md`, `ISEF_RESEARCH_PLAN.md`,
  `CLINICAL_PAPER.md` (canonical living paper + PDF renderer),
  `REFERENCES.md` (full citations)

## 12. Git / PR state

- `main`: PR #1 merged (squash `459c970`) — the 22-dim + PK + NCCR arc.
- `devin/novel-extensions`: **PR #2 OPEN, UNMERGED by user constraint** —
  everything since: c2c, L, VST, posterior UQ, PGx, optimal control,
  screening, reactivation, real-data validation, packaging, Findings
  page, clinical paper, figure polish. Latest `80df021`.
- Merge style when confirmed: `gh pr merge --squash` (matches repo's
  linear history).

## 13. How to run

```bash
cd ~/repos/BKPyV-VCM
.venv/bin/python -m pytest tests/ -x -q          # suite
.venv/bin/python scripts/validate_bkpyv.py        # gate 1
.venv/bin/python scripts/validate_isef.py         # gate 2
.venv/bin/python scripts/benchmark_viral_load.py  # 6 signatures
.venv/bin/python scripts/fit_patient_series.py    # ~13 min: real fits
.venv/bin/python scripts/generate_figures.py      # composite figure
.venv/bin/python -m streamlit run src/vcm/ui/streamlit_app.py
```

## 14. What's next (queue, from HANDOFF.md)

- Thompson longitudinal urine+plasma series — THE next dataset; the fit
  pipeline is built and waiting (`compare_patient_series.py`,
  `fit_patient_series.py`).
- Fit with reservoir-release events enabled (mechanism-pointer experiment).
- Agnoprotein latency refinement; Bayesian fit of the virtual cohort.
- Merge PR #2 after Abhinav's confirm + a final audit pass.

## 15. Key references

Funk 2008 (model + Fig 3 data) · Funk 2006 (t½) · Gosert 2008 (rr u/p) ·
Hirsch 2016 (sirolimus mechanism) · AST IDCOP 2019 + Kotton 2024
(thresholds) · Needham 2024 (S-phase gate; Thompson lab) · CPIC 2015
(CYP3A5) · Myers 2025 (prior patient-scale ODE — positions novelty) ·
GSE317012 (scRNA cohort). Full list: `docs/REFERENCES.md`.
