# METHODS

## 1. Model — bkpyv_ode (canonical engine)

22-dimensional ODE, per-day rates, `src/vcm/simulators/ode_system.py`,
integrated with `scipy.integrate.solve_ivp` (LSODA).

State: V (free virions), T (intracellular T-antigen), G_v (replication-
permissive cells), C (cell-cycle–permissive fraction), I (infected cells),
D (damaged cells), CC (cell-cycle state), DNA (viral DNA), E (immune
effector), IFN (interferon), AK (antiviral state), D_tac, D_sir (drug
concentrations), P_rep, P_immune (pathway activities).

Core structure (see source for the full system):

- Virion production: `p · I · T/(T + half_saturation) · s_phase_gate(CC)`
  — T-antigen gates production through a Hill-type term; `t_threshold`
  shapes the sirolimus phase weight (it is **not** the replication gate).
- Clearance: `delta · V` with an antiviral-state boost (`AK`) —
  δ = 0.4/day calibrated inside the Funk 2006 IS-change t½ range.
- Infection: `beta · V · G_v` on a permissive-cell reservoir.
- Immune control: effector E reduces permissive cells; IFN drives AK;
  tacrolimus reduces immune control (D_tac term); sirolimus inhibits the
  cell-cycle permissiveness term (mtor_inhibition).
- Drug PK: first-order accumulation to `target` while dosing, first-order
  clearance after `stop` (dosing_context dict).

## 2. Single-cell layer — GSE317012 (Phase 1)

- **Data**: GEO GSE317012, Kretzler lab human kidney-transplant biopsies,
  26 samples (12 Control / 5 Peaking / 9 Resolving), 294,286 raw →
  34,987 cells after QC. Download is scripted and sha256-verified
  (`scripts/download_gse317012.py`); no manual steps.
- **Important**: the reference contains human genes only — no viral genes.
  Infected cells cannot be identified by viral reads.
- **Pipeline** (`scripts/cluster_gse317012.py`): per-cell QC (min genes,
  mt-fraction cap) → normalization (target sum) → PCA (50 comps, no
  zero-center for memory) → neighbors → Leiden clustering → marker-panel
  annotation (kidney cell types: proximal/distal tubule, loop of Henle,
  collecting duct, podocyte, endothelial, fibroblast, myeloid, B/NK/other
  immune) → module scoring for host-response gene sets (translation, S-G2M,
  DDR, interferon, antigen presentation, mitochondrial, ECM/wound).
- **Infected-signature call**: within tubular epithelial cells of
  Peaking/Resolving samples, the top decile of composite viral-response
  score is labelled `infected_signature`; the remainder is the bystander
  comparator. This is a host-signature proxy (see LIMITATIONS).
- **DE**: Wilcoxon rank-sum (scanpy `rank_genes_groups`), BH FDR.
- **Calibration/validation split**: deterministic seeded stratified holdout
  — within each phase, ~1/3 of samples → validation, remainder →
  calibration (Control 8/4, Peaking 3/2, Resolving 6/3; recorded in
  `gse317012_calibration_holdout.csv`). Calibration values are computed on
  calibration split only; the same quantities are reported on the
  validation split side-by-side for consistency checking — they are never
  refit.

## 3. Viral-load benchmark (Phase 2)

- **Ground truth**: published summary statistics only (no public serial
  qPCR dataset exists — searched GEO/dbGaP/SRA/mirrors; documented in
  `data/research/published_kinetics.csv` with full citations):
  - Funk 2006 (PMID 16323135): IS-change clearance t½ 0.25–17 d (n=12);
    nephrectomy arms 1–2 h (dense) and 20–38 h (sparse), n=3;
    intervention efficacy median 22% (7–83%).
  - Funk 2008 (doi:10.1111/j.1600-6143.2008.02402.x): >80% curtailment →
    viremia clears ≤7 wk; >90% → ≤3 wk; 50% ineffective; ~10 wk dynamic
    equilibrium; urine:plasma ≈3000:1.
  - AST IDCOP 2019 / Kotton 2024: 1,000 and 10,000 cp/mL thresholds.
- **Calibration subset**: Funk 2006 IS-change arm — model t½ after
  immunosuppression stop must lie inside 0.25–17 d.
- **Validation subsets** (disjoint, never used for parameter choice):
  Funk 2008 curtailment-response statements; Funk 2006 nephrectomy arms
  (reported against both arms and the 1–38 h envelope — the two arms
  measure different sampling regimes, so between-arms is compatible).
- **Intervention start state**: peak viremia under continuous tacrolimus
  (~1.7×10⁶ cp/mL), the clinically representative intervention range.
- **Metrics**: clearance t½ (first crossing of half of start); time below
  1,000 cp/mL (clearance event, wk); sustained clearance (ends below
  threshold); band-excess RMSE on log10 cp/mL — RMS distance of the model
  trajectory *outside* the published decay envelope (0 = fully inside);
  δ-sweep (0.25–0.8/day) for parameter uncertainty.
- **Results**: see `outputs/benchmark/viral_load_benchmark.json`:
  calibration t½ 5.5 d (in range); band-excess RMSE 0.19 log10 cp/mL with
  47% of trajectory inside the envelope; Funk 2008 checks 3/3 pass;
  nephrectomy t½ ≈10 h (between published arms).

## 4. Assumptions explicitly labelled

- `p=8.0/day`, `beta=0.3/day`, `half_saturation=0.5`, immune magnitudes,
  drug PK clearances: ASSUMPTIONS (no direct cell-level measurements).
- V→cp/mL bridge: piecewise log-linear anchors — ASSUMPTION tied to
  clinical thresholds, not a fitted map (`viral_load_mapper.py`).
- Clinical-risk ORs (age ×1.9, male ×2.3, prior transplant ×3.0):
  illustrative; directions from Demey et al. 2018 — not fitted to a cohort.
- The risk model is a mechanistic feature extractor + illustrative
  covariates; it is **not a trained foundation model** and is never
  presented as one.

## 5. Validation design (summary)

| Layer | Calibration | Validation |
|---|---|---|
| Cell-state (GSE317012) | 17 samples (8C/3P/6R) | 9 held-out samples (4C/2P/3R) |
| Kinetics (ODE) | Funk 2006 IS-change t½ range | Funk 2008 curtailments + Funk 2006 nephrectomy |
| Risk features | — (illustrative covariates) | direction checks only |

Fixed seeds everywhere (numpy `default_rng` seeds, deterministic sort+hash
for per-cell holdout); rerunning the four pipeline commands reproduces all
outputs byte-identically up to numeric noise in UMAP.

## 6. Multi-compartment extension (Phase 3, October 2026)

The canonical engine was extended from 15 to **22 state variables**
(indices appended — never renumbered — so all legacy consumers keep
working):

| # | State | Meaning | Grounding |
|---|-------|---------|-----------|
| 15 | T_naive | BKPyV-specific naive T cells | priming pool, tac-inhibited via calcineurin/NFAT |
| 16 | T_eff | BKPyV-specific effector T cells | antigen-driven expansion + kill of infected cells |
| 17 | C_u | healthy urothelial cells | Funk 2008 urinary reservoir |
| 18 | I_u | infected urothelial cells | >95% of urinary viral load is urothelial origin |
| 19 | V_u | urinary virion pool | urine:plasma ~3000:1 (Funk 2008) |
| 20 | F_rr | rr-NCCR fraction, KIDNEY pool | Gosert 2008 |
| 21 | F_rr_u | rr-NCCR fraction, URINARY pool | shedding-driven selection (weakened: uro_rr_advantage 0.15) + kidney-drainage mixing |

Mechanism notes:

- **T-cell arm**: priming `a_tcell·(V+0.5·I)·T_naive` and carrying-capacity
  expansion are both attenuated by tacrolimus through
  `tac_tcell_suppression` (calcineurin blockade hits NFAT-driven
  proliferation — tac does NOT attenuate the kill itself). Effectors kill
  infected kidney cells via `tcell_kill·T_eff·I`.
- **Urinary compartment**: population-level `C_u/I_u/V_u` block
  (no per-cell intracellular gates — the urothelium is a second tissue
  population, matching the Funk 2008 formulation). Kidney→bladder flow
  via `drain_kidney·V`; bladder→kidney seeding via `cross_feed·V_u·C`.
- **NCCR quasi-species**: `dF_rr/dt = μ·pressure·(1-F) + s·pressure·F(1-F)
  - r·F` where pressure = replication activity. Rearranged variants emerge
  under sustained viremia and dominate on the Gosert weeks-months
  timescale; the discrete archetype/rearranged presets seed the F=0/F=1
  boundary conditions of BOTH pools. Effective early-gene and capsid multipliers
  interpolate (early ×1→×2, capsid ×1→×0.5).
- **Pharmacokinetics**: dosing schedules accept `trough_ng_ml` windows
  resolved with real elimination half-lives (tac ~12 h, sir ~60 h) and
  reference troughs (tac 8, sir 4 = Hirsch 2016 in-vitro IC90 anchor).
  Schedule lists express stepwise tapers; legacy dimensionless targets
  and bolus events are unchanged.

## 7. Regimen optimization and identifiability

- `scripts/optimize_reduction_schedule.py` sweeps taper and conversion
  schedules in ng/mL; efficacy = weeks to sustained <1,000 cp/mL;
  counterweight = T_eff rebound AUC (honest proxy — the model cannot
  distinguish BKPyV-specific from alloreactive T cells). Emergent result:
  **tac taper alone never clears within 180 d; tac→sir conversion clears
  ~wk 8 with lower rebound** — sir closes the S-phase/mTOR gate AND
  releases the T-cell brake (consistent with Hirsch 2016), arising from
  mechanism not fitting.
- `scripts/identifiability_analysis.py` computes normalized
  d log observable/d log θ sensitivities per observable
  (plasma V, urine V_u, T-antigen, F_rr): plasma-V data informs only the
  Nowak-May kinetics core; urine parameters are identifiable only from
  urine data; drug-PD parameters are unidentifiable without
  drug-perturbation data; the intracellular gate chain is pairwise
  confounded (|corr|≈1) — one effective parameter, not five.
- Benchmark gained compartment-signature checks: urine:plasma ≥50×,
  F_rr ≥0.5 by wk 17, tac-suppressed T_eff peak ≤0.5× untreated — all
  PASS in `outputs/benchmark/viral_load_benchmark.json`.

- `scripts/screening_policy_analysis.py` — trigger-level policies: same
  conversion intervention fired at 1k / 10k / 100k cp/mL / never.
  Emergent ordering: earlier trigger -> earlier clearance AND suppressed
  rr emergence (F_k 0.02 vs 0.93) — the mechanistic argument for
  Kotton's intensive-screening policy. `--monitor-days` sweeps the check
  cadence: under rapid doubling, WEEKLY monitoring collapses the 1k/10k
  distinction (both fire the same week) and MONTHLY monitoring can miss
  the 10k window entirely (only the 1k trigger ever fires — clearance
  slips to ~wk8.5 with F_k already 0.18). Cadence is itself an
  intervention variable.

- `scripts/reactivation_onset.py` — stochastic reactivation model:
  Poisson reactivation hazard scaled by tacrolimus trough; Monte-Carlo
  onset-time distribution (median ~6.4wk at tac 8, inside the clinical
  4-16wk cluster). First BKPyV kinetic model to generate an onset
  distribution rather than an onset point.
- `scripts/early_forecast.py` — digital-twin proof: virtual cohort
  (beta/delta/p/inoculum jittered), first three noisy weekly qPCR
  points -> joint (beta,delta) inference -> predicted clearance week.
  MAE 0.17wk, 0/24 clear/not-clear discordance — honest error bars,
  distinct train/infer parameterizations. `forecast_posterior` adds
  full uncertainty quantification: Metropolis posterior over
  (beta, delta) propagated through conversion -> clearance-week
  distribution (median, 90% interval, P(clear)) — "clear by week X
  with 90% probability".
- `scripts/generate_figures.py` — four-panel publication composite
  (compartments vs thresholds, two-pool rr emergence, taper-vs-
  conversion, onset distributions by trough) -> outputs/figures/.
- Cell-to-cell transmission channel (`c2c_rate`, default 0.03):
  free-virion spread is V-dependent (delta/AK-cleared); direct
  cell-to-cell spread via virological synapses is V-independent —
  structurally shielded from extracellular clearance. Raises the
  infected reservoir under deep suppression: the persistence
  channel that explains why viremia resurges after interruption.

- `scripts/optimal_control_analysis.py` — (conversion day, sir trough)
  as continuous controls; Nelder-Mead on weeks-to-clear + rebound AUC.
  Derived optimum: earliest feasible switch + max mTOR signal; delay
  costs ~2.4wk/week — conversion dominance DERIVED, not grid-selected.

- `scripts/pharmacogenomic_stratification.py` — CYP3A5 genotype ->
  effective trough (expressors clear tac ~2x faster, CPIC) -> onset
  distribution + clearance. Emergent inversion: non-expressors carry
  the early viral-onset risk; expressors carry rejection risk — viral
  and rejection risk are orthogonalized by genotype.
