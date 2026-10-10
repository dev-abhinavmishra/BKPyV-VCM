# STATE — repository audit (Phase 0)

Audit date: 2026-10-05. Branch: `devin/finish-vcm`.

## What exists

| Area | Path | Status |
|---|---|---|
| Canonical engine | `src/vcm/simulators/ode_system.py` | Working — **22-dim ODE**, per-day rates (+T-cell arm, urinary compartment, NCCR quasi-species) |
| ODE wrapper | `src/vcm/simulators/bkpyv_ode_simulator.py` | Working |
| Legacy simulator | `src/vcm/simulators/bkpyv_simulator.py` | Working, **superseded for analysis** (see below) |
| Plugin/parameters | `src/vcm/plugins/transplant/bk_polyomavirus/` | Working |
| Clinical bridge | `src/vcm/clinical/viral_load_mapper.py` | Working — assumption-labelled V→cp/mL anchors |
| Risk prediction | `src/vcm/clinical/risk_prediction.py` | Working — illustrative ORs, not fitted |
| Dashboard | `src/vcm/ui/streamlit_app.py` | Working — 11 pages (+💉 Regimen Design) |
| Review bundle | `src/vcm/viz/bkpyv_review.py` | Working |
| GSE317012 biopsy pipeline | `scripts/analyze_gse317012.py` | Working — 26 samples, pathway scores |
| GSE317012 cell-level pipeline | `scripts/cluster_gse317012.py` | **Added this session** — QC/cluster/DE/holdout |
| GEO download | `scripts/download_gse317012.py` | **Added this session** — 1.9 GB, sha256-verified |
| Viral-load benchmark | `scripts/benchmark_viral_load.py` | **Added this session** — now incl. compartment-signature checks (urine:plasma, F_rr, T_eff) |
| Regimen optimizer | `scripts/optimize_reduction_schedule.py` | **Added Oct 2026** — ng/mL taper/conversion sweep vs T_eff rebound |
| Identifiability | `scripts/identifiability_analysis.py` | **Added Oct 2026** — per-observable sensitivity ranks |
| Validators | `validate_bkpyv.py`, `validate_isef.py` | Working — **fixed this session** |

## What runs

- `pytest tests/` — **208 tests pass** (was 184 passing before this session;
  +14 new: benchmark metrics, sha256 helper, holdout split, e2e smoke).
- `validate_bkpyv.py` — exits 0 (after fix, see below).
- `validate_isef.py` — 4/4 pattern checks PASS, "EXCELLENT" (after fix).
- `streamlit run src/vcm/ui/streamlit_app.py` — launches, all 10 pages render
  (verified in browser and via `AppTest`).
- `scripts/benchmark_viral_load.py` — overall PASS (see docs/METHODS.md).

## What was broken, and the fix

1. **`validate_bkpyv.py`** — passed only `perturbations[0]` (the drug) to the
   simulator, so the infection perturbation never ran; crashed on `None`
   results. Fix: pass the full perturbation list; guard failed lookups.
2. **`validate_isef.py`** — ran pattern checks through the legacy discrete
   `bkpyv_specific` simulator, which does not reproduce the mechanism
   directions (sirolimus suppressed V only ~30% vs observed; tacrolimus
   direction inverted). Repointed checks at the canonical `bkpyv_ode`
   engine: baseline peakV=1.924, tacrolimus=2.180, sirolimus=1.118 —
   directions match the published biology.
3. **T-antigen threshold check** — knob was wired to `t_threshold`, which
   only shapes the sirolimus `phase_weight` and does not gate replication.
   Remapped to `half_saturation`, the actual T→production gate
   (`T/(T+half_saturation)` in `ode_system.py`).
4. **GSE317012 download** — previously manual/partial (scripted but not
   committed, hash check case-sensitive). New one-command script with
   `.part` resume + sha256 verification; all 26 GSM samples restored.

## ODE state variables (bkpyv_ode, 22 dims, per-day)

| # | State | Meaning | Source/basis |
|---|-------|---------|--------------|
| 0 | V | free virions | standard viral-kinetics compartment |
| 1 | T | intracellular T-antigen | BKPyV early protein |
| 2 | G_v | replication-permissive cells | viral reservoir |
| 3 | C | cell-cycle–permissive fraction | S-phase coupling (Weissbach 2024) |
| 4 | I | infected cells | standard compartment |
| 5 | D | damaged/apoptotic cells | pytopathic conversion |
| 6 | CC | cell-cycle state | host S-G2M |
| 7 | DNA | viral DNA | genome replication |
| 8 | E | immune effector level | adaptive control |
| 9 | IFN | interferon | innate antiviral |
| 10 | AK | antiviral state | IFN-driven refractory fraction |
| 11 | D_tac | tacrolimus concentration | dosing context |
| 12 | D_sir | sirolimus concentration | dosing context |
| 13 | P_rep | replication pathway activity | module score |
| 14 | P_immune | immune pathway activity | module score |
| 15 | T_naive | BKPyV-specific naive T cells | tac-inhibited priming (NFAT) |
| 16 | T_eff | BKPyV-specific effector T cells | antigen-driven expansion/kill |
| 17 | C_u | healthy urothelial cells | Funk 2008 urinary reservoir |
| 18 | I_u | infected urothelial cells | >95% of urinary load |
| 19 | V_u | urinary virion pool | urine:plasma ~3000:1 |
| 20 | F_rr | rearranged-NCCR fraction | Gosert 2008 emergence dynamics |

## Parameters with sources

| Parameter | Value | Source/status |
|---|---|---|
| `beta` infection rate | 0.3/day | literature-informed ASSUMPTION |
| `delta` viral clearance | 0.4/day | **calibrated** to Funk 2006 PMID 16323135 (t½ 0.25–17 d) |
| `p` production | 8.0/day | ASSUMPTION |
| `t_threshold` | 0.5 | mechanistic gate shape |
| `half_saturation` | 0.5 | ASSUMPTION (T→production gate) |
| `s_phase_bonus` | 2.0 | Weissbach 2024 J Virol (S-phase coupling) |
| `innate_immune_suppression` | 0.5 | direction cited, magnitude ASSUMPTION |
| `mtor_inhibition` | 0.5 | Hirsch 2016 Am J Transplant |
| tacrolimus/sirolimus clearance | in code | ASSUMPTION within published PK bounds |
| clinical ORs (age/sex/prior tx) | 1.9/2.3/3.0 | illustrative; directions from Demey 2018 |
| V→cp/mL bridge anchors | 0/0.02/0.2/1/3/5 → 0/1e2/1e3/1e4/1e6/1e7 | ASSUMPTION anchored on AST IDCOP 2019 + Kotton 2024 thresholds |
