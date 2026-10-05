# Paper Outline — BKPyV Mechanistic Cell Model

1. **Introduction**
   - BKPyV nephropathy burden in kidney transplantation; the
     immunosuppression trade-off.
   - Gap: drug→intracellular-mechanism link is not directly observable;
     motivates a mechanistic model.
2. **Methods**
   - ODE architecture (per-day rates; TAg → replication → capsid/production;
     immune-control + drug terms) — `src/vcm/`, `configs/bkpyv_ode_*.yaml`.
   - Assumption-labelled copies/mL bridge (anchors: AST IDCOP 2019,
     Kotton 2024) — `clinical_validation_report.json`.
   - Harness + internal-consistency checks (`harness_check_report.json`,
     `isef_validation_report.json`) — 4 ODE + 4 legacy configs pass.
   - GSE317012 reanalysis methods (sparse per-sample QC; 34,987/294,286
     cells; marker-proxy fractions; five pathway scores; sample-level
     contrasts) — `docs/GSE317012_ANALYSIS.md`.
3. **Results**
   - Directional ordering reproduced: tacrolimus > untreated > sirolimus
     (peaks 6.77M / 2.43M / 31k copies/mL; `clinical_validation_report.json`).
   - NCCR archetype vs rearranged trade-off (`outputs/bkpyv_review/`).
   - Drug-switch scenario: early tac→siro switch collapses week-52 viremia
     (`outputs/clinical/drug_switch_results.json`).
   - GSE317012 descriptive trends: immune-proxy fraction +0.187 at peaking
     (n=5 vs 12; MWU p=0.16 — not significant; descriptive only).
4. **Discussion**
   - Hypothesis-generator value; what a wet-lab test could check.
   - Limitations: phenomenological parameters, no patient cohort, small
     peaking n, marker-proxy fractions are not deconvolution.
5. **Conclusion** — mechanistic framing generates testable questions;
   clinical use explicitly out of scope.
