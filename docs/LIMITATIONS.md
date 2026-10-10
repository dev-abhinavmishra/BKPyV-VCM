# LIMITATIONS — what this model cannot claim

Ordered roughly by severity. Each item states the boundary plainly.

## Data limits

1. **No patient-level viral-load data.** No public, de-identified serial
   plasma BKPyV qPCR dataset exists (GEO, dbGaP, SRA, Dryad/figshare checked
   2026-10-05). The kinetic benchmark is against *published summary
   statistics* (ranges, medians, threshold-response statements), not fitted
   curves. The model can claim "consistent with published ranges" — never
   "fitted to patients."
2. **GSE317012 has no viral reads.** The reference contains human genes only;
   infected cells are identified by a host-response *signature proxy*
   (top-decile composite score). `infected_signature` ≠ measured infection.
   It cannot report viral load, T-antigen level, or true infected fraction —
   only cells whose host program looks infection-associated.
3. **Biopsies, not transplant time series.** GSE317012 samples are
   cross-sectional biopsies labeled by clinicopathologic phase
   (Control/Peaking/Resolving), not longitudinal patient courses.
4. **No urine data.** The model addresses plasma kinetics only; the
   ~3000× urine:plasma ratio (Funk 2008) is cited but not modeled.

## Model limits

5. **Dimensionless viral state.** `V` is a normalized virtual load bridged
   to copies/mL by an *assumption* anchored on clinical thresholds
   (1,000/10,000 cp/mL). Absolute cp/mL predictions are bridge-dependent;
   only relative dynamics and threshold-crossing statements are meaningful.
6. **Most parameters are unfitted.** `p`, `beta`, `half_saturation`, immune
   magnitudes, drug PK clearances are assumptions or direction-consistent
   choices. Only `delta` is calibrated (inside Funk 2006's t½ range).
   Reported "parameter uncertainty" is a δ-sweep, not a posterior.
7. **The model clears faster than real patients.** Predicted time-to-
   clearance after curtailment (~1–1.4 wk) is faster than the published
   bounds (≤3/≤7 wk) — the model is directionally right and inside the
   published ≤-bounds, but its absolute timescale runs fast. Do not quote
   the model's clearance times as patient timescales.
8. **No pharmacodynamic dose-response.** Tacrolimus/sirolimus are binary
   dosing contexts with first-order PK envelopes — the model says
   *direction* (tac ↑, sir ↓) and crude timing, not a ng/mL-response curve.
9. **No T-cell/BK-specific immunity.** Immune control is a lumped effector;
   there is no BKPyV-specific T-cell compartment, so the model cannot speak
   to T-cell–directed therapy or immune monitoring claims.
10. **Nephrectomy result sits between published arms** (≈10–19 h vs 1–2 h
    dense-sampling and 20–38 h sparse-sampling arms). This is compatible
    with, but not a confirmation of, either arm.

## Dashboard/risk limits

11. **Risk covariates are illustrative.** age/sex/prior-transplant ORs are
    direction-only literature values, not fitted; the risk page is a
    hypothesis generator, not a predictor.
12. **Not a clinical tool.** Every page labels outputs as research
    hypotheses; nothing here should be read as patient-specific or
    clinical-decision-grade.

## Reproducibility limits

13. **Pipeline requires a ~1.9 GB download** (GEO FTP). Clean-clone
    reproduction is scripted but bandwidth-dependent; the single-cell page
    degrades to instructions when artifacts are absent.
14. **UMAP geometry is stochastic at the ~1% level** even with fixed seeds
    (BLAS/threading); cluster identities and counts are deterministic but
    embedding coordinates may differ in the last decimal.

## Multi-compartment extension (22-dim) caveats

- **Urothelial block is population-level** — no per-cell TAg/S-phase gates
  (the tissue-level Funk 2008 formulation); intracellular mechanism lives
  only in the kidney block.
- **F_rr is resolved per pool** (kidney F_rr vs urinary F_rr_u) — the
  plasma>urine rr enrichment is emergent (~3x vs Gosert's ~5x ratio,
  magnitude reported not forced). Unverified prediction: under conversion
  the urinary pool becomes rr-dominated while the kidney pool stays
  archetype — clinically unmeasured, flagged as hypothesis.
- **PK layer is a maintenance-trough approximation** (first-order
  approach to target trough at elimination t½); peak-trough oscillation
  is averaged out — appropriate for week-scale regimen questions, not
  dose-timing questions.
- **Rejection proxy is T_eff rebound AUC** — the model cannot separate
  BKPyV-specific from alloreactive immunity, so rebound is an honest
  *index*, not a predicted biopsy-proven rejection rate.
- **Identifiability is limited by design**: plasma-V data alone cannot
  learn urine-compartment, T-cell-source, or drug-PD parameters
  (sensitivity ~0); patient fitting must be restricted to the
  plasma-visible parameter subset or use urine+plasma series together.
- Emergent multipliers interpolate between archetype/rearranged
  endpoints — explicit `nccr_*_expression_multiplier` config overrides
  pin the F_rr=0 endpoint, not a fixed genotype.
