# ISEF Project Overview: BKPyV Virtual Cell Model

## Project Description

A mechanistic virtual cell model of BK polyomavirus (BKPyV)-infected renal tubular epithelial cells that simulates differential drug effects on viral replication kinetics and predicts BKPyV-associated nephropathy risk.

## Research Question

Can a mechanistic virtual cell model of BKPyV-infected renal tubular epithelial cells simulate the differential effects of tacrolimus vs. sirolimus on viral replication kinetics and improve risk prediction of BKPyV-associated nephropathy?

## Hypothesis

I hypothesized that a mechanistic virtual cell model incorporating drug-specific mechanisms informed by Hirsch 2016 (implemented as immune-control weakening for tacrolimus and reduced mTOR/S-phase permissiveness for sirolimus) and host cell cycle dynamics would reproduce the differential directional effects of immunosuppressants on BKPyV replication kinetics, and that virtual cell-derived viral load features would improve clinical risk prediction beyond traditional covariates.

## Methods Summary

- **Mechanistic Model Development**: Built a BKPyV-specific simulator incorporating drug mechanisms from Hirsch et al. (AJT 2016), cell-cycle coupling from single-cell transcriptomics (Weissbach et al., JVI 2024; Needham et al., PLoS Pathogens 2024), and viral replication dynamics using a system of ordinary differential equations.

- **Clinical Translation**: Developed an assumption-labelled anchor bridge mapping the model's dimensionless viral-load variable to indicative plasma copies/mL, anchored at consensus thresholds (AST IDCOP 2019; Kotton et al. 2024, Second International Consensus Guidelines). This is a qualitative visualisation step, not a calibration to patient data (no patient cohort was available).

- **Risk Prediction Model**: Trained logistic regression models on 500 synthetic patient cohorts using clinical covariates (age, sex, prior transplant, diabetes, tacrolimus use, HLA mismatch, donor age) and VCM-derived features (peak viral load, weeks above thresholds, area under curve, time to peak).

- **Validation**: Performed internal consistency checks (drug-effect directions, S-phase gate, persistence) and 5-fold cross-validation on the synthetic cohort.

## Key Results (v2, regenerated after ODE/bridge rebuild; supersedes all earlier numbers)

**Directional drug effects (ODE model, 52-week horizon, infection at day 21)**:
- Infection (no drug): peak ≈ 2.4×10⁶ indicative copies/mL; model sustains non-zero viral load at week 52
- Tacrolimus: peak ≈ 6.8×10⁶ copies/mL (≈2.8× higher than no-drug), higher week-52 burden — matches the in-vitro direction of Hirsch 2016 and the clinical association of tacrolimus with higher BKPyV risk
- Sirolimus: peak ≈ 3.1×10⁴ copies/mL, ≈78× lower than no-drug — matches the in-vitro inhibition direction (Hirsch 2016)
- Internal consistency checks pass: S-phase gate (disabling host-DNA coupling collapses T antigen), persistence under baseline, correct ordering across drugs

**Clinical bridge timing (assumption-labelled)**:
- Crossing of the 1,000 copies/mL screening anchor depends on the anchor placement, which is an assumption, not a measured conversion; we therefore report direction and burden rather than validated "weeks post-transplant" timing.

**Risk Prediction Performance (synthetic cohort, honestly a null result)**:
- Baseline clinical model AUC: 0.927 ± 0.012
- VCM-enhanced model AUC: 0.925 ± 0.011
- VCM features did not significantly improve prediction in synthetic cohort (honestly documented); note the synthetic cohort outcome was generated with a dependence on VCM features, so even this null is optimistic — the meaningful claim is that adding model-derived features did not inflate performance

## Limitations

- **Synthetic Data**: Risk prediction trained on synthetic outcomes (500 simulated patients) rather than real clinical cohorts; VCM-outcome relationships may not reflect reality. Outcome generation in the synthetic cohort itself depends on VCM features, which limits what the AUC comparison can show.

- **Model Simplifications**: The clinical bridge is an assumption-labelled anchor interpolation (piecewise log-linear), not a fitted mapping to measured viremia. Immune response is a coarse three-variable (E/IFN/AK) module; drug pharmacokinetics are dimensionless dosing intensities rather than plasma ng/mL concentrations.

- **Cell-context caveat (from expert review)**: the drug-effect magnitudes used at the mTOR axis are in-vitro, cell-type-dependent values; immortalized cell lines (e.g., HEK293) have dysregulated cell cycles and produce at least ~1 log lower BKPyV titers than primary RPTECs. Single-cell findings (Weissbach 2024) come from primary RPTECs in culture and are not automatically kidney-level kinetics.

- **Validation Scope**: Consistency checks are internal (direction/magnitude sanity); no prospective or retrospective patient-cohort validation has been performed. The model should be read as a mechanistic hypothesis generator.

## Future Work

- **Clinical Validation**: Prospectively validate VCM predictions against real transplant cohorts with measured BKPyV viral loads and outcomes. A collaboration channel offering longitudinal first-year **urine and plasma BKPyV series** has been opened (domain group behind the S-phase work); `scripts/compare_patient_series.py` is the ready-made comparison path.

- **NCCR rearrangement as a model axis**: the rearranged-NCCR variant (TAg overexpression, reduced capsid expression) is the clinically important form; planned experiments compare archetype vs rearranged under matched host permissiveness and map which factor combinations predict reactivation vs persistence.

- **Model Expansion**: Incorporate explicit immune response dynamics (T-cell activation, cytokine signaling), host genetics (HLA typing), and drug pharmacokinetics.

- **Risk Model Enhancement**: Train on real clinical data, explore non-linear models (random forest, gradient boosting), and integrate time-dependent covariates for dynamic risk prediction.

## Data Sources and Citations

**Mechanistic Research**:
- Hirsch HH, Yakhontova K, Lu M, Manzetti J. "BK Polyomavirus Replication in Renal Tubular Epithelial Cells Is Inhibited by Sirolimus, but Activated by Tacrolimus Through a Pathway Involving FKBP-12." *American Journal of Transplantation* 2016;16(3):821-832. PMID 26639422. (In-vitro drug mechanisms in primary RPTECs: sirolimus IC90 ~4 ng/mL via mTOR; tacrolimus activates replication via FKBP-12)

- Weissbach FH, Follonier OM, Schmid S, Leuzinger K, Schmid M, Hirsch HH. "Single-cell RNA-sequencing of BK polyomavirus replication in primary human renal proximal tubular epithelial cells identifies specific transcriptome signatures and a novel mitochondrial stress pattern." *Journal of Virology* 2024;98(12):e01382-24. PMID 39513696. (Cell-cycle coupling, DDR, mitochondrial stress; earlier drafts of this project mis-attributed this paper to "Jang JY" and cited unrelated GEO accessions - corrected here)

- Needham JM, Perritt SE, Thompson SR. "Single-cell analysis reveals host S phase drives large T antigen expression during BK polyomavirus infection." *PLoS Pathogens* 2024;20(12):e1012663. (Basis for the S-phase gate: host S phase precedes T-antigen accumulation)

- Funk GA, Steiger J, Hirsch HH. "Rapid dynamics of polyomavirus type BK in renal transplant recipients." *Journal of Infectious Diseases* 2006;193(1):80-87. PMID 16323135. (Patient-scale BKPyV clearance kinetics: half-lives 6 h-17 d after immunosuppression change; faster 1-2 h / 20-38 h phases were measured after allograft nephrectomy, i.e., source removal)

- Myers N, et al. "Modeling BK Virus Infection in Renal Transplant Recipients." *Viruses* 2025;17(1):50. PMC11768487. (Patient-scale ODE reference model; this project's contribution is the cell-mechanism coupling layer)

**Clinical thresholds / guidelines**:
- Hirsch HH, Randhawa PS; AST Infectious Diseases Community of Practice. "BK polyomavirus in solid organ transplantation." *Clinical Transplantation* 2019;33(9):e13528. (Screening schedule; 1,000 copies/mL sustained = probable PyVAN; >10,000 = presumptive PyVAN)

- Kotton CN, et al. "The Second International Consensus Guidelines on BK polyomavirus in kidney transplantation." *Transplantation* 2024;108(9):1834-1866. (Current consensus management thresholds)

**Risk factors / risk models**:
- Demey B, et al. "Risk factors for BK polyomavirus infection and BKPyV-associated nephropathy..." *Journal of Clinical Virology* 2018;109:6-12. PMID 30343190. (Systematic review; tacrolimus, male sex, older age, prior rejection as associated factors). An earlier draft cited "Fang et al. 2022 meta-analysis, PMC9428263" with specific ORs; that PMC ID is a different paper (a single-center dynamic-prediction study), so those ORs have been removed.

- Fang Y, et al. "Dynamic risk prediction of BK polyomavirus reactivation after renal transplantation." *Frontiers in Immunology* 2022;13:971531. PMC9428263. (Dynamic Cox-model risk prediction; source of the mistakenly-referenced PMC ID)

- Yamauchi K, et al. "Development of a risk prediction model for BK polyomavirus-associated nephropathy after kidney transplantation." *Renal Failure* 2025. (Integer-based risk score with age, sex, prior transplant; reported AUC 0.65-0.68)

**Single-cell data**: GEO accession GSE317012 (BKPyV-infected kidney-biopsy scRNA-seq, JCI Insight 2026) plus microarray series GSE47199/GSE75693/GSE72925 are BKPyV-relevant and were reanalyzed through the project's data pipeline — per-sample QC, marker-proxy fractions, pathway scores, and condition contrasts at the sample level (26 biopsies; 12 Control / 5 Peaking / 9 Resolving; 34,987 of 294,286 cells retained after stringent QC). See [GSE317012_ANALYSIS.md](GSE317012_ANALYSIS.md). **Model parameters remain phenomenological** — the reanalysis provides descriptive trend support only; no coefficient was updated from these data.

**Software and Algorithms**:
- Python 3.10+, NumPy, SciPy, scikit-learn, matplotlib
- ODE integration via scipy.solve_ivp (LSODA, rtol=1e-6)
- Piecewise log-linear anchor bridge for copies/mL (exact by construction)
- Logistic regression with L2 regularization (sklearn)
- Stratified 5-fold cross-validation