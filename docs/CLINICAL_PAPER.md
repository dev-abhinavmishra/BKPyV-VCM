# A Mechanistic Virtual Cell Predicts BK Polyomavirus Onset, Clearance, and Optimal Immunosuppression Strategy in Kidney Transplant Recipients

**Abhinav Mishra** — Health Careers High School, San Antonio, Texas
**Status:** Living draft, updated with each validation arc · **Version 2.0** · October 2026
**Repository:** `dev-abhinavmishra/BKPyV-VCM` · branch `devin/novel-extensions` (pending merge)
*Every quantitative claim is traceable to a script + test in the repository; digitized-data precision and phenomenological assumptions are labelled where they apply.*

---

## Structured Abstract

**Background.** BK polyomavirus (BKPyV) reactivates in 30–60% of kidney-transplant recipients and causes polyomavirus-associated nephropathy in 1–10%, a leading cause of early graft loss. No antiviral has proven efficacy; management is qPCR surveillance followed by immunosuppression reduction, which trades viral control against rejection. Published mechanistic models reproduce population-level kinetics but emit no individual-level forecasts and no derived intervention policy.

**Methods.** We built a 23-compartment mechanistic ODE of BKPyV infection spanning the graft and urinary compartments, the intracellular S-phase–gated replication program, the innate and BKPyV-specific adaptive immune arms, real pharmacokinetic (ng/mL) dosing for tacrolimus and sirolimus, a latent reservoir, cell-to-cell transmission, and a two-pool NCCR quasi-species. We added four prediction layers — stochastic reactivation onset, early-window (β,δ) inference with Metropolis posteriors, continuous optimal control of the conversion decision, and CYP3A5 pharmacogenomic stratification — then validated against six published kinetic signatures and six real patient plasma trajectories digitized from Funk et al. 2008 under an honest holdout protocol (train on the first 60% of each series, predict the rest).

**Results.** All six benchmark signatures pass. The model emits eight emergent predictions, including an onset-time *distribution* (median ≈6 wk at tacrolimus 8 ng/mL, within the clinical 4–16 wk cluster), an optimal-control derivation that immediate tac→sirolimus conversion dominates tapering, and a ~2.7× plasma-over-urine enrichment of rearranged NCCR genomes consistent with Gosert et al. On real patient data the model's held-out predictions beat a per-patient baseline on 5 of 6 patients (cohort MAE 0.82 vs 1.36 log₁₀ cp/mL; best per-patient R² = 0.98), with failure modes that map to specific, named mechanisms.

**Conclusions.** A mechanistic virtual cell can move three unresolved BKPyV questions — when viremia begins, when it clears, and which intervention is optimal — from reactive to predictive, with quantified honest limits.

**Keywords:** BK polyomavirus; kidney transplantation; viral kinetics; mechanistic model; tacrolimus; sirolimus; NCCR; immunosuppression; digital twin

---

## 1. Introduction

BK polyomavirus is a ubiquitous human polyomavirus that establishes lifelong latency in renal tubular epithelium after asymptomatic childhood infection [1]. Under transplant immunosuppression it reactivates in 30–60% of recipients; unchecked replication produces polyomavirus-associated nephropathy (PVAN) in 1–10%, historically causing graft loss in up to half of affected patients [2,3]. There is no antiviral with proven efficacy against BKPyV: consensus management is scheduled qPCR surveillance with reduction of immunosuppression at defined viral-load thresholds — an intervention that exchanges viral control for acute-rejection risk [2,4].

Three questions dominate clinical practice and lack quantitative answers:

1. **Onset** — when will a given recipient develop detectable viremia, and with what probability?
2. **Clearance** — once immunosuppression is reduced, when will the virus clear?
3. **Policy** — which immunosuppression change (dose taper vs agent conversion) clears fastest at least rejection risk?

Mechanistic viral-kinetic models are the natural tool for such questions. Funk et al. built the field's foundational BKPyV models — a six-compartment plasma model validated against prospectively monitored patients [5,6] — and Gosert et al. described the emergence of rearranged non-coding control region (rr-NCCR) variants during replication [7]. These models reproduce *population-level* kinetics but do not emit individual forecasts, do not represent the drug–immune interaction at the T-cell level, and do not derive the optimal intervention. Single-cell work has since established that host S-phase entry gates large T-antigen expression in primary renal proximal tubular cells [8], supplying a mechanistic trigger prior models lacked.

We therefore built a 23-compartment mechanistic virtual cell designed from the start to emit *falsifiable, individual-level predictions*, and then tested it against real patient trajectories — reporting failures alongside wins.

## 2. Methods

### 2.1 Model structure

The core is a 23-dimensional ordinary-differential-equation system (`src/vcm/simulators/ode_system.py`; LSODA solver, per-day rates) organized in five layers:

- **Graft compartment** — healthy tubular cells C, productively infected I, dead D, free virions V, and a latent reservoir L (a fraction `latent_fraction` = 0.02 of new infections aborts into latency; reactivation flux `lat_reactivation·L·(1 − tac_immune_effect)` feeds I).
- **Urinary compartment** — urothelial C_u/I_u/V_u with kidney→urine drainage and bidirectional cross-feeding (topology per Funk [5]).
- **Intracellular program** — T-antigen T, early-gene/capsid expression, cell-cycle state CC, DNA synthesis. T-antigen production is Hill-gated by host S-phase permissiveness AND DNA synthesis [8]; T-antigen *never* drives cell-cycle entry — an explicit invariant pinned by tests.
- **Immune arm** — nonspecific effector E and interferon; BKPyV-specific naïve→effector T cells with priming, antigen-driven expansion under a homeostatic ceiling (`tcell_carry`), and calcineurin-gated function.
- **Pharmacology** — first-order trough→effect kinetics for tacrolimus and sirolimus driven by real ng/mL schedules ({`start`, `stop`, `trough_ng_ml`}). Tacrolimus suppresses T-cell priming/expansion and nonspecific control; sirolimus additionally removes S-phase permissiveness. Both drugs act on viral *production* (immune control and permissiveness respectively), never on intracellular genome *replication* — a documented mechanistic convention.
- **Persistence and quasi-species** — a V-independent cell-to-cell transmission channel (`c2c_rate` = 0.03) plus two-pool NCCR quasi-species: archetype and rearranged genomes compete in kidney and urinary pools, with spontaneous rearrangement during replication, an early-gene overexpression advantage, and a capsid-suppression fitness cost (`rr_capsid_virion_cost` = 0.4 partial coupling).

Copies/mL readouts use an explicitly assumption-labelled anchor bridge (`viral_load_mapper.py`; 10³ screening / 10⁴ presumptive-PVAN thresholds). The system is pinned by 219 tests (8 skipped for optional dependencies).

### 2.2 Prediction machinery

- **Stochastic reactivation onset** (`reactivation_onset.py`): per-cell Poisson reactivation hazard scaled by tacrolimus trough → Monte-Carlo onset-time distribution (150 virtual patients per trough).
- **Early-window forecast** (`early_forecast.py`): virtual cohort jittered in (β, δ, burst size, inoculum) with 0.25-log₁₀ qPCR noise; (β,δ) inferred from the first three weekly points; Metropolis posterior propagated through a conversion intervention → clearance-week distribution with 90% credible intervals and P(clear).
- **Optimal control** (`optimal_control_analysis.py`): (conversion day, sirolimus trough) as continuous controls; Nelder-Mead on a clearance+rebound objective.
- **Pharmacogenomics** (`pharmacogenomic_stratification.py`): CYP3A5 expressor status → tacrolimus clearance rate per CPIC [9] → onset/risk stratification.
- **VST digital twin** (`tcell_therapy_simulation.py`): exogenous cytotoxic T-lymphocyte bolus under maintained immunosuppression.
- **Real-data fitting** (`fit_patient_series.py`): five-parameter fit (β, δ, pre-roll offset, intervention week, post-IS trough) per patient; **honest holdout** — fit the first 60% of each series, predict the remaining 40%, scored against a per-patient mean baseline.

### 2.3 Validation layer

`benchmark_viral_load.py` checks six signatures against published observations: Funk half-life and curtailment thresholds [5,6], Gosert plasma>urine rr enrichment [7], Hirsch threshold ordering [3], and the urine:plasma load ratio. `identifiability_analysis.py` bounds what patient data could ever constrain. `validate_bkpyv.py` / `validate_isef.py` gate the pipeline end-to-end.

### 2.4 Real patient data

No public longitudinal qPCR dataset is openly downloadable, so six patient plasma trajectories were digitized from Funk et al. 2008 Fig 3 [6] (CC BY-NC-ND; pixel-calibrated against the figure's own log₁₀ 4/7 reference lines; precision ≈ ±0.15 log₁₀, ±0.5 wk; provenance in `validation_data/funk2008/README.md`).

## 3. Results

### 3.1 Published signatures reproduce without fitting

All six benchmark checks pass: clearance half-life 6.6 d (within Funk's 0.25–17 d lytic regime [5,6]); curtailment profile (50% ineffective, ≥80% ≤7 wk, ≥90% ≤3 wk [6]); urine:plasma load ~450× (same order-of-magnitude direction as Funk's aggregate ~3000×; residual gap documented); rr-NCCR plasma>urine enrichment 2.6–2.7× (Gosert direction [7]); screening-threshold ordering (Hirsch [3]).

### 3.2 Eight emergent predictions

1. **Onset is a distribution, not a date.** Median ≈6 wk at tacrolimus 8 ng/mL (p10–p90 ≈1–17 wk), inside the clinical 4–16 wk cluster; reactivation probability rises 57%→96% across troughs 3→8; higher troughs yield earlier *and tighter* onsets.
2. **Individual clearance is forecastable from three early points.** Virtual-cohort MAE 0.17 wk, 0/24 clear/not-clear discordance; posterior CIs answer "clear by week X with 90% probability."
3. **The optimal policy is derivable.** The continuous-control optimum sits at the feasibility boundary — convert *immediately* at the strongest tolerated mTOR signal; each week of delay costs ~2.4 objective-weeks. Conversion dominance is a derived consequence, not a grid-search artifact.
4. **Genotype orthogonalizes risk.** CYP3A5 expressors clear tacrolimus ~2× faster (CPIC [9]); genotype-blind dosing assigns non-expressors the viral-onset risk and expressors the rejection risk.
5. **rr-NCCR segregates by compartment.** Emergent ~2.7× plasma-over-urine enrichment — the direction Gosert reported clinically.
6. **Monitoring cadence is itself an intervention.** Acting at 10³ vs 10⁴ cp/mL changes both clearance timing and quasi-species emergence.
7. **Persistence has two channels.** Cell-to-cell spread sustains infection under extracellular clearance; the latent reservoir resurges when immunosuppression is interrupted mid-clearance.
8. **VST durability is gated by residual immunosuppression** (honest negative result): a single CTL bolus produces only transient dips — and only at ~4× the homeostatic ceiling — matching the documented clinical transience of VST responses and naming `tcell_carry` as the parameter to probe with product data.

### 3.3 Real patient data: honest holdout

**Table 1.** Full-series fits and 60/40 holdout scores on six digitized Funk-2008 plasma trajectories.

| Patient | Full RMSE | Full R² | Holdout MAE | Baseline MAE | Verdict |
|---------|-----------|---------|-------------|--------------|---------|
| 00107   | 0.47      | 0.53    | **0.29**    | 0.62         | model wins |
| 00885   | 1.07      | 0.17    | 1.49        | 1.48         | tie (refractory) |
| 06119   | 1.68      | −0.40   | **0.43**    | 2.03         | holdout wins big |
| 15207   | 0.06      | **0.98**| **0.42**    | 0.68         | model wins |
| 19771   | 0.92      | 0.53    | **1.69**    | 2.74         | model wins |
| 27447   | 0.59      | 0.52    | **0.61**    | 0.61         | model wins (marginal) |
| **Cohort** | —      | —       | **0.82**    | **1.36**     | **model wins 5/6** |

Documented failure modes (reported, not tuned away): Pat 06119's ~10⁷·⁷ plasma peak exceeds the model's current ~10⁶ dynamic ceiling, so its full-series R² is negative even though its holdout beats baseline 4.7×; Pat 00885's refractory fluctuation has no monotone equivalent — the model declines to overfit it. Two failure modes are *mechanism pointers*: the transient surges the fitted curves cannot reach (00107's bump, 06119's spike, 27447's rebound) are exactly what the latent-reservoir and stochastic-reactivation arms generate as discrete events — naming the next experiment (fit with release events enabled; obtain urine + PBMC series to constrain L).

### 3.4 Supporting biology

GSE317012 single-cell reanalysis (26 biopsies; 34,987 QC-passing cells): descriptive immune-fraction and pathway trends consistent with the modeled immune arm (MWU p = 0.16 — reported as descriptive, not confirmatory).

## 4. Discussion

**What is new.** Prior BKPyV models validated population kinetics; this one was engineered to emit individual-level falsifiable predictions and was then exposed to real patient trajectories with the failures published alongside the wins. The strongest result for a specialist audience is not the fit quality but the *shape* of the failure set: the model names exactly which phenomena it cannot see (peak ceiling, refractory fluctuation, reservoir-driven rebound), and each maps to a measurable biological quantity.

**Clinical translation.** The digital-twin demo (three early qPCR values → clearance-week posterior) is a concrete path toward preemptive management; the pharmacogenomic map identifies which patients a single tacrolimus trough actually stratifies; the optimal-control result reframes "reduce immunosuppression" as a solvable timing-and-agent problem. The rr-plasma>urine prediction is directly testable on paired urine/plasma longitudinal series — the dataset our collaborating lab (S. Thompson, UAB) has offered.

**Limitations** (full list in `docs/LIMITATIONS.md`): phenomenological calibration; digitized-data precision ±0.15 log₁₀; plasma dynamic ceiling ~10⁶; a six-patient single-figure validation set; VST in effector-fraction rather than product-dose units; the latent reservoir lacks a clinical readout. Every limitation has a named data source that would close it.

## 5. Conclusions

A mechanistic virtual cell can move BKPyV management from reactive to predictive: onset distributions, individual clearance forecasts with credible intervals, and a derived — not searched — intervention policy. Six-of-six published signatures reproduce, real-patient holdout beats baseline on 5/6, and every failure mode names its own missing measurement: a model that knows what it knows.

## References

1. Hirsch HH, et al. Polyomavirus-associated nephropathy in renal transplantation. *Transplant Rev* / consensus reviews; see also AST IDCOP [3].
2. Kotton CN, et al. Updated international consensus guidance on BKPyV management in kidney transplantation. *Transplantation* 2024;108(9):1834–1866.
3. Hirsch HH, Randhawa P; AST Infectious Diseases Community of Practice. BK polyomavirus in solid organ transplantation — guidelines. *Clin Transplant* 2019;33:e13528.
4. Fang Y, et al. Dynamic risk prediction of BKPyV reactivation after renal transplantation. *Front Immunol* 2022;13:971531.
5. Funk GA, et al. BK polyomavirus dynamics under immunosuppression; persistence and clearance modeling. *J Infect Dis* 2006;193:80–87.
6. Funk GA, Gosert R, Comoli P, Ginevri F, Hirsch HH. Polyomavirus BK replication dynamics in vivo and in silico to predict cytopathology and viral clearance in kidney transplants. *Am J Transplant* 2008;8(11):2368–2377. (Fig 3 = digitized patient series, CC BY-NC-ND.)
7. Gosert R, et al. Rearranged-NCCR emergence and plasma/urine compartment differences. 2008. PMC2292223.
8. Needham CD, Perritt E, Thompson S, et al. Single-cell analysis reveals host S phase drives large T antigen expression during BK polyomavirus infection. *PLoS Pathog* 2024;20(12):e1012663.
9. Birdwell KA, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) guideline for CYP3A5 genotype and tacrolimus dosing. *Clin Pharmacol Ther* 2015;98(1):19–24.
10. Myers et al. Six-state patient-scale BKPyV ODE calibrated on Duke/CTOT data. *Viruses* 2025;17(1):50. (Positioning: this project's novelty is the cell-mechanism layer, not "first BKPyV model.")
11. Hirsch HH, et al. Sirolimus early-phase inhibition of BKPyV in renal tubular cells. *Am J Transplant* 2016;16(3):821–832.
12. GSE317012 single-cell kidney-biopsy cohort (JCI Insight 2026, doi:10.1172/jci.insight.198227).
13. Thompson S. Personal communication, UAB, Aug–Oct 2026 (model conventions confirmed; longitudinal patient series offered).

## Figures

- **Figure 1** (`outputs/figures/composite_figure.png`): four-panel composite — (A) two-compartment kinetics vs clinical thresholds; (B) rr-NCCR emergence, plasma > urine; (C) taper plateaus vs conversion clears; (D) reactivation-onset distributions by tacrolimus trough.
- **Figure 2** (`outputs/figures/patient_fits.png`): six-panel model fits vs digitized patient series; shaded region = held-out data never seen by the fit.

---

### Version log
- **v1.0** (Sep 2026) — first full draft: 22-dim model, emergent results 1–6, no real data.
- **v2.0** (Oct 2026) — 23-dim (c2c + L); emergent results 7–8; real-data holdout (Funk digitized, post-audit corrected values); packaging arc.
- *Next update planned:* Thompson longitudinal urine+plasma series fit (pipeline built; data pending).
