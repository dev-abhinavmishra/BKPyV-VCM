# Presentation Points — Judge-Facing Talking Points

## The 30-second pitch

- BKPyV damages kidney grafts under immunosuppression; the drug→virus
  mechanism inside the cell is invisible to clinicians.
- I built a mechanistic ODE of intracellular BKPyV (host S-phase gate,
  T-antigen, replication, capsid production, immune control) and reanalyzed
  real biopsy single-cell data for descriptive context.
- It is a hypothesis generator, not a clinical predictor — every clinical
  number is explicitly assumption-labelled.

## Numbers you can quote (each traced to a fresh output)

- **184 tests pass** — `outputs/test_results.txt`.
- **5 plugins / 8 configs verified** — `vcm list-plugins`;
  `harness_check_report.json`.
- Drug ordering: tacrolimus peak 6.77M > untreated 2.43M > sirolimus 31k
  copies/mL — `clinical_validation_report.json` (assumption-labelled
  bridge, not patient calibration).
- Drug switch tac→siro at week 4/8 collapses week-52 viremia to ~0 —
  `outputs/clinical/drug_switch_results.json`.
- NCCR rearranged vs archetype: production 3.01 vs 1.93 peak rate —
  `outputs/bkpyv_review/scenario_summary.csv`.
- GSE317012: 26 biopsies (12 Control / 5 Peaking / 9 Resolving),
  34,987/294,286 cells post-QC — `data/processed/gse317012_pathway_scores.csv`.
- Immune-proxy fraction trend +0.187 at peaking (p=0.16, descriptive) —
  `data/processed/gse317012_contrasts.csv`.

## If asked about rigor

- Validators fail nonzero on regression (`validate_*.py` exit codes in
  sprint evidence).
- QC retention is stated honestly (11.9%; MT-heavy biopsy tissue).
- Model parameters are phenomenological/literature-anchored — the scRNA
  reanalysis informed NO coefficient (documented freeze rationale).
- The synthetic cohort is labelled synthetic; it exercises plumbing only.

## If asked about scope limits

- No patient-specific prediction; copies/mL are bridge-labelled, not
  calibrated.
- Microarray series inventoried as metadata-only (no value tables present);
  raw-CEL analysis explicitly out of scope.
- Cell-level scRNA stats never treat cells as independent patients —
  contrasts are sample-level across the 26 biopsies.
