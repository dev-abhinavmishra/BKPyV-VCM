# ISEF Abstract

BK polyomavirus (BKPyV) reactivation causes graft-damaging nephropathy in
immunosuppressed kidney-transplant recipients, yet the intracellular
mechanisms linking immunosuppressive drugs to viral production remain
opaque. This project builds a mechanistic ordinary-differential-equation
model of BKPyV infection in renal tubular epithelial cells: viral T-antigen
accumulation gated by host S-phase entry, genome replication, capsid
production, and drug-modulated host permissiveness. The simulator
(`bkpyv_ode`, per-day rates) reproduces the literature direction of
immunosuppression effects — tacrolimus raises modeled peak viremia while
sirolimus lowers it (clinical-bridge peaks: 6.77M vs 31k copies/mL vs 2.43M
untreated; `clinical_validation_report.json`). Archetype versus rearranged
NCCR genome scenarios expose an early-gene/capsid trade-off
(`outputs/bkpyv_review/scenario_summary.csv`). To ground the model in real
tissue biology, the kidney-biopsy single-cell series GSE317012 was
reanalyzed sample-by-sample (26 biopsies: 12 control, 5 peaking, 9
resolving; 34,987 of 294,286 cells passing stringent QC;
`data/processed/gse317012_pathway_scores.csv`), showing descriptive
condition-level trends in immune-fraction and pathway scores consistent
with the model's mechanisms. The model is a hypothesis generator — a
testable map from drug mechanism to intracellular viral dynamics — not a
patient-level predictor. All parameters remain phenomenological or
literature-anchored; scRNA findings informed no coefficient. This frames
testable questions for timing and choice of immunosuppression.

*Word count target: <=250 (ISEF 2027 rule 11: maximum 250-word, one-page
abstract — societyforscience.org/isef/international-rules/rules-for-all-projects/).*
