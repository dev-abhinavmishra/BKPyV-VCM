# ISEF Research Plan — BKPyV Mechanistic Cell Model

## Rationale

BK polyomavirus (BKPyV) reactivation under immunosuppression causes
BKPyV-associated nephropathy in up to ~10% of kidney-transplant recipients
and is a leading cause of graft loss. Clinical management trades off
rejection prevention against viral control, but the intracellular mechanism
by which specific immunosuppressants (tacrolimus vs sirolimus) shift viral
production is not directly observable in patients. A mechanistic
computational model can make those hidden dynamics explicit and generate
testable hypotheses about drug timing and choice.

## Research Question

How do host-cell state (cell-cycle phase), viral genome architecture
(archetype vs rearranged NCCR), and immunosuppressive drug mechanism jointly
determine intracellular BKPyV production dynamics, and can a mechanistic ODE
model reproduce the qualitative clinical ordering observed in the
literature?

## Hypothesis and Expected Outcomes

Hypothesis: immunosuppressants differ in viral outcomes primarily through
distinct mechanistic entry points — tacrolimus through weaker T-cell immune
targeting (higher production) and sirolimus through an mTOR-dependent
S-phase gate on host permissiveness (lower production). Expected: the model
orders tacrolimus > untreated > sirolimus on peak viral production, shows a
graded early-window drug-timing effect, and exhibits an early-gene/capsid
trade-off under NCCR rearrangement.

## Materials

- This repository: `src/vcm/` package, canonical `bkpyv_ode` simulator,
  YAML experiment configs (`configs/bkpyv_ode_*.yaml`)
- Public data: GEO GSE317012 (BKPyV kidney-biopsy scRNA-seq, JCI Insight
  2026), microarray series GSE47199 / GSE75693 / GSE72925 (series-matrix
  metadata; no processed value tables present)
- Python >=3.10 with `.venv`; pinned dependencies in `pyproject.toml`
- No wet-lab work; no human participants recruited for this project

## Procedures

1. Implement the mechanistic ODE (T-antigen → early gene → replication →
   capsid/production; immune-control and drug terms) with per-day rates.
2. Run baseline/infection/tacrolimus/sirolimus/NCCR-variant configs;
   verify trajectories are finite, non-negative, and directionally ordered
   (`validate_bkpyv.py`, `validate_isef.py`, `validate_clinical.py`).
3. Map model viral load to indicative copies/mL through the
   assumption-labelled piecewise log-linear bridge anchored to consensus
   screening/treatment thresholds (AST IDCOP 2019; Kotton 2024).
4. Reanalyze GSE317012: traversal-safe extraction, authoritative series-
   matrix metadata, sparse per-sample QC (34,987/294,286 cells retained),
   marker-proxy fractions, five curated pathway scores, and sample-level
   contrasts (`scripts/analyze_gse317012.py`; `docs/GSE317012_ANALYSIS.md`).
5. Compare model directions against literature-reported directions only —
   never patient-level calibration.

## Risk and Safety

Computational-only project; no wet-lab, no human participants recruited, no
vertebrate animals. Public de-identified GEO data only. **OPEN ITEM for the
student**: whether secondary analysis of public biopsy-derived datasets
triggers SRC human-participants review must be determined with the regional
fair SRC / adult sponsor before submission — flagged OPEN, not resolved
here.

## Data Analysis

- ODE trajectories: peak, timing, threshold-crossing weeks
  (`clinical_validation_report.json`).
- Sensitivity sweep over key parameters (`outputs/sensitivity/`).
- scRNA reanalysis: sample-level summary + Mann-Whitney/Cliff's-delta
  contrasts; descriptive trends only (peaking n=5 is underpowered — no
  significance claimed); explicit note that no model coefficient was
  updated from GSE317012.
- Synthetic cohort (`data/processed/synthetic_patient_cohort.csv`,
  labelled synthetic) exercises the risk-model plumbing only.

## Bibliography

See [REFERENCES.md](REFERENCES.md) — Needham 2024, Hirsch 2016, Funk 2006,
Weissbach 2024, Kotton 2024, JCI Insight 2026 (DOI
10.1172/jci.insight.198227), AST IDCOP 2019, Fang 2022, Yamauchi 2025, ISEF
2027 rules.

## Acknowledgements and AI-Usage Disclosure

ISEF 2027 rules (Checklist 1A and the international rules book) require
disclosure of assistance, including AI tools, used in the project.

**Honest status**: this repository was developed with AI coding assistance
(including the Devin agent) under the student's direction. The authoritative
record of human direction, AI prompts, and tool usage must be assembled by
the student from actual project history — a real prompt/usage log is NOT
reproduced here and must not be fabricated for the form. Placeholder for the
student to complete with their true disclosure text:

> [STUDENT TO COMPLETE: describe AI-tool usage truthfully per the 2027
> checklist instructions, and attach the genuine usage record.]
