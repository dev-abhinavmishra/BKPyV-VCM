# Presentation Points — judge-facing narrative

## The one-sentence pitch
A mechanistic virtual cell that doesn't just *simulate* BK polyomavirus —
it *predicts*: when viremia starts, when a patient clears after
intervention, and which intervention is optimal — validated against real
transplant-patient trajectories with its failures published, not hidden.

## The three-slide story arc
1. **Problem**: no effective antiviral exists; management is "watch the
   viral load rise, reduce immunosuppression, hope" — and every reduction
   trades viral control against graft rejection.
2. **Model**: 23-compartment mechanistic ODE — intracellular replication
   gated by host S-phase, T-cell immune control, plasma–urine coupling,
   latency + cell-to-cell persistence, real PK drug dosing.
3. **Predictions → validation**: eight falsifiable signatures; all six
   published benchmarks pass; real-patient holdout beats baseline
   (0.91 vs 1.36 log10) — with two named, instructive failures.

## Numbers that stop judges (all verified, all in-repo)
- Onset is a **distribution**: median 6.4 wk at tac 8 (clinical: 4–16 wk).
- Clearance forecast from **three early measurements**: ~1.2-day median
  error + 90% credible interval → "clear by week X, 90% confidence."
- Optimal policy is **derived**, not searched: immediate tac→sir
  conversion; delay costs ~2.4 wk/week.
- Real data: **R² 0.96** on patient 15207; cohort holdout MAE **0.91 vs
  1.36** baseline on digitized Funk-2008 series.
- Honest negative: **VST fails** under maintained tacrolimus except
  transiently — matching why clinical VST responses are often transient.
- Genotype map: CYP3A5 expressors vs non-expressors carry *opposite*
  risks under genotype-blind dosing.

## Questions judges ask — and the answers
- *"Is this just curve-fitting?"* No — predictions were emitted before
  being checked; holdout validation uses data the fit never saw;
  failure cases are published with structural causes.
- *"Where's the real data?"* Six digitized patient series (Funk 2008,
  CC license) in `validation_data/`; a Thompson longitudinal series is
  the designed next dataset — the pipeline (`fit_patient_series.py`)
  is built and waiting.
- *"What can't it do?"* Documented: plasma dynamic ceiling (~6 log),
  refractory-fluctuation patients, rebound timing without a reservoir
  readout — each limitation names the data that closes it.
- *"Why compartments?"* The kidney/urine split isn't cosmetic — it
  reproduces the clinical plasma>urine NCCR enrichment (2.7×) as an
  *emergent* property.

## Stage checklist
- Composite figure: `outputs/figures/composite_figure.png` (4-panel).
- Streamlit "Findings" page: live posterior demo + onset explorer —
  the 60-second judge demo.
- Poster story order: problem → model schematic → emergent predictions
  → real-data validation → named limits (the honesty panel wins judges).
