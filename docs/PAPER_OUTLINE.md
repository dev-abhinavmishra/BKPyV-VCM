# Paper Outline — BKPyV Mechanistic Virtual Cell

**The full draft lives in `docs/PAPER_DRAFT.md`** (competition draft,
all numbers verified against `devin/novel-extensions` runs). This file
keeps the section map; edit the draft, not this.

1. **Introduction** — PVAN burden; reactive-vs-predictive gap; three
   unanswered clinical questions (onset, clearance, policy); what prior
   mechanistic models (Funk 2008) did and didn't emit.
2. **Methods** — 23-dim ODE (kidney+urinary compartments, T-cell arm,
   PK dosing, c2c + latent channels); prediction machinery (onset MC,
   posterior forecast, optimal control, genotype map, real-data fit,
   VST twin); validation layer (benchmark, identifiability, validators).
3. **Results** — 6/6 published signatures pass; 8 emergent predictions;
   real-data validation table (cohort holdout MAE 0.91 vs 1.36);
   supporting scRNA trends (descriptive).
4. **Discussion** — novelty = individual falsifiable predictions +
   published failure modes; translation path; limitations each tied to
   a named closing dataset.
5. **Conclusion** — predictive in-silico framework for preemptive,
   personalized BKPyV management.
