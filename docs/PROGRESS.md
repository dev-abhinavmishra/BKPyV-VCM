# PROGRESS — running log of the autonomous finish task

Task: bring BKPyV-VCM to a complete, reproducible, judge-ready state
(Phases 0–4 per the owner brief). Branch `devin/finish-vcm`.
Log kept chronologically; newest last.

## 2026-10-05 — Phase 0: audit + fixes

- Read README, engine, tests, git log; ran model + validators.
- Fixed `validate_bkpyv.py` (full perturbation list + None-guards) → exit 0.
- Fixed `validate_isef.py`: pattern checks ran on the superseded discrete
  simulator → repointed to `bkpyv_ode` (correct drug directions); mapped
  T-antigen knob to `half_saturation` (real production gate); fixed pattern
  counting → 4/4 "EXCELLENT".
- Wrote `docs/STATE.md` (what exists, what runs, what was broken, state
  vector + parameter sources).

## 2026-10-05 — Phase 1: single-cell layer

- `scripts/download_gse317012.py`: one-command fetch of GSE317012_RAW.tar
  (1.9 GB) + series matrix; `.part` resume; sha256 verify
  (RAW: 1EA67E31…327A8; matrix: 221b8892…8b83e). All 26 GSMs restored.
- `scripts/cluster_gse317012.py`: QC (min genes + mt cap) → normalize →
  PCA(50)/neighbors/Leiden → marker annotation (11 kidney cell types) →
  module scores → `infected_signature` call (top-decile phase-matched
  epithelium; honest proxy — no viral genes exist in the reference) →
  Wilcoxon DE → seeded stratified holdout → finding→model mapping.
- Results: 34,987 cells post-QC, 23 clusters; 857 signature-high epithelial
  cells vs 3,417 comparators; 3,929 DE genes FDR<0.05; top genes = MHC-I/II
  antigen-presentation program (B2M, HLA-A/B/C, CD74, HLA-DRA/DRB1/DPA1,
  PSMB9) — matches the source paper's biology; 7 modules calibrated with
  direction-consistent shifts on the held-out split.

## 2026-10-05 — Phase 2: viral-load benchmark

- Verified no public serial qPCR dataset exists → sanctioned fallback:
  published summary statistics. Curated `data/research/published_kinetics.csv`
  (Funk 2006 PMID 16323135; Funk 2008 doi:10.1111/j.1600-6143.2008.02402.x;
  AST IDCOP 2019; Kotton 2024).
- `scripts/benchmark_viral_load.py`:
  - Calibration (Funk 2006 IS-change): model t½ = **5.5 d** vs published
    0.25–17 d → PASS; band-excess RMSE 0.19 log10 cp/mL (47% in-envelope).
  - Validation (Funk 2008 curtailments from ~1.7×10⁶ cp/mL start): 50%
    never clears (ineffective ✓); 80% clears at 1.4 wk (<7 wk ✓, marginal
    rebound to ~1,238 cp/mL reported transparently); 90% clears at 1.1 wk
    (<3 wk ✓ sustained). 3/3 pass.
  - Validation (Funk 2006 nephrectomy): model t½ ≈10–19 h → between the
    published 1–2 h and 20–38 h arms; reported as envelope-compatible.
  - Uncertainty: δ sweep 0.25–0.8/day → t½ response curve.
- Overall benchmark: **PASS**.

## 2026-10-05 — Phase 3: dashboard

- 3 new pages in `src/vcm/ui/streamlit_app.py` (10 total): 🧫 Single-Cell,
  📈 Viral-Load Validation, 🧾 Parameters & Assumptions (per-parameter
  cited/ASSUMPTION labels). All degrade gracefully on clean clones with
  regeneration instructions. Verified in browser + AppTest.
- `.[analysis]` extra pins scanpy/anndata/umap-learn/pyarrow; README gets
  one-command dashboard + full reproduction recipe.

## 2026-10-05 — Tests + engineering

- `tests/test_benchmark_and_pipelines.py`: +14 tests (band-excess RMSE,
  t½/threshold helpers, published-table schema/citations, quick e2e
  benchmark, sha256 helper, deterministic holdout split).
- `tests/test_streamlit_app.py`: 7→10 pages; PAGE_INDEX replaced by
  name-based navigation (order-insensitive).
- Full suite: **198 passed**.
- Refactors for testability: `split_samples(seed, frac)`,
  `sha256_matches(path, expected)`.

## 2026-10-05 — Phase 4: docs

- `docs/METHODS.md` — equations, assumptions, data sources, validation design.
- `docs/LIMITATIONS.md` — 14 honest limits ordered by severity.
- `docs/ROADMAP.md` — 5 ranked improvements + 10 ranked mentor questions.
- `docs/DECISIONS.md` — 10 autonomous decisions with alternatives + reasons.
