# DECISIONS — choices made autonomously (Phase 0–4)

Each entry: decision, alternatives considered, reason. Ordered chronologically.

## D1 — Canonical engine is `bkpyv_ode`; legacy simulator deprecated for analysis
- **Decision**: All validation, benchmark, and dashboard analysis uses the
  15-dim ODE (`src/vcm/simulators/ode_system.py`). The legacy discrete
  `bkpyv_specific` simulator stays in the repo (not deleted, per rules) but
  is documented as superseded.
- **Alternatives**: (a) keep validating against the legacy sim; (b) fix the
  legacy sim's drug directions.
- **Reason**: The legacy sim inverted/suppressed the drug directions the
  published biology requires; the ODE already reproduces them (baseline
  peakV 1.924 vs tacrolimus 2.180 vs sirolimus 1.118). Rewiring the legacy
  engine would duplicate the ODE. This was also the repo's own documented
  direction (README calls the ODE the canonical engine).

## D2 — Infected cells identified by host-signature proxy, not viral reads
- **Decision**: In GSE317012, "infected" tubular epithelial cells = top-decile
  composite host viral-response score within phase-matched epithelium;
  label is `infected_signature` everywhere (never "infected" alone).
- **Alternatives**: (a) map BKPyV reads — impossible, reference has no viral
  genes (verified: features.tsv.gz contains human gene symbols only);
  (b) call whole Peaking-phase epithelium "infected" — too coarse, mixes
  bystanders; (c) skip the infected/bystander analysis — the task requires it.
- **Reason**: Signature-proxy is the only honest route on this data. All
  outputs and the dashboard state the proxy limitation explicitly. The DE
  result (MHC-I/II antigen-presentation program: B2M, HLA-A/B/C, CD74,
  HLA-DRA/DRB1/DPA1, PSMB9) matches the source paper's biology, supporting
  the proxy's face validity.

## D3 — Phase 2 uses published summary statistics, not patient-level data
- **Decision**: Benchmark the ODE against summary statistics from Funk 2006
  (PMID 16323135) and Funk 2008 (doi:10.1111/j.1600-6143.2008.02402.x),
  curated into `data/research/published_kinetics.csv`.
- **Alternatives**: (a) find a public serial qPCR dataset — searched GEO,
  dbGaP, SRA, Dryad/figshare mirrors: none is public + de-identified;
  (b) digitize published curves — the Funk papers report summary statistics
  in text, not extractable patient trajectories; (c) synthetic data —
  forbidden (do not fabricate).
- **Reason**: Task explicitly sanctions this fallback ("benchmark against
  published summary curves and say so explicitly"). The limitation is stated
  in the report JSON, the dashboard page, and docs/LIMITATIONS.md.

## D4 — Calibration = verify δ containment; validation = disjoint subsets
- **Decision**: Calibration subset = Funk 2006 IS-change arm: check that
  model's post-intervention clearance t½ (5.5 d) lands inside the published
  0.25–17 d range (δ=0.4/day was already chosen inside it — documented as
  calibration, not refit). Validation subsets = Funk 2008 curtailment
  statements + Funk 2006 nephrectomy arms — never used for parameter choice.
- **Alternatives**: (a) refit δ to the published midpoint — pointless, the
  "ground truth" is a range not a point; (b) fit all parameters — no
  patient-level data to fit against.
- **Reason**: With range-valued ground truth, containment is the correct
  criterion; refitting to a fabricated midpoint would be calibration theater.
  For the cell-state layer, a seeded stratified holdout (≈1/3 samples per
  phase) is the calibration/validation boundary, recorded in
  `gse317012_calibration_holdout.csv`.

## D5 — Intervention experiments start from peak viremia, not equilibrium
- **Decision**: Curtailment/nephrectomy benchmarks start from the
  peak-viremia state under tacrolimus (≈1.7×10⁶ cp/mL via the bridge),
  not the t→∞ equilibrium (≈6×10³).
- **Alternatives**: equilibrium start — unrealistically low; clinical
  interventions happen at 10⁵–10⁷ cp/mL.
- **Reason**: Peak state is the only model state in the clinically
  representative intervention range. Reported transparently in the report's
  `starting_state` block.

## D6 — New Streamlit content added as pages in the existing app
- **Decision**: Single-Cell, Viral-Load Validation, and Parameters &
  Assumptions added as sidebar pages of `streamlit_app.py`; they read
  precomputed artifacts and degrade gracefully (instruction banner) on a
  clean clone.
- **Alternatives**: (a) separate app — violates "launch with one command";
  (b) live compute in pages — slow (scanpy pipeline is minutes, 1.9 GB data).
- **Reason**: One-command dashboard stays fast; regeneration instructions
  are shown in place of missing artifacts.

## D7 — Data/outputs stay gitignored; reproducibility is via scripts
- **Decision**: Keep the repo's convention — `data/`, `outputs/`, `*.csv`
  untracked; the scripts are the reproducible path (README documents the
  4 commands). Only code, tests, docs, and the embedded published-kinetics
  table are versioned.
- **Alternatives**: `git add -f` the CSVs — violates repo convention and
  pollutes diffs with generated artifacts.
- **Reason**: "Reproducible from a clean clone" is satisfied by scripts +
  pinned extras; committing generated outputs would be less reproducible,
  not more.

## D8 — δ sweep as the benchmark's uncertainty statement
- **Decision**: Parameter uncertainty reported as pure-clearance t½ vs δ
  sweep (0.25–0.8/day); global OAT analysis delegated to existing
  `scripts/sensitivity_analysis.py`.
- **Alternatives**: (a) full Sobol indices — no fitted posterior to
  propagate; (b) report nothing — task requires uncertainty statement.
- **Reason**: δ is the single calibrated parameter for this benchmark;
  sweeping it answers "how sensitive is the validation claim to the
  calibrated value" directly and cheaply.

## D9 — `split_samples`/`sha256_matches` refactored to be testable
- **Decision**: `split_samples` accepts optional `seed`/`frac` (defaults
  preserve old behavior); `sha256_matches` extracted in the downloader.
- **Alternatives**: test via subprocess — slower, more brittle.
- **Reason**: "tests for every new module" is cleaner with injectable seeds.

## D10 — Nephrectomy check reported as envelope, not hard pass/fail
- **Decision**: Model nephrectomy t½ (≈10–19 h depending on start state) is
  reported against both published arms (1–2 h dense, 20–38 h sparse) and the
  1–38 h union envelope; it does not gate `overall_pass`.
- **Alternatives**: hard-fail against either arm — the two arms disagree by
  ~20× because they measure different things (pure plasma virion clearance
  vs effective dynamics with residual production); a single number can't be
  inside both.
- **Reason**: Honest reporting of a between-arms result is better than
  forcing a pass/fail onto a sampling-regime-dependent quantity.
