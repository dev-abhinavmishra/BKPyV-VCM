# Orchestration Brief — Virtual Cell Model (VCM) to "Finished" State

**To:** GPT-6-Astra (planner/orchestrator)
**From:** Project analysis pass, 2026-09-19
**Repo root:** `C:\Users\pmsma\Downloads\MacBook\virtual-cell-model` (Windows; project was originally developed on macOS — expect `python3` vs `python` and path-separator drift in docs/scripts)

---

## 0. YOUR ROLE — READ FIRST

You are the **planner and orchestrator only**. You do NOT write, edit, or refactor code, configs, tests, docs, notebooks, or scripts yourself, and you do NOT run implementation commands. Your entire job is:

1. Verify the current state (read files, inspect git status — read-only commands are allowed).
2. Produce a written execution plan (`PLAN.md` at repo root is the only file you may author).
3. Decompose the plan into concrete, self-contained task briefs and dispatch **SWE-2 Max subagents** to execute them. Every file creation/edit, every command run, every commit is performed by a subagent, never by you.
4. Review each subagent's output against the acceptance criteria in §6, and redispatch with corrections until each task passes.
5. Run the final acceptance sweep (you may dispatch a verification subagent for this) and report done/not-done per criterion.

### Preflight (do this before anything else; report and STOP if it fails)

- **Read scope**: confirm you can read files under the repo root. If reads fail with a scope/allowed-roots error, request that `C:\Users\pmsma\Downloads\MacBook\virtual-cell-model` be added as an allowed workspace root. If it still fails after the request, STOP and report.
- **Delegation tool**: use whatever subagent/delegation mechanism your environment provides (e.g., `traycer_a2a` or a generic `run_subagent`-style tool), selecting model **SWE-2 Max** for every implementation task. If the delegation tool is not configured OR SWE-2 Max cannot be selected, STOP and report — do **not** substitute another model, do **not** implement anything yourself.
- **Browser/web tool**: needed only for the GEO metadata fetch in T7b. If unavailable, flag it and proceed with the rest — do not block the whole engagement on it.

This is a **one-and-done engagement**: when you finish, the project must be in its final, finished state per §6. Do not produce a partial plan and stop. Do not leave "future work" that is inside the scope defined here.

---

## 1. WHAT THIS PROJECT IS

**Virtual Cell Model (VCM)** — a modular Python platform for mechanistic cell simulation, built as a **Regeneron ISEF 2027** science-fair project (student: Abhinav Mishra, Health Careers High School, San Antonio, TX; category: Computational Biology / Cellular & Molecular Biology). ISEF 2027 is in Los Angeles, May 2027. NOTE: some docs still say "submit to ISEF in 2026" — that is stale; the target is 2027 everywhere.

The flagship module simulates **BK polyomavirus (BKPyV) infection of renal tubular epithelial cells** under immunosuppression (tacrolimus vs sirolimus), grounded in: Hirsch et al. AJT 2016 (in-vitro drug mechanisms), Needham et al. PLoS Pathog 2024 (host S-phase gate precedes large T-antigen), Weissbach et al. J Virol 2024 (single-cell transcriptome signatures), Funk et al. J Infect Dis 2006 (patient clearance kinetics), AST IDCOP 2019 + Kotton 2024 Second International Consensus thresholds.

**Scientific posture (must not regress):** the model is a *mechanistic hypothesis generator*. All clinical-facing copies/mL numbers flow through an explicitly **assumption-labelled anchor bridge** — NOT a calibration to patient data (no longitudinal patient cohort is in the repo yet). Validation = internal consistency checks only. Never fabricate patient data, never remove caveat labels, never upgrade the honesty wording into a validation claim. If real data analysis produces quantitative results, mark their provenance precisely (which dataset, which cell context, in-vivo vs in-vitro).

---

## 2. COLLABORATION CONTEXT (from the student's email thread with Dr. Sunnie Thompson, UAB — the group behind Needham 2024)

Sent 2026-08-08, replied 2026-08-17. This exchange already drove the model card's "domain-expert review notes," and these facts constrain all downstream work:

- **She confirmed the S-phase gate architecture**: no initial host S phase → no TAg → no replication. TAg is NOT the trigger for cell-cycle re-entry (it isn't expressed in G0/G1 when entry would need to happen); the actual driver of re-entry in terminally differentiated RPTECs is **unknown — her lab has a manuscript under review on it**. Never write docs claiming TAg drives `dCC/dt`.
- **Terminology ruling**: *replication* = intracellular genome copying; *production* = system-level virion output. Tacrolimus promotes **production** by weakening immune targeting of infected cells; sirolimus may prevent **production** by blocking S-phase entry/TAg expression. Drug wording everywhere must follow this.
- **Cell-context caveat**: mTOR-axis magnitudes are in-vitro and cell-type-dependent; immortalized lines (HEK293) yield ≥1 log lower titers than primary RPTECs despite constant S phase. Direction-validated, not magnitude-validated.
- **Single-cell ≠ kidney**: RPTEC culture findings are priors, not in-vivo kinetics; her lab is testing carryover to human kidneys.
- **She offered real longitudinal first-year urine AND plasma BKPyV patient series** for model comparison — `scripts/compare_patient_series.py` is built for exactly this and stays data-gated until the CSV arrives.
- **She proposed the collaboration axis**: archetype (transmitted/persistent form) vs **rearranged NCCR** (TAg overexpression, reduced capsid expression — the clinically important form). She wants a model of which factor combinations predict **reactivation vs persistence** in a transplanted kidney. The repo already has `configs/bkpyv_ode_infection_rearranged.yaml` and the `nccr_early_expression_multiplier`/`nccr_capsid_expression_multiplier` parameters — this axis is the project's headline direction; prioritize work that strengthens it.
- She is warm and responsive ("I would be excited to talk"). Deliverables should be drafted so the student can send her an update — accuracy and honest caveats matter more than polish.

---

## 3. VERIFIED CURRENT STATE

### Working / complete
- `src/vcm/` package: `core/models.py` (pydantic domain models), `plugins/` (5 plugins: `bacteria.minimal_cell`, `virus_host.sars_cov2`, `mammalian.immune_tcell`, `mammalian.cancer_cell`, `transplant.bk_polyomavirus`), `simulators/` (mechanistic, ml_based, hybrid, legacy discrete `bkpyv_simulator` [hour units], canonical `bkpyv_ode_simulator` + `ode_system.py` [15-var ODE, day units]), `clinical/` (viral_load_mapper, calibration, risk_prediction, thresholds, external_validation), `analysis/sensitivity.py`, `evaluation/arc_metrics.py`, `experiments/runner.py`, `cli/main.py` (click: run/compare/list-plugins/inspect), `viz/` (plots, bkpyv_review, sensitivity_plots), `ui/` (dashboard.py, streamlit_app.py ~890 lines), `data/loader.py`.
- `configs/` ~20 YAML configs (bkpyv ODE set + legacy + other plugins) + `configs/archive/`.
- `scripts/` maintained pipeline (`generate_bkpyv_review_bundle.py`, `sensitivity_analysis.py`, `simulate_drug_switch.py`, `compare_patient_series.py` [data-gated], `generate_clinical_figure.py`, `generate_risk_dataset.py`, `generate_risk_figures.py`, `save_risk_results.py`, `reorganize_gse317012.py`, `check_output_files.py`, `validate_ode_implementation.py`) + `scripts/attic/` (retired scripts incl. `analyze_gse317012_scrna.py`, `process_gse317012_singlecell.py`, `process_gse75693.py/.R`, `compute_gse75693_pathway_scores.py`, `create_probe_mapping.R`).
- Root validators: `validate_bkpyv.py` (config harness), `validate_isef.py` (mechanism checks), `validate_clinical.py`. On-disk reports all pass.
- `tests/` 16 test files (~163 tests). Canonical docs: `README.md`, `ARCHITECTURE.md`, `docs/BKPYV_MODEL_CARD.md`, `docs/ISEF_PROJECT_OVERVIEW.md`, `docs/ODE_IMPLEMENTATION.md`. `docs/archive/` holds superseded docs with ARCHIVED banners. `.github/workflows/ci.yml` (pytest matrix 3.11–3.13; ruff non-blocking).
- `data/research/extracted_parameters.json` — extracted clinical anchors from 3 papers (Kotton 2024, Utah risk model, Frontiers 2025) with a cross-paper synthesis (median peak ~190,000 copies/mL [IQR 51,400–833,000], median time-to->10k = 104 days [IQR 77–167], urine precedes plasma ~6 wk, tac trough 5–7 ng/mL, expected 10-fold decline ~4 wk post-intervention). Audit note flags abstract-only values as provisional — usable as *reasonableness anchors*, not calibration.

### Data restored on 2026-09-19 (the audit's "unrelated datasets" verdict was WRONG — all are BKPyV-relevant)
- `data/GSE317012_RAW.tar` (1.9 GB) — **the highest-value asset**: 10x-format scRNA-seq of human kidney allograft biopsies, 26 samples × filtered+raw (52 `matrix.mtx.gz`), conditions = uninfected control / **peaking** viremia / **resolving** viremia. Source publication: *JCI Insight* 2026, "Single-cell profiling reveals epithelial and immune responses in BK polyomavirus–infected human kidney biopsies" (doi:10.1172/jci.insight.198227). In-vivo human data mapping directly onto the model's infection trajectory phases.
- `data/GSE47199_RAW.tar` (245 MB) — 57 Affymetrix HuGene-1_0-ST `.CEL.gz` arrays, blood + biopsy from BKV viremia/nephropathy patients.
- `data/GSE75693_RAW.tar` (637 MB) — 79 `.CEL` arrays, urine-cell transcriptome for renal transplant injury (incl. BKVN).
- `GSE72925_RAW.tar` (1.3 GB, still at repo root — move under `data/`) — 168 biopsy `.CEL` arrays (PVAN/TCMR/stable/IFTA), Sigdel et al., Transplantation 2016 (PMID 27140517).
- **Consequence**: the "retired" attic scripts were shelved only because data was absent; the documented #1 next step ("extract quantitative pathway activities from GSE317012, replace heuristic coefficients with data-derived values") is now executable. Also fix `docs/ISEF_PROJECT_OVERVIEW.md` line ~86 which calls these "unrelated" — correct it to: parameters are still phenomenological; datasets are BKPyV-relevant and now analyzed in `data/` pipeline.

### Verified problems (confirmed on this machine today)
1. **Broken local environment → test failures.** `python -m pytest tests/ -q`: **154 passed, 9 failed** — all import failures from broken site-packages, not code bugs: 8 in `tests/test_bkpyv_review.py` (`matplotlib` → `No module named 'dateutil.rrule'`), 1 in `tests/test_runner_and_units.py` (`rich` → `No module named 'rich._emoji_codes'`). pytest itself was missing (just installed); `pydantic` shows version `None` in `pip list` (corrupt metadata). **No venv**; everything is global Python 3.13. Fix = `.venv`, `pip install -e ".[dev,ui]"`, clean reinstalls, re-run to 0 failures.
2. **Entire post-audit overhaul is uncommitted** (~80 modified/deleted/untracked paths). Last commit `ab9dd88`. `origin` IS configured: `https://github.com/website-deployer/BKPyV-VCM.git` — no push attempted; pushing requires credentials and explicit user approval (see §6 non-goals).
3. **Streamlit Comparison page is a placeholder** (`st_comparison_page()`, streamlit_app.py ~line 705). Implement it (compare stored simulation results) or remove cleanly.
4. **`src/vcm/evaluation/` missing `__init__.py`** (works only via namespace fallback; inconsistent with `analysis/`).
5. **Doc/code drift:** README structure lists `src/vcm/utils/` and `data/mock/*.json` files that don't exist (`data/mock/` is empty); `judges_mentors_guide.md` still has pre-review wording (tacrolimus "activates replication via FKBP-12", hard "24h" sirolimus window, unverified "2.0–2.3× risk" claims); `technical_documentation.md` says Python 3.8+ with ad-hoc pip installs; `python3`/Mac paths throughout; "ISEF 2026" → 2027; the "unrelated datasets" claim.
6. **Unverified surfaces:** `notebooks/bkpyv_pathway_simulation.ipynb` and `examples/basic_usage.py` written against older APIs — execute or update.
7. **Repo hygiene:** `.DS_Store`, `__pycache__/`, `.pytest_cache/`, `.ruff_cache/`, empty `.claude/worktrees/`; stale pre-ODE artifacts mixed into `outputs/`.
8. **Missing ISEF writing deliverables** — and per the published 2027 rules, the Research Plan must contain: Rationale, Research Question/Hypothesis/Expected Outcomes, Materials, Procedures, Risk and Safety, Data Analysis, Bibliography.

### External dependency (NOT a defect — do not fake it)
Thompson's patient series hasn't arrived. `compare_patient_series.py` stays data-gated; note in the plan that final clinical comparison happens when data lands.

---

## 4. INVARIANTS THAT MUST NOT BE BROKEN

From `ARCHITECTURE.md` post-audit invariants + Thompson review — enforce in every subagent brief:
1. `bkpyv_ode` stays canonical, registered in `ExperimentRunner.simulator_registry`; `bkpyv_specific` is back-compat only.
2. Time units explicit: ODE rates per **day**; configs declare `time_unit:`.
3. `simulator_parameters:`/`simulator_params:` must reach simulator constructors; YAML `1e-6` written as `1.0e-6`.
4. No unconditional success banners in `validate_*.py` — nonzero exit on failure.
5. The copies/mL bridge stays assumption-labelled on all clinical-facing outputs.
6. **Never make TAg drive `dCC/dt`** — the S-phase gate is host-driven only.
7. Terminology: *replication* = intracellular genome copying; *production* = system-level virion output; drug effects are on **production**.
8. Data-derived claims cite their dataset + context (in-vivo biopsy vs in-vitro RPTEC vs microarray cohort). No laundering phenomenological parameters into "measured" ones.

---

## 5. REQUIRED EXECUTION SHAPE

Produce `PLAN.md` (your one authored artifact) then dispatch SWE-2 Max subagents; You should only plan. Suggested decomposition — reorder/merge as needed, but every item must land:

- **T1 Environment & green tests.** Create `.venv`; `pip install -e ".[dev,ui]"`; clean-reinstall broken packages (matplotlib, python-dateutil, rich, pydantic); `python -m pytest tests/ -q` → **0 failures**; write `outputs/test_results.txt`; run `ruff check src/ tests/` + `black --check src/` — fix or document.
- **T2 Repo hygiene & commit.** Add `src/vcm/evaluation/__init__.py`; remove `.DS_Store`, `__pycache__`, `.pytest_cache`, `.ruff_cache`, empty `.claude/worktrees/`; move `GSE72925_RAW.tar` → `data/`; decide `data/mock/` (regenerate the 3 documented mock JSONs or fix README); clean commits in repo style. **Do NOT delete the tars, PDFs, or anything under `data/`/`outputs/` without flagging the user first.** Treat silence as "keep".
- **T3 Streamlit completion.** Implement or remove the placeholder Comparison page; smoke-test all pages headless.
- **T4 Docs consistency pass.** Align `judges_mentors_guide.md`, `technical_documentation.md`, `QUICK_START.md`, `README.md` with canon; **ISEF 2026→2027**; fix the "unrelated datasets" claim; Windows-correct commands; keep archived docs archived.
- **T5 Examples & notebook.** Run `examples/basic_usage.py`; `nbconvert --execute` the notebook; fix to current API.
- **T6 ISEF writing deliverables** (in `docs/`): `ISEF_ABSTRACT.md` (≤250 words, official abstract format), `ISEF_RESEARCH_PLAN.md` (all 8 sections per the 2027 rules: Rationale / RQ-Hypothesis-Expected Outcomes / Materials / Procedures / Risk & Safety / Data Analysis / Bibliography — note human-data rules: using only de-identified public GEO datasets + collaborator-provided de-identified series keeps it non-human-subjects, but flag for the student's SRC), `REFERENCES.md`, `PAPER_OUTLINE.md`, `PRESENTATION_POINTS.md`. All quoted numbers regenerated from current outputs.
- **T7 Real-data integration (the big one).** Sub-steps:
  a. Extract `data/GSE317012_RAW.tar` → `data/raw/GSE317012_RAW/`; run `scripts/reorganize_gse317012.py` to group 10x triplets per GSM.
  b. Fetch the GSE317012 series-matrix/sample metadata for **true condition labels** (control/peaking/resolving) — do NOT reuse the attic script's filename-substring guessing.
  c. Build a maintained `scripts/analyze_gse317012.py` (scanpy or a minimal scipy/`mmread` path — planner's call; if scanpy, add a `[data]` extra to pyproject rather than bloating core deps): QC → per-sample tubular-epithelial + immune fractions → pathway scores (S/G2M, DDR, mitochondrial, translation, antigen-presentation gene sets already listed in `docs/bkpyv_research_to_model_map.md`) per condition → fold-changes peaking vs control and resolving vs peaking → `data/processed/gse317012_pathway_scores.csv` + a short `docs/GSE317012_ANALYSIS.md` with methods and the JCI Insight 2026 citation.
  d. If in-vivo scores clearly quantify a coefficient the model carries (e.g., `translation_enhancement`, `mitochondrial_importance`, innate/antigen-presentation suppression), update `parameters.py`/ODE defaults **with provenance comments** ("in-vivo biopsy scRNA-seq, GSE317012/JCI Insight 2026, peaking-vs-control log2FC") — and add a test asserting the value. If results are too noisy to justify changes, document that honestly instead.
  e. Microarray tars (GSE47199/75693/72925): prefer GEO **series-matrix processed data** over raw-CEL R pipelines to avoid an R dependency; a summary-level analysis (PVAN/viremia vs stable signature tables) is sufficient — do not build a full R pipeline unless trivial. If effort balloons, deliver the inventory + processed-data summary and mark raw-CEL analysis out-of-scope with rationale.
- **T8 Regenerate outputs & final verification.** Re-run review bundle, all three `validate_*.py`, sensitivity, drug-switch, `compare_patient_series.py` (graceful no-data exit); confirm CLI `list-plugins`/`run`/`compare`/`inspect`; write `FINAL_STATUS.md` with the verification transcript summary.

Parallelize where safe (T4/T5/T6 independent; T1 gates T3/T7/T8; T2 sequencing is your call). Each subagent brief must include: repo root, the §4 invariants, the §1 honesty posture, exact acceptance commands, and "report files changed + command outputs".

---

## 6. DEFINITION OF DONE (acceptance criteria)

Finished only when ALL hold on this machine:

1. `python -m pytest tests/ -q` → 0 failures in `.venv`.
2. `python -m vcm list-plugins` shows 5 plugins; `run` + `compare` on `transplant.bk_polyomavirus` ODE configs succeed.
3. `validate_bkpyv.py`, `validate_isef.py`, `validate_clinical.py` exit 0; JSON reports regenerate passing.
4. `scripts/generate_bkpyv_review_bundle.py` regenerates `outputs/bkpyv_review/` (4 PNGs + 3 CSVs + manifest).
5. Streamlit launches; no placeholder pages.
6. `git status` clean; logical commit history.
7. No doc contradicts the model card on drug mechanism, time units, or validation posture; README tree matches reality; all docs say ISEF **2027**.
8. Notebook + example execute against current API.
9. The five `docs/ISEF_*`/`REFERENCES`/`PAPER_OUTLINE`/`PRESENTATION_POINTS` deliverables exist, quote only regenerated numbers, and the Research Plan has all 8 required sections.
10. GSE317012 is extracted + organized; `data/processed/gse317012_pathway_scores.csv` and `docs/GSE317012_ANALYSIS.md` exist with true condition labels and the JCI Insight 2026 citation; any parameter updates carry provenance comments + a test.
11. `FINAL_STATUS.md` exists.

**Explicit non-goals:** GitHub push to `origin` (`website-deployer/BKPyV-VCM`) — only if the user explicitly approves AND credentials exist; otherwise leave local. The Thompson patient-data comparison (awaiting data). New biological features beyond the data-integration path. Any weakening of caveat language. Claiming microarray cohorts were "single-cell" or that biopsy data is in-vitro.

**Escalate, don't guess:** if a subagent finds the ODE science itself wrong, a task needs destructive deletion/credentials, or the GEO metadata can't be obtained, surface it to the user instead of improvising.

**Literature anchors for REFERENCES.md** (correct attributions — the repo previously mis-cited some): Hirsch 2016 AJT 16(3):821-832 (PMID 26639422); Needham/Perritt/Thompson PLoS Pathog 2024;20(12):e1012663; Weissbach J Virol 2024;98(12):e01382-24 (NOT "Jang JY"); Funk J Infect Dis 2006;193:80-87 (PMID 16323135); Kotton Transplantation 2024;108(9):1834-1866 (PMID 38605438); Hirsch & Randhawa AST IDCOP Clin Transplant 2019;33:e13528; Myers et al. Viruses 2025;17(1):50 (PMC11768487 — 6-state patient-scale ODE calibrated on Duke/CTOT data; cite to position novelty as the cell-mechanism layer, not "first BKPyV model"); Demey J Clin Virol 2018;109:6-12 (PMID 30343190); Fang Front Immunol 2022;13:971531 (PMC9428263); Yamauchi Ren Fail 2025; JCI Insight 2026 biopsy scRNA-seq (doi:10.1172/jci.insight.198227, GSE317012); JCI Insight urinary transcriptome benchmarking (doi:10.1172/jci.insight.201060); JVI 2025 IGKC+ PT paper (doi:10.1128/jvi.01394-25); Sigdel Transplantation 2016 PVAN signature (PMID 27140517, GSE72925); Gosert 2008 rr-NCCR emergence (PMC2292223); NCCR propagation/Sp1-TFBS papers (PMC3249478; jvi.03625-14; jvi.01008-16); RPTE-hTERT archetype persistence model (PMC8546605); Lauterbach-Rivière J Med Virol 2025 TNF-α/NCCR (doi:10.1002/jmv.70210); Thompson email = "personal communication, S. Thompson lab, UAB, Aug–Sep 2026".

Begin now: verify state read-only → write `PLAN.md` → dispatch SWE-2 Max subagents → iterate to §6 → report final status.
