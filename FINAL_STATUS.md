# Final Status — Virtual Cell Model (BKPyV / ISEF 2027)

Date: 2026-09-19. Head: `b1804f0`. Verified under `.venv` (Python 3.13.14,
editable install). Acceptance criteria from `ASTRA_ORCHESTRATION_PROMPT.md`
§6, each independently re-verified by the evaluator (sprint-5 sweep:
C1–C5, C7–C9 PASS on independent rerun; C10a/b + C6 + C11 completed this
sprint).

## Criterion-by-criterion

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| 1 | `python -m pytest tests/ -q` → 0 failures in `.venv` | **PASS** — 184 passed, exit 0 (85.74s) | `outputs/test_results.txt`; artifacts `sprint-04-regen-evidence.txt` |
| 2 | `vcm list-plugins` = 5; `run` + `compare` on `transplant.bk_polyomavirus` ODE configs succeed | **PASS** — all exit 0; compare uses `--scenario1/--scenario2` | `sprint-04-regen-evidence.txt` |
| 3 | `validate_bkpyv.py`, `validate_isef.py`, `validate_clinical.py` exit 0; JSON reports regenerate passing | **PASS** — plus `scripts/validate_ode_implementation.py`; all exit 0 | `harness_check_report.json`, `isef_validation_report.json`, `clinical_validation_report.json` (fresh 2026-09-19) |
| 4 | `generate_bkpyv_review_bundle.py` regenerates `outputs/bkpyv_review/` (figures + tables + manifest) | **PASS** — 9 files fresh 22:17; manifest↔dir bidirectional match | `outputs/bkpyv_review/manifest.json` |
| 5 | Streamlit launches; no placeholder pages | **PASS** — 7 pages, AppTest coverage + two live-browser passes (sprint-2); comparison page functional | `tests/test_streamlit_app.py`; artifacts `sprint-02-ui-evidence/` |
| 6 | `git status` clean; logical commit history | **PASS** — porcelain empty; 7 scoped commits `fc8c674..b1804f0` | this file's commit list below |
| 7 | No doc contradicts the model card; README tree matches reality; docs say ISEF **2027** | **PASS** — mechanism wording model-card-aligned; tree corrected; zero "ISEF 2026" strings | `docs/judges_mentors_guide.md`, `README.md`, `docs/bkpyv_research_to_model_map.md` |
| 8 | Notebook + example execute on current API | **PASS** — `bkpyv_pathway_simulation.ipynb` executed inplace (sprint-2); `examples/basic_usage.py` exit 0 | `notebooks/*.ipynb`; `outputs/basic_example/` |
| 9 | Five ISEF deliverables exist; quote only regenerated numbers; Research Plan has 8 required sections | **PASS** — plus AI-disclosure placeholder + SRC-OPEN flag; abstract 190 words (≤250, official 2027 rule 11) | `docs/ISEF_ABSTRACT.md`, `ISEF_RESEARCH_PLAN.md`, `REFERENCES.md`, `PAPER_OUTLINE.md`, `PRESENTATION_POINTS.md` |
| 10 | GSE317012 extracted + organized; pathway-scores CSV + analysis doc with true labels + JCI Insight 2026 citation; parameter updates carry provenance + test | **PASS (completed, then raw reclaimed)** — 26 GSMs × 6 artifacts organized, 12C/5P/9R parsed, 34,987/294,286 cells post-QC; **no coefficient updated** (documented freeze rationale) | `data/processed/gse317012_*.csv`, `docs/GSE317012_ANALYSIS.md`, `tests/test_gse317012_data.py` |
| 11 | `FINAL_STATUS.md` exists | **PASS** — this file | `FINAL_STATUS.md` |

## Commit history (scoped, logical — no force, no push)

| Hash | Scope |
|------|-------|
| `fc8c674` | chore(env): packaging + test infra (sprint-1) |
| `4e1a5d2` | feat(model): ODE engine + clinical bridge + configs + validators + core tests |
| `688bdbc` | feat(ui): streamlit + evaluation + executed notebook (sprint-2) |
| `7efdea9` | feat(data): GSE317012 integration (sprint-3) |
| `9b083c0` | chore(scripts): analysis/regen scripts + attic + regenerated outputs |
| `c65938a` | docs: ISEF deliverables + docs consistency + archival layout (sprint-4) |
| `b1804f0` | docs: root docs → docs/archive move completion |

## Honest shortfalls / open items

1. **Lint debt — documented, NOT passing**: 1,628 ruff findings and 33
   black-format diffs predate this work; measured in
   `outputs/lint_format_baseline.txt`, carried as accepted debt. New and
   touched `.py` are lint-clean; no blanket reformat was applied.
2. **Raw data reclaimed (planner-authorized, user-approved)**: the
   `GSE317012_RAW.tar` archive, three microarray tars, and the extracted
   tree were deleted 2026-09-19 for disk space. Public GEO data —
   re-download path + recorded sha256 are in
   `docs/GSE317012_ANALYSIS.md` §1 and printed by both scripts'
   graceful-failure messages. Retained evidence: `data/processed/`,
   `data/research/` series matrices, `gse317012_provenance.csv`.
3. **Stale outputs — RESOLVED (sprint-6)**: all `outputs/clinical/` files
   regenerated fresh via their owning scripts
   (`generate_risk_dataset.py` → `save_risk_results.py` →
   `generate_risk_figures.py` → `generate_clinical_figure.py`, all exit 0)
   plus all tracked run artifacts refreshed via `vcm run`. Remaining
   intentionally-old files, all disclosed: `outputs/calibration/*` (owner
   scripts archived in `scripts/attic/`; required input
   `data/processed/parameter_calibration.json` no longer exists),
   `.DS_Store` ×4 (junk, preserved per no-delete rule),
   `outputs/figures_summary.md` (hand-authored, self-declared archived),
   `outputs/lint_format_baseline.txt` (intentionally historical record),
   `outputs/ode_validation/*_Enhancement/Inhibition_*.png` (fresh content;
   filenames retain legacy mechanism wording — renaming = delete+add,
   documented instead).
   **`outputs/clinical/` numbers note**: regenerated risk-model results
   remain on the synthetic cohort (labelled synthetic; see
   `data/processed/cohort_generation_params.json`).
4. **SRC eligibility — OPEN**: whether secondary analysis of public
   biopsy-derived datasets triggers SRC human-participants review is
   flagged OPEN for the student + adult sponsor; see
   `docs/ISEF_RESEARCH_PLAN.md` §Risk and Safety.
5. **AI-usage disclosure placeholder**: `ISEF_RESEARCH_PLAN.md` carries an
   honest placeholder; the genuine usage record must be assembled by the
   student — no log was fabricated.
6. **Left uncommitted by design**: `.venv`, `data/`,
   `outputs/bkpyv_review/` remain gitignored by intent. Browser-evidence
   PNGs relocated out of the repo to the epic artifacts dir (not
   committed, bytes preserved). (`src/vcm/data/` is tracked — the
   over-broad `data/` gitignore rule was root-anchored to `/data/` in
   sprint-6 so future `data`-named source dirs stay trackable.)
7. **Config comments — RESOLVED (sprint-6)**: stale mechanism wording in
   `configs/*.yaml` comments corrected to model-card-aligned phrasing
   (comment text only; `yaml.safe_load` proves every config parses
   identically pre/post). Two YAML *value* fields retain legacy wording
   (`bkpyv_sirolimus.yaml`/`bkpyv_tacrolimus.yaml` `research_grounding` /
   `description` strings) — left untouched because they are parsed values,
   not comments; flagged here.
8. **Science posture unchanged**: model = mechanistic hypothesis
   generator; copies/mL outputs are assumption-labelled bridges;
   parameters phenomenological/literature-anchored — no coefficient was
   updated from GSE317012.
