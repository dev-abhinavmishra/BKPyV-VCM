# GSE317012 Reanalysis — Per-Sample Pathway Scores and Condition Contrasts

Reanalysis of the single-cell RNA-seq kidney-biopsy series **GSE317012**
(in vivo renal biopsies from BKPyV-infected transplant recipients,
JCI Insight 2026, DOI `10.1172/jci.insight.198227`) plus an inventory of
the three supporting microarray series. All quantitative outputs below are
computed from the actual downloaded source files; nothing is a fixture.

## 1. Provenance

| Artifact | Source | Recorded |
|---|---|---|
| `data/GSE317012_RAW.tar` (156 members, 1,808.1 MB) — **deleted 2026-09-19** | local archive, since reclaimed | sha256 `1EA67E31A190AA05F39A7664721BC0CAF3D219997DFE67C2FC7920C4854327A8` (unchanged by pipeline — verified before and after; hash remains the re-download check) |

> **Storage note (2026-09-19):** the raw archive `data/GSE317012_RAW.tar`
> and the extracted tree `data/raw/GSE317012_RAW/` were deleted under
> user-authorized storage reclamation (public GEO data). The retained
> evidence is this document plus `data/processed/gse317012_*` outputs and
> `data/processed/gse317012_provenance.csv`. To reproduce: re-download
> `GSE317012_RAW.tar` from
> `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/suppl/GSE317012_RAW.tar`,
> verify the sha256 above, then run
> `python scripts/reorganize_gse317012.py data/raw/GSE317012_RAW --extract data/GSE317012_RAW.tar --move --verify`
> and `python scripts/analyze_gse317012.py`. |
| `data/research/gse317012_series_matrix.txt.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/matrix/GSE317012_series_matrix.txt.gz` | fetched 2026-09-20T02:39:30Z, sha256 `221b8892bfb3e4277e467ecc8aad13acf4f4a4348d2b3e119ffcc6db1c58b83e` |
| `data/research/GSE47199_series_matrix.txt.gz` | NCBI GEO FTP HTTPS | fetched 2026-09-20T02:43:37Z, sha256 `d91c6825f0bc84f9b94830b2bc541f691a7e9d7589c1111e447683d0c00ccbc0` |
| `data/research/GSE75693_series_matrix.txt.gz` | NCBI GEO FTP HTTPS | fetched 2026-09-20T02:43:39Z, sha256 `d24457487e321e0dc658d9ba3e2a218b49fb7696fee4590ae6eb40e3278a64f9` |
| `data/research/GSE72925_series_matrix.txt.gz` | NCBI GEO FTP HTTPS | fetched 2026-09-20T02:43:45Z, sha256 `eb6079101597aa637ec1b3e0d5d0b04cec250244f969adc6027f2f5d4eff0d08` |

Every input and output file SHA256 is recorded in
`data/processed/gse317012_provenance.csv`. Web content was treated as
untrusted data: each download records URL, UTC timestamp, and checksum.

## 2. Extraction and organization

`data/GSE317012_RAW.tar` was extracted to `data/raw/GSE317012_RAW/` with a
traversal-safe extractor (`scripts/reorganize_gse317012.py::extract_tar`)
that rejects `..`/absolute member paths, skips identical existing files,
and refuses to overwrite foreign content. The organizer then **moved** the
flat files into `GSM*` directories, so no duplicate flat copy is retained.

Each of the 26 GSM directories contains **six** artifacts: the filtered
10x triplet (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`) and the
raw triplet (`raw_matrix.mtx.gz`, `raw_barcodes.tsv.gz`,
`raw_features.tsv.gz`). A pre-existing organizer bug collapsed the `raw_`
prefix into the filtered kinds and raised `Duplicate … file` — the fixed
regex keeps `raw_` as part of the kind key (`tests/test_gse317012_data.py`
captures the regression).

## 3. Authoritative sample metadata

Phase labels were parsed **only** from the downloaded series matrix
(`!Sample_characteristics_ch1`), never from filenames or hardcoded lists.
All 26 archive GSMs cross-check exactly against the 26 matrix GSMs.

Parsed phase counts: **12 Control, 5 Peaking, 9 Resolving** — consistent
with the expected composition reported by the study. Phase strings were
keyword-derived from the actual characteristics text; any unparsed sample
would surface as `Unclassified` (none did).

Sample metadata written to `data/processed/gse317012_sample_metadata.csv`.

## 4. Methods

- **Input**: per-sample **filtered** 10x matrices only (raw triplets are
  preserved but not analyzed here).
- **Representation**: `scipy.sparse` CSR throughout; per-sample streaming;
  **no whole-matrix densification**.
- **QC** (`min_genes=200`, `max_mt_frac=0.20`, per-cell): kept
  **34,987 / 294,286** barcodes (11.9%). Exclusion accounting is fully
  transparent: 189,313 cells failed `genes>=200`; 259,100 exceeded the 20%
  mitochondrial fraction (overlapping sets). The dominant high-MT
  exclusions reflect genuinely stressed biopsy cells — e.g. GSM9463812 has
  median 125 genes/cell and median MT fraction 0.969, while the largest
  sample GSM9463837 has median 958 genes and MT 0.11. Retained cells are
  therefore the higher-quality subset; this selection is a documented
  caveat, not hidden data loss.
- **Normalization**: library size to 10k per retained cell, then `log1p`.
- **Marker-proxy fractions** (clearly labelled proxies, not deconvolution):
  fraction of retained cells expressing ≥1 marker — epithelial
  {LRP2, CUBN, UMOD, SLC34A1}, immune {PTPRC}.
- **Pathway scores** = mean normalized expression over each curated gene
  set, averaged across retained cells. Gene sets come from
  `docs/bkpyv_research_to_model_map.md` where enumerated (S/G2M:
  CLSPN/TOP2A/MKI67; DDR: BRCA1/BRCA2/PRKDC/FANCI/MMS22L; mitochondrial:
  MT-ND4/MT-CO1/MT-CYB/MT-ATP6) and are canonical sets where the map names
  a pathway without enumerating it (translation: RPL/RPS panel; antigen
  presentation: MHC-I/II + TAP1/TAP2 + B2M). **Missing-marker report:** all
  curated genes were present in all 26 samples (missing_markers empty).
- **Statistics are sample-level**: the 26 biopsies are the replicates.
  Contrasts use mean difference + Cliff's delta + two-sided Mann-Whitney U
  across samples — never fold-change on centered scores, and cells are
  never treated as independent patients.
- **Resource bounds**: full run ~264–331 s, peak Python allocation
  ~540 MB, 34,987 cells analyzed.

## 5. Results

Per-sample rows: `data/processed/gse317012_pathway_scores.csv`
(26 rows; columns include `cells_pre_qc`, `cells_post_qc`,
`cells_excluded`, `barcodes_in_file`, `frac_epithelial_proxy`,
`frac_immune_proxy`, `score_s_g2m`, `score_ddr`, `score_mito`,
`score_translation`, `score_antigen_presentation`, `missing_markers`).

Contrasts: `data/processed/gse317012_contrasts.csv`. Headline rows:

| Variable | Peaking−Control (mean diff) | MWU p | Resolving−Peaking | MWU p |
|---|---|---|---|---|
| frac_immune_proxy | +0.187 | 0.160 | −0.190 | 0.190 |
| score_antigen_presentation | +0.170 | 0.442 | −0.164 | 0.298 |
| score_translation | +0.224 | 0.234 | −0.221 | 0.240 |
| score_s_g2m | −0.001 | 0.130 | −0.002 | 0.438 |
| score_ddr | +0.001 | 0.574 | −0.007 | 0.699 |
| score_mito | −0.017 | 0.879 | +0.245 | 0.438 |
| frac_epithelial_proxy | +0.055 | 0.574 | +0.078 | 0.699 |

**Interpretation is deliberately limited**: with n=5 peaking samples, no
contrast reaches conventional significance, and none should be reported as
such. Directions are nonetheless coherent with the JVI mechanism map —
immune-proxy fraction and antigen-presentation/translation scores trend up
at the peaking stage and recede while resolving — reported as descriptive
trends only.

## 6. Microarray inventory (`data/processed/microarray_inventory.csv`)

| Series | Local raw tar | Series matrix | GSMs | Platform | Context | Expression values in matrix |
|---|---|---|---|---|---|---|
| GSE47199 | 57 `.CEL.gz` | downloaded | 57 | GPL6244 (HuGene-1_0-ST) | biopsy, blood, kidney | **No — metadata only** |
| GSE75693 | 79 `.CEL.gz` | downloaded | 79 | GPL570 | biopsy, kidney | **No — metadata only** |
| GSE72925 | none local | downloaded | 168 | GPL570 | biopsy, kidney | **No — metadata only** |

Per-GSM summaries: `data/processed/{gse}_series_summary.csv`.
These are **microarrays — never labelled single-cell**. GEO series matrices
for all three series carry sample metadata but **no processed expression
value table**; raw CEL-level analysis is explicitly out of scope for this
sprint (inventory + processed-metadata summary delivered instead).

## 7. Parameter decision: NO CHANGE

No biological coefficient or ODE default was modified. The scRNA contrasts
above are observational associations across biopsy phases; they do not
identify kinetic parameters (rates, half-lives, dose-response slopes) of
the BKPyV ODE. Mapping e.g. a +0.19 immune-fraction trend onto
`innate_immune_suppression_factor` would be an unjustified leap —
correlation across conditions is not a rate constant. The dataset instead
informs *directions* that the model already encodes qualitatively
(immune/antigen-presentation/translation programs elevated at peak,
receding on resolution), which is noted here without altering parameters.

## 8. Reproduce

```powershell
.\.venv\Scripts\python.exe scripts/reorganize_gse317012.py `
  data/raw/GSE317012_RAW --extract data/GSE317012_RAW.tar --move --verify
.\.venv\Scripts\python.exe scripts/analyze_gse317012.py
```

Determinism: sorted inputs, fixed QC constants, no RNG — repeated runs
produce byte-identical CSVs (covered by
`tests/test_gse317012_data.py::test_output_schema_and_determinism`).

## 9. Caveats

- In-vivo biopsy context: phases describe biopsy disease state, not a
  controlled perturbation; causal claims are out of scope.
- Marker fractions are proxies, not deconvolution.
- This analysis does not identify patients and does not use or infer any
  patient-level clinical endpoint; it is hypothesis-generation context for
  the mechanistic model, not calibration evidence.
