# References

All sources cited across this repository's model, data pipeline, and ISEF
documentation. Attributions verified against the authoritative list in
`ASTRA_ORCHESTRATION_PROMPT.md` §6.

## Primary mechanistic / single-cell evidence

1. **Needham CD, Perritt E, Thompson S, et al.** Single-cell analysis of BK
   polyomavirus infection in primary renal proximal tubule epithelial
   cells. *PLoS Pathogens* 2024;20(12):e1012663. — Host S-phase gating of
   T-antigen accumulation; archetype vs rearranged NCCR biology.
2. **JCI Insight 2026.** BKPyV kidney-biopsy single-cell study underlying
   GEO series GSE317012. doi:10.1172/jci.insight.198227.
3. **JCI Insight** urinary-transcriptome benchmarking study.
   doi:10.1172/jci.insight.201060.
4. **Weissbach FH, et al.** BKPyV biology informing the model's
   qualitative mechanism priors. *Journal of Virology*
   2024;98(12):e01382-24. (Author: Weissbach — not "Jang JY".)
5. **JVI 2025** IGKC+ proximal-tubule paper. doi:10.1128/jvi.01394-25.
6. **Hirsch HH, et al.** BK polyomavirus replication in renal tubular
   cells: sirolimus early-phase inhibition (in-vitro IC90 ≈ 4 ng/mL) and
   the tacrolimus/FKBP-12 mechanism motif. *American Journal of
   Transplantation* 2016;16(3):821-832. PMID 26639422.
7. **Funk GA, et al.** BK polyomavirus dynamics under immunosuppression;
   persistence and clearance modeling. *Journal of Infectious Diseases*
   2006;193:80-87. PMID 16323135.

## Clinical / epidemiological anchors

8. **Kotton CN, et al.** Updated international consensus guidance on BKPyV
   management in kidney transplantation. *Transplantation*
   2024;108(9):1834-1866. PMID 38605438. — Screening/treatment copies/mL
   threshold context for the assumption-labelled bridge.
9. **Hirsch HH, Randhawa P; AST Infectious Diseases Community of
   Practice.** BK polyomavirus in solid organ transplantation —
   guidelines. *Clinical Transplantation* 2019;33:e13528. — Consensus
   viral-load thresholds used as bridge anchors.
10. **Demey B, et al.** *Journal of Clinical Virology* 2018;109:6-12.
    PMID 30343190.
11. **Myers et al.** Six-state patient-scale BKPyV ODE calibrated on
    Duke/CTOT data. *Viruses* 2025;17(1):50. PMC11768487. — Cited to
    position this project's novelty as the **cell-mechanism layer**, not
    "first BKPyV model".
12. **Fang Y, et al.** Dynamic risk prediction of BK polyomavirus
    reactivation after renal transplantation. *Frontiers in Immunology*
    2022;13:971531. PMC9428263. — Source of cohort-model odds ratios
    (tacrolimus ≈2.3, prior transplant ≈2.1, male ≈1.6, diabetes ≈1.3).
13. **Yamauchi K, et al.** Integer-based BKPyVAN risk score (age, sex,
    prior transplant; AUC 0.65–0.68). *Renal Failure* 2025.
14. **Sigdel TK, et al.** PVAN transcriptional signature (GSE72925).
    *Transplantation* 2016. PMID 27140517.
15. **Gosert R, et al.** Rearranged-NCCR emergence. 2008. PMC2292223.
16. **NCCR propagation / Sp1-TFBS papers.** PMC3249478;
    doi:10.1128/jvi.03625-14; doi:10.1128/jvi.01008-16.
17. **RPTE-hTERT archetype persistence model.** PMC8546605.
18. **Lauterbach-Rivière L, et al.** TNF-α/NCCR. *Journal of Medical
    Virology* 2025. doi:10.1002/jmv.70210.
19. **Thompson S.** Personal communication, S. Thompson lab, UAB,
    Aug–Sep 2026 (correspondence on patient-series context; no patient
    data received — comparison script exits gracefully awaiting data).

## Public datasets

20. **GSE317012** — BKPyV-infected kidney-biopsy scRNA-seq (reanalyzed
    here; see `docs/GSE317012_ANALYSIS.md`). NCBI GEO series matrix:
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317012/matrix/
21. **GSE47199, GSE75693, GSE72925** — BKPyV-related microarray series
    (GPL6244/GPL570; series-matrix metadata inventoried; processed value
    tables not present in the matrices). NCBI GEO FTP.

## Rules / process

22. **Regeneron ISEF 2027 International Rules for Pre-College Science
    Research.** Society for Science.
    https://www.societyforscience.org/isef/international-rules/
    (incl. rules-for-all-projects rule 11: maximum 250-word abstract;
    checklist/AI-disclosure requirements).
