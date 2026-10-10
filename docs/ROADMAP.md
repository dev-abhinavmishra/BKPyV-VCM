# ROADMAP

## Next five improvements (ranked by scientific value)

1. **Patient-level viral-load fitting.** The single biggest gap is the
   absence of serial qPCR data. Options: request the Funk 2008 dataset
   (518 day-matched urine/plasma pairs, n=223) from the authors; the BENEFIT/
   directed-therapy trial data sometimes surface via data-use agreements;
   or a transplant-centre mentor's de-identified qPCR pulls. With patient
   curves, replace range-containment with real fitting + posterior
   uncertainty (mixed-effects or Bayesian hierarchical).
2. **Cell-level infection readout.** Add a reference containing BKPyV genes
   (or use a dataset that does — e.g., scRNA of BKPyVAN biopsies with viral
   alignment, or the Weissbach 2024 single-cell data) so infected calls are
   measured, not signature-inferred. Then calibrate `I`, `T`, `DNA` states
   against observed infected-cell fractions.
3. **BKPyV-specific adaptive immunity.** Add a T-cell compartment
   (activation → expansion → effector → memory) so the model can speak to
   T-cell monitoring and adoptive T-cell therapy — the clinically live
   intervention space (e.g., Kotton 2024 guidelines discuss VSTs).
4. **Urine compartment + biopsy linkage.** Funk 2008's central result is
   tubular→urothelial cross-feeding with urine loads ~3000× plasma. A urine
   state variable would let the model test screening-by-urine and connect
   cellular dynamics to the most available clinical measurement.
5. **PK/PD dose-response.** Replace binary drug contexts with concentration-
   response (tac trough → immune-control reduction; sir level → mTOR
   inhibition), enabling dose-optimization questions ("what tac reduction
   clears viremia by week N?").

## Questions for a transplant ID / nephrology mentor
(ranked by how much they would change the model)

1. **What does "viremia clearance" mean operationally?** One negative qPCR,
   two consecutive below threshold, or below 1,000 cp/mL sustained? — decides
   the benchmark's clearance criterion (currently: below 1,000 cp/mL).
2. **Realistic viremia at intervention?** I assumed 10⁵–10⁷ cp/mL as the
   range where IS reduction is attempted. Is that where clinicians actually
   act, or earlier (10³–10⁴)?
3. **How is IS reduction dosed?** Step change, taper, or drug switch — and
   over what timescale? Determines whether my instant-stop intervention is
   a reasonable idealization.
4. **Is there de-identified serial qPCR data I can use?** Even 10–20
   patients' longitudinal loads would convert the whole Phase-2 exercise
   from benchmarking to fitting.
5. **Clinical urine:plasma sampling patterns?** Are they paired in practice?
   Would a urine-first model answer a real screening question?
6. **How much do real clearances vary with biopsy-proven PyVAN vs
   DNAemia-only?** The model treats clearance as one compartment; disease
   stage may split it.
7. **What intervention efficacies are realistic?** Funk 2006's median 22%
   — does that match what leflunomide/cidofovir/IS-reduction achieve in
   current practice, or is it dated?
8. **Would sirolimus conversion be modelled as sir-added or tac-removed?**
   The model supports either; clinically which is more relevant?
9. **Cell-mediated vs antibody rejection confounding?** Peaking biopsies
   often overlap with rejection — does the infected-signature risk being a
   rejection signature instead?
10. **BKPyV genotype / NCCR rearrangement** — clinically meaningful enough
    to keep the plugin's NCCR variant layer, or is that detail without
    usable parameter data?
