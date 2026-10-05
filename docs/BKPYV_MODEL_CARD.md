# BKPyV virtual cell model card

## What changed after immunology review

The model now makes four distinctions explicit:

1. **S phase is a gate.** Host cell-cycle entry and DNA synthesis must be present before large T antigen can accumulate in the ODE model. This encodes the primary renal tubular cell observation that T antigen is not assumed to be the initiator of cell-cycle re-entry.
2. **Replication is not production.** `intracellular_replication_flux` tracks genome-copying activity, while `viral_production_rate` and clinical viral load represent broader system-level output. The UI uses “production” when discussing drug effects.
3. **Tacrolimus acts through immune control.** It reduces immune-pathway control of infected-cell persistence; it is not a direct multiplier on the genome-copying term.
4. **NCCR is a hypothesis-testing axis.** `archetype` is the presumed persistent form. `rearranged` increases early-gene bias and lowers the capsid expression multiplier as a scenario motivated by the reviewer’s comment. These coefficients are not clinical estimates.

## Interpretation boundaries

The model is mechanistic and qualitative-to-semiquantitative. It is not a validated patient-specific predictor, treatment recommender, or substitute for plasma/urine BKPyV testing. Most coefficients remain phenomenological until calibrated against primary RPTEC experiments and longitudinal transplant cohorts.

The model also keeps primary-cell and in-vivo contexts separate: renal epithelial cell culture results are not automatically treated as kidney-level kinetics. The clinical mapper is an assumption-labelled visualization bridge, not a calibration claim.

**Engine conventions (since 2026-09):** the ODE engine (`bkpyv_ode`) is canonical; all rates are per **day**. Drug inputs are dimensionless dosing intensities (driven continuously while their perturbation window is active), not ng/mL pharmacokinetics. With default parameters the no-drug infection sustains a low plateau rather than self-clearing; sirolimus attenuates and tacrolimus elevates. The legacy discrete simulator remains for backward compatibility (hourly steps; drug effects are direct multipliers there).

## Suggested experiments

- Compare archetype and rearranged NCCR under identical host-cell permissiveness.
- Sweep the S-phase gate and measure time to T-antigen onset.
- Compare tacrolimus effects on immune-control index versus intracellular replication flux.
- Fit only the clearly marked parameters to urine/plasma trajectories, reserving the rest for sensitivity analysis.

## Domain-expert review notes (September 2026)

Reviewed by the group behind the S-phase single-cell finding (Needham/Perritt/Thompson). Confirmed and corrected:

1. **The gate is the right architecture.** "If there is no initial host S phase then the virus never expresses T antigen." Critically, T antigen is NOT the trigger for cell-cycle entry (it is not expressed in G0/G1 when entry would need to happen) — the driver of cell-cycle re-entry in terminally differentiated RPTECs is unknown. The ODE model matches this: CC/DNA dynamics depend on host state only, never on viral variables; do not let future edits make T antigen drive `dCC/dt`.
2. **Drug wording matters.** We describe drug effects on viral *production* (system-level virion output). Tacrolimus promotes production by weakening immune targeting of infected cells; sirolimus may prevent production by restricting S-phase entry/T-antigen expression. *Replication* is reserved for intracellular genome copying.
3. **Cell-context caution on the mTOR axis.** Drug-inhibition magnitudes depend on cell type: immortalized lines (e.g. HEK293) have dysregulated cell cycles and yield at least ~1 log lower BKPyV titers than primary RPTECs, despite constant S phase. Model parameters anchored to Hirsch 2016 (primary RPTECs) are direction-, not magnitude-, validated.
4. **Single-cell ≠ kidney.** The RPTEC single-cell findings are biological priors, not in-vivo kinetics; collaborators are testing whether they carry over to human kidneys.
5. **Archetype vs rearranged NCCR is the axis the field cares about.** The rearranged form (TAg overexpression, reduced capsid expression) is the clinically important one; the archetype is the transmitted/persistent form. The model's NCCR scenario axis exists for this question, and collaborators offered longitudinal patient urine/plasma data for comparison.

## Literature anchors

- The S-phase gate is motivated by the single-cell study showing that host S phase drives large T-antigen expression during BKPyV infection: [PLOS Pathogens](https://journals.plos.org/plospathogens/article?id=10.1371%2Fjournal.ppat.1012663).
- The NCCR axis is motivated by primary work showing that rearranged NCCR is a major determinant of replication in RPTE cells: [Efficient Propagation of Archetype BK and JC Polyomaviruses](https://pmc.ncbi.nlm.nih.gov/articles/PMC3249478/).
- The archetype/rearranged distinction and its disease association are documented in [Polyomavirus BK with rearranged NCCR emerge in vivo](https://pmc.ncbi.nlm.nih.gov/articles/PMC2292223/).
- Patient-scale BKPyV mathematical modeling exists, so this project’s novelty claim should focus on the cell-mechanism coupling rather than “the first BKPyV model”: [Modeling BK Virus Infection in Renal Transplant Recipients](https://pmc.ncbi.nlm.nih.gov/articles/PMC11768487/).
