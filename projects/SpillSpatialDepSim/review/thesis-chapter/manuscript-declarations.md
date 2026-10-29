# Manuscript Abstract and Declarations — Preserved, Not in the Chapter

Preserved per `chapter2-integration-plan.md`'s Assumption 1: dropped from
the Chapter 2 body for now, kept here verbatim (LaTeX escaping cleaned up
only — no wording changed) so re-inserting any of it later is a paste, not
a rewrite. Source: `source/Revisions_V2c.tex`, lines 126–136 (abstract) and
672–693 (Declarations).

## Abstract

**Background**: Investigators designing cluster randomized trials desire
insights directing treatment assignment methodology for studies involving
spillover effects and spatial dependence.

**Methods & Simulation**: Treatment assignment strategies including simple
random sampling (SRS) and block stratified sampling (BSS) are defined and
spatial autoregressive modeling is applied with consideration for spillover
effects and spatial dependence for estimation of intervention effects. A
simulation study is carried out comparing SRS and BSS sampling methods on
spatial grids of varying sizes. A range of spillover effects and levels of
spatial dependence were considered for estimation of the intervention
effect via a spatial autoregressive (SAR) model.

**Results**: Findings of an extensive simulation study comparing simple
random sampling and block stratification methods indicate that randomly
selected treatment assignments result in best case reduced Mean Squared
Error (MSE) when estimating intervention effects, but block stratified
treatment assignments generally lead to less variation in MSE among a series
of treatment combinations, indicating that a block stratified treatment
arrangement will not achieve the minimal level of estimation error, but it
remains robust across a range of selected parameters. While SRS achieves
lower average MSE in certain scenarios, both designs exhibit substantial
estimation bias when treatment and spillover indicators are collinear.
Variation in performance across all possible treatment combinations and the
worst-case allocation risk this entails motivates consideration of BSS as a
more consistent alternative design.

**Conclusion**: The relationship between spillover effects and mean squared
errors (MSE) of intervention effect estimation is apparent. The MSE for the
intervention effect, which is the average MSE over each of $N$ simulation
iterations, is minimized for some combinations of random sampling treatment
assignment, but block stratified assignment generally reduces variance among
combinations of possible treatment arrangements. In short, SRS may achieve
lower average MSE in some scenarios, but BSS bounds worst-case allocation
risk through lower variation in estimation error across treatment
combinations; the preferred strategy depends on the anticipated spillover
mechanism, with BSS most advantageous when spillover is expected to flow
primarily into control clusters.

## Declarations

### Ethics Approval and Consent to Participate

Not applicable. This study involves simulation data only and does not
require ethics approval or consent.

### Consent for Publication

Not applicable.

### Availability of Data and Materials

All data supporting the findings of this study are generated via simulation
for the included simulation study. The simulation data and R code are
available upon request and with permission of the authors.

### Competing Interests

The authors declare that they have no competing interests.

### Funding

Dr. Van Deinse is supported by the National Institute of Mental Health
(K01MH129619). The trauma-informed supervision approach referenced in this
manuscript was funded by the Bureau of Justice Assistance
15PBJA-24-GG-01896-PRJH and awarded to the North Carolina Department of
Adult Correction. The funding bodies had no role in the design of the
study, data collection, analysis, or interpretation, nor in writing the
manuscript.

### Authors' Contributions

A.W. conceptualized and designed the simulation study, wrote the manuscript
text, and prepared Figures 1–7 and 10. T.V. provided Figures 8 and 9 and
contributed advisement for Section 4 of the manuscript. F.-C.L. provided
methodological guidance and critical review throughout. All authors
reviewed and approved the final manuscript.

### Acknowledgements

The authors thank the editor and reviewers for their thorough and
constructive feedback, which has substantially improved the manuscript.

## Where these could go later

Amber's dissertation (`templates/source/dissertation-resources-amber/`)
keeps two of these statements, not all seven: `Funding` and `Data
Availability Statement`, as unnumbered sections placed after the chapter's
Discussion (after §3.6, before Chapter 4), appearing in the table of
contents without numbers. She drops ethics, consent, competing interests,
and author contributions entirely.

Project 1's simulation code is public in the `SpatialCRT` repo, so a data
availability statement is genuinely applicable and will likely be wanted
for the thesis even if omitted from the prelim. If reinstated, the
manuscript's own "available upon request" wording above would need
updating to reflect the code's actual public location — that is a
substantive edit, not something to make silently if this file is ever
pasted back in.
