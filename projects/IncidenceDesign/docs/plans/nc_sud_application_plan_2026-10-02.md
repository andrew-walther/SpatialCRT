# NC community-college SUD application study

Author: Andrew Walther
Date: 2026-10-02
Status: continuous SAR approach approved; remaining design decisions under interview;
implementation has not started.

## Purpose and place in the project

Assess how the treatment assignment designs estimate τ on the 58 North Carolina
community-college service-area clusters, using their observed 2018–2021 SUD
incidence. This is a simulation-based application of the design study to real
geography and baseline incidence. There are no observed intervention outcomes
from which to estimate an actual intervention effect.

**Author clarification (October 2):** use the observed NC SUD incidence surfaces;
do not generate synthetic incidence for this application. Each year's observed
rates are fixed inputs. Simulation concerns treatment allocations and outcomes
under known intervention/spillover parameters, enabling estimation-error assessment.
Resampling incidence or introducing a latent/noisy incidence surface is not part
of the primary application.

**SRS interpretation and placement agreed; application choices remain open.** The author wants Project 1's practical
allocation-consistency argument represented: BSS offered reasonably good accuracy
while limiting poor-allocation downside in its studied settings. Project 2 must
assess that criterion for its own designs, not equate low pooled mean MSE with
uniform allocation quality or assume the earlier BSS recommendation transfers.
The completed benchmark and selected-setting allocation-risk pilot inform this
discussion. On October 2 the author approved retaining that pilot as supporting
evidence and using the NC application as the next confirmation, without first
extending the grid comparison to all ten incidence surfaces. The NC comparison
must include SRS and carry this question forward.

**Writing decisions (October 2 follow-up):** the author approved this Chapter 2
linking sentence:

> Chapter 2 examined block-stratified allocation as a practical way to obtain
> reasonably good estimation while limiting poor-allocation risk in its studied
> settings; this chapter evaluates both criteria across a broader set of designs
> under heterogeneous incidence.

The body does not need a direct numerical comparison with Project 1. A brief
introduction/discussion reference may explain the extension and inclusion of SRS
and BSS (Checkerboard) alongside the other designs. An appendix comparison is
optional only if relevant and necessary; none is commissioned. The author
approved retaining Checkerboard in Project 2's main design comparison, with
detailed identification diagnostics in the appendix. Describe the NC graph-based
adaptation explicitly; it is not a literal rectangular checkerboard.
Chapter 3 follows Chapter 2 in the prelim/thesis, whereas CTJ must stand alone
and cite the verified Project 1 publication to establish the extension. This
citation is required, not merely optional. No manuscript text has yet been changed.

Agree this document before changing the simulation. Then implement, verify, run,
interpret results, and revise the written application sections and exhibits.
The application can proceed before the full CTJ rewrite; final recommendations
should consider it alongside the grid study and allocation-risk pilot.

**Author-confirmed endpoint (October 2):** complete the application and integrate
its results into full drafts of the thesis chapter/appendix and CTJ manuscript/SI.
Build and settle the thesis chapter first, then trim/reorganize it for CTJ. A plan,
pilot or proposal-only application section is not the final deliverable.

Read with `continuation_prompt_2026-10-02.md`, `TODO.md`,
`allocation_risk_findings_2026-09-27.md`, and `application/AGENTS.md` / README.

## What is already settled or available

- All 100 counties map to exactly one of 58 service areas. Cluster-level incidence
  and ages 18–64 person-years are available for each year and the pooled period.
- Corrected numerator: 23,523 deaths over 2018–2021. The earlier 21,147 deaths are
  for reconciliation, not the application numerator; preserve the settled source
  choice. No additional source adjudication or contact is needed for this study.
- Cached queen and rook W use legal boundaries, row standardization and named
  cluster IDs. Queen is primary; the existing boundary decision includes links
  across water. Reuse these weights rather than rebuilding different geography.
- Read-only inspection on October 2 found 58 finite cluster rates for each of the
  four years and the pooled period; cached W is 58×58. Inputs must still be checked
  by cluster name when implementation begins.
- Eight existing irregular-map design adaptations plus complete-randomization
  SRS (29 treated) exist in `application/code/application_designs.R`.
- Existing application outputs are from synthetic incidence and older methods;
  they are not real-SUD design-performance results. The old runner suppresses
  fitting warnings, uses only both-arms spillover, and lacks revised seed/manifest
  safeguards. It should not simply be rerun and treated as current.

## Author decisions to settle

1. **Outcome model and units — revised by author clarification October 2.**
   The primary application simulates a continuous education outcome: observed
   SUD incidence informs allocation, but does not enter the outcome baseline
   (β = 0). Do not assume an incidence–baseline-education relationship. The known
   structural treatment effect τ is constant across clusters within a scenario,
   independent of their SUD rates; spatial propagation and spillover remain.
   The author also accepted retaining an incidence-related baseline (β = 1) as
   an explicitly hypothetical sensitivity for continuity with the grid model.
   This is not a measured incidence–education relationship. Neither analysis
   estimates an observed education effect or establishes a reduction in SUD.
   Do not label simulation-scale τ as deaths prevented per 100,000.
   **Pending:** sensitivity X scale and extent, primary fitted covariates, and
   regional incidence summary for saturation. No primary outcome X scaling is
   required. This supersedes the earlier matched-primary β = 1 proposal.
2. **Which incidence informs allocation?** Recommended primary: use each year's
   observed surface for that year's allocations, evaluating four fixed planning
   settings. Only the β = 1 sensitivity uses incidence in simulated outcomes.
   An earlier-year/later-year incidence sensitivity must have a distinct purpose
   under this revised primary model; the former proposal to generate primary
   outcomes from the later year's incidence no longer applies automatically.
3. **Adaptations and treatment budget.** Recommended: retain existing adaptations
   initially and report their treatment counts and treated population shares.
   SRS and incidence quartiles enforce 29 treated; graph coloring, rounded spatial
   blocks and region saturation may differ, and Isolation Buffer generally treats
   fewer. A universal 29-cluster budget would require new rules for some designs
   and is a separate scientific choice, not an automatic correction.
4. **Manuscript decisions settled.** Keep SRS framing conditional on spillover
   regime, retain the grid pilot as supporting evidence and use the NC application
   as the next confirmation. Retain Checkerboard in the main Project 2 comparison,
   with detailed diagnostics in the appendix. Use the approved Chapter 2 link and
   a verified Project 1 citation in standalone CTJ prose. No direct cross-project
   numerical comparison is required.

## Recommended primary analysis; remaining details require agreement

- Four fixed yearly incidence surfaces (2018–2021), reported separately. Pooled
  incidence is descriptive context; do not treat five overlapping surfaces as
  independent replicates. An average across years must be labeled as such.
- Continuous cluster-level SAR DGP:

      Y = (I − ρW)^−1 [τZ + S(Z) + ε]
      S(Z) = γWZ                 (both arms)
      S(Z) = γ(1 − Z)WZ         (control-only)
      ε ~ N(0, I); τ = 1 is proposed for the primary comparison.

  The accepted hypothetical sensitivity adds βX with β = 1. X is fixed observed
  incidence on a scale still to be agreed. The grid study already used βX, both
  before and after the September revision; τ was constant within a scenario in
  that study too. βX changes baseline outcomes, not the treatment-effect coefficient.

- **Allocation inputs:** retain observed rates. High Incidence Focus, Balanced
  Quartiles and Balanced Halves can compute ranks directly; converting rates to
  numbers in [0,1] adds no allocation requirement. Keep existing tie handling.
  Incidence-guided saturation also needs a regional summary: the existing code
  averages normalized cluster ranks. Ranking regions by mean observed rate can
  give a different order. Agree that rule explicitly, including whether regional
  rates are cluster averages or deaths divided by person-years.
- **Sensitivity outcome-model input X:** separately agree its scale relative to β and σ.
  One proposal is `(average rank(rate) − 0.5)/58`, matching the grid study's
  Poisson mode; the older application uses `rank(rate)/58`. This represents
  incidence ordering rather than absolute rate gaps. Raw or commonly rescaled
  rates are also possible, with an agreed β/σ calibration. Rank-based allocation
  does not itself determine this outcome-model choice. Do not silently change
  regional allocation rules when changing X's numerical representation.
- Fit existing validated SAR ML with true spillover covariate as primary. Fit the
  same engine omitting spillover as a clearly labeled sensitivity if included in
  the agreed scope. τ is the structural treatment coefficient; spatially propagated
  total impact is a different quantity and is not the primary estimand here.
  Proposal to settle in the implementation plan: omit incidence X from the
  primary fitted model, which has no incidence-related baseline; include X in
  the matched β = 1 sensitivity. Incidence adjustment in the primary analysis
  would be a separate choice, not implied by using incidence for allocation.
- Proposed starting grid: queen, ρ = {0, 0.2, 0.5}, γ = {0.5, 0.8}, both spillover
  regimes, all nine designs: 432 design/year/parameter blocks. Rook is a planned
  sensitivity, with its extent agreed after the queen pilot; a full τ sweep is
  optional and not required to establish the application.
- Proposed replication: smoke verifies mechanics; pilot measures speed and
  Monte Carlo precision; starting production request is 100 assignment draws and
  100 outcomes per unique allocation. Refine selected allocation-risk comparisons
  to 400 outcomes if necessary. Report achieved precision before choosing further
  replication, rather than assuming these counts resolve tail uncertainty.

## Design adaptations and comparison fairness

- Use one fixed four-region partition across years, selected by the existing
  compactness/population/cluster/county-balance criterion. Proposed reference
  population is the pooled period's person-years divided by four. Freeze and save
  the partition before evaluating performance; never select it by achieved MSE.
- Rank region means using the agreed incidence signal separately each year for
  incidence-guided saturation; use the same partition for plain saturation.
- Spatial blocks are an adaptation of grid 2×2 blocking, not literal 2×2 cells.
  Record their formation rules, incidence inputs and random seeds. Agree whether
  they remain year-specific or are frozen before running the comparison.
- Graph coloring on the irregular contiguity graph is an adapted interspersed
  assignment, not a rectangular checkerboard with guaranteed WZ = 1 − Z.
  Record its actual exposure structure, rank and treatment count. Preserve the
  implemented support unless the author approves a new randomization rule.
- Report treatment fraction, treated population share, incidence balance,
  treatment/spillover overlap and identification diagnostics beside MSE. Population
  imbalance does not silently turn τ into a population-weighted estimand.

## Performance and uncertainty

Primary: mean squared τ error, bias, empirical SD, 95% coverage, valid-fit counts,
warnings/aliasing and Monte Carlo SEs; compare each design with matched-setting
SRS using ratios/differences. Do not hide partial failures in `na.rm` summaries.

Secondary: allocation-specific MSE distribution, corrected within-setting SD,
q90 and worst-decile mean, retaining assignment frequencies and finite-outcome
uncertainty. Reuse the verified covariance correction for cached duplicates and
split-half diagnostics. Sampled maximum is supplementary and sample-dependent.

The years are fixed application settings, not four independent draws from a
population of incidence surfaces. Monte Carlo uncertainty concerns simulated
allocations/outcomes; it does not quantify uncertainty in observed rates. The
primary continuous model does not automatically give smaller populations larger
error variance. A population-dependent count/rate model or uncertainty in planning
incidence is a separately agreed sensitivity.

## Implementation and verification sequence

1. **Freeze approved study decisions and verify inputs.** Check cluster identities,
   year coverage, rates/denominators, and W ordering/row sums. Save parameters,
   input/source hashes, partition/block definitions and software versions.
2. **Implement a focused revised runner under `application/code/`.** Reuse data,
   weights, existing design helpers and the validated ML engine. Key-seed every
   assignment/noise draw; keep X fixed within a year/scenario; preserve outcome
   prefixes and assignment frequencies. Use new outputs under
   `application/results/real_sud_rev_20261002/`; do not resume old chunks.
3. **Verify before production.** Test actual treatment rules, named joins, DGP
   construction in both regimes, reproducible order/worker behavior, uncertainty
   arithmetic, warning/failure handling, and lean-versus-lagsarlm agreement on a
   small selection of the irregular-map fits. Smoke then pilot, with explicit
   completeness and Monte Carlo precision checks, before a full run.
4. **Interpret and assemble exhibits.** Report yearly results by regime, SRS
   comparisons, coverage/bias and design constraints. Write a source-linked results
   extract so every manuscript number can be traced to a current output.
5. **Revise writing after results.** Chapter body: methods, cluster map,
   substantive yearly/regime comparisons and conclusions. Keep relevant findings
   in the body, following recorded advisor guidance; use the appendix for exhaustive
   breakdowns, mechanical adaptation details and precision diagnostics. Resolve
   remaining author-confirmation comments and review figure/table placement.
   Derive CTJ
   application text from the agreed chapter; detailed material goes to CTJ SI.
   Render and review in fresh passes. Commit logical steps, update README/AGENTS
   and ROADMAP/TODO; push only with explicit permission for these new commits.

The manuscript-unification plan records a November 16, 2026 committee deadline
and approximately November 1 chapter/prelim target; confirm these with the author.
Chapter/prelim completion has priority, with CTJ derivation scheduled afterward.
The application plan does not override that schedule.

## Project 1 connection and publication boundaries

Use Project 1's allocation-consistency question to motivate evaluating both mean
and allocation risk. The new study can improve mean and upper-tail performance in
some settings without making its predecessor's preferred allocation universally
optimal. Preserve the accepted paper's scope; the historical exhibit discrepancy
remains unresolved and is not a prerequisite for this application.

Dissertation references may say Chapter 2. CTJ must cite the verified BMC paper
and stand alone. Project 1 remains read-only. Chapter 3 is edited only in
IncidenceDesign; current rules authorize its established sync hook into the prelim.
Restricted county-level source and derived files stay ignored. Only authorized
cluster-level results and aggregate statistics are candidates for tracking.
