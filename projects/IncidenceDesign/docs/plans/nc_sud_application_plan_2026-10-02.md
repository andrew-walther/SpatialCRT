# NC community-college SUD application study

**Verified completion checkpoint (2026-10-02):** the observed-incidence application
and focused tail refinement are complete. Production: 1,248 reporting rows,
936 distinct sources, 7,598,000 independent outcome fits; separate tail refinement
adds 2,280,000 outcomes on the same allocations. All settings pass completeness
and the approved mean-MSE/coverage Monte Carlo gates; no warnings, aliases,
boundary fits or failures. Coverage remains below 95% and allocation tails remain
uncertain. See [completed findings](nc_sud_application_findings_2026-10-02.md).
Chapter/appendix integration and derived CTJ/SI are the active remaining work.

Earlier execution/pending statements below are historical and superseded by this checkpoint.

Author: Andrew Walther
Date: 2026-10-02
Status: study scope and concrete implementation plan approved October 2;
independent alignment review complete; implementation/verification underway.
The HTML visual companion is an additional approved deliverable. Production
findings and full thesis/CTJ drafts remain unfinished.

Read the [independent review](nc_sud_independent_review_2026-10-02.md) and
[concrete implementation plan](nc_sud_implementation_plan_2026-10-02.md).

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

**SRS interpretation and placement agreed.** The author wants Project 1's practical
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
   The accepted sensitivity uses `(average rank(rate) − 0.5)/58` with β = 1,
   on the same queen settings; its fitted model includes X. The primary fitted
   model includes treatment and true spillover without X. No primary outcome X scaling is
   required. This supersedes the earlier matched-primary β = 1 proposal.
2. **Which incidence informs allocation?** Recommended primary: use each year's
   observed surface for that year's allocations, evaluating four fixed planning
   settings. Only the β = 1 sensitivity uses incidence in simulated outcomes.
   An earlier-year/later-year incidence sensitivity must have a distinct purpose
   under this revised primary model; the former proposal to generate primary
   outcomes from the later year's incidence no longer applies automatically.
3. **Adaptations and treatment budget — approved October 2.** Retain existing
   adaptations and report their treatment counts and treated population shares,
   with the Balanced Halves correction below.
   SRS and incidence quartiles enforce 29 treated; graph coloring, rounded spatial
   blocks and region saturation may differ, and Isolation Buffer generally treats
   fewer. A universal 29-cluster budget would require new rules for some designs
   and is a separate scientific choice, not an automatic correction.
   Read-only count check (seeded 100-draw preview, 2018 rates): queen graph
   assignment treats 27; Isolation Buffer treats 13–19 (mean 16.02); High
   Incidence Focus, Balanced Quartiles and SRS treat exactly 29. Current Balanced
   Halves currently treats 28 because R rounds 29/2 = 14.5 to 14 in each half.
   **Author-approved correction:** treat 14 and 15 in the two halves, randomly
   choosing which gets 15, to retain exactly half of 58 overall. Preview counts
   are not production results; this correction has not yet been implemented.
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
      ε ~ N(0, I); τ = 1 is accepted for the primary comparison.

  The accepted hypothetical sensitivity adds βX with β = 1. X is fixed observed
  incidence on a scale still to be agreed. The grid study already used βX, both
  before and after the September revision; τ was constant within a scenario in
  that study too. βX changes baseline outcomes, not the treatment-effect coefficient.

- **Allocation inputs:** retain observed rates. High Incidence Focus, Balanced
  Quartiles and Balanced Halves can compute ranks directly; converting rates to
  numbers in [0,1] adds no allocation requirement. Keep existing tie handling.
  **Regional summary approved October 2:** calculate each year's cluster rates
  as deaths/population aged 18–64 × 100,000; assign average ranks to tied rates
  across all 58 clusters; average these ranks with equal cluster weights within
  each fixed region. Rank regions by that mean and assign 80%, 60%, 40%, 20%
  saturation from highest to lowest. Break regional-mean ties randomly per draw,
  retaining the existing rule. Affine scaling of ranks to approximately 0–1
  gives the same regional ordering. This follows the grid Poisson implementation,
  which averages its rank-normalized per-capita incidence covariate. Explain this
  procedure explicitly in the application methods; it is not a pooled regional
  death rate and does not average raw counts.

  **Focused sensitivity accepted with study scope; implementation pending:** compare this
  primary regional rule with (a) mean cluster rate and (b) summed regional
  deaths/summed regional population × 100,000. Freeze the same partition and all
  other study choices. First report yearly summary values, regional orderings and
  assigned saturation levels. If an alternative has identical ordering and tie
  structure, verify identical allocation support and reuse the corresponding
  performance results with explicit provenance. Simulate changed rules for
  incidence-guided saturation only, using matched random streams and reporting
  differences in MSE, coverage, allocation risk and treatment/population shares.
  Do not choose the primary rule after inspecting performance. At the proposed
  queen grid, two fully distinct alternatives would add at most 96 blocks
  (2 summaries × 4 years × 3 rho × 2 gamma × 2 regimes); scope remains subject
  to the concrete implementation-plan approval. No full factorial crossing with
  every other sensitivity is implied.
- **Accepted sensitivity X:** `(average rank(rate) − 0.5)/58`, matching the grid
  Poisson mode, with β = 1 and residual SD = 1. This preserves incidence ordering,
  not absolute rate gaps; it is hypothetical rather than an estimated education
  relationship. It does not change the approved regional allocation rule.
- Fit existing validated SAR ML with true spillover covariate. τ is the structural treatment coefficient; spatially propagated
  total impact is a different quantity and is not the primary estimand here.
  Accepted fits: primary `Y ~ Z + Spill`; hypothetical β = 1 sensitivity
  `Y ~ Z + Spill + X`, both with SAR dependence. The reviewer will assess how to
  explain differences as sensitivity to the matched baseline/adjustment
  specification, not a pure β effect. The reviewer found no need for an additional
  β = 0/X-adjusted bridge sweep.
- Accepted primary grid: queen, τ = 1, residual SD = 1, ρ = {0, 0.2, 0.5},
  γ = {0.5, 0.8}, both spillover regimes, all nine designs: 432
  design/year/parameter blocks. Rook sensitivity uses ρ = {0, 0.5}, γ = {0.5, 0.8}
  on all four years/nine designs, both regimes (288 blocks). The β = 1 sensitivity
  uses the same queen settings. Regional-summary alternatives apply to Design 8
  under the primary education model. Geographic/allocation sensitivities change
  one feature at a time; the baseline sensitivity deliberately changes the
  baseline and corresponding nuisance adjustment together.
- Accepted starting replication/precision targets: 100 assignment draws and 100
  outcomes per unique allocation; increase where needed to target mean-MSE
  Monte Carlo SE ≤5% of estimated MSE and coverage Monte Carlo SE ≤0.01.
  Fixed/near-fixed designs need additional outcomes, not redundant allocation
  draws treated as independent evidence. Use split-half tail diagnostics and
  targeted 400-outcome refinement, retaining remaining tail uncertainty.
  Independent review requires proven singleton designs to start production at
  1,000 outcomes, then use empirical precision gates. High Incidence Focus has
  one eligible allocation in every observed year (58 distinct rates), even though
  its generic helper metadata is not deterministic. For stochastic supports,
  use joint allocation/outcome SEs with duplicate-cache covariance for MSE and
  coverage; increase draws or outcomes according to the dominant uncertainty.
  Prespecified tiers and gates are in the implementation plan. Smoke and pilot
  precede production. No replication rule depends on a preferred ranking.

## Design adaptations and comparison fairness

- **Approved October 2:** use one fixed four-region partition across years,
  selected by the existing compactness/population/cluster/county-balance criterion.
  Reference population is the pooled period's person-years divided by four. Freeze and save
  the partition before evaluating performance; never select it by achieved MSE.
- Rank region means using the agreed incidence signal separately each year for
  incidence-guided saturation; use the same partition for plain saturation.
- Spatial blocks are an adaptation of grid 2×2 blocking, not literal 2×2 cells.
  **Author-approved October 2:** retain the existing location/incidence-rank/
  population grouping method. Form blocks once per year and freeze them across
  that year's allocation/outcome simulations and scenarios. Record block
  membership, input ranks, populations and reproducible seeds. Their year-specific
  formation differs from the four-region map, which stays fixed across years.
- Graph coloring on the irregular contiguity graph is an adapted interspersed
  assignment, not a rectangular checkerboard with guaranteed WZ = 1 − Z.
  Record its actual exposure structure, rank and treatment count. Preserve the
  implemented support unless the author approves a new randomization rule.
- Report treatment fraction, treated population share, incidence balance,
  treatment/spillover overlap and identification diagnostics beside MSE. Population
  imbalance does not silently turn τ into a population-weighted estimand.
- Verify cross-year primary equivalence for Designs 1/3/4/9 with fixed W/regions
  and β = 0. Reuse equivalent performance with explicit source IDs while computing
  annual population/incidence diagnostics separately; shared rows are not
  independent yearly confirmations. Do not assume equivalent β = 1 outcomes.

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

**Schedule superseded by author October 2:** disregard the earlier deadline dates.
Complete this project as soon as possible, finishing the application and chapter/
appendix promptly, then derive CTJ and its supplement for submission soon afterward.
CTJ remains a required endpoint; do not defer it to a historical deadline.

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
