# NC SUD implementation and completion plan — 2026-10-02

**Verified completion checkpoint (2026-10-02):** the observed-incidence application
and focused tail refinement are complete. Production: 1,248 reporting rows,
936 computational sources, 7,598,000 independent outcome fits; separate tail refinement
adds 2,280,000 outcomes on the same allocations. All settings pass completeness
and the approved mean-MSE/coverage Monte Carlo gates; no warnings, aliases,
boundary fits or failures. Coverage remains below 95% and allocation tails remain
uncertain. See [completed findings](nc_sud_application_findings_2026-10-02.md).
Full chapter/appendix and derived CTJ/SI drafts are complete and verified; author review and submission metadata remain. See [manuscript completion findings](manuscript_completion_findings_2026-10-02.md).

Earlier execution/pending statements below are historical and superseded by this checkpoint.

Status: author approved implementation and uninterrupted FAST execution, with
an HTML visual companion added as a deliverable. Implementation and verification
are underway; production results and manuscripts are not complete. Authority:
`nc_sud_application_plan_2026-10-02.md`; review:
`nc_sud_independent_review_2026-10-02.md`.

## Study being implemented

- Fixed observed 2018–2021 yearly SUD rates for 58 named service-area clusters;
  no synthetic incidence or incidence resampling.
- Primary education model: Y = (I − ρW)^−1[τZ + S(Z) + ε], τ = 1, ε ~ N(0,I).
  Fit SAR ML with intercept, Z and true regime-specific Spill. S(Z) is γWZ
  for both arms or γ(1−Z)WZ for control-only. No incidence term in primary outcomes
  or fitted model; τ is constant across clusters within a scenario.
- Queen ρ = {0,0.2,0.5}, γ = {0.5,0.8}, both regimes, four years, eight designs
  plus SRS: 432 reporting blocks. Named cached W remains authoritative.
- Matched baseline/adjustment sensitivity: queen same settings; add βX, β = 1,
  X = (average rank(rate)−0.5)/58; include X in the fitted matrix. Do not attribute
  differences solely to β or call this an observed incidence–education relationship.
- Rook sensitivity: primary education specification at ρ = {0,0.5}, γ = {0.5,0.8},
  both regimes/four years/nine designs: 288 reporting blocks.
- Regional-summary sensitivity: Design 8 under the primary queen settings;
  alternative mean cluster rates and pooled regional deaths/population. Verify
  ordering/tie structure first; identical allocation rules reuse results with
  explicit provenance. At most 96 alternative reporting blocks before reuse.
  No crossing all sensitivities or additional grid/estimator sweeps.

## 1. Validate and freeze inputs and designs

Add `application/code/real_sud_setup.R`: strict year/cluster-name joins; check
finite rates, positive denominators, totals, named W ordering, adjacency and row
sums. Use cached legal-boundary geography for coordinates/maps and verify it
against the cached cluster identities/weights. Do not rebuild different weights.

Save one four-region map chosen by the existing geographic/balance criterion,
using pooled person-years/4. Save one block map per year using existing standardized
location/incidence-rank/log-population grouping. Freeze both before performance
evaluation. Record definitions, seeds, hashes, coordinates and diagnostic tables.

Surgically update `application_designs.R` to accept frozen block membership and
explicit regional summaries, and to assign 14/15 in Balanced Halves with a random
choice of the half receiving 15. Preserve other design support and tie rules.
Primary Design 8 averages cluster ranks within regions and assigns 80/60/40/20%
from highest to lowest regional mean, randomizing regional ties per draw.

Verify: exact 29 totals for SRS/High Incidence Focus/Quartiles/Halves; half-specific
14/15 rule; region saturation and buffer independence; canonical graph order;
block membership invariant across scenarios/draws; sensitivity only changes its
specified regional summary. Detect proven singleton support and retain sampled
duplicate frequencies for every other support.

## 2. Build a separate reproducible runner and uncertainty summaries

Add focused files under `application/code/`: `real_sud_simulation.R` for nested
allocation/outcome simulations and explicit model matrices;
`run_real_sud_study.R` for smoke/pilot/production orchestration;
`real_sud_summary.R` for performance, allocation risk and traceable exhibits.
Reuse the unchanged validated `fit_one_lag_model()` and verified allocation-risk
arithmetic; keep grid modules 01–05 and completed results unchanged. Document
every new/modified function and follow application script-header conventions.

Key-seed every assignment and fixed-size outcome batch. Preserve prefixes when
extending allocation/outcome counts; record unique assignments and draw frequencies.
Use independent noise streams across design supports so SRS difference MC SEs
do not silently omit covariance. If paired noise is used for summary-method
sensitivity, compute its paired uncertainty explicitly. Derive equivalence keys
from actual support and outcome/model inputs, not merely labels.

For primary Designs 1/3/4/9, prove cross-year support/outcome equivalence and reuse
performance with shared source IDs. Continue reporting annual treatment population
shares/incidence diagnostics separately. Shared rows are not independent yearly
replicates. Do not apply this reuse to β = 1 without checking the full model inputs.

Write new manifests/checkpoints under
`application/results/real_sud_rev_20261002/`: input/code/package/BLAS/parameter
hashes, RNG configuration and frozen designs. Refuse incompatible resumes.
Preserve older synthetic outputs. Keep individual fits/checkpoints ignored;
only permitted cluster-level/aggregate outputs are candidates for tracking.

Retain errors, warnings, bounds, aliased columns, failed/nonfinite SEs and fit
statuses. Test structural τ identification using the actual fitted matrix and
residual treatment variation. Report physical independent-fit counts, unique
allocations and frequency-weighted estimands separately. Incomplete fits yield
explicit incomplete summaries; never hide failures with `na.rm`.

Verify: exact DGP calculations for both regimes; singleton versus sampled-unique
distinction; correct primary/sensitivity matrices; cache frequency weighting and
joint MC SE for MSE/coverage/bias; serial/worker/order/prefix reproducibility;
manifest rejection; unchanged protected grid hashes.

## 3. Tests, smoke, pilot and production with precision gates

Add `application/tests/test_real_sud_study.R`; retain and run the existing data
tests and applicable simulation/allocation-risk behavior suites. Test actual
statistical calculations and assignments, not only output shapes. Cross-check
selected queen/rook primary and β = 1 fits against `spatialreg::lagsarlm` (τ and
SE within the validated tolerance), including high-incidence and graph cases.
Check the exact buffer invariant: Z ⊙ WZ = 0 makes both regimes' responses match
for the same noise and β specification. Surface convergence warnings.

Run smoke, then a focused pilot spanning years, all designs, regimes, neighbor
types and model matrices. Benchmark runtime; inspect identification, warning/failure
rates and precision before production. Use local workers with one BLAS thread each;
workers process bounded blocks, with fit batches/checkpoints to limit memory.
No Longleaf submission is implied.

Production starts at J = 100 allocation draws and R = 100 outcomes per unique
allocation for stochastic supports. Proven singleton supports start R = 1,000
(one eligible allocation), not 100 duplicated draws treated as 10,000 outcomes.
Gate each reporting block on relative MCSE(MSE) ≤5% and MCSE(coverage) ≤0.01.
Use verified joint variance with cache covariance for both targets.

Prespecify refinement tiers: J = 100/200/400 and R = 100/200/400/1,000/2,000/4,000.
Increase J when allocation uncertainty dominates and R when outcome uncertainty
dominates; singleton supports increase R only. Preserve old fits/draw prefixes.
Stop refinement on achieved precision, not desired ranks. If targets remain unmet
at these tiers, flag the setting and discuss a focused extension before claiming
verification complete. Record achieved counts/SEs/status in a completeness report.

Allocation-risk comparisons report estimated q90/worst-decile means, signed
variance corrections and split-half/cross-selected diagnostics. Use R400 or a
higher already-required tier for selected saturation/SRS/balanced/buffer contrasts
where outcome noise affects interpretation. More outcomes do not replace more
tail allocations. For close strong superiority claims, choose a final fixed tier
from pilot diagnostics and run fresh seeded confirmation, or retain qualified
descriptive conclusions. No guarantee of exact true tails is required or claimed.

Verify: all planned reporting blocks accounted for, shared provenance identified,
precision gates met or explicitly unresolved, fit failures/identification issues
shown, old outputs preserved, no restricted county data staged. Review current
findings before publishing recommendations; no design is assumed to win.

## 4. Complete results and the thesis chapter/appendix

Produce yearly/regime MSE, bias, SD, coverage and SRS comparisons; treatment counts,
population shares and incidence balance; exposure/identification diagnostics;
allocation-risk and sensitivity results with achieved MC precision. Write a
source-linked results/exhibit extract and an explicit body/appendix exhibit inventory.

Revise only `paper/dissertation_chapter/Dissertation_Chapter.qmd` and its local
assets. Explain the approved education application, both incidence roles in the
grid study, changed β/τ notation from Chapter 2, and matched model sensitivity.
Use the approved Chapter 2 extension link; include Checkerboard/SRS in the main
Project 2 comparison, with detailed diagnostics in the appendix. Keep grid-pilot
and full-study/application evidence separate; no numerical Project 1 comparison
or audit reopening is needed.

Resolve author-confirmation markers through evidence or focused author questions;
do not invent intervention/source/IRB details. Keep substantive findings in the
body and reproducibility/mechanical detail in the appendix. Render and conduct
fresh numerical, methodological and visual review; trace every number/exhibit to
current outputs. Use only the approved Chapter 3 sync hook for bios-dissertation.

Verify: complete application methods/results and appendix, supported statements,
all unresolved author items identified/resolved, correct floats/references/numbering,
clean render and verified sync. Update README/AGENTS/TODO/ROADMAP at milestones.

## 5. Derive and finish CTJ and supplementary material promptly

Derive the full standalone CTJ account from the agreed chapter; preserve omitted
material in the supplement. Verify live journal requirements and the Project 1
publication before citing it. Replace all April/synthetic numbers and exhibits;
write the required abstract, references, declarations, data/code availability
and submission checklist. Obtain author-only details when needed; do not fabricate
declarations. Keep restricted sources private and identify data-use approval items
before actual publication/submission.

Render and review CTJ/SI in fresh passes and perform a consolidated chapter/CTJ/SI
consistency check. The LaTeX manuscript is a multi-file SAGE project; use its existing
build workflow and preserve source rather than assume standalone-editor support.
The author wants ASAP completion; no historical deadline gates CTJ derivation.

Verify: full chapter/appendix and full CTJ/SI drafts, current traced numbers,
standalone framing with verified Project 1 citation, live requirements satisfied,
and remaining author/submission actions clearly identified. Actual submission is
not authorized by implementation-plan approval. Commit logical steps with specific
paths and no AI co-author lines; never push without explicit permission.
