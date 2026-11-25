# Independent review of the incidence coefficient sensitivity study

Andrew Walther · 2026-11-25

The saved run supports the main comparison. The spillover regime still determines which group of designs has the lowest intervention-effect MSE. Small shifts with the incidence coefficient change some exact rankings. Estimating the incidence coefficient is much less precise under High Incidence Focus than under Simple Random Sampling (SRS). Rook Checkerboard requires a separate interpretation because its fitted intervention coefficient absorbs spillover.

## Grid, fits, and reproduction

I reviewed the approved specification, runner, summary script, and saved smoke, pilot, and full outputs. I used the stored data for the checks below. β denotes the intervention coefficient and τ denotes the incidence coefficient.

- The full outputs have 288,000 scenario rows and 2,880,000 surface rows. At each of the five τ values there are 28,800 scenario rows and 288,000 surface rows per coefficient, with both estimators represented. The design and parameter counts agree with the production grid. Every saved surface has 25 attempted and valid fits per coefficient. The smoke and pilot used three workers. The full run used five workers, and the recorded times by τ sum to about 146 minutes, compared with the planning estimate of 75–100 minutes on ten workers.
- I independently matched the saved τ = 1 intervention rows to the dated production files by scenario and surface keys. The keys were unique and complete in each estimator. All stored production columns matched exactly for 14,400 scenarios and 144,000 surface rows per estimator. This verifies the production intervention metrics, including the Wald columns. The new incidence and corrected interval columns have no counterpart in that check.
- The code uses the production incidence, assignment, and error keys, which omit τ. It holds those inputs fixed across the sweep and changes only the X term in the outcome. It reuses the two production model matrices and extracts Z and X from each fit. The corrected interval uses 95 residual degrees of freedom for a full-rank oracle fit and 96 after dropping one oracle column or in a full-rank non-oracle fit. The per-fit handling gives the few aliased Isolation Buffer fits 96 degrees of freedom. The earlier interval script 26 assigned those fits 95. The new study follows its approved specification.
- The 7,200-row performance file keeps neighbors, incidence configurations, regimes and estimators separate. The MSE ratio divides two means within the same settings. The summary also includes the β = 1 slice. The eight PDFs cover the four specified plots by neighbor type. The queen rejection PDF has the five incidence configurations, two regimes, canonical design labels, and no title.

## Intervention-effect comparisons

Across the 20 oracle panels averaging over the full β grid, the MSE-minimizing design changes somewhere along the τ sweep in six panels. Five are switches between Saturation Quadrants and Incidence-Guided Saturation under control-only spillover. The sixth is a switch between Balanced Halves and SRS with rook contiguity, spatial incidence at ρₓ = 0.5, and both-arms spillover. Across the four comparisons with τ = 1 there are nine crossings of the SRS reference, all close to parity. The largest deviation from the reference among those crossings is Isolation Buffer with rook contiguity, spatial incidence at ρₓ = 0.2, and control-only spillover. Its MSE ratio changes from 0.991 at τ = 1 to 1.001 at τ = 2.

The paired root mean square change in β̂ from τ = 1 ranges from 0.001 to 0.067 across oracle groups at the other τ values. The largest change is for Checkerboard with queen contiguity, spatial incidence at ρₓ = 0.2, and control-only spillover at τ = 0. In that group the MSE ratio is about nine at both τ values. The paired root mean square change measures how far the individual estimates move. It can be small even when two designs exchange ranks.

As a check on the size of one rank change, I averaged MSE over the grid within each saved surface for independent uniform incidence, queen contiguity and control-only spillover. The table gives Incidence-Guided Saturation MSE minus Saturation Quadrants MSE and the Monte Carlo standard error calculated as the standard deviation of ten paired surface differences divided by √10. All values are multiplied by 1,000.

| τ | Paired MSE gap | Monte Carlo SE |
|---:|---:|---:|
| 0 | −1.889 | 5.871 |
| 1 | 1.525 | 5.607 |
| 2 | 4.927 | 5.365 |

Pairing the gaps again across τ gives a change of −3.414 (SE 0.749) at τ = 0 versus 1 and +3.402 (SE 0.577) at τ = 2 versus 1 on the same scale. In this panel the shift in the relative gap is detectable because the same surfaces and draws are used. The absolute size of the shift is small compared with the MSE of about 0.15 for either design, and the rank at a given τ is uncertain. Describe stable regime-level conclusions and acknowledge the near-tie changes.

The stability of β̂ has a direct model explanation. If ρ is held at its true value, multiplication by I − ρW transforms the paired outcome difference into a multiple of X. Since X is in the regression, the extra term can change its coefficient without changing the Z coefficient. Estimating ρ allows the fitted Z coefficient to change. The same conditional algebra applies to the non-oracle model, which also includes X. Its omission of spillover can still bias β̂ or τ̂.

## Incidence coefficient estimates and power

For independent uniform incidence, queen contiguity and both-arms spillover, I reproduced the oracle High Incidence Focus summary at τ = 1 by averaging the saved surface metrics. Mean model SE was 0.668, mean scenario SD was 0.683, MSE was 0.467, bias was −0.010 (surface Monte Carlo SE about 0.011), and corrected coverage was 94.62%. SRS had mean model SE 0.336, mean scenario SD 0.355, and MSE 0.127. Two-sided power was 29.34% for High Incidence Focus and 81.10% for SRS. Their mean absolute treatment–incidence correlations were 0.866 and 0.079, respectively.

For a uniform X and a treatment indicator equal to one above its median, the population correlation is √3/2. In the simplified regression with an intercept and Z, residualizing X on Z therefore retains one quarter of its original variance. The resulting coefficient SE doubles. With n = 100 and residual SD 1, the corresponding benchmarks are √(12/100) = 0.346 before adjustment and 0.693 after adjustment. This explains the scale of the observed contrast. Actual precision depends on the variation left in X after adjustment for the full design matrix and on estimating ρ. Treat the alignment plot as descriptive, especially for the spatial and Poisson incidence surfaces.

Across all oracle panels, null rejection was 4.29–6.17% and coverage across τ was 93.83–95.75%. Power ranged from 10.14–31.68% at τ = 0.5 and 79.14–100% at τ = 2. Small bias is a fair description of magnitude. A zero-bias claim would be too strong. For example, with independent uniform incidence, rook contiguity, both-arms spillover and τ = 2, Balanced Halves had bias −0.016 with a surface Monte Carlo SE of 0.004. This is one of many cells. The individual intervals use an approximate corrected critical value and estimate the coefficient of prior incidence in the specified spatial model. Interpret τ as an adjusted association on the structural outcome scale before spatial feedback. Its estimation after a trial still depends on the spillover specification and sufficient independent variation in X.

## Fit diagnostics and the rook Checkerboard coefficient

I counted 36,000,000 attempted model fits by using one coefficient row per estimator and scenario. There were no errors, unavailable estimates, unavailable SEs, or warnings other than alias messages. The warning CSV contains 1,000,100 occurrences of “Aliased variables found: Spill,” all in the oracle. Rook Checkerboard accounts for 1,000,000, which is every one of its oracle fits. Rook Isolation Buffer accounts for the remaining 100. At each τ those occurred once in each of 20 scenarios, corresponding to four distinct allocations repeated at five intervention coefficients. The rook Checkerboard and Isolation Buffer counts include repetitions across the incidence-coefficient sweep. The saved scenario diagnostics count aliases, and the surface rows retain coefficient denominators and missingness. There are no per-surface alias or warning counters.

Under rook Checkerboard every neighbor of a treated cluster is untreated and every neighbor of an untreated cluster is treated. With row-standardized weights, WZ = 1 − Z. The spillover term is ψ(1 − Z) in either regime, since (1 − Z)² = 1 − Z. Dropping Spill leaves an intercept shifted by ψ and a Z coefficient equal to β − ψ. X remains estimable. In the independent uniform τ = 1 oracle rows, the mean bias relative to β is −0.651 with both-arms spillover and −0.661 with control-only spillover. Coverage is 25.65% and 25.02%, respectively. Keep these metrics as diagnostics of the design’s inability to identify β separately from spillover. Label the rook Checkerboard intervention curve accordingly so readers see why a complete numerical fit can still fail to estimate β.

## Changes before reporting or reusing the workflow

1. Report which near-ties change rank and use paired surface gaps when an exact ordering matters. The current CSV gives uncertainty for each design and for its change from τ = 1. It does not give uncertainty for a paired contrast between two designs or for a design-to-SRS ratio. The saved surfaces support those calculations.
2. Include one line identifying the aliased coefficient and reporting the rook Checkerboard and Isolation Buffer counts. If the workflow is reused for a design with substantial selective fit loss, publish the surface-level failure and alias counters and distinguish estimates available for point summaries from estimates with SEs available for intervals. The present code requires finite estimates and SEs for both and will reject an undefined grouped MSE. The checks showed complete denominators in this run.
3. Keep the production reproduction report with the outputs. For future reruns, tie the assembly check to the actual saved results and record the summary script and plotting dependency versions as well. The checker compares τ = 1 before the remaining fits. The summary script subsequently trusts the final line of the report, which can also exist beside stale result files if a later sweep stops before publishing. This was not a problem in the saved run and can be addressed by recording output hashes in a completion record published after all results.

## Resolution

The findings now report the near-tie rank changes, the paired saturation comparison, and the meaning of the aliased rook Checkerboard coefficient. The rook intervention plot labels that curve as aliased. The summary provenance file records the versions and checksums of the inputs and plotting code used for this run. The independent reproduction check above verified the saved output rows again after assembly.
