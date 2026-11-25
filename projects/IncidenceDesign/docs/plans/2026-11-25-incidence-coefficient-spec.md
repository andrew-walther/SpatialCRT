# Incidence coefficient sensitivity study

Author: Andrew Walther. Drafted and approved on 2026-11-25.

## Simulation grid

Sweep the incidence coefficient τ (code `beta`) over **0, 0.25, 0.5, 1, and 2**. With X on [0, 1] and residual SD 1, τ is the change in the outcome before spatial feedback across the full incidence range, in residual SD units. The intervention effect β (code `tau`) stays on the production grid of 0.8, 1, 1.5, 2, and 3.

The full run uses all nine designs, both neighbor types, both spillover regimes, and all production settings of incidence, outcome spatial dependence, and spillover magnitude. Each scenario has ten surfaces and 25 draws per surface. Reuse the production keys for surfaces, assignments, and noise so that comparisons across τ are paired. Keep incidence configurations and neighbor types separate in the summaries. Queen contiguity and the oracle fit are primary.

## Fits and summaries

In `code/31_incidence_coefficient.R`, source `05_run_simulation.R` with `sim_define_only <- TRUE`. Wrap `fit_tau_models()` locally so each of its two `fit_sar_lag()` calls supplies the named `X` coefficient and its `rest.se` alongside the existing `Z` estimate and SE. Use the same model matrices and fit calls, including the oracle's regime-specific spillover column. Reuse 05's keyed draws and outcome arithmetic. Extend the unit loop locally to summarize both coefficients from those fits.

- Save scenario and surface results for β̂ and τ̂, with mean estimate, bias, empirical SD, MSE, mean model-based SE, coverage, and rejection rates. Use surface means for Monte Carlo uncertainty, following 05. Record valid denominators by coefficient and retain every expected scenario and surface key. Count failed fits, aliases, warned fits and warning messages, and report which coefficients or SEs are unavailable. Retain and flag estimable coefficients from aliased fits.
- Use the 95% interval from 26 for each coefficient, with half-width `qt(.975, N - p) × SE × sqrt(N / (N - p))`. Here N = 100 and p counts retained regression coefficients plus ρ. Use 95 degrees of freedom for a full-rank oracle fit and 96 for a full-rank non-oracle fit or an oracle fit with one aliased column. Keep β̂'s production Wald metrics and positive-direction rejection rule for exact comparison. For τ̂, use exclusion of zero in either direction to estimate type I error at τ = 0 and power at the four positive values.
- Measure changes in β̂ with paired differences against τ = 1. Report their root mean square and compare MSE relative to Simple Random Sampling (SRS) across τ within each regime, using the ratio of mean MSEs as in 15. Check whether design rankings and the direction of the SRS comparisons change. For τ̂, compare mean model SE and empirical SD across τ and designs, and relate mean model SE to the absolute within-assignment correlation between Z and X. Treat the SE benchmarks of roughly 0.35 for weak alignment and roughly twice that for High Incidence Focus as checks for the uniform-incidence setting. The spread of X and adjustment for spillover and ρ can change both benchmarks.

## Figures

1. β̂ MSE relative to SRS versus τ, by design and spillover regime, with a reference line at one.
2. τ̂ rejection rate versus τ, with a 5% reference at the null and panels by spillover regime.
3. τ̂ mean model SE versus τ and versus mean absolute Z–X correlation, by design and regime. Use the latter comparison at τ = 1.

Use canonical design labels from `00_design_names.R`, separate queen and rook panels or files, and omit in-figure titles. Use the oracle fit for the figures and include the non-oracle results in the summary tables.

## Run order and verification

1. After approval, run 05's `test` profile as the smoke check at all five τ values. Verify complete keys and counts, the X coefficient extraction, and failure and alias reporting. This profile uses six fits per scenario and serves as a software check.
2. Run 05's pilot profile at all five values. At τ = 1, match both estimators' β̂ scenario and surface rows to the corresponding production subset from `20260927_000502`. Require unique and complete keys and exact equality of stored numerical metrics, including missing-value patterns. Check oracle τ̂ bias against surface-level Monte Carlo uncertainty and report the fit diagnostics. Stop at any reproduction mismatch. Report and pause after the pilot.
3. Start the full run with τ = 1. Repeat the exact check against the entire production grid before running the other four values. Write results and reproduction reports under `results/incidence_coef_2026-11/`, with separate folders for smoke, pilot, and full.

Pin BLAS to one thread before R starts. The smoke and pilot runs used three forked workers to leave capacity for concurrent work. The author approved five workers for the full sweep. The pilot took 6.9 minutes, and the full sweep took 145.8 minutes. Save checkpoints by profile, τ, and production work unit. Check the manifest against the run parameters, code hashes, R and package versions, and BLAS before resuming. A `.gitignore` in the new results folder covers the checkpoint directories.
