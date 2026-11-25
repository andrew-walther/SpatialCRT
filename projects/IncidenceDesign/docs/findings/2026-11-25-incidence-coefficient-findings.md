# Incidence coefficient sensitivity on the spatial grid

Author: Andrew Walther. 2026-11-25.

## Method

We varied the incidence coefficient τ on X over 0, 0.25, 0.5, 1, and 2. In the R code this is `beta`. The intervention coefficient β on Z is `tau` in the code and took its five production values in every sweep. X ranges from 0 to 1 and the residual SD is 1. A value of τ = 1 therefore contributes one residual SD across the full incidence range before spatial feedback. τ is a regression coefficient and can take either sign.

The 100-cluster grid used nine designs, queen and rook neighbors, both spillover regimes, and all production settings of incidence, outcome spatial correlation, and spillover magnitude. Each incidence configuration had ten surfaces, with 25 assignment and outcome draws per surface. The keyed incidence, assignment, and noise draws matched across τ. The oracle spatial lag model included Z, the true form of spillover, and X. The non-oracle sensitivity fit included Z and X. Both supplied the coefficient and fitted standard error for X from the same model used for Z.

We used the corrected 95% interval with half-width `qt(.975, N - p) × SE × sqrt(N / (N - p))`, where p counts the retained regression coefficients and the spatial lag parameter. The degrees of freedom reflect each fit's retained columns, including the rare Isolation Buffer aliases. Power for τ was the proportion of corrected intervals that excluded zero in either direction. We compared treatment MSE with the mean MSE for Simple Random Sampling (SRS) in the same setting. The summaries average scenarios equally within each design, incidence configuration, neighbor type and regime. The figures average over the five β values, and the table also gives β = 1 separately. Monte Carlo standard errors for summary metrics use the ten surface means within each configuration.

## Verification

The full sweep took 145.8 minutes with five forked workers and one BLAS thread per worker. At τ = 1, all stored β̂ metrics from the September 27 run matched exactly for 28,800 scenario-estimator rows and 288,000 surface-estimator rows. The scenario and surface summaries cover 36 million model fits, counting each estimator once per draw. All fits returned both estimates and their standard errors. There were 1,000,100 spillover aliases in oracle rook fits across the five τ values. Of those, 1,000,000 were Checkerboard and 100 were Isolation Buffer. For rook Checkerboard, WZ = 1 − Z gives a retained Z coefficient of β − ψ in both regimes. Its intervention MSE and coverage still use β as the target, so they reflect this identification failure. The X coefficient remains identified. The fit log contained zero errors and zero other warnings. [Reproduction check](../../results/incidence_coef_2026-11/full/reproduction_check.txt) · [Fit messages](../../results/incidence_coef_2026-11/full/fit_messages.csv)

## Intervention effect β̂

Changes in β̂ were small relative to its sampling variation. Averaging over the intervention grid, the largest root mean square change from the paired τ = 1 estimates was 0.067. It occurred for queen Checkerboard with control-only spillover and spatial incidence at ρₓ = 0.2, where the mean scenario SD was 1.44. The largest proportional change in a design's mean MSE was 3.4% over the full intervention grid and 3.7% in the β = 1 subset. For queen High Incidence Focus with control-only spillover, setting τ to zero lowered its MSE by 2.5–3.4% across incidence configurations. Incidence-Guided Saturation gained 1.7–2.2% in those settings.

The spillover regime still determined which designs led the MSE comparison. With control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest mean MSE in every oracle comparison. The best-to-SRS ratios ranged from 0.617 to 0.715 for queen and 0.760 to 0.823 for rook. With both-arms spillover, a balanced design led all queen comparisons. The best-to-SRS ratios were 0.961–1.000 for queen and 0.950–1.000 for rook. The design with the lowest mean MSE changed across τ in six of the 20 configuration-by-neighbor-by-regime comparisons. In five of these, the two saturation designs exchanged first place. Nine design ratios switched sides of one relative to τ = 1, and all of those ratios remained within 0.01 of one. [Queen MSE figure](../../results/incidence_coef_2026-11/full/summary/intervention_mse_vs_srs_queen.pdf) · [Rook MSE figure](../../results/incidence_coef_2026-11/full/summary/intervention_mse_vs_srs_rook.pdf)

At the true spatial lag parameter, adding Δτ to the data-generating incidence coefficient adds only ΔτX to the spatially filtered outcome. Both fits already adjust for X. The change can also affect the fitted spatial lag parameter, which helps explain why paired β̂ values can move.

The paired simulation resolves some small changes. For queen uniform incidence with control-only spillover, the mean MSE gap between Incidence-Guided Saturation and Saturation Quadrants decreased by 0.0034 at τ = 0 relative to τ = 1 and increased by 0.0034 at τ = 2. The paired Monte Carlo SEs for these changes were 0.00075 and 0.00058. The individual design gaps at each τ were smaller than their own Monte Carlo SEs. [Independent review](../reviews/2026-11-25-incidence-coefficient-review.md)

## Incidence coefficient τ̂

The largest absolute oracle bias among the summaries was 0.023. The corrected intervals covered τ in about 93.8–95.8% of draws across the incidence configurations, designs and regimes. At τ = 0, two-sided rejection ranged from 4.3% to 6.2%. In the non-oracle fit, bias reached 0.137 at τ = 0. Omitting the spillover covariate can shift the fitted association with incidence.

The model-based SE depended much more on design than on the value of τ. Across the oracle summaries, the largest change in mean model SE from τ = 1 was 0.0054. High Incidence Focus had the largest mean model SE and MSE in every comparison. Checkerboard and the two balanced designs often had the smallest MSE. On the queen grid with uniform incidence, mean model SE was 0.667–0.670 for High Incidence Focus and 0.336–0.337 for SRS. Their rejection rates below span the two spillover regimes in that setting. [Queen rejection figure](../../results/incidence_coef_2026-11/full/summary/incidence_rejection_queen.pdf) · [Rook rejection figure](../../results/incidence_coef_2026-11/full/summary/incidence_rejection_rook.pdf)

| True τ | High Incidence Focus | Simple Random Sampling |
|---:|---:|---:|
| 0 | 4.62–5.43% | 4.95–5.76% |
| 0.25 | 6.34–6.35% | 10.71–11.54% |
| 0.5 | 11.05–11.13% | 29.79–29.94% |
| 1 | 29.34–30.30% | 81.10–81.92% |
| 2 | 81.32–81.51% | 99.98–100.00% |

For uniform X and residual SD 1, the simple regression approximation gives an SE of about √(12/100) = 0.35 when treatment is unrelated to incidence. Splitting treatment at the median leaves about one-quarter of the original variance of X within treatment groups, which doubles the approximation. In the simulations, the mean absolute Z–X correlation for High Incidence Focus ranged from 0.856 to 0.866. For Incidence-Guided Saturation it ranged from 0.097 to 0.119, with mean model SE of 0.339–0.348 at τ = 1. This comparison describes alignment with the individual cluster covariate. Adjustment for spillover and the fitted spatial lag also affect precision. [Queen alignment figure](../../results/incidence_coef_2026-11/full/summary/incidence_se_vs_alignment_queen.pdf) · [Rook alignment figure](../../results/incidence_coef_2026-11/full/summary/incidence_se_vs_alignment_rook.pdf)

## Interpretation and limits

Once trial outcomes are observed, τ̂ describes the association of prior incidence with the later outcome conditional on the treatment and spillover terms in the spatial model. Its interpretation is on the scale of the covariate after scaling X and before spatial feedback. A design suited to estimating β can still give a wide interval for τ if treatment closely follows incidence. The oracle results rely on the correct spillover form and the grid's covariate and outcome models. The simulation used ten surfaces per configuration, a single cluster count, and nonnegative values of τ. The τ = 0 grid fits still included X. The primary North Carolina case-study model omits X, and the grid results give no empirical estimate of the association for that trial.

The author can decide whether the incidence estimation subsection belongs in the curated or full Supplementary Information and whether Chapter 3 should carry the same subsection or a shorter appendix note. The most compact exhibit would pair the queen rejection figure with the range of changes in β̂ MSE. The full results retain both neighbor types and both estimators. [Performance data](../../results/incidence_coef_2026-11/full/summary/performance.csv) · [Runner guide](../../results/incidence_coef_2026-11/README.md)
