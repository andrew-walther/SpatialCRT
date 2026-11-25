# Incidence-Guided Saturation sensitivity study: specification

Author: Andrew Walther. The initial study ran on 2026-11-24.
**Status.** Andrew approved the Phase 3a and 3b specification on 2026-11-25.
The smoke run, full grid run, summaries, and figure are complete.

## Question

The main study tested one version of Incidence-Guided Saturation: four 5×5 quadrants
with saturations 80/60/40/20% assigned by the rank of each quadrant's mean incidence.
Do versions with other numbers of regions or saturation sets perform similarly, relative
to Simple Random Sampling (SRS) and to the original design? Committee question Q-16
in `presentation/prelim-qa-prep.md` asks this directly.

## Designs

Every variant treats exactly 50 of the 100 clusters on the 10×10 grid.

| Factor | Levels | Notes |
|---|---|---|
| Regions R | 2 (5×10 halves), 4 (5×5 quadrants), 10 (2×5 rectangles), 25 (2×2 blocks) | Regions are equal in size. |
| Spread s | 0.1, 0.2, 0.3, 0.4, 0.5 | The target saturation is sat_r = 0.5 + s(1 − 2(r − 1)/(R − 1)), where r is the rank from 1 (highest saturation) to R. R = 4 and s = 0.3 give the original 80/60/40/20%. |
| Guidance | Guided, unguided | Guided variants assign saturations by the rank of each region's mean incidence, breaking ties at random on each draw. Unguided twins assign the same saturations in random order on each draw. |

- Counts per region use largest-remainder rounding of size × saturation. The total is
  exactly 50, and equal remainders are broken at random on each draw.
- At s = 0.5 the extreme regions are fully treated and fully control. For four quadrants
  the target levels are 100%, 67%, 33%, and 0% to the nearest percentage point. For two
  halves, one entire half is treated. Rounding can also yield fully treated and fully
  control regions at smaller spreads, as with 25 blocks at s = 0.4.
- The R = 4, s = 0.3 pair consists of the original guided design (Design 8) and
  Saturation Quadrants (Design 3).
- The original 30 variants keep IDs 101–130 and their seed keys. At s = 0.5, IDs 131–134
  are guided and IDs 135–138 are unguided, each in order of R = 2, 4, 10, 25. Each new ID
  has its own seed key. The full study has 38 variants plus Designs 3, 8 and 9 (SRS).
- Each region has a distinct target saturation. The study is grid only.

The 100/0 setting is a boundary case for sensitivity analysis. The expectation is that
under control-only spillover the MSE gain from wider saturations levels off or reverses,
while under both-arms spillover MSE worsens. Compare 100/0 with 90/10 for each region
count and guidance setting to see whether the results follow that expectation.

## Simulation and verification

Use the full production grid of five incidence configurations with ten surfaces each,
queen and rook contiguity, ρ ∈ {0, 0.01, 0.2, 0.5}, four spillover magnitudes, both spillover
regimes, and five intervention effect sizes. Each scenario has 25 draws per surface, yielding
250 datasets for the oracle and non-oracle spatial lag fits. In this specification β is
the intervention effect on Z (code `tau`), and τ is the incidence coefficient on X
(code `beta`).

Production code is sourced by `code/29_igs_sensitivity.R`, which wraps `get_designs()`.
Every design in a block shares the same outcome noise, which makes the comparisons paired.

Smoke-run all eight new variants before fitting the full grid. For every new assignment
draw, check the rank of [1, Z, WZ, X]. In the control-only regime, also check the rank of the
fitted oracle matrix with (1 − Z)WZ as the spillover column. Report the number of draws
checked and any rank deficiencies by variant, together with aliased fits, failed fits,
and warnings. Endpoint saturations are allowed. Surface any incomplete results or warnings.

Verify the stored 2026-11-24 results against their original manifest and use them for
Designs 3, 8, 9 and variants 101–130 if the checks pass. Fit only IDs 131–138 on that path.
Check exact agreement for all retained scenario and surface rows against 2026-11-24 and
for the stored production metrics of Designs 3, 8 and 9 against 2026-09-27. Stop on any
reproduction discrepancy and report whether the run used the stored results.

The 2026-11-24 run took 78.5 minutes. At that fit rate, the estimated time was about
20 minutes for the eight new variants or 100 minutes for a fresh fit of all 41 designs.
The 2026-11-25 run used forked `mclapply` and pinned BLAS to one thread per worker.
It started with ten workers. After about 44 minutes, Andrew approved a restart with
four workers to reduce memory pressure. The restart reused 20 completed boundary
checkpoints and finished the last 20 units in 8.2 minutes.

## Contrasts

The primary results use queen contiguity and the oracle fit, pooled over the five effect
sizes and reported by spillover regime.

1. Each design's MSE relative to SRS, MSE_d / MSE_SRS, with the paired comparison over
   the 50 configuration-by-surface units used in `code/15`. Report the t₄₉ interval and
   the Wilcoxon p-value.
2. Each guided variant's MSE relative to the original design (Design 8), with the same
   paired comparison.

Secondary outcomes include the paired comparison of each guided design with its unguided
twin, t + df coverage, bias, and burden share treated (the share of incidence in treated
clusters, as in `code/22`). Report rook and non-oracle results separately.

Apply Holm within each family for each combination of estimator, contiguity type, and
spillover regime, and separately for the labeled summaries pooled over regimes. The
family sizes are 40 comparisons with SRS (previously 32), 19 guided variants compared
with Design 8 (previously 15), and 20 guided-versus-twin comparisons including Design 8
versus Design 3 (previously 16). For the original 30 variants, report the adjusted
p-values from both family sizes for every applicable contrast. A twin contrast has one
p-value per pair. Identify any comparisons that cross the 0.05 threshold and verify that
the original unadjusted p-values agree with 2026-11-24.

The R = 25, s = 0.1 designs both assign two of four clusters to treatment in every
block, so the pair serves as a null control. With ten surfaces per configuration the
tests detect differences that are large relative to surface-to-surface variation.

## Outputs and figure

The outputs will inform the planned SI subsection and Q-16 backup slide. This rerun produces
the analysis and figure files in `results/igs_sensitivity_2026-11/`. The summary script
`code/30_igs_sensitivity_summary.R` writes to `results/igs_sensitivity_2026-11/full/summary/`.

Label the figure's x-axis with the highest and lowest target saturations as 60/40, 70/30,
80/20, 90/10, and 100/0. Keep the circles and explain in the legend that they mark the
original Incidence-Guided Saturation design and Saturation Quadrants. Use true hyphens
and minus signs and keep the text legible at manuscript size. Save PDF and PNG, with the
title and notes reserved for the caption.
