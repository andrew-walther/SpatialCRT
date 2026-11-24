# Incidence-Guided Saturation sensitivity study: specification

Author: Andrew Walther. Drafted 2026-11-24 and run overnight the same day.
**Status: provisional.** The choices below were made so the run could go overnight; every
one is open for the author's review, and a rerun with different choices costs about an hour.

## Question

The study tested one version of Incidence-Guided Saturation Quadrants: four 5×5 quadrants
with saturations 80/60/40/20% assigned by the rank of each quadrant's mean incidence. We do
not know whether versions with other numbers of regions or other saturation sets perform
similarly, relative to Simple Random Sampling (SRS) and relative to the original version.
Committee question Q-16 in `presentation/prelim-qa-prep.md` asks this directly.

## Designs

Every variant treats exactly 50 of the 100 grid clusters.

| Factor | Levels | Notes |
|---|---|---|
| Regions R | 2 (5×10 halves), 4 (5×5 quadrants), 10 (2×5 rectangles), 25 (2×2 blocks) | Equal-size partitions only, so 50% can be treated exactly. 9 or 16 regions would need unequal regions. |
| Spread s | 0.1, 0.2, 0.3, 0.4 | sat_r = 0.5 + s(1 − 2(r − 1)/(R − 1)) for incidence rank r = 1 (highest) to R. R = 4, s = 0.3 gives the original {0.8, 0.6, 0.4, 0.2}. |
| Guidance | guided, unguided | Guided: saturations by the rank of the region's mean incidence (ties at random per draw, as in the original design). Unguided twin: the same saturations in random order per draw (as in Saturation Quadrants). |

- Counts per region use largest-remainder rounding of size × saturation, so the total is
  exactly 50; equal remainders are broken at random per draw.
- s = 0.4 with R = 25 gives saturations from 0.9 to 0.1 on blocks of four, so some blocks
  are fully treated or fully control. This is allowed; within-region contrast is not
  required for the spatial lag model.
- The (R = 4, s = 0.3) pair is not duplicated: it is the original design (Design 8, guided)
  and Saturation Quadrants (Design 3, unguided).
- 30 variants (IDs 101–130, listed in `igs_variants.csv`) plus Designs 3, 8 and 9 (SRS).

## Simulation

The full production grid, unchanged: five incidence configurations × 10 surfaces, queen
and rook, ρ ∈ {0, 0.01, 0.2, 0.5}, four spillover magnitudes, both spillover regimes, the
five effect sizes, 25 draws per surface (250 fits per scenario), oracle and non-oracle
spatial lag fits. Production code is sourced, not edited (`code/29_igs_sensitivity.R`
wraps `get_designs()` the way `code/26` wraps `summarize_fits()`). Variant draws have their
own seed keys. Designs 3, 8 and 9 reuse their production keys, so they must reproduce the
2026-09-27 run exactly (checked; the run stops otherwise), and every design in a block sees
the same noise, which makes the comparisons paired.

## Contrasts (named in advance)

Primary, queen, oracle, pooled over the five effect sizes, by spillover regime:

1. Each variant's MSE relative to SRS, MSE_d / MSE_SRS, with the paired test over the 50
   (configuration × surface) units used in `code/15` (t49 interval; Holm across variants).
2. Each guided variant's MSE relative to the original design (Design 8), same paired test.

Secondary: guided minus its unguided twin (paired); t + df coverage; bias; burden share
treated (share of incidence in treated clusters, as in `code/22`). Rook is reported beside
queen; nothing pools the two.

Multiplicity (as run): the 32 designs other than SRS against SRS, the 15 other guided versions
against the original, and the 16 guided-versus-twin pairs, Holm within each family and regime.
The (R = 25, s = 0.1) version turned out to be 2×2 Blocking exactly, so its guided and unguided
forms are one design; it serves as a null control (findings report). With 10 surfaces per configuration the tests detect only
differences that are large relative to surface-to-surface variation; the point is whether
the conclusions change, not to rank every variant.

## Where results go

The SI (one short subsection and one figure or table) and a backup slide for Q-16. Not the
CTJ main text. Summaries: `code/30_igs_sensitivity_summary.R` →
`results/igs_sensitivity_2026-11/full/summary/`.

## Questions for the author

- Are these regions and spreads the right ones? Alternatives: unequal 3×3 or 4×4 regions,
  or sets that are not linear in rank (for example 90/70/30/10).
- Should the study also run on the North Carolina map, where regions follow the
  at-least-two-colleges rule? Not done here.
- Is pooling over the five effect sizes right for the headline, as in `code/15`?
