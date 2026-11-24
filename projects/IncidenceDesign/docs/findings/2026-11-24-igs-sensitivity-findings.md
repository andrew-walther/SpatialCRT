# Incidence-Guided Saturation sensitivity study: method and findings

Author: Andrew Walther. Run and summarized 2026-11-24. Spec (provisional, awaiting the author's
review): `docs/plans/2026-11-24-igs-sensitivity-spec.md`. Code: `code/29_igs_sensitivity.R`,
`code/29_igs_variant_designs.R`, `code/30_igs_sensitivity_summary.R`. Outputs:
`results/igs_sensitivity_2026-11/full/` (summary tables and figure in `summary/`). Independent
review: `docs/reviews/2026-11-24-igs-sensitivity-review.md` (code correct; this version applies its
corrections to the report).

## Summary

- Varying the number of regions and the spread of the saturation set changed Incidence-Guided
  Saturation's standing against SRS a great deal. Under control-only spillover, two or four large
  regions with a spread of at least 70/30 had 16% to 45% lower MSE than SRS (pooled over effect
  sizes), while ten regions needed a spread of 80/20 or more and 25 small regions (2×2 blocks) were
  never better than SRS.
- The same settings pulled the other way under both-arms spillover. Wide spreads cost accuracy
  there (two halves at 90/10 had almost twice SRS's MSE). Flat sets were close to SRS with two,
  four or ten regions, but not with 25 blocks.
- The original design (four 5×5 quadrants, 80/60/40/20%) sits between these: 33% lower MSE than
  SRS under control-only spillover and 9% higher under both-arms spillover. Two halves at 70/30
  behaved much like it (32% lower and 5% higher).
- Incidence guidance itself changed accuracy little. Under queen contiguity each guided version
  was within 10% of its unguided twin, the same regions and saturations in random order (up to 14%
  in the rook and non-oracle checks), and every difference that survived Holm adjustment was
  against guidance. The differences were all in variance, not bias. Accuracy depended mainly on
  the regions and the saturation set. What guidance changed was where treatment went: guided
  versions treated up to 58% of the total incidence, against 50% for their twins and for SRS.
- With the oracle fit, coverage of the t + df interval was 94% to 95% for every version, under
  queen and rook contiguity. The design ranking under control-only spillover was the same under
  rook and for the fit without the spillover term (whose coverage was poor, as in the main study).

## Method

**Question.** The main study tested one version of incidence guidance. Do versions with other
numbers of regions or other saturation sets perform similarly relative to SRS and to the original
version? This is committee question Q-16 (`presentation/prelim-qa-prep.md`).

**Designs.** Every version treats exactly 50 of the 100 grid clusters.

| Factor | Levels |
|---|---|
| Number of regions R | 2 (5×10 halves), 4 (5×5 quadrants), 10 (2×5 rectangles), 25 (2×2 blocks). Equal-size partitions only, so exactly half can be treated. |
| Spread s of the saturation set | 0.1, 0.2, 0.3, 0.4. Saturations run linearly with incidence rank r from 0.5 + s (highest-incidence region) to 0.5 − s (lowest): sat_r = 0.5 + s(1 − 2(r − 1)/(R − 1)). R = 4, s = 0.3 gives the original 80/60/40/20%. |
| Guidance | Guided: saturations follow the rank of each region's mean incidence (ties at random, as in the original). Unguided twin: the same saturations in random order for each draw (as in Saturation Quadrants). |

Treated counts use largest-remainder rounding of region size × saturation so the total is exactly
50, with equal remainders broken at random. With R = 25 and s = 0.4, a block's saturation of 0.9 or
0.1 rounds to all or none of its four clusters. With R = 25 and s = 0.1 the rounding treats exactly
two of four clusters in every block on every draw, so that version is 2×2 Blocking exactly, and its
guided and unguided forms are the same design. It is kept as a null control for the twin
comparisons. Because the saturations are linear in rank, their standard deviation at a fixed
spread falls as R grows, s√((R + 1)/(3(R − 1))): s for two regions, 0.75s for four, 0.64s for ten
and 0.60s for 25 (before rounding). The number of regions therefore also changes the saturation
contrast and the region shape (the halves and the 2×5 rectangles are elongated). The (R = 4, s = 0.3) pair is the production
Incidence-Guided Saturation design (Design 8) and Saturation Quadrants (Design 3). That gives 30 new
versions plus Designs 3, 8 and 9 (SRS), 33 in all.

**Simulation.** The full grid benchmark was rerun unchanged except for the designs. That means five
incidence configurations with ten surfaces each, queen and rook contiguity, ρ ∈ {0, 0.01, 0.2, 0.5},
ψ ∈ {0.5, 0.6, 0.7, 0.8}, both spillover regimes and β ∈ {0.8, 1, 1.5, 2, 3}. Each scenario had 25
draws per surface, 250 simulated datasets, each fit with the oracle and non-oracle spatial lag
models. That is 105,600 scenario rows, 13.2 million datasets and 26.4 million fits. Production code was sourced, not edited. The new designs have
their own seed keys, and every design in a block saw the same outcome noise, so comparisons are
paired. Designs 3, 8 and 9 reused their production keys and reproduced the 2026-09-27 production run
exactly (9,600 rows, zero difference in every metric). No fit failed, none was aliased, and none
gave a warning. The run took 78.5 minutes on 10 cores.

**Contrasts** (named in advance). The primary contrast is each version's MSE relative to SRS, as
in `code/15`. It is computed as the paired ratio 1 + D̄ / MSE_SRS, where D̄ is the mean of the
differences over the 50 independent configuration-by-surface units, with a t49 interval and a
Wilcoxon test Holm-adjusted across the 32 designs other than SRS (the spec said 30 variants; the
family also holds Designs 3 and 8). The second contrast is each of the 15 other guided versions
against the original design. The third is each of the 16 guided versions, the original included,
against its unguided twin (Design 3 for the original). Each family is adjusted within a regime. Values
are for queen contiguity and the oracle fit, pooled over the five effect sizes as in `code/15`,
by spillover regime. Rook and the non-oracle fit are reported as checks. Burden share treated is
Σ X_j z_j / Σ X_j (as `code/22`), from the actual allocations in one block per configuration.

## Findings

Figure: `results/igs_sensitivity_2026-11/full/summary/fig_igs_re_vs_srs.pdf` (MSE relative to SRS
by spread, one line per number of regions, guided and unguided side by side, both regimes; the
original design and Saturation Quadrants are circled). Full tables: `igs_summary.txt` and
`igs_contrasts.csv` in the same folder.

### MSE relative to SRS, queen, oracle fit (guided versions)

Ratios below 1 mean lower MSE than SRS. An asterisk marks a Holm-adjusted Wilcoxon p ≥ 0.05, a
difference from SRS that the 50 units do not detect. The figure's bars are unadjusted t49
intervals, so the two can disagree near the line (for example 10 regions at 70/30 under
control-only spillover).

| Regions | Control only: 60/40 | 70/30 | 80/20 | 90/10 | Both arms: 60/40 | 70/30 | 80/20 | 90/10 |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.899 | 0.676 | 0.555 | 0.551 | 0.977* | 1.051* | 1.276 | 1.940 |
| 4 | 1.040* | 0.837 | **0.667** | 0.588 | 1.005* | 1.018* | **1.088** | 1.218 |
| 10 | 1.267 | 1.041* | 0.887 | 0.767 | 1.037* | 1.003* | 1.020* | 1.077 |
| 25 | 2.674 | 1.514 | 1.273 | 1.048* | 1.463 | 1.196 | 1.140 | 1.100 |

SRS MSE pooled over the five effect sizes was 0.227 under control-only spillover and 0.043 under
both-arms spillover. The original design is in bold. The unguided twins gave nearly the same
table, for example 0.668 and 1.047 for Saturation Quadrants.

**Regions and spread under control-only spillover.** Fewer, larger regions and a wider spread both
lowered MSE relative to SRS. Two halves reached 0.555 at 80/20 and 0.551 at 90/10, and four
quadrants reached 0.588 at 90/10. With ten regions the advantage needed a spread of 80/20 or
more. With 25 blocks of four clusters no spread beat SRS. At 60/40 that version is 2×2
Blocking exactly, and its ratio (2.674) differs from the production 2×2 Blocking ratio (2.579, same
slice) only by Monte Carlo error. Part of the effect of R is contrast: at a fixed spread more
regions means less variation between saturations. But contrast is not the whole story. Twenty-five
blocks at 90/10 have about the same saturation standard deviation as the original four quadrants
(0.24 against 0.22) and still had 1.048 times SRS's MSE against 0.667. Small regions intersperse
treatment locally, which the main study found costly when spillover reaches only controls.

**Under both-arms spillover the spread worked the other way.** Flat sets (60/40) were close to SRS
for two, four and ten regions (25 blocks at 60/40, which is 2×2 Blocking, had 1.463). Wider spreads raised MSE, most for two halves (1.940 at 90/10).
Twenty-five small blocks were worse than SRS at every spread.

**The tradeoff.** No version was best under both regimes. Wider contrast between saturations helps
when spillover reaches only controls and hurts when it reaches both arms. Among the versions tested,
two halves at 70/30 and the original four quadrants at 80/20 gave similar results: 32% and 33% lower
MSE than SRS under control-only spillover, and 5% and 9% higher under both-arms spillover. Steeper
sets traded a larger control-only gain for a larger both-arms loss. This matches the main study's
regime-based guidance and adds that the size of the gain depends on the saturation contrast. It
does not single out one recommended version. The size of the gains depends on the effect size,
because the results pool five values of β: the original design's ratio to SRS under control-only
spillover ran from 0.618 at β = 0.8 to 0.718 at β = 3. The ranking of the 32 designs was stable
across effect sizes.

**Against the original design.** Under control-only spillover only two halves at 80/20 or 90/10
(0.833, 0.826) and four quadrants at 90/10 (0.882) had lower MSE than the original; two halves at
70/30 was level with it (1.013), and every other guided version was worse. Under both-arms
spillover every version with two to ten regions and a spread of 70/30 or less, and every ten-region
version, had lower MSE than the original, which itself was 9% above SRS there.

**Guidance against its unguided twin.** Under queen contiguity with the oracle fit, guided-to-twin
MSE ratios ran from 0.98 to 1.08 under control-only spillover and from 0.96 to 1.10 under both-arms
spillover. The null control, which is the same design twice, gave 1.036 [1.007, 1.065] (unadjusted
p = 0.031), so differences of 3% to 4% are within Monte Carlo error. After Holm adjustment the
differences that remained all went against guidance: under control-only spillover R = 25 at 90/10
(1.076); under both-arms spillover R = 25 at 70/30 to 90/10 (about 1.09) and R = 10 at 90/10
(1.076). Rook and the non-oracle fit showed more and larger penalties of the same kind, up to 1.111
and 1.137. In every twin pair the MSE ratio equals the variance ratio to within 0.006, so the
difference is precision, not bias. That is consistent with the grid's τ = 1: incidence shifts the
outcome, and treatment that follows incidence is partly confounded with it, which costs precision
when the regions are small.

**Where treatment went.** Guided versions treated a larger share of total incidence than their
twins (0.500) and SRS (0.500). The share rose with the spread and, at wide spreads, with the number of regions, from 0.507 (two
halves at 60/40) to 0.577 (25 blocks at 90/10). The exception is 25 blocks at 60/40 (0.500): with
four clusters per block, rounding treats two in almost every block, so guidance has nothing to move. The original treated 0.527, between
0.522 and 0.530 across the five configurations. For comparison, the stored North Carolina result
(methods extension A1, frozen regions) is a rate share of 0.511 for Incidence-Guided Saturation
against 0.504 for plain saturation (`results/grid_burden_share/grid_burden_share_log.txt`). This
study did not run the North Carolina map.

**Inference.** With the oracle fit, t + df coverage was 0.940 to 0.953 for every version, regime
and contiguity type. Bias stayed small but grew with the spread under control-only spillover,
reaching 0.086 for two halves at 90/10 (0.087 for its unguided twin).

**Checks.** Under rook contiguity the control-only ordering was the same, with smaller gains (the
original 0.810 against SRS); under both-arms spillover the rook ordering agreed less closely with
queen. With the non-oracle fit, which omits the spillover term, the pattern of ratios to SRS was the
same under both regimes (for example, two halves at 90/10: 0.642 under control-only spillover,
2.364 under both-arms spillover), but its coverage was poor, as in the main study.

## What this changes

- The answer to Q-16 is no longer "future work". The 80/60/40/20% quadrant design is one point on
  a range, and its advantage under control-only spillover depends on large regions and contrasting
  saturations. A flatter set or many small regions loses it, and a steeper set increases it at a
  cost under both-arms spillover. The number of regions, the saturation contrast and the region
  shape move together in this design, so the study shows what combinations work, not which of the
  three drives the result.
- It supports the main text's framing that incidence guidance roughly matches the accuracy of
  unguided saturation while targeting high-burden regions. On the grid, guidance changed MSE by up
  to 10% under queen contiguity, always as a precision cost when it was detectable, and moved more
  treatment to higher-incidence regions.

## Limitations of this study

- Grid only. The North Carolina map, with irregular regions and the at-least-two-colleges rule, was
  not run.
- Equal-size partitions and saturation sets linear in rank only. Unequal regions (3×3, 4×4) and
  non-linear sets (for example 90/70/30/10) were not tested.
- Pooled over the five effect sizes, as `code/15` does. Larger effects raise every design's MSE
  through the error in ρ̂ (SI Section 6), so pooling weights the larger effects more.
- Ten surfaces per configuration; the paired tests detect only differences that are large
  relative to surface-to-surface variation, and the null control shows that twin differences of
  3% to 4% can arise by chance.
- The number of regions is confounded with the saturation contrast and the region shape (see
  Method). A follow-up could plot each version's ratio to SRS against its realized saturation
  standard deviation, or against a measure of how clustered treatment is.

## Proposed placement (for the author)

- SI: a methods subsection on the variants (after the design specifications) and a results
  section with the figure and the table above (the SI review's outline already reserves both).
- CTJ: optional, one sentence in Results with an SI pointer, if the word budget allows (about 35
  words; the recommended cut set leaves about 115). Draft: "In a sensitivity analysis that varied
  the number of regions and the saturation set, the control-only advantage over SRS required large
  regions and contrasting saturations, and incidence guidance changed MSE by at most 10% relative
  to the same design without guidance (SI Section X)."
- Chapter 3: an appendix subsection mirroring the SI.
- Prelim backup slide for Q-16: the figure, with one line on the tradeoff.
