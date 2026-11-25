# Incidence-Guided Saturation sensitivity study

Author: Andrew Walther. Initial run 2026-11-24. Boundary extension 2026-11-25.
The [specification](../plans/2026-11-24-igs-sensitivity-spec.md) records the Phase 3a and 3b
decisions. The [output folder](../../results/igs_sensitivity_2026-11/full/) contains the
simulation summaries and identification checks. The independent review recomputed the
contrasts and checked the figures and reproduction against the stored runs.
[Independent review](../reviews/2026-11-25-igs-sensitivity-rerun-review.md)

## Method

The main study's Incidence-Guided Saturation design uses four 5×5 quadrants with target
saturations of 80%, 60%, 40%, and 20%. We varied the number of regions and the distance
between the highest and lowest target saturations to answer committee question Q-16.
The four equal-size partitions have 2 halves, 4 quadrants, 10 rectangles, or 25 blocks.
The target for rank r in R regions is sat_r = 0.5 + s(1 − 2(r − 1)/(R − 1)), with s
from 0.1 to 0.5 in steps of 0.1. Guided designs assign these targets in order of mean
regional incidence, from highest to lowest, with random tie-breaking. Unguided twins
randomly permute the same targets. Largest-remainder rounding sets the total to 50 of
the 100 grid clusters. When R = 25 and s = 0.1, every block has two treated clusters,
so both versions are 2×2 Blocking and serve as a null control for the twin comparison.

The s = 0.5 extension adds a guided and an unguided variant for each partition (IDs
131–138). The four-quadrant targets are 100%, 67%, 33%, and 0% to the nearest percentage
point, which gives treated counts of 25, 17, 8, and 0. With two halves, one half is
fully treated and the other fully control. When the halves have distinct mean incidence,
the guided assignment is fixed on that surface. Each draw within a block still has its own
outcome noise vector. The extension brings the study to 38 variants plus
Saturation Quadrants (Design 3), the original guided design (Design 8), and Simple Random
Sampling (SRS, Design 9).

The full grid crosses five incidence configurations with ten surfaces each, queen and rook
contiguity, ρ ∈ {0, 0.01, 0.2, 0.5}, ψ ∈ {0.5, 0.6, 0.7, 0.8}, both spillover regimes,
and β ∈ {0.8, 1, 1.5, 2, 3}. The intervention effect β on Z is `tau` in code, and the
incidence coefficient τ on X is `beta` and equals 1. Each scenario has 25 assignments
per surface. Each dataset is fit with the oracle and non-oracle spatial lag models.
New variants have their own
assignment seed keys, and every design in a grid block shares the same outcome noise.
The combined results have 131,200 scenario rows, 1,312,000 surface rows, and 32.8 million
model fits across the two runs.

We report ratios of mean squared error (MSE) with queen contiguity and the oracle fit, pooled over the
five incidence configurations and the five intervention effect sizes within each spillover
regime. Each test uses 50 paired configuration-by-surface differences in MSE. The figure scales the 95% t₄₉ interval for the mean paired difference by the
comparator MSE and adds 1 to express it on the ratio scale. The two-sided Wilcoxon tests
use Holm correction separately for each estimator, contiguity type, and spillover slice.
The expanded families have 40 comparisons with SRS (was 32), 19 with the original guided
design (was 15), and 20 guided-versus-twin comparisons (was 16). The summary also keeps
labeled results pooled over spillover regimes.

## Reproduction and identification

The full runner verified the November 24 checkpoint manifest and reused Designs 3, 8,
and 9 and variants 101–130. All 105,600 original scenario rows and 1,056,000 original
surface rows matched exactly. Designs 3, 8, and 9 also matched all 26 stored columns
in the 9,600 September 27 production rows. The smoke run used ten forked workers.
The full run also started with ten. Andrew approved resuming with four after the first
20 of 40 units because memory pressure was slowing the machine. The same manifest
validated the saved boundary checkpoints, and the last 20 units took 8.2 minutes.
Both runs pinned BLAS to one thread per worker.

The smoke run checked 960 assignment draws and 3,840 fits. In the full run, each new
variant contributed 80,000 assignment draws and 800,000 fits across the two estimators.
All 640,000 draws had rank 4 for [1, Z, WZ, X] and for the oracle design matrix with
the spillover term used in that regime. The 6.4 million new fits had zero aliases,
failures, or warnings. The identification table breaks the counts down by variant,
contiguity, regime, and estimator. Draw counts are shared across estimators.
[Identification table](../../results/igs_sensitivity_2026-11/full/igs_identification.csv)

## MSE relative to SRS

Values below 1 indicate lower MSE. An asterisk marks Holm-adjusted p ≥ 0.05 against
SRS. Bold marks the original guided design. SRS MSE was 0.227 under control-only
spillover and 0.043 under both-arms spillover. The figure shows guided and unguided
variants, the paired intervals, and the SRS reference at 1. The circles mark the
original guided design and Saturation Quadrants. The two rows use different vertical scales.
[PDF](../../results/igs_sensitivity_2026-11/full/summary/fig_igs_re_vs_srs.pdf)
[PNG](../../results/igs_sensitivity_2026-11/full/summary/fig_igs_re_vs_srs.png)

### Guided designs, queen contiguity and oracle fit

| Regions | Control only 60/40 | 70/30 | 80/20 | 90/10 | 100/0 | Both arms 60/40 | 70/30 | 80/20 | 90/10 | 100/0 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.899 | 0.676 | 0.555 | 0.551 | 0.876 | 0.977\* | 1.051\* | 1.276 | 1.940 | 8.126 |
| 4 | 1.040\* | 0.837 | **0.667** | 0.588 | 0.594 | 1.005\* | 1.018\* | **1.088** | 1.218 | 1.632 |
| 10 | 1.267 | 1.041\* | 0.887 | 0.767 | 0.678 | 1.037\* | 1.003\* | 1.020\* | 1.077 | 1.225 |
| 25 | 2.674 | 1.514 | 1.273 | 1.048\* | 0.862 | 1.463 | 1.196 | 1.140 | 1.100 | 1.152 |

### What changes from 90/10 to 100/0

Under control-only spillover the direction depended on region count. MSE rose for the
two halves and changed little for quadrants, as expected for the larger regions.
MSE continued to fall with ten or 25 regions in both guidance settings. Those increases
for two halves and decreases for ten and 25 regions occurred in all five incidence
configurations. At 100/0, both 25-block variants had detectably lower MSE than SRS for
the first time in the pooled queen/oracle comparison.
For guided halves, MSE relative to SRS rose from 0.551 to 0.876, and bias in the
estimated intervention effect rose from 0.086 to 0.161 in units of the outcome.

Under both-arms spillover, all eight MSE ratios rose in the pooled queen/oracle
results. Guided halves rose from 1.940 to 8.126 times SRS's MSE. The unguided halves
rose from 1.933 to 8.129. The smallest increase was for the 25-block unguided variant,
from 1.009 to 1.013, with an increase in three of five incidence configurations. The
other seven increases occurred in all five configurations. The comparisons between
90/10 and 100/0 are descriptive.

With one half fully treated, WZ and Z differ only for the 20 clusters along the split.
The two columns have correlation 0.941 under queen weights and 0.972 under rook weights.
This makes intervention and spillover effects harder to distinguish in the both-arms
model even though the fitted matrix has full rank.

| Regions | Control only, unguided 90/10 → 100/0 | Both arms, unguided 90/10 → 100/0 |
|---:|---:|---:|
| 2 | 0.540 → 0.846 | 1.933 → 8.129 |
| 4 | 0.570 → 0.562 | 1.210 → 1.579 |
| 10 | 0.733 → 0.626 | 1.001 → 1.185 |
| 25 | 0.974 → 0.786 | 1.009 → 1.013 |

The original guided design remained a middle case, with ratios of 0.667 under
control-only spillover and 1.088 under both-arms spillover. Two guided halves at 70/30
gave 0.676 and 1.051. At 100/0, quadrants were the only new guided design with a
detectable MSE reduction relative to the original under control-only spillover
(ratio 0.892, Holm p < 0.001).

### Incidence guidance, targeting, and coverage

At 100/0, the guided-to-twin MSE ratios were 1.035 to 1.097 under control-only
spillover. All four differences survived Holm adjustment. Under both-arms spillover,
the 25-block guided design had 1.137 times its twin's MSE (Holm p < 0.001). The other
three guided designs were within 3.4% of their twins. The range through 90/10 had been
0.955 to 1.095 across both regimes, so the boundary extends the upper end of that range.
Every detectable queen/oracle twin difference favored the unguided design.

The guided 100/0 designs treated 53.6% to 60.1% of total incidence, averaged over
configurations in the queen allocation block used for the burden check. Their twins
treated about 50%. The 25-block design reached 60.1%, compared with 57.7% at 90/10.
The original design treated 52.7%. [Burden-share table](../../results/igs_sensitivity_2026-11/full/summary/igs_burden_share.csv)

The rook oracle results showed the same directions from 90/10 to 100/0. The two halves
had slightly higher MSE than SRS at 100/0 under rook control-only spillover. With the
non-oracle fit, MSE fell for the four-region variants in that regime with both contiguity
types. Under both-arms spillover, all eight MSE ratios rose in each of the four pooled
estimator-by-contiguity comparisons.

Oracle t + df coverage at 100/0 ranged from 0.932 to 0.953 under queen contiguity
and from 0.937 to 0.951 under rook contiguity. Coverage near 0.932 occurred for
the two halves under control-only spillover. With rook weights, both-arms MSE also
rose sharply for the 100/0 halves, to about 17.3 times SRS's MSE. The non-oracle
fit also showed the large both-arms penalty for the halves, where coverage dropped
as low as 0.607 under rook contiguity.

## Holm adjustment for the original variants

The original MSE results are identical. The original unadjusted p-values and the
Holm adjustment within the original families agree with November 24 at the precision
of the archived CSV in Git revision `a7527fd`. Three of 756 original comparisons crossed 0.05, all in rook
sensitivity slices. The original queen/oracle decisions stayed the same.

| Fit | Spillover slice | Original comparison | Holm p before | Holm p after |
|---|---|---|---:|---:|
| Oracle | Control only | Unguided halves at 60/40 vs SRS (116 vs 9) | 0.04124 | 0.05499 |
| Non-oracle | Control only | Guided quadrants at 90/10 vs twin (107 vs 122) | 0.04677 | 0.05067 |
| Non-oracle | Regimes pooled | Guided 10 regions at 90/10 vs twin (111 vs 126) | 0.04398 | 0.05131 |

The [Holm comparison CSV](../../results/igs_sensitivity_2026-11/full/summary/igs_holm_comparison.csv)
contains both adjusted p-values for every original comparison and marks threshold
crossings. The 756 rows include the Design 3 and 8 tests wherever those designs were
part of the original families, with one row per twin comparison. The
[full contrasts](../../results/igs_sensitivity_2026-11/full/summary/igs_contrasts.csv)
include the new variants and all reported slices.

## Limits

- The study covers equal-size partitions of the grid with targets spaced linearly in rank.
  The North Carolina map and repeated target-level variants were outside the study.
- Region count changes region size, shape, and the variation in target saturations.
  At a fixed s, the population standard deviation of the target levels is
  s√((R + 1)/(3(R − 1))) before rounding. The comparisons apply to each combination
  of partition and saturation. The two-region partition uses one fixed orientation.
- At 100/0 the extreme regions have a single treatment arm. Full matrix rank establishes
  numerical identification for the simulated assignments. Accuracy and coverage can
  still deteriorate, as they did for the halves.
- Results pool over the five effect sizes and five incidence configurations. The paired
  tests use ten surfaces per configuration, so detectability depends on the variability
  across these surfaces. The null-control twin pair had an unadjusted p-value of 0.031
  under queen/oracle control-only spillover on November 24. Its two designs have the
  same assignment distribution.
