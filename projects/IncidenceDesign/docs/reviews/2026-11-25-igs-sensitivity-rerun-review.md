# Independent review of the Incidence-Guided Saturation boundary rerun

Author: Andrew Walther

Date: 2026-11-25

From 90/10 to 100/0, mean squared error (MSE) moved in different directions by region count under control-only spillover. Under both-arms spillover, MSE rose for every region count and both guidance settings in the summaries pooled over incidence configurations within that regime. The November 24 results and September production benchmarks reproduce exactly. Three earlier comparisons lose Holm significance after the new variants enter their families. All three use rook contiguity.

## Variant IDs and randomization

The original variant table matches its November 24 version exactly. IDs 131–134 are guided and IDs 135–138 are unguided, in order of 2, 4, 10 and 25 regions within each group. Saturation Quadrants and the original Incidence-Guided Saturation design retain production IDs 3 and 8 at four regions and 80/20 target extrema. Simple Random Sampling (SRS) retains ID 9.

The assignment seed key includes the design ID and surface as well as the incidence configuration, contiguity, outcome spatial dependence, spillover magnitude and regime. I enumerated all 131,200 assignment key strings for the 41 designs. Both the strings and their integer seeds were unique. The production runner gives every design in a block the same incidence surfaces and outcome noise, then uses each design's own assignments to generate the outcomes. The five effect sizes and both estimators reuse those assignments.

## Treatment counts and identification

Largest-remainder rounding yields exactly 50 treated clusters per 100-cluster assignment. At 100/0, two regions yield counts 50 and 0. The four-region targets of 100%, 66⅔%, 33⅓% and 0% yield realized counts of 25, 17, 8 and 0 per quadrant. The 100/0 axis label describes the target extrema. With 25 small regions, some allocations already have fully treated and fully control regions at 90/10.

The smoke table counts 960 unique boundary draws and 3,840 valid model fits. The full checkpoint diagnostics count 640,000 unique boundary draws, or 80,000 per variant and 320,000 in the control-only regime. Each draw receives the rank check for [1, Z, WZ, X]. The control-only checks also use [1, Z, (1 − Z)WZ, X]. Both runs report zero rank deficiencies, failed fits, aliases and other warnings. The identification CSV repeats the same draws and rank results on the two estimator rows, so the draw count uses one row per contiguity, regime and design.

I regenerated 8,000 of the full-run assignments across both weights matrices, both spillover regimes and all eight boundary IDs. Every assignment was binary, treated exactly 50 clusters, and had at least one fully treated and one fully control region. Both rank checks returned four wherever applicable. Across the full combined scenario file, every row had 250 valid fits and a mean treated count of 50. Every surface row had 25 valid fits.

## Manifest, reproduction and totals

The current code and environment reproduced both saved manifest objects exactly. All 40 November 24 unit checkpoints and all 40 boundary unit checkpoints are present. I compared the retained 105,600 scenario rows and 1,056,000 surface rows directly with the November 24 Git blobs using exact R equality, including the empty warning data frame. I also keyed and sorted the saved September production rows for Designs 3, 8 and 9. All 9,600 rows across the 26 columns recorded in the production files matched exactly.

The initial 10-worker run logged units 1–20 as complete before halting. The authorized four-worker resume logged units 21–40 as complete. A diff of the temporary script against the repository runner showed the worker-count override and corresponding assertion as its only changes. The resume's 8.2-minute timer covers the resumed pass, including reading prior checkpoints and fitting the remaining units. The log provides no explicit cause for the initial halt.

| Scope | Scenario rows | Surface rows | Model fits across both estimators |
|---|---:|---:|---:|
| 33 retained designs | 105,600 | 1,056,000 | 26,400,000 |
| 8 boundary variants | 25,600 | 256,000 | 6,400,000 |
| 41 designs combined | 131,200 | 1,312,000 | 32,800,000 |

The scenario total per design is 5 incidence configurations × 2 contiguity types × 4 outcome spatial-dependence settings × 4 spillover magnitudes × 2 regimes × 5 effect sizes × 2 estimators = 3,200. I verified the per-design row counts and unique scenario and surface keys. Each scenario row represents 250 model fits. This review was limited to the stored fit summaries, checkpoint diagnostics and independently regenerated assignments.

## Paired contrasts and Holm adjustment

I recalculated the mean MSE and paired differences over the 50 configuration-by-surface units in each of 12 slices. The stored Wilcoxon p-values, MSE ratios and paired t₄₉ interval endpoints matched all 948 recalculations within numerical tolerance. The paired interval rescales the interval for the mean MSE difference by the observed comparator MSE. The p-values come from a separate signed-rank test.

Within each estimator, contiguity and regime, the Holm families increased from 32 to 40 against SRS, 15 to 19 against the original guided design, and 16 to 20 guided-versus-twin pairs. The same family sizes apply to the separately adjusted summaries pooled over regimes. For all 756 retained comparisons, the raw p-values and the Holm p-values from the smaller families matched the November 24 Git archive at stored CSV precision. A fresh application of Holm to each family matched both the before and after values. Three adjusted p-values moved from below 0.05 to above 0.05.

| Estimator | Regime | Rook comparison | Holm before | Holm after |
|---|---|---|---:|---:|
| Oracle | Control only | Unguided two-region 60/40 vs SRS (116 vs 9) | 0.041240 | 0.054987 |
| Non-oracle | Control only | Guided vs unguided four-region 90/10 (107 vs 122) | 0.046772 | 0.050669 |
| Non-oracle | Pooled over regimes | Guided vs unguided ten-region 90/10 (111 vs 126) | 0.043981 | 0.051311 |

The 25-region 60/40 null-control designs both assign two of four clusters to treatment per block. Separate assignment seeds give each version its own realization. For the queen oracle control-only slice, the guided-to-twin ratio was 1.036, the raw p-value was 0.0313 and the adjusted p-value was 0.350. The null pair had no comparison below 0.05 after Holm in the 12 slices.

The summary script checks the historical contrasts directly on the first boundary run. On later runs, it compares the raw and original-family adjusted p-values with columns carried forward in the boundary comparison CSV. Keep the November 24 Git revision as the independent reference when repeating the historical check.

## The 100/0 boundary

The table gives MSE relative to SRS in the queen oracle analysis, pooled over the five incidence configurations and five effect sizes within each regime. Each cell follows the same design from 90/10 to 100/0. The comparisons across the two spreads are descriptive. The specified tests compare designs with SRS, the original design or their twin.

| Regions | Control only, guided | Control only, unguided | Both arms, guided | Both arms, unguided |
|---:|---:|---:|---:|---:|
| 2 | 0.551 → 0.876 | 0.540 → 0.846 | 1.940 → 8.126 | 1.933 → 8.129 |
| 4 | 0.588 → 0.594 | 0.570 → 0.562 | 1.218 → 1.632 | 1.210 → 1.579 |
| 10 | 0.767 → 0.678 | 0.733 → 0.626 | 1.077 → 1.225 | 1.001 → 1.185 |
| 25 | 1.048 → 0.862 | 0.974 → 0.786 | 1.100 → 1.152 | 1.009 → 1.013 |

- Under control-only spillover, MSE for two halves rose 59% with guidance and 57% without it. Four quadrants were nearly flat, with a 1.1% increase guided and a 1.4% decrease unguided. MSE fell 12–15% with ten regions and 18–19% with 25. The original expectation of a plateau or reversal applies to the two- and four-region oracle results. The rook oracle results show the same directions, and its two-region boundary MSE was slightly above SRS for both guidance settings.
- Under both-arms spillover, pooled MSE rose from 90/10 to 100/0 across every region count and guidance setting in all four estimator-by-contiguity slices. For the queen oracle fit it was roughly fourfold in two halves and 0.4% for 25 unguided blocks. In those 25-block queen oracle results, MSE rose within three incidence configurations and fell within the two spatial-uniform configurations. For the rook oracle fit, the two-half boundary MSE exceeded 17 times the SRS MSE in both guidance settings. In fully treated and fully control halves, WZ tracks Z away from the dividing line, which explains how the design can have full rank and still yield imprecise estimates.
- At 100/0 under control-only spillover, guided MSE was 3.5–9.7% higher than its unguided twin for each queen oracle region count, and all four Holm-adjusted p-values were below 0.05. Under both-arms spillover, only the 25-region queen oracle twin comparison was below 0.05 after adjustment. The allocation benefit of guidance is distinct from estimation accuracy. At 100/0 the guided designs treated 53.6%, 54.5%, 56.8% and 60.1% of simulated incidence as the number of regions increased. Unguided designs treated about 50% in each case. These shares come from the specified queen burden calculation at outcome spatial dependence zero, spillover magnitude 0.5 and control-only spillover, averaged across five incidence configurations.
- The non-oracle control-only MSE fell from 90/10 to 100/0 at four regions under both contiguity types. Include this exception when describing the oracle plateau. At 100/0 the queen non-oracle t + df coverage was about 68% for 25 regions under control-only spillover and 72–73% for two regions under both-arms spillover.

For the supporting text, qualify the expected reversal by region count, name the spillover regime in each conclusion, and describe 100/0 as a boundary case on the grid. Use the paired tests for their specified comparisons and the table for the direction and size of the change between spreads.

## Figure

The PDF and PNG have the five requested x-axis labels from 60/40 to 100/0 and one open circle at the four-region 80/20 point in each panel. The legend identifies the original design and Saturation Quadrants. The guided column identifies which original design the first phrase refers to. The caption can spell out both names and state that the axis shows target saturation extrema, the row scales differ, and the intervals are paired t₄₉ intervals.

The figure has no title or note. The rebuilt 6.5 × 5.5 inch PDF embeds Arial and uses 8 pt labels and 9 pt axis titles at export size. Extracted PDF text contains the true hyphen in “Incidence-guided.” The plotted labels and values need no negative signs. The 1950 × 1650 pixel PNG shows the same axes, circles and legend.
