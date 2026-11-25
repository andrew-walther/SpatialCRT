# Phase 2c sign-off: SI correctness and rule fixes

Author: Andrew Walther. Drafted 2026-11-25 from `docs/reviews/2026-11-24-si-review.md`, checked
against the current SI source (`supplementary_info/Supplementary_Information.tex`, after 2a/2b).
"SI n" is a line of that source today. Reply by ID: accept, reject, or edit. Nothing is applied
before sign-off. New numbers come from `results/srs_benchmark/srs_benchmark_summary.txt`
(written by `code/15_srs_benchmark.R`), section 3, "Paired surface-level tests vs SRS (tau = 1)",
which uses the same 50 configuration-by-surface units, averaged over both regimes, as the
adjacent-pair chain.

## Batch 1: decision, rule conflicts, stale statements, SRS in the tests

### D1. Ratio decimals (Decision)

Recommendation: **keep two decimals for ratios to SRS** (tables and figures), three only where
two would hide a difference, as now. Reason: you decided this on 2026-11-15 (DECISIONS.md,
"Ratios to SRS keep their current precision"), main Table 2 uses two, and a ratio is not on the
MSE scale, so the 3-decimal rule's purpose (never printing 0.000) does not arise. To remove the
conflict, add to the AGENTS.md MSE-precision rule: "Ratios to SRS are the exception: two
decimals (author, 2026-11-15)." No SI or CTJ number changes under this option. The alternative
(three decimals) would change main Table 2, SI Tables S19 and S25–S28, and Figure S16.

### R1. SI 1774, heading of 12.1 (rule: no single recommended design)

Current: `Robustness across blocks and the recommended design`
Proposed: `Robustness across blocks`
(The internal `\label` stays, so no reference breaks.)

### R2. SI 1720, non-oracle section

Current: "These results do not change the recommendation. Incidence-Guided Saturation also had the lowest non-oracle MSE under both neighbor types…"
Proposed: "These results do not change the comparison of designs. Incidence-Guided Saturation also had the lowest non-oracle MSE under both neighbor types…"

### R3. SI 1812–1814, end of 12.1

Current: "Among the six designs in the detailed comparison, Incidence-Guided Saturation combines low pooled MSE with comparatively favorable control-only performance. The SRS comparison in Section 12.2 qualifies this. A saturation design is a reasonable choice when control-only spillover is plausible, and SRS remains competitive under both-arms spillover."
Proposed: "No single design was best under both spillover regimes. Under control-only spillover, Incidence-Guided Saturation had the lowest MSE of the six designs, and it had the lowest worst-case block MSE (Table S15). Under both-arms spillover, SRS and Balanced Quartiles were slightly more accurate (Section 12.2). A saturation design is a reasonable choice when control-only spillover is plausible, and SRS remains competitive under both-arms spillover."
Check: both-arms paired differences from SRS are +0.0035 (Incidence-Guided Saturation) and −0.0002 (Balanced Quartiles, not significant), section 3b of the summary.

### R4. SI 1912, case-study accuracy

Current: "High Incidence Focus is year-sensitive under control-only spillover and cannot be recommended uniformly."
Proposed: "Under control-only spillover, the accuracy of High Incidence Focus depended on the planning year."

### S1. SI 233–236, the two variants (stale: main Table 2 now has six designs and SRS)

Current: "Both variants remain in the main-text nine-design benchmark and in the case study, and Section 7.4 reports the comparison of all eight designs at the primary scenario, β = 1."
Proposed: "Both variants remain in the case study, and Section 7.4 reports the comparison of all eight designs and SRS at the primary scenario, β = 1."

### S2. SI 1315–1317, Section 7.4 (stale: Table S12 prints three decimals)

Current: "The MSEs and coverages of the six retained designs are those of Table S7, given here to four and three decimal places."
Proposed: "The MSEs and coverages of the six retained designs are those of Table S7."

### S3. Table S3 caption (stale: the benchmark has nine strategies with SRS)

Current: "Parameter grid of the original eight-design simulation. The crossed parameters form 5 × 5 × 2 × 4 × 4 × 2 × 8 = 12,800 scenarios, 9,600 of them for the six designs of the detailed comparison; the settings in the lower panel are the same in every scenario."
Proposed: "Parameter grid of the grid simulation. The crossed parameters form 5 × 5 × 2 × 4 × 4 × 2 × 8 = 12,800 scenarios for the eight designs, 9,600 of them for the six designs of the detailed comparison. SRS, run afterward on the same grid, adds 1,600 scenarios (14,400 in all). The settings in the lower panel are the same in every scenario."

### T1. SI 451, methods for the 50-unit tests (rule: SRS in every design test)

Current: "These tests compared each pair of designs adjacent in the mean-MSE order, and each design with SRS within each spillover regime, with Holm's adjustment over the eight designs compared with SRS."
Proposed: "These tests compared each pair adjacent in the mean-MSE order, with SRS in its place in that order, and each design with SRS within each spillover regime, with Holm's adjustment over the eight designs compared with SRS."

### T2. SI 1397–1410, queen adjacent-pair chain (SRS falls between Incidence-Guided Saturation, 0.091, and Balanced Quartiles, 0.131)

Current: "We also compared each pair of designs adjacent in the mean-MSE order over the 50 independent configuration-by-surface units. Every adjacent pair differed significantly. The mean paired differences in MSE, with t49 95% intervals, were −0.040 (−0.043, −0.037) for Incidence-Guided Saturation against Balanced Quartiles, −0.029 (−0.033, −0.026) for Balanced Quartiles against Isolation Buffer, […], with Wilcoxon p ≤ 4.9×10⁻⁹ in every case. The better design of each pair had the lower MSE in all 50 units, except High Incidence Focus against 2×2 Blocking, where it did so in 48 (96%)."
Proposed: "We also compared each pair adjacent in the mean-MSE order, with SRS in its place, over the 50 independent configuration-by-surface units. Every adjacent pair differed significantly except SRS and Balanced Quartiles. The mean paired differences in MSE, with t49 95% intervals, were −0.038 (−0.040, −0.035) for Incidence-Guided Saturation against SRS, −0.002 (−0.005, 0.001) for SRS against Balanced Quartiles (Wilcoxon p = 0.15; SRS lower in 29 of the 50 units), −0.029 (−0.033, −0.026) for Balanced Quartiles against Isolation Buffer, […unchanged…], with Wilcoxon p ≤ 4.9×10⁻⁹ in every other case. The better design of each of those pairs had the lower MSE in all 50 units, except High Incidence Focus against 2×2 Blocking, where it did so in 48 (96%)."
Note: the Incidence-Guided Saturation vs Balanced Quartiles difference (−0.040) drops out because the two are no longer adjacent; the rest of the paragraph is unchanged.

### T3. SI 1590–1601, rook adjacent pairs (rook order: Incidence-Guided Saturation 0.072, Balanced Quartiles 0.084, SRS 0.085, 2×2 Blocking 0.128)

Current: "…In the surface-level comparison, Incidence-Guided Saturation had a lower MSE than Balanced Quartiles in 48 of the 50 units (mean difference −0.012; −0.014, −0.010; Wilcoxon p = 9.9×10⁻¹⁰). 2×2 Blocking and Isolation Buffer were close: […] Every other adjacent pair differed in all 50 units."
Proposed: insert after the Wilcoxon p = 9.9×10⁻¹⁰ sentence: "Balanced Quartiles and SRS were not separated. Balanced Quartiles had the lower MSE in 32 of the 50 units (mean difference −0.001; −0.003, 0.000; Wilcoxon p = 0.11)." and end with "Every other adjacent pair, including SRS against 2×2 Blocking, differed in all 50 units."
(The colon in "were close:" is left for the prose pass, 2g.)

## Batch 2: "identified" → "identifiable" (FCL-028, M-08)

Only the identifiability sense changes. "SUD deaths were identified from death certificates"
(SI 641) and "worst-decile assignments were identified" (SI 855, 2069, 2121) keep "identified".

| ID | SI line(s) | Current | Proposed |
|---|---|---|---|
| I1 | 322 | "treatment and spillover cannot be separately identified;" | "treatment and spillover are not separately identifiable;" |
| I2 | 324 | "treatment is identified but control-only spillover is weakly separated" | "treatment is identifiable but control-only spillover is weakly separated" |
| I3 | 341 | "so β is not identified: the fitting routine drops" | "so β is not identifiable: the fitting routine drops" |
| I4 | 352–353 | "β is identified, but imprecisely." | "β is identifiable but imprecisely estimated." |
| I5 | Table S7 caption (1035) | "Under rook contiguity Checkerboard is not identified (WZ = 1 − Z)" | "Under rook contiguity, β is not identifiable for Checkerboard (WZ = 1 − Z)" |
| I6 | captions at 1116, 1196, 1254, 1324, 1481, 1631, 2293 | "β is not identified for Checkerboard" | "β is not identifiable for Checkerboard" |
| I7 | captions at 1660, 1675, 1689 | "(β not identified)" | "(β not identifiable)" |
| I8 | Table S14 caption (1735) | "β is not identified for Checkerboard (WZ = 1 − Z)" | "β is not identifiable for Checkerboard (WZ = 1 − Z)" |
| I9 | prose at 1426, 1531, 1542, 1603, 1652 | "β is not identified" (five places) | "β is not identifiable" |
| I10 | 1922 | "Checkerboard's treatment effect is identified on this irregular graph […] This differs from the exact non-identification under rook contiguity" | "Checkerboard's treatment effect is identifiable on this irregular graph […] This differs from the exact non-identifiability under rook contiguity" |
| I11 | 2405 | "Isolation Buffer is never identified, […] no strategy was identified in every fit. The saturation designs were identified most often there (in 89.7% and 91.7% of fits) and, among identified fits, had MSE…" | "Isolation Buffer is never identifiable, […] no strategy was identifiable in every fit. The saturation designs were identifiable most often there (in 89.7% and 91.7% of fits) and, among identifiable fits, had MSE…" |
| I12 | Table S28 caption (2410) | "n.i.: not identified in every fit" | "n.i.: not identifiable in every fit" |

Chapter 3: same change needed wherever the appendix says "identified" in this sense (2i list).

## Batch 3: figure files (through the producing scripts; no hand edits)

Root cause of the minus signs: the scripts write ASCII hyphens, but R's default `pdf()` device
prints them with the minus glyph. String edits alone cannot fix it; saving through
`save_thesis_pdf()` (`code/00_thesis_figure_device.R`, cairo via `pdftocairo`) does, as the
thesis-chapter copies already show.

| ID | SI figure | Problem (confirmed by text extraction) | Proposed fix | Producer |
|---|---|---|---|---|
| F1 | S18 `fig:nc-setting-ratios` | Row and facet labels "rho=0.5, gamma=0.8" | Labels "ρ = 0.5, ψ = 0.8" (plotmath) | `application/code/real_sud_summary.R:118` |
| F2 | S19 `fig:nc-mse-detail` | Axis "Mean MSE across six rho/gamma settings (log scale)" | "Mean MSE across six ρ/ψ settings (log scale)" | `real_sud_summary.R:124` |
| F3 | S22 `fig:nc-risk` | Legend "Estimated q90", "Estimated worst decile" | "Estimated 90th percentile", "Estimated worst-decile mean" | `real_sud_summary.R:159, 162` |
| F4 | S5 `fig:nc-examples` | Panel titles "1. Checkerboard (n/58)" … "9. Simple Random Sampling" | Names only: "Checkerboard (n/58)" | `application/code/render_real_sud_companion.R:47` (also feeds the offline companion HTML) |
| F5 | S10 `fig:tau-mse` | In-figure axis note "(band: ± average scenario−level Monte Carlo SE)"; the caption already explains the band | Axis "MSE"; save via `save_thesis_pdf()` | `code/12_six_design_statistical_comparisons.R:330` |
| F6 | S11 `fig:tau-coverage` | Axis "Coverage (+/− 1 Monte Carlo SE)"; the caption does not mention the band | Axis "Coverage"; add to the caption "Shaded bands show ±1 times the average scenario-level Monte Carlo standard error of coverage." (sentence moved from SI 1535, so the prose pass can drop it there) | `code/06_visualizations.R:949` via `code/14…R:589` (an argument from 14, so 06's default is untouched) |
| F7 | S7, S13 `fig:coverage`, `fig:coverage-rook` | One-panel facet strip "Queen Contiguity"/"Rook Contiguity" reads as an in-figure title | Drop the strip (caption names the neighbor type) | `code/06_visualizations.R:270–271` via `code/12` (argument from 12) |
| F8 | S17 `fig:nc-pop`, S24 `fig:power-n`, and every figure above | Hyphens printed as minus signs ("working−age", "college−level", "scenario−level") | Save through `save_thesis_pdf()` | `real_sud_summary.R`, `code/23_methods_extension_figures.R`, `code/12`, `code/14` |

Reruns after sign-off: `real_sud_summary.R` (stored results and the restricted `setup.rds`;
rewrites every file in `real_sud_rev_20261002/exhibits/`, numbers unchanged), then
`application/code/manuscript_maps.R`; `code/12`, `code/14` and `code/23` from stored results;
then `bash paper/tools/copy_figures.sh`. Each rerun is checked by diffing the CSVs it rewrites
against git (they must be byte-identical or differ only in timestamps) and by re-extracting the
figure text. `code/01`–`05` are not touched. Chapter 3 figures come from `code/28` and already
print true hyphens; F1–F4 labels also appear in chapter copies of the same files, listed for 2i.

## Decisions

- Batch 1 (author, 2026-11-25): D1, R1–R4, S1–S3, T1–T3 accepted as proposed. D1 recorded in
  DECISIONS.md and the AGENTS.md MSE-precision rule (Dr. Lin: no more than three decimals in the
  main manuscript). Applied with `apply_2c.py R1 R2 R3 R4 S1 S2 S3 T1 T2 T3`.
- Batch 2 (author, 2026-11-25): I1–I12 accepted as proposed. Applied with `apply_2c.py I1 I2 I3 I4
  I5 I6+I8 I7 I9 I10 I11 I12`; the four remaining "identified" are the non-identifiability uses
  (deaths identified from certificates; worst-decile assignments identified).
