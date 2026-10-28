# CTJ v5 statistician review, round 3 (2026-10-27)

Read-only review of `CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`. Placeholders (`\rerun`, `\hl`) were ignored. I re-checked every main-text, abstract and caption number against `manuscript_regime_means.csv`, `manuscript_paired_ratios.csv`, `interval_comparison.csv`, `grid_configuration_means.csv`, `annual_srs_ratios.csv`, `refined_tail_comparison.csv`, `yearly_sensitivity_means.csv`, `sensitivity_setting_matched_comparison.csv`, `grid_tdf_scenarios.rds` and SI S11. **All numbers match.** The values checked were 0.0907/0.1285, the τ = 3 values 0.1596/0.1822, non-oracle bias and t + df coverage of −0.12 to −0.40 and 0.54–0.91, 21.1–29.9%, 20.7–43.4%, fourfold and 2.5-fold (per-setting maximum 2.52), 35–37 colleges, and the 22–24%, 27%, 27 and 29 treated counts. The body is 3,490 words (texcount). No sentence starts with "Because", and the only prose colons are the ones the author approved. Every round-2 finding was taken except NIT 13.

## BLOCKING

None.

## SHOULD

1. **Methods, last sentence: "rejecting τ = 0 when the 95% interval lay above zero".** The power code (`methods_ext_participation.R`, l. 233) rejects when |t|/infl > t₀.₉₇₅, which is a two-sided test. SI S3 also reports a size of 4.4–5.6%, which a one-sided 2.5% test could not produce. *Fix (−1 word):* "when the 95% interval excluded zero".

2. **Methods, application replication: "To avoid overstating worst-decile error by selecting on noise, the worst-decile assignments were chosen with one half of these outcomes and evaluated with the other half."** As written, this makes the split-half estimate the reported quantity. The Results headline (20.7–43.4%) and Figure 3 use the full-sample `Worst10_Mean_Estimated`, and the split-half is reported only as a check ("the pattern held when…"). *Fix (+1):* "To check that selecting on noise did not overstate worst-decile error, we also chose the worst decile with one half of these outcomes and evaluated it with the other."

3. **Table 2 caption: "Its 95% interval is the paired design-minus-SRS interval for the mean difference … divided by the SRS MSE".** Taken literally, this interval would be centred near 0. The code (`19_manuscript_exhibits.R`, l. 40–41) computes 1 + CI/MSE_SRS. *Fix (captions are not counted):* "Its 95% interval is one plus the paired … interval … divided by the SRS MSE".

## NIT

4. **Methods, queen Checkerboard: "the spillover covariate differs from treatment only at the grid boundary".** In the interior the covariate is 0.5γ(1 − Z). It clearly differs from Z there, but it is a linear function of Z. *Fix (+2):* "is a linear function of treatment except at the grid boundary".
5. **Results, NC ¶3: "about 3.74 times its MSE".** "About" with two decimals reads oddly. *Fix:* "about 3.7 times".
6. **Results, NC allocation risk: "the two outcome halves selected largely different worst deciles".** This is not true for Isolation Buffer, whose both-arms Jaccard overlap is 0.54–0.82 (other designs 0.15–0.29). *Fix (+3):* "…largely different worst deciles for most designs". Fund it from 8.
7. **Discussion: "the low both-arms error of High Incidence Focus in the application".** Results never states this value (0.92–0.95 × SRS). It appears only in Figure 4. *Fix (+2):* add "(Figure 4)" after "application".
8. **Results, last ¶: "Alternative incidence summaries reordered some saturation regions, and the prespecified mean-rank rule was kept."** The sentence does not say whether the conclusion survived. It did: Incidence-Guided Saturation's control-only MSE was 0.216–0.244, still below SRS's 0.308. *Fix (+2):* replace ", and the prespecified mean-rank rule was kept" with ", but Incidence-Guided Saturation stayed below SRS under control-only spillover".
9. **Continued work: "The expected effect size and outcome variability also need to be established. They would also need to assess…".** "They" now refers back to "effect size and outcome variability". *Fix (word-neutral):* "Investigators should also assess…".
10. **Methods, application: "Checkerboard became a deterministic graph version".** Table 1 says "Deterministic map version". Use "map" in both places.
11. **Introduction ¶3: "we compared simple random assignment with block stratified sampling".** Two sentences later the same comparator is named "simple random sampling (SRS)". *Fix:* "simple random sampling" (word-neutral).
12. **Abstract and Results coverage (94.2–94.9% and 92.7–93.6%) vs the SI "Interval calibration" paragraph (94.1–94.9% and 92.6–93.6%, and DECISIONS 2026-10-26).** The SI range covers all nine strategies, and the main text covers six plus SRS. Both are correct, but a referee cross-checking will see a mismatch. *Fix (SI only):* add "for the nine strategies" in the SI paragraph.
13. **Abstract: "τ = 1" still has no units** (carried over from round 2). Optional, because the abstract is at 407–411 of 425 words: "τ = 1 (one noise standard deviation)".

Word balance for 1, 2, 4, 6, 7, 8 and 9 is about +9. To fund it, drop "an MCSE of mean MSE at most 5% of its value and an MCSE of coverage at most 0.01" from Results ¶ "Every application setting…" and replace it with "(Methods)" (−17). That cut is optional if the author wants the targets named in both places.

**Verdict:** apart from the `\rerun` placeholders, the draft is ready for the author's sign-off once SHOULD 1–3 (one word, one sentence and one caption phrase) are applied.
