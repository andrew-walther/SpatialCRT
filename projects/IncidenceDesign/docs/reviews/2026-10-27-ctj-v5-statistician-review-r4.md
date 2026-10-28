# CTJ v5 statistician review, round 4 (2026-10-27)

This was a read-only confirmation pass on `CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`; `\hl` placeholders were ignored. Every number in the abstract, Results and captions matches its source: `manuscript_regime_means.csv`, `manuscript_paired_ratios.csv`, `interval_comparison.csv`, `grid_configuration_means.csv` (all non-HIF spreads are under 10%), `annual_srs_ratios.csv`, `refined_tail_comparison.csv`, `grid_tdf_scenarios.rds` (the τ = 2/3 values), `tab_nc_guidance.tex`/`tab_nc_power.tex` (22% → 5%, 35–37), and `surfaces10/downside_by_design.csv` (0.1912/0.3462 and 0.0583/0.0544). The new surfaces10 wording (Methods ¶ "To estimate allocation risk…" and the Figure 3 caption) is consistent with the findings doc. The IRB sentence is consistent with Ethical considerations and Chapter 3. SI S11 is the participation/power section, as cited. No sentence starts with "Because".

## BLOCKING

None.

## SHOULD

1. **Grid allocation-risk spread is mostly Monte Carlo noise, and the text does not say so (Results ¶ "Figure 3 shows…", Figure 3 caption).** On the grid each assignment has 100 outcomes. The median share of across-assignment variance due to outcome noise is 0.63 for Incidence-Guided Saturation and 0.26 for SRS under control-only spillover, and 0.92–0.93 under both-arms spillover (findings doc, tail-stability table). The grid bars in Figure 3, especially for both arms, therefore show mostly noise. The grid both-arms worst-decile ratio (1.07) simply mirrors the mean ratio (1.07). Unlike the NC analysis, the grid has no split-half check. SI Figure S10 even says noise "widens the bars slightly", which a referee comparing the two would reject. The control-only direction is safe, because noise inflates IGS more than SRS.
   *Fix (word-neutral):* replace "Which assignments fell in the worst decile was uncertain, because outcome noise still moved assignments across the decile boundary." (18 words) with "With 100 outcomes per assignment, outcome noise accounted for much of this spread, so worst-decile membership was uncertain." (18 words). Captions are not counted, so also add to the Figure 3 caption: "Grid assignments have 100 simulated outcomes each, so Monte Carlo noise accounts for much of the grid spread, most of it under both-arms spillover." Change "slightly" in SI Figure S10 to match.

## Fact check on an author decision (not re-litigated)

2. **The one-sided power premise does not match the code.** `application/code/methods_ext_participation.R` l. 224–233 computes `tt <- abs(res[,1]/se)` and rejects when `tt/infl > ct`. That is a two-sided |t| test at 5%, not "interval above zero". At τ > 0 the power numbers are practically unaffected (lower-tail rejections are negligible), but the decision rests on a premise the code does not support. If the author wants text and code to agree, "lay above zero" → "excluded zero" saves one word.

## NIT

3. **Methods ¶ "To estimate allocation risk on the grid…": "100 assignments per surface".** Assignments are drawn per surface *and* parameter setting (3,600 blocks, J = 100 each). *Fix:* "per surface" → "per setting" (word-neutral).
4. **Results ¶ power: "fell from 22% with 58 colleges".** The text sits two paragraphs after "21.1–29.9% lower" for the same 58 colleges. The 22% comes from the power run (ρ = 0.5, γ = 0.8τ, Methods), and a reader may read it as inconsistent. *Fix:* "from 22% with 58 colleges" → "from 22% at full participation". This adds 1 word; fund it by changing "an MCSE of mean MSE at most 5% of its value and an MCSE of coverage at most 0.01" in Results ¶1 to "(Methods)", which also removes a duplication.
5. **Carry-over (round-3 NIT 8, not taken): "Alternative incidence summaries reordered some saturation regions, and the prespecified mean-rank rule was kept."** The sentence still does not say whether the conclusion held. It did: IGS control-only MSE was 0.216–0.244 against 0.308 for SRS. *Fix:* replace ", and the prespecified mean-rank rule was kept" (7 words) with ", but Incidence-Guided Saturation still beat SRS under control-only spillover" (9 words). Fund it from the cut in 4.

**Verdict:** the numbers are clean and there is nothing blocking. Apply SHOULD 1 (one sentence plus a caption sentence) and the draft is ready for the author's sign-off.
