# CTJ v5 clinician/epidemiologist review, round 3 (2026-10-27)

Reader: an epidemiologist or trialist who reads *Clinical Trials* but is not a spatial statistician. Source: `CTJ_Reading.pdf` (9 pages), quoted from `CTJ_Manuscript.tex` and `CTJ_Exhibits.tex`. I left alone the highlighted placeholders, the kept colons, "Spatial blocking", the unstated pilot regime and the reason for using service areas.

**Verdict:** Ready for the author's sign-off once S1–S5 are fixed. Nothing blocks. The motivation, comparisons and findings now read clearly to a non-spatial reader.

## BLOCKING

None.

## SHOULD

**S1. Continued work, "The expected effect size and outcome variability also need to be established. They would also need to assess whether geographic contiguity..."** The key planning result never reaches the Discussion. Under control-only spillover, no design reached 80% power at τ = 1 even with all 58 colleges, and that is the regime where Incidence-Guided Saturation helps. A trialist who reads only the Discussion will miss it. Two smaller problems: "They" seems to point back to "effect size and outcome variability", and "also" appears twice. Rewrite: "Under control-only spillover even all 58 colleges gave less than 80% power at τ = 1, so the expected effect size and outcome variability must be established first. Investigators should also assess whether geographic contiguity reflects how professionals actually communicate." The rewrite adds about 13 words. Pay for it by cutting "A practical choice should also compare recruitment and exposure consequences alongside MSE." (12 words, which the next sentence repeats) and "These shares describe geographic reach, not the enrollment or cost of healthcare professionals." (13 words, which repeat Methods line 196).

**S2. Figure 3 caption.** The caption does not say how the settings are combined. A reader sees one row per design but is told about 50 units × ρ × γ, or "each year 2018–2021". According to `code/21_downside_by_design.R`, the percentiles and means are computed within each setting and then averaged. Without that sentence, "a design with one possible assignment per setting appears as a circle only" is confusing, because the grid's single assignment differs across surfaces. Add: "Percentiles and means are computed within each setting and then averaged across settings." Pay for it by cutting "SRS is shown in orange."

**S3. Application, "an odd half (such as 1.5) was rounded to the nearest even whole number".** 1.5 is not odd, so the reader stops here. Rewrite: "with fractional halves rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), so a one-cluster group was never treated." The length is about the same.

**S4. Results 1, "a small but detectable disadvantage. No design was detectably more accurate than SRS (paired interval excluding 1)".** "Detectable" is used before it is defined. Move the parenthesis to the first use: "...a small but detectable disadvantage (paired interval excluding 1). No design was detectably more accurate than SRS under both-arms spillover". This costs no words.

**S5. Abstract Results, last sentence.** Two consecutive sentences open with "In the (North Carolina) application", and "t-based ... degrees-of-freedom-corrected standard error" is heavy for this audience. Rewrite: "Across application designs, 95% t-based confidence intervals with a small-sample correction covered the true effect 94.2–94.9% of the time, against 92.7–93.6% for standard Wald intervals." This is shorter and keeps the numbers and the claim.

## NIT

- **N1. Methods, power sentence, "random subsets of 15 to 58 participating colleges ... spillover γ = 0.8τ".** "Colleges" appears here for the first time, after the text has said "clusters" throughout. Also, "0.8τ" differs from γ = 0.8 everywhere else. Write "15 to 58 of the 58 colleges (clusters)" and "spillover magnitude γ = 0.8τ, proportional to the effect".
- **N2. Spurious precision.** "about 3.74 times its MSE" should be "3.74 times", because Figure 4 shows 3.74 in every year. "only about 22–24%" should be "only 22–24%".
- **N3. Results 1, omitted-spillover paragraph, "coverage 0.54 to 0.91".** Every other coverage value is a percentage ("94–95%"). Write "54–91%".
- **N4. Table 1.** The High Incidence Focus grid cell lacks "(50 treated)", which every other grid cell has. "about 21 treated under queen weights" uses an undefined term, so write "under queen (edge or corner) neighbors".
- **N5. Equation 1.** *I* is never defined. Add "*I* is the identity matrix" to the gloss.
- **N6. Methods, estimation paragraph (about 250 words).** Start a new paragraph at "With rook contiguity, Checkerboard gives *WZ* = 1 − *Z*". The next sentence ("With queen contiguity ... because ... so ...", 37 words) chains three clauses. Split it at "treated." This costs no words.
- **N7. Results, "High Incidence Focus treated 29 clusters but reached only 22–24%".** "Only" needs a reference point. If the SI gives SRS's typical population share, cite it here.
- **N8. Results 1, Isolation Buffer.** Table 2 shows the same Isolation Buffer MSE (0.16) in both regimes but ratios of 3.80 against 0.74. The text explains that its two regimes are identical, but it never says the ratio flips because SRS improves under both-arms spillover. One clause would prevent a misreading.
- **N9. Introduction, "These priorities need not favor separating treatment from neighboring exposure."** This is cryptic. Write "None of these goals ensures that treatment can be separated from neighboring exposure." (+3 words)
- **N10. Abstract word count (411).** Check it against the journal's abstract limit. I am not certain what that limit is.

## Author's v5 concerns

Resolved: the abstract framing, criteria, allocation-versus-outcome distinction, Methods flow, τ = 1 and precision sentence; the Introduction's colons, the region and checkerboard explanation and identifiability; the Methods order, citations, β = 1, the base rate, the surfaces, Monte Carlo error, the risk-equation terms and the replication rule; no "Because" openers; the Results' ranges, run totals, regime and transfer sentences and the budget fragment; the endpoint and staffing; Figure 2 versus Figure 3, the Figure 2 intervals and the SRS "50 of 100" note.

Partly resolved: plain explanations of the rules (S3), defining terms where they are used (S4), the planning implications in the Discussion (S1), and Figure 3 standing alone (S2).
