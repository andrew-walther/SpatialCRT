# CTJ v5 statistician review, round 2 (2026-10-27)

Reviewer stance: independent Clinical Trials biostatistician. Read-only. Placeholders (`\rerun`, `\hl`) ignored.
Every number in the abstract, Results and captions was checked against `manuscript_regime_means.csv`, `manuscript_paired_ratios.csv`, `interval_comparison.csv`, `grid_configuration_means.csv`, `annual_srs_ratios.csv`, `refined_tail_comparison.csv`, `grid_tdf_scenarios.rds` (non-oracle and τ sweep) and SI Table `tab:nc-guidance`. **All checked numbers match.** No sentence starts with "Because". There is no "Chapter 2" and no DIM.

## BLOCKING

None.

## SHOULD

1. **Results, NC allocation risk: "most of the variance in estimated risk came from Monte Carlo error (SI)".** This is true for Incidence-Guided Saturation (61–86%) but not as a general statement. For SRS it is 45%, for 2×2 Blocking 35–43% and for Isolation Buffer 7%, and the SI restricts the claim to Incidence-Guided Saturation. *Fix (+3 words):* "most of the variance in Incidence-Guided Saturation's estimated risk came from…". Fund it from fix 9.

2. **Methods, allocation risk: "the allocation risk $m(a)=E_\epsilon[\cdot]$" vs abstract "error under an unfavorable realized assignment (allocation risk)" and Figure 3 "assignment-specific MSE".** One term names both the per-assignment MSE and its downside summary. The same sentence says the worst-decile mean summarizes the "spread" of $m(a)$, but it summarizes the upper tail. $\hat\tau$ is never formally introduced. *Fix (word-neutral):* write "the assignment-specific MSE $m(a)=\ldots$". Then: "we summarized allocation risk, the upper tail of $m(a)$, by the mean of the worst 10%…". Change "We estimated $\tau$ by" to "We estimated $\tau$, giving $\hat\tau$, by" (+2) and drop "in validation" (−2).

3. **Table 2 caption: "Ratio to SRS (95% CI) … its 95% interval comes from the paired design-minus-SRS differences".** The code computes $1+\mathrm{CI}(\bar d)/\mathrm{MSE}_{SRS}$, which treats the SRS MSE as fixed. That is a valid test of ratio = 1, but it is not a CI for the ratio. A referee will ask. *Fix (+4):* "its 95% interval is the paired difference interval divided by the SRS MSE, so an interval excluding 1…".

4. **Methods, last sentence, and Results: "35 to 37 colleges", "fell from 22% with 58 colleges to 5% with 29".** The power study used settings that the main text never states: spillover $\gamma=0.8\tau$, $\rho\in\{0,0.5\}$, random volunteering, and exposure counted over all geographic neighbors (SI S11). Its 22% at 58 colleges therefore differs from the primary 21.1–29.9%, and the paper doesn't say why. The test's direction is also missing: the code rejects when $t>t_{0.975}$, which is a one-sided 2.5% test. *Fix (+14):* "…to estimate the number of colleges needed for 80% power, with spillover $\gamma=0.8\tau$ and $\rho\in\{0,0.5\}$, rejecting $\tau=0$ when the 95% interval lay above zero (SI Section S11)." Fund it from fix 10.

5. **Abstract: "Incidence-Guided Saturation" is never described.** A reader of the abstract alone cannot tell what the recommended design does. *Fix (+11, the abstract is about 400 of 425 words):* "Incidence-Guided Saturation, which treats larger shares of clusters in higher-incidence regions, had…".

6. **Results: "The control-only advantage of Incidence-Guided Saturation narrowed at larger effects."** In a linear model MSE does not depend on $\tau$, so referees will ask why it does here. Isolation Buffer stays at 0.1595–0.1596 while Incidence-Guided Saturation rises from 0.1311 to 0.1822. *Fix (word-neutral):* replace "(MSE 0.1596 and 0.1604)" with "because estimation error grows with $\tau$ when treated clusters have treated neighbors, which Isolation Buffer never has (SI)". Keep the τ = 3 numbers.

7. **Discussion: "Isolation Buffer's protection against poor control-only assignments is useful only when…".** The main text never reports Isolation Buffer's worst-decile result. It appears only in Figure 3 and the SI. *Fix (+3):* add "(Figure 3, SI)".

8. **Methods, application: "where 14 groups of nearby clusters were formed each year … so a one-cluster group was never treated."** This sentence runs about 70 words with four clauses. "Rounded to the nearest even number" is R's round-half-to-even, which non-R readers will find odd. *Fix (word-neutral):* split it after "log-population". Then: "A random half of each group was treated, and an odd half (such as 1.5) was rounded to the nearest even whole number, so a one-cluster group was never treated."

9. **Limitations: "The $t$-based interval kept coverage at 94–95% on the grid and 94–96% in the application…".** This sentence is a positive result, not a limitation, and it repeats Results. *Fix:* cut it (about −27 words), or replace it with "Coverage was checked only under the assumed SAR model." (−17 words).

10. **Continued work: "Under control-only spillover no design reached 80% power at $\tau=1$ with all 58 colleges, so the expected effect size…".** This repeats Results. *Fix (−17):* "The expected effect size and outcome variability also need to be established."

## NIT

11. Methods: "Incidence followed five configurations: independent…" is a colon in prose. Change it to "The five incidence configurations were independent…" (+1).
12. Table 1 has semicolons in cells ("ties are broken at random"; "until none remain; no two treated…"; "(blue, treated; white, control)") and uses "Deterministic graph version", which is jargon. Use periods or commas, and say "Deterministic map version".
13. Abstract: "queen contiguity" is undefined, and "$\tau=1$" has no units. Write "queen contiguity (shared edge or corner)" and "$\tau=1$ (one noise standard deviation)" (+8).
14. Results ¶1: "the five incidence configurations, which changed MSE by less than 10%". The antecedent of "which" is ambiguous. Use "…configurations. Configuration changed each design's MSE by less than 10%, except…" (word-neutral).
15. Methods, grid allocation risk: it does not say whether the worst decile is taken per surface and setting and then averaged, or pooled. Settle this when the rerun text is finalized.
16. Figure 3 caption (about 115 words) contains fragments ("Simulation study, all 50 grid surfaces at…"). Merge them into one sentence: "Grid, all 50 surfaces at…; application, …". The SI caption of the same figure still says "focused grid pilot … 40 … settings" and should be updated with the rerun.
17. Title page: "Main-text word count: 3489. Abstract word count: 411". Recount both after edits (current about 3,495).
18. Introduction ¶4: "These strategies express different priorities, and they need not provide the same information…" is abstract and slightly stuffy. Change to "These priorities can conflict with separating treatment from neighboring exposure." (−8).
19. Results, NC: "Rook and queen contiguity at matched parameter values led to the same broad conclusions" is vague. Name one: "the same regime-specific ranking; rook raised MSE by 11–29%" (+6, optional).

Word balance if fixes 1–8 and 9–10 (cut versions) are all taken: roughly −10 words. That leaves room for 11, 13 and 19.

## Author's v5 concerns

| Concern | Status | Note |
|---|---|---|
| Abstract spillover framing | RESOLVED | |
| Abstract two criteria (MSE, allocation risk) | RESOLVED | The Methods term is inconsistent (SHOULD 2) |
| Incidence allocates, not outcome (β = 0/1) | RESOLVED | |
| "Simulated continuous education outcomes" | RESOLVED | |
| Last abstract Methods sentence | RESOLVED | |
| τ = 1 | RESOLVED | Units missing (NIT 13) |
| "Precision targets … no failures" in abstract | RESOLVED | Removed |
| Conclusions: minimize error, budget/population, no "interval calibration" | RESOLVED | |
| Intro colons | RESOLVED | |
| "Larger study region" | RESOLVED | |
| Explain alternating allocation | RESOLVED | |
| Interspersion → identifiability | RESOLVED | |
| Stuffy language, shorter application intro | PARTLY | NIT 18 |
| Model and estimation before simulation | RESOLVED | |
| Cite queen/rook; primary vs sensitivity | RESOLVED | |
| Why β = 1 on grid | RESOLVED | |
| Cite 35/100,000 | RESOLVED | |
| Ten independent surfaces | RESOLVED | |
| Explain/cite Monte Carlo error | RESOLVED | |
| Choppy sentences, colons, ambiguity | PARTLY | NIT 11, 14 |
| Define every allocation-risk term | PARTLY | SHOULD 2 |
| "Pilot" language | RESOLVED | Only the real "education pilot" remains |
| No "Because" openers | RESOLVED | |
| Application section plain and clear | PARTLY | SHOULD 8, NIT 12 |
| Clarify replication "extension" | RESOLVED | |
| Results: define terms, name precision targets | RESOLVED | |
| What MSE ranges span | RESOLVED | |
| Drop run totals | RESOLVED | |
| "Broad regime distinction" / transfer sentence | RESOLVED | |
| No "grid pilot" | RESOLVED | Main text; the SI caption still has it (NIT 16) |
| "Budgets differed (…)" fragment | RESOLVED | |
| Key results only | PARTLY | Power and τ-trend detail are still dense; SHOULD 4, 6 |
| Plain-language Limitations/Conclusions | PARTLY | SHOULD 9 |
| Continued work: endpoint, staffing | RESOLVED | |
| Figure 2 vs 3 distinguished; Figure 3 caption length | PARTLY | Distinguished; caption still long (NIT 16) |
| Figure 2 shows CIs | RESOLVED | Per caption; graphic not inspected |
| Table 2: SRS 50 of 100 | RESOLVED | |
| Overall narrative flow, reproducibility | PARTLY | Power settings missing (SHOULD 4) |
