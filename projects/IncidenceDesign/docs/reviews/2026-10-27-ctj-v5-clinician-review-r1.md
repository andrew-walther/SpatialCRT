# CTJ v5 reading draft: clinician/epidemiologist review (round 1)

Reviewer lens: trialist who reads *Clinical Trials* but is not a spatial statistician. Source read: `CTJ_Reading.pdf` (9 pp.) with `CTJ_Manuscript.tex`/`CTJ_Exhibits.tex`. Yellow placeholders ignored. No numbers or claims changed.

Overall, the motivation, the two criteria (average error and allocation risk) and the "incidence guides allocation, not the outcome" framing now come through clearly. The main remaining obstacles for a clinical reader are what the two spillover regimes mean in practice, a few terms used as if known ("identified", "Wald", "saturation"), and two places where a reader could draw the wrong planning conclusion.

## BLOCKING

**B1. Discussion, "the low both-arms error of High Incidence Focus came from a fixed assignment".** Table 2 shows High Incidence Focus with the second-highest both-arms error on the grid (ratio 3.82). The low error is from the application only (Figure 4, 0.92–0.95), and the sentence does not say so. A reader checking Table 2 will think the Discussion contradicts the Results.
Rewrite: "...but the low both-arms error of High Incidence Focus in the application came from a fixed assignment..." (+3 words).

**B2. Power result missing from the planning message (Results, "no design reached 80% power at that effect with all 58 colleges").** The Discussion recommends saturation under control-only spillover, yet in that regime no design reaches 80% power with all 58 colleges. A trialist reading the abstract, Discussion and Conclusions would assume the pilot is feasible as framed. This is a placement problem, not a new claim.
Rewrite (Continued work, after "set cluster and participant budgets"): "Under control-only spillover no design reached 80% power for a one-standard-deviation effect with all 58 colleges, so the expected effect size and outcome variance also need to be established." (+28 words. Fund it with the cuts in S10.)

## SHOULD

**S1. Spillover regimes are never given a practical meaning (Methods, "S(Z) = γWZ when exposure reaches clusters in both arms").** The headline finding is "the regime determines the best design", and the Discussion asks investigators to consider "plausible spillover regimes", but nothing tells a clinician how to judge which regime applies to an education intervention. Add one plain sentence in the Methods: "Control-only spillover describes, for example, trained professionals sharing materials with colleagues in untreated neighboring areas, and both-arms spillover describes neighbors in either arm being affected, such as treated areas reinforcing each other." (+38 words. The examples are illustrative, so confirm with Simpson.)

**S2. "Identified/identification" used 13 times, never defined.** For example, "Coverage was 94–95% for every design whose effect was identified". At first use (Methods, "treatment and spillover cannot be separated"), add "...cannot be separated, so τ is not identified (cannot be estimated)". Elsewhere, "identified design" could read "design that can estimate τ".

**S3. "Wald intervals" are undefined (Abstract, Results, Limitations).** In the Methods CI sentence, add "...rather than the usual Wald interval τ̂ ± 1.96 SE". In the abstract, "Wald intervals" could read "standard normal-based intervals".

**S4. Units of τ (Methods, "The direct treatment effect τ is the coefficient on treatment...").** A clinician does not learn until the power paragraph that τ = 1 is one standard deviation of the outcome noise. Add to the ε definition: "ε is independent standard normal noise, so τ is in standard-deviation units". This also lets you drop the "not an average total effect" construction, which follows the "not X but Y" pattern and is heavy for this audience.

**S5. "Saturation" is the key design and is never defined in plain words before the abstract Conclusions or the Introduction use it (Abstract, "saturation designs had the lowest estimation error"; Intro, "set local treatment saturation").** Intro rewrite: "or vary the share of clusters treated in each region (saturation)". Abstract and Conclusions say "saturation designs" in the plural while the main text shows one saturation design. Either say "incidence-guided saturation" or keep the plural and name the SI variant.

**S6. Design names drift from the agreed set.** "spatial blocking" (twice) in the application text, while Table 1 and Figures 3–4 say "2×2 Blocking". The SI variants appear as "saturation with randomly ordered quadrants and balancing within incidence halves", not Saturation Quadrants and Balanced Halves. The main text writes "incidence-guided saturation" in lower case, but "Incidence-Guided Saturation" in the pooled-ranking sentence and the exhibits. Pick one form.

**S7. Table 1 vs Methods contradiction on 2×2 Blocking in NC.** Table 1 says "14 fixed groups formed by k-means", while the Methods say "14 groups of nearby clusters were formed each year". Rewrite in Table 1: "14 groups formed each year by k-means clustering...".

**S8. "The control-only accuracy of High Incidence Focus ranged from 0.70 to 1.40 times that of SRS".** These are MSE ratios, so "accuracy 0.70 times" reads backwards. Rewrite: "The control-only MSE of High Incidence Focus ranged from 0.70 to 1.40 times that of SRS".

**S9. Pivot sentences and clichés (skill table).**
- "Allocation risk in the application told the same story (Figure 3)." Rewrite: "As on the grid, incidence-guided saturation had lower worst-decile error than SRS under control-only spillover (Figure 3): its worst-decile mean was..." and merge with the next sentence.
- "Isolation Buffer behaves differently. No treated cluster has a treated neighbor, so..." Rewrite: "Isolation Buffer never gives a treated cluster a treated neighbor, so only control clusters are ever exposed and the two regimes define the same outcome model; its small..." (use a period if semicolons are to be avoided).
- "The NC map kept the regime distinction but changed several grid results." This is still vague (author note). Rewrite: "On the NC map the best design again depended on the regime, but several grid results changed."

**S10. Results paragraphs that cut against the main message or repeat exhibits (word budget for B2/S1).**
- "Pooled over both regimes, incidence-guided saturation had a mean MSE of 0.0907 against 0.1285 for SRS, but separating the regimes was essential." Leading with pooled numbers and then disowning them is confusing. Start at "Under control-only spillover its MSE was 0.1355..." (−24 words). Pooled values can stay in the SI.
- The pooled-ranking paragraph ("Pooled across both regimes, the MSE ordering ... need not hold within each regime") again leads with the result that matters least. Consider keeping only its τ = 2/τ = 3 sentence, or moving the whole paragraph to the SI.
- The budgets paragraph ("Treatment budgets differed by design and year (SI). Checkerboard treated 27 clusters...") repeats Table 1. Keep only the population-share sentence, which is new information.
- "All 1,248 application settings completed without estimation failures..." restates the Methods targets. One clause suffices: "Every application setting met both precision targets without estimation failures."

**S11. Sensitivity outcomes unstated.** "Alternative incidence summaries reordered some saturation regions, and the prespecified mean-rank rule was kept." This does not say whether conclusions changed, which is the author's abstract concern again. Rewrite: "...reordered some saturation regions without changing [the conclusions / X], and the prespecified mean-rank rule was kept." (Fill from the SI.)

**S12. "although the baseline and the adjustment were changed together".** The "although" hides the point that their effects cannot be separated. Rewrite: "Aligning treatment with incidence leaves little information about τ after adjustment. The baseline and the adjustment were added together, so their separate contributions are not distinguished."

**S13. "To avoid selection bias, the worst-decile assignments were identified with one half..."** To trialists, "selection bias" means participant selection. Rewrite: "To avoid overstating worst-decile error by selecting on noise, the worst-decile...".

**S14. Application setting unexplained.** Nowhere does the paper say why the clusters are community-college service areas. Readers will assume the colleges deliver the education. Add a short clause at first mention, e.g. "58 community-college service-area clusters, the units through which the education would be delivered" (confirm with the authors).

**S15. Figure 3 caption is not self-contained and does not set itself apart from Figure 2 and Table 2.** "Settings are described in the Methods." Rewrite the opening: "Allocation risk. Figure 2 and Table 2 report average MSE, and this figure shows how MSE varies across individual assignments, relative to SRS's mean MSE in the same setting (queen contiguity, τ = 1; grid ρ ∈ {0, 0.5}, γ ∈ {0.5, 0.8}; application γ = 0.8)." This replaces the colon and the pointer sentence.

**S16. Limitations, "These design comparisons complement spatial analysis rather than replace it."** This is vague, and "not X but Y" again. Rewrite: "A good allocation does not remove the need for an analysis that models spillover."

**S17. Simulation incidence list (Methods, "Incidence followed five configurations, independent Uniform(0, 1) incidence, spatial incidence at two...").** Removing the colon makes the list read as an appositive of six items. A colon is clearly best here: "Incidence followed five configurations: independent...". Also say in four words why counts were rank-scaled ("to share one scale").

## NIT

- Methods, "q90": defined but never used again. Delete "by its 90th percentile (q90) and" (−8 words).
- Intro, "Spatial interference challenges the assumption... Rather, adjacent populations...": "Rather" does not follow "challenges". Use "In practice, adjacent populations...".
- Intro, "without comparing numerical results between the two studies" reads defensively. Cut it or move it to the Discussion.
- Intro, "maximal interspersion" and Table 1, "Fixed spatial interspersion". Use "strict alternation".
- Methods, "spatialreg::lagsarlm" is code syntax. Use "the lagsarlm function of the R package spatialreg".
- Methods, "configuration-by-surface units" appears 3 times. Define once: "the 50 incidence surfaces (five configurations × ten surfaces)".
- Methods, "The observed rates were used as fixed inputs without synthetic incidence". Rewrite: "...used as fixed inputs rather than simulated".
- Results, "(0.30761)", "0.07638": five significant figures against four on the grid. Use four.
- Results, "Isolation Buffer 15.89 on average" vs Table 1 "about 16". Use "about 16".
- Results, "including a one-cluster group that the rounding rule always assigned to control": the rule (round down?) is not stated.
- Results, "Checkerboard is identified... Isolation Buffer is close... Spatial blocking is close": present tense among past-tense results. Use past.
- Results, "pooled bias −0.12 to −0.40": say "across designs".
- Discussion, "The completed grid and observed-incidence application": drop "completed".
- Discussion, "incidence-balanced designs": name Balanced Quartiles.
- Conclusions, "favorable worst-case error": the paper measures the worst decile, not the worst case. Use "lower worst-decile error".
- Discussion, "Plausible spillover regimes should therefore be considered": "therefore" does not follow from the High Incidence Focus paragraph. Drop it.
- Heading, "Application: Sudden unexpected death in North Carolina": the author asked to drop colons in this heading. Use "Sudden unexpected death application".

## Exhibits

Figure 1 is helpful and its caption links year shifts to Figure 4 well. Panel B could note that region ranks change by year. Table 1 is the best orientation device in the paper. Figure 2 vs Table 2 is now clear (the Figure 2 caption points to the paired comparison). Figure 3 is not yet clearly set apart (S15). Figure 4's caption stands alone.

## Author concerns (v5_author_notes.md)

| Concern | Status | Evidence |
|---|---|---|
| Abstract spillover framing | RESOLVED | "treatment effects can spill over and contaminate clusters that neighbor treated clusters" |
| Abstract two criteria | RESOLVED | "average estimation error and the error under an unfavorable realized assignment (allocation risk)" |
| Incidence allocates, not outcome | RESOLVED | Abstract "a continuous education outcome that does not depend on incidence"; Methods β = 0/β = 1 sentence |
| "simulated continuous education outcomes" unclear | PARTLY | Still unstated what the outcome is (e.g. a knowledge score) |
| Last Methods sentence flow | RESOLVED | Now two sentences with verbs; queen/rook still undefined in abstract |
| τ = 1 | RESOLVED | "treatment effect τ = 1" |
| Precision-target sentence | RESOLVED | Removed from abstract |
| Conclusions: budget/population, no "interval calibration" | RESOLVED | "number of clusters treated and the share of the population"; "nominal coverage" |
| Intro colons | RESOLVED | Colon gone; "Rather" now awkward (NIT) |
| "larger study region" | RESOLVED | Present |
| Explain alternating allocation | RESOLVED | "block stratified sampling, an alternating allocation that keeps clusters of the same arm apart" |
| Interspersion = identifiability | PARTLY | Checkerboard example given; "identified" never defined (S2) |
| Stuffy language; shorter application intro | PARTLY | Application intro still 5 sentences; "interspersion", "saturation" |
| Model before simulation | RESOLVED | Model → designs → simulation → application |
| Cite queen/rook; primary vs sensitivity | RESOLVED | Ref 10 and explicit primary/sensitivity |
| Why β = 1 on grid | RESOLVED | "the harder case... treatment can then line up with baseline risk" |
| Cite 35/100,000 | RESOLVED | Ref 13 |
| Ten surfaces | PARTLY | Sharing explained; "surface" never glossed as an incidence map |
| Monte Carlo error | RESOLVED | Defined with Morris et al. |
| Choppy/colon/semicolon prose | PARTLY | Pivots remain (S9) |
| Allocation-risk equation terms | RESOLVED | a, X, ε, m(a), E_ε all defined |
| "pilot" language | RESOLVED | Only the real education pilot remains |
| No "Because" openers | RESOLVED | None found |
| Application heading/colons/rules | PARTLY | Heading colon; rounding rule unstated |
| Replication extension | PARTLY | "increased replication until..." gives no step size or maximum |
| Results choppy; terms defined; targets named | PARTLY | Targets named; "identified", "Wald" undefined |
| MSE ranges | PARTLY | Most ranges say what they span; "−0.12 to −0.40" does not |
| Run totals | PARTLY | Millions gone; "14,400" and "1,248" remain |
| "broad regime distinction"; ranking transfer | PARTLY | Now "kept the regime distinction" (S9); transfer sentence clear |
| No "grid pilot" | RESOLVED | Absent |
| "Budgets differed (...)" | RESOLVED | Now a sentence, but it repeats Table 1 (S10) |
| Key results only | PARTLY | Pooled-ranking and budget paragraphs (S10) |
| Plain Limitations/Conclusions | PARTLY | "complement spatial analysis", "Wald", "worst-case" |
| Endpoint and staffing resources | RESOLVED | "define the education endpoint... specify the staff and other resources" |
| Figure 2 vs 3 distinct; Fig 3 caption length | PARTLY | Length fixed; distinction still missing (S15) |
| Figure 2 CIs | RESOLVED | "95% Monte Carlo intervals... often narrower than the symbols" |
| Table 2 SRS 50 of 100 | RESOLVED | In caption |
| Overall flow/reproducibility | PARTLY | Much improved; S1, S2 and S14 gaps remain |
