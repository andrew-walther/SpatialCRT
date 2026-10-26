# Talk outline: allocating a statewide SUD-prevention education trial

20-minute talk, Thursday 2026-10-29, for clinicians and epidemiologists. Draft 1
(2026-10-25), built from Chapter 3 and the methods-extension results. It is the
first version of prelim Part II, so each slide is tagged **prelim** (carries over
unchanged), **technical** (needs a more technical version for biostatisticians) or
**talk-only**.

The visual version is `outline.html` (built from this file by
`tools/build_outline_html.py`; figures in `figures/`). Edit this file, then rebuild.

## Open questions for Andrew

- Venue, audience size and how much of the 20 minutes is Q&A? The timing below
  assumes about 17 minutes of talk and 3 of questions.
- Slide format: Quarto revealjs at 1920×1080 with `unc-gillings.css`, like the prelim
  test deck, or another format?
- Which feedback do you most want from this group? Slide 19 lists draft questions.
- Show the NC incidence map by year (slide 10)? This group may approve the data use,
  so the full map is assumed fine.
- Use the plain names "spillover to controls only" and "spillover to everyone" on
  slides, or keep "control-only" and "both-arms"?

## Timing

| Section | Slides | Minutes |
|---|---|---|
| 1. Motivation and background | 1–4 | 3 |
| 2. Models | 5–6 | 2 |
| 3. Simulation study and designs | 7–9 | 3.5 |
| 4. NC application | 10–16 | 7 |
| 5. Conclusions, limitations, next steps | 17–19 | 2.5 |
| Backup | B1–B8 | as needed |

## 1. Motivation and background

### 1. Where we treat matters when an intervention spills over
- Figure: none (title slide)
- Slide text: Title, authors (Walther, Habib, Simpson, Lin), UNC Biostatistics, date.
- Notes: One sentence of purpose. We are planning how to choose which community colleges receive an SUD-prevention education program first, and the choice affects how well a trial can measure the program's effect.
- Time: 0.5
- Carry-over: talk-only

### 2. Neighbors share more than a border
- Figure: none (two short examples with icons)
- Slide text: Bed nets in northern Ghana lowered child mortality in nearby untreated villages. School deworming in Kenya helped untreated children near treated schools. When treatment reaches neighbors, a simple treated-versus-control comparison can understate the effect.
- Notes: Spillover breaks the usual cluster-trial assumption that one cluster's outcome ignores its neighbors' assignments. Cite Binka 1998 and Miguel and Kremer 2004 on the slide footer.
- Time: 1
- Carry-over: prelim

### 3. Analysis methods exist, design guidance is thin
- Figure: none (a two-column contrast: analysis-stage tools vs design-stage choices)
- Slide text: Spatial models can separate direct and spillover effects after the trial. Few studies say which clusters to treat before the trial starts. Our earlier small-grid study (BMC Medical Research Methodology) compared block-stratified and simple random allocation.
- Notes: Position the project as design-stage. Keep the literature to one slide for this audience; the prelim version can expand it.
- Time: 0.75
- Carry-over: technical
- Prelim note: Expand to a method/year/author table ending with this project.

### 4. The question: which colleges should get the program first?
- Figure: community_college_service_area_clusters.png
- Slide text: 58 community college service areas cover North Carolina's 100 counties. A single-county pilot is moving to a statewide trial. Sudden unexpected death (SUD) incidence differs across the state.
- Notes: Introduce Dr. Simpson's program, the education endpoint (an assessment of trained first responders and clinicians, averaged by college) and why incidence might guide allocation. The trial measures education, not SUD deaths.
- Time: 0.75
- Carry-over: prelim

## 2. Models

### 5. Spillover can reach only the controls, or everyone
- Figure: none (schematic of a treated college and its neighbors under the two regimes; to draw)
- Slide text: Control-only spillover means untreated colleges benefit from treated neighbors. Both-arms spillover means every college responds to its neighbors' treatment. The best design depends on which one is true.
- Notes: This is the single most important idea in the talk; return to it in every results slide. Exposure is the share of a college's neighbors that are treated.
- Time: 1
- Carry-over: prelim

### 6. We estimate the effect with a model that knows about neighbors
- Figure: none (one equation in words, with an optional small equation)
- Slide text: Outcome = direct effect of treatment + spillover from treated neighbors + spatial correlation + noise. A spatial lag model estimates the direct effect. We judge each design by its average error (MSE) and its error under an unlucky allocation.
- Notes: Keep the math to a sentence for this audience. Mention that the model includes the spillover term, which is a best case. Define MSE as average squared distance between the estimate and the truth.
- Time: 1
- Carry-over: technical
- Prelim note: Show the equation Y = (I − ρW)^{-1}[τZ + S(Z) + ε] and both spillover forms.

## 3. Simulation study and designs

### 7. Nine ways to choose which clusters to treat
- Figure: design_assignments_panels_wide.png
- Slide text: Simple random sampling (SRS) is the default. Others balance incidence (Balanced Quartiles), target it (High Incidence Focus), separate treated clusters (Checkerboard, Isolation Buffer, 2×2 Blocking) or vary the treated share by region (Saturation Quadrants, Incidence-Guided Saturation).
- Notes: Walk through two or three panels only. Saturation designs treat 80, 60, 40 and 20 percent of four regions; the incidence-guided version gives the highest share to the highest-incidence region.
- Time: 1.25
- Carry-over: prelim

### 8. We tested them on a 100-cluster grid under many conditions
- Figure: none (small parameter table)
- Slide text: Three ways incidence varies across space, four strengths of spatial correlation, four spillover sizes, two spillover regimes, five effect sizes. 250 simulated trials per scenario.
- Notes: One sentence on why a grid first. It lets us vary conditions systematically before moving to the real map.
- Time: 0.75
- Carry-over: technical
- Prelim note: Add queen versus rook contiguity, the oracle and non-oracle fits and the seed design.

### 9. On the grid, the spillover regime decides which design wins
- Figure: fig_vs_srs_by_regime_queen.png
- Slide text: Spillover to controls only means the saturation designs had about 37% lower error than SRS. Spillover to everyone means SRS was as good as any design; the saturation designs were 4–8% worse on average, with similar worst-case allocations.
- Notes: Each point is a design's error relative to SRS (below 1 is better), averaged over 50 independent units; filled points are significant after multiple-comparison adjustment. Mention that SRS's worst allocations are much worse than its average under control-only spillover (slide 13). Checkerboard and 2×2 Blocking look fine under one regime and fail badly under the other. Alternative figure if a simpler view is wanted is grid_accuracy_srs.png.
- Time: 1.5
- Carry-over: prelim

## 4. NC application

### 10. SUD incidence in NC is uneven, and the pattern shifts year to year
- Figure: observed_years.png
- Slide text: SUD deaths per 100,000 adults aged 18–64 by college service area, 2018–2021. Each year is a separate planning input.
- Notes: 23,523 deaths over 25.6 million person-years, from Habib's corrected county counts (Gan et al. case definition). Incidence guides who is treated; it does not enter the education outcome.
- Time: 1
- Carry-over: prelim

### 11. Each design was adapted to the real map before we looked at results
- Figure: design_examples.png
- Slide text: Four fixed regions for the saturation designs. A graph version of Checkerboard. Treated budgets differ, and we report them rather than force them equal.
- Notes: Mention that regions were chosen for compactness and balance without looking at outcomes. regions.png is an alternative figure.
- Time: 0.75
- Carry-over: prelim

### 12. With all 58 colleges, saturation beats SRS when spillover reaches only controls
- Figure: nc_annual_accuracy.png
- Slide text: Spillover to controls only means Incidence-Guided Saturation had 21–30% lower error than SRS, depending on the year. Spillover to everyone means SRS and several designs were about equal.
- Notes: Each cell is a design's error divided by SRS's in the same year and regime; blue is better than SRS. Same message as the grid.
- Time: 1
- Carry-over: prelim

### 13. Saturation also protects against an unlucky allocation
- Figure: allocation_distribution_by_design.png
- Slide text: A trial gets one allocation. Under control-only spillover, most saturation allocations beat the average SRS allocation, and the worst 10% were 21–43% better than SRS's worst 10%.
- Notes: The worst-decile estimates are uncertain because only ten sampled allocations form the worst decile. Under both-arms spillover the tail differences are small and noisy.
- Time: 1
- Carry-over: prelim

### 14. About 36 colleges detect a one-standard-deviation effect under spillover to everyone
- Figure: fig_colleges_needed.png
- Slide text: Spillover to everyone means 35–37 colleges for 80% power with any design except Checkerboard and Isolation Buffer. Spillover to controls only means even 58 colleges are not enough at that effect; at 1.5 standard deviations the saturation designs need 28 colleges, against 32 for SRS.
- Notes: Effects are in college-level standard deviations. For 25 trained professionals per college and an intraclass correlation of 0.05, a one-SD college effect is about 0.3 SD for an individual; d = 0.3 gives power 0.71 with 29 colleges and 0.96 with 58 under spillover to everyone. These numbers use a small-sample t test that keeps the false-positive rate at 5%.
- Time: 1.25
- Carry-over: prelim

### 15. If fewer colleges join, saturation's advantage shrinks
- Figure: fig_power_vs_n.png
- Slide text: With 58 colleges the saturation designs have 22% lower error than SRS under control-only spillover. With 29 colleges, 4–5%. With 20, about 1%. Under control-only spillover, adding colleges past about 40 helps little.
- Notes: When nearly every college joins, nearly every control has a treated neighbor, so spillover is hard to separate from treatment. Non-participating neighbors restore some contrast. If high-incidence colleges are more likely to volunteer, the advantage at 20–29 colleges is 2% or less. Ask the audience how many colleges they expect to join.
- Time: 1
- Carry-over: prelim

### 16. Incidence guidance redirects treatment only where regions differ
- Figure: fig_spread_vs_gain.png
- Slide text: On the NC map, Incidence-Guided Saturation treats colleges with 51% of the incidence share, against 50% for plain saturation. On the grid, where regions differ more, the gain is 2–3 percentage points and grows with the regional spread.
- Notes: NC's four regions have similar average incidence. Death share runs the other way because the highest-rate colleges are small (High Incidence Focus reaches 62% of rate share but 35% of deaths). Guidance costs no accuracy, so it is a reasonable ethical choice, but it should not be oversold on this map.
- Time: 1
- Carry-over: prelim

## 5. Conclusions, limitations, next steps

### 17. The best design depends mainly on how spillover travels
- Figure: none (take-home layout: claim, highlight box, three bullets)
- Slide text: Highlight box: choose the design for the spillover you expect. Bullets: saturation when spillover reaches mainly controls; SRS or balanced designs when it reaches everyone; 35–37 colleges detect a one-SD effect under spillover to everyone.
- Notes: No single recommended design. Incidence-Guided Saturation matches plain saturation's accuracy while sending more treatment to high-burden regions, when regions differ.
- Time: 1
- Carry-over: prelim

### 18. What these results do not tell us
- Figure: none (two columns: limitations, future work)
- Slide text: Limitations: the model knows the spillover form; planning incidence is treated as exact; contiguity stands in for professional networks. Future work: a future-incidence endpoint, noisy incidence, effect heterogeneity.
- Notes: Mention that the results held when we fitted a CAR model or a Bayesian model (backup B4) and that a small-sample interval fixes the slight undercoverage (backup B2).
- Time: 0.75
- Carry-over: prelim

### 19. What we would like from you
- Figure: none
- Slide text: Do you expect spillover to reach mainly untreated colleges, or everyone? How many colleges are likely to join, and would high-incidence colleges be more likely to? How many professionals per college, and how variable is the education score within a college?
- Notes: These three answers decide the design and the sample size. Thank Dr. Simpson, Ashkan Habib and Dr. Lin.
- Time: 0.75
- Carry-over: talk-only

## Backup

### B1. Formal tests agree with the descriptive ranking, with SRS included
- Figure: app_cd_diagram_rank_queen.png
- Slide text: Friedman test across six designs and SRS on 160 grid blocks, χ² = 560, df 6. Incidence-Guided Saturation ranks first; SRS and Balanced Quartiles cannot be separated.
- Notes: Question it answers is "Are the design differences statistically significant?"
- Carry-over: technical

### B2. A small-sample t interval restores nominal coverage
- Figure: none (table from Appendix A9, coverage by number of colleges)
- Slide text: Wald tests reject a true null 11–13% of the time with 15 colleges and 6–7% with 58. The corrected t test rejects 4.4–5.6% at every size, with coverage 94–96%.
- Notes: Answers "Why was coverage below 95%, and can it be fixed?"
- Carry-over: technical

### B3. From a college-level effect to trained professionals
- Figure: none (table tab:nc-translation)
- Slide text: σ_c = √(ICC + (1 − ICC)/m). With m = 25 and ICC = 0.05, σ_c ≈ 0.30.
- Notes: Answers "What does an effect of one standard deviation mean in practice?"
- Carry-over: prelim

### B4. The guidance holds under CAR and Bayesian models
- Figure: none (table tab:nc-bayes)
- Slide text: CAR fitting or CAR data leave saturation at 0.76–0.82 of SRS's error under control-only spillover. A Bayesian spatial lag model reproduces the maximum likelihood results. Only a Bayesian CAR fitted to strongly correlated spatial lag data loses the saturation advantage.
- Notes: Answers "Did you try a Bayesian model?" Cite Leroux et al. 2000, Lee 2013 (CARBayes) and LeSage and Pace 2009 on the slide.
- Carry-over: technical

### B5. Rate share and death share pull apart in NC
- Figure: none (table of rate share and death share by design)
- Slide text: High Incidence Focus treats 62% of the rate share but 35% of deaths, because high-rate colleges are small.
- Notes: Answers "Does incidence-guided allocation reach more deaths?"
- Carry-over: prelim

### B6. Checkerboard fails for different reasons on the grid and the map
- Figure: fig_mse_by_design_ranked_6design_queen.png
- Slide text: Under queen contiguity, Checkerboard's estimates are unbiased but imprecise under control-only spillover. Under rook contiguity, the effect is not identified.
- Notes: Answers "Why does strict alternation do so badly?"
- Carry-over: technical

### B7. Binary spillover makes Isolation Buffer unidentifiable
- Figure: none (table tab:nc-forms)
- Slide text: If any treated neighbor counts the same as all, every Isolation Buffer control is exposed, so the effect cannot be separated.
- Notes: Answers "What if spillover is all-or-nothing?"
- Carry-over: technical

### B8. Unequal college sizes cost power equally for every design
- Figure: none
- Slide text: Letting the number of trained professionals grow with college size lowers power by about 0.1 for every design and changes no design's standing.
- Notes: Answers "What if colleges train different numbers of professionals?"
- Carry-over: prelim
