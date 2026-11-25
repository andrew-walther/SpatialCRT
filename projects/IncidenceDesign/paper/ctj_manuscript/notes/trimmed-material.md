# CTJ trimmed material

Author: Andrew Walther. Text removed from `CTJ_Manuscript.tex` for length, with where its content
now lives. Never delete entries; the 1h review package for Dr. Lin cites this file.

## Round 6 (applied 2026-11-25; decisions in `revisions/2026-11-24-round6-word-limit/decisions.md`)

Base: V6 (`4386bcb`, tag `ctj-v6`). "SI" means `supplementary_info/Supplementary_Information.tex`
at the time of the cut; SI section numbers will change in Phase 2 (item 2b links them to labels).

| ID | Section | Removed text | Where the content lives now |
|---|---|---|---|
| T14 | Introduction | "We do not compare numerical results between the two studies." | Nowhere; the rule is kept by making no such comparison (AGENTS 2026-10-02) |
| T12 | Methods, estimation | "Our own ML routine reproduced the estimates of the \texttt{lagsarlm} function in the R package \texttt{spatialreg} to within $10^{-6}$." | SI Section 3 (Estimation and identification), which gives the tighter $9.7\times10^{-8}$ agreement |
| C10 | Methods, estimation | "(SI Section 15.2)" after "chosen by calibration" | Same pointer two sentences later in the same paragraph |
| C-16 | Methods, estimation | "as the HC1 covariance estimator does" | Replaced by "a standard small-sample correction" (MacKinnon cite kept); Chapter 3 keeps the HC1 wording |
| C4 | Methods, designs | "SRS randomly selects 50 of the 100 clusters for treatment, and … except Isolation Buffer, whose number treated depended on the neighbor graph." | SRS: Introduction, Table 1 and Table 2 caption. Isolation Buffer counts: Table 1 |
| T4 | Methods, designs | "The designs differ in how they handle incidence and spatial exposure. Balanced Quartiles spreads treatment evenly across ranked incidence groups, and High Incidence Focus aligns treatment with incidence. $2\times2$ Blocking imposes the same local proportion of treated clusters everywhere, while saturation deliberately varies that proportion between regions. Isolation Buffer prevents treated neighbors, and Checkerboard imposes strong interspersion." | Table 1 (category and rule of each design); saturation's varying local share is in the Discussion (T3 text) |
| T11 | Methods, simulation study | "All configurations put incidence on the $(0,1)$ scale, and one value of $\tau$ therefore sets the same shift between the lowest and highest possible incidence in every configuration, although the spread between those limits differs somewhat by mechanism." | Shortened into the rank-scaling sentence; the "spread differs somewhat by mechanism" caveat (committee C-31) stays in Chapter 3 and this file |
| N3 | Methods, simulation study | "Every performance metric is estimated from a finite number of simulated datasets and therefore carries simulation (Monte Carlo) error. Its Monte Carlo standard error (MCSE) measures this uncertainty and shrinks as more datasets are simulated." | Rewritten as one sentence; the definition is kept |
| N5 | Methods, simulation study | "The value may need adjusting for a particular outcome." | This file; the τ sensitivity analysis (Phase 4, SI) will cover other values |
| T5 | Methods, case-study data | "obtained under University of North Carolina at Chapel Hill Institutional Review Board protocol 14-2036" and "The corrected counts include heart-failure deaths that had previously been excluded as presumed non-sudden," | IRB protocol: Ethical considerations (CTJ end matter). Heart-failure correction: SI Section 5.1 |
| C1 | Methods, case-study outcome | "The fit included an intercept, treatment and the true spillover term but not incidence." | Methods, estimation ("and, in the grid study, incidence"); the τ = 1 sensitivity sentence |
| T6 | Methods, case-study adaptations | "Checkerboard became a deterministic rule that visits clusters from most to fewest neighbors and gives each cluster the arm held by fewer of its already-assigned neighbors." | Table 1, NC column (same rule) |
| T6b | Methods, case-study adaptations | "with fractional halves rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), so a one-cluster group was never treated" | SI Section 5.2, sentence edited to state the rule (approved 2026-11-25) |
| N7 | Methods, budgets | "Treating half the clusters need not reach half the population, so" | This file; the reporting of population share makes the point |
| N6 | Methods, case-study settings | "Designs whose assignments ignore incidence, such as SRS, have the same possible assignments every year, so their settings for different years shared one simulation." | This file (simulation bookkeeping; the author judged it unnecessary) |
| T2 | Methods, replication | "for each of the same assignments, so that each assignment's risk was estimated precisely enough to rank. As a check against overstating worst-decile error by selecting on noise," | SI Section 5.4 (split-half check and 400 outcomes) |
| T2b | Results, case study | "an MCSE of mean MSE at most 5\% of its value and an MCSE of coverage at most 0.01" | Methods, replication paragraph (precision targets) |
| T7 | Results, case-study allocation risk | "so which assignments belonged to it was uncertain in every setting" and "the two outcome halves selected largely different worst deciles for Incidence-Guided Saturation" | SI Section 14 (tail refinement flags uncertain worst-decile membership in all 112 settings) |
| C5 | Results, case-study allocation risk | "Checkerboard and High Incidence Focus each have a single possible assignment for these inputs, so they show no variation across assignments but can still have substantial error." | Figure 3 caption ("appears as a circle only") and Discussion, implications |
| C8 | Results, budgets | "Treatment budgets differed by design and year (Table 1 and SI Table S5)." | Methods, budgets; SI Table S5 cited in the next sentence |
| N12 | Results, case-study sensitivity | "Alternative incidence summaries reordered some saturation regions, and the prespecified mean-rank rule was kept." | SI Section 5.2 (regional orderings table) and Section 14; pointer added in Methods (N12b) |
| T13 | Discussion | "In the case study, its advantage over SRS fell from 22\% with 58 colleges to 5\% with 29." | Results, colleges-needed paragraph (same numbers) |
| T3 | Discussion | "Incidence can also target treatment directly, but … rather than from randomization …" and "and is consistent with its control-only advantage" | Merged paragraph keeps the High Incidence Focus and regime points |
| C6 | Limitations | "and omitting the spillover term produced biased estimates and poor coverage" | Results, grid (omitted-spillover sentences) |
| N13 | Continued work | "Isolation Buffer treats fewer clusters and High Incidence Focus reaches only 22--24\% of the working-age population, so these results do not identify the most efficient design for a fixed implementation cost or number of participants." | Results, budgets (reach figures); Limitations (not comparisons at equal cost) |
| C7 | Continued work | "Under control-only spillover even all 58 colleges gave less than 80\% power at $\beta=1$," | Results, colleges-needed paragraph |
