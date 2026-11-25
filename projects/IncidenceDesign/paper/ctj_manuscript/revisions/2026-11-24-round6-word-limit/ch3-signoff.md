# Round 6, Phase 1e: Chapter 3 carry-over (sign-off list)

Author: Andrew Walther. Prepared 2026-11-25 from `ch3-carryover-candidates.md`, checked against
`paper/thesis_chapter/Thesis_Chapter.qmd` (lines as of `d7f7ed5`). Nothing applied yet.

Not needed, because Chapter 3 already has the content: R2 (line 309 gives each design's MSE range
across the five configurations), R3 / Q-01 (binary spillover, lines 87 and 462), HW-19 (NC data only
in the case study, line 91). Most CTJ HW edits do not apply because Chapter 3 words those
sentences differently. C-16 and T11's caveat stay as they are in Chapter 3 (decided in 1b).

Applied to `Thesis_Chapter.qmd` (the hook regenerates the `prelim/` copy), then merged by hand
into `bios-dissertation/thesis/chapter3-project2-incidence/project2-incidence.qmd` (thesis/README.md).
The Lee (2004) entry goes into `bios-dissertation/prelim/references.bib` and
`bios-dissertation/thesis/references.bib`. The frozen `prelim-submitted-2026-11-15/` is not touched.

## Content items

| ID | Line | Current | Proposed |
|---|---|---|---|
| CH-1 (A1) | 148 | The 95\% Wald interval uses $q = z_{0.975} = 1.96$ without the factor, | The 95\% Wald interval [@lee_asymptotic_2004] uses $q = z_{0.975} = 1.96$ without the factor, |
| CH-2 (TAU-1) | 78 | The incidence coefficient $\tau$ sets how strongly cluster-level incidence shifts the outcome. Its values and the reason for them are given in @sec-simulation. | The incidence coefficient $\tau$ sets how strongly cluster-level incidence shifts the outcome, holding treatment and spillover fixed. With $\tau = 0$, incidence guides only the allocation, as in the case study. With $\tau > 0$, clusters with higher incidence have higher outcomes, as on the grid. Its values and the reason for them are given in @sec-simulation. |
| CH-3 (N5) | 162 | The incidence coefficient $\tau$ is set by the investigator according to whether incidence is expected to shift the outcome, and it is not estimated at the design stage. On the grid we set $\tau = 1$. Incidence lies between 0 and 1, and this value therefore raises the outcome, before spatial spread, by one noise standard deviation from the lowest to the highest incidence, the same size as the primary intervention effect $\beta = 1$. | The investigator specifies $\tau$ according to whether incidence is expected to shift the outcome. It is not estimated at the design stage, and we set $\tau = 1$ in the grid simulation. Incidence ranges from 0 to 1, and this value raises the outcome by one noise standard deviation from the lowest to the highest incidence before spatial spread. That shift equals the primary intervention effect $\beta = 1$. |
| CH-4 (T6b) | 416 | Half of each group, rounded, is treated at random. Treatment budgets, the numbers of treated clusters, therefore differ across designs. | Half of each group, rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), is treated at random. Therefore, treatment budgets, the numbers of treated clusters, differ across designs. |

## "Therefore"/"However" at sentence start (author style, 2026-11-25)

Two pairs would put "Therefore," on consecutive sentences (TS-12/13, TS-29/30); the second of each
is reworded without it.

| ID | Line | Current | Proposed |
|---|---|---|---|
| TS-01 | 51 | The literature therefore offers much more guidance… | Therefore, the literature offers much more guidance… |
| TS-02 (TH-1) | 78 | The intervention effect $\beta$ is therefore in units of the noise standard deviation. | Therefore, the intervention effect $\beta$ is in units of the noise standard deviation. |
| TS-03 | 87 | …regardless of its own assignment, and neighboring treatment therefore adds to a treated cluster's own effect. | …regardless of its own assignment. Therefore, neighboring treatment adds to a treated cluster's own effect. |
| TS-04 | 95 | Every mechanism puts incidence on the $(0, 1)$ scale, and one value of $\tau$ therefore sets the same outcome shift… | Every mechanism puts incidence on the $(0, 1)$ scale. Therefore, one value of $\tau$ sets the same outcome shift… |
| TS-05 | 146 | Its results are therefore a best case for each design. | Therefore, its results are a best case for each design. |
| TS-06 | 152 | A trial realizes only one allocation, and we therefore judge each design on two criteria: the error of the estimate averaged over the allocations the design can produce, and how that error varies from one allocation to another. | A trial draws only one allocation. Therefore, we judge each design on two criteria, the error of the estimate averaged over the allocations the design can produce and how that error varies from one allocation to another. |
| TS-08 | 168 | Each metric is estimated from a finite number of simulated datasets and therefore carries Monte Carlo error. | Each metric is estimated from a finite number of simulated datasets. Therefore, it carries Monte Carlo error. |
| TS-09 | 241 | The two regimes in Equation \ref{eq:spillover} therefore define the same outcome model and the same expected MSE. | Therefore, the two regimes in Equation \ref{eq:spillover} define the same outcome model and the same expected MSE. |
| TS-10 | 279 | The oracle model includes the true regime-specific spillover covariate and therefore has no omitted-variable bias. | The oracle model includes the true regime-specific spillover covariate. Therefore, it has no omitted-variable bias. |
| TS-11 | 348 | The five values of $\beta$ share their assignments and noise, and these comparisons are therefore not independent confirmations of one another. | The five values of $\beta$ share their assignments and noise. Therefore, these comparisons are not independent confirmations of one another. |
| TS-12 | 370 | The grid study did not use split-half selection, and its worst-decile means are therefore biased upward for every design, … | The grid study did not use split-half selection. Therefore, its worst-decile means are biased upward for every design, … |
| TS-13 | 370 | We therefore use the NC case study in @sec-application as a setting-specific assessment instead of extrapolating the size of the tail protection from these settings. | Instead of extrapolating the size of the tail protection from these settings, we use the NC case study in @sec-application as a setting-specific assessment. |
| TS-14 | 378 | Both average estimation error and the risk of a poor allocation therefore matter. | Therefore, both average estimation error and the risk of a poor allocation matter. |
| TS-16 | 466 | We therefore repeated the case study for subsets… | Therefore, we repeated the case study for subsets… |
| TS-17 | 554 | This design treats the highest-incidence clusters, and adjusting for incidence therefore leaves less variation… | This design treats the highest-incidence clusters. Therefore, adjusting for incidence leaves less variation… |
| TS-18 | 554 | …it gives the same set of possible assignments and therefore the same results. | …it gives the same set of possible assignments and the same results. |
| TS-19 | 560 | Under the spatial lag model the outcome mean includes the spatial multiplier, and the CAR coefficient on treatment therefore estimates a different quantity from $\beta$ [@wall_close_2004]. | Under the spatial lag model the outcome mean includes the spatial multiplier. Therefore, the CAR coefficient on treatment estimates a different quantity from $\beta$ [@wall_close_2004]. |
| TS-20 | 581 | The colleges-needed numbers above therefore come from a test with the nominal size. | Therefore, the colleges-needed numbers above come from a test with the nominal size. |
| TS-21 | 593 | …for every design whose $\beta$ was identifiable, and the designs therefore differed in precision rather than in coverage. | …for every design whose $\beta$ was identifiable. Therefore, the designs differed in precision rather than in coverage. |
| TS-22 | 599 | …(Table \ref{tab:downside}), and they are therefore the designs we recommend over SRS. | …(Table \ref{tab:downside}). Therefore, they are the designs we recommend over SRS. |
| TS-23 | 599 | SRS and the saturation designs are therefore both reasonable choices. | Therefore, SRS and the saturation designs are both reasonable choices. |
| TS-24 | 607 | Their accuracy advantage over SRS under control-only spillover therefore appears to come from… | Therefore, their accuracy advantage over SRS under control-only spillover appears to come from… |
| TS-25 | 613 | It includes the exact spillover exposure term, and its results are therefore a best case for each design. | It includes the exact spillover exposure term. Therefore, its results are a best case for each design. |
| TS-26 | 617 | Their observed rates therefore differ in precision, … | Therefore, their observed rates differ in precision, … |
| TS-27 | 625 | Uncertainty in their measurement and changes between planning and implementation are therefore not propagated. | Therefore, uncertainty in their measurement and changes between planning and implementation are not propagated. |
| TS-28 | 714 | Results from different code or settings are therefore never combined. | Therefore, results from different code or settings are never combined. |
| TS-29 | 855 | The college mean therefore has standard deviation… | Therefore, the college mean has standard deviation… |
| TS-30 | 855 | A college-level effect of $\beta$ standard deviations is therefore an individual-level standardized effect $d = \beta\sigma_c$. | A college-level effect of $\beta$ standard deviations then corresponds to an individual-level standardized effect $d = \beta\sigma_c$. |
| TS-31 | 869 | The variance due to outcome noise therefore involves squared weights… | Therefore, the variance due to outcome noise involves squared weights… |
| TS-32 | 904 | Retaining Incidence-Guided Saturation as the representative of the saturation mechanism therefore does not rest on performance. | Therefore, retaining Incidence-Guided Saturation as the representative of the saturation mechanism does not rest on performance. |
| TS-33 | 961 | We therefore also ran paired Wilcoxon signed-rank tests… | Therefore, we also ran paired Wilcoxon signed-rank tests… |
| TS-34 | 1279 | Their rows therefore repeat and are not independent annual confirmation. | Therefore, their rows repeat and are not independent annual confirmation. |
| TS-35 | 1509 | They therefore describe the consistency of the ordering rather than independent evidence. | Therefore, they describe the consistency of the ordering rather than independent evidence. |

TS-07 (line 162) is covered by CH-3 and TS-15 (line 416) by CH-4.
