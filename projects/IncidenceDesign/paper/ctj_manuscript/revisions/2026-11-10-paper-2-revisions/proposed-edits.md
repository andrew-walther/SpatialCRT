# CTJ first revision: proposed edits for sign-off

Author: Andrew Walther. Prepared 2026-11-15 (Phase 4a of
`docs/plans/2026-11-14-chapter3-review-session-plan.md`).

Line numbers are in `manuscript_ctj/CTJ_Manuscript.tex` (C), `CTJ_Exhibits.tex` (E) and
`supplementary_info/Supplementary_Information.tex` (S) as of SpatialCRT `74f8d48`. "T-" refers
to the triage rows in the plan; "FCL-" to Dr. Lin's annotations.

**Decisions received 2026-11-15** (Andrew): all accepted except the edits recorded in the Decision column. **How to answer (original).** Mark each item **A** (accept), **E** (accept with your edit) or **D**
(decline), or write "accept all except ...". Proposed text is written in the new notation
(item P-01): β is the intervention effect, τ the incidence coefficient, ψ the spillover
magnitude. Nothing below changes a number except P-30 and P-31, which correct or add values
from the stored result files named there.

Word count: the main text is at 3,498 of 3,500 words. Most items add words; the limit is
deferred until after the proposal (Decision C).

## P. Notation (applied first, as its own commit)

| ID | Rows | Change | Decision |
|---|---|---|---|
| P-01 | T-11, M-01, D-1, FCL-030 | One-pass swap in C, E and S: τ → β (intervention effect), β → τ (incidence coefficient), γ → ψ (spillover magnitude), including hats, subscripts and derived forms (τ − γ → β − ψ, γ = 0.8τ → ψ = 0.8β, γ = κτ → ψ = κβ). The spillover term is written ψ S(**Z**) with S(**Z**) = **WZ** (both arms) or (**1** − **Z**) ∘ **WZ** (control only). Checked by grep afterward: no τ left on **Z**, no β left on **X**, no γ anywhere || A |
| P-02 | T-11 | S 336–338: the earlier study's "β̂ ≈ β − ψ in its notation" now uses the same symbols as this paper. Replace with "the same consequence of exact collinearity between treatment and spillover exposure" (drop "in its notation") || A |
| P-03 | T-11 | Figure label strings only (no analysis code): "True τ" → "True β" in `code/14_manuscript_supplement_figures.R` (lines 249, 252, 585) and `code/12_six_design_statistical_comparisons.R` (252, 266, 297, 322); legend γ → ψ in `code/19_manuscript_exhibits.R:62`. Rebuild only the figures those lines draw, from stored results; script and result names (`MLE_tau_sweep` and so on) stay as they are. If any rebuild needs a simulation rerun, I stop and report instead || A |
| P-04 | T-45, FCL-056 | E Table 2 caption "the true spillover term" → "the true spillover covariate S(**Z**), with ψ estimated". The code fits the column ψ**WZ** with coefficient 1; fitting **WZ** and estimating ψ gives the same β̂ and SE, so the description is exact || A |

## I. Introduction and abstract

| ID | Rows | Current (C line) | Proposed | Decision |
|---|---|---|---|---|
| I-01 | T-01, FCL-013 | 141 "Spatial interference challenges the assumption that ..." | "Spatial interference violates the assumption that a cluster's outcome depends only on its own treatment." || A |
| I-02 | T-02, FCL-008, I-01 | 141 "... very different information about the treatment effect." (no citation) | Add `\cite{baird_optimal_2018,leung_rateoptimal_2022,walther_spatialcrt_2026}` at the end (all three already in the bib) || A |
| I-03 | T-03, FCL-001, FCL-009 | 143 first sentence | Dr. Lin's wording: "Design-stage methods, including covariate-constrained randomization, randomized saturation, and independent-set and spatial cluster designs, serve different objectives.\cite{...}" (citation at the end) || A |
| I-04 | T-04, FCL-002, -003, -004, -010 | 143 "In a strict checkerboard with edge neighbors, for example, a cluster's exposure to treated neighbors is fully determined by its own treatment." | "For example, in a checkerboard design with neighbors defined by shared edges, every neighbor of a treated cluster is a control and every neighbor of a control is treated. A cluster's share of treated neighbors is then fixed by its own arm, and the treatment effect cannot be separated from the spillover effect.\cite{walther_spatialcrt_2026}" || A |
| I-05 | T-05 to T-07, FCL-005, -006, -007, -011 | 145 whole paragraph | "In previous work on small grids, we compared simple random sampling (SRS), meaning complete randomization of a fixed number of clusters, with block stratified sampling, an alternating allocation in which no two clusters that share an edge are in the same arm; we call this design Checkerboard here.\cite{walther_spatialcrt_2026} In that study, we judged each design by its average estimation accuracy and by the risk that the single allocation a trial draws gives a poor estimate. We extend the question to heterogeneous incidence, a larger study region, more designs and an irregular observed geography. In this paper, we compare six designs and SRS on the mean squared error (MSE) of the treatment-effect estimate and on how that error varies across possible assignments. We do not compare numerical results between the two studies." || A |
| I-06 | T-08, FCL-012 | 147 "Our application supports planning an education intervention ..." | "In this paper, we also plan an education pilot trial for healthcare professionals involved in sudden unexpected death (SUD) prevention, using a case study based on simulation with observed incidence." Rest of paragraph unchanged except "application" → "case study" (P-01 does not touch it; see M-14) || E: not planning the pilot; recommending a design for a statewide extension of a pilot study. Use: "In this paper, we also use a case study, a simulation based on observed incidence, to recommend a design for a statewide extension of a pilot study of education for healthcare professionals involved in sudden unexpected death (SUD) prevention." Abstract (I-08) and C 248 heading follow the same framing |
| I-07 | T-09, FCL-014 | Abstract "the error under an unfavorable realized assignment (allocation risk)" | "how large the error can be for the particular assignment a trial happens to draw (allocation risk)" || A |
| I-08 | T-10, FCL-015 | Abstract "education pilot" | "education pilot trial" || A |
| I-09 | T-37, D-3 | Abstract Methods, Results, Conclusions: "application" (4 uses) | "case study"; "The application used observed ..." → "The case study used observed ..." || A |
| I-10 | T-59, CC-01 | Abstract "covered the true effect 94.2--94.9\% ... 92.7--93.6\%" | "94.1--94.9\% ... 92.6--93.6\%" (see M-30) || D: the 94.2–94.9 / 92.7–93.6 range is the six designs and SRS, which is what the main text presents. Numbers stay; the sentence says "across the six designs and SRS" |
| I-11 | T-15, M-05 | Abstract "added an incidence-related baseline to the outcome" | "let incidence shift the outcome" || A |

## M. Methods

| ID | Rows | Current (C line) | Proposed | Decision |
|---|---|---|---|---|
| M-01 | T-13, M-03, T-17, FCL-017 | 151 queen/rook sentences and "the $i$th element of **WZ**" | "Clusters are linked by a row-standardized spatial weights matrix **W** with elements $w_{ij} = 1/d_i$ when clusters $i$ and $j$ are neighbors and $w_{ij}=0$ otherwise, where $d_i$ is the number of neighbors of cluster $i$. For a treatment indicator vector **Z**, $(\mathbf{WZ})_i$ is then the proportion of cluster $i$'s neighbors that are treated." The queen/rook sentences move to Simulation study (M-10) || A |
| M-02 | T-12, T-64, M-02, FCL-016, -018, -019, -020 | 152–163 two-line equation in the reading layout; "(so τ is in units of its standard deviation)" | One-line equation in both layouts, without the error distribution: $\mathbf{Y}=(\mathbf{I}-\rho\mathbf{W})^{-1}[\beta\mathbf{Z}+\psi\mathbf{S}(\mathbf{Z})+\tau\mathbf{X}+\boldsymbol{\varepsilon}]$, cited to `ord_estimation_1975`. Prose: "where **Y** is the vector of cluster outcomes, **X** is cluster-level incidence, **I** is the identity matrix and ρ is the spatial dependence coefficient. The noise $\boldsymbol\varepsilon$ is independent standard normal, and the intervention effect β is therefore in units of its standard deviation." || A |
| M-03 | T-14, M-04 | 163 "cluster incidence" | "cluster-level incidence" (in M-02) || A |
| M-04 | T-15, T-16, M-05, M-06, FCL-021, -022 | 163 "The coefficient β sets how strongly incidence shifts the baseline outcome. We set β=1 in the grid study ... Investigators can choose β to suit their own outcome." | "The incidence coefficient τ sets how strongly incidence shifts the outcome." The value and its reason move to Simulation study (M-11) || A |
| M-05 | T-18, M-13, FCL-023 | 165 "We estimated τ, giving τ̂, by maximum likelihood (ML), fitting ..." | "We estimated β by maximum likelihood (ML).\cite{bivand_comparing_2015} The fitted SAR lag model included an intercept, treatment, the true regime-specific spillover covariate S(**Z**) with ψ estimated, and, in the grid study, incidence." || A |
| M-06 | T-19, M-10, D-4, FCL-024 | 165 interval sentence | "We report 95\% confidence intervals $\hat\beta \pm q\,\widehat{\mathrm{SE}}\sqrt{N/(N-p)}$, where $\widehat{\mathrm{SE}}$ is the ML standard error, $N$ the number of clusters, $p$ the number of estimated coefficients including ρ, and $q = t_{0.975,\,N-p}$. The factor $\sqrt{N/(N-p)}$ replaces the ML variance divisor $N$ with $N-p$,\cite{mackinnon_heteroskedasticity_1985} and the $t$ quantile allows for the small number of clusters.\cite{leyrat_cluster_2018} In SI Section 15, this interval covered β more reliably than the Wald interval, which uses $q = z_{0.975} = 1.96$ without the factor." (two new bib entries, verified on Crossref) || A |
| M-07 | T-20, T-21, M-12, FCL-025, -026, -027 | 165 "Accuracy was measured by MSE, ..." | "Accuracy was measured by the mean squared error (MSE), the average of $(\hat\beta-\beta)^2$ over simulated trials, and coverage was measured as the proportion of 95\% confidence intervals containing β." || A |
| M-08 | T-22, M-08, FCL-028 | 167 "τ is not identified (cannot be estimated) ... Checkerboard is identified but can be imprecise" | "With rook contiguity, Checkerboard gives **WZ** = **1** − **Z**, so β is not identifiable and the treatment coefficient estimates β − ψ. With queen contiguity, β is identifiable under Checkerboard but can be estimated imprecisely: ..." (rest unchanged) || A |
| M-09 | T-23, T-24, M-09, M-11, FCL-014, FCL-029 | 169 m(a) definition; last sentence "A design with only one possible assignment has no spread ..." | Before the formula: "For each design, we also asked how large the MSE can be for a particular assignment." Then the formula as now (with β). Delete the last sentence || A |
| M-10 | T-17, M-07 | 180 start | Insert: "Neighbors were defined by queen contiguity (clusters sharing an edge or a corner) in the primary analysis and by rook contiguity (sharing an edge) in a sensitivity analysis.\cite{bivand_applied_2008}" || A |
| M-11 | T-16, D-2, Q-02, FCL-047, FCL-053 | 180 (new sentences after the parameter list) | "The incidence coefficient τ is set by the investigator according to whether incidence is expected to shift the outcome, and it is not estimated at the design stage. On the grid we set τ = 1. Incidence lies between 0 and 1, so this value raises the outcome by one noise standard deviation from the lowest to the highest incidence, the same size as the primary intervention effect β = 1. With incidence shifting the outcome, treatment that follows incidence is partly confounded with it, which is the harder case for designs that use incidence. The value may need adjusting for a particular outcome." || A |
| M-12 | T-25, FCL-031 | 180 "both spillover regimes and" | "both spillover regimes, and" || A |
| M-13 | T-26, T-28, T-30, S-02, S-04, FCL-032 to -036, FCL-042 | 180 configuration sentences | "An incidence configuration is the model that generates incidence over the grid; an incidence surface is one draw from it. The five configurations were independent Uniform(0, 1) incidence; spatial incidence $\mathbf{X}=\Phi\{(\mathbf{I}-\rho_X\mathbf{W})^{-1}\mathbf{u}\}$, with **u** independent standard normal, Φ the standard normal distribution function and $\rho_X = 0.2$ or 0.5; and Poisson incidence, in which death counts were drawn from Poisson distributions with log-rates $\log(35/100{,}000) + (\mathbf{I}-\rho_X\mathbf{W})^{-1}\mathbf{u}$ and 100,000 people per cluster, at the same two values of $\rho_X$ ..." then the existing rank-scaling sentence, followed by "All configurations put incidence on the (0, 1) scale, so one value of τ means the same shift in every configuration." || A |
| M-14 | T-29, S-01 | 180 "These are general simulation settings rather than a calibration to the application." | "These are general simulation settings; the North Carolina data enter only the case study." || E: "These are general simulation settings; the North Carolina data only enter the case study." |
| M-15 | T-30, FCL-035 | 180 "for surface $k$ the assignments, outcomes and fitted model all used the same incidence $\mathbf{X}_k$" | "Within surface $k$, the assignments, the simulated outcomes and the fitted model all used the same incidence $\mathbf{X}_k$." || A |
| M-16 | T-31, FCL-041, W-03 | 182 "and so carries simulation (Monte Carlo) error" | "and therefore carries simulation (Monte Carlo) error" || A |
| M-17 | T-32, FCL-043 | 184 "To estimate allocation risk on the grid, we used ..." | Keep in Methods (it describes a study design, not a result). Reword: "We also ran a separate allocation-risk study on the grid, using ..." || A |
| M-18 | T-37, D-3, FCL-048 | 186 heading "Application to sudden unexpected death in North Carolina"; "application" in Methods and Results (192, 194, 198, 200, 219, 220, 229, 252, 276) | Heading: "Case study: sudden unexpected death in North Carolina". "application" → "case study" throughout. Design stays in Methods and results in Results (meeting decision) || A |
| M-19 | T-38, A-02, FCL-045 | 187 "The observed rates were used as fixed inputs rather than simulated" | "Observed incidence and geography were fixed inputs, and outcomes were simulated (Figure~\ref{fig:nc-map})." || A |
| M-20 | T-39, FCL-046 | 192 "Equation~\ref{eq:grid}" | "Equation (\ref{eq:grid})" (SAGE Vancouver numbers equations in parentheses) || A |
| M-21 | T-16 | 192 "with β=0, τ=1 ... β=1 served as a hypothetical sensitivity analysis" | After P-01: "with τ = 0, β = 1 ... The education outcome is not expected to depend on SUD incidence, so incidence guided allocation only, and τ = 1 served as a hypothetical sensitivity analysis." || A |
| M-22 | T-34, T-35, A-03, FCL-039, FCL-040, FCL-054 | 194 Checkerboard and blocking sentences | Checkerboard: "Checkerboard became a deterministic rule that visits clusters from most to fewest neighbors and gives each cluster the arm held by fewer of its already-assigned neighbors." Regions: "Incidence-Guided Saturation used four fixed geographic regions. We ran $k$-means clustering on the cluster coordinates 250 times and kept the partition that best balanced compactness, population, number of clusters and number of counties." Blocking: "$2\times2$ Blocking became spatial blocking. Each year, $k$-means clustering on four standardized variables (location, that year's incidence rank and log population) formed 14 groups of similar nearby clusters, and a random half of each group was treated ..." || A |
| M-23 | T-15, I-11 | 198 "added the hypothetical incidence baseline (β=1)" | "let incidence shift the outcome (τ = 1)"; same at 235 ("Adding the hypothetical incidence baseline" → "Letting incidence shift the outcome") || A |
| M-24 | T-67, CC-09 | 200 "rejecting τ=0 when the 95\% interval excluded zero" | "rejecting β = 0 when the 95\% $t$-based interval excluded zero (a two-sided test at level 0.05; the grid power in the SI is one-sided at level 0.025)" || A |

## R. Results and Discussion

| ID | Rows | Current (C line) | Proposed | Decision |
|---|---|---|---|---|
| R-01 | T-41, FCL-049 | 204 "Incidence configuration changed each design's MSE by less than 10\% ..." second sentence | Move to the end of the paragraph || A |
| R-02 | T-42, FCL-050, FCL-051 | 204 "Pooled over both regimes, ... 0.091 against 0.128 for SRS, but separating the regimes was essential." | Delete the pooled sentence (not in Table 2). Open with "The two spillover regimes gave opposite results." || A |
| R-03 | T-22, FCL-028 | 204 "whose effect was identified, while rook Checkerboard was not identified" | "whose effect was identifiable, while under rook contiguity Checkerboard's was not identifiable" || A |
| R-04 | T-43, FCL-052 | "(SI)" at 211, 213, 217, 229 and "in the SI" at 222 | Cite the SI section or table each time (for example "(SI Table S12)"); exact labels filled from the SI at apply time || A |
| R-05 | T-65, CC-07, T-68 | 215 "(bias $-0.12$ to $-0.40$ and coverage 54--91\% ..." | "(bias $-0.123$ to $-0.400$ and $t$-based coverage 54--91\% ..." (S 1700–1703) || A |
| R-06 | T-60, CC-02 | 217 "a median 63\% ... (SI)" | Keep the number (`results/allocation_risk/surfaces10/conditional_summary.csv`: 0.628 control only, 0.924 both arms); add both medians to the SI (S-02) so the pointer holds || A |
| R-07 | T-40, FCL-060 | 219 heading "SUD application: Incidence-Guided Saturation is strong under control-only spillover" | "SUD case study: the saturation designs had the lowest error under control-only spillover" (Saturation Quadrants did as well) || E: heading "SUD case study: Incidence-Guided Saturation had lower error than SRS under control-only spillover" (no Saturation Quadrants mention; "lowest" is not claimed because Saturation Quadrants is not in the main text) |
| R-08 | T-59, CC-01 | 222 "94.2--94.9\%, compared with 92.7--93.6\%" | "94.1--94.9\%, compared with 92.6--93.6\%" (`a5_stage1_coverage_by_design_regime.csv`, τ = 0, all nine strategies; the old range left out Saturation Quadrants and Balanced Halves) || E: numbers stay (six designs and SRS); reword to "Across the six designs and SRS and both regimes, averaged over years, ..." |
| R-09 | T-48, FCL-059 | 239 Discussion (new short paragraph after the first) | "Incidence-Guided Saturation was not best in every setting. Under both-arms spillover its MSE was 8\% above that of SRS on the grid. Under control-only spillover its advantage narrowed as β grew, and at β = 3 Isolation Buffer had lower MSE. In the case study, its advantage over SRS fell from 22\% with 58 colleges to 5\% with 29." (all existing results; no new analysis) || A |
| R-10 | T-61, CC-03 | 276 "results and cluster treatment allocations for each planning year" | "treatment budgets, regional orderings and results for each planning year, an example assignment" || A |
| R-11 | T-15 | 235 "The baseline and the adjustment were added together" | "The incidence shift and the adjustment were added together" || E: C 235 becomes "Letting incidence shift the simulated outcome (τ = 1) and adding incidence as a covariate in the fitted model increased High Incidence Focus MSE ... The two changes were made together, so their separate contributions are not distinguished." |

## E. Exhibits

| ID | Rows | Current (E line) | Proposed | Decision |
|---|---|---|---|---|
| E-01 | T-33, FCL-037, FCL-038 | 21 row label "Incidence surface" | "Example incidence surface (not a design)"; keep "cluster" (FCL-038 declined: "region" means the four saturation groups) || A |
| E-02 | T-34, FCL-039 | Table 1 NC Checkerboard cell | Same wording as M-22 || A |
| E-03 | T-35, FCL-040, FCL-054 | Table 1 NC Incidence-Guided Saturation and $2\times2$ Blocking cells | Same wording as M-22, shortened || A |
| E-04 | T-44, FCL-055 | 38, 73 "outcome-dependence values ρ" | "values of the spatial dependence coefficient ρ" || A |
| E-05 | T-46, FCL-057 | 55 SRS "N/A" (two cells) | "Reference" || A |
| E-06 | T-37 | 89 "Application mean MSE" | "Case study mean MSE" || A |

## S. Supplementary information (only items this round touches)

| ID | Rows | Change | Decision |
|---|---|---|---|
| S-01 | T-72 | S 57–60 outcome equation: cite `ord_estimation_1975,lesage_introduction_2009`; S 414–418 MC SE formulas: cite `morris_using_2019` || A |
| S-02 | T-60, CC-02 | S 1890–1892 after "outcome noise accounts for much of the estimated variation": add "(a median 63\% of the variance of Incidence-Guided Saturation's estimated risk under control-only spillover and 92\% under both arms)" || A |
| S-03 | T-62, CC-04 | S 818: $n$ → $N$ ("$N$ the number of participating colleges, $N = 58$ in the main case study") || A |
| S-04 | T-67, CC-09 | S 818 "and we test τ = 0 against the same t distribution" → "and we test β = 0 two-sided at level 0.05, rejecting when the 95\% interval excludes zero. Grid power (Section ...) is one-sided at level 0.025." || A |
| S-05 | T-63, CC-05 | m for professionals per college (S 857) → $n_c$; generic metric m in MC SE formulas (S 416–423, 579–606) → $\theta$ || A |
| S-06 | T-65, CC-07 | S 1604 "SD between 0.20 and 0.21", S 1698 "−0.12 to −0.40", S 2348 two-decimal bias values → three decimals from the stored summaries || A |
| S-07 | T-70, T-71 | S 823 `rate` → `\text{rate}`; S 342 "**WZ** = 1/2" → "$(\mathbf{WZ})_i = 1/2$" || A |
| S-08 | T-05 | S ~158: confirm the single Checkerboard = block stratified sampling statement; no change if already there (consistency check: it is) || A |
| S-09 | T-37 | "application" → "case study" in headings and prose || A |

## Not changed (for the record)

- T-36 (FCL-044) and T-47 (FCL-058): already stated in C 239 and 204; reply only.
- Ratios to SRS stay as they are (Decision H).
- CTJ word limit, full SI prose review: after the proposal (Decision C).
