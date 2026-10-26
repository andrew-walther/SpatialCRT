# Talk outline: treatment allocation for spatial cluster randomized trials with spillover

Draft 2 (2026-10-25), built from the corrected Chapter 2
(`review/thesis-chapter/`) and the corrected results in `review/`. It serves two
talks: **Part I of the prelim** (Monday 2026-12-07; about 20 minutes per project,
so about an hour of talk in a 90-minute exam; committee of statisticians plus
Dr. Van Deinse) and **Dr. Lin's conference talk** in December (statisticians).
Each slide is tagged **both** (same slide in both talks), **prelim** (prelim
only) or **conference** (conference only). Prelim and conference notes say how a
shared slide changes.

The visual version is `outline.html` (built from this file by
`tools/build_outline_html.py`; figures in `figures/`). Most figures are the
manuscript's own (3×3 and 3×4 in their corrected versions). Three are new, made
by `tools/make_figures.R` from the same results: the allocation spread (slide
14), the collinearity diagram (slide 15) and the cross-grid summary (B2). Edit
this file, then rebuild. The deck is drafted here and copied into
`bios-dissertation/slides/prelim/part-1/` once a working draft has been reviewed.

## Open questions for Andrew

- **Slides 14 and B2 use a new view of the existing results.** They plot every
  allocation's MSE and count how many random allocations BSS beats. The numbers
  come from the same per-allocation results as the chapter's box plots
  (`figures/bss_rank_among_srs.csv`), but the chapter never states the
  percentages. Keep them on slides only, or add a sentence to the chapter?
- **3×3 and 3×4 stay in backup** (decided 2026-10-22). With 20 minutes there is
  room to move B1 (the corrected 3×4 figure) into the results as one "same
  pattern on larger grids" slide. Move it up?
- **Why BSS error falls when spatial dependence is added** (2×4: 0.262 to 0.165
  at ψ = 0.5) has not been checked. Slide 16 raises it. Look into it before the
  prelim, or answer "not examined"?
- **Speaker notes** are talking points for now. The kickoff prompt asks for a
  full word-for-word script per slide; draft those once the slide list is
  settled?
- **Conference length.** The outline assumes 20 minutes for Dr. Lin's talk too.
  What slot does the conference give?
- **Status line on slide 1:** "Published in BMC Medical Research Methodology
  (2026)". Add "correction submitted" once the memo goes in?

## Timing

| Section | Slides | Minutes |
|---|---|---|
| 1. Motivation and background | 1–6 | 4 |
| 2. Design and model | 7–11 | 5 |
| 3. Simulation results | 12–16 | 6 |
| 4. Application | 17–19 | 2.25 |
| 5. Choosing a design, take-home, limitations | 20–22 | 2.75 |
| Conference close | 23 | (conference only) |
| Backup | B1–B11 | as needed |

About 20 minutes of main slides.

## 1. Motivation and background

### 1. Where treatment is assigned changes how well a cluster trial measures its effect
- Figure: none (Part I divider, then title)
- Slide text: Part I · SpatialCRT: Spillover. Walther, Van Deinse, Lin. Published in BMC Medical Research Methodology (2026).
- Notes: One sentence of purpose. When an intervention is delivered to whole places, the places we treat and the places we leave as controls are neighbors, and that changes what a trial can estimate.
- Time: 0.25
- Carry-over: both
- Conference note: Use a full title slide (title, authors, affiliations, venue) instead of the prelim divider.

### 2. When an intervention is delivered to whole communities, trials randomize clusters, not people
- Figure: none (simple cluster schematic to draw: a region split into areas, some treated)
- Slide text: Community programs, policies and supervision practices are delivered to places. A cluster randomized trial (CRT) assigns whole areas to intervention or control. When the clusters are geographic, they have neighbors.
- Notes: Define a CRT for Dr. Van Deinse in one sentence. Our clusters later are North Carolina judicial districts. Everything that follows comes from the fact that clusters share borders.
- Time: 0.75
- Carry-over: both
- Conference note: Can be cut to a sentence on slide 3 for a statistics audience.

### 3. Treatment can reach untreated neighbors
- Figure: none (icon-style map to draw: treated village, arrow into a neighboring control village)
- Slide text: Bed nets in Ghana lowered child mortality in nearby untreated villages (Binka et al., 1998). School deworming in Kenya helped untreated children near treated schools (Miguel and Kremer, 2004). Spillover makes controls look better, so the measured effect is too small.
- Notes: Spillover breaks the usual assumption that a cluster's outcome ignores its neighbors' assignments.
- Time: 0.75
- Carry-over: both

### 4. Neighboring outcomes are also correlated for reasons unrelated to treatment
- Figure: none (two-panel cartoon to draw: shared shading across neighbors)
- Slide text: Nearby areas share conditions, so their outcomes move together (Tobler's first law). This spatial dependence violates independence and inflates Type I error. It exists with or without treatment.
- Notes: Contrast with slide 3: spillover is caused by the treatment, spatial dependence is background. The model keeps them as separate parameters (ψ and ρ, slide 10).
- Time: 0.75
- Carry-over: both

### 5. Spatial methods correct for these effects after the data are in, not before
- Figure: none (two-column contrast: analysis stage vs design stage)
- Slide text: Analysis stage: spatial covariates, spatial models, distance bands (reviewed by Jarvis et al., 2017). Design stage: which clusters to treat before the trial starts. Little guidance exists.
- Notes: Jarvis et al. classify the post-hoc approaches and set aside buffer or separation designs from their primary review. That is the gap this project fills.
- Time: 1
- Carry-over: both
- Conference note: Expand into a short method/year/author table ending with this study.

### 6. The question: does the allocation chosen before a trial change how well we estimate the effect?
- Figure: none (text slide)
- Slide text: Compare simple random sampling (SRS) with block stratified sampling (BSS). Measure the error of the intervention effect estimate by its mean squared error (MSE). Vary spillover size, where spillover goes, and spatial dependence.
- Notes: State the aim in one breath. Slide 8 shows the two designs.
- Time: 0.5
- Carry-over: both

## 2. Design and model

### 7. Neighbors can share an edge, a corner, or lie within a distance
- Figure: RookQueen2ndContiguityColored.png
- Slide text: Rook neighbors share an edge. Queen neighbors share an edge or a corner. Distance-based neighbors lie within a fixed radius. This study uses rook neighbors.
- Notes: The neighbor definition decides who can receive spillover and who is correlated. Rook keeps the grids simple; Project 2 uses queen as its primary definition.
- Time: 0.75
- Carry-over: both

### 8. Block stratified sampling keeps intervention clusters from sharing an edge
- Figure: GridsCombined.png
- Slide text: SRS: any balanced allocation. BSS: no two intervention clusters share an edge, and no two controls do where the grid allows. On the 2×4 and 3×4 grids that leaves only the two checkerboards.
- Notes: Top row BSS, bottom row SRS, for the three grids. On the 3×3 grid the arms cannot be balanced, so the control condition is relaxed and six allocations qualify.
- Time: 1
- Carry-over: both

### 9. Spillover may reach only control clusters, or both arms
- Figure: TrtSpillComparisonV2.png
- Slide text: Control only: spillover flows from intervention clusters into neighboring controls. Both arms: it also reaches neighboring intervention clusters. Exposure is binary: a cluster either borders a treated cluster or it does not.
- Notes: This distinction drives every result that follows, so name the two regimes the same way on every results slide.
- Time: 0.75
- Carry-over: both

### 10. A spatial lag model separates the intervention effect, spillover and spatial dependence
- Figure: none (annotated equation)
- Slide text: y_i = α + β x_i + ψ z_i + ρ Σ_j w_ij y_j + ε_i. β: intervention effect (the target). ψ: spillover, with z_i = 1 if cluster i is exposed. ρ: spatial dependence through the row-standardized rook weights w_ij. Fitted by maximum likelihood.
- Notes: Gloss each symbol in words. The model includes the true spillover term, so the comparison is a best case for estimation; only the allocation varies.
- Time: 1.25
- Carry-over: both
- Prelim note: Keep to the one equation; the reduced form and likelihood go to backup B4.
- Conference note: Show the matrix form and the reduced form y = (I − ρW)⁻¹(α + Xβ + Zψ + ε) on the slide.

### 11. We simulate every possible allocation on three small grids
- Figure: SampledPoints.png
- Slide text: 2×4: 70 allocations, 2 BSS. 3×3: 126, 6 BSS. 3×4: 924, 2 BSS. β = 1, ψ from 0.5 to 0.8, ρ = 0 or 0.01, two spillover regimes. 20 subjects per cluster; 50 replications per allocation (10 on 3×4). Each design is judged by the mean of its allocations' MSE and by how much that MSE varies across allocations.
- Notes: The figure shows simulated subjects on the 2×4 grid; only cluster membership enters the fit. Enumerating every allocation shows the whole range of error a trial could face, not just its average. MSE = bias² + variance.
- Time: 1.25
- Carry-over: both

## 3. Simulation results

### 12. With spillover into controls, block stratified sampling wins when outcomes are spatially correlated
- Figure: MSE_Estimates_2x4.png (highlight the top row)
- Slide text: Spillover into controls, ρ = 0.01: BSS 0.165 vs SRS 0.304 (ψ = 0.5). Without spatial dependence SRS is lower on average (0.173 vs 0.262). SRS error varies far more across allocations (SD of MSE 0.41 vs 0.03 at ρ = 0.01).
- Notes: 2×4 grid, the primary results. Read the top row of the figure. After the post-publication correction, the 3×3 and 3×4 grids show the same pattern (3×4 in all 16 scenarios, 3×3 in 12 of 16; backup B1).
- Time: 1.25
- Carry-over: both

### 13. With spillover into both arms, simple random sampling has the lower average error
- Figure: MSE_Estimates_2x4.png (highlight the bottom row)
- Slide text: SRS lower in every scenario: 0.095 vs 0.165 at ρ = 0.01, 0.008 vs 0.262 at ρ = 0 (ψ = 0.5). The same holds on the 3×3 and 3×4 grids.
- Notes: Same figure, bottom row. SRS still varies more across allocations; the reason BSS does worse here is on slide 15.
- Time: 1
- Carry-over: both

### 14. Randomizing rolls the dice: block stratified sampling beats most random draws when spillover reaches controls
- Figure: fig_allocation_spread_2x4.png
- Slide text: Each gray dot is one random allocation. With spillover into controls and spatial dependence, BSS beats 57% of them, and the worst random allocation has an MSE of 2.39, against at most 0.19 for BSS. When spillover reaches both arms, most random allocations do better.
- Notes: This is the risk argument. A trial runs one allocation. Under the regime we expect for probation, BSS is better than a typical draw and never in the long tail. In the bottom rows BSS sits near the top of the error range, for the reason on the next slide. Similar counts on the other grids: 55–57% at ψ = 0.5 under spillover into controls with dependence (backup B2).
- Time: 1.5
- Carry-over: both

### 15. In a checkerboard, spillover exposure is exactly the opposite of treatment
- Figure: fig_collinearity_grid.png
- Slide text: Every control in a checkerboard borders an intervention cluster, so z = 1 − x. The model cannot separate β from ψ, and the estimate absorbs the spillover: bias near −ψ. Random allocations that leave some controls unexposed avoid this.
- Notes: The left panel is the BSS allocation on 2×4; the right is a random allocation with unexposed controls. This is why BSS results are identical under the two spillover regimes on 2×4 and 3×4, and why the 2×4 BSS bias is −0.35 to −0.80. On slide 14's rows without spatial dependence, the cluster of random allocations near 0.25 is the collinear ones (ψ² = 0.25 at ψ = 0.5). Never call this a "paradox"; the published "checkerboard paradox" came from a coding error and is retired.
- Time: 1.25
- Carry-over: both
- Conference note: Add the algebra: with z = 1 − x the design matrix loses rank, so β̂ estimates β − ψ.

### 16. Larger spillover raises error for both designs; spatial dependence raises it for SRS
- Figure: MSE_Estimates_2x4.png (point to the ψ colors and the two columns)
- Slide text: MSE rises with ψ in every setting. Adding spatial dependence (ρ = 0.01) raised SRS error in every scenario. For BSS it lowered error on the 2×4 and 3×4 checkerboards (0.262 to 0.165 on 2×4) and raised it on 3×3.
- Notes: The ψ trend is the expected one: more spillover, more to confuse with the treatment effect. The BSS response to ρ was not examined in the paper; say so if asked.
- Time: 1
- Carry-over: both

## 4. Application

### 17. North Carolina probation trials assign judicial districts, a natural set of spatial clusters
- Figure: NC_DCS_Map_NoCounties.png
- Slide text: 30 judicial districts in 4 judicial divisions. Case 1: trauma-informed and engaged supervision, assigned by district.
- Notes: Probation and Parole Officers supervise people by county and are managed by Chief Probation/Parole Officers and Judicial District Managers. A manager can oversee both intervention and control counties, so spillover into control areas is plausible.
- Time: 0.75
- Carry-over: both

### 18. A second study crosses the four LME/MCO coverage areas
- Figure: NC_MCO_Map.png
- Slide text: Case 2: an implementation strategy for specialty mental health probation. 4 LME/MCOs, one serving each county. Providers can work across counties and MCOs, so spillover is likely.
- Notes: MCO-level policies and provider networks are a confounder, which is why Case 2 stratifies by MCO area.
- Time: 0.75
- Carry-over: both

### 19. When spillover flows into controls, we recommend a best-effort block stratified allocation
- Figure: BSS_Assignments_NC_Map_Applications_Text_Labels.png
- Slide text: Stratify by judicial division (Case 1) or MCO area (Case 2) and keep same-arm districts apart where the map allows. Real borders are irregular, so perfect separation is not always possible.
- Notes: The recommendation is conditional on the expected spillover mechanism. Document how closely the chosen allocation approximates the no-shared-edge rule.
- Time: 0.75
- Carry-over: both

## 5. Choosing a design, take-home and limitations

### 20. Match the design to the study's expected spillover and spatial correlation
- Figure: none (2×2 decision table: spillover into controls only vs both arms, by spatial dependence present vs absent)
- Slide text: Controls only, spatially correlated outcomes: BSS (lower MSE on every grid, e.g. 0.165 vs 0.304). Controls only, no dependence: SRS on average on 2×4 and 3×4, BSS on 3×3; check allocations for collinearity. Both arms, either case: SRS (e.g. 0.095 vs 0.165), or a buffer-zone design.
- Notes: This is the practical guidance: investigators choose by the characteristics of their study. All cells use ψ = 0.5 values; the pattern holds across ψ.
- Time: 1
- Carry-over: both

### 21. Block stratified sampling gives up some average accuracy for protection against a bad allocation
- Figure: none (take-home layout: claim title, highlight box, three bolded bullets)
- Slide text: Highlight box: choose the allocation for the spillover you expect. Bullets: **Risk vs average:** SRS can win on average but can also draw a very poor allocation. **Regime matters:** BSS helps most when spillover flows into controls and outcomes are spatially correlated. **State the mechanism:** justify the design by the anticipated spillover direction.
- Notes: One sentence per bullet. Close by saying where the project stands: published, with a correction to the 3×3 and 3×4 summaries in progress that leaves the conclusions unchanged.
- Time: 1
- Carry-over: both

### 22. These results come from small regular grids with a simple spillover model
- Figure: none (two columns: limitations, future work)
- Slide text: Limitations: treatment and spillover can be collinear; binary rook spillover; small regular grids; no intra-cluster correlation; results specific to the spatial lag model. Future work: distance-decaying and additive spillover, buffer-zone designs, irregular geographies, model misspecification.
- Notes: The collinearity limitation is the important one (slide 15). Buffer zones restore identifiability but need grids large enough to drop boundary controls (backup B8).
- Time: 0.75
- Carry-over: both
- Prelim note: End with the bridge to Part II: larger grids, heterogeneous incidence, queen neighbors and more designs.

### 23. Thank you
- Figure: none
- Slide text: Co-authors Tonya Van Deinse and Feng-Chang Lin. Funding as in the article's Declarations. Questions.
- Notes: Conference only; the prelim has one shared thank-you slide at the end.
- Carry-over: conference

## Backup

### B1. After the correction, the 3×3 and 3×4 grids follow the 2×4 pattern
- Figure: MSE_Estimates_3x4_corrected.png
- Slide text: 3×4 matches 2×4 in all 16 scenarios. 3×3 matches in 12 of 16; the exception is spillover into controls without spatial dependence, where BSS is lower. MSE_Estimates_3x3_corrected.png is the 3×3 version.
- Notes: Use when asked whether the 2×4 results generalize.
- Carry-over: both

### B2. Where block stratified sampling falls among random allocations, every grid and setting
- Figure: fig_bss_rank_all_grids.png
- Slide text: Share of random allocations with higher MSE than BSS, for each grid, spillover regime, ρ and ψ. Spillover into controls with dependence: 40–60%. Spillover into both arms: 0–29%.
- Notes: Use when asked whether slide 14 is specific to one grid or one ψ. Each point is one ψ value.
- Carry-over: both

### B3. Full result tables
- Figure: none (re-typeset the three tables from the corrected chapter)
- Slide text: MSE, SD of MSE across allocations and mean bias for SRS and BSS, every scenario, all three grids.
- Notes: Use when asked for exact values. The SD column is the spread of MSE across allocations, not the SD of β̂ across replications.
- Carry-over: both

### B4. Reduced form and likelihood
- Figure: none (equations)
- Slide text: y = (I − ρW)⁻¹(α + Xβ + Zψ + ε). Outcomes were generated from the reduced form and fitted by maximum likelihood.
- Notes: Use when asked how outcomes were generated and how the model is estimated.
- Carry-over: both

### B5. Spatial weights matrix
- Figure: none (small grid with its row-standardized W, to draw)
- Slide text: w_ij = 1/|N(i)| for each rook neighbor j of cluster i, 0 otherwise, so Wy is the average outcome of a cluster's neighbors.
- Notes: Use when asked how W is defined and why it is row-standardized.
- Carry-over: both

### B6. Spillover and spatial dependence are different parameters
- Figure: none
- Slide text: ψ is caused by the treatment crossing a boundary; ρ is background correlation that exists without any treatment. Both use the same neighbor structure, so they are hard to separate when both are large.
- Notes: Use when asked "aren't spillover and spatial autocorrelation the same thing?"
- Carry-over: both

### B7. BSS has too few allocations for randomization inference, so inference is model-based
- Figure: none
- Slide text: BSS admits 2–6 allocations, so no permutation distribution. Maximum likelihood conditions on the realized allocation and gives a valid sampling distribution either way.
- Notes: Use when asked "is BSS still randomized, and how do you do inference?" Also covers how BSS differs from covariate-constrained randomization (spatial isolation vs covariate balance).
- Carry-over: both

### B8. Buffer-zone designs could restore identifiability
- Figure: none (diagram to draw: a treated cluster, a dropped ring, distant controls)
- Slide text: Exclude controls that border intervention clusters. Needs grids large enough that interior controls remain.
- Notes: Use when asked what to do about the collinearity problem.
- Carry-over: both

### B9. The coefficients are structural, not causal estimands
- Figure: none
- Slide text: β and ψ are spatial lag model coefficients. A potential-outcome contrast is a function of all parameters and the spatial multiplier.
- Notes: Use when asked which causal estimand β corresponds to.
- Carry-over: both

### B10. Simulation algorithm, step by step
- Figure: none (flowchart to draw from Appendix A, Simulation Algorithm)
- Slide text: Define the scenario, enumerate allocations, build the grid and rook neighbors, flag BSS allocations, simulate and fit each allocation N times, aggregate, summarize by design.
- Notes: Use when asked how many replications were run and how BSS allocations were identified.
- Carry-over: both

### B11. A post-publication check corrected the 3×3 and 3×4 summaries
- Figure: none
- Slide text: Two coding errors: the summary used the 2×4 BSS allocations for every grid, and a spillover loop stopped at district 8. The 2×4 results were unaffected, and the conclusions held.
- Notes: Use if a committee member has read the published tables. Update once the correction is submitted.
- Carry-over: prelim
