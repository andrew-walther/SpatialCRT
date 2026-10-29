# Talk outline: sampling design for spatial cluster randomized trials with spillover

Draft 3 (2026-10-28). The deck is built (`spillover_talk.qmd`), and the slide entries
below are generated from it by `tools/outline_from_deck.py`, so they always describe
the deck as rendered. Edit the deck, re-render, then regenerate this file and
`outline.html`. Everything above the marker line is written by hand.

The deck serves two talks: **Part I of the prelim** (Monday 2026-12-07; about 20
minutes per part; committee of four statisticians and Dr. Van Deinse) and **Dr. Lin's
conference talk** in December (statisticians). Tags: **both** (same slide in both
talks), **prelim** (prelim only), **conference** (standalone or conference only; the
prelim drops it). Prelim and conference notes say how a slide changes.

## Decisions (author, 2026-10-28)

- The corrected 3×3 and 3×4 results are simply the results. The slides never mention
  the correction, and there is no correction backup slide or status line.
- Framing: BSS is a good design, not the best or the worst, against rolling the dice on
  one SRS allocation. The results lead with where BSS falls among the random
  allocations, and the deck says plainly that most random allocations beat BSS when
  spillover reaches both arms (checkerboard collinearity). Per-scenario "optimal
  allocation" counts are left out. "Checkerboard paradox" is not used.
- The BSS-rank percentages appear on slides only; Chapter 2 is not changed.
- The 3×3 and 3×4 grids get one main slide (both grids side by side); their full
  figures and tables are in backup.
- Background carries the prelim's shared motivation (trials of places, spillover,
  spatial dependence) and Part I's own prior work and gap (analysis-stage methods vs
  design-stage options, from the literature review). No standalone literature review.
- Notation is Part I's own (β intervention effect, ψ spillover, ρ spatial dependence,
  cluster-level outcome). Part II defines its own; no harmonizing, no bridge slide.
- The results are presented as published. What ρ = 0.01 means on the fitted scale
  (about 0.5) is explained in a backup slide; see
  `../review/notes/2026-10-28-rho-scale-note.md`.
- The paper reference is on the last Discussion slide. The standalone deck ends with a
  thank-you slide (Dr. Lin, Dr. Van Deinse; NIMH K01MH129619 and BJA
  15PBJA-24-GG-01896-PRJH). In the prelim these acknowledgements move to the single
  closing slide after Part III.
- Every number on the slides and in the script is computed at render time from
  `results/` and `review/` (`tools/load_results.R`), not typed by hand.

## Timing

| Section | Main slides | Minutes |
|---|---|---|
| Standalone frame (title, outline) | 2 | 0.5 |
| Background | 5 | 4.5 |
| Methods | 5 | 4.5 |
| Simulation results | 5 | 5.5 |
| Application | 3 | 3 |
| Discussion | 3 | 2.5 |
| Thank you (standalone only) | 1 | 0.25 |
| Backup | 15 | as needed |

20 minutes on the slides the prelim keeps, 20.75 with the standalone frame and thank-you slide. The script's
`[time | running]` markers give the per-slide timing.

## Shorter conference slots

The conference slot is not known yet; the deck assumes about 20 minutes.

- **15–20 minutes:** merge "NC MCO areas" into "NC probation districts", and say the
  neighbor definition in one sentence instead of showing "Neighbors". Saves about
  1.25 minutes.
- **10–15 minutes:** also cut "Trials of places" (one sentence on "Spillover"), fold
  "Spatial dependence" into "Spillover", reduce "The design gap" to its bottom card, and
  cut "Control & intervention clusters" and "Larger grids" (one sentence each on
  "Choosing a design"). Keep the model, "Where BSS falls", "Why: collinearity",
  "Choosing a design" and the take-home in every length. Saves about 5.5 minutes.

## Decided after review rounds 1–2 (author, 2026-10-28)

- **NC example map:** keep the published map. The script presents it as a suggested,
  best-effort allocation for the state's irregular geography and says further
  restrictions could guarantee that no intervention district borders another.
- **Why spillover into controls:** the proposed studies assume a consistent treatment
  effect, so treated districts gain nothing more from other treated districts, while
  nearby controls can pick the program up. In Case 2, a trained provider who serves
  several counties carries it into control counties too.
- **Who BSS beats:** framed as downside risk chosen before the trial. Many random draws
  that keep β identified do better than a checkerboard, but which one can only be known
  after the fact, and the worst draw on every grid passes the same check (backup
  "Downside risk").
- **Scale of ρ:** stays in backup, with a prepared answer in its notes about the
  subject-level generating weights and the cluster-level model in the paper.
- **ψ is not the target:** the project estimates β; that ψ is not estimable under
  the checkerboards matters only through β̂ ≈ β − ψ.

- **Refining the NC map:** the Recommendation slide now says a refined allocation could
  also leave a few control districts with no intervention neighbor (author approved
  rewording, 2026-10-28).
- **Larger grids:** the 3×3/3×4 slide shows every allocation (boxes) and each BSS
  allocation (diamonds), so the spread of random draws is visible.

## Open questions for Andrew

<!-- SLIDES (generated) -->

## 1. Background

### 1. Title
- Figure: none
- Slide text: Sampling Design for Spatial Cluster Randomized Trials with Spillover and Spatial Dependence Andrew J Walther Department of Biostatistics, University of North Carolina at Chapel Hill PhD Preliminary Oral Exam · December 7, 2026
- Notes: Thank you all for being here. This first part is about a design question that comes up whenever a trial treats whole places, namely which places to treat when the treatment can reach the neighbors.
- Time: 0.25
- Carry-over: conference
- Prelim note: Dropped; the prelim has its own title slide and opens Part I with its divider.

### 2. Outline
- Figure: none
- Slide text: • Background Cluster trials of places, spillover, spatial dependence, and what design guidance exists • Methods Neighbors, the two allocation designs, the two spillover types, the spatial lag model and the simulation • Simulation Results Where block stratified sampling falls among random allocations, and why • Application Probation trials across North Carolina judicial districts • Discussion Choosing a design, …
- Notes: I'll start with why geography complicates a cluster trial, then the design and model, the simulation results, an application to probation trials in North Carolina, and what I'd recommend.
- Time: 0.25
- Carry-over: conference
- Prelim note: Dropped; the prelim has one agenda for all three parts.

### 3. When a program is delivered to whole communities, trials randomize places, and places have neighbors
- Figure: none
- Slide text: A cluster randomized trial (CRT) assigns whole groups, here geographic areas, to intervention or control It is the natural design when a program, policy or supervision practice is delivered to an area rather than to one person Geographic clusters share borders, so what happens in one area can depend on its neighbors Two separate mechanisms create that dependence, spillover and spatial dependence, and both affect …
- Notes: When a program is delivered to a whole community, such as a supervision practice in a probation office, we can't randomize individual people, so a cluster randomized trial assigns whole areas to intervention or control. When those clusters are geographic, they have neighbors.
- Time: 1
- Carry-over: both
- Prelim note: Serves as the prelim's shared motivation for cluster trials of places. Part II builds on it.
- Conference note: Cut first for a statistics audience in a 10–15 minute slot (keep one sentence on the next slide).

### 4. Treatment can spill over into neighboring untreated areas
- Figure: none
- Slide text: Bed nets, northern Ghana Child mortality was lower in control villages near villages given insecticide-treated bed nets than in more distant control villages. School deworming, western Kenya Untreated children attending schools near treated schools had fewer worm infections. Spillover is a treatment effect that crosses a cluster boundary A standard CRT analysis assumes a cluster's outcome depends only on its own …
- Notes: Two classic trials show what spillover looks like. In northern Ghana, child mortality was lower in control villages near villages that received treated bed nets than in control villages farther away.
- Time: 0.75
- Carry-over: both
- Prelim note: Shared motivation; Part II refers back to it rather than repeating the examples.

### 5. Neighboring outcomes are also correlated for reasons that have nothing to do with treatment
- Figure: none
- Slide text: Spatial dependence ($\rho$) Neighbors have similar outcomes because they share conditions, with or without treatment Spillover ($\psi$) The treated area (dark) raises outcomes in its neighbors (amber); without treatment there is none Both run through the same neighbors, so they are easy to confuse, and the model keeps them as separate parameters Spatial dependence breaks the independence a standard analysis …
- Notes: The second mechanism is spatial dependence. On the left, neighboring areas have similar outcomes because they share demographics, labor markets, providers and local practice, whether or not anyone is treated.
- Time: 0.75
- Carry-over: both
- Prelim note: Shared motivation; Part II refers back to it.
- Conference note: In a 10–15 minute slot, fold into one sentence on "Spillover".

### 6. Spatial methods mostly adjust after the data are in, and guidance on which clusters to treat is thin
- Figure: none
- Slide text: Before the trial (design) Choose which clusters to treat Restricted randomization balances covariates, and recent versions spread the selected clusters apart Buffer zones drop areas near the boundary between arms Saturation designs vary the treated share across groups → During the trial Treatment reaches neighbors, and neighboring outcomes move together Spillover ($\psi$) Spatial dependence ($\rho$) → After the …
- Notes: A spatial cluster trial has three stages, and most spatial methods work in the last one, adjusting for exposure or modeling space after the data are in. At the design stage there are restricted randomization, buffer zones and saturation designs, but few studies compare simple allocation rules when spillover and spatial dependence occur together.
- Time: 0.75
- Carry-over: both
- Prelim note: Part I's own prior work and gap. Part II adds heterogeneous incidence to it.
- Conference note: In a 15–20 minute slot, keep; in a 10–15 minute slot, reduce to the bottom card.

### 7. Does the allocation chosen before the trial change how well we estimate the effect?
- Figure: fig_designs.png
- Slide text: Dark cells are intervention clusters. The three grids studied here. Simple random sampling (SRS) Any allocation with half the clusters treated (4 of 9 on the 3×3 grid), drawn at random. It can put all treated clusters together, as in the bottom row. Block stratified sampling (BSS) No two intervention clusters share an edge, and no two controls where the grid allows. That leaves the two checkerboards on the 2×4 …
- Notes: So the question is whether the allocation chosen before the trial, meaning the list of which clusters get the intervention, changes how well we can estimate the effect. Here are the two designs on the three grids we study, with block stratified allocations on top and one random allocation below.
- Time: 1.25
- Carry-over: both


## 2. Methods

### 8. This study defines neighbors as clusters that share an edge
- Figure: RookQueen2ndContiguityColored.png
- Slide text: Rook neighbors share an edge, queen neighbors share an edge or a corner, and higher-order neighbors are neighbors of neighbors Distance-based neighbors lie within a fixed radius We use rook neighbors throughout, the same edge rule that defines BSS, and they decide both who can receive spillover and who is spatially correlated
- Notes: Everything depends on how we define a neighbor. Rook neighbors share an edge, queen neighbors share an edge or a corner, and you can also define neighbors by distance or as neighbors of neighbors.
- Time: 0.5
- Carry-over: both
- Conference note: In a 10–15 minute slot, cut and say "rook neighbors" on the designs slide.

### 9. Spillover may reach only control clusters, or control and intervention clusters
- Figure: none
- Slide text: Control Clusters Only Control & Intervention Clusters Intervention cluster ($x_i = 1$) Control cluster ($x_i = 0$) Exposed to spillover ($z_i = 1$), unmarked $z_i = 0$ Exposure is binary. A cluster is exposed if at least one rook neighbor is an intervention cluster (and, for Control Clusters Only, it is itself a control), however many treated neighbors it has Control Clusters Only is classic contamination; under …
- Notes: Spillover can travel in two ways, and this turns out to drive every result. Dark squares are intervention clusters, and here the left half of a two-by-four grid is treated.
- Time: 1
- Carry-over: both

### 10. A spatial lag model separates the intervention effect, spillover and spatial dependence
- Figure: none
- Slide text: $$ y_i \;=\; \alpha \;+\; \beta\, x_i \;+\; \psi\, z_i \;+\; \rho \sum_{j \in \mathcal{N}(i)} w_{ij}\, y_j \;+\; \varepsilon_i, \qquad \varepsilon_i \overset{iid}{\sim} N(0, \sigma^2) $$ Control Clusters Only $z_i = 1$ if $x_i = 0$ and $x_j = 1$ for some $j \in \mathcal{N}(i)$ Control & Intervention Clusters $z_i = 1$ if $x_j = 1$ for some $j \in \mathcal{N}(i)$ $y_i$ Cluster-level outcome of cluster $i$ $\alpha$ …
- Notes: This is the spatial lag model behind every result, and every symbol is defined below. A cluster's outcome is a baseline alpha, plus the intervention effect beta if the cluster is treated, plus the spillover effect psi if it is exposed.
- Time: 1.25
- Carry-over: both
- Prelim note: Part I defines its own model and notation in full; Part II defines its own (no harmonizing).
- Conference note: For a statistics audience, can add the reduced form from backup.

### 11. Outcomes are generated with y on both sides, so effects ripple through the neighbors
- Figure: none
- Slide text: $$ \mathbf{y} = \rho\mathbf{W}\mathbf{y} + \boldsymbol\alpha + \mathbf{X}\beta + \mathbf{Z}\psi + \boldsymbol\varepsilon \quad\Longrightarrow\quad \mathbf{y} = (\mathbf{I} - \rho\mathbf{W})^{-1}\left(\boldsymbol\alpha + \mathbf{X}\beta + \mathbf{Z}\psi + \boldsymbol\varepsilon\right) $$ $$ (\mathbf{I} - \rho\mathbf{W})^{-1} = \mathbf{I} + \rho\mathbf{W} + \rho^2\mathbf{W}^2 + \cdots $$ The matrix form stacks the …
- Notes: The outcome appears on both sides of the model, so we can't draw it directly, and solving for it gives the reduced form on the right. Expanding the inverse shows what it does, passing each cluster's outcome to its neighbors, then to their neighbors, and so on, with rho controlling how fast it fades.
- Time: 0.5
- Carry-over: both

### 12. We enumerate every possible allocation on three small grids and compare their error
- Figure: SampledPoints_trim.png
- Slide text: One replication on the 2×4 grid, 160 subjects placed at random. Only cluster membership enters the fit. | Grid | Clusters | Allocations (SRS) | BSS allocations | Replications | | 2×4 | 8 | 70 | 2 | 50 | | 3×3 | 9 | 126 | 6 | 50 | | 3×4 | 12 | 924 | 2 | 10 | 16 scenarios per grid: $\psi \in \{0.5, 0.6, 0.7, 0.8\}$ × $\rho \in \{0, 0.01\}$ × two spillover types Outcomes generated from the model's reduced form, …
- Notes: Small grids let us enumerate every possible allocation, from 70 on the two-by-four grid to 924 on the three-by-four grid. Each grid has sixteen scenarios, four spillover sizes, two levels of spatial dependence and the two spillover types.
- Time: 1.25
- Carry-over: both


## 3. Simulation Results

### 13. With spillover into controls and correlated outcomes, block stratified sampling beats about half of random draws and avoids the worst ones
- Figure: fig_spread_2x4.png
- Slide text: Each dot is one allocation's MSE over its replications (2×4 grid, ψ = 0.5). Further left is better. Each row has its own scale. 57% of random allocations have higher MSE than BSS with spillover into controls and $\rho = 0.01$ (2×4, $\psi = 0.5$) The worst tenth of random draws average MSE 1.20, and the worst is 2.39; neither BSS allocation exceeds 0.19 Across all grids and $\psi$ in that setting, BSS beats 40–60% …
- Notes: I want to lead with the view that matters most for a trial, which runs one allocation. Each grey dot is one random allocation on the two-by-four grid with psi at one-half, and the diamonds are the two block stratified allocations.
- Time: 1.75
- Carry-over: both

### 14. With spillover into controls only, block stratified sampling has lower error when outcomes are spatially correlated
- Figure: fig_mse_2x4_controls.png
- Slide text: With ρ = 0.01, BSS 0.165 vs SRS 0.304 ($\psi = 0.5$), lower at every $\psi$, with the gap closing by $\psi = 0.8$ (0.569 vs 0.573) With ρ = 0, SRS is lower on average, 0.173 vs 0.262 Spread across allocations, SD of MSE 0.41 for random draws; BSS is effectively one allocation here (two mirror-image checkerboards)
- Notes: Now the averages, for spillover into controls only. Boxes show the random allocations and diamonds the block stratified mean, across the four spillover sizes.
- Time: 1
- Carry-over: both

### 15. With spillover into both arms, simple random sampling has the lower average error
- Figure: fig_mse_2x4_both.png
- Slide text: SRS is lower in all 24 of 24 such scenarios across the three grids, for example 0.095 vs 0.165 with $\rho = 0.01$ and 0.008 vs 0.262 with $\rho = 0$ (2×4, $\psi = 0.5$) BSS results are identical under both spillover types on every grid, since no treated cluster has a treated neighbor
- Notes: When spillover reaches both arms, the picture flips. Simple random sampling has the lower average error in all 24 of these scenarios across the three grids, for example 0.095 against 0.165 with spatial dependence on the two-by-four grid.
- Time: 0.5
- Carry-over: both
- Conference note: In a 10–15 minute slot, cut; "Where BSS falls" and "Why: collinearity" carry the point.

### 16. With spillover into controls only, an allocation that exposes every control makes exposure the exact opposite of treatment
- Figure: fig_collinearity_grid.png
- Slide text: Under Control Clusters Only, if every control is exposed then $z_i = 1 - x_i$ for every cluster, $[\mathbf{1}, x, z]$ loses rank, and $\beta$ and $\psi$ cannot be separated. `lagsarlm` drops $z$, and $\hat\beta$ estimates $\beta - \psi$ At $\rho = 0$ the 2×4 BSS bias is $-\psi$ (-0.50 to -0.80); with spatial dependence the lag term absorbs part of it (-0.35 to -0.72) Control Clusters Only. Every checkerboard and …
- Notes: Here is the reason. With spillover into controls only, whenever every control borders an intervention cluster, the spillover indicator is exactly one minus the treatment indicator, so the model can't tell the intervention effect from the spillover effect.
- Time: 1.25
- Carry-over: both

### 17. The same pattern holds on the 3×3 and 3×4 grids
- Figure: fig_mse_larger_grids.png
- Slide text: The better design matches the 2×4 grid in 16 of 16 scenarios on 3×4 and 12 of 16 on 3×3. The exception (3×3, spillover into controls, $\rho = 0$) has five BSS allocations that leave a control unexposed Random draws spread widely. On 3×4 with spillover into controls and $\rho = 0.01$ they reach MSE 2.57, while both BSS allocations stay at or below 0.18 ($\psi = 0.5$)
- Notes: The larger grids tell the same story. On the three-by-four grid, which also has only the two checkerboards, the better design matches the two-by-four grid in all 16 scenarios.
- Time: 1
- Carry-over: both
- Conference note: In a 10–15 minute slot, cut to one sentence on "Choosing a design".


## 4. Application

### 18. North Carolina probation trials assign counties or judicial districts, natural spatial clusters
- Figure: NC_DCS_Map_trim.png
- Slide text: 30 judicial districts in 4 divisions (NC Division of Community Supervision) Chief probation/parole officers and district managers oversee districts, often spanning several counties Case 1, trauma-informed supervision. Training for officers and chiefs on specialty caseloads, randomized by county or district Neighbors share courts, labor markets and providers, and leaders and providers can carry the program across arms
- Notes: These questions came from two probation studies with the North Carolina Department of Adult Correction. Supervision is organized into thirty judicial districts in four divisions, and chiefs and district managers often span several units and counties.
- Time: 1
- Carry-over: both

### 19. A second study crosses the four behavioral health management areas
- Figure: NC_MCO_Map_trim.png
- Slide text: Each county shaded by the LME/MCO that manages its behavioral health services. Case 2, specialty mental health probation. Cross-training and joint treatment-team meetings between officers and behavioral health providers Each county is served by one of 4 LME/MCOs (Local Management Entity/Managed Care Organizations) Large providers work across counties and MCO areas, so practices can travel into control counties …
- Notes: The second study tests a strategy for specialty mental health probation, built on cross-training and joint team meetings between officers and behavioral health providers. Behavioral health care is managed by four regional managed care organizations, and each county is served by one of them.
- Time: 0.75
- Carry-over: both
- Conference note: In a 15–20 minute slot, can merge into "NC probation districts"; cut in a 10–15 minute slot.

### 20. For both studies we recommend a best-effort block stratified allocation
- Figure: BSS_Assignments_NC_Map_Applications_Text_Labels.png
- Slide text: An example allocation of judicial districts (Case 1). Red, intervention; blue, control. Case 2 would do the same with counties within each MCO area. Stratify by judicial division (Case 1) or MCO area (Case 2), then keep same-arm districts apart where the map allows Why controls only? With a consistent treatment effect, treated districts gain nothing more from other treated districts, but nearby controls can pick …
- Notes: For both studies we recommend block stratified sampling, adapted to the map, with intervention in red and control in blue. The recommendation rests on the proposed studies' assumption that, with a consistent treatment effect, treated districts gain nothing more from other treated districts, while nearby control districts can pick the program up.
- Time: 1.25
- Carry-over: both


## 5. Discussion

### 21. Match the design to the spillover and spatial correlation you expect
- Figure: none
- Slide text: [table] 2×4 values at ψ = 0.5; the direction holds across ψ. For probation districts, which share courts, providers and labor markets, the first column is the realistic one.
- Notes: Putting the results together, the design should follow from what investigators expect about spillover and spatial correlation. If spillover flows into controls and neighboring outcomes move together, which is the realistic case for probation districts, block stratified sampling had the lower error and avoided the worst draws.
- Time: 1
- Carry-over: both

### 22. Block stratified sampling is a good design, not the best or the worst, against rolling the dice on one random draw
- Figure: none
- Slide text: Choose the allocation for the spillover you expect, and know what a single random draw can cost Risk, not just the average SRS can win on average, but one trial runs one allocation, and some random draws estimate the effect poorly. The spillover type matters BSS helps most when spillover flows into controls and outcomes are spatially correlated, and most random allocations beat it when spillover reaches both …
- Notes: The take-home is that block stratified sampling is a good design, not the best and not the worst, and the comparison that matters is with rolling the dice on a single random draw. Random allocation can win on average, but a trial runs one allocation, and some random draws estimate the effect poorly.
- Time: 1
- Carry-over: both
- Conference note: Keep in every length.

### 23. These results come from small regular grids and a simple spillover model
- Figure: none
- Slide text: Limitations Treatment and spillover are exactly collinear whenever every control is exposed, which includes every checkerboard Binary rook spillover, with no decay over distance and no added effect from more treated neighbors Small regular grids of equal-area clusters No intra-cluster correlation; results specific to the spatial lag model with the spillover term included Next steps Distance-decaying and additive …
- Notes: The most important limitation is the collinearity, which arises whenever every control is exposed. Spillover is binary and limited to rook neighbors, the grids are small and regular, and the model has no intra-cluster correlation and is fitted with the spillover term included.
- Time: 0.5
- Carry-over: both
- Prelim note: Keeps the paper reference. Add one closing sentence pointing to Part II (larger grids, heterogeneous incidence, queen neighbors, more designs).

### 24. Thank you!
- Figure: none
- Slide text: Advisor Feng-Chang Lin Collaborator Tonya Van Deinse Funding National Institute of Mental Health, K01MH129619 (Dr. Van Deinse) Bureau of Justice Assistance, 15PBJA-24-GG-01896-PRJH (NC Department of Adult Correction) Paper Walther et al., *BMC Medical Research Methodology* (2026) Questions welcome awalther@unc.edu
- Notes: Thank you to my advisor, Dr. Lin, and to Dr. Van Deinse, whose probation studies motivated this work. I'm happy to take questions.
- Time: 0.25
- Carry-over: conference
- Prelim note: Dropped; the acknowledgements and funding move to the single prelim closing slide after Part III, grouped with Parts II and III.


## Backup

### B1. 3×3 grid, every allocation
- Figure: MSE_Estimates_3x3_corrected.png
- Slide text: Rows are spillover types (Control Only on top), columns are $\rho$. "Block" is the six BSS allocations and "Random" every allocation, colored by $\psi$ The BSS spread comes from its six different allocations. One, the full checkerboard, is collinear; the other five leave a control unexposed and estimate $\beta$ almost exactly With spillover into controls and $\rho = 0.01$, random allocations reach MSE 1.20 at …
- Notes: Use when asked for the per-allocation distributions on the 3×3 grid. This is the chapter's figure.
- Carry-over: both

### B2. 3×4 grid, every allocation
- Figure: MSE_Estimates_3x4_corrected.png
- Slide text: Rows are spillover types (Control Only on top), columns are $\rho$. "Block" is the two checkerboards and "Random" all 924 allocations, colored by $\psi$ With spillover into controls and $\rho = 0.01$, random allocations reach MSE 2.57 at $\psi = 0.5$, while both checkerboards stay at or below 0.18 With spillover into both arms, most random allocations sit near zero and the checkerboards do not
- Notes: Use when asked for the per-allocation distributions on the 3×4 grid. This is the chapter's figure.
- Carry-over: both

### B3. Where block stratified sampling falls among random allocations, every grid and setting
- Figure: fig_bss_rank_all_grids.png
- Slide text: One point per ψ; points that share a value overlap. With ρ = 0 and spillover into controls, most draws expose every control and estimate β − ψ like BSS, so the share there splits near-ties by simulation noise.
- Notes: Use when asked whether the 57% on the 2×4 grid is special. With spillover into controls and ρ = 0.01, BSS beats 40 to 60 percent on every grid; with spillover into both arms, at most 29%. The share counts random allocations with MSE above the BSS mean; the random set includes the BSS allocations, as in the paper.
- Carry-over: both

### B4. A good random draw can beat the checkerboard, but it cannot be picked before the trial
- Figure: none
- Slide text: Spillover into controls, $\rho = 0.01$, $\psi = 0.5$. The allocation is chosen before any outcome is seen, so the relevant comparison is with the whole spread of random draws 2×4: BSS MSE 0.165; the worst tenth of random draws average 1.20, about 7 times as large (3×3: 0.147 vs 0.84; 3×4: 0.133 vs 0.77) Random draws that leave some control unexposed keep $\beta$ and $\psi$ separate, and many of them beat BSS (on …
- Notes: Use when asked "against the random allocations that keep β identified, how does BSS do?" Answer plainly that many of them do better than a strict checkerboard, but nobody can tell in advance which ones, and the worst random draw passes the same check. BSS is chosen for its bounded downside, not because it is the best allocation.
- Carry-over: both

### B5. Full results, 2×4 grid
- Figure: none
- Slide text: [table] MSE, SD of MSE across allocations, and mean bias of $\hat\beta$. Shading marks the design with lower MSE. SD here is the spread of MSE across allocations, unscaled (the chapter reports it ×10). The two 2×4 checkerboards are mirror images, so their difference (0.144 vs 0.186 at ρ = 0.01, ψ = 0.5) is a guide to the Monte Carlo error.
- Notes: Use when asked for exact values. Generated from review/tables at render time.
- Carry-over: both

### B6. Full results, 3×3 grid
- Figure: none
- Slide text: [table] As for the 2×4 table. The 3×3 grid has six BSS allocations, not all checkerboards, so BSS MSE also varies across allocations.
- Notes: Use when asked for exact values on the 3×3 grid.
- Carry-over: both

### B7. Full results, 3×4 grid
- Figure: none
- Slide text: [table] As for the 2×4 table. 10 replications per allocation on this grid.
- Notes: Use when asked for exact values on the 3×4 grid.
- Carry-over: both

### B8. Outcomes come from the reduced form, and the model is fitted by maximum likelihood
- Figure: none
- Slide text: $$ \mathbf{y} = \rho \mathbf{W}\mathbf{y} + \boldsymbol\alpha + \mathbf{X}\beta + \mathbf{Z}\psi + \boldsymbol\varepsilon \;\Longrightarrow\; \mathbf{y} = (\mathbf{I} - \rho\mathbf{W})^{-1}\left(\boldsymbol\alpha + \mathbf{X}\beta + \mathbf{Z}\psi + \boldsymbol\varepsilon\right), \;\; \boldsymbol\varepsilon \sim N(\mathbf{0}, \sigma^2\mathbf{I}) $$ $$ \ell(\alpha, \beta, \psi, \rho, \sigma^2) = …
- Notes: Use when asked how outcomes were generated and how the model was estimated. The Jacobian term log|I − ρW| is what distinguishes this from least squares on the lagged outcome.
- Carry-over: both

### B9. Row-normalized rook weights make the lag term a neighbor average
- Figure: none
- Slide text: Cells numbered row by row for this example; cluster 2 and its rook neighbors (amber). $$ w_{ij} = \frac{c_{ij}}{\sum_k c_{ik}}, \qquad c_{ij} = \begin{cases} 1 & \text{edge shared} \\ 0 & \text{otherwise} \end{cases} $$ Cluster 2 has rook neighbors 1, 3 and 6, so $w_{21} = w_{23} = w_{26} = 1/3$ Each row sums to 1, so $\sum_j w_{ij} y_j$ is the average neighbor outcome and $\rho$ is comparable across clusters …
- Notes: Use when asked how W is defined and why it is row-standardized. Corner clusters have two neighbors and edge clusters three on the 2×4 grid. For how the simulation's generating weights relate to this, see the next slide.
- Carry-over: both

### B10. How large is ρ = 0.01?
- Figure: none
- Slide text: In the simulation code, outcomes were generated for the subjects with subject-level binary weights, linking each subject to every subject in a rook-neighboring cluster, about 40 to 60 links The model is fitted with the row-normalized version of the same neighbor structure On that fitted scale, $\rho = 0.01$ in the generating process corresponds to roughly $0.01 \times (40\text{ to }60) \approx 0.4$ to $0.6$, and …
- Notes: Use when asked "isn't ρ = 0.01 essentially zero?", why adding it changes the error so much, or how the subject-level and cluster-level weights relate. On the scale of the fitted model it is about one half. If asked about the weights, say plainly that outcomes were simulated for subjects within clusters, with a subject-level weights matrix linking each subject to every subject in the neighboring clusters, and the model was fitted to the subjects. Every subject in a cluster shares x, z and neighbors, so the cluster's signal is the subjects' common mean, which is why the paper writes the model at the cluster level with row-normalized weights. The two scales of ρ differ, which this slide shows. Results are reported as published; review/notes/2026-10-28-rho-scale-note.md records a quick check of how the generating weights affect the comparison, which would be the first thing to examine in follow-up work.
- Carry-over: both

### B11. Spillover and spatial dependence are different parameters
- Figure: none
- Slide text: Spillover $\psi$ Caused by treatment crossing a boundary. A coefficient on the exposure indicator $z_i$, of direct scientific interest, and zero if no one is treated. Spatial dependence $\rho$ Background correlation from shared conditions, present with or without treatment. A nuisance parameter that must be modeled for valid inference. Both use the same rook neighbors, so they are hard to separate when both are …
- Notes: Use when asked "aren't spillover and spatial autocorrelation the same thing?"
- Carry-over: both

### B12. BSS admits too few allocations for randomization inference, so inference is model-based
- Figure: none
- Slide text: BSS admits 2–6 allocations on these grids, so it gives no usable permutation distribution With spillover into controls and correlated outcomes, BSS rules out the clustered allocations that give the worst errors; with spillover into both arms and no spatial dependence, the checkerboards are themselves the worst allocations Maximum likelihood conditions on the realized allocation and gives $\hat\beta$ a sampling …
- Notes: Use when asked "is BSS still randomized?" or "how would you do inference?" Be direct that the checkerboard's estimate is of β − ψ, which is the collinearity on the main slide.
- Carry-over: both

### B13. Buffer-zone designs could restore identifiability
- Figure: none
- Slide text: Intervention Excluded buffer Unexposed control Exclude, or analyze separately, controls that border intervention clusters Then the remaining controls have $z_i = 0$ and $x_i = 0$, and the contrast with treated clusters estimates $\beta$ without collinearity Needs grids large enough that interior controls remain, which these small grids are not Trades the number of usable clusters against separation (McCann et …
- Notes: Use when asked what to do about the collinearity problem.
- Carry-over: both

### B14. β and ψ are structural coefficients, not potential-outcome contrasts
- Figure: none
- Slide text: In the reduced form, a change in one cluster's treatment reaches every cluster through $(\mathbf{I} - \rho\mathbf{W})^{-1}$ Write $Y_i(a, \mathbf{a}_{-i})$ for cluster $i$'s potential outcome when it has assignment $a$ and the other clusters have $\mathbf{a}_{-i}$. An indirect effect such as $E[Y_i(0, \mathbf{1})] - E[Y_i(0, \mathbf{0})]$, the focal cluster held at control while all others are treated versus …
- Notes: Use when asked which causal estimand β corresponds to.
- Carry-over: both

### B15. Simulation algorithm, step by step
- Figure: none
- Slide text: Fix the scenario (grid, $\psi$, $\rho$, spillover type), with $\alpha = 0.2$, $\beta = 1$, $\sigma = 0.1$ Enumerate every allocation with $\lfloor K/2 \rfloor$ of the $K$ clusters treated (70, 126, 924) and flag the BSS ones (2, 6, 2) For each allocation and each of $N$ replications (50; 10 on 3×4): place 160, 180 or 240 subjects at random over the grid (about 20 per cluster), assign $x$ and $z$ from the cluster, …
- Notes: Use when asked how many replications were run or how the BSS allocations were identified.
- Carry-over: both
