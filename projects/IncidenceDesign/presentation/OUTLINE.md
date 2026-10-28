# Talk outline: Sampling Design for Spatial Cluster Randomized Trials Under Heterogeneous Incidence

About 20 minutes of content at the SUDDEN project meeting, Thursday 2026-10-29 (a
one-hour meeting with about 40 minutes of open discussion; questions may come during the
talk, so ±a few minutes is fine). Audience: the clinicians and epidemiologists of the
larger SUD study. Goal: share the current state of the project and collect feedback we
have not thought of. Draft 3 (2026-10-28), built from the V5-revised CTJ manuscript
(commit 00270ae plus the 2026-10-28 declarations fix), its SI, Chapter 3 and the
methods-extension results. Each slide is tagged **prelim** (carries over to prelim
Part II), **technical** (needs a more technical version for biostatisticians) or
**talk-only**.

The deck is `sudden_talk_2026-10-29.qmd` (rendered to `sudden_talk_2026-10-29.html`;
press S for the speaker view with the read-aloud script, one sentence per bullet).
The slide entries below are generated from the deck by
`python3 presentation/tools/outline_from_deck.py`; `outline.html` is the visual
version (`python3 presentation/tools/build_outline_html.py`). Files are listed in
`README.md`.

## Structure and timing

Section names appear in the footer of every slide in that section (Background, Methods,
Simulation Results, Application, Discussion, Backup); there are no section-divider
slides. Timings in the speaker notes come from each slide's word count at 140 words per
minute.

| Section | Slides | Minutes |
|---|---|---|
| Title and outline | 1–2 | about 1 |
| 1. Background | 3–7 | about 5 |
| 2. Methods | 8–11 | about 4 |
| 3. Simulation Results | 12–13 | about 1.5 |
| 4. Application | 14–23 | about 8.5 |
| 5. Discussion | 24–26 | about 2 |
| Main body | 1–26 | about 21.5 at 140 words/min (3,040 words) |
| Backup | B1–B17 | as needed |

## Author decisions shaping Draft 3 (2026-10-28; reasons in `DECISIONS.md`)

- Title is the manuscript title; only Andrew on the title slide; Dr. Lin (advisor),
  Dr. Simpson and Ashkan Habib (collaborators) credited on the thank-you slide, with UNC
  IRB protocol 14-2036 and "manuscript in preparation for Clinical Trials".
- Background covers SUD's relevance, Dr. Simpson's Lenoir County pilot idea and its
  progression to a statewide trial, and how the 100 counties became 58 college service
  areas. Flowcharts show SUD case identification and the application study pipeline.
- The model slide shows the full equation with every term defined, read in words, and
  explains how each trial is simulated and estimated.
- Take-home: Incidence-Guided Saturation is a reasonable middle ground with limited
  downside across both spillover regimes; the final choice still depends on the
  expected spillover and participation (refines the no-single-recommended-design rule).
- Saturation Quadrants, the targeting slide and the questions for the group are backup
  only; the questions are raised as needed, some to Dr. Simpson.
- Regime names stay "controls only" and "both arms", glossed as untreated colleges gain /
  every college gains.
- Script is read aloud, one sentence per bullet.
- An outline (agenda) slide after the title lists the five sections.
- Case identification follows Gan et al.'s documented steps (age, place, manner, cancer);
  the statewide data's additional exclusion rules are mentioned generically, to clarify
  verbally. The pilot's trainees (students and staff in workforce-development and
  disaster-response programs, possibly EMS and healthcare trainees) and its process
  measures (taking medicine, timely clinic visits) are flagged for Dr. Simpson to confirm.

## How Draft 1 shaped the deck

Draft 1 (2026-10-25, in git history) set the order (motivation, models, designs, grid,
NC application, conclusions), most claim titles, the plain regime names, the audience
questions and most backup topics. Kept from its notes: ask how many colleges will join;
guidance costs no accuracy, a defensible equity choice; the colleges-needed counts use a
test with 5% size. Changed since: numbers re-checked against the V5 manuscript; the
model slide gained the equation and glossary; budgets/reach and effect-size translation
became main slides; the grid result is a table; Draft 1's `design_examples.png` (in-figure
title, numbered designs) was replaced by `nc_design_examples.png`; the allocation-risk
figure was redrawn at slide size (`deck_allocation_distribution.png`).

## Open items

- Who is assessed in control colleges (eligible EMS and healthcare professionals in every
  service area, trained or not)? To confirm with Dr. Simpson.
- Before the prelim, re-review the content and move it to the prelim deck's format
  (`bios-dissertation/slides/prelim/part-2/`), only after Andrew approves.

<!-- SLIDES (generated) -->

## 1. Background

### 1. Title
- Figure: none
- Slide text: Sampling Design for Spatial Cluster Randomized Trials Under Heterogeneous Incidence Andrew J Walther Department of Biostatistics, University of North Carolina at Chapel Hill SUDDEN project meeting · October 29, 2026
- Notes: Good morning, and thank you for having me. I'm Andrew Walther, a PhD student in biostatistics working with Dr. Lin.
- Time: 0.5
- Carry-over: talk-only

### 2. Outline
- Figure: none
- Slide text: • Background Sudden unexpected death, Lenoir County pilot, spillover, and community college service area clusters • Methods Spillover variants, spatial autoregressive outcome model, sampling design evaluation, candidate designs vs. simple random sampling • Simulation Results Design intervention estimation accuracy on a simulated 100-cluster region • Application North Carolina SUD data, application of designs to …
- Notes: Here is the plan. I'll start with why sudden death matters and where the pilot idea came from, then the methods, a short simulation section, and most of the time on the North Carolina application, before the take-home and next steps.
- Time: 0.25
- Carry-over: prelim

### 3. Sudden unexpected death takes a large and uneven toll on working-age adults
- Figure: none
- Slide text: Nationally Sudden unexpected death (SUD): an unexpected natural death outside the hospital in a working-age adult About 49 per 100,000 men and 22 per 100,000 women aged 20–64 each year (US estimate from SUDDEN adjudication in Wake County) About 2 million years of life lost, 10% of all-cause years of life lost at these ages, and about $51 billion in lifetime productivity loss Mirzaei et al., Medical Care 2019; Gan …
- Notes: Sudden unexpected death is why this project exists, an unexpected natural death outside the hospital in a working-age adult. Dr. Simpson's group estimated about 49 per 100,000 men and 22 per 100,000 women aged 20 to 64 die this way each year, about two million years of life lost.
- Time: 0.75
- Carry-over: prelim

### 4. Dr. Simpson's Lenoir County pilot is the starting point for a statewide trial
- Figure: none
- Slide text: SUD as a social problem SUDDEN findings show sudden death tracks access to care, income and rurality, so prevention has to reach communities, not only hospitals → Lenoir County pilot idea Teach SUD prevention through Lenoir Community College's workforce-development and disaster-response training (Lenoir, Greene and Jones counties), where Dr. Simpson notes SUD rates spike after flooding → Trainees (to confirm) …
- Notes: SUDDEN showed that sudden death tracks access to care, income and rurality, so prevention has to reach communities, not only hospitals. Dr. Simpson's pilot idea is to teach SUD prevention through Lenoir Community College's workforce-development and disaster-response training, which serves Lenoir, Greene and Jones counties.
- Time: 1.5
- Carry-over: prelim

### 5. When a program reaches neighbors, a simple comparison can miss its effect
- Figure: none
- Slide text: Bed nets, northern Ghana Children living near households with insecticide-treated bed nets had lower mortality, even when their own household had none. School deworming, western Kenya Untreated children attending schools near treated schools had fewer worm infections. Spillover means that treating one group changes outcomes in the groups next to it A trial that compares treated with untreated clusters assumes …
- Notes: Two classic examples show why geography matters. In northern Ghana, children living near households with treated bed nets had lower mortality, even without a net of their own.
- Time: 0.75
- Carry-over: prelim

### 6. Analysis methods exist, but guidance on which clusters to treat is thin
- Figure: none
- Slide text: After the trial (analysis) Spatial models separate a program's direct effect from the effect of treated neighbors They need a trial designed to give them that contrast Before the trial (design) Which clusters should be treated? Options serve different goals (balance a covariate, vary the treated share by region, keep treated clusters apart) Few studies compare them under spillover and uneven incidence together …
- Notes: Statisticians have good tools for analyzing spillover after a trial, but only if the design gave them enough contrast to work with, and that is decided when we choose which clusters to treat. Several allocation ideas exist, and few studies compare them when spillover and uneven incidence occur together.
- Time: 0.5
- Carry-over: technical

### 7. The trial's clusters are North Carolina's 58 community college service areas
- Figure: cc_service_areas_cropped.png
- Slide text: 100 counties Each assigned to its primary serving college (state list of colleges and the counties they serve) 58 service areas 31 colleges serve one county, 18 serve two, 9 serve three to seven Shared counties Bertie and Northampton (hatched) also lie in Roanoke-Chowan's network; each is assigned to its primary college Neighbors Two service areas are neighbors if they share a border or a corner (legal …
- Notes: To set up the trial, we grouped the state's 100 counties into the service areas of the 58 community colleges, using the state's list of which college serves which county. Each color is one college's service area, and those are our clusters.
- Time: 1
- Carry-over: prelim


## 2. Methods

### 8. Spillover can reach only the untreated colleges, or every college
- Figure: none
- Slide text: Spillover to controls only untreated colleges gain Spillover to both arms treated and untreated colleges gain Treated college Control college Gains from treated neighbors (more treated neighbors, larger gain) Controls only is classic contamination, in which untreated colleges learn from treated neighbors Both arms means treated colleges also gain from treated neighbors, for example through shared trainings Which …
- Notes: This is the most important idea in the talk. Spillover can work in two ways.
- Time: 1
- Carry-over: prelim

### 9. One outcome model links treatment, spillover and neighbors
- Figure: none
- Slide text: $$ Y = (I - \rho W)^{-1}\left[\,\tau Z + S(Z) + \beta X + \epsilon\,\right], \qquad \epsilon \sim N(0, I) $$ Spillover to both arms $S(Z) = \gamma\, WZ$ Spillover to controls only $S(Z) = \gamma\,(1 - Z) \circ WZ$ $Y$ Outcome of each cluster (in the application, a college's mean education score) $Z$ Treatment indicator, 1 if the cluster is treated and 0 if control $W$ Neighbor matrix. Under queen contiguity two …
- Notes: You don't need the algebra, since every symbol is defined below, so I'll say it in words. A college's score Y is the direct effect of its own treatment, tau times Z, plus spillover from treated neighbors, S of Z, plus a possible effect of incidence, beta times X, plus noise.
- Time: 1
- Carry-over: prelim

### 10. We judge each design by its average error, its intervals and its worst allocations
- Figure: none
- Slide text: We estimate the program's effect $\tau$ with a spatial model that allows for spillover (a spatial lag model fitted by maximum likelihood, including the true spillover term), and repeat the trial thousands of times in simulation. Mean squared error (MSE) The average of $(\hat\tau - \tau)^2$ over simulated trials. It combines bias and imprecision. Lower is better. Coverage The share of 95% confidence intervals that …
- Notes: We judge each design three ways. We simulate the trial thousands of times and estimate the effect each time with a spatial model that includes the spillover.
- Time: 0.75
- Carry-over: prelim

### 11. Six designs and simple random sampling choose treated clusters differently
- Figure: design_assignments_panels_wide.png
- Slide text: One example assignment per design on the 100-cluster grid (blue, treated; white, control). Top-left panel is the incidence surface used by every example (darker is higher). Percentages give each quadrant's treated share; numbers in Balanced Quartiles give each cluster's incidence quartile. Target incidence High Incidence Focus · Balance incidence Balanced Quartiles · Separate or alternate Checkerboard, Isolation …
- Notes: Here are the designs on a simulated 10 by 10 grid. The top-left panel is the incidence map, darker meaning higher, and each other panel is one example allocation, with blue treated.
- Time: 1
- Carry-over: prelim


## 3. Simulation Results

### 12. The grid study varied every condition we cannot pin down in advance
- Figure: none
- Slide text: Condition Values Study region 10 × 10 grid of 100 clusters; queen neighbors (rook, edges only, as a sensitivity) How incidence varies across space Independent, spatially smooth (2 strengths), Poisson death counts (2 strengths) Spatial dependence of outcomes, $\rho$ 0, 0.01, 0.2, 0.5 Spillover magnitude, $\gamma$ 0.5, 0.6, 0.7, 0.8 Spillover regime Controls only, both arms Treatment effect, $\tau$ 0.8, 1, 1.5, 2, …
- Notes: We first tested the designs on this grid, where we could vary every condition we cannot pin down in advance, from the incidence pattern to the spillover regime, with 250 simulated trials per scenario. I'll show an effect of one standard deviation.
- Time: 0.5
- Carry-over: technical

### 13. On the grid, the spillover regime decides which design wins
- Figure: none
- Slide text: Design Spillover to controls only Spillover to both arms Incidence-Guided Saturation 0.63 1.08 Isolation Buffer 0.74 3.80 Balanced Quartiles 1.02 0.99 High Incidence Focus 1.53 3.82 2×2 Blocking 2.68 1.47 Checkerboard 9.84 1.47 Simple Random Sampling (MSE) 0.215 0.042 MSE relative to SRS, grid, queen neighbors, $\tau = 1$; blue lower than SRS, orange higher, grey not distinguishable from SRS (paired 95% interval …
- Notes: The two regimes tell different stories. With spillover to controls only, Incidence-Guided Saturation has 37 percent lower error than random sampling, and Isolation Buffer is also better.
- Time: 0.75
- Carry-over: prelim


## 4. Application

### 14. SUD deaths were identified from death certificates with a validated algorithm
- Figure: none
- Slide text: 412,514 North Carolina deaths, 2018–2021 → Exclude ages under 18 or over 64 ↓ 111,665 deaths among adults aged 18–64 → Exclude deaths in a hospital, hospice or nursing home ↓ Out-of-hospital deaths → Exclude suicide, homicide and accidents ↓ Natural out-of-hospital deaths → Exclude expected deaths from cancer ↓ 23,523 presumed sudden unexpected deaths The algorithm (Gan et al. 2019) was validated against SUDDEN …
- Notes: For the North Carolina data, Ashkan identified sudden deaths from death certificates with the algorithm from Gan and colleagues, which was validated against SUDDEN's physician adjudication. Of about 412,000 deaths in the state from 2018 to 2021, about 112,000 were in adults aged 18 to 64.
- Time: 0.75
- Carry-over: prelim

### 15. We built the application study in six steps, fixing every rule before looking at results
- Figure: none
- Slide text: 1 · Data County-year SUD deaths and SEER population aged 18–64, 2018–2021 → 2 · Clusters 100 counties into 58 college service areas; neighbors share a border or a corner → 3 · Incidence Yearly SUD rate per 100,000 for each service area; each year a separate planning input → 4 · Designs Each rule adapted to the map: four fixed regions, spatial blocks, a map Checkerboard → 5 · Simulate Education scores from the …
- Notes: Here is how we built the application. We started from the county counts and population, grouped counties into the 58 service areas, and computed each area's SUD rate for each year.
- Time: 0.75
- Carry-over: prelim

### 16. SUD incidence is uneven across North Carolina and shifts from year to year
- Figure: nc_planning_map.png
- Slide text: Panel A, SUD deaths per 100,000 adults aged 18–64, pooled 2018–2021, by college service area. Panel B, the four fixed regions used by Incidence-Guided Saturation. Panel C, yearly rates on the same color scale. 23,523 deaths over 25.6 million person-years (ages 18–64), corrected county death-certificate counts (Gan et al. case definition) Each year is a separate planning input; incidence decides who is treated, …
- Notes: Panel A shows the resulting SUD rates per 100,000 adults aged 18 to 64, pooled over 2018 to 2021, for each service area. They range from about 45 to 195 per 100,000, and the highest tend to be in service areas with small populations.
- Time: 0.75
- Carry-over: prelim

### 17. Each design was adapted to the map before we looked at any results
- Figure: nc_design_examples.png
- Slide text: One example allocation per design on the 58-college map, 2018 incidence (Isolation Buffer treats 14 here, about 16 on average). Checkerboard becomes an alternating pattern on the irregular map; 2×2 Blocking becomes spatial blocking of 14 groups of nearby colleges; Incidence-Guided Saturation uses the four fixed regions. Most designs treat 29 of 58 colleges. Checkerboard treats 27, spatial blocking 27–31, and …
- Notes: Each design was translated to the real map, with every rule fixed before we looked at results. Here is one example allocation of each for 2018, with blue colleges treated.
- Time: 0.75
- Carry-over: prelim

### 18. Treating half the colleges does not always reach half the population
- Figure: none
- Slide text: Design Colleges treated (of 58) 2018 2019 2020 2021 Incidence-Guided Saturation 29 53% 45% 48% 52% Balanced Quartiles 29 49% 50% 51% 50% Isolation Buffer about 16 27% 27% 27% 27% High Incidence Focus 29 22% 24% 23% 22% 2×2 Blocking (spatial blocks) 28 / 31 / 27 / 29 50% 55% 48% 50% Checkerboard 27 50% 50% 50% 51% Simple Random Sampling 29 50% 50% 50% 50% Average share of North Carolina's adults aged 18–64 living …
- Notes: Treating half the colleges does not mean reaching half the people. This table shows the share of the state's working-age adults living in treated service areas.
- Time: 0.5
- Carry-over: prelim

### 19. With 58 colleges, Incidence-Guided Saturation beats SRS under spillover to controls only
- Figure: nc_annual_accuracy_six.png
- Slide text: Mean MSE divided by SRS's in the same year and regime (each cell averages six settings of $\rho$ and $\gamma$). Blue is lower error than SRS, orange higher. Years share one map and several designs share simulations across years, so years are not independent replications. Both arms: SRS, Balanced Quartiles and Incidence-Guided Saturation are within about 3% of each other Checkerboard is about 3 times SRS's error …
- Notes: Here is accuracy on the real map with all 58 colleges, by year, with spillover to both arms on the left, controls only on the right, and blue better than random sampling. Under controls only, Incidence-Guided Saturation, the bottom row, has 21 to 30 percent lower error than random sampling in every year, the same message as the grid.
- Time: 1.25
- Carry-over: prelim

### 20. Under controls-only spillover, Incidence-Guided Saturation also lowers the risk of a bad draw
- Figure: deck_allocation_distribution.png
- Slide text: Each bar is the spread of error across one design's possible allocations, relative to SRS's average (thin line 5th–95th percentile, thick bar 25th–75th, circle the mean); top, grid; bottom, NC. Dark grey marks SRS. Under controls only, Incidence-Guided Saturation's worst 10% of allocations had 21–43% lower error than SRS's worst 10% in NC (grid 0.19 vs 0.35); under both arms the differences are small and mostly …
- Notes: A real trial gets one allocation, so we also asked how bad an unlucky draw can be. Each bar is the spread of error across a design's possible allocations, relative to the average random allocation, with the grid on top and North Carolina below.
- Time: 0.75
- Carry-over: prelim

### 21. About 36 colleges give 80% power at one SD when spillover reaches every college
- Figure: fig_colleges_needed.png
- Slide text: Colleges needed for 80% power (corrected $t$ test, 5% false-positive rate), random volunteering, $\rho = 0.5$, spillover from a fully treated neighborhood equal to 80% of the direct effect. Effects are in college-level standard deviations (SD). Bars, mean over 2018–2021; whiskers, range; grey, not reached with 58 colleges in at least one year. Saturation Quadrants and Balanced Halves are backup variants. Counts …
- Notes: How many colleges do we need? This counts participating colleges, about half treated and half control, needed for an 80 percent chance of detecting the effect, with effects in college-level standard deviations, roughly a third of an SD for an individual professional.
- Time: 0.75
- Carry-over: prelim

### 22. If fewer colleges join, Incidence-Guided Saturation's advantage shrinks
- Figure: fig_power_vs_n.png
- Slide text: Participating colleges Incidence-Guided Saturation error vs SRS (controls only) 58 22% lower 40 9% lower 29 5% lower 20 about 1% lower Under controls only, SRS power stalls past about 40 colleges (Incidence-Guided Saturation keeps rising slowly), since nearly every control then has treated neighbors If higher-incidence colleges are more likely to volunteer, the advantage is 2% or less at 20–29 colleges Figure, …
- Notes: Not every college will join. Under spillover to both arms, on the left, every added college helps.
- Time: 0.75
- Carry-over: prelim

### 23. A one-SD college effect is a modest effect for each trained professional
- Figure: none
- Slide text: A college score averages $m$ professionals. With intraclass correlation ICC (share of score variance between colleges), the college-level noise SD ($\epsilon$), in individual SD units, is $$\sigma_c = \sqrt{\text{ICC} + (1 - \text{ICC})/m}$$ so a college effect of $\tau$ SDs is an individual effect of $d = \tau\,\sigma_c$. With 25 professionals per college and ICC = 0.05 , $\sigma_c \approx 0.30$, so a one-SD …
- Notes: What does one college-level standard deviation mean in practice? The college score averages the professionals assessed there, and averaging reduces noise.
- Time: 0.75
- Carry-over: prelim


## 5. Discussion

### 24. Incidence-Guided Saturation was a reasonable middle ground with limited downside
- Figure: none
- Slide text: Incidence-Guided Saturation gained the most under controls-only spillover, at a small cost under both arms Spillover to controls only (untreated colleges gain) Lowest or near-lowest error with all colleges taking part (21–30% below SRS in NC, 37% on the grid), a lower risk of an unlucky allocation, and the most treatment in the highest-rate region (a small shift on the NC map, where regions are similar). Even …
- Notes: The goal is a design that does reasonably well whichever kind of spillover is true, with limited downside. Incidence-Guided Saturation gained the most when spillover reaches untreated colleges, at a small cost when it reaches everyone, while Balanced Quartiles stayed close to random sampling in both.
- Time: 1
- Carry-over: prelim

### 25. Key assumptions remain, and the next step is to define the endpoint and size the trial
- Figure: none
- Slide text: Limitations Our analysis assumed we know how spillover works; if the analysis ignores spillover, every design underestimates the effect Planning incidence was treated as exact and fixed Neighboring service areas stand in for how professionals actually communicate Designs treat different numbers of colleges, so this is not a comparison at equal cost This evaluates trial designs, not the program's effect on SUD …
- Notes: A few limits. Our analysis assumed we know how spillover works, and if the analysis ignores spillover, every design underestimates the effect.
- Time: 0.75
- Carry-over: prelim

### 26. Thank you!
- Figure: none
- Slide text: Advisor Feng-Chang Lin Collaborators Ross J. Simpson, Jr. · Ashkan Habib Data Corrected NC county SUD counts (A. Habib), UNC IRB protocol 14-2036; SEER population estimates Papers Design study manuscript in preparation for *Clinical Trials* Earlier work, Walther et al., *BMC Medical Research Methodology* (2026) Questions welcome awalther@unc.edu
- Notes: Thank you to my advisor, Dr. Lin, to Dr. Simpson for the program that motivated this work, and to Ashkan Habib for preparing the county death data. I'm happy to take questions.
- Time: 0.25
- Carry-over: talk-only


## Backup

### B1. Questions we could put to the group
- Figure: none
- Slide text: How does the program spread? Mainly to untreated colleges nearby, or to every college, treated ones included? Is a shared border the right idea of "neighbor"? What is the education outcome , when would it be measured, and what improvement would be worth detecting? How many colleges and professionals? How many colleges would join (and would higher-incidence ones be keener), how many professionals each would train, …
- Notes: Use to steer discussion if it stalls, or direct individual questions to Dr. Simpson. These answers decide the design and the sample size: the spillover regime, the endpoint and the effect worth detecting, participation and professionals per college, and whether to target rates or death counts.
- Carry-over: technical

### B2. Incidence guidance sends more treatment to high-burden regions only where regions differ
- Figure: none
- Slide text: Design (NC, mean over years) Treated share of total incidence (rate-weighted) Share of SUD deaths in treated areas High Incidence Focus 62.2% 35.4% Incidence-Guided Saturation 51.1% 49.2% Saturation Quadrants (same rule, random region order) 50.4% 49.6% North Carolina's four regions have similar average incidence (mean ranks 24–34 of 58), so guidance shifts the rate share by less than one point here On the grid, …
- Notes: Use when asked whether incidence guidance actually reaches the burden. NC's four regions have similar mean incidence ranks (24–34 of 58), so guidance shifts the treated rate share by under one point; on the grid the gain was 2.1–3.1 points with the same accuracy as random-order saturation. High Incidence Focus reaches 62% of the rate share but only 35% of deaths, since high-rate areas are small. Source: Chapter 3 application results.
- Carry-over: technical

### B3. Full rules for each design on the grid and the NC map
- Figure: none
- Slide text: Design Grid (100 clusters) NC (58 colleges) Checkerboard Alternating pattern; no two clusters sharing an edge in the same arm (50 treated) Deterministic map version: visits colleges from most to fewest neighbors, gives each the arm less common among its assigned neighbors (27) High Incidence Focus The 50 highest-incidence clusters, ties broken at random The 29 highest observed rates Isolation Buffer Add clusters …
- Notes: Use when asked how a design works on the real map, or about the two variants. Saturation Quadrants and Balanced Halves performed like Incidence-Guided Saturation and Balanced Quartiles, so the main slides show six designs plus SRS. Ties in incidence are broken at random for every draw.
- Carry-over: technical

### B4. Grid benchmark with absolute error and coverage
- Figure: none
- Slide text: Design Both arms Controls only MSE Coverage Ratio to SRS (95% CI) MSE Coverage Ratio to SRS (95% CI) Checkerboard 0.0624 0.948 1.47 (1.44, 1.51) 2.1114 0.951 9.84 (9.65, 10.03) High Incidence Focus 0.1617 0.950 3.82 (3.71, 3.92) 0.3292 0.945 1.53 (1.45, 1.62) Isolation Buffer 0.1608 0.948 3.80 (3.72, 3.87) 0.1595 0.950 0.74 (0.72, 0.77) 2×2 Blocking 0.0623 0.947 1.47 (1.44, 1.50) 0.5761 0.950 2.68 (2.62, 2.74) …
- Notes: Use when asked for the absolute numbers behind the grid results slide, or whether differences exceed simulation error. Source: CTJ Table 2 (results/srs_benchmark/manuscript_regime_means.csv and manuscript_paired_ratios.csv). Pooled over regimes, Incidence-Guided Saturation 0.0907 vs SRS 0.1285.
- Carry-over: technical

### B5. The grid ordering holds in every incidence configuration
- Figure: grid_accuracy_srs.png
- Slide text: Grid mean MSE at $\tau = 1$, queen neighbors, by incidence configuration (parentheses give incidence dependence $\rho_X$), log scale. Horizontal lines are 95% Monte Carlo intervals over ten incidence surfaces.
- Notes: Use when asked whether the way incidence varies changes the conclusions. It moved each design's MSE by less than 10%, except High Incidence Focus, whose allocation follows incidence directly. Source: CTJ Figure 2.
- Carry-over: technical

### B6. Under rook neighbors, Checkerboard cannot estimate the effect
- Figure: fig_vs_srs_by_regime_rook.png
- Slide text: Rook neighbors share an edge only. Under rook, Checkerboard gives $WZ = 1 - Z$, so exposure cannot be separated from treatment, $\tau$ is not identified, and the estimate centers on $\tau - \gamma$ (coverage 0.16) Under queen, Checkerboard is identified but imprecise, since interior clusters all have half their neighbors treated Other designs keep the regime-specific pattern, with differences. Under controls only …
- Notes: Use when asked about the neighbor definition or why strict alternation fails. In NC, rook and queen at matched settings gave the same broad conclusions. Figure: SI rook comparison with SRS (fig_vs_srs_by_regime_rook).
- Carry-over: technical

### B7. Mean and worst-decile error relative to SRS
- Figure: none
- Slide text: Design Grid NC application Both arms Controls only Both arms Controls only Mean Worst 10% Mean Worst 10% Mean Worst 10% Mean Worst 10% Checkerboard 1.41 1.11 9.94 6.22 2.68 2.20 3.02 1.88 High Incidence Focus 3.81 3.05 1.51 0.96 0.93 0.76 1.03 0.64 Saturation Quadrants 1.06 1.06 0.62 0.54 0.99 0.94 0.73 0.64 Isolation Buffer 3.80 3.97 0.75 0.62 3.69 4.75 0.93 0.91 2×2 Blocking 1.44 1.57 2.69 2.64 1.00 1.02 1.15 …
- Notes: Use when asked for the numbers behind the allocation-risk slide. Each worst decile has only ten sampled allocations, so membership is uncertain; split-half checks (select with half the outcomes, evaluate with the other half) kept the control-only pattern. Under both arms most of the spread is outcome noise. Source: Chapter 3 tab:downside (results/allocation_risk/summary/surfaces10/downside_by_design.csv).
- Carry-over: technical

### B8. Population reached by each design, by year
- Figure: population_shares.png
- Slide text: Share of NC adults aged 18–64 living in treated service areas; points are means over sampled allocations and lines show the spread across allocations. Includes the two variants.
- Notes: Use when asked how much of the state each design would reach, or how variable that is from draw to draw. SRS spans roughly 31% to 69% across draws; High Incidence Focus is fixed at 22–24%. Source: Chapter 3 tab:nc-budgets and SI population-share figure.
- Carry-over: technical

### B9. If incidence did shift education scores ($\beta = 1$)
- Figure: none
- Slide text: A hypothetical check adds incidence rank to the outcome ($\beta = 1$) and adjusted for in the analysis High Incidence Focus error rose about 3.7–4.4-fold under both arms and 1.4–2.5-fold under controls only Aligning treatment with incidence leaves little independent information about $\tau$ once incidence is adjusted for Randomized designs were little affected; the regime-specific guidance is unchanged Not an …
- Notes: Use when asked "what if high-incidence colleges also score differently?" This is the main reason High Incidence Focus is not recommended even though its single allocation looked fine in the primary NC model. Source: CTJ Results; Chapter 3 Allocation downside and sensitivity.
- Carry-over: technical

### B10. Leaving the spillover term out of the analysis biases every design
- Figure: none
- Slide text: Grid sensitivity: fit $Y \sim Z$ (plus incidence) without the spillover term Every identified design underestimated $\tau$: bias $-0.12$ to $-0.40$, coverage 54–91%, at $\tau = 1$ (queen and rook) Some designs had *lower* MSE than under the correct model, since less variance outweighed the bias, so low MSE alone does not mean valid inference Incidence-Guided Saturation still had the lowest MSE under both neighbor …
- Notes: Use when asked "what if we don't know the spillover form?" We did not study a spillover term of the wrong form (wrong network or distance decay); that is listed as future work. Source: CTJ Results; Chapter 3 non-oracle sensitivity.
- Carry-over: technical

### B11. The corrected interval restores 95% coverage
- Figure: none
- Slide text: Design Grid queen Wald / $t$ NC both Wald / $t$ NC controls Wald / $t$ Incidence-Guided Saturation 93.7 / 94.7 93.5 / 94.9 92.7 / 94.2 Balanced Quartiles 93.6 / 94.6 93.3 / 94.8 92.9 / 94.4 Isolation Buffer 93.9 / 94.9 93.2 / 94.7 93.1 / 94.7 High Incidence Focus 93.8 / 94.8 93.4 / 94.9 93.1 / 94.6 2×2 Blocking 93.9 / 94.8 93.2 / 94.7 92.9 / 94.4 Checkerboard 94.1 / 95.0 93.6 / 94.9 93.4 / 94.8 Simple Random …
- Notes: Use when asked why the intervals are not the usual Wald intervals. Maximum likelihood variance divides by N and the normal reference ignores estimated parameters, so with 58 or fewer clusters Wald intervals are too narrow. Bootstrap and profile-likelihood intervals did not meet the 94–96% rule. Sources: results/grid_tdf_2026-10/interval_comparison.csv; Chapter 3 tab:nc-calibration.
- Carry-over: technical

### B12. Power for individual effects by training size and intraclass correlation
- Figure: none
- Slide text: $m$ ICC $\sigma_c$ Both arms Controls only 29 colleges 58 colleges 29 colleges 58 colleges $d$=0.2 $d$=0.3 $d$=0.2 $d$=0.3 $d$=0.2 $d$=0.3 $d$=0.2 $d$=0.3 10 0.01 0.330 0.33 0.62 0.61 0.90 0.22 0.41 0.24 0.45 10 0.05 0.381 0.27 0.51 0.50 0.82 0.18 0.33 0.20 0.36 10 0.10 0.436 0.22 0.41 0.40 0.71 0.15 0.27 0.17 0.30 25 0.01 0.223 0.61 0.91 0.90 1.00 0.41 0.71 0.44 0.74 25 0.05 0.297 0.40 0.71 0.69 0.96 0.26 0.49 …
- Notes: Use when asked about other training sizes or ICCs. Scaling the effect and noise together leaves the test statistic unchanged when β = 0, so these follow from the college-level simulation (checked by rerunning one setting with both halved). Source: Chapter 3 tab:nc-translation.
- Carry-over: technical

### B13. Other spatial models give the same design guidance
- Figure: none
- Slide text: CAR (conditional autoregressive) models: generating or fitting with CAR instead of the spatial lag model left the saturation designs at 0.76–0.82 of SRS's error under controls only, and Checkerboard and Isolation Buffer worst under both arms Bayesian spatial lag model: posterior means correlated 0.999 with maximum likelihood, credible-interval coverage 94.6–95.3%, same design ordering Exception: a Bayesian CAR …
- Notes: Use when asked "did you try a Bayesian model?" or "does this depend on the spatial lag model?" 58 colleges, β = 0. Source: Chapter 3 Appendix A9 (tab:nc-car, tab:nc-bayes).
- Carry-over: technical

### B14. Other forms of spillover
- Figure: none
- Slide text: Design Both arms Controls only Share Count Any treated Share Count Any treated Checkerboard 2.97 1.78 1.07 3.09 1.41 n.i. High Incidence Focus 0.91 0.92 0.96 1.03 1.04 n.i. Isolation Buffer 3.82 3.84 n.i. 0.95 1.16 n.i. 2×2 Blocking 1.04 1.03 1.00 1.17 1.07 n.i. Balanced Quartiles 1.01 1.01 1.00 0.99 0.99 n.i. Incidence-Guided Saturation 1.01 1.03 1.04 0.75 0.79 n.i. MSE relative to SRS, 58 colleges, $\beta = 0$. …
- Notes: Use when asked "what if spillover is all-or-nothing?" The count form keeps the saturation advantage under controls only. Source: Chapter 3 tab:nc-forms.
- Carry-over: technical

### B15. Formal tests agree with the descriptive ranking
- Figure: none
- Slide text: Grid, queen, $\tau = 1$, pooled over spillover regimes: Friedman test across the six designs and SRS over 160 blocks, $\chi^2 = 560.3$, df = 6. Incidence-Guided Saturation ranked best, significantly ahead of every other design and SRS (Nemenyi). SRS and Balanced Quartiles could not be distinguished Adjacent designs in the MSE order differed significantly over the 50 independent configuration-by-surface units …
- Notes: Use when asked "are the differences statistically significant?" Source: Chapter 3 Results (design performance under queen contiguity; Appendix A9 formal tests).
- Carry-over: technical

### B16. Guidance gains grow with how much regions differ
- Figure: fig_spread_vs_gain.png
- Slide text: Each point is one simulated grid incidence surface, plotting the spread of the regions' mean incidence ranks against the extra treated share of incidence from incidence guidance over random-order saturation. NC's four planning years (open triangles) sit at the low-spread end and below the grid trend.
- Notes: Use when asked "why does guidance barely change targeting in NC?" On the grid the gain was 2.1–3.1 percentage points, correlation with spread 0.84–0.99. NC regions' mean ranks lie between 24 and 34 of 58 in every year. Source: Chapter 3 Appendix A9 (fig:spread-gain).
- Carry-over: technical

### B17. At larger effects Isolation Buffer catches up under controls only
- Figure: none
- Slide text: Grid, controls only: Incidence-Guided Saturation's advantage narrows as $\tau$ grows At $\tau = 2$, Isolation Buffer and Incidence-Guided Saturation are level (MSE 0.160 and 0.160); at $\tau = 3$ Isolation Buffer is lower (0.160 vs 0.182) Most designs carry error in estimating $\rho$ into $\hat\tau$ in proportion to $\tau$; Isolation Buffer largely avoids this since no treated cluster has a treated neighbor …
- Notes: Use when asked whether the ranking depends on effect size. Under queen the order by mean MSE among the six designs was otherwise the same at every τ from 0.8 to 3. Source: CTJ Results; Chapter 3 Stability across treatment effect sizes.
- Carry-over: technical
