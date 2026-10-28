# SUDDEN talk (2026-10-29): independent statistician review, round 3

**Date:** 2026-10-28
**Deck reviewed:** `presentation/sudden_talk_2026-10-29.qmd` (1,183 lines, slides and `{.notes}`), figures opened: `nc_annual_accuracy_six.png`, `deck_allocation_distribution.png`, `fig_colleges_needed.png`; script `presentation/tools/deck_allocation_figure.R`
**Sources checked against:** `paper/thesis_chapter/Thesis_Chapter.qmd` (437, 476–492, 516, 528, 534–538, 630–656, 1248, 1523–1529), `CTJ_Manuscript.tex` (165, 187, 215, 246), `application/AGENTS.md` (real-data section), `application/notes/SUD_Meeting_Notes.md`, `application/data/NC_Colleges_and_Service_Areas.csv`, `NC_Community_Colleges_by_County.csv`, `paper/SpatialCRT_IncidenceDesign.bib` (Mirzaei 2019 abstract, Gan 2019, Habib 2026), `results/srs_benchmark/manuscript_paired_ratios.csv`, `manuscript_regime_means.csv`, `results/allocation_risk/summary/surfaces10/allocation_distribution_by_design.csv`, `results/methods_ext_figures/tables/tab_nc_power.tex`, `tab_nc_guidance.tex`, `tables_log.txt`, `ROADMAP.md` (1006–1010), `paper/archive/pre_completion_20261002/Dissertation_Chapter.qmd` (407), `DECISIONS.md` (2026-10-28), round-2 review.

**Verdict:** Blocking issues: 0 (3 major, 14 minor)

The round-2 major (power target missing from the colleges-needed title and take-home) is fixed, and so are most round-2 minors (× in the figure, "power stalls" wording, σ_c defined as the noise SD, "roughly fourfold", the allocation-risk both-arms wording, the formal-tests unit statement, the rook colon, the β = 1 grammar). The results numbers on the new and changed slides all trace to sources. The new allocation figure reproduces its source CSV. The problems are in the new background material and in how the take-home is worded. The case-identification flowchart states exclusion rules that the current manuscripts deliberately do not state. The take-home headline generalizes beyond the evidence. And the main body now claims targeting without the "small gain on the NC map" qualifier the standing rules require.

---

## BLOCKING

None.

---

## MAJOR

### M1. The case-identification flowchart states upstream rules that the manuscripts decline to state
- **Slide:** "SUD deaths were identified from death certificates with a validated algorithm" (lines 486–501; notes 508–518).
- **Text:** "412,514 North Carolina deaths ... 111,665 deaths among adults aged 18–64 → Exclude deaths in hospitals, hospice, nursing homes and other institutions ... Exclude expected deaths (accident, homicide, suicide, overdose, COVID-19) ... Exclude end-stage chronic disease (cancer, cirrhosis, renal failure)".
- **Problem:** The current chapter says the opposite of providing this detail: "We do not infer undocumented upstream exclusion rules or zip-code allocation from the received county panel" (Thesis_Chapter.qmd 437). The CTJ says only "the case-adjudication algorithm of Gan et al." (187). The last place these counts and exclusions appeared, the archived chapter (`pre_completion_20261002/Dissertation_Chapter.qmd` 407), tagged both as `NEEDS-AUTHOR-CONFIRMATION`. It said 412,514 and 111,665 "appear only in the AGENTS.md attribution note ... confirm against Habib (2026)", and that the exclusion list was undocumented. The deck now adds more detail than that draft had (hospice, nursing homes, homicide, suicide, cirrhosis, renal failure) and an order of steps. Nothing in the repository documents that detail. Two further points:
  - (a) 412,514 and 111,665 belong to the derivation behind Habib's 21,147 (ROADMAP 1008). The slide presents them as the path to the corrected 23,523, and that needs Habib's confirmation.
  - (b) `application/AGENTS.md` notes that Habib's free-text regex "can't match 'KIDNEY FAILURE' or 'LIVER FAILURE'". So "renal failure" and "cirrhosis" may not be excluded the way the box says, at least in the free-text step.
  
  Habib and Dr. Simpson will be in the room. An error here costs credibility on the data slide that the clinicians know best.
- **Fix:** Before the talk, have Ashkan Habib confirm the counts, the step order and the exclusion lists (send him the slide). If he can't confirm in time, cut the flowchart to what is documented: "NC death certificates, 2018–2021, adults aged 18–64 → Gan et al. algorithm (out-of-hospital, unexpected, natural deaths; details in Habib 2026) → 23,523 presumed SUD deaths (corrected counts that keep heart-failure deaths; 21,147 in Habib 2026)". Put "Counts and exclusions from Habib (2026), unpublished" in a source line under the flowchart.

### M2. The take-home headline overgeneralizes and leaves out the design with the smallest downside
- **Slide:** Take-home (line 799; notes 820–822).
- **Text:** "Several designs did well in one spillover regime and poorly in the other; Incidence-Guided Saturation did reasonably well in both". Notes: "Incidence-Guided Saturation came closest."
- **Problem:** On the grid, only Isolation Buffer did well in one regime and poorly in the other (0.74 / 3.80). Checkerboard, 2×2 Blocking and High Incidence Focus were worse than SRS in both regimes (`manuscript_paired_ratios.csv`). In NC, only High Incidence Focus switches (0.92–0.95 under both arms, 0.70–1.40 under controls only). "Several" is not supported. Balanced Quartiles is the more important omission. It was within 2% of SRS in both regimes on the grid (1.02 / 0.99) and in NC (0.92–1.00), so its downside is smaller than saturation's grid both-arms cost (1.08, paired CI 1.06–1.11, significantly above 1). The middle-ground framing holds only if the trade is stated: saturation buys a 21–37% gain under controls only for an 8% grid cost (about 0–3% in NC) under both arms. Without that, the headline reads as "the design that does well everywhere", which edges toward the single recommended design the 2026-10-14 rule forbids.
- **Fix:** "Isolation Buffer did well in one spillover regime and poorly in the other. Balanced Quartiles matched SRS in both but never beat it. Incidence-Guided Saturation gained the most under controls only and gave up little under both arms." Notes 822: "Incidence-Guided Saturation offered the best trade, a large gain under controls only for a small cost under both arms."

### M3. The main body claims targeting with no NC qualifier
- **Slides:** Take-home controls-only card (line 805): "...and the most treatment in the highest-rate region". NC accuracy notes (671): "Saturation sends more of the program to higher-burden regions while keeping randomization in every region."
- **Problem:** The targeting slide is now in backup. That leaves an unqualified claim in the main body, in front of an audience that asked "where should prevention start". In NC the four regions' mean ranks lie between 24 and 34 of 58. Incidence-Guided Saturation treats 51.1% of the summed rate share against 50.4% for Saturation Quadrants, and 49.2% of deaths, which is slightly below half (chapter 538). The AGENTS manuscript-conventions rule says to frame guidance as targeting high-burden regions "when burden differs between regions, with a small gain on the NC map". Both lines are literally true, but a clinician will hear that the design reaches more of the burden, and in NC it barely does.
- **Fix:** Card: "...and more treatment in the highest-rate region (a small shift in NC, where regions differ little)". Notes 671: "Saturation gives the highest-rate region the most treatment while keeping randomization in every region, although in North Carolina the regions differ little, so the shift is small."

---

## MINOR

1. **Service-area counts don't match the clusters** (Service-area slide, line 216; notes 232). "30 colleges serve one county, 18 serve two, 10 serve three or more" is the public list (`NC_Colleges_and_Service_Areas.csv`: 30/18/7/2/1), where Roanoke-Chowan serves three counties. Under the analysis assignment (Hertford alone to Roanoke-Chowan, chapter 437; `NC_Community_Colleges_by_County.csv` primary colleges) the 58 clusters are **31 single-county, 18 two-county, 9 three or more (up to seven)**. The box is headed "58 service areas" and sits beside the cluster map. Fix: "31 serve one county, 18 serve two, 9 serve three to seven", and notes "Thirty-one colleges serve a single county".
2. **Unsupported reason for the rate gap** (SUD burden notes, line 93). "These algorithm-based rates run higher than the adjudicated national figure, since the algorithm is broader than a physician panel." No source gives this explanation. The comparisons also differ in other ways: ages 18–64 against 20–64, NC against Wake County extrapolated to the US, 2018–2021 (including COVID years) against 2013–2015, and heart-failure deaths kept. Dr. Simpson's group will know. Fix: "These algorithm-based rates run higher than the adjudicated national figure, partly because a death-certificate algorithm is broader than physician adjudication and partly because the years, ages and region differ."
3. **Mirzaei figures** (SUD burden slide, line 69; notes 89). The abstract says "the loss of 2 million years of life" and "Lifetime productivity loss ... about $51 billion". It does not say "a year". Fix: check the full text for the annual framing, or drop "a year" and write "about $51 billion in lifetime productivity lost".
4. **Rural claim not attributed** (SUD burden slide, line 79; notes 95). "Rates are higher in rural counties" sits in the NC 2018–2021 column, but the source is Gan et al. 2019 (earlier data). The deck's own data support "highest in small-population service areas" (chapter 538). Fix: "Rates are higher in rural counties (Gan et al. 2019) and differ across the state".
5. **Pilot slide embellishments** (line 118; notes 132, 136). "Neighboring colleges share staff, trainings and patients" and "recognize, report and respond" have no source. The meeting notes say colleges "talk to each other" and have satellite campuses, and that the goal is to "get people to the hospital and share treatments". DECISIONS 2026-10-28 confirms only the trainee population. Colleges don't have patients. Fix: "neighboring colleges talk to each other and share satellite campuses", and "so they can help people at risk get timely care". Also, the meeting notes give the process measures as medication adherence and timely clinic visits, while the education score is the trial endpoint the authors chose. Consider "process measures (timely care, medication adherence) and, for this trial, an education score".
6. **Algorithm citation** (Case identification, line 499). Gan et al. 2019 is a JGIM research letter on rural incidence. ROADMAP 1008 attributes the validated algorithm to `gan_factors_2019` and `nanavati_sudden_2014` together. Fix: "(Gan et al. 2019; Nanavati et al. 2014)", or confirm with Dr. Simpson which paper validates it.
7. **"Several designs" companion bullet** (Take-home, line 814). "High Incidence Focus varied by year" is true only under controls only (round-2 m3/m11 carry-over). Under both arms in NC it was steadily 5–8% below SRS. Fix: "High Incidence Focus varied by year under controls only and is not randomized".
8. **"21–37% below SRS" has no participation condition** (Take-home, line 805). The range holds with all 58 colleges (NC 21–30%) or the full grid (37%). At 29 colleges it is 5% (`tab_nc_guidance.tex`). The second bullet says this, but the number will be quoted alone. Fix: "(21–37% below SRS with all colleges taking part)".
9. **Power range mixes designs** (Effect-size fig-note line 778; notes 793; take-home 805, 824). "Under controls only, power is about 0.5–0.65 even with all 58" sits under an SRS table. For SRS at m = 25, ICC = 0.05, d = 0.3 it is 0.53. The 0.65 is Incidence-Guided Saturation's 0.64 (chapter 652). Fix: fig-note "Under controls only, with all 58, power is 0.53 with SRS and 0.64 with Incidence-Guided Saturation". In the take-home IGS card, "only about 64%".
10. **Lower-case "saturation" as a design name** (titles at lines 641 and 674, line 720; bullets 455, 472). With Saturation Quadrants in backup the meaning is clear, but AGENTS requires one design-name set "in text, tables and figures alike". Fix: use "Incidence-Guided Saturation" in the titles ("...Incidence-Guided Saturation beats SRS when...", "...also lowers the risk..."). "Saturation creates regions..." (476) is a general mechanism and can stay.
11. **Prose colons.** Line 101 "SUDDEN findings: sudden death tracks..." → "SUDDEN found that sudden death tracks access to care, income and rurality, so...". Backup notes line 911 "the sample size: the spillover regime..." → "the sample size, namely the spillover regime...". Label colons ("**Both arms:**", "Grid:", "Share:") are fine.
12. **Dangling promise** (Spillover regimes notes, line 301). "I'll ask for your view at the end." The questions slide is now backup. Fix: "I'd value your view on this in the discussion", or plan to open the backup questions slide.
13. **"Shared borders"** (Pipeline step 2, line 526). Queen neighbors also include corners (4 corner-only pairs). Fix: "neighbors share a border or corner".
14. **Allocation-risk notes** (line 689). "Saturation's whole bar sits to the left." In NC the thin line reaches 1.02 (Q95). Fix: "Saturation's thick bar sits well to the left, and its long tail is much shorter than random sampling's."

Also noted, not a finding: `fig_colleges_needed.png` now prints "2×2 Blocking" correctly but shows the two backup variants (Saturation Quadrants, Balanced Halves) without a caption remark. Its label "> 58 in 1/4 yr, else 58+" for Isolation Buffer reads oddly. A fig-note phrase such as "(includes the two backup variants)" would help.

---

## Notation, framing and style: checked

- **Notation.** All symbols are defined. τ is "Direct effect of a cluster's own treatment before outcomes spread between neighbors through $(I-\rho W)^{-1}$ (not a total effect), in standard deviations of the noise ε" (318). β: grid 1, NC 0 primary, 1 only as a hypothetical sensitivity (323; backup 1035–1041). X is defined for both grid and NC (324). σ_c is now the college-level noise SD (761). The two S(Z) forms match the CTJ. m(a) is defined.
- **No-single-design rule.** The take-home keeps "the right choice still depends on the spillover we expect" (815, notes 830), and the grid slide says "there is no single best design". These are consistent with DECISIONS 2026-10-28 once M2 and M3 are fixed.
- **Style.** No DIM, Claude or Codex. No sentence starts with "Because". Design names are exact except the lower-case shorthand (m10). IRB 14-2036 and "in preparation for *Clinical Trials*" are author-approved (DECISIONS 2026-10-28).

## Allocation figure against its CSV

`deck_allocation_distribution.png` is drawn only from `results/allocation_risk/summary/surfaces10/allocation_distribution_by_design.csv`. The script asserts 7 designs × 2 studies × 2 regimes = 28 rows, and the CSV has them. Every plotted value I checked matches the CSV:
- Grid, controls only: IGS mean 0.647, Q05 0.466, Q95 0.856; SRS Q95 1.52; Checkerboard 9.95; 2×2 Blocking mean 2.69, Q95 4.04.
- Grid, both arms: Isolation Buffer mean 3.80, Q95 4.91.
- NC, controls only: IGS mean 0.750, Q25–Q75 0.637–0.823, Q95 1.02; SRS Q95 1.44; Checkerboard 3.02; HIF 1.03.
- NC, both arms: Isolation Buffer mean 3.69, Q95 5.36; HIF 0.93.

Design order, names (with ×), the log axis and the dark-grey SRS all match the fig-note and the speaker notes.

---

## Number spot-checks (42)

| # | Deck number | Slide | Source | Match |
|--:|---|---|---|:-:|
| 1 | 49 / 22 per 100,000, ages 20–64 | SUD burden | Mirzaei abstract 49.3 / 21.7 | yes |
| 2 | ~2 million YLL, 10% | SUD burden | Mirzaei abstract | yes ("a year" unverified, m3) |
| 3 | $51 billion | SUD burden | Mirzaei abstract (lifetime) | yes (m3) |
| 4 | 23,523 deaths | SUD burden, case ID, NC incidence | application/AGENTS.md; CTJ 187 | yes |
| 5 | 84 → 101, rising each year | SUD burden | AGENTS 83.9/85.9/97.1/100.5 | yes |
| 6 | 25.6 million person-years | Case ID, NC incidence | 25,594,321 | yes |
| 7 | 21,147 earlier counts | Case ID | AGENTS | yes |
| 8 | 412,514 / 111,665 | Case ID | ROADMAP 1008; archived chapter (NEEDS-AUTHOR-CONFIRMATION) | unconfirmed (M1) |
| 9 | Lenoir, Greene, Jones | Pilot | colleges CSV | yes |
| 10 | 58 colleges, 100 counties | several | county CSV | yes |
| 11 | 30 / 18 / 10, up to seven | Service areas | public list yes; clusters 31/18/9 | no (m1) |
| 12 | Bertie, Northampton shared with Roanoke-Chowan | Service areas | county CSV; chapter 437 | yes |
| 13 | Border or corner, legal boundaries across water | Service areas | AGENTS (cb = FALSE, queen) | yes |
| 14 | IGS 0.63 / 1.08 | Grid results | paired ratios 0.631 / 1.083 | yes |
| 15 | Isolation Buffer 0.74 / 3.80 | Grid results | 0.743 / 3.796 | yes |
| 16 | BQ 1.02 / 0.99; HIF 1.53 / 3.82 | Grid results | 1.023 / 0.995; 1.534 / 3.818 | yes |
| 17 | 2×2 2.68 / 1.47; Checkerboard 9.84 / 1.47 | Grid results | paired ratios | yes |
| 18 | SRS MSE 0.215 / 0.042; IGS 0.136 | Grid results | regime means 0.2146 / 0.0424; 0.1355 | yes |
| 19 | Coverage 94–95% | Grid results | regime means 0.942–0.951 | yes |
| 20 | 250 trials per scenario | Grid design | chapter 145 | yes |
| 21 | 100 allocations × 100 outcomes | Pipeline | chapter 1248 | yes |
| 22 | 7.6 million trials, no failures | Pipeline | chapter 528 (7,598,000) | yes |
| 23 | 15 to 58 colleges | Pipeline | chapter 516 | yes |
| 24 | ρ 0/0.2/0.5, γ 0.5/0.8 (six settings) | Pipeline, NC accuracy | heat-map note; chapter 544 | yes |
| 25 | Rates about 45 to 195 | NC incidence notes | AGENTS 45.3–194.9 | yes |
| 26 | Isolation Buffer 14 here, ~16 average | NC allocations | chapter 487 (15.89) | yes |
| 27 | Checkerboard 27; blocking 27–31 | NC allocations | chapter 484, 488 | yes |
| 28 | Reach IGS 53/45/48/52; HIF 22/24/23/22 | Budgets | chapter 491, 485 | yes |
| 29 | Isolation Buffer 27%; 2×2 50/55/48/50 | Budgets | chapter 487–488 | yes |
| 30 | IGS 21–30% lower (0.70–0.79) | NC accuracy | heat map | yes |
| 31 | SRS/BQ/IGS within ~3% under both arms | NC accuracy | heat map 0.97–1.01, BQ 0.98–1.00 | yes |
| 32 | Checkerboard ~3×; Isolation Buffer 3.7× | NC accuracy | heat map 2.83/2.99; 3.74 | yes |
| 33 | HIF 0.70–1.40; 0.92–0.95 | NC accuracy | heat map | yes |
| 34 | β = 1 HIF roughly fourfold | NC accuracy notes | chapter 630 (3.68–4.44) | yes |
| 35 | Worst 10% 21–43% lower; grid 0.19 vs 0.35 | Allocation risk | chapter 575; round-2 check | yes |
| 36 | 35–37 colleges, 1 SD, both arms | Colleges needed | tab_nc_power | yes |
| 37 | IGS 28 vs SRS 32 at 1.5 SD | Colleges needed | tab_nc_power | yes |
| 38 | 22 / 9 / 5 / ~1% lower at 58/40/29/20 | Partial participation | tab_nc_guidance IGS 0.78/0.91/0.95/0.99 | yes |
| 39 | ≤2% with incidence-weighted volunteering | Partial participation | tables_log 0.99–1.00 | yes |
| 40 | Power 0.71 / 0.96 at m = 25 | Effect size | chapter 652 | yes |
| 41 | Controls-only power 0.5–0.65 at 58 | Effect size, take-home | SRS 0.53, IGS 0.64 (chapter 652) | partly (m9) |
| 42 | "Within 8% on the grid and 3% in NC" | Take-home | 1.083; NC 0.97–1.01 | yes |
