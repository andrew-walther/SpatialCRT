# SUDDEN talk (2026-10-29): independent statistician review, round 4

**Date:** 2026-10-28
**Deck reviewed:** `presentation/sudden_talk_2026-10-29.qmd` (1,197 lines, slides and `{.notes}`). Figures opened: `cc_service_areas_cropped.png`, `fig_colleges_needed.png`.
**Sources checked against:** Gan et al. 2019 (JGIM 34:815–7, full text and Figure 1), `paper/SpatialCRT_IncidenceDesign.bib` (Mirzaei 2019 abstract, Habib 2026, Nanavati 2014), `CTJ_Manuscript.tex` 187, `Thesis_Chapter.qmd` 437 and 652, `application/AGENTS.md` (real-data section), `application/README.md` (open questions for Habib), `application/notes/SUD_Meeting_Notes.md`, `application/data/NC_Community_Colleges_by_County.csv`, `ROADMAP.md` 1011–1015, `paper/archive/pre_completion_20261002/Dissertation_Chapter.qmd` 407, `application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv` and `refined_tail_comparison.csv`, `results/srs_benchmark/manuscript_paired_ratios.csv`, `results/methods_ext_figures/tables/tab_nc_power.tex`, `tab_nc_guidance.tex`, `tables_log.txt`, `DECISIONS.md` (2026-10-28), round-3 statistician review.

**Verdict:** Blocking issues: 0 (2 major, 12 minor)

The r3 take-home and targeting majors (M2, M3) are fixed, and the case-identification slide now follows Gan et al.'s documented steps in the right order. Every results number I checked still traces to its source, and the new Gan figures (83,509 → 4,278; rural 100 vs urban 68 per 100,000 in 2014) match the paper. What remains is on the case-identification slide. Two parentheticals still name specific statewide exclusions that no repository source documents, against the 2026-10-28 decision to mention them only generically. The slide also says the algorithm was validated against SUDDEN adjudication of EMS, medical-examiner and clinical records, which neither manuscript nor Gan et al. states.

---

## r3 MAJOR fixes

| r3 item | Status | Evidence |
|---|---|---|
| M1 flowchart states undocumented upstream rules | **Partly fixed** | Steps now follow Gan Figure 1 (age → hospital, hospice or nursing home → suicide, homicide or accident → cancer). Side bullet and notes mention extra rules generically. But the boxes at lines 511 and 513 still name "overdose and COVID-19" and "other end-stage disease" (see M1 below). The 412,514 and 111,665 counts are still unconfirmed in any repository source except ROADMAP 1013 (minor 1). |
| M2 take-home overgeneralizes, omits Balanced Quartiles | **Fixed** | Line 812 is now "gained the most under controls-only spillover, at a small cost under both arms". Line 827 names Balanced Quartiles. Notes 834 state the trade. One wording point remains (minor 4). |
| M3 targeting claim without NC qualifier | **Fixed** | Card 818 says "(a small shift on the NC map, where regions are similar)". Notes 685 and 837 carry the qualifier too. |

---

## BLOCKING

None.

---

## MAJOR

### M1. The flowchart still names specific statewide exclusions that no source documents
- **Slide:** "SUD deaths were identified from death certificates with a validated algorithm" (lines 511, 513).
- **Text:** "Exclude suicide, homicide and accidents (statewide data also exclude overdose and COVID-19)" and "Exclude expected deaths from cancer (statewide data also exclude other end-stage disease)".
- **Problem:** DECISIONS 2026-10-28 says the statewide data's extra exclusions are "mentioned only generically". The side bullet (521) and notes (532) do that, but the two box parentheticals don't. Their content comes from the archived chapter's exclusion list, which was tagged `NEEDS-AUTHOR-CONFIRMATION` and called undocumented. `application/README.md` lists "Should the free-text filter match KIDNEY / LIVER FAILURE?" as still open, so "other end-stage disease" may not be applied as the box says. Overdose deaths are mostly certified with manner "accident", so "also exclude overdose" reads as redundant to clinicians and invites a question the speaker can't answer. Habib and Dr. Simpson will be in the room.
- **Fix:** Drop both parentheticals and keep Gan's steps as they are: "Exclude suicide, homicide and accidents" and "Exclude deaths from cancer". Add one box note or a side-bullet clause such as "Statewide counts apply a few further exclusion rules (Habib 2026)", which the slide already has at 521. Another option is to have Habib confirm the two lists before the talk.

### M2. "Validated against SUDDEN physician adjudication of EMS, medical examiner and clinical records" has no source
- **Slide:** Case identification, title (line 502) and side bullet (519). Notes 528: "validated against SUDDEN's physician adjudication".
- **Problem:** Gan et al. 2019 describes a stepwise screen of death certificates and reports no validation. The CTJ (187) and chapter (437) say only "the case-adjudication algorithm of Gan et al." The archived chapter flagged the validation citation as "unverifiable from the bib entry or the literature review". The "EMS, medical examiner and clinical records" detail appears nowhere in the repository. It may come from SUDDEN methodology (Nanavati 2014 screens EMS-attended deaths), but that paper validates the SUDDEN screening, not this death-certificate algorithm. The claim describes the audience's own method, so an error here is costly. If Dr. Simpson confirms it, downgrade this to minor.
- **Fix:** Title: "SUD deaths were identified from death certificates with the SUDDEN group's algorithm". Bullet: "The algorithm (Gan et al. 2019) screens certificates for premature, unexpected, natural out-of-hospital deaths, following the SUDDEN case definition (Nanavati et al. 2014)". Notes 528: "...with the algorithm from Gan and colleagues, built on SUDDEN's case definition." Or confirm the validation claim and its citation with Dr. Simpson before the talk.

---

## MINOR

1. **412,514 and 111,665 are still unconfirmed** (Case ID, 507–509; notes 530). The only source is ROADMAP 1013, which attributes them to Habib (2026). The archived chapter tagged them `NEEDS-AUTHOR-CONFIRMATION` and nothing since has confirmed them. They are upstream of the heart-failure change, so they are valid for 23,523 if they are right. Fix: ask Habib to confirm them, and add "(counts, Habib 2026)" to the side bullet.
2. **"By county of residence"** (Case ID, 521; notes 536). The meeting notes (2026-05-08) say deaths were pooled to counties by death-certificate zip code, and the chapter declines to infer "zip-code allocation". Fix: "tallied by county and year", or confirm with Habib.
3. **Burden slide source line misses Habib** (line 88 area; the NC column at 94–95). 23,523 and 84 → 101 come from the corrected counts, but the source line cites only Mirzaei and Gan. Fix: add a source line under the NC column: "Corrected county counts, A. Habib (Habib 2026, with heart-failure deaths kept); Gan et al. 2019 for 2014."
4. **"Balanced Quartiles matched SRS in both regimes but never beat it"** (Take-home, 827; notes 834). This holds on the grid (1.02 / 0.99). In NC under controls only, Balanced Quartiles was below SRS in all four years (0.92, 0.94, 0.95, 0.995; `yearly_primary_design_means.csv`), and in the allocation-risk table its mean was 0.95. Fix: "Balanced Quartiles stayed close to SRS in both regimes (within 2% on the grid, up to 8% lower in NC)".
5. **Reason given for the rate gap** (Burden notes, 110). "...which used a stricter physician review, different years and a different region." This is an improvement on r3, but the reason is still asserted without a source, and the ages differ too (20–64 against 18–64). Fix: "...partly because physician review is stricter than a certificate algorithm, and partly because the years, ages and region differ."
6. **Prose colon on the burden card** (line 96). "Higher in rural counties: in 2014, about 100 per 100,000..." Fix: "Higher in rural counties, about 100 per 100,000 in highly rural NC counties against 68 in highly urban ones in 2014 (Gan et al.)". Gan reports means of county rates, so "on average" is also accurate.
7. **"Neighboring colleges share staff and trainings"** (Pilot, 135; r3 m5 carry-over). The meeting notes support "talk to each other" and "satellite campuses", not shared staff. Fix: "neighboring colleges talk to each other and share satellite campuses, so training can spread".
8. **Pilot notes state the trainees as fact** (150). "The trainees are students and staff in those programs", but the slide box says "Trainees (to confirm)". Also, the flooding spike is recorded for Lenoir County, while the slide (120) attaches it to all three counties. Fix: "The trainees would be students and staff in those programs..." and "...(Lenoir, Greene and Jones counties); Lenoir County's SUD rates spike after flooding".
9. **Lower-case "saturation" as a design name in the notes** (492, 703, 705; r3 m10 and m14 carry-over). Line 703, "Saturation's whole bar sits to the left", is also inaccurate, since the NC thin line reaches 1.02 (r3 m14). Fix: "Incidence-Guided Saturation's thick bar sits well to the left, and its tail is much shorter than random sampling's"; and "Incidence-Guided Saturation is about 8 percent worse".
10. **Prose colon in backup notes** (925, r3 m11 carry-over). "...decide the design and the sample size: the spillover regime..." Fix: "...decide the design and the sample size, namely the spillover regime...".
11. **Outline wording** (66). "Design intervention estimation efficacy on a simulated 100 cluster region" is jargon, and "100 cluster" needs a hyphen. Fix: "How accurately each design estimates the effect on a simulated 100-cluster region".
12. **"Shared borders" in limitation notes** (874), while neighbors are border-or-corner everywhere else. Fix: "used neighboring service areas as a stand-in...". Also, the colleges-needed fig-note (712) still doesn't say the figure shows the two backup variants, and the Isolation Buffer label "> 58 in 1/4 yr, else 58+" reads oddly (r3 note). Add "(includes the two backup variants)".

---

## Notation, framing and style: checked

- **Notation.** Fully defined. Y, Z, W (queen, row-standardized), τ (direct effect before propagation, in noise SD), ρ, S(Z) in both forms, (I − ρW)⁻¹ with I, γ, β (grid 1; NC 0 primary, 1 hypothetical), X (grid and NC), ε, ∘, m(a), σ_c, ICC, m, d. Backup: N and p in the corrected interval, ρ_X glossed as incidence dependence.
- **Overclaiming.** The take-home now states the trade and the NC targeting qualifier and keeps "the right choice still depends on the spillover we expect", which is consistent with the 2026-10-14 and 2026-10-28 rules. The remaining overclaims are on case identification (M1, M2) and the "never beat it" wording (m4).
- **Style.** No DIM, Claude or Codex. No sentence starts with "Because". Design names on slides match the AGENTS set exactly. "Saturation Quadrants" appears only in backup and fig-notes. The figure alt text "four saturation regions" is generic, not a design name. Prose colons remain at m6 and m10. Label colons ("Both arms:", "Grid:", "Share:", "Variant:") are fine.
- **Timing.** Running totals add up to 22 minutes, within the "about 20 ± a few" in DECISIONS.

---

## Number spot-checks (30)

| # | Deck statement | Slide | Source | Match |
|--:|---|---|---|:-:|
| 1 | 83,509 certificates → 4,278 SUD, 2014 | Case ID | Gan Fig. 1 | yes |
| 2 | Step order: age, location (hospital, hospice, nursing home), manner (suicide, homicide, accident), cancer | Case ID | Gan Fig. 1 | yes ("expected" before cancer is the deck's gloss) |
| 3 | Rural 100 vs urban 68 per 100,000, 2014 | Burden | Gan Results (county-mean incidence) | yes |
| 4 | SUD definition: unexpected natural out-of-hospital death, working age | Burden | Gan Methods | yes |
| 5 | 49 / 22 per 100,000, ages 20–64, Wake County extrapolated | Burden | Mirzaei abstract 49.3 / 21.7 | yes |
| 6 | 2 million YLL, 10%, $51 billion lifetime | Burden | Mirzaei abstract | yes |
| 7 | 23,523; 84 → 101 rising each year | Burden | AGENTS 83.9/85.9/97.1/100.5 | yes |
| 8 | 25.6 million person-years | Case ID | 25,594,321 | yes |
| 9 | 21,147 when heart-failure deaths excluded | Case ID | AGENTS, README | yes |
| 10 | 412,514 / 111,665 | Case ID | ROADMAP 1013 only | unconfirmed (m1) |
| 11 | Validated against SUDDEN adjudication of EMS, ME, clinical records | Case ID | none | no (M2) |
| 12 | Overdose, COVID-19, other end-stage disease excluded | Case ID | archived NEEDS-CONFIRMATION; README open item | no (M1) |
| 13 | 31 / 18 / 9 colleges serve one / two / three to seven | Service areas | county CSV: 31, 18, 6+2+1, max 7 | yes |
| 14 | Bertie, Northampton secondary to Roanoke-Chowan | Service areas | county CSV | yes |
| 15 | Lenoir CC is number 29 on the map | Service-area notes | figure legend | yes |
| 16 | Lenoir serves Lenoir, Greene, Jones | Pilot | county CSV | yes |
| 17 | Process measures: medicine, timely clinic visits; credit; floods | Pilot | meeting notes 2026-04-27 | yes |
| 18 | IGS 0.63 / 1.08; BQ 1.02 / 0.99 (grey) | Grid results | paired ratios 0.631 / 1.083; BQ CI 0.995–1.050 and 0.966–1.023 | yes |
| 19 | NC IGS controls only 21–30% lower | NC accuracy | 0.701–0.789 | yes |
| 20 | NC both arms SRS / BQ / IGS within ~3% | NC accuracy | IGS 0.972–1.006; BQ 0.976–0.998 | yes |
| 21 | Checkerboard ~3×; Isolation Buffer 3.7× both arms | NC accuracy | 2.83 / 2.99; 3.74 | yes |
| 22 | HIF 0.70–1.40 controls only | NC accuracy | 0.697–1.398 | yes |
| 23 | Reach IGS 53/45/48/52; HIF 22/24/23/22; Isolation Buffer 27; 2×2 50/55/48/50 with 28/31/27/29 treated | Budgets | pop shares, Mean_Treated | yes |
| 24 | Worst 10% 21–43% lower in NC | Allocation risk | IGS/SRS worst-decile 0.566–0.793 over 8 settings | yes |
| 25 | 35–37 colleges at 1 SD both arms; IGS 28 vs SRS 32 at 1.5 SD | Colleges needed | tab_nc_power (HIF 35, 2×2 37; 28 / 32) | yes |
| 26 | No design reaches 80% at 1 SD controls only with 58 | Colleges needed | MDE at 58: IGS 1.25, SRS 1.45, all > 1 | yes |
| 27 | 22 / 9 / 5 / ~1% lower at 58 / 40 / 29 / 20 | Partial participation | tab_nc_guidance 0.78 / 0.91 / 0.95 / 0.99 | yes |
| 28 | ≤ 2% with incidence-weighted volunteering | Partial participation | tables_log IGS 0.99–1.00 | yes |
| 29 | σ_c ≈ 0.30; power 0.51/0.82, 0.71/0.96, 0.80/0.98 | Effect size | √(0.05 + 0.95/25) = 0.297; translation table | yes |
| 30 | Controls-only power with 58: 0.53 (SRS) to 0.64 (IGS) | Effect size, take-home | chapter 652 | yes |
