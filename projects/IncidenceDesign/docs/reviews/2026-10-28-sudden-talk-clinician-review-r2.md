# SUDDEN talk (2026-10-29) clinician/epidemiologist review, round 2

**Date:** 2026-10-28
**Deck:** `presentation/sudden_talk_2026-10-29.qmd` (revised after round 1), figures in `presentation/figures/`
**Reviewer stance:** a clinician or epidemiologist on the SUDDEN team, new to this deck, comfortable with trial design but not spatial statistics. Read-only review; the deck was not edited.

**Blocking issues: 0**

(2 major, 11 minor)

The revision fixes the round-1 problems. The script now fits comfortably in 20 minutes. "Participating colleges" is defined, the bullet columns match the figure, and the rollout, rate instability and "why not treat the highest-need areas" points are all in the script. The questions now ask for the endpoint and the effect size. What remains is mostly about the practical message the group will take away. The deck tells them about 36 colleges are enough, but it does not say clearly enough that under classic contamination even all 58 colleges are underpowered for the same effect. The outcome definition ("trained" professionals) also does not work for control colleges.

---

## Round-1 clinician findings: status

| Finding | Status | Evidence in the revised deck |
|---|---|---|
| B1 Script runs 23–26 min | **Fixed** | 2,520 words (see Timing), 16.8–19.4 min; headers total 18.5 min |
| M1 Equation slide loses the room | **Partly fixed (author decision honoured)** | Slide kept by author request. Script cut to 166 words in 1.25 min (133 wpm) and reads the model in words. See minor m1 for what is still hard |
| M2 "Colleges needed" ambiguous | **Fixed** | "Counts are **participating** colleges, about half treated"; table headers "29 participating colleges"; "Participating colleges" column |
| M3 "36 colleges" / "one SD" lose their assumptions | **Fixed** | Take-home card adds "about 0.3 SD per professional with 25 trained per college), under our assumptions"; fig-note states 80% spillover in words; script previews the SD translation |
| M4 Partial-participation slide dense, counterintuitive | **Mostly fixed** | Script explains the fall near 58; fig-note reconciles the 22% and 21–30% figures. Four-panel figure kept (acceptable) |
| M5 Bullet columns reversed against the heat map | **Fixed** | Left column now "Both arms", matching the figure |
| M6 "Why not the highest-need areas?" unanswered | **Fixed** | NC accuracy script, "Treating the highest-rate colleges is the intuitive choice, but it is not randomized..." |
| M7 "First" implies rollout, not discussed | **Fixed** | Planning-slide script on the two-arm, one-time-point model; delayed start in question 3 and Next steps |
| M8 Questions omit endpoint and effect size | **Fixed** | Question 2 asks for outcome, timing and a worthwhile improvement; ICC asked through pilot scores on the effect-size slide |
| M9 Small-area rate instability | **Fixed** | NC incidence script, "Part of that is small numbers..." |

Round-1 minor items: fixed 1 (card only, title unchanged), 3 (partly), 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17. Not fixed: 8 (no specific college named on the yearly maps), 16 (figure labels, see m6).

---

## Timing

Words in the main-body `.notes` blocks, title through "Thank you!", with the bracketed time headers excluded: **2,520**. That is **19.4 min at 130 wpm** and **16.8 min at 150 wpm**, against the headers' 18.5 min, which leaves 1–3 minutes of slack inside 20. No slide is over-packed. The densest are budgets-and-reach (83 words in 0.5 min, 166 wpm) and allocation risk (115 in 0.75, 153 wpm), and both are fine. No timing change is needed. The slack can absorb the two MAJOR fixes below, which add about 40 words.

| Slide | Words | Budget (min) | wpm |
|---|---:|---:|---:|
| Title | 94 | 0.75 | 125 |
| Spillover examples | 94 | 0.75 | 125 |
| Design gap | 72 | 0.5 | 144 |
| NC question | 129 | 1.0 | 129 |
| Two spillover regimes | 125 | 1.0 | 125 |
| Outcome model | 166 | 1.25 | 133 |
| Performance measures | 102 | 0.75 | 136 |
| The designs | 138 | 1.0 | 138 |
| Grid study design | 66 | 0.5 | 132 |
| Grid results | 148 | 1.0 | 148 |
| NC incidence | 146 | 1.0 | 146 |
| NC allocations | 93 | 0.75 | 124 |
| Budgets and reach | 83 | 0.5 | 166 |
| NC accuracy | 159 | 1.25 | 127 |
| Allocation risk | 115 | 0.75 | 153 |
| Colleges needed | 160 | 1.25 | 128 |
| Partial participation | 114 | 0.75 | 152 |
| Effect size in practice | 109 | 0.75 | 145 |
| Targeting burden | 92 | 0.75 | 123 |
| Take-home | 102 | 0.75 | 136 |
| Limitations | 102 | 0.75 | 136 |
| Questions | 79 | 0.5 | 158 |
| Thank you | 32 | 0.25 | 128 |

---

## BLOCKING

None.

---

## MAJOR

### M1. The take-home leaves out that under contamination even all 58 colleges are underpowered (Colleges needed; Effect size; Take-home)

**Text.** Colleges-needed title: "About 36 participating colleges suffice if spillover reaches every college". Take-home, controls-only card: "**The saturation designs** ... gave the lowest or near-lowest error and a lower risk of an unlucky allocation..." (nothing on power). Navy card: "about **36 participating colleges** detect a one-SD college effect".

**Problem.** The single number the group will carry away is "about 36 colleges". The spillover slide calls controls-only spillover "classic contamination", and for an education program spread through colleagues many in this room will think it is the more likely regime. Under that regime the deck's own figures say that no design reaches 80% power for a one-SD (about 0.3 SD per professional) effect even with all 58 colleges. `fig_power_vs_n.png` shows about 0.52 for SRS and 0.63 for the saturation designs at 58, and the backup translation gives 0.53 for SRS at m = 25, ICC = 0.05, d = 0.3. This appears only in one bullet on the colleges-needed slide and in "Power is lower under controls only (backup)". It is the most important planning fact for the pilot. It means that if contamination is expected, the trial needs a larger effect, more professionals per college, or an outcome with less between-college variation, and that the design choice alone will not fix it. Without it the group may plan a 36-college trial that cannot detect the effect it is designed for.

**Fix.**
- Add to the controls-only take-home card: "Even with all 58 colleges, power for a one-SD effect is only about 50–65%; detecting it needs a larger effect or more professionals trained per college".
- Script, take-home, after the saturation sentence: "But under this kind of spillover, even all 58 colleges give only about a 50 to 65 percent chance of detecting a one-SD effect, so the size of the effect and how many professionals each college trains matter as much as the design."
- Effect-size slide script, replace "Under controls only it is lower." with "Under contamination it drops to about 50 percent even with all 58 colleges."
- Optionally retitle the colleges-needed slide "About 36 participating colleges suffice under spillover to every college; under contamination, 58 are not enough for a one-SD effect".

### M2. "Trained healthcare professionals" cannot be the outcome population in control colleges (Planning slide; Effect-size slide)

**Text.** Planning slide: "The outcome is an **education score** of trained healthcare professionals, averaged by college". Effect-size script: "The college score averages its trained professionals". Table header "Professionals per college", with the question "how many professionals each college would train".

**Problem.** Control colleges do not deliver the program before the outcome is measured (the planning script says so). So whose score is measured there? A clinician will spot this at once, and it is exactly the endpoint question (question 2) the talk wants the group to answer. As worded, the outcome sounds like it exists only in treated colleges, and that undercuts the whole two-arm comparison. It also blurs what m means in the power translation. Is it professionals trained, or professionals assessed?

**Fix.**
- Planning slide bullet: "The outcome is an **education score** of eligible healthcare professionals in each service area (trained or not), averaged by college".
- Effect-size script: "The college score averages the professionals assessed at each college" and "how many professionals each college would assess". In question 3 say "how many professionals each would train and test".
- Be ready to say that, in the model, control colleges are assessed at the same time with no training yet.

---

## MINOR

**m1. Outcome model slide (author kept it; digestibility).** The 166-word script reads the model in words and is at a natural pace, which is a clear improvement. Two things still trip a clinical listener.
(a) There is no anchor. The audience reads an 11-row glossary while listening. Open with "You don't need the algebra. Three pieces matter: the program's own effect, spillover from treated neighbors, and neighbors' scores moving together."
(b) "Tau ... is the effect of a college's own treatment before that spreading, not the total effect" invites "so what effect would a college actually see?" with no answer. Either add "the total effect a college sees is larger, because its neighbors' gains feed back into its score" or drop the clause from the spoken script and leave it in the glossary.
The τ glossary row ("before outcomes spread between neighbors through $(I-\rho W)^{-1}$ (not a total effect), in standard deviations of the noise $\epsilon$") is the densest line on the slide. It could be shortened to "Effect of a cluster's own treatment, in SD units, before spread to neighbors; what we estimate".

**m2. Grid study script, "I'll show an effect of one."** A clinician hears "an effect of one" with no unit. Say "I'll show a treatment effect of one standard deviation."

**m3. Grid study script, "Pooled over the two regimes, the ordering was the same at every effect size, and the one change at large effects is in the backup slides."** This contradicts itself ("same at every size" and "one change") and is a statistician's point. Cut the sentence. The backup slide covers it if asked.

**m4. Performance measures slide.** "(a spatial lag model fitted by maximum likelihood, including the true spillover term)" and the $m(a)$ notation remain. "True spillover term" sounds as if the real spillover is known, which the limitations slide then walks back. Suggested on-slide sentence: "We estimate the program's effect $\tau$ with a spatial model that allows for spillover, and repeat the trial thousands of times in simulation." Allocation-risk card: "A trial gets **one** allocation. We report how bad the worst 10% of allocations are."

**m5. Saturation Quadrants and Balanced Halves appear without introduction.** The designs slide shows six designs plus SRS. Saturation Quadrants first appears on the targeting slide and the take-home ("the saturation designs (Incidence-Guided Saturation and Saturation Quadrants)"). Both variants also appear in `fig_colleges_needed.png` and `fig_power_vs_n.png`. Add one sentence to the designs script: "Two variants, Saturation Quadrants, the same rule with regions in random order, and Balanced Halves, behave like their parents and appear in some later charts."

**m6. Figure labels (round-1 minor 16, not fixed).** `fig_colleges_needed.png` prints "2x2 Blocking" (letter x), where the deck uses "2×2 Blocking". The power and colleges-needed figures say "Both-arms spillover"/"Control-only spillover", the other figures say "Spillover reaches both arms/controls only", and the slide text says "both arms"/"controls only". These are small, but this audience reads labels literally. The grey-bar labels "> 58 in 1/4 yr, else 44+" are cryptic. In the script, say "grey means 58 colleges were not enough in at least one year".

**m7. NC incidence slide, Panel B.** The regions are numbered 1–4, but nothing says which gets 80% treated. Since that is the point of the panel for this audience, say in the script which region is highest in a typical year, or label the regions by name (west, northwest, northeast, southeast).

**m8. Allocation-risk figure colour.** SRS is the orange/amber bar, while on the heat maps and the grid table orange means "worse than SRS". The script's "Random sampling, the orange bar" is fine, but a viewer carrying the earlier colour key may misread it. Either recolour SRS grey in the figure or say "the reference bar in orange" once.

**m9. NC accuracy script, "its error rose up to fourfold".** The backup says 3.7–4.4-fold. Say "about fourfold" or "up to about four and a half times".

**m10. Design-gap title (round-1 minor 1, partly fixed).** "Analysis methods exist, but guidance on which clusters to treat is thin" is still statistician phrasing. Suggested: "We know how to analyze spillover, but not how best to choose who is treated".

**m11. Backup colon in prose.** Rook backup: "Under rook, Checkerboard gives $WZ = 1 - Z$: exposure equals $1 - Z$, so it cannot be separated..." → "Under rook, Checkerboard gives $WZ = 1 - Z$, so exposure equals $1 - Z$ and cannot be separated...". The main-body slides use colons only as bullet labels ("**Both arms:**"), which is acceptable. No colons were found in the scripts. No sentence starts with "Because".

---

## Framing and practical meaning (checked, no further change needed beyond M1–M2)

- Education outcome, not deaths: stated on the planning slide (amber card and script), on the NC incidence slide and in the β glossary row. Clear.
- Design evaluation, not program evidence: stated in the opening script and as the last limitation. Clear.
- Maps, budgets and reach: the budgets-and-reach table and the rates-versus-deaths contrast (62% of rate share against 35% of deaths) remain the most useful planning facts. The NC allocation maps are readable, and the fig-note now reconciles 14 against about 16 for Isolation Buffer.
- Effect size: the translation to about 0.3 SD per professional is clear and appropriately hedged. M1 asks that the controls-only consequence be carried into the take-home.
- Design names are consistent with the canonical set on slides and in tables. "Spatial blocking" for 2×2 Blocking on the NC map is explained on the allocations slide and labelled "2×2 Blocking (spatial blocks)" in the budgets table, which is acceptable.
- The IRB protocol number on the thank-you slide matches `DECISIONS.md`.
