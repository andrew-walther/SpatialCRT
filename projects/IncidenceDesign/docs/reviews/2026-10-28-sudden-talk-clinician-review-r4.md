# SUDDEN talk (2026-10-29) clinician/epidemiologist review, round 4

**Date:** 2026-10-28
**Deck:** `presentation/sudden_talk_2026-10-29.qmd`, figures in `presentation/figures/` (all main-body PNGs viewed)
**Background checked against:** `application/notes/SUD_Meeting_Notes.md`, `application/AGENTS.md`, `DECISIONS.md` (2026-10-28 talk entry), `application/data/NC_Colleges_and_Service_Areas.csv`, CTJ Methods/Results and Chapter 3 (non-oracle claim), round 3 review
**Reviewer stance:** a clinician or epidemiologist on the SUDDEN team, new to this deck. Read-only review; the deck was not edited.

**Blocking issues: 0**

(2 major, 14 minor)

The deck is in good shape for this room. All four round-3 MAJORs are fixed. The pilot slide now sticks to what Dr. Simpson said and flags the trainees as "to confirm". Process measures are credited to him, and the education score is clearly our assumption. The case-identification chart separates non-natural from expected deaths, and both application slides and the take-home now say the targeting gain on the NC map is small. Two things remain at MAJOR level. The script has grown by about 150 words since round 3 and now runs past the window at a slow pace. The new outline slide uses the word "efficacy", which this audience hears as "does the program work".

---

## Round-3 MAJOR findings: status

| Finding | Status | Evidence |
|---|---|---|
| M1 Pilot slide says trainees are EMS/healthcare professionals | **Fixed** | Flow box "Trainees (to confirm) Students and staff in those programs, possibly including EMS and healthcare trainees; students earn credit..."; script "For our model we assumed they include EMS and other healthcare trainees, and that is one of the points I'd most like clarified." Matches the notes (students/staff, credit, outreach beyond transport). Residual wording drift later in the deck, see m1 |
| M2 Education score attributed to the pilot | **Fixed** | Bullet "...process measures (taking medicine, timely clinic visits), not SUD rates; for the design study we assume a trainee education score"; script "Dr. Simpson expects early effects in process measures..." Slide bullet still states it as fact, see m2 |
| M3 "Expected deaths" category error; undocumented steps | **Fixed (category); partly open (provenance)** | Now "Exclude suicide, homicide and accidents" → "Natural out-of-hospital deaths" → "Exclude expected deaths from cancer". Gan comparison numbers and two definitional details still worth confirming with Ashkan, see m6 |
| M4 Burden-targeting overstated | **Fixed** | NC accuracy script "...though on this map the regions are similar, so the shift is small"; take-home card "(a small shift on the NC map, where regions are similar)" and script "a small shift on this map" |

Round-3 minors: m2 take-home sentence split (fixed), m3 "patients" removed (fixed), m6 neighbor definition (fixed, "share a border or a corner" in pipeline step 2), m7 heart-failure note moved to the side bullet (fixed), m10 "(corrected counts)" (fixed), m11 flooding clause (added). Still open: m1 timing cuts (now MAJOR, see M1), m4 statewide-plan attribution (largely fixed by "After the pilot, the plan is to share..."), m5 Lenoir not marked on the map (m8 here), m8 shorthand names (m11 here), m9 variants and figure labels (m12 here).

---

## Timing

Main-body `.notes` words, title through "Thank you!", time headers and bullet markers excluded: **3,044** (round 3: 2,892).

| Pace | Minutes |
|---|---:|
| 150 wpm | 20.3 |
| 140 wpm | 21.7 |
| 130 wpm | 23.4 |

Time headers sum to 22.0 minutes. Slides budgeted faster than read-aloud pace: Outline (41 words / 0.25 min, 164 wpm), Budgets and reach (83 / 0.5, 166), Title (77 / 0.5, 154), Partial participation (114 / 0.75, 152), Pilot (225 / 1.5, 150).

---

## BLOCKING

None.

---

## MAJOR

### M1. The script now runs 23.4 minutes at 130 wpm, outside the ±2-minute window (all main-body slides)

**Problem.** The additions since round 3 (outline, longer pilot script, M4 caveats) brought the script to 3,044 words. At 140 wpm it is 21.7 minutes, just inside the window. At 130 wpm, a likely pace for a careful first delivery to a non-statistical room, it is 23.4. The headers (22.0) already exceed the "~21 minutes" target and understate five slides. Overrunning eats into the 40-minute discussion that is the point of the meeting. None of the round-3 cuts were taken.

**Fix.** Cut about 200 words to reach about 2,850 (20.4 min at 140 wpm, 21.9 at 130). None of these lose content, since each repeats the slide or another line.
- Pilot: drop "Here is how the application came about." (7). Merge "Colleges talk to each other, so the training can spill over." into the preceding sentence as "..., and since colleges talk to each other, the training can spill over" (saves 3; avoid starting with "Because").
- Outline: replace the 37-word sentence with "I'll cover the background and pilot, our methods, a short simulation section, and then the North Carolina application, which takes the most time." (saves 13).
- Burden: drop "These death-certificate rates run higher than the adjudicated national figure, which used a stricter physician review, different years and a different region." (24) and keep it ready for questions. Ashkan and Dr. Simpson know why the numbers differ; if the author wants it said, keep it and cut elsewhere.
- Case identification: drop "Ashkan tallied them by county of residence and year, and we used SEER population estimates for the denominators." (19, on the slide).
- Outcome model: drop "W records who neighbors whom, and W times Z is the share of a college's neighbors that are treated." (19, in the glossary).
- Grid results: drop "Each number is a design's error divided by that of simple random sampling, so blue, below one, is better than random and orange is worse." (26). The performance-measures slide already explained ratios; say instead "Blue is better than random, orange worse." (7, saves 19).
- Colleges needed: drop "These counts assume moderate spatial dependence, spillover at 80 percent of the effect and colleges volunteering at random, and they use a test that holds the false-positive rate at 5 percent." (31, in the fig-note; answer if asked).
- Effect size: drop "The table shows power for that individual effect with random sampling under spillover to both arms." (17).
- Take-home: drop "To pull this together, the goal is a design that does reasonably well whichever kind of spillover is true, with limited downside." (22, the slide title says it).
- Then raise the Outline and Budgets headers to 0.3 and 0.6 so the running totals match the read-aloud time.

### M2. The outline calls the simulation results "Design intervention estimation efficacy" (Outline)

**Text.** "Simulation Results — Design intervention estimation efficacy on a simulated 100 cluster region".

**Problem.** To clinicians and epidemiologists "intervention efficacy" means whether the program works. The talk's core framing, repeated on the title script, the service-area card and the last limitation, is that this is a design evaluation and not evidence of effect. The outline is the first content slide after the title, so the word lands before the framing has been set. The phrase also stacks four nouns and is hard to parse. I recognise the outline wording is author-specified, so the fix keeps its structure and changes only the noun.

**Fix.** "How accurately each design estimates the intervention effect, on a simulated 100-cluster region". If the author prefers to keep his phrasing, at minimum replace "efficacy" with "accuracy" and hyphenate "100-cluster".

---

## MINOR

**m1. "Trainees" on the pilot slide become "professionals" everywhere after it (Service-area script; Effect-size slide; Budgets script; Limitations).** The pilot slide now says the trainees are students and staff, possibly including EMS and healthcare trainees. From the service-area slide on, the deck says "the professionals assessed in each college", "A one-SD college effect is a modest effect for each trained professional", "how many professionals each college trains" and "how professionals actually communicate". A clinician will notice students turning into professionals two slides later. The working assumption (DECISIONS 2026-10-28, EMS and healthcare professionals) is fine, but bridge it once. Service-area script, change to "an education score of the trainees assessed in each college, which in our model are healthcare and EMS professionals, averaged within the college". Elsewhere "trainee" or "person trained" works without changing any number.

**m2. Process measures stated as fact on the slide (Pilot).** Bullet "Early effects show in **process measures**..." The notes say Dr. Simpson expects changes in medication adherence and timely clinic visits; nothing has been observed yet. Change to "Early effects are expected in **process measures**...". The script already says "expects".

**m3. Flooding claim presented as data (Pilot).** Flow box "...where SUD rates spike after flooding". This comes from Dr. Simpson's remark in the meeting, not from an analysis in the deck. Script "Lenoir County floods regularly, and SUD rates spike after flooding" is fine said aloud by him, less so asserted to his team as a finding. Use "where Dr. Simpson has seen SUD rise after floods" on the box, and "Dr. Simpson tells us SUD rises after floods there" in the script.

**m4. "neighboring colleges share staff and trainings" is our assumption (Pilot, right bullet).** The notes say only "schools talk to each other". The flow box already uses that wording ("colleges talk to each other"). Change the bullet to "...and neighboring colleges talk to each other, so training can spill over", or mark it as our assumption.

**m5. Pilot flow still mixes steps and attributes (Pilot).** "Trainees (to confirm)" is a property of the Lenoir pilot, not a step after it (round-3 m3). Low priority now that the content is right; folding the trainee box into the Lenoir box would also free width for the five boxes.

**m6. Case-identification details to confirm with Ashkan before the talk (Case identification).**
- The comparison "Gan et al. (2014 data, 83,509 certificates → 4,278 SUD)" is not recorded anywhere in the project (no hit in ROADMAP, the bib, the chapter or the notes). Add the page or table of Gan et al. to the speaker note so it can be defended if asked.
- "statewide data also exclude overdose" sits beside "accidents". Unintentional drug overdoses are usually coded as accidental manner of death (ICD-10 X40–X44), so the team may ask what the separate overdose rule adds. A word from Ashkan on whether it targets intentional/undetermined poisonings would settle it.
- "Exclude deaths in a hospital, hospice or nursing home". Out-of-hospital SUD definitions often keep emergency-department and dead-on-arrival deaths and exclude only inpatient deaths. I could not check Gan's coding from the project files. If Gan's screen works that way, say "hospital inpatient".
- `application/AGENTS.md` notes that Habib's free-text regex cannot match "KIDNEY FAILURE" or "LIVER FAILURE", so "statewide data also exclude other end-stage disease" is only partly applied. That is fine with "a few further exclusion rules" in the script; just be ready for the question.
- The DECISIONS entry says the statewide extra exclusions should be "mentioned only generically", but the slide names overdose, COVID-19, other end-stage disease and the heart-failure change. Either is defensible; check that this level of detail is intended.

**m7. Long sentences to read aloud.** Eleven script sentences run 32–37 words: the outline sentence (37), the NC burden sentence (36), "Dr. Simpson expects early effects..." (36), "When almost everyone is in..." (35), the outcome-model Y sentence and the "two spillover forms" sentence (34 each), "We first tested the designs on this grid..." (34), the Gan rural/urban sentence (34), "Statisticians have good tools..." (33), the case-identification exclusion sentence (33), "No design reaches 80 percent power..." (33) and "In the table, Incidence-Guided Saturation's edge..." (32). Each splits at its first ", and" or ", so". The M1 cuts remove two of them.

**m8. Lenoir CC is not marked on the service-area map (Service-area clusters).** The script says "The pilot site, Lenoir Community College, is number 29 on the map." Number 29 is a small label among similar greens (27, 28, 3 and 30 are all green). A bold outline or call-out on 29 would let the room find it at once. The legend still sorts as text (1, 10, 11 ... 2, 20), round-3 m5.

**m9. Allocation-risk script overstates SRS's tail (Allocation risk).** "Random sampling ... has a long right tail, so some random allocations are much worse than average" and "Saturation's whole bar sits to the left." In the NC controls-only panel, SRS's 95th percentile is about 1.5 times its mean, and Incidence-Guided Saturation's whisker reaches about 1.0, overlapping SRS's box. Say "Random sampling's bar reaches further right, so its unlucky allocations are worse" and "Incidence-Guided Saturation's bar sits further left."

**m10. Checkerboard power loss understated (Partial participation script).** "...which is why random sampling and Checkerboard even lose a little power near 58." In `fig_power_vs_n.png` under controls only, Checkerboard falls from about 0.35 to 0.17 at one SD and from 0.65 to 0.34 at 1.5 SD. Say "random sampling loses a little power near 58, and Checkerboard loses half of it."

**m11. Shorthand design names in the script.** "saturation is about 8 percent worse" (grid), "Balanced Quartiles and saturation are within about 3 percent", "saturation and random sampling differ little" and "Saturation's whole bar" (allocation risk), "for the saturation design" (pipeline). The take-home script still says "fixed targeting swung from year to year" without naming High Incidence Focus. Once Saturation Quadrants appears in the power figures, "saturation" is ambiguous. Use "Incidence-Guided Saturation" and "High Incidence Focus".

**m12. Variants and figure labels (Colleges needed; Partial participation).** DECISIONS 2026-10-28 says Saturation Quadrants stays in the backup slides, but it appears in `fig_colleges_needed.png` and `fig_power_vs_n.png`, and Balanced Halves appears in the colleges-needed figure without any introduction in the main body. Either regenerate the deck PNGs without the variants or add one script clause ("two variants of these designs also appear and behave like their parents"). The figures' "Both-arms spillover"/"Control-only spillover" differ from the deck's "Spillover to both arms/controls only", and "> 58 in 1/4 yr, else 44+" will puzzle this audience; say "grey means 58 colleges were not enough in at least one year".

**m13. High Incidence Focus is blue under both arms on the NC heat map (NC accuracy).** Its ratios are 0.92–0.95 under both arms, better than every randomized design, and the script does not mention it. Someone will ask "so why not just treat the worst-hit colleges?". The script answers it under controls only ("not randomized ... if high-burden areas also differ ... error rose roughly fourfold"), which also covers this case. Move that answer up a line so it is clearly about both regimes, and make "In our check where they did" concrete: "When we let incidence affect the score as well, its error rose roughly fourfold."

**m14. Take-home card mixes grid and NC numbers (Take-home).** "Lowest or near-lowest error with all colleges taking part (21–37% below SRS)". The 37% is the 100-cluster grid and 21–30% is NC. Write "(21–30% below SRS in NC, 37% on the grid)", matching the both-arms card's "within 8% on the grid and 3% in NC".

---

## Pilot and SUDDEN work: accuracy check

- Lenoir CC serves Lenoir, Greene and Jones counties (state list, `NC_Colleges_and_Service_Areas.csv`). Correct.
- Workforce development and disaster response; outreach beyond transport to hospital; students earn credit; flooding; pilot first, then share with other counties; colleges talk to each other. All match the 2026-04-27 notes (subject to m2–m4).
- "SUDDEN findings show sudden death tracks access to care, income and rurality". Consistent with the notes ("social problem ... access/money"; tract rates vary by income and urban/rural). "Tracks" is appropriately softer than "driven by".
- Burden slide figures attributed to Mirzaei et al. and Gan et al.; corrected-count rates labelled. Correct.
- The algorithm is credited to Gan et al. and Nanavati et al. and to SUDDEN adjudication, and Ashkan is credited for the counts. Correct.
- IRB 14-2036 on the thank-you slide is author-supplied (DECISIONS 2026-10-28), not invented.
- The limitation "if the analysis ignores spillover, every design underestimates the effect" is supported (CTJ Results, bias −0.12 to −0.40 for every identified design).

## Framing (checked)

- Education outcome, not deaths: service-area amber card, pipeline step 5, NC incidence bullet, β glossary, budgets fig-note. Clear.
- Design evaluation, not program evidence: title script, last limitation. Clear, apart from M2 on the outline.
- Take-home: "reasonable middle ground with limited downside" matches the DECISIONS framing. Controls-only power shortfall and the small NC targeting gain are both stated.
- Style: the script has no colons; slide colons are label colons only ("**Sudden unexpected death (SUD):**", "**Both arms:**"). No sentence starts with "Because". Canonical design names used on all slides; shorthand only in the script (m11).
