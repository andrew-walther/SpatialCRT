# SUDDEN talk (2026-10-29) clinician/epidemiologist review, round 3

**Date:** 2026-10-28
**Deck:** `presentation/sudden_talk_2026-10-29.qmd` (restructured after round 2), figures in `presentation/figures/`
**Background checked against:** `application/notes/SUD_Meeting_Notes.md`, `application/AGENTS.md`, `application/README.md`, round 1 and round 2 reviews
**Reviewer stance:** a clinician or epidemiologist on the SUDDEN team (including the people who built and validated the death-certificate algorithm), new to this deck. Read-only review; the deck was not edited.

**Blocking issues: 0**

(4 major, 11 minor)

The new background works. The burden slide, the pilot slide and the service-area slide give this room the reason for the project, where it came from and what a cluster is, and the education-outcome and design-not-efficacy framing is now stated three times in the right places. The round-2 MAJORs are fixed (all-58 underpowered under contamination is now in the take-home; the outcome is "professionals assessed in each college", with the control-college question on the backup slide). What remains is about how the deck describes the SUDDEN team's own work. The pilot slide puts words in Dr. Simpson's mouth about who is trained and what the process measures are, and the case-identification flowchart calls accidents, homicides and suicides "expected deaths", in front of the people who wrote that definition. One application message also overstates how much incidence guidance targets burden on the NC map.

---

## Round-2 clinician findings: status

| Finding | Status | Evidence |
|---|---|---|
| M1 Take-home omits that 58 colleges are underpowered under contamination | **Fixed** | Take-home controls-only card and script: "Even with all 58 colleges, power for a one-SD effect is only about 50–65%" |
| M2 "Trained professionals" cannot be the control-college outcome | **Mostly fixed** | Service-area script: "education score of the professionals assessed in each college"; effect-size script "averages the professionals assessed there"; backup question 5 asks Dr. Simpson. The effect-size table header still reads "Professionals per college" with "how many professionals each college trains" (acceptable) |
| m1 Outcome model digestibility | Partly fixed (author kept slide) | Script opens "You don't need the algebra"; τ row still the densest line |
| m2, m3 Grid-study script | Fixed | "an effect of one standard deviation"; contradictory sentence replaced |
| m4 Performance-measures wording | Not changed | Still "including the true spillover term" and $m(a)$ (acceptable for this room) |
| m5 Variants not introduced | Not fixed | See m9 |
| m6 Figure labels | Not fixed | See m9 |
| m8 SRS colour on allocation figure | Fixed | New `deck_allocation_distribution.png` shows SRS dark grey; script says so |
| m9 "up to fourfold" | Fixed | "roughly fourfold" |
| m11 Backup colon | Fixed | Rook backup now uses commas |

---

## Timing

Main-body `.notes` words, title through "Thank you!", time headers and bullet markers excluded: **2,892**.

| Pace | Minutes |
|---|---:|
| 150 wpm | 19.3 |
| 140 wpm | 20.7 |
| 130 wpm | 22.2 |

The time headers sum to 21.0 minutes. The talk fits the ±2-minute window at a normal pace but runs just over it at 130 wpm. Several slides are budgeted faster than a read-aloud pace allows, so the headers understate the real time:

| Slide | Words | Budget (min) | wpm |
|---|---:|---:|---:|
| Budgets and reach | 83 | 0.5 | 166 |
| Partial participation | 121 | 0.75 | 161 |
| Effect size in practice | 120 | 0.75 | 160 |
| Case identification | 117 | 0.75 | 156 |
| From pilot to statewide | 151 | 1.0 | 151 |
| Take-home | 149 | 1.0 | 149 |
| Grid results | 148 | 1.0 | 148 |
| Outcome model | 175 | 1.25 | 140 |

See m1 for about 150 words of cuts that bring the script to about 19.6 minutes at 140 wpm.

---

## BLOCKING

None.

---

## MAJOR

### M1. The pilot slide says the trainees are EMS and healthcare professionals; the meeting notes say college students and staff (From pilot to statewide)

**Text.** Flow box "**Who is trained** EMS and other healthcare professionals in the college's service area, so they recognize, report and respond to people at risk". Script: "The training would go to EMS and other healthcare professionals in the service area, so they can better recognize, report and respond to people at risk."

**Problem.** The author's notes from Dr. Simpson (2026-04-27) describe the pilot as working through workforce development and disaster-response teaching, with the goal to "empower CC students/staff to get people to the hospital and share treatments". The notes add that students earn credit for participating and that colleges recruit from underserved, high-risk areas. Nothing names EMS or practising healthcare professionals as the trainees. The "healthcare professionals" population comes from our own trial model (Chapter 3 Application), not from the pilot. This slide is presented as a description of Dr. Simpson's pilot, to his team, so it is the most likely thing in the talk to be corrected from the floor. It matters downstream too. Who is trained decides the endpoint, the professionals-per-college numbers and the ICC on the effect-size slide.

**Fix.**
- Flow box: "**Who is trained** Students and staff in the college's workforce-development and disaster-response programs (in our model, healthcare and EMS trainees), who learn to recognize people at risk and get them to care".
- Script: "The training would reach students and staff in those programs, so they can recognize people at risk and get them to care. For the trial model we assumed the trainees are healthcare and EMS professionals, and that is something I'd like you to correct if it's wrong."
- If the author has newer information from Dr. Simpson that the trainees are EMS and healthcare professionals, keep the wording but cite it in the notes.

### M2. The pilot slide attributes the education score to the pilot and merges it with patient-level process measures (From pilot to statewide)

**Text.** Bullet "Its effect is first visible in **process measures** (an education score, timely care), not in SUD rates". Script: "...with a college's service area as the unit and an education score, not death rates, as the first measure of effect."

**Problem.** Dr. Simpson's process measures in the notes are "medication adherence, timely clinic visits when needed", which are measured in patients or the community. The education score is a trainee-level endpoint that we chose for the design study. The bullet puts both in one parenthesis as if the pilot defined them, and the script says the education score is the pilot's first measure. A clinician will notice that a trainee test score and community medication adherence are different endpoints with different units of analysis, sample sizes and ICCs. Since the endpoint is the first question the talk wants the group to answer, the deck should not pre-answer it in Dr. Simpson's name.

**Fix.**
- Bullet: "Its effect would show first in **process measures** such as medication adherence and timely clinic visits, not in SUD rates. For the design study we assume a trainee **education score** as the endpoint".
- Script: "Dr. Simpson expects the first changes in process measures like medication adherence and timely clinic visits, not death rates. For the design work we assumed an education score of the people trained, which is one of the things I'd like your view on."

### M3. The case-identification flowchart calls non-natural deaths "expected", and its steps are not documented in the project (Case identification)

**Text.** Exclusion box "Exclude expected deaths (accident, homicide, suicide, overdose, COVID-19)", followed by the main box "Natural, unexpected deaths" and then "Exclude end-stage chronic disease (cancer, cirrhosis, renal failure)". Script: "...expected deaths such as accidents, suicides, overdoses and COVID-19, and deaths from end-stage chronic disease such as cancer."

**Problem.**
1. Wrong category. Accidents, homicides, suicides and overdoses are non-natural deaths, not expected ones. The expected deaths in this case definition are the end-stage chronic-disease deaths in the next box. As drawn, the box labelled "Natural, unexpected deaths" follows only the removal of non-natural deaths, so it is neither yet. The SUDDEN team wrote this definition and will catch the swap.
2. Unconfirmed content. The 412,514 and 111,665 totals come from Habib's paper (ROADMAP 2026-09-24 attribution note). The exclusion list and its order are not in the application code, `application/README.md` or the current Chapter 3. The archived 2026-10-02 chapter draft carried a NEEDS-AUTHOR-CONFIRMATION flag on exactly this list, and `application/AGENTS.md` records that the free-text regex in Habib's script "can't match 'KIDNEY FAILURE' or 'LIVER FAILURE'", which bears on the "cirrhosis, renal failure" box. COVID-19 exclusion was presumably Habib's addition, since the Gan et al. algorithm predates COVID. Showing a step-by-step flow in front of its authors invites "that's not the order we used".

**Fix.**
- Relabel: "Exclude non-natural deaths (accident, homicide, suicide, overdose) and COVID-19" → main box "Natural out-of-hospital deaths" → "Exclude expected deaths from end-stage chronic disease (cancer, cirrhosis, renal failure)" → "**23,523** presumed sudden unexpected deaths".
- Script: "...then removes non-natural deaths such as accidents, suicides and overdoses, and COVID-19, and then expected deaths from end-stage chronic disease such as cancer."
- Send the flowchart to Ashkan before the talk to confirm the steps and order. If he has the counts at each step, add them; if not, add a fig-note "Steps summarize the case definition; counts at intermediate steps not shown."

### M4. The application slides say saturation sends more of the program to high-burden regions without saying the NC gain is small (NC accuracy; Take-home)

**Text.** NC accuracy script: "Saturation sends more of the program to higher-burden regions while keeping randomization in every region." Take-home card and script: "...and the most treatment in the highest-rate region".

**Problem.** For this audience, reaching the burden is the main reason to care about incidence guidance. The backup "Targeting burden" slide shows that on the NC map the four regions have similar average incidence, so Incidence-Guided Saturation puts 51.1% of the rate share in treated areas against 50.4% for the same rule in random region order, and 49.2% of SUD deaths. A clinician who hears "the most treatment in the highest-rate region" will reasonably assume the program reaches materially more of the burden, and will feel misled if the backup is shown in discussion. The author's framing rule (AGENTS.md, 2026-10-21) is to say guidance matches Saturation Quadrants' accuracy while targeting high-burden regions when burden differs between regions, with a small gain on the NC map. The main body never says the gain is small.

**Fix.**
- NC accuracy script, after the saturation sentence: "On this map the four regions differ little in average incidence, so the extra targeting is small, under one percentage point."
- Take-home card: "...and the most treatment in the highest-rate region (a small shift on the NC map, where regions are similar)". In the script, add "though on this map that shift is small."
- These add about 25 words; take them from the m1 cuts.

---

## MINOR

**m1. Timing (all main-body slides).** 2,892 words is 20.7 minutes at 140 wpm and 22.2 at 130. Cuts that lose nothing, about 150 words in all:
- Pilot script, drop "Here is how the application came about." and merge the "Neighboring colleges share..." sentence into the statewide sentence.
- Case-identification script, drop the last sentence ("Ashkan tallied them by county of residence and year...") since the slide shows it.
- Pipeline script, drop "For each setting we drew 100 allocations and 100 sets of scores for each." (on the slide).
- Take-home script, see m2.
- Effect-size script, drop "The table shows power for that individual effect with random sampling under spillover to both arms." and start with "With 25 per college...".
- Adjust the time headers on budgets (0.5 → 0.6) and partial participation and effect size (0.75 → 0.85) so the running totals reflect the real pace.

**m2. Long sentences to read aloud.** The take-home sentence "When spillover reaches mainly untreated colleges, it had error 21 to 37 percent below random sampling, a lower risk of an unlucky allocation, and the most treatment in the highest-rate region, though even all 58 colleges give only about 50 to 65 percent power for a one-SD effect." is 48 words, the key message, and too long to land. Split it: "When spillover reaches mainly untreated colleges, it had 21 to 37 percent less error than random sampling and fewer unlucky allocations. But even with all 58 colleges, power for a one-SD effect is only 50 to 65 percent." Eleven other script sentences run 32–36 words (for example, the burden-slide NC sentence and "When almost everyone is in, ..."). Each splits naturally at its ", and" or ", so".

**m3. Pilot flowchart mixes a sequence with attributes (From pilot to statewide).** Arrows imply "this leads to that", but "Who is trained" is a property of the pilot, not a step after it, and "SUD as a social problem" is a finding. Make it four steps (SUDDEN finding → Lenoir pilot idea, with who is trained inside that box → statewide extension → this project). Also, "neighboring colleges share staff, trainings and patients" (bullet and script) is odd, since colleges do not have patients. Say "neighboring service areas share hospitals, EMS agencies and instructors".

**m4. "Statewide extension" reads as Dr. Simpson's plan.** The notes say "pilot study first, then share findings with all other counties" and "Hope is that other schools adopt the initiative". A randomized statewide rollout is our proposal. Script: "The natural next step, and what this project is designed for, is to offer it through all 58 colleges after the pilot..." The effect-size slide's "pilot scores would tell us the second" fits this ordering well.

**m5. Service-area map (Service-area clusters).** The legend sorts as text (1, 10, 11, ..., 2, 20), 58 colours make neighbours such as 3, 29 and 30 hard to tell apart, and Lenoir CC (29), named on the previous slide, is not marked. Highlight Lenoir CC's service area (bold outline or call-out) and say "Lenoir's service area is here" in the script. Sorting the legend numerically or dropping it in favour of the numbered labels would also help.

**m6. Neighbor definition is inconsistent between slides.** Service-area slide: "share a border or a corner". Pipeline step 2: "neighbors from shared borders". Use "shared borders or corners" in step 2.

**m7. Case-identification details.** "heart-failure deaths are kept in the corrected counts" sits inside an exclusion box, which reads as if it were excluded. Move it to the right-hand bullet that already explains 21,147 against 23,523. The intermediate main boxes ("Out-of-hospital deaths", "Natural, unexpected deaths") have no counts, which an epidemiology audience expects in a CONSORT-style chart (see M3).

**m8. Shorthand design names.** "saturation" alone (grid card "saturation was 8% higher", slide titles "saturation beats SRS", "saturation also lowers the risk", several script lines) is ambiguous once Saturation Quadrants appears in the power figures. "fixed targeting swung from year to year" in the take-home script names no design. Use "Incidence-Guided Saturation" in titles and cards, and "High Incidence Focus" in the take-home script. "Spatial blocking" for 2×2 Blocking is explained and acceptable.

**m9. Carry-over from round 2 (m5, m6).** Saturation Quadrants and Balanced Halves appear in `fig_colleges_needed.png` and `fig_power_vs_n.png` but are introduced only in a fig-note. Add to the designs script: "Two variants, Saturation Quadrants with regions in random order, and Balanced Halves, behave like their parents and appear in some later charts." The colleges-needed figure still says "Both-arms spillover"/"Control-only spillover" and "> 58 in 1/4 yr, else 44+". In the script, say "grey means 58 colleges were not enough in at least one year".

**m10. Burden slide rates differ from Ashkan's paper.** The card "84 → 101" uses the corrected counts; Habib (2026) reports 75 → 91. Ashkan will likely be in the room. Add "(corrected counts)" to the card, matching the case-identification slide.

**m11. Optional motivation.** The notes record that Lenoir County floods regularly and SUD rates spike after flooding. One clause on the pilot slide ("disaster-response training, in a county where SUD spikes after floods") explains why disaster response is the vehicle and shows the team their context was heard.

---

## Flowcharts as diagrams

- **Pilot → statewide (5 boxes).** Readable, but see m3 (mixes steps and attributes) and M1 (trainee description).
- **Service-area steps (4 stacked boxes).** Clear. The amber card on the outcome is well placed.
- **Case identification.** Layout is a standard CONSORT form and easy to follow. Content issues in M3 and m7.
- **Application pipeline (6 boxes).** Correct and in the right order, and the "fixing every rule before looking at results" title is a good message for clinicians. Only m6.
- **Spillover schematic.** Checked cell by cell. Every control in the left panel has a treated queen neighbour, and in the right panel every cell except the bottom-right treated cell (no treated neighbour) gains. Correct.

## Framing (checked)

- Education outcome, not deaths: service-area amber card, pipeline step 5, NC incidence bullet, β glossary. Clear (subject to M2 on whose endpoint it is).
- Design evaluation, not program evidence: title script and last limitation. Clear.
- Take-home: the softened "reasonable middle ground with limited downside" is accurate to the results and suits this audience. Only M4 needs adding.
- Style: no colons in the script, label colons only on slides ("**Both arms:**"), acceptable. No sentence starts with "Because". Design names otherwise match the canonical set.
