# SUDDEN talk (2026-10-29) clinician/epidemiologist review, round 1

**Date:** 2026-10-28
**Deck:** `presentation/sudden_talk_2026-10-29.qmd` (Draft 1), figures in `presentation/figures/`
**Reviewer stance:** a clinician or epidemiologist on the SUDDEN team, comfortable with trial design but not spatial statistics. Read-only review; the deck was not edited.

**Blocking issues: 1**

Overall, the story is sound and well suited to this room. The planning slide states plainly that the trial measures education and not deaths, the spillover schematic is clear (and correct, as I checked cell by cell), the budgets/reach slide and the rates-versus-deaths point are useful for this group, and the limitations slide says that this is a design evaluation and not evidence the program works. The main problems are length (the script runs about 4 to 6 minutes over), one methods slide that will lose the room, and a few places where planning numbers ("36 colleges", "one SD") could be taken out of context. The closing questions are close to right but leave out the endpoint and the effect size that would matter.

---

## BLOCKING

### B1. The script is about 25 minutes, not 20 (whole main body)

Word counts of the `.notes` blocks (title to "Thank you!", time headers excluded) total **3,445 words**. At 130–150 words per minute that is **23.0–26.5 minutes**, against the 20.5 minutes in the time headers. Several slides are budgeted well above a natural speaking rate.

| Slide | Words | Budget | Implied rate (wpm) |
|---|---:|---:|---:|
| Planning question (NC map) | 186 | 1.0 | 186 |
| Outcome model (equation) | 263 | 1.5 | 175 |
| Grid results | 208 | 1.5 | 139 |
| Partial participation | 183 | 1.0 | 183 |
| One-SD effect in practice | 169 | 0.75 | 225 |
| Targeting burden | 109 | 0.5 | 218 |
| Limitations and next steps | 161 | 0.75 | 215 |
| Spillover examples | 145 | 0.75 | 193 |
| Design gap | 140 | 0.5 | 280 |

A 4–6 minute overrun eats the discussion period, which this talk exists for, or forces rushing through the application slides, which the audience cares about most.

**Fix (target 2,700–2,850 words).**
- Move "One outcome model links treatment, spillover and neighbors" to backup and replace it with two spoken sentences on the spillover slide (see M1). Saves about 230 words.
- Merge "The grid study varied every condition..." into the grid results slide as one spoken sentence ("We tried every combination of incidence pattern, spillover strength and spatial dependence, 250 simulated trials each"). Saves about 70 words and a slide.
- Cut the design-gap script to about 70 words. The literature point ("few studies compare them under spillover and uneven incidence together") is for statisticians.
- Trim the partial-participation, effect-size, targeting and limitations scripts to about 140, 120, 80 and 110 words respectively (suggested cuts under M4, M5 and the minor items).
- Then re-set the [min | running] headers from the actual word counts at 140 wpm.

---

## MAJOR

### M1. The equation slide will lose a clinical audience ("One outcome model links treatment, spillover and neighbors")

Text: the matrix equation $Y = (I - \rho W)^{-1}[\tau Z + S(Z) + \beta X + \epsilon]$, two spillover forms with an elementwise product, and an 11-row symbol glossary, for 1.5 minutes.

Problem. This is the densest slide in the deck and the only one that needs matrix algebra. "For those who like equations" signals that most of the room can tune out, and they will for 90 seconds. The two ideas a clinician needs here are already on the spillover slide (treatment effect, spillover from treated neighbors) plus one sentence: "outcomes of neighboring colleges also move together" and "SUD incidence guides allocation but is assumed not to change the education score".

Fix. Move the slide to backup ("Backup · The outcome model"). On the spillover slide add one bullet and say it aloud:

> - The analysis estimates the program's direct effect while allowing for spillover and for neighboring colleges' scores moving together

Keep the β point (education outcome, β = 0) where it already lives on the planning slide and in the β = 1 backup.

### M2. "Colleges needed" and "29 colleges / 58 colleges" are ambiguous (Colleges needed; Partial participation; One-SD effect slides)

Text: "About 36 colleges detect a one-SD effect if spillover reaches every college"; table headers "29 colleges", "58 colleges"; "Colleges joining 58 / 40 / 29 / 20".

Problem. Earlier slides use "29 of 58" to mean colleges *treated*. A clinician will read "36 colleges" as 36 colleges *receiving the program*, when it means 36 *participating* colleges, about half treated and half control, with the rest of the state not in the trial. On the effect-size slide, "29 colleges" in the power table is easy to read as "29 treated out of 58".

Fix. Say it once, explicitly, on the colleges-needed slide and keep the wording everywhere after:

> - "Colleges needed" counts participating colleges, about half of them treated and half control

Rename table headers to "29 participating colleges" / "58 participating colleges", and the partial-participation column to "Participating colleges".

### M3. "About 36 colleges" and "one SD" will travel without their assumptions (Colleges needed; Take-home)

Text: take-home card "about **36 colleges** detect a one-SD college effect"; title "About 36 colleges detect a one-SD effect if spillover reaches every college".

Problem. This is the number the group will remember and repeat. It holds only under stated assumptions (ρ = 0.5, spillover equal to 80% of the direct effect, colleges volunteering at random, a one-SD *college-level* effect). "One SD" also sounds large to a clinician (Cohen's d = 1), and the translation to 0.3 SD per professional comes two slides later.

Fix.
- Take-home card: "about **36 participating colleges** detect a one-SD college-level effect (about 0.3 SD per professional with 25 trained per college), under our assumptions".
- On the colleges-needed slide say, before the numbers: "One college-level SD is about a third of a standard deviation for an individual professional; I'll unpack that in two slides."
- In the fig-note, state the spillover assumption in words: "spillover from a fully treated neighborhood equal to 80% of the direct effect".

### M4. The partial-participation slide is dense and its key result is counterintuitive ("If fewer colleges join, the saturation advantage shrinks")

Text: a four-panel, four-line power figure; table; bullets "Under controls only, power rises slowly past about 40 colleges" and "Non-participating neighbors already create the contrast that saturation otherwise supplies"; 183 words in 1.0 minute.

Problems.
1. In the controls-only panel, SRS and Checkerboard power *falls* as more colleges join (SRS from about 0.54 at 45 to 0.52 at 58 at 1 SD; Checkerboard from 0.35 to 0.17). A clinician will notice "more colleges, less power" and the script does not explain it except obliquely.
2. "Non-participating neighbors already create the contrast" assumes knowledge of identification.
3. Two different numbers for the same claim: this slide says saturation is "22% lower" with all 58, the NC accuracy slide says "21–30% lower". They are different settings (ρ = 0.5 with proportional spillover, versus averages over six ρ/γ settings), but the audience will not know that.

Fix. Lead with the planning message and drop the mechanism unless asked:

> Script: "Not every college will join. Under spillover to both arms, every added college helps. Under spillover to controls only, the gain stalls past about 40 colleges, because when almost everyone is in, almost every control college sits next to treated ones and the spillover blurs into the treatment. Saturation's edge is largest when most colleges join, about 20 percent less error with all 58, and almost gone with 20. So how many colleges will join is a real design question."

Add to the fig-note "ρ = 0.5, spillover 80% of the effect; the 21–30% on the earlier slide averages six settings". Consider showing only the top row (1 SD) of the figure.

### M5. The NC accuracy slide puts the bullets in the opposite order to the figure ("With all 58 colleges, saturation beats SRS when spillover reaches only controls")

Text: left bullet column starts "**Controls only:**", right column "**Both arms:**". In `nc_annual_accuracy_six.png` the left half is "Spillover reaches both arms" and the right half is "Spillover reaches controls only".

Problem. The eye reads left bullets with the left heat map, and that pairs the controls-only claim ("21–30% lower") with the both-arms panel, where saturation sits at 0.97–1.01. The script itself has the halves right.

Fix. Swap the two bullet columns so the left column is "Both arms" and the right is "Controls only". Also check the same order elsewhere (colleges-needed and power figures put both arms on the left, so keep that convention throughout).

### M6. The deck never answers "why not just give it to the highest-need areas first?" in the main body

Text: High Incidence Focus is described as "isn't really randomized" (NC allocations script); take-home bullet "fixed targeting (High Incidence Focus) varied by year and assumes incidence is unrelated to the outcome".

Problem. This is the first question a clinical team will ask, and it has an equity side. The real answer (if high-burden areas also differ in baseline education or staffing, treatment and burden become confounded, and High Incidence Focus error rose 1.4–4.4-fold in the β = 1 check) is only in backup. "Assumes incidence is unrelated to the outcome" is too compressed to land.

Fix. Add one sentence to the NC accuracy or take-home script:

> "Treating the highest-rate colleges is the intuitive choice, but it isn't randomized. If high-burden areas also differ in their staff or training to begin with, we can't separate that from the program, and in our check where they did, its error rose up to fourfold. Saturation gets more of the program to higher-burden regions while keeping randomization in every region."

### M7. "Which colleges receive the program first" implies a rollout, but the deck never discusses it (planning slide; questions)

Text: title "The planning question is which colleges receive the program first".

Problem. "First" implies control colleges get the program later. Clinicians will ask whether this is a delayed-start (waitlist or stepped-wedge) trial, whether control colleges can be withheld at all, and when the education score is measured. The deck models a parallel two-arm trial at one time point and does not say so.

Fix. Add to the planning-slide script: "We model a simple two-arm comparison at one time point, as if control colleges get the program after the outcome is measured. A staggered rollout is possible and would change the design." Add it as a sub-question to question 2 (see M8), or list it in Next steps.

### M8. The closing questions omit the endpoint and the effect size that would matter ("What we would like from you")

Text: four questions on spread, participation, professionals per college and score variation, and rates versus deaths.

Problem. The four questions are good and correctly prioritized (spillover regime first). But the two inputs that most change sample size and that only this group can supply are missing: what the education outcome is (knowledge test? when measured?), and what size of improvement would be worth detecting. Next steps lists "Define the education endpoint" but the questions do not ask it. The rollout and willingness-to-be-randomized question (M7) is also missing.

Fix. Keep four questions by merging and adding:

1. **How does the program spread?** Mainly to untreated colleges nearby, or to every college, treated ones included? Is a shared border the right idea of "neighbor"?
2. **What is the education outcome**, when would it be measured, and what improvement would be worth detecting?
3. **How many colleges and professionals?** How many colleges would join (and would higher-incidence ones be keener), how many professionals each would train, and would control colleges accept a delayed start?
4. **Target rates or deaths?** Should allocation favor high-rate areas or the areas where most deaths occur?

Question 3 on score variation (ICC) is hard for clinicians to answer directly. Ask instead whether single-county pilot scores exist, from which the speaker can estimate it.

### M9. Small-area rate instability is not mentioned, and this audience will raise it ("SUD incidence is uneven across North Carolina...")

Text: script "the highest rates tend to be in service areas with small populations" and "a college that is high in one year may not be high the next".

Problem. Epidemiologists will immediately read the year-to-year shifts in small service areas as partly random noise in small counts, not real change. That reinforces the case against High Incidence Focus, but the deck presents it as "shifts from year to year" without saying why. Someone may also ask why the rates were not smoothed or pooled.

Fix. Add to the script: "Part of that year-to-year movement is just small numbers. A few deaths swing the rate in a small service area. That is one reason we treated each year as its own planning scenario rather than trusting one year's ranking." Be ready to say that pooled or smoothed rates are an option for the real plan.

---

## MINOR

1. **Design gap slide, title and card.** "guidance on which clusters to treat is thin" and "block stratified allocation" are statistician phrasing. Suggest the title "We know how to analyze spillover, but not how best to choose who is treated" and the card wording "compared random and stratified allocation on small grids".
2. **Title.** "Sampling Design for Spatial Cluster Randomized Trials Under Heterogeneous Incidence" is opaque to this audience. Keep it, but say in the opening script "In plain terms, which colleges should get the program first when death rates differ across the state". (The script already nearly says this.)
3. **Performance measures slide.** "We use a $t$-based interval with a small-sample correction", "spatial lag model fitted by maximum likelihood" and the $m(a)$ notation add nothing for clinicians. Replace the allocation-risk card text with "A trial gets **one** allocation. We also ask how bad the worst 10% of allocations are." Drop "$t$-based ..." from the coverage card (it is in backup). Replace the opening sentence with "We estimate the program's effect with a spatial model that allows for spillover, and repeat the trial thousands of times in simulation."
4. **Planning slide, who is trained.** Slide says "trained healthcare professionals"; script says "trained first responders and clinicians". Use one term, and confirm with Dr. Simpson which is right, since this audience knows the program.
5. **Grid results card.** "Checkerboard and 2×2 Blocking are fine in one regime and poor in the other" conflicts with the table (1.47 times SRS under both arms is not "fine") and with the take-home ("did not improve on SRS in either regime"). Suggested: "Checkerboard and 2×2 Blocking were worse than SRS in both regimes, and far worse under controls only."
6. **Take-home, "Strict alternation (Checkerboard, 2×2 Blocking)".** 2×2 Blocking is not alternation, and on the NC map it became spatial blocking at about 1.0 times SRS under both arms. Suggested: "Checkerboard never improved on SRS and was about three times worse on the NC map; 2×2 Blocking was no better than SRS."
7. **NC allocations slide.** Slide says Isolation Buffer treats "about 16"; the example map is labeled "14 of 58 treated". Say "about 16 on average (14 in this example)" in the fig-note.
8. **NC incidence slide, panel C.** The four yearly maps are small. Point to one or two specific colleges that change rank (for example, the bright northeastern area in 2020) rather than asking the audience to compare four small maps.
9. **Allocation-risk slide.** "21–43% lower error than SRS's worst 10% in NC" cannot be read off the backup table, whose NC controls-only worst-10% entry for Incidence-Guided Saturation is 0.68 and whose note says entries are divided by the "SRS average". Check that the slide range and the backup table use the same denominator (SRS worst decile versus SRS mean) and say which in the fig-note, or someone who opens the backup will see an apparent contradiction.
10. **Allocation-risk figure.** The x-axis label "MSE relative to simple random sampling (log scale)" and the box-and-whisker encoding take time. In the script, name one comparison only ("look at the top row in the right-hand panels, saturation's whole bar sits left of SRS's average").
11. **Targeting slide.** The column "Share of summed college SUD rates treated" is hard to parse. Rename it "Treated share of total incidence (rate-weighted)" and say in the script "how much of the high-rate burden is in treated areas". The take-home also claims saturation gives "the largest treated share in the highest-burden region", which is true by rule but in NC the regions barely differ (mean ranks 24–34). Soften to "puts the most treatment in the region with the highest rates".
12. **Effect-size slide.** The formula $\sigma_c = \sqrt{\text{ICC} + (1-\text{ICC})/m}$ can stay for epidemiologists, but the script should not walk through it. The fig-note sentence on controls-only power is dense; move it to the script or backup.
13. **Spillover-slide legend.** "in proportion to the share of its neighbors that are treated", while every circle is the same size. Add "(circle size not to scale)" or drop the proportionality phrase from the legend (the script already says it).
14. **Design evaluation, not program evidence.** Stated well on the limitations slide, but late. Add one clause to the opening script: "This is about designing the trial, not evidence yet on whether the program works."
15. **Limitations slide, first bullet.** "The analysis model knew the form of spillover; leaving it out biased every design" is a statistician's point. Suggested: "Our analysis assumed we know how spillover works; if the analysis ignores spillover, every design underestimates the effect."
16. **Design-name consistency.** The deck uses the canonical names in tables, but the scripts and bullets use the shorthand "saturation" ("saturation was 8% higher", "saturation designs"), "spatial blocking" and "random sampling" for SRS. Spoken shorthand is fine; on slides use "Incidence-Guided Saturation" or "the two saturation designs" explicitly. `fig_colleges_needed.png` prints "2x2 Blocking" (letter x) where the deck uses "2×2 Blocking", and the figures say "Control-only spillover"/"Both-arms spillover" while slide text says "controls only"/"both arms".
17. **Colons in prose.** Slide text with avoidable colons:
    - Spillover slide: "**Controls only** is classic contamination: untreated colleges learn from treated neighbors" → "**Controls only** is classic contamination, in which untreated colleges learn from treated neighbors". "**Both arms**: treated colleges also gain..." → "**Both arms** means treated colleges also gain...".
    - Effect-size slide: "Two numbers from this group would sharpen these estimates: **how many**..." → "Two numbers from this group would sharpen these estimates, **how many professionals each college would train** and **how much scores vary between colleges**".
    - Take-home card: "**Saturation** (Incidence-Guided or random-order): lowest..." → "**Saturation designs** gave the lowest or near-lowest error, ...".
    - Limitations script: "for us to do together: define the endpoint" → "for us to do together. We need to define the endpoint...".
    - Glossary labels ("Queen contiguity: ...", "Spatial dependence: ...") are acceptable as definitions but go away if M1 is adopted. Bullet labels such as "**Controls only:**" and "**Both arms, 1 SD:**" are acceptable label style.
    - No sentence starts with "Because" (checked).

---

## What works well (keep)

- The bed net and deworming examples are familiar to epidemiologists and make spillover concrete in under a minute.
- The amber card "The trial measures education, not SUD deaths. Incidence guides allocation only." is exactly the framing this audience needs.
- The budgets-and-reach slide (High Incidence Focus treats half the colleges but reaches under a quarter of working-age adults) and the rates-versus-deaths contrast (62% of rate share, 35% of deaths) are the most useful planning facts in the deck. Expect discussion here.
- The effect-size translation to individual professionals and the explicit request for training numbers and score variation turn abstract power into something the group can act on.
- The backup set anticipates the likely questions (β = 1, omitted spillover, Bayesian models, all-or-nothing spillover).
