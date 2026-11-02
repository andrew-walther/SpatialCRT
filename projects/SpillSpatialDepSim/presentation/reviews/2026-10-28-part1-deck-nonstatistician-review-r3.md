# Part I deck review (non-statistician, criminal-justice/behavioral-health reader), round 3

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.revealjs.md (Draft 1 with round-1 and round-2 revisions) and the figures it shows in presentation/figures/

**Verdict:** Blocking issues: 1 (3 major, 9 minor)

The deck now reads well for someone without statistics training. Probation is named in the first minute. MSE gets a plain definition and an anchor ("off by about 0.41 on a true effect of one"). The "about a third too low" sentence is spoken on both the collinearity and the Recommendation slides, and the scripts have no colons and no sentences that start with "Because". What still trips me is the last mile. The take-home asks an agency to look at which controls its allocation exposes, but no slide says what to do with the answer. The Recommendation slide then suggests tightening the allocation in the direction that makes the problem worse.

---

## Round-2 BLOCKING and MAJOR items

| Item | Status | What changed / what remains |
|---|---|---|
| B1. Recommendation map contradicts the script and seems to fail the deck's own check | Map settled by author decision; wording partly fixed | The wrong "east is well separated" line is gone, and the map is now framed as a suggested best-effort allocation, which is the author's call. The "confirm that some control districts border no intervention district" instruction was dropped, but the take-home still says "look at which controls the chosen allocation exposes". A new line ("further restrictions could guarantee that no intervention district borders another") points away from the fix. See B1 below. |
| M1. Checkerboard bias in its best setting never said aloud | Fixed | Collinearity script: "even in its best setting the checkerboard underestimates the effect by about a third". Recommendation script: "its value is avoiding the bad draws, not being unbiased". |
| M2. "Controls only" rationale fits Case 1, not Case 2 | Settled by author decision | The Case 2 provider channel into controls is now on the slide and in the script, and the "reinforce each other → SRS" condition is stated. One wording point remains (Minor 3). |
| M3. Leadership channel doesn't match the randomized unit; division stratification opens a channel | Partly fixed | "Training for officers and chiefs" and "chief probation/parole officers" were added. The unit is still "county or district", and the cross-arm leader at district level (the division director) isn't named. See M2 below. |
| M4. Main talk runs over 20 minutes | Partly fixed | The thank-you is now "standalone only". With the title, outline and thank-you dropped in the prelim, the timed slides sum to 20.0 min (about 2,475 words, roughly 124 wpm). That's on target but leaves no slack. See M3 below. |

---

## BLOCKING

### B1. The take-home asks for a check without saying what to do with the result, and the Recommendation slide points the other way
- **Slides:** "For both studies we recommend a best-effort block stratified allocation" and the Take-home "State the mechanism" card
- **Text:** Recommendation bullet: "further restrictions could guarantee that no intervention district borders another". Take-home card: "look at which controls the chosen allocation exposes before the trial starts". Script: "look on the map at which controls their allocation exposes."
- **Problem:** The collinearity slide taught me that trouble comes when *every control* touches a treated area, and that a tighter checkerboard makes that more likely. The Recommendation slide then suggests restrictions that would make the NC map *more* like a pure checkerboard. As a probation reader I hear "separate the red districts even more", which is the opposite of what the collinearity slide implied. Then the take-home tells me to look at the exposed controls and stops. If I look at the example map, nearly every blue district appears to touch a red one (checked by eye, not with polygons). I can't tell whether that's fine, since the decision table says "BSS" for the realistic column regardless, or a problem that needs fixing. An agency director would leave not knowing what to do after looking.
- **Fix (wording only, map unchanged):**
  1. Replace the "further restrictions" bullet and script sentence with something like: "This is a suggested, best-effort allocation for an irregular map. Separating treated districts further would also leave fewer controls with no treated neighbor, so the trade-off has to be judged map by map."
  2. Give the take-home card a consequence: "…and look at which controls the allocation exposes. If every control touches an intervention area, expect the estimate to understate the effect, and report it that way. If the map allows it, keep a few controls with no treated neighbor." Say one sentence of this in the take-home script.
  3. One spoken line on the Recommendation slide to link the map to the "a third too low" sentence: "On this map most control districts border an intervention district, so the same caution applies here." The author should confirm this with the adjacency list before saying it.

---

## MAJOR

### M1. The larger-grids figure labels two different panels "ρ = 0.0" and "ρ = 0"
- **Slide:** "The same pattern holds on the 3×3 and 3×4 grids" (fig_mse_larger_grids.png)
- **Text (facet labels):** "3×3 grid · ρ = 0.0", "3×3 grid · ρ = 0", "3×4 grid · ρ = 0.0", "3×4 grid · ρ = 0"
- **Problem:** The 0.01 label is cut off, or rounded to one decimal place. To me those read as the same setting twice. With different lines in each pair, I can't work out which panel is the one with spatial dependence, so I can't check the "16 of 16 / 12 of 16" claim against the picture.
- **Fix:** In presentation/tools/make_figures.R, use the same facet labels as the 2×4 figures: "Spatial dependence (ρ = 0.01)" and "No spatial dependence (ρ = 0)", or at least format ρ to two decimal places.

### M2. Which unit is randomized, and which leader actually crosses arms, is still unclear
- **Slide:** "North Carolina probation trials assign judicial districts…"
- **Text:** title "…assign judicial districts"; body "randomized by county or district"; "chief probation/parole officers and judicial district managers oversee districts and often span several units and counties"; script "A manager who likes the training can carry it into control areas they oversee"
- **Problem:** The title says districts, the bullet says county or district, Case 2 randomizes counties, and the map shows districts. If districts are randomized, which is what the recommendation shows, a judicial district manager oversees a single cluster and can't carry training into a control *district*. The leader who spans both arms at that level is the Judicial Division Director. Stratifying by division guarantees every director oversees both arms. Also, "CPPOs … oversee districts" isn't quite right. Chiefs lead units of officers, and JDMs lead districts (manuscript, Background of Probation). Anyone who works with NC DCS will notice.
- **Fix:** Change the title to "…assign counties or judicial districts". Change the chiefs bullet to "Chiefs lead units of officers, district managers lead districts, and division directors lead the four divisions." Change the channel bullet to "If counties are randomized, a district manager can carry the training into control counties. If districts are randomized, division directors and neighboring offices are the channel." Mirror this in the script.

### M3. No slack at 20.0 minutes, and the densest slides are the likeliest to run long
- **Slides:** Where BSS falls (223 words, 1.75 min), Why: collinearity (175 words, 1.25 min), Recommendation (158 words, 1 min), Limitations [running 20.5]
- **Problem:** The prelim version lands at exactly 20.0 min. The Recommendation slide runs about 158 words a minute, the fastest pace in the deck, on the slide where the audience also has to study a map. If B1's sentences are added, it gets longer still.
- **Fix:** Give the Recommendation slide 1.25 min. Take the time from Limitations by speaking only its first bullet (collinearity) and the citation, about 0.25 min. Or cut "On this grid, without spatial dependence or with spillover into both arms…" from the dot-plot script, since the next two slides cover it.

---

## MINOR

1. **Slide title colon.** Recommendation: "Why controls only: with a consistent treatment effect…". Fix: "Controls only, because with a consistent treatment effect…"
2. **"Consistent treatment effect" sounds like jargon.** Recommendation slide and script. To me it could mean a steady dose, or a statistical property. Fix: "the program works the same wherever it is delivered, so a treated district gains nothing extra from treated neighbors".
3. **Case 2 bullet mixes counties and districts.** Text: "In Case 2, a provider … carrying it into control counties as well. If treated districts were instead expected to reinforce each other, SRS would be the better choice". Fix: "If treated areas were instead…". Also, the map note says "Case 2 would do the same with counties", but bullet 1 says "keep same-arm districts apart". Fix: "keep same-arm units apart".
4. **Square-root axis not mentioned.** Where BSS falls. The axis title says "(square-root scale)", but neither the note nor the script does, and the uneven tick spacing (0.50, 1.00, 1.50 bunching right) looks odd to a lay reader. Fix: add to the note "The scale is stretched near zero so small errors can be told apart."
5. **Collinearity slide is technical on the slide face.** Text: "$[\mathbf{1}, x, z]$ loses rank… `lagsarlm` drops $z$". The model slide is allowed to be technical, but this is the slide a probation audience most needs. The script is clear, but the slide isn't. Fix: lead the left column with "If every control is exposed, exposure is just 'not treated', so the effects can't be separated", and move the rank and `lagsarlm` wording to a smaller second line.
6. **Moran's I and Geary's C named but not used.** Spatial dependence slide. These are unexplained terms that never come back. Fix: drop them from the slide, or say "standard measures of spatial clustering".
7. **Negative numbers read aloud.** Collinearity script: "about -0.35 at psi one-half". Fix: write "about minus 0.35" so the speaker doesn't stumble.
8. **Design-gap slide is dense for 0.75 min.** "multiple-membership models", "restricted randomization", "saturation designs" all appear on one slide, and the bottom card is small type. The script handles it. Fix: no change needed if time is tight. Otherwise drop the "multiple-membership" clause from the slide.
9. **Long spoken sentences (all 33–38 words).** None are hard, but these two will run out of breath:
   - Collinearity: "On the two-by-four grid that is a bias of minus psi without spatial dependence, and with it about -0.35 at psi one-half, so even in its best setting the checkerboard underestimates the effect by about a third." Split after "psi one-half".
   - Recommendation: "The recommendation rests on the proposed studies' assumption that, with a consistent treatment effect, treated districts gain nothing more from other treated districts, while nearby control districts can pick the program up." Split it: "The recommendation rests on what the proposed studies assume. A treated district gains nothing extra from treated neighbors, but nearby control districts can pick the program up."

**Pacing summary:** The pace is about 124 words a minute overall, which is comfortable. The per-slide timings are realistic except the Recommendation slide (M3). Motivation and the dot plot now work cold for a non-statistician. The remaining risk to the take-home is B1. An agency needs one sentence on what to do after it looks at its map.
