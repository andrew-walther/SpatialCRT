# Part I deck review (non-statistician, criminal-justice/behavioral-health reader), round 4

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.revealjs.md (Draft 1 with round-1 to round-3 revisions) and the figures it shows in presentation/figures/ (fig_spread_2x4, fig_mse_2x4_controls, fig_mse_2x4_both, fig_collinearity_grid, fig_mse_larger_grids, BSS_Assignments_NC_Map_Applications_Text_Labels)

**Verdict:** Blocking issues: 0 (3 major, 9 minor)

I can follow this talk from start to finish without losing the thread. The round-3 blocker is mostly resolved. The Recommendation slide no longer suggests tightening the map, and it now says a refined allocation could leave some controls with no intervention neighbor. The take-home card also now says what an all-exposed map does to the estimate. What is still missing is one spoken sentence that tells an agency what to do after it looks at its map. That sentence is on the card in vague words ("absorb part of the spillover") and is not in the script. Two other items remain: the Recommendation slide is still the most rushed slide, and the leadership channel doesn't quite fit district randomization. None of these would make me misunderstand the take-home, so I don't count any as blocking.

A correction to my round-3 M2. I wrote that "CPPOs … oversee districts" was wrong. The manuscript (Background of Probation, line 542) says exactly that CPPOs and JDMs oversee districts, so the deck matches its source. I withdraw that part.

---

## Round-3 BLOCKING and MAJOR items

| Item | Status | What changed / what remains |
|---|---|---|
| B1. Take-home asks for a check without saying what to do; Recommendation pointed the wrong way | Mostly fixed (map settled by author decision) | "Further restrictions could guarantee…" is gone. The bullet and script now say "a refined allocation could also leave a few control districts with no intervention neighbor, which keeps the program and its spillover separable". The Recommendation script ties the checkerboard to "underestimates the effect by about a third". The take-home card adds "If every control borders an intervention cluster, expect the estimate to absorb part of the spillover." What remains: the card's wording is vague, and the take-home script still stops at "look on the map". Downgraded to M1 below. |
| M1. Larger-grids figure labels "ρ = 0.0" and "ρ = 0" | Fixed | Facets now read "3×3, ρ = 0.01", "3×3, ρ = 0", "3×4, ρ = 0.01", "3×4, ρ = 0". |
| M2. Randomized unit and cross-arm leader unclear | Partly fixed | The slide bullet now says "Leaders can carry the training…", which avoids naming the wrong role. The title still says "assign judicial districts", the bullet says "by county or district", and the script still says "A manager … can carry it into control areas they oversee". See M2 below. My CPPO point is withdrawn (see above). |
| M3. No slack at 20.0 min; Recommendation fastest slide | Not fixed | The prelim version is still exactly 20.0 min (2,409 words, about 120 wpm). The Recommendation slide is 152 words in 1 min, and Limitations is 77 words in 0.5 min (154 wpm). See M3 below. |

---

## BLOCKING

None.

---

## MAJOR

### M1. The take-home still doesn't say aloud what to do after looking at the map
- **Slide:** Take-home, "State the mechanism" card and script
- **Text:** Card: "look at which controls the chosen allocation exposes. If every control borders an intervention cluster, expect the estimate to absorb part of the spillover." Script: "So investigators should state the spillover they expect, choose the design to match, and look on the map at which controls their allocation exposes."
- **Problem:** The audience hears the instruction to look but not what to conclude. On the card, "absorb part of the spillover" doesn't tell me the direction. As a probation reader I don't know whether the program will look better or worse than it is. The deck already said it plainly twice: "underestimates the effect by about a third". An agency director should leave with the action: keep a few unexposed controls if the map allows, and otherwise expect and report an understated effect.
- **Fix:** Card: "…If every control borders an intervention area, expect the estimate to understate the effect. Where the map allows, keep a few controls with no intervention neighbor." Add one script sentence after the last bullet: "If every control borders a trained area, expect the estimate to understate the training's effect, so keep a few controls clear of it where the map allows, or say so when you report."

### M2. The leadership channel fits county randomization, not district randomization
- **Slide:** "North Carolina probation trials assign judicial districts, a natural set of spatial clusters"
- **Text:** Title "…assign judicial districts"; bullet "randomized by county or district"; script "A manager who likes the training can carry it into control areas they oversee, which is spillover into controls."
- **Problem:** The recommended Case 1 map randomizes districts. According to the manuscript, chiefs and JDMs span units and counties *within* a district, so under district randomization their channel stays inside one cluster. A probation person will ask how a manager reaches a control *district*. The answer (division-level leadership, shared providers and courts, officers who talk across neighboring offices) isn't given. Case 2 randomizes counties, so the title's "judicial districts" fits only half the application.
- **Fix (wording only):** Title: "North Carolina probation trials assign counties or judicial districts, natural spatial clusters". Script: "If counties are randomized, a manager who likes the training can carry it into control counties they oversee. If districts are randomized, the route is through division leadership and neighboring offices that share courts and providers." The manuscript names the Judicial Division Director, so the author can use that title if they want to.

### M3. The Recommendation slide is still the most rushed, and there is easy slack elsewhere
- **Slides:** Recommendation (152 words, 1 min, four dense bullets plus a map); "With spillover into both arms…" (60 words, 0.75 min, about 80 wpm); Limitations (77 words, 0.5 min)
- **Problem:** The audience has to read a 30-district map and four bullets while the speaker covers 152 words at about 150 wpm. This is the slide a probation agency cares about most. Meanwhile, the both-arms slide has a quarter-minute of unused time.
- **Fix:** Move 0.25 min from the both-arms slide (0.75 to 0.5) to the Recommendation slide (1 to 1.25), and update the running totals. The total stays at 20.0. Splitting the 40-word map sentence (Minor 9) also helps.

---

## MINOR

1. **Slide colon (carried).** Recommendation bullet: "Why controls only: with a consistent treatment effect…". Fix: "Controls only, because with a consistent treatment effect…"
2. **"Consistent treatment effect" (carried).** Recommendation bullet and script. To me this could mean a steady dose or a statistical property. Fix: "the program works the same wherever it is delivered, so treated districts gain nothing extra from treated neighbors".
3. **Counties versus districts in the Case 2 bullet (carried).** Text: "In Case 2, a provider … into control counties as well. If treated districts were instead expected to reinforce each other…"; first bullet "keep same-arm districts apart". Fix: "If treated areas were instead…" and "keep same-arm units apart".
4. **NC slide bullet is grammatically broken.** Text: "Some channels follow the organization chart, so borders are a proxy, but the one a design can act on". Fix: "Some channels follow the organization chart, so borders are only a proxy, but they are the one thing a design can act on."
5. **Moran's I named but never used (carried).** Spatial dependence slide and script: "it is usually measured with Moran's I". Fix: "measured with standard indices of spatial clustering". Or drop it from the script and keep it on the slide.
6. **Square-root axis not explained (carried).** Where BSS falls; the axis label says "(square-root scale)", and 0.00 to 0.50 takes up half the top row. Fix: add to the note "The scale is stretched near zero so small errors can be told apart."
7. **Negative number read aloud (carried).** Collinearity script: "about -0.35 at psi one-half". Fix: "about minus 0.35".
8. **Controls-only ρ = 0 panel can mislead.** fig_mse_2x4_controls, right panel. The SRS box median sits right at the diamond, so the picture says "SRS typical = BSS" while the speaker says SRS is lower on average (0.173 vs 0.262). The script does explain it ("69 percent … also leave every control exposed"), but "The boxes on the right reach the diamonds" undersells it. Fix: "On the right, most random draws sit right at the diamond, because 69 percent of them also leave every control exposed. The few that don't are near zero and pull the random average down."
9. **Long spoken sentences.** All are 32 to 40 words. These three will run out of breath:
   - Recommendation (40 words): "This map is a suggested, best-effort allocation … so the program and its spillover can be told apart." Split it: "This map is a suggested, best-effort allocation for the state's irregular geography, stratified by division. A refined version could also leave a few control districts with no intervention neighbor, so the program and its spillover can be told apart."
   - Collinearity (37 words, carried): split after "psi one-half".
   - Recommendation (32 words, carried): "The recommendation rests on what the proposed studies assume. A treated district gains nothing extra from treated neighbors, but nearby control districts can pick the program up."

   Also on the decision table, the top-left cell has two parentheticals in a row: "(near ties at ψ = 0.8) (2×4: 0.165 vs 0.304)". Merge them into one.

**Script read aloud:** The script has no colons, and no sentences start with "Because" ("That is because…" and "…, because…" are fine). The new Recommendation sentence "To be clear, a pure checkerboard still underestimates…" reads well and does the job round 3 asked for.

**Pacing:** 20.0 min in the prelim (title, outline and thank-you dropped), about 120 wpm overall. The fastest slides are Limitations (154 wpm), Recommendation (152), Question (144) and Spatial dependence (141). Only the Recommendation slide needs a change (M3).

**Take-home for a probation agency:** The deck now says to stratify by division (or MCO area), keep treated districts apart where the map allows, use the realistic "outcomes move together" column, and accept a reasonable error to avoid a bad draw. With M1's sentence, the agency would also know what to do once it looks at its map.
