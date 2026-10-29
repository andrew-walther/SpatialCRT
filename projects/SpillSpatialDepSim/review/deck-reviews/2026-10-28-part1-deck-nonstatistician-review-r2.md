# Part I deck review (non-statistician, criminal-justice/behavioral-health reader), round 2

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.revealjs.md (Draft 1 with round-1 revisions, 2026-10-28) and the figures it shows in presentation/figures/

**Verdict:** Blocking issues: 1 (4 major, 8 minor)

This draft is much easier to follow than round 1. Probation now shows up in the first minute. MSE is explained in plain words where it first appears. The dot plot has an orientation line, and the decision table uses words I can map onto my own setting. The one place I'd lose the thread is the recommendation map. The script describes it in a way the picture doesn't support, and the map seems to fail the check the talk asks agencies to run.

---

## Round-1 BLOCKING and MAJOR items

| Item | Status | What changed / what remains |
|---|---|---|
| B1. Recommends the design just shown to be broken | Partly fixed | The Recommendation script now ties BSS to the setting where it did best ("beat about half of random draws and never landed in the tail"), and the deck asks investigators to keep some controls unexposed. Still missing: the deck never says out loud that even in its best setting the checkerboard underestimates the effect by about a third (bias −0.35 on β = 1). See M1 below. |
| B2. Agency can't place itself in the table or run the check | Fixed | (a) Explaining ρ = 0.01 in a backup slide is the author's decision, and the backup slide handles it well. (b) The columns now read "Neighboring outcomes move together / independent", and the note says the first column is the realistic one for probation. (c) The check is now an action: "confirm that some controls border no intervention cluster." |
| M1. Probation appears only at minute 16 | Fixed | Opening bullet and script name the two NCDAC trials. The question slide links back ("what the training does … its spread to nearby offices"). |
| M2. Case 1 spillover channel vs. rook neighbors | Partly fixed | Added "Some channels follow the organization chart, so borders are a proxy." Still unresolved: which unit is randomized ("county or district"), and so which leader's reach actually crosses arms. See M3 below. |
| M3. Why "controls only" is expected isn't argued | Partly fixed | The sentence about trained districts already getting the full program, plus what would change the call, are both good for Case 1. That reasoning is about officer training, though, and doesn't carry over to Case 2. See M2 below. |
| M4. Recommendation map not walked through | Partly fixed, and the walkthrough adds a new problem | The map note and a walkthrough line were added. The walkthrough says the east is well separated, and the picture shows the opposite. See B1 below. |
| M5. Dot plot hard to read cold, row 2 unexplained | Fixed | "Further left is better… each row on its own scale" is both on the slide and in the script. Row 2 is explained. The 29% figure moved to backup. The slide gets 1.75 min. |
| M6. Controls-only box plot reads as bars; spread quoted for both panels | Fixed in script | The spread figure is now limited to ρ = 0.01. The right panel is explained (69% of random draws are collinear and match the checkerboard). The figure itself still looks like bars, but the script covers it. |
| M7. MSE not explained at first use | Mostly fixed | The plain-words definition is on the slide and in the script. Still no anchor for what a given MSE means. See Minor 1. |
| M8. Application too short | Partly fixed | The Application now runs 2.5 min (was 2.25), taken from "The design gap". But the main talk now runs 20.5 min. See M4 below. |

---

## BLOCKING

### B1. The recommendation map contradicts the script, and it seems to fail the deck's own check
- **Slide:** "For both studies we recommend a best-effort block stratified allocation" (BSS_Assignments_NC_Map_Applications_Text_Labels.png)
- **Text (script):** "On this map, districts in the east are kept apart well, but elsewhere irregular borders put two intervention districts side by side" and "Investigators should also confirm that some control districts border no intervention district"
- **Problem:** On the figure, the east is where intervention districts cluster most. Red 7, 3 and 4 touch each other, 4 touches 12 and 13, and 12 touches 13. That makes a red chain across Divisions 1 and 2. The west (30 / 28–29 / 24 / 25 / 27 / 26) looks better separated. Worse, as far as I can tell from the picture, every blue district borders at least one red one: 2, 5, 6, 8, 9, 11, 14, 16, 17, 18, 20, 22, 25, 26, 28 and 29 all appear to touch a red district. If so, this example allocation is exactly the "every control exposed" case the collinearity slide warns about. The audience is told "check this before the trial" while looking at a map that would fail the check, and nobody says what to do when it fails. For me that breaks the take-home. The map comes from a hand-coded list in `code/NCTrtMaps.Rmd` (treated: 1, 3, 4, 7, 10, 12, 13, 15, 19A, 19B, 21, 23, 24, 27, 30). I checked it only by eye, so it should be confirmed with `spdep::poly2nb` on the district polygons.
- **Fix:** (1) Run the adjacency check on the district polygons. (2) If at least one control borders no intervention district, name it on the slide ("e.g., district N borders no intervention district") and fix the east/west wording to match the picture. (3) If no control qualifies, say so honestly and use it as the teaching moment: "This published example fails the check, since every control touches a red district. In practice we'd swap assignments within a division until at least one control in each division has no treated neighbor." Add a one-line "if the check fails" step to the "State the mechanism" take-home card, such as "if it fails, re-draw within strata until it passes". Don't edit the published figure. Put the new allocation, or a highlighted unexposed control, in the talk version only, if the author approves.

---

## MAJOR

### M1. The checkerboard's bias in its best setting is never said out loud on the main slides
- **Slides:** "When every control borders a treated cluster…" and "For both studies we recommend…"
- **Text:** "At ρ = 0 the 2×4 BSS bias is −ψ (−0.50 to −0.80); with spatial dependence the lag term absorbs part of it (−0.35 to −0.72)"
- **Problem:** The bias numbers are on the collinearity slide, but the script only says "somewhat less with it". An agency director will ask: "In the setting you recommend, is my estimate still about a third too low?" At ψ = 0.5 and ρ = 0.01 the answer is yes (−0.35 on 1.0). The deck's framing is good, not best, against rolling the dice. That framing can survive this, but only if the speaker says it first.
- **Fix:** Add one spoken sentence on the Recommendation slide: "To be clear, a pure checkerboard still underestimates the effect, by about a third in this setting. It wins because random draws can be much worse, not because it is unbiased. Keeping some controls unexposed is what removes that bias."

### M2. The "spillover into controls only" rationale fits Case 1 but not Case 2
- **Slide:** Recommendation, "Why controls only: trained districts already receive the full program, so the main risk is untrained staff picking it up"
- **Problem:** Case 2 is about cross-training and joint treatment-team meetings between officers and behavioral health providers. The Case 2 slide itself says large providers work across counties and MCO areas. A provider who sits in team meetings in two treated counties gets a double dose and carries practice between them. That is spillover among treated units, the case where the deck says SRS is the better choice. A behavioral-health implementation researcher will spot this right away, since the recommendation is "for both studies".
- **Fix:** Give Case 2 its own sentence. For example: "For the mental health study, the providers are the channel. If the same providers serve several intervention counties, treated counties reinforce each other, and the case for random allocation gets stronger. Mapping which providers serve which counties before randomizing tells you which case you are in." If the author still recommends BSS for Case 2, say it is a closer call than Case 1.

### M3. The leadership spillover channel doesn't match the unit being randomized, and stratifying by division opens a new channel
- **Slide:** "North Carolina probation trials assign judicial districts…"
- **Text:** "Officers supervise by county; chief officers and judicial district managers often span several counties" / "randomized by county or district" / "Leaders can carry the training into control areas they oversee"
- **Problem:** If the trial randomizes districts, which is what the recommendation map shows, the judicial district manager channel stays inside one cluster and can't reach a control. The channel crosses arms only if counties are randomized. At the district level, the leader who oversees both arms is the Judicial Division Director. Stratifying by division guarantees that every division director oversees both intervention and control districts. A probation audience will notice that the stratification creates a contamination channel through the division office. Also, the manuscript says training goes to cohorts of PPOs and CPPOs. Trained chiefs are the people most likely to "carry it" into a control unit, and the slide mentions only officers.
- **Fix:** Pick the unit for Case 1 on this slide ("randomized by district, as in the map"), then name the channel that fits it: "Division directors oversee both arms by design, and neighboring districts share courts, providers and staff transfers." If county randomization stays an option, add: "If counties are randomized, a district manager can oversee both arms." Change "Training for officers" to "Training for officers and their chiefs."

### M4. The main talk runs over 20 minutes with no slack
- **Slides:** Limitations "[0.5 min | running 20.5]", then Thank you "[0.25 min]"
- **Problem:** The take-home lands at exactly 20.0, Limitations adds 0.5, and the thanks adds another 0.25. That's 20.75 minutes before any slip. The script is about 2,600 words, roughly 126 words per minute, which is comfortable per slide but leaves no room. The likeliest overruns are the Model slide (1.5 min, 40-word sentences) and the dot plot (1.75 min).
- **Fix:** Fold the one Limitations point that matters (collinearity) into the take-home script, and show Limitations silently for 10–15 seconds while saying "the limitations and next steps are here, and the paper is in BMC MRM." Or trim the Model script by cutting the "We sample subjects within clusters…" sentence, which the slide's own note already covers, and the "structural rather than potential-outcomes" sentence, which the backup covers. Either way, aim for 19.5 minutes.

---

## MINOR

1. **No anchor for MSE size.** Slides "Where BSS falls", "Control clusters only". I still can't tell whether 0.165 is good. Fix: say once, "an MSE of 0.16 means the estimate is typically off by about 0.4, on a true effect of 1." This also sets up M1.
2. **Decision table, controls only with independent outcomes.** Text: "SRS is lower on average on the 2×4 and 3×4 grids, BSS on the 3×3". An agency can't tell which grid it is. Fix: "Either. What decides it is whether the allocation leaves some controls unexposed, so run the check." Move the grid detail to the script.
3. **Unit inconsistency on the opening slide.** Text: "two planned probation trials in North Carolina that randomize districts". Case 1 says county or district, and Case 2 randomizes counties. Fix: "that randomize counties or districts".
4. **Why not randomize officers.** NC slide script: "Randomizing individual officers isn't feasible, so it randomizes counties or districts." A probation listener knows why, but the committee may not. Fix: "Officers in the same office share cases and supervisors, so randomizing officers would contaminate within the office."
5. **Chief role description.** Text: "chief officers and judicial district managers often span several counties". In NC a CPPO leads a unit of officers, usually inside one county office. Spanning counties happens mainly in rural multi-county districts, and the JDM is the district-level leader. Fix: "Chiefs lead units of officers; judicial district managers lead each district, which in rural areas spans several counties." Use the title "chief probation/parole officers" once.
6. **The name "block stratified" is still unexplained, and "stratified" means something else to social scientists.** The question slide card. Fix: add, as a spoken phrase, "'block' for the checkerboard of separated blocks; in the application it is also stratified by division or MCO area."
7. **A remaining colon in a slide title.** "The question: does the allocation chosen before the trial change how well we estimate the effect?" Fix: "Does the allocation chosen before the trial change how well we estimate the effect?"
8. **Long sentences to split when read aloud.** No sentence starts with "Because", and "since" is down to two uses, which is good. These four still run past about 35 words.
   - Model: "A cluster's outcome is a baseline alpha, plus the intervention effect beta if the cluster is treated, plus the spillover effect psi if it is exposed, plus rho times…" Fine for the statisticians, but pause after "exposed".
   - Where BSS falls: "Each grey dot is one random allocation… and further left is better, with each row on its own scale." Split it: "Each grey dot is one random allocation, and the diamonds are block stratified. Further left is better, and each row has its own scale."
   - Collinearity: "The software drops the spillover term, and the estimate comes out near beta minus psi, which on the two-by-four grid is a bias of minus psi…" Split after "beta minus psi". Then use M1's sentence with a number.
   - Control clusters only: "The boxes on the right stretch up to the diamonds, since most random allocations, 69 percent of them on this grid, also leave every control exposed…" Change to "The boxes on the right reach the diamonds. That's because 69 percent of random allocations on this grid also leave every control exposed, and they do just as badly."

**Pacing summary:** the per-slide timings are realistic for a mixed audience except the Model slide, which is dense for 1.5 minutes. The only real pacing problem is the 20.5–20.75 minute total (M4). The Application now gets enough time, provided the map walkthrough is corrected (B1) rather than lengthened.
