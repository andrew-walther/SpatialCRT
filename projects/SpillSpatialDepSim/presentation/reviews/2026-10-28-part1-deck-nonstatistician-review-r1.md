# Part I deck review (non-statistician, criminal-justice/behavioral-health reader), round 1

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.revealjs.md (Draft 1, 2026-10-28) and the figures it shows in presentation/figures/

**Verdict:** Blocking issues: 2 (8 major, 9 minor)

Overall, the opening works. Bed nets, deworming and "a supervision practice in a probation office" told me why a trial of places has to worry about neighbors. I followed the two designs and the checkerboard picture. I lose the thread at two points. The deck recommends the design whose main flaw it has just shown, without reconciling the two. And the take-home turns on a ρ value that sounds like zero, plus a "collinearity check" an agency could not carry out.

---

## BLOCKING

### B1. The deck recommends the design it just showed to be broken, and never says why that's still fine
- **Slides:** "In a checkerboard, spillover exposure is exactly the opposite of treatment" → "For both studies we recommend a best-effort block stratified allocation" → Limitations ("Treatment and spillover can be exactly collinear under BSS")
- **Text:** "the model can't tell the intervention effect from the spillover effect, and the estimate absorbs the spillover, coming out near beta minus psi" … then "For both studies we recommend block stratified sampling" … then "The most important [limitation] is the collinearity".
- **Problem:** As a non-statistician I hear three things in a row. The checkerboard can't separate the program from its spillover. Use it anyway. Its biggest weakness is that it can't separate the program from its spillover. Nothing tells me why this flaw doesn't sink the recommendation in the controls-only, spatially correlated case. The backup table makes this sharper, since even where BSS wins (2×4, ρ = 0.01, ψ = 0.5) its bias is −0.35 on a true effect of 1.0. That's roughly a one-third underestimate, and an agency will ask about it.
- **Fix:** Add one or two spoken sentences on the collinearity slide or the Recommendation slide that bridge the two. For example: "So why recommend it? When spillover reaches only controls and outcomes are correlated, the checkerboard still lands closer to the truth than most random draws, and it never lands in the bad tail. And on a real map, best-effort separation leaves some control districts with no treated neighbor, which gives back some of the information the checkerboard loses." Check that the second clause is defensible before saying it, since the paper doesn't test irregular maps. If it isn't, say plainly that BSS trades some bias for protection against a very bad draw.

### B2. A probation agency can't place itself in the decision table or carry out the recommended check
- **Slides:** "Match the design to the spillover and spatial correlation you expect", and the Take-home card "State the mechanism"
- **Text:** column headers "Spatially correlated outcomes (ρ = 0.01)" / "No spatial dependence (ρ = 0)"; "check any allocation for treatment–exposure collinearity before the trial"
- **Problem:** (a) To a lay ear, ρ = 0.01 sounds like essentially no correlation, so a 0.01 correlation appearing to flip the answer looks fragile or odd. The fact that it is moderate on the fitted scale (about 0.5) appears only in a backup slide. (b) An agency doesn't know its ρ. It needs a plain-language cue for which column applies. (c) "Treatment–exposure collinearity" isn't something a district manager or study coordinator could check.
- **Fix:** (a) On the decision slide say once, aloud: "ρ = 0.01 sounds tiny, but on the scale of the fitted model it is moderate correlation, about one half." This keeps the published values and stops the misreading. (b) Relabel the columns, or add a subtitle, in plain words, for example "Neighboring districts' outcomes tend to move together (shared courts, providers, labor markets)" versus "Neighboring districts behave independently". Then say that for NC probation the first column is the realistic one. (c) Restate the check as an action: "Before randomizing, look at the map and confirm that at least some control districts don't border any intervention district. If every control touches a treated district, the trial can't separate the program from its spillover."

---

## MAJOR

### M1. Probation, the motivating case, doesn't appear until minute 16
- **Slide:** "When a program is delivered to whole communities, trials randomize places…"
- **Text:** "such as a supervision practice in a probation office"
- **Problem:** The real stakes (two NCDAC studies that need an allocation) arrive only in the Application. For 15 minutes the talk is about abstract grids, and I can't connect the results to a decision anyone faces.
- **Fix:** Add 2–3 spoken sentences at the opening, for example: "This came from two planned probation trials with the NC Department of Adult Correction. They randomize districts, a manager or a treatment provider can carry the practice next door, and the investigators asked us which districts to treat." Return to it in one line on the question slide.

### M2. The Case 1 spillover channel doesn't match the model's definition of a neighbor
- **Slide:** "North Carolina probation trials assign judicial districts…"
- **Text:** "Probation officers supervise by county; chiefs and district managers oversee several counties" / "A manager can oversee both intervention and control counties, so spillover into controls is plausible"
- **Problem:** A probation researcher would raise two points. (1) If the trial randomizes *districts*, the judicial district manager's channel stays inside one cluster. It only crosses arms if the trial randomizes *counties*, but the slide says "county or district". (2) The channels described (a shared manager, large providers working across counties and MCO regions) follow the organization chart and the provider network, not shared borders. The model's spillover runs only to rook (shared-edge) neighbors. The audience may wonder whether separating districts on the map actually blocks those channels.
- **Fix:** Say which unit is randomized in Case 1 and which channel crosses arms at that unit. Add one sentence acknowledging that some channels are administrative rather than geographic. Borders are the proxy the design can act on, and the limitation is that management and provider networks don't follow borders exactly. If CPPOs supervise units within a county rather than several counties (the manuscript says CPPOs and JDMs "oversee districts"), correct the slide's "chiefs … oversee several counties".

### M3. Why "spillover into controls only" is the expected case for probation isn't argued
- **Slides:** Recommendation, "The recommendation rests on the expected mechanism, spillover mainly into control districts"
- **Problem:** The whole recommendation flips on this assumption, but the talk only asserts it. For an officer-training intervention, treated districts can plausibly reinforce each other through shared trainers, regional booster sessions and peer coaching among trained officers. A probation researcher will ask why that wouldn't count as spillover into treated districts.
- **Fix:** Give one sentence of substantive reasoning. For example: "Trained districts already receive the full program, so contact with another trained district adds little. The main risk is untrained staff picking up the practice." Then say what would change the call, such as joint training cohorts or booster sessions that span treated districts.

### M4. The recommendation map shows many same-arm borders, and the script doesn't walk through it
- **Slide:** "For both studies we recommend a best-effort block stratified allocation" (BSS_Assignments_NC_Map…png)
- **Problem:** On the map, several red districts touch each other (3–4–7, 12–13, 23–24–30) and several blue ones do too. Without a guided look, the audience sees a map that doesn't look "kept apart" and wonders whether BSS was applied at all. The map is also of judicial districts, while Case 2 randomizes counties within MCO areas, so it illustrates only Case 1.
- **Fix:** Spend 15 seconds on the map. Point to one division where separation works, point to one place where borders force two red districts together, and say "this is what 'document how close you come' means." Add to the note: "Example for Case 1 (districts); Case 2 would do the same with counties within each MCO area."

### M5. The "Where BSS falls" dot plot is hard to read cold, and its second row goes unexplained
- **Slide:** "A trial runs one allocation: with spillover into controls, block stratified sampling beats most random draws" (fig_spread_2x4.png)
- **Problem:** Each row has its own x-axis range (0–2.5, 0–0.25, 0–1.0), on a square-root scale. A lay viewer compares dot positions across rows and gets it wrong. The script skips row 2 (controls only, ρ = 0), where both diamonds sit to the right of almost every random allocation. That's the most striking pattern on the slide, and I'd be puzzled. 1.5 minutes is too short for a mixed audience to get oriented and also absorb 57%, 2.39 vs 0.19, 40–60% and 29%.
- **Fix:** Add a spoken orientation line: "Further left is better. Each row has its own scale, so compare diamonds to the dots in the same row only." Add one sentence for row 2: "Without spatial dependence, the checkerboard is worse than almost every random draw, the collinearity I'll show shortly." Consider leaving the "at most 29%" figure to the Discussion slide.

### M6. The controls-only box plot reads as a bar chart at ρ = 0, and the spread number is quoted for both panels
- **Slide:** "With spillover into controls only, block stratified sampling has lower error when outcomes are spatially correlated" (fig_mse_2x4_controls.png)
- **Text (script):** "In both panels, the random allocations spread far more, with a standard deviation of MSE of 0.41 against 0.03."
- **Problem:** In the right panel, the SRS boxes run from 0 to about the BSS diamond, with the median at the top, so they look like bars that reach the diamond. A lay viewer concludes SRS equals BSS there, which contradicts "SRS lower on average, 0.173 vs 0.262". The 0.41 vs 0.03 figure is the ρ = 0.01 value. At ρ = 0 the table gives 0.117 vs 0.009, so "in both panels" doesn't match the number quoted.
- **Fix:** In the script: "On the right, most random allocations sit near zero and a group sits near the checkerboard, so the box stretches all the way up, but the random average is lower." Change "In both panels … 0.41 against 0.03" to "With spatial dependence, the spread is 0.41 against 0.03, and it is similarly lopsided without it." Consider also marking the SRS mean in the figure.

### M7. MSE isn't explained in plain words when it first appears
- **Slide:** "The question: does the allocation chosen before the trial change how well we estimate the effect?"
- **Text:** "We compare the two by the mean squared error (MSE) of the intervention effect estimate"
- **Problem:** The plain meaning arrives three slides later, and in technical form ("bias squared plus variance"). Every result after that is an MSE number, so I'm reading numbers without a sense of what's good.
- **Fix:** At first use, say: "Mean squared error is how far off the estimate tends to be, combining being systematically off and being off by chance. Smaller is better, and zero is perfect." In the results, give one anchor, such as "an MSE of 0.25 means the estimate is typically off by about half the true effect."

### M8. The Application gets about 2.25 minutes, too little for the part that grounds the take-home
- **Slides:** the three Application slides ([0.75], [0.5], [0.75])
- **Problem:** The MCO slide gets 30 seconds for a new acronym, a new map and a stratification rationale. For a mixed committee this is where relevance gets established.
- **Fix:** Take 30–45 seconds from "Spatial methods mostly adjust after the data are in…" (eight bullets and seven citations in one minute) and give it to the MCO and Recommendation slides, mainly for M4.

---

## MINOR

1. **Outline**, "Neighbors, the two allocation designs…". *Allocation* and *balanced allocation* are used before they're defined. At first use, say "an allocation is the list of which districts get the program; balanced means equal numbers in each arm."
2. **"The question" slide.** The name *block stratified* is never motivated, and to a social scientist "stratified" means something else (strata by covariates). Add one phrase, for example "'block' because treated areas are separated into a checkerboard-like pattern."
3. **Spatial dependence slide**, "Measured by Moran's I or Geary's C, and modeled through a spatial weights matrix". This is fine on the slide for statisticians, but the script could drop it or add "a weights matrix just records who neighbors whom".
4. **MCO slide**, "4 LME/MCOs … one per county". This reads as four organizations each serving one county. Use "each county is served by one of 4 LME/MCOs".
5. **Case 1 slide.** The outcome is never named. Add "outcomes are officer knowledge and attitudes and supervisee outcomes" and the caseload examples (domestic violence, sex offense, gang) so it sounds like real probation work.
6. **Colons in slide text**: the title "A trial runs one allocation: with spillover into controls, …", the bullet "Exposure is **binary**: a cluster is exposed if…", and the cards "ρ = 0.01:" / "ρ = 0:". For the title, consider "A trial runs one allocation, and with spillover into controls…". The bullet can use a comma or period. The cards are fine as labels.
7. **Long or stiff spoken sentences** (no sentence starts with "Because"):
   - Background: "Two different things can then tie neighboring areas together, spillover and spatial dependence, and the rest of this talk is about how the choice of which areas to treat changes our ability to estimate the effect when both are present." Split it: "Two things tie neighbors together, spillover and spatial dependence. This talk asks how the choice of which areas to treat changes what we can learn when both are present."
   - Where BSS falls: "The bottom two rows are the honest counterpoint, since when spillover reaches both arms most random allocations do better, and at most 29% are worse than block stratified sampling, for a reason I'll show in two slides." Drop "honest counterpoint" and split: "The bottom rows are the other side. When spillover reaches both arms, most random allocations do better. I'll show why in two slides."
   - NC slide: "A manager who likes the training can carry it into counties assigned to control, which is spillover into controls, and neighboring counties share courts, labor markets and providers, which is spatial dependence." Split into two sentences.
   - "since" is used as a causal connector about ten times ("since a trial runs one allocation", "since it is caused by the treatment", "since the organizations' own policies…"). Aloud it sounds formal, so vary it with "so", or split the sentence.
8. **Model slide**, "beta is the structural intervention effect rather than a potential-outcomes contrast". This is fine for the statisticians but opaque to me. It's acceptable as is. Optionally add "in plain terms, the effect of being treated, holding exposure fixed".
9. **"The design gap" slide.** The boxed sentence ("Even a model with both own treatment and neighboring treatment cannot separate their effects if the allocation gives little independent variation in exposure…") is the one idea a lay listener needs from this slide. Give it a probation example, for example "if every control district borders a trained district, we can't tell training from proximity to training."

**Pacing summary:** the total is 20.0 minutes with no slack. The slides at risk of running long with a mixed audience are "Where BSS falls" (1.5 → about 2 min), the controls-only box plot (1 → about 1.25 min with the M6 orientation) and the Recommendation map (0.75 → about 1 min with M4). Recover the time from "The design gap" (1 → 0.5 min) and by shortening the Limitations script.
