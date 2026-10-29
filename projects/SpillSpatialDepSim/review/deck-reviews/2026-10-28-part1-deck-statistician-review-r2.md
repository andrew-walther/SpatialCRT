# Part I deck review (independent biostatistician), round 2

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.qmd, as knitted in presentation/spillover_talk.revealjs.md, with every figure it references in presentation/figures/
Sources checked against: review/manuscript/Revisions_V2c_corrected.tex; review/tables/beta_mse_table_{2x4,3x3,3x4}_corrected.csv; per-allocation CSVs in results/beta_mse, results/rho_mse, results/psi_mse and review/results (loaded read-only through presentation/tools/load_results.R); code/SpatialSim_NC_DOC.Rmd, SpatialSim_3x3.Rmd, SpatialSim_3x4.Rmd and NCTrtMaps.Rmd; review/notes/2026-10-28-rho-scale-note.md; round-1 report review/deck-reviews/2026-10-28-part1-deck-statistician-review-r1.md. Check scripts are in the session scratchpad and write nothing to the repo.

**Verdict:** Blocking issues: 2 (5 major, 10 minor)

Every number shown is reproduced from the saved results: 57% (40 of 70, the same at all four ψ), 2.39, 0.19 (BSS 0.186 and 0.144), 40–60% (40.2–60.3), 29%, 0.165/0.304, 0.569/0.573, 0.173/0.262, SD 0.41/0.03, 0.095/0.165, 0.008/0.262, bias −0.50 to −0.80 (ρ = 0) and −0.35 to −0.72 (ρ = 0.01), 48/57/488 and 2/1/2 collinear allocations, 69%, 24 of 24, 12 of 12, 16 of 16, 12 of 16, 38 of 48, 3×3 BSS 0.246 vs ≤ 0.0008, median ρ̂ 0.50, BSS ρ̂ 0.89/0.85, worst ρ = 0 draw 0.268, and MSE rising with ψ in all 12 cells. Running times add to 20.5 min through Limitations (20.75 with the thank-you). No script sentence opens with "Because", every script bullet is one sentence, and the script has no colons.

---

## Round-1 BLOCKING and MAJOR items

| r1 item | Status | Note |
|---|---|---|
| B1 Collinearity not special to BSS | Fixed | Counts per grid and type on the slide, limitation reworded. The new first bullet overstates the condition (see M3). |
| B2 "Beats most random draws" title and row 2 | Fixed | Retitled to "about half", and the script now covers row 2. The new bullet "Without spatial dependence ... most random allocations do better" is wrong on 3×4 (see M2). |
| B3 β̂ from a singular design, inference backup | Fixed | Proviso and `lagsarlm` behaviour stated, Model script reworded. It omits the 3×3 strict checkerboard (m7). |
| M1 β̂ ≈ β − ψ only at ρ = 0 | Fixed | |
| M2 Recommending BSS plus a check BSS fails | Partly fixed | The check is concrete on Recommendation and Take-home. But the NC example allocation appears to fail that check (B2), and the "Either, with a check" cell does not say that 2×4/3×4 BSS always fails it (m9). |
| M3 3×3 BSS, only one collinear | Fixed | On "Larger grids", with numbers. |
| M4 ρ scale, BSS ρ̂ and tail at ρ = 0.01 | Fixed | Glossary marker and both bullets added. Two wording problems follow from it (B1, m6). |
| M5 "Lower at every ψ" within MC error | Fixed | "gap closing by ψ = 0.8 (0.569 vs 0.573)". See M4 for the size of Monte Carlo error. |

---

## BLOCKING

**B1. "Limitations" slide and script: "correctly specified"**
Text (slide): "No intra-cluster correlation; results specific to a correctly specified spatial lag model". Script: "...the model has no intra-cluster correlation and is fitted as the true model."
Problem: The deck's own backup, "How large is ρ = 0.01?", says outcomes were generated with subject-level binary weights and fitted with row-normalized weights. At ρ = 0.01 the fitted model is therefore not the generating model, because the row sums are 40 vs 60, so the two W are not proportional. If the backup is pulled up after this slide, the committee sees the speaker contradict himself. This is a wording fix within the author's decision (results as published, weights issue in backup only). It does not reopen that decision. The corrected manuscript (line 332) uses the same phrase, but the deck does not have to repeat it.
Fix: Slide: "No intra-cluster correlation; results specific to the spatial lag model with the true x and z terms." Script: "...and the model has no intra-cluster correlation and includes the true treatment and spillover terms." Neither line then asserts that the weights are correct.

**B2. "Recommendation" slide: the example map and the check it recommends**
Text: Bullet "Real borders are irregular, so separation is partial; document how close the allocation comes, and confirm that some control districts border no intervention district". Script: "On this map, districts in the east are kept apart well, but elsewhere irregular borders put two intervention districts side by side."
Problem: (a) By eye, every blue (control) district in BSS_Assignments_NC_Map_Applications_Text_Labels.png touches at least one red district (for example 2–1/3, 5–4/13, 8–7/4, 16–12/13, 11–10/12/19B, 20–19A/19B, 25–24/23/27, 28–24/30, 29–30/27, 17–23/21). If so, the example allocation is collinear under control-only spillover, the very setting the slide uses to justify BSS, and it fails the check in the next breath. (b) The east is where intervention districts touch most visibly (3–4, 3–7, 4–12/13, 12–13), so "in the east ... kept apart well" reads as the opposite of the figure. The map is a hand-picked list (NCTrtMaps.Rmd, `treatment_districts_list`), so neither claim was computed. I could not verify adjacency offline because the shapes come from `tigris`.
Fix: Run `spdep::poly2nb()` on the district polygons from NCTrtMaps.Rmd and count controls with no intervention neighbor and intervention–intervention adjacencies. If no control is unexposed, swap two districts so that at least one or two controls are unexposed, regenerate the map, and say so in the script ("districts 28 and 29 border no intervention district"). Either way, replace the "east" sentence with the computed fact, for example "Four pairs of intervention districts still share a border, all in Divisions 1 and 3."

---

## MAJOR

**M1. "Where BSS falls": who BSS beats (readiness; no backup covers it)**
Text: "57% of random allocations have higher MSE than BSS with spillover into controls and ρ = 0.01"; "The worst random draw has MSE 2.39".
Problem: The 57% comes almost entirely from collinear random allocations. At 2×4, ψ = 0.5, ρ = 0.01, the 22 identified random allocations (some control unexposed) have median MSE 0.094 (IQR 0.048–0.131), and BSS (0.165) beats only 2 of them. On 3×4 it beats 151 of 436 identified allocations (35%). The tail is not collinearity either. The 2.39 draw is combo 1 ({1,2,3,4}, one half of the grid treated), which is identified and has bias +1.53, and this happens only at ρ = 0.01. Since the deck's own advice is "confirm some controls border no intervention cluster", a statistician will ask "against the random allocations that pass your check, how does BSS do?" Right now the answer is "worse than 20 of 22", and the speaker has nothing prepared. This does not challenge the framing. It is a question the framing invites.
Fix: Add one backup slide, "BSS against random allocations that keep β identified", with the two numbers above and the identity of the 2.39 draw. Prepare a one-sentence answer, for example: "Among random draws that leave a control unexposed, most do better than the checkerboard. The case for BSS is protection from the tail and from the 69% of draws that share its collinearity, and on the 3×3 grid the identified BSS allocations do best of all."

**M2. "Where BSS falls", third bullet**
Text: "Without spatial dependence, or with spillover into both arms, most random allocations do better".
Problem: The bullet follows "Across all grids and ψ in that setting", so it reads as covering all grids. On 3×4, Control Clusters Only, ρ = 0, BSS beats 52.7% of allocations at every ψ, and the backup rank figure shows those four squares right of the 50% line. On 2×4 at ψ = 0.8 the share is 47%, which is only just "most". The backup the speaker would pull up contradicts the slide.
Fix: "On this grid, without spatial dependence or with spillover into both arms, most random allocations do better." Or keep the scope general and say "with spillover into both arms, 71–100% of random allocations do better; without spatial dependence, BSS ties the large group of collinear draws."

**M3. "Why: collinearity" title and first bullet**
Text: Title "When every control borders a treated cluster, spillover exposure is exactly the opposite of treatment". Bullet "If every control is exposed, $z_i = 1 - x_i$ for every cluster, $[\mathbf 1, x, z]$ loses rank, and β and ψ cannot be separated".
Problem: The implication holds only under Control Clusters Only. Under Control & Intervention it also requires that no treated cluster borders another, since a treated cluster with a treated neighbor has z = 1. The title states the condition as general, and the right column then says that under Control & Intervention only checkerboards qualify, which a listener has to reconcile. Separately, $[\mathbf 1, x, z]$ also loses rank under Control & Intervention whenever z ≡ 1 (every cluster has a treated neighbor). That is 16 of 70 (2×4), 16 of 126 (3×3) and 126 of 924 (3×4) allocations. There ψ is aliased with α but β is still identified, so "loses rank" and "β and ψ cannot be separated" are not interchangeable.
Fix: Title "With spillover into controls only, any allocation that exposes every control makes exposure the exact opposite of treatment". Bullet "Under Control Clusters Only, if every control is exposed, $z = 1 - x$ ...". For backup readiness, add one line to "Reduced form and likelihood": "Under Control & Intervention, z ≡ 1 in 16/16/126 allocations; ψ is then aliased with α, and β is unaffected."

**M4. SD of MSE for BSS on 2×4 and 3×4 is Monte Carlo noise ("Control clusters only" card, "Larger grids" bullet, table backups)**
Text: "SD of MSE 0.41 for SRS vs 0.03 for BSS"; "BSS error varies less across allocations in 38 of 48 scenarios".
Problem: On 2×4 and 3×4 the two BSS allocations are complementary checkerboards ({1,4,5,8} and {2,3,6,7}; combos 313 and 612 on 3×4), which map onto each other by reflection. Cluster sizes are exchangeable, so the two allocations have the same MSE in distribution, and the BSS "SD across allocations" on these grids measures only Monte Carlo error. Twenty-four of the 38 rest on this comparison. A statistician will ask what an SD of two numbers means. The gap also measures Monte Carlo error. At 2×4, ψ = 0.5, ρ = 0.01, the two equivalent allocations give 0.186 and 0.144, a difference of 0.04 (about 25%) at N = 50, and at 3×4 (N = 10) 0.179 and 0.088. That bears directly on the 0.569 vs 0.573 comparison and on near-tie rank shares.
Fix: On "Control clusters only", say "BSS is effectively one allocation on this grid (the two checkerboards are mirror images), so its spread is near zero by construction, against an SD of 0.41 across random draws." On "Larger grids", either drop "38 of 48" or restrict it to the 3×3 grid, the only one where BSS has distinct allocations. Prepare a backup line: "Equivalent checkerboards differ by 0.04 at N = 50, a guide to the Monte Carlo error."

**M5. Backup "Inference", second bullet**
Text: "That is by design: BSS removes the high-error allocations SRS can draw, bounding the worst case".
Problem: This is true only with control-only spillover at ρ = 0.01. With Control & Intervention at ρ = 0, the checkerboards are the worst allocations. On 2×4 the SRS maximum is 0.268, which is BSS combo 28, and on 3×4 it is 0.242, also BSS. With Control Clusters Only at ρ = 0, 2×4 BSS beats 1.4% of draws. Stated generally, the claim contradicts the main results.
Fix: "BSS rules out the clustered allocations that produce the worst errors when spillover goes into controls and outcomes are spatially correlated. With spillover into both arms, the checkerboards are themselves the worst allocations."

---

## MINOR

**m1. Backup "Simulation algorithm", step 3.** Text: "sample 20 subjects per cluster". The code draws `points_per_iteration` (160/180/240) points uniformly over the whole grid (`sample(1:grid_width-1, ...)`), so cluster sizes are Multinomial with mean 20. The Model slide already says "about 20 per cluster". Fix: "sample 160, 180 or 240 subjects over the grid (about 20 per cluster)". Step 2 says "every balanced allocation (70, 126, 924)", but 3×3 is 4 vs 5. Fix: "every allocation with ⌊K/2⌋ clusters treated".

**m2. Timing.** "NC probation districts" runs at 201 wpm (151 words in 0.75 min). "Recommendation" (151 wpm) and "Limitations" (152 wpm) are also fast against a deck average of about 115. Total is 20.5 min plus 0.25 for the thank-you. Fix: give the NC slide 1 min and cut its script to about 120 words (drop the "Randomizing individual officers isn't feasible" sentence, which the background already covered), and trim one sentence from Limitations so the total stays near 20.

**m3. fig_mse_larger_grids.png.** The x-axis labels of adjacent facets collide ("0.80.5"). Fix: `scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8), labels = c("0.5", "", "0.7", ""))` or add `panel.spacing.x`.

**m4. Model slide glossary, ε.** Text: "Independent normal error (σ = 0.1)" on a model written for the cluster-level outcome $y_i$. σ = 0.1 is the subject-level error SD. A cluster mean of about 20 subjects has error SD near 0.1/√20 ≈ 0.022. This is wording only, within the cluster-level decision. Fix: "Independent normal error; subject-level SD σ = 0.1".

**m5. Backup "Reduced form and likelihood".** The same **W** appears in the generating equation and in the likelihood, which contradicts the ρ backup. Fix: add "(the simulation's generating W is described in the ρ backup)" under the first equation.

**m6. Backup "How large is ρ = 0.01?", fourth bullet.** Text: "the lag term absorbs part of the dropped spillover, which is why BSS bias shrinks with spatial dependence". This is plausible, since median ρ̂ is 0.60 for collinear and 0.38 for identified 2×4 allocations, and BSS ρ̂ is about 0 at ρ = 0. But the rho note's rerun did not reproduce the BSS value (0.253 vs 0.165) and calls the checkerboard fits numerically unstable. Fix: "...which is consistent with BSS bias shrinking with spatial dependence".

**m7. Backup "Inference", fourth bullet.** Text: "For the 2×4 and 3×4 checkerboards it does not". The 3×3 strict checkerboard {2,4,6,8} is also collinear (MSE 0.246). Fix: "For the checkerboards (both BSS allocations on 2×4 and 3×4, one of six on 3×3) it does not".

**m8. Colons in slide text.** Title "The question: does the allocation...", the collinearity labels "Control Clusters Only:" and "Control & Intervention Clusters:", inference backup "That is by design:", algorithm "Fix the scenario:", and three in the ρ backup ("binary weights: each subject", "0.85: the lag term", "dependence: the worst"). Fix: for example, title "Does the allocation chosen before the trial change how well we estimate the effect?", and replace the colons with commas or periods.

**m9. "Choosing a design", ρ = 0 control-only cell.** Text: "Either, with a check ... make sure some controls border no intervention cluster". On 2×4 and 3×4, BSS always fails that check, so in this cell "either" means "SRS, or BSS only where the map leaves a control unexposed". Fix: "SRS, or BSS only if some controls border no intervention cluster (never true of a strict checkerboard)".

**m10. Carried over from r1 m7 and m14.** (a) The NC map uses red/blue while every other schematic uses navy/grey/amber, so "red = intervention" is the one colour key the audience has to relearn. Fix: regenerate in the deck palette or say "red, intervention" aloud. (b) The DOI 10.1186/s12874-026-02968-0 still appears in no repo file. Check it against the publisher before the talk.
