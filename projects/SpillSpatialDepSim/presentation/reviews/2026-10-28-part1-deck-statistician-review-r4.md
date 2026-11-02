# Part I deck review (independent biostatistician), round 4

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.qmd, as knitted in presentation/spillover_talk.revealjs.md (knitted 23:50, in sync with the source: same slide titles and script bullets), with every PNG it references in presentation/figures/ opened (13, including the redrawn fig_mse_larger_grids.png)
Sources: review/manuscript/Revisions_V2c_corrected.tex; review/tables/beta_mse_table_{2x4,3x3,3x4}_corrected.csv; per-allocation CSVs in results/beta_mse, results/psi_mse and review/results, read through presentation/tools/load_results.R (read-only, its two cross-checks pass); collinearity flags from the frozen rook lists in code/SpatialSim_*.Rmd; review/notes/2026-10-28-rho-scale-note.md; presentation/tools/make_figures.R; rounds 1–3 in review/deck-reviews/. Check scripts ran from the session scratchpad and wrote nothing to the repo.

**Verdict:** Blocking issues: 1 (2 major, 15 minor)

Recomputed and confirmed this round: 57% (57.1 at every ψ on 2×4), 40–60% (40.2–60.3), worst tenth 1.20 (7 allocations) and worst 2.39 on 2×4, BSS 0.144/0.186, 0.165 vs 0.304, 0.173 vs 0.262, SD of MSE 0.41, 0.569 vs 0.573, 0.095 vs 0.165, 0.008 vs 0.262, 24 of 24, 48/57/488 and 2/1/2 collinear, 16/16/126 with z ≡ 1, 69% (48/70), 3×3 0.147 vs 0.84 and 3×4 0.133 vs 0.77, 3×4 worst 2.57 and BSS ≤ 0.18 (0.179), 2 of 22 and 151 of 436, worst draws identified on every grid, 16 of 16 (3×4) and 12 of 16 (3×3), and the ρ backup's 0.27 (a checkerboard, 0.268). Timing is 20.0 min on content slides (20.5 with title and outline). Every script bullet is one sentence (the thank-you bullet is one sentence; the checker's false hit is "Dr. Lin"), and the script has no colons.

---

## Round-3 BLOCKING and MAJOR items

| r3 item | Status | Note |
|---|---|---|
| B1 Clipped facet strips on fig_mse_larger_grids.png | Fixed | Redrawn with one strip per column ("3×3, ρ = 0.01", "3×3, ρ = 0", ...), all readable. The new figure shows every allocation (boxes) and each BSS allocation (diamonds) and agrees with the slide text: five 3×3 diamonds near zero at ρ = 0 under Control Clusters Only, both 3×4 diamonds ≤ 0.18 against a tail near 2.57. |
| M1 Recommendation bullet named the wrong improvement | Fixed | Now "A refined allocation could also leave a few control districts with no intervention neighbor, which keeps the program and its spillover separable"; the script matches. Map kept by author decision. |
| M2 ρ backup misidentified the 0.27 | Fixed | Now names the 2.39 allocation, its 0.007 at ρ = 0, and that 0.27 is a checkerboard. Verified. |

Round-3 minors: m3, m4 and m12 fixed; m2 and m11 partly fixed; m1, m5–m10 and m13 not fixed; m14 (DOI) not verifiable from the repo. The open ones are carried below.

---

## BLOCKING

**B1. Backup slide title says "Corrected"**
Slide: first backup. Text: "Backup · Corrected 3×3 and 3×4 results, every allocation".
Problem: The author's decision is that results are presented as published with no mention of a correction. This title, shown on screen whenever the slide is pulled up, invites "corrected from what?" in front of the committee and contradicts the rest of the deck. (The image file names also say `_corrected`, but those are not visible.)
Fix: "Backup · 3×3 and 3×4 results, every allocation". Leave the file names alone or rename them; only the title is visible.

---

## MAJOR

**M1. Inference backup, second bullet**
Text: "With spillover into controls and correlated outcomes, BSS rules out the clustered allocations that give the worst errors; with spillover into both arms, the checkerboards are themselves the worst allocations".
Problem: True at ρ = 0 only. Under Control & Intervention Clusters with ρ = 0.01, random allocations exceed the worse checkerboard's MSE in 3 of 70 (2×4) and 4 of 126 (3×3) allocations at every ψ, and in 203 of 924 on 3×4 at ψ = 0.5 (35 at ψ = 0.8). The deck's own backup figure (fig_bss_rank_all_grids.png) shows BSS beating 6–29% of draws in that row, so a committee member can catch the contradiction on the next slide back.
Fix: "with spillover into both arms, the checkerboards are the worst allocations without spatial dependence and sit among the worst with it (worse than 71–94% of random draws)".

**M2. The "beats X% of random draws" measure counts near-ties, and the backup figure shows it**
Slides: "Where BSS falls" third bullet ("On this grid, without spatial dependence or with spillover into both arms, most random allocations do better"), its script ("most random draws do better"), and backup "BSS rank, all grids" (figure and notes).
Problem: Without spatial dependence, the collinear allocations all estimate β − ψ and their MSEs are near-identical (2×4, ψ = 0.5: 0.244 to 0.268, against BSS 0.255 and 0.268). Whether one counts as "better" or "worse" than BSS is decided by Monte Carlo noise. As a result the backup figure shows BSS beating 53% of random draws on 3×4 (Control Clusters Only, ρ = 0, every ψ; 487 of the 488 beaten are collinear ties, 10 replications) and 47% on 2×4 at ψ = 0.8, even though SRS has the much lower mean in both (3×4 0.133 vs 0.226). The backup notes cover only ρ = 0.01 and both arms, so the answer is not prepared. On the main slide, "most random allocations do better" at ρ = 0 is numerically true (69 of 70 below the BSS mean) but 47 of those 69 are ties. Only 22 do clearly better. At ρ = 0.01 the collinear allocations spread out (0.012 to 0.776), so the 57% headline is not a tie artefact, but 38 of the 40 beaten draws are collinear allocations. The Downside-risk backup already says so through "2 of 22".
Fix: (a) Add to the rank backup's notes: "Without spatial dependence, the collinear allocations tie with the checkerboards up to simulation noise, so this share is not meaningful in the ρ = 0 rows. For example, on 3×4 the 53% is ties, and SRS has the lower mean (0.133 vs 0.226)." Consider greying out or annotating the two ρ = 0 rows in the figure. (b) Main-slide bullet: "On this grid, without spatial dependence, a third of random draws do far better and the rest tie with it; with spillover into both arms, most do better." Adjust the matching script sentence.

---

## MINOR

**m1. Timing on two slides (carried).** "Recommendation" is 147 words in 1 min and "Limitations" 74 words in 0.5 min (148 wpm), against about 115 elsewhere. Fix: drop the script's "In the mental health study, for example, ..." sentence (it repeats the slide bullet), and give Limitations 0.75 min, taking 0.25 from "Larger grids" (129 words in 1 min, also fast).

**m2. Near-ties on "Larger grids" (partly fixed).** "The better design matches the 2×4 grid in **16 of 16** scenarios on 3×4" counts ψ = 0.8, ρ = 0.01, Control Clusters Only (0.472 vs 0.480, N = 10) as a win. Fix: add "(near-tied at ψ = 0.8)" as the decision table now does.

**m3. Limitations, "No intra-cluster correlation" (carried).** Every subject in a cluster shares the random lag term ρ(Wy), so outcomes within a cluster are correlated when ρ ≠ 0; there is no cluster random effect. Fix: "No cluster-level random effect" on the slide and in the script.

**m4. ρ backup, "about 40 to 60 links" (carried).** Holds for 2×4. The 3×3 centre and 3×4 interior clusters have four rook neighbors (about 80 links, about 0.8 on the fitted scale). Fix: add "on the 2×4 grid, up to about 80 (0.8) for interior clusters of the larger grids".

**m5. fig_mse_2x4_controls.png and fig_mse_2x4_both.png compare a mean with a median (carried).** Cards quote SRS means (0.304, 0.173). The boxes show medians (about 0.24 at ρ = 0.01, 0.25 at ρ = 0) and the diamonds BSS means. Fix: add an SRS-mean marker (open grey circle, "SRS mean" in the legend).

**m6. Neighbors slide vs figure (carried).** Panel (c) of RookQueen2ndContiguityColored.png is a distance radius, not second-order contiguity, and the alt text says "second-order contiguity". Fix: "(a) Rook ... (b) queen ..." in the first bullet, "(c) Distance-based neighbors lie within a fixed radius" in the second, and drop "higher-order neighbors are neighbors of neighbors" or say it is not pictured. Correct the alt text.

**m7. SampledPoints.png legend (carried).** It numbers district 1 at top right and 2 at bottom left, which conflicts with the row-major numbering on the weights backup. Fix: crop the legend; the note says only "colored by cluster".

**m8. "NC probation districts" (carried).** "30 judicial districts". Both maps show 31 areas (19A, 19B). Fix: "30 judicial districts (district 19 split into 19A and 19B)".

**m9. Colons in slide text (partly fixed).** Text: "Why controls only: with a consistent treatment effect" (Recommendation), "(2×4: 0.165 vs 0.304)" and two more in the decision table, "(2×4: MSE 2.39 ...; 3×4: 2.57)" and "(3×3: ...; 3×4: ...)" in Downside risk, "(50; 10 on 3×4): place" in the algorithm. Fix: "Controls only, since with a consistent treatment effect ..."; replace "2×4:" with "on 2×4,".

**m10. Model glossary, α (carried).** "Baseline effect (0.2 in the simulation)". α is an intercept. Fix: "Intercept, the baseline mean outcome (0.2)".

**m11. Model glossary, σ.** "Independent normal error (subject-level SD σ = 0.1)" sits beside an equation for the cluster-level outcome y_i, whose error SD would be about σ/√20. Fix (within the cluster-level wording decision): "Independent normal error (σ = 0.1 for each subject's outcome)".

**m12. Decision table, top-left cell.** "Lower MSE in 12 of 12 scenarios across grids (near ties at ψ = 0.8) (2×4: 0.165 vs 0.304)". Two parentheticals in a row. Fix: "Lower MSE in 12 of 12 scenarios across grids, near-tied at ψ = 0.8 (0.165 vs 0.304 on 2×4)".

**m13. Backup 3×3/3×4 figures use different conventions.** ρ = 0 is on the left (the main figure has ρ = 0.01 left), rows read "Intervention & Control" (deck: "Control & Intervention Clusters"), the legend uses capital Ψ, and the x-axis reads "Block / Random". Fix: none needed if the note says so; add "ρ = 0 on the left here; 'Intervention & Control' is the deck's Control & Intervention Clusters" to the speaker note.

**m14. NC_MCO_Map.png legend truncated.** "Partners Health Manageme...". Fix: re-export wider, or name it in the bullet ("Alliance, Partners, Trillium and Vaya").

**m15. DOI (carried).** 10.1186/s12874-026-02968-0 is in references.bib and the memo, but nothing in the repo shows it came from the publisher. Confirm it on the BMC page before the talk.
