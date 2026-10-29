# Part I deck review (independent biostatistician), round 3

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.qmd, as knitted in presentation/spillover_talk.revealjs.md, with every PNG it references in presentation/figures/ (all 13 opened)
Sources: review/manuscript/Revisions_V2c_corrected.tex; review/tables/beta_mse_table_{2x4,3x3,3x4}_corrected.csv; per-allocation CSVs in results/beta_mse, results/rho_mse, results/psi_mse and review/results (read through presentation/tools/load_results.R, read-only); code/SpatialSim_NC_DOC.Rmd (rook list, lines 310–319); review/notes/2026-10-28-rho-scale-note.md; presentation/tools/make_figures.R; rounds 1 and 2 in review/deck-reviews/. Check scripts were run from the session scratchpad and write nothing to the repo.

**Verdict:** Blocking issues: 1 (2 major, 14 minor)

Recomputed and confirmed: 57% (57.1 at every ψ on 2×4; 57.4 if the BSS allocations are left out of the denominator), worst tenth 1.20 (7 allocations), worst 2.39 (combo {1,2,3,4}, a 2×2 block under the code's rook list, bias +1.53, identified), BSS 0.144/0.186, 40–60% (40.2–60.3), at most 29% (28.9) with spillover into both arms, 2 of 22 and 151 of 436 identified draws beaten, 3×3 0.147 vs 0.84, 3×4 0.133 vs 0.77 and worst 2.57 (identified), 48/57/488 and 2/1/2 collinear, 16/16/126 with z ≡ 1, the six 3×3 BSS allocations are exactly the six independent 4-sets and 3×4 has exactly two, 3×3 BSS 0.246 vs ≤ 0.0008, BSS bias −0.50 to −0.80 and −0.35 to −0.72, median ρ̂ 0.50, BSS ρ̂ 0.89/0.85. Timing is 20.0 min on content slides (20.5 with title and outline, 20.75 with the thank-you). No script bullet opens with "Because", every bullet is one sentence, and the script has no colons.

---

## Round-2 BLOCKING and MAJOR items

| r2 item | Status | Note |
|---|---|---|
| B1 Limitations "correctly specified" | Fixed | Slide and script now say "with the spillover term included", which no longer asserts the weights are right. |
| B2 NC example map and the check it recommends | Settled by author decision, wording partly fixed | Map kept and described as a suggested best-effort allocation; the "east" sentence is gone. The check sentence was removed rather than computed, and the replacement bullet points the improvement the wrong way (see M1). |
| M1 Who BSS beats among identified draws | Fixed | New backup "Downside risk" carries 2 of 22, 151 of 436, the 2.39 draw and a prepared answer. |
| M2 "Most random allocations do better" scope | Fixed | Now "On this grid, ...". |
| M3 Collinearity condition stated generally | Fixed | Title and bullet scoped to Control Clusters Only; the z ≡ 1 aliasing line is in the reduced-form backup. |
| M4 BSS SD on 2×4/3×4 is Monte Carlo noise | Fixed | Card says "effectively one allocation (two mirror-image checkerboards)", "38 of 48" dropped, table note gives 0.144 vs 0.186 as a Monte Carlo guide. |
| M5 Inference backup "bounding the worst case" | Fixed | Now scoped by spillover type. |

Round-2 minors, briefly: m1 fixed in the algorithm backup but the Model slide now disagrees with it (m4 below); m2 the NC slide is fixed, Recommendation and Limitations are still fast (m1 below); m3 the x-axis collision is fixed but the facet strips are now clipped (B1 below); m4–m7 and m9 fixed; m8 mostly fixed (m11 below); m10a fixed by caption and script; m10b the DOI now appears in review/thesis-chapter/references.bib and the correction memo, but still confirm it with the publisher.

---

## BLOCKING

**B1. "The same pattern holds on the 3×3 and 3×4 grids", figure fig_mse_larger_grids.png**
Text (facet strips as rendered): "3×3 grid · ρ = 0.0", "3×3 grid · ρ = 0", "3×4 grid · ρ = 0.0", "3×4 grid · ρ = 0".
Problem: The strip labels are clipped. make_figures.R (line 85) writes "ρ = 0.01", but the panel width cuts the last digit, so the main results slide shows two panels labelled ρ = 0.0 beside two labelled ρ = 0 with very different results. The first thing a committee member reads off this slide is a wrong label for the factor the talk is about.
Fix: Shorten or wrap the strips, for example `sprintf("%s\nρ = %s", grid, ...)` (two-line strips), or `theme(strip.text = element_text(size = rel(0.8)))`, or put ρ on the column strip only once with `facet_grid(row ~ grid + rho)`. Re-render and open the PNG to confirm "0.01" is visible.

---

## MAJOR

**M1. "Recommendation" slide, third bullet, and the example map**
Text: "This is a suggested, best-effort allocation for an irregular map; further restrictions could guarantee that no intervention district borders another". Script: "...further restrictions could guarantee that no intervention district borders another." and "a pure checkerboard still underestimates the effect in this setting, by about a third".
Problem: The map stays (author decision), but this bullet now names the wrong improvement. The deck's own result is that bias comes from every control being exposed, not from intervention districts touching. More separation of intervention districts spreads treatment across the map and tends to expose more controls; the strict checkerboard is exactly the collinear case. By eye, every control district on the map borders at least one intervention district (2, 5, 6, 8, 9, 11, 14, 16, 17, 18, 20, 22, 25, 26, 28, 29 each touch red), so under the spillover-into-controls assumption the slide relies on, this allocation would estimate β − ψ, like the checkerboard. The Take-home card then tells investigators to "look at which controls the chosen allocation exposes". A statistician will ask whether the example passes that check, and the honest answer is no. I verified adjacency only by eye (the shapes come from `tigris`).
Fix: Replace the bullet with "A suggested, best-effort allocation; before the trial, check which controls it exposes and adjust so that a few control districts border no intervention district". In the script, replace the "further restrictions" sentence with "Before the trial we would check the map so that a few control districts border no intervention district, which keeps the effect estimable." Prepare a one-line answer: "As drawn, every control borders an intervention district, so we would swap a pair before using it." Ideally run `spdep::poly2nb()` on the NCTrtMaps.Rmd polygons once so the answer is computed rather than by eye.

**M2. Backup "How large is ρ = 0.01?", fifth bullet**
Text: "The long tail of random draws exists only with dependence. The worst 2×4 random draw has MSE 2.39 at ρ = 0.01 and 0.27 at ρ = 0".
Problem: The 0.27 is not a random draw. At ρ = 0 (Control Clusters Only, ψ = 0.5) the maximum, 0.268, is combo 28, one of the two BSS checkerboards; the worst non-BSS allocation is 0.259. The allocation that is worst at ρ = 0.01 ({1,2,3,4}, MSE 2.39, bias +1.53) has MSE 0.007 at ρ = 0. Pulled up in answer to a question, the bullet misidentifies the allocation, and the correct fact, one allocation going from 0.007 to 2.39 when ρ moves from 0 to 0.01, is the one a statistician will want to probe. This is wording within the author's decision on the weights; it does not reopen it.
Fix: "The long tail appears only at ρ = 0.01. The allocation with the worst error there (half the grid treated as one block, MSE 2.39) has MSE 0.007 at ρ = 0, and at ρ = 0 the largest error of any allocation is a checkerboard's (0.27)." Add to the notes that this is the allocation the rho note's follow-up check would examine first.

---

## MINOR

**m1. Timing on two slides.** "Recommendation" runs at 147 wpm (147 words in 1 min) and "Limitations" at 148 wpm (74 words in 0.5 min), against a deck average of about 112. Fix: the M1 rewrite can drop the "In the mental health study, for example, ..." sentence (it repeats the slide); give Limitations 0.75 min and take 0.25 from "Larger grids" (90 words in 1 min).

**m2. Near-ties counted as wins.** "Choosing a design": "Lower MSE in 12 of 12 scenarios across grids"; "Larger grids": "agrees with the 2×4 direction in 16 of 16". At ψ = 0.8 the margins are 0.569 vs 0.573 (2×4, N = 50) and 0.472 vs 0.480 (3×4, N = 10), well inside the Monte Carlo error the deck itself cites (0.144 vs 0.186 for mirror-image checkerboards). Fix: "Lower MSE in 12 of 12 scenarios, near-tied at ψ = 0.8".

**m3. Recommendation script, "by about a third".** Text: "a pure checkerboard still underestimates the effect in this setting, by about a third". That holds at ψ = 0.5 (−0.35); at ψ = 0.8 it is −0.72. Fix: "by about a third at ψ = 0.5, and more as spillover grows".

**m4. Model slide vs algorithm backup.** Model: "Subjects are sampled within clusters (160 per replication on the 2×4 grid, about 20 per cluster)". Algorithm step 3: "place 160, 180 or 240 subjects at random over the grid". The second is what the code does. Fix: Model slide "Subjects are placed at random over the grid (160 per replication on the 2×4 grid, about 20 per cluster)".

**m5. Limitations, "No intra-cluster correlation".** In the generator, every subject in a cluster shares the same lag term ρ(Wy), so at ρ = 0.01 subjects within a cluster are correlated; there is just no cluster random effect. Fix: "No cluster-level random effect" on the slide and in the script.

**m6. ρ backup, first and third bullets.** "about 40 to 60 links" and "0.01 × (40 to 60) ≈ 0.4 to 0.6" hold for the 2×4 grid. The 3×3 centre and 3×4 interior clusters have four rook neighbors (about 80 links, about 0.8 on the fitted scale). Fix: add "on the 2×4 grid; up to about 80 links (0.8) for interior clusters of the larger grids".

**m7. fig_mse_2x4_controls.png vs the cards.** The cards quote SRS means (0.304, 0.173), but the figure shows SRS boxes whose medians differ from those means (about 0.24 at ρ = 0.01; about 0.25 at ρ = 0, where Q1 is near 0 so the box reads as a bar from zero). Diamonds are BSS means, so the figure compares a mean with a median. Fix: add an SRS-mean marker (for example a grey open circle, "SRS mean" in the legend) to this figure and fig_mse_2x4_both.png.

**m8. "Neighbors" slide figure.** RookQueen2ndContiguityColored.png shows (a) rook, (b) queen and (c) a distance radius, not second-order contiguity, which the alt text and the bullet "higher-order neighbors are neighbors of neighbors" imply. Panels are not mapped to words on the slide. Fix: first bullet "(a) Rook neighbors share an edge, (b) queen neighbors an edge or a corner"; second bullet "(c) Distance-based neighbors lie within a fixed radius"; correct the alt text.

**m9. SampledPoints.png legend.** The legend numbers districts 1 (top right), 2 (bottom left), 3 (bottom right) ..., which contradicts the code's rook list, where 1 neighbors 2 and 3. Nobody will ask unless they read the legend. Fix: drop the legend (the slide only needs "colored by cluster").

**m10. "NC probation districts", first bullet.** Text: "30 judicial districts". Both maps show 31 areas (19A and 19B). Fix: "30 judicial districts (district 19 split into 19A and 19B)".

**m11. Colons left in slide text.** "Why controls only: with a consistent treatment effect"; Downside risk "(2×4: ...", "(3×3: ...", "(2×4: MSE 2.39"; Choosing a design "(2×4: 0.165 vs 0.304)" and two more; algorithm step 3 "(50; 10 on 3×4): place". Fix the prose one at least: "Controls only, since with a consistent treatment effect treated districts gain nothing more from other treated districts". The parenthetical "2×4:" labels can become "on 2×4,".

**m12. "Simulation design", N.** The formulas use N, but the slide never says what N is (the table column is "Replications"). Fix: "β estimated in each of N replications (table)". If asked: the saved `variance` column uses the N − 1 divisor, so bias² + variance in the CSVs differs from MSE by variance/N; the slide's 1/N definition is the one that makes the identity exact.

**m13. Model glossary, α.** Text: "Baseline effect (0.2 in the simulation)". α is an intercept, not an effect. Fix: "Intercept, the baseline mean outcome (0.2)".

**m14. DOI (carried over).** 10.1186/s12874-026-02968-0 now appears in review/thesis-chapter/references.bib and the correction memo, but no file shows it came from the publisher. Confirm it on the BMC page before the talk.
