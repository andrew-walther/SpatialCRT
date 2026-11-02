# Part I deck review (independent biostatistician), round 1

Date: 2026-10-28
Deck reviewed: presentation/spillover_talk.qmd, as knitted in presentation/spillover_talk.revealjs.md (Draft 1), with the figures in presentation/figures/
Sources checked against: review/manuscript/Revisions_V2c_corrected.tex; review/tables/beta_mse_table_{2x4,3x3,3x4}_corrected.csv; per-allocation CSVs in results/beta_mse, results/psi_mse, results/rho_mse and review/results (loaded through presentation/tools/load_results.R); code/SpatialSim_NC_DOC.Rmd (generation, fit, compute_metrics); review/notes/2026-10-28-rho-scale-note.md. The check scripts live in the session scratchpad and are read-only.

**Verdict:** Blocking issues: 3 (5 major, 14 minor)

Every number on the main slides and in the three backup tables matches the corrected tables and per-allocation files: 57%, 2.39, 0.19 (BSS max 0.186), 40–60% (40.2–60.3), 29% (28.9), 0.165/0.304, 0.173/0.262, SD 0.41/0.03, 0.095/0.165, 0.008/0.262, bias −0.80 to −0.35, 24 of 24, 12 of 12, 16 of 16, 12 of 16, 38 of 48, median ρ̂ 0.50, and MSE rising with ψ in all 12 grid × type × ρ cells. Timing markers add up to 20.0 minutes. No script sentence opens with "Because", and each bullet is one sentence. The problems are in interpretation and wording.

---

## BLOCKING

**B1. Slide "Why: collinearity" (also "Limitations", "Two designs")**
Text: "Every control in a checkerboard borders an intervention cluster, so $z_i = 1 - x_i$ for every cluster" / "Random allocations that leave some controls unexposed keep the information needed to separate β from ψ" / Limitation "Treatment and spillover can be exactly collinear under BSS".
Problem: Under Control Clusters Only, $z = 1 - x$ holds whenever every control borders a treated cluster, whether or not treated clusters touch each other. Enumerating the allocations gives this many with exact collinearity under Control Clusters Only: 2×4 **48 of 70**, 3×3 **57 of 126**, 3×4 **488 of 924**. Under Control & Intervention the only collinear allocations are the checkerboards (2, 1 and 2). The "Random" 2×4 example in GridsCombined.png (treated 1, 5, 6, 8) is itself collinear under Control Clusters Only. The data show it too. At ρ = 0 the median random 2×4 allocation has MSE 0.249 ≈ ψ², the same as BSS. A committee member who works one example will catch this, and as written the slide implies collinearity is special to BSS.
Fix: State the general condition, "Under control-only spillover, any allocation in which every control borders a treated cluster has $z = 1 - x$. That is 69% of random 2×4 allocations, and every checkerboard. Under spillover into both arms, only checkerboards." Then give the mechanism for the both-arms result plainly: there, BSS holds the only collinear allocations, which is why SRS wins in all 24. Reword the limitation to "Treatment and spillover can be exactly collinear when every control is exposed, which includes every checkerboard."

**B2. Slide "Where BSS falls" (title and script)**
Text: Title "with spillover into controls, block stratified sampling beats most random draws". Card: "Across all grids and ψ in that setting, BSS beats 40–60% of random draws".
Problem: The title drops the ρ = 0.01 condition. The second row of the slide's own figure is Control Clusters Only with ρ = 0, where BSS beats 1.4% of 2×4 draws at ψ = 0.5 and 1–47% across ψ. "Most" also overstates the claim, since 3 of the 12 values are below 50% (3×4: 47%, 44%, 40%). The script covers rows 1, 3 and 4 and says nothing about row 2.
Fix: Retitle to something like "A trial runs one allocation: with spillover into controls and spatially correlated outcomes, block stratified sampling beats about half of random draws and avoids the tail". Add one script sentence on row 2, for example "Without spatial dependence, most random draws do better, and the cluster of draws at the same error as block stratified sampling are allocations that also expose every control."

**B3. Backup "Inference", plus the "Model" script**
Text (backup): "Maximum likelihood conditions on the realized allocation and gives $\hat\beta$ a sampling distribution from the outcome error, whether the allocation came from SRS or BSS". Script (Model): "We fit it by maximum likelihood with the true spillover term included, so this is a best case for estimation".
Problem: For the BSS allocations on 2×4 and 3×4, and for every collinear random allocation, $[\mathbf 1, x, z]$ is rank-deficient. `lagsarlm` drops the aliased spillover column. ψ̂ is NA for combos 28 and 43 in results/psi_mse/psi_mse_results_TrtNoSpill_2x4_A02_B1_P05_R00.csv, and $\hat\beta$ is then an estimate of β − ψ. So the BSS fit does not include the true spillover term, and model-based inference for β under BSS on these grids is not available. The collinearity slide contradicts the inference backup, and a statistician will ask "how did you get a β̂ at all from a singular design?"
Fix: In the inference backup, add "provided the allocation leaves β and ψ identified (not every control exposed). For the 2×4 and 3×4 checkerboards it does not, so `lagsarlm` drops the spillover term and $\hat\beta$ estimates β − ψ." In the Model script, change the "best case" sentence to "...with the true spillover term included wherever the allocation allows it to be estimated". Add one line to the collinearity slide or a backup saying what the software did, for example "lagsarlm drops the aliased z; ψ̂ is not estimated".

---

## MAJOR

**M1. Slide "Why: collinearity"**
Text: "$\hat\beta$ absorbs the spillover, $\hat\beta \approx \beta - \psi$, so 2×4 BSS bias runs from **-0.80** to **-0.35**".
Problem: β̂ ≈ β − ψ holds at ρ = 0, where bias is −0.496, −0.597, −0.697 and −0.798. At ρ = 0.01 the bias is −0.349 to −0.716, and −0.35 at ψ = 0.5 is not −ψ. The "so" links a range to a formula that only half the range satisfies. The manuscript itself says "bias near −ψ when ρ = 0".
Fix: "At ρ = 0 the bias is −ψ (−0.50 to −0.80). With spatial dependence the lag term absorbs part of it (−0.35 to −0.72)."

**M2. Slides "Choosing a design", "Take-home", "Recommendation"**
Text: "Spillover into control and intervention clusters → Simple random sampling" and "check any allocation for treatment–exposure collinearity before the trial".
Problem: The 2×4 and 3×4 BSS allocations always fail that check, so the deck recommends BSS and a check that BSS fails. The data also put the problem in full exposure of controls, not in spatial separation. On 3×3, the five BSS allocations that leave a control unexposed have MSE ≈ 0.0005 at ρ = 0 under both spillover types, against an SRS mean of 0.002–0.005. Only the full checkerboard (combo 77) has 0.246. Expect the question "so is the problem BSS or the checkerboard?"
Fix: Make the check concrete and tie it to the recommendation, for example "Choose a spatially separated allocation that still leaves some controls with no treated neighbor (z ≠ 1 − x)." Note on the Recommendation slide that irregular real maps make a perfect checkerboard unlikely, which is what keeps the best-effort NC allocation identified. That point is already half there in "perfect separation is not always possible".

**M3. Slides "Two designs" and "Larger grids" (3×3 relaxation)**
Text: "On the 3×3 grid ... only intervention clusters must be separated, and 6 allocations qualify" / "its 6 BSS allocations are not all checkerboards".
Problem: The six allocations are correct: {1,3,5,7}, {1,3,5,9}, {1,3,7,9}, {1,5,7,9}, {2,4,6,8} and {3,5,7,9}, and the indices match the frozen code. But exactly one, {2,4,6,8}, the strict BSS allocation, is collinear, and it is the 3×3 example drawn in GridsCombined (c). It alone drives the 3×3 BSS mean and SD (ρ = 0, ψ = 0.5: 0.246 vs ≈0.0005 for the other five). That is also why 3×3 BSS SD exceeds SRS SD in all 8 Control & Intervention scenarios, which accounts for 8 of the 10 exceptions in "38 of 48".
Fix: Say "One of the six (the full checkerboard, shown) is collinear. The other five each leave one control unexposed and estimate β almost perfectly." Optionally say in the script that the 10 SD exceptions are 3×3 Control & Intervention (8) and 3×4 Control & Intervention at ρ = 0 at ψ = 0.5 and 0.6 (2).

**M4. Model slide glossary and the "Scale of rho" backup (wording only, not the decision to present as published)**
Text: "$w_{ij}$ Row-normalized weight ... ρ Spatial dependence coefficient (0 or 0.01)".
Problem: Read together, these two lines describe a generator the code did not use. The backup then contradicts the main slide. Two facts belong in that backup and are not there yet: (a) at ρ = 0.01 the BSS fits give ρ̂ ≈ 0.89 and 0.85, against a median of 0.50 across allocations, so in the non-identified BSS fit the lag term soaks up part of the omitted spillover, which is why BSS bias shrinks from −0.50 to −0.35; (b) the "worst random draw, 2.39" (combo 1, bias +1.53) exists only at ρ = 0.01. At ρ = 0 the maximum is 0.268. The tail-risk argument therefore rests on the setting the rho note flags.
Fix: Add a small marker on the Model slide's ρ entry, for example "(0 or 0.01 on the generating scale; see backup)". Add bullets (a) and (b) to the rho backup so the speaker is not surprised when asked why BSS improves with spatial dependence.

**M5. Slides "Control clusters only", "Choosing a design"**
Text: "lower for every ψ"; "Lower MSE in 12 of 12 scenarios across grids".
Problem: At ψ = 0.8 the gaps are 0.569 vs 0.573 (2×4) and 0.472 vs 0.480 (3×4). These are within Monte Carlo error with N = 50 or 10 replications, and the BSS rank on 3×4 falls from 55% to 40% as ψ rises. "Every ψ" invites a question about Monte Carlo error that the deck cannot answer, since it reports no Monte Carlo SE.
Fix: "Lower at every ψ, with the gap closing by ψ = 0.8 (0.569 vs 0.573)". Or keep "12 of 12" and add "narrowing as ψ grows" in the script.

---

## MINOR

**m1. "Simulation design".** Text: "$\text{MSE} = \frac{1}{N}\sum_r(\hat\beta_r-\beta)^2 = \text{bias}^2 + \text{variance}$". Bias and variance are not defined, and `compute_metrics` uses `var()` (N − 1), so the identity holds only with a 1/N variance. Fix: define them on the slide, "bias = $\bar{\hat\beta} - \beta$, variance = $\frac1N\sum_r(\hat\beta_r-\bar{\hat\beta})^2$". Optionally note that the tables report the N − 1 variance.

**m2. "Simulation design" / "Model".** Text: "20 simulated subjects per cluster" and the limitation "equal cluster sizes". The code draws `points_per_iteration <- 160` points uniformly over the whole grid and assigns them to districts, so cluster sizes are random with an expected size of 20. The manuscript also says counts are "not necessarily identical". Fix: "about 20 per cluster (160 per grid)", and change the limitation to "equal-area clusters".

**m3. "Control clusters only" script.** Text: "In both panels, the random allocations spread far more, with a standard deviation of MSE of 0.41 against 0.03". Those values are for the ρ = 0.01 panel. At ρ = 0 they are 0.117 vs 0.009. Fix: "With spatial dependence, ... 0.41 against 0.03."

**m4. fig_mse_2x4_controls.png, right panel.** The SRS boxes read as bars, because Q1 is near 0 and the median and Q3 sit at ψ² (the collinear majority), and the two panels have different y scales. Fix: add one script sentence ("without spatial dependence, the typical random draw is as bad as BSS because it also exposes every control"), or show the ρ = 0 panel as dots as in fig_spread.

**m5. fig_spread_2x4.png.** Each row has its own square-root axis, and the tick formats are mixed (0.0 vs 0.00). Fix: use the same number of decimals in every row and say "each row has its own scale" in the note under the figure.

**m6. fig_bss_rank_all_grids.png (backup).** Identical values overplot: 2×4 Control Clusters Only ρ = 0.01 is 57.1% at all four ψ, 3×3 ρ = 0 is 45.2%, and Control & Intervention ρ = 0 overlaps too, so "one point per ψ" shows one point. Fix: add a small vertical dodge, or note "overlapping points share a value".

**m7. Colour and numbering consistency.** Controls are pale grey with amber exposure circles on "Two spillover types", Carolina blue for exposed and grey for unexposed in fig_collinearity_grid, and blue (with red for intervention) in GridsCombined. GridsCombined's 2×4 uses the column-major numbering (7 5 3 1 / 8 6 4 2), while the weights backup uses row-major ("Cluster 2 has rook neighbors 1, 3 and 6"), which is wrong under GridsCombined's numbering. Fix: one palette for intervention, exposed control and unexposed control across schematics, and drop the cell numbers from GridsCombined or match them to the backup.

**m8. GridsCombined.png.** It carries in-figure "(a)–(f)" and "Block/Random" labels that the slide does not use, against the clean-figure rule. Fix: crop the labels or regenerate without them, since the slide note already says top = block stratified and bottom = random.

**m9. "The question".** Text: SRS is "Any balanced allocation of clusters". The 3×3 grid uses 4 vs 5. Fix: "balanced (4 vs 5 on the 3×3 grid)". The script for this slide runs at about 194 wpm for 0.5 min, so trim one sentence or give it 0.75 min and take the time from "Control & intervention clusters" (about 62 wpm).

**m10. "Control & intervention clusters".** Text: "The 2×4 BSS results are identical under both spillover types". This holds for every BSS allocation on every grid, since no treated cluster has a treated neighbor and the two z definitions coincide. The 3×3 tables show it too. Fix: "BSS results are identical under both spillover types on every grid".

**m11. Backup "Reduced form and likelihood".** n is undefined, and the fit is at the subject level (n = 160 with a subject-level W). Fix: "n = number of observations".

**m12. Backup "Estimands".** "$E[Y(0,\mathbf 1)]$" notation is undefined. "Following Halloran and colleagues, a CRT's intervention contrast is a total effect" does not appear in the corrected manuscript, and is loose: with full coverage in treated clusters, the Halloran–Struchiner overall and total effects coincide. Fix: define the notation in a few words ("focal cluster at control; all others treated vs none") and either cite precisely (Halloran & Struchiner 1991; Hudgens & Halloran 2008, which the paper cites) or drop the last bullet.

**m13. "Choosing a design".** Text: "Simple random sampling, or a buffer-zone design". Buffer zones were not evaluated here, and the buffer backup says these grids are too small for them. Fix: add "(not evaluated here)".

**m14. "Limitations" citation and style.** The DOI 10.1186/s12874-026-02968-0 is not in any file in the repo, so check it against the publisher before the talk. Colons in slide text: "Exposure is **binary**: a cluster..." and the "ρ = 0.01:" / "ρ = 0:" card labels. Only one script colon remains, "Answer:" in the rho backup notes. Fix: "Exposure is binary, so a cluster..." and drop the colon after "Answer".
