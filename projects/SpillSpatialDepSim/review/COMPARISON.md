# Published vs. corrected results (2026-10-22)

What changes when the two errors in `README.md` are fixed. Everything here
is built by `code/build_tables_figures.R` from
`tables/beta_mse_table_<grid>_<published|corrected>.csv`. Full tables in the
paper's units, with published values in brackets wherever they changed, are
in `tables/corrected_tables_paper_units.md`. Figures are in `figures/`.

## Checks that passed

- **The rerun code reproduces the published CSVs.** In bug mode (`1:8`
  loop) it matches the published Control Only CSVs to 3e-9 (3×3, 126
  allocations) and 1.2e-8 (3×4, 924 allocations).
- **The table builder reproduces the paper.** In "published" mode it gives
  all three published tables, and the 2×4 figure matches the published one.
- **The corrected Control Only reruns are complete.** All 16 have full row
  counts (126 or 924) and no NA estimates.
- **The model-fit warnings are expected.** lagsarlm warns "Aliased
  variables found: z" for allocations where z = 1 − x. The original
  scripts raised the same warning.

## Results by grid (lower mean MSE; 8 scenarios per spillover type)

| Grid | Control Only, ρ = 0.01 | Control Only, ρ = 0 | Intervention & Control (both ρ) | BSS SD of MSE lower than SRS |
|---|---|---|---|---|
| 2×4 (unchanged) | BSS 4/4 | SRS 4/4 | SRS 8/8 | 16/16 |
| 3×3 published | BSS 4/4 | BSS 4/4 | **BSS 8/8** | 16/16 |
| 3×3 corrected | BSS 4/4 | BSS 4/4 | **SRS 8/8** | 8/16 (Control Only only) |
| 3×4 published | BSS 4/4 | BSS 4/4 | SRS at ρ = 0.01, **BSS at ρ = 0** | 12/16 |
| 3×4 corrected | BSS 4/4 | **SRS 4/4** | **SRS 8/8** | 14/16 |

Example values (ρ, ψ = 0.5; MSE, SRS vs. BSS):

| | 2×4 | 3×3 published → corrected | 3×4 published → corrected |
|---|---|---|---|
| Control Only, ρ = 0.01 | 0.304 vs 0.165 | 0.234 vs 0.041 → **0.274 vs 0.147** | 0.120 vs 0.035 → **0.244 vs 0.133** |
| Control Only, ρ = 0 | 0.173 vs 0.262 | 0.092 vs 0.001 → **0.114 vs 0.042** | 0.034 vs 0.0004 → **0.133 vs 0.226** |
| Int. & Control, ρ = 0.01 | 0.095 vs 0.165 | 0.103 vs 0.039 → **0.103 vs 0.147** | 0.112 vs 0.223 → **0.112 vs 0.133** |
| Int. & Control, ρ = 0 | 0.008 vs 0.262 | 0.002 vs 0.0005 → **0.002 vs 0.042** | 0.001 vs 0.0003 → **0.001 vs 0.226** |

**The main finding.** With both errors fixed, **3×4 shows the same pattern
as 2×4 in every one of its 16 scenarios**, and 3×3 matches 2×4 in 12 of
16. In all three grids:
- Intervention & Control: SRS has lower mean MSE in every scenario.
- Control Only with spatial dependence (ρ = 0.01): BSS has lower mean MSE
  in every scenario.
- Control Only, ρ = 0: SRS wins in 2×4 and 3×4, BSS wins in 3×3. This is
  the only place where the grids disagree.

The published 3×3 and 3×4 tables showed BSS winning almost everywhere.
That apparent result came from the wrong allocations.

**Why the 3×4 numbers moved so much:**
- Under Control Only, every allocation in which all control clusters
  border a treated cluster has z = 1 − x exactly. In those allocations β̂
  absorbs the spillover, so bias ≈ −ψ.
- Such allocations are common: 57 of 126 in 3×3 (including BSS combo 77)
  and 488 of 924 in 3×4 (including both checkerboards).
- The `1:8` bug hid most of this collinearity in the 3×4 SRS allocations.
  Fixing it roughly doubles SRS Control Only MSE (0.120 → 0.244 at
  ρ = 0.01, ψ = 0.5) and changes the mean SRS bias from +0.03 to about
  −0.19 to −0.42.

## Which of the paper's conclusions hold

| Conclusion (Chapter 2 / paper) | 2×4 | Corrected 3×3 | Corrected 3×4 | Verdict |
|---|---|---|---|---|
| Choice between SRS and BSS is a risk-management trade-off: BSS has stable error, SRS can do better on average but varies widely | ✓ | ✓ for Control Only; under Int. & Control BSS is neither lower on average nor more stable | ✓ (BSS SD lower in 14/16) | **Holds**, with the 3×3 Int. & Control exception |
| BSS is preferable when spillover flows into control clusters, especially with spatial dependence | ✓ (ρ = 0.01) | ✓ (both ρ) | ✓ (ρ = 0.01) | **Holds, and now consistent across all three grids** |
| SRS often has lower MSE when spillover reaches both arms | ✓ | ✓ (8/8) | ✓ (8/8) | **Holds, and stronger** (published 3×3 said the opposite) |
| Spatial dependence (ρ = 0.01) raises MSE for both designs | SRS ✓; BSS ✗ (0.262 → 0.165) | ✓ both | SRS ✓ (e.g. 0.133 → 0.244); BSS ✗ (0.226 → 0.133) | **Holds for SRS everywhere; for BSS only in 3×3.** Already false for 2×4 BSS as published. Why BSS error falls with ρ for the checkerboards has not been checked |
| MSE rises with ψ | ✓ | ✓ | ✓ | **Holds** |
| SRS has lower bias, BSS lower variance (Discussion) | ✓ | ✗ (BSS lower \|bias\| in 12/16) | ✓ SRS lower \|bias\| 16/16 | **Holds for 2×4 and 3×4, not 3×3** |
| BSS "consistently exhibits minimal MSE" (Discussion, M2.5.1) | ✗ (4/16) | ✗ (8/16) | ✗ (4/16) | **Does not hold** in any grid, published or corrected |
| 7 allocations are each optimal for some scenario; BSS never the single best | 7, BSS never | 5 (was 7), BSS never | 6 (was 7), BSS never | **"BSS never best" holds; "7" holds only for 2×4** |
| 3×3: BSS better in MSE, SD and bias under both spillover types | n/a | ✗ under Int. & Control | n/a | **Does not hold** |
| 3×4: BSS "dramatically outperforms" SRS under Control Only at all ρ | n/a | n/a | ✗ at ρ = 0 (0.226 vs 0.133) | **Does not hold at ρ = 0** |
| "Checkerboard paradox" (BSS backfires only under Int. & Control with ρ > 0, through diagonal spillover) | n/a | n/a | ✗ | **Does not hold.** The real effect is the −ψ collinearity bias, under both spillover types and at both ρ |
| Recommendation for NC probation: BSS when spillover flows mainly into controls | ✓ | ✓ | ✓ (with ρ > 0) | **Holds** |

Net: the paper's general conclusions and its applied recommendation hold.
They are more consistent across grids after the correction than before.
What fails are the grid-specific statements for 3×3 and 3×4, the
"consistently minimal MSE" sentence (already false for 2×4), the count of
7 optimal allocations for 3×3 and 3×4, and the checkerboard paradox.

## Next steps

Nothing below has been done. Items 3–6 touch Chapter 2, the paper or the
slides, and need Andrew's explicit approval.

1. **Advisor discussion.** Share this file and `README.md` with Dr. Lin.
   Key message: the general conclusions and the recommendation hold; the
   3×3 and 3×4 tables, figures and grid-specific text were wrong. Decide
   whether a journal correction is warranted. A BMC correction notice
   would replace two tables and two figures and revise a few paragraphs.
2. **Record decisions in bios-dissertation `DECISIONS.md`.** These are the
   two decisions from the side chat (retire "checkerboard paradox"; 2×4
   primary in the prelim talk), plus whatever is decided in item 1.
3. **Prelim slides.**
   - Lead with 2×4, and use the collinearity-bias slide in place of slot
     I.9 (already noted in `slides/prelim/OUTLINE.md`).
   - The outline currently says to use the published 3×3 and 3×4 numbers
     in backup unless updated ones are provided. These corrected numbers
     are now available. Swapping them in would let the backup slides say
     "the same pattern holds on 3×4". That is a stronger message than the
     published tables give, and it is defensible.
4. **Chapter 2 (when approved).**
   - Replace the 3×3 and 3×4 tables and figures with the corrected ones.
   - Rewrite the 3×3 and 3×4 results paragraphs.
   - Remove the checkerboard paradox (Results, Discussion,
     Recommendations), keeping the Limitations account of collinearity.
   - Fix the "consistently minimal MSE" and "7 allocations" sentences.
   - Add a one-line note that the 3×3 and 3×4 summaries were recomputed
     after publication.
   - Record the exception to the "ported from the paper" rule in
     `DECISIONS.md`.
5. **Fix the source scripts (when approved).** Correct `SpatialSim_3x3.Rmd:455`
   and `SpatialSim_3x4.Rmd:463`, and the 3×3 and 3×4 indices in
   `SimEstimateAnalysisAll.Rmd`. Then decide whether the corrected CSVs in
   `review/results/` replace `results/beta_mse/`, or whether `results/`
   stays frozen as the published record. Keeping it frozen is the safer
   choice.
6. **Not yet rerun.** The α, ψ and ρ CSVs (`results/{alpha,psi,rho}_mse/`)
   for 3×3 and 3×4 Control Only have the same loop bug. The paper reports
   only β̂, so they matter only if those parameters are reported anywhere.
   `code/rerun_trtnospill.R` can save them with a small change.
7. **`code/UnifiedSpatialSim/`** was written after the published results.
   Check that it builds z with the full district range and the correct BSS
   indices before relying on it. Its `build_valid_block_indices()` computes
   the indices programmatically, which avoids the hardcoding problem.
