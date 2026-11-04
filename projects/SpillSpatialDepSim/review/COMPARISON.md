# Published vs. corrected results (2026-10-22)

What changes when the two errors in `README.md` are fixed. Everything here
is built by `code/build_tables_figures.R` from
`tables/beta_mse_table_<grid>_<published|corrected>.csv`. Full tables in the
paper's units, with published values in brackets wherever they changed, are
in `tables/corrected_tables_paper_units.md`. Figures are in `figures/`.

## Checks that passed

- **One stray published file regenerated (2026-10-25).** For BSS
  allocations the two spillover types give identical data, and the
  published CSVs agree exactly on the BSS rows in 31 of 32 3×3/3×4
  scenarios. The exception, `TrtSpill_3x3_A02_B1_P08_R001.csv`, differed
  from Control Only in every row, so it came from a different random
  stream. `code/rerun_trtnospill.R` in `trtspill` mode reproduces another
  published TrtSpill file (3×3, ψ = 0.5, ρ = 0.01) to 3e-9, and its rerun of
  the stray scenario matches Control Only exactly on all six BSS rows. The
  corrected tables use the rerun (`results/beta_mse_TrtSpill_rerun/`). It
  moves that 3×3 SRS row slightly (MSE 0.168 → 0.167, SD×10 2.149 → 2.133,
  bias 0.101 → 0.102) and the number of distinct optimal 3×3 allocations
  from 5 to 4. No comparison between designs changes.
- **The optimal-allocation count reproduces the paper.** On the published
  inputs, `code/build_tables_figures.R` finds 7 distinct minimum-MSE
  allocations per grid, as the paper states. Corrected: 7, 4 and 6
  (`tables/optimal_combos_<grid>_<version>.csv`). BSS is never the single
  best allocation in any grid.
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
| 7 allocations are each optimal for some scenario; BSS never the single best | 7, BSS never | 4 (was 7), BSS never | 6 (was 7), BSS never | **"BSS never best" holds; "7" holds only for 2×4** |
| 3×3: BSS better in MSE, SD and bias under both spillover types | n/a | ✗ under Int. & Control | n/a | **Does not hold** |
| 3×4: BSS "dramatically outperforms" SRS under Control Only at all ρ | n/a | n/a | ✗ at ρ = 0 (0.226 vs 0.133) | **Does not hold at ρ = 0** |
| "Checkerboard paradox" (BSS backfires only under Int. & Control with ρ > 0, through diagonal spillover) | n/a | n/a | ✗ | **Does not hold.** The real effect is the −ψ collinearity bias, under both spillover types and at both ρ |
| Recommendation for NC probation: BSS when spillover flows mainly into controls | ✓ | ✓ | ✓ (with ρ > 0) | **Holds** |

Net: the paper's general conclusions and its applied recommendation hold.
They are more consistent across grids after the correction than before.
What fails are the grid-specific statements for 3×3 and 3×4, the
"consistently minimal MSE" sentence (already false for 2×4), the count of
7 optimal allocations for 3×3 and 3×4, and the checkerboard paradox.

## Next steps (status as of 2026-11-02)

Only item 1 is still open. The decisions are D1–D7 in
`notes/2026-10-25-correction-plan.md` and in bios-dissertation
`DECISIONS.md` (2026-10-25).

1. **Advisor discussion — open.** Whether to send a correction to BMC is
   to be decided with Dr. Lin. The corrected manuscript, the latexdiff
   markup and the memo are ready in `manuscript/`, and nothing has been
   sent. A four-page before/after brief for that meeting is in
   `change-brief/` (2026-11-04).
2. **Decisions recorded — done.** bios-dissertation `DECISIONS.md`
   (2026-10-25) retires "checkerboard paradox", makes 2×4 primary in the
   talk, and records the Chapter 2 exception to the "ported from the
   paper" rule.
3. **Slides — done.** The Part I deck (`../presentation/`) uses the
   corrected 3×3 and 3×4 results as the results, leads with 2×4, and has
   a collinearity slide. The slides do not mention the correction.
4. **Chapter 2 — done.** Corrected and copied into bios-dissertation on
   2026-10-28. The uncorrected chapter is kept in
   `project1-spillover/draft-uncorrected/`.
5. **Source scripts — decided: no change (D7).** `code/` and `results/`
   stay frozen as the published record. The corrected CSVs live only in
   `review/results/`.
6. **α, ψ and ρ CSVs — decided: not rerun.** The paper reports only β̂,
   so the 3×3 and 3×4 Control Only CSVs for those parameters (which carry
   the loop bug) are left as they are.
7. **`code/UnifiedSpatialSim/` — checked.** It loops over every district
   (`seq_len(n_dist)`), so it does not have the `1:8` bug. It computes BSS
   with the strict rule (no two adjacent districts share a status), which
   gives 2 allocations on 2×4 and 3×4 but only 1 on 3×3; the paper's 3×3
   BSS used the relaxed rule (no two treated districts adjacent; 6
   allocations). It was written after the paper and has not been run.
