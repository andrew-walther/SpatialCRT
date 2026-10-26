# Chapter 2 / Project 1: spillover definition vs. 3×4 results

Read-only check, 2026-10-22. Nothing in either repo was edited. Check scripts
and their outputs are in the session scratchpad (not in either repo):
`compare.py`, `zcheck.py`, `claims.py`, `repro3x4.R`, `repro3x3.R`,
`fullrun.R`.

Paths: `SSDS` = `SpatialCRT/projects/SpillSpatialDepSim`; `CH2` =
`bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.qmd`;
`V2c` = `SSDS/paper/Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/Revisions_V2c.tex`
(the published source). Every passage and table quoted below is also in
V2c; line numbers for V2c are given where they matter.

## Summary

1. **The simulation agrees with the text's Rook definition.** In the
   result CSVs, the two 3×4 checkerboard allocations (combos 313 and 612)
   give bit-identical β̂ results under Control Only and Intervention &
   Control, as the definition predicts. The same holds for the six 3×3 BSS
   combos and the two 2×4 BSS combos.
2. **The 3×3 and 3×4 tables report the wrong allocations as BSS.** Their
   BSS columns are computed from combos **28 and 43**. Those are the 2×4
   checkerboard indices, and they were reused for the other two grids. In
   the 3×4 grid, combo 28 treats {1,2,3,4,11,12} and combo 43 treats
   {1,2,3,5,8,12}. Neither is a checkerboard, and both place treated
   clusters next to each other. That is the whole explanation of 0.346 vs.
   2.229: under Intervention & Control the treated clusters adjacent to
   other treated clusters get z = 1, and under Control Only they don't. The
   published Figures (3×3, 3×4) show the same wrong allocations.
3. **A second, smaller code bug.** The Control Only spillover loop in the
   3×3 and 3×4 scripts runs `for (i in 1:8)`, so treated districts 9 (3×3)
   and 9–12 (3×4) never send spillover. This does not affect any BSS
   allocation, but it changes z for 19 of the 126 3×3 allocations and 602
   of the 924 3×4 allocations, so the SRS Control Only rows of those two
   tables are also off. Exact reproduction confirms that the published
   CSVs were generated with this loop.
4. **With the correct BSS allocations, several conclusions change.** In
   3×4 under Intervention & Control, where no Control Only bug applies,
   BSS has higher mean MSE than SRS in all 8 scenarios. The 3×4 Control
   Only comparison cannot be settled until the SRS rows are rerun with the
   loop fixed. In the one scenario rerun, the fix doubles SRS MSE (0.120 →
   0.244), which puts BSS (0.133) ahead. The
   "checkerboard paradox" does not exist under the paper's own model: it
   is a product of the index error. Details and corrected numbers follow.

**Correction needed:** yes, to the published paper's 3×3 and 3×4 tables,
both figures, and the text built on them (3×3 and 3×4 results, the
checkerboard paradox in Discussion and Recommendations, and the summary
claims). The 2×4 table is correct. Whether this goes to the journal as a
correction is Andrew's decision (record it in `DECISIONS.md`).

---

## 1. How the 3×4 BSS allocation is generated

- All 924 allocations come from `combn(1:12, 6)`
  (`SSDS/code/SpatialSim_3x4.Rmd:366-368`). The districts are numbered
  row-major from the top left (`:398-413`) and the Rook neighbour list is at
  `:351-364`.
- `is_valid_assignment_full` (`:227-252`) keeps the allocations in which no
  treated cell has a treated Rook neighbour and no control cell has a
  control Rook neighbour. The note at `:200` and `:344` records the result:
  **313 and 612**, the two checkerboards. I recomputed this independently
  (`zcheck.py`): 313 = treated {1,3,6,8,9,11}, 612 = treated
  {2,4,5,7,10,12}. These are the only allocations that pass either the full
  or the intervention-only separation rule.
- So the 3×4 BSS allocation is a strict checkerboard, as the text says
  (CH2:510-521).

## 2. How each spillover type's indicator is computed

`SSDS/code/SpatialSim_3x4.Rmd:455-473` (the 3×3 and 2×4 scripts have the
same structure):

- **Control Only (case 1), active code:** `points$spillover = 1` if the
  point's district is a Rook neighbour of a treated district `i` **and** the
  point is in a control district, looping `for (i in 1:8)` (`:463-467`).
- **Intervention & Control (case 2), commented out:** the same rule without
  the control-only condition, looping `for (i in 1:12)` (`:470-473`). The
  TrtSpill runs were produced by switching which block is commented. The
  script was not saved in that state.
- In both cases adjacency is Rook only (the `rook_neighbors` list). There is
  no diagonal (Queen) adjacency and no distance term. This matches the
  text's definition (CH2:381-383, CH2:547-549) **except** for the `1:8`
  loop bound in case 1:
  - 2×4 (`SpatialSim_NC_DOC.Rmd:412`): `1:8` is correct (8 districts).
  - 3×3 (`SpatialSim_3x3.Rmd:455`): `1:8` skips district 9. Should be `1:9`.
  - 3×4 (`SpatialSim_3x4.Rmd:463`): `1:8` skips districts 9–12. Should be
    `1:12`.

**Exact reproduction** (`repro3x4.R`, `repro3x3.R`): the RNG was advanced
exactly as in the original loop (seed 2024; per replicate, two `sample()`
calls then one `rnorm()`), and selected allocations were refit with
`lagsarlm` under three rules. ρ = 0.01, ψ = 0.5; values are β̂ MSE per
allocation.

| Grid | Combo | sim, case 1 with `1:8` | sim, case 1 correct | sim, case 2 | CSV TrtNoSpill | CSV TrtSpill |
|---|---|---|---|---|---|---|
| 3×4 | 28 | 0.041440 | 0.025904 | 0.371322 | 0.041440 | 0.371322 |
| 3×4 | 43 | 0.027822 | 0.063029 | 0.074407 | 0.027822 | 0.074407 |
| 3×4 | 313 (BSS) | 0.178555 | 0.178555 | 0.178555 | 0.178555 | 0.178555 |
| 3×3 | 6 | 0.023429 | 0.007193 | 0.100428 | 0.023429 | 0.100428 |
| 3×3 | 11 | 0.039758 | 0.402643 | 0.026572 | 0.039758 | 0.026572 |
| 3×3 | 28 (BSS) | 0.037425 | 0.037425 | 0.037425 | 0.037425 | 0.037425 |

The published CSVs match the simulation to six decimals. The TrtNoSpill
files were made with the `1:8` loop, and the TrtSpill files with the
correct case-2 loop. For a BSS allocation all three rules give the same
data, as the Rook definition predicts.

## 3. Why 0.346 vs. 2.229

The tables were built from the per-allocation CSVs
(`SSDS/results/beta_mse/*.csv`) by `SSDS/code/SimEstimateAnalysisAll.Rmd`.
Its summary (`:1224-1235`) computes `mean(MSE)`, `sd(MSE)` and `mean(bias)`
across allocations, which are the tables' MSE, SD and Bias columns. For
BSS it subsets the CSV with **`[c(28,43),]` for all three grids**, e.g.
`SimEstimateAnalysisAll.Rmd:905` (3×4) and every other `_block <-` line: 16
per grid, 48 in total. The companion script
`SimEstimateAnalysisMeans.Rmd:636-723` uses the correct indices
({28,30,35,51,77,106} for 3×3, {313,612} for 3×4), but the tables don't
come from it.

`compare.py` rebuilds every cell of the three tables from the CSVs:

| Table | SRS columns vs. CSV (all combos) | BSS columns vs. correct BSS combos | BSS columns vs. combos 28, 43 |
|---|---|---|---|
| 2×4 (CH2:885-912) | 16/16 rows match | 16/16 match (28, 43 *are* the 2×4 BSS) | 16/16 |
| 3×3 (CH2:944-970) | 16/16 match | **0/16** | 15/16 (one SD cell off by 0.003 from rounding) |
| 3×4 (CH2:1013-1043) | 16/16 match | **0/16** | 16/16 |

So the 0.346 vs. 2.229 gap is **neither a property of the design nor a
mismatch between the code's definition and the text's**. It comes from
reporting two non-checkerboard allocations as BSS. With the correct
allocations (313, 612), Control Only and Intervention & Control are
identical: MSE×10 = 1.332, SD×100 = 6.416, bias = −0.273 at ρ = 0.01,
ψ = 0.5.

### Corrected BSS columns (β̂; SRS columns as published)

The SRS Control Only rows in 3×3 and 3×4 still contain the `1:8` bug. Their
corrected values for one scenario are in the next subsection.

**3×3** (BSS = combos 28, 30, 35, 51, 77, 106; MSE unscaled, SD ×10)

| Spillover | ρ | ψ | SRS MSE | BSS MSE pub → correct | BSS bias pub → correct | BSS SD pub → correct |
|---|---|---|---|---|---|---|
| Control Only | 0.01 | 0.5 | 0.234 | 0.041 → 0.147 | −0.012 → −0.054 | 0.044 → 1.799 |
| Control Only | 0.01 | 0.8 | 0.372 | 0.069 → 0.240 | −0.025 → −0.083 | 0.145 → 2.839 |
| Control Only | 0 | 0.5 | 0.092 | 0.001 → 0.042 | 0.003 → −0.078 | 0.002 → 1.002 |
| Control Only | 0 | 0.8 | 0.235 | 0.001 → 0.106 | 0.003 → −0.130 | 0.002 → 2.577 |
| Int. & Control | 0.01 | 0.5 | 0.103 | 0.039 → 0.147 | −0.019 → −0.054 | 0.022 → 1.799 |
| Int. & Control | 0.01 | 0.8 | 0.168 | 0.063 → 0.242 | −0.019 → −0.083 | 0.007 → 2.844 |
| Int. & Control | 0 | 0.5 | 0.002 | 0.001 → 0.042 | 0.004 → −0.078 | 0.004 → 1.002 |
| Int. & Control | 0 | 0.8 | 0.005 | 0.001 → 0.106 | 0.004 → −0.130 | 0.004 → 2.577 |

**3×4** (BSS = combos 313, 612; MSE ×10, SD ×100, as in the published table)

| Spillover | ρ | ψ | SRS MSE | BSS MSE×10 pub → correct | BSS bias pub → correct | BSS SD×100 pub → correct |
|---|---|---|---|---|---|---|
| Control Only | 0.01 | 0.5 | 0.120 | 0.346 → 1.332 | −0.152 → −0.273 | 0.963 → 6.416 |
| Control Only | 0.01 | 0.8 | 0.183 | 0.484 → 4.721 | −0.193 → −0.627 | 0.561 → 16.429 |
| Control Only | 0 | 0.5 | 0.034 | 0.004 → 2.255 | −0.006 → −0.463 | 0.022 → 2.274 |
| Control Only | 0 | 0.8 | 0.087 | 0.003 → 6.137 | −0.006 → −0.782 | 0.016 → 1.659 |
| Int. & Control | 0.01 | 0.5 | 0.112 | 2.229 → 1.332 | −0.415 → −0.273 | 20.995 → 6.416 |
| Int. & Control | 0.01 | 0.8 | 0.173 | 3.548 → 4.721 | −0.525 → −0.627 | 33.060 → 16.429 |
| Int. & Control | 0 | 0.5 | 0.001 | 0.003 → 2.255 | −0.003 → −0.463 | 0.004 → 2.274 |
| Int. & Control | 0 | 0.8 | 0.002 | 0.003 → 6.137 | −0.002 → −0.782 | 0.005 → 1.659 |

All 16 rows per grid are in `corrected_tables.md` in the scratchpad.

Why the corrected 3×4 BSS bias is about −ψ: in both checkerboards every
control cluster borders a treated cluster, so z = 1 − x exactly
(`zcheck.py`) and β is not identifiable. β̂ ≈ β − ψ. The chapter's own
Limitations section predicts exactly this (CH2:1361-1364, "an exactly
collinear allocation yields β̂ ≈ β − ψ"). The published BSS bias of
−0.006 at ρ = 0 contradicts that section. The corrected values agree with
it. The same holds for 3×3 combo 77, the one fully separated 3×3
allocation. The other five 3×3 BSS combos are not collinear.

### Size of the `1:8` bug in the SRS Control Only rows

Full rerun of one scenario (Control Only, ρ = 0.01, ψ = 0.5) with both
loop bounds (`fullrun.R`):

- **3×3:** with the `1:8` loop, the rerun reproduces all 126 CSV MSEs to
  within 3e-9, which confirms the bug was in the published run. With the
  loop corrected, the SRS row becomes mean MSE **0.274** (published 0.234),
  SD of MSE 0.263 (published 0.259) and bias **−0.089** (published
  −0.054). BSS (0.147) still has the lower MSE in this scenario.
- **3×4:** with the `1:8` loop, the rerun reproduces all 924 CSV MSEs to
  within 1.2e-8. With the loop corrected, the SRS row becomes mean MSE
  **0.244** (published 0.120), SD of MSE 0.244 (published 0.132) and bias
  **−0.190** (published +0.033). The bug therefore roughly halved the
  published 3×4 SRS Control Only MSE. With the fix, BSS (0.133) has the
  lower MSE in this scenario, which reverses the sign of the comparison
  with uncorrected SRS (0.120).

**Consequence:** every comparison below that pairs corrected BSS with a
Control Only SRS row in 3×3 or 3×4 uses an SRS value that is still wrong.
Those comparisons are provisional until the 16 Control Only scenarios are
rerun. The Intervention & Control rows and the 2×4 rows are not affected.

## 4. The other flagged contradictions (ch2-review.md)

**Flag 1: "restrictions on the geographic distance over which spillover
can occur"** (CH2:251-253; V2c:187). **The text is wrong.** The code has
no distance parameter. The only factor that varies is which clusters may
receive spillover (control only vs. control and intervention), always by
Rook adjacency (`SpatialSim_3x4.Rmd:455-473`). The Simulation Study
already describes it correctly. Correct the Aims sentence (the review's
"restrictions on which clusters can receive spillover" works).

**Flag 3: bias definition, "the parameter estimate is quite close to the
true value"** (CH2:809-811; V2c:359). **Not a text-vs-table
contradiction; the wording is imprecise.** The code computes bias as
`mean(estimates) - true_value` per allocation (`SpatialSim_3x4.Rmd:592`),
and the tables average that across allocations. "Small bias means that the
estimates are, on average, close to the true value" matches the code. This
is optional wording, not a correction.

**Flag 4: "SRS … lower bias … BSS … lower variance"** (CH2:1272-1278;
V2c:602; repeated in the Recommendations and the Discussion summary).
**Neither the published tables nor the corrected ones support it as a
general claim.** Counts over the 16 rows per grid (`claims.py`):

| Grid | SRS \|bias\| smaller (pub / correct) | BSS SD smaller (pub / correct) | BSS mean MSE smaller (pub / correct) |
|---|---|---|---|
| 2×4 | 16 / 16 | 16 / 16 | 4 / 4 |
| 3×3 | 0 / 4 | 16 / 8 | 16 / 8 |
| 3×4 | 10 / 16 | 12 / 14 | 12 / **0** |

The "correct" counts use corrected BSS with the published SRS rows. For
the 3×3 and 3×4 Control Only rows (8 of each 16), those SRS values still
contain the `1:8` bug (Section 3), so those counts are provisional.

The claim holds for 2×4. With the corrected numbers it also holds for 3×4,
but it fails for 3×3, where BSS has the smaller bias under Control Only.
Also, the tables' "SD" is the SD across allocations of each allocation's
MSE (`SimEstimateAnalysisAll.Rmd:1228`). It is not the SD of β̂ across
replicates, so "lower variance" and "SD of BSS estimates" (CH2:855) describe
a quantity the table does not report. **The text needs correcting.** Once
the tables are fixed, state the pattern by grid and spillover type, and
describe the SD column as variation in MSE across allocations.

**M2.5.1: "consistently exhibit minimal MSE … will consistently obtain
strong estimates"** (CH2:1309-1318; V2c:608). **The text is wrong, and the
correction does not rescue it.** In the rows that need no rerun, BSS has
lower mean MSE than SRS in 4 of 16 (2×4), and in 0 of 8 Intervention &
Control rows in both 3×3 and 3×4. The Control Only rows for 3×3 and 3×4
wait on the rerun; the one rerun scenario in each favours BSS. The
review's proposed "usually low MSE" is therefore not supported either. The parts of
the paragraph that are true: in each grid exactly **7** distinct
allocations are the minimum-MSE allocation for at least one of the 16
scenarios, and none of them is a BSS allocation (`claims.py`). That also
confirms flag 6, which the review left unchecked. The defensible claim is
narrower: BSS removes the worst SRS allocations in some regimes (lower SD
of MSE across allocations in 38 of 48 rows). It does not give low MSE in
general: under Intervention & Control its mean MSE is higher than SRS in
every 3×3 and 3×4 scenario.

## 5. Other passages that depend on the wrong tables

- **3×3 results** (CH2:919-934): "BSS consistently yielded significantly
  lower MSE, SD, and bias" in both spillover types. This is false for
  Intervention & Control once corrected (BSS 0.147 vs. SRS 0.103 at
  ρ = 0.01, ψ = 0.5; 0.042 vs. 0.002 at ρ = 0). Separately, the text's
  "0.004 for ψ=0.5" at ρ = 0.01 does not match even the published table
  (0.041), and neither does "around 0.0004" (0.001).
- **3×4 results** (CH2:976-1007): every BSS number and the description of
  the "inversion" are products of combos 28 and 43. Corrected, under
  Intervention & Control, BSS mean MSE is higher than SRS in all 8
  scenarios and far higher at ρ = 0 (0.226 vs. 0.001), because both
  checkerboards are exactly collinear. The Control Only comparison waits
  on the SRS rerun. At ρ = 0.01, ψ = 0.5 it favours BSS (0.133 vs. 0.244
  corrected). At ρ = 0, BSS is 0.226 against a published, still-buggy
  SRS 0.034.
- **Checkerboard paradox** (CH2:997-1001, 1159-1160, 1319-1335; V2c:481,
  486, 571, 611). Under the stated Rook rule, spillover cannot reach the
  "diagonally adjacent" treated clusters (CH2:1325), so for a strict
  checkerboard the two spillover types produce identical data. The CSVs
  confirm this. The paradox has no basis in the model and should be
  removed, along with the qualification in the Recommendations built on
  it (CH2:1150-1164). A real 3×4 caveat does exist, but it is a different
  one: identifiability is lost because z = 1 − x.
- **Flag 2 / line 417** (CH2:416-420; V2c:234): the intended statement
  (identical results under BSS for both spillover types) is **true** for
  all three grids once the correct allocations are used.
- **Discussion, spillover types** (CH2:1290-1300): "When spillover was
  restricted to control clusters only, the MSE … was higher for SRS
  compared to BSS". Corrected, this holds for 3×3 but not for 3×4, and it
  holds for 2×4 only at ρ = 0.01.
- **Abstract / prelim abstract:** not checked here. Any claim about the
  checkerboard paradox or about BSS in 3×4 should be traced through
  `prelim/abstract/notes/abstract-fact-sheet.md`.

## 6. What needs correcting

| Item | Published paper (V2c) | Chapter 2 |
|---|---|---|
| 2×4 table and text | no change | no change |
| 3×3 table, BSS columns | wrong allocations, needs replacing | same |
| 3×4 table, BSS columns | wrong allocations, needs replacing | same |
| 3×3 and 3×4 SRS Control Only rows | `1:8` bug, needs a rerun with the corrected loop | same |
| Figures for 3×3 and 3×4 (BSS boxes) | same index error: the BSS boxes sit near 0 at ρ = 0 | same |
| 3×3 and 3×4 results text | rewrite from corrected numbers | same |
| Checkerboard paradox (Results, Discussion, Recommendations) | remove; replace with the collinearity caveat | same |
| Flag 1 (Aims, "geographic distance") | text wrong | same |
| Flag 4 (bias/variance trade-off) | text over-generalises | same |
| M2.5.1 ("consistently minimal MSE") | text wrong | same |
| Flag 3 (bias definition) | wording only | optional |

The 3×3 and 3×4 Control Only scenarios (16 CSVs per parameter) need a
rerun with the `1:8` loop fixed before the SRS rows can be final. The
TrtSpill CSVs and all BSS rows are correct as they stand; only the indices
used to summarise them were wrong. For scale: one 3×3 scenario took
about 30 minutes on a laptop with the scratchpad script, and one 3×4
scenario took longer.

## 7. Re-verification after Andrew's feedback (2026-10-22)

Andrew asked whether the review mistook which script generated which grid.
It did not. Each grid's per-allocation CSVs come from its own script, and
the checks above use those scripts:

| Grid | Simulation script | Case-1 (Control Only) loop | Correct bound |
|---|---|---|---|
| 2×4 | `SpatialSim_NC_DOC.Rmd:412` | `for (i in 1:8)` | 1:8 ✓ |
| 3×3 | `SpatialSim_3x3.Rmd:455` | `for (i in 1:8)` | 1:9 ✗ |
| 3×4 | `SpatialSim_3x4.Rmd:463` | `for (i in 1:8)` | 1:12 ✗ |

In each script the neighbour/treatment columns just above the loop use the
full range (`1:9` at `SpatialSim_3x3.Rmd:448`, `1:12` at
`SpatialSim_3x4.Rmd:456`), and so does the commented case-2 loop (`:462`,
`:470`). Only the case-1 loop kept the 2×4 bound. Git history: in the first
3×4 commit (d3d74e5, 2025-04-30) case 1 was commented out with `1:8` and
case 2 was active with `1:12`. Commit 73eb744 (2025-06-19, "alter TrtSpill
setting") switched to case 1 without changing the bound. The exact
reproduction (Section 2) used each grid's own script logic: 3×3 with 180
points, 4 treated, 3-column numbering; 3×4 with 240 points, 6 treated,
4-column numbering. It reproduces those grids' CSVs to 1e-8.

The per-grid simulation scripts do contain the correct BSS indices (3×3
`:203`, `:344`; 3×4 `:200`, `:344`), but only in notes and plotting code.
They write per-allocation CSVs and compute no SRS/BSS summary. The summary
statistics in the tables come from a separate analysis step. Two analysis
scripts read those CSVs:

- `SimEstimateAnalysisMeans.Rmd`: correct indices per grid
  (`:636-723` for 3×4).
- `SimEstimateAnalysisAll.Rmd`: `c(28,43)` for every grid, from its first
  commit (7bce9c0, 2025-06-19). It never contained 313. Its
  `mean(MSE)`/`sd(MSE)`/`mean(bias)` summary has the same three columns as
  the tables.

Whichever script was run, the published numbers decide the question. All
96 BSS cells of the 3×3 and 3×4 tables (32 rows × MSE, SD, bias) equal the
summary over CSV rows 28 and 43 of that grid's own CSVs, and none equals the
summary over the correct rows (Section 3). Those numbers have been in the
paper since the original submission (`paper/Walther_Submission/Walther_SpatialCRT_LaTeX/main.tex:1008,1119,1129`).

**Checkerboard paradox vs. the −ψ bias.** The −ψ bias from z = 1 − x is
real, but it differs from the paradox the paper describes in three ways:

- The paper attributes the paradox to spillover reaching diagonal treated
  neighbours (CH2:1325); the −ψ bias comes from collinearity.
- The paper ties the paradox to Intervention & Control with ρ > 0. The
  −ψ bias occurs under both spillover types, because they give identical
  data for a checkerboard, and at both ρ values.
- The published table shows BSS bias of −0.006 at ρ = 0. With the
  correct allocations it is −0.463.

The Limitations paragraph (CH2:1355-1376) describes the real effect
correctly.

**Flag 1.** The sentence (CH2:249-253) lists "variation of restrictions on
the geographic distance over which spillover can occur" among the factors
the simulations "explore". Read as an example of how spillover could be
modelled, it is harmless. Read as a list of varied factors, which is how it
is written, it describes a factor that was not varied. This is a wording
call for Andrew, not an error in results.

## 8. Corrected results (2026-10-22, later)

All 16 Control Only scenarios for 3×3 and 3×4 were rerun with the loop
bound fixed, and all three tables and figures were rebuilt. This was done
in the isolated folder `SpatialCRT/projects/SpillSpatialDepSim/review/`.
Its `COMPARISON.md` is the current account and **supersedes the
provisional counts in Sections 3–5 of this note.** Main points:

- **3×4 now shows the same pattern as 2×4 in all 16 scenarios:**
  - Intervention & Control: SRS has lower MSE.
  - Control Only with ρ = 0.01: BSS has lower MSE.
  - Control Only with ρ = 0: SRS has lower MSE.
- **3×3 agrees except under Control Only with ρ = 0**, where BSS also has
  lower MSE.
- **The general conclusions and the NC recommendation hold.** What fails:
  - the grid-specific 3×3 and 3×4 text,
  - the checkerboard paradox,
  - "consistently minimal MSE",
  - the "7 optimal allocations" count for 3×3 (now 4; 5 before one stray TrtSpill file was regenerated on 2026-10-25) and 3×4 (now 6).
- **The bias/variance claim (flag 4)** holds for 2×4 and 3×4, but not 3×3.
