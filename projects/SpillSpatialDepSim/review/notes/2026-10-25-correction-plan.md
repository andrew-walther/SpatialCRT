# Plan: carrying the 3×3/3×4 correction through (2026-10-25)

Status: **approved by Andrew on 2026-10-25 and implemented the same day**
(see "Decisions" and "What was done" at the end). Nothing is committed or
sent. The corrected Chapter 2 is not yet copied into bios-dissertation.

Inputs read: `review/README.md`, `review/COMPARISON.md`,
`review/tables/corrected_tables_paper_units.md`,
`review/notes/bss-spillover-check.md`, the accepted source
`paper/Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/Revisions_V2c.tex`
(below "V2c"), and Chapter 2
`bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.qmd`
(below "CH2"). Every number below comes from those files or from a check
named next to it.

Ground rules carried from Andrew's instructions: frozen `code/` and
`results/` are never overwritten; the published paper, the bios-dissertation
Chapter 2 and every other bios-dissertation file stay unchanged until Andrew
approves each integration; FAST mode; commit only paths this work changes;
no push without confirmation.

---

## 0. Two things found while planning

**0a. One published CSV comes from a different random stream.** (Resolved
2026-10-25 by option (ii): the scenario was rerun; see `COMPARISON.md`. The
rerun also moved that 3×3 SRS row slightly and the 3×3 optimal-allocation
count from 5 to 4.) For BSS
allocations the two spillover types give identical data, so the Control
Only and Intervention & Control CSVs should agree exactly on the BSS rows.
They do in 31 of the 32 3×3/3×4 scenarios (max abs difference 0). The
exception is 3×3, ψ = 0.8, ρ = 0.01, where all six BSS combos differ by up
to 0.0036 in per-allocation MSE. The corrected Control Only rerun matches
the published Control Only CSV on those combos, so the odd file is the
published `results/beta_mse/beta_mse_results_TrtSpill_3x3_A02_B1_P08_R001.csv`.
Effect on the table: corrected BSS MSE is 0.240 (Control Only) vs 0.242
(Int. & Control) for that one row, SD×10 2.839 vs 2.844. No comparison
changes (SRS 0.168 is lower either way). See open decision D6.
Check used: an inline Rscript comparing BSS rows across the two CSV sets
(rerun it in step 2.1 and save it as `review/code/check_bss_identity.R`).

**0b. The evidence note has its Section 8 pasted twice** (fixed 2026-10-25 in
both copies) (lines 621–640 and
642–661 of the copy). The copy is verbatim; the duplication is in the
bios-dissertation original too. Leave both as they are unless Andrew wants
the original tidied.

---

## Task 1. Ripple into Chapters 1 and 3 (read-only check)

Searched for Project 1/Chapter 2 references, 3×3/3×4, block stratified
performance, "checkerboard", "paradox", collinearity, β − ψ, and
"nearly constant" style claims.

### Chapter 1 (`bios-dissertation/prelim/literature-review/draft/literature-review-draft.qmd`)

| Line | What it says | Change? |
|---|---|---|
| 112–116 | Chapter 2 compares SRS with BSS; examines spillover magnitude, direction, spatial dependence | **No.** Describes the design, not results |
| 455–465 | Exposure identity Wz = d1 + (c−d)z; an alternating checkerboard under rook gives c = 0, d = 1 and exact collinearity | **No.** This is the correct mechanism, and it now matches the corrected Chapter 2 better than the published one did |
| 476 | SRS can produce imbalance with few clusters (cites Raab, Hayes) | **No.** General literature claim |

No mention of the paradox, 3×3/3×4 numbers or BSS "minimal MSE". Chapter 1
needs no change. `trimmed-material.md` and the generated `.tex` were not
searched separately (the `.tex` is rendered from the `.qmd`).

### Chapter 3 (`IncidenceDesign/paper/thesis_chapter/Thesis_Chapter.qmd`; bios-dissertation copy untouched)

| Line | What it says | Change? |
|---|---|---|
| 54 | Ch2 examined BSS "as a practical way to obtain reasonably good estimation while limiting poor-allocation risk in its studied settings" | **No.** Still true after correction (BSS SD of MSE lower in 16/16, 8/16, 14/16 rows for 2×4/3×3/3×4) |
| 64 | Rook contiguity "continues Chapter 2, whose grids used rook neighbors" | No |
| 82 | Ch2 assumed uniform baseline risk | No |
| 90, 112, 125 | Checkerboard = strict BSS of Ch2; 3×3 admitted relaxed allocations; on 2×4 and 3×4 the strict allocations are the two checkerboards | **No.** Correct, and it is the allocation the corrected tables now use (313, 612) |
| 666 | Ch2 compared SRS and BSS on average MSE and spread of MSE | No |
| 672 | Checkerboard under rook: τ not identified, estimates centred on τ − γ | No. Consistent with the corrected 3×4 BSS bias ≈ −ψ |
| 682 | "Relation to Chapter 2": no numerical cross-chapter comparison | No |
| 1207 | "Chapter 2 noted the same consequence of exact collinearity … β̂ ≈ β − ψ" | **No.** Refers to Ch2 Limitations, which stays. After correction Ch2's 3×4 BSS rows actually show it (−0.463 at ρ = 0, ψ = 0.5), which strengthens the link |

Chapter 3 never cites a Chapter 2 number or the paradox. **No change.**
Memory note `feedback_chapter_cross_links` (link via β̂ ≈ β − ψ, "checkerboard
paradox" ≠ rook collinearity) is already satisfied.

### Abstract

| File:line | What it says | Change? |
|---|---|---|
| `prelim/abstract/abstract-text.qmd:12` | "Simple random sampling can achieve lower average mean squared error, but its error varies widely" | No |
| same, :12 | "**Block stratified sampling yields nearly constant error across allocations**, which limits the worst case" | **Should change.** True for 2×4 (BSS SD lower 16/16) and largely 3×4 (14/16), not for corrected 3×3, where BSS SD of MSE is lower in only 8/16 rows and exceeds SRS under Int. & Control (e.g. SD×10 1.002 vs 0.219 at ρ = 0, ψ = 0.5). Proposed: "Block stratified sampling usually varies much less in error across allocations, which limits the worst case." **Upload deadline Mon 2026-11-02** (D3) |
| same, :12 | "Its advantage is largest when spillover flows mainly into control clusters" | No. Holds in all three grids (BSS lower in all ρ = 0.01 Control Only rows) |
| same, :14 | Chapter 3 includes "the block stratified allocation from Chapter 2" | No |
| `prelim/abstract/notes/abstract-fact-sheet.md:33` | "BSS gives near-zero variation in MSE across allocations" sourced to V2c l.126–136 | **Should change** with the abstract; add the corrected counts and point the source at `review/COMPARISON.md` |
| fact sheet :34 | Collinearity bias β̂ ≈ β − ψ | No |
| fact sheet :97 | Same claim as :33 | **Should change** with :33 |
| `prelim/abstract/notes/abstract-extended-versionB.qmd:42` | "nearly constant error across allocations and **performs best** when spillover flows mainly into control clusters" | **Should change if Version B is used.** Same qualifier; "performs best" only holds with ρ > 0 in 2×4 and 3×4, so "has its largest advantage" is safer |

### Project 2 CTJ manuscript (`IncidenceDesign/paper/ctj_manuscript/`)

| File:line | What it says | Change? |
|---|---|---|
| `manuscript_ctj/CTJ_Manuscript.tex:144` | Previous work compared SRS and BSS on small grids, weighing accuracy against poor-allocation risk; cites `walther_spatialcrt_2026` | **No** |
| `supplementary_info/Supplementary_Information.tex:83` | The earlier small-grid study assumed uniform risk | No |
| SI :166–170, :205, :218 | Checkerboard = strict BSS of the earlier study; 3×3 relaxed; 2×4/3×4 checkerboards | No (correct) |
| SI :1813–1814 | Earlier study noted β̂ ≈ β − ψ | No |
| `CTJ_Reading.tex`, `CTJ_Exhibits.tex`, cover letter | No Project 1 result claims | No |
| bib entry `walther_spatialcrt_2026` | Cites the original article | **Optional.** If a Correction is published, add it to the reference or a note (COPE practice is to cite the corrected article; nothing else in CTJ depends on the wrong numbers) |

**Verification for Task 1:** rerun the same `grep` sweep after any
integration and confirm no hit for "paradox" outside notes and frozen
copies; confirm the abstract word count still fits the limit.

---

## Task 2. Corrected Project 1 manuscript + correction memo (`review/manuscript/`)

### 2.1 Set up

1. Copy the V2c folder's `Revisions_V2c.tex`, `references.bib`,
   `Revisions_V2c.bbl` and `figures/` to `review/manuscript/`. Name the
   working copy `Revisions_V2c_corrected.tex`; keep an untouched
   `Revisions_V2c_original.tex` beside it for the diff.
   → verify: `diff` of the original copy against V2c is empty; the copy
   compiles with `/Library/TeX/texbin/latexmk -pdf` to the same page count
   as `Revisions_V2c.pdf`.
2. Save the check from 0a as `review/code/check_bss_identity.R`.
   → verify: it prints 0 for 31 scenarios and 0.0036 for 3×3 ψ = 0.8 ρ = 0.01.

### 2.2 Tables (V2c:444–476 and V2c:500–532)

Replace every cell with `review/tables/beta_mse_table_3x3_corrected.csv`
and `..._3x4_corrected.csv`, same scaling as published (3×3 SD×10; 3×4 SRS
SD×100, BSS MSE×10 and SD×100). Generate the LaTeX rows with a small
script (`review/code/make_latex_tables.R`) rather than typing them.
→ verify: the script in "published" mode regenerates V2c:453–471 and
V2c:509–527 character-for-character (proves the formatting path); every
corrected cell equals `corrected_tables_paper_units.md`. Caption
unchanged except possibly a "(corrected)" note, per D2.

### 2.3 Figures (V2c:438 and V2c:494)

Replace `figures/BetaMSEresults/MSE_Estimates_3x3.png` and `_3x4.png` in
the copy with corrected versions. The `review/figures/*_corrected.png`
files are 1050×675; the published ones are 700×432 and may differ in theme.
Rebuild the corrected 3×3/3×4 figures with the same plotting code, size and
theme as the published ones (adapt `build_tables_figures.R` if needed).
→ verify: the same code in "published" mode produces images visually
indistinguishable from the published 3×3/3×4 PNGs (side-by-side check);
the 2×4 figure stays the original file.

### 2.4 Text the corrections invalidate (rewrite only these)

| # | V2c line | CH2 line | Passage | Why | Planned change |
|---|---|---|---|---|---|
| T1 | 379 | 871–876 | "spatial dependence (ρ = 0.01) consistently increased MSE for both methods … a pattern that persists across all grid sizes" | False for BSS in 2×4 (0.262 → 0.165) and corrected 3×4 (0.226 → 0.133); holds for SRS everywhere and BSS only in 3×3 | Say it raised SRS MSE in every grid; for BSS it depended on the grid. Pre-existing in the 2×4 text, see D1 |
| T2 | 428–432 | 916–935 | 3×3 results paragraphs | BSS no longer better under Int. & Control; all quoted numbers change (and "0.004", "0.0004", "0.259" never matched the published table) | Rewrite from the corrected table: Control Only BSS lower MSE in all 8 rows (e.g. 0.147 vs 0.274 at ρ = 0.01, ψ = 0.5; 0.042 vs 0.114 at ρ = 0); Int. & Control SRS lower in all 8 (0.103 vs 0.147; 0.002 vs 0.042). Note BSS SD of MSE is not smaller here: 6 BSS allocations, including collinear combo 77 |
| T3 | 482–488 | 974–1008 | 3×4 results, incl. "dramatically outperforms … across all ρ" and the "inversion"/paradox | Wrong allocations | Rewrite: same pattern as 2×4 in all 16 rows. Control Only ρ = 0.01 BSS lower (0.133 vs 0.244); Control Only ρ = 0 and all Int. & Control SRS lower (0.133 vs 0.226; 0.112 vs 0.133; 0.001 vs 0.226). Both checkerboards have z = 1 − x, so BSS bias ≈ −ψ (−0.463 at ψ = 0.5, ρ = 0) and the two spillover types give identical BSS rows |
| T4 | 571 | 1150–1164 | Recommendations: qualification "spatial separation enforced by BSS can amplify rather than reduce bias (the 'checkerboard paradox')" | Paradox does not exist | Keep the conditional recommendation (it holds). Replace the paradox clause with: under bidirectional spillover SRS had lower average MSE in every grid, and an allocation in which every control borders a treated cluster loses identifiability (β̂ ≈ β − ψ); keep the buffer-zone pointer |
| T5 | 602 | 1268–1280 | "randomly selected treatment assignments lead to reduced average error (lower bias) … block stratified … consistently result in reduced variation (lower variance)" | Holds for 2×4 and 3×4 (SRS lower \|bias\| 16/16 each), not 3×3 (BSS lower \|bias\| 12/16, BSS SD lower 8/16) | Add "on the 2×4 and 3×4 grids"; one clause on 3×3 |
| T6 | 604 | 1280–1290 | ρ "generally increased MSE for both sampling methods … consistency of the findings across … 2×4, 3×3, and 3×4" | Same as T1 | Same fix as T1; keep the consistency sentence, now true for the main SRS-vs-BSS pattern |
| T7 | 606 | 1290–1300 | "restricted to control clusters only, the MSE … was higher for SRS compared to BSS" | Holds with ρ = 0.01 in all grids; at ρ = 0 only in 3×3 | Add "when spatial dependence was present" and the ρ = 0 exception |
| T8 | 608 | 1301–1318 | "7 different treatment assignment combinations" for each grid; BSS "consistently exhibit minimal MSE … will consistently obtain strong estimates" | 7 holds only for 2×4 (3×3 now 4, 3×4 now 6); "minimal MSE" false in every grid (BSS lower mean MSE in 4/16, 8/16, 4/16 rows) | Give 7, 4, 6; keep "BSS never the single best"; replace the "minimal MSE" sentence with the defensible claim (BSS removes the worst SRS allocations in some regimes, mainly Control Only with ρ > 0). Recompute 5 and 6 in `review/` first (step 4, item A6) |
| T9 | 610–611 | 1319–1336 | Checkerboard paradox paragraph | Mechanism cannot occur under rook; artefact of combos 28/43 | **Delete.** Replace with two or three sentences on the real effect: in both 3×4 checkerboards every control borders a treated cluster, so z = 1 − x and β̂ ≈ β − ψ under both spillover types and both ρ; point to Limitations. The term "checkerboard paradox" is retired |
| T10 | 613 | 1338–1345 | Summary: "unbiased average estimate (favored by SRS) versus … low-variance estimates (achieved with BSS)" | Same as T5 | Light qualifier ("in most settings") |
| T11 | 669 | 1549–1552 | Future: "bidirectional-spillover-with-dependence scenarios where BSS exhibited amplified bias" | Built on the paradox | Replace with "allocations in which treatment and spillover exposure are collinear" |
| T12 | 618 | 1366–1370 | Limitations: "because the standard deviation of the block-stratified estimator is near zero, this error is almost entirely bias" | Still true for the 2×4 and 3×4 checkerboards; not for the six 3×3 BSS allocations, whose MSE varies (SD×10 up to 2.8) | Change "the block-stratified estimator" to "the checkerboard allocations". Optional, see D1 |
| T13 | 132, 135 | (no abstract in CH2) | Abstract: BSS gives "minimal variation in MSE" / "near-zero variation in estimation error across treatment combinations" | Fails for corrected 3×3 | Minimal qualifier ("in most settings" or "on the 2×4 and 3×4 grids"). Changing a published abstract is a bigger step; see D1 |

Left verbatim (checked, still true): Abstract claims about SRS average MSE,
collinearity bias and the spillover-mechanism recommendation; the 2×4
section except T1; Methods V2c:230, 234 (identical BSS results under both
spillover types, now true in all grids), 273, 339, 342; Limitations
V2c:617–624 apart from T12; Appendix C (it describes the intended algorithm,
which is the corrected one); the application and NC map.

Not part of the correction but flagged in the evidence note (wording only,
D1): Aims "geographic distance over which spillover can occur" (V2c:187),
bias definition wording (V2c:359), "grouped bar figures" for box plots
(V2c:366), "SD of BSS estimates" when the column is SD of per-allocation
MSE (V2c:375, 377).

→ verify for 2.4: (a) `latexdiff` original vs corrected
(`/Library/TeX/texbin/latexdiff`) shows changes only in the tables,
figure files and the passages T1–T13 that Andrew approves; (b) every number
in new text is traced to a table cell or `COMPARISON.md` in a checklist
`review/manuscript/number_trace.md`; (c) `grep -i paradox` on the corrected
tex returns nothing; (d) the PDF compiles without undefined references;
(e) a human-writer pass on the new sentences only.

### 2.5 Correction memo (`review/manuscript/correction_notice.md` + cover email draft)

**BMC's current policy** (BMC editorial policies, "Corrections, Retractions
and Matters Arising", read 2026-10-25 at
https://link.springer.com/brands/bmc/editorial-policies#corrections+and+retractions,
formerly biomedcentral.com/getpublished/editorial-policies):
- Corrections are published as a Correction note "bidirectionally linked to
  the original article", and "any alterations to the original article will
  be described in the note".
- A Correction is used for errors that "affect the accuracy of the article,
  but do not substantially invalidate the conclusion", at the Editors'
  discretion. Retraction is for when the conclusion is "substantially
  undermined". Our case fits the first: the general conclusions and the
  recommendation hold (`COMPARISON.md`).
- Authors who find errors should contact the journal through the contact
  details on its website.

**House format**, from published BMC Medical Research Methodology
corrections (e.g. doi:10.1186/s12874-022-01587-9, doi:10.1186/s12874-023-02030-3):
title "Correction: <original title>"; opens "Following publication of the
original article [1], the authors identified …"; states each error and
the corrected text/table/figure; closes "The original article has been
corrected." Reference [1] is the original (doi:10.1186/s12874-026-02968-0).

**Planned contents:**
1. Correction notice (≈300–500 words plus the two replacement tables and
   figures): what was wrong (the two errors, in plain terms; 2×4
   unaffected); how it was found (a post-publication re-check of the
   simulation code while adapting the article for a dissertation chapter;
   exact reproduction of the published CSVs confirmed both errors); what
   changed (Tables 2 and 3, Figures 6 and 7; the 3×3/3×4 results text; the
   checkerboard-paradox paragraph and the sentences T4–T11); that the
   conclusions hold, with the one-line reason (corrected 3×4 matches 2×4 in
   all 16 scenarios; 3×3 in 12 of 16).
2. Before/after table of each changed sentence (as BMC asks that
   alterations be described).
3. Cover email to the MRM editorial office from the corresponding author,
   draft only. Sending is Andrew's/Dr. Lin's action.
4. Optional: offer the corrected code on request, consistent with the
   article's Data Availability statement (V2c:681).
→ verify: every statement in the notice maps to a row of `COMPARISON.md`
or a cell of the corrected tables; the notice lists exactly the changes in
the latexdiff; word count; human-writer pass. Andrew and Dr. Lin read it
before anything goes to the journal (D5).

---

## Task 3. Corrected Chapter 2 (`review/thesis-chapter/`)

1. Copy `project1-spillover-draft.qmd` and `figures/` from
   `bios-dissertation/prelim/project-proposals/project1-spillover/draft/`.
   Resolve its two external dependencies so the copy renders standalone:
   `bios-prelim.cls` is a symlink to
   `bios-dissertation/templates/skeleton/prelim/bios-prelim.cls` (copy the
   file) and `bibliography: ../../../references.bib` (copy the master bib
   and point the path at the copy). Keep an untouched
   `project1-spillover-draft_original.qmd` for diffing.
   → verify: the unedited copy renders with RStudio's Quarto 1.10.18
   (`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`)
   and pdflatex to the same page count as the bios-dissertation PDF.
2. Apply the same table, figure and text changes as Task 2 (CH2 column of
   the T-table; CH2 tables at :944–970 and :1013–1043, figures at :937 and
   :1010). Same sentences, so Chapter 2 and the corrected paper stay in
   step; Ch2 keeps its own markdown/Quarto syntax.
3. Add one line where the 3×3 and 3×4 results begin (or as a footnote to
   the two table captions), e.g. "The 3×3 and 3×4 summaries were recomputed
   after publication to correct a coding error; see Section …". Final
   wording per D4.
   → verify: `diff` original vs corrected `.qmd` shows only the planned
   changes; render passes with no unresolved `@tbl-`/`@fig-` refs; the
   2×4 table and figure are byte-identical to the original; grep finds no
   "paradox".
4. Do **not** touch bios-dissertation. Integration (copying back, recording
   the exception to the "ported verbatim" rule in `DECISIONS.md`, and the
   `check-prose-fidelity.py` check, which will now flag the intended
   differences) waits for Andrew's approval.

Note: `prelim/full-draft/notes/human-writer-review/ch2/ch2-applied.qmd` and
`ch2-edits.md` hold the report-only human-writer edits (ROADMAP Gate 2).
This plan starts from the draft `.qmd`, not `ch2-applied.qmd`. If Andrew
later accepts any of those edits, several touch the same passages (flags
1, 3, 4, M2.5.1), so do the correction first.

---

## Task 4. Everything else that uses the 3×3/3×4 numbers, figures or claims

Status key: **must fix** = wrong and in use / about to be used;
**should fix** = wrong or stale but not on a deadline;
**optional**; **no change**.

### SpillSpatialDepSim (this repo)

| # | Artifact | Finding | Status | Plan → verification |
|---|---|---|---|---|
| A1 | `code/SpatialSim_3x3.Rmd:455`, `SpatialSim_3x4.Rmd:463` | `for (i in 1:8)` loop bug | no change (frozen) | Frozen as the published record. `review/code/rerun_trtnospill.R` holds the fix. Optional later: a fixed copy under `review/code/` if D7 says so |
| A2 | `code/SimEstimateAnalysisAll.Rmd` | `c(28,43)` BSS for every grid; built the published tables/box plots; inline HTML only, no saved outputs | no change (frozen) | Superseded by `review/code/build_tables_figures.R` |
| A3 | `code/SimEstimateAnalysisMeans.Rmd` | Correct BSS indices (e.g. :333) but reads the buggy 3×3/3×4 Control Only CSVs, so its SRS Control Only means are wrong; inline output only | no change (frozen) | Note in AGENTS.md (A10) |
| A4 | `results/beta_mse/` 3×3/3×4 TrtNoSpill CSVs | Loop bug | no change (frozen) | Corrected copies already in `review/results/beta_mse_TrtNoSpill_fixed/` |
| A5 | `results/{alpha,psi,rho}_mse/` 3×3/3×4 TrtNoSpill (48 CSVs) | Same loop bug, not rerun | **open item, not needed** (Andrew 2026-10-25: never reported) | Paper reports only β̂, and nothing found uses these. If wanted: extend `rerun_trtnospill.R` to save α̂, ψ̂, ρ̂ (same runs; ≈30 min per 3×3 scenario, longer for 3×4, run in background) → verify bug mode reproduces the published α/ψ/ρ CSVs to ~1e-8 first |
| A6 | `code/OptimalTreatmentAssignmentCombos.Rmd` | Reads `results/beta_mse/` (buggy for 3×3/3×4 Control Only); source of the "7 optimal allocations" | **must fix (in review/)** | New `review/code/optimal_combos.R` reading corrected Control Only + published TrtSpill CSVs. → verify: on the published CSVs it returns 7 for each grid (reproduces the paper); on corrected inputs 7/4/6 (5 for 3×3 before the 0a rerun); BSS never the single best |
| A7 | `code/CombinedSamplingGrids.Rmd`, `NCTrtMaps.Rmd` | No BSS indices, CSVs or MSE values; grid and map pictures only | no change | — |
| A8 | `code/UnifiedSpatialSim/` | Neither bug: case-1 loop is `seq_len(n_dist)` (`02_simulation_core.R:96`), BSS indices from `build_valid_block_indices()` (`04_run_simulation.R:64`). Never run (no `results/sim_results_unified.rds`). Uses row-major numbering for 2×4 too and its own RNG, so it will not reproduce the published CSVs number-for-number | **should check** | → verify `build_valid_block_indices()` returns 2/6/2 allocations matching the paper's sets (mapping 2×4 indices through the numbering change), and that its 3×3 rule admits the relaxed allocations. Do not run the full simulation unless D7 asks for it |
| A9 | `paper/Manuscript Figures/ConsolidatedMSE/MSE_Estimates_3x3.png`, `_3x4.png`; `paper/Submission_Figures/6_*`, `7_*`; V2c `figures/BetaMSEresults/MSE_Estimates_3x3/3x4.png` | Published wrong figures | no change (historical record) | Corrected versions live in `review/manuscript/figures/` |
| A10 | `paper/Manuscript Figures/BetaMSEtrt*3x3.png`, `*3x4.png` (also in V2c `figures/`, unused by V2c) | Older per-type figures from earlier drafts; likely carry the loop bug | no change | Mention in README as superseded |
| A11 | `README.md:54–75` Key Findings | "BSS gives nearly constant MSE across combinations" fails for corrected 3×3; the post-publication paragraph is already there | **should fix** | Add the 3×3 qualifier; point to `review/PLAN.md` and the corrected manuscript once built |
| A12 | `AGENTS.md:72–80` | Describes the review; doesn't say `SimEstimateAnalysisMeans.Rmd` also uses buggy CSVs, or list the new `review/` subfolders | **should fix** | Add both after Tasks 2–3 (global rule: update README/AGENTS after code changes) |
| A13 | `review/README.md` | Says the full evidence lives in bios-dissertation | **should fix** | Point to `review/notes/bss-spillover-check.md` and add the new folders |
| A14 | SpatialCRT root `README.md`/`AGENTS.md` | Only grid sizes and the IncidenceDesign comparison | no change | — |

### IncidenceDesign

| # | Artifact | Status |
|---|---|---|
| A15 | Thesis chapter and CTJ manuscript | no change (Task 1). No IncidenceDesign code reads SpillSpatialDepSim results (grep of `.R/.Rmd/.qmd` for "SpillSpatialDepSim" found nothing) |

### bios-dissertation (all edits need Andrew's approval; listed, not done)

| # | Artifact | Finding | Status |
|---|---|---|---|
| B1 | `prelim/abstract/abstract-text.qmd:12` + fact sheet :33, :97 | "nearly constant error" (Task 1) | **must fix before 2026-11-02 upload** (D3) |
| B2 | `prelim/abstract/notes/abstract-extended-versionB.qmd:42` | Same, plus "performs best" | should fix if used |
| B3 | Chapter 2 draft (`project1-spillover-draft.qmd/.tex/.pdf`) | Task 3 | **must fix** before the document goes to the committee (≈3 weeks before the 2026-12-07 exam), pending D4 |
| B4 | `prelim/full-draft/_build/body.qmd`, `prelim/full-draft/prelim.tex` | Generated from the chapters; contain the paradox | rebuild after B3 (no hand edits) |
| B5 | `slides/prelim/OUTLINE.md` | Note at :67–90 already retires the paradox and puts 3×3/3×4 in backup; P1-R2 (:233) and P1-R3 (:236) quote published 3×3/3×4 numbers ("BSS 0.001 vs SRS 0.092–0.235", "BSS 0.0004 vs SRS 0.034", "0.223 vs 0.112"); P1-K5 (:281) says use current Ch2 tables "unless updated results … are provided"; open question 7 (:659) and 8 (:660) rest on the published 3×3 pattern and the paradox | **should fix**: swap in corrected numbers (they now show the same pattern on 3×4, a stronger backup slide), close questions 7–8, keep the 2×4-primary decision |
| B6 | Dr. Lin's Project 1 talk deck | Does not exist yet; planned to reuse the prelim part files (`slides/prelim/notes/slide-planning-kickoff-prompt.md:139, 259`) | no change now; inherits B5 |
| B7 | `DECISIONS.md` | No entry for the 2026-10-22 decisions (retire "checkerboard paradox"; 2×4 primary in talk) or for D1–D7 below; the 2026-09-03 "ported from the accepted paper" rule (:420–430) needs an exception once B3 is approved | **must fix** (record decisions) |
| B8 | `ROADMAP.md` | :66 "Chapter 2: done"; :25 "Chapter 2 edits are report-only (adopting any breaks fidelity)"; :58 publication status | **should fix**: add the correction work item, the abstract deadline link, and the journal-correction status |
| B9 | `templates/REPO-TEMPLATE.md` lessons list (:100 onward; 25 lessons; lesson 15 covers porting accepted papers) | Nothing on analysis-code reuse | **should fix**: add a lesson, e.g. "When one script is cloned per grid or setting, grep for hard-coded sizes and indices (loop bounds, allocation IDs) and check summary scripts against each setting's own values; reproduce a published table cell from raw outputs before porting it" |
| B10 | `manuscripts/project1-spillover/project_1_paper.tex` and `prelim/project-proposals/project1-spillover/source/Revisions_V2c.tex` | Reference copies of the accepted paper | no change (record of what was published); optionally add a README line pointing to the correction once published |
| B11 | `prelim/full-draft/notes/human-writer-review/ch2/*`, `literature-review/notes/claude-committee-review-2026-10-02.md`, `project-proposals/notes/chapter2-*.md`, `human-writer-review-plan.md` | Historical notes that quote the paradox or published 3×3/3×4 numbers | no change (dated notes); optionally a one-line pointer in `ch2-review.md` to the evidence note |
| B12 | `prelim/project-proposals/project1-spillover/README.md` | Not checked line by line | check during integration |

---

## Order of work (after approval)

1. Task 2.1–2.3 and A6 (tables, figures, optimal-combo recount), since the
   text in 2.4 and Task 3 depends on them.
2. Task 2.4 text, then Task 3 with the same sentences.
3. Task 2.5 memo.
4. A11–A13 docs; commit in `SpatialCRT` (only `review/` and the two docs),
   one commit per step.
5. Hand Andrew the bios-dissertation list (B1–B12) as a separate approved
   step; B1 first because of the 2026-11-02 deadline.

Rough time: steps 1–4 about one working day; A5 adds a few hours of
background compute if chosen.

---

## Open decisions for Andrew (answered 2026-10-25; see "Decisions" below)

- **D1. Scope of text changes in the corrected paper.** (a) Only what the
  two errors invalidate (T2–T4, T7–T9, T11 and the tables/figures), or
  (b) also the pre-existing statements that are false in the published
  2×4 results or across grids (T1/T6 ρ claim, T5/T10 bias–variance, T8
  "minimal MSE", T12, T13 abstract), or (c) (b) plus the wording flags
  (Aims "geographic distance", bias definition, "grouped bar", SD wording).
  Recommendation: (b). A correction should not leave known-false
  results statements, and each of these is short; (c) is copy-editing a
  journal will likely decline.
- **D2. Table/figure captions.** Mark the replaced tables and figures as
  corrected in the corrected manuscript, or leave captions as published
  (the notice explains it)? Recommendation: leave captions; the Chapter 2
  note covers the thesis.
- **D3. Abstract fix before the 2026-11-02 upload** ("nearly constant
  error" → "usually varies much less in error"). Approve wording, or keep
  as is because the abstract summarises the published paper?
- **D4. Chapter 2 note wording and placement**: one sentence at the start
  of the 3×3 results, or a footnote on the two table captions?
- **D5. Journal correction.** Undecided (Andrew + Dr. Lin). The memo is
  drafted either way; nothing is sent from here.
- **D6. The odd 3×3 ψ = 0.8, ρ = 0.01 TrtSpill CSV (0a).** Options:
  (i) leave it (row differs by 0.002, no conclusion changes), (ii) rerun
  that one scenario so the two spillover types match on BSS (≈30 min), or
  (iii) report the BSS rows from the Control Only run for both types.
  Recommendation: (ii), and mention it in the memo only if the published
  cell changes.
- **D7. Source fixes.** Keep `code/` and `results/` frozen as the
  published record (current rule) and keep all fixed code in `review/`?
  Recommendation: yes. Also: rerun α/ψ/ρ (A5)? Recommendation: no, nothing
  reports them.

---

## Decisions (Andrew, 2026-10-25)

- **D1:** (b). Fix what the two errors break and the results statements that
  are false in the published text (T1–T13). Wording-only flags left alone.
- **D2:** no "corrected" marks on captions.
- **D3:** abstract changed from "nearly constant error" to "much less
  variation in error"; fact sheet and Version B to match.
- **D4:** two sentences in the body (end of the paragraph introducing the
  results), no footnote.
- **D5:** prepare both a corrected manuscript and a memo; decide on
  submission with Dr. Lin. Andrew notes the article is accepted but may not be
  in final published form, so the memo asks for replacement files if that is
  still possible and a Correction note otherwise.
- **D6:** (ii), rerun the stray scenario.
- **D7:** keep `code/` and `results/` frozen; fixed code in `review/`. Do not
  rerun α/ψ/ρ.
- Memo: plain text letter like the cover letters, focused on the code bugs
  and the corrected 3×3/3×4 tables and figures, stating clearly that the
  conclusions do not change; it does not list every sentence edit (the
  latexdiff markup does).
- Plan file moved here from `review/PLAN.md`.

## What was done (2026-10-25)

- `code/rerun_trtnospill.R` gained a `trtspill` mode. It reproduces a published
  TrtSpill file to 3e-9; the stray 3×3 ψ = 0.8, ρ = 0.01 file was rerun.
- `code/build_tables_figures.R` now also writes optimal-allocation counts
  (published inputs give 7/7/7, as in the paper; corrected 7/4/6) and
  700×432 figures.
- `manuscript/`: corrected manuscript (46 pages, compiles cleanly), latexdiff
  markup, correction memo (`.md` and `.pdf`). See its README.
- `thesis-chapter/`: corrected Chapter 2 (52 pages, renders cleanly), flat
  layout (`draft/`, `references.bib`, `manuscript-declarations.md`, README)
  after Andrew's feedback. See its README for the copy-back steps.
- Correction memo dated 2026-10-25.
- bios-dissertation (uncommitted): abstract wording plus re-rendered
  `abstract.pdf` (still 2 pages), fact sheet, Version B, abstract README,
  evidence-note de-duplication, `DECISIONS.md` entry, `ROADMAP.md`,
  `slides/prelim/OUTLINE.md` (corrected numbers, P1-R2/R3/R4 rewritten), and
  template lessons 15 and 30.
- SpatialCRT docs: `README.md`, `AGENTS.md`, `review/README.md`,
  `review/COMPARISON.md`, `paper/Manuscript Revisions/AGENTS.md`.
