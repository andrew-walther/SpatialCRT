# review/ — post-publication check of the 3×3 and 3×4 results (2026-10-22)

Isolated sandbox. Nothing here overwrites `code/`, `results/` or the paper;
those are read-only inputs. The accepted paper and Chapter 2 of the
dissertation are unchanged.

## Why this exists

A check of Chapter 2 found two errors that affect the 3×3 and 3×4 results
but not the 2×4 results:

1. **BSS indices.** `code/SimEstimateAnalysisAll.Rmd`, which produced the
   published tables and box plots, takes BSS as combos `c(28, 43)` for every
   grid. Those are the 2×4 checkerboards. The correct sets are
   {28, 30, 35, 51, 77, 106} (3×3) and {313, 612} (3×4), as listed in
   `SpatialSim_3x3.Rmd:203` and `SpatialSim_3x4.Rmd:200`.
2. **Control Only loop bound.** The case-1 spillover loop is `for (i in 1:8)`
   in `SpatialSim_3x3.Rmd:455` and `SpatialSim_3x4.Rmd:463` (should be 1:9
   and 1:12). Treated districts 9 (3×3) and 9–12 (3×4) never pass spillover
   to control neighbours. BSS allocations are unaffected; SRS Control Only
   averages are.

Full evidence: `notes/bss-spillover-check.md` (copy of the bios-dissertation
note). The plan for carrying the correction through, with Andrew's decisions
D1–D7, is `notes/2026-10-25-correction-plan.md`.

## Contents

| Path | What |
|---|---|
| `code/rerun_trtnospill.R` | Reruns one Control Only scenario of the 3×3 or 3×4 script, same seed and RNG order; `bug` mode keeps `1:8`, `fixed` mode uses the full bound; `trtspill` mode reruns an Intervention & Control scenario |
| `code/build_tables_figures.R` | Rebuilds the three tables and figures, `published` and `corrected` versions, the optimal-allocation counts, and paper-size (700×432) figures |
| `code/make_latex_tables.py` | Writes the 3×3 and 3×4 table cells into the manuscript or Chapter 2 from `tables/` |
| `code/text_corrections.py` | The text edits (T1–T14, N1; T14 fixes the funder name in the Declarations), applied to the manuscript `.tex` and the Chapter 2 `.qmd` |
| `results/beta_mse_TrtSpill_rerun/` | Rerun of the one published TrtSpill CSV that came from a different random stream (3×3, ψ = 0.8, ρ = 0.01) |
| `manuscript/` | Corrected manuscript, latexdiff markup, correction memo (see its README) |
| `thesis-chapter/` | Corrected Chapter 2, ready to copy into bios-dissertation (see its README) |
| `notes/` | Evidence note, the correction plan, and `2026-10-28-rho-scale-note.md` (what ρ = 0.01 means on the fitted scale; documentation only, results unchanged) |
| `code/rho_scale_diagnostic.R` | Quick 2×4 re-run behind that note |
| `results/validation/` | `bug`-mode reruns and one `trtspill`-mode rerun; match the published CSVs (max abs MSE difference 3e-9 for 3×3, 1.2e-8 for 3×4) |
| `results/beta_mse_TrtNoSpill_fixed/` | 16 corrected Control Only CSVs (3×3, 3×4 × 4 ψ × 2 ρ) |
| `tables/` | `beta_mse_table_<grid>_<published|corrected>.csv` (unscaled; SD = SD of per-allocation MSE); `optimal_combos_<grid>_<version>.csv` |
| `figures/` | `MSE_Estimates_<grid>_<published|corrected>.png`; `paper_size/` at the published 700×432 |
| `COMPARISON.md` | Published vs corrected, grid by grid |

## How to run

From `review/`:

```bash
xargs -P 8 -L 1 sh -c 'Rscript code/rerun_trtnospill.R $0 $1 $2 fixed' < logs/jobs.txt
Rscript code/rerun_trtnospill.R 3x3 0.8 0.01 trtspill
Rscript code/build_tables_figures.R
```
