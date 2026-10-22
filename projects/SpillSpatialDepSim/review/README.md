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

Full evidence: `bios-dissertation/prelim/full-draft/notes/human-writer-review/ch2/bss-spillover-check.md`.

## Contents

| Path | What |
|---|---|
| `code/rerun_trtnospill.R` | Reruns one Control Only scenario of the 3×3 or 3×4 script, same seed and RNG order; `bug` mode keeps `1:8`, `fixed` mode uses the full bound |
| `code/build_tables_figures.R` | Rebuilds the three tables and figures, `published` and `corrected` versions |
| `results/validation/` | `bug`-mode reruns; match the published CSVs (max abs MSE difference 3e-9 for 3×3, 1.2e-8 for 3×4) |
| `results/beta_mse_TrtNoSpill_fixed/` | 16 corrected Control Only CSVs (3×3, 3×4 × 4 ψ × 2 ρ) |
| `tables/` | `beta_mse_table_<grid>_<published|corrected>.csv` (unscaled; SD = SD of per-allocation MSE) |
| `figures/` | `MSE_Estimates_<grid>_<published|corrected>.png` |
| `COMPARISON.md` | Published vs corrected, grid by grid |

## How to run

From `review/`:

```bash
xargs -P 8 -L 1 sh -c 'Rscript code/rerun_trtnospill.R $0 $1 $2 fixed' < logs/jobs.txt
Rscript code/build_tables_figures.R
```
