# AGENTS.md — AI Session Context for SpillSpatialDepSim

> For human-readable overview, see [README.md](README.md).
> For cross-project AI context, see [../../AGENTS.md](../../AGENTS.md).

---

## What This Project Is

**SpillSpatialDepSim** evaluates treatment assignment strategies for spatial cluster
randomized trials (CRTs) in the presence of **spillover effects** and **spatial
dependence**. Applied context: North Carolina Department of Adult Correction (formerly
NC DOC) probation interventions assigned across judicial districts.

**Status: final.** Accepted at *BMC Medical Research Methodology*. Final LaTeX source and
PDF: `paper/Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/`
(`Revisions_V2c.tex/.pdf`). It is Chapter 2 of the dissertation (ported into
`bios-dissertation/prelim/project-proposals/project1-spillover/`). Treat code, results,
and paper as read-only unless the author asks otherwise.

**Core question:** Does block-stratified vs. random treatment assignment affect
estimation quality for intervention effects when spillover is present?

---

## Research Design

- **Grid sizes:** 2×4 (8 districts), 3×3 (9 districts), 3×4 (12 districts)
- **Model:** Spatial lag (SAR): `lagsarlm(response ~ intervention + spillover, listw = W_listw)`
- **True parameters:** alpha=0.2, beta=1.0, sd=0.1
- **Spillover types:** TrtNoSpill (case 1: spills to control neighbors only) and TrtSpill (case 2: spills to all neighbors)
- **Parameter space:** psi ∈ {0.5,0.6,0.7,0.8}, rho ∈ {0.00,0.01}
- **Total scenarios:** 3 grids × 4 psi × 2 rho × 2 spillover = 48 per parameter
- **MSE metrics:** bias, bias_abs, variance, MSE — per treatment combination

---

## File Guide

### Original Scripts (FROZEN — paper-reproducible, do not modify)

| File | Description |
|------|-------------|
| `code/SpatialSim_NC_DOC.Rmd` | 2×4 grid, 70 combos, N=50 reps |
| `code/SpatialSim_3x3.Rmd` | 3×3 grid, 126 combos, N=50 reps |
| `code/SpatialSim_3x4.Rmd` | 3×4 grid, 924 combos, N=10 reps |
| `code/SimEstimateAnalysisAll.Rmd` | Box plots of MSE distributions (by grid/psi/rho) |
| `code/SimEstimateAnalysisMeans.Rmd` | Mean MSE + CI plots per scenario group |
| `code/OptimalTreatmentAssignmentCombos.Rmd` | Finds min-MSE assignment per scenario |
| `code/CombinedSamplingGrids.Rmd` | Multi-grid visual comparison |
| `code/NCTrtMaps.Rmd` | NC judicial district map visualizations |

**Note on original scripts:** `SpatialSim_NC_DOC.Rmd` uses an anomalous column-major
district numbering scheme (district 1 = top-right corner). The 3×3 and 3×4 scripts
use standard row-major (district 1 = top-left). The block-stratification indices are
hardcoded: {28,43} for 2×4, {28,30,35,51,77,106} for 3×3, {313,612} for 3×4.

**Note on original analysis scripts:** `SimEstimateAnalysisAll.Rmd` and
`SimEstimateAnalysisMeans.Rmd` only load from `beta_mse/` — not alpha/psi/rho.
They produce inline HTML plots only (no saved figure files).

### Unified Scripts (`code/UnifiedSpatialSim/`)

| File | Description |
|------|-------------|
| `code/UnifiedSpatialSim/01_grid_setup.R` | Programmatic grid geometry + rook neighbors |
| `code/UnifiedSpatialSim/02_simulation_core.R` | `run_scenario()` engine, `compute_metrics()` |
| `code/UnifiedSpatialSim/03_analysis_functions.R` | Plot + analysis functions |
| `code/UnifiedSpatialSim/04_run_simulation.R` | Parallel orchestrator (`mclapply`), saves RDS |
| `code/UnifiedSpatialSim/05_analysis_report.Rmd` | Unified report: all grids + parameters |

### Talk deck (`presentation/`, 2026-10-28)

`presentation/spillover_talk.qmd` is the ~20-minute deck for prelim Part I
(2026-12-07) and the seed of Dr. Lin's conference talk: Quarto revealjs, renders
to one self-contained HTML, full speaker script with timing, 13 backup slides.
Every number is computed at render time from `results/` and `review/` through
`presentation/tools/load_results.R` (author rule, 2026-10-28: build documents from
saved results, never hand-typed numbers). `OUTLINE.md` slide entries are
generated from the deck (`tools/outline_from_deck.py`); `presentation/README.md`
lists the files, build commands and the steps to adapt the deck into
`bios-dissertation/slides/prelim/part-1/`, which happens only after Andrew
approves the deck. Reviews: `review/deck-reviews/`. Framing the author wants:
BSS is a good, not best or worst, design, against rolling the dice with SRS;
per-scenario "optimal allocation" counts are not the point; the slides never
mention the correction.

### Post-publication review (`review/`, 2026-10-22)

A check found two errors affecting the published 3×3 and 3×4 results (2×4
unaffected):
- `SimEstimateAnalysisAll.Rmd` uses BSS combos `c(28,43)` for every grid.
- The Control Only loop is `1:8` in `SpatialSim_3x3.Rmd:455` and
  `SpatialSim_3x4.Rmd:463`.

Corrected reruns, tables, figures and a published-vs-corrected comparison
are in `review/` (`README.md`, `COMPARISON.md`). As of 2026-10-25 it also
holds a corrected manuscript, a correction memo to BMC (not sent), and a
corrected Chapter 2 (copied into bios-dissertation on 2026-10-28, with the
uncorrected chapter kept in `project1-spillover/draft-uncorrected/`; the
corrected manuscript is also in `bios-dissertation/manuscripts/project1-spillover/`
as `project_1_paper_revised.*`); the plan and
decisions are in `review/notes/2026-10-25-correction-plan.md`. The term
"checkerboard paradox" is retired. `SimEstimateAnalysisMeans.Rmd` has the
correct BSS indices but reads the buggy 3×3/3×4 Control Only CSVs, and
`OptimalTreatmentAssignmentCombos.Rmd` does too (its "7" is 7/4/6 once
corrected). The α, ψ and ρ CSVs for 3×3/3×4 Control Only carry the loop
bug and were not rerun (never reported). The frozen scripts and
`results/` are unchanged. Do not edit them, or the paper, without the
author's approval.

**Scale of ρ (2026-10-28, documentation only).** The frozen scripts generate
outcomes with a binary subject-level weights matrix (about 40–60 neighbors per
subject) but fit `lagsarlm` with the row-standardized version, so ρ = 0.01 in the
generator is about 0.5 on the fitted scale (median ρ̂ ≈ 0.50). A quick re-run
suggests the BSS advantage under Control Only with ρ = 0.01 depends on that
difference. The author decided to keep the published results as they are; see
`review/notes/2026-10-28-rho-scale-note.md`.

---

## Results Directory

```
results/
  alpha_mse/  beta_mse/  psi_mse/  rho_mse/   # 48 CSVs each (192 total)
  SpatialSim_NC_DOC_TrtSpill.RData             # 2x4 full run (TrtSpill)
  SpatialSim_NC_DOC_noTrtSpill.RData           # 2x4 full run (TrtNoSpill)
  sim_results_unified.rds                      # (created by UnifiedSpatialSim/04_run_simulation.R)
```

CSV naming convention: `{param}_mse_results_{TrtSpill|TrtNoSpill}_{grid}_A02_B1_P{psi*10}_R{rho*100}.csv`

---

## Data

- `data/NC_Judicial_District_Assignments.csv` — NC county→judicial district mapping

---

## Reproduce Paper Results

```r
setwd("projects/SpillSpatialDepSim/code")
rmarkdown::render("SpatialSim_NC_DOC.Rmd")   # 2x4 grid
rmarkdown::render("SpatialSim_3x3.Rmd")       # 3x3 grid
rmarkdown::render("SimEstimateAnalysisAll.Rmd")  # Analysis (loads from results/)
```

## Run Unified Simulation

```r
setwd("projects/SpillSpatialDepSim/code/UnifiedSpatialSim")
source("04_run_simulation.R")  # all 48 scenarios, parallel
# Output: ../../results/sim_results_unified.rds
rmarkdown::render("05_analysis_report.Rmd")
```

---

## Relationship to IncidenceDesign

SpillSpatialDepSim is the **applied predecessor** of IncidenceDesign. Key differences:

| Aspect | SpillSpatialDepSim | IncidenceDesign |
|--------|-------------------|-----------------|
| Grid | 2×4/3×3/3×4 (8–12 districts) | 10×10 (100 clusters) |
| Estimand | alpha, beta, psi, rho | tau (direct treatment effect) |
| Incidence | Single mode | 3 modes (iid, spatial, Poisson), matched surfaces |
| Designs | Applied NC judicial-district configurations | 8 systematic designs (6 in the manuscripts) |
| Estimators | SAR lagsarlm | Oracle ML spatial lag (`fit_sar_lag()`), non-oracle sensitivity |
| Status | **Accepted** (BMC MRM) | Simulation revised 2026-09-24; chapter rewritten; observed NC SUD application completed 2026-10-02 |

Cross-reference: `here("projects", "SpillSpatialDepSim", "results")` from
IncidenceDesign code.

---

## Packages Required

```r
install.packages(c("sf", "spdep", "spatialreg", "dplyr", "ggplot2",
                   "readr", "tibble", "tidyr", "parallel"))
```
