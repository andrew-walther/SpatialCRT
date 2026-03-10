# CLAUDE.md — AI Session Context for SpillSpatialDepSim

> For human-readable overview, see [README.md](README.md).
> For cross-project AI context, see [../../CLAUDE.md](../../CLAUDE.md).

---

## What This Project Is

**SpillSpatialDepSim** evaluates treatment assignment strategies for spatial cluster
randomized trials (CRTs) in the presence of **spillover effects** and **spatial
dependence**. Applied context: NC Department of Corrections training intervention
across judicial districts, with SUD/recidivism outcomes.

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
| Incidence | Single mode | 3 modes (iid, spatial, Poisson) |
| Designs | Applied NC DOC configurations | 6 systematic designs |
| Estimators | SAR lagsarlm | DIM + MLE (lagsarlm oracle) |
| Status | **Complete** (original + UnifiedSpatialSim) | Complete |

Cross-reference: `here("projects", "SpillSpatialDepSim", "results")` from
IncidenceDesign code.

---

## Packages Required

```r
install.packages(c("sf", "spdep", "spatialreg", "dplyr", "ggplot2",
                   "readr", "tibble", "tidyr", "parallel"))
```
