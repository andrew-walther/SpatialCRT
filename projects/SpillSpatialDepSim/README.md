# SpillSpatialDepSim

**Evaluating treatment assignment strategies for spatial CRTs with spillover and spatial dependence.**

Applied context: NC Department of Corrections (NC DOC) training intervention across
judicial districts, where outcomes (SUD recidivism) exhibit spatial correlation and
treated districts may spill over into adjacent districts.

---

## Research Question

Does block-stratified treatment assignment — which ensures treated and control districts
are spatially separated — improve estimation quality compared to random assignment when
spillover effects are present?

---

## Methods

**Model:** Spatial autoregressive (SAR) lag model via `lagsarlm`:
```
Y = rho * W * Y + alpha * Intervention + psi * Spillover + beta * X + epsilon
```

**Three grid sizes** for robustness:

| Grid | Districts | Treatment | Combos | Reps/scenario |
|------|-----------|-----------|--------|---------------|
| 2×4  | 8         | 4         | C(8,4) = 70  | 50 |
| 3×3  | 9         | 4         | C(9,4) = 126 | 50 |
| 3×4  | 12        | 6         | C(12,6) = 924 | 10 |

**Parameter space:**

| Parameter | Values | Description |
|-----------|--------|-------------|
| alpha | 0.2 (fixed) | Intervention effect |
| beta | 1.0 (fixed) | Covariate effect |
| psi | 0.5, 0.6, 0.7, 0.8 | Spillover magnitude |
| rho | 0.00, 0.01 | Spatial autocorrelation |
| TrtSpill | FALSE, TRUE | Spillover reaches treatment districts? |

**Block stratification:** Valid block-stratified assignments are those where no two
adjacent districts have the same treatment status (rook contiguity). MSE is computed
separately for Random vs. Block assignment groups.

---

## Key Findings

Across grid sizes and parameter configurations:
- Block-stratified assignments do not consistently outperform random assignments
- The optimal treatment combination depends on the specific parameter configuration
- Higher spatial autocorrelation (rho=0.01 vs 0.00) affects MSE for rho estimation
  but has modest effects on beta/psi estimation

See `paper/SpatialCRT_Manuscript_V2.pdf` for the full manuscript.

---

## File Guide

### Original Scripts (frozen, paper-reproducible)

| File | Description |
|------|-------------|
| `code/SpatialSim_NC_DOC.Rmd` | 2×4 grid simulation |
| `code/SpatialSim_3x3.Rmd` | 3×3 grid simulation |
| `code/SpatialSim_3x4.Rmd` | 3×4 grid simulation |
| `code/SimEstimateAnalysisAll.Rmd` | MSE box plots across all treatment combinations |
| `code/SimEstimateAnalysisMeans.Rmd` | Mean MSE ± CI per scenario group |
| `code/OptimalTreatmentAssignmentCombos.Rmd` | Identify minimum-MSE assignment per scenario |
| `code/CombinedSamplingGrids.Rmd` | Visual comparison across grid sizes |
| `code/NCTrtMaps.Rmd` | NC judicial district map visualization |

### Unified Scripts (planned — Plan B)

| File | Description |
|------|-------------|
| `code/01_grid_setup.R` | Programmatic polygon + rook neighbor construction |
| `code/02_simulation_core.R` | `run_scenario()`, `compute_metrics()`, both spillover types |
| `code/03_analysis_functions.R` | Reusable plot and summary functions |
| `code/04_run_simulation.R` | Parallel orchestrator (`mclapply`), all 48 scenarios |
| `code/05_analysis_report.Rmd` | Combined report across all grids and parameters |

### Data

| File | Description |
|------|-------------|
| `data/NC_Judicial_District_Assignments.csv` | County → judicial district mapping for NC |

### Results

```
results/
  alpha_mse/  beta_mse/  psi_mse/  rho_mse/   # 48 CSVs each (192 total)
  SpatialSim_NC_DOC_TrtSpill.RData             # 2x4 full simulation (TrtSpill=TRUE)
  SpatialSim_NC_DOC_noTrtSpill.RData           # 2x4 full simulation (TrtSpill=FALSE)
  sim_results_unified.rds                      # (created by 04_run_simulation.R)
```

---

## Reproduce Paper Results

```r
setwd("projects/SpillSpatialDepSim/code")

# Run simulations (each saves CSVs to ../results/{param}_mse/)
rmarkdown::render("SpatialSim_NC_DOC.Rmd")
rmarkdown::render("SpatialSim_3x3.Rmd")
rmarkdown::render("SpatialSim_3x4.Rmd")

# Generate analysis plots
rmarkdown::render("SimEstimateAnalysisAll.Rmd")
rmarkdown::render("SimEstimateAnalysisMeans.Rmd")
rmarkdown::render("OptimalTreatmentAssignmentCombos.Rmd")
```

## Run Unified Simulation (Plan B)

```r
setwd("projects/SpillSpatialDepSim/code")
source("04_run_simulation.R")  # parallel, all 48 scenarios
# Output: ../results/sim_results_unified.rds (~71,680 rows)
rmarkdown::render("05_analysis_report.Rmd")
```

---

## Packages Required

```r
install.packages(c("sf", "spdep", "spatialreg", "dplyr", "ggplot2",
                   "readr", "tibble", "tidyr", "parallel"))
```

---

## Relationship to IncidenceDesign

This project is the **applied predecessor** of `projects/IncidenceDesign/`. It
established the simulation framework (SAR model, spillover types, block stratification)
that was subsequently extended in IncidenceDesign to a larger 10×10 grid with
heterogeneous outcome incidence and a direct comparison of 6 treatment assignment designs.
