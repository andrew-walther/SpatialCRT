# CLAUDE.md — AI Session Context for SpatialCRT

> Cross-project orchestrator. For project-level detail, see:
> - `projects/SpillSpatialDepSim/CLAUDE.md`
> - `projects/IncidenceDesign/CLAUDE.md`

---

## Repository Overview

**SpatialCRT** evaluates treatment assignment designs for **Spatial Cluster Randomized
Trials (CRTs)** where outcomes exhibit spatial heterogeneity and spillover is present.
Application domain: NC law enforcement / SUD prevention policy.

---

## Two Active Projects

| | SpillSpatialDepSim | IncidenceDesign |
|-|--------------------|-----------------|
| **Location** | `projects/SpillSpatialDepSim/` | `projects/IncidenceDesign/` |
| **Grid** | 2×4 / 3×3 / 3×4 (8–12 districts) | 10×10 (100 clusters) |
| **Estimand** | alpha, beta, psi, rho | tau (direct treatment effect) |
| **Question** | Block vs. random assignment with spillover | Which design minimizes MSE across incidence modes? |
| **Status** | **Complete** (original + UnifiedSpatialSim scripts) | **Complete** |
| **Entry point** | `code/SpatialSim_NC_DOC.Rmd` | `code/05_run_simulation.R` |

### How the Projects Relate

SpillSpatialDepSim is the **applied predecessor**: it established the simulation
framework (SAR model, spillover types, block stratification logic) that IncidenceDesign
extended to a larger grid with heterogeneous outcome incidence and 8 formal design
strategies.

---

## Research Focus (IncidenceDesign — PRIMARY)

**Primary question: which treatment assignment design minimizes MSE for tau?**
The 8 designs are: Checkerboard (1), High Incidence Focus (2), Saturation Quadrants (3),
Isolation Buffer (4), 2x2 Blocking (5), Balanced Quartiles (6), Balanced Halves (7),
Incidence-Guided Saturation Quadrants (8).

**MLE (lagsarlm oracle) is the primary estimator.** DIM is a naive baseline only.

Key results (2,560-scenario MLE run, 2026-03-22):
- **Best designs: Design 8 (Incidence-Guided Saturation Quadrants) ≈ Design 3 (Saturation Quadrants)** under MLE (MSE 0.079 vs 0.080)
- **Worst design: Design 1 (Checkerboard)** — MSE 0.744, coverage ~55%
- MLE coverage ~0.94 for all designs except D1; DIM coverage ~0.72

---

## Repository Structure

```
SpatialCRT/
  CLAUDE.md                    # This file
  README.md                    # Human-facing overview
  SpatialCRT.Rproj             # Single .Rproj at root
  projects/
    SpillSpatialDepSim/        # Project 1 (applied, NC DOC context)
      CLAUDE.md  README.md
      code/      data/  results/  paper/
    IncidenceDesign/            # Project 2 (systematic design study)
      CLAUDE.md  README.md
      code/      results/  paper/
  archive/                     # Legacy/exploratory (not maintained)
    README.md
    PreliminarySpatialSim/
    SpatialSim_Unified/
    OutcomeIncidenceDesign_Legacy/
```

---

## Cross-Project Path Reference

From IncidenceDesign code, reference SpillSpatialDepSim results via:
```r
here::here("projects", "SpillSpatialDepSim", "results")
# or relative: ../../SpillSpatialDepSim/results/
```

---

## Getting Started

```r
# Open the project
# File > Open Project > SpatialCRT.Rproj

# IncidenceDesign — load completed results
setwd("projects/IncidenceDesign/code")
source("06_visualizations.R")
mle_results <- load_latest_results(estimation_mode = "MLE_combined")

# SpillSpatialDepSim — reproduce paper results
setwd("projects/SpillSpatialDepSim/code")
rmarkdown::render("SpatialSim_NC_DOC.Rmd")
rmarkdown::render("SimEstimateAnalysisAll.Rmd")
```

---

## Shared Packages

```r
install.packages(c("sf", "spdep", "spatialreg", "dplyr", "tidyr",
                   "ggplot2", "viridis", "rmarkdown", "knitr",
                   "digest", "parallel", "here"))
```
