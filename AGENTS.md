# AGENTS.md — AI Session Context for SpatialCRT

> Cross-project orchestrator. For project-level detail, see:
> - `projects/SpillSpatialDepSim/AGENTS.md`
> - `projects/IncidenceDesign/AGENTS.md`

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
| **Status** | **Complete** (original + UnifiedSpatialSim scripts) | **Simulation revised + re-run 2026-09-24**; real NC SUD data aggregated to 58 clusters 2026-09-25; manuscripts being rewritten |
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

**Oracle ML spatial-lag estimator is primary** (validated lean engine `fit_sar_lag()` ≡
`lagsarlm`); a non-oracle fit is a sensitivity analysis. DIM is a pre-revision naive
baseline only.

Key results (2026-09 revision, full re-run 2026-09-24; 12,800 scenarios × 250 fits per
estimator; queen primary; τ = 1; 6 manuscript designs):
- **Best design: Incidence-Guided Saturation Quadrants** (queen MSE 0.091, rook 0.072), then
  Balanced Quartiles (0.131 / 0.084)
- **Worst: Checkerboard** (queen MSE 1.09). Under rook its τ is not identified
  (WZ = 1 − Z): coverage 0.15, estimates kept but flagged
- Coverage ≈ 0.94 for every other design; under queen the rank order is identical at every τ ∈ {0.8, …, 3.0} (under rook, 2x2 Blocking and Isolation Buffer swap at τ = 0.8)
- April 2026 numbers (MSE 0.079 / 0.802, etc.) are superseded; see
  `projects/IncidenceDesign/results/archive/pre_revision_20260924/README.md`
- Manuscripts (CTJ, SI, chapter) still carry the April numbers until they are rewritten
  (manuscript plan step 0 remainder / step 5)

---

## Repository Structure

```
SpatialCRT/
  AGENTS.md                    # This file (CLAUDE.md imports it)
  README.md                    # Human-facing overview
  SpatialCRT.Rproj             # Single .Rproj at root
  projects/
    SpillSpatialDepSim/        # Project 1 (applied, NC DOC context)
      AGENTS.md  README.md
      code/      data/  results/  paper/
    IncidenceDesign/            # Project 2 (systematic design study)
      AGENTS.md  README.md
      code/      results/
      paper/
        report/                  # Unified project report (HTML + PDF)
        manuscript/              # Modular Quarto manuscript (child sections)
        section_drafts/          # Archival LaTeX drafts + bibliography
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
