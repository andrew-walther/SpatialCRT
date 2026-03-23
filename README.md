# SpatialCRT

**Evaluating treatment assignment designs for Spatial Cluster Randomized Trials (CRTs)
under spillover and outcome incidence heterogeneity.**

Application domain: NC law enforcement training interventions for SUD prevention.
Core question: *which cluster assignment strategy minimizes estimation error for the
direct treatment effect when outcomes are spatially correlated and spillover is present?*

> For AI session context and quick technical reference, see [CLAUDE.md](CLAUDE.md).

---

## Projects at a Glance

| | SpillSpatialDepSim | IncidenceDesign |
|-|--------------------|-----------------|
| **Grid** | 2×4 / 3×3 / 3×4 (8–12 districts) | 10×10 (100 clusters) |
| **Estimand** | alpha, beta, psi, rho | tau (direct treatment effect) |
| **Designs** | Applied NC DOC configurations | 8 systematic designs |
| **Estimators** | SAR lagsarlm | DIM + MLE (lagsarlm oracle) |
| **Spillover** | TrtNoSpill + TrtSpill | rook + queen contiguity |
| **Status** | Original complete; unified scripts planned | **Complete (2,560 scenarios)** |
| **Report** | `paper/SpatialCRT_Manuscript_V2.pdf` | `paper/report/IncidenceSpatialCRT_Report.pdf` |

---

## Repository Structure

```
SpatialCRT/
  CLAUDE.md                    # AI session context (cross-project)
  README.md                    # This file
  SpatialCRT.Rproj             # Single RStudio project at root
  projects/
    SpillSpatialDepSim/        # Project 1: NC DOC applied simulation
      CLAUDE.md  README.md
      code/      data/  results/  paper/
    IncidenceDesign/            # Project 2: Systematic design study (PRIMARY)
      CLAUDE.md  README.md
      code/      results/  paper/
  archive/                     # Legacy / exploratory work (preserved, not maintained)
    README.md
    PreliminarySpatialSim/
    SpatialSim_Unified/
    OutcomeIncidenceDesign_Legacy/
```

---

## How the Projects Relate

**SpillSpatialDepSim** is the applied predecessor. It evaluated treatment assignments
for the NC DOC context (specific district grid, spillover model, SAR estimation) and
formed the methodological foundation.

**IncidenceDesign** extends the framework to a larger 10×10 grid and systematically
compares 8 formal treatment assignment designs across 3 incidence modes (iid Uniform,
Spatial, Poisson) — addressing the generalizability question: does the optimal design
depend on how outcomes are distributed? A unified project report and modular manuscript
are in development under `paper/`.

---

## IncidenceDesign Key Findings (8-design MLE sweep, 2,560 scenarios)

| Design | Coverage (MLE) | MSE (MLE) |
|--------|---------------|-----------|
| **8 — Incidence-Guided Saturation Quadrants** | **~0.94** | **~0.079** |
| **3 — Saturation Quadrants** | **~0.94** | **~0.080** |
| 6 — Balanced Quartiles | ~0.94 | ~0.10 |
| 2 — High Incidence Focus | ~0.94 | ~0.13 |
| 7 — Balanced Halves | ~0.94 | ~0.13 |
| 4 — Isolation Buffer | ~0.94 | ~0.15 |
| 5 — 2x2 Blocking | ~0.94 | ~0.23 |
| 1 — Checkerboard | ~0.55 | ~0.744 |

**MLE (spatialreg `lagsarlm`) achieves near-nominal coverage (~94%) vs. DIM (~72%)**
by properly modeling spatial autocorrelation and spillover. MLE is the primary estimator;
DIM is a naive baseline. Checkerboard is consistently worst — its perfect spatial
alternation creates systematic confounding when incidence is spatially structured.

*Results reported per incidence mode; table above averaged for reference.
See `paper/report/IncidenceSpatialCRT_Report.pdf` for the comprehensive analysis.*

---

## Getting Started

```r
# Open project
# File > Open Project > SpatialCRT.Rproj

# IncidenceDesign — load completed results and generate plots
setwd("projects/IncidenceDesign/code")
source("06_visualizations.R")
mle_results <- load_latest_results(estimation_mode = "MLE_combined")
run_all_visualizations(estimation_mode = "MLE_combined")

# SpillSpatialDepSim — reproduce applied results
setwd("projects/SpillSpatialDepSim/code")
rmarkdown::render("SpatialSim_NC_DOC.Rmd")
rmarkdown::render("SimEstimateAnalysisAll.Rmd")
```

---

## Packages

```r
install.packages(c(
  # Spatial modeling
  "sf", "spdep", "spatialreg",
  # Core simulation
  "dplyr", "tidyr", "digest", "parallel",
  # Visualization
  "ggplot2", "viridis",
  # Reporting
  "rmarkdown", "knitr", "kableExtra", "tinytex",
  # Cross-project paths
  "here"
))
tinytex::install_tinytex()  # if LaTeX not installed
```

---

## References

- LeSage, J. & Pace, R.K. (2009). *Introduction to Spatial Econometrics*. CRC Press.
- Bivand, R. et al. — `spdep` and `spatialreg` R packages.
- Manuscript: `projects/SpillSpatialDepSim/paper/SpatialCRT_Manuscript_V2.pdf`
- Comprehensive IncidenceDesign report: `projects/IncidenceDesign/paper/report/IncidenceSpatialCRT_Report.pdf`
- Manuscript (in progress): `projects/IncidenceDesign/paper/manuscript/IncidenceSpatialCRT_Manuscript.pdf`
