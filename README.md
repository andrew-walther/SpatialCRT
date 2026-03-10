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
| **Designs** | Applied NC DOC configurations | 6 systematic designs |
| **Estimators** | SAR lagsarlm | DIM + MLE (lagsarlm oracle) |
| **Spillover** | TrtNoSpill + TrtSpill | rook + queen contiguity |
| **Status** | Original complete; unified scripts planned | **Complete (1,920 scenarios)** |
| **Paper** | `paper/SpatialCRT_Manuscript_V2.pdf` | `results/07_results_summary.pdf` |

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
compares 6 formal treatment assignment designs across 3 incidence modes (iid Uniform,
Spatial, Poisson) — addressing the generalizability question: does the optimal design
depend on how outcomes are distributed?

---

## IncidenceDesign Key Findings

| Design | Coverage (MLE) | MSE (MLE) |
|--------|---------------|-----------|
| 1 — Checkerboard | ~0.55 | ~0.80 |
| 2 — High Incidence Focus | ~0.94 | ~0.13 |
| **3 — Saturation Quadrants** | **~0.94** | **~0.08** |
| 4 — Isolation Buffer | ~0.94 | ~0.15 |
| 5 — 2x2 Blocking | ~0.94 | ~0.23 |
| 7 — Balanced Quartiles | ~0.94 | ~0.10 |

**MLE (spatialreg `lagsarlm`) achieves near-nominal coverage (~94%) vs. DIM (~72%)**
by properly modeling spatial autocorrelation and spillover. MLE is the primary estimator;
DIM is a naive baseline. Checkerboard is consistently worst despite maximizing spatial
separation — it creates systematic imbalance when incidence is spatially structured.

*Results reported per incidence mode; table above averaged for reference.*

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
- Comprehensive IncidenceDesign report: `projects/IncidenceDesign/results/07_results_summary.pdf`
