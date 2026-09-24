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

## IncidenceDesign Key Findings (revised simulation, 2026-09-24; τ = 1, oracle ML)

| Design | Queen MSE | Queen coverage | Rook MSE | Rook coverage |
|--------|-----------|----------------|----------|---------------|
| Saturation Quadrants | 0.090 | 0.94 | 0.071 | 0.94 |
| Incidence-Guided Saturation Quadrants | 0.091 | 0.94 | 0.072 | 0.94 |
| Balanced Halves | 0.130 | 0.94 | 0.084 | 0.94 |
| Balanced Quartiles | 0.131 | 0.94 | 0.084 | 0.94 |
| Isolation Buffer | 0.160 | 0.94 | 0.133 | 0.94 |
| High Incidence Focus | 0.246 | 0.94 | 0.206 | 0.94 |
| 2x2 Blocking | 0.319 | 0.94 | 0.128 | 0.94 |
| Checkerboard | 1.087 | 0.94 | 0.480 | **0.15** |

Queen contiguity is primary, and rook is the sensitivity case. Each pair of adjacent rows
above the Isolation Buffer (the two saturation designs, the two balanced designs) is
statistically indistinguishable, which is why the manuscripts carry 6 designs.

- **Coverage** is near-nominal for every design except Checkerboard under rook. There the
  spillover term is exactly 1 − Z, so τ is not identified.
- **Checkerboard's queen MSE** is almost all variance: the spillover term is ½ for every
  interior cluster, leaving it almost no independent variation.
- **Estimator:** oracle ML spatial lag (`fit_sar_lag()`, validated against `lagsarlm`), with
  a non-oracle sensitivity fit. The April 2026 results are superseded.

*Results reported per incidence mode; table above averaged over incidence modes for reference.
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
