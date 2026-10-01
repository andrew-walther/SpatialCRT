# SpatialCRT

**Treatment assignment designs for spatial cluster randomized trials (CRTs)
when outcomes are spatially correlated and treatment spills over to
neighboring clusters.**

Two projects from Andrew Walther's PhD dissertation (Biostatistics,
UNC-Chapel Hill). The dissertation text itself lives in the separate
`bios-dissertation` repo.

> For AI session context and quick technical reference, see [AGENTS.md](AGENTS.md).

---

## Projects at a Glance

| | SpillSpatialDepSim (Project 1) | IncidenceDesign (Project 2) |
|-|--------------------|-----------------|
| **Question** | Block stratified vs. simple random assignment under spillover | Which of 8 designs minimizes MSE for τ under heterogeneous incidence? |
| **Grid** | 2×4 / 3×3 / 3×4 (8–12 clusters) | 10×10 (100 clusters) |
| **Estimand** | alpha, beta, psi, rho | τ (direct treatment effect) |
| **Estimator** | SAR (`lagsarlm`) | Oracle ML spatial lag (`fit_sar_lag()`), non-oracle sensitivity fit |
| **Application** | NC probation (judicial districts) | Sudden unexpected death (SUD) prevention across NC's 58 community college service areas |
| **Status** | **Accepted** at *BMC Medical Research Methodology* | Simulation revised and re-run 2026-09-24 (12,800 scenarios); dissertation chapter rewritten; real-data application in progress |
| **Paper** | `paper/Manuscript Revisions/Revision 2c/Walther_SpatialCRT_LaTeX_Revisions_V2c/Revisions_V2c.pdf` | `paper/dissertation_chapter/` (chapter); `paper/ctj_manuscript/` (*Clinical Trials* draft, still April 2026 numbers) |

SpillSpatialDepSim is the applied predecessor: it set up the simulation
framework (SAR model, spillover types, block stratification).
IncidenceDesign extends it to a 100-cluster grid with heterogeneous baseline
incidence and 8 formal designs.

---

## IncidenceDesign Key Findings (revised simulation, 2026-09-24; τ = 1, oracle ML)

| Design | Queen MSE | Queen coverage | Rook MSE | Rook coverage |
|
---

## Repository Structure

```
SpatialCRT/
  AGENTS.md  CLAUDE.md          # AI session context (CLAUDE.md imports AGENTS.md)
  README.md  LICENSE            # this file; MIT license
  SpatialCRT.Rproj              # single RStudio project at root
  projects/
    SpillSpatialDepSim/         # Project 1
      code/  data/  results/  paper/
    IncidenceDesign/            # Project 2
      code/  results/  docs/  longleaf_setup/
      application/              # NC community-college application (data/, code/, results/)
      paper/
        dissertation_chapter/   # chapter source + tools/ that sync it to bios-dissertation
        ctj_manuscript/         # Clinical Trials manuscript + supplement
        report/                 # project report
  archive/                      # legacy and exploratory work (not maintained)
```

Each project folder has its own `README.md` and `AGENTS.md`.

---

## Getting Started

```r
# File > Open Project > SpatialCRT.Rproj

# IncidenceDesign: load results and generate plots
setwd("projects/IncidenceDesign/code")
source("06_visualizations.R")
mle_results <- load_latest_results(estimation_mode = "MLE_combined")
run_all_visualizations(estimation_mode = "MLE_combined")

# SpillSpatialDepSim: reproduce the applied results
setwd("projects/SpillSpatialDepSim/code")
rmarkdown::render("SpatialSim_NC_DOC.Rmd")
rmarkdown::render("SimEstimateAnalysisAll.Rmd")
```

## Packages

```r
install.packages(c(
  "sf", "spdep", "spatialreg",                    # spatial modeling
  "dplyr", "tidyr", "digest", "parallel",         # simulation
  "ggplot2", "viridis",                           # plots
  "rmarkdown", "knitr", "kableExtra", "tinytex",  # reporting
  "here"                                          # cross-project paths
))
```

---

## Data

The IncidenceDesign application uses restricted county-level death-certificate
data. The source files are git-ignored and never committed (see
`.gitignore` and `projects/IncidenceDesign/application/README.md`); results
aggregated to the 58 community-college clusters are tracked. All other data
in the repo are simulated or public (geography, population).

## License

MIT (see `LICENSE`).

## References

- LeSage, J. & Pace, R.K. (2009). *Introduction to Spatial Econometrics*. CRC Press.
- Bivand, R. et al. — `spdep` and `spatialreg` R packages.
- Project report: `projects/IncidenceDesign/paper/report/IncidenceSpatialCRT_Report.pdf`
