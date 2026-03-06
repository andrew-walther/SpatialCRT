# SpatialCRT

**Simulation studies evaluating treatment assignment designs for Spatial Cluster
Randomized Trials (CRTs)** with application to Sudden Unexpected Death (SUD) in
North Carolina counties.

> For AI session context and quick technical reference, see [CLAUDE.md](CLAUDE.md).

---

## Research Question

When outcome incidence is spatially heterogeneous and spillover effects are present,
which treatment assignment strategy minimizes estimation error for the direct
intervention effect? Specifically, does the *pattern* of how you assign treatment
to geographic clusters matter — and if so, which pattern wins?

---

## Repository Structure

```
SpatialCRT/
  code/
    OutcomeIncidenceDesign/
      ModularIncidenceSim/      <- PRIMARY: Full simulation study (1,920 scenarios each)
      SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd  <- Predecessor (preserved)
    NCdocSpatialSim/            <- Applied: NC DOC context (4x2 grid, 8 districts)
    SpatialSim_Unified/         <- Legacy: Earlier unification efforts
    PreliminarySpatialSim/      <- Legacy: Exploratory analyses
    DistrictAssignments/        <- Assignment tooling
  paper/                        <- Manuscript drafts
```

---

## Primary Work: ModularIncidenceSim

Located in `code/OutcomeIncidenceDesign/ModularIncidenceSim/`.

This is a complete, modular simulation study evaluating **6 treatment assignment
strategies** across a 10x10 lattice grid (100 clusters), 3 incidence generation
modes, and 2 estimation methods. Both simulation runs (DIM and MLE) are complete
as of 2026-03-05.

### Key Findings

**Design 7 (Balanced Quartiles)** is the top-performing design under difference-in-means
estimation, achieving the lowest MSE and highest coverage. **Design 3 (Saturation
Quadrants)** is the winner under MLE. **Design 1 (Checkerboard)** is consistently
worst despite maximizing spatial separation — it creates systematic imbalance when
incidence is spatially patterned.

MLE (spatial autoregressive model) delivers dramatically better coverage (~94%)
than DIM (~72%) by properly modeling spatial autocorrelation, at the cost of
requiring model specification and greater compute.

| Design | DIM Avg MSE | DIM Coverage | MLE Avg MSE | MLE Coverage |
|--------|------------|--------------|------------|--------------|
| 1 — Checkerboard | ~0.31 | ~0.45 | ~0.80 | ~0.55 |
| 2 — High Incidence Focus | ~0.10 | ~0.81 | ~0.13 | ~0.94 |
| 3 — Saturation Quadrants | ~0.11 | ~0.80 | **~0.08** | ~0.94 |
| 4 — Isolation Buffer | ~0.16 | ~0.70 | ~0.15 | ~0.94 |
| 5 — 2x2 Blocking | ~0.13 | ~0.74 | ~0.23 | ~0.94 |
| **7 — Balanced Quartiles** | **~0.10** | **~0.82** | ~0.10 | ~0.94 |

*Averaged across all incidence modes, rho, gamma, spillover type, and neighbor type.
Results are always reported separately by incidence mode in the full analysis.*

### Running the Simulation

**Prerequisites:**
```r
install.packages(c("spdep", "spatialreg", "dplyr", "tidyr",
                   "digest", "parallel", "ggplot2", "viridis"))
```

**DIM estimation (~5-10 minutes):**
```r
setwd("code/OutcomeIncidenceDesign/ModularIncidenceSim")
# Ensure line 40 of 05_run_simulation.R reads: estimation_mode <- "DIM"
source("05_run_simulation.R")
```

**MLE estimation (~10 hours; run from terminal):**
```bash
cd code/OutcomeIncidenceDesign/ModularIncidenceSim
# On macOS, prevent sleep with: caffeinate -i -w $(pgrep -f 05_run_simulation) &
Rscript 05_run_simulation.R
```
In `05_run_simulation.R`, set `estimation_mode <- "MLE"` (line 40).
Set `n_cores > 1` (line 44) for parallel execution across incidence configs.

**Visualizations only (results already saved):**
```r
source("06_visualizations.R")
run_all_visualizations(estimation_mode = "DIM_combined")   # DIM results
run_all_visualizations(estimation_mode = "MLE_combined")   # MLE results
```

**Render the comprehensive PDF report:**
```r
rmarkdown::render("07_results_summary.Rmd",
  output_format = rmarkdown::pdf_document(toc=TRUE, latex_engine="xelatex"),
  output_file = "results/07_results_summary.pdf"
)
```

### Simulation Architecture

The study uses a **Spatial Durbin Model (SDM)** as the data generating process:

```
Y = (I - rho*W)^{-1} * [tau*Z + gamma*Spill(Z) + beta*X + epsilon]
```

where:
- `tau = 1.0` — direct treatment effect (the estimand)
- `rho` — outcome spatial autocorrelation (swept: 0.00, 0.01, 0.20, 0.50)
- `gamma` — spillover magnitude (swept: 0.5, 0.6, 0.7, 0.8)
- `W` — row-standardized spatial weight matrix (rook or queen contiguity)
- `X` — cluster-level incidence (generated under one of three modes)
- `Spill(Z)` — row-standardized mean treatment among neighbors

**Three incidence modes:**
1. **iid Uniform** — X_i ~ Uniform(0,1), independent
2. **Spatial** — SAR filter + pnorm transform; spatially correlated, marginally Uniform
3. **Poisson** — Spatially correlated log-rates → Poisson counts → rank-normalized rates

**Total:** 5 incidence configs × 2 neighbor types × 4 rho × 4 gamma × 2 spillover types × 6 designs = **1,920 scenarios**

### File Inventory

| File | Lines | Purpose |
|------|------:|---------|
| `00_mathematical_specification.Rmd` | ~495 | Complete LaTeX theory doc |
| `01_spatial_setup.R` | 54 | Regular lattice grid, weight matrices |
| `02_incidence_generation.R` | 112 | Three incidence generation modes |
| `03_designs.R` | 137 | Six treatment assignment strategies |
| `04_estimation.R` | 102 | DIM and MLE estimation with CI extraction |
| `05_run_simulation.R` | ~385 | Main orchestrator, nested parameter loop |
| `06_visualizations.R` | ~800 | All plots and summary tables (18 functions) |
| `07_results_summary.Rmd` | ~350 | Comprehensive PDF report: theory + results |

### Results Status

| Simulation | Scenarios | Bias Range | MSE Range | Coverage Range |
|------------|-----------|-----------|-----------|----------------|
| DIM (completed 2026-03-04) | 1,920 | [-0.883, 0.558] | [0.032, 0.799] | [0.00, 0.99] |
| MLE (completed 2026-03-05) | 1,920 | [-0.895, 0.454] | [0.009, 3.910] | [0.00, 1.00] |

MLE convergence failure rate = 0.0 (zero failures across all 480,000 `lagsarlm` calls).

---

## Applied Work: NCdocSpatialSim

Located in `code/NCdocSpatialSim/`.

Applied spatial simulation for the NC Department of Correction context. Uses a
4x2 grid of 8 districts to evaluate treatment assignment and estimate spatial
model parameters (alpha, beta, psi, rho).

**Entry point:** `SpatialSim_NC_DOC.Rmd`

**Key files:**
- `SpatialSim_NC_DOC.Rmd` — main simulation and estimation
- `SpatialSim_3x3.Rmd`, `SpatialSim_3x4.Rmd` — grid size variants
- `CombinedSamplingGrids.Rmd` — multi-grid comparison
- `OptimalTreatmentAssignmentCombos.Rmd` — optimal assignment analysis
- `SimEstimateAnalysisAll.Rmd` / `SimEstimateAnalysisMeans.Rmd` — estimate analysis
- `alpha_mse/`, `beta_mse/`, `psi_mse/`, `rho_mse/` — saved MSE result directories
- `SpatialSim_NC_DOC_TrtSpill.RData`, `SpatialSim_NC_DOC_noTrtSpill.RData` — saved results

---

## Project Lineage

This repository has evolved through approximately 5 stages:

1. **PreliminarySpatialSim** — early exploratory work; grid experiments
2. **NCdocSpatialSim** — applied NC DOC context with specific district structure
3. **SpatialSim_Unified** — attempt to consolidate earlier scripts into a single Rmd
4. **SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd** — monolithic ~600-line Rmd
   studying 6 designs with iid incidence; identified as predecessor to be preserved
5. **ModularIncidenceSim** — clean-room rewrite with 3 incidence modes, 7 bug fixes,
   modular R scripts, proper CI tracking, and full DIM + MLE runs (CURRENT)

---

## Technical Notes

**R packages required for ModularIncidenceSim:**
- `spdep`, `spatialreg` — spatial weight matrices, `lagsarlm()` MLE
- `dplyr`, `tidyr` — data manipulation
- `digest` — per-scenario deterministic seeding
- `parallel` — optional multi-core execution
- `ggplot2`, `viridis` — visualization
- `rmarkdown`, `knitr`, `kableExtra`, `tinytex` — PDF report rendering

**PDF rendering:**
The two Rmd files (`00_mathematical_specification.Rmd`, `07_results_summary.Rmd`)
can be rendered to PDF with `latex_engine = "xelatex"`. If LaTeX is not installed,
run `tinytex::install_tinytex()` first. Output PDFs are saved to `results/`.

**Long-running MLE advice:**
Prevent macOS sleep during 10-hour runs with `caffeinate -i -w <pid>`. Use
`n_cores > 1` in `05_run_simulation.R` to parallelize across incidence configs.

---

## References

- Mirzaei, A. et al. (2019). Sudden unexpected death rates. *(SUD base rate 35/100,000)*
- Gan, W. et al. (2019). County-level Poisson regression for SUD mortality.
- Watson, K. et al. (AHA). Census tract-level spatial analysis of SUD.
- LeSage, J. & Pace, R.K. (2009). *Introduction to Spatial Econometrics*. CRC Press.
- Bivand, R. et al. `spdep` and `spatialreg` R packages.
