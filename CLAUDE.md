# CLAUDE.md — AI Session Context for SpatialCRT

> **Read this file first** when starting a new Claude session on this project.
> For human-readable narrative, see [README.md](README.md).
> For ModularIncidenceSim-specific AI context, see
> `code/OutcomeIncidenceDesign/ModularIncidenceSim/CLAUDE.md`.

---

## What This Repository Is

A research codebase evaluating **treatment assignment designs for Spatial Cluster
Randomized Trials (CRTs)** — specifically which design strategy minimizes estimation
error for a direct treatment effect tau when outcome incidence is spatially
heterogeneous and spillover is present.

**Application domain:** Sudden Unexpected Death (SUD) in North Carolina counties.
Policy-relevant question: if you plan to intervene in some subset of counties/clusters,
how should you select which clusters receive treatment to maximize statistical efficiency?

**Active estimation target:** tau = 1.0 (direct treatment effect, SDM framework).

---

## Repository Status (as of 2026-03-05)

| Directory | Status | Description |
|-----------|--------|-------------|
| `code/OutcomeIncidenceDesign/ModularIncidenceSim/` | **ACTIVE — PRIMARY** | Full simulation study: 6 designs, 3 incidence modes, DIM + MLE, 1920 scenarios each — COMPLETE |
| `code/NCdocSpatialSim/` | **ACTIVE — APPLIED** | NC DOC applied simulation: 4x2 grid, 8 districts, specific policy context |
| `code/OutcomeIncidenceDesign/` (root files) | Legacy | Predecessor monolithic Rmd (~600 lines), preserved untouched |
| `code/SpatialSim_Unified/` | Legacy | Earlier unification attempts |
| `code/PreliminarySpatialSim/` | Legacy | Early exploratory analyses |
| `code/DistrictAssignments/` | Legacy | Assignment-specific tooling |
| `paper/` | Drafts | Manuscript in progress |

---

## Repository Map

```
SpatialCRT/
  CLAUDE.md                          # This file (AI context, top-level)
  README.md                          # Human-facing project overview
  PROJECT_CONTEXT.md                 # Early planning doc (historical, largely superseded)
  SpatialCRT.Rproj                   # RStudio project file
  code/
    OutcomeIncidenceDesign/
      ModularIncidenceSim/           # PRIMARY WORK — see below
        CLAUDE.md                    # Module-level AI context (detailed)
        README.md                    # Module-level human docs
        00_mathematical_specification.Rmd  # Theory doc (rendered HTML + PDF)
        01_spatial_setup.R           # Grid + weight matrices
        02_incidence_generation.R    # 3 incidence modes
        03_designs.R                 # 6 treatment designs
        04_estimation.R              # DIM + MLE estimation
        05_run_simulation.R          # Main orchestrator (sources 01-04)
        06_visualizations.R          # All plots + tables
        07_results_summary.Rmd       # Comprehensive results + theory PDF report
        complete_after_mle.R         # Post-completion script (already run)
        results/                     # All .rds results + .pdf figures
      SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd  # Legacy (unchanged)
      SpatialCRT_Incidence_Sim.Rmd   # Legacy
    NCdocSpatialSim/
      SpatialSim_NC_DOC.Rmd          # Main applied simulation (4x2 grid, 8 districts)
      SpatialSim_3x3.Rmd             # 3x3 grid variant
      SpatialSim_3x4.Rmd             # 3x4 grid variant
      CombinedSamplingGrids.Rmd      # Multi-grid comparison
      OptimalTreatmentAssignmentCombos.Rmd  # Optimal assignment analysis
      SimEstimateAnalysisAll.Rmd     # Simulation estimate analysis
      alpha_mse/ beta_mse/ psi_mse/ rho_mse/  # Saved MSE result files
    SpatialSim_Unified/              # Legacy unification attempts
    PreliminarySpatialSim/           # Legacy exploratory
    DistrictAssignments/             # Assignment tooling
  paper/                             # Manuscript drafts
```

---

## ModularIncidenceSim Quick Reference

**The primary simulation study.** All results complete as of 2026-03-05.

### Simulation Design
- **Grid:** 10x10 regular lattice, N=100 clusters, rook + queen contiguity
- **DGP:** Spatial Durbin Model: `Y = (I - rho*W)^{-1}[tau*Z + gamma*Spill + beta*X + eps]`
- **Estimand:** tau = 1.0 (direct treatment effect)
- **Three incidence modes:** iid Uniform, Spatial (SAR + pnorm), Poisson (count rates, rank-normalized)
- **Six designs:** Checkerboard (1), High Incidence Focus (2), Saturation Quadrants (3), Isolation Buffer (4), 2x2 Blocking (5), Balanced Quartiles (7)
- **Two estimators:** DIM (Neyman variance) and MLE (lagsarlm oracle)
- **1,920 scenarios** per estimator: 5 incidence configs x 2 nb_types x 4 rho x 4 gamma x 2 spill_types x 6 designs

### Key Results
| Estimator | Winner | Loser | Coverage (mean) | MSE range |
|-----------|--------|-------|-----------------|-----------|
| DIM | Design 7 (Balanced Quartiles) | Design 1 (Checkerboard) | 0.72 | [0.032, 0.799] |
| MLE | Design 3 (Saturation Quadrants) | Design 1 (Checkerboard) | ~0.94 | [0.009, 3.910] |

**Design 1 (Checkerboard) is consistently worst** under both estimators despite
maximizing spatial separation — it creates systematic treatment-control confounding
when incidence is spatially structured.

**MLE dramatically improves coverage** (~0.94 vs ~0.72 for DIM) by properly
accounting for spatial autocorrelation.

### Result Files
```
results/sim_results_DIM_combined_20260304_195321.rds   # 1,920 DIM scenarios
results/sim_results_MLE_combined_20260305_150742.rds   # 1,920 MLE scenarios
results/sim_results_{DIM|MLE}_{iid|spatial|poisson}_*.rds  # Per-mode splits
results/{DIM|MLE}_combined_*.pdf                        # Per-config 8-plot PDFs
results/00_mathematical_specification.pdf               # Math theory (PDF)
results/07_results_summary.pdf                          # Comprehensive report (PDF)
```

### Critical Invariants (DO NOT Violate)
1. **Results NEVER aggregated across incidence modes** (iid / spatial / poisson always separate)
2. **Incidence generated ONCE per (mode, rho_X) config** — not per (nb_type, rho) scenario
3. **Deterministic designs (1, 2) generate ONE assignment and replicate** across resamples
4. **Per-scenario seed:** `digest::digest2int(paste(inc_mode, rho_x, nb_type, rho, gamma, spill_type, d_id, sep="|"))`
5. **`base_incidence = X_matrix[, 1]`** — only first column used for design decisions

---

## NCdocSpatialSim Quick Reference

**The applied simulation.** Evaluates treatment assignments for the NC DOC context.

- **Grid:** 4x2 spatial arrangement, 8 districts
- **Entry point:** `SpatialSim_NC_DOC.Rmd`
- **Output parameters encoded in filenames:** e.g., `SpatialSim_NC_DOC_TrtSpill.RData`
- **Estimated parameters:** alpha, beta, psi, rho (see `alpha_mse/`, `beta_mse/`, `psi_mse/`, `rho_mse/` directories)
- **Comparison scripts:** `SimEstimateAnalysisAll.Rmd`, `SimEstimateAnalysisMeans.Rmd`

---

## Key Differences: ModularIncidenceSim vs NCdocSpatialSim

| Aspect | ModularIncidenceSim | NCdocSpatialSim |
|--------|--------------------|-----------------|
| Grid | 10x10 (100 clusters) | 4x2 (8 districts) |
| Estimand | tau (treatment effect) | alpha, beta, psi, rho |
| Incidence modes | 3 (iid, spatial, Poisson) | Single mode |
| Designs evaluated | 6 treatment strategies | Applied NC DOC configurations |
| Estimators | DIM + MLE (lagsarlm) | Various spatial models |
| Status | Complete, documented | Earlier / applied |
| Scale | Simulation study (1920 scenarios) | Applied analysis |

---

## How to Continue from Current State

### Review completed simulation results:
```r
setwd("code/OutcomeIncidenceDesign/ModularIncidenceSim")
source("06_visualizations.R")

# Load combined results:
dim_results <- load_latest_results(estimation_mode = "DIM_combined")
mle_results <- load_latest_results(estimation_mode = "MLE_combined")

# Per-config analysis:
cfgs <- split_by_incidence_config(dim_results)
run_standard_tables(cfgs[["iid Uniform"]], "iid Uniform")
```

### Render the PDF summary report:
```r
setwd("code/OutcomeIncidenceDesign/ModularIncidenceSim")
rmarkdown::render("07_results_summary.Rmd",
  output_format = rmarkdown::pdf_document(
    toc = TRUE, number_sections = TRUE, latex_engine = "xelatex"
  ),
  output_file = "results/07_results_summary.pdf"
)
```

### Planned extensions (not yet run):
- Non-oracle MLE: set `include_spill_covariate = FALSE` in `estimate_tau()` call in `05`
- Heterogeneous populations: set `pop_mode = "heterogeneous"` in `05` Poisson configs
- Grid sensitivity: change `grid_dim = 8` or `grid_dim = 15` in `01`
- Design 6 (Center Hotspot): add `6` to `design_ids` in `05`
- DIM vs MLE joint comparison: load both result sets, join on scenario keys

---

## Packages Required

```r
# Core simulation:
install.packages(c("spdep", "spatialreg", "dplyr", "tidyr", "digest", "parallel"))

# Visualization:
install.packages(c("ggplot2", "viridis"))

# PDF reporting:
install.packages(c("rmarkdown", "knitr", "kableExtra", "tinytex"))
tinytex::install_tinytex()  # If LaTeX not already installed
```

---

## Git History

All work lives on `main`. The primary development branch was `claude/gallant-buck`
(a git worktree), merged to main on 2026-03-05 with 34 files in a single atomic commit.
The predecessor Rmd (`SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd`) was
preserved untouched throughout.
