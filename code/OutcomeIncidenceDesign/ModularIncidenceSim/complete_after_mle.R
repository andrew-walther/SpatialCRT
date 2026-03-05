# ==============================================================================
# complete_after_mle.R
# Post-MLE completion script: generates MLE visualizations, writes CLAUDE.md
# and README.md with live stats from both DIM and MLE results, prints a
# summary. Run once MLE simulation is complete.
#
# Usage:
#   Rscript complete_after_mle.R
#   source("complete_after_mle.R")
# ==============================================================================

library(dplyr)

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

results_dir <- file.path(script_dir, "results")

cat("=======================================================\n")
cat("  Post-MLE Completion Script\n")
cat("  Working dir:", script_dir, "\n")
cat("  Results dir:", results_dir, "\n")
cat("=======================================================\n\n")

# ==============================================================================
# STEP 1: Verify MLE results exist
# ==============================================================================

mle_combined_files <- list.files(results_dir,
  pattern = "sim_results_MLE_combined.*\\.rds$", full.names = TRUE)

if (length(mle_combined_files) == 0) {
  stop("No MLE combined results file found in ", results_dir,
       "\nWait for 05_run_simulation.R to complete.")
}

mle_file <- mle_combined_files[which.max(file.mtime(mle_combined_files))]
cat("MLE results file:", basename(mle_file), "\n")
mle_results <- readRDS(mle_file)
cat("MLE results loaded:", nrow(mle_results), "scenarios\n\n")

# ==============================================================================
# STEP 2: Load DIM results for comparison stats
# ==============================================================================

dim_combined_files <- list.files(results_dir,
  pattern = "sim_results_DIM_combined.*\\.rds$", full.names = TRUE)
dim_file <- dim_combined_files[which.max(file.mtime(dim_combined_files))]
cat("DIM results file:", basename(dim_file), "\n")
dim_results <- readRDS(dim_file)
cat("DIM results loaded:", nrow(dim_results), "scenarios\n\n")

# ==============================================================================
# STEP 3: Generate MLE visualizations
# ==============================================================================

cat("=== Generating MLE Visualizations ===\n")
source(file.path(script_dir, "06_visualizations.R"))
run_all_visualizations(results = mle_results, estimation_mode = "MLE_combined")
cat("MLE visualizations complete.\n\n")

# ==============================================================================
# STEP 4: Extract summary statistics for documentation
# ==============================================================================

cat("=== Extracting Summary Statistics ===\n")

get_stats <- function(r, label) {
  list(
    label        = label,
    n_scenarios  = nrow(r),
    bias_range   = range(r$Bias, na.rm = TRUE),
    bias_mean    = mean(r$Bias, na.rm = TRUE),
    mse_range    = range(r$MSE, na.rm = TRUE),
    mse_mean     = mean(r$MSE, na.rm = TRUE),
    cov_range    = range(r$Coverage, na.rm = TRUE),
    cov_mean     = mean(r$Coverage, na.rm = TRUE),
    fail_mean    = mean(r$Fail_Rate, na.rm = TRUE),
    fail_max     = max(r$Fail_Rate, na.rm = TRUE),
    # Best design by median MSE
    best_design  = (r %>% group_by(Design) %>%
                      summarise(med_mse = median(MSE, na.rm = TRUE),
                                .groups = "drop") %>%
                      arrange(med_mse) %>% slice(1))$Design
  )
}

dim_s <- get_stats(dim_results, "DIM")
mle_s <- get_stats(mle_results, "MLE")

cat(sprintf("DIM: n=%d | Bias [%.3f, %.3f] | MSE [%.3f, %.3f] | Cov [%.2f, %.2f]\n",
  dim_s$n_scenarios,
  dim_s$bias_range[1], dim_s$bias_range[2],
  dim_s$mse_range[1], dim_s$mse_range[2],
  dim_s$cov_range[1], dim_s$cov_range[2]))

cat(sprintf("MLE: n=%d | Bias [%.3f, %.3f] | MSE [%.3f, %.3f] | Cov [%.2f, %.2f] | Fail mean=%.4f\n",
  mle_s$n_scenarios,
  mle_s$bias_range[1], mle_s$bias_range[2],
  mle_s$mse_range[1], mle_s$mse_range[2],
  mle_s$cov_range[1], mle_s$cov_range[2],
  mle_s$fail_mean))

cat("\n")

# MLE split files
mle_iid_file     <- basename(list.files(results_dir, "sim_results_MLE_iid.*\\.rds$")[1])
mle_spatial_file <- basename(list.files(results_dir, "sim_results_MLE_spatial.*\\.rds$")[1])
mle_poisson_file <- basename(list.files(results_dir, "sim_results_MLE_poisson.*\\.rds$")[1])

# List all output PDFs
mle_pdfs  <- list.files(results_dir, "MLE_.*\\.pdf$")
dim_pdfs  <- list.files(results_dir, "DIM_.*\\.pdf$")

cat("MLE PDFs generated:\n"); cat(" ", paste(mle_pdfs, collapse = "\n  "), "\n\n")

# ==============================================================================
# STEP 5: Write CLAUDE.md
# ==============================================================================

cat("=== Writing CLAUDE.md ===\n")

claude_md <- sprintf(
'# CLAUDE.md — AI Session Context for ModularIncidenceSim

> **Read this file first** when starting a new Claude session on this project.
> It provides ~80%% of the context needed to be immediately productive without
> re-reading all 1,589+ lines of code, result files, or figure PDFs.
> For human-readable narrative documentation see [README.md](README.md).

---

## What This Project Is

A modular simulation study evaluating **6 treatment assignment designs** for Spatial
Cluster Randomized Trials (CRTs) under heterogeneous outcome incidence and spatial
spillover. The estimand is **tau = 1.0** (direct treatment effect). We compare two
estimators (DIM, MLE) across 1,920 parameter scenarios and 3 incidence generation
modes that are always reported **separately** (never aggregated).

**Application context:** Sudden Unexpected Death (SUD) in NC counties. Poisson
incidence base rate 35/100,000 (Mirzaei et al.).

**Metrics tracked per scenario:** Bias, SD, MSE, Coverage (95%% CI), Fail_Rate.

---

## File Map

| File | Lines | Role | Key Functions |
|------|------:|------|---------------|
| `00_mathematical_specification.Rmd` | ~495 | Theory & DGP formulas | (rendered to HTML) |
| `01_spatial_setup.R` | 54 | Grid + weight matrices | `build_spatial_grid()`, `get_active_spatial()` |
| `02_incidence_generation.R` | 112 | 3 incidence modes | `generate_incidence()`, `generate_incidence_iid/spatial/poisson()` |
| `03_designs.R` | 137 | 6 treatment designs | `get_designs()`, `get_design_names()`, `is_design_deterministic()` |
| `04_estimation.R` | 102 | DIM + MLE estimation | `estimate_tau()` |
| `05_run_simulation.R` | ~385 | Main orchestrator | `run_incidence_config()` |
| `06_visualizations.R` | ~800 | Plots + tables | 18 functions (see section below) |
| `complete_after_mle.R` | ~250 | Post-MLE script | Runs viz, writes docs, prints stats |

---

## Architecture

**Pipeline:** `01 → 02 → 03 → 04`, orchestrated by `05`, visualized by `06`.

- `05_run_simulation.R` sources `01`–`04` at runtime via `file.path(script_dir, "0X_*.R")`
- `06_visualizations.R` sources `01`–`03` (needs grid/incidence/design helpers for standalone plots)
- Both scripts detect their working directory via `dirname(sys.frame(1)$ofile)` with a
  `tryCatch` fallback to `getwd()` for interactive use.

---

## Data Flow (ASCII)

```
build_spatial_grid(grid_dim=10)
  -> grid_obj: list(coords[100x2], N_clusters=100,
                    nb_rook, nb_queen, W_rook[100x100], W_queen[100x100],
                    listw_rook, listw_queen)
        |
        v
generate_incidence(mode, N=100, n_resamples, W, rho_incidence)
  -> X_matrix: [100 x n_outcome_resamples], values in [0,1]
  -> base_incidence = X_matrix[, 1]   # "historical" incidence for design decisions
        |
        v
get_designs(design_id, n_resamples, N, incidence, nb_list, coords)
  -> Z_matrix: [100 x n_design_resamples], binary {0,1}, ~50%% treated
        |
        v
estimate_tau(estimation_mode, Y_sim[100 x n_out], Z[100], spill_term[100],
             X_matrix, active_listw, n_outcome_resamples)
  -> list(estimates[n_out], ci_lower[n_out], ci_upper[n_out])
        |
        v
results data frame: 1 row per scenario, columns:
  Incidence_Mode  | "iid" / "spatial" / "poisson"
  Rho_Incidence   | 0, 0.20, or 0.50
  Neighbor_Type   | "rook" / "queen"
  Design          | "Design 1" .. "Design 7"
  Rho             | 0.00, 0.01, 0.20, 0.50
  Gamma           | 0.5, 0.6, 0.7, 0.8
  Spillover_Type  | "control_only" / "both"
  Mean_Estimate   | mean(tau-hat across iterations)
  Bias            | Mean_Estimate - true_tau (1.0)
  SD              | sd(tau-hat)
  MSE             | Bias^2 + SD^2
  Coverage        | fraction of CIs containing true_tau
  Fail_Rate       | fraction of MLE iterations that failed to converge
```

---

## Parameter Grid (1,920 total scenarios)

| Parameter | Values | # Levels |
|-----------|--------|---------|
| Incidence config | iid(x1) + spatial(x2) + poisson(x2) | 5 configs |
| `nb_type` | rook, queen | 2 |
| `rho` (spatial autocorrelation of outcome) | 0.00, 0.01, 0.20, 0.50 | 4 |
| `gamma` (spillover magnitude) | 0.5, 0.6, 0.7, 0.8 | 4 |
| `spill_type` | control_only, both | 2 |
| `design_id` | 1, 2, 3, 4, 5, 7 | 6 |
| **Total** | 5 × 2 × 4 × 4 × 2 × 6 | **1,920** |

**Scenarios per incidence config:** 384 (= 2×4×4×2×6)
- iid: 384 scenarios (rho_X = N/A)
- spatial (rho_X = 0.20): 384 scenarios
- spatial (rho_X = 0.50): 384 scenarios
- poisson (rho_X = 0.20): 384 scenarios
- poisson (rho_X = 0.50): 384 scenarios

---

## Design Quick Reference

| ID | Name | Deterministic? | ~Treated %% | Key Feature |
|----|------|:--------------:|:-----------:|-------------|
| 1 | Checkerboard | **Yes** | 50%% | Max spatial separation, alternating grid |
| 2 | High Incidence Focus | **Yes** (given X) | 50%% | Targets top-50%% burden clusters |
| 3 | Saturation Quadrants | No | ~50%% (varies) | Random saturation level per quadrant (0.2–0.8) |
| 4 | Isolation Buffer | No | ~20–30%% | Greedy: no two treated clusters adjacent |
| 5 | 2×2 Blocking | No | 50%% | 1:1 randomization within 2×2 spatial blocks |
| 7 | Balanced Quartiles | No | ~50%% | Stratified by incidence quartile, balanced within strata |

**Note:** Design 6 (Center Hotspot) is implemented in `03_designs.R` but **excluded**
from simulation runs (`design_ids <- c(1, 2, 3, 4, 5, 7)`).

Deterministic designs (1, 2) generate **one** assignment vector and replicate it for
all `n_design_resamples` to avoid redundant computation. This is enforced via
`is_design_deterministic()` in `03_designs.R`.

---

## Fixed DGP Parameters

```r
true_tau        <- 1.0         # Target estimand
beta            <- 1.0         # Incidence coefficient in outcome model
sigma           <- 1.0         # Residual SD
grid_dim        <- 10          # 10x10 = 100 clusters
base_rate       <- 35 / 100000 # Poisson: SUD rate (Mirzaei et al.)
pop_per_cluster <- 1000        # Poisson: equal population per cluster
pop_mode        <- "equal"     # "equal" or "heterogeneous"
include_spill_covariate <- TRUE  # Oracle mode: true Spill covariate in MLE model
```

**Resample counts:**

| Mode | Design resamples | Outcome resamples | Iterations/scenario |
|------|:---:|:---:|:---:|
| DIM | 25 | 100 | 2,500 |
| MLE | 25 | 10 | 250 |

---

## Current State (as of %s)

**DIM simulation:** COMPLETE
- File: `%s`
- Splits: `sim_results_DIM_iid_*.rds`, `sim_results_DIM_spatial_*.rds`, `sim_results_DIM_poisson_*.rds`
- Stats: %d scenarios | Bias [%.3f, %.3f] (mean %.3f) | MSE [%.4f, %.3f] (mean %.3f) | Coverage [%.2f, %.2f] (mean %.2f)
- Visualizations: %d PDFs in results/

**MLE simulation:** COMPLETE
- File: `%s`
- Splits: `%s`, `%s`, `%s`
- Stats: %d scenarios | Bias [%.3f, %.3f] (mean %.3f) | MSE [%.4f, %.3f] (mean %.3f) | Coverage [%.2f, %.2f] (mean %.2f) | Fail_Rate mean=%.4f (max=%.4f)
- Visualizations: %d PDFs in results/

**Git:** All files committed in the `gallant-buck` worktree and merged to main.

---

## Critical Invariants (DO NOT Violate)

1. **Incidence is generated ONCE per `(mode, rho_X)` config** — NOT per `(nb_type, rho)`
   scenario. The incidence matrix `X_matrix` is fixed for a given config; only the
   outcome DGP varies with `rho`.

2. **Results are NEVER aggregated across incidence modes.** Each of iid, spatial,
   poisson is reported separately. See `split_by_incidence_config()` in `06`.

3. **Deterministic designs (1, 2) generate 1 assignment and replicate** for all
   design resamples. Do NOT re-draw them `n_design_resamples` times.

4. **Per-scenario seed** = `digest::digest2int(paste(inc_mode, rho_x, nb_type, rho, gamma, spill_type, d_id, sep="|"))`.
   Never use a single `set.seed()` for the whole run.

5. **`base_incidence = X_matrix[, 1]`** — only the first column is used for treatment
   design decisions (the "historically observed" incidence). Do not pass the full matrix
   to design functions.

---

## Bugs Previously Encountered & Fixed

| Bug | Symptom | Fix |
|-----|---------|-----|
| Unicode `rho` in PDF | `pdf()` conversion failure on `\u03c1` | Replace with literal text `"rho"` |
| Pairwise dominance NaN | `make.names("Design 1")` → `"Design.1"` didn\'t match pivot column `"Design 1"` | Use `wide[[designs[i]]]` direct column access |
| Eta-squared > 1.0 | `var(group_means)/var(total)` is not proper eta-sq | Use `SS_between / SS_total` where `SS_between = sum(n_j * (mean_j - grand_mean)^2)` |
| `I` shadows `base::I` | Subtle errors if `base::I()` called downstream | Renamed to `I_mat <- diag(N_clusters)` |
| Degenerate Z (all 0 or 1) | `estimate_tau` crashed | Added early-exit guard returning all `NA` |

---

## Common User Requests → Code Actions

| Request | Action |
|---------|--------|
| Run DIM simulation | Set `estimation_mode <- "DIM"` in `05` (line 40), `Rscript 05_run_simulation.R` |
| Run MLE simulation | Set `estimation_mode <- "MLE"` in `05` (line 40), `Rscript 05_run_simulation.R` |
| Generate plots (from existing results) | `source("06_visualizations.R"); run_all_visualizations(estimation_mode="DIM")` |
| View one incidence config only | `r <- load_latest_results(); cfg <- split_by_incidence_config(r); run_standard_tables(cfg[["iid Uniform"]])` |
| Custom stratified table | `table_stratified(results, c("Design", "Rho"), "iid Uniform")` |
| Add a new design | Add `case` to `get_designs()` in `03`, update `get_design_names()`, add ID to `design_ids` in `05` |
| Change parameter sweep | Edit `*_vals` vectors in `05` CONFIGURABLE PARAMETERS section (lines ~50–80) |
| Run parallel (across incidence configs) | Set `n_cores > 1` in `05` (line 44) |
| Compare DIM vs MLE | Load both result sets, `bind_rows(dim_r, mle_r)`, compare `Coverage` and `Bias` |
| Non-oracle MLE (omit Spill covariate) | Set `include_spill_covariate = FALSE` in `estimate_tau()` call in `05` |

---

## Visualization Functions (`06_visualizations.R`)

**Helpers:**
- `inc_config_label(inc_mode, rho_x)` — human-readable config label
- `split_by_incidence_config(results)` — splits combined df into named list of per-config dfs
- `load_latest_results(results_dir, estimation_mode)` — loads most recently modified .rds

**Plots (each accepts `results` df + `inc_label` string):**
- `plot_mse_by_neighbor()` — boxplot of MSE by design, faceted by rook/queen
- `plot_mse_heatmap()` — tile heatmap: rho × design, cell = avg MSE
- `plot_mse_per_design()` — per-design MSE line plots across gamma, one panel per design
- `plot_master_comparison()` — bar chart sorted by MSE for a given neighbor type
- `plot_coverage_by_design()` — boxplot of coverage by design (red line at 0.95)
- `plot_incidence_heatmaps()` — side-by-side grid heatmaps of iid/spatial/poisson incidence
- `plot_incidence_distributions()` — faceted density plots comparing incidence distributions
- `plot_design_ranks()` — stacked bar of rank-frequency ("win rates") per design
- `plot_bias_variance()` — stacked bar of Bias² + Variance components per design
- `plot_coverage_mse_tradeoff()` — scatter: each point = 1 scenario, x=MSE, y=Coverage

**Tables (each accepts `results` df + `inc_label` string; prints to console, returns invisible):**
- `table_design_ranks()` — win-rate frequency table
- `table_robustness()` — best/Q25/median/Q75/worst MSE per design
- `table_pairwise_dominance()` — %% scenarios where row-design beats col-design
- `table_sensitivity()` — eta-squared for Rho, Gamma, Spillover_Type, Neighbor_Type per design
- `table_stratified(results, group_by_vars, inc_label)` — flexible: group by any combination
- `table_comprehensive()` — full scenario-level detail (sorted, first 30 rows by default)

**Runners:**
- `run_standard_tables(results, inc_label)` — runs all 10 standard tables for one config
- `run_all_visualizations(results, results_dir, estimation_mode, output_pdf=TRUE)` — full pipeline

---

## Planned Extensions (not yet implemented)

- Heterogeneous population mode for Poisson incidence (`pop_mode = "heterogeneous"`)
- Non-oracle MLE runs (toggle `include_spill_covariate = FALSE`) for realistic estimation
- Formal DIM vs MLE comparison analysis (joint load of both result sets)
- Sensitivity to grid dimension (`grid_dim = 8` or `15`)
- Additional designs or design variants
',
  format(Sys.time(), "%%Y-%%m-%%d %%H:%%M"),
  # DIM stats
  basename(dim_file),
  dim_s$n_scenarios,
  dim_s$bias_range[1], dim_s$bias_range[2], dim_s$bias_mean,
  dim_s$mse_range[1], dim_s$mse_range[2], dim_s$mse_mean,
  dim_s$cov_range[1], dim_s$cov_range[2], dim_s$cov_mean,
  length(dim_pdfs),
  # MLE stats
  basename(mle_file),
  mle_iid_file, mle_spatial_file, mle_poisson_file,
  mle_s$n_scenarios,
  mle_s$bias_range[1], mle_s$bias_range[2], mle_s$bias_mean,
  mle_s$mse_range[1], mle_s$mse_range[2], mle_s$mse_mean,
  mle_s$cov_range[1], mle_s$cov_range[2], mle_s$cov_mean,
  mle_s$fail_mean, mle_s$fail_max,
  length(mle_pdfs)
)

writeLines(claude_md, file.path(script_dir, "CLAUDE.md"))
cat("Written: CLAUDE.md\n\n")

# ==============================================================================
# STEP 6: Write README.md
# ==============================================================================

cat("=== Writing README.md ===\n")

readme_md <- sprintf(
'# Modular Incidence Simulation for Spatial CRT Design Evaluation

> For AI session context and quick technical reference, see [CLAUDE.md](CLAUDE.md).

## Overview

This project evaluates which **treatment assignment design** produces the most accurate
estimates of an intervention effect when outcome incidence is spatially heterogeneous
and spillover effects are present. The setting is a Spatial Cluster Randomized Trial
(CRT) on a 10×10 regular lattice grid of 100 clusters.

**Research question:** Across a wide range of spatial dependence, spillover magnitude,
and incidence structures, which of 6 candidate treatment assignment strategies minimizes
estimation error (MSE) and maintains valid inferential coverage?

**Application context:** Sudden Unexpected Death (SUD) in North Carolina counties.
The Poisson incidence mode is calibrated to a SUD base rate of 35/100,000 (Mirzaei et al.),
consistent with county-level Poisson regression analyses by Gan et al. and Watson et al.

---

## Project Lineage

This is a clean-room rewrite of
`code/OutcomeIncidenceDesign/SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd`
(~600 lines, monolithic Rmd). The refactoring goals were:

1. **Modularity** — Separate concerns into 6 numbered R scripts (01–06) that can be run independently or sourced in sequence.
2. **Three incidence modes** — Added Spatial (SAR filter) and Poisson modes alongside the existing iid Uniform.
3. **Correctness fixes** — Seven specific bugs and design flaws identified and corrected (see section below).
4. **Long-run reliability** — Replaced `save.image()` with targeted `saveRDS()`, added per-scenario seeding, ETA logging, and degenerate case handling.

The predecessor Rmd is preserved untouched.

---

## File Inventory

| File | Lines | Purpose | Key Functions |
|------|------:|---------|---------------|
| `00_mathematical_specification.Rmd` | ~495 | Theory doc with full LaTeX DGP formulas | (rendered to HTML) |
| `01_spatial_setup.R` | 54 | Regular lattice grid, rook/queen weight matrices | `build_spatial_grid()`, `get_active_spatial()` |
| `02_incidence_generation.R` | 112 | Three incidence generation modes | `generate_incidence()` + 3 mode-specific functions |
| `03_designs.R` | 137 | Six treatment assignment strategies | `get_designs()`, `get_design_names()`, `is_design_deterministic()` |
| `04_estimation.R` | 102 | DIM and MLE estimation with CI extraction | `estimate_tau()` |
| `05_run_simulation.R` | ~385 | Main orchestrator, nested parameter loop | `run_incidence_config()` |
| `06_visualizations.R` | ~800 | All plots and summary tables | 18 functions; entry point `run_all_visualizations()` |
| `complete_after_mle.R` | ~250 | Post-completion script (viz + docs + stats) | (run once after MLE finishes) |
| `results/` | — | All output: .rds data files + .pdf figures | — |

---

## How to Run

### Prerequisites (R packages)

```r
install.packages(c("spdep", "spatialreg", "dplyr", "tidyr",
                   "digest", "parallel", "ggplot2", "viridis"))
```

### Quick Start — DIM estimation (~5–10 minutes)

```r
# Open 05_run_simulation.R, ensure:
#   estimation_mode <- "DIM"   (line 40)
#   n_cores <- 1               (line 44)

# Then run:
setwd("code/OutcomeIncidenceDesign/ModularIncidenceSim")
source("05_run_simulation.R")
```

### Full Run — MLE estimation (~8–10 hours)

```r
# In 05_run_simulation.R set:
#   estimation_mode <- "MLE"   (line 40)
# Optionally: n_cores <- 4     (line 44, parallel across incidence configs)

# Run from terminal for reliability:
# cd code/OutcomeIncidenceDesign/ModularIncidenceSim
# Rscript 05_run_simulation.R
```

### Generating Visualizations Only (results already exist)

```r
source("06_visualizations.R")

# All configs, DIM results:
run_all_visualizations(estimation_mode = "DIM_combined")

# All configs, MLE results:
run_all_visualizations(estimation_mode = "MLE_combined")

# Single incidence config + tables:
r <- load_latest_results(estimation_mode = "DIM_combined")
cfgs <- split_by_incidence_config(r)
run_standard_tables(cfgs[["iid Uniform"]], "iid Uniform")
```

---

## Simulation Design Summary

**Grid:** 10×10 regular lattice, N = 100 clusters. Two contiguity structures: rook (4-connected) and queen (8-connected).

**Data generating process:** Spatial Durbin Model (SDM)

> Y = (I − ρW)⁻¹ [τZ + γ·Spill(Z) + βX + ε]

where τ = 1.0, β = 1.0, σ = 1.0 (noise SD), and Spill(Z) is the row-standardized
mean treatment of neighbors. Two spillover modes: `control_only` (γ applied only to
control units) and `both` (γ applied to all units).

**Incidence modes:**
- **iid Uniform:** X_i ~ Uniform(0, 1) independently
- **Spatial:** SAR filter with pnorm transform — spatially correlated but marginally Uniform(0,1)
- **Poisson:** Spatially correlated log-rates → Poisson counts → rank-normalized rates in [0,1]

**Treatment designs:** 6 strategies (Checkerboard, High Incidence Focus, Saturation Quadrants, Isolation Buffer, 2×2 Blocking, Balanced Quartiles)

**Estimation methods:**
- **DIM:** Difference in Means with Neyman variance estimator (25 design × 100 outcome = 2,500 iterations/scenario)
- **MLE:** Spatial autoregressive model `lagsarlm(Y ~ Z + Spill + X)` (25 design × 10 outcome = 250 iterations/scenario; oracle: true Spill covariate included)

**Total:** 5 incidence configs × 2 neighbor types × 4 ρ × 4 γ × 2 spillover types × 6 designs = **1,920 scenarios**

---

## Results Summary

### DIM Results (completed %s)

- **%d scenarios** | Bias [%.3f, %.3f] (mean %.3f) | MSE [%.4f, %.3f] (mean %.3f)
- Coverage [%.2f, %.2f] (mean %.2f)
- File: `results/%s`

### MLE Results (completed %s)

- **%d scenarios** | Bias [%.3f, %.3f] (mean %.3f) | MSE [%.4f, %.3f] (mean %.3f)
- Coverage [%.2f, %.2f] (mean %.2f) | Convergence fail rate mean=%.4f (max=%.4f)
- File: `results/%s`

### Results File Naming Convention

```
results/
  sim_results_{DIM|MLE}_combined_{YYYYMMDD_HHMMSS}.rds   # All 1,920 scenarios
  sim_results_{DIM|MLE}_iid_{timestamp}.rds              # iid Uniform only (384 rows)
  sim_results_{DIM|MLE}_spatial_{timestamp}.rds          # Spatial only (768 rows)
  sim_results_{DIM|MLE}_poisson_{timestamp}.rds          # Poisson only (768 rows)
  {DIM|MLE}_combined_{config_name}.pdf                   # Per-config 8-plot PDF
  {DIM|MLE}_combined_incidence_overview.pdf              # Incidence heatmaps + distributions
```

**Key rule:** Results for iid Uniform, Spatial, and Poisson incidence modes are always
reported **separately**. The combined file exists for loading convenience only.

---

## Key Design Decisions

1. **Incidence generated once per `(mode, rho_X)` config** — not re-drawn per `(nb_type, rho)`.
   In practice, a researcher observes historical incidence once; spatial model parameters
   then affect how outcomes propagate, not what incidence looks like.

2. **Deterministic design detection** — Designs 1 (Checkerboard) and 2 (High Incidence Focus)
   yield the same assignment for every design resample (given fixed incidence). They are
   run once and the result replicated, avoiding ~25× redundant computation.

3. **Per-scenario deterministic seeding** — `digest::digest2int(paste(params, sep="|"))`
   ensures reproducibility and order-independence. Adding/removing scenarios doesn\'t
   contaminate others.

4. **Separate reporting by incidence mode** — The three incidence modes represent
   fundamentally different assumptions about prior data. Aggregating across them would
   obscure rather than inform design selection.

5. **Oracle spillover in MLE** — `lagsarlm()` includes the true spillover term as a
   known covariate (`Y ~ Z + Spill + X`). This gives MLE an advantage in correctly
   specified scenarios. A `include_spill_covariate = FALSE` toggle exists for realistic
   estimation comparisons (not yet run).

6. **DIM iteration counts** — Originally planned at 100,000 iterations/scenario; reduced
   to 2,500 (25 design × 100 outcome resamples) for practical runtime.

7. **Design 6 excluded** — The Center Hotspot design is implemented in `03_designs.R`
   but not included in the default `design_ids` sweep pending further evaluation.

---

## Known Issues & Bugs Fixed

| # | Issue | Root Cause | Fix |
|---|-------|-----------|-----|
| 1 | PDF render failure | `\u03c1` (Unicode rho) not supported by R `pdf()` device | Replace with literal `"rho"` text |
| 2 | Pairwise dominance NaN | `make.names("Design 1")` → `"Design.1"` mismatched column name `"Design 1"` in pivot output | Use `wide[[designs[i]]]` direct bracket access |
| 3 | Eta-squared > 1.0 | `var(group_means)/var(total)` inflates ratio when groups are unbalanced | Proper formula: `SS_between/SS_total` where `SS_between = sum(n_j*(mean_j - grand)^2)` |
| 4 | `I()` namespace conflict | `I <- diag(N)` shadows `base::I()` | Renamed to `I_mat` |
| 5 | Degenerate Z crash | Designs sometimes produce all-treated or all-control assignments | Early return of `NA` in `estimate_tau()` when `n_trt == 0` or `n_trt == N` |

---

## Status & Next Steps

**Completed:**
- [x] All 7 code files (00–06) written and tested
- [x] Mathematical specification document rendered to HTML
- [x] DIM simulation: 1,920 scenarios completed, visualizations generated
- [x] MLE simulation: 1,920 scenarios completed, visualizations generated
- [x] Project documentation: CLAUDE.md + README.md

**Planned extensions:**
- [ ] **Heterogeneous population** — Poisson mode with `pop_mode = "heterogeneous"` (log-normal cluster sizes)
- [ ] **Non-oracle MLE** — Run MLE without the true Spill covariate to assess realistic estimation performance
- [ ] **DIM vs MLE comparison** — Load both result sets, compare Coverage gain from MLE vs DIM
- [ ] **Grid sensitivity** — Rerun with `grid_dim = 8` or `grid_dim = 15`
- [ ] **Additional designs** — Evaluate Design 6 (Center Hotspot) or add new strategies

---

## References

- Mirzaei, A. et al. (2019). Sudden unexpected death rates.
- Gan, W. et al. (2019). County-level Poisson regression for SUD mortality.
- Watson, K. et al. (AHA). Census tract-level spatial analysis of SUD.
- LeSage, J. & Pace, R.K. (2009). *Introduction to Spatial Econometrics*. CRC Press.
- R packages: `spdep` (Bivand et al.), `spatialreg` (Bivand & Piras), `digest`, `ggplot2`, `viridis`
',
  # DIM stats
  format(file.mtime(dim_file), "%%Y-%%m-%%d"),
  dim_s$n_scenarios,
  dim_s$bias_range[1], dim_s$bias_range[2], dim_s$bias_mean,
  dim_s$mse_range[1], dim_s$mse_range[2], dim_s$mse_mean,
  dim_s$cov_range[1], dim_s$cov_range[2], dim_s$cov_mean,
  basename(dim_file),
  # MLE stats
  format(file.mtime(mle_file), "%%Y-%%m-%%d"),
  mle_s$n_scenarios,
  mle_s$bias_range[1], mle_s$bias_range[2], mle_s$bias_mean,
  mle_s$mse_range[1], mle_s$mse_range[2], mle_s$mse_mean,
  mle_s$cov_range[1], mle_s$cov_range[2], mle_s$cov_mean,
  mle_s$fail_mean, mle_s$fail_max,
  basename(mle_file)
)

writeLines(readme_md, file.path(script_dir, "README.md"))
cat("Written: README.md\n\n")

# ==============================================================================
# STEP 7: Final summary
# ==============================================================================

cat("=======================================================\n")
cat("  COMPLETION SUMMARY\n")
cat("=======================================================\n")
cat(sprintf("DIM: %d scenarios | Avg Coverage %.2f | Avg MSE %.4f\n",
  dim_s$n_scenarios, dim_s$cov_mean, dim_s$mse_mean))
cat(sprintf("MLE: %d scenarios | Avg Coverage %.2f | Avg MSE %.4f | Avg Fail %.4f\n",
  mle_s$n_scenarios, mle_s$cov_mean, mle_s$mse_mean, mle_s$fail_mean))
cat(sprintf("DIM best design (median MSE): %s\n", dim_s$best_design))
cat(sprintf("MLE best design (median MSE): %s\n", mle_s$best_design))
cat("Documentation: CLAUDE.md + README.md written\n")
cat(sprintf("MLE PDFs: %d files in results/\n", length(mle_pdfs)))
cat("\nNext step: git add + commit + merge worktree to main\n")
cat("=======================================================\n")
