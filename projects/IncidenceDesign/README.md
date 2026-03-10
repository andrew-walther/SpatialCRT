# Modular Incidence Simulation for Spatial CRT Design Evaluation

> For AI session context and quick technical reference, see [CLAUDE.md](CLAUDE.md).

## Overview

This project evaluates which **treatment assignment design** produces the most accurate
estimates of an intervention effect when outcome incidence is spatially heterogeneous
and spillover effects are present. The setting is a Spatial Cluster Randomized Trial
(CRT) on a 10x10 regular lattice grid of 100 clusters.

**Research question:** Across a wide range of spatial dependence, spillover magnitude,
and incidence structures, which of 6 candidate treatment assignment strategies minimizes
estimation error (MSE) and maintains valid inferential coverage?

**Application context:** Sudden Unexpected Death (SUD) in North Carolina counties.
The Poisson incidence mode uses a SUD base rate of 35/100,000 (Mirzaei et al.),
consistent with county-level Poisson analyses by Gan et al. and Watson et al.

---

## Project Lineage

This is a clean-room rewrite of
`archive/OutcomeIncidenceDesign_Legacy/SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd`
(~600 lines, monolithic Rmd). Refactoring goals:

1. **Modularity** — 6 numbered R scripts (01-06) that can be run independently or sourced in sequence
2. **Three incidence modes** — Added Spatial (SAR filter) and Poisson modes alongside iid Uniform
3. **Correctness fixes** — Seven specific bugs and design flaws identified and corrected (see below)
4. **Long-run reliability** — Replaced `save.image()` with targeted `saveRDS()`, added per-scenario seeding, ETA logging, and degenerate case handling

The predecessor Rmd is preserved untouched.

---

## File Inventory

| File | Lines | Purpose | Key Functions |
|------|------:|---------|---------------|
| `00_mathematical_specification.Rmd` | ~495 | Theory doc with full LaTeX DGP formulas | (rendered to HTML) |
| `01_spatial_setup.R` | 54 | Regular lattice grid, rook/queen weight matrices | `build_spatial_grid()`, `get_active_spatial()` |
| `02_incidence_generation.R` | 112 | Three incidence generation modes | `generate_incidence()` + 3 mode functions |
| `03_designs.R` | 137 | Six treatment assignment strategies | `get_designs()`, `get_design_names()`, `is_design_deterministic()` |
| `04_estimation.R` | 102 | DIM and MLE estimation with CI extraction | `estimate_tau()` |
| `05_run_simulation.R` | ~385 | Main orchestrator, nested parameter loop | `run_incidence_config()` |
| `06_visualizations.R` | ~800 | All plots and summary tables | 18 functions; entry point `run_all_visualizations()` |
| `complete_after_mle.R` | ~250 | Post-completion script (viz + docs + stats) | Run once after MLE finishes |
| `results/` | — | All output: .rds data files + .pdf figures | — |

---

## How to Run

### Prerequisites

```r
install.packages(c("spdep", "spatialreg", "dplyr", "tidyr",
                   "digest", "parallel", "ggplot2", "viridis"))
```

### Quick Start — DIM estimation (~5-10 minutes)

```r
# In 05_run_simulation.R, ensure line 40 reads:
#   estimation_mode <- "DIM"

setwd("projects/IncidenceDesign/code")
source("05_run_simulation.R")
```

### Full Run — MLE estimation (~10 hours)

```bash
# Run from terminal to keep process independent of IDE session
# Use caffeinate to prevent sleep on macOS:
#   caffeinate -i -w $(pgrep -f 05_run_simulation) &
cd projects/IncidenceDesign/code
Rscript 05_run_simulation.R
```

In `05_run_simulation.R`, set `estimation_mode <- "MLE"` (line 40).
Optionally set `n_cores > 1` (line 44) for parallel execution across incidence configs.

### Generating Visualizations Only

```r
source("06_visualizations.R")

# All configs, most recent DIM results:
run_all_visualizations(estimation_mode = "DIM_combined")

# All configs, most recent MLE results:
run_all_visualizations(estimation_mode = "MLE_combined")

# Single incidence config with all tables:
r <- load_latest_results(estimation_mode = "DIM_combined")
cfgs <- split_by_incidence_config(r)
run_standard_tables(cfgs[["iid Uniform"]], "iid Uniform")
```

---

## Simulation Design Summary

**Grid:** 10x10 regular lattice, N = 100 clusters. Two contiguity structures:
rook (4-connected) and queen (8-connected).

**Data Generating Process — Spatial Durbin Model (SDM):**

```
Y = (I - rho*W)^{-1} * [tau*Z + gamma*Spill(Z) + beta*X + epsilon]
```

where tau = 1.0, beta = 1.0, sigma = 1.0, and `Spill(Z)` is the row-standardized
mean treatment of neighbors. Two spillover modes: `control_only` (gamma applied only
to control units) and `both` (applied to all units).

**Incidence modes:**
- **iid Uniform:** X_i ~ Uniform(0,1) independently
- **Spatial:** SAR filter + pnorm transform — spatially correlated but marginally Uniform(0,1)
- **Poisson:** Spatially correlated log-rates -> Poisson counts -> rank-normalized rates in [0,1]

**Treatment designs:** 6 strategies — Checkerboard, High Incidence Focus, Saturation
Quadrants, Isolation Buffer, 2x2 Blocking, Balanced Quartiles.

**Estimation methods:**
- **DIM:** Difference in Means with Neyman variance estimator (2,500 iterations/scenario)
- **MLE:** Spatial autoregressive model `lagsarlm(Y ~ Z + Spill + X)` using oracle Spill covariate (250 iterations/scenario)

**Total scenarios:** 5 incidence configs x 2 neighbor types x 4 rho x 4 gamma x 2 spillover types x 6 designs = **1,920 scenarios**

---

## Results Summary

### DIM (completed 2026-03-04)

- **1,920 scenarios** | Bias [-0.883, 0.558] (mean -0.200) | MSE [0.032, 0.799] (mean 0.179)
- Coverage [0.00, 0.99] (mean 0.72)
- File: `results/sim_results_DIM_combined_20260304_195321.rds`

### MLE (completed 2026-03-05)

- **1,920 scenarios** | Bias [-0.895, 0.454] | MSE [0.009, 3.910] | Coverage [0.00, 1.00]
- Convergence fail rate = 0.0 (zero failures across all 480,000 lagsarlm calls)
- Note: MLE max MSE (3.910) exceeds DIM (0.799) — occurs in high-rho scenarios where spatial lag model variance inflates with 250 iterations; median MSE is lower for MLE
- File: `results/sim_results_MLE_combined_20260305_150742.rds`

### Results File Naming Convention

```
results/
  sim_results_{DIM|MLE}_combined_{YYYYMMDD_HHMMSS}.rds  # All 1,920 scenarios
  sim_results_{DIM|MLE}_iid_{timestamp}.rds             # iid Uniform only (384 rows)
  sim_results_{DIM|MLE}_spatial_{timestamp}.rds         # Spatial only (768 rows)
  sim_results_{DIM|MLE}_poisson_{timestamp}.rds         # Poisson only (768 rows)
  {DIM|MLE}_combined_{config_name}.pdf                  # Per-config 8-plot PDF
  {DIM|MLE}_combined_incidence_overview.pdf             # Incidence heatmaps + distributions
```

**Key rule:** Results for iid Uniform, Spatial, and Poisson are always reported
**separately** — never aggregated. The combined .rds exists for loading convenience only.

---

## Key Design Decisions

1. **Incidence generated once per `(mode, rho_X)` config** — not re-drawn per `(nb_type, rho)`.
   In practice a researcher observes historical incidence once; spatial model parameters
   affect outcome propagation, not incidence itself.

2. **Deterministic design detection** — Designs 1 and 2 yield the same assignment for
   every design resample (given fixed incidence). They run once and replicate, saving ~25x compute.

3. **Per-scenario deterministic seeding** — `digest::digest2int(paste(params, sep="|"))` ensures
   reproducibility and order-independence. Adding/removing scenarios doesn't affect others.

4. **Separate reporting by incidence mode** — The three modes represent fundamentally
   different assumptions about prior data. Aggregating across them obscures rather than informs.

5. **Oracle spillover in MLE** — `lagsarlm()` includes the true Spill covariate.
   A `include_spill_covariate = FALSE` toggle exists in `estimate_tau()` for realistic
   estimation comparisons (not yet run systematically).

6. **DIM iteration counts** — Originally planned at 100,000 iterations/scenario; reduced
   to 2,500 (25 design x 100 outcome resamples) for practical runtime.

7. **Design 6 excluded** — Center Hotspot is implemented in `03_designs.R` but not
   included in the default `design_ids` sweep pending further evaluation.

---

## Known Issues & Bugs Fixed

| # | Issue | Root Cause | Fix |
|---|-------|-----------|-----|
| 1 | PDF render failure | `\u03c1` (Unicode rho) unsupported by `pdf()` device | Replace with literal `"rho"` text |
| 2 | Pairwise dominance NaN | `make.names()` mangled column names in `pivot_wider` output | Use `wide[[designs[i]]]` direct bracket access |
| 3 | Eta-squared > 1.0 | `var(group_means)/var(total)` incorrect formula | Proper: `SS_between/SS_total` |
| 4 | `I()` namespace conflict | `I <- diag(N)` shadowed `base::I()` | Renamed to `I_mat` |
| 5 | Degenerate Z crash | Designs occasionally produce all-treated or all-control | Early `NA` return in `estimate_tau()` |
| 6 | R sprintf 10k limit | `complete_after_mle.R` format string too long | Write docs via Claude Code `Write` tool |

---

## Status & Next Steps

**Completed:**
- [x] All 7 code files (00-06) written and tested
- [x] Mathematical specification document (00) rendered to HTML
- [x] DIM simulation: 1,920 scenarios, visualizations generated
- [x] MLE simulation: 1,920 scenarios, visualizations generated, zero convergence failures
- [x] Project documentation: CLAUDE.md + README.md

**Planned extensions:**
- [ ] Heterogeneous population — Poisson mode with `pop_mode = "heterogeneous"`
- [ ] Non-oracle MLE — run without the true Spill covariate for realistic comparison
- [ ] DIM vs MLE comparison — joint analysis of both result sets
- [ ] Grid sensitivity — rerun with `grid_dim = 8` or `grid_dim = 15`
- [ ] Additional designs — evaluate Design 6 (Center Hotspot) or add new strategies

---

## Relationship to SpillSpatialDepSim

IncidenceDesign extends the simulation framework from `projects/SpillSpatialDepSim/`.
SpillSpatialDepSim established the SAR model, spillover mechanics, and block stratification
approach in an applied NC DOC context (8–12 districts). IncidenceDesign asks the same
core design question at larger scale (100 clusters) with systematic variation across
incidence modes and formal design strategies.

---

## References

- Mirzaei, A. et al. (2019). Sudden unexpected death rates.
- Gan, W. et al. (2019). County-level Poisson regression for SUD mortality.
- Watson, K. et al. (AHA). Census tract-level spatial analysis of SUD.
- LeSage, J. & Pace, R.K. (2009). *Introduction to Spatial Econometrics*. CRC Press.
- R packages: `spdep` (Bivand et al.), `spatialreg` (Bivand & Piras), `digest`, `ggplot2`, `viridis`
