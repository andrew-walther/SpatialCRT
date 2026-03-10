# CLAUDE.md — AI Session Context for IncidenceDesign

> **Read this file first** when starting a new Claude session on this project.
> Provides ~80% of the context needed to be immediately productive without
> re-reading all 1,589+ lines of code, result files, or figure PDFs.
> For human-readable narrative documentation see [README.md](README.md).
> For cross-project context, see [../../CLAUDE.md](../../CLAUDE.md).

---

## What This Project Is

A modular simulation study evaluating **6 treatment assignment designs** for Spatial
Cluster Randomized Trials (CRTs) under heterogeneous outcome incidence and spatial
spillover. The estimand is **tau = 1.0** (direct treatment effect). Two estimators
(DIM, MLE) are run across 1,920 parameter scenarios with 3 incidence generation
modes that are always reported **separately** — never aggregated.

**Application context:** Sudden Unexpected Death (SUD) in NC counties.
Poisson base rate 35/100,000 (Mirzaei et al.).

**Metrics tracked per scenario:** Bias, SD, MSE, Coverage (95% CI), Fail_Rate.

---

## Research Focus and Framing (IMPORTANT)

**Primary goal: compare the 6 treatment sampling designs** under realistic spatial
conditions. The key question is which design minimizes estimation error and maintains
valid coverage — NOT which estimator is better.

**MLE is the primary estimator.** It is the methodologically appropriate choice
because spatial dependence ($\rho$) and spillover ($\gamma$) are both present in
the outcome DGP. DIM ignores both and produces systematically poor coverage (~72%).

**DIM serves as a proof-of-concept / naive baseline only.** Use DIM results to
sanity-check simulation mechanics and as a secondary comparator. All substantive
design recommendations should be based on MLE results.

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
| `06_visualizations.R` | ~800 | Plots + tables | 18 functions — see section below |
| `complete_after_mle.R` | ~250 | Post-MLE script | Runs viz, writes docs, prints stats |

---

## Architecture

**Pipeline:** `01 -> 02 -> 03 -> 04`, orchestrated by `05`, visualized by `06`.

- `05_run_simulation.R` sources `01`-`04` at runtime via `file.path(script_dir, "0X_*.R")`
- `06_visualizations.R` sources `01`-`03` (needs grid/incidence/design helpers)
- Both detect working directory via `dirname(sys.frame(1)$ofile)` with `tryCatch` fallback to `getwd()`

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
  -> Z_matrix: [100 x n_design_resamples], binary {0,1}, ~50% treated
        |
        v
estimate_tau(estimation_mode, Y_sim[100 x n_out], Z[100], spill_term[100],
             X_matrix, active_listw, n_outcome_resamples)
  -> list(estimates[n_out], ci_lower[n_out], ci_upper[n_out])
        |
        v
results data frame — 1 row per scenario:
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

| Parameter | Values | Levels |
|-----------|--------|:------:|
| Incidence config | iid (x1) + spatial (x2) + poisson (x2) | 5 |
| `nb_type` | rook, queen | 2 |
| `rho` (outcome spatial autocorrelation) | 0.00, 0.01, 0.20, 0.50 | 4 |
| `gamma` (spillover magnitude) | 0.5, 0.6, 0.7, 0.8 | 4 |
| `spill_type` | control_only, both | 2 |
| `design_id` | 1, 2, 3, 4, 5, 7 | 6 |
| **Total** | 5 x 2 x 4 x 4 x 2 x 6 | **1,920** |

Scenarios per incidence config: 384 (= 2x4x4x2x6)

---

## Design Quick Reference

| ID | Name | Deterministic? | ~Treated% | Key Feature |
|----|------|:-:|:-:|-------------|
| 1 | Checkerboard | **Yes** | 50% | Max spatial separation, alternating grid |
| 2 | High Incidence Focus | **Yes** (given X) | 50% | Targets top-50% burden clusters |
| 3 | Saturation Quadrants | No | ~50% (varies) | Random saturation per quadrant (0.2-0.8) |
| 4 | Isolation Buffer | No | ~20-30% | Greedy: no two treated clusters adjacent |
| 5 | 2x2 Blocking | No | 50% | 1:1 randomization within 2x2 spatial blocks |
| 7 | Balanced Quartiles | No | ~50% | Stratified by incidence quartile |

**Note:** Design 6 (Center Hotspot) is implemented in `03_designs.R` but excluded
from runs (`design_ids <- c(1, 2, 3, 4, 5, 7)`).

Deterministic designs (1, 2) generate **one** assignment and replicate it across all
`n_design_resamples` — enforced via `is_design_deterministic()` in `03_designs.R`.

---

## Fixed DGP Parameters

```r
true_tau        <- 1.0          # Target estimand
beta            <- 1.0          # Incidence coefficient in outcome model
sigma           <- 1.0          # Residual SD
grid_dim        <- 10           # 10x10 = 100 clusters
base_rate       <- 35 / 100000  # Poisson: SUD rate (Mirzaei et al.)
pop_per_cluster <- 1000         # Poisson: equal population per cluster
pop_mode        <- "equal"      # "equal" or "heterogeneous"
include_spill_covariate <- TRUE # Oracle mode: true Spill covariate in MLE
```

**Resample counts:**

| Mode | Design resamples | Outcome resamples | Iterations/scenario |
|------|:---:|:---:|:---:|
| DIM | 25 | 100 | 2,500 |
| MLE | 25 | 10 | 250 |

---

## Current State (as of 2026-03-05)

**DIM simulation:** COMPLETE
- File: `results/sim_results_DIM_combined_20260304_195321.rds`
- Splits: `sim_results_DIM_iid_*.rds`, `sim_results_DIM_spatial_*.rds`, `sim_results_DIM_poisson_*.rds`
- Stats: 1,920 scenarios | Bias [-0.883, 0.558] (mean -0.200) | MSE [0.032, 0.799] (mean 0.179) | Coverage [0.00, 0.99] (mean 0.72)
- Visualizations: 6 PDFs in results/ (5 per-config + overview)

**MLE simulation:** COMPLETE
- File: `results/sim_results_MLE_combined_20260305_150742.rds`
- Splits: `sim_results_MLE_iid_*.rds`, `sim_results_MLE_spatial_*.rds`, `sim_results_MLE_poisson_*.rds`
- Stats: 1,920 scenarios | Bias [-0.895, 0.454] | MSE [0.009, 3.910] | Coverage [0.00, 1.00] | Fail_Rate = 0.0 (zero convergence failures)
- Note: MLE MSE max (3.910) exceeds DIM max (0.799) — occurs in high-rho scenarios where spatial lag model variance inflates with only 250 iterations
- Visualizations: 6 PDFs in results/ (5 per-config + overview)

**Git:** Original work on branch `claude/gallant-buck`, merged to `main` 2026-03-05.
Reorganized into `projects/IncidenceDesign/` on branch `claude/dreamy-wiles`.

---

## Critical Invariants (DO NOT Violate)

1. **Incidence generated ONCE per `(mode, rho_X)` config** — NOT per `(nb_type, rho)` scenario.
   `X_matrix` is fixed for a given config; only the outcome DGP varies with `rho`.

2. **Results are NEVER aggregated across incidence modes.** Each of iid, spatial, poisson
   is always reported separately. See `split_by_incidence_config()` in `06_visualizations.R`.

3. **Deterministic designs (1, 2) generate 1 assignment and replicate** for all design
   resamples. Never re-draw them `n_design_resamples` times.

4. **Per-scenario seed** = `digest::digest2int(paste(inc_mode, rho_x, nb_type, rho, gamma, spill_type, d_id, sep="|"))`.
   Never use a single `set.seed()` for the whole run.

5. **`base_incidence = X_matrix[, 1]`** — only the first column is used for treatment
   design decisions (the "historically observed" incidence). Do not pass the full matrix.

---

## Bugs Previously Encountered & Fixed

| Bug | Symptom | Fix |
|-----|---------|-----|
| Unicode rho in PDF | `pdf()` conversion failure on `\u03c1` | Use literal text `"rho"` |
| Pairwise dominance NaN | `make.names("Design 1")` -> `"Design.1"` mismatched pivot column | Use `wide[[designs[i]]]` direct bracket access |
| Eta-squared > 1.0 | `var(group_means)/var(total)` incorrect | Use `SS_between/SS_total` where `SS_between = sum(n_j*(mean_j - grand_mean)^2)` |
| `I` shadows `base::I` | Subtle namespace errors | Renamed to `I_mat <- diag(N_clusters)` |
| Degenerate Z (all 0 or 1) | `estimate_tau()` crash | Early-exit guard returning all `NA` |
| R sprintf 10k char limit | `complete_after_mle.R` failed writing docs | Write docs directly via Claude Code `Write` tool instead |

---

## Common User Requests -> Code Actions

| Request | Action |
|---------|--------|
| Run DIM simulation | Set `estimation_mode <- "DIM"` in `05` line 40, `Rscript 05_run_simulation.R` |
| Run MLE simulation | Set `estimation_mode <- "MLE"` in `05` line 40, `Rscript 05_run_simulation.R` |
| Generate plots only | `source("06_visualizations.R"); run_all_visualizations(estimation_mode="DIM_combined")` |
| One incidence config | `r <- load_latest_results(); cfg <- split_by_incidence_config(r); run_standard_tables(cfg[["iid Uniform"]])` |
| Custom stratified table | `table_stratified(results, c("Design", "Rho"), "iid Uniform")` |
| Add a new design | Add `case` to `get_designs()` in `03`, update `get_design_names()`, add ID to `design_ids` in `05` |
| Change parameter sweep | Edit `*_vals` vectors in `05` CONFIGURABLE PARAMETERS section (lines ~50-80) |
| Run parallel | Set `n_cores > 1` in `05` line 44 |
| Compare DIM vs MLE | Load both result sets, join on scenario keys, compare Coverage and Bias columns |
| Non-oracle MLE | Set `include_spill_covariate = FALSE` in `estimate_tau()` call in `05` |

---

## Visualization Functions (`06_visualizations.R`)

**Helpers:**
- `inc_config_label(inc_mode, rho_x)` — human-readable config label string
- `split_by_incidence_config(results)` — splits combined df into named list of per-config dfs
- `load_latest_results(results_dir, estimation_mode)` — loads most recently modified .rds

**Plots** (each accepts `results` df + `inc_label` string):
- `plot_mse_by_neighbor()` — boxplot of MSE by design, faceted by rook/queen
- `plot_mse_heatmap()` — tile heatmap: rho x design, cell = avg MSE
- `plot_mse_per_design()` — per-design MSE line plots across gamma, one panel per design
- `plot_master_comparison()` — bar chart sorted by MSE for a given neighbor type
- `plot_coverage_by_design()` — boxplot of coverage by design (red dashed line at 0.95)
- `plot_incidence_heatmaps()` — side-by-side grid heatmaps of iid/spatial/poisson incidence
- `plot_incidence_distributions()` — faceted density plots comparing incidence distributions
- `plot_design_ranks()` — stacked bar of rank-frequency ("win rates") per design
- `plot_bias_variance()` — stacked bar of Bias-squared + Variance components per design
- `plot_coverage_mse_tradeoff()` — scatter: each point = 1 scenario, x=MSE, y=Coverage

**Tables** (each prints to console, returns invisible):
- `table_design_ranks()` — win-rate frequency table
- `table_robustness()` — best/Q25/median/Q75/worst MSE per design
- `table_pairwise_dominance()` — % of scenarios where row-design beats col-design
- `table_sensitivity()` — eta-squared for Rho, Gamma, Spillover_Type, Neighbor_Type per design
- `table_stratified(results, group_by_vars, inc_label)` — flexible grouping by any combination
- `table_comprehensive()` — full scenario-level detail (first 30 rows by default)

**Runners:**
- `run_standard_tables(results, inc_label)` — all 10 standard tables for one config
- `run_all_visualizations(results, results_dir, estimation_mode, output_pdf=TRUE)` — full pipeline

---

## Planned Extensions (not yet implemented)

- Heterogeneous population mode for Poisson (`pop_mode = "heterogeneous"`)
- Non-oracle MLE runs (`include_spill_covariate = FALSE`) for realistic estimation
- Formal DIM vs MLE comparison analysis (joint load of both result sets)
- Sensitivity to grid dimension (`grid_dim = 8` or `15`)
- Additional designs or design variants

---

## Relationship to SpillSpatialDepSim

IncidenceDesign extends the applied simulation framework from `projects/SpillSpatialDepSim/`.
SpillSpatialDepSim used a small grid (8–12 districts) with SAR estimation to evaluate
block vs. random assignment in an NC DOC applied context. IncidenceDesign asks the same
core question but at larger scale (100 clusters) with systematic design variation and
heterogeneous incidence modes.

Cross-reference SpillSpatialDepSim results:
```r
here::here("projects", "SpillSpatialDepSim", "results")
# or relative: ../../SpillSpatialDepSim/results/
```
