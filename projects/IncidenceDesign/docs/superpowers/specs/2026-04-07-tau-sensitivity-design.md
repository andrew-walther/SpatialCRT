# Combined Simulation Extension — Design Spec
## Tau Sensitivity + Monte Carlo SEs + Statistical Power
**Date:** 2026-04-07  
**Project:** IncidenceDesign (SpatialCRT)  
**Status:** Approved — implement as a single combined run

---

## Context

All 2,560 current simulation scenarios fix `true_tau = 1.0` and do not record Monte Carlo standard errors or statistical power. This spec bundles three related improvements into one implementation pass and one Longleaf run:

1. **Tau sensitivity** — vary τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}, expanding to 12,800 scenarios
2. **Simulation SEs** — store `N_Valid_Est` to enable exact SE computation for Bias, MSE, and Coverage
3. **Statistical Power** — record P(reject H₀: τ=0) per scenario alongside Coverage

All three require the same simulation-side changes (new columns in the results data frame). Bundling into one run avoids a separate re-run for SEs/Power on the existing 2,560-scenario grid.

Primary questions this extension answers:
1. Are the D3/D8 design rankings robust across effect sizes, or do they shift at small/large τ?
2. Which design achieves 80% power at the smallest detectable effect size?
3. Does MLE remain well-calibrated (coverage ≈ 0.95) at extreme τ?
4. How precise are the MSE/coverage estimates (i.e., are the Design 8 vs. Design 3 differences statistically meaningful)?

---

## Tau Values & Rationale

| τ | Rationale |
|---|---|
| 0.8 | Boundary case: τ ≈ γ_max, maximum spillover confounding |
| 1.0 | Current baseline — reproduces existing results exactly |
| 1.5 | Moderate effect, first power curve inflection |
| 2.0 | Larger effect, designs should approach full power |
| 3.0 | Strong effect, tests estimator behavior at high signal-to-noise |

---

## Scenario Count

```
5 tau values × 2,560 existing scenarios = 12,800 total scenarios
```

Breakdown: 5 (inc_config) × 2 (nb_type) × 4 (rho) × 4 (gamma) × 2 (spill_type) × 8 (design) × 5 (tau) = 12,800

---

## Files Changed

### 1. `code/05_run_simulation.R` — Local runner

**CONFIGURABLE PARAMETERS block (~line 65):**
```r
# Replace scalar:
# true_tau <- 1.0
# With vector:
true_tau_vals <- c(0.8, 1.0, 1.5, 2.0, 3.0)
```

**Outer loop:** Add `for (true_tau in true_tau_vals)` wrapping the existing incidence-config loop. All inner code already references `true_tau` correctly (lines 222, 249–250, 262, 264) and requires no modification.

**Seed string (~line 178):** Add `true_tau` to ensure independent randomness across tau values:
```r
seed_string <- paste(inc_mode, rho_x, nb_type, rho, gamma_val, spill_type, d_id, true_tau, sep = "|")
```

**Results data frame:** Add three new columns:
```r
True_Tau    = true_tau,
N_Valid_Est = sum(!is.na(valid_est)),   # exact count for SE computation
Power       = mean(all_ci_lower[valid_ci] > 0, na.rm = TRUE)
# Power: fraction of simulations where CI excludes zero (H0: tau=0 rejected)
```

**Derived SE formulas (computed post-hoc in `06_visualizations.R` from stored columns):**
```r
SE_Bias     = SD / sqrt(N_Valid_Est)
SE_Coverage = sqrt(Coverage * (1 - Coverage) / N_Valid_Est)   # binomial
SE_MSE      = sqrt(2*SD^4 + 4*Bias^2*SD^2) / sqrt(N_Valid_Est)  # under normality
```
These SEs enable ± error bars on all tables and plots, and are essential for interpreting the close D3 vs. D8 MSE comparison (0.080 vs. 0.079).

**Output filename:** Append `"tau_sweep"` to mode tag (e.g., `"MLE_tau_sweep"`) so files don't overwrite existing `MLE_combined` results.

**Local test workflow:** Set `true_tau_vals <- c(1.0)` to reproduce current results. Use `true_tau_vals <- c(1.0, 2.0)` with low rep counts for quick smoke test.

---

### 2. `longleaf_setup/simulation.R` — HPC per-task worker

**Parameter grid (~lines 107–127):** Add tau dimension:
```r
param_grid <- expand.grid(
  inc_config_id = 1:5,
  nb_type       = c("rook", "queen"),
  rho           = c(0.00, 0.01, 0.20, 0.50),
  gamma_val     = c(0.5, 0.6, 0.7, 0.8),
  spill_type    = c("control_only", "both"),
  design_id     = 1:8,
  true_tau      = c(0.8, 1.0, 1.5, 2.0, 3.0)  # NEW
)
# Grid now 12,800 rows
```

**Parameter extraction (~line 134):** Remove hardcoded `true_tau <- 1.0` (line 90); extract from params:
```r
true_tau <- params$true_tau
```

**Seed string:** Add `true_tau` (same as local script)

**Result data frame:** Add `True_Tau`, `N_Valid_Est`, and `Power` columns (same logic as local script)

---

### 3. `longleaf_setup/aggregate_results.R`

```r
expected_total <- 12800  # was 2560
```
Update output file naming to use `_tau_sweep` suffix.

---

### 4. `longleaf_setup/submit_array.sl`

```bash
#SBATCH --array=1-12800%100   # was 1-2560%100
```

---

### 5. `code/06_visualizations.R` — New and updated functions

**Backward-compatibility:** Add optional `tau_val = NULL` parameter to key plot/table functions. When `NULL`, uses all tau values (existing behavior preserved). When specified, filters to that tau value first.

**New functions:**

| Function | Description |
|---|---|
| `plot_mse_vs_tau(results, inc_label)` | MSE ± SE per design as lines across τ values (error bands via SE_MSE) |
| `plot_best_design_per_tau(results, inc_label)` | Faceted bar chart: for each τ, winning design and its MSE. One panel per incidence config |
| `plot_rank_trajectories_by_tau(results, inc_label)` | Design rank trajectory with τ on x-axis (rank 1 = best at top); shows whether D3/D8 stay on top |
| `plot_coverage_vs_tau(results, inc_label)` | Coverage ± SE per design across τ; flags calibration drift at extreme τ |
| `plot_power_curves(results, inc_label)` | Power vs. τ, one line per design; primary power comparison figure |
| `add_mc_ses(results)` | Helper: computes SE_Bias, SE_Coverage, SE_MSE from stored columns; returns results with SE columns appended |

---

### 6. `code/10_statistical_comparisons.R`

**`run_conditional_tests()`:** Add `True_Tau` to the stratification dimension list. Conditional Friedman tests will now be runnable stratified by τ value, answering: are D3/D8 rankings robust across all effect sizes?

**New derived analyses (computable from existing results columns without new simulation metrics):**

- **Relative efficiency:** `MSE(Design X, τ) / MSE(Best Design, τ)` per τ value. Shows whether the penalty for using a suboptimal design changes with effect size. Add as a table and line plot.

- **Bias decomposition across τ:** Plot mean `Bias` per design vs. τ. MLE consistency prediction: bias stays near zero across all τ. A bias spike at τ = 0.8 (where treatment ≈ spillover) would reveal spillover-confounding issues.

---

### 7. `code/complete_after_mle.R` — Post-processing documentation generator

This script auto-generates `CLAUDE.md` and `README.md` with live stats after a simulation run. It is NOT simulation code — it writes documentation templates.

- **Line 139** contains the doc string `"The estimand is **tau = 1.0**"` — update to reflect tau is now a swept parameter (e.g., `"tau ∈ {0.8, 1.0, 1.5, 2.0, 3.0}"`)
- **Line 260** contains a code block showing `true_tau <- 1.0` under "Fixed DGP Parameters" — move `true_tau` to a new "Swept Parameters" section listing it alongside the other sweep dimensions
- Update the summary statistics block to stratify by `True_Tau` when present in results
- Update the "Metrics tracked per scenario" doc string to include `Power` alongside Coverage

---

## Invariants Preserved

1. Existing `sim_results_MLE_combined_20260322_151030.rds` (τ = 1.0 only) is untouched — written to a separate file.
2. `04_estimation.R` unchanged — it is τ-agnostic.
3. `01_spatial_setup.R`, `02_incidence_generation.R`, `03_designs.R` unchanged.
4. Incidence generated once per `(mode, rho_X)` config — invariant holds because tau is in the outer loop, outside the incidence generation block.
5. Per-scenario digest seed now includes `true_tau`, ensuring different τ values have independent randomization.

---

## Verification Plan

1. **Smoke test (local):** Set `true_tau_vals <- c(1.0)` in updated `05_run_simulation.R`. Results should match existing `MLE_combined` file within Monte Carlo noise.
2. **Cross-tau test (local):** Set `true_tau_vals <- c(1.0, 2.0)` with 2 design resamples and 2 outcome resamples. Confirm `True_Tau` column appears, `Power` column appears, Bias at τ = 2.0 is near zero, Coverage ≈ 0.95.
3. **Power sanity check:** At τ = 3.0, all designs should have Power near 1.0. At τ = 0.8, designs should differ in Power and that difference should be informative.
4. **Longleaf dry run:** Submit a small test array (e.g., `--array=1-10`) before full 12,800-scenario run.
5. **Existing results unchanged:** Load `sim_results_MLE_combined_20260322_151030.rds` after all changes — all columns and values should be identical.

---

## Implementation Dependencies

This feature requires no external packages beyond those already used. The full implementation can be done in a single session. Suggested implementation order:

1. Update `05_run_simulation.R` + run local smoke test
2. Update `longleaf_setup/simulation.R` + `aggregate_results.R` + `submit_array.sl`
3. Update `06_visualizations.R` (new functions + tau_val parameter)
4. Update `10_statistical_comparisons.R` (new stratification + derived analyses)
5. Update `complete_after_mle.R`

**Note:** Simulation SEs and Statistical Power are bundled into this same implementation — no separate prior run needed.
