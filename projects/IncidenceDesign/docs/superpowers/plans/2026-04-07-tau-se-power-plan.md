# Combined Simulation Extension Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add tau sensitivity sweep (5 values), Monte Carlo SEs, and statistical power to the simulation pipeline — bundled into one implementation pass producing 12,800 total scenarios.

**Architecture:** Extend the existing simulation loop in `05_run_simulation.R` with an outer `for (true_tau in true_tau_vals)` loop and 3 new result columns (True_Tau, N_Valid_Est, Power). Mirror changes in the Longleaf HPC worker. Add 6 new visualization functions and 2 new statistical comparison functions, all backward-compatible with existing single-tau results.

**Tech Stack:** R, ggplot2, dplyr, tidyr, viridis, spatialreg, digest, SLURM

**Spec:** `docs/superpowers/specs/2026-04-07-tau-sensitivity-design.md`

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `code/05_run_simulation.R` | Modify | Add tau loop + 3 new columns |
| `longleaf_setup/simulation.R` | Modify | Add tau to grid + extraction + 3 new columns |
| `longleaf_setup/aggregate_results.R` | Modify | Update expected_total + tau summary |
| `longleaf_setup/submit_array.sl` | Modify | Update array size |
| `code/06_visualizations.R` | Modify | Add SE helper + 5 tau plots + update runner |
| `code/10_statistical_comparisons.R` | Modify | Add True_Tau stratum + 2 new functions |
| `code/complete_after_mle.R` | Modify | Update doc strings |

---

### Task 1: Update local simulation runner

**Files:**
- Modify: `code/05_run_simulation.R:65` (true_tau scalar)
- Modify: `code/05_run_simulation.R:106-108` (scenarios_per_config)
- Modify: `code/05_run_simulation.R:155-157` (insert tau loop)
- Modify: `code/05_run_simulation.R:177-178` (seed string)
- Modify: `code/05_run_simulation.R:253-282` (results data.frame — both branches)

- [ ] **Step 1: Replace scalar true_tau with vector**

At line 65, replace:
```r
true_tau <- 1.0
```
with:
```r
# Treatment effect values for tau sensitivity sweep
# Set to c(1.0) for single-tau backward-compatible runs
true_tau_vals <- c(0.8, 1.0, 1.5, 2.0, 3.0)
```

- [ ] **Step 2: Update scenarios_per_config calculation**

At lines 106-108, replace:
```r
scenarios_per_config <- length(neighbor_types) * length(rho_vals) *
                        length(gamma_vals) * length(spill_types) * length(design_ids)
total_scenarios <- length(inc_configs) * scenarios_per_config
```
with:
```r
scenarios_per_config <- length(true_tau_vals) * length(neighbor_types) * length(rho_vals) *
                        length(gamma_vals) * length(spill_types) * length(design_ids)
total_scenarios <- length(inc_configs) * scenarios_per_config
```

- [ ] **Step 3: Add outer tau loop inside run_incidence_config**

After line 155 (`config_times <- numeric(0)`), insert:
```r

  for (true_tau in true_tau_vals) {
    if (verbose) cat(sprintf("  --- true_tau = %.1f ---\n", true_tau))
```

After the closing brace of `for (nb_type in neighbor_types)` (line ~299), insert:
```r
  }  # end true_tau loop
```

- [ ] **Step 4: Add true_tau to seed string**

At line 178, replace:
```r
            seed_string <- paste(inc_mode, rho_x, nb_type, rho,
                                 gamma_val, spill_type, d_id, sep = "|")
```
with:
```r
            seed_string <- paste(inc_mode, rho_x, nb_type, rho,
                                 gamma_val, spill_type, d_id, true_tau, sep = "|")
```

- [ ] **Step 5: Add 3 new columns to valid-estimates branch**

In the valid-estimates branch of the results data.frame (lines 253-267), add these columns after `Fail_Rate = fail_rate,`:
```r
                True_Tau    = true_tau,
                N_Valid_Est = sum(valid_idx),
                Power       = if (sum(valid_ci) > 0) {
                  mean(all_ci_lower[valid_ci] > 0, na.rm = TRUE)
                } else NA,
```

- [ ] **Step 6: Add 3 new columns to NA/else branch**

In the else/NA branch (lines 270-281), add these columns after `Coverage = NA, Fail_Rate = fail_rate,`:
```r
                True_Tau    = true_tau,
                N_Valid_Est = 0L,
                Power       = NA,
```

- [ ] **Step 7: Commit**

```bash
git add code/05_run_simulation.R
git commit -m "Add tau sweep loop, Power, N_Valid_Est to local simulation runner

- Replace scalar true_tau=1.0 with true_tau_vals vector {0.8, 1.0, 1.5, 2.0, 3.0}
- Add outer for(true_tau) loop wrapping nb_type loop in run_incidence_config()
- Add true_tau to per-scenario digest seed string
- Add True_Tau, N_Valid_Est, Power columns to results data.frame
- Scenarios per config now multiplied by length(true_tau_vals)
- Set true_tau_vals <- c(1.0) for backward-compatible single-tau runs"
```

---

### Task 2: Update HPC worker

**Files:**
- Modify: `longleaf_setup/simulation.R:91` (remove hardcoded true_tau)
- Modify: `longleaf_setup/simulation.R:115-123` (expand.grid)
- Modify: `longleaf_setup/simulation.R:134-145` (parameter extraction + cat)
- Modify: `longleaf_setup/simulation.R:205-207` (seed string)
- Modify: `longleaf_setup/simulation.R:282-310` (results data.frame)

- [ ] **Step 1: Remove hardcoded true_tau**

Delete line 91:
```r
true_tau <- 1.0
```

- [ ] **Step 2: Add true_tau to expand.grid**

At lines 115-123, replace:
```r
param_grid <- expand.grid(
  inc_config_id = 1:5,
  nb_type       = c("rook", "queen"),
  rho           = c(0.00, 0.01, 0.20, 0.50),
  gamma_val     = c(0.5, 0.6, 0.7, 0.8),
  spill_type    = c("control_only", "both"),
  design_id     = 1:8,
  stringsAsFactors = FALSE
)
```
with:
```r
param_grid <- expand.grid(
  inc_config_id = 1:5,
  nb_type       = c("rook", "queen"),
  rho           = c(0.00, 0.01, 0.20, 0.50),
  gamma_val     = c(0.5, 0.6, 0.7, 0.8),
  spill_type    = c("control_only", "both"),
  design_id     = 1:8,
  true_tau      = c(0.8, 1.0, 1.5, 2.0, 3.0),
  stringsAsFactors = FALSE
)
```

- [ ] **Step 3: Extract true_tau from params**

After `d_id <- params$design_id` (line 143), add:
```r
true_tau   <- params$true_tau
```

Update the cat statement (line 145) to include true_tau:
```r
cat(sprintf("Scenario: %s(rho_x=%.2f) | nb=%s | rho=%.2f | gamma=%.2f | spill=%s | D%d | tau=%.1f\n",
            inc_mode, rho_x, nb_type, rho, gamma_val, spill_type, d_id, true_tau))
```

- [ ] **Step 4: Add true_tau to seed string**

At line 205, replace:
```r
scenario_seed <- digest2int(paste(inc_mode, rho_x, nb_type, rho,
                                  gamma_val, spill_type, d_id, sep = "|"))
```
with:
```r
scenario_seed <- digest2int(paste(inc_mode, rho_x, nb_type, rho,
                                  gamma_val, spill_type, d_id, true_tau, sep = "|"))
```

- [ ] **Step 5: Add 3 new columns to valid-estimates branch**

In the valid-estimates branch of the result data.frame (around line 290), add after `Fail_Rate = fail_rate,`:
```r
    True_Tau    = true_tau,
    N_Valid_Est = sum(valid_idx),
    Power       = if (sum(valid_ci) > 0) {
      mean(all_ci_lower[valid_ci] > 0, na.rm = TRUE)
    } else NA,
```

- [ ] **Step 6: Add 3 new columns to NA/else branch**

In the else/NA branch (around line 305), add after `Coverage = NA, Fail_Rate = fail_rate,`:
```r
    True_Tau    = true_tau,
    N_Valid_Est = 0L,
    Power       = NA,
```

- [ ] **Step 7: Update comments**

Update the header comment (~line 3-6) to say "12,800 scenarios" instead of "2,560 scenarios".

- [ ] **Step 8: Commit**

```bash
git add longleaf_setup/simulation.R
git commit -m "Add tau dimension to HPC simulation worker

- Remove hardcoded true_tau=1.0; extract from param_grid row
- Add true_tau to expand.grid: 2,560 -> 12,800 total scenarios
- Add true_tau to seed string for independent randomization
- Add True_Tau, N_Valid_Est, Power columns to result data.frame"
```

---

### Task 3: Update HPC aggregation and submission scripts

**Files:**
- Modify: `longleaf_setup/aggregate_results.R:42`
- Modify: `longleaf_setup/submit_array.sl:6,17`

- [ ] **Step 1: Update expected_total in aggregate_results.R**

At line 42, replace:
```r
expected_total <- 2560  # 5 configs x 2 nb x 4 rho x 4 gamma x 2 spill x 8 designs
```
with:
```r
expected_total <- 12800  # 5 configs x 2 nb x 4 rho x 4 gamma x 2 spill x 8 designs x 5 tau
```

- [ ] **Step 2: Add tau-stratified summary to aggregate_results.R**

After the existing MSE summary by Design (around line 100), add:
```r
if ("True_Tau" %in% names(all_results)) {
  cat("\nMSE by True_Tau (averaged across all designs and scenarios):\n")
  tau_summary <- all_results %>%
    group_by(True_Tau) %>%
    summarise(Mean_MSE = mean(MSE, na.rm = TRUE),
              Mean_Coverage = mean(Coverage, na.rm = TRUE),
              Mean_Power = mean(Power, na.rm = TRUE),
              .groups = "drop")
  print(as.data.frame(tau_summary))
}
```

- [ ] **Step 3: Update SLURM array size**

In `submit_array.sl`, replace line 6:
```bash
# Runs 2,560 scenarios in parallel (one per array task).
```
with:
```bash
# Runs 12,800 scenarios in parallel (one per array task).
```

Replace line 17:
```bash
#SBATCH --array=1-2560%100              # 2,560 tasks, max 100 concurrent
```
with:
```bash
#SBATCH --array=1-12800%100             # 12,800 tasks, max 100 concurrent
```

- [ ] **Step 4: Commit**

```bash
git add longleaf_setup/aggregate_results.R longleaf_setup/submit_array.sl
git commit -m "Update HPC aggregation and submission for 12,800-scenario grid

- aggregate_results.R: expected_total 2560 -> 12800, add tau-stratified summary
- submit_array.sl: SLURM array 1-2560 -> 1-12800"
```

---

### Task 4: Add Monte Carlo SE helper to visualizations

**Files:**
- Modify: `code/06_visualizations.R` (insert after line 61)

- [ ] **Step 1: Add add_mc_ses() function**

Insert after the `split_by_incidence_config()` function (after line 61):
```r

# ==============================================================================
# HELPER: Monte Carlo Standard Errors ----
# ==============================================================================

#' Compute Monte Carlo standard errors for simulation summary statistics
#'
#' Requires N_Valid_Est column (count of valid estimates per scenario).
#' Falls back to approximation from Fail_Rate if N_Valid_Est is missing.
#' SE formulas:
#'   SE(Bias)     = SD / sqrt(N)                   -- exact
#'   SE(Coverage) = sqrt(p(1-p) / N)               -- binomial proportion
#'   SE(MSE)      = sqrt(2*SD^4 + 4*Bias^2*SD^2) / sqrt(N) -- under normality
#'
#' @param results Data frame of simulation results
#' @return results with SE_Bias, SE_Coverage, SE_MSE columns appended
add_mc_ses <- function(results) {
  if (!"N_Valid_Est" %in% names(results)) {
    warning("N_Valid_Est column not found; approximating from Fail_Rate (n_iter=250)")
    results$N_Valid_Est <- round(250 * (1 - results$Fail_Rate))
  }
  results %>%
    mutate(
      SE_Bias     = SD / sqrt(pmax(N_Valid_Est, 1)),
      SE_Coverage = sqrt(Coverage * (1 - Coverage) / pmax(N_Valid_Est, 1)),
      SE_MSE      = sqrt(2 * SD^4 + 4 * Bias^2 * SD^2) / sqrt(pmax(N_Valid_Est, 1))
    )
}
```

- [ ] **Step 2: Commit**

```bash
git add code/06_visualizations.R
git commit -m "Add Monte Carlo SE computation helper to visualizations

add_mc_ses() computes SE_Bias, SE_Coverage, SE_MSE from stored
columns. Falls back to Fail_Rate approximation for old results."
```

---

### Task 5: Add 5 tau-dimension plot functions

**Files:**
- Modify: `code/06_visualizations.R` (insert before `run_all_visualizations` at line 713)

- [ ] **Step 1: Add plot_mse_vs_tau()**

Insert before `run_all_visualizations()`:
```r
# ==============================================================================
# TAU SENSITIVITY: MSE vs. Tau per Design ----
# ==============================================================================

#' Plot MSE per design as lines across tau values with 95% MC CI ribbons
#' @param results Data frame with True_Tau column
#' @param inc_label Character label for the incidence config
#' @return ggplot object (invisible NULL if True_Tau not present)
plot_mse_vs_tau <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results)) {
    warning("True_Tau column not found; skipping plot_mse_vs_tau")
    return(invisible(NULL))
  }
  results <- add_mc_ses(results)

  tau_mse <- results %>%
    group_by(Design, True_Tau) %>%
    summarise(
      Avg_MSE    = mean(MSE, na.rm = TRUE),
      Avg_SE_MSE = mean(SE_MSE, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(tau_mse, aes(x = True_Tau, y = Avg_MSE, color = Design, group = Design)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = pmax(Avg_MSE - 1.96 * Avg_SE_MSE, 0),
                    ymax = Avg_MSE + 1.96 * Avg_SE_MSE, fill = Design),
                alpha = 0.15, color = NA) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = expression("MSE vs. True Treatment Effect (" * tau * ")"),
      subtitle = paste0("Incidence: ", inc_label, " | Ribbons = 95% MC CI"),
      x = expression("True " * tau),
      y = "Average MSE"
    ) +
    theme(legend.position = "right")
}
```

- [ ] **Step 2: Add plot_best_design_per_tau()**

```r
# ==============================================================================
# TAU SENSITIVITY: Best Design per Tau ----
# ==============================================================================

#' Bar chart showing the winning design (lowest MSE) at each tau value
#' @param results Data frame with True_Tau column
#' @param inc_label Character label for the incidence config
#' @return ggplot object (invisible NULL if True_Tau not present)
plot_best_design_per_tau <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results)) return(invisible(NULL))

  best <- results %>%
    group_by(True_Tau, Design) %>%
    summarise(Avg_MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
    group_by(True_Tau) %>%
    mutate(Rank = rank(Avg_MSE, ties.method = "min")) %>%
    filter(Rank == 1) %>%
    ungroup()

  ggplot(best, aes(x = factor(True_Tau), y = Avg_MSE, fill = Design)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.85, width = 0.7) +
    geom_text(aes(label = gsub("Design ", "D", Design)),
              vjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = expression("Best Design (Lowest MSE) at Each " * tau),
      subtitle = paste0("Incidence: ", inc_label),
      x = expression("True " * tau),
      y = "Best Design's Avg MSE"
    ) +
    theme(legend.position = "right")
}
```

- [ ] **Step 3: Add plot_rank_trajectories_by_tau()**

```r
# ==============================================================================
# TAU SENSITIVITY: Rank Trajectories across Tau ----
# ==============================================================================

#' Rank trajectory plot with tau on x-axis (rank 1 = best at top)
#' Modeled after plot_rank_trajectories() in 08_design_recommendations.R
#' @param results Data frame with True_Tau column
#' @param inc_label Character label for the incidence config
#' @return ggplot object (invisible NULL if True_Tau not present)
plot_rank_trajectories_by_tau <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results)) return(invisible(NULL))

  ranked <- results %>%
    group_by(True_Tau, Design) %>%
    summarise(Avg_MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
    group_by(True_Tau) %>%
    mutate(Rank = rank(Avg_MSE, ties.method = "min")) %>%
    ungroup()

  ranked$Short_Name <- gsub("Design ", "D", ranked$Design)

  ggplot(ranked, aes(x = factor(True_Tau), y = Rank,
                     color = Design, group = Design)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(label = Short_Name), hjust = -0.15, vjust = -0.5,
              size = 2.8, show.legend = FALSE) +
    scale_y_reverse(breaks = 1:8, limits = c(8.5, 0.5)) +
    scale_color_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = expression("Design Rank Trajectories across " * tau),
      subtitle = paste0("Incidence: ", inc_label, " | Rank 1 = lowest avg MSE"),
      x = expression("True " * tau),
      y = "Rank (1 = Best)"
    ) +
    theme(legend.position = "right")
}
```

- [ ] **Step 4: Add plot_coverage_vs_tau()**

```r
# ==============================================================================
# TAU SENSITIVITY: Coverage vs. Tau per Design ----
# ==============================================================================

#' Coverage per design across tau values with 95% MC CI ribbons
#' @param results Data frame with True_Tau column
#' @param inc_label Character label for the incidence config
#' @return ggplot object (invisible NULL if True_Tau not present)
plot_coverage_vs_tau <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results)) return(invisible(NULL))
  results <- add_mc_ses(results)

  cov_tau <- results %>%
    filter(!is.na(Coverage)) %>%
    group_by(Design, True_Tau) %>%
    summarise(
      Avg_Coverage = mean(Coverage, na.rm = TRUE),
      Avg_SE_Cov   = mean(SE_Coverage, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(cov_tau, aes(x = True_Tau, y = Avg_Coverage, color = Design, group = Design)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = pmax(Avg_Coverage - 1.96 * Avg_SE_Cov, 0),
                    ymax = pmin(Avg_Coverage + 1.96 * Avg_SE_Cov, 1),
                    fill = Design),
                alpha = 0.15, color = NA) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = expression("95% CI Coverage across " * tau),
      subtitle = paste0("Incidence: ", inc_label, " | Red dashed = nominal 95%"),
      x = expression("True " * tau),
      y = "Average Coverage"
    ) +
    ylim(0, 1) +
    theme(legend.position = "right")
}
```

- [ ] **Step 5: Add plot_power_curves()**

```r
# ==============================================================================
# TAU SENSITIVITY: Statistical Power Curves ----
# ==============================================================================

#' Power (P(reject H0: tau=0)) vs tau, one line per design
#' @param results Data frame with True_Tau and Power columns
#' @param inc_label Character label for the incidence config
#' @return ggplot object (invisible NULL if columns not present)
plot_power_curves <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results) || !"Power" %in% names(results)) {
    return(invisible(NULL))
  }

  power_data <- results %>%
    filter(!is.na(Power)) %>%
    group_by(Design, True_Tau) %>%
    summarise(Avg_Power = mean(Power, na.rm = TRUE), .groups = "drop")

  ggplot(power_data, aes(x = True_Tau, y = Avg_Power, color = Design, group = Design)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    scale_color_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = expression("Statistical Power vs. True " * tau),
      subtitle = paste0("Incidence: ", inc_label,
                        " | Power = P(CI excludes 0) | Dashed = 80% threshold"),
      x = expression("True " * tau),
      y = "Average Power"
    ) +
    ylim(0, 1) +
    theme(legend.position = "right")
}
```

- [ ] **Step 6: Commit**

```bash
git add code/06_visualizations.R
git commit -m "Add 5 tau-dimension plot functions to visualizations

- plot_mse_vs_tau: MSE per design across tau with SE ribbons
- plot_best_design_per_tau: winning design at each tau
- plot_rank_trajectories_by_tau: rank stability across tau
- plot_coverage_vs_tau: coverage calibration across tau
- plot_power_curves: power vs tau per design
All return invisible(NULL) if True_Tau column absent (backward-compat)."
```

---

### Task 6: Integrate tau plots into run_all_visualizations

**Files:**
- Modify: `code/06_visualizations.R:796-797` (inside per-config loop)

- [ ] **Step 1: Add tau plot calls to per-config loop**

Inside `run_all_visualizations()`, after the `plot_coverage_mse_tradeoff` call (line 796) and before the `if (output_pdf)` check (line 798), insert:

```r

    # --- Tau sensitivity plots (only if multiple tau values present) ---
    if ("True_Tau" %in% names(config_results) &&
        length(unique(config_results$True_Tau)) > 1) {
      cat("  Plot T1: MSE vs Tau...\n")
      p_t1 <- plot_mse_vs_tau(config_results, config_name)
      if (!is.null(p_t1)) print(p_t1)

      cat("  Plot T2: Best Design per Tau...\n")
      p_t2 <- plot_best_design_per_tau(config_results, config_name)
      if (!is.null(p_t2)) print(p_t2)

      cat("  Plot T3: Rank Trajectories by Tau...\n")
      p_t3 <- plot_rank_trajectories_by_tau(config_results, config_name)
      if (!is.null(p_t3)) print(p_t3)

      cat("  Plot T4: Coverage vs Tau...\n")
      p_t4 <- plot_coverage_vs_tau(config_results, config_name)
      if (!is.null(p_t4)) print(p_t4)

      cat("  Plot T5: Power Curves...\n")
      p_t5 <- plot_power_curves(config_results, config_name)
      if (!is.null(p_t5)) print(p_t5)
    }
```

- [ ] **Step 2: Commit**

```bash
git add code/06_visualizations.R
git commit -m "Integrate tau plots into run_all_visualizations runner

Tau plots auto-activate when True_Tau column has >1 unique value.
Old single-tau results skip the block entirely."
```

---

### Task 7: Update statistical comparisons

**Files:**
- Modify: `code/10_statistical_comparisons.R:1016-1018` (cond_params)
- Modify: `code/10_statistical_comparisons.R` (add 2 new functions before generate_comparison_report)

- [ ] **Step 1: Add True_Tau to conditional test parameters**

At line 1016-1018, replace:
```r
  cond_params <- c("Rho", "Gamma", "Spillover_Type", "Neighbor_Type",
                    "Incidence_Mode")
```
with:
```r
  cond_params <- c("Rho", "Gamma", "Spillover_Type", "Neighbor_Type",
                    "Incidence_Mode", "True_Tau")
```

The `intersect()` on line 1018 ensures backward compatibility — if True_Tau is absent, it's silently dropped.

- [ ] **Step 2: Add table_relative_efficiency_by_tau() function**

Insert before `generate_comparison_report()` (~line 960):
```r
# ==============================================================================
# RELATIVE EFFICIENCY ACROSS TAU ----
# ==============================================================================

#' Compute relative efficiency: MSE(design) / MSE(best design) at each tau
#' 1.0 = best at that tau; higher = worse relative to best
#'
#' @param results Data frame with True_Tau column
#' @param inc_label Character label for the incidence config
#' @return Data frame of relative efficiency values (invisible)
table_relative_efficiency_by_tau <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results)) {
    cat("True_Tau not found; skipping relative efficiency.\n")
    return(invisible(NULL))
  }

  eff <- results %>%
    group_by(True_Tau, Design) %>%
    summarise(Avg_MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
    group_by(True_Tau) %>%
    mutate(
      Best_MSE = min(Avg_MSE, na.rm = TRUE),
      Rel_Efficiency = Avg_MSE / Best_MSE
    ) %>%
    ungroup()

  wide <- eff %>%
    select(Design, True_Tau, Rel_Efficiency) %>%
    tidyr::pivot_wider(names_from = True_Tau, values_from = Rel_Efficiency,
                       names_prefix = "tau=")

  cat(sprintf("\n=== Relative Efficiency by Tau | Incidence: %s ===\n", inc_label))
  cat("(1.0 = best at that tau; higher = worse relative to best)\n\n")
  print(as.data.frame(wide), row.names = FALSE, digits = 3)

  invisible(eff)
}
```

- [ ] **Step 3: Add plot_bias_across_tau() function**

```r
# ==============================================================================
# BIAS DECOMPOSITION ACROSS TAU ----
# ==============================================================================

#' Plot mean bias per design vs tau to test MLE consistency (bias near 0 for all tau)
#' A bias spike at tau=0.8 (where treatment ~ spillover) would reveal confounding.
#'
#' @param results Data frame with True_Tau column
#' @param inc_label Character label for the incidence config
#' @return ggplot object (invisible NULL if True_Tau not present)
plot_bias_across_tau <- function(results, inc_label = "") {
  if (!"True_Tau" %in% names(results)) return(invisible(NULL))

  bias_data <- results %>%
    group_by(Design, True_Tau) %>%
    summarise(
      Avg_Bias = mean(Bias, na.rm = TRUE),
      Avg_SD   = mean(SD, na.rm = TRUE),
      N        = n(),
      .groups = "drop"
    ) %>%
    mutate(SE_Bias = Avg_SD / sqrt(pmax(N, 1)))

  ggplot(bias_data, aes(x = True_Tau, y = Avg_Bias, color = Design, group = Design)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_ribbon(aes(ymin = Avg_Bias - 1.96 * SE_Bias,
                    ymax = Avg_Bias + 1.96 * SE_Bias,
                    fill = Design),
                alpha = 0.12, color = NA) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = expression("Bias Decomposition across " * tau),
      subtitle = paste0("Incidence: ", inc_label,
                        " | Dashed = zero bias (consistency)"),
      x = expression("True " * tau),
      y = "Average Bias"
    ) +
    theme(legend.position = "right")
}
```

- [ ] **Step 4: Commit**

```bash
git add code/10_statistical_comparisons.R
git commit -m "Add True_Tau stratification and tau-dimension analyses

- Add True_Tau to conditional Friedman test parameters
- Add table_relative_efficiency_by_tau(): MSE ratio vs best per tau
- Add plot_bias_across_tau(): MLE consistency diagnostic across tau"
```

---

### Task 8: Update documentation generator

**Files:**
- Modify: `code/complete_after_mle.R:139,141,145,260-261`

- [ ] **Step 1: Update estimand description**

At line 139, replace:
```
spillover. The estimand is **tau = 1.0** (direct treatment effect). We compare two
```
with:
```
spillover. The estimand is **tau** (direct treatment effect), swept across {0.8, 1.0, 1.5, 2.0, 3.0}. We compare two
```

- [ ] **Step 2: Update scenario count**

At line 140-141, replace:
```
estimators (DIM, MLE) across 2,560 parameter scenarios and 3 incidence generation
```
with:
```
estimators (DIM, MLE) across 12,800 parameter scenarios (5 tau values x 2,560 base grid) and 3 incidence generation
```

- [ ] **Step 3: Update metrics tracked**

At line 145, replace:
```
**Metrics tracked per scenario:** Bias, SD, MSE, Coverage (95%% CI), Fail_Rate.
```
with:
```
**Metrics tracked per scenario:** Bias, SD, MSE, Coverage (95%% CI), Power (P(CI excludes 0)), Fail_Rate, N_Valid_Est, True_Tau.
```

- [ ] **Step 4: Update Fixed DGP Parameters section**

At lines 260-261, remove `true_tau <- 1.0` from the Fixed DGP code block. Add after the closing triple-backtick (line 268):

```
## Swept Parameters

| Parameter | Values | Purpose |
|-----------|--------|---------|
| true_tau | {0.8, 1.0, 1.5, 2.0, 3.0} | Treatment effect sensitivity |

Combined with the existing 2,560-scenario base grid: 5 x 2,560 = **12,800 total scenarios**.
```

- [ ] **Step 5: Commit**

```bash
git add code/complete_after_mle.R
git commit -m "Update documentation generator for tau sweep

- Estimand description: tau=1.0 -> tau swept across {0.8,...,3.0}
- Scenario count: 2,560 -> 12,800
- Metrics: add Power, N_Valid_Est, True_Tau
- Fixed DGP: move true_tau to new Swept Parameters section"
```

---

## Verification

After all tasks are complete:

1. **Backward-compatibility check**: Load existing results `results/sim_data/sim_results_MLE_combined_20260322_151030.rds`. Run `run_all_visualizations()` — all existing plots should render without error. New tau plots should be silently skipped (True_Tau column absent).

2. **Local smoke test (single tau)**: In `05_run_simulation.R`, set `true_tau_vals <- c(1.0)`, `n_designs_resamples <- 2`, `n_outcome_resamples <- 2`, `design_ids <- c(1, 3)`. Run for one incidence config. Verify output has True_Tau, N_Valid_Est, Power columns. Verify results match expected structure.

3. **Local smoke test (multi tau)**: Set `true_tau_vals <- c(1.0, 2.0)` with same reduced params. Verify:
   - True_Tau column has both values
   - Power at tau=2.0 is higher than at tau=1.0 (stronger signal)
   - Coverage approximately 0.95 at both tau values
   - Bias near zero at both tau values

4. **Visualization test**: Load the multi-tau smoke test results. Run `plot_mse_vs_tau()`, `plot_power_curves()`, etc. Verify they render without error and show two tau values on x-axis.

5. **Existing results untouched**: Confirm `sim_results_MLE_combined_20260322_151030.rds` has not been modified (check file timestamp).
