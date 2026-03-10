# 03_analysis_functions.R
# Analysis and plotting functions for the SpillSpatialDepSim unified simulation.
# Replaces the 1,000+ lines of manual plotting across SimEstimateAnalysisAll.Rmd
# and SimEstimateAnalysisMeans.Rmd.
#
# Source this file from 05_analysis_report.Rmd. Does NOT source 01/02 — no
# simulation code here, only analysis on the saved RDS.

library(ggplot2)
library(dplyr)
library(tidyr)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

#' Load the unified simulation results RDS.
#'
#' @param results_dir Path to the results directory (default: "../results").
#' @return Data frame with columns: grid, psi, rho, trt_spill, combo_id,
#'   is_block, param, bias, bias_abs, variance, MSE.
load_unified_results <- function(results_dir = "../results") {
  rds_path <- file.path(results_dir, "sim_results_unified.rds")
  if (!file.exists(rds_path)) {
    stop("Results file not found: ", rds_path,
         "\nRun 04_run_simulation.R first.")
  }
  readRDS(rds_path)
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Human-readable label for a scenario's spillover type.
spill_label <- function(trt_spill) {
  ifelse(trt_spill, "TrtSpill (case 2)", "TrtNoSpill (case 1)")
}

#' Filter the results data frame to one combination of scenario parameters.
#'
#' @param df       Results data frame from load_unified_results().
#' @param grid_nm  Grid name string: "2x4", "3x3", or "3x4".
#' @param param_nm Parameter name: "alpha", "beta", "psi", or "rho".
#' @param spill    Logical trt_spill value (TRUE / FALSE).
#' @return Filtered data frame.
filter_results <- function(df, grid_nm = NULL, param_nm = NULL, spill = NULL) {
  if (!is.null(grid_nm)) df <- df[df$grid == grid_nm, ]
  if (!is.null(param_nm)) df <- df[df$param == param_nm, ]
  if (!is.null(spill))    df <- df[df$trt_spill == spill, ]
  df
}


# ---------------------------------------------------------------------------
# Plot 1: MSE boxplot (mirrors SimEstimateAnalysisAll.Rmd style)
# ---------------------------------------------------------------------------

#' Box plot of MSE distributions by sampling method (random vs. block),
#' faceted by rho, with psi encoded as color.
#'
#' @param df        Filtered or full results data frame.
#' @param param_nm  Parameter to plot ("alpha", "beta", "psi", "rho").
#' @param grid_nm   Grid to focus on ("2x4", "3x3", or "3x4").
#' @param spill     Logical: filter to this trt_spill value (or NULL = both).
#' @return ggplot object.
plot_mse_boxplot <- function(df, param_nm = "beta", grid_nm, spill = NULL) {
  d <- filter_results(df, grid_nm = grid_nm, param_nm = param_nm, spill = spill)
  d$method <- ifelse(d$is_block, "Block Stratification", "Random Sampling")
  d$psi_f  <- factor(paste("psi =", d$psi))
  d$rho_f  <- factor(paste("rho =", d$rho))

  ggplot(d, aes(x = method, y = MSE, fill = psi_f)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
    facet_wrap(~ rho_f, labeller = label_value) +
    scale_fill_brewer(palette = "Set2", name = "Spillover (psi)") +
    labs(
      title   = sprintf("MSE Distribution — %s grid, param = %s%s",
                        grid_nm, param_nm,
                        if (!is.null(spill)) paste0(", ", spill_label(spill)) else ""),
      x       = "Treatment Assignment Method",
      y       = "MSE"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1))
}


# ---------------------------------------------------------------------------
# Plot 2: Mean MSE + CI (mirrors SimEstimateAnalysisMeans.Rmd style)
# ---------------------------------------------------------------------------

#' Point + 95% CI plot of mean MSE per combo, grouped by sampling method,
#' colored by psi, faceted by rho x trt_spill.
#'
#' @param df        Full or filtered results data frame.
#' @param param_nm  Parameter to plot.
#' @param grid_nm   Grid to focus on.
#' @return ggplot object.
plot_mse_means <- function(df, param_nm = "beta", grid_nm) {
  d <- filter_results(df, grid_nm = grid_nm, param_nm = param_nm)
  d$method    <- ifelse(d$is_block, "Block Stratification", "Random Sampling")
  d$psi_f     <- factor(paste("psi =", d$psi))
  d$rho_f     <- factor(paste("rho =", d$rho))
  d$spill_f   <- factor(spill_label(d$trt_spill))
  d$se_mse    <- sqrt(d$variance / 1)  # 1 obs per combo (MSE already aggregated)

  ggplot(d, aes(x = method, y = MSE, color = psi_f,
                ymin = MSE - 1.96 * d$bias_abs,
                ymax = MSE + 1.96 * d$bias_abs)) +
    geom_point(position = position_jitter(width = 0.15), size = 1.5, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
                 position = position_dodge(0.5)) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 position = position_dodge(0.5), width = 0.3) +
    facet_grid(spill_f ~ rho_f) +
    scale_color_brewer(palette = "Set1", name = "Spillover (psi)") +
    labs(
      title = sprintf("Mean MSE by Assignment Method — %s grid, param = %s", grid_nm, param_nm),
      x     = "Treatment Assignment Method",
      y     = "MSE"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1))
}


# ---------------------------------------------------------------------------
# Plot 3: All grids side-by-side (facet_wrap on grid)
# ---------------------------------------------------------------------------

#' Faceted MSE box plot across all three grid sizes.
#'
#' @param df        Full results data frame.
#' @param param_nm  Parameter to plot.
#' @param spill     Logical trt_spill filter (NULL = both, combined).
#' @return ggplot object.
plot_mse_all_grids <- function(df, param_nm = "beta", spill = NULL) {
  d <- filter_results(df, param_nm = param_nm, spill = spill)
  d$method <- ifelse(d$is_block, "Block", "Random")
  d$grid_f <- factor(d$grid, levels = c("2x4", "3x3", "3x4"))
  d$psi_f  <- factor(paste("psi =", d$psi))

  ggplot(d, aes(x = method, y = MSE, fill = psi_f)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
    facet_wrap(~ grid_f, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", name = "Spillover (psi)") +
    labs(
      title = sprintf("MSE Distribution across all grids — param = %s%s",
                      param_nm,
                      if (!is.null(spill)) paste0(", ", spill_label(spill)) else ""),
      x = "Assignment Method",
      y = "MSE"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}


# ---------------------------------------------------------------------------
# Analysis 1: Optimal combos (mirrors OptimalTreatmentAssignmentCombos.Rmd)
# ---------------------------------------------------------------------------

#' Find the minimum-MSE combo for each scenario.
#'
#' @param df        Full results data frame.
#' @param param_nm  Parameter to minimize MSE for (default "beta").
#' @return Data frame: one row per scenario with grid, psi, rho, trt_spill,
#'   best_combo_id, is_block, MSE.
get_optimal_combos <- function(df, param_nm = "beta") {
  d <- df[df$param == param_nm, ]
  d |>
    group_by(grid, psi, rho, trt_spill) |>
    slice_min(order_by = MSE, n = 1, with_ties = FALSE) |>
    ungroup() |>
    rename(best_combo_id = combo_id, best_MSE = MSE) |>
    select(grid, psi, rho, trt_spill, best_combo_id, is_block, best_MSE)
}


# ---------------------------------------------------------------------------
# Summary table: Random vs. Block MSE comparison
# ---------------------------------------------------------------------------

#' Summary table comparing mean MSE for Random vs. Block Stratification.
#'
#' @param df        Full results data frame.
#' @param param_nm  Parameter to summarize.
#' @return Data frame with columns: grid, psi, rho, trt_spill, method, mean_MSE, median_MSE.
table_mse_comparison <- function(df, param_nm = "beta") {
  d <- df[df$param == param_nm, ]
  d$method <- ifelse(d$is_block, "Block Stratification", "Random Sampling")
  d |>
    group_by(grid, psi, rho, trt_spill, method) |>
    summarise(
      mean_MSE   = mean(MSE,   na.rm = TRUE),
      median_MSE = median(MSE, na.rm = TRUE),
      n_combos   = n(),
      .groups = "drop"
    ) |>
    arrange(grid, psi, rho, trt_spill, method)
}
