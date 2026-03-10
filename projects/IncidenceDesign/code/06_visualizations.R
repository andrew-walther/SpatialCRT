# ==============================================================================
# 06_visualizations.R
# All plots and tables for the heterogeneous incidence simulation study.
# Can be run independently after 05_run_simulation.R has saved results.
#
# IMPORTANT: Results for each incidence generation method (iid Uniform,
#   Spatial, Poisson) are reported SEPARATELY. They are never aggregated
#   across incidence modes because each represents a fundamentally different
#   assumption about how prior observed incidence arises.
#
# Usage:
#   source("06_visualizations.R")
#   -- or load results manually and call individual plot functions --
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# --- Source modules for helper functions ---
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)

source(file.path(script_dir, "01_spatial_setup.R"))
source(file.path(script_dir, "02_incidence_generation.R"))
source(file.path(script_dir, "03_designs.R"))

# ==============================================================================
# HELPERS
# ==============================================================================

#' Create a human-readable label for an incidence config
#' @param inc_mode Character: "iid", "spatial", or "poisson"
#' @param rho_x Numeric: rho_incidence value
#' @return Character label
inc_config_label <- function(inc_mode, rho_x) {
  if (inc_mode == "iid") return("iid Uniform")
  sprintf("%s (rho_X = %.2f)", tools::toTitleCase(inc_mode), rho_x)
}

#' Split combined results into a list of per-incidence-config data frames
#' @param results Combined results data frame
#' @return Named list of data frames, one per incidence config
split_by_incidence_config <- function(results) {
  configs <- results %>%
    distinct(Incidence_Mode, Rho_Incidence) %>%
    arrange(Incidence_Mode, Rho_Incidence)

  config_list <- list()
  for (i in seq_len(nrow(configs))) {
    mode <- configs$Incidence_Mode[i]
    rho_x <- configs$Rho_Incidence[i]
    label <- inc_config_label(mode, rho_x)
    config_list[[label]] <- results %>%
      filter(Incidence_Mode == mode, Rho_Incidence == rho_x)
  }
  config_list
}

# ==============================================================================
# LOAD RESULTS
# ==============================================================================

#' Load the most recent results file from the results directory
#' @param results_dir Path to results directory
#' @param estimation_mode Optional filter for "MLE" or "DIM"
#' @return Data frame of results
load_latest_results <- function(results_dir = file.path(dirname(script_dir), "results"),
                                estimation_mode = NULL) {
  files <- list.files(results_dir, pattern = "sim_results_.*\\.rds$",
                      full.names = TRUE)
  if (!is.null(estimation_mode)) {
    files <- files[grepl(estimation_mode, files)]
  }
  if (length(files) == 0) stop("No results files found in ", results_dir)

  # Pick the most recently modified
  latest <- files[which.max(file.mtime(files))]
  cat("Loading results from:", basename(latest), "\n")
  readRDS(latest)
}

# ==============================================================================
# PLOT 1: MSE Boxplot by Design, faceted by Neighbor Type
# ==============================================================================

plot_mse_by_neighbor <- function(results, inc_label = "") {
  ggplot(results, aes(x = reorder(Design, MSE, FUN = median), y = MSE, fill = Design)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    facet_wrap(~ Neighbor_Type,
               labeller = as_labeller(c("queen" = "Queen Contiguity",
                                        "rook" = "Rook Contiguity"))) +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Distribution of Estimation Error (MSE) by Spatial Contiguity",
      subtitle = paste0("Incidence: ", inc_label),
      x = "Sampling Design",
      y = "Mean Squared Error (Lower is Better)"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==============================================================================
# PLOT 2: MSE Heatmap by Rho x Design
# ==============================================================================

plot_mse_heatmap <- function(results, inc_label = "") {
  mse_summary <- results %>%
    group_by(Design, Rho, Neighbor_Type) %>%
    summarise(Avg_MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
    mutate(text_color = ifelse(Avg_MSE > quantile(Avg_MSE, 0.75, na.rm = TRUE),
                               "black", "white"))

  ggplot(mse_summary, aes(x = as.factor(Rho), y = Design, fill = Avg_MSE)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Avg_MSE), color = text_color),
              size = 3.5, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ Neighbor_Type,
               labeller = as_labeller(c("queen" = "Queen Contiguity",
                                        "rook" = "Rook Contiguity"))) +
    scale_fill_viridis_c(option = "magma") +
    scale_color_identity() +
    theme_minimal(base_size = 14) +
    labs(
      title = "Heatmap of Error by Spatial Dependence & Neighbor Type",
      subtitle = paste0("Incidence: ", inc_label),
      x = expression("Spatial Dependence (" * rho * ")"),
      y = "Sampling Design",
      fill = "Avg MSE"
    ) +
    theme(panel.grid = element_blank())
}

# ==============================================================================
# PLOT 3: Per-Design MSE Line Plots across Gamma
# ==============================================================================

plot_mse_per_design <- function(results, inc_label = "") {
  unique_designs <- sort(unique(results$Design))
  plots <- list()

  for (current_design in unique_designs) {
    design_data <- results %>% filter(Design == current_design)

    p <- ggplot(design_data,
                aes(x = as.factor(Gamma), y = MSE,
                    color = Spillover_Type, group = Spillover_Type)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      facet_grid(Neighbor_Type ~ paste("rho =", Rho)) +
      scale_color_viridis_d(option = "turbo",
                            labels = c("both" = "Treated & Control",
                                       "control_only" = "Control Only")) +
      theme_bw(base_size = 14) +
      labs(
        title = paste("MSE Dynamics for", current_design),
        subtitle = paste0("Incidence: ", inc_label),
        x = expression("Spillover Magnitude (" * gamma * ")"),
        y = "Mean Squared Error (MSE)",
        color = "Spillover Target"
      ) +
      theme(legend.position = "bottom",
            strip.background = element_rect(fill = "grey90"),
            strip.text = element_text(face = "bold"))

    plots[[current_design]] <- p
    print(p)
  }

  invisible(plots)
}

# ==============================================================================
# PLOT 4: Master Bar Chart Comparison
# ==============================================================================

plot_master_comparison <- function(results, nb_filter = "queen", inc_label = "") {
  comparison_data <- results %>% filter(Neighbor_Type == nb_filter)

  ggplot(comparison_data, aes(x = Design, y = MSE, fill = as.factor(Gamma))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.2) +
    facet_grid(
      Spillover_Type ~ paste("rho =", Rho),
      scales = "free_y",
      labeller = labeller(Spillover_Type = c("both" = "Spillover: Both",
                                             "control_only" = "Spillover: Control"))
    ) +
    scale_fill_viridis_d(option = "turbo",
                         name = expression("Spillover\nMagnitude (" * gamma * ")")) +
    theme_bw(base_size = 13) +
    labs(
      title = "Design Comparison: Mean Squared Error across Scenarios",
      subtitle = paste0(tools::toTitleCase(nb_filter), " Contiguity | Incidence: ", inc_label),
      x = "Sampling Design",
      y = "Mean Squared Error (MSE)"
    ) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
}

# ==============================================================================
# PLOT 5: Coverage by Design (per incidence config)
# ==============================================================================

plot_coverage_by_design <- function(results, inc_label = "") {
  if (!"Coverage" %in% names(results) || all(is.na(results$Coverage))) {
    return(invisible(NULL))
  }

  cov_data <- results %>%
    filter(!is.na(Coverage))

  ggplot(cov_data, aes(x = reorder(Design, Coverage, FUN = median),
                       y = Coverage, fill = Design)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ Neighbor_Type,
               labeller = as_labeller(c("queen" = "Queen Contiguity",
                                        "rook" = "Rook Contiguity"))) +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "95% CI Coverage Probability by Design",
      subtitle = paste0("Incidence: ", inc_label, " | Red dashed = nominal 95%"),
      x = "Sampling Design",
      y = "Coverage Probability"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1)
}

# ==============================================================================
# PLOT 6: Grid Heatmaps of Incidence Patterns (standalone, not per-config)
# ==============================================================================

plot_incidence_heatmaps <- function(grid_dim = 10) {
  grid_obj <- build_spatial_grid(grid_dim)
  N <- grid_obj$N_clusters
  coords <- grid_obj$coords

  set.seed(42)

  # Generate one realization per mode
  x_iid <- generate_incidence("iid", N, 1)
  x_spatial_low  <- generate_incidence("spatial", N, 1, grid_obj$W_queen, 0.20)
  x_spatial_high <- generate_incidence("spatial", N, 1, grid_obj$W_queen, 0.50)
  x_poisson_low  <- generate_incidence("poisson", N, 1, grid_obj$W_queen, 0.20,
                                        base_rate = 35/100000, pop_per_cluster = 1000,
                                        pop_mode = "equal")
  x_poisson_high <- generate_incidence("poisson", N, 1, grid_obj$W_queen, 0.50,
                                        base_rate = 35/100000, pop_per_cluster = 1000,
                                        pop_mode = "equal")

  df <- bind_rows(
    data.frame(coords, Incidence = x_iid[, 1],           Mode = "iid Uniform"),
    data.frame(coords, Incidence = x_spatial_low[, 1],    Mode = "Spatial (rho_X = 0.20)"),
    data.frame(coords, Incidence = x_spatial_high[, 1],   Mode = "Spatial (rho_X = 0.50)"),
    data.frame(coords, Incidence = x_poisson_low[, 1],    Mode = "Poisson (rho_X = 0.20)"),
    data.frame(coords, Incidence = x_poisson_high[, 1],   Mode = "Poisson (rho_X = 0.50)")
  )

  ggplot(df, aes(x = x, y = y, fill = Incidence)) +
    geom_tile(color = "grey50", linewidth = 0.3) +
    facet_wrap(~ Mode, ncol = 3) +
    scale_fill_viridis_c(option = "mako", direction = -1) +
    coord_fixed() +
    theme_void(base_size = 14) +
    labs(
      title = "Incidence Patterns by Generation Mode",
      subtitle = "Single realization on 10x10 grid",
      fill = "Incidence\n[0, 1]"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text = element_text(face = "bold", margin = margin(b = 5, t = 5)))
}

# ==============================================================================
# PLOT 7: Incidence Distribution Histograms (standalone, not per-config)
# ==============================================================================

plot_incidence_distributions <- function(grid_dim = 10, n_resamples = 50) {
  grid_obj <- build_spatial_grid(grid_dim)
  N <- grid_obj$N_clusters

  set.seed(42)

  x_iid <- as.vector(generate_incidence("iid", N, n_resamples))
  x_sp20 <- as.vector(generate_incidence("spatial", N, n_resamples, grid_obj$W_queen, 0.20))
  x_sp50 <- as.vector(generate_incidence("spatial", N, n_resamples, grid_obj$W_queen, 0.50))
  x_po20 <- as.vector(generate_incidence("poisson", N, n_resamples, grid_obj$W_queen, 0.20,
                                          base_rate = 35/100000, pop_per_cluster = 1000,
                                          pop_mode = "equal"))
  x_po50 <- as.vector(generate_incidence("poisson", N, n_resamples, grid_obj$W_queen, 0.50,
                                          base_rate = 35/100000, pop_per_cluster = 1000,
                                          pop_mode = "equal"))

  df <- data.frame(
    Value = c(x_iid, x_sp20, x_sp50, x_po20, x_po50),
    Mode = rep(c("iid Uniform", "Spatial (rho_X = 0.20)", "Spatial (rho_X = 0.50)",
                 "Poisson (rho_X = 0.20)", "Poisson (rho_X = 0.50)"),
               each = N * n_resamples)
  )

  ggplot(df, aes(x = Value, fill = Mode)) +
    geom_histogram(bins = 30, alpha = 0.8, color = "white", linewidth = 0.2) +
    facet_wrap(~ Mode, scales = "free_y") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Marginal Distribution of Incidence Values by Mode",
      subtitle = sprintf("Pooled across %d resamples on %dx%d grid", n_resamples, grid_dim, grid_dim),
      x = "Incidence Value [0, 1]",
      y = "Count"
    ) +
    theme(legend.position = "none")
}

# ==============================================================================
# ANALYSIS 1: Design Rank Frequency ("Win Rates")
#
# For each scenario (unique combo of Rho, Gamma, Spillover_Type, Neighbor_Type),
# rank the 6 designs by MSE. Then tally how often each design finishes in each
# rank position. A design that ranks #1 in 80% of scenarios is clearly dominant.
# ==============================================================================

table_design_ranks <- function(results, inc_label = "") {
  # Each scenario is a unique (Neighbor_Type, Rho, Gamma, Spillover_Type) combo
  ranked <- results %>%
    group_by(Neighbor_Type, Rho, Gamma, Spillover_Type) %>%
    mutate(Rank = rank(MSE, ties.method = "min")) %>%
    ungroup()

  n_scenarios <- ranked %>%
    distinct(Neighbor_Type, Rho, Gamma, Spillover_Type) %>%
    nrow()

  # Count how often each design gets each rank
  rank_counts <- ranked %>%
    group_by(Design, Rank) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(Pct = Count / n_scenarios * 100)

  # Pivot to wide format: rows = Design, cols = Rank
  rank_wide <- rank_counts %>%
    select(Design, Rank, Pct) %>%
    tidyr::pivot_wider(names_from = Rank, values_from = Pct,
                       names_prefix = "Rank_", values_fill = 0) %>%
    arrange(desc(Rank_1))

  cat(sprintf("\n=== Design Rank Frequency (%%) | Incidence: %s ===\n", inc_label))
  cat(sprintf("(%d scenarios; each cell = %% of scenarios design finishes at that rank)\n\n",
              n_scenarios))
  print(as.data.frame(rank_wide), row.names = FALSE, digits = 1)

  invisible(rank_wide)
}

plot_design_ranks <- function(results, inc_label = "") {
  ranked <- results %>%
    group_by(Neighbor_Type, Rho, Gamma, Spillover_Type) %>%
    mutate(Rank = rank(MSE, ties.method = "min")) %>%
    ungroup()

  ggplot(ranked, aes(x = as.factor(Rank), fill = Design)) +
    geom_bar(position = "dodge", color = "black", linewidth = 0.2) +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Design Rank Distribution across Scenarios",
      subtitle = paste0("Incidence: ", inc_label,
                        " | Rank 1 = lowest MSE in that scenario"),
      x = "Rank (1 = Best)",
      y = "Number of Scenarios"
    )
}

# ==============================================================================
# ANALYSIS 2: Robustness Profile
#
# For each design, show the best-case, 25th percentile, median, 75th
# percentile, and worst-case MSE. A design with low median but high worst-case
# is fragile. A design with moderate median but low worst-case is robust.
# ==============================================================================

table_robustness <- function(results, inc_label = "") {
  robust <- results %>%
    group_by(Design) %>%
    summarise(
      Best_MSE   = min(MSE, na.rm = TRUE),
      Q25_MSE    = quantile(MSE, 0.25, na.rm = TRUE),
      Median_MSE = median(MSE, na.rm = TRUE),
      Q75_MSE    = quantile(MSE, 0.75, na.rm = TRUE),
      Worst_MSE  = max(MSE, na.rm = TRUE),
      IQR_MSE    = IQR(MSE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Median_MSE)

  cat(sprintf("\n=== Robustness Profile | Incidence: %s ===\n", inc_label))
  cat("(Sorted by Median MSE; IQR = interquartile range = spread of middle 50%%)\n\n")
  print(as.data.frame(robust), row.names = FALSE, digits = 4)

  invisible(robust)
}

# ==============================================================================
# ANALYSIS 3: Pairwise Dominance Matrix
#
# For every pair of designs (A, B), compute the fraction of scenarios where
# A has lower MSE than B. If Design 7 beats Design 1 in 95% of scenarios,
# that's near-total dominance. Values near 50% mean the two designs are
# essentially interchangeable.
# ==============================================================================

table_pairwise_dominance <- function(results, inc_label = "") {
  designs <- sort(unique(results$Design))
  n_designs <- length(designs)

  # Build a wide table: one row per scenario, one column per design's MSE
  scenario_ids <- c("Neighbor_Type", "Rho", "Gamma", "Spillover_Type")
  wide <- results %>%
    select(all_of(c(scenario_ids, "Design", "MSE"))) %>%
    pivot_wider(names_from = Design, values_from = MSE,
                id_cols = all_of(scenario_ids))

  dom_matrix <- matrix(NA, nrow = n_designs, ncol = n_designs,
                       dimnames = list(designs, designs))

  for (i in seq_along(designs)) {
    for (j in seq_along(designs)) {
      if (i == j) next
      mse_i <- wide[[designs[i]]]
      mse_j <- wide[[designs[j]]]
      dom_matrix[i, j] <- mean(mse_i < mse_j, na.rm = TRUE) * 100
    }
  }

  cat(sprintf("\n=== Pairwise Dominance Matrix (%%) | Incidence: %s ===\n", inc_label))
  cat("(Cell [row, col] = %% of scenarios where ROW design has lower MSE than COL design)\n\n")
  print(round(dom_matrix, 1))

  invisible(dom_matrix)
}

# ==============================================================================
# ANALYSIS 4: Bias-Variance Decomposition
#
# MSE = Bias^2 + Variance. This plot shows the two components stacked for
# each design, revealing whether a design's error comes from systematic bias
# or from estimation variability.
# ==============================================================================

plot_bias_variance <- function(results, inc_label = "") {
  bv <- results %>%
    group_by(Design) %>%
    summarise(
      Avg_Bias_Sq = mean(Bias^2, na.rm = TRUE),
      Avg_Var     = mean(SD^2, na.rm = TRUE),
      Avg_MSE     = mean(MSE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c(Avg_Bias_Sq, Avg_Var),
                        names_to = "Component", values_to = "Value") %>%
    mutate(Component = ifelse(Component == "Avg_Bias_Sq",
                              "Bias-squared", "Variance"))

  ggplot(bv, aes(x = reorder(Design, Value, FUN = sum),
                 y = Value, fill = Component)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Bias-squared" = "#E69F00", "Variance" = "#56B4E9")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Bias-Variance Decomposition of MSE",
      subtitle = paste0("Incidence: ", inc_label),
      x = "Sampling Design",
      y = "MSE Component",
      fill = ""
    ) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==============================================================================
# ANALYSIS 5: Coverage vs MSE Tradeoff Scatter
#
# Each point is one scenario. X = MSE, Y = Coverage. Ideal designs cluster
# in the bottom-left (low MSE) and top (high coverage) region. Designs that
# reduce MSE at the cost of coverage are dangerous for inference.
# ==============================================================================

plot_coverage_mse_tradeoff <- function(results, inc_label = "") {
  ggplot(results, aes(x = MSE, y = Coverage, color = Design)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_color_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Coverage vs MSE Tradeoff",
      subtitle = paste0("Incidence: ", inc_label,
                        " | Each point = 1 scenario | Red dashed = nominal 95%"),
      x = "Mean Squared Error (MSE)",
      y = "95% CI Coverage Probability"
    ) +
    theme(legend.position = "right")
}

# ==============================================================================
# ANALYSIS 6: Parameter Sensitivity — which parameters drive MSE variation?
#
# For each design, compute the ratio of between-group variance to total
# variance when grouping by each parameter. High values mean that parameter
# strongly affects that design's performance.
# ==============================================================================

table_sensitivity <- function(results, inc_label = "") {
  params <- c("Rho", "Gamma", "Spillover_Type", "Neighbor_Type")
  designs <- sort(unique(results$Design))

  sens_rows <- list()
  idx <- 1

  for (d in designs) {
    d_data <- results %>% filter(Design == d)
    grand_mean <- mean(d_data$MSE, na.rm = TRUE)
    ss_total <- sum((d_data$MSE - grand_mean)^2, na.rm = TRUE)

    for (p in params) {
      # Proper eta-squared: SS_between / SS_total
      group_stats <- d_data %>%
        group_by(.data[[p]]) %>%
        summarise(grp_mean = mean(MSE, na.rm = TRUE),
                  grp_n = n(), .groups = "drop")
      ss_between <- sum(group_stats$grp_n * (group_stats$grp_mean - grand_mean)^2)
      eta_sq <- if (ss_total > 0) ss_between / ss_total else 0

      sens_rows[[idx]] <- data.frame(
        Design = d, Parameter = p, Eta_Squared = round(eta_sq, 4),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  sens_df <- bind_rows(sens_rows)

  # Pivot to wide: rows = Design, cols = Parameter
  sens_wide <- sens_df %>%
    pivot_wider(names_from = Parameter, values_from = Eta_Squared)

  cat(sprintf("\n=== Parameter Sensitivity (eta-squared) | Incidence: %s ===\n", inc_label))
  cat("(Higher = that parameter explains more MSE variation for that design)\n")
  cat("(Eta-sq = SS_between / SS_total; range [0, 1])\n\n")
  print(as.data.frame(sens_wide), row.names = FALSE, digits = 4)

  invisible(sens_wide)
}

# ==============================================================================
# TABLE: Flexible stratified summary
#
# Group by ANY combination of simulation parameters. Each call produces
# one table scoped to a single incidence config (inc_label) and stratified
# by the columns listed in `group_by_vars`.
#
# Valid stratification columns (within a single incidence config):
#   Design, Neighbor_Type, Rho, Gamma, Spillover_Type
#
# When working with combined results (all incidence configs together), you
# can also include: Incidence_Mode, Rho_Incidence
#
# Examples:
#   table_stratified(results, "Design")
#   table_stratified(results, c("Design", "Rho"))
#   table_stratified(results, c("Design", "Spillover_Type"))
#   table_stratified(results, c("Neighbor_Type", "Design", "Gamma"))
#   table_stratified(results, c("Design", "Rho", "Gamma", "Spillover_Type"))
# ==============================================================================

table_stratified <- function(results, group_by_vars = "Design",
                             inc_label = "", sort_by = "Avg_MSE",
                             max_rows = NULL) {
  # Validate requested grouping columns exist
  valid_cols <- intersect(group_by_vars, names(results))
  if (length(valid_cols) == 0) {
    cat("Warning: None of the requested grouping variables found in results.\n")
    return(invisible(NULL))
  }
  if (length(valid_cols) < length(group_by_vars)) {
    missing <- setdiff(group_by_vars, valid_cols)
    cat(sprintf("Note: Columns not in data (skipped): %s\n",
                paste(missing, collapse = ", ")))
  }

  summary_tbl <- results %>%
    group_by(across(all_of(valid_cols))) %>%
    summarise(
      N_Scenarios  = n(),
      Avg_Bias     = mean(Bias, na.rm = TRUE),
      Avg_SD       = mean(SD, na.rm = TRUE),
      Avg_MSE      = mean(MSE, na.rm = TRUE),
      Avg_Coverage = mean(Coverage, na.rm = TRUE),
      Avg_Fail     = mean(Fail_Rate, na.rm = TRUE),
      .groups = "drop"
    )

  # Sort
  if (sort_by %in% names(summary_tbl)) {
    summary_tbl <- summary_tbl %>% arrange(.data[[sort_by]])
  }

  # Header
  grouping_desc <- paste(valid_cols, collapse = " x ")
  header <- sprintf("\n=== Stratified by: %s", grouping_desc)
  if (nchar(inc_label) > 0) header <- paste0(header, " | Incidence: ", inc_label)
  header <- paste0(header, " ===\n")
  cat(header)
  cat(sprintf("(%d groups, sorted by %s)\n\n", nrow(summary_tbl), sort_by))

  out <- as.data.frame(summary_tbl)
  if (!is.null(max_rows) && nrow(out) > max_rows) {
    cat(sprintf("(Showing first %d of %d rows)\n\n", max_rows, nrow(out)))
    print(head(out, max_rows), row.names = FALSE, digits = 4)
  } else {
    print(out, row.names = FALSE, digits = 4)
  }

  invisible(summary_tbl)
}

# ==============================================================================
# TABLE: Full scenario-level results (per incidence config)
# ==============================================================================

table_comprehensive <- function(results, inc_label = "", max_rows = 30) {
  comprehensive <- results %>%
    arrange(Design, Neighbor_Type, Spillover_Type, Rho, Gamma) %>%
    select(Neighbor_Type, Design, Spillover_Type, Rho, Gamma,
           Bias, SD, MSE, Coverage, Fail_Rate)

  cat(sprintf("\n=== Full Scenario Results | Incidence: %s (%d rows) ===\n",
              inc_label, nrow(comprehensive)))

  out <- as.data.frame(comprehensive)
  if (nrow(out) > max_rows) {
    cat(sprintf("(Showing first %d of %d rows)\n\n", max_rows, nrow(out)))
    print(head(out, max_rows), row.names = FALSE, digits = 4)
  } else {
    print(out, row.names = FALSE, digits = 4)
  }

  invisible(comprehensive)
}

# ==============================================================================
# CONVENIENCE: Run standard set of stratified tables for one incidence config
# ==============================================================================

run_standard_tables <- function(results, inc_label = "") {
  # --- Design performance analyses ---
  # 1. Design rank frequency ("win rates")
  table_design_ranks(results, inc_label)

  # 2. Robustness profile (best / median / worst MSE)
  table_robustness(results, inc_label)

  # 3. Pairwise dominance matrix
  table_pairwise_dominance(results, inc_label)

  # 4. Parameter sensitivity (which params matter most for each design)
  table_sensitivity(results, inc_label)

  # --- Stratified summaries ---
  # 5. Overall by Design (the primary comparison)
  table_stratified(results, "Design", inc_label)

  # 6. Design x Spatial Dependence
  table_stratified(results, c("Design", "Rho"), inc_label)

  # 7. Design x Spillover Magnitude
  table_stratified(results, c("Design", "Gamma"), inc_label)

  # 8. Design x Spillover Type
  table_stratified(results, c("Design", "Spillover_Type"), inc_label)

  # 9. Design x Neighbor Type
  table_stratified(results, c("Design", "Neighbor_Type"), inc_label)

  # 10. Full scenario detail
  table_comprehensive(results, inc_label)
}

# ==============================================================================
# MAIN: Run all visualizations — per incidence config, never aggregated
# ==============================================================================

run_all_visualizations <- function(results = NULL, results_dir = NULL,
                                   estimation_mode = NULL,
                                   output_pdf = TRUE) {
  if (is.null(results)) {
    rd <- if (!is.null(results_dir)) results_dir else file.path(dirname(script_dir), "results")
    results <- load_latest_results(rd, estimation_mode)
  }

  est_label <- if (!is.null(estimation_mode)) estimation_mode else "SIM"

  # Split results by incidence config
  config_list <- split_by_incidence_config(results)

  cat(sprintf("\n=== Generating Visualizations for %d Incidence Configs ===\n",
              length(config_list)))
  cat("Configs:", paste(names(config_list), collapse = ", "), "\n\n")

  # --- Standalone plots (comparing across modes, produced once) ---
  if (output_pdf) {
    overview_file <- file.path(script_dir, "results",
                               sprintf("%s_incidence_overview.pdf", est_label))
    pdf(overview_file, width = 14, height = 10)
  }

  cat("Incidence overview: Grid Heatmaps...\n")
  print(plot_incidence_heatmaps())

  cat("Incidence overview: Distribution Histograms...\n")
  print(plot_incidence_distributions())

  if (output_pdf) {
    dev.off()
    cat(sprintf("Saved: %s\n\n", basename(overview_file)))
  }

  # --- Per-config plots and tables ---
  for (config_name in names(config_list)) {
    config_results <- config_list[[config_name]]

    # Create a filesystem-safe name for the PDF
    safe_name <- gsub("[^a-zA-Z0-9_]", "_", config_name)
    safe_name <- gsub("_+", "_", safe_name)
    safe_name <- gsub("_$", "", safe_name)

    cat(sprintf("--- [%s] ---\n", config_name))

    if (output_pdf) {
      pdf_file <- file.path(script_dir, "results",
                            sprintf("%s_%s.pdf", est_label, safe_name))
      pdf(pdf_file, width = 14, height = 10)
    }

    cat("  Plot 1: MSE by Neighbor Type...\n")
    print(plot_mse_by_neighbor(config_results, config_name))

    cat("  Plot 2: MSE Heatmap...\n")
    print(plot_mse_heatmap(config_results, config_name))

    cat("  Plot 3: Per-Design MSE Dynamics...\n")
    plot_mse_per_design(config_results, config_name)

    cat("  Plot 4a: Master Comparison (Queen)...\n")
    print(plot_master_comparison(config_results, "queen", config_name))

    cat("  Plot 4b: Master Comparison (Rook)...\n")
    print(plot_master_comparison(config_results, "rook", config_name))

    cat("  Plot 5: Coverage by Design...\n")
    p_cov <- plot_coverage_by_design(config_results, config_name)
    if (!is.null(p_cov)) print(p_cov)

    cat("  Plot 6: Design Rank Distribution...\n")
    print(plot_design_ranks(config_results, config_name))

    cat("  Plot 7: Bias-Variance Decomposition...\n")
    print(plot_bias_variance(config_results, config_name))

    cat("  Plot 8: Coverage vs MSE Tradeoff...\n")
    print(plot_coverage_mse_tradeoff(config_results, config_name))

    if (output_pdf) {
      dev.off()
      cat(sprintf("  Saved: %s\n", basename(pdf_file)))
    }

    # Tables — stratified by each simulation parameter
    run_standard_tables(config_results, config_name)

    cat("\n")
  }

  cat("=== All visualizations complete ===\n")
  invisible(results)
}
