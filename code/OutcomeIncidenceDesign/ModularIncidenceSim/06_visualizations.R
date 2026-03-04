# ==============================================================================
# 06_visualizations.R
# All plots and tables for the heterogeneous incidence simulation study.
# Can be run independently after 05_run_simulation.R has saved results.
#
# Usage:
#   source("06_visualizations.R")
#   -- or load results manually and call individual plot functions --
# ==============================================================================

library(ggplot2)
library(dplyr)
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
# LOAD RESULTS
# ==============================================================================

#' Load the most recent results file from the results directory
#' @param results_dir Path to results directory
#' @param estimation_mode Optional filter for "MLE" or "DIM"
#' @return Data frame of results
load_latest_results <- function(results_dir = file.path(script_dir, "results"),
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

plot_mse_by_neighbor <- function(results) {
  ggplot(results, aes(x = reorder(Design, MSE, FUN = median), y = MSE, fill = Design)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    facet_wrap(~ Neighbor_Type,
               labeller = as_labeller(c("queen" = "Queen Contiguity",
                                        "rook" = "Rook Contiguity"))) +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Distribution of Estimation Error (MSE) by Spatial Contiguity",
      subtitle = "Aggregated across all spatial dependence and spillover scenarios",
      x = "Sampling Design",
      y = "Mean Squared Error (Lower is Better)"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==============================================================================
# PLOT 2: MSE Heatmap by Rho x Design
# ==============================================================================

plot_mse_heatmap <- function(results) {
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
      x = expression("Spatial Dependence (" * rho * ")"),
      y = "Sampling Design",
      fill = "Avg MSE"
    ) +
    theme(panel.grid = element_blank())
}

# ==============================================================================
# PLOT 3: Per-Design MSE Line Plots across Gamma
# ==============================================================================

plot_mse_per_design <- function(results) {
  unique_designs <- sort(unique(results$Design))
  plots <- list()

  for (current_design in unique_designs) {
    design_data <- results %>% filter(Design == current_design)

    p <- ggplot(design_data,
                aes(x = as.factor(Gamma), y = MSE,
                    color = Spillover_Type, group = Spillover_Type)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      facet_grid(Neighbor_Type ~ paste("\u03c1 =", Rho)) +
      scale_color_viridis_d(option = "turbo",
                            labels = c("both" = "Treated & Control",
                                       "control_only" = "Control Only")) +
      theme_bw(base_size = 14) +
      labs(
        title = paste("MSE Dynamics for", current_design),
        subtitle = expression("Stratified by Spatial Dependence (" * rho * ") and Contiguity Type"),
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

plot_master_comparison <- function(results, nb_filter = "queen") {
  comparison_data <- results %>% filter(Neighbor_Type == nb_filter)

  ggplot(comparison_data, aes(x = Design, y = MSE, fill = as.factor(Gamma))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.2) +
    facet_grid(
      Spillover_Type ~ paste("\u03c1 =", Rho),
      scales = "free_y",
      labeller = labeller(Spillover_Type = c("both" = "Spillover: Both",
                                             "control_only" = "Spillover: Control"))
    ) +
    scale_fill_viridis_d(option = "turbo",
                         name = expression("Spillover\nMagnitude (" * gamma * ")")) +
    theme_bw(base_size = 13) +
    labs(
      title = "Design Comparison: Mean Squared Error across Scenarios",
      subtitle = paste0(tools::toTitleCase(nb_filter), " Contiguity Network"),
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
# PLOT 5 (NEW): MSE by Design, faceted by Incidence Mode
# ==============================================================================

plot_mse_by_incidence_mode <- function(results) {
  # Create readable incidence labels
  results_labeled <- results %>%
    mutate(Inc_Label = case_when(
      Incidence_Mode == "iid" ~ "iid Uniform",
      Incidence_Mode == "spatial" ~ paste0("Spatial (rho_X=", Rho_Incidence, ")"),
      Incidence_Mode == "poisson" ~ paste0("Poisson (rho_X=", Rho_Incidence, ")"),
      TRUE ~ Incidence_Mode
    ))

  ggplot(results_labeled,
         aes(x = reorder(Design, MSE, FUN = median), y = MSE, fill = Design)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.4) +
    facet_wrap(~ Inc_Label, scales = "free_y") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = "MSE by Design across Incidence Generation Modes",
      x = "Sampling Design",
      y = "Mean Squared Error"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==============================================================================
# PLOT 6 (NEW): Grid Heatmaps of Incidence Patterns
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
    data.frame(coords, Incidence = x_spatial_low[, 1],    Mode = "Spatial (rho_X=0.20)"),
    data.frame(coords, Incidence = x_spatial_high[, 1],   Mode = "Spatial (rho_X=0.50)"),
    data.frame(coords, Incidence = x_poisson_low[, 1],    Mode = "Poisson (rho_X=0.20)"),
    data.frame(coords, Incidence = x_poisson_high[, 1],   Mode = "Poisson (rho_X=0.50)")
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
# PLOT 7 (NEW): Incidence Distribution Histograms
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
    Mode = rep(c("iid Uniform", "Spatial (rho_X=0.20)", "Spatial (rho_X=0.50)",
                 "Poisson (rho_X=0.20)", "Poisson (rho_X=0.50)"),
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
# TABLE 1 (NEW): Coverage Probability by Design and Scenario
# ==============================================================================

table_coverage <- function(results) {
  if (!"Coverage" %in% names(results) || all(is.na(results$Coverage))) {
    cat("Coverage data not available in results.\n")
    return(invisible(NULL))
  }

  # Summary by Design x Incidence Mode
  cov_summary <- results %>%
    filter(!is.na(Coverage)) %>%
    group_by(Design, Incidence_Mode) %>%
    summarise(
      Mean_Coverage = mean(Coverage, na.rm = TRUE),
      Min_Coverage  = min(Coverage, na.rm = TRUE),
      Max_Coverage  = max(Coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Design, Incidence_Mode)

  cat("\n=== Coverage Probability Summary ===\n")
  cat("(Nominal = 0.95)\n\n")
  print(as.data.frame(cov_summary), row.names = FALSE, digits = 3)

  invisible(cov_summary)
}

# ==============================================================================
# TABLE 2: Summary by Design
# ==============================================================================

table_summary_by_design <- function(results) {
  summary_tbl <- results %>%
    group_by(Design) %>%
    summarise(
      Avg_Bias     = mean(Bias, na.rm = TRUE),
      Avg_SD       = mean(SD, na.rm = TRUE),
      Avg_MSE      = mean(MSE, na.rm = TRUE),
      Avg_Coverage = mean(Coverage, na.rm = TRUE),
      Avg_Fail     = mean(Fail_Rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Avg_MSE)

  cat("\n=== Summary by Design (sorted by Avg MSE) ===\n\n")
  print(as.data.frame(summary_tbl), row.names = FALSE, digits = 4)

  invisible(summary_tbl)
}

# ==============================================================================
# TABLE 3: Comprehensive Results
# ==============================================================================

table_comprehensive <- function(results) {
  comprehensive <- results %>%
    arrange(Design, Incidence_Mode, Rho_Incidence,
            Neighbor_Type, Spillover_Type, Rho, Gamma) %>%
    select(Incidence_Mode, Rho_Incidence, Neighbor_Type, Design,
           Spillover_Type, Rho, Gamma, Bias, SD, MSE, Coverage, Fail_Rate)

  cat(sprintf("\n=== Comprehensive Results (%d rows) ===\n", nrow(comprehensive)))
  cat("(Showing first 20 rows)\n\n")
  print(head(as.data.frame(comprehensive), 20), row.names = FALSE, digits = 4)

  invisible(comprehensive)
}

# ==============================================================================
# MAIN: Run all visualizations
# ==============================================================================

run_all_visualizations <- function(results = NULL, results_dir = NULL,
                                   estimation_mode = NULL) {
  if (is.null(results)) {
    rd <- if (!is.null(results_dir)) results_dir else file.path(script_dir, "results")
    results <- load_latest_results(rd, estimation_mode)
  }

  cat("\n=== Generating Visualizations ===\n\n")

  cat("Plot 1: MSE by Neighbor Type...\n")
  print(plot_mse_by_neighbor(results))

  cat("Plot 2: MSE Heatmap...\n")
  print(plot_mse_heatmap(results))

  cat("Plot 3: Per-Design MSE Dynamics...\n")
  plot_mse_per_design(results)

  cat("Plot 4: Master Comparison (Queen)...\n")
  print(plot_master_comparison(results, "queen"))

  cat("Plot 5: MSE by Incidence Mode...\n")
  print(plot_mse_by_incidence_mode(results))

  cat("Plot 6: Incidence Grid Heatmaps...\n")
  print(plot_incidence_heatmaps())

  cat("Plot 7: Incidence Distributions...\n")
  print(plot_incidence_distributions())

  cat("\n")
  table_summary_by_design(results)
  table_coverage(results)
  table_comprehensive(results)

  cat("\n=== All visualizations complete ===\n")
  invisible(results)
}
