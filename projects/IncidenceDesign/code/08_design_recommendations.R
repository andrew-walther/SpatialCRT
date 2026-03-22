# ==============================================================================
# 08_design_recommendations.R
# Personalized design recommendation system for the IncidenceDesign simulation.
#
# Answers three questions:
#   Q1: Which design is best per incidence mode?
#   Q2: Which design is best at each level of a single parameter (marginal)?
#   Q3: Which design is best for a specific parameter combination?
#
# Sources 06_visualizations.R (which sources 01-03).
#
# IMPORTANT: Results for each incidence generation method are reported
#   SEPARATELY. They are never aggregated across incidence modes.
#
# Usage:
#   source("08_design_recommendations.R")
#   run_recommendation_report(estimation_mode = "MLE_combined")
#   -- or call individual functions --
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# --- Source 06 (which sources 01-03) ---
# Detect script directory. When sourced interactively, sys.frame(1)$ofile gives
# the file path. When run via Rscript from the code/ directory, getwd() is code/.
# Normalize so script_dir_07 always points to the code/ directory as an absolute path.
script_dir_07 <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
source(file.path(script_dir_07, "06_visualizations.R"))
# Override script_dir (set by 06) with our normalized version
script_dir <- script_dir_07

# ==============================================================================
# LOCAL HELPERS
# ==============================================================================

#' Abbreviated design names for compact plot labels
#' @param design_str Character vector of design names (e.g., "Design 3")
#' @return Character vector of short names
short_design_name <- function(design_str) {
  mapping <- c(
    "Design 1" = "D1: Checkerboard",
    "Design 2" = "D2: High Inc. Focus",
    "Design 3" = "D3: Sat. Quadrants",
    "Design 4" = "D4: Isolation Buffer",
    "Design 5" = "D5: 2x2 Blocking",
    "Design 6" = "D6: Balanced Quartiles",
    "Design 7" = "D7: Balanced Halves",
    "Design 8" = "D8: Inc.-Guided Sat. Q."
  )
  ifelse(design_str %in% names(mapping), mapping[design_str], design_str)
}

#' Reorder factor within faceting groups (avoids tidytext dependency)
#' @param x Character vector to reorder
#' @param by Numeric vector to reorder by
#' @param within Character vector defining the faceting groups
#' @param fun Aggregation function (default: mean)
#' @param sep Internal separator string
#' @return Factor
reorder_within_local <- function(x, by, within, fun = mean, sep = "___") {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

#' Scale for reorder_within_local — strips separator suffix from labels
#' @param sep Internal separator string
#' @param ... Passed to scale_x_discrete
#' @return ggplot scale
scale_x_reordered_local <- function(sep = "___", ...) {
  scale_x_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# ==============================================================================
# CORE UTILITY
# ==============================================================================

#' Rank designs by average metric within groups
#'
#' The workhorse function: computes avg MSE, Coverage, Bias, SD per design
#' within each group defined by group_vars, then assigns ranks.
#'
#' @param results Data frame for one incidence config
#' @param group_vars Character vector of columns to group by (must NOT include "Design")
#' @param metric Character: "MSE" (default) — which metric to rank by (ascending)
#' @return Data frame with columns: [group_vars], Design, Avg_MSE, Avg_Coverage,
#'         Avg_Bias, Avg_SD, N_Scenarios, Rank
rank_designs_by_group <- function(results, group_vars = character(0), metric = "MSE") {
  all_group_vars <- c(group_vars, "Design")

  summary_tbl <- results %>%
    group_by(across(all_of(all_group_vars))) %>%
    summarise(
      Avg_MSE      = mean(MSE, na.rm = TRUE),
      Avg_Coverage = mean(Coverage, na.rm = TRUE),
      Avg_Bias     = mean(Bias, na.rm = TRUE),
      Avg_SD       = mean(SD, na.rm = TRUE),
      N_Scenarios  = n(),
      .groups = "drop"
    )

  rank_col <- paste0("Avg_", metric)

  if (length(group_vars) > 0) {
    summary_tbl <- summary_tbl %>%
      group_by(across(all_of(group_vars))) %>%
      mutate(Rank = rank(.data[[rank_col]], ties.method = "min")) %>%
      ungroup()
  } else {
    summary_tbl <- summary_tbl %>%
      mutate(Rank = rank(.data[[rank_col]], ties.method = "min"))
  }

  summary_tbl %>% arrange(across(all_of(c(group_vars, "Rank"))))
}

# ==============================================================================
# QUESTION 1: PER INCIDENCE MODE
# ==============================================================================

#' Table: Design rankings for each incidence config
#'
#' For each of the 5 incidence configs, ranks all 6 designs by avg MSE.
#' Shows avg MSE and avg Coverage side-by-side.
#'
#' @param results Combined results data frame (all 1,920 rows)
#' @return Invisible data frame with columns: Inc_Config, Design, Avg_MSE,
#'         Avg_Coverage, Avg_Bias, Rank
table_incidence_rankings <- function(results) {
  config_list <- split_by_incidence_config(results)
  all_ranked <- list()

  cat("\n==========================================================\n")
  cat("  DESIGN RANKINGS BY INCIDENCE MODE (MLE)\n")
  cat("  Ranked by average MSE; lower rank = better.\n")
  cat("==========================================================\n")

  for (config_name in names(config_list)) {
    config_results <- config_list[[config_name]]
    ranked <- rank_designs_by_group(config_results, character(0))
    ranked$Inc_Config <- config_name
    all_ranked[[config_name]] <- ranked

    cat(sprintf("\n--- %s ---\n", config_name))
    cat(sprintf("  %-4s  %-28s  %8s  %8s  %8s\n",
                "Rank", "Design", "Avg_MSE", "Coverage", "Bias"))
    for (i in seq_len(nrow(ranked))) {
      cat(sprintf("  %-4d  %-28s  %8.4f  %8.3f  %8.4f\n",
                  ranked$Rank[i],
                  short_design_name(ranked$Design[i]),
                  ranked$Avg_MSE[i],
                  ranked$Avg_Coverage[i],
                  ranked$Avg_Bias[i]))
    }
  }

  combined <- bind_rows(all_ranked)
  invisible(combined)
}

#' Plot: Faceted bar chart of avg MSE per design, one panel per incidence config
#'
#' @param results Combined results data frame (all 1,920 rows)
#' @return ggplot object
plot_incidence_rankings <- function(results) {
  config_list <- split_by_incidence_config(results)
  all_ranked <- list()

  for (config_name in names(config_list)) {
    ranked <- rank_designs_by_group(config_list[[config_name]], character(0))
    ranked$Inc_Config <- config_name
    all_ranked[[config_name]] <- ranked
  }

  plot_df <- bind_rows(all_ranked)
  plot_df$Short_Name <- short_design_name(plot_df$Design)

  ggplot(plot_df, aes(x = reorder_within_local(Short_Name, Avg_MSE, Inc_Config),
                      y = Avg_MSE, fill = Design)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.85) +
    geom_text(aes(label = sprintf("#%d", Rank)), vjust = -0.3,
              size = 3.5, fontface = "bold") +
    facet_wrap(~ Inc_Config, scales = "free_x", ncol = 3) +
    scale_fill_viridis_d(option = "turbo") +
    scale_x_reordered_local() +
    theme_minimal(base_size = 14) +
    labs(
      title = "Design Recommendation: Average MSE by Incidence Mode",
      subtitle = "Rank labels shown above bars (1 = best). MLE estimator.",
      x = "Sampling Design",
      y = "Average MSE (lower is better)"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
}

#' Plot: Faceted bar chart of avg Coverage per design, per incidence config
#'
#' @param results Combined results data frame (all 1,920 rows)
#' @return ggplot object
plot_incidence_coverage <- function(results) {
  config_list <- split_by_incidence_config(results)
  all_ranked <- list()

  for (config_name in names(config_list)) {
    ranked <- rank_designs_by_group(config_list[[config_name]], character(0))
    ranked$Inc_Config <- config_name
    all_ranked[[config_name]] <- ranked
  }

  plot_df <- bind_rows(all_ranked)
  plot_df$Short_Name <- short_design_name(plot_df$Design)

  ggplot(plot_df, aes(x = reorder_within_local(Short_Name, -Avg_Coverage, Inc_Config),
                      y = Avg_Coverage, fill = Design)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.85) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.6) +
    facet_wrap(~ Inc_Config, scales = "free_x", ncol = 3) +
    scale_fill_viridis_d(option = "turbo") +
    scale_x_reordered_local() +
    theme_minimal(base_size = 14) +
    labs(
      title = "Design Recommendation: Average Coverage by Incidence Mode",
      subtitle = "Red dashed = nominal 95%. MLE estimator.",
      x = "Sampling Design",
      y = "Average 95% CI Coverage"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold")) +
    ylim(0, 1)
}

# ==============================================================================
# QUESTION 2: PER INDIVIDUAL PARAMETER (MARGINAL)
# ==============================================================================

#' Table: Design rankings at each level of a single parameter
#'
#' For a given parameter (Rho, Gamma, Spillover_Type, Neighbor_Type), shows
#' the avg MSE for each design at each parameter level. Marginalizes over
#' all other parameters.
#'
#' @param results Data frame for one incidence config
#' @param param Character: column name to condition on
#' @param inc_label Character: incidence config label
#' @return Invisible data frame (wide format)
table_marginal_rankings <- function(results, param, inc_label = "") {
  ranked <- rank_designs_by_group(results, param)

  # Pivot to wide: rows = Design, cols = param levels, cells = avg MSE
  wide <- ranked %>%
    select(Design, all_of(param), Avg_MSE, Rank) %>%
    mutate(cell = sprintf("%.4f (#%d)", Avg_MSE, Rank)) %>%
    select(Design, all_of(param), cell) %>%
    pivot_wider(names_from = all_of(param), values_from = cell,
                names_prefix = paste0(param, "="))

  # Sort by overall rank (rank designs by their avg MSE across all levels)
  overall <- rank_designs_by_group(results, character(0))
  design_order <- overall %>% arrange(Rank) %>% pull(Design)
  wide <- wide %>%
    mutate(Design = factor(Design, levels = design_order)) %>%
    arrange(Design) %>%
    mutate(Design = as.character(Design))

  cat(sprintf("\n=== Marginal Rankings by %s | Incidence: %s ===\n", param, inc_label))
  cat("(Each cell = Avg_MSE (#Rank) at that parameter level, marginalized over all others)\n\n")
  print(as.data.frame(wide), row.names = FALSE)

  invisible(ranked)
}

#' Plot: Design rank trajectories across parameter levels
#'
#' Line plot with x = parameter level, y = rank (1 at top), color = design.
#' Shows how design rankings shift as a parameter changes.
#'
#' @param results Data frame for one incidence config
#' @param param Character: parameter to vary on x-axis
#' @param inc_label Character: incidence config label
#' @return ggplot object
plot_rank_trajectories <- function(results, param, inc_label = "",
                                   point_size = 3.5, line_size = 1.2,
                                   label_size = 2.8) {
  ranked <- rank_designs_by_group(results, param)
  ranked$Short_Name <- short_design_name(ranked$Design)
  ranked[[param]] <- factor(ranked[[param]])

  ggplot(ranked, aes(x = .data[[param]], y = Rank,
                     color = Design, group = Design)) +
    geom_line(linewidth = line_size, alpha = 0.8) +
    geom_point(size = point_size) +
    geom_text(aes(label = Short_Name), hjust = -0.1, vjust = -0.5,
              size = label_size, show.legend = FALSE) +
    scale_y_reverse(breaks = 1:6, limits = c(6.5, 0.5)) +
    scale_color_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Design Rank Trajectories across %s", param),
      subtitle = paste0("Incidence: ", inc_label, " | Rank 1 = lowest avg MSE"),
      x = param,
      y = "Rank (1 = Best)"
    ) +
    theme(legend.position = "right")
}

#' Plot: Conditional MSE boxplots by parameter level
#'
#' Faceted boxplot showing MSE distribution per design at each level of a
#' parameter. Reveals not just averages but full distributional shifts.
#'
#' @param results Data frame for one incidence config
#' @param param Character: parameter to condition on
#' @param inc_label Character: incidence config label
#' @return ggplot object
plot_conditional_mse <- function(results, param, inc_label = "") {
  results$param_facet <- paste0(param, " = ", results[[param]])

  ggplot(results, aes(x = reorder(Design, MSE, FUN = median),
                      y = MSE, fill = Design)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.4) +
    facet_wrap(~ param_facet, scales = "free_y") +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("MSE Distribution by Design, Conditioned on %s", param),
      subtitle = paste0("Incidence: ", inc_label, " | Designs sorted by median MSE"),
      x = "Sampling Design",
      y = "MSE"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
}

# ==============================================================================
# QUESTION 3: PER PARAMETER COMBINATION
# ==============================================================================

#' Plot: Heatmap showing the best design for each (param1, param2) combination
#'
#' Each tile shows which design achieves the lowest avg MSE at that
#' parameter combination. Color = design, text = design name + MSE.
#'
#' @param results Data frame for one incidence config
#' @param row_param Character: parameter for y-axis (e.g., "Rho")
#' @param col_param Character: parameter for x-axis (e.g., "Gamma")
#' @param inc_label Character: incidence config label
#' @return ggplot object
plot_best_design_heatmap <- function(results, row_param = "Rho",
                                     col_param = "Gamma", inc_label = "",
                                     text_size = 3, text_color = "white",
                                     compact_labels = FALSE) {
  ranked <- rank_designs_by_group(results, c(row_param, col_param))
  best <- ranked %>% filter(Rank == 1)
  best$Short_Name <- short_design_name(best$Design)

  # compact_labels = TRUE: "D3\n0.060"; FALSE: "D3: Sat. Quadrants\nMSE=0.060"
  if (compact_labels) {
    best$Tile_Label <- sprintf("%s\n%.3f",
                               gsub("Design ", "D", best$Design), best$Avg_MSE)
  } else {
    best$Tile_Label <- sprintf("%s\nMSE=%.3f", best$Short_Name, best$Avg_MSE)
  }

  ggplot(best, aes(x = factor(.data[[col_param]]),
                   y = factor(.data[[row_param]]),
                   fill = Design)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = Tile_Label),
              size = text_size, color = text_color,
              fontface = "bold", lineheight = 0.9) +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Best Design by %s x %s", row_param, col_param),
      subtitle = paste0("Incidence: ", inc_label,
                        " | Cell = best design (lowest avg MSE across other params)"),
      x = col_param,
      y = row_param
    ) +
    theme(panel.grid = element_blank(),
          legend.position = "right")
}

#' Table: Scenario-specific design ranking lookup
#'
#' Given specific parameter values, filters the results and ranks the 6 designs.
#' Prints a formatted recommendation with the winner and runner-up.
#'
#' @param results Data frame for one incidence config
#' @param rho Numeric: outcome spatial autocorrelation (optional)
#' @param gamma Numeric: spillover magnitude (optional)
#' @param spill_type Character: "control_only" or "both" (optional)
#' @param nb_type Character: "rook" or "queen" (optional)
#' @param inc_label Character: incidence config label
#' @return Invisible data frame of ranked designs
table_scenario_lookup <- function(results, rho = NULL, gamma = NULL,
                                  spill_type = NULL, nb_type = NULL,
                                  inc_label = "") {
  filtered <- results
  param_desc <- c()

  if (!is.null(rho)) {
    filtered <- filtered %>% filter(Rho == rho)
    param_desc <- c(param_desc, sprintf("Rho = %.2f", rho))
  }
  if (!is.null(gamma)) {
    filtered <- filtered %>% filter(Gamma == gamma)
    param_desc <- c(param_desc, sprintf("Gamma = %.1f", gamma))
  }
  if (!is.null(spill_type)) {
    filtered <- filtered %>% filter(Spillover_Type == spill_type)
    param_desc <- c(param_desc, sprintf("Spillover = %s", spill_type))
  }
  if (!is.null(nb_type)) {
    filtered <- filtered %>% filter(Neighbor_Type == nb_type)
    param_desc <- c(param_desc, sprintf("Neighbor = %s", nb_type))
  }

  if (nrow(filtered) == 0) {
    cat("No scenarios match the given parameters.\n")
    return(invisible(NULL))
  }

  # Rank designs on filtered data (no further grouping)
  ranked <- rank_designs_by_group(filtered, character(0))

  param_str <- if (length(param_desc) > 0) paste(param_desc, collapse = ", ") else "All"

  cat(sprintf("\n=== Scenario Lookup | Incidence: %s ===\n", inc_label))
  cat(sprintf("Parameters: %s\n", param_str))
  cat(sprintf("(%d matching scenarios per design)\n\n", ranked$N_Scenarios[1]))
  cat(sprintf("  %-4s  %-28s  %8s  %8s  %8s  %8s\n",
              "Rank", "Design", "Avg_MSE", "Coverage", "Bias", "SD"))
  for (i in seq_len(nrow(ranked))) {
    cat(sprintf("  %-4d  %-28s  %8.4f  %8.3f  %+8.4f  %8.4f\n",
                ranked$Rank[i],
                short_design_name(ranked$Design[i]),
                ranked$Avg_MSE[i],
                ranked$Avg_Coverage[i],
                ranked$Avg_Bias[i],
                ranked$Avg_SD[i]))
  }

  # Recommendation
  best <- ranked %>% filter(Rank == 1)
  second <- ranked %>% filter(Rank == 2)
  if (nrow(best) > 0 && nrow(second) > 0) {
    pct_diff <- (second$Avg_MSE[1] - best$Avg_MSE[1]) / best$Avg_MSE[1] * 100
    cat(sprintf("\n>> RECOMMENDATION: %s with MSE = %.4f, Coverage = %.3f\n",
                short_design_name(best$Design[1]), best$Avg_MSE[1], best$Avg_Coverage[1]))
    cat(sprintf(">> Runner-up: %s, MSE +%.1f%% higher\n",
                short_design_name(second$Design[1]), pct_diff))
  }

  invisible(ranked)
}

# ==============================================================================
# SUMMARY & COMMENTARY
# ==============================================================================

#' Generate programmatic commentary for one incidence config
#'
#' Computes and prints 6 key findings: overall winner, dominance %,
#' rank stability, coverage, parameter sensitivity, and practical recommendation.
#'
#' @param results Data frame for one incidence config
#' @param inc_label Character: incidence config label
#' @return Invisible character vector of commentary lines
generate_commentary <- function(results, inc_label = "") {
  lines <- character(0)
  add <- function(txt) { lines <<- c(lines, txt) }

  add(sprintf("\n=== Design Recommendation Commentary | Incidence: %s (MLE) ===\n",
              inc_label))

  # --- Finding 1: Overall winner ---
  overall <- rank_designs_by_group(results, character(0))
  best_row <- overall %>% filter(Rank == 1)
  worst_row <- overall %>% filter(Rank == max(Rank))
  best_name <- short_design_name(best_row$Design[1])
  worst_name <- short_design_name(worst_row$Design[1])

  add(sprintf("Finding 1: OVERALL WINNER"))
  add(sprintf("  %s achieves the lowest average MSE (%.4f) across",
              best_name, best_row$Avg_MSE[1]))
  add(sprintf("  all %d scenarios in this incidence config.\n",
              sum(overall$N_Scenarios) / nrow(overall)))

  # --- Finding 2: Dominance ---
  # Compute pairwise: best vs second, best vs worst
  scenario_ids <- c("Neighbor_Type", "Rho", "Gamma", "Spillover_Type")
  wide <- results %>%
    select(all_of(c(scenario_ids, "Design", "MSE"))) %>%
    pivot_wider(names_from = Design, values_from = MSE,
                id_cols = all_of(scenario_ids))

  second_row <- overall %>% filter(Rank == 2)
  best_design <- best_row$Design[1]
  second_design <- second_row$Design[1]
  worst_design <- worst_row$Design[1]

  dom_vs_second <- mean(wide[[best_design]] < wide[[second_design]], na.rm = TRUE) * 100
  dom_vs_worst <- mean(wide[[best_design]] < wide[[worst_design]], na.rm = TRUE) * 100

  add(sprintf("Finding 2: DOMINANCE"))
  add(sprintf("  %s beats the runner-up (%s) in %.1f%% of scenarios.",
              best_name, short_design_name(second_design), dom_vs_second))
  add(sprintf("  %s beats the worst performer (%s) in %.1f%% of scenarios.\n",
              best_name, worst_name, dom_vs_worst))

  # --- Finding 3: Rank stability ---
  rho_ranks <- rank_designs_by_group(results, "Rho") %>%
    filter(Design == best_design)
  gamma_ranks <- rank_designs_by_group(results, "Gamma") %>%
    filter(Design == best_design)

  rho_stable <- all(rho_ranks$Rank == 1)
  gamma_stable <- all(gamma_ranks$Rank == 1)

  add(sprintf("Finding 3: ROBUSTNESS"))
  if (rho_stable && gamma_stable) {
    add(sprintf("  %s holds rank 1 at ALL levels of Rho and ALL levels of Gamma.",
                best_name))
    add("  No rank crossover occurs -- this is the consistently best design.\n")
  } else {
    if (!rho_stable) {
      crossover_rho <- rho_ranks %>% filter(Rank != 1)
      add(sprintf("  %s loses rank 1 at Rho = %s (rank %d there).",
                  best_name,
                  paste(crossover_rho$Rho, collapse = ", "),
                  crossover_rho$Rank[1]))
    }
    if (!gamma_stable) {
      crossover_gamma <- gamma_ranks %>% filter(Rank != 1)
      add(sprintf("  %s loses rank 1 at Gamma = %s (rank %d there).",
                  best_name,
                  paste(crossover_gamma$Gamma, collapse = ", "),
                  crossover_gamma$Rank[1]))
    }
    add("")
  }

  # --- Finding 4: Coverage ---
  add(sprintf("Finding 4: COVERAGE"))
  add(sprintf("  %s achieves average coverage of %.3f (nominal = 0.95).",
              best_name, best_row$Avg_Coverage[1]))
  add(sprintf("  %s (worst) achieves coverage of %.3f.\n",
              worst_name, worst_row$Avg_Coverage[1]))

  # --- Finding 5: Sensitivity ---
  params <- c("Rho", "Gamma", "Spillover_Type", "Neighbor_Type")
  d_data <- results %>% filter(Design == best_design)
  grand_mean <- mean(d_data$MSE, na.rm = TRUE)
  ss_total <- sum((d_data$MSE - grand_mean)^2, na.rm = TRUE)

  eta_sq <- numeric(length(params))
  names(eta_sq) <- params
  for (p in params) {
    grp <- d_data %>%
      group_by(.data[[p]]) %>%
      summarise(grp_mean = mean(MSE, na.rm = TRUE), grp_n = n(), .groups = "drop")
    ss_b <- sum(grp$grp_n * (grp$grp_mean - grand_mean)^2)
    eta_sq[p] <- if (ss_total > 0) ss_b / ss_total else 0
  }
  most_influential <- names(which.max(eta_sq))

  add(sprintf("Finding 5: SENSITIVITY"))
  add(sprintf("  For %s, the most influential parameter is %s (eta-sq = %.3f),",
              best_name, most_influential, eta_sq[most_influential]))
  sorted_eta <- sort(eta_sq, decreasing = TRUE)
  if (length(sorted_eta) >= 2) {
    add(sprintf("  followed by %s (eta-sq = %.3f).",
                names(sorted_eta)[2], sorted_eta[2]))
  }
  low_impact <- names(sorted_eta[sorted_eta < 0.05])
  if (length(low_impact) > 0) {
    add(sprintf("  Low-impact parameters (eta-sq < 0.05): %s.\n",
                paste(low_impact, collapse = ", ")))
  } else {
    add("")
  }

  # --- Finding 6: Practical recommendation ---
  mse_ratio <- worst_row$Avg_MSE[1] / best_row$Avg_MSE[1]
  add(sprintf("Finding 6: PRACTICAL RECOMMENDATION"))
  add(sprintf("  Under %s incidence, use %s for treatment assignment.",
              inc_label, best_name))
  add(sprintf("  It offers the best MSE (%.1fx lower than worst performer),",
              mse_ratio))
  if (best_row$Avg_Coverage[1] >= 0.90) {
    add("  maintains near-nominal coverage, and is robust across parameter variations.")
  } else {
    add(sprintf("  though coverage (%.3f) is below nominal 0.95 -- interpret CIs cautiously.",
                best_row$Avg_Coverage[1]))
  }
  add("")

  cat(paste(lines, collapse = "\n"))
  invisible(lines)
}

#' Master orchestrator: generate the full recommendation report
#'
#' Produces a PDF with all recommendation plots and prints tables + commentary
#' to the console. Operates on MLE results by default.
#'
#' @param results Combined results (all 1,920 rows) or NULL to auto-load
#' @param results_dir Path to results directory (used if results is NULL)
#' @param estimation_mode Character: "MLE_combined" (default)
#' @param output_pdf Logical: save to PDF?
#' @return Invisible combined results data frame
run_recommendation_report <- function(results = NULL, results_dir = NULL,
                                      estimation_mode = "MLE_combined",
                                      output_pdf = TRUE) {
  # Load results if needed
  if (is.null(results)) {
    rd <- if (!is.null(results_dir)) {
      results_dir
    } else {
      file.path(dirname(script_dir_07), "results")
    }
    results <- load_latest_results(rd, estimation_mode)
  }

  est_label <- if (!is.null(estimation_mode)) estimation_mode else "MLE"

  config_list <- split_by_incidence_config(results)

  cat("\n##############################################################\n")
  cat("  PERSONALIZED DESIGN RECOMMENDATION REPORT\n")
  cat(sprintf("  Estimator: %s | Configs: %d | Scenarios: %d\n",
              est_label, length(config_list), nrow(results)))
  cat("##############################################################\n")

  # --- PDF setup ---
  if (output_pdf) {
    pdf_file <- file.path(dirname(script_dir_07), "results",
                          sprintf("%s_design_recommendations.pdf", est_label))
    pdf(pdf_file, width = 14, height = 10)
  }

  # --- Page 1: Cross-config MSE rankings ---
  cat("\n[Page 1] Cross-config MSE rankings...\n")
  print(plot_incidence_rankings(results))

  # --- Page 2: Cross-config Coverage ---
  cat("[Page 2] Cross-config Coverage...\n")
  print(plot_incidence_coverage(results))

  # --- Console: Cross-config table ---
  table_incidence_rankings(results)

  # --- Per-config analyses ---
  for (config_name in names(config_list)) {
    config_results <- config_list[[config_name]]

    cat(sprintf("\n\n--- CONFIG: %s ---\n", config_name))

    # Rank trajectories
    cat(sprintf("[Plot] Rank trajectories across Rho...\n"))
    print(plot_rank_trajectories(config_results, "Rho", config_name))

    cat(sprintf("[Plot] Rank trajectories across Gamma...\n"))
    print(plot_rank_trajectories(config_results, "Gamma", config_name))

    # Conditional MSE boxplots
    cat(sprintf("[Plot] Conditional MSE by Rho...\n"))
    print(plot_conditional_mse(config_results, "Rho", config_name))

    cat(sprintf("[Plot] Conditional MSE by Gamma...\n"))
    print(plot_conditional_mse(config_results, "Gamma", config_name))

    # Best design heatmaps
    cat(sprintf("[Plot] Best design: Rho x Gamma...\n"))
    print(plot_best_design_heatmap(config_results, "Rho", "Gamma", config_name))

    # Console: Marginal ranking tables
    table_marginal_rankings(config_results, "Rho", config_name)
    table_marginal_rankings(config_results, "Gamma", config_name)
    table_marginal_rankings(config_results, "Spillover_Type", config_name)
    table_marginal_rankings(config_results, "Neighbor_Type", config_name)

    # Console: Scenario lookup examples
    cat("\n--- Example scenario lookups ---\n")
    table_scenario_lookup(config_results, rho = 0.20, gamma = 0.7,
                          spill_type = "both", nb_type = "queen",
                          inc_label = config_name)
    table_scenario_lookup(config_results, rho = 0.50, gamma = 0.8,
                          spill_type = "control_only", nb_type = "rook",
                          inc_label = config_name)

    # Console: Commentary
    generate_commentary(config_results, config_name)
  }

  # --- Close PDF ---
  if (output_pdf) {
    dev.off()
    cat(sprintf("\n>> PDF saved: %s\n", basename(pdf_file)))
  }

  cat("\n##############################################################\n")
  cat("  RECOMMENDATION REPORT COMPLETE\n")
  cat("##############################################################\n")

  invisible(results)
}

# ==============================================================================
# VALIDATION
# ==============================================================================

#' Validate the recommendation functions for correctness
#'
#' Runs 8 assertion-based tests. Stops with error on first failure.
#'
#' @param results Combined results (all 1,920 rows) or NULL to auto-load
#' @return Invisible TRUE if all tests pass
validate_recommendations <- function(results = NULL) {
  if (is.null(results)) {
    results <- load_latest_results(
      results_dir = file.path(dirname(script_dir_07), "results"),
      estimation_mode = "MLE_combined"
    )
  }

  cat("\n=== Validation Suite for 08_design_recommendations.R ===\n\n")

  configs <- split_by_incidence_config(results)
  test_config <- configs[[1]]
  test_label <- names(configs)[1]

  # Test 1: rank_designs_by_group returns correct shape (no grouping)
  cat("Test 1: rank_designs_by_group shape (ungrouped)... ")
  ranked <- rank_designs_by_group(test_config, character(0))
  stopifnot(nrow(ranked) == 6)
  stopifnot(all(c("Design", "Avg_MSE", "Avg_Coverage", "Rank") %in% names(ranked)))
  stopifnot(all(ranked$Rank %in% 1:6))
  cat("PASS\n")

  # Test 2: rank_designs_by_group with group_vars
  cat("Test 2: rank_designs_by_group with grouping... ")
  ranked_rho <- rank_designs_by_group(test_config, "Rho")
  n_rho <- length(unique(test_config$Rho))
  stopifnot(nrow(ranked_rho) == 6 * n_rho)
  cat("PASS\n")

  # Test 3: Ranks are 1-6 within each group
  cat("Test 3: Rank range within groups... ")
  rank_check <- ranked_rho %>%
    group_by(Rho) %>%
    summarise(min_r = min(Rank), max_r = max(Rank), .groups = "drop")
  stopifnot(all(rank_check$min_r == 1))
  stopifnot(all(rank_check$max_r <= 6))
  cat("PASS\n")

  # Test 4: No cross-incidence aggregation in table_incidence_rankings
  cat("Test 4: No cross-incidence aggregation... ")
  all_ranked <- table_incidence_rankings(results)
  stopifnot("Inc_Config" %in% names(all_ranked))
  n_expected <- length(unique(paste(results$Incidence_Mode, results$Rho_Incidence)))
  stopifnot(length(unique(all_ranked$Inc_Config)) == n_expected)
  cat("PASS\n")

  # Test 5: MSE ordering matches Rank ordering
  cat("Test 5: MSE-Rank consistency... ")
  for (rho_val in unique(ranked_rho$Rho)) {
    sub <- ranked_rho %>% filter(Rho == rho_val) %>% arrange(Rank)
    diffs <- diff(sub$Avg_MSE)
    stopifnot(all(diffs >= -1e-10))
  }
  cat("PASS\n")

  # Test 6: table_scenario_lookup returns 6 rows for full specification
  cat("Test 6: Scenario lookup shape... ")
  lookup <- table_scenario_lookup(test_config, rho = 0.20, gamma = 0.7,
                                  spill_type = "both", nb_type = "queen",
                                  inc_label = test_label)
  stopifnot(nrow(lookup) == 6)
  cat("PASS\n")

  # Test 7: Plot functions return ggplot objects
  cat("Test 7: Plot return types... ")
  p1 <- plot_incidence_rankings(results)
  stopifnot(inherits(p1, "gg"))
  p2 <- plot_rank_trajectories(test_config, "Rho", test_label)
  stopifnot(inherits(p2, "gg"))
  p3 <- plot_best_design_heatmap(test_config, "Rho", "Gamma", test_label)
  stopifnot(inherits(p3, "gg"))
  p4 <- plot_incidence_coverage(results)
  stopifnot(inherits(p4, "gg"))
  p5 <- plot_conditional_mse(test_config, "Rho", test_label)
  stopifnot(inherits(p5, "gg"))
  cat("PASS\n")

  # Test 8: Commentary generation produces non-empty output
  cat("Test 8: Commentary generation... ")
  commentary <- generate_commentary(test_config, test_label)
  stopifnot(is.character(commentary))
  stopifnot(length(commentary) > 0)
  stopifnot(nchar(paste(commentary, collapse = "")) > 100)
  cat("PASS\n")

  cat("\n=== All 8 validation tests PASSED ===\n")
  invisible(TRUE)
}

#' Validate that sourcing 07 does not break existing module functions
#'
#' Checks that all functions from 06_visualizations.R remain callable and
#' produce expected output types after 07 is loaded.
#'
#' @param results Combined results (all 1,920 rows) or NULL to auto-load
#' @return Invisible TRUE if all tests pass
validate_no_side_effects <- function(results = NULL) {
  if (is.null(results)) {
    results <- load_latest_results(
      results_dir = file.path(dirname(script_dir_07), "results"),
      estimation_mode = "MLE_combined"
    )
  }

  cat("\n=== Side-Effect Validation: Existing Modules Unaffected ===\n\n")

  # Test 1: Core 06 functions still exist
  cat("Test 1: 06_visualizations.R functions exist... ")
  required_fns <- c("split_by_incidence_config", "load_latest_results",
                     "inc_config_label", "table_design_ranks",
                     "table_robustness", "table_pairwise_dominance",
                     "table_sensitivity", "table_stratified",
                     "plot_mse_by_neighbor", "plot_mse_heatmap",
                     "plot_coverage_by_design", "plot_design_ranks",
                     "plot_bias_variance", "plot_coverage_mse_tradeoff",
                     "run_all_visualizations")
  for (fn_name in required_fns) {
    stopifnot(exists(fn_name, mode = "function"))
  }
  cat("PASS\n")

  # Test 2: Results data integrity
  cat("Test 2: Results data integrity... ")
  stopifnot(nrow(results) == 1920)
  stopifnot(ncol(results) == 13)
  expected_cols <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design",
                     "Rho", "Gamma", "Spillover_Type", "Mean_Estimate",
                     "Bias", "SD", "MSE", "Coverage", "Fail_Rate")
  stopifnot(all(expected_cols %in% names(results)))
  cat("PASS\n")

  # Test 3: split_by_incidence_config still works correctly
  cat("Test 3: split_by_incidence_config... ")
  configs <- split_by_incidence_config(results)
  stopifnot(length(configs) == 5)
  stopifnot(all(sapply(configs, nrow) == 384))
  cat("PASS\n")

  # Test 4: table_design_ranks produces expected shape
  cat("Test 4: table_design_ranks output shape... ")
  test_results <- configs[[1]]
  rank_tbl <- table_design_ranks(test_results, names(configs)[1])
  stopifnot(nrow(rank_tbl) == 6)
  stopifnot("Design" %in% names(rank_tbl))
  cat("PASS\n")

  # Test 5: plot_mse_by_neighbor returns ggplot
  cat("Test 5: Existing plot functions return ggplot... ")
  p1 <- plot_mse_by_neighbor(test_results, names(configs)[1])
  stopifnot(inherits(p1, "gg"))
  p2 <- plot_coverage_mse_tradeoff(test_results, names(configs)[1])
  stopifnot(inherits(p2, "gg"))
  cat("PASS\n")

  # Test 6: 01-03 functions still available
  cat("Test 6: Module 01-03 functions exist... ")
  stopifnot(exists("build_spatial_grid", mode = "function"))
  stopifnot(exists("generate_incidence", mode = "function"))
  stopifnot(exists("get_designs", mode = "function"))
  stopifnot(exists("get_design_names", mode = "function"))
  cat("PASS\n")

  cat("\n=== All 6 side-effect tests PASSED ===\n")
  invisible(TRUE)
}
