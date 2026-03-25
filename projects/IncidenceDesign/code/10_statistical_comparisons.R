# ==============================================================================
# 10_statistical_comparisons.R
# Formal Statistical Hypothesis Testing for Treatment Assignment Design
# Comparison in Spatial Cluster Randomized Trials
#
# PURPOSE:
#   This module provides rigorous statistical tests to compare the performance
#   of 8 treatment assignment designs evaluated in the IncidenceDesign simulation
#   study. Rather than relying solely on descriptive comparisons of mean MSE,
#   this module applies non-parametric hypothesis tests to determine whether
#   observed differences in estimation accuracy are statistically significant.
#
# STATISTICAL FRAMEWORK:
#   The simulation produces MSE values for 8 designs across 320 unique parameter
#   combinations ("blocks"). Each block is defined by a unique combination of:
#     - Incidence_Mode (iid, spatial, poisson)
#     - Rho_Incidence (spatial correlation of incidence: 0, 0.20, 0.50)
#     - Neighbor_Type (rook, queen)
#     - Rho (outcome spatial autocorrelation: 0.00, 0.01, 0.20, 0.50)
#     - Gamma (spillover magnitude: 0.50, 0.60, 0.70, 0.80)
#     - Spillover_Type (control_only, both)
#
#   Within each block, all 8 designs share the same data-generating conditions,
#   creating a natural repeated-measures / paired structure. This structure is
#   exploited by the Friedman test (non-parametric repeated-measures ANOVA) and
#   pairwise Wilcoxon signed-rank tests.
#
# TESTS IMPLEMENTED:
#   1. Friedman Test (Omnibus)
#      - H0: All 8 designs have equal average MSE ranks across blocks
#      - H1: At least one design's average rank differs
#      - Non-parametric alternative to repeated-measures ANOVA
#      - Reference: Friedman (1937), JASA 32(198):675-701
#
#   2. Nemenyi Post-hoc Test
#      - Pairwise comparisons following a significant Friedman test
#      - Controls family-wise error rate across all C(8,2) = 28 comparisons
#      - Critical Difference: CD = q_alpha * sqrt(k*(k+1) / (6*N))
#        where k = 8 designs, N = number of blocks
#      - Two designs are significantly different if their average rank
#        difference exceeds CD
#      - Reference: Nemenyi (1963), Distribution-free Multiple Comparisons
#      - Application reference: Demšar (2006), JMLR 7:1-30
#
#   3. Wilcoxon Signed-Rank Test (Pairwise)
#      - H0: For designs i and j, paired MSE values have the same distribution
#      - H1: The distributions differ (two-sided)
#      - 28 pairwise tests with Holm correction for multiple comparisons
#      - Serves as sensitivity check alongside Friedman/Nemenyi
#      - Reference: Wilcoxon (1945), Biometrics Bulletin 1(6):80-83
#
# CONDITIONAL ANALYSIS:
#   All tests can be stratified by individual simulation parameters to assess
#   whether design rankings are robust or parameter-dependent. For example,
#   stratifying by Rho answers: "Do design rankings change with the level
#   of spatial autocorrelation in outcomes?"
#
# DEPENDENCIES:
#   - PMCMRplus (>= 1.9.0): Friedman test + Nemenyi post-hoc
#   - stats: wilcox.test(), p.adjust(), qtukey()
#   - ggplot2, dplyr, tidyr, viridis: data manipulation and visualization
#   - 06_visualizations.R: load_latest_results(), split_by_incidence_config()
#
# USAGE:
#   source("10_statistical_comparisons.R")
#   results <- load_latest_results(estimation_mode = "MLE_combined")
#   report  <- generate_comparison_report(results)
#
# OUTPUTS:
#   - Console: test statistics, p-values, average ranks, significance summaries
#   - PDF: results/MLE_statistical_comparisons.pdf containing:
#       * Critical Difference diagram (aggregate)
#       * MSE boxplot with significance stars
#       * Pairwise p-value heatmaps (Wilcoxon and Nemenyi)
#       * Conditional CD diagrams (per parameter level)
#
# AUTHOR: Andrew Walther
# DATE: 2026-03-25
# ==============================================================================

library(PMCMRplus)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# --- Source visualization module for data loading utilities ---
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
source(file.path(script_dir, "06_visualizations.R"))

# ==============================================================================
# HELPERS
# ==============================================================================

#' Build a comparison matrix from simulation results.
#'
#' Reshapes the results data frame into a wide matrix suitable for Friedman
#' and Wilcoxon tests: rows correspond to "blocks" (unique parameter
#' combinations) and columns correspond to the 8 designs. Each cell contains
#' the MSE (or other metric) for that design under that parameter combination.
#'
#' This is the core data structure underlying all statistical tests in this
#' module. The blocking structure ensures that comparisons are paired: within
#' each block, all 8 designs were evaluated under identical DGP conditions.
#'
#' @param results Data frame with columns: Design, MSE, and block identifier
#'   columns (Incidence_Mode, Rho_Incidence, Neighbor_Type, Rho, Gamma,
#'   Spillover_Type). Must contain exactly one row per design per block.
#' @param metric Character. Column name for the comparison metric.
#'   Default "MSE". Could also be "Coverage", "Bias", etc.
#' @return Named list with:
#'   \describe{
#'     \item{wide_matrix}{Numeric matrix (N_blocks x 8). Rows = blocks,
#'       columns = designs. Column names are design names.}
#'     \item{block_ids}{Data frame of block identifier columns, one row per
#'       block, aligned with rows of wide_matrix.}
#'     \item{design_names}{Character vector of design names in column order.}
#'   }
build_comparison_matrix <- function(results, metric = "MSE") {
  # Block variables: all simulation parameters EXCEPT Design and the metric
  block_vars <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type",
                   "Rho", "Gamma", "Spillover_Type")
  # Gracefully handle subsets that may be missing some block vars
  block_vars <- intersect(block_vars, names(results))

  designs <- sort(unique(results$Design))

  wide <- results %>%
    select(all_of(c(block_vars, "Design", metric))) %>%
    pivot_wider(names_from = Design, values_from = all_of(metric),
                id_cols = all_of(block_vars))

  block_ids <- wide %>% select(all_of(block_vars))
  wide_matrix <- as.matrix(wide[, designs])
  rownames(wide_matrix) <- NULL

  list(
    wide_matrix  = wide_matrix,
    block_ids    = block_ids,
    design_names = designs
  )
}

#' Create short design labels for compact plot annotations.
#'
#' Converts "Design 1" to "D1", "Design 8" to "D8", etc.
#'
#' @param design_name Character vector of full design names.
#' @return Character vector of abbreviated labels.
short_design_label <- function(design_name) {
  gsub("Design ", "D", design_name)
}

#' Convert a numeric p-value to significance star notation.
#'
#' Uses conventional thresholds:
#'   *** : p < 0.001
#'   **  : p < 0.01
#'   *   : p < 0.05
#'   ns  : p >= 0.05
#'
#' @param p Numeric p-value(s).
#' @return Character string(s) of significance stars.
p_to_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")))
}

# ==============================================================================
# CORE TEST FUNCTIONS
# ==============================================================================

#' Run the Friedman omnibus test for design differences.
#'
#' The Friedman test is a non-parametric alternative to one-way repeated-
#' measures ANOVA. For each block (parameter combination), the 8 designs
#' are ranked by MSE (rank 1 = lowest MSE = best). The test then evaluates
#' whether the average ranks across blocks differ significantly.
#'
#' Hypotheses:
#'   H0: All 8 designs have equal average MSE ranks (no design is
#'       systematically better or worse)
#'   H1: At least one design's average rank differs from the others
#'
#' The test statistic follows a chi-squared distribution with k-1 = 7 df
#' under H0.
#'
#' @param results Data frame of simulation results.
#' @param metric Character. Column to compare (default "MSE").
#' @return Named list with:
#'   \describe{
#'     \item{statistic}{Friedman chi-squared test statistic.}
#'     \item{p_value}{P-value from chi-squared approximation.}
#'     \item{avg_ranks}{Named numeric vector of average ranks per design,
#'       sorted ascending (rank 1 = best). Lower rank = lower MSE.}
#'     \item{n_blocks}{Number of blocks (unique parameter combinations) used.}
#'     \item{n_designs}{Number of designs compared (should be 8).}
#'     \item{test_object}{Full test result object from PMCMRplus::friedmanTest().}
#'   }
run_friedman_test <- function(results, metric = "MSE") {
  comp <- build_comparison_matrix(results, metric)
  mat <- comp$wide_matrix

  # stats::friedman.test operates on a matrix:
  # rows = blocks (subjects), cols = groups (treatments/designs)
  ft <- friedman.test(mat)

  # Compute average ranks: within each block, rank designs 1 (best) to k (worst)
  rank_matrix <- t(apply(mat, 1, rank, ties.method = "average"))
  avg_ranks <- colMeans(rank_matrix, na.rm = TRUE)
  names(avg_ranks) <- comp$design_names

  list(
    statistic   = as.numeric(ft$statistic),
    p_value     = as.numeric(ft$p.value),
    avg_ranks   = sort(avg_ranks),
    n_blocks    = nrow(mat),
    n_designs   = ncol(mat),
    test_object = ft
  )
}

#' Run Nemenyi post-hoc pairwise comparisons after Friedman test.
#'
#' The Nemenyi test compares all C(k,2) = 28 pairs of designs simultaneously,
#' controlling the family-wise error rate. Two designs are declared
#' significantly different if their average rank difference exceeds the
#' Critical Difference (CD):
#'
#'   CD = q_{alpha,k,Inf} / sqrt(2) * sqrt( k*(k+1) / (6*N) )
#'
#' where:
#'   q_{alpha,k,Inf} = critical value of the Studentized range distribution
#'   k = 8 (number of designs)
#'   N = number of blocks
#'
#' The results are visualized as a Critical Difference (CD) diagram: designs
#' are placed on a number line by their average rank, and horizontal bars
#' connect groups of designs that are NOT significantly different.
#'
#' @param results Data frame of simulation results.
#' @param metric Character. Column to compare (default "MSE").
#' @param alpha Numeric. Significance level (default 0.05).
#' @return Named list with:
#'   \describe{
#'     \item{p_matrix}{8x8 symmetric matrix of adjusted pairwise p-values.
#'       p_matrix[i,j] = p-value for testing Design i vs Design j.
#'       Diagonal is NA.}
#'     \item{sig_matrix}{8x8 logical matrix. TRUE where p < alpha.}
#'     \item{avg_ranks}{Named numeric vector of average ranks (sorted ascending).}
#'     \item{critical_diff}{Nemenyi critical difference value.}
#'     \item{alpha}{Significance level used.}
#'     \item{n_blocks}{Number of blocks used.}
#'     \item{test_object}{Full result from PMCMRplus::frdAllPairsNemenyiTest().}
#'   }
run_nemenyi_posthoc <- function(results, metric = "MSE", alpha = 0.05) {
  comp <- build_comparison_matrix(results, metric)
  mat <- comp$wide_matrix
  designs <- comp$design_names
  n_designs <- length(designs)

  # Nemenyi post-hoc test via PMCMRplus.
  # Use the y/groups/blocks interface (long format) because the matrix
  # interface triggers an internal error in frdRanks for some column names.
  block_vars <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type",
                   "Rho", "Gamma", "Spillover_Type")
  block_vars <- intersect(block_vars, names(results))
  results$block_id <- interaction(results[, block_vars])

  nt <- frdAllPairsNemenyiTest(
    y      = results[[metric]],
    groups = factor(results$Design),
    blocks = factor(results$block_id)
  )

  # PMCMRplus returns a lower-triangular (n-1) x (n-1) p-value matrix.
  # We reconstruct a full symmetric n x n matrix for convenient access.
  p_raw <- nt$p.value
  # Row/col names from PMCMRplus output (may differ from our design ordering)
  raw_designs <- c(colnames(p_raw), rownames(p_raw)[n_designs - 1])
  # Build mapping from PMCMRplus names back to our sorted design names
  p_matrix <- matrix(NA, n_designs, n_designs,
                     dimnames = list(designs, designs))
  raw_rows <- rownames(p_raw)
  raw_cols <- colnames(p_raw)
  for (ri in seq_along(raw_rows)) {
    for (ci in seq_along(raw_cols)) {
      p_val <- p_raw[ri, ci]
      if (!is.na(p_val)) {
        p_matrix[raw_rows[ri], raw_cols[ci]] <- p_val
        p_matrix[raw_cols[ci], raw_rows[ri]] <- p_val
      }
    }
  }

  # Average ranks (same computation as in run_friedman_test)
  rank_matrix <- t(apply(mat, 1, rank, ties.method = "average"))
  avg_ranks <- colMeans(rank_matrix, na.rm = TRUE)
  names(avg_ranks) <- designs

  # Nemenyi critical difference formula:
  #   CD = q_{alpha,k,Inf} / sqrt(2)  *  sqrt( k*(k+1) / (6*N) )
  k <- n_designs
  n <- nrow(mat)
  q_alpha <- qtukey(1 - alpha, k, Inf) / sqrt(2)
  cd <- q_alpha * sqrt(k * (k + 1) / (6 * n))

  sig_matrix <- p_matrix < alpha
  sig_matrix[is.na(sig_matrix)] <- FALSE

  list(
    p_matrix      = p_matrix,
    sig_matrix    = sig_matrix,
    avg_ranks     = sort(avg_ranks),
    critical_diff = cd,
    alpha         = alpha,
    n_blocks      = n,
    test_object   = nt
  )
}

#' Run pairwise Wilcoxon signed-rank tests with multiple testing correction.
#'
#' For each of the C(8,2) = 28 pairs of designs, performs a two-sided
#' Wilcoxon signed-rank test on the paired MSE values (one pair per block).
#'
#' Hypotheses (for each pair i, j):
#'   H0: The distribution of MSE_i - MSE_j is symmetric about zero
#'       (i.e., designs i and j have equivalent MSE performance)
#'   H1: The distribution is not symmetric about zero
#'       (i.e., one design systematically outperforms the other)
#'
#' P-values are adjusted for the 28 simultaneous tests using the Holm
#' step-down procedure (Holm, 1979), which controls the family-wise
#' error rate while being less conservative than Bonferroni.
#'
#' @param results Data frame of simulation results.
#' @param metric Character. Column to compare (default "MSE").
#' @param p_adjust Character. P-value adjustment method passed to
#'   stats::p.adjust(). Default "holm". Other options: "bonferroni",
#'   "BH" (Benjamini-Hochberg FDR), "none".
#' @return Named list with:
#'   \describe{
#'     \item{p_matrix}{8x8 symmetric matrix of adjusted p-values.}
#'     \item{p_matrix_raw}{8x8 symmetric matrix of unadjusted p-values.}
#'     \item{sig_matrix}{8x8 logical matrix. TRUE where adjusted p < 0.05.}
#'     \item{n_tests}{Number of pairwise tests (28 for 8 designs).}
#'     \item{p_adjust_method}{Correction method used.}
#'   }
run_pairwise_wilcoxon <- function(results, metric = "MSE", p_adjust = "holm") {
  comp <- build_comparison_matrix(results, metric)
  mat <- comp$wide_matrix
  designs <- comp$design_names
  n_designs <- length(designs)
  n_tests <- choose(n_designs, 2)

  # Collect raw p-values from all C(8,2) = 28 pairwise tests
  p_raw_vec <- numeric(n_tests)
  pair_indices <- matrix(NA, n_tests, 2)
  idx <- 1
  for (i in 1:(n_designs - 1)) {
    for (j in (i + 1):n_designs) {
      # Paired test: same block (row), different designs (columns)
      # exact = FALSE to use normal approximation (avoids ties warning)
      wt <- wilcox.test(mat[, i], mat[, j], paired = TRUE, exact = FALSE)
      p_raw_vec[idx] <- wt$p.value
      pair_indices[idx, ] <- c(i, j)
      idx <- idx + 1
    }
  }

  # Adjust all 28 p-values simultaneously
  p_adj_vec <- p.adjust(p_raw_vec, method = p_adjust)

  # Build full symmetric matrices for convenient lookup
  p_matrix_raw <- p_matrix <- matrix(NA, n_designs, n_designs,
                                      dimnames = list(designs, designs))
  for (idx in seq_len(n_tests)) {
    i <- pair_indices[idx, 1]
    j <- pair_indices[idx, 2]
    p_matrix_raw[i, j] <- p_raw_vec[idx]
    p_matrix_raw[j, i] <- p_raw_vec[idx]
    p_matrix[i, j]     <- p_adj_vec[idx]
    p_matrix[j, i]     <- p_adj_vec[idx]
  }

  sig_matrix <- p_matrix < 0.05
  sig_matrix[is.na(sig_matrix)] <- FALSE

  list(
    p_matrix        = p_matrix,
    p_matrix_raw    = p_matrix_raw,
    sig_matrix      = sig_matrix,
    n_tests         = n_tests,
    p_adjust_method = p_adjust
  )
}

# ==============================================================================
# CONDITIONAL / STRATIFIED TESTS
# ==============================================================================

#' Run hypothesis tests stratified by a single simulation parameter.
#'
#' Subsets the results by each level of stratify_by, then runs the Friedman
#' (+ Nemenyi if significant) and/or Wilcoxon tests within each subset. This
#' answers questions like: "Do design rankings change when we condition on
#' a specific level of spatial autocorrelation (Rho)?"
#'
#' The number of blocks per stratum depends on the parameter:
#'   - Rho (4 levels):          80 blocks per level
#'   - Gamma (4 levels):        80 blocks per level
#'   - Spillover_Type (2):      160 blocks per level
#'   - Neighbor_Type (2):       160 blocks per level
#'   - Incidence_Mode (3):      variable (64-128 blocks per mode)
#'
#' @param results Data frame of simulation results.
#' @param stratify_by Character. Column name to stratify by.
#' @param metric Character. Column to compare (default "MSE").
#' @param test_type Character. Which tests to run: "friedman" (Friedman +
#'   Nemenyi only), "wilcoxon" (Wilcoxon only), or "both" (default).
#' @return Named list with one entry per stratum level. Each entry contains:
#'   \describe{
#'     \item{level}{The parameter value for this stratum.}
#'     \item{n_blocks}{Number of blocks in this stratum.}
#'     \item{friedman}{Friedman test result (if test_type includes it).}
#'     \item{nemenyi}{Nemenyi post-hoc result (if Friedman was significant).}
#'     \item{wilcoxon}{Wilcoxon result (if test_type includes it).}
#'   }
run_conditional_tests <- function(results, stratify_by, metric = "MSE",
                                  test_type = "both") {
  if (!stratify_by %in% names(results)) {
    stop("Column '", stratify_by, "' not found in results")
  }

  levels <- sort(unique(results[[stratify_by]]))
  out <- list()

  for (lev in levels) {
    subset_data <- results[results[[stratify_by]] == lev, ]
    label <- paste0(stratify_by, " = ", lev)

    entry <- list(level = lev, n_blocks = NULL)

    if (test_type %in% c("friedman", "both")) {
      fr <- tryCatch(run_friedman_test(subset_data, metric), error = function(e) NULL)
      if (!is.null(fr)) {
        entry$friedman <- fr
        entry$n_blocks <- fr$n_blocks
        # Run Nemenyi post-hoc only if Friedman rejects H0
        if (fr$p_value < 0.05) {
          entry$nemenyi <- tryCatch(
            run_nemenyi_posthoc(subset_data, metric),
            error = function(e) NULL
          )
        }
      }
    }

    if (test_type %in% c("wilcoxon", "both")) {
      entry$wilcoxon <- tryCatch(
        run_pairwise_wilcoxon(subset_data, metric),
        error = function(e) NULL
      )
    }

    out[[label]] <- entry
  }

  out
}

# ==============================================================================
# VISUALIZATION: Critical Difference (CD) Diagram
# ==============================================================================

#' Plot a Critical Difference diagram.
#'
#' The CD diagram (Demšar, 2006) is the standard visualization for comparing
#' multiple methods across multiple datasets/scenarios. Designs are placed on
#' a horizontal number line by their average rank (leftward = better). Thick
#' horizontal bars ("cliques") connect groups of designs whose average rank
#' differences do NOT exceed the Nemenyi critical difference — i.e., designs
#' within a clique are statistically indistinguishable.
#'
#' A red arrow at the bottom indicates the CD scale: any two designs whose
#' rank positions are farther apart than this arrow are significantly different.
#'
#' @param nemenyi_result List. Output from run_nemenyi_posthoc().
#' @param title Character or NULL. Plot title. If NULL, auto-generated with
#'   alpha and CD values.
#' @return A ggplot object.
plot_cd_diagram <- function(nemenyi_result, title = NULL) {
  avg_ranks <- nemenyi_result$avg_ranks
  cd <- nemenyi_result$critical_diff
  designs <- names(avg_ranks)
  n <- length(designs)
  alpha <- nemenyi_result$alpha

  if (is.null(title)) {
    title <- sprintf("Critical Difference Diagram (alpha = %.2f, CD = %.3f)",
                     alpha, cd)
  }

  # Identify equivalence groups: maximal sets of designs where all pairwise
  # p-values are >= alpha (i.e., no significant differences within the group)
  groups <- find_cd_groups(avg_ranks, nemenyi_result$p_matrix, alpha)

  # Data frame for design points on the rank axis
  df_points <- data.frame(
    design    = factor(designs, levels = designs),
    avg_rank  = avg_ranks,
    label     = short_design_label(designs),
    stringsAsFactors = FALSE
  )
  # Alternate labels above/below axis for readability
  df_points$y_pos <- ifelse(seq_len(n) %% 2 == 1, 1.0, -1.0)

  # Build horizontal bars for equivalence groups
  df_bars <- data.frame(xmin = numeric(0), xmax = numeric(0), y = numeric(0))
  if (length(groups) > 0) {
    bar_y_start <- 1.8
    bar_y_step <- 0.4
    for (g in seq_along(groups)) {
      grp_ranks <- avg_ranks[groups[[g]]]
      df_bars <- rbind(df_bars, data.frame(
        xmin = min(grp_ranks) - 0.05,
        xmax = max(grp_ranks) + 0.05,
        y    = bar_y_start + (g - 1) * bar_y_step
      ))
    }
  }

  p <- ggplot() +
    # Rank axis line
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    # Tick marks from axis to label positions
    geom_segment(data = df_points,
                 aes(x = avg_rank, xend = avg_rank, y = 0, yend = y_pos * 0.3),
                 color = "grey40", linewidth = 0.4) +
    # Design points
    geom_point(data = df_points, aes(x = avg_rank, y = y_pos * 0.5),
               size = 3, color = "black") +
    # Design name labels
    geom_text(data = df_points,
              aes(x = avg_rank, y = y_pos * 0.8, label = label),
              size = 3.5, fontface = "bold") +
    # Numeric rank labels
    geom_text(data = df_points,
              aes(x = avg_rank, y = y_pos * 0.35,
                  label = sprintf("%.2f", avg_rank)),
              size = 2.5, color = "grey30") +
    # Equivalence group bars (thick segments connecting non-significant groups)
    {if (nrow(df_bars) > 0)
      geom_segment(data = df_bars,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   linewidth = 2.5, color = "grey30", lineend = "round")
    } +
    # CD scale indicator (red double-arrow at bottom)
    annotate("segment", x = 1, xend = 1 + cd, y = -2.0, yend = -2.0,
             linewidth = 1.0, color = "red3",
             arrow = arrow(ends = "both", length = unit(0.08, "inches"))) +
    annotate("text", x = 1 + cd / 2, y = -2.4,
             label = sprintf("CD = %.3f", cd),
             size = 3, color = "red3") +
    scale_x_continuous(
      breaks = 1:n,
      limits = c(0.5, n + 0.5),
      name = "Average Rank (lower = better MSE)"
    ) +
    coord_cartesian(ylim = c(-3, max(3, max(df_bars$y, 0) + 1))) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 13)
    ) +
    ggtitle(title)

  p
}

#' Find equivalence groups for the CD diagram.
#'
#' An equivalence group is a maximal set of designs where every pair within
#' the group has a non-significant p-value (p >= alpha) in the Nemenyi test.
#' These groups are drawn as horizontal bars on the CD diagram.
#'
#' Algorithm: greedy left-to-right sweep over designs sorted by average rank.
#' For each starting design, extend the group rightward as long as the new
#' design is non-significant with ALL current members. Only groups with 2+
#' members are retained. Subsets of existing groups are pruned.
#'
#' @param avg_ranks Named numeric vector of average ranks (sorted ascending).
#' @param p_matrix Symmetric p-value matrix from Nemenyi test.
#' @param alpha Significance threshold.
#' @return List of character vectors, each a group of design names.
find_cd_groups <- function(avg_ranks, p_matrix, alpha = 0.05) {
  designs <- names(avg_ranks)  # already sorted by rank
  n <- length(designs)
  groups <- list()

  for (i in seq_len(n)) {
    group <- designs[i]
    for (j in (i + 1):n) {
      if (j > n) break
      # Check: is design j non-significant with ALL current group members?
      all_nonsig <- TRUE
      for (member in group) {
        p_val <- p_matrix[member, designs[j]]
        if (!is.na(p_val) && p_val < alpha) {
          all_nonsig <- FALSE
          break
        }
      }
      if (all_nonsig) {
        group <- c(group, designs[j])
      } else {
        break
      }
    }
    # Only keep groups of 2+ that aren't subsets of existing groups
    if (length(group) >= 2) {
      is_subset <- FALSE
      for (existing in groups) {
        if (all(group %in% existing)) {
          is_subset <- TRUE
          break
        }
      }
      if (!is_subset) groups[[length(groups) + 1]] <- group
    }
  }

  groups
}

# ==============================================================================
# VISUALIZATION: MSE Boxplot with Significance Stars
# ==============================================================================

#' MSE boxplot by design with pairwise significance annotations.
#'
#' Displays the distribution of scenario-level MSE for each design as a
#' boxplot, ordered by median MSE. When Wilcoxon test results are provided,
#' significant pairwise comparisons are annotated with bracket lines and
#' significance stars above the plot.
#'
#' To avoid visual clutter, a maximum of 8 annotations are shown, prioritizing
#' comparisons involving the best-ranked and worst-ranked designs.
#'
#' @param results Data frame of simulation results.
#' @param wilcoxon_result List. Output from run_pairwise_wilcoxon(). If NULL,
#'   no significance annotations are drawn.
#' @param ref_design Character or NULL. If provided, only comparisons against
#'   this reference design are shown (useful for "Design 8 vs all" comparisons).
#' @param title Character or NULL. Plot title.
#' @return A ggplot object.
plot_mse_boxplot_with_stars <- function(results, wilcoxon_result = NULL,
                                        ref_design = NULL, title = NULL) {
  if (is.null(title)) title <- "MSE by Design with Pairwise Significance"

  # Order designs by median MSE (best on left)
  design_order <- results %>%
    group_by(Design) %>%
    summarise(med = median(MSE, na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(Design)

  results$Design <- factor(results$Design, levels = design_order)

  p <- ggplot(results, aes(x = Design, y = MSE, fill = Design)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.3, outlier.size = 1) +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = title, x = NULL, y = "MSE")

  # Add significance brackets and stars
  if (!is.null(wilcoxon_result)) {
    annot <- build_significance_annotations(
      wilcoxon_result, design_order, results, ref_design
    )
    if (nrow(annot) > 0) {
      max_mse <- max(results$MSE, na.rm = TRUE)
      for (r in seq_len(nrow(annot))) {
        y_pos <- max_mse * (1.05 + 0.06 * (r - 1))
        p <- p +
          annotate("segment",
                   x = annot$x1[r], xend = annot$x2[r],
                   y = y_pos, yend = y_pos,
                   linewidth = 0.4) +
          annotate("text",
                   x = (annot$x1[r] + annot$x2[r]) / 2,
                   y = y_pos + max_mse * 0.015,
                   label = annot$stars[r],
                   size = 3.5)
      }
    }
  }

  p
}

#' Build significance annotation coordinates for boxplot brackets.
#'
#' Selects which pairwise comparisons to display (limited to avoid clutter)
#' and computes bracket positions on the x-axis.
#'
#' @param wilcoxon_result Wilcoxon test output.
#' @param design_order Character vector of designs in display order.
#' @param results Data frame (for computing y-axis range).
#' @param ref_design If non-NULL, restrict to comparisons involving this design.
#' @return Data frame with columns: x1, x2 (integer positions), stars, p_val.
build_significance_annotations <- function(wilcoxon_result, design_order,
                                            results, ref_design = NULL) {
  p_mat <- wilcoxon_result$p_matrix
  designs <- design_order

  annots <- data.frame(x1 = integer(), x2 = integer(),
                       stars = character(), p_val = numeric(),
                       stringsAsFactors = FALSE)

  if (!is.null(ref_design)) {
    # Show all significant comparisons vs reference design
    ref_idx <- which(designs == ref_design)
    for (j in seq_along(designs)) {
      if (j == ref_idx) next
      p_val <- p_mat[ref_design, designs[j]]
      if (!is.na(p_val) && p_val < 0.05) {
        annots <- rbind(annots, data.frame(
          x1 = min(ref_idx, j), x2 = max(ref_idx, j),
          stars = p_to_stars(p_val), p_val = p_val
        ))
      }
    }
  } else {
    # Collect all significant pairs
    pairs <- data.frame(d1 = character(), d2 = character(),
                        p_val = numeric(), stringsAsFactors = FALSE)
    for (i in 1:(length(designs) - 1)) {
      for (j in (i + 1):length(designs)) {
        p_val <- p_mat[designs[i], designs[j]]
        if (!is.na(p_val) && p_val < 0.05) {
          pairs <- rbind(pairs, data.frame(
            d1 = designs[i], d2 = designs[j], p_val = p_val
          ))
        }
      }
    }
    # Prioritize comparisons involving best or worst design, limit to 8
    if (nrow(pairs) > 0) {
      pairs$involves_extreme <- pairs$d1 %in% designs[c(1, length(designs))] |
        pairs$d2 %in% designs[c(1, length(designs))]
      pairs <- pairs %>%
        arrange(desc(involves_extreme), p_val) %>%
        head(8)

      for (r in seq_len(nrow(pairs))) {
        x1 <- which(designs == pairs$d1[r])
        x2 <- which(designs == pairs$d2[r])
        annots <- rbind(annots, data.frame(
          x1 = min(x1, x2), x2 = max(x1, x2),
          stars = p_to_stars(pairs$p_val[r]), p_val = pairs$p_val[r]
        ))
      }
    }
  }

  annots
}

# ==============================================================================
# VISUALIZATION: P-value Heatmap
# ==============================================================================

#' Plot an 8x8 heatmap of pairwise p-values.
#'
#' Each cell shows the adjusted p-value for the corresponding design pair,
#' with significance stars overlaid. Color scale runs from red (highly
#' significant, p ~ 0) through yellow (p ~ 0.05) to green (non-significant,
#' p ~ 1). Diagonal cells are grey (self-comparison).
#'
#' @param pairwise_result List. Output from run_pairwise_wilcoxon() or
#'   run_nemenyi_posthoc() — any object with a $p_matrix element.
#' @param test_label Character. Label for the test type shown in title.
#' @param title Character or NULL. Plot title.
#' @return A ggplot object.
plot_pvalue_heatmap <- function(pairwise_result, test_label = "Wilcoxon",
                                 title = NULL) {
  p_mat <- pairwise_result$p_matrix
  designs <- rownames(p_mat)

  if (is.null(title)) {
    adj_label <- if (!is.null(pairwise_result$p_adjust_method)) {
      sprintf(" (%s-adjusted)", tools::toTitleCase(pairwise_result$p_adjust_method))
    } else ""
    title <- sprintf("Pairwise %s p-values%s", test_label, adj_label)
  }

  # Convert to long format for ggplot
  df <- expand.grid(Design1 = designs, Design2 = designs,
                    stringsAsFactors = FALSE)
  df$p_value <- mapply(function(d1, d2) {
    if (d1 == d2) return(NA)
    p_mat[d1, d2]
  }, df$Design1, df$Design2)

  df$sig_label <- ifelse(is.na(df$p_value), "",
                         p_to_stars(df$p_value))

  df$Design1 <- factor(df$Design1, levels = designs)
  df$Design2 <- factor(df$Design2, levels = rev(designs))

  ggplot(df, aes(x = Design1, y = Design2, fill = p_value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 4, fontface = "bold") +
    scale_fill_gradientn(
      colours = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850"),
      values  = c(0, 0.001, 0.01, 0.05, 0.10, 1.0),
      limits  = c(0, 1),
      na.value = "grey90",
      name = "p-value",
      breaks = c(0.001, 0.01, 0.05, 0.10, 0.50, 1.0),
      labels = c("0.001", "0.01", "0.05", "0.10", "0.50", "1.00")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    ) +
    labs(title = title, x = NULL, y = NULL)
}

# ==============================================================================
# VISUALIZATION: Conditional CD Diagrams
# ==============================================================================

#' Generate a list of CD diagrams, one per level of a stratification parameter.
#'
#' If the Friedman test is non-significant for a stratum, a placeholder plot
#' with a text note is returned instead of a CD diagram.
#'
#' @param conditional_results List. Output from run_conditional_tests().
#' @param stratify_by Character. Parameter name (used in plot titles).
#' @return Named list of ggplot objects.
plot_conditional_cd_diagrams <- function(conditional_results, stratify_by) {
  plots <- list()
  for (name in names(conditional_results)) {
    entry <- conditional_results[[name]]
    if (!is.null(entry$nemenyi)) {
      plots[[name]] <- plot_cd_diagram(entry$nemenyi, title = name)
    } else if (!is.null(entry$friedman)) {
      # Friedman was non-significant — show a note
      plots[[name]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = sprintf("%s\nFriedman p = %.3f (not significant)",
                                 name, entry$friedman$p_value),
                 size = 4) +
        theme_void() +
        ggtitle(name)
    }
  }
  plots
}

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

#' Summarize conditional Friedman test results into a single data frame.
#'
#' For each stratum, reports the Friedman chi-squared, p-value, number of
#' blocks, and which designs are in the top equivalence group (i.e., the
#' group containing the best-ranked design where all members are statistically
#' indistinguishable from the best).
#'
#' @param conditional_results List. Output from run_conditional_tests().
#' @return Data frame with columns: Stratum, Chi_squared, P_value,
#'   N_blocks, Best_Group.
summarize_conditional_tests <- function(conditional_results) {
  rows <- lapply(names(conditional_results), function(name) {
    entry <- conditional_results[[name]]
    fr <- entry$friedman

    best_group <- if (!is.null(entry$nemenyi)) {
      # Find the equivalence group containing the best-ranked design
      groups <- find_cd_groups(entry$nemenyi$avg_ranks,
                               entry$nemenyi$p_matrix,
                               entry$nemenyi$alpha)
      best_design <- names(entry$nemenyi$avg_ranks)[1]
      found <- NULL
      for (g in groups) {
        if (best_design %in% g) {
          found <- paste(short_design_label(g), collapse = ", ")
          break
        }
      }
      if (is.null(found)) short_design_label(best_design)
      else found
    } else {
      "All equivalent"
    }

    data.frame(
      Stratum      = name,
      Chi_squared  = if (!is.null(fr)) round(fr$statistic, 2) else NA,
      P_value      = if (!is.null(fr)) fr$p_value else NA,
      N_blocks     = if (!is.null(fr)) fr$n_blocks else NA,
      Best_Group   = best_group,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# ==============================================================================
# REPORT GENERATION
# ==============================================================================

#' Run all statistical comparisons and generate outputs.
#'
#' This is the top-level function that orchestrates the full analysis:
#'   1. Aggregate Friedman test (omnibus)
#'   2. Nemenyi post-hoc (if Friedman significant)
#'   3. Pairwise Wilcoxon signed-rank tests (Holm-corrected)
#'   4. Conditional tests stratified by each simulation parameter
#'   5. All visualizations (CD diagram, boxplot, heatmaps, conditional CDs)
#'   6. PDF output to results directory
#'
#' @param results Data frame of simulation results (combined MLE).
#' @param output_dir Character. Directory for PDF output.
#'   Default: results/ directory relative to code/.
#' @param alpha Numeric. Significance level (default 0.05).
#' @param output_pdf Logical. If TRUE (default), save plots to PDF.
#' @return Invisible list with all test results and plot objects:
#'   \describe{
#'     \item{friedman}{Aggregate Friedman test result.}
#'     \item{nemenyi}{Aggregate Nemenyi post-hoc result.}
#'     \item{wilcoxon}{Aggregate pairwise Wilcoxon result.}
#'     \item{conditional}{Named list of conditional test results per parameter.}
#'     \item{conditional_summaries}{Named list of summary tables per parameter.}
#'     \item{plots}{List of all ggplot objects.}
#'   }
generate_comparison_report <- function(results,
                                       output_dir = file.path(dirname(script_dir), "results"),
                                       alpha = 0.05,
                                       output_pdf = TRUE) {
  cat("\n", strrep("=", 70), "\n")
  cat("STATISTICAL DESIGN COMPARISONS\n")
  cat(strrep("=", 70), "\n\n")

  # -------------------------------------------------------------------------
  # 1. Aggregate Friedman Test
  # -------------------------------------------------------------------------
  cat("--- Aggregate Friedman Test ---\n")
  cat("  H0: All 8 designs have equal average MSE ranks\n")
  cat("  H1: At least one design differs\n\n")
  fr <- run_friedman_test(results)
  cat(sprintf("  Chi-squared = %.2f, df = %d, p-value = %s\n",
              fr$statistic, fr$n_designs - 1,
              format.pval(fr$p_value, digits = 3)))
  cat(sprintf("  N blocks = %d\n\n", fr$n_blocks))
  cat("  Average ranks (lower = better MSE):\n")
  for (d in names(fr$avg_ranks)) {
    cat(sprintf("    %s: %.3f\n", d, fr$avg_ranks[d]))
  }
  if (fr$p_value < alpha) {
    cat("\n  >> Friedman test SIGNIFICANT: proceeding to post-hoc comparisons\n")
  } else {
    cat("\n  >> Friedman test NOT significant: no evidence of design differences\n")
  }

  # -------------------------------------------------------------------------
  # 2. Nemenyi Post-hoc Test
  # -------------------------------------------------------------------------
  cat("\n--- Nemenyi Post-hoc Test ---\n")
  nem <- run_nemenyi_posthoc(results, alpha = alpha)
  cat(sprintf("  Critical Difference (CD) = %.3f at alpha = %.2f\n",
              nem$critical_diff, alpha))
  cat("  Designs with avg rank difference > CD are significantly different.\n")

  # Count significant pairs
  n_sig_nem <- sum(nem$sig_matrix[upper.tri(nem$sig_matrix)])
  cat(sprintf("  %d of %d pairs significantly different\n\n",
              n_sig_nem, choose(nem$n_blocks, 2)))

  # -------------------------------------------------------------------------
  # 3. Pairwise Wilcoxon Signed-Rank Tests
  # -------------------------------------------------------------------------
  cat("--- Pairwise Wilcoxon Signed-Rank Tests ---\n")
  cat(sprintf("  Correction method: %s\n", "Holm"))
  wil <- run_pairwise_wilcoxon(results)
  n_sig_wil <- sum(wil$sig_matrix[upper.tri(wil$sig_matrix)])
  cat(sprintf("  %d of %d pairs significantly different (adj. p < 0.05)\n",
              n_sig_wil, wil$n_tests))

  # -------------------------------------------------------------------------
  # 4. Conditional / Stratified Tests
  # -------------------------------------------------------------------------
  cat("\n--- Conditional Tests (single-parameter stratification) ---\n")
  cond_params <- c("Rho", "Gamma", "Spillover_Type", "Neighbor_Type",
                    "Incidence_Mode")
  cond_params <- intersect(cond_params, names(results))

  cond_results <- list()
  cond_summaries <- list()
  for (param in cond_params) {
    cat(sprintf("\n  Stratifying by %s:\n", param))
    cr <- run_conditional_tests(results, param, test_type = "both")
    cond_results[[param]] <- cr
    cond_summaries[[param]] <- summarize_conditional_tests(cr)
    for (name in names(cr)) {
      entry <- cr[[name]]
      if (!is.null(entry$friedman)) {
        sig_label <- if (entry$friedman$p_value < alpha) "***" else "ns"
        cat(sprintf("    %s: Friedman p = %s [%s] (n = %d blocks)\n",
                    name,
                    format.pval(entry$friedman$p_value, digits = 3),
                    sig_label,
                    entry$friedman$n_blocks))
      }
    }
  }

  # -------------------------------------------------------------------------
  # 5. Generate Visualizations
  # -------------------------------------------------------------------------
  cat("\n--- Generating visualizations ---\n")
  plots <- list()

  plots$cd_diagram <- plot_cd_diagram(nem)
  plots$boxplot_stars <- plot_mse_boxplot_with_stars(results, wil)
  plots$pvalue_heatmap_wilcoxon <- plot_pvalue_heatmap(wil, "Wilcoxon (Holm)")
  plots$pvalue_heatmap_nemenyi <- plot_pvalue_heatmap(nem, "Nemenyi")

  plots$conditional_cd <- list()
  for (param in names(cond_results)) {
    plots$conditional_cd[[param]] <- plot_conditional_cd_diagrams(
      cond_results[[param]], param
    )
  }

  # -------------------------------------------------------------------------
  # 6. Save to PDF
  # -------------------------------------------------------------------------
  if (output_pdf) {
    pdf_path <- file.path(output_dir, "MLE_statistical_comparisons.pdf")
    cat(sprintf("  Saving to: %s\n", pdf_path))

    grDevices::pdf(pdf_path, width = 10, height = 6)

    print(plots$cd_diagram)
    print(plots$boxplot_stars)
    print(plots$pvalue_heatmap_wilcoxon)
    print(plots$pvalue_heatmap_nemenyi)

    for (param in names(plots$conditional_cd)) {
      for (pl in plots$conditional_cd[[param]]) {
        print(pl)
      }
    }

    grDevices::dev.off()
    cat("  PDF saved successfully.\n")
  }

  cat("\n", strrep("=", 70), "\n")
  cat("STATISTICAL COMPARISONS COMPLETE\n")
  cat(strrep("=", 70), "\n")

  invisible(list(
    friedman            = fr,
    nemenyi             = nem,
    wilcoxon            = wil,
    conditional         = cond_results,
    conditional_summaries = cond_summaries,
    plots               = plots
  ))
}
