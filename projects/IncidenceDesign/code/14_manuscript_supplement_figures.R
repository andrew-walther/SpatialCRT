# ==============================================================================
# 14_manuscript_supplement_figures.R
# Data and figures for the CTJ manuscript's Supplementary Information document
# (paper/ctj_manuscript/Supplementary_Information.tex) and two main-text figure
# replacements.
#
# Three deliverables:
#   (a) A fresh 8-design consolidation summary (Friedman/Nemenyi/Wilcoxon on
#       ALL 8 designs, not the 6 retained for the manuscript) -- justifies
#       dropping Saturation Quadrants and Balanced Halves. This is the ONLY
#       8-design content anywhere in the manuscript + SI.
#   (b) Two main-text figure replacements (6 designs, full names, best->worst
#       ordering): reordered MSE bar chart, new ranked bias-variance figure.
#   (c) An SI figure suite (6 designs, full names, ranking-oriented): CD
#       diagrams, MSE/rank heatmaps, two-panel bar+p-value, tau-sweep lines.
#   (d) A 6-design application table from the NC application run.
#
# All 6-design descriptive numbers reuse the already-verified
# results/six_design_manuscript/six_design_comparison_report.rds (produced by
# 12_six_design_statistical_comparisons.R) rather than recomputing -- this
# guarantees the SI matches six_design_summary.txt exactly. Only the 8-design
# test results are computed fresh here.
#
# No "Design N" / "D1" shorthand anywhere: every Design column value below is
# the full manuscript name (get_design_names() with the "Design N: " prefix
# stripped), per the CTJ + dissertation manuscript no-shorthand convention.
#
# OUTPUTS:
#   results/eight_design_supplementary/eight_design_summary.txt
#   results/six_design_manuscript/fig_mse_by_design_6design.pdf   (overwritten, reordered)
#   results/six_design_manuscript/fig_biasvar_6design.pdf         (new)
#   results/six_design_manuscript/application_table_6design.txt   (new)
#   results/six_design_manuscript/si_figures/*.pdf                (new, 6 files)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
setwd(script_dir)
results_dir <- file.path(dirname(script_dir), "results")  # before sourcing 10, which overwrites script_dir
proj_dir    <- dirname(script_dir)

source("10_statistical_comparisons.R")  # sources 06_visualizations.R -> 01-03

six_dir  <- file.path(results_dir, "six_design_manuscript")
si_dir   <- file.path(six_dir, "si_figures")
eight_dir <- file.path(results_dir, "eight_design_supplementary")
dir.create(si_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(eight_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# FULL-NAME RELABELING (mirrors code/12_six_design_statistical_comparisons.R)
# ==============================================================================

raw_names <- get_design_names()  # e.g. "Design 1: Checkerboard"
full_name_map <- setNames(
  sub("^Design [0-9]+: ", "", raw_names),
  paste0("Design ", seq_along(raw_names))
)
stopifnot(identical(unname(full_name_map),
                     c("Checkerboard", "High Incidence Focus", "Saturation Quadrants",
                       "Isolation Buffer", "2x2 Blocking", "Balanced Quartiles",
                       "Balanced Halves", "Incidence-Guided Saturation Quadrants")))

retained_ids <- c("Design 1", "Design 2", "Design 4", "Design 5", "Design 6", "Design 8")

# Conceptual group membership, restricted to the 6 retained designs, keyed by
# full manuscript name (mirrors DESIGN_GROUPS in 00_design_names.R, which is
# keyed by short names and includes all 8).
DESIGN_GROUPS_6 <- c(
  "Checkerboard"                            = "Blocking",
  "2x2 Blocking"                            = "Blocking",
  "Isolation Buffer"                        = "Blocking",
  "High Incidence Focus"                    = "Stratified",
  "Balanced Quartiles"                      = "Stratified",
  "Incidence-Guided Saturation Quadrants"   = "Saturation"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

mle_full <- load_latest_results(results_dir = results_dir, estimation_mode = "MLE_tau_sweep")
stopifnot("True_Tau" %in% names(mle_full))

# --- 8-design (all designs, relabeled, no filter) ---
mle_full_8 <- mle_full
mle_full_8$Design <- full_name_map[mle_full_8$Design]
mle_tau1_8 <- mle_full_8[mle_full_8$True_Tau == 1.0, ]

# --- 6-design (retained designs only, relabeled) ---
mle_full_6 <- mle_full[mle_full$Design %in% retained_ids, ]
mle_full_6$Design <- full_name_map[mle_full_6$Design]
mle_tau1_6 <- mle_full_6[mle_full_6$True_Tau == 1.0, ]

# Reuse the already-verified 6-design test results (fr6/nem6/wil6/cond_tau6)
# rather than recomputing, so SI numbers match six_design_summary.txt exactly.
six_report <- readRDS(file.path(six_dir, "six_design_comparison_report.rds"))
fr6  <- six_report$friedman
nem6 <- six_report$nemenyi
wil6 <- six_report$wilcoxon

# Best-to-worst design ordering (pooled mean MSE at tau=1.0, all configs) --
# used for every reordered 6-design figure below. Matches the CTJ Table 2 order.
best_to_worst_6 <- mle_tau1_6 %>%
  group_by(Design) %>%
  summarise(Avg_MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
  arrange(Avg_MSE) %>%
  pull(Design)

wrap_labels <- function(x, width = 14) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1))
}

# ==============================================================================
# (a) FRESH 8-DESIGN CONSOLIDATION SUMMARY
#     Reuses code/12's corrected Nemenyi pair-count denominator: choose(k,2),
#     NOT nem$n_blocks (which is mislabeled upstream -- see 10's comment).
# ==============================================================================

sink(file.path(eight_dir, "eight_design_summary.txt"), split = TRUE)

cat("======================================================================\n")
cat("EIGHT-DESIGN STATISTICAL COMPARISONS (consolidation justification)\n")
cat("Purpose: justify dropping Saturation Quadrants (redundant with\n")
cat("Incidence-Guided Saturation Quadrants) and Balanced Halves (redundant\n")
cat("with Balanced Quartiles) from the 6-design manuscript comparison.\n")
cat("======================================================================\n\n")

cat("--- Primary scenario: True_Tau = 1.0, all 8 designs ---\n")
fr8  <- run_friedman_test(mle_tau1_8)
nem8 <- run_nemenyi_posthoc(mle_tau1_8)
wil8 <- run_pairwise_wilcoxon(mle_tau1_8)

cat(sprintf("Friedman chi-sq = %.2f, df = %d, p = %s, n_blocks = %d\n",
            fr8$statistic, fr8$n_designs - 1, format.pval(fr8$p_value, digits = 3),
            fr8$n_blocks))
cat("Average ranks (lower = better MSE):\n")
print(round(fr8$avg_ranks, 3))

n_designs8 <- length(fr8$avg_ranks)
n_sig_nem8 <- sum(nem8$sig_matrix[upper.tri(nem8$sig_matrix)])
cat(sprintf("\nNemenyi critical difference (alpha=0.05) = %.3f\n", nem8$critical_diff))
cat(sprintf("%d of %d pairs significantly different (Nemenyi)\n",
            n_sig_nem8, choose(n_designs8, 2)))

n_sig_wil8 <- sum(wil8$sig_matrix[upper.tri(wil8$sig_matrix)])
cat(sprintf("%d of %d pairs significantly different (Wilcoxon, Holm-adjusted)\n",
            n_sig_wil8, wil8$n_tests))

desc8 <- aggregate(cbind(MSE, Coverage) ~ Design, data = mle_tau1_8, FUN = mean)
desc8 <- desc8[order(desc8$MSE), ]
desc8$Rank <- seq_len(nrow(desc8))
cat("\nAll-8 mean MSE / Coverage / Rank (tau=1.0, primary scenario):\n")
print(desc8, row.names = FALSE)

# The two specific consolidation comparisons
cat("\n--- Consolidation pair 1: Saturation Quadrants vs Incidence-Guided Saturation Quadrants ---\n")
mse_sq  <- desc8$MSE[desc8$Design == "Saturation Quadrants"]
mse_isq <- desc8$MSE[desc8$Design == "Incidence-Guided Saturation Quadrants"]
p_nem_sq_isq <- nem8$p_matrix["Saturation Quadrants", "Incidence-Guided Saturation Quadrants"]
p_wil_sq_isq <- wil8$p_matrix["Saturation Quadrants", "Incidence-Guided Saturation Quadrants"]
cat(sprintf("MSE: Saturation Quadrants = %.4f, Incidence-Guided Saturation Quadrants = %.4f (delta = %.4f)\n",
            mse_sq, mse_isq, mse_sq - mse_isq))
cat(sprintf("Nemenyi p = %.4f | Wilcoxon (Holm) p = %s\n",
            p_nem_sq_isq, format.pval(p_wil_sq_isq, digits = 3)))

cat("\n--- Consolidation pair 2: Balanced Halves vs Balanced Quartiles ---\n")
mse_bh <- desc8$MSE[desc8$Design == "Balanced Halves"]
mse_bq <- desc8$MSE[desc8$Design == "Balanced Quartiles"]
p_nem_bh_bq <- nem8$p_matrix["Balanced Halves", "Balanced Quartiles"]
p_wil_bh_bq <- wil8$p_matrix["Balanced Halves", "Balanced Quartiles"]
cat(sprintf("MSE: Balanced Halves = %.4f, Balanced Quartiles = %.4f (delta = %.4f)\n",
            mse_bh, mse_bq, mse_bh - mse_bq))
cat(sprintf("Nemenyi p = %.4f | Wilcoxon (Holm) p = %s\n",
            p_nem_bh_bq, format.pval(p_wil_bh_bq, digits = 3)))

cat("\nFull Nemenyi p-value matrix (8 designs):\n")
print(round(nem8$p_matrix, 4))
cat("\nFull Wilcoxon (Holm-adjusted) p-value matrix (8 designs):\n")
print(round(wil8$p_matrix, 4))

sink()
cat("Summary written to", file.path(eight_dir, "eight_design_summary.txt"), "\n")

saveRDS(
  list(friedman = fr8, nemenyi = nem8, wilcoxon = wil8, descriptive = desc8),
  file.path(eight_dir, "eight_design_comparison_report.rds")
)

# ==============================================================================
# (b) MAIN-TEXT FIGURES: reordered MSE bar chart + new ranked bias-variance
# ==============================================================================

p_mse_reordered <- plot_master_comparison(mle_tau1_6, nb_filter = "queen",
                                          inc_label = "all configs, tau=1.0") +
  scale_x_discrete(limits = best_to_worst_6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(six_dir, "fig_mse_by_design_6design.pdf"), p_mse_reordered,
       width = 12, height = 7)

# SUD-relevant representative config for the bias-variance decomposition
configs_tau1_6 <- split_by_incidence_config(mle_tau1_6)
rep_label_6 <- "Poisson (rho_X = 0.20)"
results_rep_6 <- configs_tau1_6[[rep_label_6]]

p_biasvar <- plot_bias_variance(results_rep_6, rep_label_6) +
  scale_x_discrete(limits = best_to_worst_6, labels = wrap_labels) +
  labs(title = "Bias-Variance Decomposition of MSE (ranked best to worst)")
ggsave(file.path(six_dir, "fig_biasvar_6design.pdf"), p_biasvar, width = 9, height = 6)

# ==============================================================================
# (c) SI FIGURES (6 designs, full names, ranking-oriented)
# ==============================================================================

# --- SI Fig 1: CD diagram by average rank (lifted from
#     paper/report/IncidenceDesign_ProjectSummary.qmd inline code, ~lines 389-487,
#     NOT the buggy short-label plot_cd_diagram() in 10_statistical_comparisons.R) ---
{
  avg_ranks <- nem6$avg_ranks
  cd_val <- nem6$critical_diff
  alpha_val <- nem6$alpha
  designs <- names(avg_ranks)
  n <- length(designs)

  groups <- find_cd_groups(avg_ranks, nem6$p_matrix, alpha_val)

  df_pts <- data.frame(design = designs, rank = avg_ranks[designs],
                       stringsAsFactors = FALSE)
  df_pts <- df_pts[order(df_pts$rank), ]
  rownames(df_pts) <- NULL
  # Two-line label (wrapped name + rank value) -- full manuscript names are much
  # wider than the "D1"-style labels the original ProjectSummary.qmd code was
  # tuned for, so a single-line label collides with neighbors at this figure width.
  df_pts$label <- sprintf("%s\n(%.2f)", wrap_labels(df_pts$design, 16), df_pts$rank)

  y_levels <- c(1.3, -1.3, 1.9, -1.9)
  df_pts$y_label <- NA_real_
  for (i in seq_len(nrow(df_pts))) {
    far_enough <- if (i == 1) TRUE else min(df_pts$rank[i] - df_pts$rank[seq_len(i - 1)]) > 1.5
    if (far_enough) {
      df_pts$y_label[i] <- y_levels[((i - 1) %% 4) + 1]
    } else {
      nearby_idx <- which(df_pts$rank[i] - df_pts$rank[seq_len(i - 1)] < 1.5)
      used_levels <- df_pts$y_label[nearby_idx]
      available <- setdiff(y_levels, used_levels)
      df_pts$y_label[i] <- if (length(available) > 0) available[1] else y_levels[1]
    }
  }

  df_bars <- data.frame(xmin = numeric(0), xmax = numeric(0), y = numeric(0))
  if (length(groups) > 0) {
    bar_y_start <- 2.7
    bar_y_step <- 0.4
    for (g in seq_along(groups)) {
      grp_ranks <- avg_ranks[groups[[g]]]
      df_bars <- rbind(df_bars, data.frame(
        xmin = min(grp_ranks) - 0.05, xmax = max(grp_ranks) + 0.05,
        y    = bar_y_start + (g - 1) * bar_y_step))
    }
  }

  p_cd_rank <- ggplot() +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_segment(data = df_pts,
                 aes(x = rank, xend = rank, y = 0, yend = sign(y_label) * 0.25),
                 color = "grey40", linewidth = 0.3) +
    geom_point(data = df_pts, aes(x = rank, y = 0), size = 3, color = "black") +
    geom_segment(data = df_pts,
                 aes(x = rank, xend = rank, y = 0.05, yend = y_label * 0.55),
                 color = "grey50", linewidth = 0.2, linetype = "dotted") +
    geom_text(data = df_pts, aes(x = rank, y = y_label * 0.75, label = label),
              size = 2.6, fontface = "bold", hjust = 0.5, lineheight = 0.85) +
    {if (nrow(df_bars) > 0)
      geom_segment(data = df_bars, aes(x = xmin, xend = xmax, y = y, yend = y),
                   linewidth = 2.8, color = "grey35", lineend = "round")
    } +
    annotate("segment", x = 1, xend = 1 + cd_val, y = -2.9, yend = -2.9,
             linewidth = 1.0, color = "red3",
             arrow = arrow(ends = "both", length = unit(0.08, "inches"))) +
    annotate("text", x = 1 + cd_val / 2, y = -3.4,
             label = sprintf("CD = %.3f", cd_val), size = 2.9, color = "red3",
             fontface = "bold") +
    scale_x_continuous(breaks = 1:n, limits = c(0.5, n + 0.5),
                       name = "Average Rank (lower = better MSE)") +
    coord_cartesian(ylim = c(-4.0, max(4, max(df_bars$y, 0) + 1))) +
    theme_minimal(base_size = 11) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.margin = margin(5, 10, 5, 10)) +
    ggtitle("Critical Difference Diagram by Average Rank (6-design set, tau=1.0)")

  ggsave(file.path(si_dir, "si_fig_cd_diagram_rank.pdf"), p_cd_rank, width = 10, height = 5.2)
}

# --- SI Fig 2: CD diagram by mean MSE (log scale) ---
{
  mse_vals <- mle_tau1_6 %>%
    group_by(Design) %>%
    summarize(Mean_MSE = mean(MSE), .groups = "drop") %>%
    tibble::deframe()

  avg_ranks <- nem6$avg_ranks
  p_mat <- nem6$p_matrix
  designs <- names(avg_ranks)
  groups <- find_cd_groups(avg_ranks, p_mat, nem6$alpha)

  df_pts <- data.frame(design = designs, mse = mse_vals[designs], stringsAsFactors = FALSE)
  df_pts <- df_pts[order(df_pts$mse), ]
  rownames(df_pts) <- NULL
  # Two-line label (wrapped name + MSE value) -- see rank-CD comment above for why.
  df_pts$label <- sprintf("%s\n(%.4f)", wrap_labels(df_pts$design, 16), df_pts$mse)

  y_levels <- c(1.0, -1.0, 1.8, -1.8, 2.6, -2.6)
  log_mse <- log10(df_pts$mse)
  label_y <- numeric(nrow(df_pts))
  label_y[1] <- y_levels[1]
  for (i in seq_len(nrow(df_pts))[-1]) {
    nearby_idx <- which((log_mse[i] - log_mse[seq_len(i - 1)]) < 0.30)
    if (length(nearby_idx) == 0) {
      label_y[i] <- y_levels[((i - 1) %% length(y_levels)) + 1]
    } else {
      used <- label_y[nearby_idx]
      avail <- setdiff(y_levels, used)
      label_y[i] <- if (length(avail) > 0) avail[1] else y_levels[1]
    }
  }
  df_pts$label_y <- label_y

  bar_base <- 3.0
  min_half_log <- 0.03
  df_bars <- do.call(rbind, lapply(seq_along(groups), function(g) {
    grp_mse <- mse_vals[groups[[g]]]
    log_mid <- mean(log10(c(min(grp_mse), max(grp_mse))))
    log_xmin <- min(log10(min(grp_mse)), log_mid - min_half_log)
    log_xmax <- max(log10(max(grp_mse)), log_mid + min_half_log)
    data.frame(xmin = 10^log_xmin, xmax = 10^log_xmax, y = bar_base + (g - 1) * 0.18)
  }))

  p_cd_mse <- ggplot() +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_segment(data = df_pts,
                 aes(x = mse, xend = mse, y = 0, yend = sign(label_y) * 0.15),
                 color = "grey40", linewidth = 0.3) +
    geom_point(data = df_pts, aes(x = mse, y = 0), size = 3.5, color = "black") +
    geom_segment(data = df_pts,
                 aes(x = mse, xend = mse, y = sign(label_y) * 0.20, yend = label_y * 0.75),
                 color = "grey50", linewidth = 0.2, linetype = "dotted") +
    geom_text(data = df_pts, aes(x = mse, y = label_y, label = label),
              size = 2.6, fontface = "bold", hjust = 0.5, vjust = 0.5, lineheight = 0.85) +
    {if (!is.null(df_bars) && nrow(df_bars) > 0)
      geom_segment(data = df_bars, aes(x = xmin, xend = xmax, y = y, yend = y),
                   linewidth = 3, color = "grey35", lineend = "round")
    } +
    scale_x_log10(name = "Mean MSE (log scale, tau = 1.0, all incidence configurations)",
                  labels = scales::label_number(accuracy = 0.001),
                  expand = expansion(mult = c(0.18, 0.18))) +
    coord_cartesian(clip = "off",
                    ylim = c(-3.1, max(if (!is.null(df_bars) && nrow(df_bars) > 0) df_bars$y else 0,
                                       bar_base) + 0.5)) +
    theme_minimal(base_size = 10.5) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.margin = margin(3, 8, 2, 8)) +
    ggtitle("Design MSE with Statistical Equivalence Groups (6-design set, tau=1.0)")

  ggsave(file.path(si_dir, "si_fig_cd_diagram_mse.pdf"), p_cd_mse, width = 11, height = 6.5)
}

# --- SI Fig 3: MSE heatmap (design x rho x neighbor), rows best->worst ---
p_heatmap_mse <- plot_mse_heatmap(results_rep_6, rep_label_6) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "Avg MSE\n(log scale)") +
  scale_y_discrete(limits = rev(best_to_worst_6), labels = wrap_labels)
ggsave(file.path(si_dir, "si_fig_mse_heatmap.pdf"), p_heatmap_mse, width = 8, height = 4.5)

# --- SI Fig 4: average-rank heatmap (design x rho x neighbor), rows best->worst ---
rank_heat_df <- results_rep_6 %>%
  group_by(Rho, Neighbor_Type, Gamma, Spillover_Type) %>%
  mutate(Rank = rank(MSE, ties.method = "average")) %>%
  ungroup() %>%
  group_by(Design, Rho, Neighbor_Type) %>%
  summarize(Avg_Rank = mean(Rank), .groups = "drop")

p_heatmap_rank <- ggplot(rank_heat_df,
                         aes(x = as.factor(Rho),
                             y = factor(Design, levels = rev(best_to_worst_6)),
                             fill = Avg_Rank)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", Avg_Rank)), size = 3.5, fontface = "bold",
            color = "white") +
  facet_wrap(~ Neighbor_Type,
             labeller = as_labeller(c("queen" = "Queen Contiguity", "rook" = "Rook Contiguity"))) +
  scale_fill_distiller(palette = "RdYlGn", direction = -1, limits = c(1, 6),
                       name = "Avg Rank\n(1 = best)") +
  scale_y_discrete(labels = wrap_labels) +
  theme_minimal(base_size = 12) +
  labs(title = expression("Design Rank by Spatial Dependence (" * rho * ")"),
       subtitle = paste0("Incidence: ", rep_label_6), x = expression("Spatial Dependence (" * rho * ")"),
       y = NULL) +
  theme(panel.grid = element_blank())
ggsave(file.path(si_dir, "si_fig_rank_heatmap.pdf"), p_heatmap_rank, width = 8, height = 4.5)

# --- SI Fig 5: two-panel ranked mean-MSE bar (colored by group) + p-value heatmap ---
mse_bar_df <- mle_tau1_6 %>%
  group_by(Design) %>%
  summarize(Mean_MSE = mean(MSE), SE = sd(MSE) / sqrt(n()), .groups = "drop")

p_bar <- ggplot(mse_bar_df,
                aes(x = Mean_MSE, y = factor(Design, levels = rev(best_to_worst_6)))) +
  geom_col(aes(fill = DESIGN_GROUPS_6[as.character(Design)]), alpha = 0.85) +
  geom_errorbarh(aes(xmin = pmax(Mean_MSE - SE, 0.001), xmax = Mean_MSE + SE), height = 0.35) +
  scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_discrete(labels = wrap_labels) +
  scale_fill_manual(values = c("Blocking" = "#D55E00", "Stratified" = "#0072B2",
                               "Saturation" = "#009E73"), name = "Design type") +
  labs(x = "Mean MSE (log scale)", y = NULL, title = "Design Performance (6-design set)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_heat_pval <- plot_pvalue_heatmap(nem6, "Nemenyi") +
  scale_x_discrete(labels = wrap_labels) +
  scale_y_discrete(labels = wrap_labels) +
  theme(axis.text = element_text(size = 7.5), plot.title = element_text(size = 11))

p_twopanel <- p_bar + p_heat_pval + plot_layout(widths = c(1, 1.4))
ggsave(file.path(si_dir, "si_fig_performance_pvalue_twopanel.pdf"), p_twopanel,
       width = 11, height = 4.5)

# --- SI Fig 6: tau-sweep mean-MSE line + 95% coverage line (6 designs) ---
mle_full_6_leveled <- mle_full_6
mle_full_6_leveled$Design <- factor(mle_full_6_leveled$Design, levels = best_to_worst_6)

p_tau_mse <- plot_mse_vs_tau(mle_full_6_leveled) +
  labs(title = "MSE vs. True Tau by Design (6-design set)")
p_tau_cov <- plot_coverage_vs_tau(mle_full_6_leveled) +
  labs(title = "Coverage vs. True Tau by Design (6-design set)")

ggsave(file.path(si_dir, "si_fig_tau_mse.pdf"), p_tau_mse, width = 8, height = 5)
ggsave(file.path(si_dir, "si_fig_tau_coverage.pdf"), p_tau_cov, width = 8, height = 5)

# ==============================================================================
# (d) 6-DESIGN APPLICATION TABLE
# ==============================================================================

app_rds_path <- file.path(proj_dir, "application", "results", "full",
                          "application_full_results.rds")
app_results <- readRDS(app_rds_path)
app_summary <- app_results$summary_results

app_sub <- app_summary[app_summary$design_id %in% c(1, 2, 4, 5, 6, 8) &
                        app_summary$tau == 1.0, ]
app_agg <- aggregate(cbind(mse, coverage, power) ~ design_name, data = app_sub, FUN = mean)
app_agg <- app_agg[order(app_agg$mse), ]

sink(file.path(six_dir, "application_table_6design.txt"), split = TRUE)
cat("======================================================================\n")
cat("6-DESIGN APPLICATION-SCALE TABLE (NC 58-cluster community college network)\n")
cat("Source: application/results/full/application_full_results.rds$summary_results\n")
cat("Filtered to design_id in {1,2,4,5,6,8} (the 6 retained manuscript designs),\n")
cat("tau = 1.0, averaged over gamma/rho/spillover_type.\n")
cat("\n")
cat("NOTE (relabeling): the application study's design_id 1 is named\n")
cat("'Block Stratified Sampling' (application_designs.R), not 'Checkerboard' --\n")
cat("the Checkerboard concept was adapted to a block-stratified analog for the\n")
cat("irregular NC geography. Design_id 8 is similarly named 'Incidence-Guided\n")
cat("Saturation Regions' in the application code (vs. 'Incidence-Guided\n")
cat("Saturation Quadrants' in the grid study) for the same reason -- 'regions'\n")
cat("(k-means clusters) replace the literal 2x2 grid quadrants on the irregular\n")
cat("map. Both are the application study's own names, used verbatim below.\n")
cat("======================================================================\n\n")
print(app_agg, row.names = FALSE)
sink()
cat("Application table written to", file.path(six_dir, "application_table_6design.txt"), "\n")

cat("\n=== 14_manuscript_supplement_figures.R complete ===\n")
cat("8-design summary:  ", file.path(eight_dir, "eight_design_summary.txt"), "\n")
cat("Main-text figures: ", six_dir, "\n")
cat("SI figures:        ", si_dir, "\n")
cat("Application table: ", file.path(six_dir, "application_table_6design.txt"), "\n")
