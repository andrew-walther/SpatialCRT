# ============================================================
# Script: 28_thesis_chapter_figures.R
# Purpose: Re-plot the Chapter 3 figures whose text printed below 8 pt, saved at
#          their printed width with embedded fonts, from saved results only.
# Author: Andrew Walther
# Created: 2026-11-10
# Dependencies: dplyr, ggplot2, patchwork (via 10_statistical_comparisons.R);
#               00_thesis_figure_device.R
# ============================================================
#
# Run from code/:  Rscript 28_thesis_chapter_figures.R
#
# Why a separate script: the manuscript figures (CTJ and SI) come from 12, 14, 15,
# 20, 21, 23 and real_sud_summary.R and are copied into both paper folders by
# paper/tools/copy_figures.sh. They are sized for the journal page. The chapter
# prints figures at 6.5 in or narrower, so the same PDFs shrink and their text falls
# to 4-7 pt. This script rebuilds the chapter's copies at the printed width with
# every text element >= 8 pt, and writes them straight into
# paper/thesis_chapter/figures/ (copy_figures.sh no longer copies those names there).
# The CTJ and SI files are untouched.
#
# Data: re-plots only. Nothing is simulated or re-fitted; no results/ file is
# written. Inputs:
#   results/sim_data/sim_results_MLE_tau_sweep_combined_*.rds  (latest)
#   results/grid_tdf_2026-10/grid_tdf_scenarios.rds            (t + df intervals)
#   results/six_design_manuscript/six_design_comparison_report.rds (Nemenyi tests)
# Plot content (data, panels, designs, order, colours) follows the originals:
# code/14 (mse ranked, bias-variance, CD diagram, tau coverage) and code/12
# (coverage by design, combined tau sensitivity).

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
setwd(script_dir)
results_dir <- file.path(dirname(script_dir), "results")  # before sourcing 10, which overwrites script_dir
proj_dir    <- dirname(script_dir)
fig_out     <- file.path(proj_dir, "paper", "thesis_chapter", "figures")

source("10_statistical_comparisons.R")  # sources 06_visualizations.R -> 01-03
source("00_thesis_figure_device.R")

#' Replace ASCII hyphens with soft hyphens (prints a true hyphen under pdf())
#' @param x Character vector of labels
#' @return Labels with "-" replaced by U+00AD
hy <- function(x) gsub("-", "­", x, fixed = TRUE)

#' Wrap labels at `width` characters (as in 14_manuscript_supplement_figures.R)
#' @param x Character vector
#' @param width Wrap width
#' @return Character vector with "\n" line breaks
wrap_labels <- function(x, width = 14) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

MS_CLEAN <- theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                  plot.caption = element_blank())

PRINT_W <- 6.5   # \linewidth in the prelim class (in)

# Data ----

raw_names <- get_design_names()
full_name_map <- setNames(sub("^Design [0-9]+: ", "", raw_names),
                          paste0("Design ", seq_along(raw_names)))
full_name_map[["Design 8"]] <- "Incidence-Guided Saturation"
full_name_map[["Design 5"]] <- "2×2 Blocking"
retained_ids <- c("Design 1", "Design 2", "Design 4", "Design 5", "Design 6", "Design 8")

mle_full <- load_latest_results(results_dir = results_dir, estimation_mode = "MLE_tau_sweep")
mle_full <- apply_tdf_intervals(mle_full, results_dir)   # reported interval: t + df
mle6_all <- mle_full[mle_full$Design %in% retained_ids, ]
mle6_all$Design <- full_name_map[mle6_all$Design]
srs_all  <- mle_full[mle_full$Design == "Design 9", ]
six_report <- readRDS(file.path(results_dir, "six_design_manuscript",
                                "six_design_comparison_report.rds"))

tau_colour_order <- c("Incidence-Guided Saturation", "Balanced Quartiles",
                      "Isolation Buffer", "High Incidence Focus",
                      "2×2 Blocking", "Checkerboard")
tau_palette <- setNames(scales::viridis_pal(option = "D", end = 0.85)(6), tau_colour_order)

#' Save one chapter figure and report its path
#' @param p Plot object
#' @param name File name inside paper/thesis_chapter/figures
#' @param w,h Printed size in inches
#' @return Invisibly, the save method
save_fig <- function(p, name, w, h) {
  m <- save_thesis_pdf(p, file.path(fig_out, name), w, h)
  cat(sprintf("%-48s %.2f x %.2f in  (%s)\n", name, w, h, m))
  invisible(m)
}

#' Gray SRS reference layer for the tau-sensitivity (MSE) figure
#' @param srs SRS rows for one neighbor type, all tau
#' @return A geom_line layer: mean MSE by tau
srs_tau_layer <- function(srs) {
  s <- aggregate(MSE ~ True_Tau, data = srs, FUN = mean)
  geom_line(data = s, aes(x = True_Tau, y = MSE), inherit.aes = FALSE,
            linetype = "dashed", colour = "grey20", linewidth = 0.6)
}

tau_panel_list <- list()

for (nb in c("queen", "rook")) {
  mle6    <- mle6_all[mle6_all$Neighbor_Type == nb, ]
  mle6_t1 <- mle6[mle6$True_Tau == 1.0, ]
  srs     <- srs_all[srs_all$Neighbor_Type == nb, ]
  srs_t1  <- srs[srs$True_Tau == 1.0, ]
  best_to_worst <- mle6_t1 %>% group_by(Design) %>%
    summarise(m = mean(MSE, na.rm = TRUE), .groups = "drop") %>% arrange(m) %>% pull(Design)
  mle6_lev <- mle6
  mle6_lev$Design <- factor(mle6_lev$Design, levels = best_to_worst)

  ## MSE by design: Figure 3.2 (queen), Figure B.5 (rook) ----
  p_mse <- plot_master_comparison(mle6_t1, nb_filter = nb) +
    scale_x_discrete(limits = best_to_worst, labels = hy) +
    facet_grid(Spillover_Type ~ Rho, scales = "free_y",
               labeller = labeller(Spillover_Type = c("both" = "Spillover: Both",
                                                      "control_only" = "Spillover: Control"),
                                   Rho = as_labeller(function(x) paste("rho ==", x), label_parsed))) +
    geom_hline(data = srs_t1 %>%
                 group_by(Spillover_Type, Rho, Gamma) %>%
                 summarise(m = max(MSE), .groups = "drop") %>%
                 group_by(Spillover_Type, Rho) %>%
                 summarise(yint = mean(m), .groups = "drop"),
               aes(yintercept = yint), linetype = "dashed", colour = "grey40",
               linewidth = 0.5, inherit.aes = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    thesis_text_sizes() + MS_CLEAN +
    guides(fill = guide_legend(nrow = 1))
  save_fig(p_mse, sprintf("fig_mse_by_design_ranked_6design_%s.pdf", nb), PRINT_W, 6.4)

  ## Bias-variance decomposition: Figure B.4 (queen), Figure B.7 (rook) ----
  rep_label <- inc_config_label("poisson", 0.20)
  results_rep <- split_by_incidence_config(mle6_t1)[[rep_label]]
  stopifnot(rep_label == "Poisson (rho_X = 0.20)")
  p_bv <- plot_bias_variance(results_rep, rep_label) +
    scale_x_discrete(limits = best_to_worst, labels = function(x) hy(wrap_labels(x))) +
    scale_fill_manual(values = c("Bias-squared" = "#E69F00", "Variance" = "#56B4E9"),
                      labels = hy) +
    geom_hline(yintercept = with(srs_t1[srs_t1$Incidence_Mode == "poisson" &
                                          abs(srs_t1$Rho_Incidence - 0.20) < 1e-9, ],
                                 mean(Bias^2) + mean(SD^2)),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    thesis_text_sizes() + MS_CLEAN
  save_fig(p_bv, sprintf("fig_biasvar_6design_%s.pdf", nb), 0.9 * PRINT_W, 4.4)

  ## Coverage by design: Figure B.3 (queen), Figure B.6 (rook) ----
  srs_cov_median <- median(srs_t1$Coverage)
  p_cov <- plot_coverage_by_design(mle6_t1) +
    geom_hline(yintercept = srs_cov_median, linetype = "dotted", colour = "grey30",
               linewidth = 0.5) +
    scale_x_discrete(labels = hy) +
    thesis_text_sizes() + MS_CLEAN
  p_cov$layers[[1]] <- geom_boxplot(alpha = 0.5, outlier.alpha = 0.5, outlier.size = 0.8,
                                    colour = "grey35", linewidth = 0.25)
  p_cov <- p_cov + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p_cov, sprintf("fig_coverage_by_design_6design_%s.pdf", nb), 0.9 * PRINT_W, 4.4)

  ## Coverage across tau: Figure B.2a (queen), B.2b (rook) ----
  srs_cov_tau <- aggregate(Coverage ~ True_Tau, data = srs, FUN = mean)
  p_ctau <- plot_coverage_vs_tau(mle6_lev) +
    geom_line(data = srs_cov_tau, aes(x = True_Tau, y = Coverage), inherit.aes = FALSE,
              linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    scale_color_manual(values = tau_palette, breaks = best_to_worst, labels = hy) +
    scale_fill_manual(values = tau_palette, breaks = best_to_worst, labels = hy) +
    labs(x = expression("True " * tau)) +
    theme_bw(base_size = 10) + thesis_text_sizes() + MS_CLEAN +
    theme(legend.position = "right")
  save_fig(p_ctau, sprintf("app_tau_coverage_%s.pdf", nb), PRINT_W, 3.5)

  ## Critical-difference diagram by average rank: Figure B.1a (queen), B.1b (rook) ----
  nem6 <- six_report[[nb]]$nemenyi
  stopifnot("Simple Random Sampling" %in% names(nem6$avg_ranks))
  avg_ranks <- nem6$avg_ranks
  cd_val <- nem6$critical_diff
  designs <- names(avg_ranks)
  n <- length(designs)
  groups <- find_cd_groups(avg_ranks, nem6$p_matrix, nem6$alpha)
  df_pts <- data.frame(design = designs, rank = avg_ranks[designs], stringsAsFactors = FALSE)
  df_pts <- df_pts[order(df_pts$rank), ]
  rownames(df_pts) <- NULL
  df_pts$label <- hy(sprintf("%s\n(%.2f)", wrap_labels(df_pts$design, 16), df_pts$rank))
  y_levels <- c(1.3, -1.3, 2.6, -2.6)
  df_pts$y_label <- NA_real_
  for (i in seq_len(nrow(df_pts))) {
    far_enough <- if (i == 1) TRUE else min(df_pts$rank[i] - df_pts$rank[seq_len(i - 1)]) > 1.5
    if (far_enough) {
      df_pts$y_label[i] <- y_levels[((i - 1) %% 4) + 1]
    } else {
      nearby_idx <- which(df_pts$rank[i] - df_pts$rank[seq_len(i - 1)] < 1.5)
      available <- setdiff(y_levels, df_pts$y_label[nearby_idx])
      df_pts$y_label[i] <- if (length(available) > 0) available[1] else y_levels[1]
    }
  }
  df_bars <- data.frame(xmin = numeric(0), xmax = numeric(0), y = numeric(0))
  for (g in seq_along(groups)) {
    grp_ranks <- avg_ranks[groups[[g]]]
    df_bars <- rbind(df_bars, data.frame(xmin = min(grp_ranks) - 0.05,
                                         xmax = max(grp_ranks) + 0.05,
                                         y = 2.7 + (g - 1) * 0.4))
  }
  TXT_MM <- THESIS_TEXT_PT / ggplot2::.pt   # 8 pt in mm
  p_cd <- ggplot() +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_segment(data = df_pts, aes(x = rank, xend = rank, y = 0, yend = sign(y_label) * 0.25),
                 color = "grey40", linewidth = 0.3) +
    geom_point(data = df_pts, aes(x = rank, y = 0), size = 2.5, color = "black") +
    geom_segment(data = df_pts, aes(x = rank, xend = rank, y = 0.05, yend = y_label * 0.55),
                 color = "grey50", linewidth = 0.2, linetype = "dotted") +
    geom_text(data = df_pts, aes(x = rank, y = y_label * 0.75, label = label),
              size = TXT_MM, fontface = "bold", hjust = 0.5, lineheight = 0.85) +
    {if (nrow(df_bars) > 0)
      geom_segment(data = df_bars, aes(x = xmin, xend = xmax, y = y, yend = y),
                   linewidth = 2.2, color = "grey35", lineend = "round")} +
    annotate("segment", x = 1, xend = 1 + cd_val, y = -2.65, yend = -2.65,
             linewidth = 0.8, color = "red3",
             arrow = arrow(ends = "both", length = unit(0.07, "inches"))) +
    annotate("text", x = 1 + cd_val / 2, y = -3.05, label = sprintf("CD = %.3f", cd_val),
             size = TXT_MM, color = "red3", fontface = "bold") +
    scale_x_continuous(breaks = 1:n, limits = c(0.5, n + 0.5),
                       name = "Average Rank (lower = better MSE)") +
    coord_cartesian(ylim = c(-3.3, max(3.3, max(df_bars$y, 0) + 0.3))) +
    theme_minimal(base_size = 10) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.margin = margin(5, 10, 5, 10)) +
    thesis_text_sizes() + MS_CLEAN
  save_fig(p_cd, sprintf("app_cd_diagram_rank_%s.pdf", nb), PRINT_W, 3.6)

  ## Tau-sensitivity panel for the combined Figure 3.5 ----
  res_tau <- mle6
  res_tau$Design <- factor(res_tau$Design, levels = tau_colour_order)
  tau_panel_list[[nb]] <- plot_mse_vs_tau(res_tau) + srs_tau_layer(srs) +
    scale_color_manual(values = tau_palette, breaks = tau_colour_order, labels = hy) +
    scale_fill_manual(values = tau_palette, breaks = tau_colour_order, labels = hy) +
    labs(title = c(queen = "(a) Queen", rook = "(b) Rook")[[nb]],
         x = expression("True " * tau),
         y = if (nb == "queen") "Mean MSE\n(band: ± average scenario-level\nMonte Carlo SE)" else NULL,
         color = NULL, fill = NULL) +
    theme_bw(base_size = 10)
}

## Figure 3.5: MSE across tau, queen | rook ----
p_tau <- (tau_panel_list$queen | tau_panel_list$rook) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.key.size = unit(0.9, "lines"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) &
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))
p_tau <- p_tau & thesis_text_sizes() & theme(plot.title = element_text(size = 10))
save_fig(p_tau, "fig_tau_sensitivity_6design_combined.pdf", PRINT_W, 3.9)

# SRS comparison forest plot: Figure 3.3 ----
# Same plot as plot_vs_srs() in 15_srs_benchmark.R (queen, oracle, eight designs shown
# by regime), from the stored paired tests in results/srs_benchmark/srs_benchmark.rds.
bench <- readRDS(file.path(results_dir, "srs_benchmark", "srs_benchmark.rds"))
tests <- bench$oracle$queen$tests_by_regime
display <- c("Incidence-Guided Saturation Quadrants" = "Incidence-Guided Saturation",
             "2x2 Blocking" = "2×2 Blocking")
d <- tests[tests$Design %in% c(unname(full_name_map[c("Design 1", "Design 2", "Design 4", "Design 5",
                                                      "Design 6", "Design 8")]),
                               "Saturation Quadrants", "Balanced Halves",
                               "Incidence-Guided Saturation Quadrants", "2x2 Blocking",
                               "Incidence-Guided Saturation"), ]
d <- d[!duplicated(paste(d$Design, d$Spillover_Type)), ]
stopifnot(nrow(d) == 16)
d$Design <- ifelse(d$Design %in% names(display), display[d$Design], d$Design)
d$ratio <- 1 + d$mean_diff / d$MSE_SRS
d$lo    <- 1 + d$ci_lo / d$MSE_SRS
d$hi    <- 1 + d$ci_hi / d$MSE_SRS
d$lab   <- ifelse(d$ratio > 2, sprintf("%.1f", d$ratio), "")
SIG <- hy("Holm-adjusted Wilcoxon p < 0.05")
d$sig <- ifelse(d$p_holm < 0.05, SIG, "Not significant")
d$Regime <- factor(ifelse(d$Spillover_Type == "both", hy("Both-arms spillover"),
                          hy("Control-only spillover")),
                   levels = hy(c("Control-only spillover", "Both-arms spillover")))
ord <- d[d$Spillover_Type == "control_only", ]
d$Design <- factor(d$Design, levels = rev(ord$Design[order(ord$ratio)]))
p_vs <- ggplot(d, aes(x = ratio, y = Design)) +
  geom_vline(xintercept = 1, colour = "grey40", linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y") +
  geom_point(aes(shape = sig), size = 2.4, fill = "black") +
  geom_text(aes(label = lab), hjust = 1.4, size = THESIS_TEXT_PT / ggplot2::.pt) +
  scale_shape_manual(values = setNames(c(21, 1), c(SIG, "Not significant")), name = NULL) +
  scale_y_discrete(labels = hy) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 4, 10),
                labels = c("0.5", "0.75", "1", "1.5", "2", "4", "10"), limits = c(0.5, 12)) +
  facet_wrap(~ Regime) +
  labs(x = "MSE relative to Simple Random Sampling (log scale; below 1 = lower error)",
       y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.spacing = unit(1.2, "lines")) +
  thesis_text_sizes()
save_fig(p_vs, "fig_vs_srs_by_regime_queen.pdf", PRINT_W, 3.7)

# Methods-extension figures: Figures 3.10, 3.11, B.12 ----
# Same plots as 23_methods_extension_figures.R, re-drawn from its stored CSV inputs.
R_EXT <- file.path(proj_dir, "application/results/methods_ext_2026-10")
DISPLAY <- c("Plain Saturation" = "Saturation Quadrants", "2x2 Blocking" = "2×2 Blocking")
display_name <- function(x) ifelse(x %in% names(DISPLAY), DISPLAY[x], x)
REGIME_LAB <- c(both = hy("Both-arms spillover"), control_only = hy("Control-only spillover"))
TXT_MM <- THESIS_TEXT_PT / ggplot2::.pt

## Power against the number of colleges (Figure 3.11) ----
settings <- read.csv(file.path(R_EXT, "summaries_v2/a3_settings.csv"))
x <- settings[settings$Mechanism == "random" & settings$Beta == 0 & settings$Kappa == 0.8 &
                settings$Rho == 0.5 & settings$Design_ID %in% c(9, 3, 8, 1) &
                settings$Tau %in% c(1, 1.5), ]
stopifnot(nrow(x) > 0, all(table(paste(x$Design_ID, x$Regime, x$Tau, x$N)) == 4))
agg <- aggregate(Power_Tdf ~ Design + Regime + Tau + N, x,
                 function(v) c(Power = mean(v), Lo = min(v), Hi = max(v)))
pc <- data.frame(agg[, c("Design", "Regime", "Tau", "N")], agg$Power_Tdf)
pc$Design <- display_name(pc$Design)
pc$Design <- factor(pc$Design, levels = c("Simple Random Sampling", "Saturation Quadrants",
                                          "Incidence-Guided Saturation", "Checkerboard"))
pc$Regime <- factor(REGIME_LAB[pc$Regime], levels = REGIME_LAB)
pc$TauLab <- factor(sprintf("Effect of %s college-level SD", pc$Tau),
                    levels = sprintf("Effect of %s college-level SD", c(1, 1.5)))
pal4 <- c("Simple Random Sampling" = "grey25", "Saturation Quadrants" = "#009E73",
          "Incidence-Guided Saturation" = "#0072B2", "Checkerboard" = "#D55E00")
p_power <- ggplot(pc, aes(N, Power, colour = Design, fill = Design, linetype = Design, shape = Design)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 15, y = 0.86, label = "80% power", hjust = 0, size = TXT_MM, colour = "grey35") +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  facet_grid(TauLab ~ Regime) +
  scale_colour_manual(values = pal4, name = NULL, labels = hy) +
  scale_fill_manual(values = pal4, name = NULL, labels = hy) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash"), name = NULL, labels = hy) +
  scale_shape_manual(values = c(16, 17, 15, 4), name = NULL, labels = hy) +
  scale_x_continuous(breaks = c(15, 20, 29, 40, 50, 58)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Number of participating colleges", y = "Power (chance the trial detects the effect)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  thesis_text_sizes()
save_fig(p_power, "fig_power_vs_n.pdf", PRINT_W, 5.0)

## Colleges needed for 80% power (Figure 3.10) ----
#' Colleges-needed data for one regime and effect (as n80_data() in 23)
#' @param n80 a3_n80_summary.csv as a data frame
#' @param regime "both" or "control_only"
#' @param tau Effect in college-level SD units
#' @return Data frame for the bars, ranges and labels
n80_data <- function(n80, regime, tau) {
  x <- n80[n80$Kappa == 0.8 & n80$Rho == 0.5 & n80$Regime == regime & n80$Tau == tau, ]
  stopifnot(nrow(x) == 9)
  reached <- is.finite(x$N80_Tdf_Mean)
  data.frame(
    Design = display_name(x$Design),
    Panel = sprintf("%s, effect of %s SD", REGIME_LAB[[regime]], tau),
    Value = ifelse(reached, ceiling(x$N80_Tdf_Mean), 58),
    Lo = ifelse(reached, ceiling(x$N80_Tdf_Min), NA), Hi = ifelse(reached, ceiling(x$N80_Tdf_Max), NA),
    Reached = ifelse(reached, "Reached in every year", "Not reached with 58 colleges in at least one year"),
    Label = ifelse(reached, sprintf("%d", as.integer(ceiling(x$N80_Tdf_Mean))),
                   ifelse(is.finite(x$N80_Tdf_Min),
                          sprintf("> 58 in %d/4 yr, else %d+", x$N80_Tdf_Years_Above58, as.integer(ceiling(x$N80_Tdf_Min))),
                          "> 58 in 4/4 yr")),
    Label_x = ifelse(reached, ceiling(x$N80_Tdf_Max), 58) + 1.5)
}
n80 <- read.csv(file.path(R_EXT, "summaries_v2/a3_n80_summary.csv"))
nd <- rbind(n80_data(n80, "both", 1), n80_data(n80, "control_only", 1.5))
nd$Panel <- factor(nd$Panel, levels = unique(nd$Panel))
ord <- nd[nd$Panel == levels(nd$Panel)[2], ]
nd$Design <- factor(nd$Design, levels = rev(ord$Design[order(ord$Value, ord$Design)]))
p_n80 <- ggplot(nd, aes(Value, Design, fill = Reached)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = Lo, xmax = Hi), width = 0.25, orientation = "y", na.rm = TRUE) +
  geom_text(aes(x = Label_x, label = Label), hjust = 0, size = TXT_MM) +
  geom_vline(xintercept = 58, linetype = "dotted", colour = "grey40") +
  facet_wrap(~ Panel) +
  scale_y_discrete(labels = hy) +
  scale_fill_manual(values = c("Reached in every year" = "grey30",
                               "Not reached with 58 colleges in at least one year" = "grey80"), name = NULL) +
  scale_x_continuous(limits = c(0, 118), breaks = c(0, 20, 29, 40, 58)) +
  labs(x = "Colleges needed for 80% power", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  thesis_text_sizes()
save_fig(p_n80, "fig_colleges_needed.pdf", PRINT_W, 3.6)

## Spread against incidence-guidance gain (Figure B.12) ----
setup <- readRDS(file.path(proj_dir, "application/results/real_sud_rev_20261002/setup.rds"))
a1 <- read.csv(file.path(R_EXT, "stored/a1_burden_share.csv"))
nc <- do.call(rbind, lapply(c("2018", "2019", "2020", "2021"), function(y) {
  ann <- setup$annual[[y]]
  stopifnot(identical(ann$Primary_College, setup$ids), nrow(ann) == 58)
  region_rank <- tapply(ann$rank, setup$regions, mean)
  stopifnot(length(region_rank) == 4)
  rs <- a1[a1$Year == as.integer(y), ]
  data.frame(Year = y, Spread = diff(range(region_rank)) / 58,
             Gain = rs$Rate_Share_Mean[rs$Design_ID == 8] - rs$Rate_Share_Mean[rs$Design_ID == 3])
}))
grid <- read.csv(file.path(results_dir, "grid_burden_share/grid_quadrant_spread.csv"))
CONFIG_LAB <- c("iid" = "iid", "spatial rho_X=0.2" = "Spatial, 0.2", "spatial rho_X=0.5" = "Spatial, 0.5",
                "poisson rho_X=0.2" = "Poisson, 0.2", "poisson rho_X=0.5" = "Poisson, 0.5")
stopifnot(all(grid$Config %in% names(CONFIG_LAB)))
grid$Config <- factor(CONFIG_LAB[grid$Config], levels = CONFIG_LAB)
p_spread <- ggplot(grid, aes(Spread, Gain)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50",
              linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(colour = Config, shape = Config), size = 2, alpha = 0.9) +
  geom_point(data = nc, shape = 24, size = 3, colour = "black", fill = "white") +
  geom_text(data = nc, aes(label = paste("NC", Year)), hjust = -0.2, vjust = 0.5, size = TXT_MM) +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
                      name = expression("Grid incidence, " * rho[displaystyle(X)])) +
  scale_shape_manual(values = c(16, 17, 15, 3, 7), name = expression("Grid incidence, " * rho[displaystyle(X)])) +
  scale_y_continuous(labels = function(v) ifelse(abs(v) < 1e-9, "0", sprintf("%+.0f", 100 * v))) +
  labs(x = "Spread of regional mean incidence ranks (max - min) / number of clusters",
       y = hy("Treated incidence share, Incidence-Guided\nSaturation minus Saturation Quadrants (points)")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank()) +
  thesis_text_sizes()
save_fig(p_spread, "fig_spread_vs_gain.pdf", 0.9 * PRINT_W, 3.9)

# Allocation distribution by design: Figure 3.9 ----
# Same plot as 21_downside_by_design.R, from its allocation_distribution_by_design.csv.
dist <- read.csv(file.path(results_dir, "allocation_risk/summary/surfaces10",
                           "allocation_distribution_by_design.csv"))
ids <- c(8, 6, 4, 2, 5, 1, 9)
labels <- c("Incidence­Guided Saturation", "Balanced Quartiles", "Isolation Buffer",
            "High Incidence Focus", "2×2 Blocking", "Checkerboard", "Simple Random Sampling")
f <- dist[dist$Design_ID %in% ids, ]
f$Label <- factor(f$Design_ID, levels = rev(ids), labels = rev(labels))
f$Study <- factor(f$Source, c("grid_surfaces10", "nc_refined"), c("Simulation study", "Application study"))
f$Regime <- factor(f$Regime, c("both", "control_only"),
                   c("Spillover reaches both arms", "Spillover reaches controls only"))
f$Benchmark <- f$Design_ID == 9
stopifnot(!anyNA(f$Study), !anyNA(f$Label))
p_alloc <- ggplot(f, aes(y = Label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
  geom_linerange(aes(xmin = Q05, xmax = Q95, color = Benchmark), linewidth = 0.5) +
  geom_linerange(aes(xmin = Q25, xmax = Q75, color = Benchmark), linewidth = 2.6) +
  geom_point(aes(x = Mean), shape = 21, fill = "white", color = "black", size = 1.8) +
  facet_grid(Study ~ Regime) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8")) +
  scale_color_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#E69F00"), guide = "none") +
  labs(x = "MSE relative to simple random sampling (log scale)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = -90)) +
  thesis_text_sizes() + theme(strip.text = element_text(face = "bold"))
save_fig(p_alloc, "compact/allocation_distribution_by_design.pdf", PRINT_W, 4.2)

# NC application exhibits: Figures B.10, B.11 ----
# Same plots as rs_exhibits() in application/code/real_sud_summary.R, from its CSVs.
ex_dir <- file.path(proj_dir, "application/results/real_sud_rev_20261002/exhibits")
ex_labels <- c("Checkerboard", "High Incidence Focus", "Saturation Quadrants", "Isolation Buffer",
               "2×2 Blocking", "Balanced Quartiles", "Balanced Halves",
               "Incidence­Guided Saturation", "Simple Random Sampling")
comparisons <- read.csv(file.path(ex_dir, "primary_setting_srs_comparison.csv"),
                        colClasses = c(Year = "character"))
comparisons$Design_Label <- factor(comparisons$Design_ID, levels = 1:9, labels = ex_labels)
comparisons$Parameter <- paste0("rho=", comparisons$Rho, ", gamma=", comparisons$Gamma)
comparisons$Regime_Label <- ifelse(comparisons$Regime == "both", "Both arms", "Control only")
stopifnot(nrow(comparisons) == 4 * 9 * 6 * 2, all(is.finite(comparisons$MSE_Ratio_SRS)))
p_ratios <- ggplot(comparisons, aes(x = Design_Label, y = Parameter, fill = log10(MSE_Ratio_SRS))) +
  geom_tile() + facet_grid(Year ~ Regime_Label) +
  scale_fill_gradient2(low = "#167c80", mid = "white", high = "#b95549", midpoint = 0,
                       breaks = log10(c(0.5, 1, 2, 4)), labels = c("0.5", "1", "2", "4"),
                       name = "MSE relative to\nsimple random sampling") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), legend.position = "bottom") +
  thesis_text_sizes()
save_fig(p_ratios, "real_sud/primary_setting_ratios.pdf", PRINT_W, 7.6)

tail <- read.csv(file.path(ex_dir, "refined_tail_comparison.csv"), colClasses = c(Year = "character"))
a <- tail[tail$Rho == 0.5, ]
a$Design_Label <- factor(a$Design_ID, levels = 9:1, labels = rev(ex_labels))
a$Regime_Label <- ifelse(a$Regime == "both", "Both arms", "Control only")
p_tail <- ggplot(a, aes(y = Design_Label)) +
  geom_segment(aes(x = Mean_MSE, xend = Q90_Estimated, yend = Design_Label), color = "#b9cdd1") +
  geom_point(aes(x = Mean_MSE, shape = "Mean"), color = "#167c80", size = 2) +
  geom_point(aes(x = Q90_Estimated, shape = "Estimated q90"), color = "#b77812", size = 2) +
  geom_point(aes(x = Worst10_Mean_Estimated, shape = "Estimated worst decile"), color = "#9f5362", size = 2) +
  facet_grid(Year ~ Regime_Label) + scale_x_log10() +
  scale_shape_manual(values = c(Mean = 16, "Estimated q90" = 17, "Estimated worst decile" = 15), name = NULL) +
  labs(x = "Conditional allocation MSE (log scale)", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom") +
  thesis_text_sizes()
save_fig(p_tail, "real_sud/refined_allocation_risk.pdf", PRINT_W, 7.4)

# Compact grid and NC exhibits: Figures 3.4 and 3.8 ----
# Same plots as export_compact_exhibits() in 20_manuscript_figure_revision.R, from the
# CSV records that function writes (results/manuscript_exhibit_revision_20261002/).
rec_dir <- file.path(results_dir, "manuscript_exhibit_revision_20261002")
short2 <- c("Incidence\u00adGuided Saturation", "Balanced Quartiles", "Isolation Buffer",
            "High Incidence Focus", "2\u00d72 Blocking", "Checkerboard", "Simple Random Sampling")
ids2 <- c(8, 6, 4, 2, 5, 1, 9)
d2 <- read.csv(file.path(rec_dir, "grid_configuration_intervals.csv"))
d2$Label <- factor(d2$Design, levels = paste("Design", rev(ids2)), labels = rev(short2))
d2$Regime <- factor(d2$Spillover_Type, levels = c("both", "control_only"),
                    labels = c("Spillover reaches both arms", "Spillover reaches controls only"))
d2$Configuration <- factor(paste(d2$Incidence_Mode, d2$Rho_Incidence, sep = "/"),
                           levels = c("iid/0", "spatial/0.2", "spatial/0.5", "poisson/0.2", "poisson/0.5"),
                           labels = c("Independent Uniform", "Spatial (0.2)", "Spatial (0.5)",
                                      "Poisson (0.2)", "Poisson (0.5)"))
stopifnot(nrow(d2) == 70, !anyNA(d2$Configuration), !anyNA(d2$Label))
compact_theme <- theme_minimal(base_size = 10) + theme(
  panel.grid.minor = element_blank(), legend.position = "bottom",
  legend.title = element_blank(), strip.text = element_text(face = "bold"),
  plot.margin = margin(3, 5, 3, 3))
cfg_colors <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#D55E00")
dodge <- position_dodge(width = 0.7)
p_grid <- ggplot(d2, aes(MSE, Label, color = Configuration, fill = Configuration,
                         shape = Configuration)) +
  geom_linerange(aes(xmin = MC_Lo, xmax = MC_Hi), position = dodge, linewidth = 0.45) +
  geom_point(position = dodge, size = 1.7, stroke = 0.3) +
  facet_wrap(~Regime, nrow = 1) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2)) +
  scale_color_manual(values = cfg_colors) + scale_fill_manual(values = cfg_colors) +
  scale_shape_manual(values = c(21, 22, 24, 23, 25)) +
  labs(x = "Mean MSE (log scale)", y = NULL) + compact_theme +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  thesis_text_sizes() + theme(strip.text = element_text(face = "bold"),
                              legend.title = element_blank())
save_fig(p_grid, "compact/grid_accuracy_srs.pdf", PRINT_W, 3.9)

annual <- read.csv(file.path(rec_dir, "annual_srs_ratios.csv"))
stopifnot(nrow(annual) == 72, all(table(annual$Year, annual$Regime) == 9))
nm <- c("Checkerboard", "High Incidence Focus", "Saturation Quadrants", "Isolation Buffer",
        "2\u00d72 Blocking", "Balanced Quartiles", "Balanced Halves",
        "Incidence\u00adGuided Saturation", "Simple Random Sampling")
annual$Label <- factor(annual$Design_ID, levels = 9:1, labels = rev(nm))
annual$Regime_Label <- factor(annual$Regime, levels = c("both", "control_only"),
                              labels = c("Spillover reaches both arms", "Spillover reaches controls only"))
annual$Ratio <- annual$Ratio_of_Mean_MSE_SRS
annual$Log_Ratio <- log2(annual$Ratio)
span <- max(abs(annual$Log_Ratio))
shown <- annual[annual$Design_ID != 9, ]   # SRS is the reference (always 1.00)
p_nc <- ggplot(shown, aes(factor(Year), Label, fill = Log_Ratio)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Ratio)), size = THESIS_TEXT_PT / ggplot2::.pt) +
  facet_wrap(~Regime_Label, nrow = 1) +
  scale_fill_gradient2(low = "#56B4E9", mid = "#F7F7F7", high = "#D55E00",
                       midpoint = 0, limits = c(-span, span), breaks = log2(c(0.25, 0.5, 1, 2, 4)),
                       labels = c("0.25", "0.5", "1", "2", "4"),
                       name = "MSE relative to\nsimple random sampling") +
  labs(x = "Year (observed SUD incidence)", y = NULL) + compact_theme +
  theme(panel.grid = element_blank()) +
  thesis_text_sizes() + theme(strip.text = element_text(face = "bold"))
save_fig(p_nc, "compact/nc_annual_accuracy.pdf", PRINT_W, 3.5)

cat("\nDone. Figures written to", fig_out, "\n")
