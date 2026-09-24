# ==============================================================================
# 12_six_design_statistical_comparisons.R
# Re-run formal statistical comparisons (Friedman/Nemenyi/Wilcoxon) restricted
# to the 6 designs retained for the manuscript, dropping Saturation Quadrants
# (redundant with Incidence-Guided Saturation Quadrants) and Balanced Halves
# (redundant with Balanced Quartiles). The full 8-design comparison in
# results/11_statistical_comparisons_report.{html,pdf} and
# results/MLE_statistical_comparisons.pdf is left untouched as the
# supplementary-material version.
#
# Mirrors the primary-tau / full-sweep split used in
# 11_statistical_comparisons_report.qmd: aggregate Friedman/Nemenyi/Wilcoxon
# run on the primary tau=1.0 slice; the True_Tau conditional stratification
# uses the full tau-sweep set so dominance-across-tau can still be checked
# for the reduced design set.
#
# 2026-09 revision: neighbor types are analyzed separately -- queen (primary) and
# rook (Chapter 2 continuity / sensitivity; tau is not identified for Checkerboard
# under rook, whose oracle estimates are kept but flagged) -- plus a POOLED slice
# that is always labeled as such. The estimator is an argument, so the non-oracle
# sensitivity results go to their own directory.
#
# Usage (from code/):
#   Rscript 12_six_design_statistical_comparisons.R                 # oracle
#   Rscript 12_six_design_statistical_comparisons.R MLEnonoracle    # non-oracle
#
# OUTPUTS (results/six_design_manuscript/ for oracle,
#          results/six_design_manuscript_nonoracle/ for non-oracle):
#   - six_design_comparison_report.rds  list(queen, rook, pooled), each with
#       friedman/nemenyi/wilcoxon/conditional_tau/descriptive
#   - MLE_statistical_comparisons_6design_{queen,rook}.pdf
#   - six_design_summary.txt            (console summary, saved to file)
#   - fig_{mse_by_design,coverage_by_design,tau_sensitivity}_6design_{queen,rook}.pdf
# ==============================================================================

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
setwd(script_dir)
results_dir <- file.path(dirname(script_dir), "results")  # compute before sourcing 10, which overwrites script_dir

# Estimator: "MLE" (oracle, primary) or "MLEnonoracle" (sensitivity, M8)
if (!exists("est_prefix")) {
  trailing <- commandArgs(trailingOnly = TRUE)
  est_prefix <- if (length(trailing) > 0) trailing[1] else "MLE"
}
stopifnot(est_prefix %in% c("MLE", "MLEnonoracle"))
est_tag <- paste0(est_prefix, "_tau_sweep")
out_dir <- file.path(results_dir, if (est_prefix == "MLE") "six_design_manuscript" else
  "six_design_manuscript_nonoracle")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("10_statistical_comparisons.R")  # sources 06_visualizations.R -> 01-03

# ------------------------------------------------------------------------
# Map raw "Design N" values to full descriptive names (no "Design N:" prefix,
# per manuscript no-shorthand rule). Source of truth: get_design_names().
# ------------------------------------------------------------------------
raw_names <- get_design_names()  # e.g. "Design 1: Checkerboard"
full_name_map <- setNames(
  sub("^Design [0-9]+: ", "", raw_names),
  paste0("Design ", seq_along(raw_names))
)
stopifnot(identical(unname(full_name_map),
                     c("Checkerboard", "High Incidence Focus", "Saturation Quadrants",
                       "Isolation Buffer", "2x2 Blocking", "Balanced Quartiles",
                       "Balanced Halves", "Incidence-Guided Saturation Quadrants")))

# 6 retained designs: drop "Design 3" (Saturation Quadrants) and "Design 7" (Balanced Halves)
retained_ids <- c("Design 1", "Design 2", "Design 4", "Design 5", "Design 6", "Design 8")

ROOK_NOTE <- "Rook: tau not identified for Checkerboard (WZ = 1 - Z); its estimates are kept but flagged."

# ------------------------------------------------------------------------
# Load tau-sweep results for this estimator, relabel Design, keep 6 designs
# ------------------------------------------------------------------------
mle_full <- load_latest_results(results_dir = results_dir, estimation_mode = est_tag)
stopifnot("True_Tau" %in% names(mle_full))

mle_full_6 <- mle_full[mle_full$Design %in% retained_ids, ]
mle_full_6$Design <- full_name_map[mle_full_6$Design]

# Integrity: every block (all parameters except Design) has all 6 designs, once
block_cols <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Rho", "Gamma",
                "Spillover_Type", "True_Tau")
per_block <- table(do.call(paste, mle_full_6[, block_cols]))
stopifnot(all(per_block == 6))
cat(sprintf("6-design rows: %d (%d blocks x 6 designs; %d tau levels)\n",
            nrow(mle_full_6), length(per_block), length(unique(mle_full_6$True_Tau))))

slices <- list(
  queen  = mle_full_6[mle_full_6$Neighbor_Type == "queen", ],
  rook   = mle_full_6[mle_full_6$Neighbor_Type == "rook", ],
  pooled = mle_full_6)
slice_title <- c(queen = "QUEEN (primary)", rook = "ROOK (sensitivity / Chapter 2 continuity)",
                 pooled = "POOLED over rook and queen (label as pooled wherever used)")

# ------------------------------------------------------------------------
# Tests per slice: aggregate at tau = 1.0; conditional across tau (full sweep)
# ------------------------------------------------------------------------
#' Friedman / Nemenyi / Wilcoxon (tau = 1) and the True_Tau-conditional test for one slice
#'
#' @param res Data frame of 6-design results for one neighbor slice (all tau)
#' @return List with friedman, nemenyi, wilcoxon, conditional_tau,
#'   conditional_tau_summary, descriptive, n_rows_tau1
run_slice <- function(res) {
  tau1 <- res[res$True_Tau == 1.0, ]
  fr  <- run_friedman_test(tau1)
  nem <- run_nemenyi_posthoc(tau1)
  wil <- run_pairwise_wilcoxon(tau1)
  desc <- aggregate(cbind(MSE, Coverage) ~ Design, data = tau1, FUN = mean)
  desc <- desc[order(desc$MSE), ]
  cond <- run_conditional_tests(res, "True_Tau", test_type = "both")
  # run_conditional_tests() turns errors into NULL entries; never drop a tau silently
  stopifnot(length(cond) == length(unique(res$True_Tau)),
            all(vapply(cond, function(e) !is.null(e$friedman), logical(1))))
  list(friedman = fr, nemenyi = nem, wilcoxon = wil, conditional_tau = cond,
       conditional_tau_summary = summarize_conditional_tests(cond),
       descriptive = desc, n_rows_tau1 = nrow(tau1))
}
report <- lapply(slices, run_slice)

sink(file.path(out_dir, "six_design_summary.txt"), split = TRUE)
cat("======================================================================\n")
cat("SIX-DESIGN STATISTICAL COMPARISONS (manuscript design set)\n")
cat("Estimator: ", if (est_prefix == "MLE") "oracle (Y ~ Z + Spill + X), primary" else
      "non-oracle (Y ~ Z + X), sensitivity", "\n")
cat("Retained: ", paste(unname(full_name_map[retained_ids]), collapse = ", "), "\n")
cat("Dropped:  Saturation Quadrants, Balanced Halves (redundant, see plan)\n")
cat("Per-tau Friedman tests share Z and eps across tau (common random numbers):\n")
cat("they are not independent confirmations of one another.\n")
cat("======================================================================\n")

for (sl in names(report)) {
  r <- report[[sl]]
  cat(sprintf("\n######## %s ########\n", slice_title[[sl]]))
  if (sl != "queen") cat(ROOK_NOTE, "\n")
  cat(sprintf("--- True_Tau = 1.0 (%d rows, %d blocks) ---\n", r$n_rows_tau1, r$friedman$n_blocks))
  cat(sprintf("Friedman chi-sq = %.2f, df = %d, p = %s\n", r$friedman$statistic,
              r$friedman$n_designs - 1, format.pval(r$friedman$p_value, digits = 3)))
  cat("Average ranks (lower = better MSE):\n")
  print(round(r$friedman$avg_ranks, 3))
  # choose(k, 2) pair denominator: nemenyi$n_blocks is the block count, not designs
  n_designs6 <- length(r$friedman$avg_ranks)
  cat(sprintf("\nNemenyi critical difference (alpha=0.05) = %.3f\n", r$nemenyi$critical_diff))
  cat(sprintf("%d of %d pairs significantly different (Nemenyi)\n",
              sum(r$nemenyi$sig_matrix[upper.tri(r$nemenyi$sig_matrix)]), choose(n_designs6, 2)))
  print(round(r$nemenyi$p_matrix, 4))
  cat(sprintf("\n%d of %d pairs significantly different (Wilcoxon, Holm-adjusted)\n",
              sum(r$wilcoxon$sig_matrix[upper.tri(r$wilcoxon$sig_matrix)]), r$wilcoxon$n_tests))
  cat("\nMean MSE / Coverage by design (tau=1.0):\n")
  print(r$descriptive, row.names = FALSE)
  cat("\n--- Conditional test: stratify by True_Tau (full sweep) ---\n")
  for (nm in names(r$conditional_tau)) {
    entry <- r$conditional_tau[[nm]]
    if (!is.null(entry$friedman)) {
      cat(sprintf("  %s: Friedman chi-sq = %.2f, p = %s, n_blocks = %d; rank order: %s\n",
                  nm, entry$friedman$statistic, format.pval(entry$friedman$p_value, digits = 3),
                  entry$friedman$n_blocks, paste(names(entry$friedman$avg_ranks), collapse = " < ")))
    }
  }
}
sink()
cat("Summary written to", file.path(out_dir, "six_design_summary.txt"), "\n")

saveRDS(c(report, list(estimator = est_prefix, retained_designs = unname(full_name_map[retained_ids]),
                       rook_note = ROOK_NOTE)),
        file.path(out_dir, "six_design_comparison_report.rds"))

# ------------------------------------------------------------------------
# Figures per neighbor type (queen primary, rook separate; no pooled figures).
# CD diagram keeps short_design_label() internally (supplementary-only); the
# manuscript exhibits use full design names from the Design column.
# ------------------------------------------------------------------------
for (nb in c("queen", "rook")) {
  r <- report[[nb]]
  res <- slices[[nb]]
  tau1 <- res[res$True_Tau == 1.0, ]
  sub_note <- if (nb == "rook") paste0(" [", ROOK_NOTE, "]") else ""
  pdf(file.path(out_dir, sprintf("MLE_statistical_comparisons_6design_%s.pdf", nb)),
      width = 10, height = 6)
  print(plot_cd_diagram(r$nemenyi, title = sprintf("Critical Difference Diagram (6 designs, %s, tau=1.0)", nb)))
  print(plot_mse_boxplot_with_stars(tau1, r$wilcoxon))
  print(plot_pvalue_heatmap(r$wilcoxon, "Wilcoxon (Holm)"))
  print(plot_pvalue_heatmap(r$nemenyi, "Nemenyi"))
  for (pl in plot_conditional_cd_diagrams(r$conditional_tau, "True_Tau")) print(pl)
  dev.off()

  p_mse <- plot_master_comparison(tau1, nb_filter = nb,
                                  inc_label = paste0("all configs, tau=1.0", sub_note))
  # Chapter figures (coverage, tau sensitivity): Greek tau via plotmath
  # expressions, saved with the default pdf() device (Unicode Greek fails in
  # pdf(), and cairo_pdf cannot load on this machine without XQuartz). The
  # figure note is a Greek-letter copy of ROOK_NOTE; ROOK_NOTE itself stays
  # unchanged because it is also written to the summary .txt and .rds. It goes
  # in a left-aligned plot-wide caption so it is never clipped.
  fig_rook_note <- if (nb == "rook")
    expression("Rook: " * tau * " not identified for Checkerboard (WZ = 1 - Z); its estimates are kept but flagged.") else NULL
  note_theme <- theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0))
  p_coverage <- plot_coverage_by_design(tau1) +
    labs(subtitle = expression("Incidence: all configs, " * tau * " = 1.0 | Red dashed = nominal 95%"),
         caption = fig_rook_note) +
    note_theme
  # Legibility only: coverage is tightly clustered, so the default black outlines
  # collapse the boxes into dark bars; thinner grey outlines let the fill show.
  p_coverage$layers[[1]] <- geom_boxplot(alpha = 0.5, outlier.alpha = 0.5, outlier.size = 1,
                                         colour = "grey35", linewidth = 0.25)
  # Legend order = queen MSE at tau = 1 (best to worst), same in both versions.
  # The ribbon is Mean_MSE +/- the AVERAGE of the stored scenario-level SE_MSE
  # (surface-level MC SEs) in each Design x tau cell, not the SE of the mean.
  tau_legend_order <- c("Incidence-Guided Saturation Quadrants", "Balanced Quartiles",
                        "Isolation Buffer", "High Incidence Focus", "2x2 Blocking", "Checkerboard")
  stopifnot(setequal(tau_legend_order, unique(res$Design)))
  res_tau <- res
  res_tau$Design <- factor(res_tau$Design, levels = tau_legend_order)
  p_tau <- plot_mse_vs_tau(res_tau) +
    labs(title = expression("MSE vs. True " * tau * " by Design"), x = expression("True " * tau),
         y = "Mean MSE\n(band: ± average scenario-level Monte Carlo SE)",
         caption = fig_rook_note) +
    note_theme
  ggsave(file.path(out_dir, sprintf("fig_mse_by_design_6design_%s.pdf", nb)), p_mse, width = 10, height = 7)
  ggsave(file.path(out_dir, sprintf("fig_coverage_by_design_6design_%s.pdf", nb)), p_coverage, width = 8, height = 6)
  ggsave(file.path(out_dir, sprintf("fig_tau_sensitivity_6design_%s.pdf", nb)), p_tau, width = 8, height = 6)
}

cat("\nDone. Outputs in", out_dir, "\n")
