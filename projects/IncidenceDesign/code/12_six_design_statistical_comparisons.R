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
# OUTPUTS (results/six_design_manuscript/):
#   - six_design_comparison_report.rds   (full list: friedman/nemenyi/wilcoxon/conditional)
#   - MLE_statistical_comparisons_6design.pdf (CD diagram, boxplot, heatmaps, conditional CDs)
#   - six_design_summary.txt             (console summary, saved to file)
# ==============================================================================

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
setwd(script_dir)
results_dir <- file.path(dirname(script_dir), "results")  # compute before sourcing 10, which overwrites script_dir
out_dir <- file.path(results_dir, "six_design_manuscript")
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

# ------------------------------------------------------------------------
# Load tau-sweep MLE results, relabel Design, filter to the 6 retained designs
# ------------------------------------------------------------------------
mle_full <- load_latest_results(results_dir = results_dir, estimation_mode = "MLE_tau_sweep")
stopifnot("True_Tau" %in% names(mle_full))

mle_full_6 <- mle_full[mle_full$Design %in% retained_ids, ]
mle_full_6$Design <- full_name_map[mle_full_6$Design]

mle_tau1_6 <- mle_full_6[mle_full_6$True_Tau == 1.0, ]

cat(sprintf("Full 6-design tau-sweep rows: %d (expect %d = 8*5*2*4*4*2*6/8)\n",
            nrow(mle_full_6), 5 * 2 * 4 * 4 * 2 * 6))
cat(sprintf("Primary tau=1.0 6-design rows: %d (expect %d)\n",
            nrow(mle_tau1_6), 2 * 4 * 4 * 2 * 6))

# ------------------------------------------------------------------------
# Primary aggregate tests at tau = 1.0 (mirrors 11_statistical_comparisons_report.qmd)
# ------------------------------------------------------------------------
sink(file.path(out_dir, "six_design_summary.txt"), split = TRUE)

cat("======================================================================\n")
cat("SIX-DESIGN STATISTICAL COMPARISONS (manuscript design set)\n")
cat("Retained: ", paste(unname(full_name_map[retained_ids]), collapse = ", "), "\n")
cat("Dropped:  Saturation Quadrants, Balanced Halves (redundant, see plan)\n")
cat("======================================================================\n\n")

cat("--- Primary scenario: True_Tau = 1.0 ---\n")
fr6  <- run_friedman_test(mle_tau1_6)
nem6 <- run_nemenyi_posthoc(mle_tau1_6)
wil6 <- run_pairwise_wilcoxon(mle_tau1_6)

cat(sprintf("Friedman chi-sq = %.2f, df = %d, p = %s, n_blocks = %d\n",
            fr6$statistic, fr6$n_designs - 1, format.pval(fr6$p_value, digits = 3),
            fr6$n_blocks))
cat("Average ranks (lower = better MSE):\n")
print(round(fr6$avg_ranks, 3))

cat(sprintf("\nNemenyi critical difference (alpha=0.05) = %.3f\n", nem6$critical_diff))
n_sig_nem <- sum(nem6$sig_matrix[upper.tri(nem6$sig_matrix)])
# NOTE: nem6$n_blocks is mislabeled upstream in run_nemenyi_posthoc() (10_statistical_comparisons.R) --
# it holds the scenario-block count (320), not the design count. Use n_designs directly for the
# pairwise-comparison denominator to avoid the "13 of 51040 pairs" nonsense this bug produces.
n_designs6 <- length(fr6$avg_ranks)
cat(sprintf("%d of %d pairs significantly different (Nemenyi)\n",
            n_sig_nem, choose(n_designs6, 2)))
cat("Nemenyi p-value matrix:\n")
print(round(nem6$p_matrix, 4))

n_sig_wil <- sum(wil6$sig_matrix[upper.tri(wil6$sig_matrix)])
cat(sprintf("\n%d of %d pairs significantly different (Wilcoxon, Holm-adjusted)\n",
            n_sig_wil, wil6$n_tests))

# Descriptive MSE/Coverage per design at tau=1.0 for manuscript table
desc6 <- aggregate(cbind(MSE, Coverage) ~ Design, data = mle_tau1_6, FUN = mean)
desc6 <- desc6[order(desc6$MSE), ]
cat("\nMean MSE / Coverage by design (tau=1.0, primary scenario):\n")
print(desc6, row.names = FALSE)

# ------------------------------------------------------------------------
# Conditional test across True_Tau (full sweep, 6-design set) -- confirms
# dominance holds across all tau levels for the reduced design set
# ------------------------------------------------------------------------
cat("\n--- Conditional test: stratify by True_Tau (full tau-sweep, 6 designs) ---\n")
cond_tau6 <- run_conditional_tests(mle_full_6, "True_Tau", test_type = "both")
cond_tau6_summary <- summarize_conditional_tests(cond_tau6)
for (nm in names(cond_tau6)) {
  entry <- cond_tau6[[nm]]
  if (!is.null(entry$friedman)) {
    cat(sprintf("  True_Tau = %s: Friedman chi-sq = %.2f, p = %s, n_blocks = %d\n",
                nm, entry$friedman$statistic,
                format.pval(entry$friedman$p_value, digits = 3),
                entry$friedman$n_blocks))
  }
}

sink()
cat("Summary written to", file.path(out_dir, "six_design_summary.txt"), "\n")

# ------------------------------------------------------------------------
# Figures with full design-name labels (no "D1"/"D3" shorthand) for
# manuscript exhibits. CD diagram keeps short_design_label() internally
# (space-constrained convention) but is supplementary-only.
# ------------------------------------------------------------------------
pdf(file.path(out_dir, "MLE_statistical_comparisons_6design.pdf"), width = 10, height = 6)
print(plot_cd_diagram(nem6, title = "Critical Difference Diagram (6-design set, tau=1.0)"))
print(plot_mse_boxplot_with_stars(mle_tau1_6, wil6))
print(plot_pvalue_heatmap(wil6, "Wilcoxon (Holm)"))
print(plot_pvalue_heatmap(nem6, "Nemenyi"))
cond_plots <- plot_conditional_cd_diagrams(cond_tau6, "True_Tau")
for (pl in cond_plots) print(pl)
dev.off()

saveRDS(
  list(friedman = fr6, nemenyi = nem6, wilcoxon = wil6,
       conditional_tau = cond_tau6, conditional_tau_summary = cond_tau6_summary,
       descriptive = desc6, retained_designs = unname(full_name_map[retained_ids])),
  file.path(out_dir, "six_design_comparison_report.rds")
)

# ------------------------------------------------------------------------
# Manuscript exhibits: MSE bar chart, coverage boxplot, tau-sensitivity
# line plot -- all use the Design column directly as axis/legend labels,
# so full design names appear automatically (no short_design_label() call
# involved, unlike the CD diagram above).
# ------------------------------------------------------------------------
p_mse     <- plot_master_comparison(mle_tau1_6, nb_filter = "queen", inc_label = "all configs, tau=1.0")
p_coverage <- plot_coverage_by_design(mle_tau1_6, inc_label = "all configs, tau=1.0")
p_tau     <- plot_mse_vs_tau(mle_full_6)

ggsave(file.path(out_dir, "fig_mse_by_design_6design.pdf"), p_mse, width = 10, height = 7)
ggsave(file.path(out_dir, "fig_coverage_by_design_6design.pdf"), p_coverage, width = 8, height = 6)
ggsave(file.path(out_dir, "fig_tau_sensitivity_6design.pdf"), p_tau, width = 8, height = 6)

cat("\nDone. Outputs in", out_dir, "\n")
