# ============================================================
# Script: 19_manuscript_exhibits.R
# Purpose: Export the current grid benchmark and a legible CTJ grid figure.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing visualization modules (ggplot2, dplyr, tidyr, viridis)
# ============================================================
# Run from the project root: Rscript code/19_manuscript_exhibits.R
# This reads completed results; it neither simulates nor changes checkpoints.

# Read verified current outputs ----
script_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- dirname(dirname(normalizePath(script_file)))
source(file.path(project_dir, "code", "06_visualizations.R"))
source_file <- file.path(project_dir, "results", "sim_data",
                         "sim_results_MLE_tau_sweep_combined_20260927_000502.rds")
r <- readRDS(source_file)
stopifnot(nrow(r) == 14400, all(r$Estimator == "oracle"),
          all(r$N_Valid_Est == 250), length(unique(r$Design)) == 9)
a <- r[r$True_Tau == 1 & r$Neighbor_Type == "queen", ]
benchmark <- aggregate(cbind(MSE, Bias, Coverage) ~ Design + Spillover_Type,
                       a, mean)
stopifnot(nrow(benchmark) == 18,
          all(table(a$Design, a$Spillover_Type) == 80))
write.csv(benchmark, file.path(project_dir, "results", "srs_benchmark",
                               "manuscript_regime_means.csv"), row.names = FALSE)

# Journal figure: retain the chapter's detailed six-design scope ----
ids <- c(8, 6, 4, 2, 5, 1)
labels <- c("Guided saturation", "Balanced quartiles", "Isolation buffer",
            "High incidence focus", "Spatial blocking", "Checkerboard")
d <- a[a$Design %in% paste("Design", ids), ]
d$Design <- factor(d$Design, levels = paste("Design", ids), labels = labels)
# Five incidence configurations overlay within a design/gamma group. Equal-height
# errorbar endpoints mark each value, preserving the chapter's bar interpretation.
p <- plot_master_comparison(d, "queen") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = MSE, ymax = MSE),
                        position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7, linewidth = 0.2) +
  ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL, fill = expression(gamma)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 9),
                 legend.position = "bottom", strip.text = ggplot2::element_text(size = 11))
stopifnot(nrow(d) == 960, all(table(d$Design) == 160))
ggplot2::ggsave(file.path(project_dir, "paper", "ctj_manuscript", "figures",
                         "fig_mse_by_design_6design_queen_journal.pdf"),
                p, width = 10.5, height = 7)
writeLines(c(paste("Source:", basename(source_file)),
             paste("Source MD5:", unname(tools::md5sum(source_file))),
             "Benchmark: tau=1, queen, 80 scenarios/design/regime; 18 rows.",
             "Figure: 960 current scenario values, detailed six-design subset.",
             "No simulation outputs or checkpoints changed."),
           file.path(project_dir, "results", "srs_benchmark", "manuscript_exhibit_sources.txt"))
cat("PASS: current nine-design benchmark and six-design journal exhibit exported.\n")
