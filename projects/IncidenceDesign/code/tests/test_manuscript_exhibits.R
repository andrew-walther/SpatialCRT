# ============================================================
# Script: test_manuscript_exhibits.R
# Purpose: Verify compact exhibit aggregation, exact decomposition and matched ratios.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-02
# Dependencies: base R
# ============================================================
source("code/20_manuscript_figure_revision.R")
# With n=250, bias=0 and MSE=1 imply sample variance=250/249.
# A biased second configuration must remain separate, as must the regimes.
x <- data.frame(Neighbor_Type = c("queen", "queen", "queen", "rook"),
  True_Tau = 1, Estimator = "oracle", N_Valid_Est = 250,
  MSE = c(1, 5, 2, 99), Bias = c(0, 2, 1, 0),
  SD = sqrt(c(1, 1, 1, 99) * 250/249), Design = "Design 9",
  Incidence_Mode = c("iid", "spatial", "iid", "iid"),
  Rho_Incidence = c(0, 0.2, 0, 0), Spillover_Type = c("both", "both", "control_only", "both"))
got <- compact_grid_means(x)
stopifnot(nrow(got) == 3, setequal(got$MSE, c(1, 5, 2)),
          all(abs(got$Error_Variance - 1) < 1e-10),
          abs(got$Squared_Bias[got$Incidence_Mode == "spatial"] - 4) < 1e-10)
# Corrupting an MSE must fail, rather than produce a plausible wrong graph.
bad <- x; bad$MSE[1] <- 9
stopifnot(inherits(try(compact_grid_means(bad), silent = TRUE), "try-error"))
# Deliberately permuted rows and distinct denominators test named matching.
a <- data.frame(Year = c(2021, 2018, 2021, 2018), Regime = "both",
  Design_ID = c(8, 9, 9, 8), Mean_MSE = c(2, 10, 4, 20),
  Ratio_of_Mean_MSE_SRS = c(0.5, 1, 1, 2))
stopifnot(identical(annual_srs_ratios(a)$Ratio, c(0.5, 1, 1, 2)))
bad <- a; bad$Year[3] <- 2018
stopifnot(inherits(try(annual_srs_ratios(bad), silent = TRUE), "try-error"))
r <- readRDS("results/sim_data/sim_results_MLE_tau_sweep_combined_20260927_000502.rds")
m <- compact_grid_means(r)
stopifnot(nrow(m) == 90, max(abs(m$MSE - m$Squared_Bias - m$Error_Variance)) < 1e-10)
# Verify every exported mean directly against the original 16 scenario rows.
exported <- read.csv("results/manuscript_exhibit_revision_20261002/grid_configuration_means.csv")
for (i in seq_len(nrow(exported))) {
  e <- exported[i, ]
  rows <- r[r$Neighbor_Type == "queen" & r$True_Tau == 1 &
              r$Design == e$Design & r$Incidence_Mode == e$Incidence_Mode &
              r$Rho_Incidence == e$Rho_Incidence & r$Spillover_Type == e$Spillover_Type, ]
  stopifnot(nrow(rows) == 16, abs(e$MSE - mean(rows$MSE)) < 1e-10,
            abs(e$Squared_Bias - mean(rows$Bias^2)) < 1e-10,
            abs(e$Error_Variance - mean(rows$SD^2)*249/250) < 1e-10)
}
# Added manuscript details come from the long-form chapter, with fresh checks
# against current outputs rather than relying on rounded historical prose.
manuscript <- paste(readLines("paper/ctj_manuscript/CTJ_Manuscript.tex"), collapse = "\n")
for (tau in c(2, 3)) {
  values <- vapply(c("Design 4", "Design 8"), function(design) {
    mean(r$MSE[r$Neighbor_Type == "queen" & r$True_Tau == tau &
                 r$Spillover_Type == "control_only" & r$Design == design])
  }, numeric(1))
  stopifnot(all(vapply(sprintf("%.4f", values), grepl, logical(1),
                       x = manuscript, fixed = TRUE)))
}
nonoracle <- readRDS("results/sim_data/sim_results_MLEnonoracle_tau_sweep_combined_20260927_000502.rds")
nonoracle <- nonoracle[nonoracle$True_Tau == 1 &
  nonoracle$Design %in% paste("Design", c(1, 2, 4, 5, 6, 8)) &
  !(nonoracle$Neighbor_Type == "rook" & nonoracle$Design == "Design 1"), ]
pooled <- aggregate(cbind(Bias, Coverage) ~ Design + Neighbor_Type, nonoracle, mean)
stopifnot(identical(sprintf("%.2f", range(pooled$Bias)), c("-0.40", "-0.12")),
          identical(sprintf("%.2f", range(pooled$Coverage)), c("0.51", "0.90")))
decomposition <- read.csv("results/manuscript_exhibit_revision_20261002/grid_bias_variance.csv")
stopifnot(min(decomposition$Variance_Percent) > 97)
cat("PASS: exact error decomposition, configuration/regime separation, named annual SRS denominators, and wrong-value fixtures.\n")
cat("PASS: expanded manuscript effect-size crossover, omitted-spillover sensitivity and variance-share claims match current outputs.\n")
