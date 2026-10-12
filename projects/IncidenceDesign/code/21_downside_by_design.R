# ============================================================
# Script: 21_downside_by_design.R
# Purpose: Compare every design with SRS on mean allocation MSE and on
#          allocation downside (worst-10% mean), by spillover regime, for the
#          grid pilot and the refined NC application corners.
# Author: Andrew Walther (drafted with Claude Code)
# Created: 2026-10-11
# Dependencies: base R only
# ============================================================
# Run from the project root: Rscript code/21_downside_by_design.R
# Reads completed aggregate outputs only; it neither simulates nor changes
# checkpoints. Output: results/allocation_risk/summary/downside_by_design.csv
#
# For allocation a, m(a) = E_eps[(tau_hat - tau)^2 | a, X] is its conditional
# MSE. Each source row already summarizes the sampled allocations of one
# design: their mean m(a) and the mean of the worst 10% of m(a).
#   Grid pilot: queen, tau = 1, oracle fit, averaged over 40 fixed
#     surface/corner units (as reported in the chapter).
#   NC: refined queen corners rho in {0, 0.5}, gamma = 0.8, education model,
#     mean-rank regions; averaged with equal weight over 2018-2021 and both
#     corners (8 rows per design and regime).
# Ratios are (design average) / (SRS average) within the same source and regime,
# so a value below 1 means lower error than SRS.

# Read sources ----
grid <- read.csv("results/allocation_risk/summary/pilot/pooled_fixed_units.csv",
                 stringsAsFactors = FALSE)
nc <- read.csv(paste0("application/results/real_sud_rev_20261002/exhibits/",
                      "refined_tail_comparison.csv"), stringsAsFactors = FALSE)

grid <- grid[grid$Neighbor_Type == "queen" & grid$Estimator == "oracle" &
               grid$True_Tau == 1, ]
grid$Design_ID <- as.integer(sub("Design ", "", grid$Design))
nc <- nc[nc$Neighbor == "queen" & nc$Model == "education" &
           nc$Summary == "mean_rank" & nc$Gamma == 0.8, ]

# Fail loud if the sources are not the expected shape ----
stopifnot(nrow(grid) == 18, setequal(grid$Design_ID, 1:9),
          all(grid$N_Fixed_Units == 40))
nc_counts <- table(nc$Design_ID, nc$Regime)
stopifnot(nrow(nc_counts) == 9, all(nc_counts == 8))

# Summarize: average mean MSE and worst-10% mean, then ratio to SRS (ID 9) ----
summarize_source <- function(ids, regimes, mean_mse, worst10, source) {
  out <- aggregate(cbind(Mean_MSE = mean_mse, Worst10 = worst10) ~
                     Design_ID + Regime,
                   data = data.frame(Design_ID = ids, Regime = regimes,
                                     mean_mse = mean_mse, worst10 = worst10),
                   FUN = mean)
  srs <- out[out$Design_ID == 9, c("Regime", "Mean_MSE", "Worst10")]
  names(srs) <- c("Regime", "SRS_Mean_MSE", "SRS_Worst10")
  out <- merge(out, srs, by = "Regime")
  out$Mean_Ratio_SRS <- out$Mean_MSE / out$SRS_Mean_MSE
  out$Worst10_Ratio_SRS <- out$Worst10 / out$SRS_Worst10
  out$Source <- source
  out
}

result <- rbind(
  summarize_source(grid$Design_ID, grid$Spillover_Type,
                   grid$Mean_Allocation_MSE, grid$Average_Conditional_Worst10,
                   "grid_pilot"),
  summarize_source(nc$Design_ID, nc$Regime, nc$Mean_MSE,
                   nc$Worst10_Mean_Estimated, "nc_refined")
)
result <- result[order(result$Source, result$Regime, result$Design_ID),
                 c("Source", "Regime", "Design_ID", "Mean_MSE", "Worst10",
                   "Mean_Ratio_SRS", "Worst10_Ratio_SRS")]

# Write ----
out_file <- "results/allocation_risk/summary/downside_by_design.csv"
write.csv(result, out_file, row.names = FALSE)
cat("Wrote", nrow(result), "rows to", out_file, "\n")
