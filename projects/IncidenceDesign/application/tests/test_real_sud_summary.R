# ============================================================
# Script: test_real_sud_summary.R
# Purpose: Verify matched SRS contrasts and their Monte Carlo uncertainty.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud summary helpers
# ============================================================
test_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
real_sud_summary_define_only <- TRUE
source(file.path(test_dir, "..", "code", "real_sud_summary.R"))
p <- data.frame(Year = rep(c("2018", "2019"), each = 2), Model = "education",
  Neighbor = "queen", Summary = "mean_rank", Rho = 0, Gamma = 0.5, Regime = "both",
  Design_ID = rep(c(3, 9), 2), Mean_MSE = c(2, 4, 3, 6),
  SE_Mean_MSE_Joint = c(0.2, 0.3, 0.4, 0.5), Bias = 0, Coverage = 0.93,
  Source_ID = letters[1:4])
out <- rs_srs_comparison(p[c(4, 1, 2, 3), ])
stopifnot(identical(out$Source_ID, p$Source_ID[c(4, 1, 2, 3)]),
  all(out$MSE_Ratio_SRS[out$Design_ID == 3] == 0.5),
  all(out$MSE_Difference_SRS[out$Design_ID == 9] == 0),
  all(out$MCSE_Difference_SRS[out$Design_ID == 9] == 0),
  abs(out$MCSE_Difference_SRS[out$Source_ID == "a"]^2 - (0.2^2 + 0.3^2)) < 1e-12)
stopifnot(inherits(tryCatch(rs_srs_comparison(p[-4, ]), error = identity), "error"))
rook <- p; rook$Neighbor <- "rook"; rook$Mean_MSE <- rook$Mean_MSE * 2
rook$Source_ID <- paste0("rook-", rook$Source_ID)
matched <- rs_primary_reference(rbind(p, rook))
stopifnot(all(matched$MSE_Ratio_Matched_Primary[matched$Neighbor == "rook"] == 2),
  all(matched$MSE_Ratio_Matched_Primary[matched$Neighbor == "queen"] == 1))
stopifnot(inherits(tryCatch(rs_primary_reference(rbind(p[-4, ], rook)), error = identity), "error"))
cat("PASS: year-matched SRS ratios, independent-design MC contrasts, self-reference covariance and missing-reference rejection.\n")
cat("PASS: sensitivities match primary parameter settings and reject missing references.\n")
