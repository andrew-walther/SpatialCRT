# ============================================================
# Script: methods_ext_bayes_long_tdf.R
# Purpose: Maximum likelihood coverage in the long-chain Bayesian comparison
#   (a6bayeslong) under the corrected t interval reported throughout Chapter 3,
#   beside the Wald coverage it replaces in the text.
# Author: Andrew Walther
# Created: 2026-10-30
# Dependencies: methods_ext_sensitivity_summary.R (mz_setting)
# ============================================================
# Run from the project root:  Rscript application/code/methods_ext_bayes_long_tdf.R
# Uses the stored per-fit estimates and SEs; no refitting. Writes
# application/results/methods_ext_2026-10/bayes_long_tdf.log.

.mb_dir <- "application/code"
suppressMessages(source(file.path(.mb_dir, "methods_ext_sensitivity_summary.R")))

# Interval ----
# tau-hat +/- t(0.975, n - p) x SE x sqrt(n / (n - p)), n = 58 colleges,
# p = 3 coefficients (1, Z, Spill) + rho = 4, so df = 54 (as in A5 Stage 2).
f <- readRDS(file.path(.me_root, "a6bayeslong", "fits.rds"))
k <- 3; tc <- stats::qt(0.975, 58 - k - 1); fb <- sqrt(58 / (58 - k - 1))
f$T_lo <- f$SAR_tau - tc * fb * f$SAR_se
f$T_hi <- f$SAR_tau + tc * fb * f$SAR_se

# Coverage by regime, averaged over settings (as methods_ext_bayes_long_summary.R) ----
by <- c("Year", "Rho", "Regime", "Design_ID")
w <- aggregate(Coverage ~ Regime, mz_setting(f, "SAR", by), mean)
t <- aggregate(Coverage ~ Regime, mz_setting(f, "SAR", by, "T_lo", "T_hi"), mean)
out <- sprintf("rho = 0.5, SAR ML, %s: Wald coverage %.3f; corrected t (df 54) coverage %.3f",
               w$Regime, w$Coverage, t$Coverage[match(w$Regime, t$Regime)])
writeLines(out); writeLines(out, file.path(.me_root, "bayes_long_tdf.log"))
