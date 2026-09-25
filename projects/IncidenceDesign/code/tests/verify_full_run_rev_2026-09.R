# ============================================================
# Script: verify_full_run_rev_2026-09.R
# Purpose: B5 integrity checks on the full 2026-09 run (plan B5 "verify" line).
# Author: Andrew Walther
# Created: 2026-09-24
# Dependencies: base R
# ============================================================
#
# Run from code/:  Rscript tests/verify_full_run_rev_2026-09.R
# Checks: 12,800 rows per estimator; N_Valid_Est = 250; aliasing and rank-deficiency
# flags only for Checkerboard x rook; lists every non-aliasing warning. Writes
# results/sim_data/full_run_verification.txt; exits non-zero on failure.

sim_dir <- file.path("..", "results", "sim_data")
latest <- function(pattern) {
  f <- list.files(sim_dir, pattern = pattern, full.names = TRUE)
  stopifnot(length(f) > 0)
  f[which.max(file.mtime(f))]
}
files <- c(oracle = latest("^sim_results_MLE_tau_sweep_combined_.*\\.rds$"),
           nonoracle = latest("^sim_results_MLEnonoracle_tau_sweep_combined_.*\\.rds$"))
warn_files <- c(oracle = latest("^warnings_MLE_tau_sweep_.*\\.csv$"),
                nonoracle = latest("^warnings_MLEnonoracle_tau_sweep_.*\\.csv$"))

out <- c("FULL-RUN VERIFICATION (2026-09 revision)", "")
all_ok <- TRUE
for (est in names(files)) {
  r <- readRDS(files[[est]])
  cb_rook <- r$Design == "Design 1" & r$Neighbor_Type == "rook"
  chk <- c(
    "12,800 scenario rows" = nrow(r) == 12800,
    "one row per scenario key" = !anyDuplicated(r[, c("Incidence_Mode", "Rho_Incidence",
      "Neighbor_Type", "Design", "Rho", "Gamma", "Spillover_Type", "True_Tau")]),
    "Estimator column matches file" = all(r$Estimator == est),
    "N_Valid_Est = 250 everywhere" = all(r$N_Valid_Est == 250),
    "Fail_Rate = 0 everywhere" = all(r$Fail_Rate == 0),
    "N_Surfaces = 10 everywhere" = all(r$N_Surfaces == 10),
    # Rule (spec section 7): an oracle fit is aliased iff its Z draw has
    # rank([1, Z, WZ]) < 3. Checkerboard x rook is always rank-deficient. Isolation
    # Buffer x rook can rarely draw the exact checkerboard (a maximal independent
    # set equal to one colour class); no other design/nb can.
    "Checkerboard x rook: all 250 fits flagged" =
      all(r$Z_WZ_rank_deficient[cb_rook]) && (est != "oracle" || all(r$N_Aliased[cb_rook] == 250)),
    "rank deficiency only Checkerboard or Isolation Buffer x rook" =
      !any(r$Z_WZ_rank_deficient & !(r$Neighbor_Type == "rook" & r$Design %in% c("Design 1", "Design 4"))),
    "oracle: N_Aliased > 0 iff Z_WZ_rank_deficient; non-oracle never aliased" = if (est == "oracle")
      identical(r$N_Aliased > 0, r$Z_WZ_rank_deficient) else all(r$N_Aliased == 0),
    "fixed-size designs treat 50" = all(r$Mean_Treated[!r$Design %in% "Design 4"] == 50),
    "no NA in core metrics" = !anyNA(r[, c("Bias", "SD", "MSE", "Coverage", "Power", "SE_MSE")])
  )
  all_ok <- all_ok && all(chk)
  ib <- r$Z_WZ_rank_deficient & r$Design == "Design 4"
  # With no warnings at all, 05 writes a zero-column CSV whose only line is ""
  wl <- readLines(warn_files[[est]], warn = FALSE)
  w <- if (all(wl %in% c("", "\"\""))) data.frame() else read.csv(warn_files[[est]])
  out <- c(out, sprintf("== %s: %s", est, basename(files[[est]])),
           sprintf("  %-66s %s", names(chk), ifelse(chk, "PASS", "FAIL")),
           sprintf("  Isolation Buffer x rook scenarios with a checkerboard draw: %d (aliased fits: %d)",
                   sum(ib), sum(r$N_Aliased[ib])),
           sprintf("  non-aliasing warnings: %d scenario-message rows, %d fits in total",
                   nrow(w), if (nrow(w)) sum(w$Count) else 0L))
  if (nrow(w)) {
    agg <- aggregate(Count ~ Message, data = w, FUN = sum)
    out <- c(out, sprintf("    [%d fits] %s", agg$Count, agg$Message))
  }
  out <- c(out, "")
}
out <- c(out, sprintf("RESULT: %s", if (all_ok) "PASS" else "FAIL"))
writeLines(out, file.path(sim_dir, "full_run_verification.txt"))
cat(out, sep = "\n")
if (!all_ok) stop("Full-run verification FAILED")
