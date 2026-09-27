# ============================================================
# Script: verify_srs_reproduction.R
# Purpose: Check that adding the Simple Random Sampling benchmark (Design 9) left
#          every Design 1-8 result of the 2026-09-24 full run exactly unchanged.
# Author: Andrew Walther
# Created: 2026-09-26
# Dependencies: base R
# ============================================================
#
# Run from code/:  Rscript tests/verify_srs_reproduction.R
#
# Why this must hold: every draw is key-seeded (spec section 4). The Z key includes
# the design id d, and eps and X are keyed per block and per config, so adding d = 9
# changes no other design's assignments, noise or surfaces. Any difference in the
# Design 1-8 rows is a bug. Compares the scenario files and the surface files for both
# estimators with identical() after sorting by key. Writes
# results/srs_reproduction_check.txt; exits non-zero on any failure.

sim_dir <- file.path("..", "results", "sim_data")
old_stamp <- "20260924_025509"   # pre-SRS full run
latest <- function(pattern) {
  f <- list.files(sim_dir, pattern = pattern, full.names = TRUE)
  f <- f[!grepl(old_stamp, f)]
  stopifnot(length(f) > 0)
  f[which.max(file.mtime(f))]
}

#' Sort rows by every non-numeric-metric key column so row order can't matter
#'
#' @param df Data frame (scenario or surface results)
#' @param keys Character vector of key columns present in df
#' @return df sorted by keys, row names reset
sort_by_keys <- function(df, keys) {
  df <- df[do.call(order, unname(as.list(df[keys]))), , drop = FALSE]
  rownames(df) <- NULL
  df
}

pairs <- list(
  scen_oracle    = c("^sim_results_MLE_tau_sweep_combined_",
                     sprintf("sim_results_MLE_tau_sweep_combined_%s.rds", old_stamp)),
  scen_nonoracle = c("^sim_results_MLEnonoracle_tau_sweep_combined_",
                     sprintf("sim_results_MLEnonoracle_tau_sweep_combined_%s.rds", old_stamp)),
  surf_oracle    = c("^surface_results_MLE_tau_sweep_",
                     sprintf("surface_results_MLE_tau_sweep_%s.rds", old_stamp)),
  surf_nonoracle = c("^surface_results_MLEnonoracle_tau_sweep_",
                     sprintf("surface_results_MLEnonoracle_tau_sweep_%s.rds", old_stamp))
)

out <- c("SRS REPRODUCTION CHECK: Designs 1-8, new run vs 2026-09-24 run", "")
all_ok <- TRUE
for (nm in names(pairs)) {
  new <- readRDS(latest(pairs[[nm]][1]))
  old <- readRDS(file.path(sim_dir, pairs[[nm]][2]))
  new18 <- new[new$Design != "Design 9", , drop = FALSE]
  keys <- intersect(c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho",
                      "Gamma", "Spillover_Type", "True_Tau", "Estimator", "Surface"),
                    names(old))
  same_cols <- identical(sort(names(new18)), sort(names(old)))
  a <- sort_by_keys(new18[, names(old)], keys)
  b <- sort_by_keys(old, keys)
  same <- same_cols && identical(a, b)
  n9 <- sum(new$Design == "Design 9")
  line <- sprintf("%-15s %s  old rows %d | new rows %d (Design 1-8: %d, Design 9: %d)",
                  nm, if (same) "PASS" else "FAIL", nrow(old), nrow(new), nrow(new18), n9)
  if (!same) {
    all_ok <- FALSE
    diffs <- if (same_cols && nrow(a) == nrow(b))
      names(old)[!vapply(names(old), function(v) identical(a[[v]], b[[v]]), logical(1))]
      else "column set or row count differs"
    line <- paste0(line, "\n    differing: ", paste(diffs, collapse = ", "))
  }
  out <- c(out, line)
}
out <- c(out, "", if (all_ok) "ALL IDENTICAL" else "MISMATCH: investigate before using the new run")
writeLines(out, file.path("..", "results", "srs_reproduction_check.txt"))
cat(out, sep = "\n")
if (!all_ok) quit(status = 1)
