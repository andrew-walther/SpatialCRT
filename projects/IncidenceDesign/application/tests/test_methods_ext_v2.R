# ============================================================
# Script: test_methods_ext_v2.R
# Purpose: Checks for the 2026-10-21 changes: the >= 2-college saturation-region
#   wrapper (V6 fix) and the t + df-corrected test (coverage decision).
# Author: Andrew Walther
# Created: 2026-10-21
# Dependencies: methods_ext_participation.R (sources methods_ext_core.R)
# ============================================================
# Run from the project root (needs setup.rds, i.e. the restricted data):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript application/tests/test_methods_ext_v2.R
# Writes application/results/methods_ext_2026-10/verification/test_methods_ext_v2.txt;
# exit status 1 on any failure.

test_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
methods_ext_define_only <- TRUE
source(file.path(test_dir, "..", "code", "methods_ext_participation.R"))
log_lines <- character(0); failed <- character(0)
check <- function(name, ok, detail = "") {
  line <- sprintf("%-78s %s %s", name, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(line, "\n"); log_lines <<- c(log_lines, line); if (!isTRUE(ok)) failed <<- c(failed, name)
}
setup <- me_load_setup()

# Region wrapper versus production make_kmeans_regions() on the same subsets ----
# Same seed and candidates: whenever production's best partition has all regions
# >= 2, the wrapper must return exactly that partition; otherwise a different one
# with all regions >= 2, flagged changed.
res <- list()
for (q in 1:400) {
  N <- c(15, 20, 25, 29)[(q - 1) %% 4 + 1]; yr <- as.character(2018 + (q %% 4))
  P <- me_subset("random", N, q, yr, setup); g <- me_geometry(P, setup, need_engine = FALSE)
  key <- list("random", N, q)
  d0 <- me_design_inputs(g, yr, setup, key)                    # production rule
  d2 <- me_design_inputs(g, yr, setup, key, min_region = 2L)   # wrapper
  ok0 <- min(tabulate(d0$regions)) >= 2
  res[[q]] <- data.frame(N = N, prod_ok = ok0, identical = identical(as.integer(d0$regions), as.integer(d2$regions)),
                         changed = d2$regions_changed, min2 = min(tabulate(d2$regions)),
                         n_reg = length(unique(d2$regions)),
                         blocks_same = identical(d0$blocks, d2$blocks), X_same = identical(d0$X, d2$X))
}
res <- do.call(rbind, res)
check("Wrapper equals production partition whenever production's has regions >= 2",
      all(res$identical[res$prod_ok]), sprintf("%d of %d subsets qualify", sum(res$prod_ok), nrow(res)))
check("Wrapper flags 'changed' exactly when production's partition had a region < 2",
      identical(res$changed, !res$prod_ok), sprintf("changed %d", sum(res$changed)))
check("Wrapper regions always >= 2 colleges, always 4 regions", all(res$min2 >= 2) && all(res$n_reg == 4))
check("Wrapper leaves 2x2 blocks and ranks unchanged (RNG reseeded per key)", all(res$blocks_same) && all(res$X_same))
by_n <- aggregate(cbind(changed) ~ N, res, mean)
check("Changed share is small and falls with N (V6 rates 8% at 15, ~1% at 20)",
      by_n$changed[by_n$N == 15] < 0.2 && by_n$changed[by_n$N == 29] == 0,
      paste(sprintf("N=%d %.3f", by_n$N, by_n$changed), collapse = "; "))

# t + df-corrected test ----
# At N = 58 with intercept, Z, Spill and rho (p = 4): half-width multiplier
# t(0.975, 54) x sqrt(58/54) = 2.0049 x 1.0364 = 2.078 (vs z = 1.960).
mult <- stats::qt(0.975, 54) * sqrt(58 / 54)
check("t + df multiplier at N = 58, p = 4 is 2.078", abs(mult - 2.078) < 5e-4, sprintf("%.4f", mult))
# Pilot task: the stored counts obey n_rej_tdf <= n_rej, n_cov_tdf >= n_cov.
cfg <- mp_config("pilotv2"); tk <- mp_tasks(cfg); tk <- tk[tk$N == 20 & tk$Design_ID == 9, ][1, ]
cfg$S <- 2L; cfg$J <- 2L; cfg$R <- 5L
r <- mp_task(tk, cfg, setup)
st <- r$stats
check("Runner stores t + df stats; rejects no more and covers no less than Wald",
      all(c("n_rej_tdf", "n_cov_tdf") %in% dimnames(st)[[4]]) &&
      all(st[, , , "n_rej_tdf"] <= st[, , , "n_rej"]) && all(st[, , , "n_cov_tdf"] >= st[, , , "n_cov"]))
check("Runner records Regions_Changed per subset", "Regions_Changed" %in% names(r$subsets))

dir.create(file.path(.me_root, "verification"), recursive = TRUE, showWarnings = FALSE)
writeLines(c(sprintf("Run %s", format(Sys.time())), log_lines, if (length(failed)) paste("FAILED:", failed)),
           file.path(.me_root, "verification", "test_methods_ext_v2.txt"))
if (length(failed)) { cat("FAILED CHECKS:\n", paste(failed, collapse = "\n"), "\n"); quit(status = 1) }
cat("All v2 checks passed.\n")
