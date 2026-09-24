# ============================================================
# Script: test_estimator_equivalence.R
# Purpose: Validate fit_sar_lag() against spatialreg::lagsarlm() on >= 5,000 fits
#          before any simulation code depends on it (simulation-revision-plan.md B3).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-24
# Dependencies: spatialreg, spdep; sources code/01-04
# ============================================================
#
# Run from code/:  Rscript tests/test_estimator_equivalence.R
# Exits non-zero if any pass criterion fails. Writes per-fit differences and a
# summary to results/estimator_validation/.
#
# Coverage: every combination of 5 incidence configs x 2 neighbor types x 4 rho x
# 2 spillover regimes x 8 designs (640 data cells), 4 outcome draws per cell, each
# fit by both engines under the oracle (Y ~ Z + Spill + X) and non-oracle
# (Y ~ Z + X) models: 640 x 4 x 2 = 5,120 fits. Includes queen Checkerboard and the
# aliased Checkerboard x rook oracle fits (WZ = 1 - Z exactly).
#
# Pass criteria (plan B3): |d tau-hat| < 1e-6, |d SE(tau-hat)| / SE < 1e-5,
# |d rho-hat| < 1e-6, |d logLik| < 1e-6, identical aliasing and identical warnings.
#
# The data-generating process here is only a source of realistic inputs: designs
# are drawn from surface X[, k] and the outcome uses the same X[, k] (revision M1),
# Poisson surfaces use 100,000 people per cluster (M2). The estimator comparison does
# not depend on those choices.

# Setup ----
suppressMessages({
  source("01_spatial_setup.R")
  source("02_incidence_generation.R")
  source("03_designs.R")
  source("04_estimation.R")
  library(spatialreg)
})
set.seed(20260924, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")

out_dir <- file.path("..", "results", "estimator_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g <- build_spatial_grid(grid_dim = 10)
N <- g$N_clusters
I_mat <- diag(N)
n_surfaces <- 10
n_draws_per_cell <- 4

configs <- list(
  list(mode = "iid",     rho_x = 0),
  list(mode = "spatial", rho_x = 0.20),
  list(mode = "spatial", rho_x = 0.50),
  list(mode = "poisson", rho_x = 0.20),
  list(mode = "poisson", rho_x = 0.50)
)
X_by_config <- lapply(configs, function(cf) {
  generate_incidence(cf$mode, N, n_surfaces, W = g$W_queen,
                     rho_incidence = cf$rho_x, pop_per_cluster = 100000)
})
setups <- list(rook = sar_lag_setup(g$listw_rook),
               queen = sar_lag_setup(g$listw_queen))

cells <- expand.grid(cfg = seq_along(configs), nb = c("rook", "queen"),
                     rho = c(0, 0.01, 0.20, 0.50),
                     spill = c("control_only", "both"), d = 1:8,
                     stringsAsFactors = FALSE)

## Capture every warning an engine raises, without stopping it ----
run_capturing <- function(expr) {
  warns <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "fit_error")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warns = warns)
}

# Fit loop ----
rows <- vector("list", nrow(cells) * n_draws_per_cell * 2)
t_lagsarlm <- 0
t_lean <- 0
i_row <- 0
for (ci in seq_len(nrow(cells))) {
  cl <- cells[ci, ]
  sp <- if (cl$nb == "rook") g$W_rook else g$W_queen
  lw <- if (cl$nb == "rook") g$listw_rook else g$listw_queen
  nbl <- if (cl$nb == "rook") g$nb_rook else g$nb_queen
  k <- sample.int(n_surfaces, 1)
  X <- X_by_config[[cl$cfg]][, k]
  gamma <- sample(c(0.5, 0.6, 0.7, 0.8), 1)
  tau <- sample(c(0.8, 1.0, 1.5, 2.0, 3.0), 1)
  Z <- get_designs(cl$d, 1, N, X, nbl, g$coords)[, 1]
  WZ <- as.vector(sp %*% Z)
  spill <- if (cl$spill == "control_only") gamma * WZ * (1 - Z) else gamma * WZ
  inv <- solve(I_mat - cl$rho * sp)

  for (r in seq_len(n_draws_per_cell)) {
    # Y = (I - rho W)^{-1} (tau Z + spill + beta X + eps), beta = 1, sigma = 1
    Y <- as.vector(inv %*% (tau * Z + spill + X + rnorm(N)))
    for (model in c("oracle", "nonoracle")) {
      if (model == "oracle") {
        df <- data.frame(Y = Y, Z = Z, Spill = spill, X = X)
        fml <- Y ~ Z + Spill + X
        xm <- cbind("(Intercept)" = 1, Z = Z, Spill = spill, X = X)
      } else {
        df <- data.frame(Y = Y, Z = Z, X = X)
        fml <- Y ~ Z + X
        xm <- cbind("(Intercept)" = 1, Z = Z, X = X)
      }
      t0 <- proc.time()[["elapsed"]]
      ref <- run_capturing(lagsarlm(fml, data = df, listw = lw, quiet = TRUE))
      t1 <- proc.time()[["elapsed"]]
      new <- run_capturing(fit_sar_lag(Y, xm, setups[[cl$nb]]))
      t2 <- proc.time()[["elapsed"]]
      t_lagsarlm <- t_lagsarlm + (t1 - t0)
      t_lean <- t_lean + (t2 - t1)

      stopifnot(!inherits(ref$value, "fit_error"), !inherits(new$value, "fit_error"))
      L <- ref$value
      f <- new$value
      L_coef <- coef(L)[-1]                     # coef.Sarlm prepends rho
      common <- intersect(names(L_coef), names(f$coefficients))
      i_row <- i_row + 1
      rows[[i_row]] <- data.frame(
        cell = ci, cfg = cl$cfg, nb = cl$nb, rho_true = cl$rho, spill = cl$spill,
        design = cl$d, model = model, draw = r,
        aliased_ref = paste(names(L$aliased)[L$aliased], collapse = " "),
        aliased_new = paste(names(f$aliased)[f$aliased], collapse = " "),
        warns_ref = paste(ref$warns, collapse = " | "),
        warns_new = paste(new$warns, collapse = " | "),
        same_coef_names = identical(names(L_coef), names(f$coefficients)),
        tau_ref = L_coef[["Z"]], tau_new = f$coefficients[["Z"]],
        se_ref = L$rest.se[["Z"]], se_new = f$rest.se[["Z"]],
        rho_ref = L$rho[[1]], rho_new = f$rho,
        ll_ref = L$LL[[1]], ll_new = f$LL,
        max_dcoef = max(abs(L_coef[common] - f$coefficients[common])),
        max_rel_dse = max(abs(L$rest.se[common] - f$rest.se[common]) / L$rest.se[common]),
        d_rho_se = abs(L$rho.se - f$rho.se),
        stringsAsFactors = FALSE
      )
    }
  }
}
res <- do.call(rbind, rows[seq_len(i_row)])

# Criteria ----
res$d_tau <- abs(res$tau_ref - res$tau_new)
res$rel_d_se <- abs(res$se_ref - res$se_new) / res$se_ref
res$d_rho <- abs(res$rho_ref - res$rho_new)
res$d_ll <- abs(res$ll_ref - res$ll_new)
res$is_cb_rook_oracle <- res$design == 1 & res$nb == "rook" & res$model == "oracle"

checks <- c(
  n_fits_ge_5000        = nrow(res) >= 5000,
  tau_lt_1e6            = all(res$d_tau < 1e-6),
  se_rel_lt_1e5         = all(res$rel_d_se < 1e-5),
  rho_lt_1e6            = all(res$d_rho < 1e-6),
  loglik_lt_1e6         = all(res$d_ll < 1e-6),
  all_coef_lt_1e6       = all(res$max_dcoef < 1e-6),
  all_se_rel_lt_1e5     = all(res$max_rel_dse < 1e-5),
  identical_aliasing    = all(res$aliased_ref == res$aliased_new),
  identical_coef_names  = all(res$same_coef_names),
  identical_warnings    = all(res$warns_ref == res$warns_new),
  aliased_iff_cb_rook_oracle = all((res$aliased_ref != "") == res$is_cb_rook_oracle)
)

fmt <- function(x) formatC(x, format = "e", digits = 2)
summary_lines <- c(
  sprintf("fit_sar_lag() vs spatialreg::lagsarlm() — %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("R %s, spatialreg %s, spdep %s", getRversion(),
          packageVersion("spatialreg"), packageVersion("spdep")),
  sprintf("Fits: %d (%d oracle, %d non-oracle); aliased fits: %d (Checkerboard x rook oracle: %d)",
          nrow(res), sum(res$model == "oracle"), sum(res$model == "nonoracle"),
          sum(res$aliased_ref != ""), sum(res$is_cb_rook_oracle)),
  sprintf("Fits with any warning (lagsarlm): %d; distinct warnings: %s",
          sum(res$warns_ref != ""),
          paste(unique(res$warns_ref[res$warns_ref != ""]), collapse = " || ")),
  "",
  "Max absolute / relative differences (threshold):",
  sprintf("  tau-hat          %s  (< 1e-6)", fmt(max(res$d_tau))),
  sprintf("  SE(tau-hat) rel  %s  (< 1e-5)", fmt(max(res$rel_d_se))),
  sprintf("  rho-hat          %s  (< 1e-6)", fmt(max(res$d_rho))),
  sprintf("  logLik           %s  (< 1e-6)", fmt(max(res$d_ll))),
  sprintf("  any coefficient  %s  (< 1e-6)", fmt(max(res$max_dcoef))),
  sprintf("  any SE rel       %s  (< 1e-5)", fmt(max(res$max_rel_dse))),
  sprintf("  SE(rho-hat) abs  %s  (not a criterion)", fmt(max(res$d_rho_se))),
  "",
  sprintf("Elapsed: lagsarlm %.1f s, fit_sar_lag %.1f s, speedup %.1fx (%.2f vs %.2f ms/fit)",
          t_lagsarlm, t_lean, t_lagsarlm / t_lean,
          1000 * t_lagsarlm / nrow(res), 1000 * t_lean / nrow(res)),
  "",
  "Checks:",
  sprintf("  %-28s %s", names(checks), ifelse(checks, "PASS", "FAIL"))
)
writeLines(summary_lines, file.path(out_dir, "equivalence_summary.txt"))
write.csv(res, file.path(out_dir, "equivalence_per_fit.csv"), row.names = FALSE)
cat(summary_lines, sep = "\n")

if (!all(checks)) {
  stop("Estimator equivalence FAILED: ", paste(names(checks)[!checks], collapse = ", "))
}
