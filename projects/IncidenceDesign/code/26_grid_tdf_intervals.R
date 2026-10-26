# ============================================================
# Script: 26_grid_tdf_intervals.R
# Purpose: Re-fit the full grid benchmark exactly as 05_run_simulation.R did and
#          record coverage and power under the t + df-corrected interval beside the
#          stored Wald metrics, without touching the production outputs.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: spdep, spatialreg, dplyr, digest, parallel (via 05)
# ============================================================
#
# Run from the project root (BLAS pinned to one thread, as 05 requires):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript code/26_grid_tdf_intervals.R
#
# Why a re-fit: 05's summarize_fits() kept only scenario/surface summaries, so the
# per-fit estimates and SEs needed for another interval were never stored. Every draw
# in 05 is key-seeded (Invariant 4), so sourcing 05's functions and re-running each
# unit reproduces the production fits exactly; this script checks that by comparing
# its Wald metrics with results/sim_data/*_20260927_000502.rds (must match exactly).
#
# Interval (author decision 2026-10-21, as in the application):
#   tau-hat +/- t(0.975, n - p) x SE x sqrt(n / (n - p)),  n = 100 clusters,
#   p = regression coefficients + rho: oracle (1, Z, Spill, X) p = 5 -> df 95;
#   non-oracle (1, Z, X) p = 4 -> df 96.
# Power keeps 05's one-sided definition, now with the t + df lower bound > 0.
# Rook Checkerboard's oracle fit drops the aliased Spill column in every fit, so it
# is the non-oracle model (p = 4, df 96) and takes that df. Rook Isolation Buffer has
# one aliased fit (of 250) in 20 scenarios; those keep p = 5 (negligible).
#
# Outputs (results/grid_tdf_2026-10/; never results/sim_data/, so
# load_latest_results() is unaffected):
#   checkpoints/unit_*.rds + manifest.rds   per-unit results, manifest-checked
#   grid_tdf_scenarios.rds                  14,400 x 2 scenario rows: production keys
#                                           and Wald metrics plus Coverage_Tdf,
#                                           Power_Tdf, SE_Coverage_Tdf, SE_Power_Tdf
#   grid_tdf_surfaces.rds                   surface rows with the same additions
#   reproduction_check.txt                  Wald re-fit vs stored production metrics
#
# Parallelization: forked mclapply over 05's 40 units (n_cores = 10, ~15 min on
# this Mac); memory per worker is small. Not written for SLURM/Longleaf.

sim_define_only <- TRUE
sim_profile     <- "full"
sim_script_dir  <- normalizePath("code")
source("code/05_run_simulation.R")          # functions, parameters, grid; no run

stopifnot(all(Sys.getenv(c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS",
                           "OMP_NUM_THREADS")) == "1"))

OUT_DIR  <- "results/grid_tdf_2026-10"
CP_DIR   <- file.path(OUT_DIR, "checkpoints")
PROD_TAG <- "20260927_000502"               # production run being reproduced
dir.create(CP_DIR, recursive = TRUE, showWarnings = FALSE)

# Interval constants ----
P_EST  <- c(oracle = 5, nonoracle = 4)      # coefficients + rho
DF_EST <- N - P_EST                         # 95, 96
MULT   <- qt(0.975, DF_EST) * sqrt(N / DF_EST)   # replaces qnorm(0.975)

# Wrap 05's summarize_fits(): unchanged Wald output plus t + df columns ----
# run_unit() calls summarize_fits() once per (scenario, estimator), oracle first;
# the estimator is not an argument, so both df choices are computed and the right
# one is selected after the run by the Estimator key.
summarize_fits_wald <- summarize_fits

#' Wald summaries from 05 plus t + df-corrected coverage and power for both models
#'
#' @param tau_hat,se,surface,n_treated,true_tau As in 05's summarize_fits()
#' @return List(scen, surf) with columns Coverage_T95/T96, Power_T95/T96 and their
#'   surface-level Monte Carlo SEs appended (T95 = oracle df, T96 = non-oracle df)
summarize_fits <- function(tau_hat, se, surface, n_treated, true_tau) {
  out   <- summarize_fits_wald(tau_hat, se, surface, n_treated, true_tau)
  valid <- is.finite(tau_hat) & is.finite(se)
  ks    <- sort(unique(surface))
  for (nm in names(MULT)) {
    tag <- paste0("T", DF_EST[[nm]])
    lo  <- tau_hat - MULT[[nm]] * se
    hi  <- tau_hat + MULT[[nm]] * se
    cov <- lo <= true_tau & hi >= true_tau
    pow <- lo > 0
    cov_k <- vapply(ks, function(k) { v <- valid & surface == k
      if (any(v)) mean(cov[v]) else NA_real_ }, numeric(1))
    pow_k <- vapply(ks, function(k) { v <- valid & surface == k
      if (any(v)) mean(pow[v]) else NA_real_ }, numeric(1))
    ok_k <- out$surf$N_Valid > 0
    out$surf[[paste0("Coverage_", tag)]] <- cov_k
    out$surf[[paste0("Power_", tag)]]    <- pow_k
    out$scen[[paste0("Coverage_", tag)]] <- if (any(valid)) mean(cov[valid]) else NA_real_
    out$scen[[paste0("Power_", tag)]]    <- if (any(valid)) mean(pow[valid]) else NA_real_
    out$scen[[paste0("SE_Coverage_", tag)]] <- sd(cov_k[ok_k]) / sqrt(sum(ok_k))
    out$scen[[paste0("SE_Power_", tag)]]    <- sd(pow_k[ok_k]) / sqrt(sum(ok_k))
  }
  out
}

# Run (manifest-checked checkpoints, as in 05) ----
# The checkpoints depend on 01-05 and on the summarize_fits() wrapper above, so the
# manifest hashes those (not this whole file: post-processing below can change
# without invalidating the checkpoints).
manifest <- list(
  code_hashes = c(vapply(code_files, digest::digest, character(1), file = TRUE),
                  "summarize_fits_tdf" = digest::digest(deparse(summarize_fits))),
  packages = vapply(c("spatialreg", "spdep", "dplyr", "digest"),
                    function(p) as.character(packageVersion(p)), character(1)),
  R = R.version.string, blas = sessionInfo()$BLAS, mult = MULT)
names(manifest$code_hashes)[seq_along(code_files)] <- basename(code_files)
mf <- file.path(CP_DIR, "manifest.rds")
if (file.exists(mf)) {
  if (!identical(readRDS(mf), manifest)) {
    stop("Manifest mismatch in ", CP_DIR, "; move or delete it to start fresh.")
  }
} else {
  if (length(list.files(CP_DIR, "^unit_"))) stop("Checkpoints without a manifest in ", CP_DIR)
  saveRDS(manifest, mf)
}

units <- expand.grid(cfg_index = seq_along(inc_configs), nb_type = neighbor_types,
                     rho = rho_vals, stringsAsFactors = FALSE)
surfaces <- lapply(inc_configs, generate_surfaces)
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(units)), function(i) {
  u  <- units[i, ]
  cp <- file.path(CP_DIR, sprintf("unit_%s_rhoX%.2f_%s_rho%.2f.rds",
                  inc_configs[[u$cfg_index]]$mode, inc_configs[[u$cfg_index]]$rho_x,
                  u$nb_type, u$rho))
  if (file.exists(cp)) return(readRDS(cp))
  out <- run_unit(u$cfg_index, u$nb_type, u$rho, surfaces[[u$cfg_index]])
  saveRDS(out, cp)
  out
}, mc.cores = n_cores, mc.preschedule = FALSE)
bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
if (any(bad)) stop("Units failed: ", paste(which(bad), collapse = ", "))
cat(sprintf("Re-fit complete in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Select each estimator's df and save ----
scen_keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho",
               "Gamma", "Spillover_Type", "True_Tau", "Estimator")
scen_raw <- dplyr::bind_rows(lapply(res, `[[`, "scen"))
# Fully aliased oracle scenarios (rook Checkerboard): every fit is the non-oracle model
aliased_key <- do.call(paste, scen_raw[scen_raw$Estimator == "oracle" &
  scen_raw$N_Aliased > 0 & scen_raw$N_Aliased == scen_raw$N_Valid_Est, scen_keys])

#' Keep each row's own df: 95 for oracle fits, 96 for non-oracle and fully aliased
#' oracle fits (rook Checkerboard, where both models coincide)
#' @param d Scenario or surface rows carrying the _T95 and _T96 columns
#' @return d with *_Tdf columns in place of the _T95/_T96 pairs
pick <- function(d) {
  t95 <- d$Estimator == "oracle" & !(do.call(paste, d[scen_keys]) %in% aliased_key)
  for (m in c("Coverage", "Power", "SE_Coverage", "SE_Power")) {
    a <- paste0(m, "_T95"); b <- paste0(m, "_T96")
    if (a %in% names(d)) {
      d[[paste0(m, "_Tdf")]] <- ifelse(t95, d[[a]], d[[b]])
      d[[a]] <- NULL; d[[b]] <- NULL
    }
  }
  d
}
scen <- pick(scen_raw)
surf <- pick(dplyr::bind_rows(lapply(res, `[[`, "surf")))
saveRDS(scen, file.path(OUT_DIR, "grid_tdf_scenarios.rds"))
saveRDS(surf, file.path(OUT_DIR, "grid_tdf_surfaces.rds"))

# Reproduction check: Wald metrics must equal the stored production run ----
keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho", "Gamma",
          "Spillover_Type", "True_Tau", "Estimator")
mets <- c("Mean_Estimate", "Bias", "SD", "MSE", "Coverage", "Power", "N_Valid_Est",
          "SE_Coverage", "N_Aliased")
prod <- dplyr::bind_rows(
  readRDS(sprintf("results/sim_data/sim_results_MLE_tau_sweep_combined_%s.rds", PROD_TAG)),
  readRDS(sprintf("results/sim_data/sim_results_MLEnonoracle_tau_sweep_combined_%s.rds", PROD_TAG)))
m <- merge(prod[, c(keys, mets)], scen[, c(keys, mets)], by = keys,
           suffixes = c(".prod", ".refit"))
diffs <- vapply(mets, function(v) {
  a <- m[[paste0(v, ".prod")]]; b <- m[[paste0(v, ".refit")]]
  max(abs(a - b)[is.finite(a) | is.finite(b)], na.rm = TRUE)
}, numeric(1))
lines <- c(sprintf("Rows: production %d, re-fit %d, matched %d", nrow(prod), nrow(scen), nrow(m)),
           sprintf("max |production - re-fit|  %-14s %.3g", mets, diffs))
writeLines(lines, file.path(OUT_DIR, "reproduction_check.txt"))
cat(lines, sep = "\n")
if (nrow(m) != nrow(prod) || nrow(m) != nrow(scen) || any(diffs > 1e-12)) {
  stop("Re-fit does not reproduce the production run; see reproduction_check.txt")
}
cat("Reproduction exact; t + df metrics saved to", OUT_DIR, "\n")
