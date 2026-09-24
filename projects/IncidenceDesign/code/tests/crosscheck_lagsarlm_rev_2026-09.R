# ============================================================
# Script: crosscheck_lagsarlm_rev_2026-09.R
# Purpose: B5 1% cross-check — recompute a random 1% of full-run scenarios end to end
#          with engine = "lagsarlm" and compare with the stored (lean-engine) rows.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-24
# Dependencies: dplyr, parallel; sources code/05 (define-only)
# ============================================================
#
# Run from code/ after the full run:
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript tests/crosscheck_lagsarlm_rev_2026-09.R
# ~128 scenarios x 250 fits x 2 estimators with lagsarlm (~76 ms/fit): ~8 min on 10 cores.
#
# Why this reproduces the stored rows exactly: every draw is keyed (spec section 4), so
# restricting run_unit()'s loops to one (gamma, regime, design, tau) regenerates the same
# eps, Z and X as the full run. Only the estimation engine differs.
# Writes results/estimator_validation/crosscheck_full_run.{txt,csv}; exits non-zero on
# failure.

sim_profile <- "full"; sim_define_only <- TRUE; sim_script_dir <- normalizePath(".")
suppressMessages(source("05_run_simulation.R"))

sim_dir <- file.path("..", "results", "sim_data")
latest <- function(pattern) {
  f <- list.files(sim_dir, pattern = pattern, full.names = TRUE)
  stopifnot(length(f) > 0)
  f[which.max(file.mtime(f))]
}
stored <- bind_rows(readRDS(latest("^sim_results_MLE_tau_sweep_combined_.*\\.rds$")),
                    readRDS(latest("^sim_results_MLEnonoracle_tau_sweep_combined_.*\\.rds$")))

keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho", "Gamma",
          "Spillover_Type", "True_Tau")
scen_keys <- distinct(stored[, keys])
stopifnot(nrow(scen_keys) == 12800)
set.seed(20260924, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
pick <- scen_keys[sample.int(nrow(scen_keys), 128), ]

surfaces <- lapply(inc_configs, generate_surfaces)
cfg_of <- function(mode, rx) which(vapply(inc_configs, function(c) c$mode == mode && c$rho_x == rx, logical(1)))

recompute <- function(i) {
  p <- pick[i, ]
  # Restrict run_unit()'s loops to this one scenario (children of mclapply own their globals)
  gamma_vals    <<- p$Gamma
  spill_types   <<- p$Spillover_Type
  design_ids    <<- as.integer(sub("Design ", "", p$Design))
  true_tau_vals <<- p$True_Tau
  engine        <<- "lagsarlm"
  ci <- cfg_of(p$Incidence_Mode, p$Rho_Incidence)
  run_unit(ci, p$Neighbor_Type, p$Rho, surfaces[[ci]])$scen
}
t0 <- Sys.time()
redo <- parallel::mclapply(seq_len(nrow(pick)), recompute, mc.cores = 10, mc.preschedule = FALSE)
bad <- vapply(redo, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
if (any(bad)) stop("Cross-check units failed: ", paste(which(bad), collapse = ", "))
redo <- bind_rows(redo)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

cmp <- inner_join(stored, redo, by = c(keys, "Estimator"), suffix = c("", ".lag"))
stopifnot(nrow(cmp) == 2 * nrow(pick))
cont <- c("Mean_Estimate", "Bias", "SD", "MSE", "SE_Bias", "SE_MSE")
disc <- c("Coverage", "Power", "N_Valid_Est", "N_Aliased", "N_Warn", "Mean_Treated",
          "Z_WZ_rank_deficient")
max_diff <- vapply(cont, function(v) max(abs(cmp[[v]] - cmp[[paste0(v, ".lag")]])), numeric(1))
n_mismatch <- vapply(disc, function(v) sum(cmp[[v]] != cmp[[paste0(v, ".lag")]]), integer(1))
# Coverage/Power are fractions of 250 fits: one CI endpoint within ~1e-7 of tau (or of 0)
# could flip between engines. Report counts; any flip is shown row by row below.
ok <- all(max_diff < 1e-6) && all(n_mismatch[c("N_Valid_Est", "N_Aliased", "N_Warn",
                                               "Mean_Treated", "Z_WZ_rank_deficient")] == 0) &&
      all(abs(cmp$Coverage - cmp$Coverage.lag) <= 1 / 250) &&
      all(abs(cmp$Power - cmp$Power.lag) <= 1 / 250)

out <- c(
  sprintf("1%% lagsarlm cross-check of the full run — %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Scenarios: %d (x 2 estimators = %d rows; %d lagsarlm fits), %.1f min on 10 workers",
          nrow(pick), nrow(cmp), nrow(cmp) * 250, mins),
  "Max |stored - lagsarlm| (criterion < 1e-6):",
  sprintf("  %-14s %.2e", cont, max_diff),
  "Mismatch counts (criterion 0; Coverage/Power may differ by at most one fit):",
  sprintf("  %-20s %d", disc, n_mismatch),
  sprintf("RESULT: %s", if (ok) "PASS" else "FAIL"))
writeLines(out, file.path("..", "results", "estimator_validation", "crosscheck_full_run.txt"))
write.csv(cmp, file.path("..", "results", "estimator_validation", "crosscheck_full_run.csv"),
          row.names = FALSE)
cat(out, sep = "\n")
flips <- cmp[cmp$Coverage != cmp$Coverage.lag | cmp$Power != cmp$Power.lag, ]
if (nrow(flips) > 0) print(flips[, c(keys, "Estimator", "Coverage", "Coverage.lag", "Power", "Power.lag")])
if (!ok) stop("Cross-check FAILED")
