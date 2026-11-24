# ============================================================
# Script: 29_igs_sensitivity.R
# Purpose: Sensitivity study of Incidence-Guided Saturation Quadrants: vary the number
#          of regions and the spread of the saturation set, each with a matched
#          unguided twin, on the full 10x10 grid benchmark, beside the original design,
#          Saturation Quadrants and Simple Random Sampling.
# Author: Andrew Walther
# Created: 2026-11-24
# Dependencies: spdep, spatialreg, dplyr, digest, parallel (via 05)
# ============================================================
#
# Question (author, 2026-11-24): the study tested one version of incidence guidance
# (four 5x5 quadrants, saturations 80/60/40/20%). Do versions with other numbers of
# regions or other saturation sets perform similarly relative to SRS and to the
# original version?
#
# Run from the project root (BLAS pinned to one thread, as 05 requires):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript code/29_igs_sensitivity.R smoke    # small check, separate folder
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript code/29_igs_sensitivity.R full
# Summaries: code/30_igs_sensitivity_summary.R.
#
# Design of the study (provisional; the spec is in
# docs/plans/2026-11-24-igs-sensitivity-spec.md and awaits the author's approval):
#
#   Regions R (equal-size partitions of the 10x10 grid, so 50% can be treated exactly):
#     R = 2   two 5x10 halves (x <= 5, x > 5)
#     R = 4   four 5x5 quadrants (the original design)
#     R = 10  ten 2x5 rectangles
#     R = 25  twenty-five 2x2 blocks (the 2x2 Blocking blocks)
#   Saturation set with spread s, by incidence rank r = 1 (highest mean) .. R:
#     sat_r = 0.5 + s * (1 - 2 (r - 1) / (R - 1)),   s in {0.1, 0.2, 0.3, 0.4}
#     R = 4, s = 0.3 gives {0.8, 0.6, 0.4, 0.2}, the original set.
#   Counts: n_r = size_r * sat_r rounded by largest remainder so sum(n_r) = N/2
#     exactly; equal remainders are broken at random per draw.
#   Guided: region r gets sat_r by the rank of its mean incidence (ties at random per
#     draw, as Design 8). Unguided twin: the same saturations in random order per draw
#     (as Design 3). The pair separates the effect of using incidence from the effect
#     of the partition and the saturation set.
#   The (R = 4, s = 0.3) pair is not duplicated: it is Design 8 (guided) and Design 3
#     (unguided), which run here unchanged beside Design 9 (SRS).
#
# Reuse without editing hashed code: 05 is sourced with sim_define_only = TRUE (as in
# 26); get_designs() is wrapped so production IDs 1-9 go to 03 unchanged and the
# variant IDs 101-130 go to igs_variant_designs() in 29_igs_variant_designs.R.
# Variant draws are keyed by their own IDs through 05's draw_assignments(), so they
# never share seeds with the production designs. Designs 3, 8 and 9 reuse their production keys and must
# reproduce the 2026-09-27 run exactly; the script stops otherwise. The noise E is
# keyed by block only, so every design in a block sees the same noise (paired
# comparisons).
#
# Intervals: t + df-corrected (author decision 2026-10-26) via the 26 wrapper logic;
# Wald kept beside it.
#
# Outputs (results/igs_sensitivity_2026-11/<profile>/; never results/sim_data/):
#   checkpoints/unit_*.rds + manifest.rds  per-unit results, manifest-checked
#   igs_scenarios.rds, igs_surfaces.rds    scenario and surface rows, both estimators
#   igs_variants.csv                       variant ID -> R, s, guided, label
#   reproduction_check.txt                 Designs 3, 8, 9 vs the production run
#
# Parallelization: forked mclapply over 05's 40 units (n_cores = 10). 33 designs
# against production's 9 makes this about 3.7 times the production fit count; expect
# about an hour on this Mac. Not written for SLURM/Longleaf.

profile <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(profile)) profile <- "full"
stopifnot(profile %in% c("smoke", "full"))

sim_define_only <- TRUE
sim_profile     <- if (profile == "smoke") "test" else "full"
sim_script_dir  <- normalizePath("code")
source("code/05_run_simulation.R")          # functions, parameters, grid; no run

stopifnot(all(Sys.getenv(c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS",
                           "OMP_NUM_THREADS")) == "1"))

OUT_DIR  <- file.path("results/igs_sensitivity_2026-11", profile)
CP_DIR   <- file.path(OUT_DIR, "checkpoints")
PROD_TAG <- "20260927_000502"               # production run reproduced for 3, 8, 9
dir.create(CP_DIR, recursive = TRUE, showWarnings = FALSE)

# Variant designs (VARIANTS, region_index, treated_counts, igs_variant_designs) ----
source("code/29_igs_variant_designs.R")

# Route variant IDs to the new designs; production IDs go to 03 unchanged ----
get_designs_prod <- get_designs
get_designs <- function(design_id, n_resamples, N, incidence, nb_list, coords) {
  get_designs_with_variants(design_id, n_resamples, N, incidence, nb_list, coords,
                            prod_fun = get_designs_prod)
}

# Self-checks before any fitting ----
local({
  for (R in c(2, 4, 10, 25)) stopifnot(all(tabulate(region_index(coords, R)) == N / R))
  stopifnot(identical(treated_counts(rep(25, 4), c(0.8, 0.6, 0.4, 0.2), 50),
                      c(20, 15, 10, 5)))
  set.seed(1)
  x <- runif(N)
  for (j in seq_len(nrow(VARIANTS))) {
    zz <- igs_variant_designs(VARIANTS[j, ], 5, N, x, coords)
    stopifnot(all(colSums(zz) == N / 2))
  }
  # Guided R = 25, s = 0.4: the highest-mean block must be (nearly) fully treated
  v  <- VARIANTS[VARIANTS$R == 25 & VARIANTS$s == 0.4 & VARIANTS$guided, ]
  rg <- region_index(coords, 25)
  top <- which.max(tapply(x, rg, mean))
  zz <- igs_variant_designs(v, 20, N, x, coords)
  stopifnot(all(colSums(zz[rg == top, , drop = FALSE]) >= 3))
})

design_ids <- c(3, 8, 9, VARIANTS$id)
if (profile == "smoke") design_ids <- c(8, 9, VARIANTS$id[c(1, 8, 16, 17, 30)])
write.csv(VARIANTS, file.path(OUT_DIR, "igs_variants.csv"), row.names = FALSE)

# t + df interval wrapper (logic of 26) ----
P_EST  <- c(oracle = 5, nonoracle = 4)
DF_EST <- N - P_EST
MULT   <- qt(0.975, DF_EST) * sqrt(N / DF_EST)
summarize_fits_wald <- summarize_fits

#' Wald summaries from 05 plus t + df coverage and power for both df choices
#' @param tau_hat,se,surface,n_treated,true_tau As in 05's summarize_fits()
#' @return List(scen, surf) with Coverage/Power _T95 and _T96 columns appended
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

# Run (manifest-checked checkpoints) ----
manifest <- list(
  code_hashes = c(vapply(code_files, digest::digest, character(1), file = TRUE),
                  "summarize_fits_tdf" = digest::digest(deparse(summarize_fits)),
                  "variant_designs"    = digest::digest(list(deparse(igs_variant_designs),
                                                             deparse(treated_counts),
                                                             deparse(region_index),
                                                             VARIANTS))),
  design_ids = design_ids, profile = profile,
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
cat(sprintf("Run complete in %.1f min (%d units, %d designs)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(units),
            length(design_ids)))

# Each estimator keeps its own df (oracle 95, non-oracle 96) ----
# No variant is aliased the way rook Checkerboard is (all are random with N/2
# treated), so oracle rows take df 95; aliased fits are counted and reported.
pick <- function(d) {
  t95 <- d$Estimator == "oracle"
  for (m in c("Coverage", "Power", "SE_Coverage", "SE_Power")) {
    a <- paste0(m, "_T95"); b <- paste0(m, "_T96")
    if (a %in% names(d)) {
      d[[paste0(m, "_Tdf")]] <- ifelse(t95, d[[a]], d[[b]])
      d[[a]] <- NULL; d[[b]] <- NULL
    }
  }
  d
}
scen <- pick(dplyr::bind_rows(lapply(res, `[[`, "scen")))
surf <- pick(dplyr::bind_rows(lapply(res, `[[`, "surf")))
warn <- dplyr::bind_rows(lapply(res, `[[`, "warn"))
saveRDS(scen, file.path(OUT_DIR, "igs_scenarios.rds"))
saveRDS(surf, file.path(OUT_DIR, "igs_surfaces.rds"))
saveRDS(warn, file.path(OUT_DIR, "igs_warnings.rds"))
cat(sprintf("Fail rate max %.4f; aliased fits %d; other warnings %d\n",
            max(scen$Fail_Rate), sum(scen$N_Aliased), sum(scen$N_Warn)))

# Reproduction check: Designs 3, 8, 9 must equal the production run (full only) ----
if (profile == "full") {
  keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho", "Gamma",
            "Spillover_Type", "True_Tau", "Estimator")
  mets <- c("Mean_Estimate", "Bias", "SD", "MSE", "Coverage", "Power", "N_Valid_Est")
  prod <- dplyr::bind_rows(
    readRDS(sprintf("results/sim_data/sim_results_MLE_tau_sweep_combined_%s.rds", PROD_TAG)),
    readRDS(sprintf("results/sim_data/sim_results_MLEnonoracle_tau_sweep_combined_%s.rds",
                    PROD_TAG)))
  prod <- prod[prod$Design %in% c("Design 3", "Design 8", "Design 9"), ]
  mine <- scen[scen$Design %in% c("Design 3", "Design 8", "Design 9"), ]
  m <- merge(prod[, c(keys, mets)], mine[, c(keys, mets)], by = keys,
             suffixes = c(".prod", ".rerun"))
  diffs <- vapply(mets, function(v) {
    a <- m[[paste0(v, ".prod")]]; b <- m[[paste0(v, ".rerun")]]
    max(abs(a - b)[is.finite(a) | is.finite(b)], na.rm = TRUE)
  }, numeric(1))
  lines <- c(sprintf("Rows: production %d, rerun %d, matched %d", nrow(prod), nrow(mine),
                     nrow(m)),
             sprintf("max |production - rerun|  %-14s %.3g", mets, diffs))
  writeLines(lines, file.path(OUT_DIR, "reproduction_check.txt"))
  cat(lines, sep = "\n")
  if (nrow(m) != nrow(prod) || nrow(m) != nrow(mine) || any(diffs > 1e-12)) {
    stop("Designs 3, 8, 9 do not reproduce the production run; see reproduction_check.txt")
  }
  cat("Reproduction exact for Designs 3, 8, 9\n")
}
cat("Saved to", OUT_DIR, "\n")
