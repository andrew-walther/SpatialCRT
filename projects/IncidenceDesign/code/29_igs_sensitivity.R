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
# Study specification: docs/plans/2026-11-24-igs-sensitivity-spec.md.
#
#   Regions R (equal-size partitions of the 10x10 grid, so 50% can be treated exactly):
#     R = 2   two 5x10 halves (x <= 5, x > 5)
#     R = 4   four 5x5 quadrants (the original design)
#     R = 10  ten 2x5 rectangles
#     R = 25  twenty-five 2x2 blocks (the 2x2 Blocking blocks)
#   Saturation set with spread s, by incidence rank r = 1 (highest mean) .. R:
#     sat_r = 0.5 + s * (1 - 2 (r - 1) / (R - 1)),   s in {0.1, 0.2, 0.3, 0.4, 0.5}
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
# variant IDs 101-138 go to igs_variant_designs() in 29_igs_variant_designs.R.
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
#   checkpoints/boundary/                 per-unit results, manifest-checked
#   igs_scenarios.rds, igs_surfaces.rds    scenario and surface rows, both estimators
#   igs_variants.csv                       variant ID -> R, s, guided, label
#   igs_identification.csv                 ranks, aliasing, failures and warnings
#   reproduction_check.txt                 November 24 and September 27 checks
#
# Parallelization: forked mclapply over 05's 40 units with n_cores = 10 and one
# BLAS thread per worker. Reuse the manifest-checked November 24 fits and fit the
# eight boundary designs in about 20 minutes on this Mac.

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
BASE_CP  <- file.path(OUT_DIR, "checkpoints")
CP_DIR   <- file.path(BASE_CP, "boundary")
PROD_TAG <- "20260927_000502"               # production run reproduced for 3, 8, 9

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

design_ids <- VARIANTS$id[VARIANTS$s == 0.5]
stopifnot(identical(as.integer(design_ids), 131:138), n_cores == 10)
if (profile == "smoke") spill_types <- c("control_only", "both")

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
# Manifest and work units ----
manifest <- list(
  code_hashes = c(vapply(code_files, digest::digest, character(1), file = TRUE),
                  "summarize_fits_tdf" = digest::digest(deparse(summarize_fits)),
                  "variant_designs"    = digest::digest(list(deparse(igs_variant_designs),
                                                             deparse(treated_counts),
                                                             deparse(region_index),
                                                             ORIGINAL_VARIANTS))),
  design_ids = c(3, 8, 9, ORIGINAL_VARIANTS$id), profile = profile,
  packages = vapply(c("spatialreg", "spdep", "dplyr", "digest"),
                    function(p) as.character(packageVersion(p)), character(1)),
  R = R.version.string, blas = sessionInfo()$BLAS, mult = MULT)
names(manifest$code_hashes)[seq_along(code_files)] <- basename(code_files)

units <- expand.grid(cfg_index = seq_along(inc_configs), nb_type = neighbor_types,
                     rho = rho_vals, stringsAsFactors = FALSE)
unit_files <- with(units, sprintf("unit_%s_rhoX%.2f_%s_rho%.2f.rds",
  vapply(inc_configs[cfg_index], `[[`, character(1), "mode"),
  vapply(inc_configs[cfg_index], `[[`, numeric(1), "rho_x"), nb_type, rho))

keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho", "Gamma",
          "Spillover_Type", "True_Tau", "Estimator")
stats <- c("scen", "surf", "warn")
files <- setNames(file.path(OUT_DIR, paste0("igs_", c("scenarios", "surfaces", "warnings"),
                                           ".rds")), stats)

#' Combine per-unit fit summaries and select the t + df columns for each estimator
#' @param res List of unit results from run_unit()
#' @return Named list of scenario, surface and warning data frames
#' @examples combine_units(list(readRDS(file.path(BASE_CP, unit_files[1]))))
combine_units <- function(res) {
  setNames(lapply(stats, function(nm) {
    rows <- dplyr::bind_rows(lapply(res, `[[`, nm))
    if (nm == "warn") rows else pick(rows)
  }), stats)
}

#' Verify that a set of summaries has one row per scenario and per scenario-surface
#' @param summaries List of scenario, surface and warning data frames
#' @param ids Numeric vector of expected design IDs
#' @return Invisibly NULL, or an error identifying a missing or duplicated result
#' @examples check_rows(baseline, c(3, 8, 9, ORIGINAL_VARIANTS$id))
check_rows <- function(summaries, ids) {
  per_design <- nrow(units) * length(gamma_vals) * length(spill_types) *
    length(true_tau_vals) * 2L
  ids <- paste("Design", ids)
  for (nm in c("scen", "surf")) {
    rows <- summaries[[nm]]
    expected <- per_design * if (nm == "surf") n_surfaces else 1L
    row_keys <- if (nm == "surf") c(keys, "Surface") else keys
    counts <- table(rows$Design)
    if (!setequal(names(counts), ids) || any(counts != expected) ||
        anyNA(rows[row_keys]) || anyDuplicated(rows[row_keys])) {
      stop("Incomplete or duplicate ", nm, " rows")
    }
  }
  invisible(NULL)
}

baseline <- NULL
reproduction <- character()
if (profile == "full") {
  if (!identical(readRDS(file.path(BASE_CP, "manifest.rds")), manifest)) {
    stop("November 24 checkpoint manifest mismatch in ", BASE_CP)
  }
  baseline <- combine_units(lapply(file.path(BASE_CP, unit_files), readRDS))
  check_rows(baseline, manifest$design_ids)
  old_designs <- paste("Design", manifest$design_ids)
  for (nm in stats) {
    previous <- readRDS(files[[nm]])
    if ("Design" %in% names(previous)) {
      previous <- previous[previous$Design %in% old_designs, , drop = FALSE]
    }
    rownames(previous) <- NULL
    if (!identical(previous, baseline[[nm]])) {
      stop("November 24 checkpoint summaries differ from stored ", nm, " rows")
    }
  }

  prod <- dplyr::bind_rows(
    readRDS(sprintf("results/sim_data/sim_results_MLE_tau_sweep_combined_%s.rds", PROD_TAG)),
    readRDS(sprintf("results/sim_data/sim_results_MLEnonoracle_tau_sweep_combined_%s.rds",
                    PROD_TAG)))
  benchmarks <- paste("Design", c(3, 8, 9))
  prod <- prod[prod$Design %in% benchmarks, , drop = FALSE]
  mine <- baseline$scen[baseline$scen$Design %in% benchmarks, names(prod), drop = FALSE]
  prod <- prod[do.call(order, prod[keys]), ]
  mine <- mine[do.call(order, mine[keys]), ]
  rownames(prod) <- rownames(mine) <- NULL
  if (!identical(prod, mine)) stop("Designs 3, 8, 9 differ from September 27 production")
  reproduction <- c(
    "November 24 manifest matches. Reused Designs 3, 8, 9 and variants 101-130.",
    sprintf("November 24 exact agreement: %d scenario rows, %d surface rows, %d warning rows.",
            nrow(baseline$scen), nrow(baseline$surf), nrow(baseline$warn)),
    sprintf("September 27 exact agreement: %d scenario rows across all %d stored columns.",
            nrow(prod), ncol(prod)))
  cat(reproduction, sep = "\n")
}

# Boundary checkpoints have their own manifest. The original checkpoint files stay available.
manifest$design_ids <- design_ids
manifest$code_hashes <- c(manifest$code_hashes, vapply(
  c("code/29_igs_sensitivity.R", "code/29_igs_variant_designs.R"),
  digest::digest, character(1), file = TRUE, USE.NAMES = TRUE))
manifest$spill_types <- spill_types
dir.create(CP_DIR, recursive = TRUE, showWarnings = FALSE)
mf <- file.path(CP_DIR, "manifest.rds")
if (file.exists(mf)) {
  if (!identical(readRDS(mf), manifest)) stop("Boundary checkpoint manifest mismatch in ", CP_DIR)
} else {
  if (length(list.files(CP_DIR, "^unit_"))) stop("Checkpoints without a manifest in ", CP_DIR)
  saveRDS(manifest, mf)
}

# Identification checks use the actual assignments returned to run_unit() ----
#' Count rank deficiencies over the surface-major assignment draws of one design block
#' @param Z N x (K*J) assignment matrix from draw_assignments()
#' @param X N x K incidence surfaces in draw order
#' @param W N x N row-standardized spatial weights
#' @param spill_type "control_only" or "both"
#' @return One row with the draw count and both rank-deficiency counts
#' @examples
#' X <- matrix(seq_len(100), 100, 1)
#' Z <- matrix(rep(0:1, 50 * n_design_draw), 100, n_design_draw)
#' assignment_ranks(Z, X, diag(100), "both")
assignment_ranks <- function(Z, X, W, spill_type) {
  WZ <- W %*% Z
  surface <- rep(seq_len(ncol(X)), each = n_design_draw)
  stopifnot(length(surface) == ncol(Z), all(colSums(Z) == N / 2))
  rank_wz <- vapply(seq_len(ncol(Z)), function(j) {
    qr(cbind(1, Z[, j], WZ[, j], X[, surface[j]]))$rank
  }, integer(1))
  rank_fit <- rank_wz
  if (spill_type == "control_only") {
    spill <- (1 - Z) * WZ
    rank_fit <- vapply(seq_len(ncol(Z)), function(j) {
      qr(cbind(1, Z[, j], spill[, j], X[, surface[j]]))$rank
    }, integer(1))
  }
  data.frame(Draws = length(surface), Rank_WZ_Deficient = sum(rank_wz < 4L),
             Rank_Oracle_Deficient = sum(rank_fit < 4L))
}

surfaces <- lapply(inc_configs, generate_surfaces)
draw_assignments_prod <- draw_assignments
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(units)), function(i) {
  cp <- file.path(CP_DIR, unit_files[i])
  if (file.exists(cp)) return(readRDS(cp))
  u <- units[i, ]
  W <- get_active_spatial(grid_obj, u$nb_type)$W
  diagnostics <- list()
  # Give this worker's copy of run_unit() a local assignment wrapper. The
  # production function receives the same Z and uses the same RNG state.
  draw_assignments <- function(cfg, nb_type, rho, gamma_val, spill_type, d_id, X) {
    Z <- draw_assignments_prod(cfg, nb_type, rho, gamma_val, spill_type, d_id, X)
    diagnostics[[length(diagnostics) + 1L]] <<- cbind(
      data.frame(Neighbor_Type = nb_type, Spillover_Type = spill_type, Design = paste("Design", d_id)),
      assignment_ranks(Z, X, W, spill_type))
    Z
  }
  run_boundary_unit <- run_unit
  environment(run_boundary_unit) <- environment()
  out <- run_boundary_unit(u$cfg_index, u$nb_type, u$rho, surfaces[[u$cfg_index]])
  out$ranks <- dplyr::bind_rows(diagnostics)
  tmp <- paste0(cp, ".tmp")
  saveRDS(out, tmp)
  if (!file.rename(tmp, cp)) stop("Could not save checkpoint ", cp)
  cat(sprintf("Boundary unit %02d/%02d complete\n", i, nrow(units)))
  out
}, mc.cores = n_cores, mc.preschedule = FALSE)
bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
if (any(bad)) stop("Units failed: ", paste(which(bad), collapse = ", "))
cat(sprintf("Boundary profile took %.1f min (%d units, %d designs)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(units), length(design_ids)))

added <- combine_units(res)
check_rows(added, design_ids)
by_design <- c("Neighbor_Type", "Spillover_Type", "Design")
rank_rows <- dplyr::bind_rows(lapply(res, `[[`, "ranks"))
stopifnot(nrow(rank_rows) == nrow(units) * length(gamma_vals) * length(spill_types) *
            length(design_ids), all(rank_rows$Draws == n_surfaces * n_design_draw))
ranks <- rank_rows |>
  dplyr::group_by(dplyr::across(dplyr::all_of(by_design))) |>
  dplyr::summarise(Draws = sum(Draws), Rank_WZ_Deficient = sum(Rank_WZ_Deficient),
                   Rank_Oracle_Deficient = sum(Rank_Oracle_Deficient), .groups = "drop")
n_estimands <- length(true_tau_vals)
fits <- added$scen |>
  dplyr::group_by(dplyr::across(dplyr::all_of(c(by_design, "Estimator")))) |>
  dplyr::summarise(Valid = sum(N_Valid_Est), Aliased = sum(N_Aliased),
                   Warnings = sum(N_Warn), .groups = "drop")
identification <- merge(ranks, fits, by = by_design, all = TRUE)
identification$Fits <- identification$Draws * n_estimands
identification$Failed <- identification$Fits - identification$Valid
stopifnot(!anyNA(identification), nrow(identification) ==
            length(neighbor_types) * length(spill_types) * length(design_ids) * 2L,
          all(identification$Failed >= 0L),
          sum(identification$Warnings) == sum(added$warn[["Count"]]))
write.csv(identification, file.path(OUT_DIR, "igs_identification.csv"), row.names = FALSE)
cat(sprintf("Assignment draws checked: %d. Rank deficiencies: WZ %d, oracle %d.\n",
            sum(ranks$Draws), sum(ranks$Rank_WZ_Deficient), sum(ranks$Rank_Oracle_Deficient)))
cat(sprintf("Fit attempts: %d. Failed: %d. Aliased: %d. Other warnings: %d.\n",
            sum(identification$Fits), sum(identification$Failed), sum(identification$Aliased),
            sum(identification$Warnings)))
if (any(identification$Warnings > 0L)) {
  print(as.data.frame(added$warn), max = Inf)
}
if (any(identification[c("Rank_WZ_Deficient", "Rank_Oracle_Deficient", "Failed", "Aliased")] > 0L)) {
  stop("Identification or fit failures in new variants; see igs_identification.csv")
}

results <- added
if (profile == "full") {
  results <- setNames(lapply(stats, function(nm) dplyr::bind_rows(baseline[[nm]], added[[nm]])),
                      stats)
  check_rows(results, c(3, 8, 9, VARIANTS$id))
  for (nm in stats) {
    kept <- head(results[[nm]], nrow(baseline[[nm]]))
    if (!identical(kept, baseline[[nm]])) stop("Combining changed retained ", nm, " rows")
  }
}
for (nm in stats) saveRDS(results[[nm]], files[[nm]])
write.csv(VARIANTS, file.path(OUT_DIR, "igs_variants.csv"), row.names = FALSE)
if (profile == "full") writeLines(reproduction, file.path(OUT_DIR, "reproduction_check.txt"))
cat("Saved to", OUT_DIR, "\n")
