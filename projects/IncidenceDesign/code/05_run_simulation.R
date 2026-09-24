# ==============================================================================
# 05_run_simulation.R
# Main orchestrator (2026-09 revision): sources modules 01-04, defines parameters,
# runs the simulation over (config x neighbor type x rho) work units, and saves
# oracle and non-oracle results.
#
# Specification: docs/plans/simulation-revision-spec.md (DGP, seeds, estimand,
# dependence structure, result schema). Fix IDs M1-M8 refer to
# docs/plans/simulation-revision-plan.md.
#
# Usage (from code/), BLAS pinned to one thread (the runner refuses otherwise):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript 05_run_simulation.R pilot     # B4 pilot -> results/pilot_rev_2026-09/
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript 05_run_simulation.R full      # full grid -> results/sim_data/
#
# Parallelization: mclapply (fork) over the 40 (config x nb x rho) units, one unit
# per core at a time (mc.preschedule = FALSE). Every random draw is seeded from a
# key (spec section 4), so results do not depend on scheduling. Memory per worker is
# small (a few N x 250 matrices, N = 100). Not written for SLURM/Longleaf.
#
# DIM is not run by this script (out of scope for the revision; the legacy DIM
# path is estimate_tau() in 04).
# ==============================================================================

# --- Load Libraries ---
library(spdep)
library(spatialreg)
library(dplyr)
library(digest)
library(parallel)

# --- Source Module Scripts ---
# Detect script directory: works via Rscript --file=, source(), or interactive.
# A caller (the tests) can preset sim_script_dir instead.
script_dir <- if (exists("sim_script_dir")) sim_script_dir else local({
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("--file=", "", file_arg[1])))
  } else {
    tryCatch(
      normalizePath(dirname(sys.frame(1)$ofile)),
      error = function(e) normalizePath(getwd())
    )
  }
})

code_files <- file.path(script_dir, c("01_spatial_setup.R", "02_incidence_generation.R",
                                      "03_designs.R", "04_estimation.R",
                                      "05_run_simulation.R"))
for (f in code_files[1:4]) source(f)

# ==============================================================================
# CONFIGURABLE PARAMETERS
# ==============================================================================

# Profile: "full" (default) or "pilot" (first trailing argument). A variable named
# sim_profile set before source() takes precedence (used by the tests).
if (!exists("sim_profile")) {
  trailing <- commandArgs(trailingOnly = TRUE)
  sim_profile <- if (length(trailing) > 0) trailing[1] else "full"
}
stopifnot(sim_profile %in% c("full", "pilot", "test"))

# Parallelization: number of forked workers (1 = sequential)
if (!exists("n_cores")) n_cores <- 10

# Estimation engine: "lean" (fit_sar_lag, validated) or "lagsarlm"
if (!exists("engine")) engine <- "lean"

# Grid
grid_dim <- 10                 # 10x10 = 100 clusters

# Incidence generation
incidence_modes    <- c("iid", "spatial", "poisson")
rho_incidence_vals <- c(0.20, 0.50)  # for spatial & poisson modes

# Poisson-specific parameters (M2: 100,000 per cluster -> ~35 expected deaths)
base_rate       <- 35 / 100000   # SUD base rate (Mirzaei et al.)
pop_per_cluster <- 100000
pop_mode        <- "equal"       # "equal" or "heterogeneous"

# Spatial / spillover parameters
rho_vals       <- c(0.00, 0.01, 0.20, 0.50)
gamma_vals     <- c(0.5, 0.6, 0.7, 0.8)
spill_types    <- c("control_only", "both")
neighbor_types <- c("rook", "queen")

# Outcome parameters: tau swept; beta and sigma fixed
true_tau_vals <- c(0.8, 1.0, 1.5, 2.0, 3.0)
beta          <- 1.0
sigma         <- 1.0

# Monte Carlo structure (M1): K incidence surfaces x J design draws per surface;
# every one of the K*J fits gets its own noise column.
n_surfaces    <- 10    # K
n_design_draw <- 25    # J

# Designs to evaluate
design_ids <- 1:8

# Pilot profile (plan B4): tau = 1; all configs, nb and regimes; rho in {0, 0.5};
# gamma in {0.5, 0.8}
if (sim_profile == "pilot") {
  true_tau_vals <- 1.0
  rho_vals      <- c(0.00, 0.50)
  gamma_vals    <- c(0.5, 0.8)
}
# Test profile (code/tests): tiny grid slice for determinism checks
if (sim_profile == "test") {
  true_tau_vals <- c(1.0, 2.0)
  rho_vals      <- 0.20
  gamma_vals    <- 0.6
  spill_types   <- "both"
  n_surfaces    <- 2
  n_design_draw <- 3
  design_ids    <- c(1, 2, 6)
}

# ==============================================================================
# BUILD SPATIAL GRID, SETUPS, INCIDENCE CONFIGS
# ==============================================================================

grid_obj <- build_spatial_grid(grid_dim)
N        <- grid_obj$N_clusters
coords   <- grid_obj$coords
I_mat    <- diag(N)
lag_setups <- list(rook  = sar_lag_setup(grid_obj$listw_rook),
                   queen = sar_lag_setup(grid_obj$listw_queen))

inc_configs <- list()
for (inc_mode in incidence_modes) {
  rho_sweep <- if (inc_mode == "iid") 0 else rho_incidence_vals
  for (rho_x in rho_sweep) {
    inc_configs[[length(inc_configs) + 1]] <- list(mode = inc_mode, rho_x = rho_x)
  }
}

# ==============================================================================
# SEEDING (M4; spec section 4)
# ==============================================================================

#' Build the canonical seed key string for one random draw
#'
#' Numeric fields are formatted with sprintf("%.2f") so that, e.g., rho = 0.5 and
#' rho = 0.50 give the same key; fields are joined with "|".
#'
#' @param ... Key fields (character or numeric)
#' @return Character scalar
#' @examples
#' seed_key("eps", "spatial", 0.2, "queen", 0.5, 0.8, "both")
seed_key <- function(...) {
  fields <- lapply(list(...), function(v) if (is.numeric(v)) sprintf("%.2f", v) else v)
  paste(unlist(fields), collapse = "|")
}

#' Seed the RNG from a key, with RNG kinds fixed explicitly
#'
#' @param ... Key fields passed to seed_key()
#' @return Invisibly, the integer seed
set_seed_key <- function(...) {
  s <- digest::digest2int(seed_key(...))
  set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  invisible(s)
}

#' Generate the K incidence surfaces for one config, each from its own seed
#'
#' Surface k is keyed by ("X", mode, rho_X, k) only, so it is shared by every block
#' of the config (Invariant 1) and is identical in every run.
#'
#' @param cfg List(mode, rho_x)
#' @return N x K matrix, column k = X_k in [0, 1]
generate_surfaces <- function(cfg) {
  sapply(seq_len(n_surfaces), function(k) {
    set_seed_key("X", cfg$mode, cfg$rho_x, k)
    generate_incidence(cfg$mode, N, 1, W = grid_obj$W_queen,
                       rho_incidence = cfg$rho_x, base_rate = base_rate,
                       pop_per_cluster = pop_per_cluster, pop_mode = pop_mode)[, 1]
  })
}

# ==============================================================================
# SUMMARIES (spec section 7)
# ==============================================================================

#' Summarize one estimator's fits for one scenario
#'
#' Point estimates pool all valid fits (finite tau-hat and SE). Monte Carlo SEs come
#' from the K surface-level means: SE(m) = sd(m_1..m_K) / sqrt(K) (the surface is the
#' primary sampling unit; spec section 5). Intervals downstream use t_{K-1}.
#'
#' @param tau_hat,se Numeric vectors length K*J (fit order: surface-major)
#' @param surface Integer vector length K*J, surface index of each fit
#' @param n_treated Integer vector length K*J, sum(Z) of each fit
#' @param true_tau Numeric scalar
#' @return List with scen (one-row data frame of metrics) and surf (K-row data frame)
summarize_fits <- function(tau_hat, se, surface, n_treated, true_tau) {
  z975  <- qnorm(0.975)
  valid <- is.finite(tau_hat) & is.finite(se)
  lo    <- tau_hat - z975 * se
  hi    <- tau_hat + z975 * se
  err   <- tau_hat - true_tau
  covered  <- lo <= true_tau & hi >= true_tau
  rejected <- lo > 0                              # one-sided, as before

  surf <- do.call(rbind, lapply(sort(unique(surface)), function(k) {
    v <- valid & surface == k
    data.frame(Surface = k, N_Valid = sum(v),
               Bias = if (any(v)) mean(err[v]) else NA_real_,
               MSE = if (any(v)) mean(err[v]^2) else NA_real_,
               Coverage = if (any(v)) mean(covered[v]) else NA_real_,
               Power = if (any(v)) mean(rejected[v]) else NA_real_,
               Mean_Treated = mean(n_treated[surface == k]))
  }))
  ok_k <- surf$N_Valid > 0
  mc_se <- function(m) sd(m[ok_k]) / sqrt(sum(ok_k))

  scen <- data.frame(
    Mean_Estimate = if (any(valid)) mean(tau_hat[valid]) else NA_real_,
    Bias          = if (any(valid)) mean(err[valid]) else NA_real_,
    SD            = if (sum(valid) > 1) sd(tau_hat[valid]) else NA_real_,
    MSE           = if (any(valid)) mean(err[valid]^2) else NA_real_,
    Coverage      = if (any(valid)) mean(covered[valid]) else NA_real_,
    Fail_Rate     = 1 - mean(valid),
    N_Valid_Est   = sum(valid),
    Power         = if (any(valid)) mean(rejected[valid]) else NA_real_,
    SE_Bias       = mc_se(surf$Bias),
    SE_MSE        = mc_se(surf$MSE),
    SE_Coverage   = mc_se(surf$Coverage),
    SE_Power      = mc_se(surf$Power),
    N_Surfaces    = sum(ok_k),
    Mean_Treated  = mean(n_treated)
  )
  list(scen = scen, surf = surf)
}

#' Draw all K*J assignments of one design in one block (M1 + M4)
#'
#' Columns (k-1)*J + 1..k*J are J draws of design d given surface X_k, keyed by
#' ("Z", mode, rho_X, nb, rho, gamma, regime, k, d).
#'
#' @param cfg List(mode, rho_x)
#' @param nb_type,rho,gamma_val,spill_type Block keys
#' @param d_id Design ID
#' @param X N x K matrix of surfaces
#' @return N x (K*J) binary matrix, surface-major column order
draw_assignments <- function(cfg, nb_type, rho, gamma_val, spill_type, d_id, X) {
  nb_list <- get_active_spatial(grid_obj, nb_type)$nb
  do.call(cbind, lapply(seq_len(ncol(X)), function(k) {
    set_seed_key("Z", cfg$mode, cfg$rho_x, nb_type, rho, gamma_val, spill_type, k, d_id)
    get_designs(d_id, n_design_draw, N, X[, k], nb_list, coords)
  }))
}

# ==============================================================================
# ONE WORK UNIT: (config, nb type, rho)
# ==============================================================================

#' Run every (gamma, regime, design, tau) scenario of one (config, nb, rho) unit
#'
#' For each block (gamma, regime): draw eps (N x K*J) from the block key; for each
#' design d and surface k draw J assignments from X_k (M1); build
#' Y_kj = (I - rho W)^{-1} (tau Z_kj + S(Z_kj) + beta X_k + eps_kj) for each tau; fit
#' oracle and non-oracle models to the same Y (M8). Nothing here draws random
#' numbers except under set_seed_key(), so the unit's output is fixed by its keys.
#'
#' @param cfg_index Integer index into inc_configs
#' @param nb_type "rook" or "queen"
#' @param rho Outcome spatial autocorrelation
#' @param X N x K matrix of surfaces for this config
#' @return List with scen, surf, warn data frames (both estimators stacked, column
#'   Estimator)
run_unit <- function(cfg_index, nb_type, rho, X) {
  cfg     <- inc_configs[[cfg_index]]
  spatial <- get_active_spatial(grid_obj, nb_type)
  W       <- spatial$W
  setup   <- lag_setups[[nb_type]]
  inv_mat <- solve(I_mat - rho * W)                 # (I - rho W)^{-1}
  KJ      <- n_surfaces * n_design_draw
  surface <- rep(seq_len(n_surfaces), each = n_design_draw)   # fit f -> surface k

  scen_rows <- list(); surf_rows <- list(); warn_rows <- list()

  for (gamma_val in gamma_vals) {
    for (spill_type in spill_types) {
      # eps: one column per fit, shared by all designs, taus and both estimators
      set_seed_key("eps", cfg$mode, cfg$rho_x, nb_type, rho, gamma_val, spill_type)
      E <- matrix(rnorm(N * KJ, 0, sigma), nrow = N, ncol = KJ)

      for (d_id in design_ids) {
        Zm <- draw_assignments(cfg, nb_type, rho, gamma_val, spill_type, d_id, X)
        WZ <- W %*% Zm
        Sm <- if (spill_type == "control_only") gamma_val * WZ * (1 - Zm) else gamma_val * WZ
        Xm <- X[, surface]                                  # X_k for each fit
        n_treated <- colSums(Zm)
        # Design-level identifiability flag: rank([1, Z, WZ]) < 3 for any draw
        rank_def <- any(vapply(seq_len(KJ), function(f) {
          qr(cbind(1, Zm[, f], WZ[, f]))$rank < 3
        }, logical(1)))

        base <- inv_mat %*% (Sm + beta * Xm + E)            # tau-free part of Y
        IZ   <- inv_mat %*% Zm

        for (true_tau in true_tau_vals) {
          Y <- true_tau * IZ + base
          fits <- lapply(seq_len(KJ), function(f) {
            fit_tau_models(Y[, f], Zm[, f], Sm[, f], Xm[, f], setup,
                           listw = spatial$listw, engine = engine)
          })
          for (est in c("oracle", "nonoracle")) {
            tau_hat <- vapply(fits, function(x) x[[est]]$tau, numeric(1))
            se_hat  <- vapply(fits, function(x) x[[est]]$se, numeric(1))
            warns   <- lapply(fits, function(x) x[[est]]$warns)
            is_alias <- vapply(warns, function(w) any(startsWith(w, "Aliased variables found")),
                               logical(1))
            other   <- lapply(warns, function(w) w[!startsWith(w, "Aliased variables found")])

            keys <- data.frame(
              Incidence_Mode = cfg$mode, Rho_Incidence = cfg$rho_x,
              Neighbor_Type = nb_type, Design = paste("Design", d_id),
              Rho = rho, Gamma = gamma_val, Spillover_Type = spill_type,
              True_Tau = true_tau, Estimator = est, stringsAsFactors = FALSE)
            sm <- summarize_fits(tau_hat, se_hat, surface, n_treated, true_tau)
            scen_rows[[length(scen_rows) + 1]] <- cbind(
              keys, sm$scen,
              N_Aliased = sum(is_alias),
              N_Warn = sum(lengths(other)),
              Z_WZ_rank_deficient = rank_def)
            surf_rows[[length(surf_rows) + 1]] <- cbind(keys, sm$surf)
            msgs <- unlist(other)
            if (length(msgs) > 0) {
              tab <- table(msgs)
              warn_rows[[length(warn_rows) + 1]] <- cbind(
                keys, Message = names(tab), Count = as.integer(tab))
            }
          }
        }
      }
    }
  }
  list(scen = bind_rows(scen_rows), surf = bind_rows(surf_rows),
       warn = bind_rows(warn_rows))
}

#' Run the simulation for the current profile: BLAS guard, checkpoints/manifest,
#' parallel execution over units, integrity stop, and saving (spec sections 8-9).
#'
#' @return Invisibly, the output directory
main <- function() {
  # ==============================================================================
  # BLAS GUARD (one BLAS thread per worker; must be set before R starts)
  # ==============================================================================

  blas_path <- sessionInfo()$BLAS
  blas_vars <- c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
  if (!all(Sys.getenv(blas_vars) == "1")) {
    stop("Set ", paste0(blas_vars, "=1", collapse = " "),
         " before starting R (BLAS: ", blas_path, ")")
  }

  # ==============================================================================
  # CHECKPOINTS + MANIFEST (refuse to load on any mismatch)
  # ==============================================================================

  results_dir <- file.path(dirname(script_dir), "results")
  checkpoint_dir <- file.path(results_dir, "checkpoints", "rev_2026-09", sim_profile)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  params <- list(
    sim_profile = sim_profile, engine = engine, grid_dim = grid_dim,
    inc_configs = inc_configs, base_rate = base_rate, pop_per_cluster = pop_per_cluster,
    pop_mode = pop_mode, rho_vals = rho_vals, gamma_vals = gamma_vals,
    spill_types = spill_types, neighbor_types = neighbor_types,
    true_tau_vals = true_tau_vals, beta = beta, sigma = sigma,
    n_surfaces = n_surfaces, n_design_draw = n_design_draw, design_ids = design_ids)
  manifest <- list(
    param_hash = digest::digest(params),
    code_hashes = vapply(code_files, digest::digest, character(1), file = TRUE),
    packages = vapply(c("spatialreg", "spdep", "dplyr", "digest"),
                      function(p) as.character(packageVersion(p)), character(1)),
    R = R.version.string, blas = blas_path)
  names(manifest$code_hashes) <- basename(code_files)

  manifest_file <- file.path(checkpoint_dir, "manifest.rds")
  existing_units <- list.files(checkpoint_dir, pattern = "^unit_.*\\.rds$")
  if (file.exists(manifest_file)) {
    old <- readRDS(manifest_file)
    if (!identical(old, manifest)) {
      stop("Checkpoint manifest mismatch in ", checkpoint_dir,
           " (parameters, code, packages or BLAS changed). Refusing to load; ",
           "move or delete that directory to start fresh.")
    }
  } else if (length(existing_units) > 0) {
    stop("Checkpoints without a manifest in ", checkpoint_dir, "; refusing to load.")
  } else {
    saveRDS(manifest, manifest_file)
  }

  units <- expand.grid(cfg_index = seq_along(inc_configs), nb_type = neighbor_types,
                       rho = rho_vals, stringsAsFactors = FALSE)
  unit_file <- function(u) {
    file.path(checkpoint_dir, sprintf("unit_%s_rhoX%.2f_%s_rho%.2f.rds",
              inc_configs[[u$cfg_index]]$mode, inc_configs[[u$cfg_index]]$rho_x,
              u$nb_type, u$rho))
  }

  # ==============================================================================
  # EXECUTE
  # ==============================================================================

  n_scen <- nrow(units) * length(gamma_vals) * length(spill_types) *
            length(design_ids) * length(true_tau_vals)
  cat(sprintf("\n=== Simulation (profile: %s, engine: %s) ===\n", sim_profile, engine))
  cat(sprintf("Units: %d | scenarios per estimator: %d | fits per scenario: %d | total fits: %d\n",
              nrow(units), n_scen, n_surfaces * n_design_draw,
              2 * n_scen * n_surfaces * n_design_draw))
  cat(sprintf("Workers: %d | checkpoints: %s\n\n", n_cores, checkpoint_dir))

  global_start <- Sys.time()
  surfaces <- lapply(inc_configs, generate_surfaces)     # X per config, keyed per k

  run_or_load <- function(i) {
    u  <- units[i, ]
    cp <- unit_file(u)
    if (file.exists(cp)) return(readRDS(cp))
    t0  <- Sys.time()
    out <- run_unit(u$cfg_index, u$nb_type, u$rho, surfaces[[u$cfg_index]])
    out$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    saveRDS(out, cp)
    out
  }

  unit_results <- if (n_cores > 1) {
    mclapply(seq_len(nrow(units)), run_or_load, mc.cores = n_cores,
             mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(units)), run_or_load)
  }

  # A NULL or try-error result stops the run (never silently dropped)
  bad <- vapply(unit_results, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  if (any(bad)) {
    msgs <- vapply(which(bad), function(i) {
      sprintf("unit %d (%s): %s", i, basename(unit_file(units[i, ])),
              if (is.null(unit_results[[i]])) "NULL" else as.character(unit_results[[i]]))
    }, character(1))
    stop("Simulation units failed:\n", paste(msgs, collapse = "\n"))
  }

  scen_all <- bind_rows(lapply(unit_results, `[[`, "scen"))
  surf_all <- bind_rows(lapply(unit_results, `[[`, "surf"))
  warn_all <- bind_rows(lapply(unit_results, `[[`, "warn"))
  unit_secs <- vapply(unit_results, function(r) r$elapsed_sec %||% NA_real_, numeric(1))

  total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
  cat(sprintf("=== Complete: %.1f min wall; unit compute %.1f min total (%.0f s mean/unit) ===\n",
              total_elapsed, sum(unit_secs, na.rm = TRUE) / 60, mean(unit_secs, na.rm = TRUE)))

  # ==============================================================================
  # SAVE (spec section 8)
  # ==============================================================================

  out_dir <- if (sim_profile == "full") file.path(results_dir, "sim_data") else
    file.path(results_dir, paste0(sim_profile, "_rev_2026-09"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  est_tags <- c(oracle = "MLE_tau_sweep", nonoracle = "MLEnonoracle_tau_sweep")

  for (est in names(est_tags)) {
    tag  <- est_tags[[est]]
    scen <- scen_all[scen_all$Estimator == est, ]
    saveRDS(scen, file.path(out_dir, sprintf("sim_results_%s_combined_%s.rds", tag, timestamp_str)))
    for (mode_name in unique(scen$Incidence_Mode)) {
      saveRDS(scen[scen$Incidence_Mode == mode_name, ],
              file.path(out_dir, sprintf("sim_results_%s_%s_%s.rds", tag, mode_name, timestamp_str)))
    }
    saveRDS(surf_all[surf_all$Estimator == est, ],
            file.path(out_dir, sprintf("surface_results_%s_%s.rds", tag, timestamp_str)))
    write.csv(warn_all[if (nrow(warn_all)) warn_all$Estimator == est else 0, , drop = FALSE],
              file.path(out_dir, sprintf("warnings_%s_%s.csv", tag, timestamp_str)),
              row.names = FALSE)
    cat(sprintf("Saved %s: %d scenario rows, %d surface rows, %d warning rows\n",
                tag, nrow(scen), sum(surf_all$Estimator == est),
                if (nrow(warn_all)) sum(warn_all$Estimator == est) else 0L))
  }
  saveRDS(list(manifest = manifest, params = params, sessionInfo = sessionInfo(),
               unit_elapsed_sec = unit_secs, wall_min = total_elapsed),
          file.path(out_dir, sprintf("run_info_%s.rds", timestamp_str)))
  invisible(out_dir)
}

# Tests source this file with sim_define_only <- TRUE to get the functions only
if (!isTRUE(get0("sim_define_only", ifnotfound = FALSE))) main()
