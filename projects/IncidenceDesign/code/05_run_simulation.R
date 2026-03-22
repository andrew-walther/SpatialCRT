# ==============================================================================
# 05_run_simulation.R
# Main orchestrator: sources modules 01-04, defines parameters, runs the
# nested simulation loop, saves results.
#
# Usage:
#   source("05_run_simulation.R")            # sequential
#   Rscript 05_run_simulation.R              # sequential
#   Set n_cores > 1 before sourcing for parallel execution across
#   incidence configs (iid, spatial, poisson).
# ==============================================================================

# --- Load Libraries ---
library(spdep)
library(spatialreg)
library(dplyr)
library(digest)
library(parallel)

# --- Source Module Scripts ---
# Detect script directory: works via Rscript, source(), or interactive
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    if (exists("ofile", envir = sys.frame(1))) dirname(sys.frame(1)$ofile)
    else getwd()
  }
)

source(file.path(script_dir, "01_spatial_setup.R"))
source(file.path(script_dir, "02_incidence_generation.R"))
source(file.path(script_dir, "03_designs.R"))
source(file.path(script_dir, "04_estimation.R"))

# ==============================================================================
# CONFIGURABLE PARAMETERS
# ==============================================================================

# Estimation mode
estimation_mode <- "MLE"      # "MLE" or "DIM"

# Parallelization: number of cores for parallel incidence config processing
# Set to 1 for sequential, or detectCores() - 1 for max parallelism
n_cores <- 1

# Grid
grid_dim <- 10                 # 10x10 = 100 clusters

# Incidence generation
incidence_modes    <- c("iid", "spatial", "poisson")
rho_incidence_vals <- c(0.20, 0.50)  # for spatial & poisson modes

# Poisson-specific parameters
base_rate       <- 35 / 100000   # SUD base rate (Mirzaei et al.)
pop_per_cluster <- 1000
pop_mode        <- "equal"       # "equal" or "heterogeneous"

# Spatial / spillover parameters
rho_vals       <- c(0.00, 0.01, 0.20, 0.50)
gamma_vals     <- c(0.5, 0.6, 0.7, 0.8)
spill_types    <- c("control_only", "both")
neighbor_types <- c("rook", "queen")

# SDM true parameters
true_tau <- 1.0
beta     <- 1.0
sigma    <- 1.0

# Resample counts (mode-dependent)
if (estimation_mode == "MLE") {
  n_designs_resamples <- 25
  n_outcome_resamples <- 10
} else {
  n_designs_resamples <- 25
  n_outcome_resamples <- 100
}

# Designs to evaluate
design_ids <- c(1, 2, 3, 4, 5, 6, 7, 8)

# Oracle spillover toggle
include_spill_covariate <- TRUE

# ==============================================================================
# BUILD SPATIAL GRID
# ==============================================================================

cat("Building spatial grid...\n")
grid_obj <- build_spatial_grid(grid_dim)
N <- grid_obj$N_clusters
coords <- grid_obj$coords
I_mat <- diag(N)

# ==============================================================================
# BUILD INCIDENCE CONFIGS
# ==============================================================================

inc_configs <- list()
for (inc_mode in incidence_modes) {
  rho_sweep <- if (inc_mode == "iid") 0 else rho_incidence_vals
  for (rho_x in rho_sweep) {
    inc_configs[[length(inc_configs) + 1]] <- list(mode = inc_mode, rho_x = rho_x)
  }
}

scenarios_per_config <- length(neighbor_types) * length(rho_vals) *
                        length(gamma_vals) * length(spill_types) * length(design_ids)
total_scenarios <- length(inc_configs) * scenarios_per_config

cat(sprintf("\n=== Simulation Configuration ===\n"))
cat(sprintf("Estimation mode:      %s\n", estimation_mode))
cat(sprintf("Incidence configs:    %d\n", length(inc_configs)))
cat(sprintf("Scenarios per config: %d\n", scenarios_per_config))
cat(sprintf("Total scenarios:      %d\n", total_scenarios))
cat(sprintf("Iterations/scenario:  %d (%d design x %d outcome)\n",
            n_designs_resamples * n_outcome_resamples,
            n_designs_resamples, n_outcome_resamples))
cat(sprintf("Parallel cores:       %d\n", n_cores))
cat(sprintf("================================\n\n"))

# ==============================================================================
# CORE SIMULATION FUNCTION (processes one incidence config)
# ==============================================================================

run_incidence_config <- function(ic_index, ic_config, grid_obj, coords, I_mat,
                                  verbose = TRUE) {
  inc_mode <- ic_config$mode
  rho_x    <- ic_config$rho_x
  inc_label <- if (inc_mode == "iid") "iid" else sprintf("%s(rho_X=%.2f)", inc_mode, rho_x)

  if (verbose) cat(sprintf("\n--- [Config %d] %s ---\n", ic_index, inc_label))

  N <- grid_obj$N_clusters

  # Generate incidence ONCE for this config
  X_matrix <- generate_incidence(
    mode          = inc_mode,
    N             = N,
    n_resamples   = n_outcome_resamples,
    W             = grid_obj$W_queen,
    rho_incidence = rho_x,
    base_rate     = base_rate,
    pop_per_cluster = pop_per_cluster,
    pop_mode      = pop_mode
  )
  base_incidence <- X_matrix[, 1]

  if (verbose) {
    cat(sprintf("  X range: [%.4f, %.4f], mean: %.4f\n",
                min(X_matrix), max(X_matrix), mean(X_matrix)))
  }

  config_results <- list()
  config_idx <- 1
  config_times <- numeric(0)

  for (nb_type in neighbor_types) {
    spatial <- get_active_spatial(grid_obj, nb_type)
    W <- spatial$W
    active_listw <- spatial$listw
    active_nb <- spatial$nb

    for (rho in rho_vals) {
      inv_mat <- solve(I_mat - rho * W)
      Epsilon_matrix <- matrix(rnorm(N * n_outcome_resamples, 0, sigma),
                               nrow = N, ncol = n_outcome_resamples)
      baseline_part <- inv_mat %*% (X_matrix * beta + Epsilon_matrix)

      for (gamma_val in gamma_vals) {
        for (spill_type in spill_types) {
          for (d_id in design_ids) {

            scenario_start <- Sys.time()

            # Per-scenario deterministic seed
            seed_string <- paste(inc_mode, rho_x, nb_type, rho,
                                 gamma_val, spill_type, d_id, sep = "|")
            scenario_seed <- digest2int(seed_string)
            set.seed(scenario_seed)

            # Progress logging with ETA
            if (verbose) {
              avg_time <- if (length(config_times) > 0) mean(config_times) else NA
              eta_str <- if (!is.na(avg_time)) {
                remaining <- avg_time * (scenarios_per_config - config_idx + 1)
                if (remaining > 3600) sprintf("ETA: %.1f hrs", remaining / 3600)
                else if (remaining > 60) sprintf("ETA: %.1f min", remaining / 60)
                else sprintf("ETA: %.0f sec", remaining)
              } else {
                "ETA: calculating..."
              }
              cat(sprintf("  [%3d/%d] NB:%-5s | rho:%.2f | gam:%.2f | spill:%-12s | D%d | %s\n",
                          config_idx, scenarios_per_config,
                          nb_type, rho, gamma_val, spill_type, d_id, eta_str))
              flush.console()
            }

            # Generate treatment assignments
            if (is_design_deterministic(d_id)) {
              Z_single <- get_designs(d_id, 1, N, base_incidence, active_nb, coords)[, 1]
              Z_matrix <- matrix(rep(Z_single, n_designs_resamples),
                                 nrow = N, ncol = n_designs_resamples)
            } else {
              Z_matrix <- get_designs(d_id, n_designs_resamples, N,
                                      base_incidence, active_nb, coords)
            }

            # Collect estimates
            all_estimates <- numeric(0)
            all_ci_lower  <- numeric(0)
            all_ci_upper  <- numeric(0)

            for (d_resample in seq_len(n_designs_resamples)) {
              Z <- Z_matrix[, d_resample]
              WZ <- as.vector(W %*% Z)
              spill_term <- if (spill_type == "control_only") {
                gamma_val * WZ * (1 - Z)
              } else {
                gamma_val * WZ
              }

              trt_effect <- as.vector(inv_mat %*% (true_tau * Z + spill_term))
              Y_sim <- sweep(baseline_part, 1, trt_effect, "+")

              est_result <- estimate_tau(
                estimation_mode     = estimation_mode,
                Y_sim               = Y_sim,
                Z                   = Z,
                spill_term          = spill_term,
                X_matrix            = X_matrix,
                active_listw        = active_listw,
                n_outcome_resamples = n_outcome_resamples,
                include_spill_covariate = include_spill_covariate
              )

              all_estimates <- c(all_estimates, est_result$estimates)
              all_ci_lower  <- c(all_ci_lower,  est_result$ci_lower)
              all_ci_upper  <- c(all_ci_upper,  est_result$ci_upper)
            }

            # Aggregate
            valid_idx <- !is.na(all_estimates)
            valid_est <- all_estimates[valid_idx]
            fail_rate <- 1 - sum(valid_idx) / length(all_estimates)

            if (length(valid_est) > 0) {
              valid_ci <- valid_idx & !is.na(all_ci_lower) & !is.na(all_ci_upper)
              coverage_val <- if (sum(valid_ci) > 0) {
                mean((all_ci_lower[valid_ci] <= true_tau) &
                     (all_ci_upper[valid_ci] >= true_tau))
              } else NA

              config_results[[config_idx]] <- data.frame(
                Incidence_Mode  = inc_mode,
                Rho_Incidence   = rho_x,
                Neighbor_Type   = nb_type,
                Design          = paste("Design", d_id),
                Rho             = rho,
                Gamma           = gamma_val,
                Spillover_Type  = spill_type,
                Mean_Estimate   = mean(valid_est),
                Bias            = mean(valid_est) - true_tau,
                SD              = sd(valid_est),
                MSE             = mean((valid_est - true_tau)^2),
                Coverage        = coverage_val,
                Fail_Rate       = fail_rate,
                stringsAsFactors = FALSE
              )
            } else {
              config_results[[config_idx]] <- data.frame(
                Incidence_Mode  = inc_mode,
                Rho_Incidence   = rho_x,
                Neighbor_Type   = nb_type,
                Design          = paste("Design", d_id),
                Rho             = rho,
                Gamma           = gamma_val,
                Spillover_Type  = spill_type,
                Mean_Estimate   = NA, Bias = NA, SD = NA, MSE = NA,
                Coverage        = NA, Fail_Rate = fail_rate,
                stringsAsFactors = FALSE
              )
            }

            scenario_elapsed <- as.numeric(difftime(Sys.time(), scenario_start, units = "secs"))
            config_times <- c(config_times, scenario_elapsed)

            if (verbose) {
              mins <- floor(scenario_elapsed / 60)
              secs <- round(scenario_elapsed %% 60, 1)
              cat(sprintf("     -> %d min %.1f sec\n", mins, secs))
              flush.console()
            }

            config_idx <- config_idx + 1
          }
        }
      }
    }
  }

  config_df <- bind_rows(config_results)
  total_time <- sum(config_times)
  if (verbose) {
    cat(sprintf("  Config %s complete: %d rows, %.1f min total\n",
                inc_label, nrow(config_df), total_time / 60))
  }
  config_df
}

# ==============================================================================
# CHECKPOINTING SETUP
# ==============================================================================

results_dir <- file.path(dirname(script_dir), "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
sim_data_dir <- file.path(results_dir, "sim_data")
if (!dir.exists(sim_data_dir)) dir.create(sim_data_dir, recursive = TRUE)
checkpoint_dir <- file.path(results_dir, "checkpoints")
if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

# Helper: canonical checkpoint filename for a config
config_checkpoint_file <- function(ic_config) {
  label <- paste0(ic_config$mode, "_rhoX", sprintf("%.2f", ic_config$rho_x))
  file.path(checkpoint_dir,
            sprintf("checkpoint_%s_%s.rds", estimation_mode, label))
}

# ==============================================================================
# EXECUTE SIMULATION (sequential or parallel)
# ==============================================================================

global_start <- Sys.time()

if (n_cores > 1 && length(inc_configs) > 1) {
  cat(sprintf("Running %d incidence configs in parallel across %d cores...\n",
              length(inc_configs), min(n_cores, length(inc_configs))))

  results_by_config <- mclapply(
    seq_along(inc_configs),
    function(i) {
      cp_file <- config_checkpoint_file(inc_configs[[i]])
      if (file.exists(cp_file)) {
        cat(sprintf("  Config %d: loading checkpoint (skipping)\n", i))
        return(readRDS(cp_file))
      }
      result <- run_incidence_config(i, inc_configs[[i]], grid_obj, coords, I_mat,
                                     verbose = FALSE)
      saveRDS(result, cp_file)
      result
    },
    mc.cores = min(n_cores, length(inc_configs))
  )

  # Check for errors
  errors <- sapply(results_by_config, inherits, "try-error")
  if (any(errors)) {
    cat("WARNING: Some configs failed:\n")
    for (i in which(errors)) {
      cat(sprintf("  Config %d: %s\n", i, as.character(results_by_config[[i]])))
    }
    results_by_config <- results_by_config[!errors]
  }
} else {
  cat("Running sequentially...\n")
  results_by_config <- vector("list", length(inc_configs))
  for (i in seq_along(inc_configs)) {
    cp_file <- config_checkpoint_file(inc_configs[[i]])
    if (file.exists(cp_file)) {
      cat(sprintf("\n--- [Config %d] checkpoint found, loading and skipping ---\n", i))
      results_by_config[[i]] <- readRDS(cp_file)
    } else {
      results_by_config[[i]] <- run_incidence_config(
        i, inc_configs[[i]], grid_obj, coords, I_mat, verbose = TRUE)
      saveRDS(results_by_config[[i]], cp_file)
      cat(sprintf("  Checkpoint saved: %s\n", basename(cp_file)))
    }
  }
}

# ==============================================================================
# COMBINE & SAVE RESULTS
# ==============================================================================

granular_results <- bind_rows(results_by_config)

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\n=== Simulation Complete ===\n"))
cat(sprintf("Total time: %.1f minutes\n", total_elapsed))
cat(sprintf("Scenarios completed: %d\n", nrow(granular_results)))

# --- Save combined results ---
timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

combined_file <- file.path(sim_data_dir,
                            sprintf("sim_results_%s_combined_%s.rds",
                                    estimation_mode, timestamp_str))
saveRDS(granular_results, combined_file)
cat(sprintf("Combined results saved: %s\n", basename(combined_file)))

# --- Save split results by incidence mode ---
inc_mode_groups <- split(granular_results, granular_results$Incidence_Mode)
for (mode_name in names(inc_mode_groups)) {
  split_file <- file.path(sim_data_dir,
                           sprintf("sim_results_%s_%s_%s.rds",
                                   estimation_mode, mode_name, timestamp_str))
  saveRDS(inc_mode_groups[[mode_name]], split_file)
  cat(sprintf("  Split results saved: %s (%d rows)\n",
              basename(split_file), nrow(inc_mode_groups[[mode_name]])))
}

cat(sprintf("\nTotal files saved: %d (1 combined + %d per incidence mode)\n",
            1 + length(inc_mode_groups), length(inc_mode_groups)))
