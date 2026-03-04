# ==============================================================================
# 05_run_simulation.R
# Main orchestrator: sources modules 01-04, defines parameters, runs the
# nested simulation loop, saves results.
#
# Usage:
#   source("05_run_simulation.R")
#   -- or --
#   Rscript 05_run_simulation.R
# ==============================================================================

# --- Load Libraries ---
library(spdep)
library(spatialreg)
library(dplyr)
library(digest)

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
estimation_mode <- "DIM"      # "MLE" or "DIM"

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
  n_designs_resamples <- 100
  n_outcome_resamples <- 1000
}

# Designs to evaluate
design_ids <- c(1, 2, 3, 4, 5, 7)

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
# COMPUTE TOTAL SCENARIO COUNT & ETA TRACKER
# ==============================================================================

# Build incidence config list: (mode, rho_X) pairs
inc_configs <- list()
for (inc_mode in incidence_modes) {
  rho_sweep <- if (inc_mode == "iid") 0 else rho_incidence_vals
  for (rho_x in rho_sweep) {
    inc_configs[[length(inc_configs) + 1]] <- list(mode = inc_mode, rho_x = rho_x)
  }
}

total_scenarios <- length(inc_configs) * length(neighbor_types) * length(rho_vals) *
                   length(gamma_vals) * length(spill_types) * length(design_ids)

cat(sprintf("\n=== Simulation Configuration ===\n"))
cat(sprintf("Estimation mode:    %s\n", estimation_mode))
cat(sprintf("Incidence configs:  %d\n", length(inc_configs)))
cat(sprintf("Total scenarios:    %d\n", total_scenarios))
cat(sprintf("Iterations/scenario: %d (%d design x %d outcome)\n",
            n_designs_resamples * n_outcome_resamples,
            n_designs_resamples, n_outcome_resamples))
cat(sprintf("================================\n\n"))

# ETA tracking
scenario_times <- numeric(0)
scenario_counter <- 0
global_start <- Sys.time()

# ==============================================================================
# RESULTS STORAGE
# ==============================================================================

results_list <- vector("list", total_scenarios)
list_index <- 1

# ==============================================================================
# MAIN SIMULATION LOOP
# ==============================================================================

for (ic in seq_along(inc_configs)) {
  inc_mode <- inc_configs[[ic]]$mode
  rho_x    <- inc_configs[[ic]]$rho_x

  # --- Generate incidence ONCE per config (using queen W as reference) ---
  inc_label <- if (inc_mode == "iid") "iid" else sprintf("%s(rho_X=%.2f)", inc_mode, rho_x)
  cat(sprintf("\n--- Incidence config %d/%d: %s ---\n", ic, length(inc_configs), inc_label))

  X_matrix <- generate_incidence(
    mode         = inc_mode,
    N            = N,
    n_resamples  = n_outcome_resamples,
    W            = grid_obj$W_queen,
    rho_incidence = rho_x,
    base_rate    = base_rate,
    pop_per_cluster = pop_per_cluster,
    pop_mode     = pop_mode
  )

  # "Historical" incidence for design decisions (first column)
  base_incidence <- X_matrix[, 1]

  for (nb_type in neighbor_types) {
    spatial <- get_active_spatial(grid_obj, nb_type)
    W <- spatial$W
    active_listw <- spatial$listw
    active_nb <- spatial$nb

    for (rho in rho_vals) {
      # SDM spatial multiplier: (I - rho*W)^{-1}
      inv_mat <- solve(I_mat - rho * W)

      # Baseline part: inv_mat %*% (X * beta + epsilon)
      Epsilon_matrix <- matrix(rnorm(N * n_outcome_resamples, 0, sigma),
                               nrow = N, ncol = n_outcome_resamples)
      baseline_part <- inv_mat %*% (X_matrix * beta + Epsilon_matrix)

      for (gamma_val in gamma_vals) {
        for (spill_type in spill_types) {
          for (d_id in design_ids) {

            scenario_start <- Sys.time()
            scenario_counter <- scenario_counter + 1

            # Per-scenario deterministic seed
            seed_string <- paste(inc_mode, rho_x, nb_type, rho,
                                 gamma_val, spill_type, d_id, sep = "|")
            scenario_seed <- digest2int(seed_string)
            set.seed(scenario_seed)

            # Progress logging with ETA
            avg_time <- if (length(scenario_times) > 0) mean(scenario_times) else NA
            eta_str <- if (!is.na(avg_time)) {
              remaining <- avg_time * (total_scenarios - scenario_counter + 1)
              if (remaining > 3600) {
                sprintf("ETA: %.1f hrs", remaining / 3600)
              } else if (remaining > 60) {
                sprintf("ETA: %.1f min", remaining / 60)
              } else {
                sprintf("ETA: %.0f sec", remaining)
              }
            } else {
              "ETA: calculating..."
            }

            cat(sprintf("[%4d/%d] %s | NB:%-5s | rho:%.2f | gam:%.2f | spill:%-12s | D%d | %s\n",
                        scenario_counter, total_scenarios, inc_label,
                        nb_type, rho, gamma_val, spill_type, d_id, eta_str))
            flush.console()

            # Generate treatment assignments
            if (is_design_deterministic(d_id)) {
              Z_single <- get_designs(d_id, 1, N, base_incidence, active_nb, coords)[, 1]
              Z_matrix <- matrix(rep(Z_single, n_designs_resamples),
                                 nrow = N, ncol = n_designs_resamples)
            } else {
              Z_matrix <- get_designs(d_id, n_designs_resamples, N,
                                      base_incidence, active_nb, coords)
            }

            # Collect estimates across all (design resample x outcome resample)
            all_estimates <- numeric(0)
            all_ci_lower  <- numeric(0)
            all_ci_upper  <- numeric(0)

            for (d_resample in seq_len(n_designs_resamples)) {
              Z <- Z_matrix[, d_resample]

              # Spillover term: WZ with type mask
              WZ <- as.vector(W %*% Z)
              spill_term <- if (spill_type == "control_only") {
                gamma_val * WZ * (1 - Z)
              } else {
                gamma_val * WZ
              }

              # Outcome: Y = inv_mat %*% (tau*Z + spillover) + baseline
              trt_effect <- as.vector(inv_mat %*% (true_tau * Z + spill_term))
              Y_sim <- sweep(baseline_part, 1, trt_effect, "+")

              # Estimate tau
              est_result <- estimate_tau(
                estimation_mode    = estimation_mode,
                Y_sim              = Y_sim,
                Z                  = Z,
                spill_term         = spill_term,
                X_matrix           = X_matrix,
                active_listw       = active_listw,
                n_outcome_resamples = n_outcome_resamples,
                include_spill_covariate = include_spill_covariate
              )

              all_estimates <- c(all_estimates, est_result$estimates)
              all_ci_lower  <- c(all_ci_lower,  est_result$ci_lower)
              all_ci_upper  <- c(all_ci_upper,  est_result$ci_upper)
            }

            # Aggregate results
            valid_idx <- !is.na(all_estimates)
            valid_est <- all_estimates[valid_idx]
            fail_rate <- 1 - sum(valid_idx) / length(all_estimates)

            if (length(valid_est) > 0) {
              bias_val <- mean(valid_est) - true_tau
              sd_val   <- sd(valid_est)
              mse_val  <- mean((valid_est - true_tau)^2)

              # Coverage: proportion of CIs that contain true_tau
              valid_ci <- valid_idx & !is.na(all_ci_lower) & !is.na(all_ci_upper)
              if (sum(valid_ci) > 0) {
                covers <- (all_ci_lower[valid_ci] <= true_tau) &
                          (all_ci_upper[valid_ci] >= true_tau)
                coverage_val <- mean(covers)
              } else {
                coverage_val <- NA
              }

              results_list[[list_index]] <- data.frame(
                Incidence_Mode  = inc_mode,
                Rho_Incidence   = rho_x,
                Neighbor_Type   = nb_type,
                Design          = paste("Design", d_id),
                Rho             = rho,
                Gamma           = gamma_val,
                Spillover_Type  = spill_type,
                Mean_Estimate   = mean(valid_est),
                Bias            = bias_val,
                SD              = sd_val,
                MSE             = mse_val,
                Coverage        = coverage_val,
                Fail_Rate       = fail_rate,
                stringsAsFactors = FALSE
              )
            } else {
              results_list[[list_index]] <- data.frame(
                Incidence_Mode  = inc_mode,
                Rho_Incidence   = rho_x,
                Neighbor_Type   = nb_type,
                Design          = paste("Design", d_id),
                Rho             = rho,
                Gamma           = gamma_val,
                Spillover_Type  = spill_type,
                Mean_Estimate   = NA,
                Bias            = NA,
                SD              = NA,
                MSE             = NA,
                Coverage        = NA,
                Fail_Rate       = fail_rate,
                stringsAsFactors = FALSE
              )
            }

            # Timing
            scenario_elapsed <- as.numeric(difftime(Sys.time(), scenario_start, units = "secs"))
            scenario_times <- c(scenario_times, scenario_elapsed)
            mins <- floor(scenario_elapsed / 60)
            secs <- round(scenario_elapsed %% 60, 1)
            cat(sprintf("   -> %d min %.1f sec\n", mins, secs))
            flush.console()

            list_index <- list_index + 1
          }
        }
      }
    }
  }
}

# ==============================================================================
# COMBINE & SAVE RESULTS
# ==============================================================================

granular_results <- bind_rows(results_list)

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\n=== Simulation Complete ===\n"))
cat(sprintf("Total time: %.1f minutes\n", total_elapsed))
cat(sprintf("Scenarios completed: %d\n", nrow(granular_results)))
cat(sprintf("Results columns: %s\n", paste(names(granular_results), collapse = ", ")))

# Save results
results_dir <- file.path(script_dir, "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- file.path(results_dir,
                          sprintf("sim_results_%s_%s.rds", estimation_mode, timestamp_str))
saveRDS(granular_results, output_file)
cat(sprintf("Results saved to: %s\n", output_file))
