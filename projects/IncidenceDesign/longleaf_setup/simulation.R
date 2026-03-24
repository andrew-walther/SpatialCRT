# =============================================================================
# simulation.R — Per-scenario HPC script for Longleaf SLURM job arrays
#
# Each SLURM array task runs this script with a task_id (1-2560).
# The task_id indexes into a flat parameter grid of all 2,560 scenarios.
# Each task independently rebuilds the spatial grid, generates incidence/epsilon
# with deterministic seeding, runs one scenario, and saves one .rds file.
#
# Usage:
#   Rscript simulation.R <task_id>          # called by submit_array.sl
#   Rscript simulation.R 1                  # test locally with task 1
#
# Output:
#   results/rds/scenario_<task_id>.rds      # 1-row data.frame per scenario
#
# NOTE ON RESAMPLES (Phase 2 — increase after validation):
#   After verifying the pipeline works with the current values (25 x 10),
#   increase to e.g. 50 x 50 for production-quality precision, then also
#   increase --time in submit_array.sl from 00:30:00 to 01:00:00.
# =============================================================================

# --- Libraries ---
library(spdep)
library(spatialreg)
library(dplyr)
library(digest)

# =============================================================================
# 1. RECEIVE TASK ID
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

task_id <- if (length(args) > 0) {
  as.integer(args[1])
} else {
  message("No task ID supplied -- running locally with task_id = 1")
  1L
}

cat("=== Task ID:", task_id, "===\n")

# =============================================================================
# 2. SOURCE SIMULATION MODULES (01-04)
# =============================================================================

# Resolve the path to code/ relative to this script's location
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[grep("--file=", script_args)])
if (length(script_path) > 0) {
  script_dir <- dirname(normalizePath(script_path))
} else {
  # Fallback for interactive/source() usage
  script_dir <- normalizePath(getwd())
}

code_dir <- normalizePath(file.path(script_dir, "..", "code"))

# Verify the code directory exists
if (!dir.exists(code_dir)) {
  stop("Cannot find code directory at: ", code_dir,
       "\n  Expected ../code/ relative to this script's location.",
       "\n  Script dir: ", script_dir)
}

cat("Code directory:", code_dir, "\n")

source(file.path(code_dir, "01_spatial_setup.R"))
source(file.path(code_dir, "02_incidence_generation.R"))
source(file.path(code_dir, "03_designs.R"))
source(file.path(code_dir, "04_estimation.R"))

cat("Modules loaded.\n")

# =============================================================================
# 3. CONFIGURABLE PARAMETERS
# =============================================================================

# Estimation mode
estimation_mode <- "MLE"

# Grid
grid_dim <- 10  # 10x10 = 100 clusters

# Resample counts -- EDIT THESE FOR PRODUCTION RUNS
n_design_resamples  <- 25    # Phase 2: increase to 50
n_outcome_resamples <- 10    # Phase 2: increase to 50

# SDM true parameters
true_tau <- 1.0
beta     <- 1.0
sigma    <- 1.0

# Poisson-specific
base_rate       <- 35 / 100000
pop_per_cluster <- 1000
pop_mode        <- "equal"

# Oracle spillover
include_spill_covariate <- TRUE

# =============================================================================
# 4. BUILD PARAMETER GRID (2,560 scenarios)
# =============================================================================

# The 5 incidence configurations
inc_configs <- data.frame(
  inc_config_id = 1:5,
  inc_mode = c("iid", "spatial", "spatial", "poisson", "poisson"),
  rho_x    = c(0.00, 0.20, 0.50, 0.20, 0.50),
  stringsAsFactors = FALSE
)

# Full factorial parameter grid
param_grid <- expand.grid(
  inc_config_id = 1:5,
  nb_type       = c("rook", "queen"),
  rho           = c(0.00, 0.01, 0.20, 0.50),
  gamma_val     = c(0.5, 0.6, 0.7, 0.8),
  spill_type    = c("control_only", "both"),
  design_id     = 1:8,
  stringsAsFactors = FALSE
)

total_scenarios <- nrow(param_grid)
cat("Total scenarios in grid:", total_scenarios, "\n")

# Validate task_id
if (task_id < 1 || task_id > total_scenarios) {
  stop("task_id ", task_id, " is out of range [1, ", total_scenarios, "]")
}

# Extract this task's parameters
params <- param_grid[task_id, ]
ic     <- inc_configs[params$inc_config_id, ]

inc_mode   <- ic$inc_mode
rho_x      <- ic$rho_x
nb_type    <- params$nb_type
rho        <- params$rho
gamma_val  <- params$gamma_val
spill_type <- params$spill_type
d_id       <- params$design_id

cat(sprintf("Scenario: %s(rho_x=%.2f) | nb=%s | rho=%.2f | gamma=%.2f | spill=%s | D%d\n",
            inc_mode, rho_x, nb_type, rho, gamma_val, spill_type, d_id))

# =============================================================================
# 5. BUILD SPATIAL GRID
# =============================================================================

cat("Building spatial grid...\n")
grid_obj <- build_spatial_grid(grid_dim)
N <- grid_obj$N_clusters
coords <- grid_obj$coords
I_mat <- diag(N)

# Select active spatial structure (rook or queen)
spatial <- get_active_spatial(grid_obj, nb_type)
W <- spatial$W
active_listw <- spatial$listw
active_nb <- spatial$nb

# =============================================================================
# 6. GENERATE INCIDENCE (deterministic seed per config)
#
# All tasks sharing the same incidence config produce identical X_matrix.
# =============================================================================

inc_seed <- digest2int(paste("incidence", inc_mode, rho_x, sep = "|"))
set.seed(inc_seed)

X_matrix <- generate_incidence(
  mode          = inc_mode,
  N             = N,
  n_resamples   = n_outcome_resamples,
  W             = grid_obj$W_queen,    # always queen for incidence generation
  rho_incidence = rho_x,
  base_rate     = base_rate,
  pop_per_cluster = pop_per_cluster,
  pop_mode      = pop_mode
)
base_incidence <- X_matrix[, 1]

cat(sprintf("Incidence: range [%.4f, %.4f], mean %.4f\n",
            min(X_matrix), max(X_matrix), mean(X_matrix)))

# =============================================================================
# 7. GENERATE EPSILON MATRIX (deterministic seed per config + nb_type + rho)
#
# All tasks sharing the same (config, nb_type, rho) produce identical epsilon.
# =============================================================================

eps_seed <- digest2int(paste("epsilon", inc_mode, rho_x, nb_type, rho, sep = "|"))
set.seed(eps_seed)

Epsilon_matrix <- matrix(rnorm(N * n_outcome_resamples, 0, sigma),
                         nrow = N, ncol = n_outcome_resamples)

# =============================================================================
# 8. COMPUTE BASELINE OUTCOME COMPONENTS
# =============================================================================

inv_mat <- solve(I_mat - rho * W)
baseline_part <- inv_mat %*% (X_matrix * beta + Epsilon_matrix)

# =============================================================================
# 9. GENERATE TREATMENT ASSIGNMENTS (per-scenario seed)
# =============================================================================

scenario_seed <- digest2int(paste(inc_mode, rho_x, nb_type, rho,
                                  gamma_val, spill_type, d_id, sep = "|"))
set.seed(scenario_seed)

if (is_design_deterministic(d_id)) {
  Z_single <- get_designs(d_id, 1, N, base_incidence, active_nb, coords)[, 1]
  Z_matrix <- matrix(rep(Z_single, n_design_resamples),
                     nrow = N, ncol = n_design_resamples)
} else {
  Z_matrix <- get_designs(d_id, n_design_resamples, N,
                          base_incidence, active_nb, coords)
}

# =============================================================================
# 10. RUN SIMULATION + ESTIMATION
# =============================================================================

cat("Running estimation...\n")
start_time <- proc.time()

all_estimates <- numeric(0)
all_ci_lower  <- numeric(0)
all_ci_upper  <- numeric(0)

for (d_resample in seq_len(n_design_resamples)) {
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

elapsed <- proc.time() - start_time
cat(sprintf("Estimation done in %.1f seconds\n", elapsed["elapsed"]))

# =============================================================================
# 11. AGGREGATE INTO 1-ROW DATA.FRAME
# =============================================================================

valid_idx <- !is.na(all_estimates)
valid_est <- all_estimates[valid_idx]
fail_rate <- 1 - sum(valid_idx) / length(all_estimates)

if (length(valid_est) > 0) {
  valid_ci <- valid_idx & !is.na(all_ci_lower) & !is.na(all_ci_upper)
  coverage_val <- if (sum(valid_ci) > 0) {
    mean((all_ci_lower[valid_ci] <= true_tau) &
         (all_ci_upper[valid_ci] >= true_tau))
  } else NA

  result <- data.frame(
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
  result <- data.frame(
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

# =============================================================================
# 12. SAVE OUTPUT
# =============================================================================

dir.create("results/rds", showWarnings = FALSE, recursive = TRUE)

rds_path <- file.path("results/rds", paste0("scenario_", task_id, ".rds"))
saveRDS(result, rds_path)
cat("Saved:", rds_path, "\n")
cat("Task", task_id, "complete.\n")
