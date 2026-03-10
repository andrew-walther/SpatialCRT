# 04_run_simulation.R
# Parallel orchestrator for the SpillSpatialDepSim unified simulation.
#
# Runs all 48 scenarios (3 grids x 4 psi x 2 rho x 2 spillover types) using
# mclapply() on Mac/Linux. Each worker is an independent forked process — no
# cluster setup overhead. Workers share the pre-built grid/combo objects via
# OS copy-on-write semantics.
#
# Output: ../results/sim_results_unified.rds
#
# Usage (from terminal, recommended for long runs):
#   cd projects/SpillSpatialDepSim/code
#   caffeinate -i Rscript 04_run_simulation.R
#
# Usage (from RStudio):
#   setwd("projects/SpillSpatialDepSim/code")
#   source("04_run_simulation.R")

# --- Detect script directory (works from Rscript or RStudio source()) ---
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
results_dir <- file.path(dirname(script_dir), "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --- Source dependencies ---
source(file.path(script_dir, "01_grid_setup.R"))
source(file.path(script_dir, "02_simulation_core.R"))

library(parallel)

cat("=== SpillSpatialDepSim Unified Simulation ===\n")
cat(sprintf("Results directory: %s\n", results_dir))

# ============================================================
# CONFIGURABLE PARAMETERS
# ============================================================
psi_vals  <- c(0.5, 0.6, 0.7, 0.8)   # spillover effect values
rho_vals  <- c(0.00, 0.01)            # spatial autocorrelation values
alpha_val <- 0.2                       # true intercept
beta_val  <- 1.0                       # true intervention effect
sd_val    <- 0.1                       # residual SD

# Grid configurations: name -> (n_rows, n_cols, N_iterations)
# N matches original scripts: 50 for 2x4/3x3, 10 for 3x4 (924 combos)
grid_configs <- list(
  "2x4" = list(n_rows = 2, n_cols = 4, N = 50),
  "3x3" = list(n_rows = 3, n_cols = 3, N = 50),
  "3x4" = list(n_rows = 3, n_cols = 4, N = 10)
)

n_cores <- max(1L, detectCores() - 1L)  # leave one core free for system
cat(sprintf("Using %d cores\n\n", n_cores))

# ============================================================
# PRE-BUILD GRID OBJECTS (shared read-only across workers)
# ============================================================
cat("Building grid objects and combination sets...\n")
grids <- lapply(names(grid_configs), function(nm) {
  cfg     <- grid_configs[[nm]]
  grid    <- setup_grid(n_rows = cfg$n_rows, n_cols = cfg$n_cols)
  combos  <- combn(seq_len(grid$n_districts), grid$n_trt, simplify = FALSE)
  blk_idx <- build_valid_block_indices(grid, combos)
  cat(sprintf("  %s: %d districts, %d combos, %d valid block combos, N=%d\n",
              nm, grid$n_districts, length(combos), length(blk_idx), cfg$N))
  list(grid = grid, combos = combos, blk_idx = blk_idx, N = cfg$N)
})
names(grids) <- names(grid_configs)

# ============================================================
# BUILD FLAT SCENARIO LIST (48 total)
# ============================================================
scenarios <- expand.grid(
  grid_name = names(grid_configs),
  psi       = psi_vals,
  rho       = rho_vals,
  trt_spill = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

cat(sprintf("\nRunning %d scenarios...\n\n", nrow(scenarios)))
overall_start <- proc.time()

# ============================================================
# WORKER FUNCTION — one call = one full scenario
# ============================================================
run_one <- function(i) {
  s   <- scenarios[i, ]
  env <- grids[[s$grid_name]]

  # Deterministic per-scenario seed: reproducible regardless of run order
  set.seed(2024L + i)

  result <- tryCatch(
    run_scenario(
      grid            = env$grid,
      combinations    = env$combos,
      valid_block_idx = env$blk_idx,
      psi             = s$psi,
      rho             = s$rho,
      trt_spill       = s$trt_spill,
      N               = env$N,
      alpha           = alpha_val,
      beta            = beta_val,
      sd              = sd_val
    ),
    error = function(e) {
      cat(sprintf("[%d/%d] ERROR: %s | psi=%.1f | rho=%.2f | TrtSpill=%s — %s\n",
                  i, nrow(scenarios), s$grid_name, s$psi, s$rho, s$trt_spill, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(result)) {
    result$grid      <- s$grid_name
    result$psi       <- s$psi
    result$rho       <- s$rho
    result$trt_spill <- s$trt_spill
    cat(sprintf("[%d/%d] Done: %s | psi=%.1f | rho=%.2f | TrtSpill=%s\n",
                i, nrow(scenarios), s$grid_name, s$psi, s$rho, s$trt_spill))
  }
  result
}

# ============================================================
# EXECUTE PARALLEL
# ============================================================
results_list <- mclapply(
  seq_len(nrow(scenarios)),
  run_one,
  mc.cores = n_cores
)

# ============================================================
# COMBINE & SAVE
# ============================================================
n_failed <- sum(vapply(results_list, is.null, logical(1)))
if (n_failed > 0) {
  warning(sprintf("%d scenario(s) failed — they will be excluded from the saved results.", n_failed))
}

final <- do.call(rbind, Filter(Negate(is.null), results_list))

# Column order
col_order <- c("grid", "psi", "rho", "trt_spill", "combo_id", "is_block",
               "param", "bias", "bias_abs", "variance", "MSE")
final <- final[, col_order]

out_path <- file.path(results_dir, "sim_results_unified.rds")
saveRDS(final, out_path)

elapsed <- proc.time() - overall_start
cat(sprintf("\n=== Complete ===\n"))
cat(sprintf("Rows: %d | Scenarios: %d/%d succeeded\n",
            nrow(final), nrow(scenarios) - n_failed, nrow(scenarios)))
cat(sprintf("Saved: %s\n", out_path))
cat(sprintf("Total elapsed: %.1f minutes\n", elapsed["elapsed"] / 60))

# Quick sanity check
cat("\nRow counts by grid and param:\n")
print(table(final$grid, final$param))
