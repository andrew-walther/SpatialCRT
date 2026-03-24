# =============================================================================
# simulation.R — Template R simulation script for Longleaf job arrays
#
# This script is called by submit_array.sl with the SLURM array task ID
# as a command-line argument. Each array task runs independently with its
# own seed / parameter set.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Receive task ID from SLURM
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

task_id <- if (length(args) > 0) {
  as.integer(args[1])
} else {
  # Fallback: run locally with task_id = 1 for testing
  message("No task ID supplied — running locally with task_id = 1")
  1L
}

cat("Task ID:", task_id, "\n")

# -----------------------------------------------------------------------------
# 2. Set reproducible seed from task ID
# -----------------------------------------------------------------------------

set.seed(task_id)
cat("Seed set to:", task_id, "\n")

# -----------------------------------------------------------------------------
# 3. Define your parameter grid (optional)
#
# If your simulations vary by MORE than just seed (e.g. different sample sizes,
# effect sizes, or model settings), define a parameter grid here and use
# task_id to index into it.
#
# Example: 100 tasks covering 10 sample sizes × 10 effect sizes
# -----------------------------------------------------------------------------

# Uncomment and adapt this block if you have a parameter grid:
# param_grid <- expand.grid(
#   n       = c(50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000),
#   effect  = c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0)
# )
# params <- param_grid[task_id, ]
# n      <- params$n
# effect <- params$effect

# For seed-only variation, just use task_id as seed (already done above).

# -----------------------------------------------------------------------------
# 4. Your simulation code goes here
# -----------------------------------------------------------------------------

cat("Starting simulation...\n")
start_time <- proc.time()

# --- REPLACE THIS SECTION WITH YOUR ACTUAL SIMULATION ---

# Example: simple Monte Carlo
n_obs    <- 1000
effect   <- 0.5
x        <- rnorm(n_obs)
y        <- effect * x + rnorm(n_obs)
fit      <- lm(y ~ x)
results  <- data.frame(
  task_id   = task_id,
  seed      = task_id,
  n         = n_obs,
  estimate  = coef(fit)["x"],
  se        = sqrt(vcov(fit)["x", "x"]),
  p_value   = summary(fit)$coefficients["x", "Pr(>|t|)"]
)

# --------------------------------------------------------

elapsed <- proc.time() - start_time
cat("Simulation done in", round(elapsed["elapsed"], 1), "seconds\n")

# -----------------------------------------------------------------------------
# 5. Save outputs — use task_id in filenames to avoid collisions
# -----------------------------------------------------------------------------

# Create output directories if they don't exist (safe for parallel runs)
dir.create("results/rds",     showWarnings = FALSE, recursive = TRUE)
dir.create("results/csv",     showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# .rds — full R object, exact precision
rds_path <- file.path("results/rds", paste0("sim_", task_id, ".rds"))
saveRDS(results, rds_path)
cat("Saved RDS:", rds_path, "\n")

# CSV — human-readable, easy to load later
csv_path <- file.path("results/csv", paste0("sim_", task_id, ".csv"))
write.csv(results, csv_path, row.names = FALSE)
cat("Saved CSV:", csv_path, "\n")

# Figure — save as PNG (no display available on HPC)
fig_path <- file.path("results/figures", paste0("sim_", task_id, ".png"))
png(fig_path, width = 800, height = 600, res = 120)
plot(x, y, main = paste("Task", task_id), col = "steelblue", pch = 19, cex = 0.5)
abline(fit, col = "red", lwd = 2)
dev.off()
cat("Saved figure:", fig_path, "\n")

cat("All outputs saved. Task", task_id, "complete.\n")
