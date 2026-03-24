# =============================================================================
# aggregate_results.R — Combine outputs from all job array tasks
#
# Run this LOCALLY on your Mac (or interactively on Longleaf) after all
# array tasks have completed.
#
# Usage (local, after scp-ing results/ back):
#   Rscript R/aggregate_results.R
#
# Usage (on Longleaf interactively):
#   module add r/4.4.0
#   Rscript R/aggregate_results.R
# =============================================================================

results_dir <- "results"

# -----------------------------------------------------------------------------
# Option A: Aggregate from .rds files (preserves R types exactly)
# -----------------------------------------------------------------------------

rds_files <- list.files(
  file.path(results_dir, "rds"),
  pattern  = "\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No .rds files found in results/rds/ — check your output path.")
}

cat("Found", length(rds_files), "RDS files\n")

all_results <- do.call(rbind, lapply(rds_files, readRDS))

# Sort by task_id for clean output
all_results <- all_results[order(all_results$task_id), ]

cat("Total rows:", nrow(all_results), "\n")
print(head(all_results))

# Save combined output
saveRDS(all_results, file.path(results_dir, "all_results.rds"))
write.csv(all_results, file.path(results_dir, "all_results.csv"), row.names = FALSE)
cat("Saved combined results to results/all_results.rds and .csv\n")

# -----------------------------------------------------------------------------
# Option B: Aggregate from .csv files (simpler, no R needed to inspect)
# -----------------------------------------------------------------------------

# csv_files <- list.files(file.path(results_dir, "csv"), pattern = "\\.csv$", full.names = TRUE)
# all_results_csv <- do.call(rbind, lapply(csv_files, read.csv))

# -----------------------------------------------------------------------------
# Check for missing tasks (did any array jobs fail?)
# -----------------------------------------------------------------------------

expected_ids <- 1:100   # change to match your --array= range

completed_ids <- sort(unique(all_results$task_id))
missing_ids   <- setdiff(expected_ids, completed_ids)

if (length(missing_ids) == 0) {
  cat("All", length(expected_ids), "tasks completed successfully.\n")
} else {
  cat("WARNING:", length(missing_ids), "tasks missing:", paste(missing_ids, collapse = ", "), "\n")
  cat("Check logs/slurm-*_<task_id>.err for those task IDs.\n")
  cat("Resubmit failed tasks with: sbatch --array=",
      paste(missing_ids, collapse = ","), "submit_array.sl\n")
}

# -----------------------------------------------------------------------------
# Quick summary statistics
# -----------------------------------------------------------------------------

cat("\nSummary of estimates:\n")
print(summary(all_results$estimate))
