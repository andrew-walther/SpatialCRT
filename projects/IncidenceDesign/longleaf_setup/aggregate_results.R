# =============================================================================
# aggregate_results.R — Combine per-scenario outputs from SLURM job array
#
# Run this AFTER all 12,800 array tasks have completed.
# Produces the same output format that 06_visualizations.R expects.
# Output files use the "MLE_tau_sweep" mode tag to avoid overwriting the
# existing MLE_combined baseline (tau = 1.0 only).
#
# Usage (on Longleaf):
#   module add r/4.4.0
#   Rscript aggregate_results.R
#
# Usage (locally, after scp-ing results/ back):
#   Rscript aggregate_results.R
# =============================================================================

library(dplyr)

# =============================================================================
# 1. READ ALL PER-SCENARIO .rds FILES
# =============================================================================

rds_dir <- "results/rds"
rds_files <- list.files(rds_dir, pattern = "^scenario_.*\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No scenario_*.rds files found in ", rds_dir,
       "\n  Make sure all SLURM tasks have completed.")
}

cat("Found", length(rds_files), "scenario files\n")

all_results <- bind_rows(lapply(rds_files, readRDS))
all_results <- all_results[order(as.integer(
  sub("scenario_(\\d+)\\.rds", "\\1", basename(rds_files))
)), ]

cat("Combined data frame:", nrow(all_results), "rows\n")

# =============================================================================
# 2. CHECK FOR MISSING TASKS
# =============================================================================

expected_total <- 12800  # 5 tau x 5 configs x 2 nb x 4 rho x 4 gamma x 2 spill x 8 designs

# Extract task IDs from filenames
completed_ids <- sort(as.integer(
  sub("scenario_(\\d+)\\.rds", "\\1", basename(rds_files))
))
missing_ids <- setdiff(1:expected_total, completed_ids)

if (length(missing_ids) == 0) {
  cat("All", expected_total, "tasks completed successfully.\n")
} else {
  cat("WARNING:", length(missing_ids), "tasks missing!\n")
  if (length(missing_ids) <= 20) {
    cat("  Missing IDs:", paste(missing_ids, collapse = ", "), "\n")
  } else {
    cat("  First 20 missing:", paste(head(missing_ids, 20), collapse = ", "), "...\n")
  }
  cat("  Check logs/slurm-*_<task_id>.err for those task IDs.\n")
  cat("  Resubmit with: sbatch --array=",
      paste(missing_ids, collapse = ","), " submit_array.sl\n", sep = "")
}

# =============================================================================
# 3. SAVE COMBINED RESULTS (same format as 05_run_simulation.R output)
# =============================================================================

# Save to results/longleaf/ to keep separate from local Mac results.
# Load with: load_latest_results(results_dir = "results/longleaf/sim_data")
longleaf_dir <- file.path("..", "results", "longleaf")
sim_data_dir <- file.path(longleaf_dir, "sim_data")
reports_dir  <- file.path(longleaf_dir, "reports")
dir.create(sim_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(reports_dir,  showWarnings = FALSE, recursive = TRUE)

timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Combined file — use MLE_tau_sweep tag to preserve existing MLE_combined baseline
combined_file <- file.path(sim_data_dir,
                           sprintf("sim_results_MLE_tau_sweep_combined_%s.rds", timestamp_str))
saveRDS(all_results, combined_file)
cat("Saved combined:", combined_file, "\n")

# Per-incidence-mode splits
mode_groups <- split(all_results, all_results$Incidence_Mode)
for (mode_name in names(mode_groups)) {
  split_file <- file.path(sim_data_dir,
                          sprintf("sim_results_MLE_tau_sweep_%s_%s.rds", mode_name, timestamp_str))
  saveRDS(mode_groups[[mode_name]], split_file)
  cat(sprintf("  Split saved: %s (%d rows)\n", basename(split_file), nrow(mode_groups[[mode_name]])))
}

cat(sprintf("\nTotal files saved: %d (1 combined + %d per incidence mode)\n",
            1 + length(mode_groups), length(mode_groups)))

# =============================================================================
# 4. SUMMARY STATISTICS
# =============================================================================

cat("\n=== Summary ===\n")
cat("Scenarios:", nrow(all_results), "\n")
cat("Designs:", paste(sort(unique(all_results$Design)), collapse = ", "), "\n")
cat("Incidence modes:", paste(unique(all_results$Incidence_Mode), collapse = ", "), "\n")

cat("\nMSE by Design (averaged across all scenarios):\n")
mse_summary <- all_results %>%
  group_by(Design) %>%
  summarise(Mean_MSE = mean(MSE, na.rm = TRUE),
            Mean_Coverage = mean(Coverage, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(Mean_MSE)
print(as.data.frame(mse_summary))

cat("\nFail rate:", mean(all_results$Fail_Rate, na.rm = TRUE), "\n")

cat("\n=== Next Steps: Generate Reports ===\n")
cat("Results saved to:", normalizePath(sim_data_dir, mustWork = FALSE), "\n")
cat("Reports directory:", normalizePath(reports_dir, mustWork = FALSE), "\n")
cat("\nTo generate reports from these results, run in R:\n")
cat('  source("code/06_visualizations.R")\n')
cat('  results <- load_latest_results(results_dir = "results/longleaf/sim_data",\n')
cat('                                 estimation_mode = "MLE_tau_sweep_combined")\n')
cat('  run_all_visualizations(results=results, results_dir="results/longleaf",\n')
cat('                         estimation_mode="MLE_tau_sweep")\n')
cat('\n  source("code/08_design_recommendations.R")\n')
cat('  run_recommendation_report(results=results, estimation_mode="MLE_tau_sweep",\n')
cat('                            results_dir="results/longleaf")\n')
