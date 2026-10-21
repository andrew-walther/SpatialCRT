# ============================================================
# Script: check_a3v2_reproduction.R
# Purpose: Verify the 2026-10-21 A3 rerun (profile a3v2) against the 2026-10-20
#   A3 run: identical seeds and code paths must give bit-identical Wald and
#   size-adjusted counts, except saturation subsets whose regions changed under
#   the >= 2-college rule (V6 fix). Also checks the new t + df-corrected counts.
# Author: Andrew Walther
# Created: 2026-10-21
# Dependencies: base R
# ============================================================
# Run from the project root:  Rscript application/tests/check_a3v2_reproduction.R [old] [new]
# Exit status 1 (and "FAIL" lines) on any mismatch, so an unattended queue stops.
#
# Rules
#   - Same task list and the same cells in every task.
#   - Designs other than 3 and 8 (Plain / Incidence-Guided Saturation): every stat of
#     the old run (n_fit ... n_rej_adj) and every allocation diagnostic is identical in
#     every cell, subset and allocation, including subsets whose regions were rebuilt.
#   - Designs 3 and 8: subsets whose regions did not change are identical in
#     n_fit ... n_unident. n_rej_adj is compared only when no subset of the task
#     changed (the size-adjusted critical value pools all subsets).
#   - Every subset of a3v2 has saturation regions of >= 2 colleges.
#   - t + df test: n_rej_tdf <= n_rej and n_cov_tdf >= n_cov everywhere (a wider
#     interval can only reject less and cover more).

args <- commandArgs(TRUE)
root <- "application/results/methods_ext_2026-10/participation"
old_dir <- file.path(root, if (length(args) > 0) args[1] else "a3", "checkpoints")
new_dir <- file.path(root, if (length(args) > 1) args[2] else "a3v2", "checkpoints")
fails <- character(0); fail <- function(...) fails <<- c(fails, sprintf(...))

old_f <- sort(list.files(old_dir, "\\.rds$")); new_f <- sort(list.files(new_dir, "\\.rds$"))
if (!identical(old_f, new_f)) fail("Task lists differ: %d old vs %d new checkpoints", length(old_f), length(new_f))
common <- intersect(old_f, new_f)
base_stats <- c("n_fit", "n_est", "n_ci", "n_rej", "n_cov", "sum_err", "sum_err2",
                "n_alias", "n_bound", "n_error", "n_warn", "n_unident")
n_changed <- 0L; n_sub <- 0L; n_cmp <- 0
for (f in common) {
  o <- readRDS(file.path(old_dir, f)); n <- readRDS(file.path(new_dir, f))
  if (!identical(o$cells, n$cells)) { fail("%s: cells differ", f); next }
  id <- o$task$Design_ID
  chg <- if (is.null(n$subsets$Regions_Changed)) rep(FALSE, nrow(n$subsets)) else n$subsets$Regions_Changed
  n_changed <- n_changed + sum(chg); n_sub <- n_sub + length(chg)
  if (any(n$subsets$Min_Region_Size < 2)) fail("%s: a subset has a region below 2 colleges", f)
  # Regions_Changed is a property of the subset, recorded for every design run on it;
  # only designs 3 and 8 use the regions, so every other design must match on ALL
  # subsets (corrected 2026-10-21 07:35: the first version wrongly failed when a
  # non-saturation task's subset had rebuilt regions).
  keep <- if (id %in% c(3, 8)) which(!chg) else seq_along(chg)
  if (!(id %in% c(3, 8)) && !identical(o$allocations, n$allocations)) fail("%s: allocations differ for a non-saturation design", f)
  so <- o$stats[, keep, , base_stats, drop = FALSE]; sn <- n$stats[, keep, , base_stats, drop = FALSE]
  if (!identical(so, sn)) fail("%s: base stats differ on unchanged subsets (max |diff| %g)", f, max(abs(so - sn), na.rm = TRUE))
  n_cmp <- n_cmp + length(so)
  if ((!any(chg) || !(id %in% c(3, 8))) && !identical(o$stats[, , , "n_rej_adj"], n$stats[, , , "n_rej_adj"]))
    fail("%s: size-adjusted rejections differ with no changed subsets", f)
  if (any(n$stats[, , , "n_rej_tdf"] > n$stats[, , , "n_rej"], na.rm = TRUE)) fail("%s: t + df rejects more than Wald", f)
  if (any(n$stats[, , , "n_cov_tdf"] < n$stats[, , , "n_cov"], na.rm = TRUE)) fail("%s: t + df covers less than Wald", f)
}
cat(sprintf("a3v2 reproduction: %d tasks compared, %s stat cells identical-checked; %d of %d subset-task pairs had saturation regions rebuilt (designs 3/8)\n",
            length(common), format(n_cmp, big.mark = ","), n_changed, n_sub))
if (length(fails)) { cat("FAIL\n", paste(head(fails, 50), collapse = "\n"), "\n"); quit(status = 1) }
cat("PASS\n")
