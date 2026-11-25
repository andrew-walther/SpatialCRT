# ============================================================
# Script: resign_production_manifest.R
# Purpose: Re-record the case-study production manifest after a change to a hashed file
#          that the case study does not use, failing unless that is the only difference.
# Author: Andrew Walther
# Created: 2026-11-25
# Dependencies: the real_sud modules (via real_sud_summary.R), digest
# ============================================================
#
# Why: the production run (2026-10-02) hashed code/16_allocation_risk.R. On 2026-10-27
# (commit c5d070f) code/16 gained the grid "surfaces10" profile: comments and
# allocation_parameters() only, which the case study never calls. The hash changed, so
# rs_manifest() refused every later rs_exhibits() run. This script re-records the manifest
# only if the source hash of code/16 is the sole difference (same config, inputs, regions,
# R, packages, BLAS, RNG). The original manifest stays in git history.
#
# Run from the project root:  Rscript application/code/resign_production_manifest.R
# Then confirm: Rscript application/code/real_sud_summary.R reproduces every exhibit CSV
# byte for byte (git diff --stat application/results/real_sud_rev_20261002/exhibits/).

real_sud_summary_define_only <- TRUE
source("application/code/real_sud_summary.R")

root <- file.path(.rs_app, "results", "real_sud_rev_20261002")
production <- file.path(root, "production")
b <- readRDS(file.path(production, "results.rds"))
setup <- readRDS(file.path(root, "setup.rds"))
stored <- readRDS(file.path(production, "manifest.rds"))

# Current manifest, written to a scratch folder so nothing is overwritten yet
tmp <- tempfile("manifest_"); dir.create(tmp)
current <- rs_manifest(tmp, setup, b$config)

# Every field except the source hashes must be identical
other <- setdiff(union(names(stored), names(current)), "source")
same <- vapply(other, function(k) identical(stored[[k]], current[[k]]), logical(1))
if (!all(same)) stop("Manifest fields differ beyond source hashes: ", paste(other[!same], collapse = ", "))

# Among the source hashes, only code/16 may differ
if (!identical(names(stored$source), names(current$source))) stop("Hashed file list changed")
changed <- basename(names(stored$source))[stored$source != current$source]
if (!identical(changed, "16_allocation_risk.R"))
  stop("Unexpected changed source files: ", paste(changed, collapse = ", "))

saveRDS(current, file.path(production, "manifest.rds"))
cat("Re-recorded", file.path(production, "manifest.rds"),
    "\nOnly difference: source hash of 16_allocation_risk.R (2026-10-27 surfaces10 profile).\n")
