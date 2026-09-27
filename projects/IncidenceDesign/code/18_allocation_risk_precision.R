# ============================================================
# Script: 18_allocation_risk_precision.R
# Purpose: Extend 50 matched pilot blocks from R = 100 to R = 400 on the same assignments.
# Author: Andrew Walther
# Created: 2026-09-27
# Dependencies: unchanged IncidenceDesign16 allocation-risk runner
# ============================================================
# Keeps pilot source/checkpoints/results unchanged.
# Run from code/: VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# OMP_NUM_THREADS=1 Rscript 18_allocation_risk_precision.R
# Eight local fork workers each use single-thread BLAS. Memory per worker is one
# 100-cluster block's 100-allocation × 400-outcome records; not a SLURM runner.
allocation_define_only <- TRUE
allocation_script_dir <- local({
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("--file=", "", f[1])))
})
source(file.path(allocation_script_dir, "16_allocation_risk.R"))
stopifnot(all(Sys.getenv(c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")) == "1"))
p <- allocation_parameters("pilot")
pilot <- file.path(dirname(allocation_script_dir), "results", "allocation_risk", "pilot")
mf <- readRDS(file.path(pilot, "manifest.rds"))
model_p <- p; model_p$n_allocations <- model_p$n_outcomes <- NULL
stopifnot(identical(mf$params, model_p), identical(mf$R, R.version.string),
          identical(mf$BLAS, sessionInfo()$BLAS))
for (f in names(mf$source_hashes)) {
  stopifnot(identical(unname(mf$source_hashes[[f]]),
    digest::digest(file.path(allocation_script_dir, f), file = TRUE)))
}
for (pkg in names(mf$package_versions)) {
  stopifnot(identical(unname(mf$package_versions[[pkg]]), as.character(packageVersion(pkg))))
}
units <- expand.grid(config = seq_along(p$configs), surface = p$surfaces, rho = p$rho,
  gamma = p$gamma, regime = p$regimes, design = p$designs, stringsAsFactors = FALSE)
selected <- which(units$surface == 1 & units$rho == 0.5 & units$gamma == 0.8 &
                    units$design %in% c(3, 4, 6, 8, 9))
stopifnot(length(selected) == 50)
p$n_outcomes <- 400L
root <- file.path(dirname(pilot), "precision_R400")
cp_dir <- file.path(root, "checkpoints")
dir.create(cp_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- list(parent_manifest = mf, parameters = p, selected = units[selected, ],
                 original_block_ids = selected, driver_hash = digest::digest(
                   file.path(allocation_script_dir, "18_allocation_risk_precision.R"), file = TRUE))
dest_mf <- file.path(root, "manifest.rds")
if (file.exists(dest_mf)) stopifnot(identical(manifest, readRDS(dest_mf))) else saveRDS(manifest, dest_mf)
grid <- build_spatial_grid(p$grid_dim)
setup <- sar_lag_setup(grid$listw_queen)
surfaces <- lapply(p$configs, function(cfg) allocation_surface(cfg, 1, grid, p))
t0 <- proc.time()[[3]]
# One task extends an explicitly validated original fixed-X/scenario/design key.
results <- parallel::mclapply(selected, function(i) {
  cp <- file.path(cp_dir, sprintf("block_%04d.rds", i))
  if (file.exists(cp)) {
    cached <- readRDS(cp); cached$fits <- NULL
    return(cached)
  }
  u <- units[i, ]; cfg <- p$configs[[u$config]]
  old <- readRDS(file.path(pilot, "checkpoints", sprintf("block_%04d.rds", i)))
  key <- old$summary
  stopifnot(key$Incidence_Mode == cfg$mode, key$Rho_Incidence == cfg$rho_x,
    key$Surface == u$surface, key$Neighbor_Type == p$nb_type,
    key$Rho == u$rho, key$Gamma == u$gamma, key$Spillover_Type == u$regime,
    key$True_Tau == p$tau, key$Design == paste("Design", u$design),
    all(old$allocations$N_Outcome == 100), all(old$allocations$Complete))
  out <- allocation_block(u, cfg, surfaces[[u$config]], grid, setup, p, old)
  # Prefix equality is asserted for every allocation, not inferred from seeds.
  stopifnot(identical(names(old$fits), names(out$fits)),
    identical(old$draws, out$draws),
    all(vapply(names(old$fits), function(h) identical(old$fits[[h]], out$fits[[h]][1:100, ]), logical(1))))
  temp <- paste0(cp, ".tmp"); saveRDS(out, temp)
  if (!file.rename(temp, cp)) stop("Could not save precision checkpoint")
  cat(sprintf("precision block %d: %d unique, %.1f s\n", i, nrow(out$allocations), out$elapsed_sec))
  out$fits <- NULL; out
}, mc.cores = 8, mc.preschedule = FALSE, mc.set.seed = FALSE)
bad <- vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(bad)) stop("Precision workers failed for original block IDs: ", paste(selected[bad], collapse = ", "))
a <- do.call(rbind, lapply(results, `[[`, "allocations"))
s <- do.call(rbind, lapply(results, `[[`, "summary"))
draws <- do.call(rbind, lapply(results, `[[`, "draws"))
saveRDS(list(allocations = a, summary = s, draws = draws, params = p,
            original_block_ids = selected), file.path(root, "results.rds"))
write.csv(a, file.path(root, "allocation_metrics.csv"), row.names = FALSE)
write.csv(s, file.path(root, "conditional_summary.csv"), row.names = FALSE)
write.csv(draws, file.path(root, "draw_mapping.csv"), row.names = FALSE)
status <- list(elapsed_sec = proc.time()[[3]] - t0, complete = all(a$Complete),
               n_blocks = nrow(s), total_fits = sum(a$N_Outcome),
               extra_fits = sum(a$N_Outcome) * 0.75,
               warnings = sum(a$N_Warn), aliases = sum(a$N_Aliased))
saveRDS(status, file.path(root, "status.rds")); print(status)
if (!all(a$Complete)) stop("Precision extension had incomplete fits; diagnostics saved")
