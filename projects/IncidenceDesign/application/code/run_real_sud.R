# ============================================================
# Script: run_real_sud.R
# Purpose: Execute manifest-checked observed-incidence NC application profiles.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: real_sud_simulation.R, parallel
# ============================================================

# Entry point and shared helpers ----
.rs_run_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(p)) return(normalizePath(dirname(p)))
  }
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else stop("Run or source run_real_sud.R by path")
})
source(file.path(.rs_run_dir, "real_sud_simulation.R"))

#' Attach yearly budget diagnostics to canonical performance results
#' @param u One yearly reporting block.
#' @param b Cached canonical result; equivalent years share Source_ID.
#' @param setup Named annual populations and incidence.
#' @return Summary row and allocation rows with actual annual budget diagnostics.
#' @family real_sud
#' @seealso rs_block_key
#' @examples
#' # rs_reporting_rows(u, b, setup)
rs_reporting_rows <- function(u, b, setup) {
  a <- setup$annual[[u$Year]]; m <- b$allocations
  for (i in seq_len(nrow(m))) {
    z <- b$Z[, match(m$Allocation_ID[i], b$draws)]
    m$Population_Share[i] <- sum(a$person_years * z) / sum(a$person_years)
    m$Incidence_Difference[i] <- mean(a$rate_per_100k[z == 1]) - mean(a$rate_per_100k[z == 0])
  }
  w <- m$Frequency / sum(m$Frequency)
  s <- b$summary
  s$Mean_Treated <- sum(w * m$Treated)
  s$Min_Treated <- min(m$Treated); s$Max_Treated <- max(m$Treated)
  s$Mean_Population_Share <- sum(w * m$Population_Share)
  s$Min_Population_Share <- min(m$Population_Share); s$Max_Population_Share <- max(m$Population_Share)
  s$Mean_Incidence_Difference <- sum(w * m$Incidence_Difference)
  s$N_Unidentified_Allocations <- sum(!m$Identified)
  s$Fail_Rate_Est <- 1 - s$N_Valid_Est / s$N_Independent_Fits
  s$Fail_Rate_CI <- 1 - s$N_Valid_CI / s$N_Independent_Fits
  list(summary = cbind(u, s), allocations = cbind(u[rep(1, nrow(m)), ], Source_ID = b$key, m))
}

#' Run each distinct distribution once, without treating reused years as replication
#' @param profile smoke, pilot or production.
#' @param workers Number of forked workers; each has separate key-seeded RNG.
#' @param root Study root containing frozen setup and separate profile directories.
#' @return Reporting summary; checkpoints retain every warning and failed fit.
#' @family real_sud
#' @seealso rs_manifest, rs_run_block, rs_reporting_rows
#' @examples
#' # rs_run("smoke", workers = 2)
rs_run <- function(profile = "smoke", workers = 2L,
                   root = file.path(.rs_app, "results", "real_sud_rev_20261002")) {
  config <- rs_config(profile); setup <- rs_setup(root)
  outdir <- file.path(root, profile)
  rs_manifest(outdir, setup, config)
  cp <- file.path(outdir, "checkpoints")
  dir.create(cp, recursive = TRUE, showWarnings = FALSE)
  u <- rs_units(config)
  keys <- vapply(seq_len(nrow(u)), function(i) rs_block_key(u[i, ], setup, config), character(1))
  unique_keys <- unique(keys); canonical <- match(unique_keys, keys)
  cat(sprintf("%s: %d reporting blocks, %d unique distributions, %d workers\n", profile, nrow(u), length(canonical), workers))
  flush.console()
  # BLAS is limited by environment in the command below. On Longleaf use one
  # process per allocated core; checkpoints can occupy several GB at large tiers.
  # Worker return values are small; full fit caches stay on disk and are read once
  # per source while reporting. mc.set.seed=FALSE preserves explicit key streams.
  status <- parallel::mclapply(canonical, function(i) {
    key <- keys[i]
    b <- rs_run_block(u[i, ], setup, config, file.path(cp, paste0(key, ".rds")))
    cat(sprintf("%s %s %s %s rho=%.1f gamma=%.1f: J=%d R=%d complete=%s precision=%s\n",
      u$Year[i], u$Model[i], u$Neighbor[i], u$Design[i], u$Rho[i], u$Gamma[i],
      b$J, b$R, b$summary$Complete, b$summary$Precision_OK))
    flush.console()
    key
  }, mc.cores = workers, mc.set.seed = FALSE)
  failed <- vapply(status, inherits, logical(1), "try-error")
  if (any(failed)) stop("Worker failure(s); inspect profile log and preserved checkpoints")
  summaries <- allocations <- vector("list", nrow(u)); warning_rows <- list()
  for (key in unique_keys) {
    b <- readRDS(file.path(cp, paste0(key, ".rds")))
    for (i in which(keys == key)) {
      r <- rs_reporting_rows(u[i, ], b, setup)
      summaries[[i]] <- r$summary; allocations[[i]] <- r$allocations
    }
    for (h in names(b$fits)) {
      f <- b$fits[[h]]
      f <- f[nzchar(f$messages) | !is.finite(f$error) | !is.finite(f$se), ]
      if (nrow(f)) warning_rows[[length(warning_rows) + 1L]] <- cbind(Source_ID = key, Allocation_ID = h, f)
    }
  }
  s <- do.call(rbind, summaries); rownames(s) <- NULL
  s$Source_Reporting_Uses <- as.integer(table(keys)[s$Source_ID])
  am <- do.call(rbind, allocations); rownames(am) <- NULL
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(Source_ID = character(), Allocation_ID = character(), replicate = integer(),
      tau = numeric(), error = numeric(), se = numeric(), alias = logical(), messages = character(), identified = logical())
  stopifnot(nrow(s) == nrow(u), !anyDuplicated(s[names(u)]), setequal(s$Design_ID, config$designs))
  write.csv(s, file.path(outdir, "performance.csv"), row.names = FALSE)
  write.csv(am, file.path(outdir, "allocation_metrics.csv"), row.names = FALSE)
  write.csv(warnings, file.path(outdir, "warnings_failures.csv"), row.names = FALSE)
  saveRDS(list(summary = s, allocations = am, config = config), file.path(outdir, "results.rds"))
  cat(sprintf("Saved %d rows: %d incomplete, %d unresolved precision; %d warning/failure fit rows.\n",
    nrow(s), sum(!s$Complete), sum(!s$Precision_OK), nrow(warnings)))
  invisible(s)
}

# Source for tests without executing; use Rscript for the profile entry point.
if (!exists("real_sud_define_only") || !isTRUE(real_sud_define_only)) {
  args <- commandArgs(TRUE)
  profile <- if (length(args)) args[1] else "smoke"
  workers <- if (length(args) > 1) as.integer(args[2]) else 2L
  rs_run(profile, workers)
}
