# ============================================================
# Script: refine_real_sud_tail.R
# Purpose: Refine conditional-risk estimates at prespecified queen corners.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud modules
# ============================================================
real_sud_define_only <- TRUE
.rs_tail_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(normalizePath(dirname(f)))
  }
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
})
source(file.path(.rs_tail_dir, "run_real_sud.R"))

#' Measure selection stability using independent outcome halves
#' @param b Cached block with unique allocation metrics and draw frequencies.
#' @return Cross-selected upper-tail means and tail draw count, or NA if incomplete.
#' @family real_sud_tail
#' @seealso allocation_distribution
#' @examples
#' # rs_cross_tail(b)
rs_cross_tail <- function(b) {
  a <- b$allocations
  n <- sum(a$Frequency); top_n <- max(1L, ceiling(n / 10))
  if (!all(a$Complete)) return(data.frame(Tail_Draws = top_n,
    Worst10_Select1_Evaluate2 = NA_real_, Worst10_Select2_Evaluate1 = NA_real_))
  ix <- rep(seq_len(nrow(a)), a$Frequency)
  h1 <- a$MSE_Half1[ix]; h2 <- a$MSE_Half2[ix]
  data.frame(Tail_Draws = top_n,
    Worst10_Select1_Evaluate2 = mean(h2[order(h1, decreasing = TRUE)[seq_len(top_n)]]),
    Worst10_Select2_Evaluate1 = mean(h1[order(h2, decreasing = TRUE)[seq_len(top_n)]]))
}

#' Extend R at fixed sampled allocations, preserving the completed main study
#' @param workers Forked worker count, with explicit seeded streams.
#' @param root Study root, requiring completed production outputs.
#' @return Tail-confirmation reporting data frame, saved separately.
#' @family real_sud_tail
#' @seealso rs_block, rs_cross_tail
#' @examples
#' # rs_refine_tail(workers = 8)
rs_refine_tail <- function(workers = 8L, root = file.path(.rs_app, "results", "real_sud_rev_20261002")) {
  setup <- readRDS(file.path(root, "setup.rds"))
  main <- readRDS(file.path(root, "production", "results.rds"))
  if (any(!main$summary$Complete) || any(!main$summary$Precision_OK)) stop("Complete/precision-valid main results required before tail refinement")
  config <- main$config
  rs_manifest(file.path(root, "production"), setup, config)
  outdir <- file.path(root, "tail_confirmation")
  # The same main streams are extended: this is outcome refinement, not a new
  # independent confirmatory allocation sample. State that distinction in prose.
  tail_config <- config; tail_config$profile <- "tail_confirmation"
  tail_config$selection <- "education queen mean_rank rho 0/0.5 gamma 0.8 both regimes all years/designs"
  tail_config$minimum_R <- 400L
  rs_manifest(outdir, setup, tail_config)
  extra_path <- file.path(outdir, "tail_source_hash.rds")
  extra <- rs_hash_files(file.path(.rs_tail_dir, "refine_real_sud_tail.R"))
  if (file.exists(extra_path) && !identical(readRDS(extra_path), extra)) stop("Tail source manifest mismatch")
  if (!file.exists(extra_path)) saveRDS(extra, extra_path)
  cp <- file.path(outdir, "checkpoints"); dir.create(cp, recursive = TRUE, showWarnings = FALSE)
  u <- subset(main$summary, Model == "education" & Neighbor == "queen" & Summary == "mean_rank" & Rho %in% c(0, 0.5) & Gamma == 0.8)
  units <- u[names(rs_units(config))]; keys <- u$Source_ID
  unique_keys <- unique(keys)
  cat(sprintf("Tail refinement: %d reporting rows, %d sources, R >= 400 at fixed main-study allocations\n", nrow(u), length(unique_keys)))
  status <- parallel::mclapply(unique_keys, function(key) {
    path <- file.path(cp, paste0(key, ".rds"))
    old <- if (file.exists(path)) readRDS(path) else readRDS(file.path(root, "production", "checkpoints", paste0(key, ".rds")))
    i <- match(key, keys)
    b <- rs_block(units[i, ], setup, config, old$J, max(old$R, 400L), old)
    stopifnot(identical(old$Z, b$Z), all(vapply(names(old$fits), function(h) identical(old$fits[[h]], b$fits[[h]][seq_len(old$R), ]), logical(1))))
    saveRDS(b, path); cat(key, "J", b$J, "R", b$R, "complete", b$summary$Complete, "\n")
    key
  }, mc.cores = workers, mc.set.seed = FALSE)
  if (any(vapply(status, inherits, logical(1), "try-error"))) stop("Tail worker failure; preserve caches and inspect log")
  summaries <- allocations <- vector("list", nrow(units)); warns <- list()
  for (key in unique_keys) {
    b <- readRDS(file.path(cp, paste0(key, ".rds")))
    for (i in which(keys == key)) {
      r <- rs_reporting_rows(units[i, ], b, setup)
      summaries[[i]] <- cbind(r$summary, rs_cross_tail(b), Tail_Refinement = "Same allocations, extended independent outcome prefix")
      allocations[[i]] <- r$allocations
    }
    for (h in names(b$fits)) {
      f <- b$fits[[h]]; f <- f[nzchar(f$messages) | !is.finite(f$error) | !is.finite(f$se), ]
      if (nrow(f)) warns[[length(warns) + 1L]] <- cbind(Source_ID = key, Allocation_ID = h, f)
    }
  }
  report <- do.call(rbind, summaries); rownames(report) <- NULL
  report$Source_Reporting_Uses <- as.integer(table(keys)[report$Source_ID])
  am <- do.call(rbind, allocations); rownames(am) <- NULL
  warnings <- if (length(warns)) do.call(rbind, warns) else
    read.csv(file.path(root, "production", "warnings_failures.csv"), stringsAsFactors = FALSE)[0, ]
  write.csv(report, file.path(outdir, "performance.csv"), row.names = FALSE)
  write.csv(am, file.path(outdir, "allocation_metrics.csv"), row.names = FALSE)
  write.csv(warnings, file.path(outdir, "warnings_failures.csv"), row.names = FALSE)
  saveRDS(list(summary = report, allocations = am, config = config,
    reporting_units = units, tail_config = tail_config), file.path(outdir, "results.rds"))
  if (any(!report$Complete)) stop("Tail comparisons incomplete; inspect explicit warning/failure outputs")
  cat("Tail refinement saved separately; main study outputs preserved.\n")
  invisible(report)
}
if (!exists("real_sud_tail_define_only") || !isTRUE(real_sud_tail_define_only)) {
  args <- commandArgs(TRUE)
  rs_refine_tail(if (length(args)) as.integer(args[1]) else 8L)
}
