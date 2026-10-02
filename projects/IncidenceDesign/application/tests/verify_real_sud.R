# ============================================================
# Script: verify_real_sud.R
# Purpose: Independently reconcile application outputs with every cached source.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud modules
# ============================================================
test_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
real_sud_define_only <- TRUE
source(file.path(test_dir, "..", "code", "run_real_sud.R"))
args <- commandArgs(TRUE); profile <- if (length(args)) args[1] else "production"
root <- file.path(.rs_app, "results", "real_sud_rev_20261002")
s <- readRDS(file.path(root, "setup.rds"))
d <- file.path(root, profile); bundle <- readRDS(file.path(d, "results.rds"))
config <- bundle$config; report <- bundle$summary
u <- if (!is.null(bundle$reporting_units)) bundle$reporting_units else rs_units(config)
keys <- vapply(seq_len(nrow(u)), function(i) rs_block_key(u[i, ], s, config), character(1))
stopifnot(nrow(report) == nrow(u),
  all(vapply(names(u), function(n) identical(report[[n]], u[[n]]), logical(1))),
  identical(report$Source_ID, keys))
csv <- read.csv(file.path(d, "performance.csv"), stringsAsFactors = FALSE, colClasses = c(Year = "character"))
stopifnot(isTRUE(all.equal(csv, report, check.attributes = FALSE, tolerance = 1e-12)))
am <- read.csv(file.path(d, "allocation_metrics.csv"), stringsAsFactors = FALSE)
warnings <- read.csv(file.path(d, "warnings_failures.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(am) == nrow(bundle$allocations))
all_sources <- unique(keys); nf <- nw <- nb <- 0L
for (key in all_sources) {
  b <- readRDS(file.path(d, "checkpoints", paste0(key, ".rds")))
  stopifnot(identical(b$key, key), ncol(b$Z) == b$J, length(b$draws) == b$J,
    nrow(b$Z) == 58, all(b$Z %in% 0:1), sum(b$allocations$Frequency) == b$J,
    setequal(names(b$fits), b$allocations$Allocation_ID))
  hashes <- apply(b$Z, 2, digest::digest, algo = "xxhash64")
  stopifnot(identical(hashes, b$draws))
  for (h in names(b$fits)) {
    f <- b$fits[[h]]
    stopifnot(nrow(f) == b$R, identical(as.integer(f$replicate), seq_len(b$R)))
    m <- rs_allocation_metrics(f)
    old <- b$allocations[b$allocations$Allocation_ID == h, names(m), drop = FALSE]
    stopifnot(isTRUE(all.equal(m, old, check.attributes = FALSE, tolerance = 1e-12)))
    nf <- nf + nrow(f); nw <- nw + sum(nzchar(f$messages)); nb <- nb + m$N_Bound
  }
  recomputed <- rs_distribution(b$allocations, b$summary$Proven_Singleton)
  stopifnot(isTRUE(all.equal(recomputed, b$summary[names(recomputed)], check.attributes = FALSE, tolerance = 1e-12)))
  for (i in which(keys == key)) {
    expected <- rs_reporting_rows(u[i, ], b, s)$summary
    actual <- report[i, names(expected), drop = FALSE]
    stopifnot(isTRUE(all.equal(expected, actual, check.attributes = FALSE, tolerance = 1e-12)))
  }
}
stopifnot(nw == sum(report$N_Warn[!duplicated(keys)]),
  nf == sum(report$N_Independent_Fits[!duplicated(keys)]),
  nb == sum(report$N_Bound[!duplicated(keys)]))
if (all(report$Complete)) stopifnot(nrow(warnings) == nw)
if (profile == "production") rs_manifest(d, s, config)
if (profile == "tail_confirmation") {
  rs_manifest(d, s, bundle$tail_config)
  stopifnot(identical(readRDS(file.path(d, "tail_source_hash.rds")),
    rs_hash_files(file.path(.rs_code_dir, "refine_real_sud_tail.R"))))
}
lines <- c(paste("Profile:", profile), paste("Reporting rows:", nrow(report)),
  paste("Distinct sources:", length(all_sources)), paste("Independent outcome fits:", nf),
  paste("Incomplete rows:", sum(!report$Complete)), paste("Precision-unresolved rows:", sum(!report$Precision_OK)),
  paste("Warning fits:", nw), paste("Boundary fits:", nb),
  "PASS: reporting grid, CSV/RDS, every unique allocation/draw/fit count, conditional metrics,",
  "joint uncertainty and annual diagnostics independently reconciled to cached sources.")
writeLines(lines, file.path(d, "verification.txt"))
cat(paste(lines, collapse = "\n"), "\n")
if (any(!report$Complete)) stop("Incomplete comparisons remain; verification reconciliation passed but study is incomplete")
if (profile == "production" && any(!report$Precision_OK)) stop("Production precision targets remain unresolved")
