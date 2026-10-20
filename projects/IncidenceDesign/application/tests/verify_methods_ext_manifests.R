# ============================================================
# Script: verify_methods_ext_manifests.R
# Purpose: Check V9: the hashed production computation files are unchanged after
#   the methods-extension work (NC production and grid simulation manifests).
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: digest
# ============================================================
# Run from the project root:  Rscript application/tests/verify_methods_ext_manifests.R
# Exits with status 1 if any recorded hash differs from the file on disk.

root <- normalizePath(".")
out <- character(0); fail <- FALSE

# NC production manifest: SHA-256 of each source file, keyed by absolute path
m <- readRDS(file.path(root, "application/results/real_sud_rev_20261002/production/manifest.rds"))
for (p in names(m$source)) {
  now <- digest::digest(p, file = TRUE, algo = "sha256"); ok <- identical(now, unname(m$source[[p]]))
  out <- c(out, sprintf("%-4s NC production  %s", if (ok) "OK" else "DIFF", sub(paste0(root, "/"), "", p)))
  fail <- fail || !ok
}

# Grid simulation manifest: code hashes for 01-05 (algorithm detected from the record)
g <- readRDS(file.path(root, "results/checkpoints/rev_2026-09/full/manifest.rds"))
for (f in names(g$code_hashes)) {
  path <- file.path(root, "code", f)
  cand <- vapply(c("md5", "sha256", "xxhash64"), function(a) digest::digest(path, file = TRUE, algo = a), character(1))
  ok <- unname(g$code_hashes[[f]]) %in% cand
  out <- c(out, sprintf("%-4s grid full run   code/%s", if (ok) "OK" else "DIFF", f))
  fail <- fail || !ok
}

out <- c(sprintf("V9 run %s: %s", format(Sys.time()), if (fail) "FAIL" else "PASS"), out)
cat(out, sep = "\n")
dir.create(file.path(root, "application/results/methods_ext_2026-10/verification"), recursive = TRUE, showWarnings = FALSE)
writeLines(out, file.path(root, "application/results/methods_ext_2026-10/verification/v9_manifests.txt"))
if (fail) quit(status = 1)
