# ============================================================
# Script: build_cluster_weights.R
# Purpose: Build and save the queen and rook contiguity weights (W) for the 58 NC
#          community college clusters.
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: sf, spdep, tigris, dplyr; sources run_application_profiles.R
# ============================================================
#
# Run from anywhere:  Rscript projects/IncidenceDesign/application/code/build_cluster_weights.R
#
# The clusters are the 100 NC counties (tigris boundaries) dissolved by primary
# community college (build_nc_application_clusters()). Contiguity comes from shared
# polygon borders:
#   queen: clusters i and j are neighbors if their borders share any point
#   rook:  ... only if they share a border segment (rook is a subset of queen)
# For each type the output holds:
#   nb         spdep neighbor list
#   adjacency  binary matrix A, A_ij = 1 if j neighbors i (symmetric, zero diagonal)
#   W          row-standardized weights, W_ij = A_ij / sum_j A_ij (rows sum to 1)
# Rows and columns are named by Primary_College. Join other cluster data by name,
# not by position. The weights are public geography, so the output is tracked in git:
#   application/data/nc_cluster_weights.rds

# Locate this script: innermost source() frame (when sourced), else Rscript --file ----
code_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(normalizePath(dirname(ofile)))
  }
  file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(getwd())
})

# Provides build_nc_application_clusters() and build_application_spatial_weights()
source(file.path(code_dir, "run_application_profiles.R"))

# Legal county boundaries (cb = FALSE), not the shoreline-clipped cartographic ones
# (cb = TRUE), so counties that meet across water are neighbors. Author decision
# (2026-09-25): an education intervention can plausibly reach nearby counties across a
# river or sound. This adds 3 pairs relative to cb = TRUE: College of The Albemarle-Martin
# CC (Chowan River), Carteret-Pamlico (Neuse River), Carteret-Beaufort County CC (Pamlico
# Sound).
clusters <- build_nc_application_clusters(cb = FALSE)$clusters

# Build both contiguity types with the existing helper (zero diagonal, row-standardized) ----
weights <- lapply(c(queen = TRUE, rook = FALSE), function(queen) {
  w <- build_application_spatial_weights(clusters, queen = queen)
  W <- w$W
  dimnames(W) <- list(clusters$Primary_College, clusters$Primary_College)
  list(nb = w$nb, adjacency = (W > 0) * 1, W = W)
})

out <- c(weights, list(clusters = sf::st_drop_geometry(clusters)[c("ID", "Primary_College", "n_counties")]))
out_file <- file.path(code_dir, "..", "data", "nc_cluster_weights.rds")
saveRDS(out, out_file)

# Summary: neighbor counts, and the pairs that are neighbors under queen only ----
n_nb <- sapply(weights, function(w) rowSums(w$adjacency))
queen_only <- which(weights$queen$adjacency == 1 & weights$rook$adjacency == 0, arr.ind = TRUE)
queen_only <- queen_only[queen_only[, 1] < queen_only[, 2], , drop = FALSE]
message(sprintf("Saved %s\n  neighbors per cluster: queen %d-%d (mean %.2f), rook %d-%d (mean %.2f)",
                normalizePath(out_file),
                min(n_nb[, "queen"]), max(n_nb[, "queen"]), mean(n_nb[, "queen"]),
                min(n_nb[, "rook"]), max(n_nb[, "rook"]), mean(n_nb[, "rook"])))
message("  queen-only pairs (touch at a point): ",
        if (nrow(queen_only) == 0) "none" else
          paste(clusters$Primary_College[queen_only[, 1]], "--", clusters$Primary_College[queen_only[, 2]], collapse = "; "))
