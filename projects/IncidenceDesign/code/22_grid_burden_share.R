# ============================================================
# Script: 22_grid_burden_share.R
# Purpose: Burden share treated by design on the stored 10x10 grid incidence
#   surfaces, to test whether Incidence-Guided Saturation out-targets Plain
#   Saturation when incidence differs between quadrants (author decision 2026-10-21).
# Author: Andrew Walther
# Created: 2026-10-21
# Dependencies: 05_run_simulation.R (sourced with sim_define_only = TRUE; sources
#   01-04), digest
# ============================================================
# Run from code/ (no fitting; about a minute):
#   Rscript 22_grid_burden_share.R
# Outputs: results/grid_burden_share/ (small CSVs):
#   grid_burden_share_allocations.csv  one row per (config, nb, block, design, surface)
#   grid_burden_share_summary.csv      mean and 10th/90th percentiles by config x nb x design
#   grid_quadrant_spread.csv           per surface: quadrant mean-rank spread and the
#                                      Incidence-Guided minus Plain Saturation gain
#   grid_burden_share_log.txt
#
# Definitions (match the NC A1 analysis, methods_ext_stored.R):
#   rate share  = sum_j X_j z_j / sum_j X_j   (grid clusters have equal population,
#                 so this is also the expected death share)
#   quadrant spread = (max - min of the four quadrant mean incidence ranks) / N,
#                 ranks 1..N with N = 100 (NC: (max - min) / 58 over its 4 regions)
# The allocations are exactly the production ones: generate_surfaces() and
# draw_assignments() are called with the production keys for every queen and rook
# block, so nothing here draws outside set_seed_key(). Neighbor types are reported
# separately (queen primary); incidence configs are never pooled.

# Setup ----
sim_define_only <- TRUE
sim_profile <- "full"
sim_script_dir <- normalizePath(".")
source("05_run_simulation.R")

out_dir <- file.path(dirname(script_dir), "results", "grid_burden_share")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Quadrant ids, exactly as designs 3 and 8 in 03_designs.R
half <- coords$x[which.max(coords$x)] / 2
q_id <- ifelse(coords$x <= half & coords$y <= half, 1,
        ifelse(coords$x > half  & coords$y <= half, 2,
        ifelse(coords$x <= half & coords$y > half,  3, 4)))

#' Rate share treated for each column of an allocation matrix
#'
#' share_i = sum_j X_j Z_ji / sum_j X_j, the share of incidence in treated clusters.
#' @param Zk N x J allocation matrix for one surface.
#' @param x Incidence surface (length N).
#' @return Numeric vector of length J.
#' @examples
#' rate_share(matrix(c(1, 0, 0, 1), 2), c(3, 1))   # 0.75, 0.25
rate_share <- function(Zk, x) as.vector(crossprod(Zk, x)) / sum(x)

# Main loop ----
rows <- list(); spread <- list()
for (ci in seq_along(inc_configs)) {
  cfg <- inc_configs[[ci]]
  X <- generate_surfaces(cfg)
  if (any(X < 0)) stop("Negative incidence in config ", cfg$mode, " ", cfg$rho_x, ": share undefined")
  lab <- if (cfg$mode == "iid") "iid" else sprintf("%s rho_X=%.1f", cfg$mode, cfg$rho_x)
  for (k in seq_len(ncol(X))) {
    r <- rank(X[, k], ties.method = "average"); qm <- tapply(r, q_id, mean)
    spread[[length(spread) + 1L]] <- data.frame(Config = lab, Surface = k,
      Quadrant_Rank_Min = min(qm), Quadrant_Rank_Max = max(qm), Spread = (max(qm) - min(qm)) / N)
  }
  for (nb in neighbor_types) for (rho in rho_vals) for (gm in gamma_vals) for (sp in spill_types) {
    for (d in design_ids) {
      Z <- draw_assignments(cfg, nb, rho, gm, sp, d, X)
      for (k in seq_len(ncol(X))) {
        cols <- (k - 1) * n_design_draw + seq_len(n_design_draw)
        s <- rate_share(Z[, cols, drop = FALSE], X[, k])
        rows[[length(rows) + 1L]] <- data.frame(Config = lab, Neighbor = nb, Rho = rho, Gamma = gm,
          Regime = sp, Design_ID = d, Surface = k, Share_Mean = mean(s), Share_P10 = unname(stats::quantile(s, 0.1)),
          Share_P90 = unname(stats::quantile(s, 0.9)), Treated_Fraction = mean(colSums(Z[, cols, drop = FALSE])) / N)
      }
    }
  }
}
al <- do.call(rbind, rows); sp_tab <- do.call(rbind, spread)
names_map <- c("Checkerboard", "High Incidence Focus", "Plain Saturation", "Isolation Buffer",
               "2x2 Blocking", "Balanced Quartiles", "Balanced Halves", "Incidence-Guided Saturation",
               "Simple Random Sampling")
al$Design <- names_map[al$Design_ID]
write.csv(al, file.path(out_dir, "grid_burden_share_allocations.csv"), row.names = FALSE)

# Summaries ----
sm <- aggregate(cbind(Share_Mean, Share_P10, Share_P90, Treated_Fraction) ~ Config + Neighbor + Design_ID + Design, al, mean)
sm <- sm[order(sm$Neighbor != "queen", sm$Config, -sm$Share_Mean), ]
write.csv(sm, file.path(out_dir, "grid_burden_share_summary.csv"), row.names = FALSE)

# Per-surface gain of Incidence-Guided over Plain Saturation (queen, all blocks)
q <- al[al$Neighbor == "queen" & al$Design_ID %in% c(3, 8), ]
g <- aggregate(Share_Mean ~ Config + Surface + Design_ID, q, mean)
g <- merge(g[g$Design_ID == 8, c("Config", "Surface", "Share_Mean")], g[g$Design_ID == 3, c("Config", "Surface", "Share_Mean")],
           by = c("Config", "Surface"), suffixes = c("_IGS", "_Plain"))
g$Gain <- g$Share_Mean_IGS - g$Share_Mean_Plain
sp_tab <- merge(sp_tab, g, by = c("Config", "Surface"))
write.csv(sp_tab, file.path(out_dir, "grid_quadrant_spread.csv"), row.names = FALSE)

# Log ----
log <- c(sprintf("Run %s; production surfaces (K = %d) and allocations (J = %d per surface) for %d configs x %d neighbor types x %d blocks",
                 format(Sys.time()), n_surfaces, n_design_draw, length(inc_configs), length(neighbor_types),
                 length(rho_vals) * length(gamma_vals) * length(spill_types)))
for (nb in neighbor_types) for (cf in unique(sm$Config)) {
  d <- sm[sm$Neighbor == nb & sm$Config == cf, ]
  log <- c(log, sprintf("%s, %s: %s", nb, cf, paste(sprintf("%s %.3f", d$Design, d$Share_Mean), collapse = "; ")))
}
for (cf in unique(sp_tab$Config)) {
  d <- sp_tab[sp_tab$Config == cf, ]
  log <- c(log, sprintf("Quadrant spread %s: mean %.3f (range %.3f-%.3f); IGS - Plain gain mean %+.3f (range %+.3f to %+.3f); corr(spread, gain) %.2f",
    cf, mean(d$Spread), min(d$Spread), max(d$Spread), mean(d$Gain), min(d$Gain), max(d$Gain),
    if (stats::sd(d$Spread) > 0) stats::cor(d$Spread, d$Gain) else NA_real_))
}
log <- c(log, "NC reference (A1, frozen regions): region mean-rank spread (max - min)/58 = 0.10-0.14 per year; IGS rate share 0.511 vs Plain 0.504")
writeLines(log, file.path(out_dir, "grid_burden_share_log.txt"))
cat(log, sep = "\n")
