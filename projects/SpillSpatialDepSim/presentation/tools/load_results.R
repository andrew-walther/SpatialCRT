# ============================================================
# Script: presentation/tools/load_results.R
# Purpose: Load the Project 1 results the talk uses (per-allocation beta-hat MSE
#          and the corrected summary tables), so figures and slide numbers are
#          computed from saved results rather than typed by hand.
# Author: Andrew Walther
# Created: 2026-10-28
# Dependencies: dplyr
# ============================================================
# Sourced by presentation/tools/make_figures.R (run from the project root) and by
# the deck's setup chunk (rendered from presentation/). All inputs are read-only:
#   results/beta_mse/                          published per-allocation CSVs (frozen)
#   review/results/beta_mse_TrtNoSpill_fixed/  corrected 3x3/3x4 Control Only reruns
#   review/results/beta_mse_TrtSpill_rerun/    the one regenerated Intervention & Control file
#   review/tables/beta_mse_table_<grid>_corrected.csv   summary tables (as in Chapter 2)
# The file-choice rule matches review/code/build_tables_figures.R ("corrected").

suppressMessages(library(dplyr))

# Project root: "." from the project root, ".." from presentation/
p1_root <- if (dir.exists("presentation")) "." else ".."

# BSS allocations (indices into combn(districts, k)), as in the frozen scripts
bss_idx <- list(`2x4` = c(28, 43), `3x3` = c(28, 30, 35, 51, 77, 106), `3x4` = c(313, 612))
grid_info <- data.frame(grid = c("2x4", "3x3", "3x4"), clusters = c(8, 9, 12),
                        n_alloc = c(70, 126, 924), n_bss = c(2, 6, 2), reps = c(50, 50, 10))

spill_labels <- c(TrtNoSpill = "Control Clusters Only", TrtSpill = "Control & Intervention Clusters")

## Per-allocation results ----

#' Path to the corrected per-allocation beta-hat MSE CSV for one scenario.
#'
#' @param spill "TrtNoSpill" (Control Clusters Only) or "TrtSpill" (Control &
#'   Intervention Clusters).
#' @param grid "2x4", "3x3" or "3x4".
#' @param psi,rho Scenario values (psi in 0.5-0.8, rho 0 or 0.01).
#' @return File path; stops if no file exists.
#' @family load_results
csv_path <- function(spill, grid, psi, rho) {
  f <- sprintf("beta_mse_results_%s_%s_A02_B1_P%02d_R%s.csv", spill, grid,
               round(psi * 10), ifelse(rho == 0, "00", "001"))
  dirs <- file.path(p1_root, c(
    if (spill == "TrtNoSpill" && grid != "2x4") "review/results/beta_mse_TrtNoSpill_fixed",
    if (spill == "TrtSpill") "review/results/beta_mse_TrtSpill_rerun",
    "results/beta_mse"))
  hit <- file.path(dirs, f)[file.exists(file.path(dirs, f))]
  if (length(hit) == 0) stop("No results file for ", f)
  hit[1]
}

#' Per-allocation MSE of beta-hat for one grid, every scenario.
#'
#' Every balanced allocation is an SRS allocation (as in the paper); the BSS
#' allocations are flagged. MSE for an allocation = mean over its N replications
#' of (beta-hat - beta)^2, with beta = 1.
#' @param grid Grid label.
#' @return Data frame: grid, spill (label), Psi, Rho, combo, MSE, is_BSS.
#' @family load_results
load_grid <- function(grid) {
  out <- list()
  for (spill in names(spill_labels)) for (psi in c(0.5, 0.6, 0.7, 0.8)) for (rho in c(0.01, 0)) {
    d <- read.csv(csv_path(spill, grid, psi, rho))
    expected <- grid_info$n_alloc[grid_info$grid == grid]
    if (nrow(d) != expected) stop(sprintf("%s: %d rows, expected %d", csv_path(spill, grid, psi, rho), nrow(d), expected))
    out[[length(out) + 1]] <- data.frame(
      grid = grid, spill = unname(spill_labels[spill]), Psi = psi, Rho = rho,
      combo = seq_len(nrow(d)), MSE = d$MSE, is_BSS = seq_len(nrow(d)) %in% bss_idx[[grid]])
  }
  bind_rows(out)
}

all_alloc <- bind_rows(lapply(names(bss_idx), load_grid))

#' Where BSS falls among the random allocations, per scenario.
#'
#' pct_SRS_worse = share of all balanced allocations whose MSE exceeds the mean
#' MSE of the BSS allocations (the BSS allocations themselves are in the
#' denominator, as in the paper's SRS set).
#' @family load_results
rank_tab <- all_alloc %>%
  group_by(grid, spill, Psi, Rho) %>%
  summarise(BSS_mean = mean(MSE[is_BSS]), BSS_max = max(MSE[is_BSS]),
            SRS_mean = mean(MSE), SRS_median = median(MSE), SRS_max = max(MSE),
            pct_SRS_worse = 100 * mean(MSE > mean(MSE[is_BSS])), .groups = "drop")

## Summary tables (as reported in Chapter 2) ----

#' Corrected summary table for one grid (SD = SD of per-allocation MSE, unscaled).
#' @param grid Grid label.
#' @return Data frame with TrtSpill relabelled to the deck's spillover names.
#' @family load_results
load_table <- function(grid) {
  d <- read.csv(file.path(p1_root, sprintf("review/tables/beta_mse_table_%s_corrected.csv", grid)))
  d$grid <- grid
  d$spill <- ifelse(d$TrtSpill == "Control Only", spill_labels[["TrtNoSpill"]], spill_labels[["TrtSpill"]])
  d
}
summary_tab <- bind_rows(lapply(names(bss_idx), load_table))

## Collinearity of treatment and spillover exposure, per allocation ----

#' Rook neighbor list the simulation used, parsed from the frozen script.
#'
#' Each script defines rook_neighbors three times; the last definition is the one
#' the simulation loop uses (2x4: line 310; 3x3 and 3x4: line 351).
#' @param grid Grid label. @return Named list: cluster label -> neighbor labels.
#' @family load_results
sim_neighbors <- function(grid) {
  f <- file.path(p1_root, "code", c(`2x4` = "SpatialSim_NC_DOC.Rmd", `3x3` = "SpatialSim_3x3.Rmd",
                                    `3x4` = "SpatialSim_3x4.Rmd")[[grid]])
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  blocks <- regmatches(txt, gregexpr("rook_neighbors <- list\\((.*?)\\n\\)", txt))[[1]]
  pairs <- regmatches(blocks[length(blocks)], gregexpr('"(\\d+)" = c\\(([0-9, ]+)\\)', blocks[length(blocks)]))[[1]]
  nb <- lapply(pairs, function(p) as.integer(strsplit(sub('.*c\\(([0-9, ]+)\\)', "\\1", p), ",\\s*")[[1]]))
  names(nb) <- sub('"(\\d+)".*', "\\1", pairs)
  nb
}

#' Is z = 1 - x for every cluster (treatment and exposure exactly collinear)?
#'
#' z_i = 1 if cluster i has a treated rook neighbor and, for Control Clusters Only,
#' is itself a control. Collinear allocations leave beta and psi unidentified.
#' @param trt Treated cluster labels. @param nb Neighbor list. @param both TRUE for
#'   Control & Intervention Clusters. @return Logical.
#' @family load_results
is_collinear <- function(trt, nb, both) {
  ids <- as.integer(names(nb)); x <- as.integer(ids %in% trt)
  z <- sapply(seq_along(ids), function(k) as.integer(any(nb[[k]] %in% trt) && (both || x[k] == 0)))
  all(z == 1 - x)
}

#' Per-allocation collinearity flags for one grid.
#' @param grid Grid label. @return Data frame: grid, combo, is_BSS, collinear_CO, collinear_CI.
#' @family load_results
collin_grid <- function(grid) {
  nb <- sim_neighbors(grid); k <- grid_info$clusters[grid_info$grid == grid]
  combos <- combn(seq_len(k), floor(k / 2), simplify = FALSE)
  # all_exposed_CI: under Control & Intervention every cluster has z = 1, so psi is
  # aliased with the intercept (beta is still identified)
  all_exposed <- function(trt) all(sapply(nb, function(n) any(n %in% trt)))
  data.frame(grid = grid, combo = seq_along(combos), is_BSS = seq_along(combos) %in% bss_idx[[grid]],
             collinear_CO = sapply(combos, is_collinear, nb = nb, both = FALSE),
             collinear_CI = sapply(combos, is_collinear, nb = nb, both = TRUE),
             all_exposed_CI = sapply(combos, all_exposed))
}
collin <- bind_rows(lapply(names(bss_idx), collin_grid))

# Cross-check: on the 2x4 grid the collinear Control-Only allocations are exactly the
# ones whose psi-hat could not be estimated (NA) in the published results
psi24 <- read.csv(file.path(p1_root, "results/psi_mse/psi_mse_results_TrtNoSpill_2x4_A02_B1_P05_R00.csv"))
if (!identical(is.na(psi24$MSE), collin$collinear_CO[collin$grid == "2x4"]))
  stop("Collinearity flags do not match the NA psi-hat allocations in results/psi_mse (2x4)")

# Cross-check: the per-allocation means must reproduce the summary tables
chk <- inner_join(rank_tab, summary_tab, by = c("grid", "spill", "Psi", "Rho"))
if (nrow(chk) != 48 || max(abs(chk$SRS_mean - chk$SRS_MSE)) > 1e-9 || max(abs(chk$BSS_mean - chk$BSS_MSE)) > 1e-9)
  stop("Per-allocation results do not reproduce review/tables (", nrow(chk), " scenarios matched)")
