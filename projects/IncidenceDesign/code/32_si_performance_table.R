# ============================================================
# Script: 32_si_performance_table.R
# Purpose: Build the body rows of the supplement table "performance of all eight
#          designs and SRS" from saved results, after reproducing every cell of the
#          existing Tables S7 (tab:performance) and S12 (tab:eight-design).
# Author: Andrew Walther
# Created: 2026-11-25
# Dependencies: base R only (reads results/sim_data and results/grid_tdf_2026-10)
# Output: results/si_tables/performance_nine.tex (body rows) and performance_nine_check.txt
# ============================================================
# Run from the project root:  Rscript code/32_si_performance_table.R
# Reads saved results only; no simulation is run.

# Setup ----

res_dir <- "results"
out_dir <- file.path(res_dir, "si_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Inputs ----
# Oracle production run (September 27) and the t + df re-fit (Coverage, Power).
sim <- readRDS(file.path(res_dir, "sim_data",
                         "sim_results_MLE_tau_sweep_combined_20260927_000502.rds"))
tdf <- readRDS(file.path(res_dir, "grid_tdf_2026-10", "grid_tdf_scenarios.rds"))
# Surface-level oracle results of the same run, for the Monte Carlo SE of pooled MSE
surf <- readRDS(file.path(res_dir, "sim_data",
                          "surface_results_MLE_tau_sweep_20260927_000502.rds"))

keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design", "Rho",
          "Gamma", "Spillover_Type", "True_Tau", "Estimator")
stopifnot(all(sim$Estimator == "oracle"))
idx <- match(do.call(paste, sim[keys]), do.call(paste, tdf[keys]))
if (anyNA(idx)) stop(sum(is.na(idx)), " scenario rows have no t + df re-fit")
if (max(abs(sim$MSE - tdf$MSE[idx])) > 1e-12) stop("t + df re-fit MSE differs")
# Reported interval is t + df (AGENTS.md "Intervals everywhere")
sim$Coverage <- tdf$Coverage_Tdf[idx]
sim$Power    <- tdf$Power_Tdf[idx]

## Design labels ----
design_name <- c(`Design 1` = "Checkerboard", `Design 2` = "High Incidence Focus",
                 `Design 3` = "Saturation Quadrants", `Design 4` = "Isolation Buffer",
                 `Design 5` = "2$\\times$2 Blocking", `Design 6` = "Balanced Quartiles",
                 `Design 7` = "Balanced Halves",
                 `Design 8` = "Incidence-Guided Saturation",
                 `Design 9` = "Simple Random Sampling")
retained <- c("Design 1", "Design 2", "Design 4", "Design 5", "Design 6", "Design 8")
srs <- "Design 9"

# Functions ----

#' Average within-block MSE rank
#'
#' Blocks are config x rho x psi (Gamma) x regime (Spillover_Type) at one neighbor
#' type and beta = 1; each block ranks the designs in `ids` by MSE (1 = lowest,
#' ties averaged), and the ranks are averaged over blocks.
#' @param d Scenario rows for one neighbor type at True_Tau == 1
#' @param ids Design IDs entering the ranking
#' @return Named numeric vector of average ranks
avg_rank <- function(d, ids) {
  d <- d[d$Design %in% ids, ]
  blk <- do.call(paste, d[c("Incidence_Mode", "Rho_Incidence", "Rho", "Gamma",
                            "Spillover_Type")])
  stopifnot(all(table(blk) == length(ids)))
  r <- ave(d$MSE, blk, FUN = function(v) rank(v, ties.method = "average"))
  tapply(r, d$Design, mean)
}

#' Per-design summary for one neighbor type at beta = 1
#'
#' MSE, Bias, Coverage, Power, Mean_Treated are means over the 160 scenarios.
#' SE of MSE = sd(m_1..m_50) / sqrt(50), where m_u is the design's mean MSE over the
#' blocks of (configuration, surface) unit u (as code/15_srs_benchmark.R, Table S7).
#' @param nb "queen" or "rook"
#' @return Data frame, one row per design, with Rank7 (six retained designs + SRS)
#'   and Rank9 (all nine designs)
summarize_nb <- function(nb) {
  d <- sim[sim$Neighbor_Type == nb & sim$True_Tau == 1, ]
  stopifnot(nrow(d) == 9 * 160)
  s <- do.call(rbind, lapply(split(d, d$Design), function(x) {
    stopifnot(nrow(x) == 160)
    data.frame(Design = x$Design[1], MSE = mean(x$MSE), SE100 = 100 * mean(x$SE_MSE),
               Bias = mean(x$Bias), Coverage = mean(x$Coverage),
               Power = mean(x$Power), Treated = mean(x$Mean_Treated))
  }))
  # Monte Carlo SE of pooled MSE from the 50 independent configuration-by-surface units
  su <- surf[surf$Neighbor_Type == nb & surf$True_Tau == 1 & surf$Estimator == "oracle", ]
  unit_mean <- aggregate(MSE ~ Design + Incidence_Mode + Rho_Incidence + Surface, su, mean)
  n_units <- tapply(unit_mean$MSE, unit_mean$Design, length)
  stopifnot(all(n_units == 50))
  se <- tapply(unit_mean$MSE, unit_mean$Design, function(m) sd(m) / sqrt(length(m)))
  s$SE100 <- 100 * se[s$Design]
  s$Rank7 <- avg_rank(d, c(retained, srs))[s$Design]
  s$Rank9 <- avg_rank(d, names(design_name))[s$Design]
  s$Rank7[!(s$Design %in% c(retained, srs))] <- NA  # the two variants are outside that comparison
  s
}

fmt3 <- function(x) ifelse(round(x, 3) < 0, sprintf("$-%.3f$", abs(x)), sprintf("%.3f", x))
fmt_treated <- function(x) ifelse(abs(x - 50) < 1e-9, "50", sprintf("%.1f", x))

# Summaries ----
S <- list(queen = summarize_nb("queen"), rook = summarize_nb("rook"))
# Row order: queen MSE, applied to both panels (as in S12)
ord <- S$queen$Design[order(S$queen$MSE)]

# Table body ----
panel <- function(nb, title) {
  s <- S[[nb]][match(ord, S[[nb]]$Design), ]
  nm <- design_name[s$Design]
  if (nb == "rook") nm[s$Design == "Design 1"] <- "Checkerboard (flagged)"
  r7 <- ifelse(is.na(s$Rank7), "--", fmt3(s$Rank7))
  c(sprintf("\\multicolumn{8}{l}{\\textit{%s}} \\\\", title),
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", nm, fmt3(s$MSE), fmt3(s$SE100),
            fmt3(s$Bias), fmt3(s$Coverage), fmt3(s$Power), r7, fmt_treated(s$Treated)))
}
body <- c(panel("queen", "Queen contiguity (primary)"), "", "\\midrule",
          panel("rook", "Rook contiguity (sensitivity)"))
writeLines(body, file.path(out_dir, "performance_nine.tex"))

# Reproduction check against the printed S7 and S12 ----
si <- readLines("paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex")

#' Pull the printed rows of one labelled table's tabular body
#' @param label LaTeX label of the table
#' @return Character vector of the rows (split into queen, rook by \midrule)
printed_rows <- function(label) {
  i <- grep(paste0("\\\\label\\{", label, "\\}"), si)
  stopifnot(length(i) == 1)
  b <- grep("\\\\begin\\{tabular\\}", si); b <- b[b > i][1]
  e <- grep("\\\\end\\{tabular\\}", si);   e <- e[e > b][1]
  x <- si[(b + 1):(e - 1)]
  x <- x[grepl("\\\\\\\\$", x) & !grepl("multicolumn|^Design &", x)]
  trimws(sub("\\s*\\\\\\\\$", "", x))
}
# Split a row into (name, cells); print the stored cell strings exactly
split_row <- function(r) strsplit(r, " & ", fixed = TRUE)[[1]]

lines <- c("Reproduction check of Tables S7 (tab:performance) and S12 (tab:eight-design)",
           "Source: results/sim_data oracle 20260927_000502 + grid_tdf_scenarios.rds (t + df)",
           "Cell format: panel | design | column | printed | computed | status", "")
bad <- 0
chk <- function(panel, des, col, printed, computed) {
  ok <- identical(trimws(printed), computed)
  if (!ok) bad <<- bad + 1
  lines <<- c(lines, sprintf("%s | %s | %s | %s | %s | %s", panel, des, col, printed,
                             computed, if (ok) "OK" else "MISMATCH"))
}

## S7: six designs x 7 numeric columns x 2 panels ----
p7 <- printed_rows("tab:performance"); stopifnot(length(p7) == 12)
cols7 <- c("MSE", "SE100", "Bias", "Coverage", "Power", "Rank7", "Treated")
n7 <- 0
for (k in seq_along(p7)) {
  nb <- if (k <= 6) "queen" else "rook"
  cells <- split_row(p7[k]); nm <- cells[1]
  des <- names(design_name)[match(sub(" \\(flagged\\)", "", nm), design_name)]
  stopifnot(!is.na(des))
  s <- S[[nb]][S[[nb]]$Design == des, ]
  comp <- c(fmt3(s$MSE), fmt3(s$SE100), fmt3(s$Bias), fmt3(s$Coverage), fmt3(s$Power),
            fmt3(s$Rank7), fmt_treated(s$Treated))
  stopifnot(length(cells) == 8)
  for (j in 1:7) { chk(paste0("S7 ", nb), nm, cols7[j], cells[j + 1], comp[j]); n7 <- n7 + 1 }
}

## S12: nine designs x (MSE, Coverage, Rank9) x 2 panels ----
p12 <- printed_rows("tab:eight-design"); stopifnot(length(p12) == 18)
n12 <- 0
for (k in seq_along(p12)) {
  nb <- if (k <= 9) "queen" else "rook"
  cells <- split_row(p12[k]); nm <- cells[1]
  des <- names(design_name)[match(sub(" \\(flagged\\)", "", nm), design_name)]
  stopifnot(!is.na(des))
  s <- S[[nb]][S[[nb]]$Design == des, ]
  chk(paste0("S12 ", nb), nm, "MSE", cells[3], fmt3(s$MSE))
  chk(paste0("S12 ", nb), nm, "Coverage", cells[4], fmt3(s$Coverage))
  chk(paste0("S12 ", nb), nm, "Rank9", cells[5], fmt3(s$Rank9))
  n12 <- n12 + 3
}
lines <- c(lines, "", sprintf("Cells checked: S7 %d, S12 %d; mismatches: %d", n7, n12, bad))
writeLines(lines, file.path(out_dir, "performance_nine_check.txt"))
if (bad > 0) stop(bad, " printed cells are not reproduced; see performance_nine_check.txt")
cat(sprintf("All %d S7 and %d S12 cells reproduced exactly.\n", n7, n12))
