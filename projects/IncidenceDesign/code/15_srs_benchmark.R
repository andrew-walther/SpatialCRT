# ============================================================
# Script: 15_srs_benchmark.R
# Purpose: Tabulate the Simple Random Sampling (SRS) benchmark (Design 9) against
#          the six proposed designs and the two consolidated ones: relative
#          efficiency, share of blocks beating SRS, and paired surface-level tests.
# Author: Andrew Walther
# Created: 2026-09-26
# Dependencies: dplyr, tidyr (via 10_statistical_comparisons.R)
# ============================================================
#
# SRS is a benchmark, not a proposed design (author, 2026-09-26). It never enters
# the rank-based six-design statistics of 12-14 (Friedman, Nemenyi, average rank,
# win rates, CD diagrams), which therefore do not change. Its role is to show each
# design's value added over naive complete randomization:
#
#   relative efficiency  RE_d = MSE_d / MSE_SRS        (< 1 means d beats SRS)
#   block share          P_d  = share of blocks with MSE_d < MSE_SRS
#   paired test          D_u  = m_{d,u} - m_{SRS,u} over the 50 independent
#                        (configuration x surface) units, t49 CI + Wilcoxon,
#                        Holm-adjusted across the 8 designs within each slice
#
# Queen is primary, rook the sensitivity (tau not identified for Checkerboard under
# rook). Pooled numbers are labeled as pooled.
#
# STOP RULE (author): if SRS has lower pooled MSE than Incidence-Guided Saturation
# Quadrants (queen, tau = 1, or pooled over tau), the recommendation needs the
# author's review; the script prints STOP-RULE TRIGGERED.
#
# Usage (from code/):  Rscript 15_srs_benchmark.R
# OUTPUTS (results/srs_benchmark/):
#   srs_benchmark_summary.txt   console summary (oracle primary, non-oracle appended)
#   srs_benchmark.rds           list of the tables below, per estimator and nb

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
setwd(script_dir)
results_dir <- file.path(dirname(script_dir), "results")
source("10_statistical_comparisons.R")

# Labels ----
raw_names <- get_design_names()
full_name_map <- setNames(sub("^Design [0-9]+: ", "", raw_names),
                          paste0("Design ", seq_along(raw_names)))
SRS <- "Simple Random Sampling"
stopifnot(full_name_map[["Design 9"]] == SRS)
retained <- unname(full_name_map[c("Design 1", "Design 2", "Design 4", "Design 5",
                                   "Design 6", "Design 8")])
consolidated <- unname(full_name_map[c("Design 3", "Design 7")])
IGSQ <- "Incidence-Guided Saturation Quadrants"

relabel <- function(df) { df$Design <- full_name_map[df$Design]; df }
latest_file <- function(pattern) {
  f <- list.files(file.path(results_dir, "sim_data"), pattern = pattern, full.names = TRUE)
  stopifnot(length(f) > 0)
  f[which.max(file.mtime(f))]
}

scen <- list(
  oracle    = relabel(load_latest_results(results_dir = results_dir, estimation_mode = "MLE_tau_sweep")),
  nonoracle = relabel(load_latest_results(results_dir = results_dir, estimation_mode = "MLEnonoracle_tau_sweep")))
surf <- list(
  oracle    = relabel(readRDS(latest_file("^surface_results_MLE_tau_sweep_.*\\.rds$"))),
  nonoracle = relabel(readRDS(latest_file("^surface_results_MLEnonoracle_tau_sweep_.*\\.rds$"))))
for (e in names(scen)) stopifnot(SRS %in% scen[[e]]$Design, SRS %in% surf[[e]]$Design)

block_cols <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Rho", "Gamma",
                "Spillover_Type", "True_Tau")

# Helpers ----

#' Mean MSE per design with its Monte Carlo SE from the 50 unit means
#'
#' SE(MSE) = sd(m_1..m_50) / sqrt(50), where m_u is the design's mean MSE over the
#' blocks of (configuration, surface) unit u -- the chapter's Table 2 convention.
#'
#' @param sc Scenario rows for one slice (one nb, one tau)
#' @param su Surface rows for the same slice
#' @return Data frame: Design, MSE, SE_MSE, Bias, SD, Coverage, Power, Mean_Treated, RE
perf_table <- function(sc, su) {
  perf <- sc %>% group_by(Design) %>%
    summarise(MSE = mean(MSE), Bias = mean(Bias), SD = mean(SD), Coverage = mean(Coverage),
              Power = mean(Power), Mean_Treated = mean(Mean_Treated), .groups = "drop")
  se <- su %>% group_by(Design, Incidence_Mode, Rho_Incidence, Surface) %>%
    summarise(m = mean(MSE), .groups = "drop") %>%
    group_by(Design) %>% summarise(SE_MSE = sd(m) / sqrt(n()), n_units = n(), .groups = "drop")
  stopifnot(all(se$n_units == 50))
  out <- merge(perf, se[, c("Design", "SE_MSE")], by = "Design")
  out$RE <- out$MSE / out$MSE[out$Design == SRS]
  out[order(out$MSE), c("Design", "MSE", "SE_MSE", "RE", "Bias", "SD", "Coverage", "Power",
                        "Mean_Treated")]
}

#' Share of blocks in which each design's MSE is below SRS's
#'
#' @param sc Scenario rows for one slice
#' @return Named numeric vector (designs other than SRS)
block_share <- function(sc) {
  w <- tidyr::pivot_wider(sc[, c(block_cols, "Design", "MSE")], names_from = Design,
                          values_from = MSE)
  ds <- setdiff(unique(sc$Design), SRS)
  sapply(ds, function(d) mean(w[[d]] < w[[SRS]]))
}

#' Paired surface-level tests of every design against SRS, Holm-adjusted
#'
#' @param su Surface rows for one slice
#' @return Data frame: Design, n_units, mean_diff, ci_lo, ci_hi, wilcoxon_p,
#'   p_holm, frac_better
vs_srs_tests <- function(su) {
  ds <- setdiff(unique(su$Design), SRS)
  rows <- lapply(ds, function(d) {
    r <- run_surface_paired_test(su, d, SRS)
    data.frame(Design = d, n_units = r$n_units, mean_diff = r$mean_diff, ci_lo = r$ci[1],
               ci_hi = r$ci[2], wilcoxon_p = r$wilcoxon_p, frac_better = r$frac_a_better)
  })
  out <- do.call(rbind, rows)
  out$p_holm <- p.adjust(out$wilcoxon_p, method = "holm")
  out[order(out$mean_diff), ]
}

#' SRS MSE and each design's relative efficiency within levels of one parameter
#'
#' @param sc Scenario rows for one slice
#' @param by Column name (e.g. "Rho", "Gamma", "Spillover_Type", or "Config")
#' @return Wide data frame: level x (MSE_SRS, RE per design)
re_by <- function(sc, by) {
  m <- sc %>% group_by(.data[[by]], Design) %>% summarise(MSE = mean(MSE), .groups = "drop")
  srs <- m[m$Design == SRS, c(by, "MSE")]; names(srs)[2] <- "MSE_SRS"
  m <- merge(m[m$Design != SRS, ], srs, by = by)
  m$RE <- m$MSE / m$MSE_SRS
  w <- tidyr::pivot_wider(m[, c(by, "MSE_SRS", "Design", "RE")], names_from = Design,
                          values_from = RE)
  as.data.frame(w)
}

# Main ----
out_dir <- file.path(results_dir, "srs_benchmark")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
res <- list()
stop_rule <- character(0)

sink(file.path(out_dir, "srs_benchmark_summary.txt"), split = TRUE)
cat("======================================================================\n")
cat("SIMPLE RANDOM SAMPLING BENCHMARK (Design 9; complete randomization, 50 of 100)\n")
cat("RE = MSE_design / MSE_SRS (< 1: design beats SRS). Paired tests over 50\n")
cat("configuration x surface units, Holm-adjusted within slice. tau levels share\n")
cat("draws (common random numbers).\n")
cat("======================================================================\n")

for (est in c("oracle", "nonoracle")) {
  cat(sprintf("\n\n==================== ESTIMATOR: %s ====================\n",
              if (est == "oracle") "oracle (Y ~ Z + Spill + X), PRIMARY" else
                "non-oracle (Y ~ Z + X), sensitivity"))
  for (nb in c("queen", "rook")) {
    sc_nb <- scen[[est]][scen[[est]]$Neighbor_Type == nb, ]
    su_nb <- surf[[est]][surf[[est]]$Neighbor_Type == nb, ]
    sc1 <- sc_nb[sc_nb$True_Tau == 1, ]; su1 <- su_nb[su_nb$True_Tau == 1, ]
    sc1$Config <- mapply(inc_config_label, sc1$Incidence_Mode, sc1$Rho_Incidence)
    six <- c(retained, SRS)
    cat(sprintf("\n######## %s %s ########\n", toupper(nb),
                if (nb == "queen") "(primary)" else "(sensitivity; Checkerboard tau not identified)"))

    r <- list()
    r$perf6 <- perf_table(sc1[sc1$Design %in% six, ], su1[su1$Design %in% six, ])
    r$perf8 <- perf_table(sc1, su1)
    cat("\n-- 1. Performance at tau = 1, six proposed designs + SRS benchmark --\n")
    print(r$perf6, row.names = FALSE, digits = 3)
    cat("\n-- 1b. Consolidated designs (Appendix A1) vs SRS, tau = 1 --\n")
    print(r$perf8[r$perf8$Design %in% c(consolidated, SRS), ], row.names = FALSE, digits = 3)

    r$block_share <- block_share(sc1)
    cat(sprintf("\n-- 2. Share of the %d blocks with MSE below SRS (tau = 1) --\n",
                nrow(sc1) / length(unique(sc1$Design))))
    print(round(sort(r$block_share, decreasing = TRUE), 3))

    r$tests <- vs_srs_tests(su1)
    cat("\n-- 3. Paired surface-level tests vs SRS (tau = 1; diff = design - SRS, negative = better) --\n")
    print(r$tests, row.names = FALSE, digits = 3)

    r$by_config <- re_by(sc1, "Config")
    r$by_rho    <- re_by(sc1, "Rho")
    r$by_gamma  <- re_by(sc1, "Gamma")
    r$by_regime <- re_by(sc1, "Spillover_Type")
    r$by_tau    <- re_by(sc_nb, "True_Tau")
    for (nm in c("by_config", "by_rho", "by_gamma", "by_regime", "by_tau")) {
      cat(sprintf("\n-- 4. MSE_SRS and RE %s%s --\n", sub("by_", "by ", nm),
                  if (nm == "by_tau") " (all tau)" else " (tau = 1)"))
      print(r[[nm]][, c(1, 2, match(intersect(c(retained, consolidated), names(r[[nm]])), names(r[[nm]])))],
            row.names = FALSE, digits = 3)
    }

    blk <- sc1[sc1$Design == SRS, "MSE"]
    r$srs_block_quantiles <- quantile(blk, c(0, 0.25, 0.5, 0.75, 1))
    cat("\n-- 5. SRS block-level MSE quantiles (tau = 1; min, Q1, median, Q3, max) --\n")
    print(round(r$srs_block_quantiles, 4))

    # Stop rule (oracle only): pooled MSE of IGSQ vs SRS
    if (est == "oracle") {
      m1  <- tapply(sc1$MSE, sc1$Design, mean)
      mall <- tapply(sc_nb$MSE, sc_nb$Design, mean)
      cat(sprintf("\n-- STOP-RULE CHECK: pooled MSE, %s --\n", nb))
      cat(sprintf("   tau = 1:      IGSQ %.4f vs SRS %.4f (RE %.3f)\n", m1[IGSQ], m1[SRS], m1[IGSQ] / m1[SRS]))
      cat(sprintf("   pooled tau:   IGSQ %.4f vs SRS %.4f (RE %.3f)\n", mall[IGSQ], mall[SRS], mall[IGSQ] / mall[SRS]))
      beaten6 <- names(m1)[names(m1) %in% retained & m1 > m1[SRS]]
      cat("   proposed designs with higher pooled MSE than SRS at tau = 1:",
          if (length(beaten6)) paste(beaten6, collapse = ", ") else "none", "\n")
      if (nb == "queen" && (m1[SRS] < m1[IGSQ] || mall[SRS] < mall[IGSQ]))
        stop_rule <- c(stop_rule, sprintf("%s: SRS pooled MSE below IGSQ", nb))
    }
    res[[est]][[nb]] <- r
  }
}
cat("\n\n", if (length(stop_rule)) paste("STOP-RULE TRIGGERED:", paste(stop_rule, collapse = "; "))
    else "Stop rule not triggered: IGSQ has lower pooled MSE than SRS (queen, tau = 1 and pooled over tau).",
    "\n", sep = "")
sink()
saveRDS(c(res, list(stop_rule = stop_rule)), file.path(out_dir, "srs_benchmark.rds"))
cat("Written to", out_dir, "\n")
