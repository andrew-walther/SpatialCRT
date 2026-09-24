# ==============================================================================
# 13_dissertation_results_extract.R
# Pulls the additional real, verified numbers needed for the dissertation
# manuscript's fuller Results section (per-incidence-mode rankings, parameter
# sensitivity by Rho/Gamma, subgroup robustness, best-design-by-config), all
# restricted to the 6 manuscript designs. Reuses the relabeling logic from
# 12_six_design_statistical_comparisons.R. Read-only console summary; no PDF.
#
# 2026-09 revision: every section is reported separately for queen (primary) and
# rook (sensitivity; tau not identified for Checkerboard). Added: per-mode win
# counts, Mean_Treated, flag summary, non-oracle headline, average-rank order per
# tau, and the surface-level robustness tests (50 independent config x surface
# units). Per-tau results share Z and eps across tau (common random numbers), so
# they are not independent confirmations.
# ==============================================================================

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) normalizePath(getwd())
)
setwd(script_dir)
results_dir <- file.path(dirname(script_dir), "results")
source("10_statistical_comparisons.R")

raw_names <- get_design_names()
full_name_map <- setNames(
  sub("^Design [0-9]+: ", "", raw_names),
  paste0("Design ", seq_along(raw_names))
)
retained_ids <- c("Design 1", "Design 2", "Design 4", "Design 5", "Design 6", "Design 8")
ROOK_NOTE <- "Rook: tau not identified for Checkerboard (WZ = 1 - Z); estimates kept but flagged."

relabel6 <- function(df) {
  df <- df[df$Design %in% retained_ids, ]
  df$Design <- full_name_map[df$Design]
  df
}
mle_full_6 <- relabel6(load_latest_results(results_dir = results_dir, estimation_mode = "MLE_tau_sweep"))
non_full_6 <- relabel6(load_latest_results(results_dir = results_dir, estimation_mode = "MLEnonoracle_tau_sweep"))
surf_files <- list.files(file.path(results_dir, "sim_data"),
                         pattern = "^surface_results_MLE_tau_sweep_.*\\.rds$", full.names = TRUE)
stopifnot(length(surf_files) > 0)
surf_6 <- relabel6(readRDS(surf_files[which.max(file.mtime(surf_files))]))
cat("Surface-level results from:", basename(surf_files[which.max(file.mtime(surf_files))]), "\n")

block_cols <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Rho", "Gamma", "Spillover_Type")

out_dir <- file.path(results_dir, "six_design_manuscript")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "dissertation_results_extract.txt"), split = TRUE)

for (nb in c("queen", "rook")) {
  mle_tau1_6 <- mle_full_6[mle_full_6$True_Tau == 1.0 & mle_full_6$Neighbor_Type == nb, ]
  cat(sprintf("\n################ %s ################\n", toupper(nb)))
  if (nb == "rook") cat(ROOK_NOTE, "\n")

  cat("\n=== 1. Overall performance (tau=1.0), all incidence configs averaged ===\n")
  overall <- aggregate(cbind(MSE, Coverage, Power, Mean_Treated) ~ Design, data = mle_tau1_6, FUN = mean)
  print(overall[order(overall$MSE), ], row.names = FALSE)

  cat("\n=== 2. Design ranking by incidence configuration (tau=1.0) ===\n")
  mle_tau1_6$inc_label <- mapply(inc_config_label, mle_tau1_6$Incidence_Mode, mle_tau1_6$Rho_Incidence)
  by_inc <- aggregate(MSE ~ Design + inc_label, data = mle_tau1_6, FUN = mean)
  for (cfg in unique(by_inc$inc_label)) {
    sub <- by_inc[by_inc$inc_label == cfg, ]
    sub <- sub[order(sub$MSE), c("Design", "MSE")]
    cat(sprintf("\n-- %s --\n", cfg))
    print(sub, row.names = FALSE)
  }

  cat("\n=== 3. Parameter sensitivity: mean MSE by Design x Rho (tau=1.0) ===\n")
  by_rho <- aggregate(MSE ~ Design + Rho, data = mle_tau1_6, FUN = mean)
  print(reshape(by_rho, idvar = "Design", timevar = "Rho", direction = "wide"), row.names = FALSE)

  cat("\n=== 4. Parameter sensitivity: mean MSE by Design x Gamma (tau=1.0) ===\n")
  by_gamma <- aggregate(MSE ~ Design + Gamma, data = mle_tau1_6, FUN = mean)
  print(reshape(by_gamma, idvar = "Design", timevar = "Gamma", direction = "wide"), row.names = FALSE)

  cat("\n=== 5. Parameter sensitivity: mean MSE by Design x Spillover_Type (tau=1.0) ===\n")
  by_spill <- aggregate(MSE ~ Design + Spillover_Type, data = mle_tau1_6, FUN = mean)
  print(reshape(by_spill, idvar = "Design", timevar = "Spillover_Type", direction = "wide"), row.names = FALSE)

  cat("\n=== 6. Subgroup robustness: MSE quantiles by design (tau=1.0) ===\n")
  robust <- do.call(rbind, lapply(split(mle_tau1_6, mle_tau1_6$Design), function(d) {
    q <- quantile(d$MSE, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    data.frame(Design = d$Design[1], Best = q[1], Q25 = q[2], Median = q[3], Q75 = q[4], Worst = q[5])
  }))
  print(robust[order(robust$Median), ], row.names = FALSE)

  wide_mse <- reshape(mle_tau1_6[, c(block_cols, "Design", "MSE")],
                      idvar = block_cols, timevar = "Design", direction = "wide")
  design_cols <- grep("^MSE\\.", names(wide_mse), value = TRUE)
  # NA guard: a block with any missing MSE has no well-defined winner; count and drop it
  complete_blk <- stats::complete.cases(wide_mse[, design_cols])
  if (any(!complete_blk)) cat(sprintf("NOTE: %d blocks with missing MSE excluded from win rates\n",
                                      sum(!complete_blk)))
  wide_ok <- wide_mse[complete_blk, ]
  winners <- sub("^MSE\\.", "", design_cols[apply(wide_ok[, design_cols], 1, which.min)])

  cat(sprintf("\n=== 7. Win rate: fraction of %d blocks where design has lowest MSE (tau=1.0) ===\n",
              length(winners)))
  print(sort(table(winners) / length(winners), decreasing = TRUE))

  cat("\n=== 7b. Win counts by incidence configuration (tau=1.0) ===\n")
  wins_cfg <- table(mapply(inc_config_label, wide_ok$Incidence_Mode, wide_ok$Rho_Incidence), winners)
  print(addmargins(wins_cfg, 1))

  cat("\n=== 8. Checkerboard coverage across Rho (tau=1.0) ===\n")
  ck_cov <- aggregate(Coverage ~ Rho, data = mle_tau1_6[mle_tau1_6$Design == "Checkerboard", ], FUN = mean)
  print(ck_cov, row.names = FALSE)

  cat("\n=== 9. Average-rank order per tau (within-block MSE ranks; tau levels share draws) ===\n")
  nb_all <- mle_full_6[mle_full_6$Neighbor_Type == nb, ]
  for (tv in sort(unique(nb_all$True_Tau))) {
    fr <- run_friedman_test(nb_all[nb_all$True_Tau == tv, ])
    cat(sprintf("  tau = %.1f: %s\n", tv,
                paste(sprintf("%s (%.2f)", names(fr$avg_ranks), fr$avg_ranks), collapse = " < ")))
  }

  cat("\n=== 10. Non-oracle (Y ~ Z + X) vs oracle, tau=1.0 ===\n")
  non_tau1 <- non_full_6[non_full_6$True_Tau == 1.0 & non_full_6$Neighbor_Type == nb, ]
  cmp <- merge(aggregate(cbind(MSE, Bias, Coverage) ~ Design, data = mle_tau1_6, FUN = mean),
               aggregate(cbind(MSE, Bias, Coverage) ~ Design, data = non_tau1, FUN = mean),
               by = "Design", suffixes = c(".oracle", ".nonoracle"))
  print(cmp[order(cmp$MSE.oracle), ], row.names = FALSE, digits = 3)

  cat("\n=== 11. Surface-level robustness (paired over 5 configs x 10 surfaces, tau=1.0) ===\n")
  surf_nb <- surf_6[surf_6$Neighbor_Type == nb & surf_6$True_Tau == 1.0, ]
  ord <- overall$Design[order(overall$MSE)]
  for (i in seq_len(length(ord) - 1)) {
    r <- run_surface_paired_test(surf_nb, ord[i], ord[i + 1])
    cat(sprintf("  %-38s vs %-38s n=%d diff %+.4f [%+.4f, %+.4f] Wilcoxon p=%s; first better in %.0f%%\n",
                ord[i], ord[i + 1], r$n_units, r$mean_diff, r$ci[1], r$ci[2],
                format.pval(r$wilcoxon_p, digits = 2), 100 * r$frac_a_better))
  }
}

cat("\n################ FLAG SUMMARY (all tau, both estimators) ################\n")
flag_sum <- function(df, lab) {
  agg <- aggregate(cbind(N_Aliased, N_Warn, Z_WZ_rank_deficient, Fail_Rate) ~ Design + Neighbor_Type,
                   data = df, FUN = sum)
  cat(sprintf("-- %s: scenarios with any aliased fit = %d; total non-aliasing warnings = %d; min N_Valid_Est = %d\n",
              lab, sum(df$N_Aliased > 0), sum(df$N_Warn), min(df$N_Valid_Est)))
  print(agg[agg$N_Aliased > 0 | agg$N_Warn > 0 | agg$Z_WZ_rank_deficient > 0 | agg$Fail_Rate > 0, ],
        row.names = FALSE)
}
flag_sum(mle_full_6, "oracle")
flag_sum(non_full_6, "non-oracle")

sink()
cat("Written to", file.path(out_dir, "dissertation_results_extract.txt"), "\n")
