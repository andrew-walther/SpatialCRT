# ==============================================================================
# 13_dissertation_results_extract.R
# Pulls the additional real, verified numbers needed for the dissertation
# manuscript's fuller Results section (per-incidence-mode rankings, parameter
# sensitivity by Rho/Gamma, subgroup robustness, best-design-by-config), all
# restricted to the 6 manuscript designs. Reuses the relabeling logic from
# 12_six_design_statistical_comparisons.R. Read-only console summary; no PDF.
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

mle_full <- load_latest_results(results_dir = results_dir, estimation_mode = "MLE_tau_sweep")
mle_full_6 <- mle_full[mle_full$Design %in% retained_ids, ]
mle_full_6$Design <- full_name_map[mle_full_6$Design]
mle_tau1_6 <- mle_full_6[mle_full_6$True_Tau == 1.0, ]

out_dir <- file.path(results_dir, "six_design_manuscript")
sink(file.path(out_dir, "dissertation_results_extract.txt"), split = TRUE)

cat("=== 1. Overall performance (tau=1.0), all configs pooled ===\n")
overall <- aggregate(cbind(MSE, Coverage, Power) ~ Design, data = mle_tau1_6, FUN = mean)
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

cat("\n=== 7. Win-rate: fraction of 320 scenarios where design has lowest MSE (tau=1.0) ===\n")
wide_mse <- reshape(mle_tau1_6[, c("Incidence_Mode","Rho_Incidence","Neighbor_Type","Rho","Gamma","Spillover_Type","Design","MSE")],
                     idvar = c("Incidence_Mode","Rho_Incidence","Neighbor_Type","Rho","Gamma","Spillover_Type"),
                     timevar = "Design", direction = "wide")
design_cols <- grep("^MSE\\.", names(wide_mse), value = TRUE)
winners <- design_cols[apply(wide_mse[, design_cols], 1, which.min)]
winners <- sub("^MSE\\.", "", winners)
print(sort(table(winners) / length(winners), decreasing = TRUE))

cat("\n=== 8. Coverage by design, Checkerboard detail across Rho (tau=1.0) ===\n")
ck_cov <- aggregate(Coverage ~ Rho, data = mle_tau1_6[mle_tau1_6$Design == "Checkerboard", ], FUN = mean)
print(ck_cov, row.names = FALSE)

sink()
cat("Written to", file.path(out_dir, "dissertation_results_extract.txt"), "\n")
