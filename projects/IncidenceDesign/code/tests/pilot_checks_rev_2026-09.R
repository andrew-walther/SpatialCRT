# ============================================================
# Script: pilot_checks_rev_2026-09.R
# Purpose: B4 pilot checks for the 2026-09 simulation revision (plan B4):
#          integrity, before/after directions, 8-design consolidation, surface-level
#          robustness, non-oracle headline, projected full-run time.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-24
# Dependencies: dplyr, tidyr, PMCMRplus (via 10); sources code/10 (-> 06 -> 01-03)
# ============================================================
#
# Run from code/ after `Rscript 05_run_simulation.R pilot`:
#   Rscript tests/pilot_checks_rev_2026-09.R
# Writes results/pilot_rev_2026-09/pilot_report.txt. Queen is primary; rook is
# reported separately; anything pooled over nb is labeled "POOLED".

suppressMessages(source("10_statistical_comparisons.R"))
library(tidyr)

pilot_dir <- file.path("..", "results", "pilot_rev_2026-09")
latest <- function(pattern) {
  f <- list.files(pilot_dir, pattern = pattern, full.names = TRUE)
  stopifnot(length(f) > 0)
  f[which.max(file.mtime(f))]
}
orc  <- readRDS(latest("^sim_results_MLE_tau_sweep_combined_.*\\.rds$"))
non  <- readRDS(latest("^sim_results_MLEnonoracle_tau_sweep_combined_.*\\.rds$"))
surf <- readRDS(latest("^surface_results_MLE_tau_sweep_.*\\.rds$"))
info <- readRDS(latest("^run_info_.*\\.rds$"))

nm <- setNames(sub("^Design [0-9]+: ", "", get_design_names()), paste("Design", 1:8))
six <- paste("Design", c(1, 2, 4, 5, 6, 8))
label <- function(df) mutate(df, Design = factor(nm[Design], levels = nm))
cfg_lab <- function(mode, rx) ifelse(mode == "iid", "iid", sprintf("%s %.2f", mode, rx))

sink(file.path(pilot_dir, "pilot_report.txt"), split = TRUE)
cat("B4 PILOT REPORT — 2026-09 simulation revision\n")
cat(sprintf("Oracle file: %s\n", basename(latest("^sim_results_MLE_tau_sweep_combined_"))))
cat("Slice: tau = 1; all 5 configs x rook/queen x both regimes; rho in {0, 0.5}; gamma in {0.5, 0.8}; 8 designs\n\n")

# 1. Integrity ----
cat("== 1. Integrity ==\n")
chk <- c(
  "640 scenario rows per estimator"          = nrow(orc) == 640 && nrow(non) == 640,
  "N_Valid_Est = 250 everywhere"             = all(orc$N_Valid_Est == 250) && all(non$N_Valid_Est == 250),
  "No non-aliasing warnings (N_Warn = 0)"    = all(orc$N_Warn == 0) && all(non$N_Warn == 0),
  "Aliasing only Checkerboard x rook (oracle, all 250)" =
    all((orc$N_Aliased > 0) == (orc$Design == "Design 1" & orc$Neighbor_Type == "rook")) &&
    all(orc$N_Aliased[orc$N_Aliased > 0] == 250) && all(non$N_Aliased == 0),
  "Z_WZ_rank_deficient only Checkerboard x rook (both files)" =
    all(orc$Z_WZ_rank_deficient == (orc$Design == "Design 1" & orc$Neighbor_Type == "rook")) &&
    all(non$Z_WZ_rank_deficient == (non$Design == "Design 1" & non$Neighbor_Type == "rook")),
  "HIF and Balanced Quartiles Mean_Treated = 50" = all(orc$Mean_Treated[orc$Design %in% c("Design 2", "Design 6")] == 50)
)
for (n in names(chk)) cat(sprintf("  %-58s %s\n", n, if (chk[[n]]) "PASS" else "FAIL"))

# Poisson distinct X values per surface (regenerate the surfaces from their keys)
sim_profile <- "full"; sim_define_only <- TRUE; sim_script_dir <- normalizePath(".")
suppressMessages(source("05_run_simulation.R"))
for (ci in which(vapply(inc_configs, `[[`, "", "mode") == "poisson")) {
  Xp <- generate_surfaces(inc_configs[[ci]])
  nd <- apply(Xp, 2, function(x) length(unique(x)))
  cat(sprintf("  Poisson rhoX %.2f: distinct X values per surface median %g (range %d-%d) [target >= ~55] %s\n",
              inc_configs[[ci]]$rho_x, median(nd), min(nd), max(nd),
              if (median(nd) >= 55) "PASS" else "FAIL"))
}
mt <- label(orc) %>% filter(Estimator == "oracle") %>% group_by(Design, Neighbor_Type) %>%
  summarise(Mean_Treated = mean(Mean_Treated), .groups = "drop") %>%
  pivot_wider(names_from = Neighbor_Type, values_from = Mean_Treated)
cat("\n  Mean treated per design:\n"); print(as.data.frame(mt), row.names = FALSE)

# 2. Before/after (same slice of the April 2026 tau-sweep results) ----
cat("\n== 2. Before (April 2026) vs after, 6 manuscript designs, MSE / coverage ==\n")
old <- readRDS(file.path("..", "results", "archive", "pre_revision_20260924", "sim_data",
                         "sim_results_MLE_tau_sweep_combined_20260408_191916.rds"))
old <- old %>% filter(True_Tau == 1, Rho %in% c(0, 0.5), Gamma %in% c(0.5, 0.8))
stopifnot(nrow(old) == 640)
ba <- bind_rows(
  old %>% mutate(When = "before"),
  orc %>% mutate(When = "after")) %>%
  filter(Design %in% six) %>% label() %>%
  mutate(Config = cfg_lab(Incidence_Mode, Rho_Incidence)) %>%
  group_by(Neighbor_Type, Config, Design, When) %>%
  summarise(MSE = mean(MSE), Cov = mean(Coverage), .groups = "drop") %>%
  pivot_wider(names_from = When, values_from = c(MSE, Cov)) %>%
  group_by(Neighbor_Type, Config) %>%
  mutate(Rank_before = rank(MSE_before), Rank_after = rank(MSE_after)) %>%
  ungroup() %>%
  arrange(desc(Neighbor_Type), Config, Rank_after)
print(as.data.frame(ba %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\n  Rank of each design (queen, per config), before -> after:\n")
print(as.data.frame(ba %>% filter(Neighbor_Type == "queen") %>%
  transmute(Config, Design, change = paste(Rank_before, "->", Rank_after)) %>%
  pivot_wider(names_from = Config, values_from = change)), row.names = FALSE)

# 3. 8-design consolidation check (14(a) logic), by neighbor type ----
cat("\n== 3. 8-design consolidation (Friedman / Nemenyi / Wilcoxon-Holm) ==\n")
cons_pairs <- list(c("Saturation Quadrants", "Incidence-Guided Saturation Quadrants"),
                   c("Balanced Halves", "Balanced Quartiles"))
for (nbt in c("queen", "rook", "POOLED")) {
  d8 <- label(orc) %>% mutate(Design = as.character(Design))
  if (nbt != "POOLED") d8 <- d8 %>% filter(Neighbor_Type == nbt)
  fr <- run_friedman_test(d8); nem <- run_nemenyi_posthoc(d8); wil <- run_pairwise_wilcoxon(d8)
  cat(sprintf("\n-- %s: %d blocks; Friedman chi2 = %.1f, p = %s\n", nbt, fr$n_blocks,
              fr$statistic, format.pval(fr$p_value, digits = 3)))
  desc <- d8 %>% group_by(Design) %>%
    summarise(MSE = mean(MSE), Coverage = mean(Coverage), .groups = "drop") %>% arrange(MSE)
  desc$AvgRank <- round(fr$avg_ranks[desc$Design], 2)
  print(as.data.frame(desc %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
  for (pr in cons_pairs) {
    cat(sprintf("  %s vs %s: MSE %.4f vs %.4f | Nemenyi p = %.4f | Wilcoxon-Holm p = %s\n",
                pr[1], pr[2], desc$MSE[desc$Design == pr[1]], desc$MSE[desc$Design == pr[2]],
                nem$p_matrix[pr[1], pr[2]], format.pval(wil$p_matrix[pr[1], pr[2]], digits = 3)))
  }
}

# 4. Surface-level robustness (50 independent config x surface units), queen / rook ----
cat("\n== 4. Surface-level paired tests (units = 5 configs x 10 surfaces; MSE averaged over each unit's blocks) ==\n")
key_pairs <- list(
  c("Incidence-Guided Saturation Quadrants", "Balanced Quartiles"),
  c("Incidence-Guided Saturation Quadrants", "Saturation Quadrants"),
  c("Balanced Quartiles", "Balanced Halves"),
  c("Balanced Quartiles", "High Incidence Focus"),
  c("Isolation Buffer", "2x2 Blocking"),
  c("2x2 Blocking", "Checkerboard"))
sl <- label(surf) %>% mutate(Design = as.character(Design))
for (nbt in c("queen", "rook")) {
  cat(sprintf("-- %s\n", nbt))
  for (pr in key_pairs) {
    r <- run_surface_paired_test(sl %>% filter(Neighbor_Type == nbt), pr[1], pr[2])
    cat(sprintf("  %-38s - %-38s n=%d  diff %+.4f  95%% CI [%+.4f, %+.4f]  Wilcoxon p=%s  a better in %2.0f%%\n",
                pr[1], pr[2], r$n_units, r$mean_diff, r$ci[1], r$ci[2],
                format.pval(r$wilcoxon_p, digits = 2), 100 * r$frac_a_better))
  }
}

# 5. Non-oracle headline ----
cat("\n== 5. Non-oracle (Y ~ Z + X) vs oracle, 6 designs, queen and rook ==\n")
nv <- bind_rows(orc, non) %>% filter(Design %in% six) %>% label() %>%
  group_by(Neighbor_Type, Design, Estimator) %>%
  summarise(MSE = mean(MSE), Bias = mean(Bias), Cov = mean(Coverage), .groups = "drop") %>%
  pivot_wider(names_from = Estimator, values_from = c(MSE, Bias, Cov)) %>%
  arrange(desc(Neighbor_Type), MSE_oracle)
print(as.data.frame(nv %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# 6. Runtime ----
cat("\n== 6. Runtime ==\n")
per_fit_ms <- 1000 * sum(info$unit_elapsed_sec) / (2 * 640 * 250)
full_fits <- 2 * 12800 * 250
cat(sprintf("  Pilot: %.1f min wall (10 workers), %.1f min unit compute, %.2f ms per fit incl. DGP\n",
            info$wall_min, sum(info$unit_elapsed_sec) / 60, per_fit_ms))
cat(sprintf("  Projected full run (%s fits): %.0f CPU-min -> ~%.0f min wall on 10 workers (40 units, 4 rounds)\n",
            format(full_fits, big.mark = ","), full_fits * per_fit_ms / 60000,
            full_fits * per_fit_ms / 60000 / 10))
sink()
