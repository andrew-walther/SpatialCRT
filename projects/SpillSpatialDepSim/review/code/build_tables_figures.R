# ============================================================
# Script: review/code/build_tables_figures.R
# Purpose: Rebuild the beta-MSE summary tables and box-plot figures for the
#          2x4, 3x3 and 3x4 grids, both as published and corrected.
# Author: Andrew Walther (review written with Claude Code)
# Created: 2026-10-22
# Dependencies: dplyr, ggplot2
# ============================================================
# Two versions are built from the same code:
#   "published": original CSVs (../results/beta_mse) and BSS = combos c(28, 43) for every
#                grid, as in code/SimEstimateAnalysisAll.Rmd. Must reproduce the paper.
#   "corrected": BSS = each grid's own block-stratified combos, and for 3x3/3x4 the
#                TrtNoSpill CSVs rerun with the case-1 loop bound fixed
#                (review/results/beta_mse_TrtNoSpill_fixed). 2x4 inputs are unchanged.
# Reads ../results (read-only); writes only to review/tables and review/figures.
# Run from review/:  Rscript code/build_tables_figures.R

suppressMessages({library(dplyr); library(ggplot2)})

# Inputs ----
orig_dir  <- "../results/beta_mse"
fixed_dir <- "results/beta_mse_TrtNoSpill_fixed"
grids     <- c("2x4", "3x3", "3x4")
bss_idx   <- list(published = list(`2x4` = c(28, 43), `3x3` = c(28, 43), `3x4` = c(28, 43)),
                  corrected = list(`2x4` = c(28, 43), `3x3` = c(28, 30, 35, 51, 77, 106),
                                   `3x4` = c(313, 612)))
psis <- c(0.5, 0.6, 0.7, 0.8); rhos <- c(0.01, 0)

#' Path to the per-allocation beta-MSE CSV for one scenario.
#'
#' @param version "published" or "corrected".
#' @param spill "TrtNoSpill" or "TrtSpill".
#' @param grid "2x4", "3x3" or "3x4".
#' @param psi,rho Scenario values.
#' @return File path. Corrected 3x3/3x4 TrtNoSpill files come from the rerun.
csv_path <- function(version, spill, grid, psi, rho) {
  f <- sprintf("beta_mse_results_%s_%s_A02_B1_P%02d_R%s.csv", spill, grid,
               round(psi * 10), ifelse(rho == 0, "00", "001"))
  if (version == "corrected" && spill == "TrtNoSpill" && grid != "2x4") file.path(fixed_dir, f)
  else file.path(orig_dir, f)
}

#' Long data frame of per-allocation results, labelled SRS ("Random", all combos)
#' and BSS ("Block", the given combos), as in SimEstimateAnalysisAll.Rmd.
#'
#' @param version "published" or "corrected".
#' @param grid Grid label.
#' @return Data frame with MSE, bias, Psi, Rho, TrtSpill, Sampling.
load_grid <- function(version, grid) {
  out <- list()
  for (spill in c("TrtNoSpill", "TrtSpill")) for (psi in psis) for (rho in rhos) {
    d <- read.csv(csv_path(version, spill, grid, psi, rho))
    d <- data.frame(MSE = d$MSE, bias = d$bias, Psi = psi, Rho = rho,
                    TrtSpill = ifelse(spill == "TrtSpill", "Intervention & Control", "Control Only"))
    out[[length(out) + 1]] <- mutate(d, Sampling = "Random")
    out[[length(out) + 1]] <- mutate(d[bss_idx[[version]][[grid]], ], Sampling = "Block")
  }
  bind_rows(out)
}

#' Summary table: mean MSE, SD of per-allocation MSE, mean bias, by design.
#'
#' @param df Output of load_grid().
#' @return One row per spillover type x rho x psi, SRS and BSS side by side (unscaled).
summarise_grid <- function(df) {
  df %>% group_by(TrtSpill, Rho, Psi, Sampling) %>%
    summarise(mean_MSE = mean(MSE), sd_MSE = sd(MSE), mean_bias = mean(bias), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sampling, values_from = c(mean_MSE, sd_MSE, mean_bias)) %>%
    arrange(TrtSpill, desc(Rho), Psi) %>%
    select(TrtSpill, Rho, Psi,
           SRS_MSE = mean_MSE_Random, SRS_SD = sd_MSE_Random, SRS_bias = mean_bias_Random,
           BSS_MSE = mean_MSE_Block,  BSS_SD = sd_MSE_Block,  BSS_bias = mean_bias_Block)
}

#' Box-plot figure in the published layout (rows = spillover type, cols = rho).
#'
#' @param df Output of load_grid().
#' @return ggplot object (code copied from SimEstimateAnalysisAll.Rmd, 3x4 2x2 chunk).
plot_grid <- function(df) {
  ggplot(df, aes(x = Sampling, y = MSE, fill = as.factor(Psi))) +
    geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.75)) +
    theme_minimal(base_size = 14) +
    labs(x = "Sampling Group", y = "Mean Squared Error (MSE)", fill = expression(Psi)) +
    facet_grid(TrtSpill ~ Rho, scales = "free_y",
               labeller = label_bquote(rows = .(TrtSpill), cols = rho == .(Rho))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
          strip.text.y = element_text(size = 8, angle = 90),
          axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 10))
}

# Build ----
for (version in c("published", "corrected")) for (grid in grids) {
  df  <- load_grid(version, grid)
  tab <- summarise_grid(df)
  stopifnot(nrow(tab) == 16, !anyNA(tab))
  write.csv(tab, sprintf("tables/beta_mse_table_%s_%s.csv", grid, version), row.names = FALSE)
  ggsave(sprintf("figures/MSE_Estimates_%s_%s.png", grid, version), plot_grid(df),
         width = 7, height = 4.5, dpi = 150, bg = "white")
}
cat("Wrote", length(list.files("tables")), "tables and", length(list.files("figures")), "figures\n")
