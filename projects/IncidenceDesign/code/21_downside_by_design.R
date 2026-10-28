# ============================================================
# Script: 21_downside_by_design.R
# Purpose: Compare every design with SRS on mean allocation MSE and on
#          allocation downside (worst-10% mean), by spillover regime, for the
#          grid allocation-risk study and the refined NC application corners.
# Author: Andrew Walther
# Created: 2026-10-11
# Dependencies: base R; ggplot2 for the distribution figure
# ============================================================
# Run from the project root: Rscript code/21_downside_by_design.R [profile]
# profile is the grid run of code/16: surfaces10 (default; all ten surfaces, the
# reported study) or pilot (first two surfaces). Reads completed aggregate outputs
# only; it neither simulates nor changes checkpoints. Outputs go to
# results/allocation_risk/summary/surfaces10/ (pilot: results/allocation_risk/summary/):
# downside_by_design.csv, allocation_distribution_by_design.csv and
# allocation_distribution_by_design.pdf (CTJ Figure 3). Only the surfaces10 figure
# is copied to the CTJ and chapter figure folders.
#
# For allocation a, m(a) = E_eps[(tau_hat - tau)^2 | a, X] is its conditional
# MSE. Each source row already summarizes the sampled allocations of one
# design: their mean m(a) and the mean of the worst 10% of m(a).
#   Grid: queen, tau = 1, oracle fit, averaged over 20 x (surfaces per config)
#     fixed surface/corner units (200 for surfaces10, 40 for the pilot).
#   NC: refined queen corners rho in {0, 0.5}, gamma = 0.8, education model,
#     mean-rank regions; averaged with equal weight over 2018-2021 and both
#     corners (8 rows per design and regime).
# Ratios are (design average) / (SRS average) within the same source and regime,
# so a value below 1 means lower error than SRS.

# Grid profile and output folder ----
args <- commandArgs(TRUE); profile <- if (length(args)) args[1] else "surfaces10"
n_surfaces <- c(pilot = 2L, surfaces10 = 10L)[[profile]]
out_dir <- if (profile == "pilot") "results/allocation_risk/summary" else
  file.path("results/allocation_risk/summary", profile)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
grid_source <- paste0("grid_", profile)

# Read sources ----
grid <- read.csv(file.path("results/allocation_risk/summary", profile, "pooled_fixed_units.csv"),
                 stringsAsFactors = FALSE)
nc <- read.csv(paste0("application/results/real_sud_rev_20261002/exhibits/",
                      "refined_tail_comparison.csv"), stringsAsFactors = FALSE)

grid <- grid[grid$Neighbor_Type == "queen" & grid$Estimator == "oracle" &
               grid$True_Tau == 1, ]
grid$Design_ID <- as.integer(sub("Design ", "", grid$Design))
nc <- nc[nc$Neighbor == "queen" & nc$Model == "education" &
           nc$Summary == "mean_rank" & nc$Gamma == 0.8, ]

# Fail loud if the sources are not the expected shape ----
stopifnot(nrow(grid) == 18, setequal(grid$Design_ID, 1:9),
          all(grid$N_Fixed_Units == 20 * n_surfaces))
nc_counts <- table(nc$Design_ID, nc$Regime)
stopifnot(nrow(nc_counts) == 9, all(nc_counts == 8))

# Summarize: average mean MSE and worst-10% mean, then ratio to SRS (ID 9) ----
summarize_source <- function(ids, regimes, mean_mse, worst10, source) {
  out <- aggregate(cbind(Mean_MSE = mean_mse, Worst10 = worst10) ~
                     Design_ID + Regime,
                   data = data.frame(Design_ID = ids, Regime = regimes,
                                     mean_mse = mean_mse, worst10 = worst10),
                   FUN = mean)
  srs <- out[out$Design_ID == 9, c("Regime", "Mean_MSE", "Worst10")]
  names(srs) <- c("Regime", "SRS_Mean_MSE", "SRS_Worst10")
  out <- merge(out, srs, by = "Regime")
  out$Mean_Ratio_SRS <- out$Mean_MSE / out$SRS_Mean_MSE
  out$Worst10_Ratio_SRS <- out$Worst10 / out$SRS_Worst10
  out$Source <- source
  out
}

result <- rbind(
  summarize_source(grid$Design_ID, grid$Spillover_Type,
                   grid$Mean_Allocation_MSE, grid$Average_Conditional_Worst10,
                   grid_source),
  summarize_source(nc$Design_ID, nc$Regime, nc$Mean_MSE,
                   nc$Worst10_Mean_Estimated, "nc_refined")
)
result <- result[order(result$Source, result$Regime, result$Design_ID),
                 c("Source", "Regime", "Design_ID", "Mean_MSE", "Worst10",
                   "Mean_Ratio_SRS", "Worst10_Ratio_SRS")]

# Write ----
out_file <- file.path(out_dir, "downside_by_design.csv")
write.csv(result, out_file, row.names = FALSE)
cat("Wrote", nrow(result), "rows to", out_file, "\n")

# Allocation-to-allocation distribution (CTJ Figure 3; added 2026-10-12) ----
# Per-allocation conditional MSE estimates m(a), each divided by the
# frequency-weighted SRS mean in the same unit. Percentiles are computed within
# each unit and then averaged over units, so the bars show allocation spread
# only, not differences between units. A design with one allocation per unit has
# q05 = q95. m(a) is estimated from 100 (grid) or >= 400 (NC) outcomes, so the
# spreads include some outcome Monte Carlo noise and slightly overstate the
# true allocation spread.

#' Weighted quantile with frequency weights (how often an allocation was drawn)
#' @param x Numeric values. @param w Nonnegative weights. @param p Probabilities.
#' @return Values at the requested weighted quantiles.
weighted_quantile <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w) / sum(w)
  sapply(p, function(q) x[which(cw >= q)[1]])
}

#' Divide each allocation's m(a) by the weighted SRS (ID 9) mean in its unit
#' @param d Per-allocation rows with MSE, Frequency, Design_ID.
#' @param unit_cols Columns that identify one fixed setting.
#' @return d with Unit and Ratio columns.
relative_to_srs <- function(d, unit_cols) {
  d$Unit <- do.call(paste, c(d[unit_cols], sep = "/"))
  srs <- d[d$Design_ID == 9, ]
  ref <- tapply(srs$MSE * srs$Frequency, srs$Unit, sum) / tapply(srs$Frequency, srs$Unit, sum)
  stopifnot(all(d$Unit %in% names(ref)))
  d$Ratio <- d$MSE / ref[d$Unit]
  d
}

ga <- read.csv(file.path("results/allocation_risk", profile, "allocation_metrics.csv"),
               stringsAsFactors = FALSE)
ga <- ga[ga$Neighbor_Type == "queen" & ga$True_Tau == 1 & ga$Estimator == "oracle", ]
ga$Design_ID <- as.integer(sub("Design ", "", ga$Design))
ga$Regime <- ga$Spillover_Type
ga <- relative_to_srs(ga, c("Incidence_Mode", "Rho_Incidence", "Surface", "Rho", "Gamma", "Regime"))
ga$Source <- grid_source
na <- read.csv(paste0("application/results/real_sud_rev_20261002/tail_confirmation/",
                      "allocation_metrics.csv"), stringsAsFactors = FALSE)
na <- na[na$Model == "education" & na$Neighbor == "queen" & na$Summary == "mean_rank" &
           na$Gamma == 0.8, ]
na <- relative_to_srs(na, c("Year", "Rho", "Regime"))
na$Source <- "nc_refined"
keep <- c("Source", "Regime", "Design_ID", "Unit", "Ratio", "Frequency")
alloc <- rbind(ga[keep], na[keep])
stopifnot(all(is.finite(alloc$Ratio)), length(unique(ga$Unit)) == 40 * n_surfaces,
          length(unique(na$Unit)) == 16)
per_unit <- do.call(rbind, lapply(split(alloc, list(alloc$Source, alloc$Regime, alloc$Design_ID,
                                                    alloc$Unit), drop = TRUE), function(x) {
  q <- weighted_quantile(x$Ratio, x$Frequency, c(0.05, 0.25, 0.5, 0.75, 0.95))
  data.frame(Source = x$Source[1], Regime = x$Regime[1], Design_ID = x$Design_ID[1],
             Q05 = q[1], Q25 = q[2], Q50 = q[3], Q75 = q[4], Q95 = q[5],
             Mean = sum(x$Ratio * x$Frequency) / sum(x$Frequency))
}))
dist <- aggregate(cbind(Q05, Q25, Q50, Q75, Q95, Mean) ~ Source + Regime + Design_ID, per_unit, mean)
dist <- dist[order(dist$Source, dist$Regime, dist$Design_ID), ]
write.csv(dist, file.path(out_dir, "allocation_distribution_by_design.csv"), row.names = FALSE)

# Figure: six detailed designs and SRS; vector PDF in the default Helvetica, like Figure 4.
ids <- c(8, 6, 4, 2, 5, 1, 9)
# "\u00ad" prints as a true hyphen; R's pdf() device draws "-" as a minus sign.
# Names match CTJ Table 1 and Figure 4 (Title Case).
labels <- c("Incidence\u00adGuided Saturation", "Balanced Quartiles", "Isolation Buffer",
            "High Incidence Focus", "2\u00d72 Blocking", "Checkerboard", "Simple Random Sampling")
f <- dist[dist$Design_ID %in% ids, ]
f$Label <- factor(f$Design_ID, levels = rev(ids), labels = rev(labels))
f$Study <- factor(f$Source, c(grid_source, "nc_refined"), c("Simulation study", "Application study"))
f$Regime <- factor(f$Regime, c("both", "control_only"),
                   c("Spillover reaches both arms", "Spillover reaches controls only"))
f$Benchmark <- f$Design_ID == 9
fig <- ggplot2::ggplot(f, ggplot2::aes(y = Label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
  ggplot2::geom_linerange(ggplot2::aes(xmin = Q05, xmax = Q95, color = Benchmark), linewidth = 0.5) +
  ggplot2::geom_linerange(ggplot2::aes(xmin = Q25, xmax = Q75, color = Benchmark), linewidth = 2.6) +
  ggplot2::geom_point(ggplot2::aes(x = Mean), shape = 21, fill = "white", color = "black", size = 1.8) +
  ggplot2::facet_grid(Study ~ Regime) +
  ggplot2::scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8")) +
  ggplot2::scale_color_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#E69F00"), guide = "none") +
  ggplot2::labs(x = "MSE relative to simple random sampling (log scale)",
                y = NULL) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(face = "bold"),
                 strip.text.y = ggplot2::element_text(angle = -90))
fig_file <- file.path(out_dir, "allocation_distribution_by_design.pdf")
ggplot2::ggsave(fig_file, fig, width = 174, height = 105, units = "mm")
if (profile == "surfaces10") {
  for (folder in c("paper/ctj_manuscript/figures/compact", "paper/thesis_chapter/figures/compact")) {
    stopifnot(file.copy(fig_file, folder, overwrite = TRUE))
  }
}
cat("Wrote allocation distribution summary and figure\n")

