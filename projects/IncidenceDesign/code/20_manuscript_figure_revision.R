# ============================================================
# Script: 20_manuscript_figure_revision.R
# Purpose: Export compact, traceable manuscript figures from completed studies.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-02
# Dependencies: ggplot2, sf, grid (existing project dependencies)
# ============================================================
# Run from the project root: Rscript code/20_manuscript_figure_revision.R
# This is a read-only analysis export: no simulation, checkpoint or weights change.

# Numerical helpers ----
#' Aggregate grid accuracy without mixing incidence configurations or regimes
#' @param results Completed oracle scenario data frame.
#' @return Data frame of queen, tau=1, configuration/design/regime means.
#' @examples
#' # compact_grid_means(readRDS("results/sim_data/current.rds"))
compact_grid_means <- function(results) {
  x <- results[results$Neighbor_Type == "queen" & results$True_Tau == 1, ]
  stopifnot(all(x$Estimator == "oracle"), all(x$N_Valid_Est == 250),
            all(is.finite(x$MSE)), all(is.finite(x$SD)))
  # SD uses n-1 whereas empirical MSE uses n. This finite-sample correction
  # makes MSE = Bias^2 + ((n-1)/n)*SD^2 exact, not an approximation.
  x$Squared_Bias <- x$Bias^2
  x$Error_Variance <- (x$N_Valid_Est - 1) / x$N_Valid_Est * x$SD^2
  stopifnot(max(abs(x$MSE - x$Squared_Bias - x$Error_Variance)) < 1e-10)
  aggregate(cbind(MSE, Squared_Bias, Error_Variance) ~ Design +
              Incidence_Mode + Rho_Incidence + Spillover_Type, x, mean)
}

#' Compute annual MSE ratios using the year- and regime-matched SRS denominator
#' @param annual Authorized yearly primary design-means data frame.
#' @return Same rows with a recomputed Ratio column; input is not modified.
#' @examples
#' # annual_srs_ratios(read.csv("yearly_primary_design_means.csv"))
annual_srs_ratios <- function(annual) {
  keys <- paste(annual$Year, annual$Regime, sep = "/")
  reference <- annual[annual$Design_ID == 9, ]
  ref_keys <- paste(reference$Year, reference$Regime, sep = "/")
  stopifnot(!anyDuplicated(ref_keys), all(keys %in% ref_keys),
            all(is.finite(annual$Mean_MSE)), all(reference$Mean_MSE > 0))
  annual$Ratio <- annual$Mean_MSE / reference$Mean_MSE[match(keys, ref_keys)]
  stopifnot(max(abs(annual$Ratio - annual$Ratio_of_Mean_MSE_SRS)) < 1e-10)
  annual
}

#' Export four figures and their numerical records for both manuscript families
#' @param project_dir Absolute IncidenceDesign checkout directory.
#' @return Invisibly, the directory containing the four PDFs and three CSVs.
#' @examples
#' # export_compact_exhibits(normalizePath("."))
export_compact_exhibits <- function(project_dir) {
  out <- file.path(project_dir, "results/manuscript_exhibit_revision_20261002")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  grid_file <- file.path(project_dir, "results/sim_data",
                        "sim_results_MLE_tau_sweep_combined_20260927_000502.rds")
  grid_results <- readRDS(grid_file)
  stopifnot(nrow(grid_results) == 14400)
  m <- compact_grid_means(grid_results)
  stopifnot(nrow(m) == 90)
  write.csv(m, file.path(out, "grid_configuration_means.csv"), row.names = FALSE)
  app_root <- file.path(project_dir, "application/results/real_sud_rev_20261002")
  annual <- annual_srs_ratios(read.csv(file.path(app_root, "exhibits",
                                                "yearly_primary_design_means.csv")))
  stopifnot(nrow(annual) == 72, all(table(annual$Year, annual$Regime) == 9))
  write.csv(annual, file.path(out, "annual_srs_ratios.csv"), row.names = FALSE)

  # Plot labels are presentation only; numeric IDs remain in the exported records.
  ids <- c(8, 6, 4, 2, 5, 1)
  short <- c("Guided saturation", "Balanced quartiles", "Isolation buffer",
             "High incidence focus", "2 x 2 blocking", "Checkerboard")
  d <- m[m$Design %in% paste("Design", ids), ]
  d$Label <- factor(d$Design, levels = paste("Design", rev(ids)), labels = rev(short))
  d$Regime <- factor(d$Spillover_Type, levels = c("both", "control_only"),
                     labels = c("Both arms", "Control only"))
  configs <- c("iid/0", "spatial/0.2", "spatial/0.5", "poisson/0.2", "poisson/0.5")
  d$Configuration <- factor(paste(d$Incidence_Mode, d$Rho_Incidence, sep = "/"),
                            levels = configs,
                            labels = c("Independent uniform", "Spatial (0.2)",
                                       "Spatial (0.5)", "Poisson (0.2)", "Poisson (0.5)"))
  stopifnot(!anyNA(d$Configuration))
  theme <- ggplot2::theme_minimal(base_size = 9) + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
    legend.title = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(3, 5, 3, 3))
  colors <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#D55E00")
  p <- ggplot2::ggplot(d, ggplot2::aes(MSE, Label, color = Configuration, shape = Configuration)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 1.8) +
    ggplot2::facet_wrap(~Regime, nrow = 1) +
    ggplot2::scale_x_log10(breaks = c(0.03, 0.1, 0.3, 1, 3)) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = "Mean MSE (log scale)", y = NULL) + theme +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2),
                    shape = ggplot2::guide_legend(nrow = 2))
  ggplot2::ggsave(file.path(out, "grid_accuracy.pdf"), p, width = 174, height = 93, units = "mm")

  # Decomposition: selected configuration, 16 rho/gamma settings in each regime.
  b <- d[d$Incidence_Mode == "poisson" & d$Rho_Incidence == 0.2, ]
  b$Bias_Percent <- 100 * b$Squared_Bias / b$MSE
  b$Variance_Percent <- 100 * b$Error_Variance / b$MSE
  stopifnot(nrow(b) == 12, max(abs(b$Bias_Percent + b$Variance_Percent - 100)) < 1e-10)
  write.csv(b, file.path(out, "grid_bias_variance.csv"), row.names = FALSE)
  components <- rbind(transform(b, Component = "Squared bias", Percent = Bias_Percent),
                      transform(b, Component = "Error variance", Percent = Variance_Percent))
  components$Component <- factor(components$Component, levels = c("Squared bias", "Error variance"))
  p <- ggplot2::ggplot(components, ggplot2::aes(Percent, Label, fill = Component)) +
    ggplot2::geom_col(width = 0.65, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_text(data = b, ggplot2::aes(x = 97, y = Label,
                      label = sprintf("bias %.1f%%", Bias_Percent)), inherit.aes = FALSE,
                      hjust = 1, size = 2.6, color = "#222222") +
    ggplot2::facet_wrap(~Regime, nrow = 1) +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#B9DDE9")) +
    ggplot2::scale_x_continuous(breaks = c(0, 50, 100)) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::labs(x = "Share of mean MSE (%)", y = NULL) + theme
  ggplot2::ggsave(file.path(out, "grid_bias_variance.pdf"), p, width = 174, height = 82, units = "mm")

  # A diverging log-ratio scale is symmetric about SRS=1; printed cells retain
  # actual ratios. Saturation below SRS and very imprecise designs both remain visible.
  names <- c("Graph Checkerboard", "High incidence focus", "Plain saturation", "Isolation buffer",
             "Spatial blocking", "Balanced quartiles", "Balanced halves", "Guided saturation", "SRS")
  annual$Label <- factor(annual$Design_ID, levels = 9:1, labels = rev(names))
  annual$Regime_Label <- factor(annual$Regime, levels = c("both", "control_only"),
                                labels = c("Both arms", "Control only"))
  annual$Log_Ratio <- log2(annual$Ratio)
  span <- max(abs(annual$Log_Ratio))
  p <- ggplot2::ggplot(annual, ggplot2::aes(factor(Year), Label, fill = Log_Ratio)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Ratio)), size = 2.7) +
    ggplot2::facet_wrap(~Regime_Label, nrow = 1) +
    ggplot2::scale_fill_gradient2(low = "#56B4E9", mid = "#F7F7F7", high = "#D55E00",
      midpoint = 0, limits = c(-span, span), breaks = log2(c(0.25, 0.5, 1, 2, 4)),
      labels = c("0.25", "0.5", "1", "2", "4"), name = "MSE / SRS") +
    ggplot2::labs(x = "Observed planning year", y = NULL) + theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.title = ggplot2::element_text(size = 9))
  ggplot2::ggsave(file.path(out, "nc_annual_accuracy.pdf"), p, width = 174, height = 93, units = "mm")

  # Named joins prevent an unnoticed permutation of observed rates or regions.
  # Simplification is display only and never feeds the cached queen/rook weights.
  setup <- readRDS(file.path(app_root, "setup.rds"))
  observed <- read.csv(file.path(app_root, "exhibits/observed_cluster_inputs.csv"))
  observed <- observed[observed$period == 2018, ]
  partition <- read.csv(file.path(app_root, "frozen_partitions.csv"))
  g <- sf::st_simplify(setup$clusters, dTolerance = 250)
  stopifnot(nrow(g) == 58, !anyDuplicated(g$Primary_College),
            setequal(g$Primary_College, observed$Primary_College),
            setequal(g$Primary_College, partition$Primary_College))
  g$Rate <- observed$rate_per_100k[match(g$Primary_College, observed$Primary_College)]
  g$Region <- factor(partition$Region[match(g$Primary_College, partition$Primary_College)])
  stopifnot(all(g$Region == setup$regions), all(is.finite(g$Rate)))
  map_theme <- ggplot2::theme_void(base_size = 9) + ggplot2::theme(
    legend.position = "bottom", legend.key.height = grid::unit(2, "mm"),
    plot.title = ggplot2::element_text(size = 9, face = "bold"),
    plot.margin = ggplot2::margin(2, 2, 2, 2))
  rates <- ggplot2::ggplot(g) + ggplot2::geom_sf(ggplot2::aes(fill = Rate), color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_viridis_c(name = "SUD per 100,000", option = "C") +
    ggplot2::labs(title = "A. Observed 2018 rates") + map_theme
  regions <- ggplot2::ggplot(g) + ggplot2::geom_sf(ggplot2::aes(fill = Region), color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"), name = "Region") +
    ggplot2::labs(title = "B. Frozen saturation regions") + map_theme
  grDevices::pdf(file.path(out, "nc_planning_map.pdf"), width = 174/25.4, height = 64/25.4)
  grid::grid.newpage()
  print(rates, vp = grid::viewport(x = 0.25, width = 0.5))
  print(regions, vp = grid::viewport(x = 0.75, width = 0.5))
  grDevices::dev.off()
  for (folder in c("ctj_manuscript", "dissertation_chapter")) {
    destination <- file.path(project_dir, "paper", folder, "figures/compact")
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    selected <- if (folder == "ctj_manuscript") list.files(out, pattern = "[.]pdf$", full.names = TRUE) else
      file.path(out, "nc_annual_accuracy.pdf")
    stopifnot(all(file.copy(selected, destination, overwrite = TRUE)))
  }
  writeLines(c("Read-only export of completed oracle/primary aggregate results.",
               "PDF widths: 174 mm; fonts 9 pt (cell labels approximately 7.7 pt).",
               "Grid accuracy: six-design subset, tau=1, queen; 16 settings/configuration/regime.",
               "Bias/variance: Poisson rho_X=0.2, each regime separate; exact empirical decomposition.",
               "NC means: 6 queen parameter pairs/year/regime/design; shared years not independent.",
               "Map: observed 2018 rates, authorized 58-cluster aggregates, frozen partition."),
             file.path(out, "sources.txt"))
  invisible(out)
}

# Entry point ----
if (sys.nframe() == 0L) {
  script_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
  export_compact_exhibits(dirname(dirname(normalizePath(script_file))))
  cat("PASS: four compact exhibits and numerical records exported; simulation unchanged.\n")
}
