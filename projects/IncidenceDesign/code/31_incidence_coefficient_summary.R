# ============================================================
# Script: 31_incidence_coefficient_summary.R
# Purpose: Summarize and plot the incidence coefficient sweep after reproduction.
# Author: Andrew Walther
# Created: 2026-11-25
# Dependencies: dplyr, ggplot2
# ============================================================
#
# From the IncidenceDesign root after the full run:
#   Rscript code/31_incidence_coefficient_summary.R
# Document notation: β multiplies Z and τ multiplies X.

suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
source("code/00_design_names.R")

summary_keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Spillover_Type",
                  "Estimator", "Coefficient", "Incidence_Coef", "Design")

#' Average matched scenarios and compare each design's mean MSE to SRS
#' @param d Scenario data frame from the incidence coefficient runner.
#' @return One performance row per incidence configuration, neighbor type, regime,
#'   estimator, coefficient, true τ, and design.
#' @examples
#' # summarize_incidence(readRDS("results/incidence_coef_2026-11/full/scenarios.rds"))
summarize_incidence <- function(d) {
  out <- d |>
    group_by(across(all_of(summary_keys))) |>
    summarise(
      N_Scenarios = n(), N_Attempted = sum(N_Attempted), N_Valid_Est = sum(N_Valid_Est),
      Bias = mean(Bias), Mean_Scenario_SD = mean(SD), MSE = mean(MSE),
      Coverage_Tdf = mean(Coverage_Tdf), Power_Tdf = mean(Power_Tdf),
      Mean_Model_SE = mean(Mean_Model_SE),
      # Pool paired squared differences over valid pairs, then take the square root.
      Paired_RMS_Change = if (sum(N_Paired, na.rm = TRUE) > 0L) {
        sqrt(sum(ifelse(N_Paired > 0L, Paired_MSD * N_Paired, 0)) / sum(N_Paired))
      } else NA_real_,
      N_Paired = if (all(is.na(N_Paired))) NA_integer_ else sum(N_Paired),
      Mean_Abs_Cor_ZX = if (sum(N_Cor_ZX) > 0L) {
        sum(ifelse(N_Cor_ZX > 0L, Mean_Abs_Cor_ZX * N_Cor_ZX, 0)) / sum(N_Cor_ZX)
      } else NA_real_,
      N_Cor_ZX = sum(N_Cor_ZX), .groups = "drop"
    )
  srs_keys <- setdiff(summary_keys, "Design")
  ref <- out |>
    filter(Design == "Design 9") |>
    select(all_of(srs_keys), MSE_SRS = MSE)
  out <- left_join(out, ref, by = srs_keys, relationship = "many-to-one")
  stopifnot(all(is.finite(out$MSE)), all(is.finite(out$MSE_SRS) & out$MSE_SRS > 0))
  out |> mutate(MSE_Ratio_SRS = MSE / MSE_SRS)
}

#' Compute uncertainty across surfaces for the summaries and paired changes
#' @param d Surface data frame from the incidence coefficient runner.
#' @param n_surfaces Expected number of surfaces per incidence configuration.
#' @return Summary keys with Monte Carlo SEs and changes in MSE and model SE
#'   relative to τ = 1. The changes pair surface means at the two τ values.
#' @examples
#' # summarize_surface_mc(readRDS("results/incidence_coef_2026-11/full/surfaces.rds"), 10)
summarize_surface_mc <- function(d, n_surfaces) {
  surface_keys <- c(summary_keys, "Surface")
  metrics <- c("Bias", "MSE", "Coverage_Tdf", "Power_Tdf", "Mean_Model_SE")
  by_surface <- d |>
    group_by(across(all_of(surface_keys))) |>
    summarise(across(all_of(metrics), mean), .groups = "drop")
  pair_keys <- setdiff(surface_keys, "Incidence_Coef")
  ref <- by_surface |>
    filter(Incidence_Coef == 1) |>
    select(all_of(pair_keys), Reference_MSE = MSE, Reference_Model_SE = Mean_Model_SE)
  by_surface <- left_join(by_surface, ref, by = pair_keys, relationship = "many-to-one")
  stopifnot(all(is.finite(by_surface$Reference_MSE)), all(is.finite(by_surface$Reference_Model_SE)))
  out <- by_surface |>
    mutate(Change_MSE = MSE - Reference_MSE, Change_Model_SE = Mean_Model_SE - Reference_Model_SE) |>
    group_by(across(all_of(summary_keys))) |>
    summarise(
      N_Surfaces = n(),
      across(all_of(c(metrics, "Change_MSE", "Change_Model_SE")),
             ~ sd(.x) / sqrt(n()), .names = "MC_SE_{.col}"),
      Change_MSE = mean(Change_MSE), Change_Model_SE = mean(Change_Model_SE), .groups = "drop"
    )
  stopifnot(all(out$N_Surfaces == n_surfaces))
  out
}

#' Write the performance tables and title-free oracle figures
#' @param in_dir Profile directory with completed, checked simulation results.
#' @param out_dir Destination for the CSVs and PDFs.
#' @return Invisibly, the aggregated scenario metrics with their Monte Carlo SEs.
#' @examples
#' # write_incidence_summary("results/incidence_coef_2026-11/full")
write_incidence_summary <- function(in_dir, out_dir = file.path(in_dir, "summary")) {
  input_names <- c("reproduction_check.txt", "run_info.rds", "scenarios.rds", "surfaces.rds")
  inputs <- setNames(file.path(in_dir, input_names), input_names)
  checksums <- tools::md5sum(inputs)
  stopifnot(!anyNA(checksums))
  verification <- readLines(inputs[["reproduction_check.txt"]])
  stopifnot(identical(tail(verification, 1), "PASS: all production metrics match exactly."))
  info <- readRDS(inputs[["run_info.rds"]])
  scen <- readRDS(inputs[["scenarios.rds"]])
  surf <- readRDS(inputs[["surfaces.rds"]]) |>
    select(all_of(c(summary_keys, "Surface", "True_Tau", "Bias", "MSE",
                    "Coverage_Tdf", "Power_Tdf", "Mean_Model_SE")))
  design_names <- DESIGN_FULL_NAMES
  design_names["Design 1"] <- "Checkerboard"
  design_names["Design 5"] <- "2×2 Blocking"
  stopifnot(setequal(unique(scen$Design), names(design_names)))
  slices <- list(All_Beta = rep(TRUE, nrow(scen)), Beta_1 = scen$True_Tau == 1)
  all_results <- lapply(names(slices), function(slice) {
    sc <- scen[slices[[slice]], ]
    su <- if (slice == "All_Beta") surf else filter(surf, True_Tau == 1)
    perf <- summarize_incidence(sc)
    mc <- summarize_surface_mc(su, info$manifest$params$n_surfaces)
    stopifnot(nrow(perf) == nrow(mc))
    left_join(perf, mc, by = summary_keys, relationship = "one-to-one") |>
      mutate(Intervention_Grid = slice, Design_Name = unname(design_names[Design]))
  })
  results <- bind_rows(all_results)
  stopifnot(all(is.finite(results$N_Surfaces)))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(results, file.path(out_dir, "performance.csv"), row.names = FALSE)

  # Base PDF maps the ASCII hyphen to a minus sign in some fonts.
  plot_names <- gsub("-", "\u00ad", design_names, fixed = TRUE)
  primary <- results |>
    filter(Intervention_Grid == "All_Beta", Estimator == "oracle") |>
    mutate(
      Design_Name = factor(unname(plot_names[Design]), levels = unname(plot_names)),
      Incidence = factor(
        case_when(Incidence_Mode == "iid" ~ "Independent uniform",
                  Incidence_Mode == "spatial" ~ sprintf("Spatial, ρₓ = %.1f", Rho_Incidence),
                  Incidence_Mode == "poisson" ~ sprintf("Poisson, ρₓ = %.1f", Rho_Incidence)),
        levels = c("Independent uniform", "Spatial, ρₓ = 0.2", "Spatial, ρₓ = 0.5",
                   "Poisson, ρₓ = 0.2", "Poisson, ρₓ = 0.5")
      ),
      Regime = factor(Spillover_Type, levels = c("control_only", "both"),
                      labels = c("Control\u00adonly spillover", "Both\u00adarms spillover"))
    )
  stopifnot(!anyNA(primary$Incidence), !anyNA(primary$Regime))
  colours <- setNames(c("#A3214B", "#AD4B00", "#127187", "#6F5AA1", "#55742A",
                        "#AF7300", "#745348", "#0060A7", "#272727"), plot_names)
  shapes <- setNames(c(0, 17, 15, 5, 8, 6, 3, 18, 1), plot_names)
  design_scales <- list(scale_colour_manual(values = colours), scale_shape_manual(values = shapes))
  rook_labels <- unname(plot_names)
  rook_labels[plot_names == "Checkerboard"] <- "Checkerboard (aliased Z)"
  rook_z_scales <- list(scale_colour_manual(values = colours, labels = rook_labels,
                                            breaks = unname(plot_names)),
                        scale_shape_manual(values = shapes, labels = rook_labels,
                                           breaks = unname(plot_names)))
  appearance <- list(
    facet_grid(Incidence ~ Regime, labeller = labeller(Incidence = as_labeller(
      setNames(c('"Independent uniform"', 'paste("Spatial, ", rho[X] == 0.2)',
                 'paste("Spatial, ", rho[X] == 0.5)', 'paste("Poisson, ", rho[X] == 0.2)',
                 'paste("Poisson, ", rho[X] == 0.5)'), levels(primary$Incidence)), label_parsed))),
    theme_bw(base_size = 10),
    labs(colour = NULL, shape = NULL),
    theme(legend.position = "bottom", panel.grid.minor = element_blank()),
    guides(colour = guide_legend(ncol = 3))
  )
  tau_scale <- scale_x_continuous(breaks = sort(unique(primary$Incidence_Coef)),
                                  labels = function(x) format(x, trim = TRUE, drop0trailing = TRUE))
  tau_label <- expression("Incidence coefficient (" * tau * ")")
  se_label <- expression("Mean model SE of " * hat(tau))
  for (nb in c("queen", "rook")) {
    z <- filter(primary, Coefficient == "Z", Neighbor_Type == nb)
    x <- filter(primary, Coefficient == "X", Neighbor_Type == nb)
    x_at_one <- filter(x, Incidence_Coef == 1)
    stopifnot(nrow(z) > 0, nrow(x) > 0, nrow(x_at_one) > 0, all(z$MSE_Ratio_SRS > 0),
              all(is.finite(x$Power_Tdf)), all(is.finite(x$Mean_Model_SE)),
              all(is.finite(x_at_one$Mean_Abs_Cor_ZX)))
    beta_ratio <- ggplot(z, aes(Incidence_Coef, MSE_Ratio_SRS, colour = Design_Name,
                                shape = Design_Name)) +
      geom_hline(yintercept = 1, colour = "grey45", linetype = "dashed") +
      geom_line(linewidth = 0.45) + geom_point(size = 1.5) +
      scale_y_log10() + tau_scale + appearance +
      (if (nb == "rook") rook_z_scales else design_scales) +
      labs(x = tau_label, y = "Intervention MSE relative to SRS (log scale)")
    tau_power <- ggplot(x, aes(Incidence_Coef, Power_Tdf, colour = Design_Name,
                               shape = Design_Name)) +
      geom_hline(yintercept = 0.05, colour = "grey45", linetype = "dashed") +
      geom_line(linewidth = 0.45) + geom_point(size = 1.5) +
      scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(p) paste0(100 * p, "%"),
                         limits = c(0, 1)) +
      tau_scale + appearance + design_scales + labs(x = tau_label, y = "Two\u00adsided rejection rate")
    tau_se <- ggplot(x, aes(Incidence_Coef, Mean_Model_SE, colour = Design_Name,
                            shape = Design_Name)) +
      geom_line(linewidth = 0.45) + geom_point(size = 1.5) +
      tau_scale + appearance + design_scales + labs(x = tau_label, y = se_label)
    tau_alignment <- ggplot(x_at_one, aes(Mean_Abs_Cor_ZX, Mean_Model_SE, colour = Design_Name,
                                         shape = Design_Name)) +
      geom_point(size = 2.1) + appearance + design_scales +
      labs(x = "Mean absolute correlation of treatment and incidence", y = se_label)

    plots <- list(intervention_mse_vs_srs = beta_ratio, incidence_rejection = tau_power,
                  incidence_se = tau_se, incidence_se_vs_alignment = tau_alignment)
    for (nm in names(plots)) {
      ggsave(file.path(out_dir, sprintf("%s_%s.pdf", nm, nb)), plots[[nm]],
             width = 8.5, height = 10, units = "in", device = grDevices::pdf)
    }
  }
  if (!identical(checksums, tools::md5sum(inputs))) stop("Simulation inputs changed during summary assembly.")
  versions <- vapply(c("dplyr", "ggplot2"), function(p) as.character(packageVersion(p)), character(1))
  sources <- tools::md5sum(c("code/31_incidence_coefficient_summary.R", "code/00_design_names.R"))
  writeLines(c(R.version.string, paste(names(versions), versions), "Input and source MD5 hashes",
               paste(c(input_names, names(sources)), c(checksums, sources))),
             file.path(out_dir, "provenance.txt"))
  invisible(results)
}

if (sys.nframe() == 0L) write_incidence_summary("results/incidence_coef_2026-11/full")
