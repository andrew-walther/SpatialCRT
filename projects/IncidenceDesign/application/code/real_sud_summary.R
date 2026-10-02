# ============================================================
# Script: real_sud_summary.R
# Purpose: Extract traceable yearly comparisons and publication figures.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud modules, ggplot2
# ============================================================
real_sud_define_only <- TRUE
.rs_summary_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(normalizePath(dirname(f)))
  }
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
})
source(file.path(.rs_summary_dir, "run_real_sud.R"))

#' Compare design and SRS risks within the same planning/parameter setting
#' @param p Complete reporting table for one model, neighbor and summary rule.
#' @return Matched table with mean-MSE ratios, differences and within-setting MC SE.
#' @family real_sud_summary
#' @seealso rs_exhibits
#' @examples
#' # rs_srs_comparison(primary)
rs_srs_comparison <- function(p) {
  keys <- c("Year", "Model", "Neighbor", "Summary", "Rho", "Gamma", "Regime")
  ref <- p[p$Design_ID == 9, c(keys, "Mean_MSE", "SE_Mean_MSE_Joint", "Source_ID")]
  names(ref)[!names(ref) %in% keys] <- paste0("SRS_", names(ref)[!names(ref) %in% keys])
  p$Row_ID <- seq_len(nrow(p))
  out <- merge(p, ref, by = keys, all.x = TRUE, sort = FALSE)
  out <- out[order(out$Row_ID), ]; rownames(out) <- NULL
  if (nrow(out) != nrow(p) || anyNA(out$SRS_Source_ID)) stop("Incomplete/duplicate SRS matching")
  out$MSE_Ratio_SRS <- out$Mean_MSE / out$SRS_Mean_MSE
  out$MSE_Difference_SRS <- out$Mean_MSE - out$SRS_Mean_MSE
  # Design supports have independent allocation/noise streams within one setting.
  # This formula is not an MC SE for pooled setting means or cross-model contrasts.
  out$MCSE_Difference_SRS <- sqrt(out$SE_Mean_MSE_Joint^2 + out$SRS_SE_Mean_MSE_Joint^2)
  same <- out$Source_ID == out$SRS_Source_ID
  out$MCSE_Difference_SRS[same] <- 0
  out$MC_Lower_Difference <- out$MSE_Difference_SRS - 1.96 * out$MCSE_Difference_SRS
  out$MC_Upper_Difference <- out$MSE_Difference_SRS + 1.96 * out$MCSE_Difference_SRS
  out
}

#' Match each sensitivity to primary queen results at identical parameter values
#' @param p Complete reporting table, including the primary reference rows.
#' @return Original-order table with matched primary risks and ratios; no pooled MC SE.
#' @family real_sud_summary
#' @seealso rs_exhibits
#' @examples
#' # rs_primary_reference(all_results)
rs_primary_reference <- function(p) {
  keys <- c("Year", "Rho", "Gamma", "Regime", "Design_ID")
  ref <- subset(p, Model == "education" & Neighbor == "queen" & Summary == "mean_rank")
  ref <- ref[c(keys, "Mean_MSE", "Bias", "Coverage", "Source_ID")]
  if (anyDuplicated(ref[keys])) stop("Duplicate primary reference")
  names(ref)[!names(ref) %in% keys] <- paste0("Primary_", names(ref)[!names(ref) %in% keys])
  p$Comparison_ID <- seq_len(nrow(p))
  out <- merge(p, ref, by = keys, all.x = TRUE, sort = FALSE)
  out <- out[order(out$Comparison_ID), ]; rownames(out) <- NULL
  if (nrow(out) != nrow(p) || anyNA(out$Primary_Source_ID)) stop("Missing matched primary reference")
  out$MSE_Ratio_Matched_Primary <- out$Mean_MSE / out$Primary_Mean_MSE
  out
}

#' Save yearly descriptive tables, full setting comparisons and traced exhibits
#' @param root Study root with independently verified production outputs.
#' @return List of yearly, setting-level and optional refined-tail tables.
#' @family real_sud_summary
#' @seealso rs_srs_comparison
#' @examples
#' # rs_exhibits()
rs_exhibits <- function(root = file.path(.rs_app, "results", "real_sud_rev_20261002")) {
  production <- file.path(root, "production")
  if (!file.exists(file.path(production, "verification.txt"))) stop("Run independent production verification before extracts")
  b <- readRDS(file.path(production, "results.rds")); p <- b$summary
  if (any(!p$Complete) || any(!p$Precision_OK)) stop("Incomplete/precision-unresolved production results")
  setup <- readRDS(file.path(root, "setup.rds")); rs_manifest(production, setup, b$config)
  dest <- file.path(root, "exhibits"); dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  primary <- subset(p, Model == "education" & Neighbor == "queen" & Summary == "mean_rank")
  comparisons <- rs_srs_comparison(primary)
  write.csv(comparisons, file.path(dest, "primary_setting_srs_comparison.csv"), row.names = FALSE)
  yearly <- do.call(rbind, lapply(split(comparisons, list(comparisons$Year, comparisons$Regime, comparisons$Design_ID), drop = TRUE), function(x) {
    if (nrow(x) != 6) stop("Expected six rho/gamma pairs per year/regime/design")
    data.frame(Year = x$Year[1], Regime = x$Regime[1], Design_ID = x$Design_ID[1], Design = x$Design[1],
      Mean_MSE = mean(x$Mean_MSE), Min_MSE = min(x$Mean_MSE), Max_MSE = max(x$Mean_MSE),
      Bias = mean(x$Bias), Coverage = mean(x$Coverage), Min_Coverage = min(x$Coverage), Max_Coverage = max(x$Coverage),
      Ratio_of_Mean_MSE_SRS = mean(x$Mean_MSE) / mean(x$SRS_Mean_MSE),
      Min_Setting_Ratio_SRS = min(x$MSE_Ratio_SRS), Max_Setting_Ratio_SRS = max(x$MSE_Ratio_SRS),
      N_Point_Lower_SRS = sum(x$MSE_Difference_SRS < 0), N_MC_Upper_Below_Zero = sum(x$MC_Upper_Difference < 0),
      Mean_Treated = mean(x$Mean_Treated), Min_Treated = min(x$Min_Treated), Max_Treated = max(x$Max_Treated),
      Mean_Population_Share = mean(x$Mean_Population_Share),
      Min_Population_Share = min(x$Min_Population_Share), Max_Population_Share = max(x$Max_Population_Share),
      Max_Relative_MCSE_MSE = max(x$Relative_MCSE_MSE), Max_MCSE_Coverage = max(x$MCSE_Coverage),
      Source_IDs = paste(unique(x$Source_ID), collapse = " | "),
      Averaging = "Equal weight to six rho/gamma pairs; descriptive, shared allocation streams")
  }))
  rownames(yearly) <- NULL
  write.csv(yearly, file.path(dest, "yearly_primary_design_means.csv"), row.names = FALSE)
  matched <- rs_primary_reference(p)
  write.csv(matched, file.path(dest, "sensitivity_setting_matched_comparison.csv"), row.names = FALSE)
  sensitivity <- do.call(rbind, lapply(split(matched, list(matched$Model, matched$Neighbor, matched$Summary, matched$Year, matched$Regime, matched$Design_ID), drop = TRUE), function(x) {
    data.frame(Model = x$Model[1], Neighbor = x$Neighbor[1], Summary = x$Summary[1], Year = x$Year[1], Regime = x$Regime[1],
      Design_ID = x$Design_ID[1], Design = x$Design[1], Mean_MSE = mean(x$Mean_MSE), Bias = mean(x$Bias), Coverage = mean(x$Coverage),
      Mean_Matched_Primary_MSE = mean(x$Primary_Mean_MSE),
      Ratio_of_Matched_Mean_MSE = mean(x$Mean_MSE) / mean(x$Primary_Mean_MSE),
      N_Settings = nrow(x), Source_IDs = paste(unique(x$Source_ID), collapse = " | "),
      Primary_Source_IDs = paste(unique(x$Primary_Source_ID), collapse = " | "))
  }))
  write.csv(sensitivity, file.path(dest, "yearly_sensitivity_means.csv"), row.names = FALSE)
  labels <- unname(application_design_names(1:9, "short")); labels[1] <- "Graph Checkerboard"; labels[5] <- "Spatial Blocking"
  yearly$Design_Label <- factor(yearly$Design_ID, levels = 9:1, labels = rev(labels))
  yearly$Regime_Label <- ifelse(yearly$Regime == "both", "Both arms", "Control only")
  comparisons$Design_Label <- factor(comparisons$Design_ID, levels = 1:9, labels = labels)
  comparisons$Parameter <- paste0("rho=", comparisons$Rho, ", gamma=", comparisons$Gamma)
  comparisons$Regime_Label <- ifelse(comparisons$Regime == "both", "Both arms", "Control only")
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11))
  g <- ggplot2::ggplot(yearly, ggplot2::aes(x = Mean_MSE, y = Design_Label)) +
    ggplot2::geom_point(ggplot2::aes(color = Design_ID == 9), size = 2.5) + ggplot2::facet_grid(Year ~ Regime_Label) +
    ggplot2::scale_x_log10() + ggplot2::scale_color_manual(values = c("FALSE" = "#167c80", "TRUE" = "#b77812"), guide = "none") +
    ggplot2::labs(x = "Mean MSE across six rho/gamma settings (log scale)", y = NULL,
      title = "Primary education-outcome performance by observed planning year", subtitle = "Equal-weight descriptive parameter averages; gold identifies SRS. No pooled MC intervals.")
  ggplot2::ggsave(file.path(dest, "primary_yearly_mse.pdf"), g, width = 11, height = 10)
  ggplot2::ggsave(file.path(dest, "primary_yearly_mse.png"), g, width = 11, height = 10, dpi = 150)
  g <- ggplot2::ggplot(comparisons, ggplot2::aes(x = Design_Label, y = Parameter, fill = log10(MSE_Ratio_SRS))) +
    ggplot2::geom_tile() + ggplot2::facet_grid(Year ~ Regime_Label) +
    ggplot2::scale_fill_gradient2(low = "#167c80", mid = "white", high = "#b95549", midpoint = 0,
      breaks = log10(c(0.5, 1, 2, 4)), labels = c("0.5", "1", "2", "4"), name = "MSE / SRS") +
    ggplot2::labs(x = NULL, y = NULL, title = "Every primary queen setting compared with SRS", subtitle = "Teal: lower MSE; red: higher MSE. Setting-specific MC errors are in the source table.") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1, size = 8))
  ggplot2::ggsave(file.path(dest, "primary_setting_ratios.pdf"), g, width = 13, height = 10)
  ggplot2::ggsave(file.path(dest, "primary_setting_ratios.png"), g, width = 13, height = 10, dpi = 150)
  g <- ggplot2::ggplot(yearly[yearly$Regime == "both", ], ggplot2::aes(x = Mean_Population_Share, y = Design_Label)) +
    ggplot2::geom_segment(ggplot2::aes(x = Min_Population_Share, xend = Max_Population_Share, yend = Design_Label), color = "#b9cdd1") +
    ggplot2::geom_point(color = "#167c80") + ggplot2::facet_wrap(~Year, ncol = 2) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, color = "#59718c") +
    ggplot2::labs(x = "Treated working-age population share", y = NULL, title = "Cluster budgets and population budgets differ",
      subtitle = "Both-arms run budgets (support is regime-independent); points: mean; lines: sampled allocation range.")
  ggplot2::ggsave(file.path(dest, "population_shares.pdf"), g, width = 11, height = 7)
  ggplot2::ggsave(file.path(dest, "population_shares.png"), g, width = 11, height = 7, dpi = 150)
  tail <- NULL
  tp <- file.path(root, "tail_confirmation", "performance.csv")
  if (file.exists(tp)) {
    if (!file.exists(file.path(dirname(tp), "verification.txt"))) stop("Tail outputs need independent verification")
    tail <- read.csv(tp, stringsAsFactors = FALSE, colClasses = c(Year = "character"))
    if (any(!tail$Complete)) stop("Incomplete tail comparisons")
    write.csv(tail, file.path(dest, "refined_tail_comparison.csv"), row.names = FALSE)
    a <- tail[tail$Rho == 0.5, ]; a$Design_Label <- factor(a$Design_ID, levels = 9:1, labels = rev(labels))
    a$Regime_Label <- ifelse(a$Regime == "both", "Both arms", "Control only")
    g <- ggplot2::ggplot(a, ggplot2::aes(y = Design_Label)) +
      ggplot2::geom_segment(ggplot2::aes(x = Mean_MSE, xend = Q90_Estimated, yend = Design_Label), color = "#b9cdd1") +
      ggplot2::geom_point(ggplot2::aes(x = Mean_MSE, shape = "Mean"), color = "#167c80", size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = Q90_Estimated, shape = "Estimated q90"), color = "#b77812", size = 2) +
      ggplot2::facet_grid(Year ~ Regime_Label) + ggplot2::scale_x_log10() +
      ggplot2::scale_shape_manual(values = c(Mean = 16, "Estimated q90" = 17), name = NULL) +
      ggplot2::labs(x = "Conditional allocation MSE (log scale)", y = NULL, title = "Average accuracy and allocation downside are different criteria",
        subtitle = "rho=0.5, gamma=0.8; same sampled allocations, R >= 400. Finite-sample tails remain uncertain.")
    ggplot2::ggsave(file.path(dest, "refined_allocation_risk.pdf"), g, width = 11, height = 10)
    ggplot2::ggsave(file.path(dest, "refined_allocation_risk.png"), g, width = 11, height = 10, dpi = 150)
  }
  physical <- p[!duplicated(p$Source_ID), ]
  # Authorized aggregate inputs can be inspected without restricted county files.
  input_cols <- c("Primary_College", "period", "deaths", "person_years", "rate_per_100k", "n_counties", "rank", "X")
  write.csv(do.call(rbind, lapply(setup$annual, function(a) a[input_cols])),
    file.path(dest, "observed_cluster_inputs.csv"), row.names = FALSE)
  block_diagnostics <- do.call(rbind, lapply(as.character(2018:2021), function(y) {
    sizes <- table(setup$blocks[[y]])
    data.frame(Year = y, Block = names(sizes), N_Clusters = as.integer(sizes),
      N_Treated = as.integer(round(sizes / 2)), Singleton_Always_Control = sizes == 1)
  }))
  write.csv(block_diagnostics, file.path(dest, "spatial_block_diagnostics.csv"), row.names = FALSE)
  saveRDS(list(source = rs_hash_files(file.path(.rs_summary_dir, "real_sud_summary.R")),
    inputs = rs_hash_files(c(file.path(root, "setup.rds"), file.path(production, "results.rds"),
      if (!is.null(tail)) file.path(root, "tail_confirmation", "results.rds")))),
    file.path(dest, "exhibit_manifest.rds"))
  lines <- c("NC application exhibit extract", paste("Reporting rows:", nrow(p)), paste("Distinct sources:", nrow(physical)),
    paste("Physical independent outcome fits:", sum(physical$N_Independent_Fits)),
    paste("Maximum relative MCSE(mean MSE):", max(p$Relative_MCSE_MSE)),
    paste("Maximum MCSE(coverage):", max(p$MCSE_Coverage)),
    paste("J range:", paste(range(p$J), collapse = " to ")), paste("R range:", paste(range(p$R), collapse = " to ")),
    "Yearly tables equally average rho/gamma pairs. Reused years and shared allocation streams are not independent replications.",
    "Tail extension refines the same sampled allocations; it is not a fresh allocation confirmation.",
    "Primary yearly MSE: yearly_primary_design_means.csv, Mean_MSE; full parameter ratios: primary_setting_srs_comparison.csv.",
    "Sensitivity comparisons: yearly_sensitivity_means.csv; rook/queen contrasts match the same rho/gamma corners.",
    "No independent-error pooled or cross-model MC intervals are supplied: allocation streams may be shared.",
    "Tail plot/table: refined_tail_comparison.csv if present.")
  writeLines(lines, file.path(dest, "manuscript_extract.txt"))
  write.csv(data.frame(Exhibit = c("primary_yearly_mse", "primary_setting_ratios", "population_shares", "refined_allocation_risk"),
    Source = c("yearly_primary_design_means.csv", "primary_setting_srs_comparison.csv", "yearly_primary_design_means.csv", "refined_tail_comparison.csv"),
    Proposed_Placement = c("Chapter application body; CTJ main/SI after exhibit review", "Appendix/SI", "Body table plus appendix figure", "Chapter application body; CTJ main/SI after exhibit review")),
    file.path(dest, "exhibit_inventory.csv"), row.names = FALSE)
  cat(paste(lines, collapse = "\n"), "\n")
  invisible(list(yearly = yearly, settings = comparisons, sensitivity = sensitivity, tail = tail))
}
if (!exists("real_sud_summary_define_only") || !isTRUE(real_sud_summary_define_only)) rs_exhibits()
