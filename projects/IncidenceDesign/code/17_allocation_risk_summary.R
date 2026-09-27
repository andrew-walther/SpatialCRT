# ============================================================
# Script: 17_allocation_risk_summary.R
# Purpose: Report fixed-X conditional allocation risk and pilot precision.
# Author: Andrew Walther
# Created: 2026-09-27
# Dependencies: base R
# ============================================================
# Rscript code/17_allocation_risk_summary.R pilot (from project directory).
# Outputs are isolated under results/allocation_risk/summary/<profile>/.

# Identity and grouping ----
allocation_unit_keys <- c("Incidence_Mode", "Rho_Incidence", "Surface", "Neighbor_Type",
                          "Rho", "Gamma", "Spillover_Type", "True_Tau", "Estimator")

#' Identify fixed-X scenario units without treating them as population samples
#' @param x Data frame with the requested keys.
#' @param keys Character grouping columns.
#' @return Character key for each row.
#' @examples
#' allocation_summary_key(data.frame(x = c(1, 2)), "x")
allocation_summary_key <- function(x, keys) {
  do.call(paste, c(x[keys], sep = "|"))
}

#' Average conditional metrics equally across the stated fixed scenario units
#'
#' Mean MC SE is sqrt(sum(SE_unit^2))/B: independent design/noise streams conditional
#' on the fixed X surfaces. The signed corrected variances are retained in their
#' average; its square root estimates RMS within-setting allocation SD, not the
#' arithmetic mean SD. Averaged conditional q90 is not a pooled quantile.
#' @param s Complete block-summary data frame from script 16.
#' @param groups Character columns defining reporting strata, excluding Design.
#' @return One row per stratum/design, including MC SE and negative corrections.
#' @examples
#' # allocation_risk_average(results$summary, "Spillover_Type")
allocation_risk_average <- function(s, groups) {
  if (!all(s$Complete) || anyNA(s$Complete)) stop("Incomplete blocks: comparisons withheld")
  keys <- c(groups, "Design", "Design_Name")
  out <- lapply(split(seq_len(nrow(s)), allocation_summary_key(s, keys)), function(i) {
    z <- s[i, ]; v <- mean(z$Variance_Corrected)
    cbind(z[1, keys, drop = FALSE], data.frame(N_Fixed_Units = nrow(z),
      Mean_Allocation_MSE = mean(z$Mean_MSE),
      MCSE_Mean_MSE = sqrt(sum(z$SE_Mean_MSE_Joint^2)) / nrow(z),
      Average_Corrected_Within_Variance = v,
      RMS_Corrected_Within_SD = if (v >= 0) sqrt(v) else NA_real_,
      N_Negative_Variance = sum(z$Variance_Corrected < 0),
      Average_Conditional_Q90 = mean(z$Q90_Estimated),
      Average_Conditional_Worst10 = mean(z$Worst10_Mean_Estimated),
      Observed_Max_Sample_Dependent = max(z$Sampled_Max_Estimated),
      Average_Conditional_Coverage = mean(z$Coverage)))
  })
  do.call(rbind, out)
}

#' Compare ratios of averages and fixed-unit differences against SRS
#' @param s Block summaries, with all designs sharing identical unit keys.
#' @param groups Reporting-stratum columns.
#' @return Average risk, SRS ratios, conditional mean-difference MC SE and unit shares.
#' @examples
#' # allocation_risk_compare(results$summary, "Spillover_Type")
allocation_risk_compare <- function(s, groups) {
  a <- allocation_risk_average(s, groups); base <- s[s$Design == "Design 9", ]
  if (anyDuplicated(allocation_summary_key(base, allocation_unit_keys))) stop("Duplicate SRS units")
  b <- a[a$Design == "Design 9", ]; ix <- match(allocation_summary_key(a, groups),
                                               allocation_summary_key(b, groups))
  if (anyNA(ix)) stop("SRS missing from a reporting stratum")
  metrics <- c("Mean_Allocation_MSE", "Average_Conditional_Q90", "Average_Conditional_Worst10")
  for (metric in metrics) a[[paste0(metric, "_Ratio_To_SRS")]] <- a[[metric]] / b[[metric]][ix]
  a$Mean_Difference_To_SRS <- a$Mean_Allocation_MSE - b$Mean_Allocation_MSE[ix]
  a$MCSE_Difference_To_SRS <- sqrt(a$MCSE_Mean_MSE^2 + b$MCSE_Mean_MSE[ix]^2)
  a$MCSE_Difference_To_SRS[a$Design == "Design 9"] <- 0
  for (j in seq_len(nrow(a))) {
    keep <- s$Design == a$Design[j]
    for (g in groups) keep <- keep & s[[g]] == a[[g]][j]
    z <- s[keep, ]; m <- match(allocation_summary_key(z, allocation_unit_keys),
                              allocation_summary_key(base, allocation_unit_keys))
    if (anyNA(m) || length(m) != sum(allocation_summary_key(base, groups) ==
                                   allocation_summary_key(a[j, , drop = FALSE], groups)))
      stop("Unmatched design/SRS units")
    for (metric in c("Mean_MSE", "Q90_Estimated", "Worst10_Mean_Estimated"))
      a[j, paste0("Share_Fixed_Units_Below_SRS_", metric)] <- mean(z[[metric]] < base[[metric]][m])
  }
  a$Design_Name <- sub("^Design [0-9]+: ", "", a$Design_Name)
  rownames(a) <- NULL; a
}

#' Measure noisy-tail precision with split batches and cross-selected evaluation
#'
#' Half1-selected worst draws are evaluated with Half2 and vice versa. The selected
#' minus held-out difference diagnoses selection/noise optimism, not a deconvolved
#' true-risk tail. Frequency expansion preserves the assignment sampling measure.
#' @param allocations Unique-allocation data frame with Frequency and half MSEs.
#' @param s Matching full-data block summaries.
#' @return One diagnostic row per fixed unit/design.
#' @examples
#' # allocation_risk_precision(results$allocations, results$summary)
allocation_risk_precision <- function(allocations, s) {
  keys <- c(allocation_unit_keys, "Design"); ids <- allocation_summary_key(s, keys)
  parts <- split(seq_len(nrow(allocations)), allocation_summary_key(allocations, keys))
  if (anyDuplicated(ids) || !setequal(ids, names(parts))) stop("Allocation/block identity mismatch")
  out <- lapply(seq_len(nrow(s)), function(j) {
    a <- allocations[parts[[ids[j]]], ]; idx <- rep(seq_len(nrow(a)), a$Frequency)
    h1 <- a$MSE_Half1[idx]; h2 <- a$MSE_Half2[idx]; n <- max(1L, ceiling(length(idx) / 10))
    t1 <- order(h1, decreasing = TRUE)[seq_len(n)]; t2 <- order(h2, decreasing = TRUE)[seq_len(n)]
    selected <- (mean(h1[t1]) + mean(h2[t2])) / 2
    held <- (mean(h2[t1]) + mean(h1[t2])) / 2
    cbind(s[j, c(keys, "Design_Name", "Q90_Estimated", "Q90_Half1", "Q90_Half2",
      "Worst10_Mean_Estimated", "Worst10_Half1", "Worst10_Half2", "Tail_Half_Jaccard",
      "Tail_Half_Rank_Cor", "Tail_Boundary_Uncertain")],
      Cross_Selected_Worst10 = held, Half_Selected_Worst10 = selected,
      Selection_Noise_Optimism = selected - held,
      Q90_Half_Absolute_Difference = abs(s$Q90_Half1[j] - s$Q90_Half2[j]),
      Worst10_Half_Absolute_Difference = abs(s$Worst10_Half1[j] - s$Worst10_Half2[j]))
  })
  do.call(rbind, out)
}

#' Compare a targeted outcome-replication extension with its matching pilot units
#' @param current Extended block summaries.
#' @param pilot Original pilot block summaries.
#' @return Matched full/extended metrics and ratios; descriptive precision diagnostics.
#' @examples
#' # allocation_risk_extension(extended$summary, pilot$summary)
allocation_risk_extension <- function(current, pilot) {
  keys <- c(allocation_unit_keys, "Design")
  if (anyDuplicated(allocation_summary_key(pilot, keys))) stop("Duplicate pilot units")
  m <- match(allocation_summary_key(current, keys), allocation_summary_key(pilot, keys))
  if (anyNA(m)) stop("Extension has units missing from pilot")
  out <- current[c(keys, "Design_Name")]
  metrics <- c("Mean_MSE", "Q90_Estimated", "Worst10_Mean_Estimated", "SD_Corrected",
    "Variance_Corrected", "Variance_MC_Noise", "Variance_Noise_Fraction",
    "SE_Mean_MSE_Joint", "Tail_Half_Jaccard", "Tail_Half_Rank_Cor")
  for (metric in metrics) {
    out[[paste0(metric, "_Pilot")]] <- pilot[[metric]][m]
    out[[paste0(metric, "_Extended")]] <- current[[metric]]
    out[[paste0(metric, "_Change")]] <- current[[metric]] - pilot[[metric]][m]
    out[[paste0(metric, "_Ratio")]] <- current[[metric]] / pilot[[metric]][m]
  }
  out
}

# File runner ----
#' Write transparent conditional pilot tables and an interpretation guide
#' @param input Path to script 16 results.rds.
#' @param output Output directory within results/allocation_risk/summary/.
#' @return Invisibly, named list of summary tables.
#' @examples
#' # allocation_risk_report("../results/allocation_risk/pilot/results.rds", "../results/allocation_risk/summary/pilot")
allocation_risk_report <- function(input, output) {
  r <- readRDS(input); s <- r$summary; a <- r$allocations
  if (!all(s$Complete) || !all(a$Complete) || anyNA(c(s$Complete, a$Complete)))
    stop("Incomplete fits: summary comparisons withheld")
  if (anyDuplicated(allocation_summary_key(s, c(allocation_unit_keys, "Design"))))
    stop("Duplicate summary units")
  # Coverage is averaged within allocation frequency, then equally across fixed units.
  for (j in seq_len(nrow(s))) {
    m <- allocation_summary_key(a, c(allocation_unit_keys, "Design")) ==
      allocation_summary_key(s[j, ], c(allocation_unit_keys, "Design"))
    s$Coverage[j] <- weighted.mean(a$Coverage[m], a$Frequency[m])
  }
  common <- c("Neighbor_Type", "True_Tau", "Estimator", "Spillover_Type")
  tables <- list(pooled_fixed_units = allocation_risk_compare(s, common),
    by_incidence_config = allocation_risk_compare(s, c(common, "Incidence_Mode", "Rho_Incidence")),
    by_parameter_corner = allocation_risk_compare(s, c(common, "Rho", "Gamma")),
    by_config_corner = allocation_risk_compare(s, c(common, "Incidence_Mode", "Rho_Incidence", "Rho", "Gamma")),
    precision_by_fixed_unit = allocation_risk_precision(a, s))
  pilot_file <- file.path(dirname(dirname(input)), "pilot", "results.rds")
  if (basename(dirname(input)) != "pilot" && file.exists(pilot_file)) {
    pilot <- readRDS(pilot_file)
    tables$matched_pilot_replication_comparison <- allocation_risk_extension(s, pilot$summary)
  }
  baseline <- file.path(dirname(input), "main_run_surface_baseline.csv")
  if (file.exists(baseline)) {
    old <- read.csv(baseline); keys <- c(allocation_unit_keys, "Design")
    m <- match(allocation_summary_key(old, keys), allocation_summary_key(s, keys))
    if (anyNA(m)) stop("Main-run baseline contains unmatched units")
    old$Nested_MSE <- s$Mean_MSE[m]; old$Nested_Minus_Main <- old$Nested_MSE - old$MSE
    tables$main_run_descriptive_fixed_units <- old
    bg <- c(common, "Incidence_Mode", "Rho_Incidence", "Design")
    tables$main_run_descriptive_averages <- do.call(rbind,
      lapply(split(seq_len(nrow(old)), allocation_summary_key(old, bg)), function(i) {
        cbind(old[i[1], bg, drop = FALSE], N_Fixed_Units = length(i),
          Main_Run_Mean_MSE = mean(old$MSE[i]), Nested_Mean_MSE = mean(old$Nested_MSE[i]),
          Nested_Minus_Main = mean(old$Nested_Minus_Main[i]))
      }))
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) write.csv(tables[[nm]], file.path(output, paste0(nm, ".csv")), row.names = FALSE)
  flags <- data.frame(Complete = all(a$Complete), N_Unique_Allocations = nrow(a),
    N_Outcome_Fits = sum(a$N_Outcome), N_Aliased = sum(a$N_Aliased), N_Warn = sum(a$N_Warn),
    N_Negative_Corrected_Variance = sum(s$Variance_Corrected < 0),
    N_Structural_Singleton_Blocks = length(unique(allocation_summary_key(
      a[a$Structural_Singleton, ], c(allocation_unit_keys, "Design")))))
  write.csv(flags, file.path(output, "runtime_flags.csv"), row.names = FALSE)
  writeLines(c("Allocation-risk conditional summary", paste("Input:", normalizePath(input)),
    paste("Output slice:", basename(dirname(input)), "; actual surfaces:", paste(sort(unique(s$Surface)), collapse = ",")),
    paste("Actual configs:", length(unique(allocation_summary_key(s, c("Incidence_Mode", "Rho_Incidence")))),
      "; blocks:", nrow(s), "; outcome replicates per unique allocation:", paste(sort(unique(a$N_Outcome)), collapse = ",")),
    "Pooled tables equally average the selected fixed surfaces/configurations and rho/gamma corners within each regime.",
    "The pilot selects the first two surfaces of five configs; targeted extensions may use a strict subset, as recorded above.",
    "MC SEs condition on these fixed X surfaces, with independent block/design streams; no population inference or p-values.",
    "Ratios are ratios of averages. Unit shares are descriptive comparisons on matched fixed settings.",
    "Conditional q90 and worst10 averages concern estimated risks; they are NOT pooled quantiles or deconvolved true tails.",
    "RMS corrected within SD = sqrt(average signed corrected variance). Negative component estimates are counted and retained.",
    "Observed maximum depends on sampled allocations/outcomes; it is not an exhaustive worst case.",
    "Split-half and cross-selected tails diagnose outcome-noise sensitivity and selection optimism.",
    "Main-run differences are descriptive only: the old slice used one outcome per allocation, with 25 fits per surface.",
    paste("Warnings:", flags$N_Warn, "; aliases:", flags$N_Aliased)), file.path(output, "interpretation.txt"))
  print(flags); invisible(tables)
}

if (sys.nframe() == 0L) {
  arg <- commandArgs(TRUE); profile <- if (length(arg)) arg[1] else "pilot"
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)[1]
  project <- dirname(dirname(normalizePath(sub("--file=", "", f))))
  root <- file.path(project, "results", "allocation_risk")
  allocation_risk_report(file.path(root, profile, "results.rds"), file.path(root, "summary", profile))
}
