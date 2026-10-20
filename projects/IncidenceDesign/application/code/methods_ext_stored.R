# ============================================================
# Script: methods_ext_stored.R
# Purpose: Methods-extension analyses on the stored NC production results:
#   A0 power at tau = 1, A1 burden share treated (N = 58), A2 formal design
#   comparisons, A5 Stage 1 coverage corrections. No new simulation.
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: methods_ext_core.R, PMCMRplus, ggplot2
# ============================================================
# Run from the project root:  Rscript application/code/methods_ext_stored.R
# Reads application/results/real_sud_rev_20261002/production/{performance.csv,
# checkpoints/} (read-only). Writes application/results/methods_ext_2026-10/stored/.
# The design-test helpers (me_design_tests, me_cd_plot) are reused by the A3
# summaries (methods_ext_summaries.R) for the per-N comparisons.

.ms_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(p)) return(normalizePath(dirname(p)))
  }
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else stop("Run or source by path")
})
if (!exists("me_seed")) source(file.path(.ms_dir, "methods_ext_core.R"))
suppressPackageStartupMessages({ library(ggplot2); library(PMCMRplus) })

# Design names (manuscript convention, author 2026-10-12) ----
ME_DESIGNS <- c("1" = "Checkerboard", "2" = "High Incidence Focus", "3" = "Plain Saturation",
                "4" = "Isolation Buffer", "5" = "2x2 Blocking", "6" = "Balanced Quartiles",
                "7" = "Balanced Halves", "8" = "Incidence-Guided Saturation",
                "9" = "Simple Random Sampling")

#' Map design IDs to manuscript names
#' @param id Integer design IDs.
#' @return Character names.
#' @family methods_ext_stored
#' @examples
#' me_design_name(9)
me_design_name <- function(id) unname(ME_DESIGNS[as.character(id)])

# Formal design comparisons (adapted from code/10, not edited) ----
#' Friedman, Nemenyi and Holm-adjusted paired Wilcoxon tests over matched blocks
#'
#' Same tests and formulas as code/10_statistical_comparisons.R
#' (run_friedman_test, run_nemenyi_posthoc, run_pairwise_wilcoxon), but with the
#' NC block variables given explicitly. Within each block the designs are ranked
#' by the metric (rank 1 = lowest MSE). Nemenyi critical difference:
#'   CD = q_{alpha,k,Inf} / sqrt(2) * sqrt(k (k + 1) / (6 N)).
#' P-values are a consistency summary, not trial-level evidence: the blocks share
#' one map, the years are correlated, and allocation streams are shared.
#' @param d Data frame with Design, the metric and block columns; one row per
#'   design per block.
#' @param block_cols Block-defining columns.
#' @param metric Metric column (default Mean_MSE).
#' @return List: friedman (statistic, df, p), avg_ranks, nemenyi p matrix, cd,
#'   wilcoxon pair table (raw and Holm p, median difference), n_blocks, k.
#' @family methods_ext_stored
#' @examples
#' # me_design_tests(perf, c("Year", "Rho", "Gamma"))
me_design_tests <- function(d, block_cols, metric = "Mean_MSE") {
  d$Block <- interaction(d[block_cols], drop = TRUE)
  wide <- tapply(d[[metric]], list(d$Block, d$Design), function(x) { if (length(x) != 1) stop("Not one row per design per block"); x })
  if (anyNA(wide)) stop("Incomplete block x design matrix")
  ft <- stats::friedman.test(wide)
  rk <- t(apply(wide, 1, rank, ties.method = "average"))
  avg <- sort(colMeans(rk))
  nt <- PMCMRplus::frdAllPairsNemenyiTest(y = d[[metric]], groups = factor(d$Design), blocks = factor(d$Block))
  k <- ncol(wide); n <- nrow(wide); des <- colnames(wide)
  pm <- matrix(NA_real_, k, k, dimnames = list(des, des)); pr <- nt$p.value
  for (r in rownames(pr)) for (c in colnames(pr)) if (!is.na(pr[r, c])) pm[r, c] <- pm[c, r] <- pr[r, c]
  cd <- stats::qtukey(0.95, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * n))
  pairs <- utils::combn(des, 2)
  w <- data.frame(Design_A = pairs[1, ], Design_B = pairs[2, ], stringsAsFactors = FALSE)
  w$Median_Diff <- apply(pairs, 2, function(p) stats::median(wide[, p[1]] - wide[, p[2]]))
  w$A_Better_Blocks <- apply(pairs, 2, function(p) sum(wide[, p[1]] < wide[, p[2]]))
  w$P_Raw <- apply(pairs, 2, function(p) stats::wilcox.test(wide[, p[1]], wide[, p[2]], paired = TRUE, exact = FALSE)$p.value)
  w$P_Holm <- stats::p.adjust(w$P_Raw, "holm")
  w$Nemenyi_P <- apply(pairs, 2, function(p) pm[p[1], p[2]])
  list(friedman = c(statistic = unname(ft$statistic), df = unname(ft$parameter), p = ft$p.value),
       avg_ranks = avg, nemenyi_p = pm, cd = cd, wilcoxon = w, n_blocks = n, k = k)
}

#' Critical-difference diagram using code/10's plot_cd_diagram() unchanged
#'
#' The three helper functions are evaluated from code/10's source without
#' sourcing the whole module (which would load the grid visualization stack).
#' @param tests Output of me_design_tests().
#' @param title Plot title.
#' @return ggplot object.
#' @family methods_ext_stored
#' @examples
#' # me_cd_plot(tests, "Both-arms spillover")
me_cd_plot <- function(tests, title) {
  if (!exists("plot_cd_diagram")) {
    ex <- parse(file.path(.rs_project, "code", "10_statistical_comparisons.R"))
    for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                      as.character(e[[2]]) %in% c("plot_cd_diagram", "find_cd_groups", "short_design_label"))
      eval(e, envir = globalenv())
  }
  plot_cd_diagram(list(avg_ranks = tests$avg_ranks, critical_diff = tests$cd,
                       p_matrix = tests$nemenyi_p, alpha = 0.05), title = title)
}

# Production fit cache access ----
#' Load the primary queen mean-rank production rows (education and beta = 1)
#' @return performance.csv rows with Beta and Design_Name added.
#' @family methods_ext_stored
#' @examples
#' # p <- ms_perf()
ms_perf <- function() {
  p <- read.csv(file.path(.me_prod_root, "production", "performance.csv"), stringsAsFactors = FALSE)
  p <- p[p$Neighbor == "queen" & p$Summary == "mean_rank", ]
  p$Beta <- ifelse(p$Model == "education", 0, 1)
  p$Design_Name <- me_design_name(p$Design_ID)
  p$Year <- as.character(p$Year)
  if (nrow(p) != 864) stop("Expected 864 queen mean-rank rows, got ", nrow(p))
  p
}

#' Frequency-weighted per-fit statistic over one production source
#'
#' Each unique allocation's fits are weighted by its draw frequency, as in
#' production's rs_distribution(). Boundary/error fits are excluded and counted.
#' @param b Production checkpoint (rs_block output).
#' @param fun Function(fits data frame) returning a named numeric vector of
#'   per-allocation means.
#' @return Named vector: weighted means, plus MC SE over allocations and counts.
#' @family methods_ext_stored
#' @examples
#' # ms_weighted(b, function(f) c(Reject = mean(abs(f$tau / f$se) > 1.96)))
ms_weighted <- function(b, fun) {
  a <- b$allocations
  vals <- t(vapply(a$Allocation_ID, function(h) {
    f <- b$fits[[h]]
    bad <- grepl("interval bound|ERROR:", f$messages) | !is.finite(f$error) | !is.finite(f$se)
    c(fun(f[!bad, , drop = FALSE]), N_Excluded = sum(bad), N_Fits = nrow(f))
  }, numeric(length(fun(b$fits[[1]])) + 2)))
  w <- a$Frequency / sum(a$Frequency); J <- sum(a$Frequency)
  m <- colSums(w * vals[, seq_len(ncol(vals) - 2), drop = FALSE])
  se <- if (nrow(vals) > 1) sqrt(colSums(w * sweep(vals[, seq_len(ncol(vals) - 2), drop = FALSE], 2, m)^2) / (J - 1) / J +
                                 colSums(w^2 * vals[, seq_len(ncol(vals) - 2), drop = FALSE] * (1 - vals[, seq_len(ncol(vals) - 2), drop = FALSE]) / b$R))
        else sqrt(m * (1 - m) / b$R)
  c(m, stats::setNames(se, paste0("MCSE_", names(m))), N_Excluded = sum(vals[, "N_Excluded"] * a$Frequency),
    N_Fits_Unique = sum(vals[, "N_Fits"]))
}

# A0 and A5 Stage 1 ----
#' A0 power and A5 Stage 1 interval corrections from the stored fits
#'
#' Per fit: t = tau-hat / SE. Power: |t| > z_{0.975} (two-sided Wald, as A3), and
#' with a t_{n-p} critical value. Coverage corrections of the Wald interval
#' tau-hat +/- c x SE x f: (a) c = t_{n-p}; (b) f = sqrt(n / (n - p)) df-corrected
#' residual variance; (c) both; with n = 58 and p = mean coefficients + rho
#' (4 for the education model, 5 with X), plus (b') p = mean coefficients only.
#' MC SEs here are approximate (allocation variance / J plus binomial noise).
#' @param p Production rows from ms_perf().
#' @return Data frame with one row per production row.
#' @family methods_ext_stored
#' @examples
#' # ms_a0_a5(ms_perf())
ms_a0_a5 <- function(p) {
  cp <- file.path(.me_prod_root, "production", "checkpoints")
  out <- lapply(seq_len(nrow(p)), function(i) {
    b <- readRDS(file.path(cp, paste0(p$Source_ID[i], ".rds")))
    k <- if (p$Beta[i] == 1) 4 else 3; n <- 58; z <- stats::qnorm(0.975)
    tc <- stats::qt(0.975, n - k - 1)
    fb <- sqrt(n / (n - k - 1)); fb2 <- sqrt(n / (n - k))
    r <- ms_weighted(b, function(f) {
      t <- abs(f$tau / f$se); e <- abs(f$error)
      c(Power_Wald = mean(t > z), Power_t = mean(t > tc),
        Cov_Wald = mean(e <= z * f$se), Cov_t = mean(e <= tc * f$se),
        Cov_dfvar = mean(e <= z * fb * f$se), Cov_dfvar_both = mean(e <= tc * fb * f$se),
        Cov_dfvar_meanonly = mean(e <= z * fb2 * f$se))
    })
    cbind(p[i, c("Year", "Beta", "Rho", "Gamma", "Regime", "Design_ID", "Design_Name", "Coverage", "MCSE_Coverage")],
          as.data.frame(t(r)), t_crit = tc, df_factor = fb)
  })
  do.call(rbind, out)
}

# A1 burden share at N = 58 ----
#' Share of each year's SUD deaths and person-years in treated clusters
#'
#' Death share (plan A1) counts deaths, so it is driven by population; the rate
#' share (sum of treated colleges' rates / sum of all rates) and the treated vs
#' control mean rate measure targeting of high-incidence colleges per college, and
#' Enrichment = death share / person-year share measures deaths reached per
#' person-year treated.
#'
#' Production allocations: for a given year and design, the support key excludes
#' rho, gamma and regime, so every row of that year/design has the same draws
#' (checked); the first J = 100 draws are used.
#' @param p Production rows.
#' @param setup Frozen setup.
#' @return list(alloc = per-allocation shares, summary = mean, p10, p90 by year x design).
#' @family methods_ext_stored
#' @examples
#' # ms_a1(ms_perf(), me_load_setup())
ms_a1 <- function(p, setup) {
  cp <- file.path(.me_prod_root, "production", "checkpoints")
  pe <- p[p$Beta == 0, ]
  rows <- list()
  for (y in as.character(2018:2021)) for (id in 1:9) {
    src <- unique(pe$Source_ID[pe$Year == y & pe$Design_ID == id])
    Zs <- lapply(src, function(s) readRDS(file.path(cp, paste0(s, ".rds")))$Z[, 1:100, drop = FALSE])
    if (!all(vapply(Zs, identical, logical(1), Zs[[1]]))) stop("Allocation draws differ across settings: ", y, " design ", id)
    Z <- Zs[[1]]; a <- setup$annual[[y]]
    rows[[length(rows) + 1L]] <- data.frame(Year = y, Design_ID = id, Design_Name = me_design_name(id),
      Allocation = 1:100, Treated = colSums(Z),
      Death_Share = colSums(a$deaths * Z) / sum(a$deaths),
      PY_Share = colSums(a$person_years * Z) / sum(a$person_years),
      Rate_Share = colSums(a$rate_per_100k * Z) / sum(a$rate_per_100k),
      Treated_Rate = colSums(a$rate_per_100k * Z) / colSums(Z),
      Control_Rate = colSums(a$rate_per_100k * (1 - Z)) / colSums(1 - Z))
  }
  al <- do.call(rbind, rows)
  q <- function(x, pr) unname(stats::quantile(x, pr, type = 7))
  s <- do.call(rbind, lapply(split(al, list(al$Year, al$Design_ID)), function(d) data.frame(
    Year = d$Year[1], Design_ID = d$Design_ID[1], Design_Name = d$Design_Name[1],
    Mean_Treated_Fraction = mean(d$Treated) / 58,
    Death_Share_Mean = mean(d$Death_Share), Death_Share_P10 = q(d$Death_Share, 0.1), Death_Share_P90 = q(d$Death_Share, 0.9),
    PY_Share_Mean = mean(d$PY_Share), PY_Share_P10 = q(d$PY_Share, 0.1), PY_Share_P90 = q(d$PY_Share, 0.9),
    Death_minus_PY = mean(d$Death_Share - d$PY_Share),
    Enrichment = mean(d$Death_Share / d$PY_Share),
    Rate_Share_Mean = mean(d$Rate_Share), Rate_Share_P10 = q(d$Rate_Share, 0.1), Rate_Share_P90 = q(d$Rate_Share, 0.9),
    Treated_Rate_Mean = mean(d$Treated_Rate), Control_Rate_Mean = mean(d$Control_Rate))))
  rownames(s) <- NULL
  list(alloc = al, summary = s[order(s$Year, s$Design_ID), ])
}

# Driver ----
#' Run A0, A1, A2 and A5 Stage 1 and write all outputs
#' @return Invisibly, a list of the result tables.
#' @family methods_ext_stored
#' @examples
#' # ms_run()
ms_run <- function() {
  out <- file.path(.me_root, "stored"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  me_manifest(out, list(analysis = "stored A0/A1/A2/A5-stage1", production = "real_sud_rev_20261002"),
              file.path(.ms_dir, c("methods_ext_core.R", "methods_ext_stored.R")))
  setup <- me_load_setup(); p <- ms_perf()
  log <- character(0); say <- function(...) { l <- sprintf(...); cat(l, "\n"); log <<- c(log, l) }

  ## A0 + A5 Stage 1 ----
  r <- ms_a0_a5(p)
  say("A0/A5: %d production rows; excluded (boundary/error) fits: %d", nrow(r), sum(r$N_Excluded))
  say("A0/A5: max |stored coverage - recomputed Wald coverage| = %.2e", max(abs(r$Coverage - r$Cov_Wald)))
  write.csv(r[, c("Year", "Beta", "Rho", "Gamma", "Regime", "Design_ID", "Design_Name", "Power_Wald",
                  "MCSE_Power_Wald", "Power_t", "N_Excluded", "N_Fits_Unique")],
            file.path(out, "a0_stored_power.csv"), row.names = FALSE)
  write.csv(r[, setdiff(names(r), c("Power_Wald", "MCSE_Power_Wald", "Power_t", "MCSE_Power_t"))],
            file.path(out, "a5_stage1_coverage.csv"), row.names = FALSE)
  say("A0: tau = 1 Wald power, min %.4f, median %.4f (beta = 0); min %.4f (beta = 1)",
      min(r$Power_Wald[r$Beta == 0]), median(r$Power_Wald[r$Beta == 0]), min(r$Power_Wald[r$Beta == 1]))
  cov_cols <- c("Cov_Wald", "Cov_t", "Cov_dfvar", "Cov_dfvar_both", "Cov_dfvar_meanonly")
  dr <- aggregate(r[cov_cols], r[c("Beta", "Regime", "Design_Name")], mean)
  write.csv(dr, file.path(out, "a5_stage1_coverage_by_design_regime.csv"), row.names = FALSE)
  for (cc in cov_cols) {
    v <- dr[[cc]]; vs <- r[[cc]]
    say("A5 %-19s design x regime mean %.4f-%.4f; per setting %.4f-%.4f; rule (94-96%% every design x regime): %s",
        cc, min(v), max(v), min(vs), max(vs), if (all(v >= 0.94 & v <= 0.96)) "MET" else "NOT MET")
  }
  say("A5 Stage 1 grid: not checked; the grid caches (results/checkpoints/rev_2026-09) keep surface-level summaries, not per-fit tau-hat/SE.")

  ## A1 ----
  a1 <- ms_a1(p, setup)
  write.csv(a1$summary, file.path(out, "a1_burden_share.csv"), row.names = FALSE)
  saveRDS(a1$alloc, file.path(out, "a1_burden_share_allocations.rds"))
  am <- aggregate(cbind(Death_Share_Mean, PY_Share_Mean, Mean_Treated_Fraction, Enrichment, Rate_Share_Mean,
                        Treated_Rate_Mean, Control_Rate_Mean) ~ Design_Name, a1$summary, mean)
  for (i in seq_len(nrow(am))) say("A1 %-28s deaths %.3f, person-years %.3f, enrichment %.2f, rate share %.3f, treated/control rate %.1f/%.1f, treated %.3f (mean of years)",
                                   am$Design_Name[i], am$Death_Share_Mean[i], am$PY_Share_Mean[i], am$Enrichment[i],
                                   am$Rate_Share_Mean[i], am$Treated_Rate_Mean[i], am$Control_Rate_Mean[i], am$Mean_Treated_Fraction[i])
  s <- a1$summary
  s$Design_Name <- factor(s$Design_Name, levels = am$Design_Name[order(am$Death_Share_Mean)])
  fig <- ggplot(s, aes(x = Death_Share_Mean, y = Design_Name, colour = Year)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = Death_Share_P10, xmax = Death_Share_P90), width = 0.2, orientation = "y",
                   position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6)) +
    labs(x = "Share of the year's SUD deaths in treated clusters (mean; bars 10th-90th percentile of allocations)",
         y = NULL, colour = "Year") + theme_minimal(base_size = 11)
  ggsave(file.path(out, "a1_burden_share.pdf"), fig, width = 8, height = 4.5)
  ggsave(file.path(out, "a1_burden_share.png"), fig, width = 8, height = 4.5, dpi = 150)

  ## A2 ----
  tests_all <- list(); pdf(file.path(out, "a2_cd_plots.pdf"), width = 9, height = 4)
  for (b in c(0, 1)) for (g in c("both", "control_only")) {
    d <- p[p$Beta == b & p$Regime == g, ]; d$Design <- d$Design_Name
    tt <- me_design_tests(d, c("Year", "Rho", "Gamma"))
    lab <- sprintf("beta%d_%s", b, g); tests_all[[lab]] <- tt
    say("A2 %s: %d blocks, Friedman chi2 = %.1f (df %d), p = %.2g; best avg rank %s (%.2f); CD = %.2f",
        lab, tt$n_blocks, tt$friedman[["statistic"]], tt$friedman[["df"]], tt$friedman[["p"]],
        names(tt$avg_ranks)[1], tt$avg_ranks[1], tt$cd)
    print(me_cd_plot(tt, sprintf("NC consistency summary: beta = %d, %s spillover (24 blocks)", b, g)))
  }
  dev.off()
  write.csv(do.call(rbind, lapply(names(tests_all), function(l) data.frame(Setting = l,
    Statistic = tests_all[[l]]$friedman[["statistic"]], DF = tests_all[[l]]$friedman[["df"]],
    P = tests_all[[l]]$friedman[["p"]], N_Blocks = tests_all[[l]]$n_blocks, CD = tests_all[[l]]$cd))),
    file.path(out, "a2_tests_friedman.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(names(tests_all), function(l) data.frame(Setting = l,
    Design = names(tests_all[[l]]$avg_ranks), Avg_Rank = unname(tests_all[[l]]$avg_ranks)))),
    file.path(out, "a2_tests_avg_ranks.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(names(tests_all), function(l) cbind(Setting = l, tests_all[[l]]$wilcoxon))),
    file.path(out, "a2_tests_pairwise.csv"), row.names = FALSE)
  saveRDS(tests_all, file.path(out, "a2_tests.rds"))
  writeLines(c(sprintf("Run %s", format(Sys.time())),
               "A2 p-values are a consistency summary (one map, four correlated years, shared allocation streams), not trial-level evidence.",
               log), file.path(out, "stored_summary.txt"))
  invisible(list(a0a5 = r, a1 = a1, a2 = tests_all))
}

if (!exists("methods_ext_define_only") || !isTRUE(methods_ext_define_only)) ms_run()
