# ============================================================
# Script: methods_ext_sensitivity_summary.R
# Purpose: Summaries of A5 Stage 2 (coverage methods and adoption rule), A6
#   (CAR, reverse misspecification, Bayesian CAR; check V10) and A7 (spillover
#   forms) from the per-fit outputs of methods_ext_sensitivity.R.
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: methods_ext_stored.R (names, production rows)
# ============================================================
# Run from the project root:  Rscript application/code/methods_ext_sensitivity_summary.R
# Writes application/results/methods_ext_2026-10/sensitivity_summaries/.
# MC SE of a setting metric: sd of allocation-level means / sqrt(J) for random
# designs; binomial / iid fit-level SE for single (deterministic) allocations.

.mz_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(p)) return(normalizePath(dirname(p)))
  }
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else stop("Run or source by path")
})
methods_ext_define_only <- TRUE
source(file.path(.mz_dir, "methods_ext_stored.R"))

#' Setting-level mean and MC SE of a per-fit quantity, allocation as top unit
#' @param v Per-fit values (NA = failed fit, excluded and counted by the caller).
#' @param alloc Allocation index per fit.
#' @return c(mean, mcse).
#' @family methods_ext_sensitivity_summary
#' @examples
#' mz_mean_se(c(1, 0, 1, 1), c(1, 1, 2, 2))
mz_mean_se <- function(v, alloc) {
  ok <- is.finite(v); v <- v[ok]; a <- alloc[ok]
  if (!length(v)) return(c(NA_real_, NA_real_))
  m <- tapply(v, a, mean)
  c(mean(v), if (length(m) > 1) stats::sd(m) / sqrt(length(m)) else stats::sd(v) / sqrt(length(v)))
}

#' Summarize per-fit estimates of one estimator by setting
#' @param f Fits. @param est Column prefix ("SAR", "CAR", "BCAR"). @param by Setting columns.
#' @param lo,hi Optional interval columns (default Wald from <est>_se).
#' @return Data frame: MSE, Bias, Coverage with MC SEs and failure counts.
#' @family methods_ext_sensitivity_summary
#' @examples
#' # mz_setting(fits, "SAR", c("Year", "Rho", "Regime", "Design_ID"))
mz_setting <- function(f, est, by, lo = NULL, hi = NULL) {
  tau <- f[[paste0(est, "_tau")]]
  if (is.null(lo)) { se <- f[[paste0(est, "_se")]]; lo_v <- tau - stats::qnorm(0.975) * se; hi_v <- tau + stats::qnorm(0.975) * se }
  else { lo_v <- f[[lo]]; hi_v <- f[[hi]] }
  f$.e <- tau - 1; f$.e2 <- (tau - 1)^2; f$.c <- as.numeric(lo_v <= 1 & hi_v >= 1)
  do.call(rbind, lapply(split(f, f[by], drop = TRUE), function(d) {
    m <- mz_mean_se(d$.e2, d$Allocation); b <- mz_mean_se(d$.e, d$Allocation); c <- mz_mean_se(d$.c, d$Allocation)
    cbind(d[1, by, drop = FALSE], Estimator = est, N_Fits = nrow(d), N_Failed = sum(!is.finite(d$.e)),
          N_CI_Failed = sum(!is.finite(d$.c)), MSE = m[1], MCSE_MSE = m[2], Bias = b[1], MCSE_Bias = b[2],
          Coverage = c[1], MCSE_Coverage = c[2])
  }))
}

#' MSE ratio to SRS within blocks, then mean / range by design x regime
#' @param s Setting summary. @param block Block columns (excluding Design_ID, Regime).
#' @param group Extra grouping columns (e.g. Estimator, DGP, Form).
#' @return Data frame of mean, min and max ratio, with unidentified settings marked.
#' @family methods_ext_sensitivity_summary
#' @examples
#' # mz_ratio(s, c("Year", "Rho", "Gamma"), "Estimator")
mz_ratio <- function(s, block, group) {
  srs <- s[s$Design_ID == 9, intersect(c(group, block, "Regime", "MSE", "Share_Identified"), names(s))]
  m <- merge(s, srs, by = c(group, block, "Regime"), suffixes = c("", "_SRS"))
  m$Ratio <- m$MSE / m$MSE_SRS
  r <- do.call(rbind, lapply(split(m, m[c(group, "Regime", "Design_ID")], drop = TRUE), function(d)
    cbind(d[1, c(group, "Regime", "Design_ID")], Design = me_design_name(d$Design_ID[1]),
          Ratio_Mean = mean(d$Ratio), Ratio_Min = min(d$Ratio), Ratio_Max = max(d$Ratio),
          Blocks = nrow(d), Blocks_Unidentified = sum(!is.finite(d$Ratio)),
          Min_Share_Identified = if (is.null(d$Share_Identified)) 1 else min(d$Share_Identified),
          Min_Share_Identified_SRS = if (is.null(d$Share_Identified_SRS)) 1 else min(d$Share_Identified_SRS))))
  rownames(r) <- NULL; r
}

#' Spearman correlation of the nine-design MSE ordering between two estimators
#' @param s Setting summary with both estimators. @param block Block columns.
#' @param a,b Estimator labels.
#' @return Per-block correlations.
#' @family methods_ext_sensitivity_summary
#' @examples
#' # mz_spearman(s, c("DGP", "Year", "Rho", "Gamma", "Regime"), "SAR", "CAR")
mz_spearman <- function(s, block, a, b) {
  do.call(rbind, lapply(split(s, s[block], drop = TRUE), function(d) {
    x <- d[d$Estimator == a, ]; y <- d[d$Estimator == b, ]; y <- y[match(x$Design_ID, y$Design_ID), ]
    ok <- is.finite(x$MSE) & is.finite(y$MSE)
    cbind(d[1, block, drop = FALSE], Spearman = if (sum(ok) > 2) stats::cor(x$MSE[ok], y$MSE[ok], method = "spearman") else NA_real_,
          Best_A = me_design_name(x$Design_ID[which.min(x$MSE)]), Best_B = me_design_name(y$Design_ID[which.min(y$MSE)]))
  }))
}

#' Run all sensitivity summaries that have outputs
#' @return Invisibly, list of tables.
#' @family methods_ext_sensitivity_summary
#' @examples
#' # mz_run()
mz_run <- function(suffix = "") {
  out <- file.path(.me_root, paste0("sensitivity_summaries", suffix)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- character(0); say <- function(...) { l <- sprintf(...); cat(l, "\n"); log <<- c(log, l) }
  path <- function(a) file.path(.me_root, paste0(a, suffix), "fits.rds")
  res <- list()

  ## A5 Stage 2 ----
  if (file.exists(path("a5"))) {
    f <- readRDS(path("a5")); by <- c("Year", "Rho", "Regime", "Design_ID")
    k <- 3; tc <- stats::qt(0.975, 58 - k - 1); fb <- sqrt(58 / (58 - k - 1))
    f$S1c_lo <- f$SAR_tau - tc * fb * f$SAR_se; f$S1c_hi <- f$SAR_tau + tc * fb * f$SAR_se
    meth <- list(Wald = NULL, Stage1_t_df = c("S1c_lo", "S1c_hi"), Percentile = c("Pct_lo", "Pct_hi"),
                 Bootstrap_t = c("Bt_lo", "Bt_hi"), Profile = c("Profile_lo", "Profile_hi"))
    s <- do.call(rbind, lapply(names(meth), function(mn) { x <- mz_setting(f, "SAR", by, meth[[mn]][1], meth[[mn]][2]); x$Method <- mn; x }))
    write.csv(s, file.path(out, "a5_stage2_coverage_settings.csv"), row.names = FALSE)
    dr <- aggregate(Coverage ~ Method + Regime + Design_ID, s, mean); dr$Design <- me_design_name(dr$Design_ID)
    write.csv(dr, file.path(out, "a5_stage2_coverage_by_design_regime.csv"), row.names = FALSE)
    say("A5 Stage 2: %d fits; failed SAR fits %d; profile intervals missing %d; bootstrap refits failed %d; max |profile LL at tau-hat - LL| %.1e",
        nrow(f), sum(!is.finite(f$SAR_tau)), sum(!is.finite(f$Profile_lo) | !is.finite(f$Profile_hi)), sum(f$Boot_fail, na.rm = TRUE),
        max(f$Profile_LL_check, na.rm = TRUE))
    for (mn in names(meth)) { d <- dr[dr$Method == mn, ]; ds <- s[s$Method == mn, ]
      say("A5 %-12s design x regime %.4f-%.4f; settings %.4f-%.4f; adoption rule (94-96%% every design x regime; rankings unchanged since tau-hat is unchanged): %s",
          mn, min(d$Coverage), max(d$Coverage), min(ds$Coverage), max(ds$Coverage), if (all(d$Coverage >= 0.94 & d$Coverage <= 0.96)) "MET" else "NOT MET") }
    res$a5 <- s
  } else say("A5 Stage 2: no output")

  ## A6 CAR and reverse misspecification ----
  if (file.exists(path("a6car"))) {
    f <- readRDS(path("a6car")); by <- c("DGP", "Year", "Rho", "Gamma", "Regime", "Design_ID")
    s <- rbind(mz_setting(f, "SAR", by), mz_setting(f, "CAR", by))
    write.csv(s, file.path(out, "a6_car_settings.csv"), row.names = FALSE)
    say("A6 CAR: %d fits; SAR failed %d (boundary %d); CAR failed %d (boundary %d, errors %d, warnings %d); unidentified fits %d",
        nrow(f), sum(!is.finite(f$SAR_tau)), sum(f$SAR_bound), sum(!is.finite(f$CAR_tau)), sum(f$CAR_bound), sum(f$CAR_error), sum(f$CAR_warn), sum(!f$Identified))
    rt <- mz_ratio(s, c("Year", "Rho", "Gamma"), c("DGP", "Estimator"))
    write.csv(rt, file.path(out, "a6_car_ratio_to_srs.csv"), row.names = FALSE)
    rtr <- mz_ratio(s, c("Year", "Gamma"), c("DGP", "Estimator", "Rho"))
    write.csv(rtr, file.path(out, "a6_car_ratio_to_srs_by_rho.csv"), row.names = FALSE)
    sp <- mz_spearman(s, c("DGP", "Year", "Rho", "Gamma", "Regime"), "SAR", "CAR")
    write.csv(sp, file.path(out, "a6_car_spearman.csv"), row.names = FALSE)
    for (dg in unique(sp$DGP)) for (rg in unique(sp$Regime)) { d <- sp[sp$DGP == dg & sp$Regime == rg, ]
      say("A6 %s DGP, %s: Spearman(SAR, CAR design MSE order) mean %.2f (range %.2f-%.2f) over %d blocks; best design agrees in %d/%d",
          dg, rg, mean(d$Spearman), min(d$Spearman), max(d$Spearman), nrow(d), sum(d$Best_A == d$Best_B), nrow(d)) }
    cv <- aggregate(cbind(Coverage, Bias) ~ DGP + Estimator + Rho, s, mean)
    for (i in seq_len(nrow(cv))) say("A6 %s DGP, %s fit, rho/lambda = %.1f: mean coverage %.3f, mean bias %+.3f", cv$DGP[i], cv$Estimator[i], cv$Rho[i], cv$Coverage[i], cv$Bias[i])
    z <- f[f$DGP == "sar" & f$Rho == 0 & is.finite(f$SAR_tau) & is.finite(f$CAR_tau), ]
    say("V10 CAR vs SAR at rho = 0 (same data): mean |diff| %.4f, median %.4f, correlation %.4f, mean diff %+.4f over %d fits",
        mean(abs(z$CAR_tau - z$SAR_tau)), stats::median(abs(z$CAR_tau - z$SAR_tau)), stats::cor(z$CAR_tau, z$SAR_tau), mean(z$CAR_tau - z$SAR_tau), nrow(z))
    res$a6 <- s
  } else say("A6 CAR: no output")

  ## A6 Bayesian CAR ----
  if (file.exists(path("a6bayes"))) {
    f <- readRDS(path("a6bayes")); by <- c("Year", "Rho", "Regime", "Design_ID")
    cfgb <- readRDS(file.path(.me_root, paste0("a6bayes", suffix), "manifest.rds"))$config
    f$BCAR_converged <- is.finite(f$BCAR_ess) & f$BCAR_ess >= cfgb$ess_min & abs(f$BCAR_geweke) < 1.96
    s <- rbind(mz_setting(f, "SAR", by), mz_setting(f, "BCAR", by, "BCAR_lo", "BCAR_hi"))
    write.csv(s, file.path(out, "a6_bayes_settings.csv"), row.names = FALSE)
    say("A6 Bayes: %d MCMC fits (burn-in %d, %d iterations, thin %d); errors %d; ESS < %d: %d; |Geweke z| >= 1.96: %d; not converged (either) %d (%.1f%%); ESS median %.0f (min %.0f)",
        nrow(f), cfgb$burnin, cfgb$n_sample, cfgb$thin, sum(f$BCAR_error), cfgb$ess_min, sum(f$BCAR_ess < cfgb$ess_min, na.rm = TRUE),
        sum(abs(f$BCAR_geweke) >= 1.96, na.rm = TRUE), sum(!f$BCAR_converged), 100 * mean(!f$BCAR_converged),
        stats::median(f$BCAR_ess, na.rm = TRUE), min(f$BCAR_ess, na.rm = TRUE))
    say("A6 Bayes: |Geweke| >= 1.96 is expected in about 5%% of converged chains by chance")
    rt <- mz_ratio(s, c("Year", "Rho"), "Estimator"); write.csv(rt, file.path(out, "a6_bayes_ratio_to_srs.csv"), row.names = FALSE)
    sp <- mz_spearman(s, c("Year", "Rho", "Regime"), "SAR", "BCAR"); write.csv(sp, file.path(out, "a6_bayes_spearman.csv"), row.names = FALSE)
    for (rg in unique(sp$Regime)) { d <- sp[sp$Regime == rg, ]
      say("A6 Bayes %s: Spearman(SAR ML, Bayesian CAR design MSE order) mean %.2f (range %.2f-%.2f); best design agrees %d/%d",
          rg, mean(d$Spearman), min(d$Spearman), max(d$Spearman), sum(d$Best_A == d$Best_B), nrow(d)) }
    cv <- aggregate(cbind(Coverage, Bias, MSE) ~ Estimator + Rho, s, mean)
    for (i in seq_len(nrow(cv))) say("A6 Bayes %s rho = %.1f: coverage %.3f, bias %+.3f, MSE %.3f", cv$Estimator[i], cv$Rho[i], cv$Coverage[i], cv$Bias[i], cv$MSE[i])
    z <- f[is.finite(f$BCAR_tau) & is.finite(f$SAR_tau), ]
    say("V10 Bayes vs ML: correlation %.3f, mean |diff| %.3f, mean diff %+.3f over %d fits", stats::cor(z$BCAR_tau, z$SAR_tau),
        mean(abs(z$BCAR_tau - z$SAR_tau)), mean(z$BCAR_tau - z$SAR_tau), nrow(z))
    me_seed("v10-hand"); h <- z[sample.int(nrow(z), 5), c("Year", "Rho", "Regime", "Design_ID", "SAR_tau", "SAR_se", "BCAR_tau", "BCAR_lo", "BCAR_hi", "BCAR_ess")]
    write.csv(h, file.path(out, "v10_hand_checked_fits.csv"), row.names = FALSE)
    for (i in seq_len(nrow(h))) say("V10 hand check: %s rho=%.1f %s %s: ML %.3f (SE %.3f) vs posterior mean %.3f [%.3f, %.3f], ESS %.0f",
      h$Year[i], h$Rho[i], h$Regime[i], me_design_name(h$Design_ID[i]), h$SAR_tau[i], h$SAR_se[i], h$BCAR_tau[i], h$BCAR_lo[i], h$BCAR_hi[i], h$BCAR_ess[i])
    res$bayes <- s
  } else say("A6 Bayes: no output")

  ## A7 spillover forms ----
  if (file.exists(path("a7"))) {
    f <- readRDS(path("a7")); by <- c("Form", "Year", "Rho", "Gamma", "Regime", "Design_ID")
    s <- mz_setting(f, "SAR", by)
    idt <- aggregate(cbind(Identified, Alias = SAR_alias) ~ Form + Regime + Design_ID, f, mean)
    s <- merge(s, aggregate(cbind(Share_Identified = Identified) ~ Form + Year + Rho + Gamma + Regime + Design_ID, f, mean))
    # MSE is over identified fits only; Share_Identified travels with every ratio, and
    # designs never identified (MSE NA) are reported as not identified, not dropped.
    write.csv(s, file.path(out, "a7_forms_settings.csv"), row.names = FALSE)
    write.csv(idt, file.path(out, "a7_identification.csv"), row.names = FALSE)
    un <- idt[idt$Identified < 1, ]
    for (i in seq_len(nrow(un))) say("A7 %s %s %s: identified in %.1f%% of fits (aliased Spill %.1f%%) - reported as not identified",
                                     un$Form[i], un$Regime[i], me_design_name(un$Design_ID[i]), 100 * un$Identified[i], 100 * un$Alias[i])
    rt <- mz_ratio(s, c("Year", "Rho", "Gamma"), "Form"); write.csv(rt, file.path(out, "a7_ratio_to_srs.csv"), row.names = FALSE)
    for (fm in unique(rt$Form)) for (rg in unique(rt$Regime)) { d <- rt[rt$Form == fm & rt$Regime == rg & is.finite(rt$Ratio_Mean), ]; d <- d[order(d$Ratio_Mean), ]
      say("A7 %-12s %-12s best-to-worst (mean MSE ratio to SRS; identified share if < 1): %s", fm, rg,
          paste(sprintf("%s %.2f%s", d$Design, d$Ratio_Mean, ifelse(d$Min_Share_Identified < 1, sprintf(" [id %.2f]", d$Min_Share_Identified), "")), collapse = ", "))
      if (any(d$Min_Share_Identified_SRS < 1)) say("A7 %s %s: SRS reference itself identified in only %.1f%% of fits in some block", fm, rg, 100 * min(d$Min_Share_Identified_SRS)) }
    pr <- ms_perf(); pr <- pr[pr$Beta == 0 & pr$Rho %in% c(0, 0.5), ]
    cmp <- merge(s[s$Form == "proportional", ], pr[c("Year", "Rho", "Gamma", "Regime", "Design_ID", "Mean_MSE", "SE_Mean_MSE_Joint")])
    cmp$Z <- (cmp$MSE - cmp$Mean_MSE) / sqrt(cmp$MCSE_MSE^2 + cmp$SE_Mean_MSE_Joint^2)
    say("A7 proportional vs production MSE (V11 at N = 58): %d settings, |z| > 2 in %d, max |z| %.2f", nrow(cmp), sum(abs(cmp$Z) > 2), max(abs(cmp$Z)))
    res$a7 <- s
  } else say("A7: no output")
  writeLines(c(sprintf("Run %s", format(Sys.time())), log), file.path(out, "sensitivity_summary_log.txt"))
  invisible(res)
}

if (sys.nframe() == 0L) { args <- commandArgs(TRUE); mz_run(if (length(args)) args[1] else "") }
