# ============================================================
# Script: methods_ext_summaries.R
# Purpose: Summaries of the A3/A4 participation and power runs: setting metrics,
#   N80 and MDE with bootstrap intervals, m x ICC translation tables, the
#   guidance-holds table with A2 tests per N, A1 burden share per N, A4
#   comparison, and verification checks V1, V2 and V8.
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: methods_ext_core.R, methods_ext_stored.R (test helpers), parallel, ggplot2
# ============================================================
# Run from the project root after the a3 and a4 profiles finish:
#   Rscript application/code/methods_ext_summaries.R [workers]
# Reads participation/{a3,a4}/checkpoints (gitignored); writes
# application/results/methods_ext_2026-10/summaries/ (small CSVs, figures, text).
#
# Monte Carlo units: the participant subset is the top-level unit for N < 58
# (MC SE = sd of subset-level metrics / sqrt(50)); at N = 58 the allocation is
# (sd / sqrt(J)); a single deterministic allocation uses the fit-level SE.

.mu_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(p)) return(normalizePath(dirname(p)))
  }
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else stop("Run or source by path")
})
methods_ext_define_only <- TRUE
source(file.path(.mu_dir, "methods_ext_stored.R"))   # sources core; me_design_tests, me_design_name

# Load and reduce checkpoints ----
#' Reduce one task checkpoint to setting rows and unit-level counts
#'
#' Ratios per setting: Power = n_rej / n_ci, Power_Adj = n_rej_adj / n_ci,
#' Coverage = n_cov / n_ci, Bias = sum_err / n_est, MSE = sum_err2 / n_est.
#' Fit-level SE of MSE for a single allocation uses the normal approximation
#' Var(e^2) = 2 s^4 + 4 b^2 s^2 (b = bias, s^2 = error variance).
#' @param f Checkpoint path.
#' @return list(rows = setting data frame, units = list of unit count matrices,
#'   alloc = allocation diagnostics, subsets = subset diagnostics).
#' @family methods_ext_summaries
#' @examples
#' # mu_reduce(f)
mu_reduce <- function(f) {
  x <- readRDS(f); st <- x$stats; tk <- x$task
  S <- dim(st)[2]; J <- dim(st)[3]
  unit_dim <- if (S > 1) 2 else if (J > 1) 3 else NA
  tot <- apply(st, c(1, 4), sum)
  units <- lapply(seq_len(nrow(x$cells)), function(ci)
    if (is.na(unit_dim)) NULL else apply(st[ci, , , , drop = FALSE], c(unit_dim, 4), sum))
  se_ratio <- function(u, num, den) { if (is.null(u)) return(NA_real_); v <- u[, num] / u[, den]; v <- v[is.finite(v)]
    if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_ }
  rows <- do.call(rbind, lapply(seq_len(nrow(x$cells)), function(ci) {
    t <- tot[ci, ]; u <- units[[ci]]
    p <- t[["n_rej"]] / t[["n_ci"]]; pa <- t[["n_rej_adj"]] / t[["n_ci"]]; cv <- t[["n_cov"]] / t[["n_ci"]]
    b <- t[["sum_err"]] / t[["n_est"]]; mse <- t[["sum_err2"]] / t[["n_est"]]; s2 <- max(0, mse - b^2)
    fit_se <- function(q) sqrt(q * (1 - q) / t[["n_ci"]])
    data.frame(Mechanism = tk$Mechanism, Year = tk$Year, N = tk$N, Design_ID = tk$Design_ID,
      x$cells[ci, ], N_Subsets = S, N_Alloc = J, N_Fit = t[["n_fit"]], N_Est = t[["n_est"]], N_CI = t[["n_ci"]],
      Power = p, MCSE_Power = if (is.na(unit_dim)) fit_se(p) else se_ratio(u, "n_rej", "n_ci"),
      Power_Adj = pa, MCSE_Power_Adj = if (is.na(unit_dim)) fit_se(pa) else se_ratio(u, "n_rej_adj", "n_ci"),
      Coverage = cv, MCSE_Coverage = if (is.na(unit_dim)) fit_se(cv) else se_ratio(u, "n_cov", "n_ci"),
      Bias = b, MCSE_Bias = if (is.na(unit_dim)) sqrt(s2 / t[["n_est"]]) else se_ratio(u, "sum_err", "n_est"),
      MSE = mse, MCSE_MSE = if (is.na(unit_dim)) sqrt((2 * s2^2 + 4 * b^2 * s2) / t[["n_est"]]) else se_ratio(u, "sum_err2", "n_est"),
      N_Alias = t[["n_alias"]], N_Bound = t[["n_bound"]], N_Error = t[["n_error"]], N_Warn = t[["n_warn"]],
      N_Unidentified = t[["n_unident"]], stringsAsFactors = FALSE)
  }))
  crit <- vapply(x$t0, `[[`, numeric(1), "crit"); names(crit) <- vapply(x$t0, function(z) as.character(z$cell), "")
  rows$Crit_Adj <- NA_real_
  for (ci in seq_len(nrow(rows))) {
    c0 <- which(x$cells$Tau == 0 & x$cells$Beta == x$cells$Beta[ci] & x$cells$Rho == x$cells$Rho[ci] &
                x$cells$Regime == x$cells$Regime[ci] & x$cells$Arm == x$cells$Arm[ci])
    rows$Crit_Adj[ci] <- crit[as.character(c0)]
  }
  sb <- x$subsets
  rows$Islands_Mean <- mean(sb$Islands); rows$Islands_Max <- max(sb$Islands)
  rows$Subsets_With_Islands <- sum(sb$Islands > 0); rows$Guard_Bypass <- sum(sb$Guard_Bypass)
  rows$Engine_Errors <- sum(!is.na(sb$Engine_Error)); rows$Min_Region_Size <- min(sb$Min_Region_Size)
  rows$Subsets_Region_Lt2 <- sum(sb$Min_Region_Size < 2); rows$Subsets_Region_Lt3 <- sum(sb$Regions_Below3 > 0)
  rows$Mean_Treated <- mean(x$allocations$Treated); rows$Death_Share <- mean(x$allocations$Death_Share)
  rows$Islands_Treated <- mean(x$allocations$Islands_Treated)
  list(rows = rows, units = units, task = tk$Task, cells = x$cells,
       alloc = cbind(tk[c("Mechanism", "Year", "N", "Design_ID")], x$allocations),
       subsets = cbind(tk[c("Mechanism", "Year", "N", "Design_ID")], sb), n_fits = x$n_fits, elapsed = x$elapsed)
}

#' Duplicate the shared tau = 0 rows under each kappa
#' @param r Setting rows (Kappa NA at tau = 0).
#' @param kappas Kappa values.
#' @return Rows with Kappa filled and Shared_Tau0 flag.
#' @family methods_ext_summaries
#' @examples
#' # mu_expand_tau0(rows, c(0.5, 0.8))
mu_expand_tau0 <- function(r, kappas) {
  z <- r[is.na(r$Kappa), ]; nz <- r[!is.na(r$Kappa), ]; nz$Shared_Tau0 <- FALSE
  rbind(nz, do.call(rbind, lapply(kappas, function(k) { z$Kappa <- k; z$Shared_Tau0 <- TRUE; z })))
}

# Interpolation ----
#' First crossing of a target along a curve, by linear interpolation
#'
#' Scans the curve in order of x; returns the x where y first reaches the target,
#' interpolating linearly from the previous point. Left-censored (target reached at
#' the first x) returns -Inf; never reached returns +Inf. No smoothing: a curve that
#' later dips below the target is flagged by V8, not here.
#' @param x Increasing grid (N or effect). @param y Power values.
#' @param target Target power (0.8).
#' @return Crossing x, -Inf or Inf.
#' @family methods_ext_summaries
#' @examples
#' mu_cross(c(15, 20, 25), c(0.5, 0.7, 0.9))
mu_cross <- function(x, y, target = 0.8) {
  i <- which(y >= target)[1]
  if (is.na(i)) return(Inf)
  if (i == 1) return(-Inf)
  x[i - 1] + (target - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
}

#' Format a censored crossing for tables
#' @param v Crossing. @param lo,hi Grid limits. @param digits Rounding.
#' @return Character.
#' @family methods_ext_summaries
#' @examples
#' mu_fmt(Inf, 15, 58)
mu_fmt <- function(v, lo, hi, digits = 0) ifelse(is.infinite(v) & v > 0, paste0(">", hi),
  ifelse(is.infinite(v), paste0("<=", lo), formatC(v, format = "f", digits = digits)))

#' Power at an arbitrary effect by linear interpolation over the effect grid
#' @param tau Effect grid. @param p Power. @param t Target effects.
#' @return Interpolated power (NA beyond the grid's top; flagged by the caller).
#' @family methods_ext_summaries
#' @examples
#' mu_power_at(c(0, 1, 2), c(0.05, 0.5, 0.9), 1.5)
mu_power_at <- function(tau, p, t) ifelse(t > max(tau), NA_real_, stats::approx(tau, p, xout = pmin(t, max(tau)))$y)

# Driver ----
#' Build every A3/A4 summary table and verification check
#' @param workers Forked workers for loading checkpoints.
#' @param B Bootstrap resamples for N80/MDE intervals.
#' @param profiles Named profile folders (a3, a4); pilot can stand in for a3 to test.
#' @param out_name Output folder name under methods_ext_2026-10/.
#' @return Invisibly, the main tables.
#' @family methods_ext_summaries
#' @examples
#' # mu_run(10)
mu_run <- function(workers = 8L, B = 500L, profiles = c(a3 = "a3", a4 = "a4"), out_name = "summaries") {
  out <- file.path(.me_root, out_name); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- character(0); say <- function(...) { l <- sprintf(...); cat(l, "\n"); log <<- c(log, l) }
  me_seed("summaries-bootstrap")
  res <- list()
  for (prof in names(profiles)) {
    cp <- file.path(.me_root, "participation", profiles[[prof]], "checkpoints")
    if (!dir.exists(cp)) { say("%s: no checkpoint folder", prof); next }
    cfg <- readRDS(file.path(.me_root, "participation", profiles[[prof]], "manifest.rds"))$config
    tasks <- mp_tasks_from(cfg)
    have <- file.exists(file.path(cp, paste0(tasks, ".rds")))
    say("%s: %d / %d task checkpoints present", prof, sum(have), length(tasks))
    if (!any(have)) next
    red <- parallel::mclapply(file.path(cp, paste0(tasks[have], ".rds")), mu_reduce, mc.cores = workers)
    bad <- vapply(red, inherits, logical(1), "try-error")
    if (any(bad)) stop("Failed to reduce: ", paste(tasks[have][bad], collapse = ", "))
    rows <- mu_expand_tau0(do.call(rbind, lapply(red, `[[`, "rows")), cfg$kappa)
    rows$Design <- me_design_name(rows$Design_ID)
    say("%s: %s fits in %d tasks; total task time %.1f core-h (%.3f ms/fit)", prof,
        format(sum(sapply(red, `[[`, "n_fits")), big.mark = ","), length(red),
        sum(sapply(red, `[[`, "elapsed")) / 3600, 1000 * sum(sapply(red, `[[`, "elapsed")) / sum(sapply(red, `[[`, "n_fits")))
    say("%s failures: estimates missing %d, CI missing %d, boundary %d, errors %d, aliased %d, unidentified %d, warnings %d (of %s fits)",
        prof, sum(rows$N_Fit - rows$N_Est), sum(rows$N_Fit - rows$N_CI), sum(rows$N_Bound), sum(rows$N_Error),
        sum(rows$N_Alias), sum(rows$N_Unidentified), sum(rows$N_Warn), format(sum(rows$N_Fit), big.mark = ","))
    write.csv(rows, file.path(out, paste0(prof, "_settings.csv")), row.names = FALSE)
    units <- list(); for (r in red) for (ci in seq_along(r$units))
      units[[paste(r$task, ci, sep = "|")]] <- r$units[[ci]]
    res[[prof]] <- list(rows = rows, units = units, red = red, cfg = cfg)
    al <- do.call(rbind, lapply(red, `[[`, "alloc")); sb <- do.call(rbind, lapply(red, `[[`, "subsets"))
    saveRDS(list(alloc = al, subsets = sb), file.path(out, paste0(prof, "_allocations_subsets.rds")))
  }
  if (is.null(res$a3)) stop("No A3 results")
  a3 <- res$a3$rows; cfg <- res$a3$cfg
  main <- a3[a3$Beta == 0 & a3$Arm == "unit", ]

  ## Diagnostics: islands, guard bypass, small regions ----
  isl <- unique(main[main$Mechanism == "random", c("N", "Year", "Design_ID", "Islands_Mean", "Islands_Max",
                                                    "Subsets_With_Islands", "Guard_Bypass", "Subsets_Region_Lt2", "Subsets_Region_Lt3")])
  isl <- aggregate(cbind(Islands_Mean, Islands_Max, Subsets_With_Islands, Guard_Bypass, Subsets_Region_Lt2, Subsets_Region_Lt3) ~ N,
                   isl[isl$Design_ID == 9 & isl$Year == "2018", ], max)
  write.csv(isl, file.path(out, "a3_subset_diagnostics_random.csv"), row.names = FALSE)
  for (i in seq_len(nrow(isl))) say("Subsets N = %d: islands mean %.2f (max %d), subsets with islands %d/50, symmetric-W bypass %d, region <2: %d, <3: %d",
    isl$N[i], isl$Islands_Mean[i], isl$Islands_Max[i], isl$Subsets_With_Islands[i], isl$Guard_Bypass[i], isl$Subsets_Region_Lt2[i], isl$Subsets_Region_Lt3[i])

  ## V1 reproduction ----
  prod <- ms_perf(); prod <- prod[prod$Beta == 0 & prod$Rho %in% c(0, 0.5), ]
  v1 <- merge(main[main$N == 58 & main$Mechanism == "random" & main$Tau == 1,
                   c("Year", "Design_ID", "Design", "Rho", "Regime", "Kappa", "MSE", "MCSE_MSE", "Coverage", "MCSE_Coverage")],
              prod[c("Year", "Design_ID", "Rho", "Gamma", "Regime", "Mean_MSE", "SE_Mean_MSE_Joint", "Coverage", "MCSE_Coverage")],
              by.x = c("Year", "Design_ID", "Rho", "Regime", "Kappa"), by.y = c("Year", "Design_ID", "Rho", "Regime", "Gamma"),
              suffixes = c("_A3", "_Prod"))
  v1$Z_MSE <- (v1$MSE - v1$Mean_MSE) / sqrt(v1$MCSE_MSE^2 + v1$SE_Mean_MSE_Joint^2)
  v1$Z_Coverage <- (v1$Coverage_A3 - v1$Coverage_Prod) / sqrt(v1$MCSE_Coverage_A3^2 + v1$MCSE_Coverage_Prod^2)
  write.csv(v1, file.path(out, "v1_reproduction.csv"), row.names = FALSE)
  nz <- 2 * nrow(v1); n2 <- sum(abs(v1$Z_MSE) > 2) + sum(abs(v1$Z_Coverage) > 2); n3 <- sum(abs(c(v1$Z_MSE, v1$Z_Coverage)) > 3)
  p_binom <- stats::binom.test(n2, nz, 0.0455, alternative = "greater")$p.value
  v1_pass <- nrow(v1) == 288 && p_binom > 0.01 && n3 <= 1   # 4 years x 9 designs x 2 regimes x 2 rho x 2 kappa
  say("V1 %s: %d comparisons (MSE + coverage over %d cells); |z| > 2: %d (expected %.1f by chance, one-sided binomial p = %.3f); |z| > 3: %d; max |z| %.2f",
      if (v1_pass) "PASS" else "FAIL", nz, nrow(v1), n2, 0.0455 * nz, p_binom, n3, max(abs(c(v1$Z_MSE, v1$Z_Coverage))))

  ## V2 size ----
  v2 <- main[main$Tau == 0 & main$Kappa == cfg$kappa[1], ]   # tau = 0 rows are shared across kappa
  v2s <- aggregate(cbind(Size = Power, Size_Adj = Power_Adj, MCSE_Adj = MCSE_Power_Adj) ~ Mechanism + N + Design, v2, mean)
  v2r <- aggregate(cbind(Size_Min = Power) ~ Mechanism + N + Design, v2, min); v2x <- aggregate(cbind(Size_Max = Power) ~ Mechanism + N + Design, v2, max)
  v2s <- merge(merge(v2s, v2r), v2x); write.csv(v2s, file.path(out, "v2_size.csv"), row.names = FALSE)
  for (n in sort(unique(v2s$N))) { d <- v2s[v2s$N == n & v2s$Mechanism == "random", ]
    say("V2 N = %2d (random): Wald size by design (mean over settings) %.3f-%.3f, within 5-8%% for every design: %s; settings %.3f-%.3f; size-adjusted %.3f-%.3f",
        n, min(d$Size), max(d$Size), if (all(d$Size >= 0.05 & d$Size <= 0.08)) "yes" else "NO", min(d$Size_Min), max(d$Size_Max),
        min(d$Size_Adj), max(d$Size_Adj)) }

  ## V8 monotonicity ----
  key <- c("Mechanism", "Year", "Design_ID", "Rho", "Regime", "Kappa")
  v8 <- list()
  for (metric in c("Power", "Power_Adj")) {
    se <- paste0("MCSE_", metric)
    sp <- split(main, main[c(key, "N")], drop = TRUE)
    for (d in sp) { d <- d[order(d$Tau), ]; dd <- diff(d[[metric]]); tol <- 2 * sqrt(head(d[[se]], -1)^2 + tail(d[[se]], -1)^2)
      if (any(dd < -tol, na.rm = TRUE)) v8[[length(v8) + 1L]] <- data.frame(Metric = metric, Along = "tau", d[1, c(key, "N")], Worst_Drop = min(dd)) }
    sp <- split(main[main$Mechanism == "random", ], main[main$Mechanism == "random", c(key, "Tau")], drop = TRUE)
    for (d in sp) { d <- d[order(d$N), ]; dd <- diff(d[[metric]]); tol <- 2 * sqrt(head(d[[se]], -1)^2 + tail(d[[se]], -1)^2)
      if (any(dd < -tol, na.rm = TRUE)) v8[[length(v8) + 1L]] <- data.frame(Metric = metric, Along = "N", d[1, c(key, "Tau")], Worst_Drop = min(dd)) }
  }
  v8 <- if (length(v8)) do.call(rbind, lapply(v8, function(x) { for (m in c("N", "Tau")) if (is.null(x[[m]])) x[[m]] <- NA; x })) else data.frame()
  write.csv(v8, file.path(out, "v8_monotonicity_violations.csv"), row.names = FALSE)
  n_curves <- 2 * (nrow(unique(main[c(key, "N")])) + nrow(unique(main[main$Mechanism == "random", c(key, "Tau")])))
  say("V8: %d of %d power curves (Wald and size-adjusted, along tau and N) drop by more than 2 MC SEs (flagged, not smoothed)", nrow(v8), n_curves)

  ## N80 and MDE with bootstrap over subsets ----
  rnd <- main[main$Mechanism == "random", ]
  # Map each row to its unit counts: cells are identical across tasks of one (mechanism, N) shape.
  cell_key <- function(d) paste(d$Beta, d$Rho, d$Regime, d$Arm, d$Tau, ifelse(d$Shared_Tau0, "NA", d$Kappa))
  red_index <- list()
  for (r in res$a3$red) red_index[[r$task]] <- paste(r$cells$Beta, r$cells$Rho, r$cells$Regime, r$cells$Arm, r$cells$Tau, ifelse(is.na(r$cells$Kappa), "NA", r$cells$Kappa))
  get_units <- function(d) { tk <- sprintf("%s_%s_N%02d_D%d", d$Mechanism, d$Year, d$N, d$Design_ID)
    lapply(seq_len(nrow(d)), function(i) res$a3$units[[paste(tk[i], match(cell_key(d[i, ]), red_index[[tk[i]]]), sep = "|")]]) }
  boot_power <- function(u, adj) { if (is.null(u)) return(NULL)
    num <- if (adj) "n_rej_adj" else "n_rej"
    replicate(B, { i <- sample.int(nrow(u), replace = TRUE); sum(u[i, num]) / sum(u[i, "n_ci"]) }) }
  curves <- split(rnd[rnd$Tau > 0, ], rnd[rnd$Tau > 0, c("Year", "Design_ID", "Rho", "Regime", "Kappa", "Tau")], drop = TRUE)
  n80 <- parallel::mclapply(curves, function(d) {
    d <- d[order(d$N), ]; u <- get_units(d)
    o <- d[1, c("Year", "Design_ID", "Design", "Rho", "Regime", "Kappa", "Tau")]
    for (adj in c(FALSE, TRUE)) {
      p <- if (adj) d$Power_Adj else d$Power
      bp <- sapply(seq_len(nrow(d)), function(i) { b <- boot_power(u[[i]], adj); if (is.null(b)) rep(p[i], B) else b })
      bn <- apply(bp, 1, function(pp) mu_cross(d$N, pp))
      nm <- if (adj) "N80_Adj" else "N80"
      o[[nm]] <- mu_cross(d$N, p); o[[paste0(nm, "_lo")]] <- unname(stats::quantile(bn, 0.025, type = 1))
      o[[paste0(nm, "_hi")]] <- unname(stats::quantile(bn, 0.975, type = 1))
    }
    pick <- function(n) if (any(d$N == n)) d$Power[d$N == n] else NA_real_
    o$Power_29 <- pick(29); o$Power_58 <- pick(58)
    o
  }, mc.cores = workers)
  if (any(vapply(n80, inherits, logical(1), "try-error"))) stop("N80 failed: ", as.character(n80[vapply(n80, inherits, logical(1), "try-error")][[1]]))
  n80 <- do.call(rbind, n80)
  write.csv(n80, file.path(out, "a3_n80_by_year.csv"), row.names = FALSE)
  agg <- function(v) c(mean = if (all(is.finite(v))) mean(v) else NA_real_, min = min(v), max = max(v), n_censored_high = sum(v == Inf), n_censored_low = sum(v == -Inf))
  n80s <- do.call(rbind, lapply(split(n80, n80[c("Design_ID", "Rho", "Regime", "Kappa", "Tau")], drop = TRUE), function(d)
    cbind(d[1, c("Design_ID", "Design", "Rho", "Regime", "Kappa", "Tau")],
          t(setNames(agg(d$N80), paste0("N80_", c("Mean", "Min", "Max", "Years_Above58", "Years_AtMost15")))),
          t(setNames(agg(d$N80_Adj), paste0("N80_Adj_", c("Mean", "Min", "Max", "Years_Above58", "Years_AtMost15")))),
          N80_lo_min = min(d$N80_lo), N80_hi_max = max(d$N80_hi))))
  write.csv(n80s, file.path(out, "a3_n80_summary.csv"), row.names = FALSE)

  mde <- do.call(rbind, lapply(split(rnd, rnd[c("Year", "Design_ID", "Rho", "Regime", "Kappa", "N")], drop = TRUE), function(d) {
    d <- d[order(d$Tau), ]
    cbind(d[1, c("Year", "Design_ID", "Design", "Rho", "Regime", "Kappa", "N")],
          MDE = mu_cross(d$Tau, d$Power), MDE_Adj = mu_cross(d$Tau, d$Power_Adj))
  }))
  write.csv(mde, file.path(out, "a3_mde_by_year.csv"), row.names = FALSE)
  mdes <- aggregate(cbind(MDE, MDE_Adj) ~ Design_ID + Design + Rho + Regime + Kappa + N, mde,
                    function(v) if (all(is.finite(v))) mean(v) else Inf)
  write.csv(mdes, file.path(out, "a3_mde_summary.csv"), row.names = FALSE)

  ## Translation tables (beta = 0, random volunteering; exact by scale invariance, V3) ----
  grid <- expand.grid(m = c(10, 25, 50), ICC = c(0.01, 0.05, 0.10), d = c(0.1, 0.2, 0.3, 0.5))
  grid$Sigma_c <- sqrt(grid$ICC + (1 - grid$ICC) / grid$m); grid$Tau_c <- grid$d / grid$Sigma_c
  pw <- aggregate(cbind(Power, Power_Adj) ~ Design_ID + Design + Rho + Regime + Kappa + N + Tau, rnd, mean)
  tr <- do.call(rbind, lapply(split(pw, pw[c("Design_ID", "Rho", "Regime", "Kappa", "N")], drop = TRUE), function(d) {
    d <- d[order(d$Tau), ]
    cbind(d[1, c("Design_ID", "Design", "Rho", "Regime", "Kappa", "N")][rep(1, nrow(grid)), ], grid,
          Power = mu_power_at(d$Tau, d$Power, grid$Tau_c), Power_Adj = mu_power_at(d$Tau, d$Power_Adj, grid$Tau_c),
          Beyond_Grid = grid$Tau_c > max(d$Tau), Power_At_Grid_Top = d$Power[nrow(d)])
  }))
  rownames(tr) <- NULL
  write.csv(tr, file.path(out, "a3_translation_m_icc.csv"), row.names = FALSE)
  say("Translation: %d of %d (m, ICC, d) cells have standardized effects beyond tau = 2 (power reported as >= power at tau = 2)",
      sum(grid$Tau_c > 2), nrow(grid))

  ## Guidance holds: MSE ratio to SRS at tau = 1 by N, with A2 tests per N ----
  g1 <- main[main$Tau == 1, ]
  srs <- g1[g1$Design_ID == 9, c("Mechanism", "Year", "N", "Rho", "Regime", "Kappa", "MSE")]
  g1 <- merge(g1, srs, by = c("Mechanism", "Year", "N", "Rho", "Regime", "Kappa"), suffixes = c("", "_SRS"))
  g1$Ratio <- g1$MSE / g1$MSE_SRS
  gh <- aggregate(Ratio ~ Mechanism + N + Regime + Design_ID + Design, g1, function(v) c(Mean = mean(v), Min = min(v), Max = max(v)))
  gh <- do.call(data.frame, gh)
  write.csv(gh, file.path(out, "a3_guidance_ratio_to_srs.csv"), row.names = FALSE)
  tests <- list(); trows <- list()
  for (mech in unique(g1$Mechanism)) for (n in sort(unique(g1$N[g1$Mechanism == mech]))) for (rg in cfg$regimes) {
    d <- g1[g1$Mechanism == mech & g1$N == n & g1$Regime == rg, ]
    d$Mean_MSE <- d$MSE
    lab <- sprintf("%s_N%02d_%s", mech, n, rg)
    tt <- tryCatch(me_design_tests(d, c("Year", "Rho", "Kappa")), error = function(e) { say("A2 per N %s: test not run (%s)", lab, conditionMessage(e)); NULL })
    if (is.null(tt)) next
    tests[[lab]] <- tt
    trows[[lab]] <- data.frame(Mechanism = mech, N = n, Regime = rg, Blocks = tt$n_blocks,
      Friedman = tt$friedman[["statistic"]], P = tt$friedman[["p"]], CD = tt$cd,
      Design = names(tt$avg_ranks), Avg_Rank = unname(tt$avg_ranks))
  }
  write.csv(do.call(rbind, trows), file.path(out, "a2_tests_per_N_ranks.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(names(tests), function(l) cbind(Setting = l, tests[[l]]$wilcoxon))),
            file.path(out, "a2_tests_per_N_pairwise.csv"), row.names = FALSE)
  saveRDS(tests, file.path(out, "a2_tests_per_N.rds"))
  pdf(file.path(out, "a2_cd_plots_per_N.pdf"), width = 9, height = 4)
  for (l in names(tests)) if (grepl("N(20|29|40|58)", l)) print(me_cd_plot(tests[[l]], paste("Consistency summary:", l, "(16 blocks: year x rho x kappa)")))
  dev.off()

  ## Saturation designs at small N: restrict to subsets whose regions all have >= 2 colleges ----
  ss <- readRDS(file.path(out, "a3_allocations_subsets.rds"))$subsets
  keep_sub <- unique(ss[ss$Design_ID == 3 & ss$Min_Region_Size >= 2, c("Mechanism", "Year", "N", "Subset")])
  rs <- list()
  for (r in res$a3$red) {
    tk <- r$task; parts <- strsplit(tk, "_")[[1]]
    if (r$rows$N[1] == 58) next
    ks <- keep_sub$Subset[keep_sub$Mechanism == r$rows$Mechanism[1] & keep_sub$Year == r$rows$Year[1] & keep_sub$N == r$rows$N[1]]
    for (ci in which(r$cells$Beta == 0 & r$cells$Tau == 1)) {
      u <- r$units[[ci]]
      rs[[length(rs) + 1L]] <- data.frame(Mechanism = r$rows$Mechanism[1], Year = r$rows$Year[1], N = r$rows$N[1],
        Design_ID = r$rows$Design_ID[1], Rho = r$cells$Rho[ci], Regime = r$cells$Regime[ci], Kappa = r$cells$Kappa[ci],
        MSE_All = sum(u[, "sum_err2"]) / sum(u[, "n_est"]), MSE_Regions2 = sum(u[ks, "sum_err2"]) / sum(u[ks, "n_est"]),
        Subsets_Kept = length(ks))
    }
  }
  rs <- do.call(rbind, rs)
  rs <- merge(rs, rs[rs$Design_ID == 9, c("Mechanism", "Year", "N", "Rho", "Regime", "Kappa", "MSE_All", "MSE_Regions2")],
              by = c("Mechanism", "Year", "N", "Rho", "Regime", "Kappa"), suffixes = c("", "_SRS"))
  rs$Ratio_All <- rs$MSE_All / rs$MSE_All_SRS; rs$Ratio_Regions2 <- rs$MSE_Regions2 / rs$MSE_Regions2_SRS
  rss <- aggregate(cbind(Ratio_All, Ratio_Regions2, Subsets_Kept) ~ Mechanism + N + Regime + Design_ID, rs[rs$Design_ID %in% c(3, 8), ], mean)
  write.csv(rss, file.path(out, "a3_saturation_region_restriction.csv"), row.names = FALSE)

  ## A1 per N (share of the participants' deaths) ----
  al <- readRDS(file.path(out, "a3_allocations_subsets.rds"))$alloc
  a1n <- aggregate(cbind(Death_Share, PY_Share, Treated_Fraction = Treated / N) ~ Mechanism + N + Design_ID, al, mean)
  a1n$Enrichment <- a1n$Death_Share / a1n$PY_Share; a1n$Design <- me_design_name(a1n$Design_ID)
  write.csv(a1n, file.path(out, "a1_burden_share_per_N.csv"), row.names = FALSE)

  ## beta = 1 sensitivity beside beta = 0 ----
  b1 <- a3[a3$Arm == "unit" & a3$N %in% cfg$beta1$N & a3$Tau %in% cfg$beta1$tau & a3$Mechanism == "random", ]
  if (any(b1$Beta == 1)) {
    b1s <- aggregate(cbind(Power, Power_Adj, MSE, Coverage) ~ Beta + N + Tau + Regime + Design_ID, b1, mean)
    b1s$Design <- me_design_name(b1s$Design_ID)
    write.csv(b1s, file.path(out, "a3_beta1_vs_beta0.csv"), row.names = FALSE)
  } else say("beta = 1: no rows in this profile")

  ## A4: unequal vs equal trained professionals ----
  if (!is.null(res$a4)) {
    a4 <- res$a4$rows
    a4w <- merge(a4[a4$Arm == "equal_m", ], a4[a4$Arm == "unequal_m", ],
                 by = c("Mechanism", "Year", "N", "Design_ID", "Design", "Beta", "Rho", "Regime", "Tau", "Kappa", "Shared_Tau0"),
                 suffixes = c("_Equal", "_Unequal"))
    a4w$Power_Diff <- a4w$Power_Unequal - a4w$Power_Equal; a4w$MSE_Ratio <- a4w$MSE_Unequal / a4w$MSE_Equal
    write.csv(a4w, file.path(out, "a4_equal_vs_unequal.csv"), row.names = FALSE)
    a4s <- aggregate(cbind(Power_Equal, Power_Unequal, Power_Diff, MSE_Ratio, Size_Unequal = Power_Unequal) ~ N + Regime + Design_ID + Design + Tau, a4w, mean)
    write.csv(a4s, file.path(out, "a4_summary.csv"), row.names = FALSE)
    say("A4: unequal minus equal power at d > 0, by design (mean over settings): %s",
        paste(sprintf("%s %+.3f", me_design_name(1:9), tapply(a4w$Power_Diff[a4w$Tau > 0], a4w$Design_ID[a4w$Tau > 0], mean)), collapse = "; "))
    say("A4: MSE ratio unequal/equal by design: %s",
        paste(sprintf("%s %.3f", me_design_name(1:9), tapply(a4w$MSE_Ratio[a4w$Tau > 0], a4w$Design_ID[a4w$Tau > 0], mean)), collapse = "; "))
  }
  mu_headline(out)
  writeLines(c(sprintf("Run %s", format(Sys.time())), log), file.path(out, "summaries_log.txt"))
  invisible(list(a3 = a3, n80 = n80s, mde = mdes, guidance = gh, v1 = v1))
}

# Headline markdown tables ----
#' Write the headline N80, MDE, guidance and translation tables as markdown
#'
#' Reads the summary CSVs written by mu_run() and formats: colleges for 80% power
#' (size-adjusted and Wald, each kappa), MDE at N = 20/29/40/58, the guidance-holds
#' ratios (random and incidence-weighted) and the m x ICC translation for SRS and
#' Incidence-Guided Saturation at rho = 0.5, kappa = 0.8.
#' @param out Summaries folder.
#' @return Invisibly NULL; writes a3_headline_tables.md and a3_translation_headline.md.
#' @family methods_ext_summaries
#' @examples
#' # mu_headline(file.path(.me_root, "summaries"))
mu_headline <- function(out) {
  n <- read.csv(file.path(out, "a3_n80_summary.csv")); m <- read.csv(file.path(out, "a3_mde_summary.csv")); g <- read.csv(file.path(out, "a3_guidance_ratio_to_srs.csv"))
  ord <- c("Simple Random Sampling","High Incidence Focus","Balanced Halves","Balanced Quartiles","Plain Saturation","Incidence-Guided Saturation","2x2 Blocking","Checkerboard","Isolation Buffer")
  fm <- function(mean, mn, mx) ifelse(is.infinite(mn) & mn > 0, ">58", ifelse(is.na(mean), sprintf("%s–%s", ifelse(is.infinite(mn), ifelse(mn<0,"≤15",">58"), sprintf("%.0f",mn)), ifelse(is.infinite(mx), ifelse(mx<0,"≤15",">58"), sprintf("%.0f",mx))), ifelse(mx - mn < 1.5, sprintf("%.0f", mean), sprintf("%.0f (%.0f–%.0f)", mean, mn, mx))))
  L <- c()
  for (adj in c(TRUE, FALSE)) for (k in c(0.8, 0.5)) {
    pre <- if (adj) "N80_Adj_" else "N80_"
    L <- c(L, sprintf("\n**Colleges for 80%% power (%s), κ = %.1f.** Mean over 2018–2021 (range when years differ by > 1); `>58` = not reached with all colleges.\n", if (adj) "size-adjusted" else "Wald", k))
    cols <- expand.grid(Tau = c(1, 1.5, 2), Rho = c(0, 0.5), Regime = c("both", "control_only"), stringsAsFactors = FALSE)
    L <- c(L, paste0("| Design | ", paste(sprintf("%s ρ=%.1f τ=%.2g", ifelse(cols$Regime=="both","Both","Ctrl"), cols$Rho, cols$Tau), collapse=" | "), " |"), paste0("|", strrep("---|", nrow(cols)+1)))
    for (d in ord) { v <- sapply(seq_len(nrow(cols)), function(i) { r <- n[n$Design==d & n$Kappa==k & n$Tau==cols$Tau[i] & n$Rho==cols$Rho[i] & n$Regime==cols$Regime[i],]
      fm(r[[paste0(pre,"Mean")]], r[[paste0(pre,"Min")]], r[[paste0(pre,"Max")]]) }); L <- c(L, paste0("| ", d, " | ", paste(v, collapse=" | "), " |")) }
  }
  L <- c(L, "\n**Minimum detectable effect (size-adjusted, cluster-SD units), κ = 0.8, mean over years; `>2` = above the grid.**\n")
  cols <- expand.grid(N = c(20, 29, 40, 58), Regime = c("both","control_only"), Rho = c(0, 0.5), stringsAsFactors=FALSE)
  L <- c(L, paste0("| Design | ", paste(sprintf("%s ρ=%.1f N=%d", ifelse(cols$Regime=="both","Both","Ctrl"), cols$Rho, cols$N), collapse=" | "), " |"), paste0("|", strrep("---|", nrow(cols)+1)))
  for (d in ord) { v <- sapply(seq_len(nrow(cols)), function(i) { r <- m[m$Design==d & m$Kappa==0.8 & m$N==cols$N[i] & m$Rho==cols$Rho[i] & m$Regime==cols$Regime[i],]; if (is.infinite(r$MDE_Adj)) ">2" else sprintf("%.2f", r$MDE_Adj) }); L <- c(L, paste0("| ", d, " | ", paste(v, collapse=" | "), " |")) }
  L <- c(L, "\n**Guidance holds? Mean-MSE ratio to SRS at τ = 1 (mean over 16 blocks: year × ρ × κ; range in brackets), random volunteering.**\n")
  cols <- expand.grid(N = c(20, 29, 40, 58), Regime = c("both","control_only"), stringsAsFactors=FALSE)
  L <- c(L, paste0("| Design | ", paste(sprintf("%s N=%d", ifelse(cols$Regime=="both","Both","Ctrl"), cols$N), collapse=" | "), " |"), paste0("|", strrep("---|", nrow(cols)+1)))
  for (d in ord[-1]) { v <- sapply(seq_len(nrow(cols)), function(i) { r <- g[g$Design==d & g$Mechanism=="random" & g$N==cols$N[i] & g$Regime==cols$Regime[i],]; sprintf("%.2f [%.2f–%.2f]", r$Ratio.Mean, r$Ratio.Min, r$Ratio.Max) }); L <- c(L, paste0("| ", d, " | ", paste(v, collapse=" | "), " |")) }
  L <- c(L, "\n**Incidence-weighted (3×) volunteering, same table.**\n")
  cols <- expand.grid(N = c(20, 29), Regime = c("both","control_only"), stringsAsFactors=FALSE)
  L <- c(L, paste0("| Design | ", paste(sprintf("%s N=%d", ifelse(cols$Regime=="both","Both","Ctrl"), cols$N), collapse=" | "), " |"), paste0("|", strrep("---|", nrow(cols)+1)))
  for (d in ord[-1]) { v <- sapply(seq_len(nrow(cols)), function(i) { r <- g[g$Design==d & g$Mechanism=="incidence3x" & g$N==cols$N[i] & g$Regime==cols$Regime[i],]; sprintf("%.2f [%.2f–%.2f]", r$Ratio.Mean, r$Ratio.Min, r$Ratio.Max) }); L <- c(L, paste0("| ", d, " | ", paste(v, collapse=" | "), " |")) }
  writeLines(L, file.path(out, "a3_headline_tables.md"))
  t <- read.csv(file.path(out, "a3_translation_m_icc.csv")); t <- t[t$Rho==0.5 & t$Kappa==0.8,]
  L <- "\n**Translation (size-adjusted power), ρ = 0.5, κ = 0.8, random volunteering.** Individual effect d (in individual-SD units), m trained professionals per college, ICC; σ_c = √(ICC + (1 − ICC)/m), τ = d/σ_c. `≥x` = effect beyond the simulated τ ≤ 2, power at least x.\n"
  for (reg in c("both","control_only")) for (des in c("Simple Random Sampling","Incidence-Guided Saturation")) {
    L <- c(L, sprintf("\n*%s, %s spillover*\n", des, if (reg=="both") "both-arms" else "control-only"))
    cols <- expand.grid(N=c(20,29,40,58), d=c(0.1,0.2,0.3)); L <- c(L, paste0("| m | ICC | ", paste(sprintf("d=%.1f N=%d", cols$d, cols$N), collapse=" | "), " |"), paste0("|", strrep("---|", nrow(cols)+2)))
    for (m in c(10,25,50)) for (icc in c(0.01,0.05,0.10)) { v <- sapply(seq_len(nrow(cols)), function(i) { r <- t[t$Design==des & t$Regime==reg & t$N==cols$N[i] & t$d==cols$d[i] & t$m==m & t$ICC==icc,]
       if (r$Beyond_Grid) sprintf("≥%.2f", r$Power_At_Grid_Top) else sprintf("%.2f", r$Power_Adj) }); L <- c(L, sprintf("| %d | %.2f | %s |", m, icc, paste(v, collapse=" | "))) }
  }
  writeLines(L, file.path(out, "a3_translation_headline.md"))
  invisible(NULL)
}

#' Task names of a profile, from its stored configuration
#' @param cfg Configuration from the manifest.
#' @return Character task names.
#' @family methods_ext_summaries
#' @examples
#' # mp_tasks_from(cfg)
mp_tasks_from <- function(cfg) {
  if (!exists("mp_tasks")) source(file.path(.mu_dir, "methods_ext_participation.R"))
  mp_tasks(cfg)$Task
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  mu_run(if (length(args)) as.integer(args[1]) else 8L)
}
