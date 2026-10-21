# ============================================================
# Script: methods_ext_sensitivity.R
# Purpose: N = 58 sensitivity runs: A5 Stage 2 coverage calibration (parametric
#   bootstrap, profile likelihood), A6 estimation method (CAR, reverse
#   misspecification, Bayesian CAR) and A7 spillover forms (binary, additive count).
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: methods_ext_core.R, spatialreg, spdep, CARBayes, coda, parallel
# ============================================================
# Run from the project root (BLAS pinned to one thread):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript application/code/methods_ext_sensitivity.R <analysis> <workers> [smoke]
# analysis: a5 | a6car | a6bayes | a6bayeslong | a6bayeslag | a7 ; "smoke" runs a reduced replicate count.
# Outputs: application/results/methods_ext_2026-10/<analysis>[_smoke]/ with
#   manifest.rds, checkpoints/<year>_D<id>.rds (gitignored) and fits/summary CSVs.
# One task = (year, design) [A7: (form, year, design)]; each worker holds one task.
#
# Shared N = 58 setup: frozen production regions/blocks/ranks, queen W, designs via
# me_draw_design() (identical to production), every fit's mean = intercept + Z +
# Spill (gamma-scaled exposure, as production), key-seeded allocations and noise.

.mx_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(p)) return(normalizePath(dirname(p)))
  }
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else stop("Run or source by path")
})
if (!exists("me_seed")) source(file.path(.mx_dir, "methods_ext_core.R"))

# Configuration ----
#' Settings and replicate counts for each sensitivity analysis
#' @param analysis "a5", "a6car", "a6bayes" or "a7".
#' @param smoke Reduced replicate counts for a smoke check.
#' @return Configuration list.
#' @family methods_ext_sensitivity
#' @examples
#' mx_config("a7")
mx_config <- function(analysis, smoke = FALSE) {
  c0 <- list(analysis = analysis, smoke = smoke, years = as.character(2018:2021), designs = 1:9,
             regimes = c("both", "control_only"), tau = 1, sigma = 1)
  cfg <- switch(analysis,
    a5 = c(c0, list(rho = c(0, 0.5), gamma = 0.8, J = 25L, R = 20L, det_R = 500L, B = 499L)),
    a6car = c(c0, list(rho = c(0, 0.2, 0.5), gamma = c(0.5, 0.8), J = 25L, R = 20L, det_R = 500L,
                       dgp = c("sar", "car"))),
    a6bayes = c(c0, list(rho = c(0, 0.5), gamma = 0.8, J = 10L, R = 10L, det_R = 100L,
                         burnin = 5000L, n_sample = 25000L, thin = 5L, ess_min = 400, check_frac = 1)),
    # 2026-10-21 (author): rho = 0.5 cells again with chains 4x longer, on the same
    # data as a6bayes (same seed prefix), so convergence is compared fit by fit.
    a6bayeslong = c(c0, list(rho = 0.5, gamma = 0.8, J = 10L, R = 10L, det_R = 100L,
                             burnin = 20000L, n_sample = 100000L, thin = 20L, ess_min = 400, check_frac = 1)),
    # 2026-10-21 (author): Bayesian SAR lag (spatialreg::spBreg_lag), the correctly
    # specified Bayesian model, on the same data as a6bayes (both rho), to separate
    # "Bayesian" from "CAR" in the a6bayes results.
    a6bayeslag = c(c0, list(rho = c(0, 0.5), gamma = 0.8, J = 10L, R = 10L, det_R = 100L,
                            ndraw = 10000L, nomit = 2000L, ess_min = 400)),
    a7 = c(c0, list(rho = c(0, 0.5), gamma = c(0.5, 0.8), J = 100L, R = 100L, det_R = 1000L,
                    forms = c("proportional", "binary", "count"))),
    stop("Unknown analysis"))
  if (smoke) { cfg$years <- "2018"; cfg$designs <- c(1, 4, 8, 9); cfg$J <- 3L; cfg$R <- 4L; cfg$det_R <- 6L
    if (analysis == "a5") cfg$B <- 99L
    if (analysis == "a6bayeslag") { cfg$ndraw <- 2000L; cfg$nomit <- 500L }
    if (analysis %in% c("a6bayes", "a6bayeslong")) { cfg$burnin <- 1000L; cfg$n_sample <- 3000L; cfg$thin <- 1L } }
  cfg
}

# Fitters ----
#' Profile log-likelihood of the SAR lag model at a fixed tau
#'
#' Model y = rho W y + tau z + X b + e. At fixed tau the residual for given rho is
#' M[(y - tau z) - rho W y] with M the projector off X (intercept, Spill), so
#'   SSE(rho) = e_a - 2 rho e_b + rho^2 e_c,  e_a = |M(y - tau z)|^2,
#'   e_b = (M W y)'M(y - tau z),  e_c = |M W y|^2,
#'   l(rho) = log|I - rho W| - (n/2) log(2 pi) - (n/2) log(SSE/n) - n/2,
#' the same concentrated likelihood fit_sar_lag() maximizes (Wy is NOT offset).
#' @param tau Fixed treatment effect.
#' @param y Outcome. @param z Treatment. @param xn Nuisance matrix (intercept, Spill).
#' @param engine sar_lag_setup() list.
#' @return Maximized log-likelihood over rho (and b, sigma^2).
#' @family methods_ext_sensitivity
#' @examples
#' # mx_profile_ll(1, y, z, xm[, -2], g58$engine)
mx_profile_ll <- function(tau, y, z, xn, engine) {
  n <- engine$n; eig <- engine$eig
  e0 <- stats::lm.fit(xn, y - tau * z)$residuals
  eW <- stats::lm.fit(xn, as.vector(engine$W %*% y))$residuals
  ea <- sum(e0^2); eb <- sum(eW * e0); ec <- sum(eW^2)
  ldet <- if (is.complex(eig)) function(r) Re(sum(log(1 - r * eig))) else function(r) sum(log(1 - r * eig))
  ll <- function(r) { SSE <- ea - 2 * r * eb + r * r * ec; ldet(r) - (n / 2) * log(2 * pi) - (n / 2) * log(SSE / n) - n / 2 }
  stats::optimize(ll, engine$interval, maximum = TRUE, tol = .Machine$double.eps^0.5)$objective
}

#' Profile-likelihood (LR) 95% interval for tau
#' {tau : 2 [l(tau-hat) - l_p(tau)] <= chi2_{1,0.95}}, roots by uniroot on each side.
#' @param fit fit_sar_lag() result on the data.
#' @param y,z,xn,engine As mx_profile_ll().
#' @return c(lo, hi), NA where a root could not be bracketed.
#' @family methods_ext_sensitivity
#' @examples
#' # mx_profile_ci(fit, y, z, xn, engine)
mx_profile_ci <- function(fit, y, z, xn, engine) {
  th <- fit$coefficients[["Z"]]; se <- fit$rest.se[["Z"]]; q <- stats::qchisq(0.95, 1)
  f <- function(t) 2 * (fit$LL - mx_profile_ll(t, y, z, xn, engine)) - q
  side <- function(dir) {
    for (k in c(2, 4, 8, 16)) { b <- th + dir * k * se; if (f(b) > 0) return(stats::uniroot(f, sort(c(th, b)), tol = 1e-8)$root) }
    NA_real_
  }
  c(side(-1), side(1))
}

#' Parametric bootstrap percentile and bootstrap-t intervals for tau
#' y* = (I - rho-hat W)^{-1} (X b-hat + s-hat eps*), refit by fit_sar_lag().
#' Percentile: [q_0.025, q_0.975] of tau*. Bootstrap-t: t* = (tau* - tau-hat)/SE*,
#' interval [tau-hat - q_0.975(t*) SE, tau-hat - q_0.025(t*) SE]. Failed refits
#' (boundary rho, errors) are excluded and counted.
#' @param fit Original fit. @param xm Model matrix. @param engine Setup.
#' @param B Bootstrap replicates. @param key Seed key fields.
#' @return c(pct_lo, pct_hi, t_lo, t_hi, n_fail).
#' @family methods_ext_sensitivity
#' @examples
#' # mx_boot(fit, xm, engine, 499, list("a5-boot", 1))
mx_boot <- function(fit, xm, engine, B, key) {
  do.call(me_seed, key)
  Ainv <- solve(diag(engine$n) - fit$rho * engine$W)
  mu <- as.vector(Ainv %*% (xm %*% fit$coefficients))
  Ys <- mu + Ainv %*% matrix(stats::rnorm(engine$n * B, 0, sqrt(fit$s2)), engine$n)
  r <- t(vapply(seq_len(B), function(b) me_fit(Ys[, b], xm, engine)[c("tau", "se")], numeric(2)))
  ok <- is.finite(r[, 1]) & is.finite(r[, 2]) & r[, 2] > 0
  th <- fit$coefficients[["Z"]]; se <- fit$rest.se[["Z"]]
  ts <- (r[ok, 1] - th) / r[ok, 2]
  c(pct_lo = unname(stats::quantile(r[ok, 1], 0.025)), pct_hi = unname(stats::quantile(r[ok, 1], 0.975)),
    t_lo = th - unname(stats::quantile(ts, 0.975)) * se, t_hi = th - unname(stats::quantile(ts, 0.025)) * se,
    n_fail = sum(!ok))
}

#' Gaussian CAR fit (spautolm, ML) with symmetric similar weights
#' Var(u) = sigma^2 (I - lambda S)^{-1}, S = D^{-1/2} A D^{-1/2} (same eigenvalues as
#' the row-standardized W, so lambda is on rho's scale).
#' @param y Outcome. @param xm Model matrix. @param slw similar.listw() object.
#' @return c(tau, se, lambda, bound, error, warn).
#' @family methods_ext_sensitivity
#' @examples
#' # mx_fit_car(y, xm, slw)
mx_fit_car <- function(y, xm, slw) {
  warns <- character(0)
  df <- data.frame(y = y, xm[, -1, drop = FALSE])
  f <- withCallingHandlers(tryCatch(spatialreg::spautolm(stats::as.formula(paste("y ~", paste(colnames(xm)[-1], collapse = " + "))),
          data = df, listw = slw, family = "CAR"), error = function(e) { warns <<- c(warns, paste("ERROR:", conditionMessage(e))); NULL }),
        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(f)) return(c(tau = NA, se = NA, lambda = NA, bound = 0, error = 1, warn = 1))
  cf <- summary(f)$Coef
  bound <- min(abs(f$lambda - f$interval)) < 1e-6
  c(tau = if (bound) NA else cf["Z", 1], se = if (bound) NA else cf["Z", 2], lambda = unname(f$lambda),
    bound = bound, error = 0, warn = length(warns) > 0)
}

#' Bayesian Leroux CAR fit (CARBayes::S.CARleroux, Gaussian)
#' Posterior mean and 95% equal-tailed credible interval of the Z coefficient,
#' with its effective sample size and Geweke z (coda).
#' @param y,xm Outcome and model matrix. @param A Binary adjacency.
#' @param cfg Configuration (burn-in, samples, thinning). @param key Seed key.
#' @return c(tau, lo, hi, ess, geweke, error).
#' @family methods_ext_sensitivity
#' @examples
#' # mx_fit_bcar(y, xm, A, cfg, list("x"))
mx_fit_bcar <- function(y, xm, A, cfg, key) {
  do.call(me_seed, key)
  df <- data.frame(y = y, xm[, -1, drop = FALSE])
  f <- tryCatch(suppressMessages(CARBayes::S.CARleroux(stats::as.formula(paste("y ~", paste(colnames(xm)[-1], collapse = " + "))),
         family = "gaussian", data = df, W = A, burnin = cfg$burnin, n.sample = cfg$n_sample,
         thin = cfg$thin, verbose = FALSE)), error = function(e) NULL)
  if (is.null(f)) return(c(tau = NA, lo = NA, hi = NA, ess = NA, geweke = NA, error = 1))
  s <- coda::as.mcmc(as.numeric(f$samples$beta[, which(colnames(xm) == "Z")]))
  c(tau = mean(s), lo = unname(stats::quantile(s, 0.025)), hi = unname(stats::quantile(s, 0.975)),
    ess = unname(coda::effectiveSize(s)), geweke = unname(coda::geweke.diag(s)$z), error = 0)
}

#' Bayesian SAR lag fit (spatialreg::spBreg_lag) with convergence diagnostics
#'
#' Same model as the ML fit, y = rho W y + X b + e, estimated by MCMC with
#' spBreg_lag's default priors, cfg$ndraw draws with the first cfg$nomit discarded.
#' Returns the posterior mean and 2.5/97.5% quantiles of the Z coefficient, its
#' effective sample size and Geweke z (coda).
#' @param y,xm Outcome and model matrix. @param lw Row-standardized queen listw.
#' @param cfg Configuration (ndraw, nomit). @param key Seed key.
#' @return c(tau, lo, hi, ess, geweke, error).
#' @family methods_ext_sensitivity
#' @examples
#' # mx_fit_blag(y, xm, lw, cfg, list("x"))
mx_fit_blag <- function(y, xm, lw, cfg, key) {
  do.call(me_seed, key)
  df <- data.frame(y = y, xm[, -1, drop = FALSE])
  f <- tryCatch(suppressMessages(spatialreg::spBreg_lag(stats::as.formula(paste("y ~", paste(colnames(xm)[-1], collapse = " + "))),
         data = df, listw = lw, control = list(ndraw = cfg$ndraw, nomit = cfg$nomit))), error = function(e) NULL)
  if (is.null(f)) return(c(tau = NA, lo = NA, hi = NA, ess = NA, geweke = NA, error = 1))
  s <- coda::as.mcmc(as.numeric(f[, "Z"]))
  c(tau = mean(s), lo = unname(stats::quantile(s, 0.025)), hi = unname(stats::quantile(s, 0.975)),
    ess = unname(coda::effectiveSize(s)), geweke = unname(coda::geweke.diag(s)$z), error = 0)
}

# One task ----
#' Run all settings for one (year, design) [or (form, year, design)] task
#' @param task Task row. @param cfg Configuration. @param setup Frozen setup.
#' @return list(task, fits data frame, elapsed).
#' @family methods_ext_sensitivity
#' @examples
#' # mx_task(mx_tasks(cfg)[1, ], cfg, me_load_setup())
mx_task <- function(task, cfg, setup) {
  t0 <- proc.time()[[3]]
  g <- me_geometry(1:58, setup); year <- task$Year; id <- task$Design_ID
  di <- me_design_inputs(g, year, setup, list("random", 58, 1))
  det <- me_deterministic(id, di); J <- if (det) 1L else cfg$J; R <- if (det) cfg$det_R else cfg$R
  form <- if (is.null(task$Form)) "proportional" else task$Form
  dbar <- mean(rowSums(setup$weights$queen$adjacency))
  slw <- spatialreg::similar.listw(spdep::nb2listw(setup$weights$queen$nb, style = "W"))
  g$lw_full <- spdep::nb2listw(setup$weights$queen$nb, style = "W")   # Bayesian SAR lag (a6bayeslag)
  S <- spdep::listw2mat(slw); A01 <- setup$weights$queen$adjacency
  pref <- if (cfg$analysis %in% c("a6bayeslong", "a6bayeslag")) "a6bayes" else cfg$analysis   # same data as a6bayes
  out <- list()
  for (j in seq_len(J)) {
    me_seed(paste0(pref, "-Z"), year, id, j); z <- me_draw_design(id, di)
    me_seed(paste0(pref, "-eps"), year, id, j); E <- matrix(stats::rnorm(58 * R), 58)
    ex <- me_exposure_form(z, g, form, dbar)
    for (reg in cfg$regimes) for (gm in cfg$gamma) for (r in cfg$rho) {
      sp <- gm * ex * if (reg == "control_only") (1 - z) else 1
      xm <- cbind("(Intercept)" = 1, Z = z, Spill = sp)
      ident <- sum(stats::lm.fit(xm[, -2, drop = FALSE], z)$residuals^2) > 1e-10
      dgps <- if (cfg$analysis == "a6car") cfg$dgp else "sar"
      for (dg in dgps) {
        if (dg == "sar") { Ar <- solve(diag(58) - r * g$WP); Y <- as.vector(Ar %*% (cfg$tau * z + sp)) + Ar %*% (cfg$sigma * E) }
        else { L <- t(chol(solve(diag(58) - r * S))); Y <- (cfg$tau * z + sp) + L %*% (cfg$sigma * E) }
        sar <- t(vapply(seq_len(R), function(k) me_fit(Y[, k], xm, g$engine), numeric(6)))
        d <- data.frame(Year = year, Design_ID = id, Form = form, Regime = reg, Gamma = gm, Rho = r,
                        DGP = dg, Allocation = j, Replicate = seq_len(R), Identified = ident,
                        SAR_tau = if (ident) sar[, 1] else NA_real_, SAR_se = if (ident) sar[, 2] else NA_real_,
                        SAR_alias = sar[, 3], SAR_bound = sar[, 4], SAR_error = sar[, 5], SAR_warn = sar[, 6],
                        stringsAsFactors = FALSE)
        if (cfg$analysis == "a5") {
          ex5 <- t(vapply(seq_len(R), function(k) {
            y <- Y[, k]
            if (!ident || !is.finite(sar[k, 1])) return(rep(NA_real_, 8))
            fit <- fit_sar_lag(y, xm, g$engine)
            pc <- mx_profile_ci(fit, y, z, xm[, -2, drop = FALSE], g$engine)
            bt <- mx_boot(fit, xm, g$engine, cfg$B, list("a5-boot", year, id, j, reg, gm, r, k))
            c(pc, bt[["pct_lo"]], bt[["pct_hi"]], bt[["t_lo"]], bt[["t_hi"]], bt[["n_fail"]],
              abs(mx_profile_ll(fit$coefficients[["Z"]], y, z, xm[, -2, drop = FALSE], g$engine) - fit$LL))
          }, numeric(8)))
          colnames(ex5) <- c("Profile_lo", "Profile_hi", "Pct_lo", "Pct_hi", "Bt_lo", "Bt_hi", "Boot_fail", "Profile_LL_check")
          d <- cbind(d, ex5)
        }
        if (cfg$analysis == "a6car") {
          cr <- t(vapply(seq_len(R), function(k) mx_fit_car(Y[, k], xm, slw), numeric(6)))
          colnames(cr) <- paste0("CAR_", c("tau", "se", "lambda", "bound", "error", "warn"))
          if (!ident) cr[, 1:2] <- NA
          d <- cbind(d, cr)
        }
        if (cfg$analysis %in% c("a6bayes", "a6bayeslong")) {
          bc <- t(vapply(seq_len(R), function(k) mx_fit_bcar(Y[, k], xm, A01, cfg, list("a6-mcmc", year, id, j, reg, gm, r, k)), numeric(6)))
          colnames(bc) <- paste0("BCAR_", c("tau", "lo", "hi", "ess", "geweke", "error"))
          d <- cbind(d, bc)
        }
        if (cfg$analysis == "a6bayeslag") {
          bl <- t(vapply(seq_len(R), function(k) mx_fit_blag(Y[, k], xm, g$lw_full, cfg, list("a6-blag", year, id, j, reg, gm, r, k)), numeric(6)))
          colnames(bl) <- paste0("BLAG_", c("tau", "lo", "hi", "ess", "geweke", "error"))
          if (!ident) bl[, 1:3] <- NA
          d <- cbind(d, bl)
        }
        out[[length(out) + 1L]] <- d
      }
    }
  }
  fits <- do.call(rbind, out)
  list(task = task, fits = fits, elapsed = proc.time()[[3]] - t0, deterministic = det)
}

#' Task list for an analysis
#' @param cfg Configuration.
#' @return Data frame of tasks.
#' @family methods_ext_sensitivity
#' @examples
#' nrow(mx_tasks(mx_config("a7")))
mx_tasks <- function(cfg) {
  t <- expand.grid(Year = cfg$years, Design_ID = cfg$designs, stringsAsFactors = FALSE)
  if (cfg$analysis == "a7") t <- merge(data.frame(Form = setdiff(cfg$forms, character(0)), stringsAsFactors = FALSE), t)
  t$Task <- if (cfg$analysis == "a7") sprintf("%s_%s_D%d", t$Form, t$Year, t$Design_ID) else sprintf("%s_D%d", t$Year, t$Design_ID)
  t[order(t$Design_ID %in% c(1, 2), t$Task), ]
}

# Driver ----
#' Run one sensitivity analysis with forked workers and checkpoints, then write fits
#' @param analysis Analysis name. @param workers Workers. @param smoke Smoke profile.
#' @return Invisibly, all fits.
#' @family methods_ext_sensitivity
#' @examples
#' # mx_run("a7", 10)
mx_run <- function(analysis, workers = 2L, smoke = FALSE) {
  me_check_threads()
  cfg <- mx_config(analysis, smoke); setup <- me_load_setup()
  out <- file.path(.me_root, paste0(analysis, if (smoke) "_smoke" else ""))
  me_manifest(out, cfg, c(file.path(.mx_dir, "methods_ext_core.R"), file.path(.mx_dir, "methods_ext_sensitivity.R")))
  cp <- file.path(out, "checkpoints"); dir.create(cp, recursive = TRUE, showWarnings = FALSE)
  tasks <- mx_tasks(cfg); todo <- which(!file.exists(file.path(cp, paste0(tasks$Task, ".rds"))))
  cat(sprintf("[%s] %s: %d tasks, %d to run\n", format(Sys.time()), analysis, nrow(tasks), length(todo)))
  st <- parallel::mclapply(todo, function(i) {
    r <- tryCatch({ res <- mx_task(tasks[i, ], cfg, setup); saveRDS(res, file.path(cp, paste0(tasks$Task[i], ".rds")))
      sprintf("OK %.0fs %d rows", res$elapsed, nrow(res$fits)) }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("[%s] %s %s\n", format(Sys.time(), "%H:%M:%S"), tasks$Task[i], r)); r
  }, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  bad <- !grepl("^OK", vapply(st, as.character, character(1)))
  if (any(bad)) cat("FAILED TASKS:\n", paste(tasks$Task[todo][bad], unlist(st)[bad], collapse = "\n"), "\n")
  fits <- do.call(rbind, lapply(tasks$Task, function(t) { f <- file.path(cp, paste0(t, ".rds")); if (file.exists(f)) readRDS(f)$fits }))
  saveRDS(fits, file.path(out, "fits.rds"))
  cat(sprintf("[%s] %s done: %d fit rows, %d failed tasks\n", format(Sys.time()), analysis, nrow(fits), sum(bad)))
  invisible(fits)
}

if (!exists("methods_ext_define_only") || !isTRUE(methods_ext_define_only)) {
  args <- commandArgs(TRUE)
  mx_run(args[1], if (length(args) > 1) as.integer(args[2]) else 2L, smoke = length(args) > 2 && args[3] == "smoke")
}
