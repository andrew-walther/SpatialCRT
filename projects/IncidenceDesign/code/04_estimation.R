# ==============================================================================
# 04_estimation.R
# Estimation of the treatment effect tau.
#   - fit_tau_models(): oracle + non-oracle ML spatial-lag fits of one outcome,
#     with every warning captured (used by the 2026-09 runner, 05).
#   - fit_sar_lag() / sar_lag_setup(): lean ML engine, validated against lagsarlm.
#   - estimate_tau(): legacy DIM / lagsarlm path of the pre-2026-09 runner; kept for
#     DIM and for complete_after_mle.R / longleaf_setup, not used by 05.
# ==============================================================================

#' Estimate the treatment effect tau across outcome resamples
#'
#' @param estimation_mode Character, "DIM" or "MLE"
#' @param Y_sim N x n_outcome_resamples matrix of simulated outcomes
#' @param Z Numeric vector length N, binary treatment assignment
#' @param spill_term Numeric vector length N, spillover covariate
#' @param X_col Numeric vector length N, incidence covariate for this resample
#'        (or N x n_outcome_resamples matrix — only used column-wise for MLE)
#' @param X_matrix N x n_outcome_resamples matrix of incidence values
#' @param active_listw listw object for lagsarlm (MLE only)
#' @param n_outcome_resamples Integer, number of outcome columns
#' @param include_spill_covariate Logical, include Spill in MLE model (default TRUE)
#' @return List with:
#'   estimates - numeric vector of tau-hat values
#'   ci_lower  - numeric vector of 95% CI lower bounds
#'   ci_upper  - numeric vector of 95% CI upper bounds
estimate_tau <- function(estimation_mode, Y_sim, Z, spill_term, X_matrix,
                         active_listw, n_outcome_resamples,
                         include_spill_covariate = TRUE) {

  estimates <- numeric(n_outcome_resamples)
  ci_lower  <- numeric(n_outcome_resamples)
  ci_upper  <- numeric(n_outcome_resamples)

  N <- length(Z)
  n_trt  <- sum(Z)
  n_ctrl <- N - n_trt

  # Degenerate case: all treated or all control
  if (n_trt == 0 || n_trt == N) {
    estimates[] <- NA
    ci_lower[]  <- NA
    ci_upper[]  <- NA
    return(list(estimates = estimates, ci_lower = ci_lower, ci_upper = ci_upper))
  }

  if (estimation_mode == "DIM") {
    # --- Difference in Means (vectorized across outcome resamples) ---
    mean_Y_trt  <- colMeans(Y_sim[Z == 1, , drop = FALSE])
    mean_Y_ctrl <- colMeans(Y_sim[Z == 0, , drop = FALSE])
    estimates <- mean_Y_trt - mean_Y_ctrl

    # Neyman variance estimator for each resample
    var_trt  <- apply(Y_sim[Z == 1, , drop = FALSE], 2, var)
    var_ctrl <- apply(Y_sim[Z == 0, , drop = FALSE], 2, var)
    se_dim <- sqrt(var_trt / n_trt + var_ctrl / n_ctrl)

    ci_lower <- estimates - qnorm(0.975) * se_dim
    ci_upper <- estimates + qnorm(0.975) * se_dim

  } else if (estimation_mode == "MLE") {
    # --- Maximum Likelihood via spatialreg::lagsarlm() ---
    for (k in seq_len(n_outcome_resamples)) {
      # Build data frame for this outcome resample
      if (include_spill_covariate) {
        df <- data.frame(Y = Y_sim[, k], Z = Z, Spill = spill_term,
                         X = X_matrix[, k])
        formula_str <- Y ~ Z + Spill + X
      } else {
        df <- data.frame(Y = Y_sim[, k], Z = Z, X = X_matrix[, k])
        formula_str <- Y ~ Z + X
      }

      fit <- tryCatch({
        suppressWarnings(
          spatialreg::lagsarlm(formula_str, data = df,
                               listw = active_listw, quiet = TRUE)
        )
      }, error = function(e) NULL)

      if (!is.null(fit) && "Z" %in% names(coef(fit))) {
        # Extract point estimate
        estimates[k] <- coef(fit)["Z"]

        # Extract SE from summary for CI
        fit_summary <- tryCatch(summary(fit), error = function(e) NULL)
        if (!is.null(fit_summary)) {
          se_z <- fit_summary$Coef["Z", "Std. Error"]
          ci_lower[k] <- estimates[k] - qnorm(0.975) * se_z
          ci_upper[k] <- estimates[k] + qnorm(0.975) * se_z
        } else {
          ci_lower[k] <- NA
          ci_upper[k] <- NA
        }
      } else {
        estimates[k] <- NA
        ci_lower[k]  <- NA
        ci_upper[k]  <- NA
      }
    }
  } else {
    stop("Unknown estimation_mode: ", estimation_mode,
         ". Must be 'DIM' or 'MLE'.")
  }

  list(estimates = estimates, ci_lower = ci_lower, ci_upper = ci_upper)
}

# Lean spatial-lag ML estimator (replicates lagsarlm(method = "eigen")) ----
#
# Added 2026-09-24 (simulation revision, docs/plans/simulation-revision-plan.md B2).
# Validated against spatialreg::lagsarlm() by code/tests/test_estimator_equivalence.R
# before anything depends on it. Each step mirrors the corresponding line of
# spatialreg 1.4.3's lagsarlm(); comments name the lagsarlm quantity reproduced.

#' Precompute the per-weights-matrix quantities for fit_sar_lag()
#'
#' The spatial-lag log-likelihood needs log|I - rho W| = sum_i log(1 - rho * lambda_i)
#' over the eigenvalues lambda_i of W. These depend only on W, so they are computed
#' once per neighbor type instead of once per fit (the main source of the speedup).
#' Follows lagsarlm's eigen_setup(): eigenvalues come from spatialreg::eigenw()
#' (the path lagsarlm takes when the listw cannot be symmetrized by similarity),
#' and the search interval for rho is 1 / range(Re(lambda)) shrunk by machine eps.
#'
#' @param listw listw object (row-standardized, style "W")
#' @return List with eig (possibly complex eigenvalues), interval (length-2 search
#'   interval for rho), W (dense N x N weights matrix), n (number of units)
#' @family sar_lag
#' @seealso [fit_sar_lag()]
#' @examples
#' # g <- build_spatial_grid(10); setup <- sar_lag_setup(g$listw_queen)
sar_lag_setup <- function(listw) {
  # lagsarlm uses a symmetric similar matrix when can.be.simmed(); our mat2listw()
  # objects never qualify, and supporting that branch would be untested code.
  if (spatialreg::can.be.simmed(listw)) {
    stop("sar_lag_setup(): symmetrizable listw not supported; use lagsarlm()")
  }
  eig <- spatialreg::eigenw(listw)
  # eig.range: only exactly-real eigenvalues bound the interval (as in eigen_setup)
  eig_real <- if (is.complex(eig)) Re(eig[which(Im(eig) == 0)]) else eig
  eig_range <- 1 / range(eig_real)
  list(
    eig      = eig,
    interval = c(eig_range[1] + .Machine$double.eps,
                 eig_range[2] - .Machine$double.eps),
    W        = spdep::listw2mat(listw),
    n        = length(listw$neighbours)
  )
}

#' Fit the spatial lag model y = rho W y + X b + e by maximum likelihood
#'
#' Same estimator as spatialreg::lagsarlm(method = "eigen"): rho maximizes the
#' concentrated log-likelihood
#'   l(rho) = log|I - rho W| - (n/2) log(2 pi) - (n/2) log(s2(rho)) - n/2,
#'   s2(rho) = SSE(rho) / n,  SSE(rho) = e_a - 2 rho e_b + rho^2 e_c,
#' where e_a = e0'e0, e_b = eW'e0, e_c = eW'eW for the OLS residuals e0 of y on X
#' and eW of Wy on X. Then b = OLS of (y - rho W y) on X, and standard errors come
#' from inverting the analytic information matrix for (sigma2, rho, b), scaled by
#' s2^2 (lagsarlm's "resvar").
#'
#' Warnings are raised with lagsarlm's exact messages ("Aliased variables found: ...",
#' "rho on interval bound - results should not be used"), so a caller's
#' withCallingHandlers() treats both engines identically. If the information matrix
#' cannot be inverted, lagsarlm falls back to a numerical Hessian; this function
#' instead warns and returns NA standard errors (se_ok = FALSE).
#'
#' @param y Numeric vector length n, outcome
#' @param x Numeric n x p model matrix with column names, including "(Intercept)"
#' @param setup List from sar_lag_setup() for the weights used to generate y
#' @param tol.opt Tolerance passed to optimize() (lagsarlm default)
#' @param tol.solve Tolerance passed to solve() for the information matrix
#'   (lagsarlm default)
#' @return List with coefficients (named, aliased columns dropped), rest.se (named
#'   SEs), rho, rho.se, LL (maximized log-likelihood), s2, aliased (named logical
#'   over the original columns), se_ok (logical)
#' @family sar_lag
#' @seealso [sar_lag_setup()], spatialreg::lagsarlm
#' @examples
#' # x <- cbind("(Intercept)" = 1, Z = Z, Spill = spill, X = X)
#' # fit <- fit_sar_lag(Y, x, sar_lag_setup(g$listw_queen)); fit$coefficients["Z"]
fit_sar_lag <- function(y, x, setup,
                        tol.opt   = .Machine$double.eps^0.5,
                        tol.solve = .Machine$double.eps) {
  n   <- setup$n
  W   <- setup$W
  eig <- setup$eig
  wy  <- as.vector(W %*% y)                        # lag.listw(listw, y)

  ## Aliasing: same LINPACK QR (tol 1e-7) as lm(y ~ x - 1) ----
  aliased <- is.na(stats::lm.fit(x, y)$coefficients)
  names(aliased) <- colnames(x)
  if (any(aliased)) {
    warning("Aliased variables found: ",
            paste(names(aliased)[aliased], collapse = " "))
    x <- x[, !aliased, drop = FALSE]
  }

  ## Concentrated log-likelihood in rho (sar.lag.mixed.f) ----
  e0  <- stats::lm.fit(x, y)$residuals             # e.null
  eW  <- stats::lm.fit(x, wy)$residuals            # e.w
  e_a <- sum(e0 * e0)
  e_b <- sum(eW * e0)
  e_c <- sum(eW * eW)
  ldet <- if (is.complex(eig)) {
    function(rho) Re(sum(log(1 - rho * eig)))      # eigen_ldet, complex branch
  } else {
    function(rho) sum(log(1 - rho * eig))
  }
  loglik <- function(rho) {
    SSE <- e_a - 2 * rho * e_b + rho * rho * e_c
    s2  <- SSE / n
    ldet(rho) - (n / 2) * log(2 * pi) - (n / 2) * log(s2) - (1 / (2 * s2)) * SSE
  }
  opt <- stats::optimize(loglik, interval = setup$interval, maximum = TRUE,
                         tol = tol.opt)
  rho <- opt$maximum
  if (isTRUE(all.equal(rho, setup$interval[1])) ||
      isTRUE(all.equal(rho, setup$interval[2]))) {
    warning("rho on interval bound - results should not be used")
  }

  ## Coefficients given rho-hat (lm.lag) ----
  lag_fit <- stats::lm.fit(x, y - rho * wy)
  b   <- lag_fit$coefficients
  r   <- lag_fit$residuals
  SSE <- sum(r * r)                                # deviance(lm.lag)
  s2  <- SSE / n

  ## Analytic information matrix for (sigma2, rho, b) (lagsarlm eigen branch) ----
  O     <- (eig / (1 - rho * eig))^2
  omega <- sum(O)
  if (is.complex(omega)) omega <- Re(omega)
  A     <- solve(diag(n) - rho * W)                # (I - rho W)^{-1}
  AW    <- A %*% W
  xb    <- x %*% b
  zero  <- rbind(rep(0, length(b)))
  xtawxb <- s2 * (t(x) %*% AW %*% xb)
  V     <- s2 * (s2 * sum(diag(crossprod(AW))) + crossprod(AW %*% xb)) +
           omega * s2^2
  inf1  <- rbind(n / 2, s2 * sum(diag(AW)), t(zero))
  inf2  <- rbind(s2 * sum(diag(AW)), V, xtawxb)
  inf3  <- rbind(zero, t(xtawxb), s2 * crossprod(x))
  inf   <- cbind(inf1, inf2, inf3)
  varb  <- try(solve(inf, tol = tol.solve), silent = TRUE)

  if (inherits(varb, "try-error")) {
    warning("inversion of asymptotic covariance matrix failed; SEs set to NA")
    rest_se <- stats::setNames(rep(NA_real_, length(b)), names(b))
    rho_se  <- NA_real_
    se_ok   <- FALSE
  } else {
    varb    <- (s2^2) * varb                       # lagsarlm "resvar"
    se_all  <- sqrt(diag(varb))
    rest_se <- stats::setNames(se_all[-c(1:2)], names(b))
    rho_se  <- unname(se_all[2])
    se_ok   <- TRUE
  }

  list(coefficients = b, rest.se = rest_se, rho = rho, rho.se = rho_se,
       LL = opt$objective, s2 = s2, aliased = aliased, se_ok = se_ok)
}

# Per-fit wrapper used by the 2026-09 runner ----

#' Fit one spatial-lag model, capturing warnings and errors instead of hiding them
#'
#' Replaces the old tryCatch(suppressWarnings(lagsarlm(...))) pattern (revision
#' M6/M7): every warning is recorded, and an error yields NA estimates plus an
#' "ERROR: ..." message rather than a silent NULL.
#'
#' @param y Numeric outcome vector
#' @param xm Model matrix with an "(Intercept)" column and a "Z" column
#' @param setup List from sar_lag_setup() (lean engine)
#' @param listw listw object (lagsarlm engine only)
#' @param engine "lean" (fit_sar_lag) or "lagsarlm"
#' @return List with tau (tau-hat), se (its ML standard error), warns (character
#'   vector of warning / error messages, possibly empty)
#' @family sar_lag
#' @seealso [fit_tau_models()], [fit_sar_lag()]
fit_one_lag_model <- function(y, xm, setup, listw = NULL,
                              engine = c("lean", "lagsarlm")) {
  engine <- match.arg(engine)
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      if (engine == "lean") {
        fit_sar_lag(y, xm, setup)
      } else {
        df <- data.frame(Y = y, xm[, colnames(xm) != "(Intercept)", drop = FALSE])
        rhs <- paste(setdiff(colnames(xm), "(Intercept)"), collapse = " + ")
        spatialreg::lagsarlm(stats::as.formula(paste("Y ~", rhs)), data = df,
                             listw = listw, quiet = TRUE)
      },
      error = function(e) {
        warns <<- c(warns, paste("ERROR:", conditionMessage(e)))
        NULL
      }
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(fit)) return(list(tau = NA_real_, se = NA_real_, warns = warns))
  # lean: coefficients / rest.se are named vectors; lagsarlm: same fields on the fit
  b  <- fit$coefficients   # same field name for both engines (lagsarlm: coef() adds rho)
  se <- fit$rest.se
  list(tau = if ("Z" %in% names(b)) unname(b[["Z"]]) else NA_real_,
       se  = if (!is.null(se) && "Z" %in% names(se)) unname(se[["Z"]]) else NA_real_,
       warns = warns)
}

#' Fit the oracle and non-oracle models to one simulated outcome
#'
#' Oracle (primary):      Y ~ Z + Spill + X  (Spill = true regime-specific S(Z))
#' Non-oracle (M8):       Y ~ Z + X
#' Both use the same Y, so the two estimators share every random draw.
#'
#' @param y Numeric outcome vector length N
#' @param Z Binary treatment vector length N
#' @param spill Spillover covariate S(Z), length N
#' @param X Incidence covariate for this fit's surface, length N
#' @param setup List from sar_lag_setup()
#' @param listw listw object (lagsarlm engine only)
#' @param engine "lean" or "lagsarlm"
#' @return Named list with elements oracle and nonoracle, each as returned by
#'   fit_one_lag_model()
#' @family sar_lag
#' @seealso [fit_one_lag_model()]
fit_tau_models <- function(y, Z, spill, X, setup, listw = NULL, engine = "lean") {
  x_oracle <- cbind("(Intercept)" = 1, Z = Z, Spill = spill, X = X)
  list(
    oracle    = fit_one_lag_model(y, x_oracle, setup, listw, engine),
    nonoracle = fit_one_lag_model(y, x_oracle[, c("(Intercept)", "Z", "X")],
                                  setup, listw, engine)
  )
}
