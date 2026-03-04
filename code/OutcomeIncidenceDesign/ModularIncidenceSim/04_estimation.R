# ==============================================================================
# 04_estimation.R
# DIM and MLE estimation of the treatment effect tau, with CI extraction
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
