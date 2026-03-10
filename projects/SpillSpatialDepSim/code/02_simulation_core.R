# 02_simulation_core.R
# Simulation engine for the SpillSpatialDepSim unified simulation.
# Provides compute_metrics() and run_scenario() — the inner double loop
# (combo_index x dataset_num) extracted from SpatialSim_NC_DOC.Rmd and
# generalized across grid sizes and spillover types.
#
# Depends on: 01_grid_setup.R (grid object with rook_neighbors, assign_district, etc.)

library(spdep)
library(spatialreg)


#' Compute bias, variance, and MSE for a matrix of estimates.
#'
#' Identical to compute_metrics() in the original simulation scripts.
#'
#' @param estimates  Numeric matrix: rows = iterations, cols = combinations.
#' @param true_value Scalar: the true parameter value (e.g. beta=1.0).
#' @return Data frame (n_combos rows) with columns: bias, bias_abs, variance, MSE.
compute_metrics <- function(estimates, true_value) {
  as.data.frame(t(apply(estimates, 2, function(col) {
    bias     <- mean(col)    - true_value
    bias_abs <- abs(bias)
    variance <- var(col)
    mse      <- mean((col - true_value)^2)
    c(bias = bias, bias_abs = bias_abs, variance = variance, MSE = mse)
  })))
}


#' Run one full simulation scenario (all combos x N datasets).
#'
#' Replicates the double loop from SpatialSim_NC_DOC.Rmd and its 3x3/3x4
#' variants, parameterized over grid size and spillover type.
#'
#' @param grid          Output of setup_grid() — grid metadata + rook_neighbors
#'                      + assign_district() function.
#' @param combinations  List of integer vectors — all C(n_districts, n_trt) combos.
#' @param valid_block_idx Integer vector: indices in `combinations` that satisfy
#'                        full block stratification (output of build_valid_block_indices()).
#' @param psi           True spillover effect (0.5, 0.6, 0.7, or 0.8).
#' @param rho           True spatial autocorrelation (0.00 or 0.01).
#' @param trt_spill     Logical. FALSE = TrtNoSpill (case 1: spillover to control
#'                      neighbors only); TRUE = TrtSpill (case 2: spillover to
#'                      all neighbors of treated districts).
#' @param N             Number of datasets to simulate per combination.
#' @param alpha         True intercept (default 0.2).
#' @param beta          True intervention effect (default 1.0).
#' @param sd            Residual standard deviation (default 0.1).
#' @return Data frame with one row per (combo_id x parameter), columns:
#'   combo_id, is_block, param, bias, bias_abs, variance, MSE.
run_scenario <- function(grid, combinations, valid_block_idx,
                         psi, rho, trt_spill, N,
                         alpha = 0.2, beta = 1.0, sd = 0.1) {

  n_combos <- length(combinations)
  n_pts    <- grid$points_per_iteration
  L        <- grid$cell_length
  rn       <- grid$rook_neighbors
  n_dist   <- grid$n_districts

  # Storage matrices: rows = iterations, cols = combos
  alpha_est <- matrix(NA_real_, nrow = N, ncol = n_combos)
  beta_est  <- matrix(NA_real_, nrow = N, ncol = n_combos)
  psi_est   <- matrix(NA_real_, nrow = N, ncol = n_combos)
  rho_est   <- matrix(NA_real_, nrow = N, ncol = n_combos)

  for (combo_index in seq_len(n_combos)) {
    current_combination <- combinations[[combo_index]]

    for (dataset_num in seq_len(N)) {

      # --- Generate random points on grid ---
      points <- data.frame(
        x = sample(seq(0, grid$grid_width  - 1), n_pts, replace = TRUE),
        y = sample(seq(0, grid$grid_height - 1), n_pts, replace = TRUE)
      )

      # --- Assign districts (row-major) ---
      points$district <- grid$assign_district(points$x, points$y)

      # --- Assign treatment status ---
      points$treatment    <- ifelse(points$district %in% current_combination, "intervention", "control")
      points$intervention <- as.integer(points$treatment == "intervention")
      points$spillover    <- 0L

      # --- Add per-district neighbor / treatment columns ---
      for (i in seq_len(n_dist)) {
        points[[paste0("nb.", i)]]  <- as.integer(points$district %in% rn[[as.character(i)]])
        points[[paste0("trt.", i)]] <- as.integer(i %in% current_combination)
      }

      # --- Assign spillover ---
      if (!trt_spill) {
        # Case 1 (TrtNoSpill): spills only to control neighbors of treated districts
        for (i in seq_len(n_dist)) {
          nb_col  <- paste0("nb.", i)
          trt_col <- paste0("trt.", i)
          points$spillover <- ifelse(
            points[[nb_col]] == 1L & points[[trt_col]] == 1L & points$treatment == "control",
            1L, points$spillover
          )
        }
      } else {
        # Case 2 (TrtSpill): spills to all neighbors of treated districts
        for (i in seq_len(n_dist)) {
          nb_col  <- paste0("nb.", i)
          trt_col <- paste0("trt.", i)
          points$spillover <- ifelse(
            points[[nb_col]] == 1L & points[[trt_col]] == 1L,
            1L, points$spillover
          )
        }
      }

      # Remove helper columns
      drop_cols   <- c(grep("^nb\\.", names(points), value = TRUE),
                       grep("^trt\\.", names(points), value = TRUE))
      points      <- points[, setdiff(names(points), drop_cols)]

      # --- Build point-level W matrix (n_pts x n_pts) ---
      W <- matrix(0L, nrow = n_pts, ncol = n_pts)
      for (i in seq_len(n_pts)) {
        d_i      <- points$district[i]
        nb_dists <- rn[[as.character(d_i)]]
        for (nd in nb_dists) {
          nb_pts <- which(points$district == nd)
          W[i, nb_pts] <- 1L
        }
      }

      # --- Generate response via spatial lag model ---
      I_mat          <- diag(n_pts)
      epsilon        <- rnorm(n_pts, mean = 0, sd = sd)
      linear_resp    <- alpha + beta * points$intervention + psi * points$spillover + epsilon
      A              <- I_mat - rho * W
      points$response <- solve(A, linear_resp)

      # --- Fit spatial lag model ---
      listw <- mat2listw(W, style = "W")
      model <- tryCatch(
        lagsarlm(response ~ intervention + spillover, data = points, listw = listw),
        error = function(e) NULL
      )

      if (!is.null(model)) {
        rho_est[dataset_num, combo_index]   <- coef(model)[1]  # rho
        alpha_est[dataset_num, combo_index] <- coef(model)[2]  # intercept
        beta_est[dataset_num, combo_index]  <- coef(model)[3]  # intervention
        psi_est[dataset_num, combo_index]   <- coef(model)[4]  # spillover
      }
      # NA stays in matrix on convergence failure — compute_metrics() uses mean/var which propagate NA
    }
  }

  # --- Compute metrics per combo for each parameter ---
  metrics_alpha <- compute_metrics(alpha_est, true_value = alpha)
  metrics_beta  <- compute_metrics(beta_est,  true_value = beta)
  metrics_psi   <- compute_metrics(psi_est,   true_value = psi)
  metrics_rho   <- compute_metrics(rho_est,   true_value = rho)

  is_block <- seq_len(n_combos) %in% valid_block_idx

  make_param_df <- function(metrics, param_name) {
    data.frame(
      combo_id  = seq_len(n_combos),
      is_block  = is_block,
      param     = param_name,
      bias      = metrics$bias,
      bias_abs  = metrics$bias_abs,
      variance  = metrics$variance,
      MSE       = metrics$MSE,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_param_df(metrics_alpha, "alpha"),
    make_param_df(metrics_beta,  "beta"),
    make_param_df(metrics_psi,   "psi"),
    make_param_df(metrics_rho,   "rho")
  )
}
