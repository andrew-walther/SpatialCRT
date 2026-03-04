# ==============================================================================
# 02_incidence_generation.R
# Three incidence generation modes: iid, spatial (SAR), Poisson count-based
# All modes return N x n_resamples matrices with values in [0, 1]
# ==============================================================================

#' Mode 1: iid Uniform(0,1) incidence (naive baseline)
#'
#' @param N Integer, number of clusters
#' @param n_resamples Integer, number of incidence resamples (columns)
#' @return N x n_resamples matrix with values in [0, 1]
generate_incidence_iid <- function(N, n_resamples) {
  matrix(runif(N * n_resamples, 0, 1), nrow = N, ncol = n_resamples)
}

#' Mode 2: Spatially correlated incidence via SAR filter
#'
#' X* = (I - rho_X * W)^{-1} * u,  u ~ N(0, I)
#' X  = pnorm(X*)  ->  marginally Uniform(0,1) but spatially correlated
#'
#' @param N Integer, number of clusters
#' @param n_resamples Integer, number of incidence resamples
#' @param W N x N row-standardized spatial weight matrix
#' @param rho_incidence Numeric in [0, 1), spatial correlation for incidence
#' @return N x n_resamples matrix with values in [0, 1]
generate_incidence_spatial <- function(N, n_resamples, W, rho_incidence) {
  I_mat <- diag(N)
  inv_filter <- solve(I_mat - rho_incidence * W)

  # Generate iid standard normal innovations
  U <- matrix(rnorm(N * n_resamples), nrow = N, ncol = n_resamples)

  # Apply SAR filter -> spatially correlated normals
  X_star <- inv_filter %*% U

  # Probability integral transform -> marginally Uniform(0,1)
  pnorm(X_star)
}

#' Mode 3: Poisson count-based incidence with spatially correlated rates
#'
#' 1. Spatially correlated log-rates via SAR filter
#' 2. Poisson counts: C ~ Poisson(lambda * P)
#' 3. Observed rates: R = C / P
#' 4. Rank-normalize to [0, 1]: (rank - 0.5) / N
#'
#' @param N Integer, number of clusters
#' @param n_resamples Integer, number of incidence resamples
#' @param W N x N row-standardized spatial weight matrix
#' @param rho_incidence Numeric in [0, 1), spatial correlation for rates
#' @param base_rate Numeric, baseline incidence rate (e.g., 35/100000 for SUD)
#' @param pop_per_cluster Numeric, mean cluster population (default 1000)
#' @param pop_mode Character, "equal" or "heterogeneous"
#' @return N x n_resamples matrix with values in [0, 1]
generate_incidence_poisson <- function(N, n_resamples, W, rho_incidence,
                                       base_rate = 35 / 100000,
                                       pop_per_cluster = 1000,
                                       pop_mode = "equal") {
  I_mat <- diag(N)
  inv_filter <- solve(I_mat - rho_incidence * W)
  mu_0 <- log(base_rate)

  # Generate cluster populations
  if (pop_mode == "equal") {
    P <- rep(pop_per_cluster, N)
  } else {
    # LogNormal with mean ~ pop_per_cluster, moderate heterogeneity
    P <- round(rlnorm(N, meanlog = log(pop_per_cluster), sdlog = 0.5))
    P <- pmax(P, 100)  # floor at 100 to avoid degenerate small populations
  }

  # Output matrix
  X_matrix <- matrix(0, nrow = N, ncol = n_resamples)

  for (j in seq_len(n_resamples)) {
    # Spatially correlated log-rates
    u <- rnorm(N)
    log_lambda <- mu_0 + as.vector(inv_filter %*% u)
    lambda <- exp(log_lambda)

    # Poisson counts
    counts <- rpois(N, lambda = lambda * P)

    # Observed rates
    rates <- counts / P

    # Rank-normalize to [0, 1]
    X_matrix[, j] <- (rank(rates, ties.method = "average") - 0.5) / N
  }

  X_matrix
}

#' Unified dispatcher for incidence generation
#'
#' @param mode Character, one of "iid", "spatial", "poisson"
#' @param N Integer, number of clusters
#' @param n_resamples Integer, number of incidence resamples
#' @param W N x N spatial weight matrix (required for spatial and poisson modes)
#' @param rho_incidence Numeric, spatial correlation for incidence (ignored for iid)
#' @param ... Additional arguments passed to mode-specific generators
#' @return N x n_resamples matrix with values in [0, 1]
generate_incidence <- function(mode, N, n_resamples, W = NULL,
                               rho_incidence = 0, ...) {
  switch(mode,
    "iid"     = generate_incidence_iid(N, n_resamples),
    "spatial" = generate_incidence_spatial(N, n_resamples, W, rho_incidence),
    "poisson" = generate_incidence_poisson(N, n_resamples, W, rho_incidence, ...),
    stop("Unknown incidence mode: ", mode,
         ". Must be one of: iid, spatial, poisson")
  )
}
