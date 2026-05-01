# ==============================================================================
# Synthetic SUD incidence generation for the NC Community College application
# ==============================================================================

#' Build Queen-contiguity spatial weights for irregular application clusters
#'
#' @description
#' Converts an `sf` polygon object into a spatial econometric weights object.
#' The default is the standard application setup: Queen contiguity, zero diagonal,
#' and row-standardized weights.
#'
#' @param clusters_sf An `sf` object with one row per spatial cluster.
#' @param queen Logical. If TRUE, use Queen contiguity; if FALSE, use Rook.
#' @param diagonal_weight Numeric self-weight to place on the diagonal before
#'   row standardization. The default, 0, gives the usual `W_ii = 0`.
#' @param zero_policy Logical passed to `spdep` for possible isolated clusters.
#'
#' @return A list with `nb`, `W`, `listw`, `row_sums`, and `zero_policy`.
#' @export
build_application_spatial_weights <- function(clusters_sf,
                                              queen = TRUE,
                                              diagonal_weight = 0,
                                              zero_policy = TRUE) {
  if (!inherits(clusters_sf, "sf")) {
    stop("`clusters_sf` must be an sf object.", call. = FALSE)
  }

  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package `spdep` is required.", call. = FALSE)
  }

  # `poly2nb()` constructs the binary neighborhood graph from polygon borders.
  # Queen contiguity treats shared boundaries and shared vertices as adjacency.
  nb <- spdep::poly2nb(clusters_sf, queen = queen)

  # `style = "B"` gives the unstandardized binary adjacency matrix A, where
  # A_ij = 1 if j is a neighbor of i and A_ij = 0 otherwise.
  adjacency <- spdep::nb2mat(nb, style = "B", zero.policy = zero_policy)

  # Spatial econometric SAR/SDM models conventionally set W_ii = 0. The
  # parameter is exposed because sensitivity analyses sometimes include a
  # small self-weight before row normalization.
  diag(adjacency) <- diagonal_weight

  # Row sums define each cluster's total outgoing neighbor mass before
  # standardization. Zero rows are left as zeros under `zero_policy = TRUE`.
  row_sums <- rowSums(adjacency)

  # Row standardization creates W, where each non-isolated row sums to one.
  # This makes WZ and WY interpretable as neighbor averages.
  W <- adjacency
  nonzero_rows <- row_sums > 0
  W[nonzero_rows, ] <- W[nonzero_rows, , drop = FALSE] / row_sums[nonzero_rows]

  # `listw` is the object expected by spatialreg::lagsarlm().
  listw <- spdep::mat2listw(W, style = "W", zero.policy = zero_policy)

  list(
    nb = nb,
    W = W,
    listw = listw,
    row_sums = rowSums(W),
    zero_policy = zero_policy
  )
}

#' Simulate a latent SAR log-risk surface
#'
#' @description
#' Generates `eta = (I - rho W)^-1 (log(base_rate) + epsilon)`. The inverse
#' spatial multiplier propagates shocks through neighboring clusters and creates
#' spatially autocorrelated baseline incidence risk.
#'
#' @param W Row-standardized spatial weights matrix.
#' @param base_rate Baseline event rate per person.
#' @param rho_incidence Spatial autoregressive coefficient.
#' @param sigma Standard deviation of the Gaussian innovation.
#'
#' @return Numeric vector of cluster-specific log risks.
#' @export
simulate_sar_log_risk <- function(W,
                                  base_rate = 35 / 100000,
                                  rho_incidence = 0.20,
                                  sigma = 0.35) {
  N <- nrow(W)
  if (ncol(W) != N) {
    stop("`W` must be a square matrix.", call. = FALSE)
  }
  if (abs(rho_incidence) >= 1) {
    stop("`rho_incidence` should be inside (-1, 1) for a row-standardized W.", call. = FALSE)
  }
  if (base_rate <= 0) {
    stop("`base_rate` must be positive.", call. = FALSE)
  }

  # Identity matrix for the SAR spatial multiplier.
  I_mat <- diag(N)

  # Independent innovation term on the log-risk scale.
  epsilon <- stats::rnorm(N, mean = 0, sd = sigma)

  # The SAR multiplier `(I - rho W)^-1` injects spatial dependence. A positive
  # innovation in one cluster propagates to its Queen-contiguous neighbors and
  # then recursively through the spatial graph.
  as.vector(solve(I_mat - rho_incidence * W, log(base_rate) + epsilon))
}

#' Generate placeholder synthetic SUD incidence for application clusters
#'
#' @description
#' Appends synthetic SUD incidence fields to dissolved Community College service
#' areas. This is a placeholder generator; future true incidence data should use
#' the same output schema.
#'
#' @param clusters_sf An `sf` object with one row per application cluster.
#' @param cluster_id_col Optional name of an existing cluster ID column.
#' @param population_col Optional name of an existing population column.
#' @param base_rate Baseline event rate per person.
#' @param rho_incidence SAR coefficient for the latent log-risk surface.
#' @param sigma Standard deviation of SAR innovations.
#' @param queen Logical. If TRUE, use Queen contiguity.
#' @param diagonal_weight Numeric self-weight for W before row standardization.
#' @param seed Optional integer seed.
#'
#' @return The input `sf` object with synthetic incidence columns appended.
#' @export
generate_synthetic_sud_data <- function(clusters_sf,
                                        cluster_id_col = NULL,
                                        population_col = NULL,
                                        base_rate = 35 / 100000,
                                        rho_incidence = 0.20,
                                        sigma = 0.35,
                                        queen = TRUE,
                                        diagonal_weight = 0,
                                        seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  weights <- build_application_spatial_weights(
    clusters_sf = clusters_sf,
    queen = queen,
    diagonal_weight = diagonal_weight
  )

  N <- nrow(clusters_sf)
  cluster_id <- if (!is.null(cluster_id_col) && cluster_id_col %in% names(clusters_sf)) {
    clusters_sf[[cluster_id_col]]
  } else {
    seq_len(N)
  }

  population <- if (!is.null(population_col) && population_col %in% names(clusters_sf)) {
    clusters_sf[[population_col]]
  } else {
    rep(100000, N)
  }

  if (any(!is.finite(population)) || any(population <= 0)) {
    stop("Population values must be positive and finite.", call. = FALSE)
  }

  # Generate a spatially autocorrelated log-risk surface.
  log_risk <- simulate_sar_log_risk(
    W = weights$W,
    base_rate = base_rate,
    rho_incidence = rho_incidence,
    sigma = sigma
  )

  # Convert log-risk to positive rates and then to expected counts.
  risk_per_person <- exp(log_risk)
  expected_count <- risk_per_person * population

  # Poisson sampling maps the continuous risk surface to realistic event counts.
  sud_count <- stats::rpois(N, lambda = expected_count)
  sud_rate_per_100k <- (sud_count / population) * 100000

  # Rank-normalized incidence is the design-facing covariate used by the
  # original simulation designs.
  incidence_rank01 <- rank(sud_rate_per_100k, ties.method = "average") / N

  clusters_sf$cluster_id <- cluster_id
  clusters_sf$population <- as.numeric(population)
  clusters_sf$sud_count <- as.integer(sud_count)
  clusters_sf$sud_rate_per_100k <- as.numeric(sud_rate_per_100k)
  clusters_sf$incidence_rank01 <- as.numeric(incidence_rank01)
  clusters_sf$synthetic_log_risk <- as.numeric(log_risk)
  clusters_sf$synthetic_expected_count <- as.numeric(expected_count)
  clusters_sf$incidence_source <- "synthetic_poisson_sar_placeholder"

  clusters_sf
}

