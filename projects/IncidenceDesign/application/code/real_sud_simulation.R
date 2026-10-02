# ============================================================
# Script: real_sud_simulation.R
# Purpose: Fit nested allocation/outcome experiments with verified MC precision.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud_setup.R helpers, parallel, digest
# ============================================================
.rs_sim_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(normalizePath(dirname(f)))
  }
  stop("Source this script by path")
})
source(file.path(.rs_sim_dir, "real_sud_setup.R"))

#' Return the approved grid and small verification profiles
#' @param profile smoke, pilot or production.
#' @return Explicit model, parameter and replication list.
#' @family real_sud
#' @seealso rs_units
#' @examples
#' rs_config("smoke")
rs_config <- function(profile = c("smoke", "pilot", "production")) {
  profile <- match.arg(profile)
  list(profile = profile, tau = 1, sigma = 1, years = as.character(2018:2021),
       rho = if (profile == "production") c(0, 0.2, 0.5) else c(0, 0.5),
       gamma = if (profile == "smoke") 0.5 else c(0.5, 0.8),
       regimes = c("both", "control_only"), designs = 1:9,
       J = switch(profile, smoke = 3L, pilot = 20L, production = 100L),
       R = switch(profile, smoke = 8L, pilot = 100L, production = 100L),
       singleton_R = switch(profile, smoke = 8L, pilot = 200L, production = 1000L),
       batch = 100L, J_tiers = c(100L, 200L, 400L),
       R_tiers = c(100L, 200L, 400L, 1000L, 2000L, 4000L),
       target_mse = 0.05, target_coverage = 0.01)
}

#' Form reporting blocks, including single-feature sensitivities
#' @param config Approved configuration.
#' @return Data frame of year, model, neighbor, summary, parameter and design keys.
#' @family real_sud
#' @seealso rs_support_key
#' @examples
#' nrow(rs_units(rs_config("production")))
rs_units <- function(config) {
  #' Construct one reporting-grid component
  #' @param model Outcome model label.
  #' @param nb Neighbor type.
  #' @param summary Regional summary rule.
  #' @param designs Design IDs.
  #' @param rho Spatial dependence values.
  #' @return Reporting-grid data frame.
  #' @examples
  #' # grid("education", "queen", "mean_rank", 1:9, c(0, 0.5))
  grid <- function(model, nb, summary, designs, rho) {
    expand.grid(Year = config$years, Model = model, Neighbor = nb, Summary = summary,
                Rho = rho, Gamma = config$gamma, Regime = config$regimes,
                Design_ID = designs, stringsAsFactors = FALSE)
  }
  u <- rbind(grid("education", "queen", "mean_rank", config$designs, config$rho),
             grid("baseline_sensitivity", "queen", "mean_rank", config$designs, config$rho),
             grid("education", "rook", "mean_rank", config$designs, c(0, 0.5)),
             grid("education", "queen", "mean_rate", 8L, config$rho),
             grid("education", "queen", "pooled_rate", 8L, config$rho))
  if (config$profile == "smoke") u <- u[u$Year == "2018", ]
  rownames(u) <- NULL
  u$Design <- unname(application_design_names(u$Design_ID))
  u$Design[u$Design_ID == 1] <- "Checkerboard (graph adaptation)"
  u
}

#' Get regional summary values for the requested allocation rule
#' @param a Named-order annual observation table.
#' @param regions Frozen region membership.
#' @param method mean_rank, mean_rate or pooled_rate.
#' @return Named numeric regional values; regional ties are retained.
#' @family real_sud
#' @seealso get_application_designs
#' @examples
#' # rs_region_values(a, regions, "mean_rank")
rs_region_values <- function(a, regions, method) {
  switch(method, mean_rank = tapply(a$rank, regions, mean),
         mean_rate = tapply(a$rate_per_100k, regions, mean),
         pooled_rate = tapply(a$deaths, regions, sum) / tapply(a$person_years, regions, sum) * 1e5,
         stop("Unknown regional summary"))
}

#' Canonicalize allocation support, proving cross-year/summary equivalence
#' @param u One reporting block.
#' @param setup Frozen design/input object.
#' @return Content-derived key; nuisance X is included in outcome keys separately.
#' @family real_sud
#' @seealso rs_block_key
#' @examples
#' # rs_support_key(u, setup)
rs_support_key <- function(u, setup) {
  a <- setup$annual[[u$Year]]; id <- u$Design_ID
  signal <- if (id %in% c(2, 6, 7)) a$rank else if (id == 5) setup$blocks[[u$Year]] else NULL
  if (id == 8) {
    v <- rs_region_values(a, setup$regions, u$Summary)
    # Only order and tie groups determine region saturation, not numeric gaps.
    signal <- rank(v, ties.method = "average")
  }
  digest::digest(list(design = id, nb = u$Neighbor,
                      regions = if (id %in% c(3, 8)) setup$regions else NULL,
                      signal = signal, ids = setup$ids))
}

#' Canonical key for equivalent outcome/analysis distributions
#' @param u One reporting block.
#' @param setup Frozen inputs.
#' @param config Parameters, excluding replication counts.
#' @return Stable content hash; reused blocks retain explicit provenance.
#' @family real_sud
#' @seealso rs_support_key
#' @examples
#' # rs_block_key(u, setup, config)
rs_block_key <- function(u, setup, config) {
  digest::digest(list(support = rs_support_key(u, setup), model = u$Model,
    X = if (u$Model == "baseline_sensitivity") setup$annual[[u$Year]]$X else NULL,
    W = setup$weights[[u$Neighbor]]$W, rho = u$Rho, gamma = u$Gamma,
    regime = u$Regime, tau = config$tau, sigma = config$sigma))
}

#' Draw assignments from frozen supports with extensible per-draw seeds
#' @param u Reporting block.
#' @param setup Frozen inputs.
#' @param J Number of allocation draws.
#' @return Integer 58 by J matrix; no geography is recomputed.
#' @family real_sud
#' @seealso rs_seed
#' @examples
#' # rs_draws(u, setup, 100)
rs_draws <- function(u, setup, J) {
  a <- setup$annual[[u$Year]]; key <- rs_support_key(u, setup)
  vapply(seq_len(J), function(j) {
    rs_seed("Z", key, j)
    get_application_designs(u$Design_ID, 1L, a$X, setup$weights[[u$Neighbor]]$nb,
      setup$coords, a$person_years, a$n_counties, setup$regions,
      block_id = setup$blocks[[u$Year]],
      region_summary = rs_region_values(a, setup$regions, u$Summary))[, 1]
  }, integer(58))
}

#' Construct the agreed SAR mean and actual fitted matrix
#' @param u Reporting block.
#' @param z Binary treatment vector.
#' @param setup Frozen inputs.
#' @param config Known effect/noise settings.
#' @return Mean, multiplier, fitted matrix, identification and conditioning diagnostics.
#' @family real_sud
#' @seealso fit_one_lag_model
#' @examples
#' # rs_model(u, z, setup, config)
rs_model <- function(u, z, setup, config) {
  W <- setup$weights[[u$Neighbor]]$W
  wz <- as.vector(W %*% z)
  spill <- u$Gamma * wz * if (u$Regime == "control_only") (1 - z) else 1
  xm <- cbind("(Intercept)" = 1, Z = z, Spill = spill)
  betaX <- rep(0, length(z))
  if (u$Model == "baseline_sensitivity") {
    X <- setup$annual[[u$Year]]$X
    xm <- cbind(xm, X = X); betaX <- X
  }
  nuisance <- xm[, colnames(xm) != "Z", drop = FALSE]
  residual_Z <- stats::lm.fit(nuisance, z)$residuals
  identified <- sum(residual_Z^2) > 1e-10
  norms <- sqrt(colSums(xm^2))
  scaled <- sweep(xm, 2, pmax(norms, .Machine$double.eps), "/")
  A <- solve(diag(length(z)) - u$Rho * W)
  list(mean = as.vector(A %*% (config$tau * z + spill + betaX)), A = A,
       xm = xm, Identified = identified, Rank = qr(xm)$rank,
       Residual_Z_SS = sum(residual_Z^2), Condition = kappa(scaled, exact = TRUE),
       Treatment_Spill_Correlation = if (sd(spill) > 0) cor(z, spill) else NA_real_)
}

#' Summarize independent outcome fits for one fixed allocation
#' @param fits Data frame of errors, SEs, identification and warning statuses.
#' @return Allocation metrics, including independent-fit counts and MC variances.
#' @family real_sud
#' @seealso allocation_fit_summary, rs_distribution
#' @examples
#' # rs_allocation_metrics(fits)
rs_allocation_metrics <- function(fits) {
  invalid <- grepl("interval bound|ERROR:", fits$messages)
  # Preserve raw fit values in the cache; invalid boundary fits cannot contribute
  # to performance or a claimed complete comparison merely because they are finite.
  fits$error[invalid] <- fits$se[invalid] <- NA_real_
  a <- allocation_fit_summary(fits, tau = 1)
  a$N_Bound <- sum(grepl("interval bound", fits$messages))
  a$Error_Second_Moment <- if (a$Complete) mean(fits$error^2) else NA_real_
  a$Covered <- if (a$Complete) mean(abs(fits$error) <= qnorm(0.975) * fits$se) else NA_real_
  a
}

#' Compute joint uncertainty with cache covariance, not expanded-fit independence
#' @param a Unique allocation metrics including Frequency.
#' @param singleton Whether eligible support is proven singleton.
#' @return Frequency-weighted performance and allocation-risk diagnostics.
#' @family real_sud
#' @seealso allocation_distribution
#' @examples
#' # rs_distribution(a, FALSE)
rs_distribution <- function(a, singleton) {
  s <- allocation_distribution(a)
  s$N_Independent_Fits <- sum(a$N_Outcome)
  s$N_Valid_Est <- sum(a$N_Valid_Est)
  s$N_Valid_CI <- sum(a$N_Valid_CI)
  s$N_Aliased <- sum(a$N_Aliased); s$N_Warn <- sum(a$N_Warn)
  s$N_Bound <- sum(a$N_Bound)
  s$Proven_Singleton <- singleton
  s$Bias <- s$MCSE_Bias <- s$Coverage <- s$MCSE_Coverage <- s$SD <- s$Relative_MCSE_MSE <- NA_real_
  s$Coverage_Variance_Allocation <- s$Coverage_Variance_Noise <- NA_real_
  if (!isTRUE(s$Complete)) return(s)
  w <- a$Frequency / sum(a$Frequency); J <- sum(a$Frequency)
  for (metric in c("Bias", "Coverage")) {
    m <- a[[metric]]; v <- a[[paste0("SE_", metric)]]^2
    raw <- if (J > 1) sum(a$Frequency * (m - sum(w * m))^2) / (J - 1) else 0
    # Equivalent to corrected allocation variance/J + sum(w_i² v_i).
    joint <- if (J == 1) v else raw / J + sum(a$Frequency * (a$Frequency - 1) * v) / (J * (J - 1))
    s[[metric]] <- sum(w * m)
    s[[paste0("MCSE_", metric)]] <- sqrt(joint)
    if (metric == "Coverage") {
      noise_sample <- if (J > 1) sum(a$Frequency * (1 - w) * v) / (J - 1) else 0
      s$Coverage_Variance_Allocation <- raw - noise_sample
      s$Coverage_Variance_Noise <- sum(w^2 * v)
    }
  }
  s$SD <- sqrt(max(0, sum(w * a$MSE) - s$Bias^2))
  s$Relative_MCSE_MSE <- s$SE_Mean_MSE_Joint / s$Mean_MSE
  s
}

#' Simulate or extend one canonical block, preserving independent noise prefixes
#' @param u Reporting block.
#' @param setup Frozen inputs.
#' @param config Statistical and replication settings.
#' @param J Allocation draws.
#' @param R Outcomes per unique allocation.
#' @param old Previous matching block result, including fit cache.
#' @return Allocation/draw/fit caches, performance and identification diagnostics.
#' @family real_sud
#' @seealso rs_distribution
#' @examples
#' # rs_block(u, setup, config, 100, 100)
rs_block <- function(u, setup, config, J, R, old = NULL) {
  key <- rs_block_key(u, setup, config)
  Z <- rs_draws(u, setup, J)
  hashes <- apply(Z, 2, digest::digest, algo = "xxhash64")
  unique_hashes <- unique(hashes)
  a <- setup$annual[[u$Year]]
  singleton <- u$Design_ID == 1 || (u$Design_ID == 2 &&
    length(unique(a$rate_per_100k)) == nrow(a))
  lw <- spdep::mat2listw(setup$weights[[u$Neighbor]]$W, style = "W")
  engine_setup <- sar_lag_setup(lw)
  fits_out <- list(); metrics <- list()
  for (h in unique_hashes) {
    z <- Z[, match(h, hashes)]; model <- rs_model(u, z, setup, config)
    fits <- if (!is.null(old)) old$fits[[h]] else NULL
    if (!is.null(fits)) fits <- fits[fits$replicate <= R, , drop = FALSE]
    previous <- if (is.null(fits)) 0L else nrow(fits)
    new <- list()
    if (previous < R) {
      batches <- seq.int(previous %/% config$batch + 1L, ceiling(R / config$batch))
      for (b in batches) {
        rs_seed("eps", key, h, b)
        E <- matrix(rnorm(58 * config$batch, 0, config$sigma), nrow = 58)
        ix <- (b - 1L) * config$batch + seq_len(config$batch)
        keep <- which(ix > previous & ix <= R)
        Y <- model$mean + model$A %*% E[, keep, drop = FALSE]
        for (k in seq_along(keep)) {
          fit <- fit_one_lag_model(Y[, k], model$xm, engine_setup)
          messages <- paste(fit$warns, collapse = " | ")
          new[[length(new) + 1L]] <- data.frame(replicate = ix[keep[k]],
            tau = fit$tau, error = if (model$Identified) fit$tau - config$tau else NA_real_,
            se = if (model$Identified) fit$se else NA_real_,
            alias = any(grepl("Aliased", fit$warns)), messages = messages,
            identified = model$Identified, stringsAsFactors = FALSE)
        }
      }
    }
    if (length(new)) fits <- rbind(fits, do.call(rbind, new))
    stopifnot(nrow(fits) == R, identical(as.integer(fits$replicate), seq_len(R)))
    fits_out[[h]] <- fits
    m <- rs_allocation_metrics(fits)
    m$Allocation_ID <- h; m$Frequency <- sum(hashes == h)
    m$Identified <- model$Identified; m$Rank <- model$Rank
    m$Residual_Z_SS <- model$Residual_Z_SS; m$Condition <- model$Condition
    m$Treatment_Spill_Correlation <- model$Treatment_Spill_Correlation
    m$Treated <- sum(z); m$Population_Share <- sum(a$person_years * z) / sum(a$person_years)
    m$Incidence_Difference <- mean(a$rate_per_100k[z == 1]) - mean(a$rate_per_100k[z == 0])
    metrics[[h]] <- m
  }
  allocations <- do.call(rbind, metrics); rownames(allocations) <- NULL
  s <- rs_distribution(allocations, singleton)
  s$J <- J; s$R <- R; s$Source_ID <- key
  s$Precision_OK <- isTRUE(s$Complete) && is.finite(s$Relative_MCSE_MSE) &&
    s$Relative_MCSE_MSE <= config$target_mse && s$MCSE_Coverage <= config$target_coverage
  list(summary = s, allocations = allocations, draws = hashes, Z = Z,
       fits = fits_out, unit = u, J = J, R = R, key = key)
}

#' Run fixed refinement tiers according to the dominant uncertainty component
#' @param u Canonical block.
#' @param setup Frozen inputs.
#' @param config Profile configuration.
#' @param checkpoint Matching block cache file.
#' @return Verified block or explicit incomplete/precision-unresolved result.
#' @family real_sud
#' @seealso rs_block
#' @examples
#' # rs_run_block(u, setup, config, checkpoint)
rs_run_block <- function(u, setup, config, checkpoint) {
  old <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
  if (!is.null(old) && !identical(old$key, rs_block_key(u, setup, config))) stop("Block key mismatch")
  a <- setup$annual[[u$Year]]
  singleton <- u$Design_ID == 1 || (u$Design_ID == 2 && length(unique(a$rank)) == 58)
  J <- config$J; R <- if (singleton) config$singleton_R else config$R
  if (!is.null(old)) { J <- old$J; R <- old$R }
  repeat {
    out <- rs_block(u, setup, config, J, R, old)
    saveRDS(out, checkpoint)
    if (config$profile != "production" || !isTRUE(out$summary$Complete) || isTRUE(out$summary$Precision_OK)) break
    next_J <- config$J_tiers[config$J_tiers > J]
    next_R <- config$R_tiers[config$R_tiers > R]
    s <- out$summary
    # Choose the dominant component of the target with the largest relative miss.
    coverage_drives <- s$MCSE_Coverage / config$target_coverage >
      s$Relative_MCSE_MSE / config$target_mse
    allocation_var <- max(0, if (coverage_drives) s$Coverage_Variance_Allocation else s$Variance_Corrected) / J
    outcome_var <- if (coverage_drives) s$Coverage_Variance_Noise else s$SE_Mean_MSE_Noise^2
    increase_J <- !singleton && length(next_J) && (allocation_var >= outcome_var || !length(next_R))
    if (increase_J) J <- next_J[1] else if (length(next_R)) R <- next_R[1] else if (!singleton && length(next_J)) J <- next_J[1] else break
    old <- out
  }
  out
}
