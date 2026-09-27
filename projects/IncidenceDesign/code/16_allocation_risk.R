# ============================================================
# Script: 16_allocation_risk.R
# Purpose: Estimate conditional allocation MSE and its distribution at fixed X.
# Author: Andrew Walther
# Created: 2026-09-27
# Dependencies: spdep, spatialreg, dplyr, digest, parallel; modules 01-04
# ============================================================
# Run from code/: VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# OMP_NUM_THREADS=1 Rscript 16_allocation_risk.R pilot 8
# Profiles: smoke, benchmark, pilot. All outputs/checkpoints are isolated under
# results/allocation_risk/<profile>/. No existing simulation result is overwritten.
# Local fork workers use one BLAS thread; each worker holds only one fixed-X block
# (100 clusters, at most 100 allocation x 100 outcome fits). Not a SLURM runner.

# Load unchanged model modules ----
if (!exists("allocation_script_dir")) {
  allocation_script_dir <- local({
    f <- grep("--file=", commandArgs(FALSE), value = TRUE)
    if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else
      normalizePath(dirname(sys.frame(1)$ofile))
  })
}
for (f in c("01_spatial_setup.R", "02_incidence_generation.R", "03_designs.R",
            "04_estimation.R")) source(file.path(allocation_script_dir, f))

# Seeds and profile ----
#' Seed one allocation-risk draw using an explicit, stable key
#'
#' X deliberately uses the existing main-run key; Z and eps use a new namespace.
#' Numeric fields match the main runner's two-decimal convention. An assignment
#' draw has its own index, and noise has fixed-size batch keys, so extending J or R
#' never changes any earlier assignment or outcome. Distinct assignments/designs
#' have independent noise streams (not common random numbers).
#' @param ... Character/numeric key fields.
#' @return Invisibly, the integer RNG seed.
#' @examples
#' allocation_seed("allocation-risk-v1", "Z", "iid", 0, "queen", 0.5, 0.8, "both", 1, 9, 1)
allocation_seed <- function(...) {
  fields <- lapply(list(...), function(x) if (is.numeric(x)) sprintf("%.2f", x) else x)
  s <- digest::digest2int(paste(unlist(fields), collapse = "|"))
  set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  invisible(s)
}

#' Specify the authorized pilot or a small verification slice
#' @param profile One of smoke, benchmark, pilot.
#' @return List of fixed DGP parameters and Monte Carlo replication counts.
#' @examples
#' p <- allocation_parameters("pilot")
allocation_parameters <- function(profile = "pilot") {
  stopifnot(profile %in% c("smoke", "benchmark", "pilot"))
  p <- list(profile = profile, grid_dim = 10, nb_type = "queen", tau = 1,
            beta = 1, sigma = 1, base_rate = 35 / 100000,
            pop_per_cluster = 100000, pop_mode = "equal",
            configs = list(list(mode = "iid", rho_x = 0),
              list(mode = "spatial", rho_x = 0.2), list(mode = "spatial", rho_x = 0.5),
              list(mode = "poisson", rho_x = 0.2), list(mode = "poisson", rho_x = 0.5)),
            surfaces = 1:2, rho = c(0, 0.5), gamma = c(0.5, 0.8),
            regimes = c("control_only", "both"), designs = 1:9,
            n_allocations = 100L, n_outcomes = 100L, noise_batch_size = 25L,
            seed_namespace = "allocation-risk-v1")
  if (profile != "pilot") {
    p$configs <- p$configs[1]; p$surfaces <- 1L; p$rho <- 0.5
    p$gamma <- 0.5; p$regimes <- "both"
    p$n_allocations <- if (profile == "smoke") 4L else 10L
    p$n_outcomes <- if (profile == "smoke") 8L else 100L
  }
  p
}

#' Regenerate one existing keyed incidence surface without touching main files
#' @param cfg List with mode and rho_x.
#' @param k Surface index.
#' @param grid Grid from build_spatial_grid().
#' @param p Parameters from allocation_parameters().
#' @return Numeric length-N incidence vector identical to the main runner's X_k.
#' @examples
#' # X <- allocation_surface(p$configs[[1]], 1, build_spatial_grid(), p)
allocation_surface <- function(cfg, k, grid, p) {
  allocation_seed("X", cfg$mode, cfg$rho_x, k)
  generate_incidence(cfg$mode, grid$N_clusters, 1, grid$W_queen, cfg$rho_x,
    base_rate = p$base_rate, pop_per_cluster = p$pop_per_cluster,
    pop_mode = p$pop_mode)[, 1]
}

#' Draw an extensible sequence of assignments from the existing design rules
#' @param u One-row block data frame, including design and surface.
#' @param cfg Config list.
#' @param X Fixed incidence vector.
#' @param grid Grid object.
#' @param p Parameter list.
#' @return N x J matrix; Checkerboard has only one column.
#' @examples
#' # Z <- allocation_draws(u, cfg, X, grid, p)
allocation_draws <- function(u, cfg, X, grid, p) {
  J <- if (is_design_deterministic(u$design)) 1L else p$n_allocations
  spatial <- get_active_spatial(grid, p$nb_type)
  vapply(seq_len(J), function(j) {
    allocation_seed(p$seed_namespace, "Z", cfg$mode, cfg$rho_x, p$nb_type,
                    u$rho, u$gamma, u$regime, u$surface, u$design, j)
    get_designs(u$design, 1, length(X), X, spatial$nb, grid$coords)[, 1]
  }, numeric(length(X)))
}

#' Generate one fixed-size independent outcome-noise batch
#' @param key Block/assignment key fields, excluding batch number.
#' @param batch Integer batch index.
#' @param N Number of clusters.
#' @param p Parameters; batch size stays fixed when replication increases.
#' @return N x noise_batch_size Gaussian noise matrix.
#' @examples
#' e <- allocation_noise(list("example"), 1, 100, allocation_parameters())
allocation_noise <- function(key, batch, N, p) {
  do.call(allocation_seed, c(list(p$seed_namespace, "eps"), key, list(batch)))
  matrix(rnorm(N * p$noise_batch_size, sd = p$sigma), nrow = N)
}

# Conditional summaries ----
#' Summarize repeated oracle fits for one fixed assignment
#'
#' m(a|X) = Eε[(τ̂−τ)²]; SE_MSE² estimates Varε((τ̂−τ)²)/R.
#' Point estimates retain finite τ̂ even if its CI SE fails. Coverage has a separate
#' denominator; any failure marks the block incomplete for design comparisons.
#' @param fits Data frame with replicate, error, se, alias, messages.
#' @param tau True treatment effect.
#' @return One-row data frame of MSE, bias, coverage, MC errors and two halves.
#' @examples
#' # allocation_fit_summary(fits, 1)
allocation_fit_summary <- function(fits, tau) {
  ok <- is.finite(fits$error); ci <- ok & is.finite(fits$se) & fits$se >= 0
  e <- fits$error[ok]; l2 <- e^2; R <- nrow(fits)
  half <- fits$replicate <= floor(R / 2)
  half_mse <- function(h) if (any(ok & h)) mean(fits$error[ok & h]^2) else NA_real_
  data.frame(N_Outcome = R, N_Valid_Est = sum(ok), N_Valid_CI = sum(ci),
    MSE = if (length(e)) mean(l2) else NA_real_,
    Bias = if (length(e)) mean(e) else NA_real_,
    Coverage = if (any(ci)) mean(abs(fits$error[ci]) <= qnorm(0.975) * fits$se[ci]) else NA_real_,
    SE_MSE = if (length(e) > 1) sd(l2) / sqrt(length(e)) else NA_real_,
    SE_Bias = if (length(e) > 1) sd(e) / sqrt(length(e)) else NA_real_,
    SE_Coverage = if (sum(ci) > 1) {
      sd(as.numeric(abs(fits$error[ci]) <= qnorm(0.975) * fits$se[ci])) / sqrt(sum(ci))
    } else NA_real_,
    MSE_Half1 = half_mse(half), MSE_Half2 = half_mse(!half),
    Fail_Rate_Est = mean(!ok), Fail_Rate_CI = mean(!ci),
    N_Aliased = sum(fits$alias), N_Warn = sum(nzchar(fits$messages)),
    Complete = all(ok & ci))
}

#' Estimate allocation variance and tail metrics while retaining draw frequencies
#'
#' For n sampled draws, unique allocations i occur f_i times. Their independently
#' estimated risks have MC variance v_i. Reused cache results induce covariance
#' between duplicate draws: the noise contribution to the draw sample variance is
#' Σ_i f_i(1−f_i/n)v_i/(n−1), not simply mean(v_i). Subtract this from the raw
#' weighted sample variance. Signed negative corrections are reported verbatim;
#' their square root is NA. A single sampled unique allocation has zero sampled
#' allocation variance; population variance is zero only with proven singleton
#' support. Its conditional MSE can remain high and uncertain.
#' @param a Unique-allocation summary with Frequency, MSE, SE_MSE and halves.
#' @return One-row data frame; tails concern estimated risks and are not deconvolved.
#' @examples
#' # allocation_distribution(a)
allocation_distribution <- function(a) {
  n <- sum(a$Frequency); w <- a$Frequency / n; single <- nrow(a) == 1L
  complete <- all(a$Complete) && all(is.finite(a$MSE))
  if (!complete) return(data.frame(N_Draws = n, N_Unique = nrow(a), Complete = FALSE,
    Mean_MSE = NA_real_, Variance_Raw = NA_real_, Variance_MC_Noise = NA_real_,
    Variance_Corrected = NA_real_, SD_Raw = NA_real_, SD_Corrected = NA_real_,
    Q90_Estimated = NA_real_, Worst10_Mean_Estimated = NA_real_, Sampled_Max_Estimated = NA_real_,
    SE_Mean_MSE_Noise = NA_real_, SE_Mean_MSE_Joint = NA_real_,
    Variance_Split_Half = NA_real_, Variance_Noise_Fraction = NA_real_,
    Q90_Half1 = NA_real_, Q90_Half2 = NA_real_, Worst10_Half1 = NA_real_, Worst10_Half2 = NA_real_,
    Tail_Half_Jaccard = NA_real_, Tail_Half_Rank_Cor = NA_real_,
    Tail_Boundary_Uncertain = NA, No_Sampled_Allocation_Variation = single))
  mu <- sum(w * a$MSE)
  raw <- if (single) 0 else sum(a$Frequency * (a$MSE - mu)^2) / (n - 1)
  noise <- if (single) 0 else sum(a$Frequency * (1 - w) * a$SE_MSE^2) / (n - 1)
  corrected <- raw - noise
  draw <- rep(seq_len(nrow(a)), a$Frequency); m <- a$MSE[draw]
  top_n <- max(1L, ceiling(n * 0.1)); order_full <- order(m, decreasing = TRUE)
  top <- order_full[seq_len(top_n)]
  top1 <- order(a$MSE_Half1[draw], decreasing = TRUE)[seq_len(top_n)]
  top2 <- order(a$MSE_Half2[draw], decreasing = TRUE)[seq_len(top_n)]
  half_metrics <- function(v) {
    expanded <- v[draw]
    c(Q90 = unname(quantile(expanded, 0.9, type = 7)),
      Worst10 = mean(sort(expanded, decreasing = TRUE)[seq_len(top_n)]))
  }
  hm1 <- half_metrics(a$MSE_Half1); hm2 <- half_metrics(a$MSE_Half2)
  cov_halves <- if (single) 0 else sum(a$Frequency *
    (a$MSE_Half1 - sum(w * a$MSE_Half1)) *
    (a$MSE_Half2 - sum(w * a$MSE_Half2))) / (n - 1)
  # Joint SE includes allocation draws and cached-noise covariance. The expanded
  # nonnegative expression equals Variance_Corrected/n + Σw_i²v_i.
  joint_var <- if (n == 1) a$SE_MSE^2 else raw / n +
    sum(a$Frequency * (a$Frequency - 1) * a$SE_MSE^2) / (n * (n - 1))
  # Approximate pointwise ±1.96 MC SE intervals straddling the tail boundary signal
  # uncertain membership. This is a diagnostic, not simultaneous tail inference.
  boundary <- min(m[top]); lo <- a$MSE - 1.96 * a$SE_MSE; hi <- a$MSE + 1.96 * a$SE_MSE
  data.frame(N_Draws = n, N_Unique = nrow(a), Complete = TRUE, Mean_MSE = mu,
    Variance_Raw = raw, Variance_MC_Noise = noise, Variance_Corrected = corrected,
    SD_Raw = sqrt(raw), SD_Corrected = if (corrected >= 0) sqrt(corrected) else NA_real_,
    Q90_Estimated = unname(quantile(m, 0.9, type = 7)),
    Worst10_Mean_Estimated = mean(m[top]), Sampled_Max_Estimated = max(m),
    SE_Mean_MSE_Noise = sqrt(sum(w^2 * a$SE_MSE^2)),
    SE_Mean_MSE_Joint = sqrt(joint_var), Variance_Split_Half = cov_halves,
    Variance_Noise_Fraction = if (raw > 0) noise / raw else NA_real_,
    Q90_Half1 = hm1[["Q90"]], Q90_Half2 = hm2[["Q90"]],
    Worst10_Half1 = hm1[["Worst10"]], Worst10_Half2 = hm2[["Worst10"]],
    Tail_Half_Jaccard = if (single) 1 else length(intersect(top1, top2)) / length(union(top1, top2)),
    Tail_Half_Rank_Cor = if (single) NA_real_ else suppressWarnings(cor(a$MSE_Half1[draw], a$MSE_Half2[draw], method = "spearman")),
    Tail_Boundary_Uncertain = !single && sum(lo <= boundary & hi >= boundary) > 1,
    No_Sampled_Allocation_Variation = single)
}

# One fixed-X/design/scenario block ----
#' Fit cached unique allocations with independent noise and preserve all draws
#' @param u One-row block data frame.
#' @param cfg Incidence config.
#' @param X Fixed incidence vector.
#' @param grid Spatial grid object.
#' @param setup Validated ML setup.
#' @param p Parameters.
#' @param old Optional matching checkpoint used to extend outcomes/allocations.
#' @return List: draws, allocations, summary, fits (squared errors recoverable), elapsed.
#' @examples
#' # out <- allocation_block(u, cfg, X, grid, setup, p)
allocation_block <- function(u, cfg, X, grid, setup, p, old = NULL) {
  t0 <- proc.time()[[3]]
  Z <- allocation_draws(u, cfg, X, grid, p)
  hashes <- apply(Z, 2, digest::digest, algo = "xxhash64")
  unique_hashes <- unique(hashes)
  spatial <- get_active_spatial(grid, p$nb_type); W <- spatial$W
  inv_filter <- solve(diag(length(X)) - u$rho * W)
  keys <- data.frame(Incidence_Mode = cfg$mode, Rho_Incidence = cfg$rho_x,
    Surface = u$surface, Neighbor_Type = p$nb_type, Rho = u$rho, Gamma = u$gamma,
    Spillover_Type = u$regime, True_Tau = p$tau, Design = paste("Design", u$design),
    Design_Name = unname(get_design_names(u$design)), Estimator = "oracle")
  all_fits <- list(); rows <- list()
  for (h in unique_hashes) {
    z <- Z[, match(h, hashes)]; wz <- as.vector(W %*% z)
    spill <- u$gamma * wz * if (u$regime == "control_only") (1 - z) else 1
    xm <- cbind("(Intercept)" = 1, Z = z, Spill = spill, X = X)
    fits <- if (!is.null(old) && h %in% names(old$fits)) old$fits[[h]] else NULL
    if (!is.null(fits)) fits <- fits[fits$replicate <= p$n_outcomes, , drop = FALSE]
    previous <- if (is.null(fits)) 0L else nrow(fits)
    noise_key <- list(cfg$mode, cfg$rho_x, p$nb_type, u$rho, u$gamma,
                      u$regime, u$surface, u$design, h)
    new_rows <- list()
    if (previous < p$n_outcomes) {
      # Generate each fixed-size batch in full, selecting only requested columns.
      # This preserves prefix draws even when an extension ends inside a batch.
      batches <- unique((seq.int(previous + 1L, p$n_outcomes) - 1L) %/% p$noise_batch_size + 1L)
      for (b in batches) {
        eps <- allocation_noise(noise_key, b, length(X), p)
        Y <- inv_filter %*% (p$tau * z + spill + p$beta * X + eps)
        ids <- (b - 1L) * p$noise_batch_size + seq_len(p$noise_batch_size)
        for (col in which(ids > previous & ids <= p$n_outcomes)) {
          fit <- fit_one_lag_model(Y[, col], xm, setup, engine = "lean")
          new_rows[[length(new_rows) + 1L]] <- data.frame(replicate = ids[col],
            error = fit$tau - p$tau, se = fit$se,
            alias = any(startsWith(fit$warns, "Aliased variables found")),
            messages = paste(fit$warns, collapse = " || "))
        }
      }
      fits <- rbind(fits, do.call(rbind, new_rows))
    }
    all_fits[[h]] <- fits
    rows[[length(rows) + 1L]] <- cbind(keys, Assignment_Hash = h,
      Frequency = sum(hashes == h), Assignment = paste0(z, collapse = ""),
      Structural_Singleton = u$design == 1 || (u$design == 2 &&
        sort(X)[length(X) / 2] < sort(X)[length(X) / 2 + 1]),
      N_Treated = sum(z), Z_WZ_rank_deficient = qr(cbind(1, z, wz))$rank < 3,
      Oracle_Model_Rank = qr(xm)$rank, allocation_fit_summary(fits, p$tau))
  }
  a <- do.call(rbind, rows)
  list(draws = cbind(keys, Draw = seq_along(hashes), Assignment_Hash = hashes),
       allocations = a, summary = cbind(keys, allocation_distribution(a)),
       fits = all_fits, elapsed_sec = proc.time()[[3]] - t0)
}

# Guarded runner ----
#' Execute isolated, manifest-checked blocks; save failures before stopping
#' @param p Parameters.
#' @param workers Integer local fork workers (default 8).
#' @return Invisibly, output directory. Each saved fit retains error and SE.
#' @examples
#' # allocation_main(allocation_parameters("pilot"), workers = 8)
allocation_main <- function(p, workers = 8L) {
  stopifnot(p$n_outcomes >= 4, p$n_allocations >= 1, workers >= 1, workers <= 8)
  vars <- c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
  if (!all(Sys.getenv(vars) == "1")) stop("Set all BLAS thread variables to 1 before R starts")
  root <- file.path(dirname(allocation_script_dir), "results", "allocation_risk", p$profile)
  cp_dir <- file.path(root, "checkpoints")
  dir.create(cp_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(allocation_script_dir, c(sprintf("%02d_%s.R", 1:4,
      c("spatial_setup", "incidence_generation", "designs", "estimation")), "16_allocation_risk.R"))
  # Replication can increase under an unchanged model manifest. Every invocation
  # separately records its requested J/R and worker count in request_<time>.rds.
  model_p <- p; model_p$n_allocations <- model_p$n_outcomes <- NULL
  manifest <- list(params = model_p, param_hash = digest::digest(model_p),
    source_hashes = setNames(vapply(files, digest::digest, character(1), file = TRUE), basename(files)),
    package_versions = setNames(vapply(c("spdep", "spatialreg", "dplyr", "digest"),
      function(x) as.character(packageVersion(x)), character(1)), c("spdep", "spatialreg", "dplyr", "digest")),
    R = R.version.string, BLAS = sessionInfo()$BLAS)
  mf <- file.path(root, "manifest.rds")
  if (file.exists(mf)) {
    if (!identical(readRDS(mf), manifest)) stop("Allocation-risk manifest mismatch; refusing cached fits: ", root)
  } else {
    if (length(list.files(cp_dir))) stop("Allocation-risk checkpoints have no manifest")
    saveRDS(manifest, mf)
  }
  saveRDS(list(params = p, workers = workers, started = Sys.time()),
          file.path(root, paste0("request_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")))
  grid <- build_spatial_grid(p$grid_dim); setup <- sar_lag_setup(get_active_spatial(grid, p$nb_type)$listw)
  surfaces <- lapply(p$configs, function(cfg) lapply(p$surfaces, function(k) allocation_surface(cfg, k, grid, p)))
  units <- expand.grid(config = seq_along(p$configs), surface = p$surfaces, rho = p$rho,
    gamma = p$gamma, regime = p$regimes, design = p$designs, stringsAsFactors = FALSE)
  cat(sprintf("Allocation risk %s: %d blocks; J=%d, R=%d, %d workers\n", p$profile, nrow(units), p$n_allocations, p$n_outcomes, workers))
  t0 <- proc.time()[[3]]
  run <- function(i) {
    u <- units[i, ]; cp <- file.path(cp_dir, sprintf("block_%04d.rds", i))
    old <- if (file.exists(cp)) readRDS(cp) else NULL
    out <- allocation_block(u, p$configs[[u$config]], surfaces[[u$config]][[match(u$surface, p$surfaces)]],
                            grid, setup, p, old)
    temp <- paste0(cp, ".tmp"); saveRDS(out, temp)
    if (!file.rename(temp, cp)) stop("Could not atomically save checkpoint ", cp)
    cat(sprintf("block %d/%d: %d unique allocations, %.1f s\n", i, nrow(units), nrow(out$allocations), out$elapsed_sec))
    # Fit arrays stay in block checkpoints, keeping aggregate worker returns small.
    out$fits <- NULL; out
  }
  results <- if (workers == 1) lapply(seq_len(nrow(units)), run) else
    parallel::mclapply(seq_len(nrow(units)), run, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  bad <- vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
  if (any(bad)) stop("Allocation-risk worker failed for blocks: ", paste(which(bad), collapse = ", "))
  a <- do.call(rbind, lapply(results, `[[`, "allocations"))
  s <- do.call(rbind, lapply(results, `[[`, "summary"))
  draws <- do.call(rbind, lapply(results, `[[`, "draws"))
  write.csv(a, file.path(root, "allocation_metrics.csv"), row.names = FALSE)
  write.csv(s, file.path(root, "conditional_summary.csv"), row.names = FALSE)
  write.csv(draws, file.path(root, "draw_mapping.csv"), row.names = FALSE)
  saveRDS(list(allocations = a, summary = s, draws = draws, params = p), file.path(root, "results.rds"))
  writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
  status <- list(elapsed_sec = proc.time()[[3]] - t0, complete = all(a$Complete),
    requested_fits = sum(a$N_Outcome), valid_estimates = sum(a$N_Valid_Est),
    valid_CI = sum(a$N_Valid_CI), alias = sum(a$N_Aliased), warns = sum(a$N_Warn),
    negative_variance_corrections = sum(s$Variance_Corrected < 0, na.rm = TRUE))
  saveRDS(status, file.path(root, "status.rds")); print(status)
  if (!status$complete) stop("Partial fit failure: saved diagnostics; design summaries withheld for incomplete blocks")
  invisible(root)
}

if (!exists("allocation_define_only") || !isTRUE(allocation_define_only)) {
  args <- commandArgs(TRUE)
  allocation_main(allocation_parameters(if (length(args)) args[1] else "pilot"),
                  if (length(args) > 1) as.integer(args[2]) else 8L)
}
