# ============================================================
# Script: methods_ext_core.R
# Purpose: Shared helpers for the 2026-10-20 methods extensions (participation,
#   power, coverage and estimation sensitivities) on the 58 NC clusters.
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: real_sud_simulation.R and real_sud_setup.R (sourced, unchanged; they source
#   application_designs.R and code/16_allocation_risk.R -> code/01-04),
#   spdep, spatialreg, digest
# ============================================================
# Plan: docs/plans/methods_extensions_plan_2026-10-20.md (D1-D9).
# Nothing here edits or re-saves a production object. The frozen production
# setup (regions, yearly blocks, observed ranks) is read from setup.rds, and the
# fitter is the lean eigen ML engine fit_sar_lag() from code/04_estimation.R.

# Paths and unchanged production helpers ----
.me_code_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(normalizePath(dirname(f)))
  }
  stop("Source methods_ext_core.R by its file path")
})
source(file.path(.me_code_dir, "real_sud_simulation.R"))   # sources real_sud_setup.R; defines rs_*, allocation_seed, fit_sar_lag
.me_root <- file.path(.rs_app, "results", "methods_ext_2026-10")
.me_prod_root <- file.path(.rs_app, "results", "real_sud_rev_20261002")

#' Key-seed a methods-extension draw (new namespace; never collides with production)
#'
#' Same mechanism as production rs_seed(): numeric key fields are formatted with
#' two decimals and all RNG kinds are set explicitly. Effect sizes (tau, kappa, d)
#' are never part of a key, so all effects share allocations and noise
#' (common random numbers).
#' @param ... Key fields, e.g. "pp-subset", "random", 20, 3.
#' @return Invisibly, the integer seed.
#' @family methods_ext
#' @seealso allocation_seed, rs_seed
#' @examples
#' me_seed("pp-subset", "random", 20, 1)
me_seed <- function(...) allocation_seed("nc-methods-ext-v1", ...)

#' Load the frozen production setup and confirm it is the one production used
#'
#' The design digest (regions, yearly blocks, ids) must equal the production
#' manifest's, so N = 58 runs reuse exactly the production partitions (V1).
#' @return The setup list from real_sud_rev_20261002/setup.rds.
#' @family methods_ext
#' @examples
#' s <- me_load_setup()
me_load_setup <- function() {
  s <- readRDS(file.path(.me_prod_root, "setup.rds"))
  m <- readRDS(file.path(.me_prod_root, "production", "manifest.rds"))
  if (!identical(m$design, digest::digest(list(s$regions, s$blocks, s$ids))))
    stop("Frozen setup does not match the production manifest design digest")
  s
}

# Participation ----
#' Randomized systematic probability-proportional-to-size sampling
#'
#' Units are put in random order and a systematic sample with a uniform random
#' start is taken on the cumulative inclusion probabilities. Each unit k is then
#' selected with probability exactly pi_k (requires 0 < pi_k <= 1, sum(pi) = n).
#' @param pi Inclusion probabilities summing to the sample size.
#' @return Integer vector of selected unit indices (unsorted).
#' @family methods_ext
#' @examples
#' me_upsystematic(rep(0.5, 10))
me_upsystematic <- function(pi) {
  n <- round(sum(pi))
  if (abs(sum(pi) - n) > 1e-8 || any(pi <= 0) || any(pi > 1)) stop("Invalid inclusion probabilities")
  o <- sample.int(length(pi))
  cp <- cumsum(pi[o]); cp[length(cp)] <- n          # remove floating drift at the end
  u <- stats::runif(1)
  sel <- floor(cp - u) - floor(c(0, cp[-length(cp)]) - u)
  out <- o[sel == 1]
  if (length(out) != n) stop("Systematic sample has the wrong size")
  out
}

#' Participation weights for the incidence-weighted (3x) mechanism
#'
#' w_j = 1 + 2 (r_j - 1) / 57 for the year's incidence rank r_j (1 = lowest), so
#' the highest-burden college has 3 times the lowest one's weight. Inclusion
#' probabilities are pi_j = N w_j / sum(w), realized exactly by me_upsystematic().
#' @param rank Year's incidence ranks (1..58).
#' @param N Number of participants.
#' @return Inclusion probabilities pi_j.
#' @family methods_ext
#' @examples
#' # me_incidence_pi(setup$annual[["2018"]]$rank, 20)
me_incidence_pi <- function(rank, N) {
  w <- 1 + 2 * (rank - 1) / (length(rank) - 1)
  pi <- N * w / sum(w)
  if (any(pi > 1)) stop("Inclusion probability above 1 at N = ", N)
  pi
}

#' Draw the s-th participant subset
#'
#' Random volunteering: a simple random sample of N colleges, keyed without the
#' year or design, so every year and design uses the same subsets. Incidence-
#' weighted: pi-proportional to the year's rank weights (year-specific key).
#' @param mech "random" or "incidence3x".
#' @param N Number of participants (58 = everyone).
#' @param s Subset index.
#' @param year Year label (used only by the incidence-weighted mechanism).
#' @param setup Frozen setup.
#' @return Sorted integer indices into the 58 named clusters.
#' @family methods_ext
#' @examples
#' # me_subset("random", 20, 1, "2018", setup)
me_subset <- function(mech, N, s, year, setup) {
  if (N == 58) return(seq_len(58))
  if (mech == "random") {
    me_seed("pp-subset", "random", N, s)
    return(sort(sample.int(58, N)))
  }
  if (mech != "incidence3x") stop("Unknown mechanism: ", mech)
  me_seed("pp-subset", "incidence3x", year, N, s)
  sort(me_upsystematic(me_incidence_pi(setup$annual[[year]]$rank, N)))
}

# Participant geometry and spillover exposure ----
#' Build the participant-only spatial structures for one subset
#'
#' W_P: queen contiguity among participants, row-standardized; a participant with
#' no participating neighbor (island) has a zero row. The listw object is built
#' with zero.policy = TRUE and must pass sar_lag_setup()'s can.be.simmed() guard.
#' Exposure (rule A) uses all geographic queen neighbors: e_j = (sum of treated
#' neighbors) / (number of geographic neighbors); non-participants are untreated.
#' @param P Sorted participant indices into the 58 clusters.
#' @param setup Frozen setup.
#' @param need_engine Build the fitter setup (FALSE for design-only checks).
#' @return List: P, n, AP (binary), WP, nb (plain integer list), deg_full,
#'   islands, components, n_edges, engine (sar_lag_setup) or a failure message.
#' @family methods_ext
#' @examples
#' # g <- me_geometry(1:58, setup)
me_geometry <- function(P, setup, need_engine = TRUE) {
  Adj <- setup$weights$queen$adjacency
  AP <- Adj[P, P, drop = FALSE]
  dP <- rowSums(AP)
  WP <- AP
  WP[dP > 0, ] <- AP[dP > 0, , drop = FALSE] / dP[dP > 0]
  nb <- lapply(seq_along(P), function(i) which(AP[i, ] > 0))
  g <- list(P = P, n = length(P), AP = AP, WP = WP, nb = nb,
            deg_full = rowSums(Adj)[P], islands = sum(dP == 0),
            n_edges = sum(AP) / 2, engine = NULL, engine_error = NA_character_,
            listw_warnings = character(0))
  # Connected components of the participant graph (islands count as components)
  comp <- rep(NA_integer_, g$n); k <- 0L
  for (i in seq_len(g$n)) if (is.na(comp[i])) {
    k <- k + 1L; stack <- i
    while (length(stack)) { v <- stack[1]; stack <- stack[-1]
      if (is.na(comp[v])) { comp[v] <- k; stack <- c(stack, nb[[v]][is.na(comp[nb[[v]]])]) } }
  }
  g$components <- k
  if (!need_engine) return(g)
  if (g$n_edges == 0) { g$engine_error <- "no participant edges"; return(g) }
  lw <- withCallingHandlers(
    spdep::mat2listw(WP, style = "W", zero.policy = TRUE),
    warning = function(w) { g$listw_warnings <<- c(g$listw_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (max(abs(spdep::listw2mat(lw) - WP)) > 1e-12) stop("listw does not reproduce W_P")
  g$listw <- lw
  g$guard_bypass <- spatialreg::can.be.simmed(lw)
  g$engine <- tryCatch(me_engine_setup(lw), error = function(e) { g$engine_error <<- conditionMessage(e); NULL })
  g
}

#' Lean-engine setup for a participant listw, including symmetric W_P
#'
#' sar_lag_setup() (unchanged) refuses listw objects that can.be.simmed(), because
#' lagsarlm() then takes its symmetric-similar eigen branch. A row-standardized W_P
#' is itself symmetric when every non-island participant has the same within-subset
#' degree (e.g. all components are pairs). Similar matrices share eigenvalues, so
#' for those subsets a copy of the listw is marked non-symmetric and passed to the
#' unchanged sar_lag_setup(); V5 checks these fits against lagsarlm(), which uses
#' its own symmetric branch on the original listw. Counted per subset (Guard_Bypass).
#' @param lw Participant listw (style "W", zero.policy = TRUE).
#' @return sar_lag_setup() list.
#' @family methods_ext
#' @seealso sar_lag_setup
#' @examples
#' # me_engine_setup(g$listw)
me_engine_setup <- function(lw) {
  if (!spatialreg::can.be.simmed(lw)) return(sar_lag_setup(lw))
  W <- spdep::listw2mat(lw)
  if (!isSymmetric(unname(W))) stop("can.be.simmed() but W_P is not symmetric")
  lw2 <- lw; attr(lw2$weights, "glistsym") <- FALSE
  sar_lag_setup(lw2)
}

#' Rule-A spillover exposure for participants
#'
#' e_j = (sum over all geographic queen neighbors k of z_k) / |nb(j)|, where
#' non-participants have z_k = 0. Unscaled (the DGP multiplies by kappa*tau).
#' @param z Treatment among participants (length n).
#' @param g Participant geometry from me_geometry().
#' @return Numeric exposure vector (length n).
#' @family methods_ext
#' @examples
#' # me_exposure(z, g)
me_exposure <- function(z, g) as.vector(g$AP %*% z) / g$deg_full

#' Spillover forms for the A7 sensitivity (N = 58)
#'
#' proportional: e_j as in me_exposure(); binary: 1{any treated queen neighbor};
#' count: (number of treated queen neighbors) / dbar, dbar = mean queen degree.
#' @param z Treatment vector over the participants.
#' @param g Participant geometry.
#' @param form "proportional", "binary" or "count".
#' @param dbar Mean queen degree of the full map (count form only).
#' @return Unscaled exposure vector.
#' @family methods_ext
#' @examples
#' # me_exposure_form(z, g, "binary")
me_exposure_form <- function(z, g, form, dbar = NULL) {
  cnt <- as.vector(g$AP %*% z)
  switch(form, proportional = cnt / g$deg_full, binary = as.numeric(cnt > 0),
         count = cnt / dbar, stop("Unknown spillover form"))
}

# Designs on a participant subset ----
#' Assemble the frozen (N = 58) or rebuilt (N < 58) design inputs
#'
#' N = 58 reuses the production regions, yearly 2x2 blocks and rank-scaled X,
#' with the regional mean-rank summary, exactly as rs_draws(). For N < 58 the
#' saturation regions (4 k-means regions), the 2x2 blocking analogue (target size 4)
#' and the participant ranks are rebuilt on the subset; regions use the pooled
#' population, as in production, and blocks the year's incidence and population.
#' @param g Participant geometry.
#' @param year Year label.
#' @param setup Frozen setup.
#' @param key List of key fields identifying the subset (mechanism, N, s, year if
#'   the subset is year-specific).
#' @return List of design inputs plus region/block diagnostics.
#' @family methods_ext
#' @examples
#' # di <- me_design_inputs(g, "2018", setup, list("random", 20, 1))
me_design_inputs <- function(g, year, setup, key) {
  a <- setup$annual[[year]]; P <- g$P
  if (g$n == 58) {
    return(list(X = a$X, rank = a$rank, nb = setup$weights$queen$nb, coords = setup$coords,
      pop = a$person_years, county = a$n_counties, regions = setup$regions,
      blocks = setup$blocks[[year]], region_summary = rs_region_values(a, setup$regions, "mean_rank"),
      rates = a$rate_per_100k, frozen = TRUE))
  }
  rP <- rank(a$rate_per_100k[P], ties.method = "average")
  XP <- (rP - 0.5) / g$n
  pooled_pop <- setup$annual[["2018-2021"]]$person_years / 4
  rseed <- do.call(me_seed, c(list("pp-regions"), key))
  reg <- make_kmeans_regions(setup$coords[P, , drop = FALSE], a$n_counties[P], pooled_pop[P], XP,
                             seed = rseed)
  do.call(me_seed, c(list("pp-blocks"), key, list(year)))
  blk <- make_spatial_blocks(setup$coords[P, , drop = FALSE], XP, a$person_years[P])
  list(X = XP, rank = rP, nb = g$nb, coords = setup$coords[P, , drop = FALSE],
       pop = a$person_years[P], county = a$n_counties[P], regions = reg$region_id,
       blocks = blk, region_summary = tapply(rP, reg$region_id, mean),
       rates = a$rate_per_100k[P], frozen = FALSE)
}

#' Is a design deterministic on this subset?
#' Checkerboard always; High Incidence Focus when the participants' rates have no ties.
#' @param id Design ID (1-9).
#' @param di Design inputs.
#' @return Logical.
#' @family methods_ext
#' @examples
#' # me_deterministic(2, di)
me_deterministic <- function(id, di) id == 1 || (id == 2 && !anyDuplicated(di$rates))

#' Draw one allocation of design id on the participants
#'
#' Designs 1, 3-8 call get_application_designs() unchanged with the participant
#' inputs. High Incidence Focus (2) and SRS (9) treat floor(N/2) (approved 2026-10-20;
#' identical to production's round(N/2) at N = 58, where the production call is used).
#' @param id Design ID.
#' @param di Design inputs.
#' @return Integer 0/1 vector of length n.
#' @family methods_ext
#' @examples
#' # me_draw_design(9, di)
me_draw_design <- function(id, di) {
  n <- length(di$X)
  if (!di$frozen && id == 2) return(as.integer(application_random_tie_rank(di$X) > n - floor(n / 2)))
  if (!di$frozen && id == 9) { k <- floor(n / 2); return(sample(rep(c(1L, 0L), c(k, n - k)))) }
  get_application_designs(id, 1L, di$X, di$nb, di$coords, di$pop, di$county,
                          di$regions, block_id = di$blocks, region_summary = di$region_summary)[, 1]
}

# Fitting ----
#' Fit one outcome with the lean engine and classify the result
#'
#' Boundary-rho fits and errors are invalid (as in production); aliased columns
#' are counted but the fit is kept when Z itself is estimable.
#' @param y Outcome.
#' @param xm Model matrix (with "(Intercept)" and "Z").
#' @param engine sar_lag_setup() list.
#' @return c(tau, se, alias, bound, error, warn) as a numeric vector.
#' @family methods_ext
#' @examples
#' # me_fit(y, xm, g$engine)
me_fit <- function(y, xm, engine) {
  f <- fit_one_lag_model(y, xm, engine)
  bound <- any(grepl("interval bound", f$warns)); err <- any(grepl("^ERROR:", f$warns))
  c(tau = if (bound || err) NA_real_ else f$tau, se = if (bound || err) NA_real_ else f$se,
    alias = any(grepl("Aliased", f$warns)), bound = bound, error = err,
    warn = length(f$warns) > 0)
}

#' Hash files for a methods-extension manifest
#' @param paths File paths.
#' @return Named SHA-256 vector.
#' @family methods_ext
#' @examples
#' # me_hash(file.path(.me_code_dir, "methods_ext_core.R"))
me_hash <- function(paths) rs_hash_files(normalizePath(paths))

#' Write or check a manifest for one methods-extension output folder
#'
#' Records the configuration, the hashes of the computation files (new and
#' sourced production files), package versions and BLAS. An existing manifest
#' that differs stops execution: checkpoints are never reused across versions.
#' @param dir Output folder.
#' @param config Analysis configuration list.
#' @param files Computation files to hash.
#' @return Manifest invisibly.
#' @family methods_ext
#' @examples
#' # me_manifest(dir, config, files)
me_manifest <- function(dir, config, files) {
  prod <- file.path(.rs_project, "code", c("01_spatial_setup.R", "02_incidence_generation.R",
    "03_designs.R", "04_estimation.R", "16_allocation_risk.R"))
  prod <- c(prod, file.path(.rs_code_dir, c("real_sud_setup.R", "application_designs.R")))
  pkgs <- c("spdep", "spatialreg", "digest")
  m <- list(config = config, config_hash = digest::digest(config),
            source = me_hash(c(files, prod)),
            setup_digest = digest::digest(readRDS(file.path(.me_prod_root, "setup.rds"))),
            R = R.version.string,
            versions = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)),
            BLAS = unname(extSoftVersion()["BLAS"]),
            RNG = c("Mersenne-Twister", "Inversion", "Rejection"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(dir, "manifest.rds")
  if (file.exists(p) && !identical(readRDS(p)[names(m)], m)) stop("Manifest mismatch in ", dir, ": refuse checkpoint reuse")
  if (!file.exists(p)) saveRDS(m, p)
  invisible(m)
}

#' Refuse to run with multithreaded BLAS settings (forked workers need 1 thread)
#' @return Invisibly TRUE; stops otherwise.
#' @family methods_ext
#' @examples
#' # me_check_threads()
me_check_threads <- function() {
  v <- Sys.getenv(c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"))
  if (!all(v == "1")) stop("Set VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1")
  invisible(TRUE)
}
