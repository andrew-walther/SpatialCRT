# ============================================================
# Script: test_methods_extensions.R
# Purpose: Verification checks V3-V7 and V11 for the 2026-10-20 methods
#   extensions (scale invariance, rule A, islands in the fitter, designs on
#   subsets, informative sampling, spillover-form exposures).
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: methods_ext_core.R, methods_ext_participation.R, spatialreg, spdep
# ============================================================
# Run from the project root:
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript application/tests/test_methods_extensions.R
# V4 and V11's toy-map parts need no restricted data; the rest uses setup.rds
# (the frozen observed ranks), which exists only where the data are available.
# SKIP_LAGSARLM=1 skips the ~2 min lagsarlm part of V5.

test_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
methods_ext_define_only <- TRUE
source(file.path(test_dir, "..", "code", "methods_ext_participation.R"))
source(file.path(test_dir, "..", "code", "real_sud_simulation.R"))   # rs_config/rs_units/rs_model (unchanged)
log_lines <- character(0)

#' Assert a named behavior; record it and stop at the first failure
#' @param name Check description.
#' @param ok Logical condition, required TRUE.
#' @param detail Optional numeric/text detail printed beside the result.
#' @param stop_on_fail FALSE records a failure and continues (used for V6's region
#'   check, whose failure concerns two designs' partitions, not the runner).
#' @return Invisibly TRUE; otherwise stops.
#' @examples
#' check("arithmetic", 1 + 1 == 2)
check <- function(name, ok, detail = "", stop_on_fail = TRUE) {
  line <- sprintf("%-78s %s %s", name, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(line, "\n"); log_lines <<- c(log_lines, line)
  if (!isTRUE(ok)) { failed_checks <<- c(failed_checks, name); if (stop_on_fail) stop(name, call. = FALSE) }
  invisible(TRUE)
}
failed_checks <- character(0)

# Toy map (no restricted data) ----
# Nodes: 1 is a control college with neighbors 2, 3 (non-participants), 4 (treated),
# 5 (control). 6 participates; its only neighbor (2) does not, so it is an island.
toy_adj <- matrix(0, 6, 6); e <- rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 6), c(4, 5))
toy_adj[e] <- 1; toy_adj[e[, 2:1]] <- 1
toy <- list(weights = list(queen = list(adjacency = toy_adj)))
gt <- me_geometry(c(1, 4, 5, 6), toy, need_engine = FALSE)
zt <- c(0, 1, 0, 0)                     # participants 1, 4, 5, 6; only college 4 treated
et <- me_exposure(zt, gt)
check("V4 rule A: control with 4 neighbors (2 non-part., 1 treated, 1 control) = 0.25", isTRUE(all.equal(et[1], 0.25)))
check("V4 rule A: island exposure is 0 and its W_P row is zero", et[4] == 0 && all(gt$WP[4, ] == 0) && gt$islands == 1)
check("V4 rule A: W_P rows sum to 1 for non-islands", isTRUE(all.equal(rowSums(gt$WP)[1:3], c(1, 1, 1), check.attributes = FALSE)))
check("V11 binary form: 1{any treated neighbor}", identical(me_exposure_form(zt, gt, "binary"), c(1, 0, 1, 0)))
check("V11 count form: treated neighbors / dbar (dbar = 2)", isTRUE(all.equal(me_exposure_form(zt, gt, "count", dbar = 2), c(0.5, 0, 0.5, 0))))
check("V11 proportional form equals rule A", isTRUE(all.equal(me_exposure_form(zt, gt, "proportional"), et)))

# Real-map checks (frozen setup) ----
setup <- me_load_setup()
Wq <- setup$weights$queen$W
g58 <- me_geometry(1:58, setup)
check("N = 58 geometry: W_P equals production queen W, no islands", max(abs(g58$WP - Wq)) == 0 && g58$islands == 0)
set.seed(11); z58 <- sample(rep(0:1, 29))
check("V11 N = 58: proportional exposure equals production W z", max(abs(me_exposure(z58, g58) - as.vector(Wq %*% z58))) < 1e-14)
dbar <- mean(rowSums(setup$weights$queen$adjacency))
check("V11 N = 58: count form at full treatment equals degree / dbar", isTRUE(all.equal(me_exposure_form(rep(1, 58), g58, "count", dbar), rowSums(setup$weights$queen$adjacency) / dbar, check.attributes = FALSE)), sprintf("dbar = %.4f", dbar))

# Production DGP reproduced at N = 58 (rule A = production spillover)
cfg_prod <- rs_config("production"); u <- rs_units(cfg_prod)
for (reg in c("both", "control_only")) {
  v <- u[u$Year == "2019" & u$Model == "education" & u$Neighbor == "queen" & u$Design_ID == 9 &
         u$Rho == 0.5 & u$Gamma == 0.8 & u$Regime == reg, ][1, ]
  m <- rs_model(v, z58, setup, cfg_prod)
  e58 <- me_exposure(z58, g58); sp <- if (reg == "control_only") e58 * (1 - z58) else e58
  mine <- as.vector(solve(diag(58) - 0.5 * g58$WP) %*% (1 * z58 + 0.8 * 1 * sp))
  check(sprintf("V1/V11 N = 58 DGP mean and fitted Spill equal production (%s)", reg),
        max(abs(mine - m$mean)) < 1e-12 && max(abs(0.8 * sp - m$xm[, "Spill"])) < 1e-14)
}

# N = 58 frozen designs: same call (and draw, for the same seed) as production rs_draws
for (id in 1:9) {
  di <- me_design_inputs(g58, "2020", setup, list("random", 58, 1))
  a <- setup$annual[["2020"]]
  set.seed(100 + id); mine <- me_draw_design(id, di)
  set.seed(100 + id)
  prod <- get_application_designs(id, 1L, a$X, setup$weights$queen$nb, setup$coords, a$person_years,
    a$n_counties, setup$regions, block_id = setup$blocks[["2020"]],
    region_summary = rs_region_values(a, setup$regions, "mean_rank"))[, 1]
  check(sprintf("N = 58 design %d reuses the frozen production rule (identical draw)", id), identical(mine, prod))
}

# V3 Scale invariance ----
# beta = 0: Y(tau, sigma) = A (tau z + kappa tau e + sigma eps) is linear in (tau, sigma),
# and ML is scale-equivariant, so (0.5, 0.5) halves tau-hat and SE and keeps t.
P <- me_subset("random", 29, 3, "2018", setup); g <- me_geometry(P, setup)
di <- me_design_inputs(g, "2018", setup, list("random", 29, 3))
set.seed(5); z <- me_draw_design(9, di); sp <- me_exposure(z, g) * (1 - z)
xm <- cbind("(Intercept)" = 1, Z = z, Spill = sp); A <- solve(diag(g$n) - 0.5 * g$WP)
set.seed(6); E <- matrix(rnorm(g$n * 200), g$n)
d3 <- t(vapply(1:200, function(k) {
  f1 <- me_fit(as.vector(A %*% (1 * z + 0.8 * 1 * sp + 1 * E[, k])), xm, g$engine)
  f2 <- me_fit(as.vector(A %*% (0.5 * z + 0.8 * 0.5 * sp + 0.5 * E[, k])), xm, g$engine)
  c(f2["tau"] - 0.5 * f1["tau"], f2["se"] - 0.5 * f1["se"], f2["tau"] / f2["se"] - f1["tau"] / f1["se"])
}, numeric(3)))
check("V3 scale invariance: (0.5, 0.5) vs (1, 1) on 200 fits (tau-hat, SE, t)", max(abs(d3)) < 1e-6,
      sprintf("max |diff| = %.2e", max(abs(d3))))

# V5 Islands in the fitter ----
# Subsets with islands, including symmetric W_P (can.be.simmed() TRUE), which the
# unchanged sar_lag_setup() guard refuses; me_engine_setup() handles those and
# lagsarlm() (its own symmetric branch, original listw) is the reference.
isl <- list(); sym <- list(); s <- 0
while (length(isl) < 20 || length(sym) < 3) {
  s <- s + 1; Pi <- me_subset("random", 15, s, "2018", setup); gi <- me_geometry(Pi, setup)
  if (gi$islands > 0 && !is.null(gi$engine)) {
    if (gi$guard_bypass && length(sym) < 3) sym[[length(sym) + 1L]] <- c(gi, s = s)
    else if (!gi$guard_bypass && length(isl) < 20) isl[[length(isl) + 1L]] <- c(gi, s = s)
  }
}
n_bypass <- sum(vapply(1:2000, function(q) me_geometry(me_subset("random", 15, q, "2018", setup), setup)$guard_bypass, logical(1)))
check("V5 non-symmetric subset listw objects pass the can.be.simmed() guard",
      all(vapply(isl, function(x) !spatialreg::can.be.simmed(x$listw), logical(1))),
      sprintf("islands per subset %d-%d; symmetric W_P at N = 15: %d / 2000 subsets", min(sapply(isl, `[[`, "islands")), max(sapply(isl, `[[`, "islands")), n_bypass))
check("V5 symmetric W_P subsets fail the guard and get an engine via me_engine_setup()",
      all(vapply(sym, function(x) spatialreg::can.be.simmed(x$listw) && !is.null(x$engine), logical(1))))
if (Sys.getenv("SKIP_LAGSARLM") != "1") {
  #' lean vs lagsarlm(zero.policy = TRUE) on k fits of one subset
  #' @param gi Geometry. @param k Fits. @param seed Seed.
  #' @return Matrix of |dtau| and |dSE|.
  #' @examples
  #' # compare_subset(isl[[1]], 10, 1)
  compare_subset <- function(gi, k, seed) {
    dii <- me_design_inputs(gi, "2018", setup, list("random", 15, gi$s))
    set.seed(seed); zi <- me_draw_design(9, dii); spi <- me_exposure(zi, gi)
    xmi <- cbind("(Intercept)" = 1, Z = zi, Spill = spi); Ai <- solve(diag(gi$n) - 0.4 * gi$WP)
    t(vapply(seq_len(k), function(i) {
      y <- as.vector(Ai %*% (zi + 0.5 * spi + rnorm(gi$n)))
      lean <- fit_sar_lag(y, xmi, gi$engine)
      ref <- spatialreg::lagsarlm(Y ~ Z + Spill, data = data.frame(Y = y, Z = zi, Spill = spi),
                                  listw = gi$listw, zero.policy = TRUE, quiet = TRUE)
      c(abs(lean$coefficients[["Z"]] - coef(ref)[["Z"]]), abs(lean$rest.se[["Z"]] - ref$rest.se[["Z"]]))
    }, numeric(2)))
  }
  d1 <- do.call(rbind, lapply(seq_along(isl), function(q) compare_subset(isl[[q]], 50, q)))
  d2 <- do.call(rbind, lapply(seq_along(sym), function(q) compare_subset(sym[[q]], 100, 100 + q)))
  check("V5 fit_sar_lag matches lagsarlm(zero.policy = TRUE) with islands (1,000 fits)",
        nrow(d1) >= 1000 && max(d1[, 1]) < 1e-6, sprintf("max |dtau| = %.2e, max |dSE| = %.2e", max(d1[, 1]), max(d1[, 2])))
  check("V5 symmetric W_P subsets match lagsarlm's symmetric branch (300 fits)",
        max(d2[, 1]) < 1e-6, sprintf("max |dtau| = %.2e, max |dSE| = %.2e", max(d2[, 1]), max(d2[, 2])))
}

# V6 Designs on subsets ----
Ns <- c(15, 20, 25, 29, 30, 35, 40, 45, 50)
v6 <- list(); bad <- character(0)
for (q in 1:1000) {
  N <- Ns[(q - 1) %% length(Ns) + 1]; yr <- as.character(2018 + (q %% 4))
  Pq <- me_subset("random", N, 1000 + q, yr, setup); gq <- me_geometry(Pq, setup, need_engine = FALSE)
  dq <- me_design_inputs(gq, yr, setup, list("random", N, 1000 + q))
  me_seed("v6", q)
  Zq <- lapply(1:9, function(id) me_draw_design(id, dq))
  tab_r <- tabulate(dq$regions)
  ib <- Zq[[4]]; adj_t <- as.vector(gq$AP %*% ib)
  cb2 <- { me_seed("v6b", q); me_draw_design(1, dq) }
  # Budget expected for each rule-based saturation/blocking design
  sat_ok <- function(zz, levels) { n_r <- as.numeric(table(dq$regions)); tr <- sapply(split(zz, dq$regions), sum)
    all(vapply(seq_along(n_r), function(i) any(abs(tr[i] - round(n_r[i] * levels)) == 0), logical(1))) }
  blk_ok <- all(sapply(split(Zq[[5]], dq$blocks), function(b) sum(b) == round(length(b) / 2)))
  v6[[q]] <- data.frame(N = N, floor_ok = all(sapply(Zq[c(2, 6, 7, 9)], sum) == floor(N / 2)),
    ib_indep = sum(ib * adj_t) == 0, ib_maximal = all(adj_t[ib == 0] > 0),
    ib_islands_treated = all(ib[rowSums(gq$AP) == 0] == 1),
    n_regions = length(unique(dq$regions)), min_region = min(tab_r), below3 = sum(tab_r < 3),
    sat_ok = sat_ok(Zq[[3]], c(0.2, 0.4, 0.6, 0.8)) && sat_ok(Zq[[8]], c(0.2, 0.4, 0.6, 0.8)),
    blk_ok = blk_ok, cb_det = identical(Zq[[1]], cb2), cb_half = sum(Zq[[1]]) <= N / 2,
    islands = gq$islands, treated_cb = sum(Zq[[1]]), treated_ib = sum(ib),
    treated_sat = sum(Zq[[3]]), treated_igs = sum(Zq[[8]]), treated_blk = sum(Zq[[5]]))
}
v6 <- do.call(rbind, v6)
check("V6 fixed-budget designs (2, 6, 7, 9) treat floor(N/2) on 1,000 subsets", all(v6$floor_ok))
check("V6 Isolation Buffer: no two adjacent treated participants", all(v6$ib_indep))
check("V6 Isolation Buffer: maximal (every control has a treated participant neighbor)", all(v6$ib_maximal))
check("V6 Isolation Buffer: islands always treated", all(v6$ib_islands_treated))
check("V6 saturation designs use their regional rules; 2x2 blocking balances each block", all(v6$sat_ok) && all(v6$blk_ok))
check("V6 checkerboard deterministic given the subset (and treats <= N/2)", all(v6$cb_det) && all(v6$cb_half))
check("V6 four regions with >= 2 members each", all(v6$n_regions == 4) && all(v6$min_region >= 2),
      sprintf("min region size %d; subsets with a region < 2: %d, < 3: %d / 1000", min(v6$min_region),
              sum(v6$min_region < 2), sum(v6$below3 > 0)), stop_on_fail = FALSE)
v6_by_N <- aggregate(cbind(below2 = min_region < 2, below3 = below3 > 0, islands, treated_cb, treated_ib, treated_sat, treated_igs, treated_blk) ~ N,
                     data = v6, FUN = mean)
print(v6_by_N)

# V7 Informative sampling ----
for (N in c(20, 29)) {
  rk <- setup$annual[["2018"]]$rank; pi <- me_incidence_pi(rk, N)
  inc <- numeric(58)
  for (s in 1:20000) { ps <- me_subset("incidence3x", N, s, "2018", setup); inc[ps] <- inc[ps] + 1 }
  inc <- inc / 20000
  ext <- inc[rk == 58] / inc[rk == 1]
  top <- rk > quantile(rk, 0.9); bot <- rk <= quantile(rk, 0.1)
  dec <- mean(inc[top]) / mean(inc[bot]); dec_target <- mean(pi[top]) / mean(pi[bot])
  check(sprintf("V7 N = %d: realized inclusion probabilities match pi (max |diff|)", N),
        max(abs(inc - pi)) < 4 * sqrt(max(pi * (1 - pi)) / 20000) + 1e-3, sprintf("%.4f", max(abs(inc - pi))))
  check(sprintf("V7 N = %d: highest/lowest-burden inclusion ratio about 3", N), abs(ext - 3) < 0.25,
        sprintf("extreme %.2f; top/bottom decile %.2f (design value %.2f)", ext, dec, dec_target))
}

dir.create(file.path(.me_root, "verification"), recursive = TRUE, showWarnings = FALSE)
writeLines(c(sprintf("Run %s", format(Sys.time())), log_lines),
           file.path(.me_root, "verification", "test_methods_extensions.txt"))
write.csv(v6_by_N, file.path(.me_root, "verification", "v6_designs_by_N.csv"), row.names = FALSE)
if (length(failed_checks)) {
  cat("FAILED CHECKS:\n", paste(failed_checks, collapse = "\n"), "\n")
  writeLines(c(readLines(file.path(.me_root, "verification", "test_methods_extensions.txt")),
               paste("FAILED:", failed_checks)), file.path(.me_root, "verification", "test_methods_extensions.txt"))
  quit(status = 1)
}
cat("All methods-extension checks passed.\n")
