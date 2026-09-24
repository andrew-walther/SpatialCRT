# ============================================================
# Script: test_simulation_revision.R
# Purpose: Tests for the 2026-09 simulation revision (simulation-revision-plan.md B3).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-24
# Dependencies: spatialreg, spdep, dplyr, digest, parallel; sources code/05 (define-only)
# ============================================================
#
# Run from code/:
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript tests/test_simulation_revision.R
# Exits non-zero on the first failed check. Set SKIP_EQUIVALENCE=1 to skip the
# ~7-minute estimator-equivalence test (it is re-run otherwise; its own output goes
# to results/estimator_validation/).
#
# The surface-level robustness check (paired test over 50 config x surface units)
# needs full results; it is run by the pilot/full analysis, not here.

sim_profile <- "full"
sim_define_only <- TRUE
sim_script_dir <- normalizePath(".")   # run from code/
suppressMessages(source("05_run_simulation.R"))

check <- function(name, ok) {
  cat(sprintf("  %-66s %s\n", name, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("Check failed: ", name, call. = FALSE)
}

# 1. Designs (M3) ----
cat("\n[1] Designs with tied incidence (M3)\n")
set.seed(1)
X_tied <- sample(rep(c(0.1, 0.3, 0.5, 0.7, 0.9), each = 20))   # 5 values, heavy ties
Z2 <- get_designs(2, 1000, N, X_tied, grid_obj$nb_queen, grid_obj$coords)
check("High Incidence Focus treats exactly 50 on tied X (1,000 draws)",
      all(colSums(Z2) == 50))
check("High Incidence Focus treats every cluster above the tied median value",
      all(Z2[X_tied > 0.5, ] == 1) && all(Z2[X_tied < 0.5, ] == 0))
check("High Incidence Focus breaks the tie differently across draws",
      ncol(unique(Z2, MARGIN = 2)) > 900)

# Tie assignment vs grid position: the 20 clusters tied at the median value fill 10
# treated slots, so each must be treated with probability 0.5 whatever its grid
# index. The pre-revision rule gave probabilities of exactly 0 or 1 here.
# Bound: 0.07 ~ 4.4 binomial SDs over 1,000 draws.
tied <- which(X_tied == 0.5)
p_tied <- rowMeans(Z2[tied, ])
check("HIF tie-break independent of grid position (all tied p = 0.5 ± 0.07)",
      all(abs(p_tied - 0.5) < 0.07))

for (d in c(6, 7)) {
  Zd <- get_designs(d, 1000, N, X_tied, grid_obj$nb_queen, grid_obj$coords)
  n_strata <- if (d == 6) 4 else 2
  # Reconstruct each draw's strata is impossible after the fact, so check the
  # consequence: treated counts per incidence level are balanced in expectation,
  # and the within-tie treatment probability is flat in grid index.
  lvl_rate <- tapply(rowMeans(Zd), X_tied, mean)
  expected <- if (d == 6) 12 / 25 else 0.5
  check(sprintf("Design %d: treated per draw = %d", d, if (d == 6) 48 else 50),
        all(colSums(Zd) == if (d == 6) 48 else 50))
  check(sprintf("Design %d: treatment rate at every tied incidence level ~ %.2f (±0.03)",
                d, expected),
        all(abs(lvl_rate - expected) < 0.03))
  pr <- rowMeans(Zd)
  check(sprintf("Design %d: every cluster's treatment rate ~ %.2f ± 0.07 (no grid-position tie-break)",
                d, expected),
        all(abs(pr - expected) < 0.07))
}
# Strata are equal-sized rank strata: verify directly from the stratum rule
r <- random_tie_rank(X_tied)
check("ntile(random_tie_rank(X), 4) gives 4 strata of 25",
      all(table(dplyr::ntile(r, 4)) == 25))
check("ntile(random_tie_rank(X), 2) gives 2 strata of 50",
      all(table(dplyr::ntile(r, 2)) == 50))
check("Checkerboard is the only deterministic design",
      identical(vapply(1:8, is_design_deterministic, logical(1)), c(TRUE, rep(FALSE, 7))))

# 2. Design k is drawn from X[, k] (M1) ----
cat("\n[2] Matched surfaces (M1)\n")
X_sp <- generate_surfaces(inc_configs[[3]])          # spatial rhoX = 0.50, K = 10
Zm <- draw_assignments(inc_configs[[3]], "queen", 0.2, 0.6, "both", 2, X_sp)
cors <- sapply(seq_len(n_surfaces), function(j) sapply(seq_len(n_surfaces), function(k) {
  cor(rowMeans(Zm[, (k - 1) * n_design_draw + seq_len(n_design_draw)]), X_sp[, j])
}))                                                  # rows: design surface k, cols: X_j
check("High Incidence Focus draws for surface k track X_k (diag cor > 0.8)",
      all(diag(cors) > 0.8))
check("... and not other surfaces (|off-diag cor| < 0.4)",
      all(abs(cors[row(cors) != col(cors)]) < 0.4))
check("Surfaces are distinct columns (max |cor(X_j, X_k)| < 0.5)",
      max(abs(cor(X_sp)[upper.tri(diag(n_surfaces))])) < 0.5)
X_p <- generate_surfaces(inc_configs[[5]])           # Poisson rhoX = 0.50, P = 100,000
check("Poisson surfaces: median distinct X values per surface >= 55 (M2)",
      median(apply(X_p, 2, function(x) length(unique(x)))) >= 55)

# 3. Determinism: sequential = parallel = any order (M4) ----
cat("\n[3] Seeding (M4)\n")
check("Surfaces regenerate identically", identical(generate_surfaces(inc_configs[[5]]), X_p))
n_surfaces <- 2; n_design_draw <- 3; gamma_vals <- c(0.5, 0.8); spill_types <- "both"
true_tau_vals <- c(1, 2); design_ids <- c(1, 2, 4)
surf_small <- lapply(inc_configs, generate_surfaces)
small_units <- expand.grid(cfg_index = c(1, 4), nb_type = c("rook", "queen"),
                           rho = c(0, 0.5), stringsAsFactors = FALSE)
run_i <- function(i) {
  u <- small_units[i, ]
  run_unit(u$cfg_index, u$nb_type, u$rho, surf_small[[u$cfg_index]])
}
seq_res <- lapply(seq_len(nrow(small_units)), run_i)
par_res <- parallel::mclapply(seq_len(nrow(small_units)), run_i, mc.cores = 4,
                              mc.preschedule = FALSE)
rev_res <- rev(lapply(rev(seq_len(nrow(small_units))), run_i))
check("Sequential and parallel results identical", identical(seq_res, par_res))
check("Reverse-order results identical", identical(seq_res, rev_res))
check("Unit result does not depend on RNG state before the call", {
  set.seed(999); a <- run_i(3); runif(10); b <- run_i(3); identical(a, b)
})
s1 <- seq_res[[1]]$scen
check("Z shared across tau (same Mean_Treated at tau = 1 and 2)", {
  a <- s1[s1$True_Tau == 1, c("Design", "Gamma", "Estimator", "Mean_Treated")]
  b <- s1[s1$True_Tau == 2, c("Design", "Gamma", "Estimator", "Mean_Treated")]
  rownames(a) <- rownames(b) <- NULL
  identical(a, b)
})

# 4. Aliasing and rank-deficiency flags (M6/M7), gamma >= 0.5 ----
cat("\n[4] Aliasing flags (M6/M7)\n")
n_surfaces <- 2; n_design_draw <- 3; gamma_vals <- c(0.5, 0.8)
spill_types <- c("control_only", "both"); true_tau_vals <- 1; design_ids <- 1:8
flag <- dplyr::bind_rows(lapply(c("rook", "queen"), function(nb) {
  run_unit(2, nb, 0.2, surf_small[[2]])$scen
}))
cb_rook <- flag$Design == "Design 1" & flag$Neighbor_Type == "rook"
fits_per <- n_surfaces * n_design_draw
orc <- flag$Estimator == "oracle"
check("Oracle: aliased in 100% of Checkerboard x rook fits",
      all(flag$N_Aliased[orc & cb_rook] == fits_per))
check("Oracle: aliased in 0% of all other fits", all(flag$N_Aliased[orc & !cb_rook] == 0))
check("Non-oracle: never aliased", all(flag$N_Aliased[!orc] == 0))
check("Z_WZ_rank_deficient TRUE only for Checkerboard x rook (both estimators)",
      identical(flag$Z_WZ_rank_deficient, cb_rook))
check("No other warnings", all(flag$N_Warn == 0))
check("High Incidence Focus Mean_Treated = 50",
      all(flag$Mean_Treated[flag$Design == "Design 2"] == 50))

# 5. Known answer: correctly specified models are ~unbiased with ~nominal coverage ----
cat("\n[5] Known answer (Balanced Quartiles, queen, spatial rhoX = 0.20, rho = 0.20)\n")
n_surfaces <- 10; n_design_draw <- 25; spill_types <- "both"; true_tau_vals <- 1
design_ids <- 6
X_ka <- generate_surfaces(inc_configs[[2]])
known <- function(g) { gamma_vals <<- g; run_unit(2, "queen", 0.2, X_ka)$scen }
ka <- rbind(known(0.5)[1, ], known(0)[2, ])          # oracle at gamma 0.5; non-oracle at 0
stopifnot(ka$Estimator == c("oracle", "nonoracle"))
print(ka[, c("Estimator", "Gamma", "Bias", "SE_Bias", "Coverage", "SE_Coverage", "N_Valid_Est")],
      row.names = FALSE)
for (i in 1:2) {
  lab <- sprintf("%s, gamma = %.1f", ka$Estimator[i], ka$Gamma[i])
  check(sprintf("%s: |bias| < 3 MC SE", lab), abs(ka$Bias[i]) < 3 * ka$SE_Bias[i])
  # Coverage SE from surface means can be ~0 when every surface covers ~95%;
  # floor it at the binomial SE over all fits so the band is not degenerate.
  se_cov <- max(ka$SE_Coverage[i], sqrt(0.95 * 0.05 / ka$N_Valid_Est[i]))
  check(sprintf("%s: |coverage - 0.95| < 3 MC SE", lab), abs(ka$Coverage[i] - 0.95) < 3 * se_cov)
}

# 6. Estimator equivalence (lean engine vs lagsarlm) ----
if (Sys.getenv("SKIP_EQUIVALENCE") != "1") {
  cat("\n[6] Estimator equivalence (>= 5,000 fits; ~7 min)\n")
  source("tests/test_estimator_equivalence.R", local = new.env())
} else {
  cat("\n[6] Estimator equivalence SKIPPED (SKIP_EQUIVALENCE=1); last result in",
      "results/estimator_validation/equivalence_summary.txt\n")
}

cat("\nAll simulation-revision tests passed.\n")
