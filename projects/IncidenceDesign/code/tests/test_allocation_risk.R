# ============================================================
# Script: test_allocation_risk.R
# Purpose: Verify nested allocation-risk summaries, seeds, and actual ML fitting.
# Author: Andrew Walther
# Created: 2026-09-27
# Dependencies: spdep, spatialreg, dplyr, digest; sources code/16 (define-only)
# ============================================================
# Run from code/ with the three BLAS thread variables set to1.
allocation_define_only <- TRUE
allocation_script_dir <- normalizePath(".")
source("16_allocation_risk.R")

#' Stop on a failed behavioral assertion
#' @param name Description of intended behavior.
#' @param ok Logical scalar.
#' @return Invisibly NULL; raises error on failure.
#' @examples
#' check("identity", 1 == 1)
check <- function(name, ok) {
  cat(sprintf("%-75s %s\n", name, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("Check failed: ", name, call. = FALSE)
  invisible(NULL)
}

# Weighted duplicate cache: known arithmetic ----
a <- data.frame(Frequency = c(2, 1), MSE = c(1, 3), SE_MSE = c(0.2, 0.3),
  MSE_Half1 = c(0.9, 2.8), MSE_Half2 = c(1.1, 3.2), Complete = TRUE)
s <- allocation_distribution(a)
noise_expected <- (2 * (1 - 2 / 3) * 0.04 + (1 - 1 / 3) * 0.09) / 2
check("draw-frequency mean is5/3 (not the unweighted unique mean2)", abs(s$Mean_MSE - 5 / 3) < 1e-14)
check("raw variance matches expanded draws1,1,3", abs(s$Variance_Raw - var(c(1, 1, 3))) < 1e-14)
check("cache covariance correction matches exact frequency formula", abs(s$Variance_MC_Noise - noise_expected) < 1e-14)
check("signed corrected variance equals raw minus MC variance", abs(s$Variance_Corrected - (4 / 3 - noise_expected)) < 1e-14)
check("weighted q90 matches expanded draw quantile", s$Q90_Estimated == unname(quantile(c(1, 1, 3), 0.9)))
check("cached mean noise SE includes squared frequency weights", abs(s$SE_Mean_MSE_Noise^2 - (4 / 9 * 0.04 + 1 / 9 * 0.09)) < 1e-14)
check("joint mean SE agrees with corrected variance plus cache noise", abs(s$SE_Mean_MSE_Joint^2 - (s$Variance_Corrected / 3 + s$SE_Mean_MSE_Noise^2)) < 1e-14)
singleton <- a[1, ]; singleton$Frequency <- 100L
ss <- allocation_distribution(singleton)
check("deterministic assignment allocation variance0 despite MSE1", ss$Variance_Raw == 0 && ss$Variance_Corrected == 0 && ss$Mean_MSE == 1)
check("singleton cache mean SE retains outcome uncertainty", ss$SE_Mean_MSE_Noise == singleton$SE_MSE)
negative <- a; negative$MSE <- c(1, 1.001)
sn <- allocation_distribution(negative)
check("negative variance correction retained; SD missing instead of truncation", sn$Variance_Corrected < 0 && is.na(sn$SD_Corrected))
f <- data.frame(replicate = 1:4, error = c(-1, 1, -2, 2), se = 1, alias = FALSE, messages = "")
sf <- allocation_fit_summary(f, 1)
check("nonzero error kept even when bias is0", sf$Bias == 0 && sf$MSE == 2.5)
f$se[1] <- NA_real_; sf <- allocation_fit_summary(f, 1)
check("CI failure retains finite estimate and flags incomplete", sf$N_Valid_Est == 4 && sf$N_Valid_CI == 3 && sf$MSE == 2.5 && !sf$Complete)
bad <- a; bad$Complete[1] <- FALSE
check("incomplete block cannot silently drop an allocation", is.na(allocation_distribution(bad)$Mean_MSE))

# Seed prefixes and actual design behavior ----
p <- allocation_parameters("smoke"); grid <- build_spatial_grid(p$grid_dim)
setup <- sar_lag_setup(grid$listw_queen); cfg <- p$configs[[1]]
X <- allocation_surface(cfg, 1, grid, p)
u <- data.frame(config = 1, surface = 1, rho = 0.5, gamma = 0.5, regime = "both", design = 9)
z4 <- allocation_draws(u, cfg, X, grid, p)
p8 <- p; p8$n_allocations <- 8L
z8 <- allocation_draws(u, cfg, X, grid, p8)
check("increasing assignment reps preserves earlier draws", identical(z4, z8[, 1:4]))
allocation_seed("X", "iid", 0, 1)
check("existing keyed X surface reproduced exactly", identical(X, generate_incidence_iid(100, 1)[, 1]))
noise <- allocation_noise(list("first-allocation"), 1, 100, p)
check("noise batch reproduces independently of outside RNG state", {
  set.seed(14); runif(10); identical(noise, allocation_noise(list("first-allocation"), 1, 100, p))
})
check("distinct allocations get distinct noise streams", !identical(noise, allocation_noise(list("second-allocation"), 1, 100, p)))

# Actual oracle ML fits, extension and exact engine comparison ----
p$n_allocations <- 2L
out <- allocation_block(u, cfg, X, grid, setup, p)
p16 <- p; p16$n_outcomes <- 16L
extended <- allocation_block(u, cfg, X, grid, setup, p16, out)
fresh <- allocation_block(u, cfg, X, grid, setup, p16)
check("actual ML fits finite and nonconstant", all(out$allocations$Complete) && all(out$allocations$MSE > 0) && length(unique(out$fits[[1]]$error)) == 8)
check("outcome extension cache gives exactly fresh extended fits", identical(extended$fits, fresh$fits))
check("outcome extension retains original8 errors and SEs", identical(out$fits[[1]], extended$fits[[1]][1:8, ]))
rev_p <- p; rev_p$n_allocations <- 4L
out4 <- allocation_block(u, cfg, X, grid, setup, rev_p, out)
check("assignment extension retains existing allocation outcomes", all(vapply(names(out$fits), function(h) identical(out$fits[[h]], out4$fits[[h]]), logical(1))))
u$design <- 1
cb <- allocation_block(u, cfg, X, grid, setup, p)
check("Checkerboard drawn and fit once, allocation variance zero", nrow(cb$draws) == 1 && cb$summary$Variance_Corrected == 0 && cb$allocations$MSE > 0)
u$design <- 2
hif <- allocation_block(u, cfg, X, grid, setup, rev_p)
check("continuous-X HIF cache retains all4 draw frequencies", nrow(hif$allocations) == 1 && hif$allocations$Frequency == 4 && nrow(hif$draws) == 4 && hif$allocations$Structural_Singleton)

for (d in c(1, 6, 9)) {
  u$design <- d; z <- allocation_draws(u, cfg, X, grid, p)[, 1]
  spill <- 0.5 * as.vector(grid$W_queen %*% z)
  allocation_seed("allocation-risk-test-equivalence", d)
  y <- as.vector(solve(diag(100) - 0.5 * grid$W_queen, z + spill + X + rnorm(100)))
  xm <- cbind("(Intercept)" = 1, Z = z, Spill = spill, X = X)
  lean <- fit_one_lag_model(y, xm, setup, engine = "lean")
  reference <- fit_one_lag_model(y, xm, setup, grid$listw_queen, "lagsarlm")
  check(paste("actual ML engine equivalence design", d), abs(lean$tau - reference$tau) < 1e-6 && abs(lean$se - reference$se) < 1e-6)
}
cat("All allocation-risk tests passed.\n")
