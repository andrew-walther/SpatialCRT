# ============================================================
# Script: review/code/rho_scale_diagnostic.R
# Purpose: Check what the published rho = 0.01 means on the row-standardized scale
#          the model is fitted on, by re-running the 2x4 data-generating process.
# Author: Andrew Walther
# Created: 2026-10-28
# Dependencies: spatialreg, spdep, parallel
# Note: documentation only; it does not change any reported result. See
#       review/notes/2026-10-28-rho-scale-note.md. Runs in about 2 minutes on 8 cores.
# ============================================================
# Replicates the frozen 2x4 data-generating process (code/SpatialSim_NC_DOC.Rmd, lines 380-471):
# 160 subjects, binary subject-level W (every subject in a rook-neighboring district), outcomes
# y = (I - rho W)^{-1}(alpha + beta x + psi z + e), fit by lagsarlm with row-standardized W.
# Settings (psi = 0.5, Control Clusters Only, N = 20 per allocation, all 70 allocations):
#   A = published DGP (binary W, rho = 0.01); B = rho = 0;
#   C = row-standardized W in the DGP with rho = 0.5 (the same W as the fit).
# BSS allocations are combn indices 28 and 43, as in the frozen script.
suppressMessages({library(spatialreg); library(spdep); library(parallel)})
set.seed(20261028)
nb <- list(c(2,3), c(1,4), c(1,4,5), c(2,3,6), c(3,6,7), c(4,5,8), c(5,8), c(6,7))
combos <- combn(1:8, 4, simplify = FALSE)
alpha <- 0.2; beta <- 1; psi <- 0.5; sdv <- 0.1; N <- 20; n <- 160

one_fit <- function(trt, rho, wtype) {
  d <- sample(1:8, n, replace = TRUE)            # district of each subject (uniform points)
  x <- as.integer(d %in% trt)
  z <- as.integer(x == 0 & sapply(d, function(k) any(nb[[k]] %in% trt)))  # controls only
  Wb <- outer(d, d, function(a, b) mapply(function(a, b) b %in% nb[[a]], a, b)) * 1
  Wr <- Wb / pmax(rowSums(Wb), 1)
  Wg <- if (wtype == "binary") Wb else Wr
  y <- solve(diag(n) - rho * Wg, alpha + beta * x + psi * z + rnorm(n, 0, sdv))
  m <- lagsarlm(y ~ x + z, data = data.frame(y, x, z), listw = mat2listw(Wb, style = "W"), quiet = TRUE)
  c(beta = unname(coef(m)["x"]), rho = unname(m$rho))
}
run <- function(rho, wtype) {
  res <- mclapply(seq_along(combos), function(i) {
    est <- replicate(N, one_fit(combos[[i]], rho, wtype))
    c(mse = mean((est["beta", ] - beta)^2), rhohat = mean(est["rho", ]))
  }, mc.cores = max(1, detectCores() - 1))
  do.call(rbind, res)
}
bss <- c(28, 43)
for (s in list(list("A published: binary W, rho=0.01", 0.01, "binary"),
               list("B rho=0", 0, "binary"),
               list("C row-std W in DGP, rho=0.5", 0.5, "rowstd"))) {
  r <- run(s[[2]], s[[3]])
  cat(sprintf("%-34s SRS mean MSE %.3f | BSS mean MSE %.3f | mean rho-hat %.2f | BSS beats %.0f%% of allocations\n",
              s[[1]], mean(r[, "mse"]), mean(r[bss, "mse"]), mean(r[, "rhohat"]),
              100 * mean(r[, "mse"] > mean(r[bss, "mse"]))))
}
