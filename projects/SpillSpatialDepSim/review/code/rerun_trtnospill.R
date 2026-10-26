# ============================================================
# Script: review/code/rerun_trtnospill.R
# Purpose: Re-run one Control Only (TrtNoSpill) scenario of SpatialSim_3x3.Rmd /
#          SpatialSim_3x4.Rmd with the case-1 loop bound either as published (1:8)
#          or corrected (1:n_districts). Same seed and RNG order as the originals.
#          Mode trtspill re-runs the Intervention & Control (case 2) scenario instead,
#          used to regenerate the one published TrtSpill CSV that came from a
#          different random stream (3x3, psi = 0.8, rho = 0.01; see notes/2026-10-25-correction-plan.md).
# Usage:   Rscript rerun.R <3x3|3x4> <psi> <rho> <bug|fixed|trtspill>
# Output:  review/results/beta_mse_TrtNoSpill_fixed/ (mode fixed), review/results/validation/ (mode bug),
#          review/results/beta_mse_TrtSpill_rerun/ (mode trtspill)
# Author:  Andrew Walther (review written with Claude Code)
# Created: 2026-10-22
# Dependencies: spdep, spatialreg
# Note:    writes only inside review/; never touches code/ or results/ of the project.
# ============================================================
suppressMessages({library(spdep); library(spatialreg)})
a <- commandArgs(TRUE); g <- a[1]; psi <- as.numeric(a[2]); rho <- as.numeric(a[3]); mode <- a[4]
if (g == "3x4") { nc <- 4; k <- 6; N <- 10; ppi <- 240; gw <- 400
  nb <- list(c(2,5),c(1,3,6),c(2,4,7),c(3,8),c(1,6,9),c(2,5,7,10),c(3,6,8,11),c(4,7,12),c(5,10),c(6,9,11),c(7,10,12),c(8,11))
} else { nc <- 3; k <- 4; N <- 50; ppi <- 180; gw <- 300
  nb <- list(c(2,4),c(1,3,5),c(2,6),c(1,5,7),c(2,4,6,8),c(3,5,9),c(4,8),c(5,7,9),c(6,8)) }
n <- 3 * nc; alpha <- 0.2; beta <- 1; sdv <- 0.1
loop_max <- if (mode == "bug") 8 else n          # the only difference between the two modes
A <- matrix(0, n, n); for (i in 1:n) A[i, nb[[i]]] <- 1   # district Rook adjacency
set.seed(2024)
combos <- combn(1:n, k, simplify = FALSE)
est <- matrix(NA, N, length(combos))
for (ci in seq_along(combos)) {
  trt <- combos[[ci]]
  # z_d = 1 if control district d has a treated Rook neighbour among districts 1..loop_max
  # (mode trtspill: z_d = 1 for any district, treated or control, with a treated Rook neighbour)
  zd <- if (mode == "trtspill") sapply(1:n, function(d) as.integer(any(nb[[d]] %in% trt))) else
    sapply(1:n, function(d) as.integer(!(d %in% trt) && any(nb[[d]][nb[[d]] <= loop_max] %in% trt)))
  for (r in 1:N) {
    x <- sample(1:gw - 1, ppi, replace = TRUE); y <- sample(1:300 - 1, ppi, replace = TRUE)
    eps <- rnorm(ppi, 0, sdv)                      # same RNG order as the original scripts
    d <- (2 - (y %/% 100)) * nc + (x %/% 100) + 1  # row-major district, 1 = top-left
    xi <- as.integer(d %in% trt); z <- zd[d]
    W <- A[d, d]                                   # point i ~ point j iff their districts are Rook neighbours
    yv <- solve(diag(ppi) - rho * W, alpha + beta * xi + psi * z + eps)
    m <- lagsarlm(yv ~ xi + z, data = data.frame(yv, xi, z), listw = mat2listw(W, style = "W"))
    est[r, ci] <- coef(m)[3]                       # beta-hat (coef order: rho, intercept, xi, z)
  }
}
out <- t(apply(est, 2, function(e) c(bias = mean(e) - beta, bias.abs = abs(mean(e) - beta), variance = var(e), MSE = mean((e - beta)^2))))
rownames(out) <- seq_along(combos)
outdir <- switch(mode, bug = "results/validation", fixed = "results/beta_mse_TrtNoSpill_fixed",
                 trtspill = "results/beta_mse_TrtSpill_rerun")
dir.create(outdir, showWarnings = FALSE)
write.csv(out, sprintf("%s/beta_mse_results_%s_%s_A02_B1_P%02d_R%s.csv", outdir, ifelse(mode == "trtspill", "TrtSpill", "TrtNoSpill"), g, round(psi * 10), ifelse(rho == 0, "00", "001")))
