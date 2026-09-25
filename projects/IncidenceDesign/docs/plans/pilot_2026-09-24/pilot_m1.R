# ============================================================
# Script: pilot_m1.R
# Purpose: Scratch pilot (evidence for simulation-revision-plan.md) comparing the
#          current incidence structure against M1 matched surfaces.
# Author: Andrew Walther
# Created: 2026-09-24
# Dependencies: spatialreg, parallel; sources code/01-03
# Note: run once from a session scratchpad; the saveRDS path below points there and
#       is not part of the repo. Output as printed is in pilot_m1_output.txt.
#       This is not the plan's B4 pilot, which runs on the revised code.
# ============================================================
# Scratch pilot: current incidence structure vs. M1(A) (design and outcome see the same X).
# Slice: rho = 0.2, gamma = 0.7, tau = 1, both spill types, both nb types, 6 designs.
# Incidence: spatial rhoX = 0.5 (pop n/a) and Poisson rhoX = 0.5 at pop 1,000 (current) vs 100,000.
suppressMessages({
  setwd("/Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign/code")
  source("01_spatial_setup.R"); source("02_incidence_generation.R"); source("03_designs.R")
  library(spatialreg); library(parallel)
})
g <- build_spatial_grid(grid_dim = 10); N <- 100; I_mat <- diag(N)
rho <- 0.2; gam <- 0.7; tau <- 1; nD <- 25; nO <- 10

inc_cfgs <- list(
  list(lab = "spatial 0.5",           mode = "spatial", pop = NA),
  list(lab = "poisson 0.5, pop 1e3",  mode = "poisson", pop = 1000),
  list(lab = "poisson 0.5, pop 1e5",  mode = "poisson", pop = 1e5)
)
grid <- expand.grid(ic = seq_along(inc_cfgs), nb = c("rook", "queen"),
                    spill = c("both", "control_only"), d = c(1, 2, 4, 5, 6, 8),
                    struct = c("current", "M1A"), stringsAsFactors = FALSE)

run_one <- function(i) {
  s <- grid[i, ]; ic <- inc_cfgs[[s$ic]]
  set.seed(1000 + s$ic)  # same incidence across designs/nb/spill/struct within a config
  X <- if (ic$mode == "spatial") generate_incidence("spatial", N, nO, g$W_queen, 0.5)
       else generate_incidence("poisson", N, nO, g$W_queen, 0.5, pop_per_cluster = ic$pop)
  if (s$struct == "M1A") X <- matrix(X[, 1], N, nO)   # outcome/analysis use the design's X
  W  <- if (s$nb == "rook") g$W_rook else g$W_queen
  lw <- if (s$nb == "rook") g$listw_rook else g$listw_queen
  nbl <- if (s$nb == "rook") g$nb_rook else g$nb_queen
  inv <- solve(I_mat - rho * W)
  set.seed(2000 + s$ic + 10 * (s$nb == "rook"))
  E <- matrix(rnorm(N * nO), N, nO)
  set.seed(3000 + i)
  Zm <- get_designs(s$d, nD, N, X[, 1], nbl, g$coords)
  est <- lo <- hi <- c(); alias <- 0
  for (j in seq_len(nD)) {
    Z <- Zm[, j]; WZ <- as.vector(W %*% Z)
    sp <- if (s$spill == "control_only") gam * WZ * (1 - Z) else gam * WZ
    for (k in seq_len(nO)) {
      Y <- as.vector(inv %*% (tau * Z + sp + X[, k] + E[, k]))
      df <- data.frame(Y = Y, Z = Z, Spill = sp, X = X[, k])
      f <- tryCatch(suppressWarnings(lagsarlm(Y ~ Z + Spill + X, df, lw, quiet = TRUE)), error = function(e) NULL)
      if (is.null(f) || !"Z" %in% names(coef(f))) { est <- c(est, NA); lo <- c(lo, NA); hi <- c(hi, NA); next }
      if (!"Spill" %in% names(coef(f)) || any(is.na(coef(f)))) alias <- alias + 1
      se <- summary(f)$Coef["Z", "Std. Error"]; b <- coef(f)["Z"]
      est <- c(est, b); lo <- c(lo, b - 1.96 * se); hi <- c(hi, b + 1.96 * se)
    }
  }
  ok <- !is.na(est)
  data.frame(inc = ic$lab, nb = s$nb, spill = s$spill, design = s$d, struct = s$struct,
             MSE = mean((est[ok] - tau)^2), Bias = mean(est[ok]) - tau,
             Cov = mean(lo[ok] <= tau & hi[ok] >= tau), aliased = alias,
             corZX = mean(abs(apply(Zm, 2, cor, y = X[, 1]))))
}
t0 <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), run_one, mc.cores = 10))
saveRDS(res, "/private/tmp/claude-501/-Users-ajwalther-GithubProjects-SpatialCRT-projects-IncidenceDesign/9c7b3414-64ed-4762-917d-d48f8c5824e9/scratchpad/pilot_m1.rds")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min for", nrow(grid) * nD * nO, "fits\n")
nm <- c("1" = "Checkerboard", "2" = "HighIncFocus", "4" = "IsolBuffer", "5" = "2x2Block", "6" = "BalQuartiles", "8" = "IGSatQuad")
res$design <- nm[as.character(res$design)]
agg <- aggregate(cbind(MSE, Cov) ~ inc + struct + design, res, mean)
w <- reshape(agg, idvar = c("inc", "design"), timevar = "struct", direction = "wide")
w <- w[order(w$inc, w$MSE.M1A), ]
print(w, digits = 3, row.names = FALSE)
cat("\nAliased fits (spill term dropped) by design x nb, M1A:\n")
print(aggregate(aliased ~ design + nb, subset(res, struct == "M1A"), sum))
