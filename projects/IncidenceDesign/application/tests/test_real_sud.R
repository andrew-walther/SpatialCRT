# ============================================================
# Script: test_real_sud.R
# Purpose: Verify model meaning, allocation budgets, caching and uncertainty.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud modules, spatialreg
# ============================================================
test_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
real_sud_define_only <- TRUE
source(file.path(test_dir, "..", "code", "run_real_sud.R"))
s <- readRDS(file.path(.rs_app, "results", "real_sud_rev_20261002", "setup.rds"))
c <- rs_config("smoke"); u <- rs_units(c)

#' Assert a named scientific behavior, failing immediately if wrong
#' @param name Human-readable check.
#' @param ok Logical condition, required TRUE.
#' @return Invisibly TRUE; otherwise stops.
#' @examples
#' check("arithmetic", 1 + 1 == 2)
check <- function(name, ok) {
  cat(sprintf("%-70s %s\n", name, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop(name, call. = FALSE)
  invisible(TRUE)
}

# Fixed observations and budgets ----
check("Production has 1,248 annual reporting blocks", nrow(rs_units(rs_config("production"))) == 1248)
check("Each observed year has 58 unique rates", all(vapply(s$annual[1:4], function(a) length(unique(a$rate_per_100k)) == 58, logical(1))))
for (id in c(2, 6, 7, 9)) {
  v <- u[which(u$Design_ID == id)[1], ]; Z <- rs_draws(v, s, 100)
  check(paste("Design", id, "treats exactly 29 in every draw"), all(colSums(Z) == 29))
  if (id == 7) {
    half <- dplyr::ntile(s$annual[[v$Year]]$rank, 2)
    check("Balanced Halves alternates the extra treatment between odd halves",
      setequal(colSums(Z[half == 1, , drop = FALSE]), c(14, 15)))
  }
}
v <- u[which(u$Design_ID == 4)[1], ]; Z <- rs_draws(v, s, 20)
check("Buffer never treats adjacent clusters", all(vapply(seq_len(ncol(Z)), function(j) sum(Z[, j] * (s$weights$queen$adjacency %*% Z[, j])) == 0, logical(1))))
v <- u[which(u$Design_ID == 8)[1], ]
rv <- rs_region_values(s$annual[[v$Year]], s$regions, "mean_rank")
check("Regional primary rule equals mean of per-capita cluster ranks",
  isTRUE(all.equal(rv, tapply(rank(s$annual[[v$Year]]$deaths / s$annual[[v$Year]]$person_years), s$regions, mean))))

# Model meaning and estimator equivalence ----
v <- u[which(u$Design_ID == 9 & u$Rho == 0.5 & u$Regime == "control_only")[1], ]
z <- rs_draws(v, s, 1)[, 1]; m <- rs_model(v, z, s, c)
check("Education SAR mean has no incidence baseline",
  max(abs((diag(58) - v$Rho * s$weights$queen$W) %*% m$mean - (z + v$Gamma * as.vector(s$weights$queen$W %*% z) * (1 - z)))) < 1e-12 && !"X" %in% colnames(m$xm))
vv <- v; vv$Model <- "baseline_sensitivity"; mm <- rs_model(vv, z, s, c)
check("Matched baseline sensitivity adds X to both DGP and fit",
  max(abs((diag(58) - v$Rho * s$weights$queen$W) %*% (mm$mean - m$mean) - s$annual[[v$Year]]$X)) < 1e-12 && "X" %in% colnames(mm$xm))
vv <- v; vv$Regime <- "both"; mb <- rs_model(vv, z, s, c)
check("Both-arms spillover can affect treated clusters", any(mb$xm[z == 1, "Spill"] > 0) && all(m$xm[z == 1, "Spill"] == 0))
lw <- spdep::mat2listw(s$weights$queen$W, style = "W"); engine <- sar_lag_setup(lw)
for (model in c("education", "baseline_sensitivity")) {
  vv <- v; vv$Model <- model; mm <- rs_model(vv, z, s, c)
  rs_seed("test-equivalence", model); y <- mm$mean + as.vector(mm$A %*% rnorm(58))
  lean <- fit_one_lag_model(y, mm$xm, engine)
  ref <- fit_one_lag_model(y, mm$xm, engine, lw, "lagsarlm")
  check(paste("Lean engine matches lagsarlm:", model),
    abs(lean$tau - ref$tau) < 1e-6 && abs(lean$se - ref$se) < 1e-6)
}
# A constructed exact alias must never be counted as identified tau.
ss <- s; ss$weights$queen$W <- (matrix(1, 58, 58) - diag(58)) / 57
vv <- v; vv$Regime <- "both"
check("Exact treatment/spillover alias is flagged", !rs_model(vv, z, ss, c)$Identified)

# Reuse, prefixes and independent-fit accounting ----
b <- rs_block(v, s, c, 3, 8)
extended <- rs_block(v, s, c, 5, 13, b); fresh <- rs_block(v, s, c, 5, 13)
check("Cached J/R extension equals a fresh run exactly", identical(extended, fresh))
parallel_blocks <- parallel::mclapply(1:2, function(i) rs_block(v, s, c, 5, 13),
  mc.cores = 2, mc.set.seed = FALSE)
check("Worker count/order cannot change keyed allocations or outcomes",
  all(vapply(parallel_blocks, function(b) identical(b, fresh), logical(1))))
check("Extension preserves allocation and outcome prefixes",
  identical(b$Z, extended$Z[, 1:3, drop = FALSE]) && all(vapply(names(b$fits), function(h) identical(b$fits[[h]], extended$fits[[h]][1:8, ]), logical(1))))
vb <- u[which(u$Design_ID == 4 & u$Rho == 0.5 & u$Regime == "both")[1], ]
zb <- rs_draws(vb, s, 1)[, 1]; vc <- vb; vc$Regime <- "control_only"
check("Buffer has exactly equal SAR mean/matrix under both spillover regimes",
  identical(rs_model(vb, zb, s, c)$mean, rs_model(vc, zb, s, c)$mean) &&
  identical(rs_model(vb, zb, s, c)$xm, rs_model(vc, zb, s, c)$xm))
vv <- v; vv$Year <- "2021"
check("Primary SRS shares performance source across years", identical(rs_block_key(v, s, c), rs_block_key(vv, s, c)))
rr <- rs_reporting_rows(vv, b, s)
check("Reused performance still recomputes annual population shares",
  abs(rr$allocations$Population_Share[1] - sum(s$annual$`2021`$person_years * b$Z[, 1]) / sum(s$annual$`2021`$person_years)) < 1e-12)
one <- u[which(u$Design_ID == 1)[1], ]; bo <- rs_block(one, s, c, 10, 20)
check("Singleton 10 duplicate draws count as only 20 independent fits", bo$summary$N_Independent_Fits == 20 && bo$summary$Proven_Singleton)
check("Singleton precision is conditional R precision, not divided by sqrt(J)",
  abs(bo$summary$SE_Mean_MSE_Joint - bo$allocations$SE_MSE) < 1e-12)

# Hand-calculated duplicate covariance and invalid fits ----
f <- data.frame(replicate = 1:4, error = c(-1, 0, 1, 2), se = rep(1, 4), alias = FALSE, messages = "")
a <- rs_allocation_metrics(f); a$Frequency <- 3L
aa <- a; aa$MSE <- 4; aa$Frequency <- 1L
sdist <- rs_distribution(rbind(a, aa), FALSE)
raw <- (3 * (1.5 - 2.125)^2 + (4 - 2.125)^2) / 3
check("Duplicate-aware joint MSE variance matches hand formula",
  abs(sdist$SE_Mean_MSE_Joint^2 - (raw / 4 + 3 * 2 * a$SE_MSE^2 / 12)) < 1e-12)
f$messages[1] <- "rho on interval bound - results should not be used"
bad <- rs_allocation_metrics(f)
check("Finite boundary fit invalidates complete comparison", !bad$Complete && bad$N_Bound == 1 && bad$N_Valid_Est == 3)
tmp <- tempfile(); dir.create(tmp)
rs_manifest(tmp, s, c); cc <- c; cc$tau <- 2
check("Manifest rejects changed statistical inputs", inherits(tryCatch(rs_manifest(tmp, s, cc), error = identity), "error"))
cat("\nAll real-SUD behavioral checks passed.\n")
