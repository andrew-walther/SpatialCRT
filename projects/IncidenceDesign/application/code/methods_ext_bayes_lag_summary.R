# ============================================================
# Script: methods_ext_bayes_lag_summary.R
# Purpose: Summarize the Bayesian SAR-lag fits (a6bayeslag, spatialreg::spBreg_lag)
#   beside SAR ML and the Bayesian CAR (a6bayes at rho = 0, the long-chain
#   a6bayeslong at rho = 0.5) on the same data, to separate "Bayesian" from "CAR".
# Author: Andrew Walther
# Created: 2026-10-21
# Dependencies: methods_ext_sensitivity_summary.R (mz_setting, mz_ratio, mz_spearman)
# ============================================================
# Run from the project root:  Rscript application/code/methods_ext_bayes_lag_summary.R [suffix]
# ("_smoke" reads the *_smoke folders; the long-chain comparison is skipped if absent.)
# Writes application/results/methods_ext_2026-10/sensitivity_summaries_v2/:
#   a6_bayes_lag_settings.csv, a6_bayes_lag_ratio_to_srs.csv, a6_bayes_lag_spearman.csv,
#   a6_bayes_lag_log.txt
# Fails loud (exit 1) if the SAR ML fits differ from a6bayes's on the shared data.
# Converged: ESS(tau) >= 400 and |Geweke z| < 1.96, as for the Bayesian CAR.

.ml_dir <- local({
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else normalizePath("application/code")
})
source(file.path(.ml_dir, "methods_ext_sensitivity_summary.R"))

args <- commandArgs(TRUE); suffix <- if (length(args)) args[1] else ""
out <- file.path(.me_root, paste0("sensitivity_summaries_v2", suffix)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
log <- character(0); say <- function(...) { l <- sprintf(...); cat(l, "\n"); log <<- c(log, l) }
bail <- function(...) { say(...); writeLines(log, file.path(out, "a6_bayes_lag_log.txt")); quit(status = 1) }

lag <- readRDS(file.path(.me_root, paste0("a6bayeslag", suffix), "fits.rds"))
car <- readRDS(file.path(.me_root, paste0("a6bayes", suffix), "fits.rds"))
cfg <- readRDS(file.path(.me_root, paste0("a6bayeslag", suffix), "manifest.rds"))$config
key <- c("Year", "Design_ID", "Regime", "Gamma", "Rho", "Allocation", "Replicate")
m <- merge(lag[c(key, "SAR_tau")], car[c(key, "SAR_tau")], by = key, suffixes = c("_lag", "_car"))
if (nrow(m) != nrow(lag) || nrow(m) != nrow(car)) bail("FAIL: %d lag fits, %d CAR fits, %d matched", nrow(lag), nrow(car), nrow(m))
if (!identical(m$SAR_tau_lag, m$SAR_tau_car)) bail("FAIL: SAR ML fits differ between a6bayeslag and a6bayes (not the same data)")
say("Same data confirmed: %d fits, SAR ML estimates identical to a6bayes", nrow(m))

conv <- function(ess, gz) is.finite(ess) & ess >= cfg$ess_min & abs(gz) < 1.96
lag$Conv <- conv(lag$BLAG_ess, lag$BLAG_geweke)
say("Bayesian SAR lag (spBreg_lag, %d draws, first %d discarded): errors %d; converged %.1f%%; ESS median %.0f (min %.0f); |Geweke z| >= 1.96 in %.1f%% (about 5%% expected by chance)",
    cfg$ndraw, cfg$nomit, sum(lag$BLAG_error), 100 * mean(lag$Conv), stats::median(lag$BLAG_ess, na.rm = TRUE),
    min(lag$BLAG_ess, na.rm = TRUE), 100 * mean(abs(lag$BLAG_geweke) >= 1.96, na.rm = TRUE))
for (r in sort(unique(lag$Rho))) say("  rho = %.1f: converged %.1f%%", r, 100 * mean(lag$Conv[lag$Rho == r]))
z <- lag[is.finite(lag$BLAG_tau) & is.finite(lag$SAR_tau), ]
say("Posterior mean vs ML: correlation %.4f, mean |diff| %.4f, mean diff %+.4f over %d fits",
    stats::cor(z$BLAG_tau, z$SAR_tau), mean(abs(z$BLAG_tau - z$SAR_tau)), mean(z$BLAG_tau - z$SAR_tau), nrow(z))

# Bayesian CAR beside it: short chains at rho = 0 (converged 92-94%), long chains at rho = 0.5
bc <- car; bc$Source <- "a6bayes"
lp <- file.path(.me_root, paste0("a6bayeslong", suffix), "fits.rds")
if (file.exists(lp)) { bl <- readRDS(lp); bl$Source <- "a6bayeslong"; bc <- rbind(bc[bc$Rho != 0.5, ], bl[names(bc)]) }
bc$Conv <- conv(bc$BCAR_ess, bc$BCAR_geweke)
say("Bayesian CAR beside it: %s", paste(sprintf("rho = %.1f from %s", sort(unique(bc$Rho)), tapply(bc$Source, bc$Rho, `[`, 1)), collapse = "; "))

by <- c("Year", "Rho", "Regime", "Design_ID")
lc <- lag; lc[!lc$Conv, c("BLAG_tau", "BLAG_lo", "BLAG_hi")] <- NA
cc <- bc; cc[!cc$Conv, c("BCAR_tau", "BCAR_lo", "BCAR_hi")] <- NA
s <- rbind(mz_setting(lag, "SAR", by), mz_setting(lag, "BLAG", by, "BLAG_lo", "BLAG_hi"),
           transform(mz_setting(lc, "BLAG", by, "BLAG_lo", "BLAG_hi"), Estimator = "BLAG_converged"),
           transform(mz_setting(cc, "BCAR", by, "BCAR_lo", "BCAR_hi"), Estimator = "BCAR_converged"))
write.csv(s, file.path(out, "a6_bayes_lag_settings.csv"), row.names = FALSE)
cv <- aggregate(cbind(Coverage, Bias, MSE) ~ Estimator + Rho + Regime, s, mean)
for (i in seq_len(nrow(cv))) say("rho = %.1f, %-12s %-15s coverage %.3f, bias %+.3f, MSE %.3f", cv$Rho[i], cv$Regime[i], cv$Estimator[i], cv$Coverage[i], cv$Bias[i], cv$MSE[i])
rt <- mz_ratio(s, "Year", c("Estimator", "Rho")); write.csv(rt, file.path(out, "a6_bayes_lag_ratio_to_srs.csv"), row.names = FALSE)
for (e in unique(rt$Estimator)) for (r in sort(unique(rt$Rho))) for (rg in unique(rt$Regime)) {
  d <- rt[rt$Estimator == e & rt$Rho == r & rt$Regime == rg, ]; d <- d[order(d$Ratio_Mean), ]
  say("MSE ratio to SRS, %s, rho = %.1f, %s: %s", e, r, rg, paste(sprintf("%s %.2f", d$Design, d$Ratio_Mean), collapse = ", ")) }
sp <- do.call(rbind, lapply(c("BLAG", "BLAG_converged", "BCAR_converged"), function(b)
  cbind(Comparison = paste("SAR vs", b), mz_spearman(s, c("Year", "Rho", "Regime"), "SAR", b))))
write.csv(sp, file.path(out, "a6_bayes_lag_spearman.csv"), row.names = FALSE)
for (cmp in unique(sp$Comparison)) for (r in sort(unique(sp$Rho))) for (rg in unique(sp$Regime)) {
  d <- sp[sp$Comparison == cmp & sp$Rho == r & sp$Regime == rg, ]
  say("Spearman %s, rho = %.1f, %s: mean %.2f (range %.2f-%.2f); best design agrees %d/%d", cmp, r, rg,
      mean(d$Spearman, na.rm = TRUE), min(d$Spearman, na.rm = TRUE), max(d$Spearman, na.rm = TRUE), sum(d$Best_A == d$Best_B), nrow(d)) }
writeLines(c(sprintf("Run %s", format(Sys.time())), log), file.path(out, "a6_bayes_lag_log.txt"))
