# ============================================================
# Script: methods_ext_bayes_long_summary.R
# Purpose: Compare the 2026-10-21 long-chain Bayesian CAR rerun (a6bayeslong,
#   rho = 0.5, 20,000 burn-in + 100,000 iterations thinned by 20) with the
#   2026-10-20 run (a6bayes, 5,000 + 25,000 thinned by 5) on the same data:
#   convergence, then coverage, bias and design ordering with all fits and with
#   converged fits only.
# Author: Andrew Walther
# Created: 2026-10-21
# Dependencies: methods_ext_sensitivity_summary.R (mz_setting, mz_ratio, mz_spearman)
# ============================================================
# Run from the project root:  Rscript application/code/methods_ext_bayes_long_summary.R [suffix]
# ("_smoke" reads a6bayes_smoke / a6bayeslong_smoke.)
# Writes application/results/methods_ext_2026-10/sensitivity_summaries_v2/:
#   a6_bayes_long_settings.csv, a6_bayes_long_ratio_to_srs.csv,
#   a6_bayes_long_spearman.csv, a6_bayes_long_paired.csv, a6_bayes_long_log.txt
# Fails loud (exit 1) if the SAR fits on the shared data are not identical, which
# would mean the two runs did not see the same data.
#
# Converged: ESS(tau) >= 400 and |Geweke z| < 1.96 (as the 2026-10-20 run).

.mb_dir <- local({
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else normalizePath("application/code")
})
source(file.path(.mb_dir, "methods_ext_sensitivity_summary.R"))

args <- commandArgs(TRUE); suffix <- if (length(args)) args[1] else ""
out <- file.path(.me_root, paste0("sensitivity_summaries_v2", suffix)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
log <- character(0); say <- function(...) { l <- sprintf(...); cat(l, "\n"); log <<- c(log, l) }

old <- readRDS(file.path(.me_root, paste0("a6bayes", suffix), "fits.rds"))
new <- readRDS(file.path(.me_root, paste0("a6bayeslong", suffix), "fits.rds"))
cfg_new <- readRDS(file.path(.me_root, paste0("a6bayeslong", suffix), "manifest.rds"))$config
old <- old[old$Rho %in% unique(new$Rho), ]
key <- c("Year", "Design_ID", "Regime", "Gamma", "Rho", "Allocation", "Replicate")
m <- merge(old, new, by = key, suffixes = c("_short", "_long"))
if (nrow(m) != nrow(new) || nrow(m) != nrow(old)) {
  say("FAIL: %d short and %d long fits, %d matched", nrow(old), nrow(new), nrow(m)); writeLines(log, file.path(out, "a6_bayes_long_log.txt")); quit(status = 1)
}
same <- isTRUE(all.equal(m$SAR_tau_short, m$SAR_tau_long, tolerance = 0)) || identical(m$SAR_tau_short, m$SAR_tau_long)
if (!same) { say("FAIL: SAR fits differ between runs (max |diff| %g); the runs did not share data",
                 max(abs(m$SAR_tau_short - m$SAR_tau_long), na.rm = TRUE)); writeLines(log, file.path(out, "a6_bayes_long_log.txt")); quit(status = 1) }
say("Same data confirmed: %d fits, SAR estimates identical in both runs", nrow(m))

conv <- function(ess, gz) is.finite(ess) & ess >= cfg_new$ess_min & abs(gz) < 1.96
m$Conv_short <- conv(m$BCAR_ess_short, m$BCAR_geweke_short); m$Conv_long <- conv(m$BCAR_ess_long, m$BCAR_geweke_long)
say("Long chains (burn-in %d, %d iterations, thin %d): errors %d; converged %.1f%% (short chains %.1f%%); ESS median %.0f (short %.0f), min %.0f",
    cfg_new$burnin, cfg_new$n_sample, cfg_new$thin, sum(m$BCAR_error_long), 100 * mean(m$Conv_long), 100 * mean(m$Conv_short),
    stats::median(m$BCAR_ess_long, na.rm = TRUE), stats::median(m$BCAR_ess_short, na.rm = TRUE), min(m$BCAR_ess_long, na.rm = TRUE))
for (rg in unique(m$Regime)) { d <- m[m$Regime == rg, ]
  say("  %s: converged long %.1f%% vs short %.1f%%; ESS < 400 long %d; |Geweke| >= 1.96 long %d", rg,
      100 * mean(d$Conv_long), 100 * mean(d$Conv_short), sum(d$BCAR_ess_long < cfg_new$ess_min, na.rm = TRUE),
      sum(abs(d$BCAR_geweke_long) >= 1.96, na.rm = TRUE)) }
say("Posterior mean, long vs short chain: correlation %.3f, mean |diff| %.3f", stats::cor(m$BCAR_tau_long, m$BCAR_tau_short, use = "complete.obs"),
    mean(abs(m$BCAR_tau_long - m$BCAR_tau_short), na.rm = TRUE))
write.csv(m[c(key, "SAR_tau_short", "BCAR_tau_short", "BCAR_tau_long", "BCAR_ess_short", "BCAR_ess_long",
              "BCAR_geweke_short", "BCAR_geweke_long", "Conv_short", "Conv_long")], file.path(out, "a6_bayes_long_paired.csv"), row.names = FALSE)

# Setting metrics: SAR ML, long-chain Bayes (all fits), long-chain Bayes (converged only)
f <- new; f$BCAR_converged <- conv(f$BCAR_ess, f$BCAR_geweke)
by <- c("Year", "Rho", "Regime", "Design_ID")
fc <- f; fc[!fc$BCAR_converged, c("BCAR_tau", "BCAR_lo", "BCAR_hi")] <- NA
s <- rbind(mz_setting(f, "SAR", by), mz_setting(f, "BCAR", by, "BCAR_lo", "BCAR_hi"),
           transform(mz_setting(fc, "BCAR", by, "BCAR_lo", "BCAR_hi"), Estimator = "BCAR_converged"))
write.csv(s, file.path(out, "a6_bayes_long_settings.csv"), row.names = FALSE)
cv <- aggregate(cbind(Coverage, Bias, MSE) ~ Estimator + Regime, s, mean)
for (i in seq_len(nrow(cv))) say("rho = 0.5, %s, %s: coverage %.3f, bias %+.3f, MSE %.3f", cv$Estimator[i], cv$Regime[i], cv$Coverage[i], cv$Bias[i], cv$MSE[i])
rt <- mz_ratio(s, c("Year", "Rho"), "Estimator"); write.csv(rt, file.path(out, "a6_bayes_long_ratio_to_srs.csv"), row.names = FALSE)
for (e in unique(rt$Estimator)) for (rg in unique(rt$Regime)) { d <- rt[rt$Estimator == e & rt$Regime == rg, ]; d <- d[order(d$Ratio_Mean), ]
  say("MSE ratio to SRS, %s, %s: %s", e, rg, paste(sprintf("%s %.2f", d$Design, d$Ratio_Mean), collapse = ", ")) }
sp <- rbind(cbind(Comparison = "SAR vs BCAR", mz_spearman(s, c("Year", "Rho", "Regime"), "SAR", "BCAR")),
            cbind(Comparison = "SAR vs BCAR_converged", mz_spearman(s, c("Year", "Rho", "Regime"), "SAR", "BCAR_converged")))
write.csv(sp, file.path(out, "a6_bayes_long_spearman.csv"), row.names = FALSE)
for (cmp in unique(sp$Comparison)) for (rg in unique(sp$Regime)) { d <- sp[sp$Comparison == cmp & sp$Regime == rg, ]
  say("Spearman %s, %s: mean %.2f (range %.2f-%.2f); best design agrees %d/%d", cmp, rg, mean(d$Spearman, na.rm = TRUE),
      min(d$Spearman, na.rm = TRUE), max(d$Spearman, na.rm = TRUE), sum(d$Best_A == d$Best_B), nrow(d)) }
writeLines(c(sprintf("Run %s", format(Sys.time())), log), file.path(out, "a6_bayes_long_log.txt"))
