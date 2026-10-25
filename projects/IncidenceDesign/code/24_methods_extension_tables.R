# ============================================================
# Script: 24_methods_extension_tables.R
# Purpose: LaTeX table bodies for the NC methods extensions, written straight from
#          the stored 2026-10-21 summaries so no number is transcribed by hand.
# Author: Andrew Walther
# Created: 2026-10-25
# Dependencies: base R only
# ============================================================
#
# Run from the project root:  Rscript code/24_methods_extension_tables.R
#
# Output: results/methods_ext_figures/tables/<name>.tex, one tabular body per table
# (rows between \midrule and \bottomrule), plus tables_log.txt with the in-text
# numbers. The manuscripts paste these bodies; verify_manuscripts.py does not check
# them, so rerun this script and diff after any upstream change.
#
# Conventions (author, 2026-10-21): t + df-corrected tests, random volunteering,
# beta = 0, kappa = 0.8; four planning years averaged; "Plain Saturation" in the
# stored files is displayed as "Saturation Regions".

R_EXT   <- "application/results/methods_ext_2026-10"
OUT_DIR <- "results/methods_ext_figures/tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Display order and names (manuscript order of the nine strategies)
ORDER <- c("Checkerboard", "High Incidence Focus", "Saturation Regions", "Isolation Buffer",
           "2x2 Blocking", "Balanced Quartiles", "Balanced Halves",
           "Incidence-Guided Saturation", "Simple Random Sampling")
display <- function(x) ifelse(x == "Plain Saturation", "Saturation Regions", x)
tex_name <- function(x) sub("2x2", "$2\\\\times2$", x)

rd <- function(f) read.csv(file.path(R_EXT, f), stringsAsFactors = FALSE)
f2 <- function(v) ifelse(is.finite(v), sprintf("%.2f", v), "--")
f3 <- function(v) ifelse(is.finite(v), sprintf("%.3f", v), "--")

#' Write one tabular body: one row per design in ORDER
#'
#' @param name File stem
#' @param cells Named list of character vectors, each named by design
#' @return Invisibly the lines written
write_body <- function(name, cells) {
  lines <- vapply(ORDER, function(d) {
    vals <- vapply(cells, function(col) { stopifnot(d %in% names(col)); col[[d]] }, character(1))
    paste(tex_name(d), "&", paste(vals, collapse = " & "), "\\\\")
  }, character(1))
  writeLines(lines, file.path(OUT_DIR, paste0(name, ".tex")))
  invisible(lines)
}

#' Pick one value per design from a summary data frame
#'
#' @param df Data frame with a Design column (stored labels)
#' @param value Column to return
#' @param fmt Formatter
#' @param ... Equality filters (column = value)
#' @return Named character vector over ORDER
pick <- function(df, value, fmt, ...) {
  f <- list(...)
  for (k in names(f)) df <- df[df[[k]] == f[[k]], ]
  df$Design <- display(df$Design)
  stopifnot(setequal(df$Design, ORDER), !anyDuplicated(df$Design))
  setNames(fmt(df[[value]]), df$Design)
}

log <- character(0)

# 1. Colleges needed and smallest detectable effect ----
n80 <- rd("summaries_v2/a3_n80_summary.csv")
mde <- rd("summaries_v2/a3_mde_summary.csv")
n80_fmt <- function(v) ifelse(is.finite(v), sprintf("%.0f", v), "$>58$")
mde_fmt <- function(v) ifelse(is.finite(v), sprintf("%.2f", v), "$>2$")  # 2 = top of the effect grid
write_body("tab_nc_power", list(
  pick(n80, "N80_Tdf_Mean", n80_fmt, Kappa = 0.8, Rho = 0.5, Regime = "both", Tau = 1),
  pick(n80, "N80_Tdf_Mean", n80_fmt, Kappa = 0.8, Rho = 0.5, Regime = "control_only", Tau = 1.5),
  pick(mde, "MDE_Tdf", mde_fmt, Kappa = 0.8, Rho = 0.5, Regime = "both", N = 29),
  pick(mde, "MDE_Tdf", mde_fmt, Kappa = 0.8, Rho = 0.5, Regime = "control_only", N = 29),
  pick(mde, "MDE_Tdf", mde_fmt, Kappa = 0.8, Rho = 0.5, Regime = "both", N = 58),
  pick(mde, "MDE_Tdf", mde_fmt, Kappa = 0.8, Rho = 0.5, Regime = "control_only", N = 58)))

# 2. Interval calibration by number of colleges (random volunteering) ----
st <- rd("summaries_v2/a3_settings.csv")
sz <- rd("summaries_v2/v2_size.csv")
cal <- do.call(rbind, lapply(sort(unique(st$N[st$Mechanism == "random"])), function(n) {
  s <- sz[sz$Mechanism == "random" & sz$N == n, ]
  c1 <- st[st$Mechanism == "random" & st$N == n & st$Beta == 0 & st$Tau == 1, ]
  dr <- aggregate(Coverage_Tdf ~ Design_ID + Regime, c1, mean)
  data.frame(N = n, Wald_size = sprintf("%.1f--%.1f", 100 * min(s$Size), 100 * max(s$Size)),
             Tdf_size = sprintf("%.1f--%.1f", 100 * min(s$Size_Tdf), 100 * max(s$Size_Tdf)),
             Wald_cov = sprintf("%.1f", 100 * mean(c1$Coverage)),
             Tdf_cov = sprintf("%.1f (%.1f--%.1f)", 100 * mean(c1$Coverage_Tdf),
                               100 * min(dr$Coverage_Tdf), 100 * max(dr$Coverage_Tdf)))
}))
writeLines(sprintf("%d & %s & %s & %s & %s \\\\", cal$N, cal$Wald_size, cal$Tdf_size,
                   cal$Wald_cov, cal$Tdf_cov), file.path(OUT_DIR, "tab_nc_calibration.tex"))
log <- c(log, "Calibration by N:", capture.output(print(cal, row.names = FALSE)))
s1 <- rd("stored/a5_stage1_coverage_by_design_regime.csv")
s1 <- s1[s1$Beta == 0, ]
log <- c(log, sprintf("Stage 1 (production fits, beta = 0) by design x regime: Wald %.1f-%.1f, t %.1f-%.1f, df-var %.1f-%.1f, both %.1f-%.1f",
                      100 * min(s1$Cov_Wald), 100 * max(s1$Cov_Wald), 100 * min(s1$Cov_t), 100 * max(s1$Cov_t),
                      100 * min(s1$Cov_dfvar), 100 * max(s1$Cov_dfvar),
                      100 * min(s1$Cov_dfvar_both), 100 * max(s1$Cov_dfvar_both)))
s2 <- rd("sensitivity_summaries/a5_stage2_coverage_by_design_regime.csv")
for (m in unique(s2$Method)) log <- c(log, sprintf("Stage 2 %s: %.1f-%.1f", m,
  100 * min(s2$Coverage[s2$Method == m]), 100 * max(s2$Coverage[s2$Method == m])))

# 3. Guidance under partial participation: MSE ratio to SRS at tau = 1 ----
gd <- rd("summaries_v2/a3_guidance_ratio_to_srs.csv")
names(gd)[names(gd) == "Ratio.Mean"] <- "Ratio"
cols <- list()
for (reg in c("both", "control_only")) for (n in c(20, 29, 40, 58))
  cols[[paste(reg, n)]] <- pick(gd, "Ratio", f2, Mechanism = "random", Regime = reg, N = n)
write_body("tab_nc_guidance", cols)
for (reg in c("both", "control_only")) for (n in c(20, 29)) {
  x <- gd[gd$Mechanism == "incidence3x" & gd$Regime == reg & gd$N == n, ]
  log <- c(log, sprintf("incidence-weighted %s N=%d: %s", reg, n,
                        paste(display(x$Design), f2(x$Ratio), collapse = "; ")))
}

# 4. Translation to trained professionals (SRS) ----
tr <- rd("summaries_v2/a3_translation_m_icc.csv")
tr <- tr[tr$Design_ID == 9 & tr$Rho == 0.5 & tr$Kappa == 0.8, ]
rows <- character(0)
for (m in c(10, 25, 50)) for (icc in c(0.01, 0.05, 0.10)) {
  v <- vapply(list(c("both", 29), c("both", 58), c("control_only", 29), c("control_only", 58)),
              function(k) vapply(c(0.2, 0.3), function(d) {
                x <- tr[tr$Regime == k[1] & tr$N == as.numeric(k[2]) & tr$m == m &
                          abs(tr$ICC - icc) < 1e-9 & abs(tr$d - d) < 1e-9, ]
                stopifnot(nrow(x) == 1)
                sprintf("%.2f", x$Power_Tdf)
              }, character(1)), character(2))
  rows <- c(rows, sprintf("%d & %.2f & %.3f & %s \\\\", m, icc, sqrt(icc + (1 - icc) / m),
                          paste(v, collapse = " & ")))
}
writeLines(rows, file.path(OUT_DIR, "tab_nc_translation.tex"))

# 5. Estimation method: CAR fitting and CAR data ----
car <- rd("sensitivity_summaries/a6_car_ratio_to_srs.csv")
cols <- list()
for (reg in c("both", "control_only")) for (k in list(c("sar", "SAR"), c("sar", "CAR"), c("car", "SAR"), c("car", "CAR")))
  cols[[paste(reg, k, collapse = " ")]] <- pick(car, "Ratio_Mean", f2, DGP = k[1], Estimator = k[2], Regime = reg)
write_body("tab_nc_car", cols)

# 6. Bayesian fits at rho = 0.5 beside ML on the same data ----
bl <- rd("sensitivity_summaries_v2/a6_bayes_lag_ratio_to_srs.csv")
cols <- list()
for (reg in c("both", "control_only")) for (e in c("SAR", "BLAG", "BCAR_converged"))
  cols[[paste(reg, e)]] <- pick(bl, "Ratio_Mean", f2, Estimator = e, Rho = 0.5, Regime = reg)
write_body("tab_nc_bayes", cols)

# 7. Spillover form ----
fm <- rd("sensitivity_summaries/a7_ratio_to_srs.csv")
idt <- rd("sensitivity_summaries/a7_identification.csv")
form_col <- function(form, reg) {
  v <- pick(fm, "Ratio_Mean", f2, Form = form, Regime = reg)
  i <- idt[idt$Form == form & idt$Regime == reg, ]; i$Design <- ORDER[match(i$Design_ID, c(1:9))]
  stopifnot(nrow(i) == 9)
  v[i$Design[i$Identified < 1]] <- "n.i."
  v
}
write_body("tab_nc_forms", list(form_col("proportional", "both"), form_col("count", "both"),
                                form_col("binary", "both"), form_col("proportional", "control_only"),
                                form_col("count", "control_only"), form_col("binary", "control_only")))

# 8. Unequal trained professionals and NC tests (in-text numbers) ----
a4 <- rd("summaries/a4_summary.csv")
a4d <- a4[abs(a4$Tau - 0.3) < 1e-9, ]
for (n in c(20, 58)) for (reg in c("both", "control_only")) {
  x <- a4d[a4d$N == n & a4d$Regime == reg, ]
  log <- c(log, sprintf("A4 d=0.3 N=%d %s: mean power equal %.2f -> unequal %.2f", n, reg,
                        mean(x$Power_Equal), mean(x$Power_Unequal)))
}
log <- c(log, sprintf("A4 overall: mean power difference %.3f; MSE ratio range %.2f-%.2f",
                      mean(a4$Power_Diff[a4$Tau > 0]), min(a4$MSE_Ratio), max(a4$MSE_Ratio)))
fr <- rd("stored/a2_tests_friedman.csv")
ar <- rd("stored/a2_tests_avg_ranks.csv")
for (i in seq_len(nrow(fr))) {
  r <- ar[ar$Setting == fr$Setting[i], ]; r <- r[order(r$Avg_Rank), ]
  log <- c(log, sprintf("A2 %s: chi2 = %.1f, df = %d, p = %.1e, blocks = %d, CD = %.2f; ranks %s",
                        fr$Setting[i], fr$Statistic[i], fr$DF[i], fr$P[i], fr$N_Blocks[i], fr$CD[i],
                        paste(sprintf("%s %.2f", display(r$Design), r$Avg_Rank), collapse = ", ")))
}

writeLines(log, file.path(OUT_DIR, "tables_log.txt"))
cat(log, sep = "\n")
cat("Tables written to", OUT_DIR, "\n")
