# ============================================================
# Script: 27_interval_comparison.R
# Purpose: Coverage of the Wald and t + df-corrected 95% intervals side by side,
#          for the grid benchmark and the NC application, plus the application's
#          yearly coverage under the t + df interval (the reported interval).
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: base R
# ============================================================
#
# Run from the project root:  Rscript code/27_interval_comparison.R
#
# Inputs (read only; both intervals were computed from the same stored fits):
#   results/grid_tdf_2026-10/grid_tdf_scenarios.rds        grid, Wald + t + df (26_)
#   application/results/methods_ext_2026-10/stored/a5_stage1_coverage.csv
#       application primary rows: Cov_Wald and Cov_dfvar_both (= t + df, df = 54)
#   application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv
#       production yearly means (Wald), used only to check the Wald reproduction
#
# Outputs:
#   application/results/methods_ext_2026-10/stored/yearly_primary_coverage_tdf.csv
#       Year x Regime x Design: Coverage_Wald, Coverage_Tdf (equal weight to the six
#       rho/gamma pairs, as the production yearly means); feeds tab:nc-annual
#   results/grid_tdf_2026-10/interval_comparison.csv
#       six designs + SRS: grid (tau = 1, oracle, pooled over regimes and the five
#       incidence configurations) by neighbor type, and the application (queen,
#       beta = 0, mean over years and the six rho/gamma pairs) by regime, each Wald
#       and t + df; feeds the chapter/SI interval-comparison table
#
# Both intervals are kept so either can be reported later; the manuscripts report
# t + df (author, 2026-10-26) with Wald in the comparison table.

GRID  <- "results/grid_tdf_2026-10/grid_tdf_scenarios.rds"
A5    <- "application/results/methods_ext_2026-10/stored/a5_stage1_coverage.csv"
YEARS <- "application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv"
OUT_YEAR <- "application/results/methods_ext_2026-10/stored/yearly_primary_coverage_tdf.csv"
OUT_CMP  <- "results/grid_tdf_2026-10/interval_comparison.csv"

# Table order (chapter/SI design table order, then SRS) and design numbers
IDS   <- c(8, 6, 4, 2, 5, 1, 9)
NAMES <- c("Incidence-Guided Saturation", "Balanced Quartiles", "Isolation Buffer",
           "High Incidence Focus", "2x2 Blocking", "Checkerboard", "Simple Random Sampling")

# Application: yearly coverage under both intervals ----
a5 <- read.csv(A5)
a5 <- a5[a5$Beta == 0, ]
stopifnot(nrow(a5) == 4 * 6 * 2 * 9)                 # years x rho/gamma x regimes x designs
year <- aggregate(cbind(Coverage_Wald = Cov_Wald, Coverage_Tdf = Cov_dfvar_both) ~
                    Year + Regime + Design_ID, data = a5, FUN = mean)
stopifnot(all(aggregate(Cov_Wald ~ Year + Regime + Design_ID, a5, length)$Cov_Wald == 6))

# Fail loud unless the Wald column reproduces the production yearly means exactly
prod <- read.csv(YEARS)
key  <- function(d) paste(d$Year, d$Regime, d$Design_ID)
m <- match(key(year), key(prod))
stopifnot(!anyNA(m), nrow(year) == nrow(prod))
wald_diff <- max(abs(year$Coverage_Wald - prod$Coverage[m]))
if (wald_diff > 1e-12) stop("Wald coverage does not reproduce the production yearly means: ", wald_diff)
year$Design <- prod$Design[m]
year <- year[order(year$Year, year$Regime, year$Design_ID),
             c("Year", "Regime", "Design_ID", "Design", "Coverage_Wald", "Coverage_Tdf")]
write.csv(year, OUT_YEAR, row.names = FALSE)

# Both-interval comparison table ----
g <- readRDS(GRID)
g <- g[g$Estimator == "oracle" & g$True_Tau == 1, ]
grid_cov <- function(nb, col) {
  d <- g[g$Neighbor_Type == nb, ]
  vapply(IDS, function(i) mean(d[[col]][d$Design == paste("Design", i)]), numeric(1))
}
app_cov <- function(regime, col) {
  d <- a5[a5$Regime == regime, ]
  vapply(IDS, function(i) mean(d[[col]][d$Design_ID == i]), numeric(1))
}
cmp <- data.frame(
  Design = NAMES, Design_ID = IDS,
  Grid_Queen_Wald = grid_cov("queen", "Coverage"), Grid_Queen_Tdf = grid_cov("queen", "Coverage_Tdf"),
  Grid_Rook_Wald  = grid_cov("rook", "Coverage"),  Grid_Rook_Tdf  = grid_cov("rook", "Coverage_Tdf"),
  App_Both_Wald    = app_cov("both", "Cov_Wald"),
  App_Both_Tdf     = app_cov("both", "Cov_dfvar_both"),
  App_Control_Wald = app_cov("control_only", "Cov_Wald"),
  App_Control_Tdf  = app_cov("control_only", "Cov_dfvar_both"))
stopifnot(all(table(g$Design, g$Neighbor_Type) == 160))   # 5 configs x 4 rho x 4 gamma x 2 regimes
write.csv(cmp, OUT_CMP, row.names = FALSE)

cat(sprintf("Application yearly coverage: Wald %.3f-%.3f, t + df %.3f-%.3f (%d rows; Wald reproduced, max diff %.1e)\n",
            min(year$Coverage_Wald), max(year$Coverage_Wald), min(year$Coverage_Tdf),
            max(year$Coverage_Tdf), nrow(year), wald_diff))
print(format(cmp, digits = 3), row.names = FALSE)
