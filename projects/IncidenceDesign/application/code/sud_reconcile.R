# ============================================================
# Script: sud_reconcile.R
# Purpose: QC report comparing the SUD numerator (corrected counts in
#          final_county_sudden.csv) with the older counts behind Habib (2026).
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: dplyr (via sud_incidence.R)
# ============================================================
#
# Step 3 of the real-data aggregation (after sud_incidence.R). The numerator is the
# corrected count (deaths, 23,523). The Habib (2026) counts (deaths_habib2026, 21,147)
# are kept so the paper's published numbers stay reproducible: the paper checks run on
# that secondary column. Background and the 2026-09-25 decision: application/README.md.

#' Published statewide figures from Habib (2026) that the secondary counts must reproduce
SUD_PAPER_TARGETS <- list(
  deaths = 21147,
  person_years = 25594321,
  overall_rate = 82.6,
  yearly_rate = c(`2018` = 75.1, `2019` = 77.9, `2020` = 86.7, `2021` = 90.6)
)

#' Write the numerator reconciliation report
#'
#' @description Four sections: (1) statewide deaths and rates for the numerator
#'   (corrected counts) and the Habib (2026) counts, next to the paper; (2)
#'   pass/fail checks of the Habib (2026) counts against the paper; (3) how the two
#'   definitions differ county-year by county-year; (4) whether they rank the 58
#'   clusters the same way (Spearman, same top half, same quartile: the inputs that
#'   High Incidence Focus and Balanced Quartiles use). A mismatch with the paper
#'   raises a warning rather than an error, so a revised source file can still be
#'   aggregated and inspected.
#'
#' @param county_df County-year data.frame from `load_sud_county_data()`.
#' @param mapping County-to-college table from `get_cc_mapping_data()`.
#' @param out_file Path of the text report.
#' @return Invisibly, the named logical vector of paper checks (on deaths_habib2026).
#' @family sud_aggregation
reconcile_sud_counts <- function(county_df, mapping, out_file) {
  pooled <- paste0(min(SUD_YEARS), "-", max(SUD_YEARS))
  state_df <- transform(county_df, state = "NC")
  st <- compute_incidence(state_df, "state", "deaths")
  st_hab <- compute_incidence(state_df, "state", "deaths_habib2026")
  cl <- aggregate_sud_to_clusters(county_df, mapping, "deaths")
  cl_hab <- aggregate_sud_to_clusters(county_df, mapping, "deaths_habib2026")

  # (2) Checks against the paper, using the Habib (2026) counts ----
  yr <- st_hab$period != pooled
  checks <- c(
    deaths       = st_hab$deaths[!yr] == SUD_PAPER_TARGETS$deaths,
    person_years = st_hab$person_years[!yr] == SUD_PAPER_TARGETS$person_years,
    overall_rate = round(st_hab$rate_per_100k[!yr], 1) == SUD_PAPER_TARGETS$overall_rate,
    yearly_rates = all(round(st_hab$rate_per_100k[yr], 1) == SUD_PAPER_TARGETS$yearly_rate[st_hab$period[yr]])
  )
  if (!all(checks)) {
    warning("deaths_habib2026 no longer matches Habib (2026): ", paste(names(checks)[!checks], collapse = ", "),
            call. = FALSE)
  }

  # (4) Cluster rank agreement between the two definitions, per period ----
  agree <- sapply(unique(cl$period), function(p) {
    a <- cl$rate_per_100k[cl$period == p]
    b <- cl_hab$rate_per_100k[cl_hab$period == p]
    n <- length(a)
    c(spearman = stats::cor(a, b, method = "spearman"),
      top_half = sum((rank(-a) <= n / 2) == (rank(-b) <= n / 2)),
      quartile = sum(ceiling(4 * rank(a) / n) == ceiling(4 * rank(b) / n)),
      n = n)
  })

  paper_rate <- c(SUD_PAPER_TARGETS$yearly_rate, stats::setNames(SUD_PAPER_TARGETS$overall_rate, pooled))
  diff <- county_df$deaths - county_df$deaths_habib2026
  lines <- c(
    "NC SUD numerator reconciliation (written by run_sud_aggregation.R)",
    "Numerator: final_county_sudden.csv num_obs (Habib's corrected case filtering; heart-failure",
    "deaths no longer excluded; author decision 2026-09-25, Habib personal communication).",
    "Habib (2026): sudden_county_year.csv num_obs, the older counts behind the published paper.",
    "",
    "1. Statewide deaths and rates per 100,000, ages 18-64",
    sprintf("  %-10s %8s %8s %8s %8s %8s", "period", "deaths", "h2026", "rate", "h2026", "paper"),
    sprintf("  %-10s %8d %8d %8.1f %8.1f %8.1f", st$period, st$deaths, st_hab$deaths,
            st$rate_per_100k, st_hab$rate_per_100k, paper_rate[st$period]),
    "",
    "2. Habib (2026) counts checked against the paper",
    sprintf("  %-13s %s", names(checks), ifelse(checks, "MATCH", "MISMATCH")),
    "",
    "3. County-years, numerator minus Habib (2026)",
    sprintf("  equal %d | numerator higher %d | numerator lower %d",
            sum(diff == 0), sum(diff > 0), sum(diff < 0)),
    "",
    "4. Cluster rank agreement, numerator vs Habib (2026) (58 clusters)",
    sprintf("  %-10s %9s %9s %9s", "period", "Spearman", "top half", "quartile"),
    sprintf("  %-10s %9.3f %6d/%d %6d/%d", colnames(agree), agree["spearman", ],
            agree["top_half", ], agree["n", ], agree["quartile", ], agree["n", ])
  )
  writeLines(lines, out_file)
  invisible(checks)
}
