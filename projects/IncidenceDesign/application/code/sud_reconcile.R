# ============================================================
# Script: sud_reconcile.R
# Purpose: QC report checking the SUD numerator against Habib (2026) and against the
#          alternative counts in final_county_sudden.csv.
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: dplyr (via sud_incidence.R)
# ============================================================
#
# Step 3 of the real-data aggregation (after sud_incidence.R). Background on the two
# count definitions and the open questions for Ashkan Habib: application/README.md.

#' Published statewide figures from Habib (2026) that the numerator must reproduce
SUD_PAPER_TARGETS <- list(
  deaths = 21147,
  person_years = 25594321,
  overall_rate = 82.6,
  yearly_rate = c(`2018` = 75.1, `2019` = 77.9, `2020` = 86.7, `2021` = 90.6)
)

#' Write the numerator reconciliation report
#'
#' @description Four sections: (1) statewide deaths and rates for both count
#'   definitions next to the paper; (2) pass/fail checks against the paper;
#'   (3) how the two definitions differ county-year by county-year; (4) whether
#'   they rank the 58 clusters the same way (Spearman, same top half, same
#'   quartile: the inputs that High Incidence Focus and Balanced Quartiles use).
#'   A mismatch with the paper raises a warning rather than an error, so a revised
#'   source file can still be aggregated and inspected.
#'
#' @param county_df County-year data.frame from `load_sud_county_data()`.
#' @param mapping County-to-college table from `get_cc_mapping_data()`.
#' @param out_file Path of the text report.
#' @return Invisibly, the named logical vector of paper checks.
#' @family sud_aggregation
reconcile_sud_counts <- function(county_df, mapping, out_file) {
  pooled <- paste0(min(SUD_YEARS), "-", max(SUD_YEARS))
  state_df <- transform(county_df, state = "NC")
  st <- compute_incidence(state_df, "state", "deaths")
  st_alt <- compute_incidence(state_df, "state", "deaths_final_csv")
  cl <- aggregate_sud_to_clusters(county_df, mapping, "deaths")
  cl_alt <- aggregate_sud_to_clusters(county_df, mapping, "deaths_final_csv")

  # (2) Checks against the paper, using the methodology counts ----
  yr <- st$period != pooled
  checks <- c(
    deaths       = st$deaths[!yr] == SUD_PAPER_TARGETS$deaths,
    person_years = st$person_years[!yr] == SUD_PAPER_TARGETS$person_years,
    overall_rate = round(st$rate_per_100k[!yr], 1) == SUD_PAPER_TARGETS$overall_rate,
    yearly_rates = all(round(st$rate_per_100k[yr], 1) == SUD_PAPER_TARGETS$yearly_rate[st$period[yr]])
  )
  if (!all(checks)) {
    warning("Numerator no longer matches Habib (2026): ", paste(names(checks)[!checks], collapse = ", "),
            call. = FALSE)
  }

  # (4) Cluster rank agreement between the two definitions, per period ----
  agree <- sapply(unique(cl$period), function(p) {
    a <- cl$rate_per_100k[cl$period == p]
    b <- cl_alt$rate_per_100k[cl_alt$period == p]
    n <- length(a)
    c(spearman = stats::cor(a, b, method = "spearman"),
      top_half = sum((rank(-a) <= n / 2) == (rank(-b) <= n / 2)),
      quartile = sum(ceiling(4 * rank(a) / n) == ceiling(4 * rank(b) / n)),
      n = n)
  })

  diff <- county_df$deaths_final_csv - county_df$deaths
  lines <- c(
    "NC SUD numerator reconciliation (written by run_sud_aggregation.R)",
    "Numerator: sudden_county_year.csv. Alternative: final_county_sudden.csv num_obs.",
    "",
    "1. Statewide deaths and rates per 100,000, ages 18-64",
    sprintf("  %-10s %8s %8s %8s %8s %8s", "period", "deaths", "alt", "rate", "alt", "paper"),
    sprintf("  %-10s %8d %8d %8.1f %8.1f %8.1f", st$period, st$deaths, st_alt$deaths,
            st$rate_per_100k, st_alt$rate_per_100k,
            c(SUD_PAPER_TARGETS$yearly_rate, SUD_PAPER_TARGETS$overall_rate)[st$period]),
    "",
    "2. Numerator checks against Habib (2026)",
    sprintf("  %-13s %s", names(checks), ifelse(checks, "MATCH", "MISMATCH")),
    "",
    "3. County-years, alternative minus numerator",
    sprintf("  equal %d | alternative higher %d | alternative lower %d",
            sum(diff == 0), sum(diff > 0), sum(diff < 0)),
    "",
    "4. Cluster rank agreement, numerator vs alternative (58 clusters)",
    sprintf("  %-10s %9s %9s %9s", "period", "Spearman", "top half", "quartile"),
    sprintf("  %-10s %9.3f %6d/%d %6d/%d", colnames(agree), agree["spearman", ],
            agree["top_half", ], agree["n", ], agree["quartile", ], agree["n", ])
  )
  writeLines(lines, out_file)
  invisible(checks)
}
