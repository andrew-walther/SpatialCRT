# ============================================================
# Script: sud_load_data.R
# Purpose: Load NC county-year SUD deaths and population at risk into one table.
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: base R only
# ============================================================
#
# Step 1 of the real-data aggregation (run_sud_aggregation.R runs all steps):
#   sud_load_data.R  -> sud_incidence.R -> sud_reconcile.R
#
# Two source files, both gitignored (restricted death-certificate data; public repo):
#   final_county_sudden.csv  NUMERATOR + DENOMINATOR. num_obs = sudden unexpected
#                            out-of-hospital deaths, ages 18-64, from Habib's corrected
#                            case filtering (23,523 deaths); pop_18_64 = SEER population
#                            aged 18-64 for each county and year.
#   sudden_county_year.csv   RECONCILIATION ONLY. CORES (3-digit county FIPS), DOD_YR,
#                            num_obs = the older counts behind Habib (2026), which also
#                            excluded heart-failure deaths (21,147 deaths). Kept as
#                            deaths_habib2026 so the published numbers stay reproducible.
# Decision (author, 2026-09-25; Habib, personal communication): the corrected counts are
# the numerator. The corrected filtering no longer excludes heart-failure deaths, because
# adjudicated sudden cardiac deaths overlap non-negligibly with heart-failure patients.
# The files join on 5-digit FIPS: county_fips = "37" + 3-digit CORES.

SUD_YEARS <- 2018:2021

#' Load and join county-year SUD deaths and population aged 18-64
#'
#' @description Reads both source files, converts CORES to 5-digit FIPS, and joins
#'   on FIPS and year. Stops unless each file is a complete 100-county x 4-year
#'   panel and every row joins, so a bad input can never pass silently.
#'
#' @param habib2026_path Path to `sudden_county_year.csv` (reconciliation counts).
#' @param source_path Path to `final_county_sudden.csv` (numerator and denominator).
#' @return data.frame, one row per county-year (400 rows): county_name,
#'   county_fips, year, deaths (numerator, corrected counts), deaths_habib2026
#'   (reconciliation only), pop_18_64 (denominator).
#' @examples
#' county_df <- load_sud_county_data("data/sudden_county_year.csv",
#'                                   "data/final_county_sudden.csv")
#' @family sud_aggregation
load_sud_county_data <- function(habib2026_path, source_path) {
  counts <- utils::read.csv(habib2026_path)
  src <- utils::read.csv(source_path, colClasses = c(county_fips = "character"))
  require_columns(counts, c("CORES", "DOD_YR", "num_obs"), habib2026_path)
  require_columns(src, c("county_name", "county_fips", "year", "num_obs", "pop_18_64"), source_path)

  # Put both files on the same keys: 5-digit FIPS + integer year
  counts <- data.frame(county_fips = sprintf("37%03d", as.integer(counts$CORES)),
                       year = as.integer(counts$DOD_YR),
                       deaths_habib2026 = counts$num_obs)
  src <- data.frame(county_name = src$county_name,
                    county_fips = src$county_fips,
                    year = as.integer(src$year),
                    deaths = src$num_obs,
                    pop_18_64 = src$pop_18_64)

  check_panel(counts, "deaths_habib2026", habib2026_path)
  check_panel(src, c("deaths", "pop_18_64"), source_path)

  joined <- merge(src, counts, by = c("county_fips", "year"))
  if (nrow(joined) != nrow(counts)) {
    stop(sprintf("Only %d of %d county-years joined on FIPS + year.", nrow(joined), nrow(counts)),
         call. = FALSE)
  }

  joined <- joined[order(joined$county_name, joined$year),
                   c("county_name", "county_fips", "year", "deaths", "deaths_habib2026", "pop_18_64")]
  rownames(joined) <- NULL
  joined
}

#' Stop if a data frame lacks required columns
#'
#' @param df data.frame to check.
#' @param cols Required column names.
#' @param label File name used in the error message.
#' @return Invisibly TRUE.
#' @family sud_aggregation
require_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop("Missing column(s) in ", label, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless a table is a complete 100-county x 2018-2021 panel
#'
#' @description Checks: no duplicate county-years, exactly the years in
#'   `SUD_YEARS`, 100 counties x 4 years = 400 rows, and no missing or negative
#'   values in `value_cols`.
#'
#' @param df data.frame with county_fips and year.
#' @param value_cols Columns that must be complete and non-negative.
#' @param label File name used in the error message.
#' @return Invisibly TRUE.
#' @family sud_aggregation
check_panel <- function(df, value_cols, label) {
  key <- paste(df$county_fips, df$year)
  if (any(duplicated(key))) stop("Duplicate county-year rows in ", label, call. = FALSE)
  if (!setequal(df$year, SUD_YEARS)) stop("Unexpected years in ", label, call. = FALSE)
  if (length(unique(df$county_fips)) != 100 || nrow(df) != 100 * length(SUD_YEARS)) {
    stop(label, " is not a complete 100-county x ", length(SUD_YEARS), "-year panel (",
         nrow(df), " rows).", call. = FALSE)
  }
  for (v in value_cols) {
    if (anyNA(df[[v]]) || any(df[[v]] < 0)) stop("Missing or negative ", v, " in ", label, call. = FALSE)
  }
  invisible(TRUE)
}
