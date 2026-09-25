# ============================================================
# Script: sud_incidence.R
# Purpose: Compute SUD incidence per 100,000 (by year, total, average) for counties
#          and community college clusters.
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: dplyr
# ============================================================
#
# Step 2 of the real-data aggregation (after sud_load_data.R).
#
# Rates, with d = deaths, P = population aged 18-64, t = year, T = number of years:
#   yearly   r_t        = d_t / P_t * 1e5
#   average  r          = sum_t d_t / sum_t P_t * 1e5          per 100k person-years
#   total    r_cum      = sum_t d_t / (sum_t P_t / T) * 1e5    = T * r; deaths per 100k
#                                                                over the whole period
#   diagnostic: mean_t r_t (equal weight per year; differs from r when P varies)

#' Compute yearly and pooled incidence per 100,000 for any grouping
#'
#' @description One function serves counties, clusters and the whole state, so
#'   the rate formulas live in one place. Input rows may be finer than the
#'   grouping (e.g. counties within a cluster); deaths and population are summed
#'   within group-year first. When more than one year is present, one pooled row
#'   per group is added with period "2018-2021".
#'
#' @param df data.frame with `year`, `pop_18_64`, the deaths column and `group_vars`.
#' @param group_vars Character vector of grouping columns.
#' @param deaths_col Numerator column: "deaths" (methodology counts, default) or
#'   "deaths_final_csv" (reconciliation).
#' @return Long data.frame, one row per group and period, with columns:
#'   group_vars, period, n_years, deaths, person_years, rate_per_100k (yearly or
#'   average), rate_cumulative_per_100k (total; NA on yearly rows),
#'   mean_of_yearly_rates (diagnostic; NA on yearly rows).
#' @examples
#' compute_incidence(county_df, "county_name")
#' @family sud_aggregation
compute_incidence <- function(df, group_vars, deaths_col = "deaths") {
  if (anyNA(df[c(group_vars, "year", "pop_18_64", deaths_col)])) {
    stop("Missing values in the incidence input.", call. = FALSE)
  }
  years <- sort(unique(df$year))

  # Yearly rows: r_t = d_t / P_t ----
  yearly <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "year")))) |>
    dplyr::summarise(deaths = sum(.data[[deaths_col]]), person_years = sum(pop_18_64), .groups = "drop") |>
    dplyr::mutate(period = as.character(year), n_years = 1L,
                  rate_per_100k = deaths / person_years * 1e5,
                  rate_cumulative_per_100k = NA_real_, mean_of_yearly_rates = NA_real_)

  # Every group must have every year, or the pooled rates would mix periods
  if (any(table(interaction(yearly[group_vars], drop = TRUE)) != length(years))) {
    stop("Some groups lack data for every year.", call. = FALSE)
  }

  cols <- c(group_vars, "period", "n_years", "deaths", "person_years",
            "rate_per_100k", "rate_cumulative_per_100k", "mean_of_yearly_rates")
  if (length(years) == 1) return(as.data.frame(yearly[cols]))

  # Pooled rows: average = sum d / sum P; total = sum d / mean annual P ----
  pooled <- yearly |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(n_years = dplyr::n(), mean_of_yearly_rates = mean(rate_per_100k),
                     deaths = sum(deaths), person_years = sum(person_years), .groups = "drop") |>
    dplyr::mutate(period = paste0(min(years), "-", max(years)),
                  rate_per_100k = deaths / person_years * 1e5,
                  rate_cumulative_per_100k = deaths / (person_years / n_years) * 1e5)

  as.data.frame(dplyr::bind_rows(yearly[cols], pooled[cols]))
}

#' Aggregate county-year SUD data to community college clusters
#'
#' @description Assigns each county to its single primary college (the mapping
#'   partitions 100 counties into 58 colleges) and computes cluster incidence with
#'   `compute_incidence()`. Stops if a county is unmapped, a mapped county is
#'   missing from the data, or a county-year appears twice.
#'
#' @param county_df County-year data.frame from `load_sud_county_data()`.
#' @param mapping County-to-college table with NAME, Primary_College,
#'   Full_College_Name (from `get_cc_mapping_data()` in cc_mapping_data.R).
#' @param deaths_col Numerator column passed to `compute_incidence()`.
#' @return Long data.frame, one row per college and period: Primary_College,
#'   Full_College_Name, n_counties, county_names, then the
#'   `compute_incidence()` columns.
#' @examples
#' aggregate_sud_to_clusters(county_df, get_cc_mapping_data())
#' @family sud_aggregation
#' @seealso [compute_incidence()]
aggregate_sud_to_clusters <- function(county_df, mapping, deaths_col = "deaths") {
  if (any(duplicated(county_df[c("county_name", "year")]))) {
    stop("Duplicate county-year rows.", call. = FALSE)
  }
  unmapped <- setdiff(county_df$county_name, mapping$NAME)
  absent <- setdiff(mapping$NAME, county_df$county_name)
  if (length(unmapped) > 0) stop("Counties without a college: ", paste(unmapped, collapse = ", "), call. = FALSE)
  if (length(absent) > 0) stop("Mapped counties missing from the data: ", paste(absent, collapse = ", "), call. = FALSE)

  joined <- merge(county_df, mapping[c("NAME", "Primary_College", "Full_College_Name")],
                  by.x = "county_name", by.y = "NAME")
  incidence <- compute_incidence(joined, c("Primary_College", "Full_College_Name"), deaths_col)

  # Cluster descriptors: how many counties, and which
  info <- mapping |>
    dplyr::group_by(Primary_College) |>
    dplyr::summarise(n_counties = dplyr::n(), county_names = paste(sort(NAME), collapse = "; "))

  out <- merge(info, incidence, by = "Primary_College")
  out <- out[order(out$Primary_College, out$n_years, out$period),
             c("Primary_College", "Full_College_Name", "n_counties", "county_names",
               setdiff(names(incidence), c("Primary_College", "Full_College_Name")))]
  rownames(out) <- NULL
  out
}
