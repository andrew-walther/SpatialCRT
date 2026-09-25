# ============================================================
# Script: test_application_data.R
# Purpose: Verify the SUD incidence aggregation and the cluster contiguity weights.
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: dplyr; sources application/code/{cc_mapping_data,sud_load_data,sud_incidence,sud_reconcile}.R
# ============================================================
#
# Run from anywhere:
#   Rscript projects/IncidenceDesign/application/tests/test_application_data.R
# Stops at the first failed check.
#   [1]-[2] hand-computed fixture: always run
#   [3]     real data: numerator = final_county_sudden num_obs (23,523), plus reconciliation of
#           the secondary Habib (2026) counts with the paper; needs the gitignored sources
#   [4]     weights: reads the tracked application/data/nc_cluster_weights.rds

test_dir <- local({
  file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(getwd())
})
code_dir <- file.path(test_dir, "..", "code")
data_dir <- file.path(test_dir, "..", "data")
suppressPackageStartupMessages(library(dplyr))
for (f in c("cc_mapping_data.R", "sud_load_data.R", "sud_incidence.R", "sud_reconcile.R")) source(file.path(code_dir, f))

check <- function(name, ok) {
  cat(sprintf("  %-66s %s\n", name, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("Check failed: ", name, call. = FALSE)
}
fails <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

# [1] Rate formulas on a hand-computed fixture ----
# Counties A, B -> college X; C -> college Y. Populations change by year, so the pooled
# rate differs from the mean of yearly rates and a wrong formula would fail.
#   X: 2018 (10+5)/15000 = 100; 2019 (30+5)/25000 = 140; pooled 50/40000 = 125;
#      mean of rates 120; total 50/(40000/2) = 250
#   Y: 2018 2/1000 = 200; 2019 1/4000 = 25; pooled 3/5000 = 60; mean 112.5; total 120
cat("\n[1] Rate formulas (fixture)\n")
fixture <- data.frame(county_name = rep(c("A", "B", "C"), each = 2), year = rep(2018:2019, 3),
                      deaths = c(10, 30, 5, 5, 2, 1), pop_18_64 = c(10000, 20000, 5000, 5000, 1000, 4000))
fixture_map <- data.frame(NAME = c("A", "B", "C"), Primary_College = c("X", "X", "Y"),
                          Full_College_Name = c("X College", "X College", "Y College"))
cl <- aggregate_sud_to_clusters(fixture, fixture_map)
row <- function(college, period) cl[cl$Primary_College == college & cl$period == period, ]

check("yearly: X = 100 (2018), 140 (2019)",
      isTRUE(all.equal(c(row("X", "2018")$rate_per_100k, row("X", "2019")$rate_per_100k), c(100, 140))))
check("average: X = sum d / sum P = 125 (mean of yearly rates = 120)",
      isTRUE(all.equal(c(row("X", "2018-2019")$rate_per_100k, row("X", "2018-2019")$mean_of_yearly_rates), c(125, 120))))
check("total: X = sum d / mean annual P = 250",
      isTRUE(all.equal(row("X", "2018-2019")$rate_cumulative_per_100k, 250)))
check("Y: average 60, mean of rates 112.5, total 120; X has 2 counties",
      isTRUE(all.equal(unlist(row("Y", "2018-2019")[c("rate_per_100k", "mean_of_yearly_rates", "rate_cumulative_per_100k")]),
                       c(60, 112.5, 120), check.attributes = FALSE)) && row("X", "2018")$n_counties == 2)
check("single-year input gives yearly rows only",
      all(compute_incidence(fixture[fixture$year == 2018, ], "county_name")$period == "2018"))

# [2] Bad input stops ----
cat("\n[2] Bad input stops\n")
check("unmapped county", fails(aggregate_sud_to_clusters(rbind(fixture, transform(fixture[1, ], county_name = "D")), fixture_map)))
check("group missing a year", fails(compute_incidence(fixture[-1, ], "county_name")))

# [3] Real data: primary numerator and Habib (2026) reconciliation ----
# Numerator (author decision 2026-09-25): final_county_sudden$num_obs, Habib's corrected
# filtering. The older sudden_county_year counts (deaths_habib2026) must still reproduce
# the published paper. Expected values are checked against an independent read of the raw
# file and against hard-coded published/derived numbers, so a column swap would fail.
cat("\n[3] Real data: primary numerator (final_county_sudden) and Habib (2026) reconciliation\n")
habib_path <- file.path(data_dir, "sudden_county_year.csv")
source_path <- file.path(data_dir, "final_county_sudden.csv")
if (!file.exists(habib_path) || !file.exists(source_path)) {
  cat("  SKIPPED: gitignored source files not in application/data/\n")
} else {
  # Loader rejects an incomplete panel (in either file) and FIPS codes that do not join
  habib_raw <- read.csv(habib_path)
  src_raw <- read.csv(source_path, colClasses = c(county_fips = "character"))
  write_tmp <- function(df) { tmp <- tempfile(fileext = ".csv"); write.csv(df, tmp, row.names = FALSE); tmp }
  check("loader: incomplete Habib (2026) panel stops", fails(load_sud_county_data(write_tmp(habib_raw[-1, ]), source_path)))
  check("loader: incomplete final_county_sudden panel stops", fails(load_sud_county_data(habib_path, write_tmp(src_raw[-1, ]))))
  check("loader: FIPS that do not join stop",
        fails(load_sud_county_data(write_tmp(transform(habib_raw, CORES = replace(CORES, CORES == 1, 2))), source_path)))

  county_df <- load_sud_county_data(habib_path, source_path)

  ## Primary numerator = final_county_sudden$num_obs ----
  key_df <- paste(county_df$county_fips, county_df$year)
  key_raw <- paste(src_raw$county_fips, src_raw$year)
  check("primary deaths equal final_county_sudden num_obs in all 400 county-years",
        nrow(county_df) == 400 && identical(as.numeric(county_df$deaths), as.numeric(src_raw$num_obs[match(key_df, key_raw)])))
  check("primary: 23,523 deaths over 25,594,321 person-years",
        sum(county_df$deaths) == 23523 && sum(src_raw$num_obs) == 23523 && sum(county_df$pop_18_64) == 25594321)
  check("primary >= Habib (2026) in every county-year (2,376 more in total)",
        all(county_df$deaths >= county_df$deaths_habib2026) &&
          sum(county_df$deaths) - sum(county_df$deaths_habib2026) == 2376)

  st <- compute_incidence(transform(county_df, state = "NC"), "state")
  raw_yearly <- tapply(src_raw$num_obs, src_raw$year, sum) / tapply(src_raw$pop_18_64, src_raw$year, sum) * 1e5
  check("statewide rates use the primary counts (= raw-file d/P by year)",
        isTRUE(all.equal(st$rate_per_100k[st$n_years == 1], as.numeric(raw_yearly[st$period[st$n_years == 1]]))) &&
          isTRUE(all.equal(st$rate_per_100k[st$n_years > 1], sum(src_raw$num_obs) / sum(src_raw$pop_18_64) * 1e5)))
  check("statewide primary rates 83.9 / 85.9 / 97.1 / 100.5, overall 91.9",
        all(round(st$rate_per_100k, 1) == c(83.9, 85.9, 97.1, 100.5, 91.9)))

  mapping <- get_cc_mapping_data()
  cl_real <- aggregate_sud_to_clusters(county_df, mapping)
  sums <- aggregate(cbind(deaths, person_years) ~ period, cl_real, sum)
  check("58 clusters whose primary totals equal the state totals in every year and pooled",
        length(unique(cl_real$Primary_College)) == 58 && nrow(sums) == 5 &&
          all(sums$deaths == st$deaths[match(sums$period, st$period)]) &&
          all(sums$person_years == st$person_years[match(sums$period, st$period)]) &&
          sums$deaths[sums$period == "2018-2021"] == 23523)
  alb_fips <- src_raw$county_fips[src_raw$county_name %in% mapping$NAME[mapping$Primary_College == "College of The Albemarle"]]
  alb <- src_raw[src_raw$county_fips %in% alb_fips, ]
  got <- cl_real[cl_real$Primary_College == "College of The Albemarle" & cl_real$n_years > 1, ]
  check("College of The Albemarle (7 counties) matches a raw-file hand calculation",
        length(unique(alb$county_fips)) == 7 && got$deaths == sum(alb$num_obs) &&
          isTRUE(all.equal(got$rate_per_100k, sum(alb$num_obs) / sum(alb$pop_18_64) * 1e5)))

  ## Secondary column reproduces Habib (2026) ----
  st_hab <- compute_incidence(transform(county_df, state = "NC"), "state", "deaths_habib2026")
  check("Habib (2026) column: 21,147 deaths, = sudden_county_year num_obs",
        sum(county_df$deaths_habib2026) == 21147 && sum(habib_raw$num_obs) == 21147)
  check("Habib (2026) column: rates 75.1 / 77.9 / 86.7 / 90.6, overall 82.6",
        all(round(st_hab$rate_per_100k, 1) == c(75.1, 77.9, 86.7, 90.6, 82.6)))
  cty <- compute_incidence(county_df, "county_name", "deaths_habib2026")
  cty <- cty[cty$n_years > 1, ]
  check("Habib (2026) column: county extremes Orange 41.2, Swain 215.6 (paper: 216.0)",
        cty$county_name[which.min(cty$rate_per_100k)] == "Orange" && round(min(cty$rate_per_100k), 1) == 41.2 &&
          cty$county_name[which.max(cty$rate_per_100k)] == "Swain" && round(max(cty$rate_per_100k), 1) == 215.6)
  rep_file <- tempfile(fileext = ".txt")
  rep_checks <- reconcile_sud_counts(county_df, mapping, rep_file)
  check("reconcile_sud_counts(): all four paper checks MATCH on the Habib (2026) column",
        length(rep_checks) == 4 && all(rep_checks) && sum(grepl("MATCH$", readLines(rep_file))) == 4)
}

# [4] Contiguity weights ----
cat("\n[4] Contiguity weights (nc_cluster_weights.rds)\n")
w <- readRDS(file.path(data_dir, "nc_cluster_weights.rds"))
colleges <- sort(unique(get_cc_mapping_data()$Primary_College))
for (type in c("queen", "rook")) {
  A <- w[[type]]$adjacency
  W <- w[[type]]$W
  check(sprintf("%s: 58 x 58, rows/cols named by the 58 colleges", type),
        all(dim(W) == 58) && setequal(rownames(W), colleges) && identical(rownames(W), colnames(W)) &&
          identical(dimnames(A), dimnames(W)))
  check(sprintf("%s: adjacency binary, symmetric, zero diagonal", type),
        all(A %in% c(0, 1)) && isSymmetric(unname(A)) && all(diag(A) == 0))
  check(sprintf("%s: W = A / row sums, every row sums to 1 (no isolated cluster)", type),
        isTRUE(all.equal(W, A / rowSums(A), check.attributes = FALSE)) &&
          isTRUE(all.equal(unname(rowSums(W)), rep(1, 58))))
}
check("rook neighbors are a subset of queen neighbors", all(w$rook$adjacency <= w$queen$adjacency))
check("geography: Durham Tech-Wake Tech neighbors; A-B Tech-Brunswick CC not",
      w$queen$adjacency["Durham Tech", "Wake Tech"] == 1 && w$queen$adjacency["A-B Tech", "Brunswick CC"] == 0)
check("legal boundaries: Albemarle-Martin CC neighbors across the Chowan River",
      w$rook$adjacency["College of The Albemarle", "Martin CC"] == 1)

cat("\nAll checks passed.\n")
