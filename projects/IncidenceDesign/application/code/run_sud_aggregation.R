# ============================================================
# Script: run_sud_aggregation.R
# Purpose: Build county- and cluster-level SUD incidence per 100,000 from the source
#          data and write the derived files.
# Author: Andrew Walther
# Created: 2026-09-25
# Dependencies: dplyr; sources cc_mapping_data.R, sud_load_data.R, sud_incidence.R,
#               sud_reconcile.R
# ============================================================
#
# Run from anywhere:  Rscript projects/IncidenceDesign/application/code/run_sud_aggregation.R
#
# Inputs  (gitignored): application/data/final_county_sudden.csv (numerator: corrected
#                       counts; denominator: pop_18_64), sudden_county_year.csv (Habib
#                       2026 counts, reconciliation only)
# Outputs (gitignored): application/data/derived/
#   sud_cluster_incidence.csv      58 colleges x period {2018, ..., 2021, 2018-2021}
#   sud_county_incidence.csv       100 counties x the same periods
#   sud_reconciliation_report.txt  numerator vs the Habib (2026) counts and the paper

suppressPackageStartupMessages(library(dplyr))

# Locate this script: innermost source() frame (when sourced), else Rscript --file ----
code_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(normalizePath(dirname(ofile)))
  }
  file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(getwd())
})

for (f in c("cc_mapping_data.R", "sud_load_data.R", "sud_incidence.R", "sud_reconcile.R")) {
  source(file.path(code_dir, f))
}

data_dir <- normalizePath(file.path(code_dir, "..", "data"))
out_dir <- file.path(data_dir, "derived")
dir.create(out_dir, showWarnings = FALSE)

# Load, aggregate ----
county_df <- load_sud_county_data(file.path(data_dir, "sudden_county_year.csv"),
                                  file.path(data_dir, "final_county_sudden.csv"))
mapping <- get_cc_mapping_data()
county_inc <- compute_incidence(county_df, c("county_name", "county_fips"))
cluster_inc <- aggregate_sud_to_clusters(county_df, mapping)

# Conservation: in every period, the 58 clusters must add up to the 100 counties ----
totals <- function(d) aggregate(cbind(deaths, person_years) ~ period, d, sum)
if (!isTRUE(all.equal(totals(county_inc), totals(cluster_inc)))) {
  stop("Cluster totals do not match county totals.", call. = FALSE)
}

# Write ----
utils::write.csv(cluster_inc, file.path(out_dir, "sud_cluster_incidence.csv"), row.names = FALSE)
utils::write.csv(county_inc, file.path(out_dir, "sud_county_incidence.csv"), row.names = FALSE)
checks <- reconcile_sud_counts(county_df, mapping, file.path(out_dir, "sud_reconciliation_report.txt"))

pooled_rates <- cluster_inc$rate_per_100k[cluster_inc$n_years > 1]
message(sprintf(
  "%d clusters | %s deaths | %s person-years | pooled cluster rates %.1f-%.1f per 100k | Habib (2026) paper checks: %s\nWrote %s",
  length(pooled_rates), format(sum(county_df$deaths), big.mark = ","),
  format(sum(county_df$pop_18_64), big.mark = ","), min(pooled_rates), max(pooled_rates),
  if (all(checks)) "all match" else "MISMATCH (see report)", out_dir))
