# ============================================================
# Script: manuscript_maps.R
# Purpose: Export the NC maps and application figures used in Chapter 3 and the SI
#          as clean figures (no in-figure title, subtitle or note).
# Author: Andrew Walther
# Created: 2026-10-25
# Dependencies: ggplot2, sf, tigris, ggpattern (via the sourced map scripts)
# ============================================================
#
# Run from the project root:  Rscript application/code/manuscript_maps.R
#
# Manuscript figures carry their explanation in the LaTeX caption, so this script
# saves title-free versions of figures whose report versions have titles (author,
# 2026-10-25). It reads only frozen production inputs and outputs; nothing is
# simulated, and the hashed real_sud_setup.R / real_sud_simulation.R are untouched.
#
# Writes to paper/{thesis_chapter,ctj_manuscript}/figures/:
#   community_college_service_area_clusters.png   county-to-college map
#   real_sud/observed_years.png, regions.png, design_examples.png   NC maps
#   (design_examples.png: CTJ folder only; thesis_design_examples.R writes the chapter's)
#   real_sud/{primary_yearly_mse,primary_setting_ratios,population_shares,
#             refined_allocation_risk}.pdf        copies of the *_journal.pdf exports
#   (primary_setting_ratios and refined_allocation_risk: CTJ folder only; code/28 builds
#   the chapter's copies)
#
# Restricted data: the maps show only 58-cluster aggregate rates, as before.

rs_companion_define_only <- TRUE          # define functions, do not render the companion
source("application/code/render_real_sud_companion.R")   # rs_maps(), via run_real_sud.R
source("application/code/nc_community_college_map.R")     # get_nc_cc_data(), map builder

ROOT <- "application/results/real_sud_rev_20261002"
FOLDERS <- file.path("paper", c("thesis_chapter", "ctj_manuscript"), "figures")
stopifnot(all(dir.exists(file.path(FOLDERS, "real_sud"))))

#' Save one plot into every manuscript figure folder, failing loud on any write error
#'
#' @param p ggplot object (titles already removed)
#' @param rel Path relative to a figures/ folder, e.g. "real_sud/regions.png"
#' @param w,h Size in inches
#' @param dpi Resolution for PNG output
#' @return Invisibly, the written paths
save_all <- function(p, rel, w, h, dpi = 150) {
  paths <- file.path(FOLDERS, rel)
  # White background: theme_void() leaves the plot background transparent
  for (path in paths) ggplot2::ggsave(path, p, width = w, height = h, dpi = dpi, bg = "white")
  stopifnot(all(file.exists(paths)))
  invisible(paths)
}

no_title <- ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)

# NC maps (same builders and sizes as the report companion) ----
maps <- rs_maps(readRDS(file.path(ROOT, "setup.rds")))
save_all(maps$observed + no_title, "real_sud/observed_years.png", 11, 7)
save_all(maps$regions + no_title, "real_sud/regions.png", 10, 4)
# The chapter's design_examples.png is a thesis-width export from thesis_design_examples.R,
# so this 14 x 8 in copy goes to the CTJ folder only.
for (path in file.path("paper/ctj_manuscript/figures", "real_sud/design_examples.png"))
  ggplot2::ggsave(path, maps$examples + no_title, width = 14, height = 8, dpi = 150, bg = "white")

# County-to-college service-area map (12 x 8 in at 300 dpi, as the existing copy) ----
# tigris reads the project's cached county boundaries rather than downloading
options(tigris_use_cache = TRUE)
Sys.setenv(TIGRIS_CACHE_DIR = normalizePath("application/results/tigris_cache"))
cc_map <- generate_county_borders_map(get_nc_cc_data(), show_title = FALSE)
save_all(cc_map, "community_college_service_area_clusters.png", 12, 8, dpi = 300)

# Application summary figures: the summary step already writes title-free *_journal.pdf ----
THESIS_BUILT <- c("primary_setting_ratios", "refined_allocation_risk")  # written by code/28
for (stem in c("primary_yearly_mse", "primary_setting_ratios", "population_shares",
               "refined_allocation_risk")) {
  src <- file.path(ROOT, "exhibits", paste0(stem, "_journal.pdf"))
  stopifnot(file.exists(src))
  # The chapter's copies of these two are built by code/28_thesis_chapter_figures.R
  # (6.5 in, 8 pt text), so they go to the CTJ folder only, as in copy_figures.sh's
  # THESIS_BUILT list (fixed 2026-11-25: this loop had overwritten the chapter copies).
  folders <- if (stem %in% THESIS_BUILT) FOLDERS[basename(dirname(FOLDERS)) != "thesis_chapter"] else FOLDERS
  ok <- file.copy(src, file.path(folders, "real_sud", paste0(stem, ".pdf")), overwrite = TRUE)
  stopifnot(length(folders) >= 1, all(ok))
}
cat("Clean manuscript maps and application figures written to", paste(FOLDERS, collapse = ", "), "\n")
