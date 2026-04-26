# ============================================================
# Script: 00_design_names.R
# Purpose: Canonical design name mappings for IncidenceDesign.
#          Source this file and call apply_design_names(results)
#          to replace "Design N" column values with ordered short
#          display names. All other scripts/reports should source
#          this rather than defining their own label vectors.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-26
# Dependencies: dplyr
# ============================================================

# Full names (for figure panels, tables, and prose)
DESIGN_FULL_NAMES <- c(
  "Design 1" = "Block Stratified Sampling",
  "Design 2" = "High Incidence Focus",
  "Design 3" = "Saturation Quadrants",
  "Design 4" = "Isolation Buffer",
  "Design 5" = "2x2 Blocking",
  "Design 6" = "Balanced Quartiles",
  "Design 7" = "Balanced Halves",
  "Design 8" = "Incidence-Guided Saturation Quadrants"
)

# Short names (for plot axes and legends where space is limited)
DESIGN_SHORT_NAMES <- c(
  "Design 1" = "Block Stratified",
  "Design 2" = "Hi-Inc. Focus",
  "Design 3" = "Sat. Quadrants",
  "Design 4" = "Iso. Buffer",
  "Design 5" = "2x2 Blocking",
  "Design 6" = "Bal. Quartiles",
  "Design 7" = "Bal. Halves",
  "Design 8" = "Incidence Sat. Quad."
)

# Abbreviations (for very space-constrained contexts)
DESIGN_ABBREVS <- c(
  "Design 1" = "BSS",
  "Design 2" = "HIF",
  "Design 3" = "SatQ",
  "Design 4" = "IsoB",
  "Design 5" = "2x2B",
  "Design 6" = "BalQ",
  "Design 7" = "BalH",
  "Design 8" = "ISQ"
)

# Display order: Blocking -> Stratified -> Saturation (roughly worst -> best)
# Controls factor levels and axis/legend ordering in all plots.
DESIGN_DISPLAY_ORDER <- c(
  "Block Stratified",    # D1 — Blocking
  "2x2 Blocking",        # D5 — Blocking
  "Iso. Buffer",         # D4 — Blocking
  "Hi-Inc. Focus",       # D2 — Stratified
  "Bal. Halves",         # D7 — Stratified
  "Bal. Quartiles",      # D6 — Stratified
  "Sat. Quadrants",      # D3 — Saturation
  "Incidence Sat. Quad." # D8 — Saturation
)

# Conceptual group membership (keyed by short name)
DESIGN_GROUPS <- c(
  "Block Stratified"     = "Blocking",
  "2x2 Blocking"         = "Blocking",
  "Iso. Buffer"          = "Blocking",
  "Hi-Inc. Focus"        = "Stratified",
  "Bal. Halves"          = "Stratified",
  "Bal. Quartiles"       = "Stratified",
  "Sat. Quadrants"       = "Saturation",
  "Incidence Sat. Quad." = "Saturation"
)

#' Apply canonical short display names to the Design column
#'
#' Mutates the Design column from simulation "Design N" strings to ordered
#' factor short display names defined in DESIGN_DISPLAY_ORDER. Call this
#' immediately after loading results so all downstream code and figures
#' automatically use canonical labels.
#'
#' @param results Data frame with a "Design N" character or factor column
#' @return results with Design replaced by an ordered factor of short names
#' @examples
#' results <- apply_design_names(readRDS("sim_results.rds"))
apply_design_names <- function(results) {
  results |>
    dplyr::mutate(
      Design = factor(
        DESIGN_SHORT_NAMES[as.character(Design)],
        levels = DESIGN_DISPLAY_ORDER
      )
    )
}
