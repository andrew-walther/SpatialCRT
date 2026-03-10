# 01_grid_setup.R
# Programmatic grid geometry for SpillSpatialDepSim unified simulation.
# Replaces the hardcoded polygons, neighbor lists, and assign_district()
# functions found in SpatialSim_NC_DOC.Rmd, SpatialSim_3x3.Rmd, SpatialSim_3x4.Rmd.
#
# District numbering: row-major (top-left = 1, increases L->R then T->B).
# NOTE: The original SpatialSim_NC_DOC.Rmd used an anomalous column-major scheme.
# This unified code standardizes to row-major. Block indices are computed
# programmatically so no hardcoded values are needed.

library(sf)

#' Build a rectangular grid setup for the spatial simulation.
#'
#' @param n_rows  Number of rows of districts (e.g. 2 for 2x4)
#' @param n_cols  Number of columns of districts (e.g. 4 for 2x4)
#' @param cell_length  Side length of each square district cell (default 100 units)
#' @param points_per_district  Number of simulated subjects per district per iteration (default 20)
#' @return A named list with grid metadata, neighbor structure, and sf objects.
setup_grid <- function(n_rows, n_cols, cell_length = 100, points_per_district = 20) {
  L <- cell_length
  n_districts   <- n_rows * n_cols
  n_trt         <- floor(n_districts / 2)
  grid_width    <- n_cols * L
  grid_height   <- n_rows * L
  points_per_iteration <- n_districts * points_per_district

  # --- Rook neighbors (row-major, 1-indexed from top-left) ---
  # District d is at row r = (d-1) %/% n_cols, col c = (d-1) %% n_cols  (0-indexed)
  rook_neighbors <- vector("list", n_districts)
  names(rook_neighbors) <- as.character(seq_len(n_districts))
  for (d in seq_len(n_districts)) {
    r <- (d - 1) %/% n_cols  # 0-indexed row
    c <- (d - 1) %%  n_cols  # 0-indexed col
    nb <- integer(0)
    if (r > 0)        nb <- c(nb, d - n_cols)  # up
    if (r < n_rows-1) nb <- c(nb, d + n_cols)  # down
    if (c > 0)        nb <- c(nb, d - 1)        # left
    if (c < n_cols-1) nb <- c(nb, d + 1)        # right
    rook_neighbors[[d]] <- nb
  }

  # --- Spatial polygons (sf) ---
  # District d at (r, c) spans x: [c*L, (c+1)*L], y: [(n_rows-1-r)*L, (n_rows-r)*L]
  district_polys <- lapply(seq_len(n_districts), function(d) {
    r <- (d - 1) %/% n_cols
    c <- (d - 1) %%  n_cols
    x0 <- c * L;       x1 <- x0 + L
    y0 <- (n_rows - 1 - r) * L; y1 <- y0 + L
    st_polygon(list(matrix(c(x0,y0, x1,y0, x1,y1, x0,y1, x0,y0), ncol=2, byrow=TRUE)))
  })
  districts_sf <- st_sfc(district_polys, crs = NA_crs_)

  # Centroid label positions
  label_positions <- do.call(rbind, lapply(seq_len(n_districts), function(d) {
    r <- (d - 1) %/% n_cols
    c <- (d - 1) %%  n_cols
    data.frame(district = as.character(d),
               x = (c + 0.5) * L,
               y = (n_rows - 0.5 - r) * L)
  }))
  label_sf <- st_as_sf(label_positions, coords = c("x", "y"), crs = NA_crs_)

  # --- assign_district(x, y) function ---
  # Maps a point coordinate to a district number (row-major).
  # x ∈ [0, grid_width), y ∈ [0, grid_height).
  assign_district <- function(x, y) {
    r <- n_rows - 1L - as.integer(floor(y / L))
    r <- pmax(0L, pmin(r, n_rows - 1L))  # clamp
    c <- as.integer(floor(x / L))
    c <- pmax(0L, pmin(c, n_cols - 1L))  # clamp
    r * n_cols + c + 1L
  }

  list(
    n_rows               = n_rows,
    n_cols               = n_cols,
    n_districts          = n_districts,
    n_trt                = n_trt,
    cell_length          = L,
    grid_width           = grid_width,
    grid_height          = grid_height,
    points_per_district  = points_per_district,
    points_per_iteration = points_per_iteration,
    rook_neighbors       = rook_neighbors,
    districts_sf         = districts_sf,
    label_sf             = label_sf,
    assign_district      = assign_district
  )
}


#' Compute which combination indices satisfy full block stratification.
#'
#' A valid block-stratified assignment requires that:
#'   (a) no two adjacent districts are BOTH intervention, AND
#'   (b) no two adjacent districts are BOTH control.
#'
#' This is the same criterion used in all three original simulation scripts
#' (is_valid_assignment in NC_DOC.Rmd and is_valid_assignment_full in 3x3/3x4.Rmd).
#'
#' @param grid        Output of setup_grid()
#' @param combinations List of integer vectors (each = set of intervention district indices)
#' @return Integer vector of valid combination indices
build_valid_block_indices <- function(grid, combinations) {
  rn <- grid$rook_neighbors
  n  <- grid$n_districts

  is_valid <- function(combo) {
    intervention <- combo
    control      <- setdiff(seq_len(n), combo)

    # Check no two adjacent intervention cells
    for (d in intervention) {
      if (any(rn[[d]] %in% intervention)) return(FALSE)
    }
    # Check no two adjacent control cells
    for (d in control) {
      if (any(rn[[d]] %in% control)) return(FALSE)
    }
    TRUE
  }

  which(vapply(combinations, is_valid, logical(1)))
}
