# ==============================================================================
# 01_spatial_setup.R
# Grid construction and spatial weight matrices for the CRT simulation
# ==============================================================================

#' Build a regular square lattice grid with rook and queen contiguity structures
#'
#' @param grid_dim Integer, side length of the grid (default 10 -> 100 clusters)
#' @return A list with: coords, N_clusters, nb_rook, nb_queen,
#'         W_rook, W_queen, listw_rook, listw_queen
build_spatial_grid <- function(grid_dim = 10) {
  library(spdep)

  N_clusters <- grid_dim^2
  coords <- expand.grid(x = 1:grid_dim, y = 1:grid_dim)

  # Neighbor lists (rook = edge-sharing, queen = edge+vertex-sharing)
  nb_rook  <- cell2nb(grid_dim, grid_dim, type = "rook")
  nb_queen <- cell2nb(grid_dim, grid_dim, type = "queen")

  # Row-standardized weight matrices
  W_rook  <- nb2mat(nb_rook,  style = "W", zero.policy = TRUE)
  W_queen <- nb2mat(nb_queen, style = "W", zero.policy = TRUE)

  # listw objects for spatialreg::lagsarlm()
  listw_rook  <- mat2listw(W_rook,  style = "W")
  listw_queen <- mat2listw(W_queen, style = "W")

  list(
    coords      = coords,
    N_clusters  = N_clusters,
    grid_dim    = grid_dim,
    nb_rook     = nb_rook,
    nb_queen    = nb_queen,
    W_rook      = W_rook,
    W_queen     = W_queen,
    listw_rook  = listw_rook,
    listw_queen = listw_queen
  )
}

#' Select the active spatial structure (rook or queen) from a grid object
#'
#' @param grid_obj List returned by build_spatial_grid()
#' @param nb_type Character, "rook" or "queen"
#' @return A list with: W, nb, listw
get_active_spatial <- function(grid_obj, nb_type = c("rook", "queen")) {
  nb_type <- match.arg(nb_type)
  if (nb_type == "rook") {
    list(W = grid_obj$W_rook, nb = grid_obj$nb_rook, listw = grid_obj$listw_rook)
  } else {
    list(W = grid_obj$W_queen, nb = grid_obj$nb_queen, listw = grid_obj$listw_queen)
  }
}
