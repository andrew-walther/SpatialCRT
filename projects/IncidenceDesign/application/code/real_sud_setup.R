# ============================================================
# Script: real_sud_setup.R
# Purpose: Validate observed inputs and freeze NC design geography.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: sf, spdep, dplyr, digest; existing grid/application helpers
# ============================================================

# Paths and shared, unchanged statistical helpers ----
.rs_code_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(normalizePath(dirname(f)))
  }
  stop("Source real_sud_setup.R by its file path")
})
.rs_app <- normalizePath(file.path(.rs_code_dir, ".."))
.rs_project <- normalizePath(file.path(.rs_app, ".."))
source(file.path(.rs_code_dir, "application_designs.R"))
source(file.path(.rs_code_dir, "cc_mapping_data.R"))
allocation_define_only <- TRUE
allocation_script_dir <- file.path(.rs_project, "code")
source(file.path(allocation_script_dir, "16_allocation_risk.R"))

#' Key-seed a draw without changing prefixes when replication increases
#' @param ... Explicit character/numeric key fields; numeric fields use two decimals.
#' @return Invisibly, the seed, with all RNG kinds specified.
#' @family real_sud
#' @seealso allocation_seed
#' @examples
#' rs_seed("blocks", 2018)
rs_seed <- function(...) allocation_seed("nc-real-sud-v1", ...)

#' Hash inputs and source files by content
#' @param paths Existing absolute file paths.
#' @return Named SHA-256 character vector.
#' @family real_sud
#' @seealso rs_setup
#' @examples
#' rs_hash_files(file.path(.rs_code_dir, "application_designs.R"))
rs_hash_files <- function(paths) {
  if (any(!file.exists(paths))) stop("Missing manifest input: ", paste(paths[!file.exists(paths)], collapse = ", "))
  stats::setNames(vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256"), paths)
}

#' Load observed rates, cached legal geometry and named weights; freeze partitions
#' @param root New study output directory, never the historical synthetic directory.
#' @return List with observed annual tables, sf geometry, weights, fixed regions,
#'   yearly blocks and content hashes. No synthetic incidence is generated.
#' @family real_sud
#' @seealso get_application_designs, rs_manifest
#' @examples
#' s <- rs_setup(file.path(.rs_app, "results", "real_sud_rev_20261002"))
rs_setup <- function(root) {
  rate_path <- file.path(.rs_app, "data", "derived", "sud_cluster_incidence.csv")
  weight_path <- file.path(.rs_app, "data", "nc_cluster_weights.rds")
  weights <- readRDS(weight_path)
  ids <- rownames(weights$queen$W)
  d <- read.csv(rate_path, stringsAsFactors = FALSE)
  years <- as.character(2018:2021)
  expected <- c(years, "2018-2021")
  if (!setequal(unique(d$period), expected)) stop("Unexpected observation periods")
  annual <- lapply(expected, function(y) {
    a <- d[d$period == y, ]
    if (nrow(a) != 58 || anyDuplicated(a$Primary_College) || !setequal(a$Primary_College, ids)) stop("Invalid named panel: ", y)
    a <- a[match(ids, a$Primary_College), ]
    if (any(!is.finite(a$rate_per_100k)) || any(!is.finite(a$person_years)) || any(a$person_years <= 0)) stop("Invalid rates/denominators")
    if (max(abs(a$rate_per_100k - a$deaths / a$person_years * 1e5)) > 1e-10) stop("Rate mismatch")
    a$rank <- rank(a$rate_per_100k, ties.method = "average")
    a$X <- (a$rank - 0.5) / 58
    a
  })
  names(annual) <- expected
  if (sum(d$deaths[d$period %in% years]) != 23523 || sum(annual[[5]]$deaths) != 23523) stop("Unexpected primary count total")
  if (sum(annual[[5]]$person_years) != 25594321 ||
      !identical(as.numeric(annual[[5]]$deaths), Reduce(`+`, lapply(annual[1:4], function(a) as.numeric(a$deaths)))) ||
      !identical(as.numeric(annual[[5]]$person_years), Reduce(`+`, lapply(annual[1:4], function(a) as.numeric(a$person_years))))) stop("Pooled panel disagrees with annual sums")
  for (nb in c("queen", "rook")) {
    w <- weights[[nb]]
    if (!identical(rownames(w$W), ids) || !identical(colnames(w$W), ids) ||
        !isTRUE(all.equal(w$adjacency, t(w$adjacency))) || any(diag(w$adjacency) != 0) ||
        max(abs(rowSums(w$W) - 1)) > 1e-12 || any(w$W < 0) ||
        max(abs(spdep::nb2mat(w$nb, style = "B") - w$adjacency)) > 0 ||
        max(abs(w$W - w$adjacency / rowSums(w$adjacency))) > 1e-12) stop("Invalid named weights: ", nb)
  }
  shp <- file.path(.rs_app, "results", "tigris_cache", "tl_2024_us_county.shp")
  if (!file.exists(shp)) stop("Cached legal-boundary geometry is missing; do not substitute clipped boundaries")
  counties <- sf::st_read(shp, query = "SELECT * FROM tl_2024_us_county WHERE STATEFP = '37'", quiet = TRUE)
  counties <- sf::st_transform(counties, 32119)
  map <- get_cc_mapping_data()
  counties$Primary_College <- map$Primary_College[match(counties$NAME, map$NAME)]
  if (nrow(counties) != 100 || anyNA(counties$Primary_College)) stop("Incomplete county geography/mapping")
  clusters <- dplyr::summarise(dplyr::group_by(counties, Primary_College),
                             n_counties = dplyr::n(), .groups = "drop")
  clusters <- clusters[match(ids, clusters$Primary_College), ]
  # Recheck geometry, but cached W is the authority and is never overwritten.
  for (nb in c("queen", "rook")) {
    gnb <- spdep::poly2nb(clusters, queen = nb == "queen", row.names = ids)
    A <- spdep::nb2mat(gnb, style = "B")
    if (max(abs(A - weights[[nb]]$adjacency)) != 0) stop("Geometry disagrees with cached weights: ", nb)
  }
  points <- sf::st_point_on_surface(sf::st_geometry(clusters))
  xy <- sf::st_coordinates(points)
  coords <- data.frame(x = xy[, 1], y = xy[, 2])
  pop <- annual[[5]]$person_years / 4
  rs_seed("regions")
  region <- make_kmeans_regions(coords, clusters$n_counties, pop, annual[[5]]$X,
                               seed = digest::digest2int("nc-real-sud-v1|regions"))
  blocks <- lapply(years, function(y) {
    rs_seed("blocks", y)
    make_spatial_blocks(coords, annual[[y]]$X, annual[[y]]$person_years)
  })
  names(blocks) <- years
  geometry_paths <- Sys.glob(file.path(dirname(shp), "tl_2024_us_county.*"))
  hashes <- rs_hash_files(c(rate_path, weight_path, geometry_paths,
                           file.path(.rs_app, "data", "cc_name_crosswalk.csv")))
  out <- list(annual = annual, ids = ids, clusters = clusters, coords = coords,
              weights = weights, regions = region$region_id, region_diagnostics = region$diagnostics,
              blocks = blocks, input_hashes = hashes)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, "setup.rds")
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(digest::digest(old), digest::digest(out))) stop("Frozen setup changed; use a new study directory")
  } else saveRDS(out, path)
  write.csv(data.frame(Primary_College = ids, Region = out$regions,
                      out$blocks, check.names = FALSE), file.path(root, "frozen_partitions.csv"), row.names = FALSE)
  out
}

#' Construct a content manifest and reject incompatible checkpoint reuse
#' @param root Profile output directory.
#' @param setup Frozen setup object.
#' @param config Statistical/replication configuration.
#' @return Manifest invisibly; existing mismatches stop execution.
#' @family real_sud
#' @seealso rs_hash_files
#' @examples
#' # rs_manifest(root, setup, config)
rs_manifest <- function(root, setup, config) {
  files <- c(file.path(.rs_code_dir, c("real_sud_setup.R", "real_sud_simulation.R", "run_real_sud.R")),
             file.path(.rs_code_dir, "application_designs.R"),
             file.path(.rs_project, "code", c("01_spatial_setup.R", "02_incidence_generation.R",
                                              "03_designs.R", "04_estimation.R", "05_run_simulation.R",
                                              "16_allocation_risk.R")))
  pkgs <- c("sf", "spdep", "spatialreg", "dplyr", "digest", "ggplot2")
  m <- list(config = config, input = setup$input_hashes, source = rs_hash_files(files),
            design = digest::digest(list(setup$regions, setup$blocks, setup$ids)),
            R = R.version.string, versions = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)),
            BLAS = unname(extSoftVersion()["BLAS"]),
            RNG = c("Mersenne-Twister", "Inversion", "Rejection"))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(root, "manifest.rds")
  if (file.exists(p) && !identical(readRDS(p), m)) stop("Manifest mismatch: refuse checkpoint reuse")
  if (!file.exists(p)) saveRDS(m, p)
  invisible(m)
}
