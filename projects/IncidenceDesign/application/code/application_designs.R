# ==============================================================================
# Irregular-map treatment designs for the NC SUD application
# ==============================================================================

APPLICATION_DESIGN_FULL_NAMES <- c(
  "1" = "Block Stratified Sampling",
  "2" = "High Incidence Focus",
  "3" = "Saturation Regions",
  "4" = "Isolation Buffer",
  "5" = "2x2 Blocking",
  "6" = "Balanced Quartiles",
  "7" = "Balanced Halves",
  "8" = "Incidence-Guided Saturation Regions"
)

APPLICATION_DESIGN_SHORT_NAMES <- c(
  "1" = "Block Stratified",
  "2" = "Hi-Inc. Focus",
  "3" = "Sat. Regions",
  "4" = "Iso. Buffer",
  "5" = "2x2 Blocking",
  "6" = "Bal. Quartiles",
  "7" = "Bal. Halves",
  "8" = "Incidence Sat. Regions"
)

APPLICATION_DESIGN_DISPLAY_ORDER <- c(
  "Block Stratified",
  "2x2 Blocking",
  "Iso. Buffer",
  "Hi-Inc. Focus",
  "Bal. Halves",
  "Bal. Quartiles",
  "Sat. Regions",
  "Incidence Sat. Regions"
)

APPLICATION_DESIGN_GROUPS <- c(
  "Block Stratified" = "Blocking",
  "2x2 Blocking" = "Blocking",
  "Iso. Buffer" = "Blocking",
  "Hi-Inc. Focus" = "Stratified",
  "Bal. Halves" = "Stratified",
  "Bal. Quartiles" = "Stratified",
  "Sat. Regions" = "Saturation",
  "Incidence Sat. Regions" = "Saturation"
)

application_design_names <- function(design_id = NULL, style = c("full", "short")) {
  style <- match.arg(style)
  names <- if (style == "full") APPLICATION_DESIGN_FULL_NAMES else APPLICATION_DESIGN_SHORT_NAMES
  if (is.null(design_id)) names else names[as.character(design_id)]
}

application_design_metadata <- function(design_id) {
  short_name <- unname(application_design_names(design_id, style = "short"))
  data.frame(
    design_id = design_id,
    design_label = paste("Design", design_id),
    design_name = unname(application_design_names(design_id, style = "full")),
    design_short_name = short_name,
    design_group = unname(APPLICATION_DESIGN_GROUPS[short_name]),
    stringsAsFactors = FALSE
  )
}

application_design_is_deterministic <- function(design_id) {
  design_id %in% c(1, 2)
}

make_balanced_assignment <- function(idx, treated_fraction = 0.5) {
  n <- length(idx)
  n_treated <- max(0, min(n, round(n * treated_fraction)))
  sample(c(rep(1L, n_treated), rep(0L, n - n_treated)))
}

assign_by_region_saturation <- function(region_id, saturations) {
  z <- integer(length(region_id))
  for (region in sort(unique(region_id))) {
    idx <- which(region_id == region)
    z[idx] <- make_balanced_assignment(idx, saturations[as.character(region)])
  }
  z
}

graph_checkerboard_assignment <- function(nb_list) {
  N <- length(nb_list)
  degree <- lengths(nb_list)
  order_idx <- order(degree, decreasing = TRUE)
  color <- rep(NA_integer_, N)

  for (i in order_idx) {
    neighbors <- nb_list[[i]]
    neighbor_colors <- color[neighbors]
    conflicts_if_zero <- sum(neighbor_colors == 0L, na.rm = TRUE)
    conflicts_if_one <- sum(neighbor_colors == 1L, na.rm = TRUE)
    color[i] <- ifelse(conflicts_if_one < conflicts_if_zero, 1L, 0L)
  }

  if (sum(color) > N / 2) {
    color <- 1L - color
  }

  color
}

isolation_buffer_assignment <- function(nb_list) {
  N <- length(nb_list)
  z <- integer(N)
  available <- seq_len(N)

  while (length(available) > 0) {
    target <- if (length(available) == 1) available else sample(available, 1)
    z[target] <- 1L
    blocked <- unique(c(target, nb_list[[target]]))
    available <- setdiff(available, blocked)
  }

  z
}

make_spatial_blocks <- function(coords, incidence, population, target_block_size = 4) {
  N <- nrow(coords)
  k <- max(2, round(N / target_block_size))
  scaled <- scale(cbind(coords$x, coords$y, incidence, log1p(population)))
  stats::kmeans(scaled, centers = k, nstart = 25, iter.max = 100)$cluster
}

score_region_partition <- function(region_id, coords, county_count, population) {
  regions <- sort(unique(region_id))
  tab <- data.frame(
    region_id = regions,
    n_clusters = as.numeric(tabulate(match(region_id, regions), nbins = length(regions))),
    n_counties = as.numeric(tapply(county_count, region_id, sum)[as.character(regions)]),
    population = as.numeric(tapply(population, region_id, sum)[as.character(regions)])
  )

  pop_share <- tab$population / sum(tab$population)
  cluster_share <- tab$n_clusters / sum(tab$n_clusters)
  county_share <- tab$n_counties / sum(tab$n_counties)

  compactness <- mean(vapply(regions, function(region) {
    idx <- which(region_id == region)
    center <- colMeans(coords[idx, c("x", "y"), drop = FALSE])
    mean(sqrt((coords$x[idx] - center[1])^2 + (coords$y[idx] - center[2])^2))
  }, numeric(1)))

  3 * stats::sd(pop_share) +
    1.5 * stats::sd(cluster_share) +
    1.5 * stats::sd(county_share) +
    compactness / max(stats::dist(coords[, c("x", "y")]))
}

make_kmeans_regions <- function(coords,
                                county_count,
                                population,
                                incidence,
                                n_regions = 4,
                                n_candidates = 250,
                                seed = 20260430) {
  set.seed(seed)
  best_region <- NULL
  best_score <- Inf

  scaled_coords <- scale(coords[, c("x", "y")])
  for (candidate in seq_len(n_candidates)) {
    km <- stats::kmeans(
      scaled_coords,
      centers = n_regions,
      nstart = 1,
      iter.max = 100
    )
    score <- score_region_partition(km$cluster, coords, county_count, population)
    if (is.finite(score) && score < best_score) {
      best_score <- score
      best_region <- km$cluster
    }
  }

  region_summary <- data.frame(
    region_id = sort(unique(best_region)),
    n_clusters = as.numeric(tabulate(best_region, nbins = n_regions)),
    n_counties = as.numeric(tapply(county_count, best_region, sum)),
    population = as.numeric(tapply(population, best_region, sum)),
    mean_incidence_rank01 = as.numeric(tapply(incidence, best_region, mean))
  )
  region_summary$population_share <- region_summary$population / sum(region_summary$population)
  region_summary$score <- best_score

  list(region_id = best_region, diagnostics = region_summary)
}

get_application_designs <- function(design_id,
                                    n_resamples,
                                    incidence,
                                    nb_list,
                                    coords,
                                    population,
                                    county_count,
                                    region_id = NULL) {
  N <- length(incidence)
  mat <- matrix(0L, nrow = N, ncol = n_resamples)

  if (design_id == 1) {
    z <- graph_checkerboard_assignment(nb_list)
    mat <- matrix(rep(z, n_resamples), nrow = N, ncol = n_resamples)
  } else if (design_id == 2) {
    cutoff <- stats::median(incidence)
    z <- as.integer(incidence > cutoff)
    if (sum(z) == 0 || sum(z) == N) {
      z[order(incidence, decreasing = TRUE)[seq_len(round(N / 2))]] <- 1L
    }
    mat <- matrix(rep(z, n_resamples), nrow = N, ncol = n_resamples)
  } else if (design_id == 3) {
    stopifnot(!is.null(region_id))
    base_saturations <- c(0.20, 0.40, 0.60, 0.80)
    for (i in seq_len(n_resamples)) {
      sats <- stats::setNames(sample(base_saturations), sort(unique(region_id)))
      mat[, i] <- assign_by_region_saturation(region_id, sats)
    }
  } else if (design_id == 4) {
    for (i in seq_len(n_resamples)) {
      mat[, i] <- isolation_buffer_assignment(nb_list)
    }
  } else if (design_id == 5) {
    block_id <- make_spatial_blocks(coords, incidence, population)
    for (i in seq_len(n_resamples)) {
      z <- integer(N)
      for (block in sort(unique(block_id))) {
        idx <- which(block_id == block)
        z[idx] <- make_balanced_assignment(idx, 0.5)
      }
      mat[, i] <- z
    }
  } else if (design_id == 6) {
    quartile <- dplyr::ntile(incidence, 4)
    for (i in seq_len(n_resamples)) {
      z <- integer(N)
      for (q in sort(unique(quartile))) {
        idx <- which(quartile == q)
        z[idx] <- make_balanced_assignment(idx, 0.5)
      }
      mat[, i] <- z
    }
  } else if (design_id == 7) {
    half <- dplyr::ntile(incidence, 2)
    for (i in seq_len(n_resamples)) {
      z <- integer(N)
      for (h in sort(unique(half))) {
        idx <- which(half == h)
        z[idx] <- make_balanced_assignment(idx, 0.5)
      }
      mat[, i] <- z
    }
  } else if (design_id == 8) {
    stopifnot(!is.null(region_id))
    means <- tapply(incidence, region_id, mean)
    ranked_regions <- names(sort(means, decreasing = TRUE))
    sats <- stats::setNames(rep(NA_real_, length(means)), names(means))
    sats[ranked_regions] <- c(0.80, 0.60, 0.40, 0.20)
    for (i in seq_len(n_resamples)) {
      mat[, i] <- assign_by_region_saturation(region_id, sats)
    }
  } else {
    stop("Unknown design_id: ", design_id, call. = FALSE)
  }

  storage.mode(mat) <- "integer"
  mat
}
