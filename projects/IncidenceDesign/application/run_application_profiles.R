# ==============================================================================
# Smoke/pilot runner for the NC SUD application
# ==============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(spatialreg)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tigris)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

detect_application_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("--file=", "", file_arg[1])), mustWork = TRUE))
  }

  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(source_file)) {
    return(normalizePath(dirname(source_file), mustWork = TRUE))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

application_dir <- detect_application_dir()
if (!dir.exists(application_dir)) {
  application_dir <- normalizePath(getwd(), mustWork = TRUE)
}
if (basename(application_dir) != "application") {
  application_dir <- file.path(application_dir, "application")
}
project_dir <- normalizePath(file.path(application_dir, ".."), mustWork = TRUE)

source(file.path(application_dir, "cc_mapping_data.R"))
source(file.path(application_dir, "synthetic_sud_data.R"))
source(file.path(application_dir, "application_designs.R"))

get_application_profile <- function(profile = c("smoke", "pilot", "full")) {
  profile <- match.arg(profile)

  profiles <- list(
    smoke = list(
      profile = "smoke",
      tau_vals = 1.0,
      gamma_vals = 0.6,
      rho_vals = 0.20,
      spillover_type = "both",
      n_design_resamples = 2,
      n_outcome_resamples = 2
    ),
    pilot = list(
      profile = "pilot",
      tau_vals = c(1.0, 2.0),
      gamma_vals = c(0.6, 0.8),
      rho_vals = c(0, 0.20),
      spillover_type = "both",
      n_design_resamples = 5,
      n_outcome_resamples = 5
    ),
    full = list(
      profile = "full",
      tau_vals = c(0.8, 1.0, 1.5, 2.0, 3.0),
      gamma_vals = c(0.5, 0.6, 0.7, 0.8),
      rho_vals = c(0, 0.01, 0.20, 0.50),
      spillover_type = "both",
      n_design_resamples = 25,
      n_outcome_resamples = 10
    )
  )

  profiles[[profile]]
}

load_county_population <- function(population_path = NULL) {
  if (is.null(population_path) || !file.exists(population_path)) {
    return(NULL)
  }

  ext <- tolower(tools::file_ext(population_path))
  raw <- if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(population_path)
  } else {
    read.csv(population_path, stringsAsFactors = FALSE)
  }

  names(raw) <- make.names(names(raw))
  county_col <- names(raw)[grepl("county|name", names(raw), ignore.case = TRUE)][1]
  pop_col <- names(raw)[grepl("population|estimate|pop", names(raw), ignore.case = TRUE)][1]

  if (is.na(county_col) || is.na(pop_col)) {
    stop("Population file must include county/name and population/estimate columns.", call. = FALSE)
  }

  data.frame(
    NAME = gsub("\\s+County$", "", raw[[county_col]], ignore.case = TRUE),
    county_population = as.numeric(gsub("[^0-9.]", "", raw[[pop_col]])),
    stringsAsFactors = FALSE
  )
}

download_osbm_county_population <- function(cache_path = file.path(application_dir, "results", "osbm_county_population_2024.csv"),
                                            refresh = FALSE) {
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  url <- paste0(
    "https://demography.osbm.nc.gov/api/explore/v2.1/catalog/datasets/",
    "county-population-totals/exports/csv?lang=en&timezone=America%2FNew_York",
    "&use_labels=false&delimiter=%2C"
  )

  if (refresh || !file.exists(cache_path)) {
    ok <- tryCatch({
      utils::download.file(url, cache_path, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) {
      return(NULL)
    }
  }

  raw <- utils::read.csv(cache_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(raw) <- tolower(names(raw))

  county_col <- "county"
  value_col <- "value"
  year_col <- "year"
  variable_col <- "variable"
  vintage_col <- "vintage"
  estimate_col <- "estimateprojection"

  required <- c(county_col, value_col, year_col, variable_col, vintage_col, estimate_col)
  if (!all(required %in% names(raw))) {
    return(NULL)
  }

  pop <- raw |>
    dplyr::filter(
      .data[[year_col]] == 2024,
      .data[[variable_col]] == "Total Population",
      .data[[vintage_col]] == 2024,
      .data[[estimate_col]] == "Estimate"
    ) |>
    dplyr::transmute(
      NAME = .data[[county_col]],
      county_population = as.numeric(.data[[value_col]])
    )

  if (nrow(pop) < 100 || any(is.na(pop$county_population))) {
    return(NULL)
  }

  pop
}

build_nc_application_clusters <- function(population_path = NULL,
                                          cb = TRUE,
                                          year = NULL) {
  options(tigris_use_cache = TRUE)
  tigris_cache <- file.path(application_dir, "results", "tigris_cache")
  dir.create(tigris_cache, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TIGRIS_CACHE_DIR = tigris_cache)

  counties <- tigris::counties(
    state = "NC",
    cb = cb,
    year = year %||% 2024,
    class = "sf",
    progress_bar = FALSE
  ) |>
    sf::st_transform(32119)

  mapping <- get_cc_mapping_data()
  college_ids <- data.frame(
    Primary_College = sort(unique(mapping$Primary_College)),
    ID = seq_along(sort(unique(mapping$Primary_College)))
  )
  college_ids$Legend_Label <- paste0(college_ids$ID, ") ", college_ids$Primary_College)

  mapping <- dplyr::left_join(mapping, college_ids, by = "Primary_College")
  county_pop <- load_county_population(population_path)
  if (is.null(county_pop)) {
    county_pop <- download_osbm_county_population()
    if (is.null(county_pop)) {
      county_pop <- data.frame(
        NAME = mapping$NAME,
        county_population = 1,
        stringsAsFactors = FALSE
      )
      population_source <- "equal_county_placeholder"
    } else {
      population_source <- "NC OSBM county-population-totals API, 2024 Total Population Estimate"
    }
  } else {
    population_source <- normalizePath(population_path, mustWork = FALSE)
  }

  county_data <- counties |>
    dplyr::left_join(mapping, by = "NAME") |>
    dplyr::left_join(county_pop, by = "NAME")

  if (any(is.na(county_data$Primary_College))) {
    missing <- paste(county_data$NAME[is.na(county_data$Primary_College)], collapse = ", ")
    stop("Missing Community College mapping for counties: ", missing, call. = FALSE)
  }

  county_data$county_population[is.na(county_data$county_population)] <- 1

  clusters <- county_data |>
    dplyr::group_by(ID, Legend_Label, Primary_College) |>
    dplyr::summarise(
      n_counties = dplyr::n(),
      population = sum(county_population, na.rm = TRUE),
      county_names = paste(sort(NAME), collapse = "; "),
      geometry = sf::st_union(geometry),
      .groups = "drop"
    ) |>
    dplyr::arrange(ID)

  clusters$cluster_id <- clusters$ID
  clusters$population_source <- population_source

  list(counties = county_data, clusters = clusters)
}

get_cluster_coords <- function(clusters_sf) {
  points <- suppressWarnings(sf::st_point_on_surface(clusters_sf))
  xy <- sf::st_coordinates(points)
  data.frame(x = xy[, 1], y = xy[, 2])
}

simulate_application_outcomes <- function(W, X_matrix, Z, tau, gamma, rho, sigma = 1, beta = 1) {
  N <- nrow(W)
  I_mat <- diag(N)

  # SDM/SAR outcome multiplier: spatial dependence in outcomes is generated
  # through `(I - rho W)^-1`, matching the original simulation.
  inv_mat <- solve(I_mat - rho * W)

  # Neighbor exposure under spillover_type = "both": every cluster receives
  # the weighted average of treated neighbors, regardless of own treatment.
  WZ <- as.vector(W %*% Z)
  spill_term <- gamma * WZ

  epsilon <- matrix(stats::rnorm(N * ncol(X_matrix), 0, sigma), nrow = N)
  baseline_part <- inv_mat %*% (beta * X_matrix + epsilon)
  treatment_part <- as.vector(inv_mat %*% (tau * Z + spill_term))
  Y_sim <- sweep(baseline_part, 1, treatment_part, "+")

  list(Y_sim = Y_sim, spill_term = spill_term)
}

estimate_application_tau <- function(Y_sim, Z, spill_term, X_matrix, listw) {
  n_outcomes <- ncol(Y_sim)
  estimates <- rep(NA_real_, n_outcomes)
  se <- rep(NA_real_, n_outcomes)
  ci_lower <- rep(NA_real_, n_outcomes)
  ci_upper <- rep(NA_real_, n_outcomes)
  converged <- rep(FALSE, n_outcomes)

  if (sum(Z) == 0 || sum(Z) == length(Z)) {
    return(data.frame(
      outcome_resample = seq_len(n_outcomes),
      estimate = estimates,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      converged = converged
    ))
  }

  for (k in seq_len(n_outcomes)) {
    df <- data.frame(
      Y = Y_sim[, k],
      Z = Z,
      Spill = spill_term,
      X = X_matrix[, k]
    )

    fit <- tryCatch(
      suppressWarnings(spatialreg::lagsarlm(Y ~ Z + Spill + X, data = df, listw = listw, quiet = TRUE)),
      error = function(e) NULL
    )

    if (!is.null(fit) && "Z" %in% names(stats::coef(fit))) {
      fit_summary <- tryCatch(summary(fit), error = function(e) NULL)
      if (!is.null(fit_summary) && "Z" %in% rownames(fit_summary$Coef)) {
        estimates[k] <- unname(stats::coef(fit)["Z"])
        se[k] <- fit_summary$Coef["Z", "Std. Error"]
        ci_lower[k] <- estimates[k] - stats::qnorm(0.975) * se[k]
        ci_upper[k] <- estimates[k] + stats::qnorm(0.975) * se[k]
        converged[k] <- TRUE
      }
    }
  }

  data.frame(
    outcome_resample = seq_len(n_outcomes),
    estimate = estimates,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    converged = converged
  )
}

summarise_application_results <- function(iteration_results) {
  iteration_results |>
    dplyr::group_by(
      profile, design_id, design_label, design_name, design_short_name,
      design_group, tau, gamma, rho, spillover_type,
      n_design_resamples, n_outcome_resamples
    ) |>
    dplyr::summarise(
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate - tau, na.rm = TRUE),
      sd = stats::sd(estimate, na.rm = TRUE),
      mse = mean((estimate - tau)^2, na.rm = TRUE),
      coverage = mean(ci_lower <= tau & ci_upper >= tau, na.rm = TRUE),
      power = mean(ci_lower > 0, na.rm = TRUE),
      fail_rate = mean(!converged | is.na(estimate)),
      n_valid_est = sum(converged & !is.na(estimate)),
      .groups = "drop"
    )
}

plot_incidence_map <- function(incidence_sf, out_file) {
  p <- ggplot2::ggplot(incidence_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = sud_rate_per_100k), color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_viridis_c(name = "SUD rate\nper 100k") +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Synthetic Placeholder SUD Incidence",
      subtitle = "Poisson SAR surface over NC Community College service areas"
    )
  ggplot2::ggsave(out_file, p, width = 8, height = 5.5, dpi = 300)
  p
}

plot_region_map <- function(incidence_sf, region_id, out_file) {
  plot_sf <- incidence_sf
  plot_sf$region_id <- factor(region_id)
  p <- ggplot2::ggplot(plot_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = region_id), color = "white", linewidth = 0.25) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "K-Means Saturation Regions",
      subtitle = "Four centroid regions selected with population/count balance diagnostics",
      fill = "Region"
    )
  ggplot2::ggsave(out_file, p, width = 8, height = 5.5, dpi = 300)
  p
}

run_application_profile <- function(profile = c("smoke", "pilot", "full"),
                                    population_path = NULL,
                                    output_root = file.path(application_dir, "results"),
                                    seed = 20260430,
                                    allow_full = identical(Sys.getenv("APPLICATION_RUN_FULL"), "TRUE"),
                                    resume = TRUE) {
  profile <- match.arg(profile)
  config <- get_application_profile(profile)

  if (profile == "full" && !allow_full) {
    stop(
      "The full profile is defined but requires explicit confirmation. ",
      "Set APPLICATION_RUN_FULL=TRUE or call run_application_profile('full', allow_full = TRUE).",
      call. = FALSE
    )
  }

  out_dir <- file.path(output_root, profile)
  chunk_dir <- file.path(out_dir, "scenario_chunks")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  spatial <- build_nc_application_clusters(population_path = population_path)
  clusters <- spatial$clusters

  incidence_sf <- generate_synthetic_sud_data(
    clusters,
    cluster_id_col = "cluster_id",
    population_col = "population",
    seed = seed
  )

  weights <- build_application_spatial_weights(incidence_sf, queen = TRUE)
  coords <- get_cluster_coords(incidence_sf)
  region_obj <- make_kmeans_regions(
    coords = coords,
    county_count = incidence_sf$n_counties,
    population = incidence_sf$population,
    incidence = incidence_sf$incidence_rank01,
    seed = seed
  )

  plot_incidence_map(incidence_sf, file.path(out_dir, "figures", "synthetic_incidence_map.png"))
  plot_region_map(incidence_sf, region_obj$region_id, file.path(out_dir, "figures", "kmeans_regions_map.png"))

  N <- nrow(incidence_sf)
  X_matrix <- matrix(
    rep(incidence_sf$incidence_rank01, config$n_outcome_resamples),
    nrow = N,
    ncol = config$n_outcome_resamples
  )

  scenario_grid <- tidyr::expand_grid(
    tau = config$tau_vals,
    gamma = config$gamma_vals,
    rho = config$rho_vals,
    design_id = 1:8
  )

  iteration_results <- list()
  result_idx <- 1L
  total <- nrow(scenario_grid)

  for (scenario_idx in seq_len(total)) {
    scenario <- scenario_grid[scenario_idx, ]
    chunk_file <- file.path(
      chunk_dir,
      sprintf("scenario_%03d_tau_%s_gamma_%s_rho_%s_design_%d.csv",
              scenario_idx,
              gsub("\\.", "p", sprintf("%.2f", scenario$tau)),
              gsub("\\.", "p", sprintf("%.2f", scenario$gamma)),
              gsub("\\.", "p", sprintf("%.2f", scenario$rho)),
              scenario$design_id)
    )

    if (resume && file.exists(chunk_file)) {
      iteration_results[[result_idx]] <- utils::read.csv(chunk_file, stringsAsFactors = FALSE)
      result_idx <- result_idx + 1L
      next
    }

    message(sprintf(
      "[%s %03d/%03d] tau=%.2f gamma=%.2f rho=%.2f design=%d",
      profile, scenario_idx, total, scenario$tau, scenario$gamma,
      scenario$rho, scenario$design_id
    ))

    set.seed(seed + scenario_idx)
    Z_matrix <- get_application_designs(
      design_id = scenario$design_id,
      n_resamples = config$n_design_resamples,
      incidence = incidence_sf$incidence_rank01,
      nb_list = weights$nb,
      coords = coords,
      population = incidence_sf$population,
      county_count = incidence_sf$n_counties,
      region_id = region_obj$region_id
    )

    scenario_results <- list()
    scenario_result_idx <- 1L
    design_meta <- application_design_metadata(scenario$design_id)

    for (design_resample in seq_len(config$n_design_resamples)) {
      Z <- Z_matrix[, design_resample]
      outcome <- simulate_application_outcomes(
        W = weights$W,
        X_matrix = X_matrix,
        Z = Z,
        tau = scenario$tau,
        gamma = scenario$gamma,
        rho = scenario$rho
      )

      est <- estimate_application_tau(
        Y_sim = outcome$Y_sim,
        Z = Z,
        spill_term = outcome$spill_term,
        X_matrix = X_matrix,
        listw = weights$listw
      )

      est$profile <- profile
      est$scenario_id <- scenario_idx
      est$design_resample <- design_resample
      est$design_id <- scenario$design_id
      est$design_label <- design_meta$design_label
      est$design_name <- design_meta$design_name
      est$design_short_name <- design_meta$design_short_name
      est$design_group <- design_meta$design_group
      est$tau <- scenario$tau
      est$gamma <- scenario$gamma
      est$rho <- scenario$rho
      est$spillover_type <- config$spillover_type
      est$n_design_resamples <- config$n_design_resamples
      est$n_outcome_resamples <- config$n_outcome_resamples
      est$n_treated <- sum(Z)
      est$n_control <- length(Z) - sum(Z)

      scenario_results[[scenario_result_idx]] <- est
      scenario_result_idx <- scenario_result_idx + 1L
    }

    scenario_results <- dplyr::bind_rows(scenario_results)
    utils::write.csv(scenario_results, chunk_file, row.names = FALSE)
    iteration_results[[result_idx]] <- scenario_results
    result_idx <- result_idx + 1L
  }

  iteration_results <- dplyr::bind_rows(iteration_results)
  summary_results <- summarise_application_results(iteration_results)

  utils::write.csv(sf::st_drop_geometry(incidence_sf), file.path(out_dir, "synthetic_incidence.csv"), row.names = FALSE)
  utils::write.csv(region_obj$diagnostics, file.path(out_dir, "region_balance_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(iteration_results, file.path(out_dir, "iteration_results.csv"), row.names = FALSE)
  utils::write.csv(summary_results, file.path(out_dir, "summary_results.csv"), row.names = FALSE)
  saveRDS(
    list(
      config = config,
      clusters = incidence_sf,
      weights = weights,
      region = region_obj,
      iteration_results = iteration_results,
      summary_results = summary_results
    ),
    file.path(out_dir, paste0("application_", profile, "_results.rds"))
  )

  invisible(summary_results)
}

run_requested_application_profiles <- function(profiles = c("smoke", "pilot"),
                                               population_path = NULL,
                                               seed = 20260430,
                                               ...) {
  lapply(profiles, run_application_profile, population_path = population_path, seed = seed, ...)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  profiles <- if (length(args) > 0) strsplit(args[[1]], ",", fixed = TRUE)[[1]] else c("smoke", "pilot")
  run_requested_application_profiles(profiles = profiles)
}
