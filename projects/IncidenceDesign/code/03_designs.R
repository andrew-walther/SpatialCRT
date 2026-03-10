# ==============================================================================
# 03_designs.R
# Six treatment assignment strategies for the spatial CRT simulation
# ==============================================================================

#' Human-readable names for each design ID
#'
#' @param design_id Optional integer or vector of IDs. If NULL, returns all.
#' @return Named character vector
get_design_names <- function(design_id = NULL) {
  all_names <- c(
    "1" = "Design 1: Checkerboard",
    "2" = "Design 2: High Incidence Focus",
    "3" = "Design 3: Saturation Quadrants",
    "4" = "Design 4: Isolation Buffer",
    "5" = "Design 5: 2x2 Blocking",
    "6" = "Design 6: Center Hotspot",
    "7" = "Design 7: Balanced Quartiles"
  )
  if (is.null(design_id)) return(all_names)
  all_names[as.character(design_id)]
}

#' Check whether a design produces deterministic assignments given incidence
#'
#' Designs 1 (Checkerboard) and 2 (High Incidence Focus) are deterministic:
#' they produce the same Z vector for every resample.
#'
#' @param design_id Integer, design identifier
#' @return Logical
is_design_deterministic <- function(design_id) {
  design_id %in% c(1, 2)
}

#' Generate treatment assignment matrices for a given design
#'
#' @param design_id Integer, design identifier (1-7, excluding 6 by default)
#' @param n_resamples Integer, number of assignment vectors to generate
#' @param N Integer, number of clusters
#' @param incidence Numeric vector length N, baseline incidence values in [0,1]
#' @param nb_list Neighbor list object (from spdep), used for design 4
#' @param coords Data frame with columns x, y for cluster coordinates
#' @return N x n_resamples binary matrix of treatment assignments
get_designs <- function(design_id, n_resamples, N, incidence, nb_list, coords) {
  mat <- matrix(0, nrow = N, ncol = n_resamples)

  if (design_id == 1) {
    # Design 1: Checkerboard — alternating assignment on grid (deterministic)
    cb <- (coords$x + coords$y) %% 2
    mat <- matrix(rep(cb, n_resamples), nrow = N, ncol = n_resamples)

  } else if (design_id == 2) {
    # Design 2: High Incidence Focus — treat top-50% clusters (deterministic)
    med_val <- median(incidence)
    trt_assign <- as.integer(incidence > med_val)
    mat <- matrix(rep(trt_assign, n_resamples), nrow = N, ncol = n_resamples)

  } else if (design_id == 3) {
    # Design 3: Saturation Quadrants — varying density by spatial region
    half <- coords$x[which.max(coords$x)] / 2
    q_id <- ifelse(coords$x <= half & coords$y <= half, 1,
            ifelse(coords$x > half  & coords$y <= half, 2,
            ifelse(coords$x <= half & coords$y > half,  3, 4)))

    for (i in seq_len(n_resamples)) {
      z <- numeric(N)
      current_sats <- sample(c(0.20, 0.40, 0.60, 0.80))
      for (q in 1:4) {
        idx <- which(q_id == q)
        n_trt <- round(length(idx) * current_sats[q])
        z[idx] <- sample(c(rep(1, n_trt), rep(0, length(idx) - n_trt)))
      }
      mat[, i] <- z
    }

  } else if (design_id == 4) {
    # Design 4: Isolation Buffer — no adjacent treated clusters
    for (i in seq_len(n_resamples)) {
      z <- numeric(N)
      available <- 1:N
      while (length(available) > 0) {
        target <- if (length(available) == 1) available else sample(available, 1)
        z[target] <- 1
        neighbors <- nb_list[[target]]
        blocked <- c(target, neighbors)
        available <- setdiff(available, blocked)
      }
      mat[, i] <- z
    }

  } else if (design_id == 5) {
    # Design 5: 2x2 Blocking — 1:1 randomization within spatial blocks
    block_id <- paste((coords$x - 1) %/% 2, (coords$y - 1) %/% 2)
    unique_blocks <- unique(block_id)
    for (i in seq_len(n_resamples)) {
      z <- numeric(N)
      for (b in unique_blocks) {
        idx <- which(block_id == b)
        z[idx] <- sample(c(1, 1, 0, 0))
      }
      mat[, i] <- z
    }

  } else if (design_id == 6) {
    # Design 6: Center Hotspot — concentric rings (optional)
    dist_from_center <- pmax(abs(coords$x - (max(coords$x) + 1) / 2),
                             abs(coords$y - (max(coords$y) + 1) / 2))
    rings <- as.numeric(as.factor(dist_from_center))
    for (i in seq_len(n_resamples)) {
      z <- numeric(N)
      for (r in unique(rings)) {
        idx <- which(rings == r)
        n_half <- round(length(idx) / 2)
        z[idx] <- sample(c(rep(1, n_half), rep(0, length(idx) - n_half)))
      }
      mat[, i] <- z
    }

  } else if (design_id == 7) {
    # Design 7: Balanced Quartiles — stratified by incidence quartile
    quartiles <- dplyr::ntile(incidence, 4)
    for (i in seq_len(n_resamples)) {
      z <- numeric(N)
      for (q in 1:4) {
        idx <- which(quartiles == q)
        n_half <- round(length(idx) / 2)
        z[idx] <- sample(c(rep(1, n_half), rep(0, length(idx) - n_half)))
      }
      mat[, i] <- z
    }

  } else {
    stop("Unknown design_id: ", design_id)
  }

  mat
}
