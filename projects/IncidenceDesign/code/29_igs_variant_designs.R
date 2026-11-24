# ============================================================
# Script: 29_igs_variant_designs.R
# Purpose: Variant Incidence-Guided Saturation designs (number of regions R, spread s,
#          guided or unguided) for the sensitivity study; sourced by 29 and 30.
# Author: Andrew Walther
# Created: 2026-11-24
# Dependencies: base R; random_tie_rank() from 03_designs.R
# ============================================================
#
# Spec: docs/plans/2026-11-24-igs-sensitivity-spec.md. Defines VARIANTS, region_index(),
# treated_counts(), igs_variant_designs() and get_designs_with_variants(), a drop-in
# for 03's get_designs() that sends IDs >= 100 here and every other ID to 03 unchanged.

# Variant table ----
# ID 100 + index; the (R = 4, s = 0.3) pair is dropped because Designs 8 and 3 are it.
VARIANTS <- expand.grid(s = c(0.1, 0.2, 0.3, 0.4), R = c(2, 4, 10, 25),
                        guided = c(TRUE, FALSE))
VARIANTS <- VARIANTS[!(VARIANTS$R == 4 & VARIANTS$s == 0.3), ]
VARIANTS$id    <- 100 + seq_len(nrow(VARIANTS))
VARIANTS$label <- sprintf("%s R=%d s=%.1f", ifelse(VARIANTS$guided, "Guided", "Unguided"),
                          VARIANTS$R, VARIANTS$s)
rownames(VARIANTS) <- NULL

#' Region index of each cluster for an equal-size partition of the 10x10 grid
#'
#' @param coords Data frame with integer x, y in 1..10
#' @param R Number of regions: 2 (5x10 halves), 4 (5x5 quadrants), 10 (2x5
#'   rectangles) or 25 (2x2 blocks)
#' @return Integer vector, length nrow(coords), values 1..R, every region the same size
#' @examples region_index(expand.grid(x = 1:10, y = 1:10), 4)
region_index <- function(coords, R) {
  x0 <- coords$x - 1; y0 <- coords$y - 1
  id <- switch(as.character(R),
               "2"  = x0 %/% 5,
               "4"  = (x0 %/% 5) * 2 + y0 %/% 5,
               "10" = (x0 %/% 2) * 2 + y0 %/% 5,
               "25" = (x0 %/% 2) * 5 + y0 %/% 2,
               stop("Unsupported R: ", R))
  id + 1L
}

#' Treated counts per region that hit the target total exactly
#'
#' Largest-remainder rounding of size_r * sat_r: floor every target, then add one to
#' the regions with the largest fractional parts until the total is N/2. Equal
#' remainders are broken at random (the caller's seed), so no region is favored.
#'
#' @param sizes Integer region sizes
#' @param sats Saturations in [0, 1], same length
#' @param total Required total treated (N/2)
#' @return Integer vector of treated counts, sum = total, each within [0, size]
#' @examples treated_counts(rep(25, 4), c(0.8, 0.6, 0.4, 0.2), 50)  # 20 15 10 5
treated_counts <- function(sizes, sats, total) {
  tgt  <- sizes * sats
  n    <- floor(tgt + 1e-9)
  left <- total - sum(n)
  stopifnot(left >= 0, left <= length(n))
  if (left > 0) {
    frac <- tgt - n
    ord  <- order(-round(frac, 9), runif(length(frac)))   # random among equal remainders
    n[ord[seq_len(left)]] <- n[ord[seq_len(left)]] + 1
  }
  stopifnot(sum(n) == total, all(n >= 0), all(n <= sizes))
  n
}

#' Draw assignments for one variant design
#'
#' sat_r = 0.5 + s (1 - 2 (r - 1) / (R - 1)) for rank r. Guided: rank by the
#' region's mean incidence (highest mean -> r = 1; ties at random per draw, as in
#' Design 8). Unguided: random order per draw (as in Design 3).
#'
#' @param v One row of VARIANTS
#' @param n_resamples Number of draws (J)
#' @param N Number of clusters
#' @param incidence Numeric length-N surface X_k
#' @param coords Grid coordinates
#' @return N x n_resamples binary matrix, each column with exactly N/2 treated
igs_variant_designs <- function(v, n_resamples, N, incidence, coords) {
  reg   <- region_index(coords, v$R)
  sizes <- tabulate(reg, nbins = v$R)
  sats  <- 0.5 + v$s * (1 - 2 * (seq_len(v$R) - 1) / (v$R - 1))   # rank 1 .. R
  means <- vapply(seq_len(v$R), function(r) mean(incidence[reg == r]), numeric(1))
  mat <- matrix(0, N, n_resamples)
  for (i in seq_len(n_resamples)) {
    rank_r <- if (v$guided) v$R + 1 - random_tie_rank(means) else sample.int(v$R)
    sat_by_region <- sats[rank_r]
    n_trt <- treated_counts(sizes, sat_by_region, N %/% 2)
    z <- numeric(N)
    for (r in seq_len(v$R)) {
      idx <- which(reg == r)
      z[idx[sample.int(length(idx), n_trt[r])]] <- 1
    }
    mat[, i] <- z
  }
  mat
}

#' get_designs() with the variant IDs added
#'
#' @param design_id Production ID (1-9, sent to 03's get_designs()) or variant ID (101+)
#' @param n_resamples,N,incidence,nb_list,coords As in 03's get_designs()
#' @param prod_fun 03's get_designs(), captured before any override
#' @return N x n_resamples binary assignment matrix
get_designs_with_variants <- function(design_id, n_resamples, N, incidence, nb_list,
                                      coords, prod_fun) {
  if (design_id < 100) return(prod_fun(design_id, n_resamples, N, incidence, nb_list, coords))
  igs_variant_designs(VARIANTS[VARIANTS$id == design_id, ], n_resamples, N, incidence, coords)
}
