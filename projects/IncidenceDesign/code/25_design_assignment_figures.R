# ============================================================
# Script: 25_design_assignment_figures.R
# Purpose: Example treatment assignments for each design on one illustrative
#          incidence surface: per-design thumbnails for the CTJ hybrid Table 1,
#          a 2x4 panel figure for the SI and Chapter 3, and a wide slide version.
# Author: Andrew Walther
# Created: 2026-10-25
# Dependencies: ggplot2, patchwork, dplyr (via 03_designs.R), spdep (via 01)
# ============================================================
#
# Run from the project root:  Rscript code/25_design_assignment_figures.R
#
# Inputs: none (draws one illustrative surface and one assignment per design with
# the production functions in 01_spatial_setup.R, 02_incidence_generation.R and
# 03_designs.R; no simulation output is read).
#
# Outputs:
#   paper/ctj_manuscript/figures/compact/design_thumbs/thumb_<id>.pdf
#       one map per CTJ Table 1 row (id = design number; "incidence" = surface)
#   paper/ctj_manuscript/figures/design_assignments_panels.pdf      SI figure
#   paper/thesis_chapter/figures/design_assignments_panels.pdf      Chapter 3 figure
#   presentation/figures/design_assignments_panels_wide.{pdf,png}   slides
#
# Conventions: queen contiguity; treated clusters solid blue, controls white;
# overlays show each design's structure (2x2 blocks, incidence quartiles,
# saturation quadrants with their treated shares); no in-figure title or note
# (author, 2026-10-25); design names from the manuscript name set.
#
# The surface is illustrative, NOT one of the simulation's surfaces: a spatially
# autocorrelated surface X = Phi{(I - rho_X W_Q)^{-1} u} with rho_X = 0.5, chosen
# so that regional differences in burden (and hence High Incidence Focus and
# Incidence-Guided Saturation) are visible. One fixed seed drives every draw.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("code/01_spatial_setup.R")
source("code/02_incidence_generation.R")
source("code/03_designs.R")

pdf(NULL)  # keep ggplot from opening a stray Rplots.pdf

# Setup ----

SEED    <- 2026
RHO_X   <- 0.5
TRT_COL <- "#2166AC"   # treated fill; controls are white

OUT_THUMB <- "paper/ctj_manuscript/figures/compact/design_thumbs"
OUT_SI    <- "paper/ctj_manuscript/figures/design_assignments_panels.pdf"
OUT_CH3   <- "paper/thesis_chapter/figures/design_assignments_panels.pdf"
OUT_SLIDE <- "presentation/figures/design_assignments_panels_wide"
for (d in c(OUT_THUMB, dirname(OUT_SI), dirname(OUT_CH3), dirname(OUT_SLIDE))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

#' Replace ASCII hyphens with soft hyphens: R's pdf() draws "-" as a minus sign,
#' while "­" prints as a true hyphen (as in 23_methods_extension_figures.R)
#' @param x Character vector of labels
#' @return The labels with "-" replaced by "­"
hy <- function(x) gsub("-", "­", x, fixed = TRUE)

# Design numbers (03_designs.R) and manuscript names
DESIGN_NAME <- c("1" = "Checkerboard", "2" = "High Incidence Focus",
                 "4" = "Isolation Buffer", "5" = "2×2 Blocking",
                 "6" = "Balanced Quartiles", "8" = "Incidence-Guided Saturation",
                 "9" = "Simple Random Sampling")
CTJ_ORDER <- c(1, 2, 4, 5, 6, 8, 9)   # CTJ Table 1 row order
CH3_ORDER <- c(8, 6, 4, 2, 5, 1, 9)   # Chapter 3 / SI design table order, then SRS

# Illustrative surface and one assignment per design ----

grid <- build_spatial_grid(10)
co   <- grid$coords
N    <- grid$N_clusters

set.seed(SEED)
inc <- as.numeric(generate_incidence_spatial(N, 1, grid$W_queen, RHO_X))
Z <- sapply(CTJ_ORDER, function(d)
  get_designs(d, 1, N, inc, grid$nb_queen, co)[, 1])
colnames(Z) <- CTJ_ORDER

# Structure behind the overlays, computed the way 03_designs.R does
quad  <- with(co, ifelse(x <= 5 & y <= 5, 1, ifelse(x > 5 & y <= 5, 2,
                  ifelse(x <= 5 & y > 5, 3, 4))))
qmean <- tapply(inc, quad, mean)
qsat  <- c(80, 60, 40, 20)[rank(-qmean)]                 # highest mean -> 80%
quart <- dplyr::ntile(rank(inc, ties.method = "random"), 4)

# Fail loud if the drawn assignments do not follow their rules
n_trt <- colSums(Z)
stopifnot(all(n_trt[c("1", "2", "5", "6", "8", "9")] == 50))
stopifnot(all(sapply(which(Z[, "4"] == 1),
                     function(i) all(Z[grid$nb_queen[[i]], "4"] == 0))))
stopifnot(all(round(tapply(Z[, "8"], quad, mean) * 100) == qsat))
stopifnot(all(tapply(Z[, "6"], quart, sum) %in% c(12, 13)))
stopifnot(all(tapply(Z[, "5"], paste((co$x - 1) %/% 2, (co$y - 1) %/% 2), sum) == 2))
cat(sprintf("Treated: %s\nIsolation Buffer treated %d; quadrant saturations %s\n",
            paste(n_trt, collapse = " "), n_trt[["4"]], paste(qsat, collapse = "/")))

# Panel builders ----

base_df <- data.frame(x = co$x, y = co$y, inc = inc, quart = quart)
frame_theme <- function(legend) {
  theme_void() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.position = legend)
}

#' Map of one design's example assignment with its structural overlay
#'
#' @param d Design number (character or numeric) present in colnames(Z)
#' @param text_size Size (mm) of overlay text
#' @param lw Multiplier for overlay line widths (smaller for thumbnails)
#' @param legend Legend position passed to theme()
#' @return A ggplot object
design_map <- function(d, text_size = 2.6, lw = 1, legend = "none") {
  d  <- as.character(d)
  df <- base_df
  df$z <- factor(Z[, d], levels = c(0, 1))
  p <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = z), color = "grey60", linewidth = 0.3) +
    scale_fill_manual(values = c("0" = "white", "1" = TRT_COL),
                      labels = c("Control", "Treated"), name = NULL, drop = FALSE)
  if (d == "5") {          # 2x2 Blocking: the 25 blocks
    p <- p + geom_vline(xintercept = seq(0.5, 10.5, 2), linewidth = 0.7 * lw) +
      geom_hline(yintercept = seq(0.5, 10.5, 2), linewidth = 0.7 * lw)
  }
  if (d == "6") {          # Balanced Quartiles: each cluster's incidence quartile
    p <- p + geom_text(aes(label = quart, color = z), size = text_size * 1.15,
                       fontface = "bold") +
      scale_color_manual(values = c("0" = "grey25", "1" = "white"), guide = "none")
  }
  if (d == "8") {          # Incidence-Guided Saturation: quadrants and treated shares
    lab <- data.frame(x = c(3, 8, 3, 8), y = c(3, 3, 8, 8), lab = paste0(qsat, "%"))
    p <- p + geom_vline(xintercept = 5.5, linewidth = 1.1 * lw) +
      geom_hline(yintercept = 5.5, linewidth = 1.1 * lw) +
      geom_label(data = lab, aes(x, y, label = lab), inherit.aes = FALSE,
                 size = text_size, linewidth = 0, fill = "white", alpha = 0.7,
                 fontface = "bold")
  }
  p + coord_fixed(expand = FALSE) + frame_theme(legend)
}

#' Grayscale map of the illustrative incidence surface, quadrants dashed
#'
#' @param legend Legend position passed to theme()
#' @return A ggplot object
incidence_map <- function(legend = "none") {
  ggplot(base_df, aes(x, y)) +
    geom_tile(aes(fill = inc), color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 5.5, linewidth = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 5.5, linewidth = 0.5, linetype = "dashed") +
    scale_fill_gradient(low = "white", high = "grey15", limits = c(0, 1),
                        breaks = c(0, 0.5, 1), name = "Incidence") +
    coord_fixed(expand = FALSE) + frame_theme(legend)
}

#' 2x4 panel figure: incidence surface first, then the designs in `order`
#'
#' @param order Design numbers in panel order
#' @param title_size Panel title size (pt)
#' @param wrap Break long names onto two lines
#' @return A patchwork object
panel_figure <- function(order, title_size = 8.5, wrap = TRUE) {
  nm <- hy(DESIGN_NAME[as.character(order)])
  if (wrap) {
    nm <- sub("Incidence­Guided ", "Incidence­Guided\n", nm)
    nm <- sub("Random ", "Random\n", nm)
  }
  titled <- function(p, t) p + ggtitle(t) +
    theme(plot.title = element_text(size = title_size, face = "bold", hjust = 0.5,
                                    lineheight = 0.9, margin = margin(b = 3)))
  plots <- c(list(titled(incidence_map("bottom"), "Incidence (all panels)")),
             Map(function(d, t) titled(design_map(d, legend = "bottom"), t), order, nm))
  wrap_plots(plots, ncol = 4) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.title = element_text(size = 8),
          legend.text = element_text(size = 8), legend.key.size = unit(0.35, "cm"),
          plot.margin = margin(3, 3, 3, 3))
}

# Outputs ----

## CTJ Table 1 thumbnails ----
ggsave(file.path(OUT_THUMB, "thumb_incidence.pdf"), incidence_map(),
       width = 0.85, height = 0.85)
for (d in CTJ_ORDER) {
  ggsave(file.path(OUT_THUMB, sprintf("thumb_%d.pdf", d)),
         design_map(d, text_size = 1.6, lw = 0.55), width = 0.85, height = 0.85)
}

## SI and Chapter 3 panel figure ----
fig <- panel_figure(CH3_ORDER)
ggsave(OUT_SI, fig, width = 7.2, height = 4.2)
ggsave(OUT_CH3, fig, width = 7.2, height = 4.2)

## Slide version (wider, larger titles) ----
ggsave(paste0(OUT_SLIDE, ".pdf"), panel_figure(CH3_ORDER, title_size = 10),
       width = 11, height = 5.8)
# The PNG is rendered from the PDF: raster devices drop the soft hyphens
status <- system2("pdftoppm", c("-png", "-r", "200", "-singlefile",
                                shQuote(paste0(OUT_SLIDE, ".pdf")), shQuote(OUT_SLIDE)))
if (status != 0) stop("pdftoppm failed for ", OUT_SLIDE, ".pdf")

cat("Wrote", length(CTJ_ORDER) + 1, "thumbnails,", OUT_SI, ",", OUT_CH3, "and",
    paste0(OUT_SLIDE, ".{pdf,png}"), "\n")
