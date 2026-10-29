# ============================================================
# Script: presentation/tools/make_figures.R
# Purpose: Build the slide-sized figures for the Project 1 talk from the saved
#          per-allocation beta-hat MSE results (no hand-typed numbers).
# Author: Andrew Walther
# Created: 2026-10-25 (slide-sized versions 2026-10-28)
# Dependencies: ggplot2, dplyr, ragg
# ============================================================
# Run from the project root (SpillSpatialDepSim/):
#   Rscript presentation/tools/make_figures.R
# Inputs are read through presentation/tools/load_results.R (read-only).
# Outputs go to presentation/figures/. Figures carry no titles or notes; the
# slide title carries the claim. Sizes are for a 1920x1080 slide.

suppressMessages({library(ggplot2); library(dplyr)})
source("presentation/tools/load_results.R")

out_dir <- "presentation/figures"
navy <- "#13294B"; carolina <- "#4B9CD3"; gray <- "#7D8590"; light <- "#E1E5EA"
design_cols <- c("Simple random sampling" = gray, "Block stratified sampling" = navy)

#' Save a figure as a PNG with ragg (sharp text at slide size).
#' @param p ggplot object. @param file File name in out_dir. @param w,h Inches.
#' @return Path, invisibly.
save_png <- function(p, file, w, h) {
  path <- file.path(out_dir, file)
  ggsave(path, p, width = w, height = h, dpi = 150, bg = "white", device = ragg::agg_png)
  invisible(path)
}

theme_slide <- function(base = 19) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.text = element_text(face = "bold", hjust = 0, size = rel(0.95)),
          axis.title = element_text(size = rel(0.9)), legend.text = element_text(size = rel(0.95)))
}
rho_lab <- function(r) ifelse(r == 0, "No spatial dependence (ρ = 0)", "Spatial dependence (ρ = 0.01)")

write.csv(rank_tab, file.path(out_dir, "bss_rank_among_srs.csv"), row.names = FALSE)

# Figure 1: every 2x4 allocation at psi = 0.5, BSS marked ----
# One row per spillover type x rho; square-root x scale, free per row, because
# the four settings differ by an order of magnitude.
d1 <- all_alloc %>% filter(grid == "2x4", Psi == 0.5) %>%
  mutate(panel = factor(sprintf("%s · %s", spill, rho_lab(Rho)),
                        levels = c(sprintf("%s · %s", spill_labels[1], rho_lab(c(0.01, 0))),
                                   sprintf("%s · %s", spill_labels[2], rho_lab(c(0.01, 0))))),
         design = ifelse(is_BSS, "Block stratified sampling", "Simple random sampling"))
set.seed(1)  # vertical jitter only, for legibility
p1 <- ggplot(d1, aes(x = MSE, y = 0)) +
  geom_jitter(data = filter(d1, !is_BSS), aes(colour = design, shape = design),
              width = 0, height = 0.25, alpha = 0.7, size = 3.4) +
  geom_point(data = filter(d1, is_BSS), aes(colour = design, shape = design), fill = navy, size = 6) +
  scale_colour_manual(NULL, values = design_cols,
                      labels = c("Block stratified sampling", "Simple random sampling")) +
  scale_shape_manual(NULL, values = c("Simple random sampling" = 16, "Block stratified sampling" = 23),
                     labels = c("Block stratified sampling", "Simple random sampling")) +
  scale_x_sqrt(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.03)), labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = NULL) +
  facet_wrap(~ panel, ncol = 1, scales = "free_x") +
  labs(x = "Mean squared error of the intervention effect estimate (square-root scale)", y = NULL) +
  theme_slide(24)
save_png(p1, "fig_spread_2x4.png", 13, 8.4)

# Figures 2a/2b: 2x4 MSE by psi, one spillover type each (slide-sized redraw of the
# chapter's box-plot figure). Boxes = the 70 random allocations; diamonds = BSS mean.
mse_by_psi <- function(gr, sp) {
  d <- all_alloc %>% filter(grid == gr, spill == sp) %>%
    mutate(rho = factor(rho_lab(Rho), levels = rho_lab(c(0.01, 0))))
  b <- d %>% filter(is_BSS) %>% group_by(rho, Psi) %>% summarise(MSE = mean(MSE), .groups = "drop")
  ggplot(d, aes(x = factor(Psi), y = MSE)) +
    geom_boxplot(aes(fill = "Simple random sampling"), colour = "#4a5058", width = 0.55,
                 outlier.colour = gray, outlier.alpha = 0.6, outlier.size = 2.4) +
    geom_point(data = b, aes(fill = "Block stratified sampling"), shape = 23, size = 6.5, colour = "white", stroke = 1.2) +
    scale_fill_manual(NULL, values = c("Simple random sampling" = light, "Block stratified sampling" = navy),
                      labels = c("Block stratified sampling (mean)", "Simple random sampling (every allocation)")) +
    facet_wrap(~ rho, nrow = 1, scales = "free_y") +
    labs(x = "Spillover effect ψ", y = "MSE of the intervention effect") + theme_slide(20)
}
save_png(mse_by_psi("2x4", spill_labels[["TrtNoSpill"]]), "fig_mse_2x4_controls.png", 13, 6.2)
save_png(mse_by_psi("2x4", spill_labels[["TrtSpill"]]), "fig_mse_2x4_both.png", 13, 6.2)

# Figure 3: the 3x3 and 3x4 grids side by side, every allocation, both spillover types ----
# Boxes = all random allocations' MSE at each psi (the spread is the point: a trial runs
# one allocation); diamonds = each BSS allocation (2 on 3x4, 6 on 3x3). Square-root y
# scale, free per panel, because the random-allocation tails reach about 2.5.
lab_col <- function(grid, Rho) factor(sprintf("%s, ρ = %s", sub("x", "×", grid), ifelse(Rho == 0, "0", "0.01")),
                                      levels = c("3×3, ρ = 0.01", "3×3, ρ = 0", "3×4, ρ = 0.01", "3×4, ρ = 0"))
lab_row <- function(spill) factor(ifelse(spill == spill_labels[1], "Control Clusters Only", "Control & Intervention\nClusters"),
                                  levels = c("Control Clusters Only", "Control & Intervention\nClusters"))
d3 <- all_alloc %>% filter(grid != "2x4") %>% mutate(col = lab_col(grid, Rho), row = lab_row(spill))
set.seed(3)  # horizontal jitter of the BSS diamonds only
p3 <- ggplot(d3, aes(x = factor(Psi), y = MSE)) +
  geom_boxplot(aes(fill = "Simple random sampling"), colour = "#4a5058", width = 0.6,
               outlier.colour = gray, outlier.alpha = 0.5, outlier.size = 1.6) +
  geom_point(data = filter(d3, is_BSS), aes(fill = "Block stratified sampling"), shape = 23, size = 4.2,
             colour = "white", stroke = 0.8, position = position_jitter(width = 0.12, height = 0)) +
  scale_fill_manual(NULL, values = c("Simple random sampling" = light, "Block stratified sampling" = navy),
                    labels = c("Block stratified sampling (each allocation)", "Simple random sampling")) +
  scale_y_sqrt() + scale_x_discrete(labels = c("0.5", "", "0.7", "")) +
  facet_grid(row ~ col, scales = "free_y") +
  labs(x = "Spillover effect ψ", y = "MSE (square-root scale)") +
  theme_slide(25) + theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.x = unit(1.6, "lines"))
save_png(p3, "fig_mse_larger_grids.png", 16, 7)

# Figure 4: where BSS falls among random allocations, every grid and scenario ----
d4 <- rank_tab %>%
  mutate(setting = factor(sprintf("%s\n%s", spill, ifelse(Rho == 0, "ρ = 0", "ρ = 0.01")),
                          levels = rev(sprintf("%s\n%s", rep(spill_labels, each = 2), c("ρ = 0.01", "ρ = 0")))),
         grid = factor(grid, levels = c("2x4", "3x3", "3x4"), labels = c("2×4", "3×3", "3×4")))
p4 <- ggplot(d4, aes(x = pct_SRS_worse, y = setting)) +
  geom_vline(xintercept = 50, colour = "#b9bec4", linetype = "dashed") +
  geom_point(aes(shape = grid, colour = grid, group = interaction(grid, Psi)), size = 4.2, alpha = 0.9,
             position = position_dodge(width = 0.75)) +
  scale_shape_manual("Grid", values = c(16, 17, 15)) +
  scale_colour_manual("Grid", values = c(navy, carolina, "#7a5c00")) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Random allocations with higher MSE than block stratified sampling (%), one point per ψ", y = NULL) +
  theme_slide(18)
save_png(p4, "fig_bss_rank_all_grids.png", 13, 6.4)

# Figure 5: why a checkerboard makes treatment and spillover collinear ----
# Generic 2x4 rook grid. Spillover indicator (Control Clusters Only): z = 1 for a
# control cell with a treated rook neighbor.
spill_z <- function(trt, nr = 2, nc = 4) {
  g <- expand.grid(r = 1:nr, c = 1:nc); g$x <- trt[cbind(g$r, g$c)]
  g$z <- mapply(function(r, c, x) {
    if (x == 1) return(0L)
    nb <- rbind(c(r - 1, c), c(r + 1, c), c(r, c - 1), c(r, c + 1))
    nb <- nb[nb[, 1] >= 1 & nb[, 1] <= nr & nb[, 2] >= 1 & nb[, 2] <= nc, , drop = FALSE]
    as.integer(any(trt[nb] == 1))
  }, g$r, g$c, g$x)
  g
}
checker <- outer(1:2, 1:4, function(r, c) as.integer((r + c) %% 2 == 0))
blocked <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0), nrow = 2, byrow = TRUE)
d5 <- bind_rows(mutate(spill_z(checker), allocation = "Checkerboard (block stratified)"),
                mutate(spill_z(blocked), allocation = "Treated clusters together (one random allocation)")) %>%
  mutate(cell = case_when(x == 1 ~ "Intervention (x = 1)",
                          z == 1 ~ "Control, exposed to spillover (z = 1)",
                          TRUE ~ "Control, not exposed (z = 0)"))
p5 <- ggplot(d5, aes(x = c, y = -r, fill = cell)) +
  geom_tile(colour = "white", linewidth = 2) +
  geom_text(aes(label = sprintf("x = %d\nz = %d", x, z)), size = 7.5, lineheight = 0.9,
            colour = ifelse(d5$cell == "Intervention (x = 1)", "white", "#1d1f21")) +
  scale_fill_manual(NULL, breaks = c("Intervention (x = 1)", "Control, exposed to spillover (z = 1)",
                                     "Control, not exposed (z = 0)"),
                    values = c("Intervention (x = 1)" = navy, "Control, exposed to spillover (z = 1)" = "#FBAF3D",
                               "Control, not exposed (z = 0)" = light)) +
  facet_wrap(~ allocation, ncol = 2) + coord_equal() +
  theme_void(base_size = 24) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", margin = margin(b = 8)))
save_png(p5, "fig_collinearity_grid.png", 14, 4.6)

# Figure 6: the two designs on each grid (replaces the manuscript's GridsCombined
# panel on the slide: same palette as the other figures, no panel letters or numbers).
# Top: strict block stratified allocation; bottom: one random allocation that puts the
# treated clusters together. Generic row-major grids.
design_tiles <- function(nr, nc, trt, grid_lab, design) {
  g <- expand.grid(r = 1:nr, c = 1:nc); g$id <- (g$r - 1) * nc + g$c
  transform(g, x = as.integer(id %in% trt), grid = grid_lab, design = design)
}
d6 <- bind_rows(
  design_tiles(2, 4, c(1, 3, 6, 8), "2×4 grid", "Block\nstratified"),
  design_tiles(3, 3, c(2, 4, 6, 8), "3×3 grid", "Block\nstratified"),
  design_tiles(3, 4, c(1, 3, 6, 8, 9, 11), "3×4 grid", "Block\nstratified"),
  design_tiles(2, 4, c(1, 2, 5, 6), "2×4 grid", "One random\ndraw"),
  design_tiles(3, 3, c(1, 2, 4, 5), "3×3 grid", "One random\ndraw"),
  design_tiles(3, 4, c(1, 2, 5, 6, 9, 10), "3×4 grid", "One random\ndraw")) %>%
  mutate(arm = ifelse(x == 1, "Intervention", "Control"))
p6 <- ggplot(d6, aes(c, -r, fill = arm)) + geom_tile(colour = "white", linewidth = 2) +
  scale_fill_manual(NULL, values = c(Intervention = navy, Control = light)) +
  facet_grid(design ~ grid, switch = "y") + coord_equal() + theme_void(base_size = 30) +
  theme(legend.position = "bottom", strip.text.x = element_text(face = "bold", margin = margin(b = 6)),
        strip.text.y.left = element_text(face = "bold", angle = 0, hjust = 1, margin = margin(r = 10)))
save_png(p6, "fig_designs.png", 12, 7.4)

# Trimmed copies of the manuscript's maps and sampled-points figure ----
# The originals carry wide white margins, which shrink them on a slide. Crop to the
# non-white content plus a small pad; the image content itself is unchanged.
#' Crop white margins from a PNG.
#' @param src,dst File names in out_dir. @param pad Pixels kept around the content.
trim_png <- function(src, dst, pad = 8) {
  img <- png::readPNG(file.path(out_dir, src))
  rgb <- img[, , 1:3, drop = FALSE]
  ink <- apply(rgb, c(1, 2), function(v) any(v < 0.97))
  if (dim(img)[3] == 4) ink <- ink & img[, , 4] > 0.05
  r <- range(which(rowSums(ink) > 0)); c <- range(which(colSums(ink) > 0))
  r <- c(max(1, r[1] - pad), min(nrow(ink), r[2] + pad)); c <- c(max(1, c[1] - pad), min(ncol(ink), c[2] + pad))
  png::writePNG(img[r[1]:r[2], c[1]:c[2], , drop = FALSE], file.path(out_dir, dst))
}
trim_png("NC_DCS_Map_NoCounties.png", "NC_DCS_Map_trim.png")
trim_png("NC_MCO_Map.png", "NC_MCO_Map_trim.png")
trim_png("SampledPoints.png", "SampledPoints_trim.png")

cat("Wrote 7 figures, 3 trimmed copies and bss_rank_among_srs.csv to", out_dir, "\n")
