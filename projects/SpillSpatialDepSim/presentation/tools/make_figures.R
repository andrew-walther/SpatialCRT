# ============================================================
# Script: presentation/tools/make_figures.R
# Purpose: Build the slide-sized figures for the Project 1 talk from the saved
#          per-allocation beta-hat MSE results (no hand-typed numbers).
# Author: Andrew Walther
# Created: 2026-10-25 (slide-sized versions 2026-10-28)
# Dependencies: ggplot2, dplyr, ragg, sf, tigris (NC county shapes; downloaded once, then cached)
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

# Figures 7a/7b: every 3x3 and 3x4 allocation (backup slides) ----
# Slide-sized redraw of the chapter's box-plot figure (review/code/build_tables_figures.R,
# plot_grid()) from the same corrected data: "Random" = all balanced allocations,
# "Block" = the BSS allocations; rows = spillover type, columns = rho.
mse_all_alloc <- function(gr) {
  d <- bind_rows(mutate(filter(all_alloc, grid == gr), Sampling = "Random"),
                 mutate(filter(all_alloc, grid == gr, is_BSS), Sampling = "Block")) %>%
    mutate(Sampling = factor(Sampling, levels = c("Block", "Random")),
           spill = factor(spill, levels = spill_labels),
           rho = factor(sprintf("ρ = %s", ifelse(Rho == 0, "0", "0.01")), levels = c("ρ = 0", "ρ = 0.01")))
  ggplot(d, aes(x = Sampling, y = MSE, fill = factor(Psi))) +
    geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.75)) +
    scale_fill_manual("Spillover ψ", values = c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B")) +
    facet_grid(spill ~ rho, scales = "free_y") +
    labs(x = "Sampling group", y = "MSE of the intervention effect") +
    theme_slide(20) + theme(strip.text.y = element_text(angle = 0, hjust = 0))
}
save_png(mse_all_alloc("3x3"), "fig_mse_3x3_all.png", 15, 8)
save_png(mse_all_alloc("3x4"), "fig_mse_3x4_all.png", 15, 8)

# Figure 8: simulated subjects on the 2x4 grid ----
# Same scheme as SpatialSim_NC_DOC.Rmd: 160 subjects per dataset at integer
# coordinates drawn uniformly over the 400 x 200 grid of 100 x 100 districts;
# colour = the district each subject falls in.
set.seed(8)
pts <- data.frame(x = sample(0:399, 160, replace = TRUE), y = sample(0:199, 160, replace = TRUE)) %>%
  mutate(district = factor((y %/% 100) * 4 + x %/% 100 + 1))
cells8 <- expand.grid(cx = seq(50, 350, 100), cy = c(50, 150))
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
p8 <- ggplot() +
  geom_tile(data = cells8, aes(cx, cy), width = 100, height = 100, fill = "white", colour = "#1d1f21", linewidth = 1) +
  geom_point(data = pts, aes(x, y, colour = district), size = 3.6) +
  scale_colour_manual(values = okabe_ito, guide = "none") + coord_equal() + theme_void()
save_png(p8, "fig_sampled_points.png", 12, 6)

# Schematic figures (formerly inline SVG in the deck) ----
# Generic rook grids; no data. Palette as above, amber = spillover.
amber <- "#FBAF3D"
#' Tile data for an nr x nc grid, cells numbered row by row from the top left.
#' @param nr,nc Rows and columns. @return data.frame with r, c, id, and plot coords px, py.
grid_cells <- function(nr, nc) {
  g <- expand.grid(c = 1:nc, r = 1:nr)
  transform(g, id = (r - 1) * nc + c, px = c, py = -r)
}
#' Binary rook exposure, as in the simulation: a cluster is exposed if a rook
#' neighbour is treated; under "Control Clusters Only" only controls can be exposed.
#' @param g Output of grid_cells() with a 0/1 column x. @param controls_only Logical.
#' @return Integer vector z.
exposure <- function(g, controls_only) {
  vapply(seq_len(nrow(g)), function(i) {
    nb <- abs(g$r - g$r[i]) + abs(g$c - g$c[i]) == 1
    as.integer(any(g$x[nb] == 1) && (!controls_only || g$x[i] == 0))
  }, integer(1))
}
theme_grid <- function(base = 24) {
  theme_void(base_size = base) +
    theme(legend.position = "bottom", legend.key.size = unit(1.6, "lines"),
          strip.text = element_text(face = "bold", margin = margin(b = 8)))
}
arm_fill <- c("Intervention" = navy, "Control" = light)

# 9: a region of nine areas, four treated, arrows from treated into neighbouring controls
g9 <- grid_cells(3, 3) %>% mutate(arm = ifelse(id %in% c(1, 3, 5, 7), "Intervention area", "Control area"))
arrows9 <- data.frame(x = c(1, 2, 2, 1), y = c(-1, -2, -2, -3), dx = c(1, 0, 1, 0), dy = c(0, 1, 0, 1)) %>%
  mutate(x0 = x + 0.25 * dx, y0 = y + 0.25 * dy, x1 = x + 0.8 * dx, y1 = y + 0.8 * dy)
p9 <- ggplot() + geom_tile(data = g9, aes(px, py, fill = arm), colour = "white", linewidth = 3) +
  geom_segment(data = arrows9, aes(x = x0, y = y0, xend = x1, yend = y1, colour = "Possible spillover"),
               linewidth = 2.6, arrow = arrow(type = "closed", length = unit(0.22, "in"))) +
  scale_fill_manual(NULL, values = c("Intervention area" = navy, "Control area" = light),
                    breaks = c("Intervention area", "Control area")) +
  scale_colour_manual(NULL, values = c("Possible spillover" = amber)) +
  coord_equal() + theme_grid(36) + theme(legend.box = "vertical", legend.spacing.y = unit(0.2, "lines")) +
  guides(fill = guide_legend(override.aes = list(colour = NA)))
save_png(p9, "fig_places.png", 8, 9)

# 10a: spatial dependence, outcomes graded smoothly across neighbours, no treatment
dep_shades <- c("#1f4e79", "#2f6aa0", "#5b8fc0", "#8fb5da", "#c3d8ec")
g10a <- grid_cells(3, 3) %>% mutate(shade = dep_shades[r + c - 1])
p10a <- ggplot(g10a, aes(px, py, fill = shade)) + geom_tile(colour = "white", linewidth = 3) +
  scale_fill_identity() + coord_equal() + theme_void()
save_png(p10a, "fig_dependence.png", 6, 6)

# 10b: spillover, a treated centre raises its four rook neighbours
g10b <- grid_cells(3, 3) %>% mutate(shade = case_when(id == 5 ~ navy, id %in% c(2, 4, 6, 8) ~ amber, TRUE ~ light))
arrows10 <- data.frame(dx = c(0, 0, -1, 1), dy = c(1, -1, 0, 0)) %>%
  mutate(x0 = 2 + 0.3 * dx, y0 = -2 + 0.3 * dy, x1 = 2 + 0.78 * dx, y1 = -2 + 0.78 * dy)
p10b <- ggplot() + geom_tile(data = g10b, aes(px, py, fill = shade), colour = "white", linewidth = 3) +
  geom_segment(data = arrows10, aes(x = x0, y = y0, xend = x1, yend = y1), colour = navy,
               linewidth = 2.2, arrow = arrow(type = "closed", length = unit(0.2, "in"))) +
  scale_fill_identity() + coord_equal() + theme_void()
save_png(p10b, "fig_spillover_center.png", 6, 6)

# 11: the two spillover types on a 2x4 grid with the left half treated
g11 <- grid_cells(2, 4) %>% mutate(x = as.integer(c <= 2))
d11 <- bind_rows(mutate(g11, z = exposure(g11, TRUE), type = spill_labels[["TrtNoSpill"]]),
                 mutate(g11, z = exposure(g11, FALSE), type = spill_labels[["TrtSpill"]])) %>%
  mutate(type = factor(type, levels = spill_labels),
         arm = ifelse(x == 1, "Intervention cluster (x = 1)", "Control cluster (x = 0)"))
p11 <- ggplot(d11, aes(px, py)) +
  geom_tile(aes(fill = arm), colour = "white", linewidth = 3) +
  geom_point(data = filter(d11, z == 1), aes(colour = "Exposed to spillover (z = 1)"),
             shape = 21, fill = amber, stroke = 1.2, size = 22) +
  geom_text(data = filter(d11, z == 1), label = "z = 1", size = 5.2, fontface = "bold", colour = navy) +
  scale_fill_manual(NULL, values = c("Intervention cluster (x = 1)" = navy, "Control cluster (x = 0)" = light),
                    breaks = c("Intervention cluster (x = 1)", "Control cluster (x = 0)")) +
  scale_colour_manual(NULL, values = c("Exposed to spillover (z = 1)" = navy),
                      guide = guide_legend(override.aes = list(size = 7))) +
  facet_wrap(~ type, nrow = 1) + coord_equal() + theme_grid(32) + theme(panel.spacing.x = unit(3, "lines")) +
  guides(fill = guide_legend(override.aes = list(colour = NA)))
save_png(p11, "fig_spill_types.png", 14, 5.4)

# 12: row-normalised rook weights, cluster 2 and its neighbours on a 2x4 grid
g12 <- grid_cells(2, 4) %>%
  mutate(role = case_when(id == 2 ~ "Cluster 2", id %in% c(1, 3, 6) ~ "Rook neighbour of cluster 2",
                          TRUE ~ "Other cluster"),
         role = factor(role, levels = c("Cluster 2", "Rook neighbour of cluster 2", "Other cluster")))
p12 <- ggplot(g12, aes(px, py)) + geom_tile(aes(fill = role), colour = "white", linewidth = 3) +
  geom_text(aes(label = id, colour = role == "Cluster 2"), size = 12, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(NULL, values = c("Cluster 2" = navy, "Rook neighbour of cluster 2" = amber, "Other cluster" = light)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = navy)) +
  coord_equal() + theme_grid(28) + guides(fill = guide_legend(nrow = 2, override.aes = list(colour = NA)))
save_png(p12, "fig_weights_2x4.png", 9, 6)

# 13: buffer-zone design, one treated cluster, a ring of excluded clusters, unexposed controls
g13 <- grid_cells(4, 5) %>%
  mutate(role = case_when(r == 3 & c == 3 ~ "Intervention",
                          r >= 2 & c >= 2 & c <= 4 ~ "Excluded buffer",
                          TRUE ~ "Unexposed control"),
         role = factor(role, levels = c("Intervention", "Excluded buffer", "Unexposed control")))
p13 <- ggplot(g13, aes(px, py, fill = role, colour = role, linetype = role)) +
  geom_tile(width = 0.92, height = 0.92, linewidth = 1.1) +
  scale_fill_manual(NULL, values = c("Intervention" = navy, "Excluded buffer" = "white", "Unexposed control" = light)) +
  scale_colour_manual(NULL, values = c("Intervention" = navy, "Excluded buffer" = "#9aa3ad", "Unexposed control" = light)) +
  scale_linetype_manual(NULL, values = c("Intervention" = "solid", "Excluded buffer" = "dashed", "Unexposed control" = "solid")) +
  coord_equal() + theme_grid(34)
save_png(p13, "fig_buffer_zone.png", 9, 8)

# North Carolina maps ----
# County shapes: US Census cartographic boundaries (tigris, cached), 2023 vintage.
# County -> Division of Community Supervision district: data/NC_County_Assignments.csv,
# the county file behind the paper's maps (read off the DAC district map). County -> LME/MCO:
# data/NC_LME_MCO_Counties_2025.csv, from the NCDHHS LME/MCO directory as archived
# 2025-05-31 (four LME/MCOs, before the Vaya-Partners merger). See data/README.md.
suppressMessages({library(sf); library(tigris)})
options(tigris_use_cache = TRUE)
dcs <- read.csv("data/NC_County_Assignments.csv", colClasses = "character") %>%
  transmute(county = trimws(County), district = trimws(District), division = trimws(Division))
nc <- counties(state = "NC", cb = TRUE, year = 2023, progress_bar = FALSE) %>%
  select(county = NAME, geometry) %>%
  left_join(dcs, by = "county") %>%
  left_join(read.csv("data/NC_LME_MCO_Counties_2025.csv"), by = "county")
stopifnot(nrow(nc) == 100, !anyNA(nc$district), !anyNA(nc$lme_mco))
districts <- nc %>% group_by(division, district) %>% summarise(geometry = st_union(geometry), .groups = "drop")
divisions <- nc %>% group_by(division) %>% summarise(geometry = st_union(geometry), .groups = "drop")
lab_pt <- function(x) cbind(st_drop_geometry(x), st_coordinates(suppressWarnings(st_point_on_surface(x))))

# 14: districts coloured by judicial division (redraw of the manuscript's DCS map)
p14 <- ggplot() + geom_sf(data = divisions, aes(fill = factor(division)), colour = NA) +
  geom_sf(data = districts, fill = NA, colour = "black", linewidth = 0.6) +
  geom_text(data = lab_pt(districts), aes(X, Y, label = district), size = 6, fontface = "bold") +
  scale_fill_manual("Judicial division", values = c("1" = "#ffeb3b", "2" = "#e53935", "3" = "#43a047", "4" = "#1e88e5")) +
  theme_void(base_size = 22) + theme(legend.position = "bottom")
save_png(p14, "fig_nc_dcs_districts.png", 14, 6.2)

# 15: counties coloured by LME/MCO (redraw of the manuscript's MCO map)
p15 <- ggplot(nc) + geom_sf(aes(fill = lme_mco), colour = "#6b7178", linewidth = 0.25) +
  geom_text(data = lab_pt(nc), aes(X, Y, label = county), size = 2.6, colour = "#1d1f21") +
  scale_fill_manual("LME/MCO", values = c("Alliance Health" = "#D9A8D2", "Partners Health Management" = "#7FA6CE",
                                          "Trillium Health Resources" = "#BFE0B1", "Vaya Health" = "#C9E3F3")) +
  theme_void(base_size = 22) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
save_png(p15, "fig_nc_mco.png", 14, 6.6)

# 16: the example best-effort BSS allocation of the 31 districts (manuscript figure).
# Treated districts as hard-coded in code/NCTrtMaps.Rmd (treatment_districts_list).
treated16 <- c("1", "3", "4", "7", "10", "12", "13", "15", "19A", "19B", "21", "23", "24", "27", "30")
d16 <- districts %>% mutate(arm = ifelse(district %in% treated16, "Intervention", "Control"))
div_lab <- lab_pt(divisions) %>% mutate(label = paste("Division", division))
p16 <- ggplot() + geom_sf(data = d16, aes(fill = arm), colour = "gray60", linewidth = 0.3) +
  geom_sf(data = divisions, fill = NA, colour = "gray20", linewidth = 1.4) +
  geom_text(data = lab_pt(d16), aes(X, Y, label = district), colour = "white", size = 5.5, fontface = "bold") +
  geom_label(data = div_lab, aes(X, Y, label = label), nudge_y = 0.32, size = 6, fontface = "bold.italic",
             fill = "white", colour = "gray20", linewidth = 0) +
  scale_fill_manual("Treatment assignment", values = c("Intervention" = "red3", "Control" = "dodgerblue3"),
                    breaks = c("Intervention", "Control")) +
  theme_void(base_size = 22) + theme(legend.position = "bottom")
save_png(p16, "fig_nc_bss_assignment.png", 14, 6.4)

cat("Wrote 19 figures and bss_rank_among_srs.csv to", out_dir, "\n")
