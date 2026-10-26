# ============================================================
# Script: presentation/tools/make_figures.R
# Purpose: Build the new talk figures for Project 1 (allocation spread,
#          where BSS falls among SRS allocations, collinearity diagram)
#          from the per-allocation beta-hat MSE results.
# Author: Andrew Walther
# Created: 2026-10-25
# Dependencies: ggplot2, dplyr
# ============================================================
# Run from the project root (SpillSpatialDepSim/):
#   Rscript presentation/tools/make_figures.R
# Inputs are read-only: results/beta_mse (published 2x4; Intervention & Control
# for every grid except one rerun scenario) and the corrected reruns in
# review/results/. Outputs go to presentation/figures/. No figure titles: the
# slide title and caption carry them.

suppressMessages({library(ggplot2); library(dplyr)})

# Inputs ----
out_dir <- "presentation/figures"
bss_idx <- list(`2x4` = c(28, 43), `3x3` = c(28, 30, 35, 51, 77, 106), `3x4` = c(313, 612))
navy <- "#13294B"; gray <- "#7D8590"

#' Path to the corrected per-allocation beta-hat MSE CSV for one scenario.
#'
#' Same rule as review/code/build_tables_figures.R ("corrected" version):
#' 3x3/3x4 Control Only come from the loop-fixed rerun, the one stray
#' Intervention & Control file from its rerun, everything else from results/.
#'
#' @param spill "TrtNoSpill" (Control Only) or "TrtSpill" (Intervention & Control).
#' @param grid "2x4", "3x3" or "3x4".
#' @param psi,rho Scenario values.
#' @return File path.
csv_path <- function(spill, grid, psi, rho) {
  f <- sprintf("beta_mse_results_%s_%s_A02_B1_P%02d_R%s.csv", spill, grid,
               round(psi * 10), ifelse(rho == 0, "00", "001"))
  dirs <- c(if (spill == "TrtNoSpill" && grid != "2x4") "review/results/beta_mse_TrtNoSpill_fixed",
            if (spill == "TrtSpill") "review/results/beta_mse_TrtSpill_rerun",
            "results/beta_mse")
  hit <- file.path(dirs, f)[file.exists(file.path(dirs, f))]
  stopifnot(length(hit) >= 1)
  hit[1]
}

#' Per-allocation MSE for one grid, every scenario, labelled SRS or BSS.
#'
#' Every allocation is an SRS allocation; the BSS allocations are also flagged.
#' @param grid Grid label.
#' @return Data frame: grid, spill label, Psi, Rho, combo, MSE, is_BSS.
load_grid <- function(grid) {
  out <- list()
  for (spill in c("TrtNoSpill", "TrtSpill")) for (psi in c(0.5, 0.6, 0.7, 0.8)) for (rho in c(0.01, 0)) {
    d <- read.csv(csv_path(spill, grid, psi, rho))
    out[[length(out) + 1]] <- data.frame(
      grid = grid, Psi = psi, Rho = rho, combo = seq_len(nrow(d)), MSE = d$MSE,
      spill = ifelse(spill == "TrtSpill", "Spillover into both arms", "Spillover into controls only"),
      is_BSS = seq_len(nrow(d)) %in% bss_idx[[grid]])
  }
  bind_rows(out)
}

all_grids <- bind_rows(lapply(names(bss_idx), load_grid))

# Share of SRS allocations whose MSE exceeds the BSS mean MSE, per scenario ----
# (SRS = all balanced allocations, as in the paper; BSS mean = mean over the BSS set)
rank_tab <- all_grids %>%
  group_by(grid, spill, Psi, Rho) %>%
  summarise(BSS_mean = mean(MSE[is_BSS]), SRS_mean = mean(MSE),
            pct_SRS_worse = 100 * mean(MSE > mean(MSE[is_BSS])),
            SRS_max = max(MSE), .groups = "drop")
write.csv(rank_tab, file.path(out_dir, "bss_rank_among_srs.csv"), row.names = FALSE)

# Figure 1: every 2x4 allocation, BSS marked (psi = 0.5) ----
# One row per spillover type x rho; x on a square-root scale and free per row,
# because the four settings differ by orders of magnitude.
# Row header = setting plus the share of random allocations that BSS beats
lab1 <- rank_tab %>% filter(grid == "2x4", Psi == 0.5) %>%
  mutate(order = (spill == "Spillover into both arms") * 2 + (Rho == 0),
         panel = sprintf("%s, %s:  BSS beats %.0f%% of random allocations", spill,
                         ifelse(Rho == 0, "no spatial dependence", "spatial dependence (ρ = 0.01)"),
                         pct_SRS_worse)) %>% arrange(order)
d1 <- all_grids %>% filter(grid == "2x4", Psi == 0.5) %>%
  left_join(select(lab1, spill, Rho, panel), by = c("spill", "Rho")) %>%
  mutate(panel = factor(panel, levels = lab1$panel))
set.seed(1)  # vertical jitter only, for legibility
p1 <- ggplot(d1, aes(x = MSE, y = 0)) +
  geom_jitter(data = filter(d1, !is_BSS), aes(shape = "Simple random sampling (each allocation)"),
              width = 0, height = 0.25, colour = gray, alpha = 0.7, size = 2.6) +
  geom_point(data = filter(d1, is_BSS), aes(shape = "Block stratified sampling"),
             colour = navy, fill = navy, size = 4.2) +
  scale_shape_manual(NULL, values = c("Simple random sampling (each allocation)" = 16,
                                      "Block stratified sampling" = 23)) +
  scale_x_sqrt() + scale_y_continuous(limits = c(-0.6, 0.6), breaks = NULL) +
  facet_wrap(~ panel, ncol = 1, scales = "free_x") +
  labs(x = "Mean squared error of the intervention effect (square-root scale)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", strip.text = element_text(hjust = 0, face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "fig_allocation_spread_2x4.png"), p1, width = 9, height = 6.2, dpi = 150, bg = "white")

# Figure 2: where BSS falls among SRS allocations, every grid and scenario ----
d2 <- rank_tab %>%
  mutate(setting = factor(sprintf("%s\n%s", spill, ifelse(Rho == 0, "ρ = 0", "ρ = 0.01")),
                          levels = rev(c("Spillover into controls only\nρ = 0.01", "Spillover into controls only\nρ = 0",
                                         "Spillover into both arms\nρ = 0.01", "Spillover into both arms\nρ = 0"))),
         grid = factor(grid, levels = c("2x4", "3x3", "3x4"), labels = c("2×4", "3×3", "3×4")))
p2 <- ggplot(d2, aes(x = pct_SRS_worse, y = setting)) +
  geom_vline(xintercept = 50, colour = "#b9bec4", linetype = "dashed") +
  geom_point(aes(shape = grid), colour = navy, size = 3,
             position = position_dodge(width = 0.6)) +
  scale_shape_manual("Grid", values = c(16, 17, 15)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Random allocations with higher MSE than block stratified sampling (%)", y = NULL) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(out_dir, "fig_bss_rank_all_grids.png"), p2, width = 9, height = 5, dpi = 150, bg = "white")

# Figure 3: why a checkerboard makes treatment and spillover collinear ----
# Generic 2x4 rook grid (not tied to the simulation's combo numbering).
# Spillover indicator, Control Only: z = 1 for a control cell with a treated rook neighbour.
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
d3 <- bind_rows(mutate(spill_z(checker), allocation = "Checkerboard (block stratified)"),
                mutate(spill_z(blocked), allocation = "Clustered (a random allocation)")) %>%
  mutate(cell = case_when(x == 1 ~ "Intervention",
                          z == 1 ~ "Control, exposed to spillover (z = 1)",
                          TRUE ~ "Control, not exposed (z = 0)"))
p3 <- ggplot(d3, aes(x = c, y = -r, fill = cell)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("x=%d\nz=%d", x, z)), size = 3.6,
            colour = ifelse(d3$cell == "Intervention", "white", "#1d1f21")) +
  scale_fill_manual(NULL, breaks = c("Intervention", "Control, exposed to spillover (z = 1)",
                                     "Control, not exposed (z = 0)"),
                    values = c("Intervention" = navy,
                                     "Control, exposed to spillover (z = 1)" = "#4B9CD3",
                                     "Control, not exposed (z = 0)" = "#E1E5EA")) +
  facet_wrap(~ allocation, ncol = 2) + coord_equal() +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", margin = margin(b = 6)))
ggsave(file.path(out_dir, "fig_collinearity_grid.png"), p3, width = 9, height = 3.6, dpi = 150, bg = "white")

cat("Wrote 3 figures and bss_rank_among_srs.csv to", out_dir, "\n")
print(filter(rank_tab, Psi == 0.5), n = 24)
