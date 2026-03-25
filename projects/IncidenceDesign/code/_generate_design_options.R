# Temporary script to generate Option 1 and Option 2 design figures
# Run from projects/IncidenceDesign/

source("code/06_visualizations.R")

set.seed(2026)
grid <- build_spatial_grid(10)
coords     <- grid$coords
N_clusters <- grid$N_clusters
nb_queen   <- grid$nb_queen
dummy_incidence <- runif(N_clusters, 0, 1)

design_labels <- c(
  "1" = "1: Checkerboard",
  "2" = "2: High Incidence Focus",
  "3" = "3: Saturation Quadrants",
  "4" = "4: Isolation Buffer",
  "5" = "5: 2\u00d72 Blocking",
  "6" = "6: Balanced Quartiles",
  "7" = "7: Balanced Halves",
  "8" = "8: Guided Sat. Quadrants"
)

# Generate assignments (same seed as main figure)
# We need to track the random state for Design 3 to recover its saturation mapping
sample_designs_list <- list()
d3_sats <- NULL  # will capture Design 3's random saturations

for (d in 1:8) {
  if (d == 3) {
    # Capture the RNG state before Design 3 so we can recover its sats
    rng_before <- .Random.seed
    z <- get_designs(d, 1, N_clusters, dummy_incidence, nb_queen, coords)[, 1]
    # Re-derive what saturation was assigned: restore seed and call sample()
    .Random.seed <<- rng_before
    d3_sats <- sample(c(0.20, 0.40, 0.60, 0.80))
    # Restore the RNG state to after the design call
    # (just re-run to advance the seed properly)
    .Random.seed <<- rng_before
    z <- get_designs(d, 1, N_clusters, dummy_incidence, nb_queen, coords)[, 1]
  } else {
    z <- get_designs(d, 1, N_clusters, dummy_incidence, nb_queen, coords)[, 1]
  }

  sample_designs_list[[d]] <- data.frame(
    x = coords$x, y = coords$y,
    Incidence = dummy_incidence,
    Assignment = as.factor(z),
    Design = factor(design_labels[as.character(d)], levels = design_labels)
  )
}
sample_designs_df <- do.call(rbind, sample_designs_list)

# --- Compute strata metadata ---
med_val   <- median(dummy_incidence)
quartiles <- dplyr::ntile(dummy_incidence, 4)

# Spatial quadrant IDs
half <- 5
q_id <- ifelse(coords$x <= half & coords$y <= half, 1,
        ifelse(coords$x > half  & coords$y <= half, 2,
        ifelse(coords$x <= half & coords$y > half,  3, 4)))

# Design 8 saturation mapping
quad_means <- sapply(1:4, function(q) mean(dummy_incidence[q_id == q]))
rank_order <- order(quad_means, decreasing = TRUE)
sat_levels <- c(0.80, 0.60, 0.40, 0.20)
quad_sats  <- numeric(4)
quad_sats[rank_order] <- sat_levels

# =====================================================================
# OPTION 1: Same 8-panel figure with subtle boundary overlays
# =====================================================================

# Helper to make quadrant boundary segments for a specific facet
make_quad_lines <- function(design_label) {
  rbind(
    data.frame(x = 0.5, xend = 10.5, y = 5.5, yend = 5.5,
               Design = factor(design_label, levels = design_labels)),
    data.frame(x = 5.5, xend = 5.5, y = 0.5, yend = 10.5,
               Design = factor(design_label, levels = design_labels))
  )
}

make_block_lines <- function(design_label) {
  h <- data.frame(x = 0.5, xend = 10.5,
                  y = seq(2.5, 8.5, by = 2), yend = seq(2.5, 8.5, by = 2),
                  Design = factor(design_label, levels = design_labels))
  v <- data.frame(x = seq(2.5, 8.5, by = 2), xend = seq(2.5, 8.5, by = 2),
                  y = 0.5, yend = 10.5,
                  Design = factor(design_label, levels = design_labels))
  rbind(h, v)
}

# Quadrant lines for designs 3 and 8
quad_line_data <- rbind(
  make_quad_lines(design_labels["3"]),
  make_quad_lines(design_labels["8"])
)

# Block lines for design 5
block_line_data <- make_block_lines(design_labels["5"])

# H/L labels per cell (Design 2)
d2_stratum <- data.frame(
  x = coords$x, y = coords$y,
  stratum = ifelse(dummy_incidence > med_val, "H", "L"),
  Design = factor(design_labels["2"], levels = design_labels)
)

# Saturation labels for Design 3 (quadrant centers)
# d3_sats[q] is the saturation for quadrant q (1=SW, 2=SE, 3=NW, 4=NE)
d3_sat_label_data <- data.frame(
  x = c(3, 8, 3, 8), y = c(3, 3, 8, 8),
  label = paste0(d3_sats * 100, "%"),
  Design = factor(design_labels["3"], levels = design_labels)
)

# Saturation labels for Design 8 — include mean incidence to show link
d8_sat_label_data <- data.frame(
  x = c(3, 8, 3, 8), y = c(3, 3, 8, 8),
  label = paste0(quad_sats * 100, "%"),
  Design = factor(design_labels["8"], levels = design_labels)
)

# Mean incidence per quadrant (smaller, below the sat label) for Design 8
d8_mean_label_data <- data.frame(
  x = c(3, 8, 3, 8), y = c(3, 3, 8, 8),
  label = paste0("\u0078\u0304=", sprintf("%.2f", quad_means)),
  Design = factor(design_labels["8"], levels = design_labels)
)

# Quartile number per cell (Design 6)
d6_stratum <- data.frame(
  x = coords$x, y = coords$y, stratum = quartiles,
  Design = factor(design_labels["6"], levels = design_labels)
)

# Half label per cell (Design 7)
d7_stratum <- data.frame(
  x = coords$x, y = coords$y,
  stratum = ifelse(dummy_incidence > med_val, "H", "L"),
  Design = factor(design_labels["7"], levels = design_labels)
)

p1 <- ggplot(sample_designs_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = Incidence), color = "grey50", linewidth = 0.4) +
  geom_point(aes(shape = Assignment), fill = "white", color = "black",
             size = 2.2, stroke = 0.9) +
  # Quadrant boundaries (Designs 3, 8)
  geom_segment(data = quad_line_data,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "red", linewidth = 0.8, linetype = "dashed",
               inherit.aes = FALSE) +
  # Block boundaries (Design 5)
  geom_segment(data = block_line_data,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "red", linewidth = 0.6, linetype = "dashed",
               inherit.aes = FALSE) +
  # H/L labels in cell corner (Design 2)
  geom_text(data = d2_stratum, aes(x = x, y = y, label = stratum),
            size = 2.0, color = "red", fontface = "bold",
            nudge_x = 0.3, nudge_y = 0.3, inherit.aes = FALSE) +
  # Saturation % labels (Design 3)
  geom_label(data = d3_sat_label_data, aes(x = x, y = y, label = label),
             size = 3.5, fontface = "bold", fill = "white", alpha = 0.75,
             linewidth = 0, inherit.aes = FALSE) +
  # Saturation % labels (Design 8) — top line
  geom_label(data = d8_sat_label_data, aes(x = x, y = y, label = label),
             size = 3.2, fontface = "bold", fill = "white", alpha = 0.8,
             linewidth = 0, nudge_y = 0.8, inherit.aes = FALSE) +
  # Mean incidence per quadrant (Design 8) — bottom line
  geom_label(data = d8_mean_label_data, aes(x = x, y = y, label = label),
             size = 2.3, fill = "white", alpha = 0.8,
             linewidth = 0, nudge_y = -0.8, inherit.aes = FALSE) +
  # Quartile number in cell corner (Design 6)
  geom_text(data = d6_stratum, aes(x = x, y = y, label = stratum),
            size = 1.8, color = "red", fontface = "bold",
            nudge_x = 0.3, nudge_y = 0.3, inherit.aes = FALSE) +
  # Half label in cell corner (Design 7)
  geom_text(data = d7_stratum, aes(x = x, y = y, label = stratum),
            size = 1.8, color = "red", fontface = "bold",
            nudge_x = 0.3, nudge_y = 0.3, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Baseline\nIncidence") +
  scale_shape_manual(
    values = c("0" = 4, "1" = 21),
    labels = c("0" = "Control (X)", "1" = "Treated (\u25cf)"),
    name = "Treatment\nAssignment"
  ) +
  facet_wrap(~ Design, ncol = 4) +
  theme_void(base_size = 12) +
  theme(
    legend.position  = "right",
    strip.text       = element_text(face = "bold", size = 10, margin = margin(b = 6, t = 6)),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

ggsave("results/figures/design_samples_option1_overlays.png", p1,
       width = 14, height = 7, dpi = 300, bg = "white")
ggsave("results/figures/design_samples_option1_overlays.pdf", p1,
       width = 14, height = 7, bg = "white")
cat("Option 1 saved (PNG + PDF)\n")
