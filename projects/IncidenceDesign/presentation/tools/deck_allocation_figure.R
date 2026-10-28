# ============================================================
# Script: deck_allocation_figure.R
# Purpose: Slide-sized version of the allocation-risk figure (CTJ Figure 3) with
#          larger text, drawn from the same stored summary, for the talk deck.
# Author: Andrew Walther
# Created: 2026-10-28
# Dependencies: ggplot2; pdftoppm (poppler) for the PNG
# ============================================================
#
# Run from the project root:  Rscript presentation/tools/deck_allocation_figure.R
#
# Same data and encoding as code/21_downside_by_design.R (thin line 5th-95th
# percentile, thick bar 25th-75th, circle the mean, all relative to SRS's mean in the
# same setting), with three slide changes: base font 17 pt instead of 9, the design
# order and names unchanged, and SRS drawn in dark grey rather than orange so it does
# not read as "worse than SRS" next to the orange cells of the NC accuracy heat map.
# Reads results/allocation_risk/summary/surfaces10/allocation_distribution_by_design.csv
# only; nothing is simulated.

library(ggplot2)

SRC <- "results/allocation_risk/summary/surfaces10/allocation_distribution_by_design.csv"
OUT <- "presentation/figures/deck_allocation_distribution"

dist <- read.csv(SRC)
ids <- c(8, 6, 4, 2, 5, 1, 9)
# "­" prints a true hyphen in pdf(); "-" would print as a minus sign
labels <- c("Incidence­Guided Saturation", "Balanced Quartiles", "Isolation Buffer",
            "High Incidence Focus", "2×2 Blocking", "Checkerboard", "Simple Random Sampling")
f <- dist[dist$Design_ID %in% ids, ]
stopifnot(nrow(f) == 7 * 2 * 2)                     # 7 designs x 2 studies x 2 regimes
f$Label <- factor(f$Design_ID, levels = rev(ids), labels = rev(labels))
f$Study <- factor(f$Source, c("grid_surfaces10", "nc_refined"), c("Simulation grid", "North Carolina"))
f$Regime <- factor(f$Regime, c("both", "control_only"),
                   c("Spillover to both arms", "Spillover to controls only"))
f$Benchmark <- f$Design_ID == 9

fig <- ggplot(f, aes(y = Label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
  geom_linerange(aes(xmin = Q05, xmax = Q95, color = Benchmark), linewidth = 0.9) +
  geom_linerange(aes(xmin = Q25, xmax = Q75, color = Benchmark), linewidth = 5) +
  geom_point(aes(x = Mean), shape = 21, fill = "white", color = "black", size = 3.2, stroke = 0.9) +
  facet_grid(Study ~ Regime) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8")) +
  scale_color_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#4A4F55"), guide = "none") +
  labs(x = "Error relative to the average SRS allocation (log scale)", y = NULL) +
  theme_minimal(base_size = 17) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 17),
        strip.text.y = element_text(angle = -90),
        panel.spacing = unit(1.2, "lines"))

ggsave(paste0(OUT, ".pdf"), fig, width = 15, height = 8.2)
status <- system2("pdftoppm", c("-png", "-r", "150", "-singlefile", shQuote(paste0(OUT, ".pdf")), shQuote(OUT)))
stopifnot(status == 0, file.exists(paste0(OUT, ".png")))
file.remove(paste0(OUT, ".pdf"))
cat("Wrote", paste0(OUT, ".png"), "\n")
