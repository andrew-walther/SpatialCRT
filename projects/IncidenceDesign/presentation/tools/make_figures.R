# ============================================================
# Script: make_figures.R
# Purpose: Draw the talk deck's two diagrams as PNG files (spillover regimes and
#          SUD case identification), replacing the inline SVG/HTML versions.
# Author: Andrew Walther
# Created: 2026-11-02
# Dependencies: ggplot2
# ============================================================
#
# Run from the project root:  Rscript presentation/tools/make_figures.R
#
# Each figure is a file built from code and data so it can be reused in the
# manuscripts or the prelim deck. Legends sit inside the figure; slide captions
# stay as slide text. Colors follow the deck theme (theme/unc-gillings.css).
#
# Inputs:
#   application/results/real_sud_rev_20261002/exhibits/observed_cluster_inputs.csv
#     (tracked 58-cluster aggregate; its 2018-2021 deaths sum to the case count)
#   presentation/data/case_identification_counts.csv
#     (the two upstream counts cited from Habib 2026; no repository data holds them)
# Outputs:
#   presentation/figures/fig_spillover_regimes.png
#   presentation/figures/fig_case_identification.png

library(ggplot2)

OUT <- "presentation/figures"
stopifnot(dir.exists(OUT))

# Deck theme colors
NAVY     <- "#13294B"
CAROLINA <- "#4B9CD3"
TINT     <- "#EBF4FA"
TREATED  <- "#1f5fa6"
CONTROL  <- "#e0e7ef"
GAIN     <- "#FBAF3D"

# Spillover regimes ----

#' Share of queen neighbors that are treated, for every cell of a grid
#'
#' For cell i, share_i = (W Z)_i with W the row-standardized queen contiguity
#' matrix, the quantity the spillover term S(Z) scales by gamma.
#'
#' @param Z integer matrix of 0/1 treatment indicators (rows x columns)
#' @return numeric matrix of the same shape with each cell's treated-neighbor share
#' @examples
#' treated_neighbor_share(matrix(c(1, 0, 0, 1), 2))
treated_neighbor_share <- function(Z) {
  nr <- nrow(Z); nc <- ncol(Z)
  share <- matrix(0, nr, nc)
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    rr <- max(1, r - 1):min(nr, r + 1)
    cc <- max(1, c - 1):min(nc, c + 1)
    nb <- Z[rr, cc, drop = FALSE]
    share[r, c] <- (sum(nb) - Z[r, c]) / (length(nb) - 1)
  }
  share
}

# One illustrative 3 x 4 allocation (row 1 at the top), the same in both panels
Z <- matrix(c(1, 0, 1, 0,
              0, 1, 0, 0,
              0, 0, 0, 1), nrow = 3, byrow = TRUE)
share <- treated_neighbor_share(Z)

cells <- expand.grid(row = 1:3, col = 1:4)
cells$treated <- Z[cbind(cells$row, cells$col)] == 1
cells$share   <- share[cbind(cells$row, cells$col)]

# Who gains: S(Z) = gamma * W Z (both arms) or gamma * (1 - Z) o W Z (controls only)
regimes <- rbind(
  transform(cells, regime = "Spillover to controls only\n(untreated colleges gain)",
            gain = ifelse(treated, 0, share)),
  transform(cells, regime = "Spillover to both arms\n(treated and untreated colleges gain)",
            gain = share)
)
regimes$regime <- factor(regimes$regime, levels = unique(regimes$regime))
regimes$arm <- factor(ifelse(regimes$treated, "Treated college", "Control college"),
                      levels = c("Treated college", "Control college"))
gainers <- regimes[regimes$gain > 0, ]

p_spill <- ggplot(regimes, aes(x = col, y = -row)) +
  geom_tile(aes(fill = arm), colour = "#8a96a3", linewidth = 0.8, width = 0.92, height = 0.92) +
  geom_point(data = gainers, aes(size = gain, colour = "Gains from treated neighbors"),
             shape = 21, fill = GAIN, stroke = 0.9) +
  geom_text(data = gainers, label = "+", colour = NAVY, fontface = "bold", size = 7) +
  facet_wrap(~ regime) +
  scale_fill_manual(values = c("Treated college" = TREATED, "Control college" = CONTROL), name = NULL) +
  scale_colour_manual(values = c("Gains from treated neighbors" = NAVY), name = NULL) +
  scale_size_continuous(range = c(9, 15), guide = "none") +
  coord_equal() +
  labs(caption = "Circle size grows with the share of queen neighbors that are treated") +
  theme_void(base_size = 20) +
  theme(strip.text = element_text(face = "bold", colour = NAVY, lineheight = 1.1,
                                  margin = margin(b = 8)),
        legend.position = "bottom", legend.text = element_text(size = 17),
        plot.caption = element_text(hjust = 0.5, colour = "#505962", size = 15),
        panel.spacing = unit(3, "lines")) +
  guides(fill = guide_legend(order = 1, override.aes = list(colour = "#8a96a3")),
         colour = guide_legend(order = 2, override.aes = list(size = 8, fill = GAIN)))

ggsave(file.path(OUT, "fig_spillover_regimes.png"), p_spill, width = 13, height = 6.4, dpi = 150, bg = "white")
cat("Wrote", file.path(OUT, "fig_spillover_regimes.png"), "\n")

# Case identification ----

upstream <- read.csv("presentation/data/case_identification_counts.csv")
n_all <- upstream$count[upstream$step == "all_deaths_2018_2021"]
n_age <- upstream$count[upstream$step == "deaths_age_18_64"]
stopifnot(length(n_all) == 1, length(n_age) == 1)

clusters <- read.csv("application/results/real_sud_rev_20261002/exhibits/observed_cluster_inputs.csv")
yearly <- clusters[clusters$period %in% as.character(2018:2021), ]
stopifnot(nrow(yearly) == 58 * 4)            # 58 service areas x 4 years, fail loud otherwise
n_sud <- sum(yearly$deaths)

fmt <- function(n) format(n, big.mark = ",")
steps <- data.frame(
  y = 5:1,
  main = c(paste0(fmt(n_all), "\nNorth Carolina deaths, 2018–2021"),
           paste0(fmt(n_age), "\ndeaths among adults aged 18–64"),
           "Out-of-hospital deaths",
           "Natural out-of-hospital deaths",
           paste0(fmt(n_sud), "\npresumed sudden unexpected deaths")),
  excl = c("Exclude ages under 18\nor over 64",
           "Exclude deaths in a hospital,\nhospice or nursing home",
           "Exclude suicide, homicide\nand accidents",
           "Exclude expected deaths\nfrom cancer",
           NA),
  final = c(FALSE, FALSE, FALSE, FALSE, TRUE)
)

MAIN_X <- 0; EXCL_X <- 1.58; W_MAIN <- 1.4; W_EXCL <- 1.04; H <- 0.7
p_case <- ggplot(steps) +
  # downward arrows between main boxes
  geom_segment(data = steps[-5, ], aes(x = MAIN_X, xend = MAIN_X, y = y - H / 2, yend = y - 1 + H / 2),
               colour = CAROLINA, linewidth = 1.2, arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  # side arrows from the flow line to the exclusion boxes
  geom_segment(data = steps[!is.na(steps$excl), ],
               aes(x = MAIN_X, xend = EXCL_X - W_EXCL / 2, y = y - 0.5, yend = y - 0.5),
               colour = "#505962", linewidth = 0.9, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_tile(aes(x = MAIN_X, y = y, fill = final), width = W_MAIN, height = H,
            colour = ifelse(steps$final, NAVY, CAROLINA), linewidth = 1) +
  geom_text(aes(x = MAIN_X, y = y, label = main, colour = final), size = 7, lineheight = 0.95) +
  geom_tile(data = steps[!is.na(steps$excl), ], aes(x = EXCL_X, y = y - 0.5),
            width = W_EXCL, height = H * 0.8, fill = "#f6f6f6", colour = "#9aa3ad",
            linetype = "dashed", linewidth = 0.7) +
  geom_text(data = steps[!is.na(steps$excl), ], aes(x = EXCL_X, y = y - 0.5, label = excl),
            size = 6.6, lineheight = 0.95, colour = "#222222") +
  scale_fill_manual(values = c("FALSE" = TINT, "TRUE" = NAVY), guide = "none") +
  scale_colour_manual(values = c("FALSE" = NAVY, "TRUE" = "white"), guide = "none") +
  coord_cartesian(xlim = c(-0.72, 2.12), ylim = c(0.6, 5.4), expand = FALSE) +
  theme_void()

ggsave(file.path(OUT, "fig_case_identification.png"), p_case, width = 10, height = 7.6, dpi = 150, bg = "white")
cat("Wrote", file.path(OUT, "fig_case_identification.png"), "(", fmt(n_sud), "cases )\n")
