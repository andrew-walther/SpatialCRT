# ============================================================
# Script: nc_design_examples.R
# Purpose: Draw one example treatment allocation per design on the NC
#          58-college map (six designs + SRS, 2018), title-free, for the talk deck.
# Author: Andrew Walther
# Created: 2026-10-28
# Dependencies: ggplot2, sf (via the sourced application scripts)
# ============================================================
#
# Run from the project root:  Rscript presentation/tools/nc_design_examples.R
#
# Same allocations as the manuscript example map (application/code/manuscript_maps.R,
# rs_maps() in render_real_sud_companion.R: 2018, queen graph, first draw of each design), but restricted to the six
# main-text designs and SRS, labeled by name only (no design numbers) and with no
# in-figure title, per the slide figure conventions. Nothing is simulated beyond the
# allocation draw itself, and no application computation file is modified.
# Restricted data: the map shows only 58-cluster treatment indicators.

real_sud_define_only <- TRUE              # define functions, run no profile
source("application/code/run_real_sud.R")   # sources real_sud_simulation.R: rs_config(), rs_units(), rs_draws()

ROOT <- "application/results/real_sud_rev_20261002"
OUT  <- "presentation/figures/nc_design_examples.png"

s <- readRDS(file.path(ROOT, "setup.rds"))
g <- sf::st_simplify(s$clusters, dTolerance = 250)     # display only, as in rs_maps()
units <- rs_units(rs_config("smoke"))

# Main-text designs in the deck's order; IDs follow the production numbering
shown <- c(`8` = "Incidence-Guided Saturation", `6` = "Balanced Quartiles",
           `4` = "Isolation Buffer", `2` = "High Incidence Focus",
           `5` = "2×2 Blocking", `1` = "Checkerboard", `9` = "Simple Random Sampling")

examples <- do.call(rbind, lapply(names(shown), function(id) {
  u <- units[which(units$Design_ID == as.integer(id))[1], ]
  stopifnot(nrow(u) == 1)
  z <- rs_draws(u, s, 1)[, 1]                            # first draw, as in rs_maps()
  a <- g
  a$Design <- sprintf("%s (%d of 58 treated)", shown[[id]], sum(z))
  a$Treatment <- factor(z, levels = 0:1, labels = c("Control", "Treated"))
  a
}))
examples$Design <- factor(examples$Design, levels = unique(examples$Design))

p <- ggplot2::ggplot(examples) +
  ggplot2::geom_sf(ggplot2::aes(fill = Treatment), color = "white", linewidth = 0.12) +
  ggplot2::scale_fill_manual(values = c(Control = "#e0e7ef", Treated = "#1f5fa6"), name = NULL) +
  ggplot2::facet_wrap(~Design, ncol = 4) +
  ggplot2::theme_void(base_size = 13) +
  ggplot2::theme(legend.position = c(0.88, 0.22), legend.text = ggplot2::element_text(size = 14),
                 strip.text = ggplot2::element_text(size = 12, margin = ggplot2::margin(b = 4)))

ggplot2::ggsave(OUT, p, width = 16, height = 5.0, dpi = 150, bg = "white")
stopifnot(file.exists(OUT))
cat("Wrote", OUT, "\n")
