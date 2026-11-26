# ============================================================
# Script: thesis_design_examples.R
# Purpose: Thesis-width export of the NC example-allocation map (Chapter 3 Figure
#          fig:nc-examples): 6.5 in wide, panel and legend text >= 8 pt, 300 dpi.
# Author: Andrew Walther
# Created: 2026-11-26
# Dependencies: ggplot2, sf (via render_real_sud_companion.R)
# ============================================================
#
# Run from the project root:  Rscript application/code/thesis_design_examples.R
#
# Same plot as maps$examples in rs_maps() (no numeric panel IDs), without title or
# subtitle (the caption carries them). Reads the frozen setup.rds; nothing is
# simulated. Writes only paper/thesis_chapter/figures/real_sud/design_examples.png;
# the CTJ and SI copies (14 x 8 in) come from manuscript_maps.R.

rs_companion_define_only <- TRUE
source("application/code/render_real_sud_companion.R")

out <- "paper/thesis_chapter/figures/real_sud/design_examples.png"
stopifnot(dir.exists(dirname(out)))
maps <- rs_maps(readRDS("application/results/real_sud_rev_20261002/setup.rds"))
p <- maps$examples +
  ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
  ggplot2::theme(text = ggplot2::element_text(size = 8),
                 strip.text = ggplot2::element_text(size = 8),
                 legend.text = ggplot2::element_text(size = 8),
                 legend.title = ggplot2::element_blank(),
                 legend.position = "bottom", legend.margin = ggplot2::margin(0, 0, 0, 0),
                 plot.margin = ggplot2::margin(2, 2, 2, 2))
ggplot2::ggsave(out, p, width = 6.5, height = 4.2, dpi = 300, bg = "white")
stopifnot(file.exists(out))
cat("Wrote", out, "\n")
