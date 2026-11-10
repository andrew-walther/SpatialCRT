# ============================================================
# Script: 00_thesis_figure_device.R
# Purpose: Helpers for Chapter 3 figures: save a ggplot as a PDF with embedded
#          fonts, and set every text element to at least 8 pt.
# Author: Andrew Walther
# Created: 2026-11-10
# Dependencies: ggplot2; poppler's pdftocairo (fallback when cairo_pdf cannot load)
# ============================================================
#
# Chapter rules (2026-11-10): thesis-chapter figures are saved with cairo_pdf so
# fonts are embedded, and all figure text is at least 8 pt at the printed size.
# Figures are saved at the width they print (6.5 in = \linewidth) so the PDF is
# not rescaled in the chapter.
#
# On this machine grDevices::cairo_pdf() cannot load (R_X11.so needs XQuartz's
# libSM). save_thesis_pdf() therefore tries cairo_pdf and, if it fails, draws with
# the default pdf() device and runs `pdftocairo -pdf`, which re-writes the page
# through cairo with the fonts embedded (TrueType subsets). Vector content is kept.

THESIS_TEXT_PT <- 8  # smallest text size (pt) allowed in a chapter figure

#' Save a plot as a PDF with embedded fonts
#'
#' @param plot A ggplot or patchwork object
#' @param path Output .pdf path
#' @param width,height Size in inches (the printed size)
#' @return Invisibly, the method used: "cairo_pdf" or "pdftocairo"
#' @examples
#' # save_thesis_pdf(ggplot(mtcars, aes(wt, mpg)) + geom_point(), "t.pdf", 6.5, 4)
save_thesis_pdf <- function(plot, path, width, height) {
  if (file.exists(path)) file.remove(path)   # so a failed device cannot leave the old file
  ok <- tryCatch({
    ggplot2::ggsave(path, plot, device = grDevices::cairo_pdf,
                    width = width, height = height)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (ok && file.exists(path) && file.size(path) > 0) return(invisible("cairo_pdf"))
  tmp <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(tmp, plot, device = "pdf", width = width, height = height)
  status <- system2("pdftocairo", c("-pdf", shQuote(tmp), shQuote(path)))
  if (status != 0 || !file.exists(path)) stop("pdftocairo failed for ", path)
  invisible("pdftocairo")
}

#' Theme layer that sets every text element to at least THESIS_TEXT_PT
#'
#' Adds explicit sizes over whatever theme the plot already has (face, hjust and
#' colour are inherited). Text drawn by geom_text()/annotate() is sized in mm and
#' must be set separately: 8 pt = 8 / .pt mm.
#'
#' @param axis_title Size (pt) of axis titles
#' @return A ggplot2 theme object to add to a plot
thesis_text_sizes <- function(axis_title = 9) {
  s <- THESIS_TEXT_PT
  ggplot2::theme(
    text         = ggplot2::element_text(size = s),
    axis.text    = ggplot2::element_text(size = s),
    axis.title   = ggplot2::element_text(size = axis_title),
    strip.text   = ggplot2::element_text(size = s),
    legend.text  = ggplot2::element_text(size = s),
    legend.title = ggplot2::element_text(size = s),
    plot.title   = ggplot2::element_text(size = axis_title + 1),
    plot.subtitle = ggplot2::element_text(size = s),
    plot.caption = ggplot2::element_text(size = s))
}
