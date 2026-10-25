# ============================================================
# Script: 23_methods_extension_figures.R
# Purpose: Figures for the NC methods extensions (partial participation, power,
#          targeting) from the stored 2026-10-21 summaries; no fitting or simulation.
# Author: Andrew Walther
# Created: 2026-10-25
# Dependencies: ggplot2
# ============================================================
#
# Run from the project root:  Rscript code/23_methods_extension_figures.R
#
# Inputs (read only):
#   application/results/methods_ext_2026-10/summaries_v2/a3_settings.csv   power by N
#   application/results/methods_ext_2026-10/summaries_v2/a3_n80_summary.csv colleges for 80% power
#   application/results/methods_ext_2026-10/stored/a1_burden_share.csv     NC rate shares
#   application/results/real_sud_rev_20261002/setup.rds                    NC regions, yearly ranks
#   results/grid_burden_share/grid_quadrant_spread.csv                     grid spread and gain
#
# Outputs (results/methods_ext_figures/, each as .pdf and .png):
#   fig_power_vs_n          t + df power against the number of participating colleges
#   fig_colleges_needed     colleges needed for 80% t + df power, all designs
#   fig_spread_vs_gain      grid quadrant spread vs incidence-guidance gain, NC marked
#
# Conventions (author decisions 2026-10-21): t + df-corrected power is the headline;
# random volunteering, beta = 0, kappa = 0.8 (spillover gamma = kappa * tau); the four
# planning years are summarized by their mean (line or bar) and range (band or bar);
# "Plain Saturation" in the stored files is displayed as "Saturation Regions".

library(ggplot2)

R_EXT   <- "application/results/methods_ext_2026-10"
OUT_DIR <- "results/methods_ext_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DISPLAY <- c("Plain Saturation" = "Saturation Regions")
display_name <- function(x) ifelse(x %in% names(DISPLAY), DISPLAY[x], x)
REGIME_LAB <- c(both = "Both-arms spillover", control_only = "Control-only spillover")

#' Save one ggplot as PDF (manuscripts) and PNG (HTML outline, slides)
#'
#' @param p ggplot object
#' @param name File stem in OUT_DIR
#' @param w,h Size in inches
#' @return Invisibly, the two paths
save_both <- function(p, name, w, h) {
  paths <- file.path(OUT_DIR, paste0(name, c(".pdf", ".png")))
  ggsave(paths[1], p, width = w, height = h)
  ggsave(paths[2], p, width = w, height = h, dpi = 200)
  invisible(paths)
}

# Figure 1: power against the number of colleges ----
#' Power-vs-N data: t + df power, mean and range over the four years
#'
#' Power_N(design) = mean over years y of Power_Tdf(y, N, design); the band is the
#' minimum and maximum over years. Random volunteering, beta = 0, kappa = 0.8, rho = 0.5.
#'
#' @param settings a3_settings.csv as a data frame
#' @param designs Design_IDs to keep
#' @param taus Effects (college-level SD units) to keep
#' @return Data frame: Design, Regime, Tau, N, Power, Lo, Hi
power_curve_data <- function(settings, designs, taus) {
  x <- settings[settings$Mechanism == "random" & settings$Beta == 0 &
                  settings$Kappa == 0.8 & settings$Rho == 0.5 &
                  settings$Design_ID %in% designs & settings$Tau %in% taus, ]
  # Every (design, regime, tau, N) cell must have all four years: fail loud otherwise
  cells <- table(paste(x$Design_ID, x$Regime, x$Tau, x$N))
  stopifnot(nrow(x) > 0, all(cells == 4))
  agg <- aggregate(Power_Tdf ~ Design + Regime + Tau + N, x,
                   function(v) c(Power = mean(v), Lo = min(v), Hi = max(v)))
  out <- data.frame(agg[, c("Design", "Regime", "Tau", "N")], agg$Power_Tdf)
  out$Design <- display_name(out$Design)
  out
}

settings <- read.csv(file.path(R_EXT, "summaries_v2/a3_settings.csv"))
pc <- power_curve_data(settings, designs = c(9, 3, 8, 1), taus = c(1, 1.5))
pc$Design <- factor(pc$Design, levels = c("Simple Random Sampling", "Saturation Regions",
                                          "Incidence-Guided Saturation", "Checkerboard"))
pc$Regime <- factor(REGIME_LAB[pc$Regime], levels = REGIME_LAB)
pc$TauLab <- factor(sprintf("tau == %s", pc$Tau), levels = sprintf("tau == %s", c(1, 1.5)))

pal4 <- c("Simple Random Sampling" = "grey25", "Saturation Regions" = "#009E73",
          "Incidence-Guided Saturation" = "#0072B2", "Checkerboard" = "#D55E00")
p_power <- ggplot(pc, aes(N, Power, colour = Design, fill = Design)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3) +
  facet_grid(TauLab ~ Regime, labeller = labeller(TauLab = label_parsed)) +
  scale_colour_manual(values = pal4, name = NULL) +
  scale_fill_manual(values = pal4, name = NULL) +
  scale_x_continuous(breaks = c(15, 20, 29, 40, 50, 58)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Participating colleges", y = "Power (t test, df-corrected)",
       caption = expression("Random volunteering; " * rho * " = 0.5, " * gamma * " = 0.8" * tau *
                              "; lines: mean of 2018-2021 planning years, bands: range. Dashed: 80%.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8), plot.caption.position = "plot") +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
save_both(p_power, "fig_power_vs_n", 7, 5.2)

# Figure 2: colleges needed for 80% power ----
#' Colleges-needed data for one regime and effect
#'
#' N80 is the interpolated number of participating colleges at which t + df power
#' reaches 0.8 (a3_n80_summary.csv, monotone interpolation over N). A design that
#' does not reach 0.8 with all 58 colleges in some years has no mean; it is drawn to
#' 58 and labeled with its range instead, never dropped.
#'
#' @param n80 a3_n80_summary.csv as a data frame
#' @param regime "both" or "control_only"
#' @param tau Effect in college-level SD units
#' @return Data frame: Design, Panel, Value, Lo, Hi, Reached, Label, Label_x
n80_data <- function(n80, regime, tau) {
  x <- n80[n80$Kappa == 0.8 & n80$Rho == 0.5 & n80$Regime == regime & n80$Tau == tau, ]
  stopifnot(nrow(x) == 9)
  reached <- is.finite(x$N80_Tdf_Mean)
  data.frame(
    Design = display_name(x$Design),
    Panel = sprintf('"%s,"~tau == %s', REGIME_LAB[[regime]], tau),  # plotmath label
    Value = ifelse(reached, x$N80_Tdf_Mean, 58),
    Lo = ifelse(reached, x$N80_Tdf_Min, NA), Hi = ifelse(reached, x$N80_Tdf_Max, NA),
    Reached = ifelse(reached, "Reached in every year", "Not reached in some years (count shown)"),
    Label = ifelse(reached, sprintf("%.0f", x$N80_Tdf_Mean),
                   sprintf("> 58 (%d/4 years)", x$N80_Tdf_Years_Above58)),
    Label_x = ifelse(reached, x$N80_Tdf_Max, 58) + 1.5)
}

n80 <- read.csv(file.path(R_EXT, "summaries_v2/a3_n80_summary.csv"))
nd <- rbind(n80_data(n80, "both", 1), n80_data(n80, "control_only", 1.5))
nd$Panel <- factor(nd$Panel, levels = unique(nd$Panel))
ord <- nd[nd$Panel == levels(nd$Panel)[2], ]
nd$Design <- factor(nd$Design, levels = rev(ord$Design[order(ord$Value, ord$Design)]))
p_n80 <- ggplot(nd, aes(Value, Design, fill = Reached)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = Lo, xmax = Hi), width = 0.25, orientation = "y", na.rm = TRUE) +
  geom_text(aes(x = Label_x, label = Label), hjust = 0, size = 3) +
  geom_vline(xintercept = 58, linetype = "dotted", colour = "grey40") +
  facet_wrap(~ Panel, labeller = label_parsed) +
  scale_fill_manual(values = c("Reached in every year" = "#0072B2",
                               "Not reached in some years (count shown)" = "grey75"), name = NULL) +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 15, 29, 45, 58)) +
  labs(x = "Colleges needed for 80% power (t test, df-corrected)", y = NULL,
       caption = expression("Random volunteering; " * rho * " = 0.5, " * gamma * " = 0.8" * tau *
                              "; bars: mean of 2018-2021, whiskers: range. Dotted: all 58 colleges.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8), plot.caption.position = "plot")
save_both(p_n80, "fig_colleges_needed", 7.5, 4)

# Figure 3: quadrant spread vs incidence-guidance gain ----
#' NC reference points: one per planning year
#'
#' spread_y = (max_r - min_r of region mean incidence rank in year y) / 58 over the
#' four frozen production regions; gain_y = rate share treated by Incidence-Guided
#' Saturation minus Saturation Regions (A1, production allocations).
#'
#' @param setup Production setup.rds list (regions, annual)
#' @param a1 stored/a1_burden_share.csv as a data frame
#' @return Data frame: Year, Spread, Gain
nc_points <- function(setup, a1) {
  years <- c("2018", "2019", "2020", "2021")
  do.call(rbind, lapply(years, function(y) {
    ann <- setup$annual[[y]]
    stopifnot(identical(ann$Primary_College, setup$ids), nrow(ann) == 58)
    region_rank <- tapply(ann$rank, setup$regions, mean)
    stopifnot(length(region_rank) == 4)
    rs <- a1[a1$Year == as.integer(y), ]
    data.frame(Year = y, Spread = diff(range(region_rank)) / 58,
               Gain = rs$Rate_Share_Mean[rs$Design_ID == 8] - rs$Rate_Share_Mean[rs$Design_ID == 3])
  }))
}

grid <- read.csv("results/grid_burden_share/grid_quadrant_spread.csv")
CONFIG_LAB <- c("iid" = "iid", "spatial rho_X=0.2" = "Spatial, 0.2", "spatial rho_X=0.5" = "Spatial, 0.5",
                "poisson rho_X=0.2" = "Poisson, 0.2", "poisson rho_X=0.5" = "Poisson, 0.5")
stopifnot(all(grid$Config %in% names(CONFIG_LAB)))
grid$Config <- factor(CONFIG_LAB[grid$Config], levels = CONFIG_LAB)
setup <- readRDS("application/results/real_sud_rev_20261002/setup.rds")
a1 <- read.csv(file.path(R_EXT, "stored/a1_burden_share.csv"))
nc <- nc_points(setup, a1)
write.csv(nc, file.path(OUT_DIR, "nc_spread_gain.csv"), row.names = FALSE)

p_spread <- ggplot(grid, aes(Spread, Gain)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50",
              linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(colour = Config), size = 2, alpha = 0.85) +
  geom_point(data = nc, shape = 17, size = 3, colour = "black") +
  geom_text(data = nc, aes(label = paste("NC", Year)), hjust = -0.2, vjust = 0.5, size = 2.8) +
  scale_colour_viridis_d(end = 0.85, name = expression("Grid incidence, " * rho[X])) +
  scale_y_continuous(labels = function(v) sprintf("%+.0f", 100 * v)) +
  labs(x = "Spread of regional mean incidence ranks (max - min) / number of clusters",
       y = "Treated incidence share, guided minus plain\n(percentage points)",
       caption = paste("Grid: 10 surfaces per incidence configuration, queen, production allocations.",
                       "NC: four frozen regions, one point per planning year (triangles).")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8), plot.caption.position = "plot")
save_both(p_spread, "fig_spread_vs_gain", 7, 4.2)

cat("NC reference points (spread, gain):\n"); print(nc, row.names = FALSE)
cat("Written to", OUT_DIR, "\n")
