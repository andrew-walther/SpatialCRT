# ============================================================
# Script: 30_igs_sensitivity_summary.R
# Purpose: Summarize the Incidence-Guided Saturation sensitivity study (code/29):
#          each variant's MSE relative to SRS and to the original design, guided vs
#          unguided twins, coverage, bias and burden share treated.
# Author: Andrew Walther
# Created: 2026-11-24
# Dependencies: dplyr, tidyr, ggplot2 (via 10), digest (via 05)
# ============================================================
#
# Run from the project root after code/29 (full profile):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript code/30_igs_sensitivity_summary.R
#
# Contrasts (spec: docs/plans/2026-11-24-igs-sensitivity-spec.md), per estimator x
# neighbor type x spillover regime (and pooled over regimes, labeled "pooled"):
#   RE_SRS  = MSE_d / MSE_SRS            (< 1: lower MSE than SRS)
#   RE_orig = MSE_d / MSE_original       (original = Design 8, R = 4, s = 0.3, guided)
#   paired tests over the 50 (configuration x surface) units with
#   run_surface_paired_test() from 10 (t49 interval, Wilcoxon), Holm within each
#   family: 40 designs vs SRS, 19 guided variants vs the original, and 20
#   guided-vs-unguided pairs (including original vs Saturation Quadrants).
#   Original family sizes were 32, 15 and 16.
#   ratio = 1 + mean paired difference / MSE of the comparator (as in 15's forest plot).
# MSE values pool the five effect sizes, as in code/15; queen and rook are never pooled.
#
# Burden share treated (as code/22): share_i = sum_j X_j z_ji / sum_j X_j, from the
# actual allocations of one block per configuration (queen, rho = 0, psi = 0.5,
# control-only; allocations do not depend on the block except through their seed).
#
# Outputs: results/igs_sensitivity_2026-11/full/summary/
#   igs_contrasts.csv        one row per slice x design: performance, RE, paired tests
#   igs_holm_comparison.csv  original comparisons: raw p and Holm p before/after
#   igs_burden_share.csv     mean and 10th/90th percentile share by configuration x design
#   igs_summary.txt          headline tables (queen, oracle)
#   fig_igs_re_vs_srs.{pdf,png}  RE vs SRS by spread, lines by R, guided/unguided x regime

# Order matters: 10 sources 06, which re-sources 03 and would mask the wrapper below,
# so 10 comes first and the get_designs() override comes last.
source("code/10_statistical_comparisons.R")  # run_surface_paired_test()
sim_define_only <- TRUE
sim_profile     <- "full"
sim_script_dir  <- normalizePath("code")
source("code/05_run_simulation.R")          # grid, configs, surfaces, draw_assignments()
source("code/29_igs_variant_designs.R")
source("code/00_thesis_figure_device.R")
get_designs_prod <- get_designs
get_designs <- function(design_id, n_resamples, N, incidence, nb_list, coords) {
  get_designs_with_variants(design_id, n_resamples, N, incidence, nb_list, coords,
                            prod_fun = get_designs_prod)
}

IN_DIR  <- "results/igs_sensitivity_2026-11/full"
OUT_DIR <- file.path(IN_DIR, "summary")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(file.path(IN_DIR, "reproduction_check.txt")))

# Labels: every design gets R, s, guided ----
info <- rbind(
  data.frame(Design = "Design 3", R = 4, s = 0.3, guided = FALSE,
             label = "Saturation Quadrants (Unguided R=4 s=0.3)"),
  data.frame(Design = "Design 8", R = 4, s = 0.3, guided = TRUE,
             label = "Original IGS (Guided R=4 s=0.3)"),
  data.frame(Design = "Design 9", R = NA, s = NA, guided = NA,
             label = "Simple Random Sampling"),
  data.frame(Design = paste("Design", VARIANTS$id), R = VARIANTS$R, s = VARIANTS$s,
             guided = VARIANTS$guided, label = VARIANTS$label))
SRS <- "Design 9"; ORIG <- "Design 8"

scen <- readRDS(file.path(IN_DIR, "igs_scenarios.rds"))
surf <- readRDS(file.path(IN_DIR, "igs_surfaces.rds"))
stopifnot(setequal(unique(scen$Design), info$Design),
          setequal(unique(surf$Design), info$Design), all(is.finite(surf$MSE)))
stopifnot(all(abs(scen$Mean_Treated - 50) < 1e-9))            # exactly half treated
cat(sprintf("Scenario rows %d; max fail rate %.4f; aliased fits %d; warnings %d\n",
            nrow(scen), max(scen$Fail_Rate), sum(scen$N_Aliased), sum(scen$N_Warn)))

#' Paired surface-level comparisons of designs against one comparator, Holm-adjusted
#'
#' @param su Surface rows for one slice
#' @param designs Designs to compare
#' @param comp Comparator design, or a named vector design -> its own comparator
#' @param mse_comp Named numeric, mean MSE of each comparator in the slice
#' @return Data frame: Design, Comparator, mean_diff, ci_lo, ci_hi, ratio, ratio_lo,
#'   ratio_hi, wilcoxon_p, p_holm, frac_better, n_units
paired_vs <- function(su, designs, comp, mse_comp) {
  comps <- if (length(comp) == 1) setNames(rep(comp, length(designs)), designs) else comp
  rows <- lapply(designs, function(d) {
    cd <- comps[[d]]
    r <- run_surface_paired_test(su, d, cd)
    if (r$n_units != length(inc_configs) * n_surfaces || !is.finite(r$wilcoxon_p)) {
      stop("Incomplete paired comparison: ", d, " vs ", cd)
    }
    m <- mse_comp[[cd]]
    data.frame(Design = d, Comparator = cd, mean_diff = r$mean_diff, ci_lo = r$ci[1],
               ci_hi = r$ci[2], ratio = 1 + r$mean_diff / m, ratio_lo = 1 + r$ci[1] / m,
               ratio_hi = 1 + r$ci[2] / m, wilcoxon_p = r$wilcoxon_p,
               frac_better = r$frac_a_better, n_units = r$n_units)
  })
  out <- do.call(rbind, rows)
  out$p_holm <- p.adjust(out$wilcoxon_p, method = "holm")
  out
}

# Twin of each guided design: same R and s, unguided (Design 3 for the original)
guided_designs <- info$Design[info$guided %in% TRUE]
twin <- vapply(guided_designs, function(d) {
  i <- info[info$Design == d, ]
  info$Design[info$guided %in% FALSE & info$R %in% i$R & info$s %in% i$s]
}, character(1))

# Contrasts by slice ----
slices <- expand.grid(Estimator = c("oracle", "nonoracle"), Neighbor_Type = c("queen", "rook"),
                      Regime = c("control_only", "both", "pooled"), stringsAsFactors = FALSE)
contrasts <- do.call(rbind, lapply(seq_len(nrow(slices)), function(i) {
  sl <- slices[i, ]
  pick <- function(d) d[d$Estimator == sl$Estimator & d$Neighbor_Type == sl$Neighbor_Type &
                        (sl$Regime == "pooled" | d$Spillover_Type == sl$Regime), ]
  sc <- pick(scen); su <- pick(surf)
  perf <- sc %>% group_by(Design) %>%
    summarise(MSE = mean(MSE), Bias = mean(Bias), SD = mean(SD),
              Coverage_Tdf = mean(Coverage_Tdf), Coverage_Wald = mean(Coverage),
              Power_Tdf = mean(Power_Tdf), n_scen = n(), .groups = "drop") %>%
    as.data.frame()
  mse <- setNames(perf$MSE, perf$Design)
  perf$RE_SRS  <- perf$MSE / mse[[SRS]]
  perf$RE_orig <- perf$MSE / mse[[ORIG]]
  v_srs  <- paired_vs(su, setdiff(perf$Design, SRS), SRS, mse)
  v_orig <- paired_vs(su, setdiff(guided_designs, ORIG), ORIG, mse)
  v_twin <- paired_vs(su, guided_designs, twin, mse)
  names(v_srs)[-1]  <- paste0("vsSRS_", names(v_srs)[-1])
  names(v_orig)[-1] <- paste0("vsOrig_", names(v_orig)[-1])
  names(v_twin)[-1] <- paste0("vsTwin_", names(v_twin)[-1])
  out <- Reduce(function(a, b) merge(a, b, by = "Design", all.x = TRUE),
                list(perf, v_srs, v_orig, v_twin))
  data.frame(sl, merge(info, out, by = "Design"), row.names = NULL)
}))
# Compare the original multiplicity families before replacing their summary ----
contrast_file <- file.path(OUT_DIR, "igs_contrasts.csv")
holm_file <- file.path(OUT_DIR, "igs_holm_comparison.csv")
slice_keys <- names(slices)
original <- contrasts[contrasts$Design %in% paste("Design", c(3, 8, 9, ORIGINAL_VARIANTS$id)), ]
families <- c(vsSRS = 32L, vsOrig = 15L, vsTwin = 16L)
holm <- do.call(rbind, lapply(names(families), function(family) {
  cols <- paste0(family, c("_Comparator", "_wilcoxon_p", "_p_holm"))
  old <- original[!is.na(original[[cols[1]]]), c(slice_keys, "Design", cols)]
  names(old)[match(cols, names(old))] <- c("Comparator", "p_raw", "p_holm_after")
  old$Contrast <- family
  old <- old |>
    dplyr::group_by(dplyr::across(dplyr::all_of(slice_keys))) |>
    dplyr::mutate(n_before = dplyr::n(),
                  p_holm_before = p.adjust(p_raw, method = "holm")) |>
    dplyr::ungroup()
  old$n_after <- families[[family]] + if (family == "vsSRS") 8L else 4L
  stopifnot(nrow(old) == nrow(slices) * families[[family]], all(old$n_before == families[[family]]),
            sum(!is.na(contrasts[[cols[1]]])) == nrow(slices) * unique(old$n_after))
  old
}))
holm$crosses_0.05 <- (holm$p_holm_before < 0.05) != (holm$p_holm_after < 0.05)
comparison_keys <- c(slice_keys, "Contrast", "Design", "Comparator")
if (file.exists(holm_file)) {
  archived <- read.csv(holm_file)[c(comparison_keys, "p_raw", "p_holm_before")]
} else {
  previous <- read.csv(contrast_file)
  archived <- do.call(rbind, lapply(names(families), function(family) {
    cols <- paste0(family, c("_Comparator", "_wilcoxon_p", "_p_holm"))
    d <- previous[!is.na(previous[[cols[1]]]), c(slice_keys, "Design", cols)]
    names(d)[match(cols, names(d))] <- c("Comparator", "p_raw", "p_holm_before")
    d$Contrast <- family
    d
  }))
}
# Compare at the precision used by the archived CSV writer and reader.
p_columns <- c(comparison_keys, "p_raw", "p_holm_before")
computed_csv <- capture.output(write.csv(holm[p_columns], row.names = FALSE))
recomputed <- read.csv(text = computed_csv)
p_check <- merge(archived, recomputed, by = comparison_keys, suffixes = c("_stored", "_computed"))
stopifnot(nrow(p_check) == nrow(holm), nrow(p_check) == nrow(archived),
          !anyDuplicated(p_check[comparison_keys]))
for (nm in c("p_raw", "p_holm_before")) {
  stored <- p_check[[paste0(nm, "_stored")]]
  computed <- p_check[[paste0(nm, "_computed")]]
  if (!identical(stored, computed)) {
    stop("November 24 p-values differ at stored CSV precision: ", nm)
  }
}
write.csv(holm, holm_file, row.names = FALSE)
write.csv(contrasts, contrast_file, row.names = FALSE)

# Burden share treated ----
burden <- do.call(rbind, lapply(seq_along(inc_configs), function(ci) {
  cfg <- inc_configs[[ci]]
  X <- generate_surfaces(cfg)
  stopifnot(all(X >= 0))
  surface <- rep(seq_len(ncol(X)), each = n_design_draw)
  do.call(rbind, lapply(info$Design, function(d) {
    d_id <- as.integer(sub("Design ", "", d))
    Z <- draw_assignments(cfg, "queen", 0, 0.5, "control_only", d_id, X)
    sh <- vapply(seq_len(ncol(Z)), function(f) sum(X[, surface[f]] * Z[, f]) /
                   sum(X[, surface[f]]), numeric(1))
    data.frame(Incidence_Mode = cfg$mode, Rho_Incidence = cfg$rho_x, Design = d,
               share_mean = mean(sh), share_p10 = unname(quantile(sh, 0.1)),
               share_p90 = unname(quantile(sh, 0.9)))
  }))
}))
burden <- merge(info, burden, by = "Design")
write.csv(burden, file.path(OUT_DIR, "igs_burden_share.csv"), row.names = FALSE)

# Headline text (queen, oracle) ----
fmt <- function(d) {
  d <- d[order(is.na(d$guided), -d$guided, d$R, d$s), ]
  with(d, sprintf("%-42s MSE %.3f  RE_SRS %.3f [%.3f, %.3f] holm %.3f  RE_orig %.3f  vsTwin %s  cov_t %.3f  bias %.3f",
                  label, MSE, RE_SRS, vsSRS_ratio_lo, vsSRS_ratio_hi, vsSRS_p_holm, RE_orig,
                  ifelse(is.na(vsTwin_ratio), "   -  ",
                         sprintf("%.3f", vsTwin_ratio)), Coverage_Tdf, Bias))
}
sink(file.path(OUT_DIR, "igs_summary.txt"), split = TRUE)
cat("Incidence-Guided Saturation sensitivity study: queen, oracle fit\n")
cat("Pooled over five incidence configurations and five intervention effect sizes.\n")
cat("RE_SRS = MSE / MSE_SRS with the paired-unit ratio interval; RE_orig = MSE / MSE_original;\n")
cat("vsTwin = guided / unguided twin (paired ratio). cov_t = t + df coverage.\n")
cat(readLines(file.path(IN_DIR, "reproduction_check.txt")), sep = "\n")
cat(sprintf("Original Holm comparisons crossing 0.05: %d of %d across all slices.\n",
            sum(holm$crosses_0.05), nrow(holm)))
for (rg in c("control_only", "both", "pooled")) {
  cat(sprintf("\n== Spillover regime: %s%s ==\n", rg, if (rg == "pooled") " (pooled over regimes)" else ""))
  d <- contrasts[contrasts$Estimator == "oracle" & contrasts$Neighbor_Type == "queen" &
                   contrasts$Regime == rg, ]
  cat(fmt(d[d$Design != SRS, ]), sep = "\n")
  cat(sprintf("SRS MSE %.3f\n", d$MSE[d$Design == SRS]))
}
cat("\n== Burden share treated (queen block, mean over configurations) ==\n")
b <- burden %>% group_by(label, R, s, guided) %>% summarise(share = mean(share_mean), .groups = "drop") %>%
  arrange(is.na(guided), desc(guided), R, s)
cat(sprintf("%-42s %.3f", b$label, b$share), sep = "\n")
sink()

# Figure: RE vs SRS by spread (queen, oracle) ----
fd <- contrasts[contrasts$Estimator == "oracle" & contrasts$Neighbor_Type == "queen" &
                  contrasts$Regime != "pooled" & contrasts$Design != SRS, ]
fd$Guidance <- ifelse(fd$guided, "Incidence-guided", "Unguided")
fd$Regions  <- factor(paste(fd$R, "regions"), levels = paste(c(2, 4, 10, 25), "regions"))
fd$Regime   <- factor(ifelse(fd$Regime == "control_only", "Spillover to control clusters only",
                             "Spillover to both arms"),
                      levels = c("Spillover to control clusters only", "Spillover to both arms"))
# Apply the same horizontal offset to the line, interval, point and its circle.
fd$plot_s <- fd$s + (as.integer(fd$Regions) - 2.5) * 0.011
ring_label <- "Original design and Saturation Quadrants"
p <- ggplot(fd, aes(plot_s, vsSRS_ratio, colour = Regions)) +
  geom_hline(yintercept = 1, colour = "grey40") +
  geom_errorbar(aes(ymin = vsSRS_ratio_lo, ymax = vsSRS_ratio_hi), width = 0.008) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_point(data = fd[fd$Design %in% c(ORIG, "Design 3"), ], aes(shape = ring_label),
             size = 3.7, colour = "black") +
  facet_grid(Regime ~ Guidance, scales = "free_y") +
  scale_x_continuous(breaks = seq(0.1, 0.5, by = 0.1),
                     labels = c("60/40", "70/30", "80/20", "90/10", "100/0")) +
  scale_y_continuous(breaks = function(limits) sort(unique(c(1, scales::breaks_extended()(limits)))),
                     labels = scales::label_number(style_negative = "minus")) +
  scale_colour_viridis_d(end = 0.85) +
  scale_shape_manual(values = setNames(1, ring_label)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "Highest/lowest saturation (%)", y = "MSE relative to SRS",
       colour = NULL, shape = NULL) +
  theme_bw(base_size = 10) + thesis_text_sizes() +
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text = element_text(colour = "black"))
# Quartz embeds fonts and preserves the hyphen glyph on macOS when Cairo is unavailable.
pdf_path <- file.path(OUT_DIR, "fig_igs_re_vs_srs.pdf")
if (capabilities("aqua")) {
  grDevices::quartz(type = "pdf", file = pdf_path, width = 6.5, height = 5.5)
  print(p)
  invisible(grDevices::dev.off())
} else {
  save_thesis_pdf(p, pdf_path, width = 6.5, height = 5.5)
}
ggsave(file.path(OUT_DIR, "fig_igs_re_vs_srs.png"), p, width = 6.5, height = 5.5, dpi = 300)
cat("Saved to", OUT_DIR, "\n")
