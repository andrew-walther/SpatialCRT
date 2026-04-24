# Modular Incidence Simulation for Spatial CRT Design Evaluation

> For AI session context and quick technical reference, see [CLAUDE.md](CLAUDE.md).

## Overview

This project evaluates which **treatment assignment design** produces the most accurate
estimates of an intervention effect when outcome incidence is spatially heterogeneous
and spillover effects are present. The setting is a Spatial Cluster Randomized Trial
(CRT) on a 10x10 regular lattice grid of 100 clusters.

**Research question:** Across a wide range of spatial dependence, spillover magnitude,
and incidence structures, which of 8 candidate treatment assignment strategies minimizes
estimation error (MSE) and maintains valid inferential coverage?

**Application context:** Sudden Unexpected Death (SUD) in North Carolina counties.
The Poisson incidence mode uses a SUD base rate of 35/100,000 (Mirzaei et al.),
consistent with county-level Poisson analyses by Gan et al. and Watson et al.

---

## Project Lineage

This is a clean-room rewrite of
`archive/OutcomeIncidenceDesign_Legacy/SpatialCRT_Incidence_TreatmentAssignment_Simulation.Rmd`
(~600 lines, monolithic Rmd). Refactoring goals:

1. **Modularity** — numbered R scripts (01-09) that can be run independently or sourced in sequence
2. **Three incidence modes** — Added Spatial (SAR filter) and Poisson modes alongside iid Uniform
3. **Correctness fixes** — Seven specific bugs and design flaws identified and corrected (see below)
4. **Long-run reliability** — Replaced `save.image()` with targeted `saveRDS()`, added per-scenario seeding, ETA logging, and degenerate case handling

The predecessor Rmd is preserved untouched.

---

## File Inventory

| File | Lines | Purpose | Key Functions |
|------|------:|---------|---------------|
| `00_mathematical_specification.Rmd` | ~495 | Theory doc with full LaTeX DGP formulas | (rendered to HTML) |
| `01_spatial_setup.R` | 54 | Regular lattice grid, rook/queen weight matrices | `build_spatial_grid()`, `get_active_spatial()` |
| `02_incidence_generation.R` | 112 | Three incidence generation modes | `generate_incidence()` + 3 mode functions |
| `03_designs.R` | ~162 | Eight treatment assignment strategies | `get_designs()`, `get_design_names()`, `is_design_deterministic()` |
| `04_estimation.R` | 102 | DIM and MLE estimation with CI extraction | `estimate_tau()` |
| `05_run_simulation.R` | ~385 | Main orchestrator, nested parameter loop | `run_incidence_config()` |
| `06_visualizations.R` | ~800 | All plots and summary tables | 18 functions; entry point `run_all_visualizations()` |
| `07_results_summary.Rmd` | ~800 | Rendered results report (HTML/PDF) | Knitted summary of all MLE findings |
| `08_design_recommendations.R` | ~560 | Personalized design recommendations | `run_recommendation_report()`, `table_scenario_lookup()`, `generate_commentary()` |
| `09_MLE_design_recommendation_report.Rmd` | ~984 | Companion narrative PDF report | Knitted to `results/MLE_design_recommendation_report.pdf` |
| `10_statistical_comparisons.R` | ~1050 | Formal statistical hypothesis tests | `run_friedman_test()`, `run_nemenyi_posthoc()`, `run_pairwise_wilcoxon()`, `run_conditional_tests()`, `plot_cd_diagram()`, `plot_mse_boxplot_with_stars()`, `plot_pvalue_heatmap()` |
| `11_statistical_comparisons_report.qmd` | ~543 | Statistical comparisons narrative report | Rendered to `results/11_statistical_comparisons_report.{html,pdf}` |
| `complete_after_mle.R` | ~250 | Post-completion script (viz + docs + stats) | Run once after MLE finishes |

---

## How to Run

### Prerequisites

```r
install.packages(c("spdep", "spatialreg", "dplyr", "tidyr",
                   "digest", "parallel", "ggplot2", "viridis"))
```

### Quick Start — DIM estimation (~5-10 minutes)

```r
# In 05_run_simulation.R, ensure line 40 reads:
#   estimation_mode <- "DIM"

setwd("projects/IncidenceDesign/code")
source("05_run_simulation.R")
```

### Full Run — MLE estimation (~10 hours)

```bash
# Run from terminal to keep process independent of IDE session
# Use caffeinate to prevent sleep on macOS:
#   caffeinate -i -w $(pgrep -f 05_run_simulation) &
cd projects/IncidenceDesign/code
Rscript 05_run_simulation.R
```

In `05_run_simulation.R`, set `estimation_mode <- "MLE"` (line 40).
Optionally set `n_cores > 1` (line 44) for parallel execution across incidence configs.

### Generating Visualizations Only

```r
source("06_visualizations.R")

# All configs, most recent MLE results (primary estimator):
run_all_visualizations(estimation_mode = "MLE_tau_sweep")

# Single incidence config with all tables:
r <- load_latest_results(estimation_mode = "MLE_tau_sweep")
cfgs <- split_by_incidence_config(r)
run_standard_tables(cfgs[["iid Uniform"]], "iid Uniform")
```

### Generating Design Recommendations

```r
source("08_design_recommendations.R")

# Full recommendation report (PDF + console output):
run_recommendation_report(estimation_mode = "MLE_tau_sweep")

# Specific scenario lookup:
r <- load_latest_results(estimation_mode = "MLE_tau_sweep")
cfgs <- split_by_incidence_config(r)
table_scenario_lookup(cfgs[["iid Uniform"]], rho = 0.2, gamma = 0.7,
                      spill_type = "both", nb_type = "queen")
```

---

## Simulation Design Summary

**Grid:** 10x10 regular lattice, N = 100 clusters. Two contiguity structures:
rook (4-connected) and queen (8-connected).

**Data Generating Process — Spatial Durbin Model (SDM):**

```
Y = (I - rho*W)^{-1} * [tau*Z + gamma*Spill(Z) + beta*X + epsilon]
```

where tau = 1.0, beta = 1.0, sigma = 1.0, and `Spill(Z)` is the row-standardized
mean treatment of neighbors. Two spillover modes: `control_only` (gamma applied only
to control units) and `both` (applied to all units).

**Incidence modes:**
- **iid Uniform:** X_i ~ Uniform(0,1) independently
- **Spatial:** SAR filter + pnorm transform — spatially correlated but marginally Uniform(0,1)
- **Poisson:** Spatially correlated log-rates -> Poisson counts -> rank-normalized rates in [0,1]

**Treatment designs:** 8 strategies — Checkerboard, High Incidence Focus, Saturation
Quadrants, Isolation Buffer, 2x2 Blocking, Balanced Quartiles, Balanced Halves,
Incidence-Guided Saturation Quadrants.

<p align="center">
  <img src="results/figures/design_samples_8panel.png" alt="Sample treatment assignments for each of the 8 designs on a 10x10 lattice with iid Uniform baseline incidence" width="75%">
  <br><em>Sample treatment assignments for each design applied to a single realization of iid Uniform(0,1) baseline incidence. Tile color indicates baseline incidence (darker = higher). Circles = treated, crosses = control.</em>
</p>

**Estimation:** MLE via spatial autoregressive model `lagsarlm(Y ~ Z + Spill + X)` using oracle Spill covariate (250 iterations/scenario).

**True tau (tau-sweep):** τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0} — swept to assess design robustness across effect sizes and estimate power curves.

**Total scenarios:** 5 tau × 5 incidence configs × 2 neighbor types × 4 rho × 4 gamma × 2 spillover types × 8 designs = **12,800 scenarios**
(Baseline τ=1.0 slice = 2,560 scenarios, completed 2026-03-22)

---

## Results Summary

### Tau-sweep (COMPLETE — 2026-04-08, 12,800 scenarios)

- **12,800 scenarios** across τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0} | Fail_Rate = 0.0 | N_Valid_Est = 250 (all)
- Primary scenario (τ=1.0): **Best: D8** MSE=0.079 ≈ **D3** MSE=0.080 | **Worst: D1** MSE=0.802, coverage ~55%
- D3/D8 dominance holds across **all τ levels** (Friedman p < 2.2×10⁻¹⁶ at each τ)
- Power curves: D3/D8 reach 80% power already at τ=0.8 (smallest tested effect); D1 requires τ≥2.0
- Files: `results/sim_data/sim_results_MLE_tau_sweep_combined_20260408_191916.rds` (12,800 rows)
- All reports regenerated with tau sensitivity sections (2026-04-08)

### MLE baseline (tau=1.0, archived — superseded by tau-sweep)

- `results/sim_data/sim_results_MLE_combined_20260322_151030.rds` (2,560 rows, preserved)
- Pre-sweep deliverables: `results/archive/pre_tau_sweep_20260408/`

### Statistical Comparisons (updated for tau-sweep, 2026-04-08)

- **Friedman test by tau level:** χ²=1052.77 (τ=0.8) through χ²=743.41 (τ=3.0), all p < 2.2×10⁻¹⁶
- **Top equivalence group:** Design 3 and Design 8 — not significantly different from each other, significantly better than all others at all tau levels
- **Design 1 (Checkerboard) is significantly worse than all alternatives** at every τ
- Full report (with tau-stratified tests + power analysis): `results/11_statistical_comparisons_report.pdf`

### Results Directory Structure

```
results/
  MLE_tau_sweep_design_recommendations.pdf    # PRIMARY — figures/tables PDF
  MLE_tau_sweep_incidence_overview.pdf        # Incidence heatmaps + distributions
  09_MLE_design_recommendation_report.{html,pdf}  # Narrative rec report with tau sensitivity
  11_statistical_comparisons_report.{html,pdf}    # PRIMARY — formal hypothesis testing
  00_mathematical_specification.pdf           # Theory document
  07_results_summary.pdf                      # Results summary (with tau sensitivity section)
  sim_data/
    sim_results_MLE_tau_sweep_combined_20260408_191916.rds   # PRIMARY — 12,800 rows
    sim_results_MLE_tau_sweep_{iid,spatial,poisson}_{ts}.rds # Per-incidence splits
    sim_results_MLE_combined_20260322_151030.rds             # Archived baseline (2,560 rows)
  mle_per_config/
    MLE_tau_sweep_{config_name}.pdf             # Per-config 8-plot PDF (5 configs, tau=1.0)
    MLE_tau_sweep_{config_name}_tau_sensitivity.pdf  # Per-config tau sensitivity PDFs (5)
    MLE_tau_sweep_incidence_overview.pdf        # Incidence heatmaps + distributions
  archive/
    pre_tau_sweep_20260408/                     # All pre-sweep deliverables (archived)
  figures/
    design_samples_8panel.{png,pdf}            # 8-panel design illustration (clean)
    design_samples_option1_overlays.{png,pdf}  # 8-panel with saturation % annotations
  archive/
    test_plots.pdf                             # Dev artifacts
    completion_log.txt
```

**Key rule:** Results for iid Uniform, Spatial, and Poisson are always reported
**separately** — never aggregated. The combined .rds exists for loading convenience only.
`load_latest_results()` automatically detects the `sim_data/` subdirectory.

### Paper Directory Structure

```
paper/
  spatialCRT.bib                              # Zotero bibliography export
  SpatialCRT_IncidenceDesign_Presentation.qmd # Presentation slides (scaffold)
  report/
    IncidenceSpatialCRT_Report.qmd            # Unified project report (sources R modules)
    IncidenceSpatialCRT_Report.pdf            # Rendered PDF
    IncidenceSpatialCRT_Report.html           # Rendered HTML
    IncidenceDesign_ProjectSummary.qmd        # Brief project summary report
    IncidenceDesign_ProjectSummary.pdf        # Rendered PDF summary
    IncidenceDesign_ProjectSummary.html       # Rendered HTML summary
  manuscript/
    IncidenceSpatialCRT_Manuscript.qmd        # Master document (includes child sections)
    _abstract.qmd                             # Structured abstract
    _introduction.qmd                         # Lit review, motivation, aims
    _methods.qmd                              # Spatial structure, SDM, 8 designs, MLE
    _simulation.qmd                           # Parameter space, metrics, results (WIP)
    _application.qmd                          # SUD in NC (skeleton)
    _discussion.qmd                           # Relevance, limitations, future (skeleton)
  section_drafts/                             # Archival LaTeX drafts (Gemini, 6-design era)
    abstract.tex, introduction.tex, methods.tex, simulation.tex
    application.tex, disussion.tex, main.tex
    references.bib                            # Canonical bibliography (113K, Zotero)
    figures/SamplingGridExamples_crop.png
```

The **unified report** (`paper/report/`) consolidates all code-side reports into a single
end-to-end reference (50+ pages): spatial setup, DGP, 8 design illustrations, estimation
methods, simulation design, all MLE results, full statistical comparisons (Section 10 with
Friedman/Nemenyi/Wilcoxon tests and conditional CD diagrams), and design recommendations.

The **project summary report** (`paper/report/IncidenceDesign_ProjectSummary.{qmd,html,pdf}`)
is a shorter companion document that presents the same core findings in a more concise format
for quick review and sharing.

The **manuscript** (`paper/manuscript/`) uses a modular Quarto structure with
`{{< include >}}` directives. Each section is independently editable. The LaTeX
originals in `section_drafts/` are preserved as archival reference.

Figures available in `paper/manuscript/figures/`:
- `design_samples_8panel.{png,pdf}` — clean 2×4 design panel
- `design_samples_option1_overlays.{png,pdf}` — annotated version with saturation % labels

---

## Key Design Decisions

1. **Incidence generated once per `(mode, rho_X)` config** — not re-drawn per `(nb_type, rho)`.
   In practice a researcher observes historical incidence once; spatial model parameters
   affect outcome propagation, not incidence itself.

2. **Deterministic design detection** — Designs 1 and 2 yield the same assignment for
   every design resample (given fixed incidence). They run once and replicate, saving ~25x compute.

3. **Per-scenario deterministic seeding** — `digest::digest2int(paste(params, sep="|"))` ensures
   reproducibility and order-independence. Adding/removing scenarios doesn't affect others.

4. **Separate reporting by incidence mode** — The three modes represent fundamentally
   different assumptions about prior data. Aggregating across them obscures rather than informs.

5. **Oracle spillover in MLE** — `lagsarlm()` includes the true Spill covariate.
   A `include_spill_covariate = FALSE` toggle exists in `estimate_tau()` for realistic
   estimation comparisons (not yet run systematically).

6. **DIM iteration counts** — Originally planned at 100,000 iterations/scenario; reduced
   to 2,500 (25 design x 100 outcome resamples) for practical runtime.

7. **All 8 designs included** — Design IDs 1–8 are all included in the default
   `design_ids` sweep in `05_run_simulation.R`.

---

## Known Issues & Bugs Fixed

| # | Issue | Root Cause | Fix |
|---|-------|-----------|-----|
| 1 | PDF render failure | `\u03c1` (Unicode rho) unsupported by `pdf()` device | Replace with literal `"rho"` text |
| 2 | Pairwise dominance NaN | `make.names()` mangled column names in `pivot_wider` output | Use `wide[[designs[i]]]` direct bracket access |
| 3 | Eta-squared > 1.0 | `var(group_means)/var(total)` incorrect formula | Proper: `SS_between/SS_total` |
| 4 | `I()` namespace conflict | `I <- diag(N)` shadowed `base::I()` | Renamed to `I_mat` |
| 5 | Degenerate Z crash | Designs occasionally produce all-treated or all-control | Early `NA` return in `estimate_tau()` |
| 6 | R sprintf 10k limit | `complete_after_mle.R` format string too long | Write docs via Claude Code `Write` tool |

---

## Status & Next Steps

**Completed:**
- [x] All code files (00–11) written and tested
- [x] Mathematical specification document (00) rendered to HTML
- [x] DIM simulation: 1,920 scenarios (prior 6-design sweep), visualizations generated (baseline)
- [x] MLE simulation: 1,920 scenarios (prior 6-design sweep), visualizations generated, zero convergence failures
- [x] Design recommendations module (08) with validation tests
- [x] Narrative PDF report (09): `results/MLE_design_recommendation_report.pdf`
- [x] Results directory reorganized: `sim_data/`, `mle_per_config/`, `dim/`, `archive/`
- [x] Project documentation: CLAUDE.md + README.md
- [x] Design set expanded to 8 designs (Balanced Halves, Incidence-Guided Saturation Quadrants added)
- [x] Re-run MLE simulation: 2,560 scenarios covering all 8 designs (1–8) — completed 2026-03-22
- [x] Unified project report consolidating all code-side reports — completed 2026-03-23
- [x] Modular manuscript framework with converted LaTeX section drafts — completed 2026-03-23
- [x] Statistical comparisons module (10) + report (11): Friedman/Nemenyi/Wilcoxon tests — completed 2026-03-25
- [x] Comprehensive report (`paper/report/`) expanded with full statistical section, design figures
- [x] **Tau-sweep simulation: 12,800 scenarios across τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}** — completed 2026-04-08
- [x] Power metric added (P(reject H₀: τ=0)) and Monte Carlo SEs via delta-method (N_Valid_Est, SE_MSE)
- [x] All reports regenerated with tau sensitivity sections (MSE vs τ, power curves, coverage, rank stability)
- [x] Statistical comparisons updated for tau-sweep; conditional Friedman/Nemenyi across all τ levels
- [x] Pre-sweep deliverables archived: `results/archive/pre_tau_sweep_20260408/`

---

## Open To-Dos / Future Work Roadmap

The core simulation is complete. The following extensions would strengthen the study:

### Simulation Extensions

| Priority | Extension | Description | Notes |
|----------|-----------|-------------|-------|
| High | **Non-oracle MLE** | Re-run MLE without true Spill covariate (`include_spill_covariate = FALSE`) | Toggle already exists in `estimate_tau()`; reveals realistic vs. oracle performance gap |
| Medium | **Heterogeneous population** | Poisson mode with `pop_mode = "heterogeneous"` (unequal cluster sizes) | `pop_mode` parameter exists in `05`; needs `generate_incidence_poisson()` extension |
| Medium | **Grid sensitivity** | Rerun with `grid_dim = 8` (64 clusters) or `grid_dim = 15` (225 clusters) | Tests whether D3/D8 dominance holds at different spatial scales |
| Low | **DIM tau-sweep** | Re-run DIM across τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0} to compute DIM power curves | Low priority — DIM is confirmed naive baseline; MLE results are primary |
| Low | **Heterogeneous beta** | Vary `beta` (incidence coefficient) across incidence modes | Would test robustness to signal strength of incidence covariate |

### Manuscript Development

| Priority | Task | Description |
|----------|------|-------------|
| High | **Simulation results section** | Populate `paper/manuscript/_simulation.qmd` with live R code chunks (figures + tables from `06_visualizations.R` and `08_design_recommendations.R`) |
| High | **Prose review** | Careful editing of `_introduction.qmd` and `_methods.qmd` (converted from LaTeX, not yet polished) |
| Medium | **Application section** | Develop `_application.qmd` — SUD in NC counties context, policy implications, connection to NC DOC work in SpillSpatialDepSim |
| Medium | **Discussion section** | Develop `_discussion.qmd` — limitations (oracle MLE, 50% treated, homogeneous population), future work, recommendations for practitioners |
| Low | **Presentation slides** | Expand `SpatialCRT_IncidenceDesign_Presentation.qmd` scaffold into full conference slides |

### Code Quality / Infrastructure

| Priority | Task | Description |
|----------|------|-------------|
| Low | **`add_mc_ses()` doc fix** | The Roxygen comment says N_Valid_Est=0 produces NA but actually produces Inf; minor documentation inaccuracy |
| Low | **`06_visualizations.R` comment** | Clarify that `_combined_` preference in `load_latest_results()` applies per estimation-mode (not globally) — no behavior change needed, just comment clarity |
| Low | **DIM/MLE join note** | The compare-table in `07_results_summary.Rmd` left-joins DIM (6 designs) onto MLE (8 designs); D7/D8 DIM columns show NA. Caption should mention this asymmetry explicitly. |

---

## Relationship to SpillSpatialDepSim

IncidenceDesign extends the simulation framework from `projects/SpillSpatialDepSim/`.
SpillSpatialDepSim established the SAR model, spillover mechanics, and block stratification
approach in an applied NC DOC context (8–12 districts). IncidenceDesign asks the same
core design question at larger scale (100 clusters) with systematic variation across
incidence modes and formal design strategies.

---

## References

- Mirzaei, A. et al. (2019). Sudden unexpected death rates.
- Gan, W. et al. (2019). County-level Poisson regression for SUD mortality.
- Watson, K. et al. (AHA). Census tract-level spatial analysis of SUD.
- LeSage, J. & Pace, R.K. (2009). *Introduction to Spatial Econometrics*. CRC Press.
- R packages: `spdep` (Bivand et al.), `spatialreg` (Bivand & Piras), `digest`, `ggplot2`, `viridis`
