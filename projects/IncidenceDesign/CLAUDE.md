# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Provides ~80% of the context needed to be immediately productive without
> re-reading all code, result files, or figure PDFs.
> For human-readable narrative documentation see [README.md](README.md).
> For cross-project context, see [../../CLAUDE.md](../../CLAUDE.md).

---

## What This Project Is

A modular simulation study evaluating **8 treatment assignment designs** for Spatial
Cluster Randomized Trials (CRTs) under heterogeneous outcome incidence and spatial
spillover. **tau is swept across {0.8, 1.0, 1.5, 2.0, 3.0}** (direct treatment effect)
for 12,800 total scenarios. Two estimators (DIM, MLE) with 3 incidence generation
modes that are always reported **separately** — never aggregated.

**Application context:** Sudden Unexpected Death (SUD) in NC counties.
Poisson base rate 35/100,000 (Mirzaei et al.).

**Metrics tracked per scenario:** Bias, SD, MSE, Coverage (95% CI), Fail_Rate,
N_Valid_Est (Monte Carlo SE denominator), Power (P(reject H₀: τ=0)).

---

## Research Focus and Framing (IMPORTANT)

**Primary goal: compare the 8 treatment sampling designs** under realistic spatial
conditions. The key question is which design minimizes estimation error and maintains
valid coverage — NOT which estimator is better.

**MLE is the primary estimator.** It is the methodologically appropriate choice
because spatial dependence ($\rho$) and spillover ($\gamma$) are both present in
the outcome DGP. DIM ignores both and produces systematically poor coverage (~72%).

**DIM serves as a proof-of-concept / naive baseline only.** Use DIM results to
sanity-check simulation mechanics and as a secondary comparator. All substantive
design recommendations should be based on MLE results.

**Manuscript prose rule (both `paper/ctj_manuscript/` and `paper/dissertation_chapter/`):
no DIM-vs-MLE comparison, ever.** DIM was always an internal validation baseline,
never a scientific comparator — MLE (`lagsarlm`) is simply presented as *the*
estimator. This is a stricter rule than the codebase/report language above, which is
allowed to mention DIM as a baseline; manuscript prose should not mention DIM at all.
See `feedback_no_dim_comparison_no_shorthand` in the Claude Code memory system for
the reasoning behind this distinction.

---

## File Map

| File | Lines | Role | Key Functions |
|------|------:|------|---------------|
| `00_mathematical_specification.Rmd` | ~495 | Theory & DGP formulas | (rendered to HTML) |
| `01_spatial_setup.R` | 54 | Grid + weight matrices | `build_spatial_grid()`, `get_active_spatial()` |
| `02_incidence_generation.R` | 112 | 3 incidence modes | `generate_incidence()`, `generate_incidence_iid/spatial/poisson()` |
| `03_designs.R` | ~162 | 8 treatment designs | `get_designs()`, `get_design_names()`, `is_design_deterministic()` |
| `04_estimation.R` | 102 | DIM + MLE estimation | `estimate_tau()` |
| `05_run_simulation.R` | ~385 | Main orchestrator | `run_incidence_config()` |
| `06_visualizations.R` | ~800 | Plots + tables | 18 functions — see section below |
| `07_results_summary.Rmd` | ~800 | Rendered results report | Knitted HTML/PDF summary |
| `08_design_recommendations.R` | ~891 | Personalized design recs | `run_recommendation_report()`, `table_scenario_lookup()`, `generate_commentary()` |
| `09_MLE_design_recommendation_report.Rmd` | ~984 | Companion narrative report | Knitted to `results/MLE_design_recommendation_report.pdf` |
| `10_statistical_comparisons.R` | ~1050 | Formal hypothesis tests on design MSE | `run_friedman_test()`, `run_nemenyi_posthoc()`, `run_pairwise_wilcoxon()`, `run_conditional_tests()`, `plot_cd_diagram()`, `plot_mse_boxplot_with_stars()`, `plot_pvalue_heatmap()`, `plot_conditional_cd_diagrams()`, `generate_comparison_report()` |
| `11_statistical_comparisons_report.qmd` | ~543 | Narrative statistical comparisons report | Renders to `results/11_statistical_comparisons_report.{html,pdf}` |
| `12_six_design_statistical_comparisons.R` | ~180 | Re-run Friedman/Nemenyi/Wilcoxon on the 6 manuscript designs + named-label figures | Outputs to `results/six_design_manuscript/` |
| `13_dissertation_results_extract.R` | ~75 | Pulls per-incidence-mode / per-parameter / robustness numbers for the dissertation chapter's fuller Results section | Outputs `results/six_design_manuscript/dissertation_results_extract.txt` |
| `14_manuscript_supplement_figures.R` | ~400 | Fresh 8-design consolidation test + CTJ SI figure suite + reordered/new main-text figures + application table | Outputs `results/eight_design_supplementary/`, `results/six_design_manuscript/si_figures/`, reordered `fig_mse_by_design_6design.pdf`, new `fig_biasvar_6design.pdf`, `application_table_6design.txt` |
| `complete_after_mle.R` | ~685 | Post-MLE script | Runs viz, writes docs, prints stats |
| **paper/report/** | | | |
| `IncidenceSpatialCRT_Report.qmd` | ~1100 | Unified project report (50+ pages, 8-design) | Covers spatial setup → DGP → 8 designs → estimation → simulation → all MLE results → full statistical comparisons (Section 10) → design recommendations |
| **paper/archive_manuscript/** | | | |
| (retired) | — | Byte-identical snapshot of the pre-2026-07-02 `paper/manuscript/` + `paper/section_drafts/` content | Reference/fact-check source only — do not edit in place; both current manuscripts were written fresh, not derived from this |
| **paper/ctj_manuscript/** | | | |
| `CTJ_Manuscript.tex` | ~220 | *Clinical Trials* (SAGE) submission draft, 6 designs, condensed, exactly 6 exhibits (2 tables + 4 figures) | Compiles via `sagej.cls` — see Current State section for `TEXINPUTS`/`BSTINPUTS` setup |
| `Supplementary_Information.tex` | ~400 | SI: reproducibility/seeding (S1), parameter grid (S2), metric formulas (S3), estimation model (S4), full 8-design table (S5), 8-design consolidation justification (S6, only 8-design section), 6-design rankings/sensitivity/robustness/application (S7-S11) | Plain `article` class, S-prefixed numbering, no bibtex needed |
| **paper/dissertation_chapter/** | | | |
| `Dissertation_Chapter.qmd` | ~330 | Longer-form dissertation chapter, 6 designs, full detail, no length ceiling | Quarto → simple double-spaced `article`-class PDF matching Project 1's format |
| **paper/** | | | |
| `spatialCRT.bib` | ~86K | Base bibliography (Zotero export) | Included in `IncidenceDesign_shared.bib` |
| `IncidenceDesign_shared.bib` | ~99K, 79 entries | Single shared bibliography for both manuscripts | `spatialCRT.bib` + SUD/SUDDEN citations + design-theory citations |
| `SAGE_Journal_Template/` | — | `sagej.cls`, `SageH.bst`, `SageV.bst` | CTJ manuscript render dependency |
| `SpatialCRT_IncidenceDesign_Presentation.qmd` | ~26 | Presentation template | Quarto revealjs slides (scaffold) |

---

## Architecture

**Pipeline:** `01 -> 02 -> 03 -> 04`, orchestrated by `05`, visualized by `06`, recommended by `08`.

- `05_run_simulation.R` sources `01`-`04` at runtime via `file.path(script_dir, "0X_*.R")`
- `06_visualizations.R` sources `01`-`03` (needs grid/incidence/design helpers)
- `08_design_recommendations.R` sources `06` (which sources `01`-`03`)
- All detect working directory via `normalizePath(dirname(sys.frame(1)$ofile))` with `tryCatch` fallback to `normalizePath(getwd())`

---

## Data Flow (ASCII)

```
build_spatial_grid(grid_dim=10)
  -> grid_obj: list(coords[100x2], N_clusters=100,
                    nb_rook, nb_queen, W_rook[100x100], W_queen[100x100],
                    listw_rook, listw_queen)
        |
        v
generate_incidence(mode, N=100, n_resamples, W, rho_incidence)
  -> X_matrix: [100 x n_outcome_resamples], values in [0,1]
  -> base_incidence = X_matrix[, 1]   # "historical" incidence for design decisions
        |
        v
get_designs(design_id, n_resamples, N, incidence, nb_list, coords)
  -> Z_matrix: [100 x n_design_resamples], binary {0,1}, ~50% treated
        |
        v
estimate_tau(estimation_mode, Y_sim[100 x n_out], Z[100], spill_term[100],
             X_matrix, active_listw, n_outcome_resamples)
  -> list(estimates[n_out], ci_lower[n_out], ci_upper[n_out])
        |
        v
results data frame — 1 row per scenario:
  Incidence_Mode  | "iid" / "spatial" / "poisson"
  Rho_Incidence   | 0, 0.20, or 0.50
  Neighbor_Type   | "rook" / "queen"
  Design          | "Design 1" .. "Design 8"
  Rho             | 0.00, 0.01, 0.20, 0.50
  Gamma           | 0.5, 0.6, 0.7, 0.8
  Spillover_Type  | "control_only" / "both"
  True_Tau        | 0.8, 1.0, 1.5, 2.0, or 3.0 (swept parameter)
  Mean_Estimate   | mean(tau-hat across iterations)
  Bias            | Mean_Estimate - true_tau
  SD              | sd(tau-hat)
  MSE             | Bias^2 + SD^2
  Coverage        | fraction of CIs containing true_tau
  Fail_Rate       | fraction of MLE iterations that failed to converge
  N_Valid_Est     | count of non-NA estimates (denominator for MC SEs)
  Power           | fraction of CIs excluding zero (P(reject H0: tau=0))
```

---

## Parameter Grid (12,800 total scenarios — tau-sweep)

| Parameter | Values | Levels |
|-----------|--------|:------:|
| `true_tau` (treatment effect) | 0.8, 1.0, 1.5, 2.0, 3.0 | **5** |
| Incidence config | iid (x1) + spatial (x2) + poisson (x2) | 5 |
| `nb_type` | rook, queen | 2 |
| `rho` (outcome spatial autocorrelation) | 0.00, 0.01, 0.20, 0.50 | 4 |
| `gamma` (spillover magnitude) | 0.5, 0.6, 0.7, 0.8 | 4 |
| `spill_type` | control_only, both | 2 |
| `design_id` | 1, 2, 3, 4, 5, 6, 7, 8 | 8 |
| **Total** | 5 x 5 x 2 x 4 x 4 x 2 x 8 | **12,800** |

Scenarios per (tau × incidence config): 512 (= 2x4x4x2x8).
Baseline tau=1.0 slice (2,560 scenarios) reproduces the existing MLE_combined results.

---

## Design Quick Reference

| ID | Name | Deterministic? | ~Treated% | Key Feature |
|----|------|:-:|:-:|-------------|
| 1 | Checkerboard | **Yes** | 50% | Max spatial separation, alternating grid |
| 2 | High Incidence Focus | **Yes** (given X) | 50% | Targets top-50% burden clusters |
| 3 | Saturation Quadrants | No | ~50% (varies) | Random saturation per quadrant (0.2-0.8) |
| 4 | Isolation Buffer | No | ~20-30% | Greedy: no two treated clusters adjacent |
| 5 | 2x2 Blocking | No | 50% | 1:1 randomization within 2x2 spatial blocks |
| 6 | Balanced Quartiles | No | ~50% | Stratified by incidence quartile |
| 7 | Balanced Halves | No | ~50% | 2-strata median split, balanced within each half |
| 8 | Incidence-Guided Saturation Quadrants | No | ~50% | Saturation by quadrant avg incidence rank |

Deterministic designs (1, 2) generate **one** assignment and replicate it across all
`n_design_resamples` — enforced via `is_design_deterministic()` in `03_designs.R`.

---

## Swept Parameters

```r
true_tau_vals <- c(0.8, 1.0, 1.5, 2.0, 3.0)  # Direct treatment effect
```

## Fixed DGP Parameters

```r
beta            <- 1.0          # Incidence coefficient in outcome model
sigma           <- 1.0          # Residual SD
grid_dim        <- 10           # 10x10 = 100 clusters
base_rate       <- 35 / 100000  # Poisson: SUD rate (Mirzaei et al.)
pop_per_cluster <- 1000         # Poisson: equal population per cluster
pop_mode        <- "equal"      # "equal" or "heterogeneous"
include_spill_covariate <- TRUE # Oracle mode: true Spill covariate in MLE
```

**Resample counts:**

| Mode | Design resamples | Outcome resamples | Iterations/scenario |
|------|:---:|:---:|:---:|
| DIM | 25 | 100 | 2,500 |
| MLE | 25 | 10 | 250 |

---

## Current State (as of 2026-07-03)

**CTJ Supplementary Information:** COMPLETE
- `paper/ctj_manuscript/Supplementary_Information.tex` (+ compiled PDF, 14 pages,
  S-prefixed numbering: Tables S1-S12, Figures S1-S7) fulfills the three items
  the main text defers to "online supplementary material": reproducibility/
  seeding/code pipeline (S1), full parameter grid + resample counts (S2),
  metric formulas (S3), estimation model incl. non-oracle variant (S4), the
  full 8-design table (S5), and the complete eight-design comparison (S6) —
  the only section anywhere in either document that presents all 8 designs.
  Sections S7-S11 cover the 6 retained designs only: per-incidence-mode
  rankings, Rho/Gamma/spillover-regime sensitivity (with CD-diagram, heatmap,
  and two-panel figures), tau sensitivity, robustness/win-rate, and a
  6-design application-scale table.
- New `code/14_manuscript_supplement_figures.R` generates all of this: a
  fresh 8-design Friedman/Nemenyi/Wilcoxon run (`choose(8,2)=28` Nemenyi
  denominator, `results/eight_design_supplementary/`) confirming the main
  text's consolidation claims (Nemenyi p=1.0000 for Saturation Quadrants vs.
  Incidence-Guided Saturation Quadrants; MSE 0.1024 vs 0.1095, Wilcoxon
  p=1.34e-6, for Balanced Quartiles vs. Balanced Halves); the SI figure suite
  (clean CD diagrams lifted from `paper/report/IncidenceDesign_ProjectSummary.qmd`'s
  inline code, not the buggy short-label `plot_cd_diagram()` in
  `10_statistical_comparisons.R`); a reordered (best-to-worst) main-text
  `fig_mse_by_design_6design.pdf`; and a new ranked bias-variance
  `fig_biasvar_6design.pdf`.
- `CTJ_Manuscript.tex` updated: swapped in the reordered Figure 1, added the
  new bias-variance Figure 2 (explains *why* designs differ — Checkerboard's
  MSE is nearly all variance from Z/WZ collinearity), added an SI pointer
  sentence. Now exactly 6 exhibits (2 tables + 4 figures, at the CTJ cap);
  body ~2,540 words incl. captions (limit 3,500), abstract 400 words (limit 425).
- Reviewed by a fresh agent against the plan's decisions (6-design scope
  except S6, no DIM-vs-MLE comparison, no design shorthand, numeric accuracy,
  S-numbering) — one gap found and fixed (two generated CD-diagram figures
  were missing from the SI; added to S7 as Figures S1-S2).
- **Not yet done:** the application-table naming wrinkle extends beyond the
  single relabeling the original plan anticipated — the application study
  also calls design 8 "Incidence-Guided Saturation Regions" (not "...
  Quadrants"), both noted with a footnote in SI Section S11; consolidated
  user review pass on the SI still pending.

## Prior State (as of 2026-07-02)

**Manuscripts (CTJ + dissertation chapter):** First full drafts COMPLETE
- Two documents, both written fresh (not edits of the retired `paper/archive_manuscript/`
  files): a short manuscript targeting *Clinical Trials* (SAGE) and a longer-form
  dissertation chapter, sharing a bibliography and verified facts but not prose files.
- `paper/ctj_manuscript/CTJ_Manuscript.tex` — compiles via `sagej.cls` (Sage Vancouver
  style). 7 pages, ~2,255-word body (limit 3,500), ~398-word abstract (limit 425), 6
  exhibits (limit 6). Requires TinyTeX on `PATH` plus
  `TEXINPUTS`/`BSTINPUTS=".:../SAGE_Journal_Template:"` to find `sagej.cls`/`SageV.bst`.
- `paper/dissertation_chapter/Dissertation_Chapter.qmd` — Quarto, renders to a
  simple double-spaced `article`-class PDF (`linestretch: 2`, `unsrt` numbered
  citations) matching the format of the SpillSpatialDepSim (Project 1) manuscript,
  per user direction — not the UNC Graduate School template yet (to be applied in a
  later pass once all dissertation chapters are ready to merge). 28 pages, no length
  ceiling, full theory + full simulation detail + a fuller Results section (per-
  incidence-mode rankings, Rho/Gamma/spillover-regime sensitivity, subgroup
  robustness/win-rate) than the CTJ version. Bibliography must be referenced via the
  `shared-refs.bib` symlink in that directory, not the `../IncidenceDesign_shared.bib`
  path directly — Quarto's pandoc→LaTeX conversion mishandles underscores in bib
  filenames referenced from YAML.
- **Both manuscripts present 6 designs, not 8** (Checkerboard, High Incidence Focus,
  Isolation Buffer, 2x2 Blocking, Balanced Quartiles, Incidence-Guided Saturation
  Quadrants), dropping Saturation Quadrants and Balanced Halves as statistically
  redundant with a retained design (see `code/12_six_design_statistical_comparisons.R`
  and `results/six_design_manuscript/`). **This is unrelated to the older "6-design
  era" mentioned elsewhere in this file** (the pre-Balanced-Halves/pre-Incidence-
  Guided-Saturation-Quadrants phase, before those two designs existed) — the overlap
  in count is coincidental.
- Shared bibliography: `paper/IncidenceDesign_shared.bib` (79 entries) — base
  `paper/spatialCRT.bib` plus 13 new SUD/SUDDEN citations (converted from
  `application/Bibliography_Sudden_Death.html`), 1 Project-1 carryover citation, and
  6 new 2001–2025 design-theory citations, all verified against CrossRef/arXiv.
- Both manuscripts reviewed by fresh agents against the plan's decisions (6 designs,
  no DIM comparison, no design shorthand, accurate Project-1/meeting-notes grounding,
  numeric accuracy against `results/six_design_manuscript/`) — all checks passed
  after two small fixes (a table-layout bug and a swapped rho_X label).
- **Not yet done:** actual journal/advisor submission (waiting on Ashkan Habib's real
  SUDDEN-derived county dataset — current Application-section numbers are an
  explicitly labeled simulated placeholder), UNC dissertation template formatting,
  and a consolidated user review/revision pass on both documents together.

## Prior Simulation State (as of 2026-04-08)

**Tau-sweep simulation:** COMPLETE — 12,800 scenarios, τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}
- Data: `results/sim_data/sim_results_MLE_tau_sweep_combined_20260408_191916.rds`
- Splits: `results/sim_data/sim_results_MLE_tau_sweep_{iid|spatial|poisson}_20260408_191916.rds`
- Stats: 12,800 scenarios | 2,560 per tau level | N_Valid_Est = 250 (all converged) | Fail_Rate = 0.0
- Primary scenario (tau=1.0): D8 MSE=0.079, D3 MSE=0.080 | Worst: D1 MSE=0.802, coverage=55%
- All reports regenerated and integrated with tau sensitivity sections (2026-04-08)

**MLE simulation (baseline tau=1.0):** ARCHIVED — superseded by tau-sweep
- Data preserved: `results/sim_data/sim_results_MLE_combined_20260322_151030.rds`
- Pre-sweep deliverables archived: `results/archive/pre_tau_sweep_20260408/`

**Statistical comparisons (tau-sweep, primary tau=1.0):** UPDATED
- χ² by tau level: τ=0.8: 1052.77, τ=1.0: 1091.90, τ=1.5: 964.05, τ=2.0: 909.95, τ=3.0: 743.41
- All tau levels p < 2.2×10⁻¹⁶ — D3/D8 dominance holds across all effect sizes
- Full report: `results/11_statistical_comparisons_report.{html,pdf}`

**Comprehensive report:** `paper/report/IncidenceSpatialCRT_Report.{html,pdf}` — 50+ pages
- Now includes: Tau Sensitivity section (MSE vs τ, power curves, coverage, rank stability)
- Monte Carlo SEs section (N_Valid_Est, SE_MSE per design at primary tau=1.0)

**DIM simulation:** COMPLETE for prior 6-design sweep (naive baseline only)
- Data: `results/sim_data/sim_results_DIM_combined_20260304_195321.rds`
- Coverage systematically ~72% — use MLE for all substantive analysis.

**Results directory layout:**
```
results/
  MLE_tau_sweep_design_recommendations.pdf    # PRIMARY figures/tables PDF
  MLE_tau_sweep_incidence_overview.pdf        # Incidence overview
  11_statistical_comparisons_report.{html,pdf} # Formal hypothesis testing + tau-strata
  00_mathematical_specification.pdf
  07_results_summary.pdf
  09_MLE_design_recommendation_report.{html,pdf}
  sim_data/          # All .rds files (load_latest_results() auto-detects this)
  mle_per_config/    # MLE per-config PDFs (tau_sweep_* + tau_sweep_*_tau_sensitivity)
  figures/           # design_samples_8panel + design_samples_option1_overlays
  six_design_manuscript/  # 6-design stats re-run + named-label figures (both manuscripts pull from here)
  archive/
    pre_tau_sweep_20260408/  # Archived pre-sweep deliverables
paper/
  report/
    IncidenceSpatialCRT_Report.{qmd,html,pdf}  # Unified report (now 50+ pages with tau section)
  archive_manuscript/  # Retired paper/manuscript/ + paper/section_drafts/ content (reference only,
                        # byte-identical to last commit before this location was retired 2026-07-02)
  ctj_manuscript/
    CTJ_Manuscript.{tex,pdf}  # Clinical Trials (SAGE) submission draft
    figures/
  dissertation_chapter/
    Dissertation_Chapter.{qmd,pdf}  # Longer-form dissertation chapter draft
    figures/
    shared-refs.bib -> ../IncidenceDesign_shared.bib  # symlink, see note above
  IncidenceDesign_shared.bib  # Single bibliography shared by both manuscripts
  SAGE_Journal_Template/      # sagej.cls, SageH.bst, SageV.bst (CTJ render dependency)
```

**Git:** Original work on `claude/gallant-buck` → `main` 2026-03-05.
Reorganized into `projects/IncidenceDesign/` on `claude/dreamy-wiles`.
Statistical comparisons + report expansions on `claude/elated-lederberg`.
Tau-sweep + MC SEs + Power implementation on `main` 2026-04-08.
Tau-sweep results integrated, all reports regenerated on `main` 2026-04-08.
CTJ + dissertation-chapter manuscript first drafts on `main` 2026-07-02.
CTJ Supplementary Information + ranked figures on `main` 2026-07-03.

---

## Critical Invariants (DO NOT Violate)

1. **Incidence generated ONCE per `(mode, rho_X)` config** — NOT per `(nb_type, rho)` scenario.
   `X_matrix` is fixed for a given config; only the outcome DGP varies with `rho`.

2. **Results are NEVER aggregated across incidence modes.** Each of iid, spatial, poisson
   is always reported separately. See `split_by_incidence_config()` in `06_visualizations.R`.

3. **Deterministic designs (1, 2) generate 1 assignment and replicate** for all design
   resamples. Never re-draw them `n_design_resamples` times.

4. **Per-scenario seed** = `digest::digest2int(paste(inc_mode, rho_x, nb_type, rho, gamma, spill_type, d_id, sep="|"))`.
   Never use a single `set.seed()` for the whole run.

5. **`base_incidence = X_matrix[, 1]`** — only the first column is used for treatment
   design decisions (the "historically observed" incidence). Do not pass the full matrix.

---

## Bugs Previously Encountered & Fixed

| Bug | Symptom | Fix |
|-----|---------|-----|
| Unicode rho in PDF | `pdf()` conversion failure on `\u03c1` | Use literal text `"rho"` |
| Pairwise dominance NaN | `make.names("Design 1")` -> `"Design.1"` mismatched pivot column | Use `wide[[designs[i]]]` direct bracket access |
| Eta-squared > 1.0 | `var(group_means)/var(total)` incorrect | Use `SS_between/SS_total` where `SS_between = sum(n_j*(mean_j - grand_mean)^2)` |
| `I` shadows `base::I` | Subtle namespace errors | Renamed to `I_mat <- diag(N_clusters)` |
| Degenerate Z (all 0 or 1) | `estimate_tau()` crash | Early-exit guard returning all `NA` |
| R sprintf 10k char limit | `complete_after_mle.R` failed writing docs | Write docs directly via Claude Code `Write` tool instead |

---

## Common User Requests -> Code Actions

| Request | Action |
|---------|--------|
| Run DIM simulation | Set `estimation_mode <- "DIM"` in `05` line 40, `Rscript 05_run_simulation.R` |
| Run MLE simulation | Set `estimation_mode <- "MLE"` in `05` line 40, `Rscript 05_run_simulation.R` |
| Generate MLE plots only | `source("06_visualizations.R"); run_all_visualizations(estimation_mode="MLE_combined")` |
| One incidence config | `r <- load_latest_results(estimation_mode="MLE_combined"); cfg <- split_by_incidence_config(r); run_standard_tables(cfg[["iid Uniform"]])` |
| Custom stratified table | `table_stratified(results, c("Design", "Rho"), "iid Uniform")` |
| Add a new design | Add `case` to `get_designs()` in `03`, update `get_design_names()`, add ID to `design_ids` in `05` |
| Change parameter sweep | Edit `*_vals` vectors in `05` CONFIGURABLE PARAMETERS section (lines ~50-80) |
| Run parallel | Set `n_cores > 1` in `05` line 44 |
| Compare DIM vs MLE | Load both result sets, join on scenario keys, compare Coverage and Bias columns |
| Non-oracle MLE | Set `include_spill_covariate = FALSE` in `estimate_tau()` call in `05` |
| Full recommendation report | `source("08_design_recommendations.R"); run_recommendation_report(estimation_mode="MLE_combined")` |
| Rankings per incidence mode | `table_incidence_rankings(results)` or `plot_incidence_rankings(results)` |
| Rankings for one parameter | `table_marginal_rankings(results, "Rho", "iid Uniform")` |
| Rank trajectory plot | `plot_rank_trajectories(results, "Gamma", "iid Uniform")` |
| Best design heatmap | `plot_best_design_heatmap(results, "Rho", "Gamma", "iid Uniform")` |
| Specific scenario lookup | `table_scenario_lookup(results, rho=0.2, gamma=0.7, spill_type="both", nb_type="queen")` |
| Validate recommendations | `validate_recommendations()` then `validate_no_side_effects()` |

---

## Visualization Functions (`06_visualizations.R`)

**Helpers:**
- `inc_config_label(inc_mode, rho_x)` — human-readable config label string
- `split_by_incidence_config(results)` — splits combined df into named list of per-config dfs
- `load_latest_results(results_dir, estimation_mode)` — loads most recently modified .rds; auto-detects `sim_data/` subdirectory if present

**Plots** (each accepts `results` df + `inc_label` string):
- `plot_mse_by_neighbor()` — boxplot of MSE by design, faceted by rook/queen
- `plot_mse_heatmap()` — tile heatmap: rho x design, cell = avg MSE
- `plot_mse_per_design()` — per-design MSE line plots across gamma, one panel per design
- `plot_master_comparison()` — bar chart sorted by MSE for a given neighbor type
- `plot_coverage_by_design()` — boxplot of coverage by design (red dashed line at 0.95)
- `plot_incidence_heatmaps()` — side-by-side grid heatmaps of iid/spatial/poisson incidence
- `plot_incidence_distributions()` — faceted density plots comparing incidence distributions
- `plot_design_ranks()` — stacked bar of rank-frequency ("win rates") per design
- `plot_bias_variance()` — stacked bar of Bias-squared + Variance components per design
- `plot_coverage_mse_tradeoff()` — scatter: each point = 1 scenario, x=MSE, y=Coverage

**Tables** (each prints to console, returns invisible):
- `table_design_ranks()` — win-rate frequency table
- `table_robustness()` — best/Q25/median/Q75/worst MSE per design
- `table_pairwise_dominance()` — % of scenarios where row-design beats col-design
- `table_sensitivity()` — eta-squared for Rho, Gamma, Spillover_Type, Neighbor_Type per design
- `table_stratified(results, group_by_vars, inc_label)` — flexible grouping by any combination
- `table_comprehensive()` — full scenario-level detail (first 30 rows by default)

**Runners:**
- `run_standard_tables(results, inc_label)` — all 10 standard tables for one config
- `run_all_visualizations(results, results_dir, estimation_mode, output_pdf=TRUE)` — full pipeline

---

## Design Recommendation Functions (`08_design_recommendations.R`)

Sources `06_visualizations.R`. Answers three personalization questions.
**Always run on MLE results only** (`estimation_mode = "MLE_combined"`).

**Core utility:**
- `rank_designs_by_group(results, group_vars, metric)` — rank designs by avg MSE within groups

**Q1 — Per incidence mode:**
- `table_incidence_rankings(results)` — console table: design rank + MSE + Coverage per config
- `plot_incidence_rankings(results)` — faceted bar chart: avg MSE per design per config
- `plot_incidence_coverage(results)` — faceted bar chart: avg Coverage per design per config

**Q2 — Per parameter level (marginal):**
- `table_marginal_rankings(results, param, inc_label)` — wide table: param levels × designs, cells = MSE (#Rank)
- `plot_rank_trajectories(results, param, inc_label)` — line plot of rank vs parameter level (1 = best at top)
- `plot_conditional_mse(results, param, inc_label)` — faceted boxplot: MSE distribution per design at each param level

**Q3 — Per parameter combination:**
- `plot_best_design_heatmap(results, row_param, col_param, inc_label)` — tile heatmap: winning design per (row, col) combo
- `table_scenario_lookup(results, rho, gamma, spill_type, nb_type, inc_label)` — filter to specific params, rank designs, print recommendation

**Summary:**
- `generate_commentary(results, inc_label)` — 6-finding programmatic narrative: winner, dominance %, stability, coverage, sensitivity, recommendation
- `run_recommendation_report(results, estimation_mode, output_pdf)` — master orchestrator → `results/MLE_combined_design_recommendations.pdf`

**Validation:**
- `validate_recommendations(results)` — 8 unit tests for new functions
- `validate_no_side_effects(results)` — 6 integration tests verifying existing modules unaffected

---

## Extensions & Future Work Roadmap

### Completed (as of 2026-04-08)
- **Tau-sweep:** COMPLETE — 12,800 scenarios, τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}, all reports integrated
- **True_Tau dimension in Quarto reports:** COMPLETE — all 6 report files updated
- **Power metric:** COMPLETE — P(reject H₀: τ=0) tracked per scenario
- **Monte Carlo SEs:** COMPLETE — N_Valid_Est, SE_MSE columns; delta-method SE_MSE = √(2SD⁴ + 4Bias²SD²)/√N

### Simulation Extensions (open)

| Priority | Extension | Implementation Note |
|----------|-----------|---------------------|
| High | **Non-oracle MLE** (`include_spill_covariate = FALSE`) | Toggle already in `estimate_tau()`; run `05_run_simulation.R` with `include_spill_covariate = FALSE` |
| Medium | **Heterogeneous population** Poisson mode | `pop_mode = "heterogeneous"` in `05`; extend `generate_incidence_poisson()` for unequal cluster sizes |
| Medium | **Grid sensitivity** | Change `grid_dim` to 8 or 15 in `05`; tests D3/D8 dominance at different spatial scales |
| Low | **DIM tau-sweep** | Re-run DIM across all τ levels for power curve comparison (DIM is confirmed naive baseline) |
| Low | **Heterogeneous beta** | Vary `beta` coefficient across simulation configs |

### Manuscript Development (open)

First full drafts of both manuscripts, plus the CTJ Supplementary Information,
are complete (see Current State above: `paper/ctj_manuscript/`,
`paper/dissertation_chapter/`). Remaining work:

| Priority | Task | File |
|----------|------|------|
| High | Consolidated user review/revision pass on all three documents together (CTJ main text, CTJ SI, dissertation chapter) | All |
| High | Replace placeholder Application-section numbers with real SUDDEN-derived data once finalized | `paper/ctj_manuscript/CTJ_Manuscript.tex`, `paper/ctj_manuscript/Supplementary_Information.tex` (Section S11), `paper/dissertation_chapter/Dissertation_Chapter.qmd` |
| Medium | Apply UNC Graduate School dissertation template/formatting | `paper/dissertation_chapter/Dissertation_Chapter.qmd` |
| Medium | Fix the pre-existing `plot_cd_diagram()` label-collision bug (designs with adjacent ranks overlap regardless of image width) — currently worked around by omitting the CD diagram from the dissertation chapter's inline exhibits | `code/10_statistical_comparisons.R` |
| Low | Expand presentation scaffold into full conference slides | `SpatialCRT_IncidenceDesign_Presentation.qmd` |

### Known Minor Issues (low-priority cleanup)

- `add_mc_ses()` Roxygen doc: says N_Valid_Est=0 → NA but actually produces Inf
- `load_latest_results()` comment (line ~87 of `06_visualizations.R`): clarify `_combined_` preference applies per estimation-mode, not globally
- `07_results_summary.Rmd` compare-table caption: should note that DIM only ran 6 designs (D7/D8 NAs are expected)
- `IncidenceSpatialCRT_Report.qmd` caption/text alignment: MC SEs table uses tau=1.0 slice (2,560 rows), not full 12,800 — prose now correctly clarifies this distinction

---

## Relationship to SpillSpatialDepSim

IncidenceDesign extends the applied simulation framework from `projects/SpillSpatialDepSim/`.
SpillSpatialDepSim used a small grid (8–12 districts) with SAR estimation to evaluate
block vs. random assignment in an NC DOC applied context. IncidenceDesign asks the same
core question but at larger scale (100 clusters) with systematic design variation and
heterogeneous incidence modes.

Cross-reference SpillSpatialDepSim results:
```r
here::here("projects", "SpillSpatialDepSim", "results")
# or relative: ../../SpillSpatialDepSim/results/
```
