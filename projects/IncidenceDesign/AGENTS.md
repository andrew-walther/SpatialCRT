# AGENTS.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Provides ~80% of the context needed to be immediately productive without
> re-reading all code, result files, or figure PDFs.
> For human-readable narrative documentation see [README.md](README.md).
> For cross-project context, see [../../AGENTS.md](../../AGENTS.md).
> Progress checkpoints (dated status logs) and future work live in
> [ROADMAP.md](ROADMAP.md), not here. Keep this file to reference material and
> standing rules.

---

## Standing Author Rules

- **Reference rule (author-confirmed 2026-09-27):** dissertation prose may refer directly
  to "Chapter 2". CTJ must use self-contained wording (e.g., "In previous work...")
  with a citation to the Project 1 BMC Medical Research Methodology paper, never
  "Chapter 2" or "Project 1" as reader-facing references. Describe findings within
  their studied conditions; verify final bibliographic details before adding the citation.
- **Project 1 is final.** Treat `projects/SpillSpatialDepSim/` as read-only. Changes
  reach `bios-dissertation` only through the Chapter 3 sync (below) unless the author
  asks otherwise.
- **Allocation-risk numbers** (selected-surface/corner pilot, `results/allocation_risk/`)
  are never substituted for full-study numbers.
- **Where to write things down:** key findings go in `docs/plans/` (open items:
  `docs/plans/TODO.md`); dated progress checkpoints go in the progress log of
  `ROADMAP.md`. Restricted SUD data details: `application/AGENTS.md`.

---

## What This Project Is

A modular simulation study evaluating **8 treatment assignment designs** for Spatial
Cluster Randomized Trials (CRTs) under heterogeneous outcome incidence and spatial
spillover. **tau is swept across {0.8, 1.0, 1.5, 2.0, 3.0}** (direct treatment effect)
for 12,800 total scenarios. Since the 2026-09 revision, two ML estimators are fit to
the same data: oracle (primary) and non-oracle (sensitivity). DIM is a pre-revision
baseline only. The 3 incidence generation modes and the 2 neighbor types are reported
separately, with queen primary.

**Application context:** Sudden Unexpected Death (SUD) in NC counties.
Poisson base rate 35/100,000 (Mirzaei et al.).

**Metrics tracked per scenario:** Bias, SD, MSE, Coverage (95% CI), Fail_Rate,
N_Valid_Est, Power (P(reject H₀: τ=0)), surface-level MC SEs (SE_Bias/SE_MSE/SE_Coverage/
SE_Power), Mean_Treated, and flags N_Aliased / N_Warn / Z_WZ_rank_deficient.

---

## Research Focus and Framing (IMPORTANT)

**Primary goal: compare the treatment sampling designs** (6 in the manuscripts, 8 in the
reports) under realistic spatial conditions. The key question is which designs estimate
the treatment effect well under heterogeneous incidence (low MSE, valid coverage). It is
NOT which estimator is better, and NOT whether knowing prior incidence is necessary
(user framing, 2026-09-24).

**The oracle ML spatial-lag estimator is primary** (Y ~ Z + Spill + X, true spillover
covariate). It's the methodologically appropriate choice because spatial dependence ($\rho$)
and spillover ($\gamma$) are both present in the outcome DGP. The non-oracle ML fit
(Y ~ Z + X) is a sensitivity analysis. DIM ignores both and produced ~72% coverage in the
(pre-revision) April runs.

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
| `06_visualizations.R` | ~800 | Plots + tables | 18 plot/table/runner functions (see the script) |
| `07_results_summary.Rmd` | ~800 | Rendered results report | Knitted HTML/PDF summary |
| `08_design_recommendations.R` | ~891 | Personalized design recs | `run_recommendation_report()`, `table_scenario_lookup()`, `generate_commentary()` |
| `09_MLE_design_recommendation_report.Rmd` | ~984 | Companion narrative report | Knitted to `results/MLE_design_recommendation_report.pdf` |
| `10_statistical_comparisons.R` | ~1050 | Formal hypothesis tests on design MSE | `run_friedman_test()`, `run_nemenyi_posthoc()`, `run_pairwise_wilcoxon()`, `run_conditional_tests()`, `plot_cd_diagram()`, `plot_mse_boxplot_with_stars()`, `plot_pvalue_heatmap()`, `plot_conditional_cd_diagrams()`, `generate_comparison_report()` |
| `11_statistical_comparisons_report.qmd` | ~543 | Narrative statistical comparisons report | Renders to `results/11_statistical_comparisons_report.{html,pdf}` |
| `12_six_design_statistical_comparisons.R` | ~180 | Re-run Friedman/Nemenyi/Wilcoxon on the 6 manuscript designs + named-label figures | Outputs to `results/six_design_manuscript/` |
| `13_dissertation_results_extract.R` | ~75 | Pulls per-incidence-mode / per-parameter / robustness numbers for the dissertation chapter's fuller Results section | Outputs `results/six_design_manuscript/dissertation_results_extract.txt` |
| `14_manuscript_supplement_figures.R` | ~400 | Fresh 8-design consolidation test + CTJ SI figure suite + reordered/new main-text figures + application table | Outputs `results/eight_design_supplementary/`, `results/six_design_manuscript/si_figures/`, reordered `fig_mse_by_design_6design.pdf`, new `fig_biasvar_6design.pdf`, `application_table_6design.txt` |
| `complete_after_mle.R` | ~685 | Post-MLE script | Runs viz, writes docs, prints stats |
| **application/** | | | |
| `application/code/run_sud_aggregation.R` | ~65 | Entry script: real SUD data → county + 58-cluster rates per 100k (by year, total, average) → `application/data/derived/` | Sources the three `sud_*.R` files below |
| `application/code/sud_load_data.R` | ~105 | Load and join numerator + denominator; complete-panel checks | `load_sud_county_data()` |
| `application/code/sud_incidence.R` | ~110 | Rate formulas; county → college aggregation | `compute_incidence()`, `aggregate_sud_to_clusters()` |
| `application/code/sud_reconcile.R` | ~100 | QC report vs Habib (2026) and the 23,523 counts | `reconcile_sud_counts()` |
| `application/code/build_cluster_weights.R` | ~70 | Queen + rook contiguity (nb, binary A, row-standardized W) for the 58 clusters → tracked `application/data/nc_cluster_weights.rds` | Reuses `build_nc_application_clusters()`, `build_application_spatial_weights()` |
| `application/tests/test_application_data.R` | ~130 | 23 checks: rate formulas (fixture), bad input, real data vs Habib (2026), weights | Real-data part skips if the gitignored sources are absent |
| **paper/report/** | | | |
| `IncidenceSpatialCRT_Report.qmd` | ~1100 | Unified project report (50+ pages, 8-design) | Covers spatial setup → DGP → 8 designs → estimation → simulation → all MLE results → full statistical comparisons (Section 10) → design recommendations |
| **paper/archive_manuscript/** | | | |
| (retired) | — | Byte-identical snapshot of the pre-2026-07-02 `paper/manuscript/` + `paper/section_drafts/` content | Reference/fact-check source only — do not edit in place; both current manuscripts were written fresh, not derived from this |
| **paper/ctj_manuscript/** | | | |
| `CTJ_Manuscript.tex` | ~220 | *Clinical Trials* (SAGE) submission draft, 6 designs, condensed, exactly 6 exhibits (2 tables + 4 figures) | Compiles via `sagej.cls` — see Build gotchas for `TEXINPUTS`/`BSTINPUTS` setup |
| `Supplementary_Information.tex` | ~400 | SI: reproducibility/seeding (S1), parameter grid (S2), metric formulas (S3), estimation model (S4), full 8-design table (S5), 8-design consolidation justification (S6, only 8-design section), 6-design rankings/sensitivity/robustness/application (S7-S11) | Plain `article` class, S-prefixed numbering, no bibtex needed |
| **paper/dissertation_chapter/** | | | |
| `Dissertation_Chapter.qmd` | ~330 | Longer-form dissertation chapter, 6 designs, full detail, no length ceiling | Quarto → simple double-spaced `article`-class PDF matching Project 1's format |
| **paper/** | | | |
| `SpatialCRT_IncidenceDesign.bib` | ~99K, 81 entries | Single shared bibliography for both manuscripts (renamed from `IncidenceDesign_shared.bib` 2026-09-24; `baird`/`leung` entries synced to `bios-dissertation/prelim/references.bib`; `habib_temporal_2026` added 2026-09-24 = Ashkan Habib's unpublished BIOS master's paper, source of the SUD county data) | Zotero base + SUD/SUDDEN citations + design-theory citations |
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
generate_surfaces(cfg)            # 05; X_k keyed by ("X", mode, rho_X, k)
  -> X: [100 x K=10], values in [0,1]; shared by every block of the config
        |
        v
draw_assignments(cfg, nb, rho, gamma, regime, d, X)   # 05 -> get_designs() per surface
  -> Zm: [100 x K*J = 250]; columns (k-1)*25+1..k*25 drawn from X[, k]  (M1)
        |
        v
Y_kj = (I - rho W)^{-1} (tau Z_kj + S(Z_kj) + beta X_k + eps_kj)   # eps: 250 columns per block
        |
        v
fit_tau_models(y, Z, spill, X_k, setup)   # 04; oracle + non-oracle via fit_sar_lag()
  -> list(oracle = list(tau, se, warns), nonoracle = ...)
        |
        v
summarize_fits()                          # 05; metrics + surface-level MC SEs
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
  Estimator       | "oracle" / "nonoracle" (one file per estimator)
  Mean_Estimate   | mean(tau-hat over the 250 fits)
  Bias            | Mean_Estimate - true_tau
  SD              | sd(tau-hat)
  MSE             | mean((tau-hat - true_tau)^2)
  Coverage        | fraction of CIs containing true_tau
  Fail_Rate       | fraction of fits with non-finite tau-hat or SE
  N_Valid_Est     | count of valid fits (250)
  Power           | fraction of CIs with lower bound > 0 (one-sided)
  SE_Bias/SE_MSE/SE_Coverage/SE_Power | sd of 10 surface means / sqrt(10)  (t9 intervals)
  N_Surfaces      | surfaces with valid fits (10)
  Mean_Treated    | mean sum(Z)
  N_Aliased       | fits where the engine dropped an aliased column
  N_Warn          | all other warnings (logged in warnings_*.csv)
  Z_WZ_rank_deficient | rank([1, Z, WZ]) < 3 for any draw
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

Scenarios per (tau × incidence config): 512 (= 2x4x4x2x8). Each scenario has 250 fits per
estimator (K = 10 surfaces × J = 25 design draws). The `pilot` profile of 05 runs τ = 1,
ρ ∈ {0, 0.5}, γ ∈ {0.5, 0.8}.

---

## Design Quick Reference

| ID | Name | Deterministic? | ~Treated% | Key Feature |
|----|------|:-:|:-:|-------------|
| 1 | Checkerboard | **Yes** | 50 | Alternating (x + y) mod 2 grid (maximal interspersion; under rook WZ = 1 − Z) |
| 2 | High Incidence Focus | No (random ties) | exactly 50 | Treat the 50 highest-incidence clusters |
| 3 | Saturation Quadrants | No | 50 | Random saturation {0.2, 0.4, 0.6, 0.8} per 5×5 quadrant |
| 4 | Isolation Buffer | No | rook ≈ 38, queen ≈ 22 | Greedy random maximal independent set (rarely the exact checkerboard under rook) |
| 5 | 2x2 Blocking | No | 50 | 2 of 4 within each 2×2 block |
| 6 | Balanced Quartiles | No | exactly 50 | Equal rank quartiles; 12/12/13/13 treated in random order |
| 7 | Balanced Halves | No | 50 | Exact rank halves, 25 treated in each |
| 8 | Incidence-Guided Saturation Quadrants | No | 50 | Saturations {0.8…0.2} by rank of quadrant mean incidence |

Only Checkerboard is deterministic (`is_design_deterministic()` in `03_designs.R`). Ties in
incidence are broken at random for every draw (`random_tie_rank()`).

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
pop_per_cluster <- 100000       # Poisson: equal population per cluster (M2; ≈35 deaths/yr)
pop_mode        <- "equal"      # "equal" or "heterogeneous"
n_surfaces      <- 10           # K incidence surfaces per config
n_design_draw   <- 25           # J design draws per surface; 250 fits per scenario
```
Both the oracle and the non-oracle model are fit to every simulated outcome. The pre-revision
DIM baseline used 25 design × 100 outcome resamples, and wasn't re-run.

---

## Current State (summary; dated detail in ROADMAP.md "Progress log")

- **Simulation:** revised and re-run 2026-09-24 (12,800 scenarios × 250 fits × oracle and
  non-oracle ML). These are the current numbers; the April 2026 numbers are superseded
  (`results/archive/pre_revision_20260924/`).
- **Dissertation chapter:** rewritten on the revised run (Phase D) and synced to
  `bios-dissertation` as Chapter 3. Inline `NEEDS-AUTHOR-CONFIRMATION` items remain in the
  Application section.
- **Allocation-risk pilot:** complete; `docs/plans/allocation_risk_findings_2026-09-27.md`.
- **CTJ manuscript + SI:** still carry the April numbers until manuscript step 5
  (`docs/plans/manuscript-unification-plan.md`).
- **Next:** Simple Random Sampling benchmark design (before the CTJ); the design comparison
  on the real SUD data (aggregated to the 58 clusters 2026-09-25); CTJ derivation.

## Dissertation Chapter 3 Sync

- `paper/dissertation_chapter/Dissertation_Chapter.qmd` is the only file anyone edits
  (author decision 2026-09-25).
- `paper/dissertation_chapter/tools/sync_to_prelim.sh` generates
  `bios-dissertation/prelim/project-proposals/project2-incidence/draft/project2-incidence-draft.qmd`
  (CHAPTER 3 / APPENDIX B). It checks every citekey against the master bib, failing loud,
  renders with RStudio's Quarto 1.10, and commits only `project2-incidence/` paths there.
  It never pushes and skips PDF-only re-renders.
- A `post-commit` hook (install with `tools/install_hooks.sh`) runs it when a commit touches
  the chapter, its figures, `tools/` or the bib. Its log is `tools/sync.log` (gitignored).
- The prelim YAML header is `tools/prelim_header.yml` (`colorlinks: false`, so links print
  black). Cross-references are `\label`/`\ref`, so numbers read 3.x / B.x in the prelim.

## Build gotchas

- **CTJ manuscript:** needs TinyTeX on `PATH` plus `TEXINPUTS`/`BSTINPUTS=".:../SAGE_Journal_Template:"` to find `sagej.cls`/`SageV.bst`.
- **Dissertation chapter bibliography:** reference it via the `shared-refs.bib` symlink in `paper/dissertation_chapter/`, not `../SpatialCRT_IncidenceDesign.bib` directly — Quarto's pandoc→LaTeX conversion mishandles underscores in bib filenames referenced from YAML.
- Dated progress logs moved to `ROADMAP.md` ("Progress log") on 2026-10-01; logs from
  2026-04-08 to 2026-09-04 are in `git log -p AGENTS.md`.

---

## Critical Invariants (DO NOT Violate)

Revised 2026-09-24 for the simulation revision (spec: `docs/plans/simulation-revision-spec.md`).

1. **Incidence surfaces are generated per `(mode, rho_X, k)`, k = 1..10, and shared by every
   block of that config** — NOT regenerated per `(nb_type, rho, gamma, regime)` block. Only the
   outcome DGP, the design draws and the noise vary across blocks.

2. **Never pool silently.** Incidence modes and neighbor types are reported separately:
   queen is primary, rook is the sensitivity / Chapter 2 continuity case. Pooled numbers
   appear only when labeled "pooled", with splits wherever a conclusion changes. See
   `split_by_incidence_config()` in 06 and the queen/rook/pooled slices in 12-14.

3. **Checkerboard is the only deterministic design.** It uses one assignment for every
   draw, but each of its 250 fits still gets its own noise column. Every other design,
   including High Incidence Focus, is re-drawn for each of the J = 25 draws per surface.

4. **Every random draw is key-seeded** with `set_seed_key()` (05): X by `("X", mode, rho_X, k)`,
   Z by `("Z", mode, rho_X, nb, rho, gamma, regime, k, d)`, eps by
   `("eps", mode, rho_X, nb, rho, gamma, regime)`. Numeric fields use `sprintf("%.2f")`, and RNG
   kinds are set explicitly. Never draw outside `set_seed_key()` in the runner.

5. **Matched surfaces (M1):** designs for surface k are drawn from `X[, k]`, and the outcome
   and the analysis use the same `X[, k]`. (This reverses the old `base_incidence = X_matrix[, 1]` rule.)

6. **Checkpoints are manifest-checked** (`results/checkpoints/rev_2026-09/<profile>/manifest.rds`:
   parameter hash, code hashes for 01-05, package versions, BLAS). If anything changed, the
   runner refuses to load. Never copy or reuse checkpoints across code versions.

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
| Design/outcome incidence mismatch (fixed 2026-09) | Designs used `X[,1]`; outcome replicate k used `X[,k]` (cor ≈ 0) | Matched surfaces (M1) |
| Poisson ties (fixed 2026-09) | 1,000 per cluster → ~5 distinct X; High Incidence Focus treated 14–50; `ntile()` broke ties by grid position | 100,000 per cluster; `random_tie_rank()` per draw; exact N/2 |
| Unreproducible seeds (fixed 2026-09) | X and eps drawn before the per-scenario seed; forked workers reseeded | Key-based `set_seed_key()` for every draw |
| Overstated MC precision (fixed 2026-09) | Deterministic designs copied 25× against 10 noise draws | One noise column per fit; SEs from surface means |
| Silent aliasing (fixed 2026-09) | Checkerboard × rook: WZ = 1 − Z; `lagsarlm` dropped Spill, and `suppressWarnings` hid it | `withCallingHandlers`; `N_Aliased`, `Z_WZ_rank_deficient` |
| Stale checkpoints (fixed 2026-09) | Old checkpoints loaded silently on re-run | Manifest check refuses mismatches |
| `lagsarlm` slowness | ~76 ms/fit, ~94% of it an unconditional `gc()` on exit | `fit_sar_lag()` (validated, ~0.6 ms/fit) |
| Rank-trajectory axis (fixed 2026-09) | `plot_rank_trajectories()` y-limits hard-coded to 6 silently dropped ranks 7–8 | Limits from `max(Rank)` |
| Commentary pivot (fixed 2026-09) | `generate_commentary()` pivot lacked `True_Tau` → list-columns on tau-sweep files | `True_Tau` added to pivot keys |
| Naive MC SEs (fixed 2026-09) | `add_mc_ses()` overwrote stored surface-level SEs | Returns stored SEs when present |

---

## Common User Requests -> Code Actions

| Request | Action |
|---------|--------|
| Run the simulation (oracle + non-oracle) | From `code/`: `VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript 05_run_simulation.R full` (or `pilot`); ~14 min on 10 workers |
| Run the tests | `VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript tests/test_simulation_revision.R` (`SKIP_EQUIVALENCE=1` skips the ~7 min lagsarlm check); then `Rscript tests/verify_full_run_rev_2026-09.R` |
| Run DIM simulation | Not supported by the revised 05 (DIM is pre-revision only; legacy path `estimate_tau()` in 04) |
| Generate MLE plots only | `source("06_visualizations.R"); run_all_visualizations(results_dir = normalizePath("../results"), estimation_mode = "MLE_tau_sweep")` |
| One incidence config | `r <- load_latest_results(estimation_mode="MLE_tau_sweep"); cfg <- split_by_incidence_config(r); run_standard_tables(cfg[["iid Uniform"]])` |
| Custom stratified table | `table_stratified(results, c("Design", "Rho"), "iid Uniform")` |
| Add a new design | Add `case` to `get_designs()` in `03`, update `get_design_names()`, add ID to `design_ids` in `05` |
| Change parameter sweep | Edit `*_vals` vectors in `05` CONFIGURABLE PARAMETERS section (the checkpoint manifest then forces a fresh run) |
| Run parallel | `n_cores` in `05` (default 10; forked `mclapply` over 40 units) |
| Non-oracle results | `load_latest_results(estimation_mode = "MLEnonoracle_tau_sweep")`; 6-design tests: `Rscript 12_six_design_statistical_comparisons.R MLEnonoracle` |
| Full recommendation report | `source("08_design_recommendations.R"); run_recommendation_report(estimation_mode="MLE_tau_sweep")` |
| Rankings per incidence mode | `table_incidence_rankings(results)` or `plot_incidence_rankings(results)` |
| Rankings for one parameter | `table_marginal_rankings(results, "Rho", "iid Uniform")` |
| Rank trajectory plot | `plot_rank_trajectories(results, "Gamma", "iid Uniform")` |
| Best design heatmap | `plot_best_design_heatmap(results, "Rho", "Gamma", "iid Uniform")` |
| Specific scenario lookup | `table_scenario_lookup(results, rho=0.2, gamma=0.7, spill_type="both", nb_type="queen")` |
| Validate recommendations | `validate_recommendations()` then `validate_no_side_effects()` |

---

## Known Minor Issues (low-priority cleanup)

- `02_incidence_generation.R`: `generate_incidence_poisson()` still defaults `pop_per_cluster = 1000` (and its roxygen says so). 05 always passes 100,000, so results aren't affected. It was left unchanged after the full run so the run manifest's code hash still matches.
- `DESIGN_FULL_NAMES` in `00_design_names.R` (used only by the two design-sample figures) labels Design 1 "Block Stratified Sampling (Checkerboard)" (user decision 2026-09-24: Chapter 2's name plus this study's). The manuscripts, 03 and 12–14 say "Checkerboard", and `DESIGN_SHORT_NAMES` still says "Block Stratified". This is not a bug to fix: the author uses "Block Stratified Sampling" and "Checkerboard" interchangeably, and the chapter says so once (Phase D).
- `load_latest_results()` comment (line ~87 of `06_visualizations.R`): clarify `_combined_` preference applies per estimation-mode, not globally
- `07_results_summary.Rmd` compare-table caption: should note that DIM only ran 6 designs (D7/D8 NAs are expected)
- `IncidenceSpatialCRT_Report.qmd` caption/text alignment: MC SEs table uses tau=1.0 slice (2,560 rows), not full 12,800 — prose now correctly clarifies this distinction

---

## Relationship to SpillSpatialDepSim

IncidenceDesign extends the applied simulation framework from `projects/SpillSpatialDepSim/`.
SpillSpatialDepSim used a small grid (8–12 districts) with SAR estimation to evaluate
block vs. random assignment for NC Department of Adult Correction probation interventions. IncidenceDesign asks the same
core question but at larger scale (100 clusters) with systematic design variation and
heterogeneous incidence modes.

Cross-reference SpillSpatialDepSim results:
```r
here::here("projects", "SpillSpatialDepSim", "results")
# or relative: ../../SpillSpatialDepSim/results/
```
