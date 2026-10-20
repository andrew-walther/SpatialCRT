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

Reasons and history for each rule: `DECISIONS.md`.

- **Application (approved 2026-10-02):** continue the approved plan without
  re-interviewing settled choices. Design authority:
  `docs/plans/nc_sud_application_plan_2026-10-02.md` (its status separates settled from
  proposed choices; never treat proposals as approved). Smoke/pilot outputs are
  verification stages, not conclusions.
- **Application outcome (2026-10-02):** the primary continuous education outcome has no
  SUD-incidence baseline term (β = 0); observed rates inform allocation only. β = 1 is an
  explicitly hypothetical sensitivity. τ is constant within a scenario.
- **Methods extensions (author, 2026-10-20):** Y is the measured education endpoint, not
  future SUD incidence (no count-outcome model); partial participation and power run on
  NC only, with the fast fitter and new output folders, never editing hashed production
  code. Plan: `docs/plans/methods_extensions_plan_2026-10-20.md`.
- **CTJ authors:** Walther, Habib, Simpson, Lin, in that order; Lin corresponding. Never
  invent an IRB protocol number or claim approval/exemption.
- **Manuscript conventions (author, 2026-10-12):** define every abbreviation at first
  body use and avoid internal pipeline jargon; one design-name set everywhere
  (Checkerboard, High Incidence Focus, Plain Saturation, Isolation Buffer, 2×2 Blocking,
  Balanced Quartiles, Balanced Halves, Incidence-Guided Saturation, Simple Random Sampling);
  frame incidence guidance as matching plain saturation's accuracy while targeting
  high-burden regions (do not call it "no benefit"); no repository hyperlink in manuscripts;
  no Codex/Claude mentions in code headers or documents (AI disclosure is a later TODO).
- **References to Project 1 (2026-09-27):** the chapter may say "Chapter 2". CTJ must be
  self-contained ("In previous work...") and cite the BMC Medical Research Methodology
  paper, never "Chapter 2" or "Project 1".
- **Checkerboard placement (2026-10-02):** keep it in the main design comparison, with
  identification diagnostics in the appendix; no numerical cross-project comparison.
- **No single recommended design (author, 2026-10-14):** report the best designs by
  scenario (mainly spillover regime).
- **Revisions reach every version (author, 2026-10-19):** CTJ edits go in
  `CTJ_Manuscript.tex`/`CTJ_Exhibits.tex` (the reading layout follows); rebuild both
  `CTJ_Manuscript.pdf` and `CTJ_Reading.pdf`. After any edit to the CTJ paper or Chapter 3,
  check whether the other needs it and say so (content edits usually do; journal-structure
  edits usually do not). The author's annotated copies live in the gitignored
  `manuscript_ctj/revision_markup/` as `CTJ_Reading_Markup_V<n>.pdf`; never edit or overwrite them.
- **Project 1 is final.** Treat `projects/SpillSpatialDepSim/` as read-only. Changes
  reach `bios-dissertation` only through the Chapter 3 sync (below) unless the author
  asks otherwise.
- **Allocation-risk numbers** (selected-surface/corner pilot, `results/allocation_risk/`)
  are never substituted for full-study numbers.
- **Where to write things down:** key findings go in `docs/plans/` (open items:
  `docs/plans/TODO.md`); dated progress checkpoints go in the progress log of
  `ROADMAP.md`; author decisions and their reasons go in `DECISIONS.md`, with only a
  one-line rule here. Restricted SUD data details: `application/AGENTS.md`.

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

**Primary goal: compare the treatment sampling designs** (eight strategies plus SRS in the full benchmark/application; six in detailed grid rankings) under realistic spatial conditions. The key question is which designs estimate
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

**Manuscript prose rule (both `paper/ctj_manuscript/` and `paper/thesis_chapter/`):
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
| `19_manuscript_exhibits.R` | ~60 | Current nine-design grid regime CSV, provenance and six-design CTJ vector figure | Reads verified September 27 oracle outputs only |
| `21_downside_by_design.R` | ~165 | Per-design mean and worst-decile conditional MSE as ratios to SRS (Table `tab:downside`) and the allocation-distribution figure (CTJ Figure 3, chapter/SI `fig:alloc-dist`) | Reads completed aggregate CSVs; writes `results/allocation_risk/summary/downside_by_design.csv` |
| `paper/tools/verify_manuscripts.py` | ~290 | CSV-to-table and chapter-to-supplement reconciliation (incl. merged `tab:nc-annual`); assets/references/citations/exhibit count; `--chapter`, `--pending-si` for candidate copies | Stdlib; deliberately wrong-number fixture must fail |
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
| **paper/archive/** | | | |
| `pre_unification_20260702/` (retired) | — | Byte-identical snapshot of the pre-2026-07-02 `paper/manuscript/` + `paper/section_drafts/` content | Reference/fact-check source only — do not edit in place; both current manuscripts were written fresh, not derived from this |
| `pre_completion_20261002/`, `manuscript_provenance_20261002.json` | — | Sources before the Oct 2 completion pass; Oct 2 SHA-256 manifest | Frozen provenance; do not reuse numbers |
| **paper/ctj_manuscript/** | | | |
| `manuscript_ctj/CTJ_Manuscript.tex` | ~220 | *Clinical Trials* (SAGE) submission draft, six designs + SRS in main text (variants in SI), exactly 6 exhibits (2 tables + 4 figures) | Compiles via `sagej.cls` — see Build gotchas |
| `manuscript_ctj/CTJ_Reading.tex`, `CTJ_Exhibits.tex` | wrapper/shared exhibits | Reading wrapper selects the master's two-column branch; six shared exhibit macros provide identical captions, cells and graphics in both layouts | Build reading wrapper or submission master, not exhibit file |
| `supplementary_info/Supplementary_Information.tex` | full draft | Chapter-derived detailed methods/results, full application and seven appendix components (S1–S10); all 22 source tables retained | Plain article class, S-prefixed exhibits, shared bib |
| `cover_letter/`, `figures/` (incl. `Figure_Legends.md`), `submission_checklist.md` | — | Cover letter draft; shared main/SI figures written by the R scripts | Figure paths in both TeX sources are `../figures/...` |
| **paper/thesis_chapter/** | | | |
| `Thesis_Chapter.qmd` | ~330 | Longer-form dissertation chapter, full long-form chapter/appendix, no length ceiling | Quarto → simple double-spaced `article`-class PDF matching Project 1's format |
| **paper/** | | | |
| `SpatialCRT_IncidenceDesign.bib` | ~99K, 81 entries | Single shared bibliography for both manuscripts (renamed from `IncidenceDesign_shared.bib` 2026-09-24; `baird`/`leung` entries synced to `bios-dissertation/prelim/references.bib`; `habib_temporal_2026` added 2026-09-24 = Ashkan Habib's unpublished BIOS master's paper, source of the SUD county data) | Zotero base + SUD/SUDDEN citations + design-theory citations |
| `SAGE_Journal_Template/` | — | `sagej.cls`, `SageH.bst`, `SageV.bst` | CTJ manuscript render dependency |

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

## Current State

Status and dated checkpoints: `ROADMAP.md` (Progress Log). Open items:
`docs/plans/TODO.md`. Numerical/exhibit checks after any manuscript edit:
`Rscript code/tests/test_manuscript_exhibits.R` and
`python3 paper/tools/verify_manuscripts.py`.

## Dissertation Chapter 3 Sync

- `paper/thesis_chapter/Thesis_Chapter.qmd` is the only file anyone edits
  (author decision 2026-09-25).
- `paper/thesis_chapter/tools/sync_to_prelim.sh` generates
  `bios-dissertation/prelim/project-proposals/project2-incidence/draft/project2-incidence-draft.qmd`
  (CHAPTER 3 / APPENDIX B). It checks every citekey against the master bib, failing loud,
  renders with RStudio's Quarto 1.10, and commits only `project2-incidence/` paths there.
  It never pushes and skips PDF-only re-renders.
- `paper/thesis_chapter/tools/sync_ctj_to_bios.sh` (same post-commit hook, triggered by
  `paper/ctj_manuscript/manuscript_ctj/`) copies the committed `CTJ_Manuscript.tex`,
  `CTJ_Exhibits.tex` and `CTJ_Reading.pdf` to `bios-dissertation/manuscripts/project2-incidence/`
  as `project_2_paper.tex`, `project_2_paper_exhibits.tex` and `project_2_paper.pdf`, and
  commits them there. No render; it warns if the `.tex` was committed after the reading PDF.
- A `post-commit` hook (install with `tools/install_hooks.sh`) runs it when a commit touches
  the chapter, its figures, `tools/` or the bib. Its log is `tools/sync.log` (gitignored).
- The prelim YAML header is `tools/prelim_header.yml` (`colorlinks: false`, so links print
  black). Cross-references are `\label`/`\ref`, so numbers read 3.x / B.x in the prelim.

- **Length review (adopted 2026-10-13):** cut or moved text is in
  `notes/trimmed-material.md`. The verifier's `--chapter` and `--pending-si` options
  check any future candidate copy.

## Build gotchas

- **CTJ layouts:** `CTJ_Manuscript.tex` is the master for both layouts;
  `CTJ_Reading.tex` only selects the two-column reading layout. Edit captions and cells in
  `CTJ_Exhibits.tex`; never duplicate manuscript content in the wrapper. Main figures are
  `figures/compact/` exports from `code/20_manuscript_figure_revision.R`. The reading layout
  uses the `balance` package (`\balance` before the bibliography) and no forced page breaks.
- **CTJ manuscript/SI:** need TinyTeX on `PATH`; run `latexmk <file>.tex` inside
  `paper/ctj_manuscript/manuscript_ctj/` or `.../supplementary_info/`. Symlinks there (`sagej.cls`,
  `SageV.bst`, `SpatialCRT_IncidenceDesign.bib`) replace the old `TEXINPUTS`/`BSTINPUTS` setup,
  and each folder's `latexmkrc` sends intermediates to its gitignored `.build/`.
  A plain `pdflatex` build also works but leaves intermediates beside the sources;
  `git clean -fX -- paper/ctj_manuscript` removes them (it only touches gitignored files).
- **Dissertation chapter bibliography:** reference it via the `shared-refs.bib` symlink in `paper/thesis_chapter/`, not `../SpatialCRT_IncidenceDesign.bib` directly — Quarto's pandoc→LaTeX conversion mishandles underscores in bib filenames referenced from YAML.
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

Full table: `README.md` ("Known Issues & Bugs Fixed").

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

## Known Minor Issues

Low-priority cleanup items: `docs/plans/TODO.md` ("Known minor code issues").

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
