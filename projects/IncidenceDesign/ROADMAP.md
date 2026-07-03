# IncidenceDesign — Future Work Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the IncidenceDesign project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-07-03:** Tau-sweep simulation complete (8 designs × 12,800 MLE
> scenarios across τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}), power metric and Monte Carlo SEs added,
> statistical comparisons done at every τ level. Two manuscripts drafted from scratch — a
> *Clinical Trials* (SAGE) submission and a longer-form dissertation chapter — plus a
> completed Supplementary Information document for the CTJ submission. The application
> study (NC's 58 Community College service areas) has a complete synthetic-incidence run;
> real SUDDEN-derived data is still pending. The modular Quarto manuscript framework
> referenced by older items below (`paper/manuscript/_application.qmd` etc.) was retired
> 2026-07-02 in favor of the two manuscripts written fresh — see `paper/archive_manuscript/`.

---

## How to Use This Document

Each item follows this format:

```
- [ ] **Title** `[Priority: High/Medium/Low]` `[Effort: Small/Medium/Large]`
  Description of what this involves and why it matters.
  *Notes: dependencies, caveats, or relevant files*
```

**Priority:** How much this would strengthen the work (High = core contribution; Medium = meaningful addition; Low = nice-to-have)

**Effort:** Rough implementation cost (Small = hours; Medium = days; Large = weeks or HPC run)

Move completed items to the [Completed](#-completed) section at the bottom.

---

## 🔬 Statistical Rigor

- [ ] **Sensitivity to spatial weight matrix misspecification** `[Priority: Medium]` `[Effort: Medium]`
  The simulation uses a known W (rook or queen). Real analysts specify W with uncertainty.
  Test design robustness when the analyst's assumed W differs from the true DGP W
  (e.g., DGP uses queen, analyst uses rook; or distance-based vs. contiguity-based).
  *Notes: MLE (`lagsarlm`) is particularly sensitive to W misspecification — important to characterize.*

- [ ] **Randomization-based inference as an alternative estimator** `[Priority: Medium]` `[Effort: Large]`
  All CIs currently come from MLE's model-based SEs. Randomization/permutation inference
  is increasingly preferred in CRT settings (no distributional assumptions). Compare design
  rankings under RI vs. MLE. Particularly relevant for the DIM estimator and small-N settings.

---

## 🌍 Applicability & Generalizability

- [ ] **Replace synthetic application incidence with real SUDDEN-derived data** `[Priority: High]` `[Effort: Medium]`
  `application/` runs the full design comparison on NC's actual 58 Community College service
  areas, but all current incidence numbers and maps are a synthetic Poisson SAR placeholder.
  Real SUD incidence (~100,000 NC death certificates, 2018-2021, classified via the SUDDEN
  algorithm, pooled to county level with SEER population denominators) is being finalized
  with an epidemiology co-investigator under an active IRB protocol.
  *Notes: same schema as the synthetic run — see `application/README.md` "Outstanding Work".
  Affects both manuscripts' Application sections and the CTJ Supplementary Information (S11).*

- [ ] **Add cluster-level baseline covariates to the DGP** `[Priority: Medium]` `[Effort: Medium]`
  Current DGP generates outcomes from spatial structure alone. Adding covariates
  (e.g., poverty rate, population density) and testing whether covariate-adaptive designs
  (Design 8) outperform non-adaptive ones more when covariates are predictive would
  directly inform the SUD/NC application.
  *Notes: the SDM model in `02_incidence_generation.R` / `04_estimation.R` would need covariate terms.*

- [ ] **Heterogeneous (spatially-varying) treatment effects** `[Priority: Medium]` `[Effort: Large]`
  Currently tau is a global scalar. If treatment effects vary spatially (intervention works
  better in high-incidence areas), the "best" design shifts. This connects directly to
  the saturation/incidence-guided designs, which implicitly assume high-incidence areas matter more.

- [ ] **Partial compliance and attrition sensitivity** `[Priority: Low]` `[Effort: Medium]`
  CRTs rarely achieve full compliance. A sensitivity analysis with 10–20% non-compliance
  (random or spatially clustered) would make design recommendations more practically defensible.
  Particularly relevant for the NC law enforcement application context.

---

## ⚙️ Infrastructure & Scalability

- [ ] **Increase replications per scenario for tighter simulation SEs, on Longleaf if needed** `[Priority: Medium]` `[Effort: Large]`
  The tau-sweep (12,800 scenarios) already varies `true_tau`, so this item now covers only
  increasing `n_design_resamples`/`n_outcome_resamples` for tighter Monte Carlo SEs — most
  useful for coverage estimates and Poisson/rare-event incidence scenarios. Current
  `SE_MSE`/`SE_Coverage` (delta-method, `add_mc_ses()`) are already small and uniform across
  designs at the current rep counts, so this is a precision refinement, not a correctness gap.
  *Notes: HPC scripts in `longleaf_setup/` are ready for a per-scenario SLURM array if needed.*

- [ ] **Vary number of clusters (25, 50, 100)** `[Priority: Medium]` `[Effort: Medium]`
  All results are for N = 100 clusters (10×10 grid). Design rankings may differ substantially
  at 25 clusters (5×5) or 50 clusters (realistic for many CRT settings). Design 8 in particular
  relies on an accurate incidence surface, which is harder to estimate from fewer clusters.
  *Notes: requires changes to `01_spatial_setup.R` and `03_designs.R` to support non-10×10 grids.*

---

## 📄 Manuscript & Dissemination

> The items below superseded the pre-2026-07-02 modular-Quarto-manuscript items
> (`_application.qmd`, `_discussion.qmd`, target-journal selection, etc.), which are now
> obsolete — target journal (*Clinical Trials*, SAGE) is chosen, both the Application and
> Discussion sections are fully written in `CTJ_Manuscript.tex` and `Dissertation_Chapter.qmd`,
> and manuscript figures are finalized (see README "Manuscript Development" for the current
> figure inventory). See `paper/archive_manuscript/` for the retired framework.

- [ ] **Consolidated review of all three manuscript documents** `[Priority: High]` `[Effort: Medium]`
  CTJ main text, CTJ Supplementary Information, and the dissertation chapter all have complete
  first drafts (the SI as of 2026-07-03) but have not yet had one combined user review/revision
  pass together.

- [ ] **Deliberate figure-list review across main text + SI** `[Priority: High]` `[Effort: Small]`
  The coverage and tau-sensitivity figures were merged into one 2-panel figure on 2026-07-03
  specifically to fit the new NC incidence map within the CTJ's 6-exhibit cap — a quick fit,
  not a considered final selection. Review the full figure list (both documents) and decide
  what to keep, drop, or further combine.
  *Notes: see `code/14_manuscript_supplement_figures.R` for how main-text/SI figures are generated.*

- [ ] **Real SUDDEN data + maps in the Application sections** `[Priority: High]` `[Effort: Medium]`
  Both manuscripts' Application sections, plus SI Section S11, currently use an explicitly
  labeled synthetic placeholder (incidence numbers *and* maps). Replace once the real
  SUDDEN-derived NC county dataset is finalized.
  *Notes: duplicate of the Applicability & Generalizability item above — tracked in both
  places since it blocks manuscript finalization specifically.*

- [ ] **UNC Graduate School dissertation template** `[Priority: Medium]` `[Effort: Small]`
  `Dissertation_Chapter.qmd` currently renders to a simple double-spaced `article`-class PDF
  matching the SpillSpatialDepSim (Project 1) manuscript format, per user direction. Apply the
  actual UNC template once all dissertation chapters are ready to merge.

- [ ] **Fix the `plot_cd_diagram()` label-collision bug upstream** `[Priority: Medium]` `[Effort: Small]`
  `10_statistical_comparisons.R`'s `plot_cd_diagram()` has a label-collision bug for designs
  with adjacent ranks, independent of image width. Currently worked around in two different
  ways (omitted from the dissertation chapter; re-implemented cleanly, with better label
  staggering, inline in `paper/report/IncidenceDesign_ProjectSummary.qmd` and lifted from there
  into `code/14_manuscript_supplement_figures.R` for the CTJ SI) rather than fixed at the source.

- [ ] **Presentation slides** `[Priority: Low]` `[Effort: Medium]`
  Expand `SpatialCRT_IncidenceDesign_Presentation.qmd` scaffold into full conference slides.

---

## ✅ Completed

- [x] **Build modular 8-design simulation pipeline** — Numbered R scripts (01–11) replacing monolithic Rmd; three incidence modes (iid Uniform, Spatial SAR, Poisson). *(2026-03)*
- [x] **Run full 2,560-scenario MLE simulation** — 8 designs × 4 gamma levels × 4 rho levels × 4 incidence configs × 2 spillover types × `n_sim` reps. Best: D8 (MSE 0.079), Worst: D1 (MSE 0.744). *(2026-03-22)*
- [x] **Formal statistical comparisons** — Friedman test, Nemenyi post-hoc, pairwise Wilcoxon, CD diagrams; rendered to `11_statistical_comparisons_report.{html,pdf}`. *(2026-03-25)*
- [x] **Comprehensive unified project report** — 50+ page `IncidenceSpatialCRT_Report.qmd` covering full pipeline through design recommendations. *(2026-03-23)*
- [x] **Modular Quarto manuscript framework** — Master + child sections (abstract, intro, methods, simulation, application skeleton, discussion skeleton). Retired 2026-07-02 in favor of two manuscripts written fresh; preserved as `paper/archive_manuscript/`. *(2026-03-23)*
- [x] **HPC setup for Longleaf** — SLURM job array scripts for per-scenario parallelization (2,560 tasks), ready to scale. *(2026-03-23)*
- [x] **8-panel design sample figures** — Generated and integrated into manuscript and README. *(2026-03-24)*
- [x] **Tau-sweep: vary `true_tau` across the full parameter grid** — 12,800 scenarios across τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}; design ranking (D3/D8 dominance) stable at every τ level (all conditional Friedman p < 2.2×10⁻¹⁶). *(2026-04-08)*
- [x] **Statistical power added as a primary metric** — P(reject H₀: τ=0) tracked per scenario; power curves computed as a function of τ per design. *(2026-04-08)*
- [x] **Monte Carlo SEs for MSE/coverage estimates** — Delta-method `SE_MSE`, `SE_Coverage`, `SE_Bias` via `add_mc_ses()`; confirmed small and uniform across designs at current rep counts. *(2026-04-08)*
- [x] **NC application study (irregular geometry)** — Full design-comparison pipeline re-implemented for NC's actual 58 Community College service areas (`application/`), all 8 designs adapted, 640-scenario/160,000-fit synthetic-incidence run complete. Real SUDDEN-derived data still pending (tracked above). *(2026-06, ongoing)*
- [x] **Two manuscripts drafted from scratch** — Short *Clinical Trials* (SAGE) submission (`paper/ctj_manuscript/`) and longer-form dissertation chapter (`paper/dissertation_chapter/`), both presenting 6 of the 8 designs (2 dropped as statistically redundant), no DIM-vs-MLE comparison, full design names throughout. *(2026-07-02)*
- [x] **CTJ Supplementary Information drafted and reviewed** — Fulfills the main text's 3 deferred "online supplementary material" items (reproducibility/code, full 8-design comparison, metric formulas/parameter grid); reviewed by a fresh agent against the plan's hard constraints. *(2026-07-03)*
- [x] **NC application maps added to both manuscript documents** — 3 maps from `application/` added to the SI (Figures S8-S10) and main text (Figure 4); coverage + tau-sensitivity figures merged into one to preserve the CTJ's 6-exhibit cap. *(2026-07-03)*
