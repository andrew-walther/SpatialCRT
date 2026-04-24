# IncidenceDesign — Future Work Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the IncidenceDesign project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-03-26:** Core simulation complete (8 designs × 2,560 MLE scenarios),
> statistical comparisons done, manuscript framework and comprehensive report in place.

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

- [ ] **Vary `true_tau` across the parameter grid** `[Priority: High]` `[Effort: Large]`
  All 2,560 scenarios currently fix tau = 1.0. Design rankings may shift near the null (tau ≈ 0),
  at small effects (tau = 0.25, 0.5), or large effects (tau = 2.0+). Essential for generalizable
  design recommendations and for computing power curves.
  *Notes: primary scale-up candidate for Longleaf. Expands grid to ~5× current size. See `longleaf_setup/`.*

- [ ] **Report simulation standard errors for MSE and coverage estimates** `[Priority: High]` `[Effort: Small]`
  Monte Carlo estimates of MSE and coverage carry simulation uncertainty. Reporting SE alongside
  point estimates (or confidence bands on figures) strengthens credibility of design comparisons,
  especially for close calls (e.g., Design 8 vs. Design 3, MSE 0.079 vs. 0.080).
  *Notes: can be added directly to `07_results_summary.qmd` and `11_statistical_comparisons_report.qmd`.*

- [ ] **Add statistical power as a primary metric** `[Priority: High]` `[Effort: Medium]`
  Currently tracking bias, SD, MSE, and coverage — but not power (P(reject H₀ | tau ≠ 0)).
  Power is what study planners actually need. Requires varying tau (see above) to produce
  power curves as a function of effect size per design.
  *Notes: add to `04_estimation.R` output and `05_run_simulation.R` metrics tracked.*

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

- [ ] **Test on irregular geometries (NC county adjacency network)** `[Priority: Medium]` `[Effort: Medium]`
  The 10×10 regular lattice is analytically clean but doesn't reflect real geographic
  structures. Running the simulation on NC's actual county adjacency graph (irregular shapes,
  water boundaries, urban/rural clustering) would directly strengthen the SUD application framing.
  *Notes: `sf` + `spdep` can build the NC county weights matrix from a shapefile. See `01_spatial_setup.R`.*

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

- [ ] **Longleaf scale-up: expanded parameter grid (vary tau + more reps)** `[Priority: High]` `[Effort: Large]`
  The HPC scripts in `longleaf_setup/` are ready for a per-scenario SLURM array (2,560 tasks).
  The next major run should expand to include `true_tau` variation and/or increase replications
  per scenario for tighter simulation SEs. Estimate: ~5–10× more compute than current run.
  *Notes: `longleaf_setup/submit_array.sl`, `simulation.R`, `aggregate_results.R` are all ready.*

- [ ] **Vary number of clusters (25, 50, 100)** `[Priority: Medium]` `[Effort: Medium]`
  All results are for N = 100 clusters (10×10 grid). Design rankings may differ substantially
  at 25 clusters (5×5) or 50 clusters (realistic for many CRT settings). Design 8 in particular
  relies on an accurate incidence surface, which is harder to estimate from fewer clusters.
  *Notes: requires changes to `01_spatial_setup.R` and `03_designs.R` to support non-10×10 grids.*

- [ ] **Increase replications per scenario for tighter simulation SEs** `[Priority: Medium]` `[Effort: Large]`
  More Monte Carlo reps reduce simulation noise, especially important for coverage estimates
  (which require many reps to distinguish, say, 0.92 from 0.94) and for Poisson/rare-event
  incidence scenarios.
  *Notes: natural companion to the Longleaf scale-up item above.*

---

## 📄 Manuscript & Dissemination

- [ ] **Update the full-length report to match the summary-report rendering fixes** `[Priority: Low]` `[Effort: Small]`
  The brief `IncidenceDesign_ProjectSummary.qmd` has been updated with safer table formatting and
  multi-format output handling. Mirror those fixes in `IncidenceSpatialCRT_Report.qmd` later so the
  full report renders cleanly in both HTML and PDF with consistent table widths and math rendering.
  *Notes: lower priority than manuscript completion; defer until the next full-report refresh.*

- [ ] **Complete `_application.qmd` (SUD in NC section)** `[Priority: High]` `[Effort: Medium]`
  Currently a skeleton (~15 lines). Needs: description of the NC SUD context, county-level
  data summary, justification for Poisson incidence mode, and a worked design recommendation
  for a hypothetical NC-scale trial.
  *Notes: see `paper/manuscript/_application.qmd`; Poisson base rate 35/100,000 from Mirzaei et al.*

- [ ] **Complete `_discussion.qmd`** `[Priority: High]` `[Effort: Medium]`
  Currently a skeleton (~30 lines). Needs: interpretation of Design 8 vs. Design 3 result,
  limitations (regular grid, fixed tau, model-based inference), connections to broader CRT
  design literature, and directions for future work (can draw directly from this ROADMAP).
  *Notes: see `paper/manuscript/_discussion.qmd`.*

- [ ] **Polish and finalize manuscript figures** `[Priority: Medium]` `[Effort: Small]`
  Placeholder figure calls exist in `_simulation.qmd`. Finalize which figures go into the
  manuscript (vs. supplementary) and ensure captions, labels, and color schemes are
  publication-ready.
  *Notes: 8-panel design figure (`design_samples_8panel.png`) is already integrated.*

- [ ] **Select target journal and adapt manuscript formatting** `[Priority: Low]` `[Effort: Small]`
  Candidate journals: *Statistics in Medicine*, *Clinical Trials*, *Spatial and Spatio-temporal
  Epidemiology*, *American Journal of Epidemiology*. Each has different length/format requirements.
  *Notes: Quarto supports journal-specific templates (`quarto use template ...`).*

---

## ✅ Completed

- [x] **Build modular 8-design simulation pipeline** — Numbered R scripts (01–11) replacing monolithic Rmd; three incidence modes (iid Uniform, Spatial SAR, Poisson). *(2026-03)*
- [x] **Run full 2,560-scenario MLE simulation** — 8 designs × 4 gamma levels × 4 rho levels × 4 incidence configs × 2 spillover types × `n_sim` reps. Best: D8 (MSE 0.079), Worst: D1 (MSE 0.744). *(2026-03-22)*
- [x] **Formal statistical comparisons** — Friedman test, Nemenyi post-hoc, pairwise Wilcoxon, CD diagrams; rendered to `11_statistical_comparisons_report.{html,pdf}`. *(2026-03-25)*
- [x] **Comprehensive unified project report** — 50+ page `IncidenceSpatialCRT_Report.qmd` covering full pipeline through design recommendations. *(2026-03-23)*
- [x] **Modular Quarto manuscript framework** — Master + child sections (abstract, intro, methods, simulation, application skeleton, discussion skeleton). *(2026-03-23)*
- [x] **HPC setup for Longleaf** — SLURM job array scripts for per-scenario parallelization (2,560 tasks), ready to scale. *(2026-03-23)*
- [x] **8-panel design sample figures** — Generated and integrated into manuscript and README. *(2026-03-24)*
