# IncidenceDesign — Future Work Roadmap

> **Current status (2026-10-25):** observed application (2026-10-02) and full chapter/appendix, CTJ/SI drafts complete. Methods extensions A0–A7 ran 2026-10-20/21 and the author decisions on them were recorded 2026-10-21. Steps 7–8 of `docs/plans/2026-10-20-methods-extensions-plan.md` were drafted overnight 2026-10-25 (figures, SRS in tests, chapter/SI/CTJ text, talk outline); next is the author's review of those drafts. A 20-minute application talk is due 2026-10-29 (`docs/prompts/2026-10-25-talk-kickoff-prompt.md`). See current progress log and `docs/TODO.md`. Historical statuses below are retained by date.
>
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
> real SUDDEN-derived data is still pending *(superseded 2026-10-02: the application now uses observed 2018–2021 SUD data)*. The modular Quarto manuscript framework
> referenced by older items below (`paper/manuscript/_application.qmd` etc.) was retired
> 2026-07-02 in favor of the two manuscripts written fresh — see `paper/archive/pre_unification_20260702/`.
>
> **Status as of 2026-09-24 (later): simulation revised and re-run; downstream regeneration
> (Phase C) in progress.**
> - The M1–M8 fixes are implemented and tested. The pilot directions matched the scratch pilot.
> - Full run: 12,800 scenarios × 250 fits × 2 estimators in 13.9 min, using the validated lean
>   ML engine. It's verified, and a 1% `lagsarlm` cross-check passed (max difference 6e-9).
> - Headline (queen, τ = 1): Incidence-Guided Saturation Quadrants best (MSE 0.091), then
>   Balanced Quartiles (0.131). Checkerboard is worst (1.09). Coverage ≈ 0.94 for every design
>   except Checkerboard × rook (0.15, τ not identified).
> - Decisions: keep the 6 designs; Balanced Quartiles treats exactly 50.
> - Framing: which designs work under heterogeneous incidence, not whether incidence
>   knowledge is needed.
> - April 2026 numbers are superseded, including those in the Completed section below.
>   The CTJ/SI still carry them until manuscript step 5.
>
> Earlier the same day — **simulation revision planned; current results are provisional.**
> A review of the dissertation chapter traced several results to simulation oversights:
> - Designs saw a different incidence surface than the outcome model used.
> - Poisson incidence at 1,000 people per cluster produced heavy ties.
> - Seeding was not reproducible.
> - Deterministic designs were copied 25×.
> - Checkerboard under rook contiguity was silently aliased.
>
> A pilot suggests the top and bottom designs are stable, while High Incidence Focus changes
> substantially. Fix, test, re-run, and regenerate every downstream result before the
> manuscripts' Methods/Results are rewritten:
> [`docs/plans/2026-09-24-simulation-revision-plan.md`](docs/plans/2026-09-24-simulation-revision-plan.md),
> which is step 0.5 of `docs/plans/2026-09-18-manuscript-unification-plan.md`.
> The 2026-04-08 numbers below and in the Completed section are superseded once the re-run lands.

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

- [x] **Replace synthetic application incidence with real SUDDEN-derived data** *(2026-10-02)*
  Done: the application now uses observed NC county SUD incidence (2018–2021, SEER
  denominators; restricted source files gitignored) and the observed-incidence production run
  is complete and verified. The synthetic Poisson SAR run is kept as historical. See
  `application/README.md` and `docs/findings/2026-10-02-nc-sud-application-findings.md`.

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
> Discussion sections are fully written in `CTJ_Manuscript.tex` and `Thesis_Chapter.qmd`,
> and manuscript figures are finalized (see README "Manuscript Development" for the current
> figure inventory). See `paper/archive/pre_unification_20260702/` for the retired framework.

- [x] **Consolidated review of all three manuscript documents** `[Priority: High]` `[Effort: Medium]`
  CTJ main text, CTJ Supplementary Information, and the dissertation chapter all have complete
  first drafts (the SI as of 2026-07-03) but have not yet had one combined user review/revision
  pass together.
  *Done 2026-10-12: four-track comprehensive review, approved fixes applied
  (`docs/reviews/2026-10-12-ctj-comprehensive-review-findings.md`).*

- [x] **Deliberate figure-list review across main text + SI** `[Priority: High]` `[Effort: Small]`
  The coverage and tau-sensitivity figures were merged into one 2-panel figure on 2026-07-03
  specifically to fit the new NC incidence map within the CTJ's 6-exhibit cap — a quick fit,
  not a considered final selection. Review the full figure list (both documents) and decide
  what to keep, drop, or further combine.
  *Notes: see `code/14_manuscript_supplement_figures.R` for how main-text/SI figures are generated.*
  *Done 2026-10-02: approved six-exhibit revision (`docs/findings/2026-10-02-ctj-exhibit-revision-findings.md`).*

- [x] **Real SUDDEN data + maps in the Application sections** *(2026-10-02)*
  Done: the chapter, CTJ manuscript and SI Application sections use the observed-incidence
  results (`docs/findings/2026-10-02-manuscript-completion-findings.md`). Duplicate of the
  Applicability & Generalizability item above.

- [ ] **UNC Graduate School dissertation template** `[Priority: Medium]` `[Effort: Small]`
  `Thesis_Chapter.qmd` currently renders to a simple double-spaced `article`-class PDF
  matching the SpillSpatialDepSim (Project 1) manuscript format, per user direction. Apply the
  actual UNC template once all dissertation chapters are ready to merge.

- [ ] **Fix the `plot_cd_diagram()` label-collision bug upstream** `[Priority: Medium]` `[Effort: Small]`
  `10_statistical_comparisons.R`'s `plot_cd_diagram()` has a label-collision bug for designs
  with adjacent ranks, independent of image width. Currently worked around in two different
  ways (omitted from the dissertation chapter; re-implemented cleanly, with better label
  staggering, inline in `paper/report/IncidenceDesign_ProjectSummary.qmd` and lifted from there
  into `code/14_manuscript_supplement_figures.R` for the CTJ SI) rather than fixed at the source.

- [ ] **Presentation slides** `[Priority: Low]` `[Effort: Medium]`
  Slides are built in `bios-dissertation/slides/` (prelim deck in per-part pieces; a Project 1
  talk would be trimmed from its prelim section). The unused scaffold here was deleted 2026-10-13.

## 2026-09 Simulation Revision (step 0.5)

- [x] **Revision M1–M8** — matched surfaces, 100,000 per Poisson cluster, random ties with exact
  N/2, key-based seeds, one noise column per fit, aliasing flags, non-oracle sensitivity estimator,
  surface-level MC SEs, manifest-checked checkpoints, lean validated ML engine. Spec:
  `docs/plans/2026-09-24-simulation-revision-spec.md`. *(2026-09-24)*
- [x] **Full re-run + verification + 1% lagsarlm cross-check** *(2026-09-24)*
- [x] **Downstream regeneration (Phase C)** — 12/13/14 and the 06/08/10 PDFs, reports
  00/07/11 (`516e127`), 09 and the project Report/ProjectSummary (`6e3aaa0`), and the docs
  (`fe4d65f`) regenerated. *(2026-09-24)*
- [x] **Chapter Methods/Results rewrite (Phase D)** — robust vs. fragile; queen primary; rook
  sensitivity including the Checkerboard-rook anomaly; non-oracle sensitivity. Chapter complete
  2026-10-02 (`docs/findings/2026-10-02-manuscript-completion-findings.md`). Instead of one
  recommended design, the Discussion gives guidance by spillover scenario (added 2026-10-11).
  *(2026-10-02)*

### Future work surfaced by the revision
- [ ] **Latent-risk / noisy-snapshot DGP** `[Priority: Medium]` — designs see a noisy snapshot of a
  latent risk surface, and the outcome depends on the latent surface (not the observed one).
- [ ] **Heterogeneous populations in the simulation** `[Priority: Medium]` — `pop_mode = "heterogeneous"`
  exists in 02 but is unused; the application will use real county/CC populations.
- [ ] **Port the lean engine to the application** `[Priority: Low]` — the application still uses `lagsarlm`.
- [ ] **Poisson ρX = 0 config** `[Priority: Low]` — offered separately; not in the grid.

---

## ✅ Completed

- [x] **Build modular 8-design simulation pipeline** — Numbered R scripts (01–11) replacing monolithic Rmd; three incidence modes (iid Uniform, Spatial SAR, Poisson). *(2026-03)*
- [x] **Run full 2,560-scenario MLE simulation** *(SUPERSEDED 2026-09-24 by the revision)* — 8 designs × 4 gamma levels × 4 rho levels × 4 incidence configs × 2 spillover types × `n_sim` reps. Best: D8 (MSE 0.079), Worst: D1 (MSE 0.744). *(2026-03-22)*
- [x] **Formal statistical comparisons** — Friedman test, Nemenyi post-hoc, pairwise Wilcoxon, CD diagrams; rendered to `11_statistical_comparisons_report.{html,pdf}`. *(2026-03-25)*
- [x] **Comprehensive unified project report** — 50+ page `IncidenceSpatialCRT_Report.qmd` covering full pipeline through design recommendations. *(2026-03-23)*
- [x] **Modular Quarto manuscript framework** — Master + child sections (abstract, intro, methods, simulation, application skeleton, discussion skeleton). Retired 2026-07-02 in favor of two manuscripts written fresh; preserved as `paper/archive/pre_unification_20260702/`. *(2026-03-23)*
- [x] **HPC setup for Longleaf** — SLURM job array scripts for per-scenario parallelization (2,560 tasks), ready to scale. *(2026-03-23)*
- [x] **8-panel design sample figures** — Generated and integrated into manuscript and README. *(2026-03-24)*
- [x] **Tau-sweep: vary `true_tau` across the full parameter grid** *(results SUPERSEDED 2026-09-24; the sweep is kept in the revised run)* — 12,800 scenarios across τ ∈ {0.8, 1.0, 1.5, 2.0, 3.0}; design ranking (D3/D8 dominance) stable at every τ level (all conditional Friedman p < 2.2×10⁻¹⁶). *(2026-04-08)*
- [x] **Statistical power added as a primary metric** — P(reject H₀: τ=0) tracked per scenario; power curves computed as a function of τ per design. *(2026-04-08)*
- [x] **Monte Carlo SEs for MSE/coverage estimates** *(SUPERSEDED 2026-09-24: SEs now come from the 10 surface-level means, stored by 05)* — Delta-method `SE_MSE`, `SE_Coverage`, `SE_Bias` via `add_mc_ses()`; confirmed small and uniform across designs at current rep counts. *(2026-04-08)*
- [x] **NC application study (irregular geometry)** — Full design-comparison pipeline re-implemented for NC's actual 58 Community College service areas (`application/`), all 8 designs adapted, 640-scenario/160,000-fit synthetic-incidence run complete. Observed 2018–2021 SUD-incidence run completed 2026-10-02 (`docs/findings/2026-10-02-nc-sud-application-findings.md`). *(2026-06 to 2026-10-02)*
- [x] **Two manuscripts drafted from scratch** — Short *Clinical Trials* (SAGE) submission (`paper/ctj_manuscript/`) and longer-form dissertation chapter (`paper/thesis_chapter/`), both presenting 6 of the 8 designs (2 dropped as statistically redundant), no DIM-vs-MLE comparison, full design names throughout. *(2026-07-02)* *(Design count history: 6 designs at first, expanded to 8, cut back to 6 for the CTJ main text because some were similar, then the SRS benchmark was added. Current: CTJ main text 6 + SRS; thesis chapter 8 + SRS, with 6 in the detailed grid rankings.)*
- [x] **CTJ Supplementary Information drafted and reviewed** — Fulfills the main text's 3 deferred "online supplementary material" items (reproducibility/code, full 8-design comparison, metric formulas/parameter grid); reviewed by a fresh agent against the plan's hard constraints. *(2026-07-03)*
- [x] **NC application maps added to both manuscript documents** — 3 maps from `application/` added to the SI (Figures S8-S10) and main text (Figure 4); coverage + tau-sensitivity figures merged into one to preserve the CTJ's 6-exhibit cap. *(2026-07-03)*

---

## Open Work Moved from AGENTS.md (2026-10-01)

Moved verbatim from `AGENTS.md` "Extensions & Future Work Roadmap" (done rows removed
2026-10-01).

### Simulation Extensions (open)

| Priority | Extension | Implementation Note |
|----------|-----------|---------------------|
| Medium | **Heterogeneous population** Poisson mode | `pop_mode = "heterogeneous"` in `05`; extend `generate_incidence_poisson()` for unequal cluster sizes |
| Medium | **Grid sensitivity** | Change `grid_dim` to 8 or 15 in `05`; tests D3/D8 dominance at different spatial scales |
| Low | **DIM tau-sweep** | Re-run DIM across all τ levels for power curve comparison (DIM is confirmed naive baseline) |
| Low | **Heterogeneous beta** | Vary `beta` coefficient across simulation configs |

### Manuscript Development (open)

First full drafts of both manuscripts, plus the CTJ Supplementary Information,
are complete (see Current State above: `paper/ctj_manuscript/`,
`paper/thesis_chapter/`). Remaining work:

| Priority | Task | File |
|----------|------|------|
| **NEXT UP** | Consolidated author/Claude review of completed study and full drafts; finish author-only ethics/data-use and submission metadata | `docs/prompts/2026-10-02-claude-comprehensive-review-prompt.md`; `paper/ctj_manuscript/submission_checklist.md` |
| **High** | Coauthor approval, mandatory reviewer corrections if any, then approved CTJ submission | Full chapter, CTJ and supplement |
| Medium | Final UNC dissertation-template integration when chapters are combined | Canonical Chapter 3 source and approved sync |
| Low | Project 2 prelim/conference slides | `bios-dissertation/slides/` |

---

## Progress Log (moved from AGENTS.md 2026-10-01)

Dated status checkpoints, newest first, moved verbatim from `AGENTS.md`. Add new
checkpoints here, not in `AGENTS.md`.

### 2026-10-26 — application t + df coverage; both-interval table

- `code/27_interval_comparison.R`: application yearly t + df coverage (Wald column
  reproduces the production means exactly) and a Wald vs t + df comparison CSV.
- Application coverage switched to t + df in the yearly table (chapter, SI), Results
  text and the CTJ abstract/Results; new `tab:interval-comparison` (chapter A9, SI S11),
  checked by the verifier. CTJ body 3,489 words, abstract 411. Removed a duplicated
  Checkerboard sentence in the SI application paragraph.

### 2026-10-26 — grid t + df intervals; Table 1 categories

- `code/26_grid_tdf_intervals.R` re-fit the full grid (16 min, 10 workers) and
  reproduced every stored Wald metric exactly; t + df coverage/power now feed scripts
  12/13/14/19 through `apply_tdf_intervals()` (06). Updated every grid coverage and
  power number in the CTJ (Methods, Results, Limitations, Table 2), SI and Chapter 3
  (five tables, four figures, prose, metric formulas). CTJ body 3,486 words.
- Table 1 shows each design's category in italics under its name.

### 2026-10-26 — pre-review cleanup

- Removed "personal communication" for Habib's counts (CTJ, SI, chapter); CTJ body
  recounted at 3,477 (title page, checklist, cover letter updated).
- New review brief `docs/prompts/2026-10-26-claude-comprehensive-review-prompt.md`;
  stale TODO entries closed. Grid coverage with the t + df interval needs a re-fit
  (no per-fit grid output stored) and awaits the author's go-ahead.

### 2026-10-25 — hybrid Table 1, new design figure, flushend reading layout

- Reading layout uses `flushend` instead of `balance`, so References follows the
  declarations on page 8 again.
- New `code/25_design_assignment_figures.R`: CTJ Table 1 example maps, SI/Chapter 3
  panel figure (`design_assignments_panels.pdf`) and a wide slide copy in
  `presentation/figures/`. Table 1 is now a hybrid table (reading layout still 9 pages);
  the verifier counts its eight rows and checks every map file exists.
- The old `design_samples_6panel.pdf` (CTJ and chapter figure folders) is no longer
  referenced by any manuscript; `presentation/OUTLINE.md` still points to its PNG.

### 2026-10-25 (afternoon) — author feedback on the overnight drafts

- "Saturation Quadrants" restored as the one name in text, tables and figures (reversing
  "Saturation Regions"); the grid figures and the design-sample figure were regenerated so
  every manuscript figure matches.
- Manuscript figures lose in-figure titles and notes (captions hold them) and print true
  hyphens (soft-hyphen labels; PNGs rendered from the PDFs). The power figure was relabeled.
- Colleges needed are whole colleges, rounded up (35 to 37 under both-arms spillover; 28
  for the saturation designs against 29 to 32 at 1.5 SD under control-only).
- The SRS comparison now reports allocation spread beside the mean difference.
- All remaining titled manuscript figures cleaned: bias-variance and coverage-by-design
  (codes 12/14), the NC maps, service-area map and application figures
  (`application/code/manuscript_maps.R`, title-free `_journal` PDFs); captions extended.

### 2026-10-25 (overnight) — methods-extension steps 7–8: SRS in tests, figures, manuscripts, talk outline

- SRS joins every grid hypothesis test (author decision, `DECISIONS.md`): `code/12` (six
  designs + SRS), `code/14` (eight + SRS, CD diagrams), `code/15` (each design vs SRS
  within regime, ratio forest plot). Queen Friedman chi2 = 560.33 on 6 df; SRS and
  Balanced Quartiles not separated; under both-arms spillover the saturation designs
  have 4-8% higher MSE than SRS (significant), under control-only about 37% lower.
- New figures (`code/23`): power vs participating colleges, colleges needed for 80%
  power, grid spread vs incidence-guidance gain with NC points. Tables from the stored
  CSVs (`code/24`). A figure reviewer's fixes applied (ratio axis, grayscale, labels).
- Chapter 3: application methods (partial participation, sigma_c and ICC, corrected t
  interval, further sensitivity analyses), results (realistic participation, targeting,
  calibration, model robustness), discussion and limitations, new Appendix A9 with
  model specifications (CAR, Bayesian CAR, Bayesian SAR) and seven tables. SI mirrors it
  (S10 subsection). CTJ: about 95 words (methods, results, limitations, abstract);
  body 3,483 of 3,500, abstract 414 of 425. "Saturation Regions" in text and tables.
- New citations verified at doi.org and added to both bibs (Besag 1974, Leroux 2000,
  Lee 2013, Geweke 1992). Chapter 3 sync ran (bios-dissertation f0647c8).
- Talk outline draft: `presentation/OUTLINE.md` and `outline.html` (19 + 8 slides).
- Open: author review of the drafts; older figures still print "Saturation Quadrants".

### 2026-10-25 — docs folder reorganized by type

- `docs/plans/` held every project document. They now sit in `docs/plans/`, `findings/`,
  `reviews/`, `prompts/` and `meetings/`, named `YYYY-MM-DD-topic-type.md`. Each folder
  has a README giving every document's purpose and status. `TODO.md` is now `docs/TODO.md`.
- References were updated everywhere except hashed code (`code/04`, `05`,
  `application/code/`), `CTJ_Manuscript.tex` and frozen archives. `docs/README.md` maps
  their old paths to the new ones. Links resolve, and `verify_manuscripts.py` and
  `test_manuscript_exhibits.R` pass.

### 2026-10-21 — affiliation, Gan citation and IRB protocol in all versions

- Habib → Department of Epidemiology (CTJ, Chapter 3); Gan et al. (2019) cited for the
  SUD case-adjudication algorithm (CTJ, SI, Chapter 3); IRB 14-2036 stated in CTJ Methods
  (determination still a placeholder) and Chapter 3. Rebuilt CTJ_Manuscript, CTJ_Reading,
  SI and Thesis_Chapter PDFs; `verify_manuscripts.py` and `test_manuscript_exhibits.R` pass.
- Awaiting Simpson: IRB disclosures/permissions and feedback on the application's fit
  to his pilot idea (`docs/TODO.md`).

### 2026-10-20 — methods extensions (post-advisor-meeting analyses), steps 1–6

- Implemented the approved plan (`docs/plans/2026-10-20-methods-extensions-plan.md`)
  in new files `application/code/methods_ext_*.R` and tests; production files unchanged (V9).
- Ran A0–A2 and A5 Stage 1 on stored fits; A3 (301.8M fits) and A4 (98.5M) overnight
  locally, all 576 tasks OK, no failed fits; then A7, A5 Stage 2, A6 CAR/reverse and
  Bayesian CAR. Findings, failed checks (V2 at N <= 35, V6 small regions, V8 along N,
  V5 guard deviation) and open author decisions:
  `docs/findings/2026-10-20-methods-extensions-findings.md`.
- No manuscript, chapter, CTJ, SI or slide edits (steps 7–8 await approval). No push.

### Status text moved from AGENTS.md (2026-10-14)

Moved verbatim when `AGENTS.md` was trimmed to standing rules and reference material.

**Verified completion checkpoint (2026-10-02):** the observed-incidence application
and focused tail refinement are complete. Production: 1,248 reporting rows,
936 computational sources, 7,598,000 independent outcome fits; separate tail refinement
adds 2,280,000 outcomes on the same allocations. All settings pass completeness
and the approved mean-MSE/coverage Monte Carlo gates; no warnings, aliases,
boundary fits or failures. Coverage remains below 95% and allocation tails remain
uncertain. See [completed findings](docs/findings/2026-10-02-nc-sud-application-findings.md).
Full chapter/appendix and derived CTJ/SI drafts are complete and verified; author review and submission metadata remain. 2026-10-11: chapter, SI (Table S12) and CTJ add a per-design mean/worst-decile comparison with SRS and scenario guidance (see ROADMAP progress log). See [manuscript completion findings](docs/findings/2026-10-02-manuscript-completion-findings.md).

**CTJ layout authority:** `CTJ_Manuscript.tex` remains the master for both layouts;
`CTJ_Reading.tex` only selects the two-column author preview. Edit captions/cells
in `CTJ_Exhibits.tex`, never duplicate manuscript content in the wrapper. Current
PDFs: submission 21 pp, reading 8 pp, chapter 68 pp, supplement 45 pp (2026-10-13, after the Chapter 3 length review). The author
approved original-title restoration and the six-exhibit revision. Main figures
are `figures/compact/` exports from `code/20_manuscript_figure_revision.R`; old
parameter/tail/budget graphics remain supporting material. Author order is
Walther, Habib, Simpson, Lin, with Lin corresponding. See
[exhibit revision](docs/findings/2026-10-02-ctj-exhibit-revision-findings.md) and
[layout review](docs/findings/2026-10-02-ctj-layout-review.md) for builds and checks.

**Full prelim/reference-flow checkpoint:** bios-dissertation's
`prelim/full-draft/prelim.pdf` is the 192-page working assembly (2026-10-13); Chapter 3 has 41 body
and 27 appendix pages after the length review. The master explicitly includes the author-requested
List of Tables (26 entries) and List of Figures (35 entries) even in prelim mode.
Rebuild there with `prelim/full-draft/tools/build_prelim.py` after
individual chapter syncs; source ownership and hooks are unchanged. Chapter 4
continuing work now uses narrative aims/expected evidence without a dated
timeline or advisor-meeting references. CTJ reading
uses the `balance` package (`\balance` before the bibliography) and no forced page breaks;
references begin on page 7 and occupy both columns on page 8. The review layout
alone loads `lastpage`. See
[checkpoint/teach-back](docs/findings/2026-10-02-prelim-assembly-ctj-reference-review.md).

**Current State (summary), as it stood in AGENTS.md:**

- Original revised grid: 12,800 eight-design scenarios; SRS benchmark extends to
  14,400 per estimator. September 24/27 current outputs, no checkpoint edits.
- Observed NC education application and tail extension: complete/verified;
  1,248 main/144 tail reporting rows, 9,878,000 distinct outcomes across stages.
- Full thesis chapter/appendix and derived CTJ/SI: complete, independently reviewed,
  rendered and numerically checked. Old application comments resolved by focused
  rewriting or agreed implementation, not by claiming unsupported facts confirmed.
- CTJ: standalone preceding-work citation verified; 3,165 prose/heading words
  (3,314 including end declarations)/279 structured
  abstract words, six exhibits, full detailed supplement and submission checklist.
- `code/20_manuscript_figure_revision.R` exports the compact figures and plot-data
  CSVs without simulation; `code/19_manuscript_exhibits.R` retains the prior detailed
  export. Run `Rscript code/tests/test_manuscript_exhibits.R` and
  `python3 paper/tools/verify_manuscripts.py` for numerical/exhibit consistency.
- Chapter checkpoint `58c1c82` synced through the approved hook to bios-dissertation
  `fc63326`; final manuscript checkpoint `a47b0ae` resynced to `6b07b95`
  (72 pp, 31 citekeys); editorial revision `3fc8385` resynced to `fa49f47`
  (72 pp, one compact figure copied); no push.
- Comprehensive four-track review completed and approved fixes applied (2026-10-12):
  `docs/reviews/2026-10-12-ctj-comprehensive-review-findings.md`. CTJ body 3,280 words
  (3,425 with declarations; limit 3,500), abstract 317.
- Advisor-meeting report: `application/report/IncidenceDesign_Application_Study_Report.{qmd,pdf}`
  (reads tracked aggregates only; inline numbers computed from CSVs).
- Next: author review; IRB statement, AI disclosure, coauthor approval (docs/TODO.md).
- Full figure selection, code teach-back, source mapping and live journal check:
  `docs/findings/2026-10-02-manuscript-completion-findings.md`.

- **Length review (adopted 2026-10-13):** the chapter is the reviewed concise version
  (prelim body 41 pp, was 45.5; Appendix B 27 pp). Detail moved to Appendix B; nothing
  numeric was dropped. Cut or moved text: `notes/trimmed-material.md`. The four annual NC
  tables are one merged `tab:nc-annual` in both chapter and SI. The verifier's `--chapter`
  and `--pending-si` options remain for checking any future candidate copy.

### Advisor meeting with Dr. Lin (2026-10-14)

- Reviewed the CTJ draft and application report. Main feedback: the outcome should
  be correlated with the incidence that guides allocation; test partial
  college participation, sample size/power and estimation-method sensitivity; group
  the designs in Table 1. Prelim: Project 1/2/3 structure, brief literature review.
- Record: `docs/meetings/2026-10-14-advisor-meeting.md`; action items in
  `docs/TODO.md` ("Advisor meeting 2026-10-14 follow-ups"). Only quick checks
  so far: equal-weight vs pooled rates make no meaningful difference, and the CTJ
  words were recounted. The outcome–incidence decision (possibly reporting β = 0 and
  β = 1 both) comes first.

### CTJ headings, abstract and advisor-meeting report polish (2026-10-14)

- CTJ subheadings rewritten to be informative (Methods: simulation study, estimating
  treatment effect and allocation risk, NC sudden unexpected death application; Results
  headings state the finding); Discussion split into Implications for trial design,
  Limitations, Continued work for the education pilot, Conclusions (author wording).
- Structured abstract rewritten for clarity; every number unchanged. Counts: body 3,458
  with declarations (limit 3,500), abstract 393 (limit 425). PDFs: submission 22 pp,
  reading 8 pp. Both verifiers pass. Chapter and SI unchanged.
- Application report (`application/report/`): black links, layout tightened (non-splitting
  tables, floating/resized figures), Appendix A removed, Figures 7-8 redrawn in-report
  from the same CSVs (shared exhibits unchanged); 12 pp.
- CTJ ethics placeholder (Methods; assumes not-human-subjects, IRB number/date blank)
  and coauthor-approval placeholder sentence added (DECISIONS.md); offset by a
  word-neutral Methods trim; Acknowledgements AI disclosure draft added. Body 3,326;
  3,526 with declarations; author: the 3,500 limit applies to the body.
  Placeholders highlighted yellow; AI disclosure states AI wrote code from author-specified
  methods, no AI-generated data, all code and results validated by the authors.
- Monte Carlo error and Monte Carlo standard error (MCSE) now defined at first use in the
  CTJ, SI and chapter (finite number of simulated datasets); later CTJ uses say MCSE or
  explain the source of error. "Treatment budgets" defined at first use as numbers of
  treated clusters (budget and allocation wording otherwise kept, author). Supplementary
  material statement rewritten. CTJ body 3,362, abstract 394; SI 45 pp, chapter 68 pp.

### CTJ and paper/ folder cleanup (2026-10-13)

- Build intermediates: a local `latexmkrc` sends latexmk's `.aux/.log/.bbl/...` to a hidden,
  gitignored `.build/`; symlinks to `sagej.cls`, `SageV.bst` and the shared bib replace the
  `TEXINPUTS`/`BSTINPUTS` setup, so plain `pdflatex` builds also work. Leftovers from plain
  builds: `git clean -fX -- paper/ctj_manuscript` (gitignored files only). Commit `8838311`.
- `paper/ctj_manuscript/` reorganized with `git mv`: `manuscript_ctj/` (CTJ_Manuscript,
  CTJ_Reading, CTJ_Exhibits), `supplementary_info/`, `cover_letter/`; `Figure_Legends.md`
  moved into `figures/`, which stays shared. TeX figure paths are now `../figures/...`.
  Verifier, exhibit test, companion-report links and current docs updated; dated plan docs
  keep their historical paths. Rebuilt PDFs matched the committed text (21/8/45 pp); the
  committed PDFs were kept. Commit `7d6558f`; both commits pushed.
- `paper/` tidied (`35f16d9`): `archive_manuscript/` renamed `archive/` with the April-July
  manuscript grouped into `pre_unification_20260702/` and the Oct 2 provenance JSON moved
  in; empty `manuscript/`, `.gitkeep` and the default presentation scaffold removed (slides
  are built in `bios-dissertation/slides/`); chapter figures: two current MSE figures kept,
  three July pre-revision figures removed (copies in `results/archive/pre_revision_20260924/`).
- Chapter folder and files renamed (`96a7b68`): `paper/thesis_chapter/Thesis_Chapter.{qmd,pdf}`
  (was `dissertation_chapter/Dissertation_Chapter`). Sync tools, scripts 20/21, verifier and
  current docs updated; hook reinstalled; its sync (bios-dissertation `4c20a57`) changed only
  path comments (73 pp). Dated logs and plans keep the old names.

### Concise Chapter 3 adopted (2026-10-13)

- Independent read-only review of the concise copy: 3 meaning drifts, 6 warnings, 7 notes;
  all fixed (queen-only scope, replication-schedule wording, regional-summary antecedent,
  tau/regime qualifiers, non-oracle rook contrast, citation placement, B.7 ordering).
  The merged annual table matched the four old tables cell for cell.
- CTJ Figure 2 (grid accuracy vs SRS, by configuration and regime) added to the SRS
  benchmark subsection at the author's request.
- Adopted: the copy replaced `Dissertation_Chapter.qmd`; SI's four annual tables replaced
  by `tab:nc-annual`; verifier passes without `--pending-si`. Hook synced to bios-dissertation
  `f4e36b0` (draft 76 -> 73 pp, 29 citekeys); full prelim rebuilt `81802a0` (194 -> 192 pp;
  Ch3 body 70-110, Appendix B 133-159; 26 tables, 35 figures); prelim tests pass.
  Standalone chapter 72 -> 68 pp; SI 47 -> 45 pp. The CTJ is unchanged and cites no SI table numbers.
- Both repositories pushed (SpatialCRT `f69a626`, bios-dissertation `81802a0`).
- Reference-list check (bios-dissertation): entries are longer than Eric Zhang's prelim because
  APA 7 prints DOIs and up to 20 authors (same font and spacing; Amber Young's dissertation
  matches ours). A first-six-then-et-al. variant saved 31 lines and no pages, so it was reverted.

### Concise Chapter 3 copy for length review (2026-10-12)

- Author-approved length review (`docs/plans/2026-10-12-chapter3-length-review-plan.md`)
  carried out on a **separate copy**, `paper/dissertation_chapter/Dissertation_Chapter_concise.qmd`
  (branch `ch3-concise`), awaiting author approval as the replacement; the canonical
  chapter, the hook sync, bios-dissertation and the SI are unchanged.
- Prelim-format preview (scratch build): Ch3 body 70-115 -> 70-109 (about 45.5 -> 39.5 pp;
  Ch2 body is 42 pp); Appendix B 25 -> 27 pp. Five batches: intro/methods trims
  (detail to B.1/B.2), merged queen results with test detail to B.3, rho/gamma/Checkerboard
  coverage tables to a new B.4, condensed rook/non-oracle/tau text, application procedures
  to B.8, four annual NC tables merged into one 36-row `tab:nc-annual`, two B.8 figures cut
  (in SI), Discussion restatements condensed, queen coverage and bias-variance figures to B.5.
  NC maps kept at full size (author). Cut text: `paper/dissertation_chapter/notes/trimmed-material.md`.
- Checks: no number in the copy that is not in the canonical chapter; verifier
  (`--chapter ... --pending-si tab:nc-annual`) and prelim assembly tests pass; no `??`.
- On adoption: copy over `Dissertation_Chapter.qmd`, swap the SI's four `tab:nc-full-*`
  tables for `tab:nc-annual`, rerun the verifier without `--pending-si`, let the hook sync,
  rebuild the full prelim, update page counts here, in AGENTS.md and in bios-dissertation.

### Comprehensive review and advisor-meeting report (2026-10-12)

- Four read-only review tracks (numbers, live journal compliance/clarity, chapter vs
  CTJ/SI, code + tests) found no numeric or code errors; findings and author decisions:
  `docs/reviews/2026-10-12-ctj-comprehensive-review-findings.md`.
- Approved fixes applied to CTJ, SI, chapter, exhibits, bib and cover letter:
  abbreviations defined at first body use; pipeline jargon rewritten; one design-name
  set across text, tables and figures (code/20 and code/21 labels regenerated); positive,
  accurate incidence-guidance framing in SI/chapter; SI data/ethics wording aligned with
  the CTJ (no repository link); two chapter overstatements qualified; Cai bib note and
  Baird DOI fixed; Lin telephone added; Codex/Claude mentions removed from code headers
  and documents. Habib citation unchanged (author).
- Counts: CTJ body 3,280 words, 3,425 with declarations (limit 3,500); abstract 317.
  PDFs: submission 21 pp, reading 8 pp, chapter 72 pp, SI 47 pp. Both verifiers pass.
- New advisor-meeting report: `application/report/IncidenceDesign_Application_Study_Report.{qmd,pdf}`.
- Open (TODO): IRB statement, generative-AI disclosure, coauthor approval.

### CTJ Figures 2 and 3 replaced (2026-10-12)

- Figure 2 adds a simple random sampling row (`grid_accuracy_srs.pdf`, from code/20).
- Figure 3 now shows the distribution of allocation-specific MSE relative to the SRS
  mean (5th-95th and 25th-75th percentiles within each setting, averaged over settings;
  `allocation_distribution_by_design.pdf`, from code/21). It is also in the chapter
  and SI (`fig:alloc-dist`). The previous bias/variance figure is dropped from the CTJ;
  SI Figure S4 already shows that decomposition, and the text cites it.
- Both new figures are vector PDFs in R's default Helvetica, like Figure 4. A 600-dpi
  Verdana PNG version was tried and replaced because its text looked soft in the PDFs.
- Counts: CTJ body 3,165 words, 3,314 with end declarations; abstract 279. PDFs:
  submission 20 pp, reading 8 pp, chapter 69 pp, SI 47 pp.

### CTJ exhibit and terminology revision (2026-10-12)

- CTJ main text and exhibits now show the six detailed designs and simple random
  sampling; the two related variants (plain saturation, balanced halves) are
  reported in the SI and chapter only. Table 1 rewritten with separate grid and NC
  rules; Table 2 caption defines MSE, coverage and the oracle fit; Figure 4 uses
  `nc_annual_accuracy_six.pdf`.
- Figure 1A shows pooled 2018-2021 rates (deaths / person-years); "frozen" is now
  "fixed" in all prose and figures. Tables and figures in the CTJ, SI and chapter label
  SRS as "Simple Random Sampling" (prose keeps the abbreviation). `code/20` also
  regenerates `figures/real_sud/regions.png`; application figure text relabeled in
  `application/code/real_sud_summary.R` (CSVs unchanged).
- `verify_manuscripts.py` checks the seven-row CTJ grid table and the new Table 1.
- Counts: CTJ body 3,177 words, 3,326 with end declarations (3,365 with math units);
  abstract 279. PDFs: submission 20 pp, reading 8 pp, chapter 68 pp, SI 46 pp.

### Per-design downside comparison and count update (2026-10-11/12)

- Added `code/21_downside_by_design.R` and
  `results/allocation_risk/summary/downside_by_design.csv`: each design's mean and
  worst-decile conditional MSE as ratios to SRS (grid pilot; refined NC corners).
- Chapter and SI: new Table `tab:downside` (SI Table S12) with a results paragraph.
  The chapter Discussion now leads with mean-plus-spread criteria and gives guidance
  by spillover scenario. The CTJ condenses the same point. The CTJ abstract
  conclusion describes incidence-guided saturation's performance. Fixed the CTJ
  typo (plain saturation 0.22329, was 0.32159). `verify_manuscripts.py` checks
  the new table against its CSV (commit `d1be46d`).
- Watson 2026 / Cai 2024 bib entries now match the verified bios-dissertation
  master bib (`b20179b`).
- Counts (texcount fragments, same method as before): CTJ body 3,304 prose/heading
  words, 3,470 with end declarations (3,508 also counting 38 math units); abstract
  285. PDFs: submission 21 pp, reading 8 pp, chapter 68 pp, SI 45 pp.
- 2026-10-12 author edits to the CTJ: removed the ethics placeholder, the telephone
  placeholder and the generative-assistance acknowledgement (acknowledgements
  omitted for now); added Walther's ORCID and an author-supplied contributions
  statement (AW, AH, RS, FL); Methods credits Habib for preparing the data. Data
  availability: simulation data/code on request, NC death-certificate data restricted,
  no repository link. Counts: body 3,263, 3,412 with declarations (3,450 with math
  units); abstract 285.

### Approved compact exhibit revision and author order (2026-10-02)

- Restored *Sampling Design for Spatial Cluster Randomized Trials Under Heterogeneous
  Incidence* in both CTJ layouts, SI and cover letter. Author order now Walther,
  Habib, Simpson, Lin; Lin remains corresponding.
- Six main exhibits: all-nine-design allocation rules and current grid MSE/coverage;
  observed NC rates/frozen-region map, compact simulation accuracy/decomposition,
  and nine-design annual NC MSE/SRS heatmap. Numbering follows first citation.
- Chapter uses the compact annual comparison; its substantive tail/budget findings
  remain in the body, detailed absolute/tail displays move to Appendix A7. SI mirrors
  this change. No simulation or analytical recommendation changed; old assets retained.
- New read-only export and intent-based checks verify exact empirical decomposition,
  configuration/regime separation and year-matched denominators; all manuscript
  tables/figures/citations checked. Counts: 3,215 body/283 abstract; reading 8 pp,
  review 21 pp, chapter 66 pp, SI 44 pp. Fresh rendered-page review and hook sync
  recorded in `docs/findings/2026-10-02-ctj-exhibit-revision-findings.md`.
- On author feedback about the 3,500-word cap, expanded CTJ from the chapter to
  3,215 body prose/heading words (3,381 with end declarations); abstract 283.
  Added mechanisms, effect-size crossover, omitted-spillover sensitivity, allocation
  risk interpretation and NC trial-planning implications; added numbers checked.
- Checkpoint `3fc8385` synced through the approved hook to bios-dissertation
  `fa49f47` (72 pp, 31 citekeys, one figure copied; Chapter 3 paths only).
- Next: author full draft review and requested Claude review; institutional metadata
  remain explicit placeholders. Nothing pushed or submitted.

### Paired CTJ reading and submission layouts (2026-10-02)

- Author approved building both layouts to assess the eventual article appearance.
  Created an 8-page two-column reading preview with embedded shared exhibits;
  retained the 18-page double-spaced submission version with exhibits at the end.
- One master owns prose, abstract and metadata; a small wrapper selects reading
  layout, and six shared exhibit definitions prevent drift in captions/cells/assets.
- Scientific body/abstract unchanged, all numerical/table/citation checks passed;
  both builds and all reading/submission pages visually reviewed. No new simulations.
- Clarified word-count convention: 2,184 prose/heading words; prior 2,217 also
  counted 33 math expressions. Abstract remains 283 words.
- Title restoration and less repetitive exhibit selection are proposed, pending
  author review. Preview charts deliberately retain current verified figures.
  Source/provenance/build details and teach-back:
  `docs/findings/2026-10-02-ctj-layout-review.md`. No push or submission.

### Complete chapter/appendix and derived CTJ/SI drafts (2026-10-02)

- Completed the observed-incidence application and same-allocation risk extension,
  then integrated full methods/results into the canonical Chapter 3 and appendix.
- Derived a standalone full CTJ article and detailed supplement from that source;
  replaced April/synthetic results. CTJ has 2,217 body words/283 structured-abstract
  words and six exhibits. Live journal requirements and preceding BMC citation verified.
- Independent reviewer cleared scientific wording after three scoped corrections;
  36 grid/20 primary MSE cells in all three sources, 360 annual cells/72 budget cells/
  32 matched-sensitivity cells and all 22 supplement tables reconciled.
- Fresh chapter/CTJ/SI renders and visual passes complete. Updated offline HTML,
  figure exports, declarations, submission checklist/cover letter, teach-back and
  comprehensive ready-to-paste Claude review prompt. Claude review remains future work.
- Author confirmed current authors/Lin corresponding and no funding/conflicts;
  author-approved IRB/data-use placeholder remains. Telephone/ORCID/contributions,
  coauthor approval and submission attestations remain author-side work.
- Project 1 unchanged/read-only; no audit or resolution claim. Restricted source
  data/caches remain ignored, older outputs preserved. No push or submission.
- Chapter checkpoint `58c1c82` synced successfully through the existing hook to
  bios-dissertation `fc63326` (72 pp, seven new assets, 31 citekeys). Final
  manuscript checkpoint `a47b0ae` resynced the coordinate wording to `6b07b95`
  (72 pp, 31 citekeys); no push.
  Detailed evidence: `docs/findings/2026-10-02-manuscript-completion-findings.md`.

### Completed observed-incidence application (2026-10-02)

- Completed production and the prespecified allocation-risk refinement; all
  1,248 main and 144 tail reporting settings verified. Combined distinct
  outcome effort is 9,878,000 fits, with copied tail prefixes counted once.
- Independent reviewer checked matching, actual estimator equivalence and all
  96 tail prefixes. Behavioral suites, independent cache recomputation and the
  saved companion's JavaScript checks pass.
- Primary control-only results favor saturation; both-arms SRS remains competitive.
  Undercoverage, budget differences and finite-tail uncertainty remain explicit.
- Traced publication exhibits and the production/refined HTML companion are built.
  Findings: `docs/findings/2026-10-02-nc-sud-application-findings.md`.
- Next required deliverables: full Chapter 3/appendix integration, fresh rendering
  and approved sync, then full standalone CTJ/SI and comprehensive Claude prompt.
  Author declarations confirmed except the authorized IRB/data-use placeholder.

### Observed-incidence implementation and visual companion (2026-10-02)

- Author approved the concrete plan, added a continuously updated HTML companion
  and requested a comprehensive Claude review prompt after the completed analyses
  and full thesis/CTJ drafts. Continue here after context compaction; no repeated
  approval needed. Old deadlines are superseded by ASAP completion.
- Setup validates named observations, pooled/annual cluster totals, cached
  nb/adjacency/W consistency and legal geometry. Regions are frozen across years;
  spatial blocks are frozen within each year. Balanced Halves now treats 29.
- New nested runner preserves allocation/outcome prefixes, handles duplicates
  and singleton precision, withholds complete comparisons on invalid boundary
  fits, and reports shared yearly provenance plus separate annual budgets.
- Behavioral suite and existing application data tests passed. Smoke: 116 rows,
  no incomplete settings or warning/failure fit rows; deliberately insufficient
  replication for conclusions. Larger pilot completed 928 rows/688 sources,
  no incomplete rows or warning records, and 524 expected precision misses.
  Production is running with approved refinement. Production verification, tail review,
  chapter integration and CTJ derivation remain unfinished.
- Offline companion saved with model flow, observed-year maps, example assignments
  and interactive performance filters. Maps reviewed as local images; in-app
  browser blocks file URLs, so interactive browser rendering is not yet verified.
- Details: `docs/findings/2026-10-02-nc-sud-implementation-findings.md`.

### SRS interpretation interview follow-up (2026-10-02)

- Author accepted the proposed primary queen grid, one-feature-at-a-time
  sensitivities (rank-based β = 1 baseline, regional summaries, rook corners),
  and starting Monte Carlo targets, while requesting an independent reviewer
  to check alignment with earlier simulations before implementation. Review
  completed read-only; no production code or runs yet. Verdict and safeguards are
  in `docs/reviews/2026-10-02-nc-sud-independent-review.md`: ready for implementation
  planning, no extra scientific sweep. Concrete plan saved as
  `docs/plans/2026-10-02-nc-sud-implementation-plan.md`, awaiting author approval.
- Reviewer verified singleton High Incidence Focus in all four observed years;
  start singleton production at R = 1,000 and refine using empirical MC targets.
  Stochastic designs need duplicate-aware joint MSE/coverage uncertainty and
  allocation-versus-outcome refinement. Primary Designs 1/3/4/9 are cross-year
  invariant with fixed W/regions and β = 0; disclose reused performance and keep
  annual population/incidence diagnostics separate. The β = 1 analysis is a
  matched baseline/adjustment sensitivity, not an isolated β comparison.
- Author superseded the date-based schedule: finish Project 2 ASAP, then derive
  CTJ/SI promptly for submission. Chapter-first hierarchy remains; disregard
  earlier November dates as scheduling constraints for this work.

- Author approved retaining the irregular-map adaptations and reporting their
  differing treatment counts/population shares. Spatial blocks use existing
  location/incidence-rank/population grouping, formed once per year and frozen
  within that year. Balanced Halves will treat 14+15, randomly choosing which
  half gets 15, to avoid the existing 28-treated rounding outcome. No code yet;
  remaining model/sensitivity/precision details and implementation plan pending.

- Author approved one fixed four-region map across all four yearly analyses,
  selected by the existing compactness and population/cluster/county-balance
  criterion using pooled person-years divided by four. Save before performance
  evaluation; yearly incidence changes saturation ranking, not region membership.

- Regional allocation summary approved: observed per-capita cluster rates →
  average ranks across 58 clusters → equal-weight mean rank within each fixed
  region → 80/60/40/20% saturation, with existing regional tie randomization.
  Explain explicitly in the application methods. This follows the grid Poisson
  version's mean of rank-normalized incidence, rather than mean raw rates.
- Author raised summary-method sensitivity. Added a focused proposal comparing
  mean cluster rates and pooled regional death/population rates on the same
  partition, first checking yearly orderings and then simulating changed
  incidence-guided rules. Scope awaits implementation-plan approval; no code/run.

- Outcome clarification: author accepted an education-focused NC primary model
  with observed incidence used for allocation, β = 0 in the outcome baseline,
  and an explicitly hypothetical β = 1 baseline sensitivity. Scale/extent and
  fitted covariates remain to settle before implementation. τ is constant within
  each scenario in both versions; βX does not make τ incidence-dependent.
- Read-only source check confirms both April and revised grid simulations, and
  the older synthetic application, included βX in outcomes and adjusted for X.
  Chapter 3 already states the revised grid equation and heterogeneous-baseline
  extension; the later manuscript pass should make the notation change from
  Chapter 2 and the distinct NC primary model explicit. No existing grid results
  or manuscript text changed.

- Author approved retaining the completed allocation-risk grid pilot as supporting
  evidence and using the real-SUD NC application as the next confirmation, without
  first expanding the grid comparison to all ten incidence surfaces.
- Approved the Chapter 2 linking sentence, recorded verbatim in
  `docs/plans/2026-10-02-nc-sud-application-plan.md`. The author does not need a
  direct cross-project numerical comparison in the body. A brief extension
  reference and inclusion of SRS/BSS can appear in the introduction/discussion;
  an appendix comparison is optional if necessary. In the follow-up the author
  approved retaining Checkerboard in Project 2's main comparison, with detailed
  identification diagnostics in the appendix. CTJ must stand alone and cite the
  verified Project 1 publication as the work being extended.
- Read-only checks passed: application data test suite; unique named cluster sets
  in every year and pooled period match cached queen/rook weights; finite rates,
  positive denominators and row-standardized weights verified. Restricted source
  and derived files remain ignored. No code or manuscript changes.

### Application design and integrated continuation (2026-10-02)

- Created `docs/plans/2026-10-02-nc-sud-application-plan.md` and an integrated
  `docs/prompts/2026-10-02-continuation-prompt.md`, as requested by the author. They combine
  Project 1 continuity, open manuscript decisions, and application execution.
- Author chose the existing continuous SAR approach, with observed NC SUD
  incidence informing allocation. Clarified after the author's follow-up:
  rank-based allocation does not require rate normalization. The outcome-model
  covariate scale and the regional saturation summary are separate choices.
  Remaining study details await agreement; no application code or manuscripts
  changed.
- Read-only input check: 58 finite cluster rates in each of 2018–2021 and the
  pooled period; named queen/rook weights available. Existing runner uses only
  both-arms spillover, suppresses fit warnings and has old cache/seeding behavior.
- September-session commits through `e4599ca` were pushed after explicit approval.
  New October work has no push authorization. Current Chapter 3 sync rules allow
  the established hook to update the prelim; Project 1 remains read-only.
- Next: finish the study interview and approve the concrete implementation plan,
  then implement/verify/run the revised application before replacing old results.

### Continuation note (2026-09-27): SRS and Chapter 2 decisions

Read `docs/findings/2026-09-27-ch2-check-findings.md` with `docs/TODO.md` before
continuing manuscript work. It records the restricted Chapter 2 reproduction and
subsequent review of the accepted PDF, including the allocation-level figures/tables.
Current code and stored CSV agree at 3×4 checkerboard MSE 0.22555, whereas the accepted
exhibit reports 0.0004; the exhibit's historical inputs have not been traced. This is
a limited reproduction discrepancy, not a settled explanation of cross-study results.
The author regards Project 1 as final: keep SpillSpatialDepSim and bios-dissertation
read-only. SRS framing, Checkerboard placement and the Chapter 2 link remain open;
no chapter/CTJ edits were made. Save future key findings in this project's `docs/plans/`.

**Author-confirmed reference rule (2026-09-27):** dissertation prose may refer directly
to "Chapter 2". CTJ must use self-contained wording (e.g., "In previous work...")
with a citation to the Project 1 BMC Medical Research Methodology paper, never
"Chapter 2" or "Project 1" as reader-facing references. Describe findings within
their studied conditions; verify final bibliographic details before adding the citation.

**New analysis requested (2026-09-27):** assess allocation-specific MSE variation
and upper tails for every candidate design and SRS before final recommendations.
Existing runs use one noise draw per allocation draw and cannot isolate this risk
from outcome noise. **The author subsequently authorized the necessary pilot/analysis
and subagents.** The nested pilot is complete; see
`docs/findings/2026-09-27-allocation-risk-findings.md` for results, precision checks and
reproduction commands. The older
"pending approval" statements in the chronological findings note predate this authorization.

**Allocation-risk pilot (2026-09-27):** new scripts 16–18 and isolated
`results/allocation_risk/` outputs; existing modules/main outputs unchanged. Queen,
oracle, tau=1; first two surfaces of five configurations, rho=0/0.5,
gamma=0.5/0.8, both regimes, all nine designs. Pilot: 100 assignment draws × 100
outcomes per unique allocation, 5,617,600 successful fits, no warnings/aliases.
IGSQ/SRS mean-MSE ratio = 0.649 under control-only, 1.070 under both-arms.
Control-only estimated worst-decile ratio = 0.555; finite-outcome noise requires
reading this with the targeted R=400 check in the findings note. Plain Saturation
Quadrants also performs strongly; incidence guidance is not shown uniquely optimal.
The R=400 check completed 50 selected blocks (1.5 million additional fits, no
failures/warnings/aliases), preserving the regime distinction; its control-only
IGSQ/SRS mean and estimated worst-decile ratios are 0.658 and 0.527. Individual
tail membership remains noisy. Full pilot plus refinement: 7,117,600 distinct fits.
Never substitute these selected-surface/corner numbers for the full-study numbers.
Manuscript framing/Checkerboard placement remain author decisions; revised real-SUD
application remains pending. Project 1 and bios-dissertation remain read-only.

### Current State (as of 2026-09-25): real SUD data aggregated

**The real NC SUD data are aggregated to the 58 community college clusters.** The design
comparison on them is the next plan. Plan of record:
`~/.claude/plans/read-the-prompt-at-gentle-willow.md`.

- **SUD data details** (sources, numerator decision 23,523 vs 21,147, rate formulas, outputs,
  contiguity weights, code layout, Habib script notes): see [application/AGENTS.md](application/AGENTS.md).
- **Chapter Application text (updated 2026-09-25, commit eb4db79):** framed as a proposal.
  Data preparation (aggregation, yearly rates, weights, checks) is reported as done; the design
  application is planned. A hidden HTML-comment TODO stub lists the planned exhibits. Open
  `NEEDS-AUTHOR-CONFIRMATION` comments sit inline, including a possible error: the text says
  "zip-code pooling", but the data are keyed by county code. The application runner itself is
  still pre-revision (`lagsarlm`, old seeding).

#### Phase D chapter rewrite (2026-09-24/25): steps 0–4 and 7 done; chapter synced to the prelim
- `paper/dissertation_chapter/Dissertation_Chapter.qmd` is rewritten on the revised run:
  - Methods (from the spec), Results (robust vs. fragile; queen primary; rook sensitivity;
    non-oracle; τ stability), Discussion, and a brief Chapter 2 link through its β̂ ≈ β − ψ
    collinearity result.
  - Static Table 2; queen/rook figure pairs (Figures 1–4, including bias–variance).
  - Renders at 40 pp.
- Six fresh-reviewer passes. Passes 4–6 had no ERROR outside the Application.
  `docs/prompts/2026-09-23-step0-chapter-review-prompt.md` is the current review brief; it holds the
  settled decisions.
- **Finding to remember:** under queen, the control-only lead of Incidence-Guided Saturation
  Quadrants holds only up to τ = 1.5. Isolation Buffer matches it at τ = 2 and passes it at
  τ = 3. The recommendation rests on lowest pooled MSE and first/second in each regime at
  every τ.
- **Step 3 draft for author review:** `paper/dissertation_chapter/notes/step3_exhibit_fate.md`.
  It covers body vs. appendix, the appendix outline A1–A6 (A5 = oracle vs non-oracle table; A6 = technical notes on
  anticipated committee questions), and a new deliverable: Q&A backup slides plus a written
  companion doc for the full prelim.
- **Figure code:** `code/12`/`14` got cosmetic fixes (plotmath τ/ρ, captions, legend order).
  Numeric outputs are unchanged.
- **Bib:** `paper/SpatialCRT_IncidenceDesign.bib` is synced to the master. Printing fields
  absent from the master were dropped, so the CTJ loses publisher cities (revisit at step 5).
  The five test citations were added to `bios-dissertation/prelim/references.bib`
  (local commit b1e90b3, not pushed).
- **Step 3 applied (author-approved, 2026-09-25):**
  - body adds Table 3 (MSE by configuration), Figure 1 (design samples) and Figure 6 (service areas);
  - the appendix has A1 (eight designs), A2 (simulation mechanics), A3 (CD diagrams, coverage
    across τ) and A4 (technical notes from a full re-run, `notes/a5_technical_notes.md`);
  - the final fresh review had no ERRORs, and its WARNs are fixed. Standalone render: 56 pp.
- **Step 7 / sync (author decision: the SpatialCRT `.qmd` is the only file edited):**
  - `paper/dissertation_chapter/tools/sync_to_prelim.sh` generates
    `bios-dissertation/prelim/project-proposals/project2-incidence/draft/project2-incidence-draft.qmd`
    (CHAPTER 3 / APPENDIX B). It checks every citekey against the master bib, failing loud,
    renders, and commits only `project2-incidence/` paths there. It never pushes.
  - A `post-commit` hook (install with `tools/install_hooks.sh`) runs it when a commit touches
    the chapter, its figures, `tools/` or the bib. Its log is `tools/sync.log` (gitignored).
  - Cross-references are `\label`/`\ref`, so numbers read 3.x / B.x in the prelim.
  - Prelim render (step 4): body 41 pp, appendix 13, references 3; 59 pp with the TOC preview.
  - The class (`bios-prelim.cls`) handles both pandoc `CSLReferences` forms since 2026-09-25,
    so the header no longer patches it. `tools/prelim_header.yml` sets `colorlinks: false`
    (2026-10-01): Quarto 1.10 otherwise prints links blue.
- **Author review of the prelim render (2026-09-25), applied (commit 02a70d9):**
  - queen figures in the body; the rook MSE/coverage/bias-variance figures are in Appendix A4;
    the τ figure is a single side-by-side queen|rook plot;
  - floats are `[!htb]` with smaller figures, so none lands alone on a page or after the
    chapter end;
  - A2 has no script or function names;
  - the technical notes (now A6) are one formal subsection, and the regime-gap note moved to
    the Q&A notes;
  - SUD counts are 23,523, with the heart-failure inclusion explained.
  - Standalone 53 pp; prelim 56 pp.
  - The design palette ends at viridis `end = 0.85` so Checkerboard is visible; the colours
    are consistent across the τ and coverage figures.
  - The sync skips PDF-only re-renders.
  - The Chapter 2 audit is paused by the author. Chapter 2 is accepted and final: any critique
    of it stays light.
- **Pending:**
  - author confirmation of the inline `NEEDS-AUTHOR-CONFIRMATION` items;
  - a Chapter 2 audit (separate session), because Chapter 2's 3×4/3×3 block-stratified numbers
    contradict its own β̂ ≈ β − ψ.
- **TO DO before the CTJ (author, 2026-09-25): add Simple Random Sampling as a benchmark design.**
  There is currently no naive baseline: Checkerboard continues Chapter 2 but isn't a generic
  benchmark. Plan it in a dedicated session. Open choices:
  - complete randomization with exactly 50 treated (recommended) vs Bernoulli;
  - SRS as a 7th ranked design vs a separate reference;
  - a relative-efficiency column (MSE / MSE_SRS).

  Key-seeded draws (spec §4: the Z key includes d; ε and X are shared per block/config) mean
  existing designs' results should be reproduced exactly by the rerun, which makes a free
  verification check. It touches:
  - `03`, `05` (re-run, ~14 min);
  - summaries `12`–`14` (rank-based statistics change: Friedman, Nemenyi, win rates, CD
    diagrams);
  - chapter numbers and figures;
  - the application study.
- **Next:** CTJ derivation (manuscript step 5); integrate the application results as they land.

### Prior State (as of 2026-09-24)

**Simulation revision (step 0.5): re-run COMPLETE; downstream regeneration (Phase C) in progress.**
Plan: `docs/plans/2026-09-24-simulation-revision-plan.md`; method authority:
`docs/plans/2026-09-24-simulation-revision-spec.md`. All April 2026 numbers are superseded (archived in
`results/archive/pre_revision_20260924/`, with a README).

- **What changed:** M1–M8. Matched surfaces; 100,000 people per Poisson cluster; random
  tie-breaking, with High Incidence Focus and Balanced Quartiles treating exactly 50 (the
  latter a user decision after the pilot); key-based seeds; one noise column per fit;
  aliasing flagged (estimates kept); a non-oracle estimator (Y ~ Z + X) as sensitivity; MC SEs
  from 10 surface means; manifest-checked checkpoints; the lean `fit_sar_lag()` engine
  (validated on 5,120 fits, max |Δτ̂| 9.7e-8, ~120× faster).
- **Full run:** `results/sim_data/sim_results_{MLE,MLEnonoracle}_tau_sweep_combined_20260924_025509.rds`.
  12,800 scenarios × 250 fits per estimator (6.4M fits, 13.9 min on 10 workers).
  `full_run_verification.txt`: all checks pass; zero non-aliasing warnings. Aliasing hits
  every Checkerboard × rook fit, plus 20 Isolation Buffer × rook scenarios where one draw
  happened to be the exact checkerboard. The 1% `lagsarlm` cross-check is in
  `results/estimator_validation/crosscheck_full_run.txt`.
- **Headline (oracle, τ = 1, 6 designs):**

  | Design | Queen MSE | Rook MSE |
  |---|---|---|
  | Incidence-Guided Saturation Quadrants | 0.091 | 0.072 |
  | Balanced Quartiles | 0.131 | 0.084 |
  | Isolation Buffer | 0.160 | 0.133 |
  | High Incidence Focus | 0.245 | 0.206 |
  | 2x2 Blocking | 0.319 | 0.128 |
  | Checkerboard | 1.087 | 0.480 |

  - Coverage is ≈0.94 for every design except Checkerboard × rook (0.15; τ not identified).
  - Under queen the rank order is the same at every τ (under rook, 2x2 Blocking and Isolation Buffer swap at τ = 0.8).
  - Every adjacent pair is significant on 50 independent config × surface units.
  - Saturation Quadrants and Balanced Halves remain indistinguishable from their retained
    counterparts, so the 6-design set is kept (user decision, 2026-09-24).
  - Non-oracle: bias −0.12 to −0.40 and coverage 0.51–0.90 for all designs, but lower
    variance where Z and WZ are collinear (queen Checkerboard MSE 1.09 → 0.11).
- **Research framing (user, 2026-09-24):** the question is which designs estimate τ well
  under heterogeneous incidence, NOT whether knowing incidence helps.
- **Outputs regenerated:**
  - `results/six_design_manuscript/{,queen/,rook/}`, `six_design_manuscript_nonoracle/`,
    `eight_design_supplementary/`
  - `results/MLE_tau_sweep_*.pdf`, `results/MLE_statistical_comparisons.pdf`, `mle_per_config/`
  - `results/figures/design_samples_8panel.*`
- **Not yet updated:** the CTJ manuscript and SI still carry the April numbers. They stay
  untouched until manuscript step 5 ("numbers superseded"). The application is not re-run
  (14(d)'s table is labeled STALE).
- **SUD data attribution (2026-09-24):** the Application passages of the CTJ, SI (S11), and
  chapter now credit Habib's master's paper (`habib_temporal_2026`) for the county-level
  data: working-age (18–64) sudden unexpected out-of-hospital deaths, 2018–2021, 21,147 of
  412,514 NC deaths (111,665 working-age), via the SUDDEN-validated algorithm
  (`gan_factors_2019`, `nanavati_sudden_2014`); we aggregate them to the 58 CC clusters. This
  replaced the wrong "~100,000 death certificates" and "epidemiology co-investigator" wording.
- (Superseded: Phase C and Phase D are done; see Current State above. The single
  recommended design for investigators is still the author's decision to make.)

### 2026-10-02 — full prelim assembly and CTJ reference flow

- Approved full-length prelim assembled in bios-dissertation: 177 pages, four
  chapter bodies (24/42/42/18), appendices A/B (4/25), one combined bibliography.
  Existing source prose/math/captions preserved; five assembly tests pass.
- Author retains Chapter 3’s full appendix pending content review. New initial
  front abstract is 299 words. Current literature-review working revisions are
  included without staging or overwriting them.
- CTJ reading references flow from page 7 across both final-page columns; eight
  pages retained. Submission remains 21 pages with identical extracted text.
  Numerical/exhibit checks and visual review pass; no scientific result changes.
- Findings/teach-back: `docs/findings/2026-10-02-prelim-assembly-ctj-reference-review.md`.
  No pushes; full-draft review remains.

### 2026-10-02 — full prelim exhibit lists

- Added the author-requested List of Tables (29 entries) and List of Figures
  (35 entries) in bios-dissertation, including appendix exhibits and contents
  entries. Full prelim now 187 pages; chapter-body lengths and printed numbers
  unchanged. All 64 page references checked, all ten list pages visually
  reviewed, five assembly tests and manifest PDF hash check pass.
- Front-matter master edits only; no chapter/template/scientific or CTJ changes.
  Pending literature-review/bib revisions remain unstaged; no push.

### 2026-10-02 — Chapter 4 narrative continuing-work revision

- At author request, removed the fixed-date timeline from canonical SSBMF
  chapter source, retained cut material in multiomicsGEP notes, verified hook
  sync and rebuilt the full prelim. Planned advances are described directly
  in prose informed by the October 2 plan, without meeting references.
- Full prelim remains 187 pages; Chapter 4 body 18 pages. Lists now have
  28 tables/35 figures and all 63 references verified. Five assembly tests,
  PDF hash and visual review pass; other chapters/appendices/reference page
  text unchanged. No Project 2 scientific or CTJ changes; no push.
