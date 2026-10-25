# Ready-to-paste Claude review prompt

The review below has not yet been performed by Claude. Paste it with this project
open. It requests assessment first, not unilateral new studies or publication.

---

Please conduct a comprehensive independent review of the completed Project 2
application study and full thesis/CTJ drafts in:
`/Users/ajwalther/GithubProjects/SpatialCRT/projects/IncidenceDesign`.

The intended endpoint is completed analyses and full thesis chapter/appendix plus
a standalone Clinical Trials (CTJ) manuscript/supplement derived from that chapter.
Plans and pilots are intermediate stages. These full drafts now exist; assess
whether they meet the scientific goals and whether anything material is still
missing. Earlier deadline dates are superseded: finish promptly, with CTJ still
required. Start read-only; do not change code, documents, checkpoints, or analytical
choices until you present the review and the author authorizes corrections.
Do not push, submit, send messages, or edit Project 1/bios-dissertation.

Read in order:
1. Root `../../AGENTS.md`, project `AGENTS.md`, `ROADMAP.md`, `docs/plans/TODO.md`.
2. `docs/plans/ch2_check_findings_2026-09-27.md` and
   `allocation_risk_findings_2026-09-27.md`.
3. `docs/plans/nc_sud_application_plan_2026-10-02.md`,
   `nc_sud_implementation_plan_2026-10-02.md`,
   `nc_sud_independent_review_2026-10-02.md`,
   `nc_sud_implementation_findings_2026-10-02.md`,
   `nc_sud_application_findings_2026-10-02.md`.
4. `docs/plans/manuscript-unification-plan.md` and
   `manuscript_completion_findings_2026-10-02.md`,
   `ctj_layout_review_2026-10-02.md` and
   `ctj_exhibit_revision_findings_2026-10-02.md`.
5. `application/AGENTS.md`, `application/README.md` and project `README.md`.
6. The canonical `paper/dissertation_chapter/Dissertation_Chapter.qmd` and its PDF,
   `paper/ctj_manuscript/CTJ_Manuscript.tex`/PDF, shared `CTJ_Exhibits.tex`,
   the `CTJ_Reading.tex`/PDF layout preview and
   `Supplementary_Information.tex`/PDF, plus submission checklist, figure legends,
   cover letter and updated manuscript declarations.
7. The saved `application/report/real_sud_companion.html` and current output
   manifests/verification/exhibits; check Git status before interpreting files.

Protect restricted data. County source files and `application/data/derived/`
are ignored and must never be committed, uploaded or copied into the review.
Report aggregate diagnostics only; authorized 58-cluster aggregates/maps may be
reviewed. Fit caches are ignored. Work in the existing checkout. Project 1 is
accepted/final/read-only: do not reopen its numerical audit or claim its historical
reproduction discrepancy has been resolved. Bibliographic verification alone is
appropriate. The verified preceding-paper DOI is 10.1186/s12874-026-02968-0.

Review the implementation against the original grid mechanisms:
- Read `code/01_spatial_setup.R` through `05_run_simulation.R`, the validated
  `fit_sar_lag()` engine in `code/04_estimation.R`, grid/SRS/risk verification scripts, and
  `code/03_designs.R` in detail. Then inspect `application/code/application_designs.R`,
  `real_sud_setup.R`, `real_sud_simulation.R`, `run_real_sud.R`,
  `refine_real_sud_tail.R`, `real_sud_summary.R`, and the companion renderer.
- Confirm observed 2018–2021 incidence is fixed on 58 named clusters; no primary
  synthetic incidence is generated. Primary education outcomes are
  Y=(I−ρW)⁻¹[τZ+S(Z)+ε], β=0, τ=1, σ=1, fitted with intercept/Z/true Spill and no X.
  Incidence directs allocation only; it is not education ability or an assumed
  moderator of the intervention. τ is a structural coefficient, not a total
  spatial impact. Both spillover regimes are correctly generated and fitted.
- Distinguish the original grid β=1 model, primary application β=0 model and
  hypothetical β=1/X-adjusted sensitivity. The latter changes baseline and
  adjustment jointly; no empirical incidence–knowledge link is claimed.
- Queen is primary. Rook uses matched parameter corners, not comparison with
  unconditioned broader queen averages. Named cached legal-boundary weights,
  cluster order and raw/standardized matrices must agree. Geographic adjacency is
  an exposure proxy, not a measured professional communication network.
- Check every one of eight designs plus SRS: graph Checkerboard rather than
  rectangular coloring; greedy random independent-set buffer; rank/tie rules;
  29-treatment SRS/HIF/balanced assignments; 14+15 Balanced Halves; frozen four
  saturation regions; 14 frozen yearly spatial blocks including singleton rounding.
  Adaptations retain mechanisms as closely as possible; do not silently equalize budgets.
- Application guided saturation uses annual per-capita cluster rates → average
  ranks → equally weighted regional mean ranks → 80/60/40/20%. Original grid
  Design 8 averages X directly; mean rank corresponds only to rank-scaled Poisson X,
  not continuous iid/spatial grid X. Mean raw-rate and pooled death/population
  summary sensitivities retain regions and do not select the primary by results.
- Confirm frozen regions, yearly blocks and reference population were chosen
  before outcome performance, annual settings kept separate, treatment counts and
  population shares reported alongside estimation. Population shares are geographic
  reach, not participant enrollment or equal implementation cost.

Audit precision and numerical evidence:
- Production: 1,248 reporting rows, 936 computational sources and 7,598,000 independent
  outcomes. Refined tail: 144 rows, 96 sources, 3,060,000 outcomes including 780,000
  copied prefixes and 2,280,000 new; combined distinct outcomes 9,878,000.
  Computational-source count is not necessarily a count of distinct distributions.
- Check all completeness/precision gates, finite coefficients/SEs, fitted-matrix
  rank/residual treatment information, warnings, aliases, boundary fits and failures.
  A zero warning/failure claim applies to the application, not the intentionally
  aliased rectangular-rook Checkerboard grid cases.
- Assess randomization-key reproducibility, manifest refusal, matching observed
  inputs, independent outcome replication, duplicate-frequency covariance, singleton
  support handling, preserved main prefixes and explicitly shared yearly sources.
  Reused annual performance is not independent confirmation.
- Assess bias/MSE/coverage and MCSE formulas against correct conditional expectations,
  not only output shape. Main stochastic settings use J100/R100, singleton R1000,
  two hypothetical-sensitivity HIF sources R2000. Maximum relative mean-MSE MCSE
  0.0487115 and coverage MCSE 0.00891687 satisfy approved gates. Coverage about
  92.6–93.6% remains below nominal; precision is not interval calibration.
- Allocation-specific m(a)=Eε[(τhat−τ)² | a,X] is different from realized squared
  errors and scenario-level MSE. Verify quantiles/worst-decile summaries, signed
  noise-corrected variance, singleton qualifications, cross-selected outcome halves
  and their finite-outcome interpretation. Same-allocation R400 refinement is not
  fresh allocation sampling. J100 yields only ten sampled worst-decile draws;
  all 112 non-singleton tail reporting rows flag uncertain boundary membership.
- Test the conclusions rather than assuming the preferred design wins: control-only
  plain saturation average MSE is 27.4% lower than SRS; guided reduction is 21.1–29.9%
  by year. Both-arms SRS remains competitive; no uniqueness/equivalence or additional
  incidence-guidance benefit is established. Conditional deterministic HIF accuracy
  does not demonstrate randomized causal identification.
- Confirm grid allocation-risk pilot values remain separate supporting evidence;
  fair representation of preceding BSS work includes reasonable estimation and
  limitation of poor-allocation downside within its studied settings. No direct
  numerical cross-project comparison is required.

Read current outputs below `application/results/real_sud_rev_20261002/`, including
production/tail verification records and all `exhibits/` CSVs. Run proportionate
read-only checks as needed. Existing commands:
`Rscript application/tests/test_application_data.R`,
`Rscript application/tests/test_real_sud.R`,
`Rscript application/tests/test_real_sud_tail.R`,
`Rscript application/tests/test_real_sud_summary.R`,
`node application/tests/test_companion.mjs`,
`Rscript code/tests/test_manuscript_exhibits.R`,
`python3 paper/tools/verify_manuscripts.py`.
Set VECLIB_MAXIMUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before R.
The exhaustive cache verification and full ML cross-check were already run;
inspect their records before repeating expensive work. Do not rerun production,
change manifests or start new simulations without author approval.
The DOM fixture checks the actual HTML script's behavior, not browser appearance;
in-app local-file browsing was blocked, so no visual-browser check was claimed.

Review all document sources AND rendered PDFs:
- Does Chapter 3 clearly extend Chapter 2 while explaining β/τ/model differences,
  and retain substantive full-study/application findings in its body? Are appendix
  details appropriately supporting? Checkerboard stays in main comparisons and
  detailed diagnostics in the appendix. Distinguish the detailed six-design ranking
  subset from full nine-design benchmark/application.
- Were unsupported old NEEDS-AUTHOR comments resolved through removal/rewording or
  actual author decisions, rather than falsely claiming verification? Check the
  recorded resolution table and corrected numerator/Habib provenance.
- Was CTJ trimmed/reorganized from the completed chapter, without a competing
  account or stale April/synthetic application claims? Verify all 22 supplement
  tables and source references, citations, figure values/captions, annual/sensitivity
  parameter matching and every material prose number. No DIM discussion in manuscripts.
- CTJ/SI must stand alone: no reader-facing “Chapter 2” or “Project 1”. Verify the
  BMC publication citation against its official publisher page, without reopening
  accepted-study findings or inventing volume/page details.
- Verify live journal requirements: original-research body ≤3,500 words, structured
  abstract ≤425 with appropriate headings, six main exhibits, keywords, review
  formatting/title page, legends, reference style, declarations and submission parts.
  Draft counts are 3,215 prose/heading words; 3,381 including end declarations
  (3,419 also counting 38 mathematical-expression units)/283 abstract; two tables/four figures. Old plan's
  250-word unstructured-abstract instruction was superseded by the live check.
- Author confirmed Walther, Habib, Simpson and Lin in that order with Lin corresponding, currently
  no funding/conflicts. IRB/data-use/consent documentation is unavailable and the
  author explicitly approved a placeholder. Do not fabricate an approval, waiver,
  protocol or permission. Telephone/ORCID/contributions/coauthor consent/final
  attestations remain submission tasks. The implementation is not an actual trial.
- Assess the Codex assistance-disclosure draft against current Sage policy and
  author responsibility. Do not treat a future Claude review as already performed.
- Verify the approved six-exhibit revision: all-nine-design allocation rules and
  grid MSE/coverage tables; observed 2018/frozen-region map, compact six-design
  grid accuracy/decomposition, and all-nine-design annual NC ratio heatmap.
  Numbering follows first citation (map first). Check the read-only export
  `code/20_manuscript_figure_revision.R` and plot-data CSVs: configuration/regime
  separation, exact MSE=Bias²+(249/250)SD², and year/regime-matched SRS denominators.
  Do not confuse ratio-of-means with mean setting-specific ratios. Chapter body
  retains substantive tails/budgets; detailed displays are in Appendix A7/SI.
- Inspect figures/tables/math/citations/page flow at readable size. Build success
  alone is insufficient. Confirm approved Chapter 3 sync occurred through its hook,
  no direct bios-dissertation edits, old outputs preserved, restricted sources not
  tracked and no unapproved push. Read-only `git status/log` is appropriate.

Return:
1. A concise verdict: completed scientific/drafting goal, mandatory blockers,
   remaining author-only submission tasks, and whether an additional analysis is
   actually necessary. Do not equate IRB placeholders with incomplete simulations.
2. Findings by severity, each with file/line/output evidence, why it matters and
   the smallest verifiable correction. Separate scientific correctness from
   numerical provenance, presentation and optional extensions.
3. A cross-document consistency/number-and-exhibit audit and reviewer-check record.
4. Any genuinely unresolved decision requiring the author's input. Do not reopen
   settled preferences or speculate about improvements as if mandatory.
5. Prioritized correction plan for author approval, with checks that would establish
   completion. If no mandatory scientific/drafting issue remains, say so plainly.
