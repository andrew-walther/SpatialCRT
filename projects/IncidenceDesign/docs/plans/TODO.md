# IncidenceDesign to-dos — completed study and drafts

Updated 2026-10-13. The author superseded earlier deadline dates and requested ASAP
completion with CTJ remaining a required deliverable. Historical interview and
execution to-dos are preserved in Git history; settled choices are not reopened.

## Completed

- [x] SRS benchmark added and verified without altering original eight-design results.
- [x] Average MSE and allocation-specific downside interpreted together; grid pilot
  distinguished from the full study and finite tail uncertainty acknowledged.
- [x] Checkerboard retained in main comparisons; technical diagnostics in appendix.
- [x] Approved conceptual Chapter 2 link; CTJ standalone preceding-work citation.
  Project 1 is final/read-only; historical numerical discrepancy remains unresolved.
- [x] Author-approved education application: observed 2018–2021 inputs, 58 clusters,
  β=0 primary, constant τ, no synthetic incidence, all eight designs and SRS,
  queen primary/both spillover regimes/yearly reporting, geographic adaptations,
  unequal budgets and focused sensitivities.
- [x] Independent alignment review and concrete implementation approval.
- [x] Implementation, behavioral/ML checks, manifests, named cached weights,
  reproducible seeds and explicit warnings/alias/failure accounting.
- [x] Full production: 1,248 rows/936 computational sources/7,598,000 independent
  outcome fits; all complete and mean/coverage MC gates passed; no invalid fits.
- [x] Separate same-allocation tail extension: 144 rows/96 sources, 2,280,000
  additional outcomes; tail-membership uncertainty remains explicitly reported.
- [x] Offline HTML companion: development, model, maps, mean/tail performance,
  budgets, yearly filters and links to completed drafts; DOM behavior tested.
- [x] Full thesis chapter and appendix with number/exhibit provenance, resolved
  former application comments and fresh PDF review; approved Chapter 3 sync.
- [x] Full CTJ manuscript/supplement derived from the chapter; old results replaced,
  live journal requirements checked, confirmed declarations and explicit ethics
  placeholder, figure selection, legends and cover-letter draft.
- [x] Independent scientific source review, aggregate-to-table verification,
  fresh document builds and visual review.
- [x] Requested comprehensive Claude review prompt saved.
- [x] Author-requested paired CTJ layouts: shared-source two-column reading preview
  and double-spaced submission draft, rendered and numerically checked.

## Remaining author/submission work

- [x] Assemble full-length prelim in bios-dissertation, retaining the full
  Chapter 3 appendix; fix CTJ reading bibliography flow/last-page columns.
  See [assembly/reference review](prelim_assembly_ctj_reference_review_2026-10-02.md).
- [x] Add author-requested List of Tables and List of Figures to the full prelim;
  verify all 63 current exhibit page references and preserve main-matter numbering.
- [x] Remove Chapter 4's fixed-date timeline in canonical multiomicsGEP source,
  sync and rebuild the prelim; frame continuing work as methodology and expected evidence.
- [ ] Author review of full prelim and initial front abstract. Chapter 4's later
  substantive revision will follow its October 2 plan and verified new results.
  Rebuild the assembly after later chapter syncs.

- [x] Implement approved title restoration and six-exhibit revision in both CTJ
  layouts; carry applicable clarity changes into the chapter and derived supplement.
  See [exhibit review](ctj_exhibit_revision_findings_2026-10-02.md).
- [x] Chapter 3 length review (2026-10-12/13): concise version independently reviewed
  and adopted (prelim body 45.5 to 41 pp; Appendix B 27 pp; standalone chapter 68 pp;
  SI 45 pp with merged annual table); synced, full prelim rebuilt (192 pp), pushed.
  See [plan/status](chapter3_length_review_plan_2026-10-12.md); cut text in
  `paper/thesis_chapter/notes/trimmed-material.md`.
- [x] Reference-list length check: strict APA 7 with DOIs and full author lists kept
  (a six-author et-al. variant saved no pages; reverted).
- [ ] Author full review of the revised (concise) chapter, both CTJ layouts and supplement.
- [x] Comprehensive review (four tracks: numbers, journal compliance/clarity,
  chapter consistency, code) completed 2026-10-12; approved fixes applied. See
  [findings](ctj_comprehensive_review_findings_2026-10-12.md).
- [x] Application study report for the advisor meeting:
  `application/report/IncidenceDesign_Application_Study_Report.{qmd,pdf}`.
- [ ] **Ethics/IRB statement** (CTJ requires one in Methods): placeholder added
  2026-10-14 (assumes not-human-subjects; `[IRB #XX-XXXX, date]`). Confirm the
  determination with Habib/Simpson, fill in number/date (also on the title page),
  and confirm permission to publish 58-cluster aggregates and maps.
- [ ] **Generative-AI disclosure** (Sage policy requires one for AI-assisted code,
  code checking and visualization): draft in CTJ Acknowledgements (2026-10-14);
  author to confirm tool list and wording.
- [ ] Coauthor ORCIDs if the submission system requests them (Walther's is on the
  title page). Lin's telephone (919) 966-9431 added from his UNC profile.
- [ ] Confirm contribution wording; obtain all coauthors' final approval and
  submission attestations; confirm the pilot description with Simpson.
- [ ] Submit CTJ with the author-completed
  [checklist](../../paper/ctj_manuscript/submission_checklist.md)/cover letter.
- [x] Completion work is pushed: `main` matches `origin/main` (verified 2026-10-13,
  SpatialCRT `f67ab4b`, bios-dissertation `7d10596`). Future pushes still require explicit author authorization.

No additional analysis was required by the final independent source review.
Actual pilot endpoint/enrollment/budget/network specification and interval
calibration belong to planning the implemented trial; this simulation evaluation
is complete and does not claim observed intervention effectiveness.

Findings: [application](nc_sud_application_findings_2026-10-02.md),
[manuscripts and teach-back](manuscript_completion_findings_2026-10-02.md).
Commands and safeguards: project/application README and AGENTS.

## Known minor code issues (low priority; moved from AGENTS.md 2026-10-14)

- `02_incidence_generation.R`: `generate_incidence_poisson()` still defaults `pop_per_cluster = 1000` (and its roxygen says so). 05 always passes 100,000, so results aren't affected. It was left unchanged after the full run so the run manifest's code hash still matches.
- `DESIGN_FULL_NAMES` in `00_design_names.R` (used only by the two design-sample figures) labels Design 1 "Block Stratified Sampling (Checkerboard)" (user decision 2026-09-24: Chapter 2's name plus this study's). The manuscripts, 03 and 12–14 say "Checkerboard", and `DESIGN_SHORT_NAMES` still says "Block Stratified". This is not a bug to fix: the author uses "Block Stratified Sampling" and "Checkerboard" interchangeably, and the chapter says so once (Phase D).
- `load_latest_results()` comment (line ~87 of `06_visualizations.R`): clarify `_combined_` preference applies per estimation-mode, not globally
- `07_results_summary.Rmd` compare-table caption: should note that DIM only ran 6 designs (D7/D8 NAs are expected)
- `IncidenceSpatialCRT_Report.qmd` caption/text alignment: MC SEs table uses tau=1.0 slice (2,560 rows), not full 12,800 — prose now correctly clarifies this distinction
