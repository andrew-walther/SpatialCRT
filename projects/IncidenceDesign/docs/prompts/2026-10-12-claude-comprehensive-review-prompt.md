# Comprehensive review brief — current state 2026-10-12

Supersedes `claude_comprehensive_review_prompt_2026-10-02.md`, which predates the
2026-10-12 CTJ revision (commits `22b7765`..`2b142f9`). Used to drive four parallel
read-only review tracks before the advisor meeting of 2026-10-14.

## Current state being reviewed

- CTJ (`paper/ctj_manuscript/`): `CTJ_Manuscript.tex` is the master for both layouts;
  `CTJ_Reading.tex` only selects the two-column branch; captions/cells live in
  `CTJ_Exhibits.tex`. Six main exhibits (2 tables + 4 figures). Main text now shows
  the six detailed designs plus Simple Random Sampling (SRS); the two related variants
  (plain saturation, ID 3; balanced halves, ID 7) are in the supplement and chapter only.
  "Frozen" partitions are now called "fixed" everywhere.
- 2026-10-12 changes: ethics and telephone placeholders removed; generative-AI
  acknowledgement removed; Walther ORCID added; author contributions added; data
  availability now "simulation data and code on request; NC death-certificate data
  restricted"; Figures 2–4 rebuilt as vector PDFs (Fig 2 adds SRS; Fig 3 shows
  allocation MSE distributions; annual heatmaps omit the all-1.00 SRS row).
- Supplement `Supplementary_Information.tex`; chapter `paper/dissertation_chapter/Dissertation_Chapter.qmd`
  (canonical source; synced to bios-dissertation by a post-commit hook).
- Figure scripts: `code/20_manuscript_figure_revision.R`, `code/21_downside_by_design.R` (newest).

## Guardrails (all tracks)

- **Read-only.** Report findings; do not edit, commit, push or re-run simulations.
- Restricted data: never open `application/data/*.csv` or `application/data/derived/`.
  Aggregate exhibits under `application/results/real_sud_rev_20261002/exhibits/` are fine.
- Project 1 (`projects/SpillSpatialDepSim/`) is final/read-only; do not reopen its audit.
  CTJ must cite it as the BMC paper (DOI 10.1186/s12874-026-02968-0), never "Chapter 2"/"Project 1".
- Manuscript prose: no DIM mention at all; descriptive design names, not "Design k".
- Grid allocation-risk pilot numbers (`results/allocation_risk/`) must not stand in for
  full-study numbers; label them as pilot/selected-surface.
- Claims must be regime-specific (control-only vs both-arms spillover); coverage
  (~92.6–93.6%) is below nominal and must not be called nominal; tail membership is uncertain.
- Application primary outcome: β = 0, τ = 1, σ = 1, queen primary, fit with intercept/Z/true
  Spill, no X; β = 1 X-adjusted is a sensitivity only.

## Reference findings

`docs/plans/nc_sud_application_findings_2026-10-02.md`,
`docs/plans/manuscript_completion_findings_2026-10-02.md`,
`docs/plans/ctj_exhibit_revision_findings_2026-10-02.md`,
`paper/ctj_manuscript/submission_checklist.md`, project `AGENTS.md`.

## Output format

Findings tagged **ERROR** (wrong number/claim, broken build/ref, rule violation),
**WARN** (likely confusing, inconsistent or non-compliant), **NOTE** (polish), each with
`file:line`, evidence (quote ≤ 1 line or computed value vs source value) and a proposed fix.
State explicitly what was checked and found clean.

## Tracks

1. Numbers & internal consistency (CTJ, exhibits, supplement vs aggregate CSVs).
2. Live CTJ journal compliance + first-reader clarity + rendered-PDF visual check.
3. Dissertation chapter ↔ revised CTJ/SI consistency.
4. Code audit of grid pipeline, application pipeline and figure scripts; cheap test suites.
