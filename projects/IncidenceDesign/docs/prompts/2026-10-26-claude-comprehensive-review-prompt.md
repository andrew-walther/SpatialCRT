# Comprehensive review brief — current state 2026-10-26

Supersedes `2026-10-12-claude-comprehensive-review-prompt.md`. Since then the
methods extensions (partial participation, power, interval calibration, targeting)
were written into all three documents, Table 1 became a hybrid table with example
maps, and several author decisions changed the wording rules below. Same four
read-only tracks and output format as before.

## Current state being reviewed

- CTJ (`paper/ctj_manuscript/manuscript_ctj/`): `CTJ_Manuscript.tex` is the master for
  both layouts; `CTJ_Reading.tex` only selects the two-column branch; captions, cells
  and example maps live in `CTJ_Exhibits.tex`. Six main exhibits (2 tables +
  4 figures). Body 3,486 of 3,500 words (texcount text + headers, Introduction through
  Discussion), so any proposed addition must name an equal cut.
- Supplement `paper/ctj_manuscript/supplementary_info/Supplementary_Information.tex`;
  chapter `paper/thesis_chapter/Thesis_Chapter.qmd` (canonical; synced to
  bios-dissertation by a post-commit hook).
- New since 2026-10-12 (review these closely):
  - Methods extensions: CTJ Methods (college-subset power sentence), Results (power
    paragraph), Limitations (interval calibration); SI Section S11; chapter
    "Planning for realistic participation" and Appendix A9 (`app:nc-ext`). Sources:
    `docs/findings/2026-10-21-methods-extensions-report.md`,
    `application/results/methods_ext_2026-10/`, `code/23_*` and `code/24_*`.
  - Hybrid Table 1 (example map per design, incidence row) from
    `code/25_design_assignment_figures.R`; SI/chapter design figure
    `design_assignments_panels.pdf` from the same script (illustrative surface,
    ρ_X = 0.5, not a simulated surface).
  - Habib's affiliation (Epidemiology), Gan et al. (2019) cited for the SUD
    case-adjudication algorithm, IRB protocol 14-2036 stated; the counts are called
    corrected, with no "personal communication".
  - Reading layout uses `flushend` (References follows the declarations).

## Guardrails (all tracks)

- **Read-only.** Report findings; do not edit, commit, push or re-run simulations.
- Restricted data: never open `application/data/*.csv` or `application/data/derived/`.
  Aggregate exhibits and summaries under `application/results/` are fine.
- Project 1 (`projects/SpillSpatialDepSim/`) is final/read-only. CTJ cites it as the
  BMC paper (DOI 10.1186/s12874-026-02968-0), never "Chapter 2"/"Project 1".
- Manuscript prose: no DIM mention; the one design-name set in `AGENTS.md`
  ("Saturation Quadrants" restored 2026-10-25); no in-figure titles or notes.
- Grid allocation-risk pilot numbers (`results/allocation_risk/`) must not stand in for
  full-study numbers.
- Claims must be regime-specific; tail membership is uncertain; incidence guidance is
  framed as matching Saturation Quadrants' accuracy while targeting high-burden regions.
- Intervals: the grid reports the t + df-corrected interval throughout (coverage
  0.946–0.952 for identified designs; exact re-fit in `code/26_grid_tdf_intervals.R`,
  applied by `apply_tdf_intervals()` in 06). The application reports Wald (92.6–93.6%)
  in its main tables with t + df (94–96%) beside it; whether to switch those tables is
  an open author decision, so flag inconsistent wording, not that choice itself.
- Known placeholders (do not report as errors): IRB determination/consent wording and
  coauthor approval (highlighted), pending Dr. Simpson and the coauthors.
- Application primary outcome: β = 0, τ = 1, σ = 1, queen primary, fit with
  intercept/Z/true Spill, no X; β = 1 X-adjusted is a sensitivity only.

## Reference findings

`docs/findings/2026-10-21-methods-extensions-report.md`,
`docs/findings/2026-10-02-nc-sud-application-findings.md`,
`docs/findings/2026-10-02-manuscript-completion-findings.md`,
`paper/ctj_manuscript/submission_checklist.md`, `DECISIONS.md`, project `AGENTS.md`.
Run `python3 paper/tools/verify_manuscripts.py` and
`Rscript code/tests/test_manuscript_exhibits.R` first; both passed on 2026-10-26.

## Output format

Findings tagged **ERROR** (wrong number/claim, broken build/ref, rule violation),
**WARN** (likely confusing, inconsistent or non-compliant), **NOTE** (polish), each with
`file:line`, evidence (quote ≤ 1 line or computed value vs source value) and a proposed fix.
State explicitly what was checked and found clean.

## Tracks

1. Numbers & internal consistency (CTJ, exhibits, supplement vs aggregate CSVs,
   including the methods-extension summaries).
2. Live CTJ journal compliance + first-reader clarity + rendered-PDF visual check
   (both layouts; the hybrid Table 1 on one page).
3. Dissertation chapter ↔ CTJ/SI consistency.
4. Code audit of grid pipeline, application pipeline, methods-extension and figure
   scripts (including `code/25_design_assignment_figures.R`); cheap test suites.
