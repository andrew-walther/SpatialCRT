# Manuscript numbers from saved results: plan

Date: 2026-11-02. Status: **plan agreed in outline by the author; not started.** The author
asked to hold implementation for now. When it starts, run the stages in order, each
committed on its own, with the CTJ first.

Goal: the CTJ paper, the SI and Chapter 3 take their table bodies and their numbers in the
text from saved results through code, as the figures already do (author goal, 2026-11-02;
`docs/TODO.md`). A re-run or revision then changes every document by re-running one
generator, not by hand edits.

Inventory behind this plan: the session notes of 2026-11-02 (summarized in the next section).

## Where things stand

- **Figures** are already built by scripts and copied with `paper/tools/copy_figures.sh`.
- **Tables.** Chapter and SI share 28 table labels; the CTJ main text has two tables.
  - Generated, then pasted, with no check on the copy (7): nc-power, nc-translation,
    nc-calibration, nc-guidance, nc-car, nc-bayes, nc-forms (`code/24` writes them).
  - Checked against CSVs by `paper/tools/verify_manuscripts.py` (7): grid-srs (all three
    documents), nc-year, nc-annual, nc-budgets, nc-summary-sens, downside,
    interval-comparison.
  - Not checked (9, about 410 numbers): performance, by-config, regime, robustness,
    eight-design, rho, checkerboard-coverage, gamma, oracle-nonoracle. No CSV exists for
    these; they were most likely copied from `code/13`'s text extract.
  - Small or fixed (5): designs, excluded-designs, parameter-grid, seed-keys, nc-order.
- **Numbers in the text** (decimals and percentages outside tables): CTJ about 110, SI
  about 700, chapter about 740. Some are fixed settings (τ values, saturations, cluster
  counts). Only four CTJ values are checked (`code/tests/test_manuscript_exhibits.R`).
- All result sources needed are tracked in git (2026-09-27 simulation RDS files,
  `results/grid_tdf_2026-10/`, `results/srs_benchmark/` CSVs, methods-extension tables,
  `application/results/real_sud_rev_20261002/exhibits/`).

## Rules for every stage

- Never edit the hashed simulation code `code/01`–`05`.
- Numbers must not change. Any generated value that differs from the printed one is
  reported to the author before anything is switched (rounding, or a table that missed the
  t + df change of 2026-10-26).
- `verify_manuscripts.py` and the R tests pass after every stage.
- CTJ edits go in `CTJ_Manuscript.tex`/`CTJ_Exhibits.tex`, so both layouts follow; rebuild
  `CTJ_Manuscript.pdf` and `CTJ_Reading.pdf`.
- Fixed settings stay as literals; only numbers computed from results become macros.
- Macro names read like the quantity (for example `\gridIGSCtrlMSE`), and
  `paper/generated/numbers.csv` lists each macro with its value, source file and formula.

## Stage 0. Pipeline review and cleanup (before any generator)

Builds on the 2026-11-02 cleanup (`37bf50b`, `c7d62fb`).

1. Map the pipeline: for each manuscript exhibit and number group, the script that makes
   it, the results file it reads and that file's inputs. Write it as a table in
   `docs/findings/`. → check: every table label and figure in the three documents has a row.
2. Code review of `code/06`–`27`, `application/code/` and `paper/tools/` for simplification:
   repeated helpers (for example the relabel and `apply_tdf_intervals()` steps repeated in
   12–14), overlapping aggregations in 19–21, dead branches and unused outputs. Report
   findings first; change only what the author approves.
3. Apply approved simplifications. → check: each changed producer rewrites its outputs
   unchanged (CSV numerically identical, figures rebuilt and compared), tests pass.
4. Results and application folders: list superseded or unused outputs (tracked and
   untracked caches) for the author; archive or delete only on approval. Untracked files
   cannot be recovered once deleted.

## Stage 1. Generator and the CTJ

1. New `code/28_manuscript_numbers.R` reads only saved results and writes
   `paper/generated/tables/<label>.tex` (tabular bodies), `paper/generated/numbers.tex`
   (`\newcommand` per number) and `numbers.csv`. It reuses `apply_tdf_intervals()` from 06.
   → check: each generated table equals the current printed table cell by cell.
2. CTJ: `tab:grid-srs` body and every results number in `CTJ_Manuscript.tex` and the
   `CTJ_Exhibits.tex` captions become `\input`/macros.
3. `verify_manuscripts.py` expands `\input` and macros before its comparisons and fails
   on an undefined macro or a stale table; `test_manuscript_exhibits.R` reads the
   expanded text.
4. Submission package (the submission layout: one column, double-spaced, exhibits at the
   end, which is `CTJ_Manuscript.tex` built on its own; the two-column `CTJ_Reading` layout
   is for reading only): a script writes `CTJ_Submission.tex` from the master with the
   exhibits, tables and numbers inlined, and zips it with the figure files, `sagej.cls`,
   `SageV.bst`, the `.bib` and the rendered PDF (the author expects SAGE to want such a zip,
   2026-11-02). → check: the zip unpacked in an empty folder compiles, and its PDF text
   matches `CTJ_Manuscript.pdf`.
5. `sync_ctj_to_bios.sh` copies the flattened text so bios-dissertation stays
   self-contained.
6. Rebuild both layouts; `pdftotext` output matches the previous PDFs except build dates.
   Commit.

## Stage 2. Chapter 3

1. Chapter tables → `\input` of the same generated bodies; results numbers in the text →
   the same macros, loaded through `header-includes` (pandoc passes unknown commands
   through as raw LaTeX in PDF output).
2. `prelim_transform.py` inlines the tables and expands the macros, so the prelim copy in
   bios-dissertation keeps plain numbers. The post-commit hook syncs it.
3. Render with RStudio's Quarto 1.10; text matches the previous PDF. Commit.

## Stage 3. SI

1. Same switch as the chapter; most macros are shared with it.
2. The verifier's chapter-to-SI table check compares the expanded bodies.
3. Rebuild; text matches. Commit.

## Stage 4 (optional). One build command

A small `Makefile` for the deliverable end only: saved results → `code/28` and the figure
producers → `copy_figures.sh` → latexmk/Quarto → verifier. It reruns only the steps whose
inputs changed. It never runs the simulation or the application runs, which are
manifest-checked and expensive. The R `targets` package is the alternative; it is heavier
and a new dependency.

## Open decisions for the author

- Whether to add the Stage 4 `Makefile`.
- Which Stage 0 simplifications to apply, once the review is reported.
- Whether to merge `code/` and `application/code/` (author undecided, 2026-11-02).
  Recommendation: keep them apart and tidy the boundary instead. Reasons: they are two
  studies with different inputs (simulated grid vs restricted county data with its own
  `application/AGENTS.md` and gitignored sources); the application manifests store file
  hashes keyed by full path (`rs_hash_files()` in `application/code/real_sud_setup.R`), so
  moving production files would break the manifest checks of the stored runs; `code/01`–`05`
  must not change. Tidying instead: name the one real crossing (`code/16` is sourced by
  `real_sud_setup.R`) in both AGENTS files, give `application/code/` a clear entry-point
  order in its README, and let the Stage 0 pipeline map and the Stage 4 `Makefile` span
  both folders.
- Bibliography in the submission: settled, the `.bib` goes in the zip (see Stage 1).
- Numbered application scripts (author prefers numbered scripts, 2026-11-02). Renaming
  hashed production files (`real_sud_setup.R`, `real_sud_simulation.R`, `run_real_sud.R`,
  `application_designs.R`, the `methods_ext_*` computation files) breaks the stored-run
  manifest checks, so the choice is: (a) a numbered run order in `application/README.md`
  and numbered names only for scripts outside any manifest; or (b) rename everything and
  re-record the manifests, which needs a check that the stored runs still verify. Stage 0
  lists which files are in a manifest.

## Simulation vs application model code (checked 2026-11-02)

The author's concern with keeping the folders apart is that the two studies might fit
different models. They do not:

- **Fitter.** Both call `fit_one_lag_model()` from `code/04_estimation.R` (the lean eigen
  ML engine `fit_sar_lag()`): the grid through 05, the application through
  `code/16_allocation_risk.R` and `application/code/real_sud_simulation.R`. The
  application manifests hash `code/01`–`05` and `16`, so a change to the fitter would stop
  the application runs from loading.
- **Outcome model.** Same SAR lag form, Y = (I − ρW)⁻¹(τZ + S(Z) + βX + ε); the
  application builds it in its own `rs_model()` because its W comes from the 58-college
  map and its primary model has β = 0 with no X term (author decision 2026-10-02). Those
  differences are by design, not separate model code.
- **Interval.** Both use τ̂ ± t(0.975, n − p) × SE × √(n / (n − p)), with p the
  coefficients plus ρ (grid n = 100, p = 5; application n = 58, p = 4, or 5 with X), but
  written twice: `code/26_grid_tdf_intervals.R` and `ms_a0_a5()` in
  `application/code/methods_ext_stored.R`. Stage 0 candidate: one shared helper, with a
  test that both studies' stored values are unchanged.

## Afterwards

- AGENTS.md standing rule (numbers come from `code/28`), reason in DECISIONS.md, TODO item
  closed, ROADMAP entry.
