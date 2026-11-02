# Project cleanup review (2026-11-02)

A read-only review of the IncidenceDesign folder for stale code, duplicated or unused
files and out-of-date status documents, followed by the cleanup the author approved.
Commits `c7d62fb` to `230ea6b` on `main` (the deck speaker-script edits in `fd10271` came
just before). Details per commit: `git log --stat c7d62fb^..230ea6b -- .`.

## Findings

**Hazards**
- `code/complete_after_mle.R` overwrote `CLAUDE.md` and `README.md` with April 2026
  content if run. The 2026-09-24 revision plan had already slated it for deletion.
- `.claude/agents/phase-d-*.md` still loaded as agents but pointed at the retired
  `paper/dissertation_chapter/Dissertation_Chapter.qmd` and the superseded
  `20260924_025509` run.
- The `.gitignore` rule `paper/**/*.tex` silently ignored any new `.tex` file under
  `paper/` (the live sources were force-tracked).
- README and AGENTS named run `20260924_025509` (12,800 rows) as primary. The production
  run is `20260927_000502` with 14,400 rows (nine designs including Simple Random
  Sampling), verified in `results/sim_data/full_run_verification.txt`; Designs 1–8 are
  identical between the two runs (`results/srs_reproduction_check.txt`).

**Unused or duplicated material**
- Synthetic-incidence application outputs and the proposal report, superseded by the
  observed-incidence run of 2026-10-02.
- Figures in the CTJ, chapter and deck figure folders that no source includes.
- Byte-identical tracked copies, the largest set being `results/archive/pre_revision_20260924`
  against `pre_tau_sweep_20260408`.
- Removing tracked files does not shrink the Git history.
- Large gitignored local caches (methods-extension checkpoints and fits, allocation-risk
  caches, the tigris cache) include clearly superseded runs: the first A3 participation
  checkpoints (replaced by `a3v2`), the two-surface allocation-risk `pilot/` and
  `precision_R400/`, smoke folders. These were listed for the author, not deleted, because
  untracked files cannot be recovered.

**Talk deck** (`presentation/sudden_talk_2026-10-29.qmd`)
- The spillover-regime schematic and the case-identification flowchart were inline
  SVG/HTML, and the flowchart carried hand-typed counts.
- Many PNG files in `presentation/figures/` were unused leftovers from the first draft
  (ten PNGs are now in `presentation/figures/archive/`). Every used PNG is built by a
  script.

**Status documents**
- README, AGENTS, ROADMAP, `docs/TODO.md`, `application/README.md`, `application/AGENTS.md`
  and `longleaf_setup/README.md` carried stale run tags, missing paths, superseded
  synthetic-run sections and stacked dated banners.

## Done

- Deleted `code/complete_after_mle.R`, `code/_generate_design_options.R` and
  `.claude/agents/phase-d-*.md`.
- Moved the synthetic runs and proposal report to `application/archive/synthetic_2026-07/`
  (with a README); `code/14` reads its STALE application table from there.
- Moved unused figures to `archive/` inside `paper/ctj_manuscript/figures/`,
  `paper/thesis_chapter/figures/` and `presentation/figures/`, each with a README.
- Moved `run_overnight.sh` into `longleaf_setup/`; removed the empty `results/longleaf/`.
- Replaced ten byte-identical files in `results/archive/pre_revision_20260924` with
  relative symlinks to `pre_tau_sweep_20260408`.
- Narrowed the `.gitignore` rule to explicit `.tex` by-products.
- Codes 12, 14 and 15 print "2×2 Blocking" in figures (the CSVs the verifier reads keep
  "2x2 Blocking"); figures regenerated, SI and Chapter 3 PDFs rebuilt.
- `code/24` reads `summaries_v2/a4_summary.csv` (byte-identical to the old copy).
- New `paper/tools/copy_figures.sh` copies the 20 figures that live under `results/` into
  both manuscript figure folders, a step that was previously undocumented.
- New `presentation/tools/make_figures.R` draws the two deck diagrams as PNG files. The
  two upstream case counts are recorded once in `presentation/data/case_identification_counts.csv`.
- Status documents refreshed (this review's companion edit). New standing rule in AGENTS:
  figures are files built by code.

## Left open

- Confirm the two upstream case-identification counts with Habib (`docs/TODO.md`).
- Gitignored superseded caches remain on disk until the author decides to delete them.
- Byte-identical figure copies between the CTJ and chapter folders remain; both folders
  are referenced by their own sources.
- Manuscript numbers are still hand-entered and checked by
  `paper/tools/verify_manuscripts.py`; generating them from results is a `docs/TODO.md` item.
