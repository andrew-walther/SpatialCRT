# AGENTS.md — SpatialCRT

Cross-project orientation for agent sessions (Claude Code, Codex). `CLAUDE.md`
only imports this file; edit it here. Project-level detail:
- `projects/SpillSpatialDepSim/AGENTS.md`
- `projects/IncidenceDesign/AGENTS.md` (and `application/AGENTS.md`)

## Projects

| | SpillSpatialDepSim (Project 1) | IncidenceDesign (Project 2, primary) |
|-|--------------------|-----------------|
| **Location** | `projects/SpillSpatialDepSim/` | `projects/IncidenceDesign/` |
| **Question** | Block stratified vs. simple random assignment under spillover | Which of 8 designs minimizes MSE for τ under heterogeneous incidence? |
| **Status** | **Accepted** at *BMC Medical Research Methodology* (final source: `paper/Manuscript Revisions/Revision 2c/`) | Simulation revised + re-run 2026-09-24; chapter rewritten (Phase D); CTJ manuscript + SI still carry April 2026 numbers (manuscript plan step 5); real-data application next |
| **Entry point** | `code/SpatialSim_NC_DOC.Rmd` | `code/05_run_simulation.R` |

IncidenceDesign's current results and design list are in `README.md` and
`projects/IncidenceDesign/AGENTS.md`. The April 2026 numbers are superseded
(`projects/IncidenceDesign/results/archive/pre_revision_20260924/README.md`).

## Rules

- **Writing:** no mannered prose. Say what you mean in literal terms, in
  prose, docs, and comments alike.
- **Restricted data.** This repo is public. The NC sudden-unexpected-death
  source files (`application/data/final_county_sudden.csv`,
  `sudden_county_year.csv`, `data/derived/`) are git-ignored and must never
  be committed. Results aggregated to the 58 community-college clusters may
  be tracked, and aggregate statistics in docs are acceptable (author
  decision 2026-10-01).
- **Dissertation Chapter 3 sync.**
  `projects/IncidenceDesign/paper/dissertation_chapter/Dissertation_Chapter.qmd`
  is the only copy anyone edits. A post-commit hook
  (`dissertation_chapter/tools/`, installed by `install_hooks.sh`) runs
  `sync_to_prelim.sh`, which renders it with RStudio's Quarto 1.10 and
  commits the result in `~/GithubProjects/bios-dissertation`. The prelim
  YAML header is `tools/prelim_header.yml` (sets `colorlinks: false`, so
  links print black). The sync never pushes.
- **Commits:** stage specific paths, never `git commit -a`.
- License: MIT (`LICENSE`).

## Structure

```
projects/
  SpillSpatialDepSim/   code/ data/ results/ paper/
  IncidenceDesign/      code/ results/ docs/ longleaf_setup/
    application/        data/ (restricted files ignored) code/ results/ report/
    paper/              dissertation_chapter/ ctj_manuscript/ report/ archive_manuscript/
archive/                legacy work, not maintained
```

Cross-project paths from IncidenceDesign code:
`here::here("projects", "SpillSpatialDepSim", "results")`.

## Getting started

```r
# File > Open Project > SpatialCRT.Rproj
setwd("projects/IncidenceDesign/code")
source("06_visualizations.R")
mle_results <- load_latest_results(estimation_mode = "MLE_combined")
```

Shared packages: `sf`, `spdep`, `spatialreg`, `dplyr`, `tidyr`, `ggplot2`,
`viridis`, `rmarkdown`, `knitr`, `digest`, `parallel`, `here`.

## Planned

Split into separate SpatialCRT-Spillover and SpatialCRT-Incidence repos
after the CTJ submission.
