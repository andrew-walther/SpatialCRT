# Approved CTJ exhibit revision and chapter review — 2026-10-02

The author approved restoring the original title, revising the six main exhibits
for readable two-column publication, and reviewing the chapter for corresponding
clarity changes. During implementation the author specified the order **Andrew
Walther, Ashkan Habib, Ross Joseph Simpson, Jr., Feng-Chang Lin**; Lin remains
corresponding. No new scientific design choice or simulation was introduced. The author then
requested fuller CTJ content within the 3,500-word cap. The master now restores
substantive chapter material on design mechanisms, conditional effect-size
rankings, the spillover-omission sensitivity, allocation-risk interpretation,
NC adaptations and protocol implications.

Title: *Sampling Design for Spatial Cluster Randomized Trials Under Heterogeneous
Incidence*. Both CTJ versions use this title and author block from one master.
The supplement and cover letter use the restored title as well.

## Exhibit selection and placement

| Exhibit | Content and source | Role |
|---|---|---|
| Table 1 | All eight strategies plus SRS; canonical grid rules and approved NC adaptations | Explain the compared allocations concisely |
| Table 2 | Nine-design queen, τ=1 MSE/coverage benchmark; `results/srs_benchmark/manuscript_regime_means.csv` | Include omitted plain saturation, Balanced Halves and SRS explicitly |
| Figure 1 | Observed 2018 cluster rates and frozen saturation regions | Show application inputs and geography; all-year maps remain in SI |
| Figure 2 | Six-design queen, τ=1 MSE by all five incidence configurations and two regimes | Replace eight dense overlaid-bar panels with two compact colored point panels |
| Figure 3 | Squared bias versus empirical variance for Poisson ρ_X=0.2, separately by regime | Explain the accuracy mechanism without implying this configuration represents every setting |
| Figure 4 | Nine-design annual NC mean-MSE/SRS ratios, separately by regime | Show year and design differences in 72 clearly labeled cells |

Figure numbering follows first citation: the application input map is introduced
in Methods before the grid results. This changes the provisional ordering of
the approved list, without changing its content. Tables retain separate counters.
Both submission and reading versions have exactly two tables and four figures.

All four graphics are vector PDFs exported at 174 mm, with 9-point base labels
and approximately 7.4–7.7-point printed component/cell labels. These multi-panel
figures span both reading columns; forcing them into a single column would make
the design names and ratio cells too small. The actual reading display is about
174 mm and the review display about 160 mm; both are visually checked. Color is
reinforced by shapes in the grid comparison and printed numerical values in the
annual heatmap. The map uses authorized cluster-level aggregates, not county
source records or a synthetic surface.

The double-spaced submission draft retains separate exhibit sheets after
references, following the [Clinical Trials instructions](https://journals.sagepub.com/author-instructions/ctj)
checked earlier this date. The two-column reading preview embeds the shared
exhibits in the body before references. It remains an author preview, not a
publisher proof. The supplied Sage class is unchanged.

## Chapter and supplement review

The chapter body now uses the compact annual ratio heatmap. Its annual absolute
MSE table, substantive tail findings, tail uncertainty, coverage limitation,
treatment counts/population shares and sensitivity interpretation remain in the
body. The detailed absolute-MSE and allocation-tail figures are retained in
Appendix A7. The derived supplement mirrors the heatmap and moves the two
detailed figures to its application-detail section. All 22 chapter/SI table
bodies remain equivalent; no values or analytical recommendation changed.

The full parameter-specific grid graphics, all-year incidence maps, allocation
examples and detailed diagnostics remain appropriate supporting material in the
longer chapter/SI. They were not given the CTJ six-exhibit limit. Existing figure
assets remain available; previously rendered draft PDFs are recoverable from Git
before this checkpoint. No core simulation outputs, cached weights, computational
manifests, restricted county inputs or Project 1 files changed.

## Export and verification teach-back

### `compact_grid_means()` in `code/20_manuscript_figure_revision.R`

**Purpose:** produce compact summaries of the completed grid experiment without
mixing incidence configurations or spillover regimes. **Logic:** select queen and
τ=1, require valid oracle fits, calculate squared bias and empirical error variance,
check their sum against MSE, then average by design/configuration/regime. **Choice:**
variance is `(n-1)/n * SD^2`, not uncorrected sample variance: MSE uses denominator
n while SD uses n−1. With n=250 the correction is 249/250. This exact displayed
decomposition does not change the simulation estimates. **Inputs/outputs:** the
completed 14,400-row oracle data frame produces 90 means (nine designs, five
configurations, two regimes), each from 16 settings. **Gotchas:** failed fits or
changed denominators must not be silently accepted; figures are descriptive means,
not pooled MC confidence intervals. **R notes:** base `aggregate()` retains explicit
grouping columns; no new statistical dependency is used.

### `annual_srs_ratios()`

**Purpose:** make the annual NC heatmap compare like with like. **Logic:** key each
row by year/regime, require one SRS denominator for each key, divide each design's
mean MSE by that matched SRS mean, and check against the approved aggregate ratio.
**Choice:** ratio of descriptive means matches the annual tables; it is not a mean
of setting-specific ratios. **Inputs/outputs:** a 72-row primary annual table gains
a `Ratio` column. **Gotchas:** missing/duplicate SRS keys, nonfinite MSE or incorrect
precomputed ratios fail; shared annual sources are not independent evidence.
**R notes:** named `match()` avoids positional or reordered-year mistakes.

### `export_compact_exhibits()` and entry point

**Purpose:** render the four selected displays with traceable numerical records.
**Logic:** read completed results, call the two numerical helpers, select the
labeled six-design grid subset, draw separate-regime comparisons/decompositions,
draw the nine-design annual heatmap, and join observed 2018 rates/region IDs to
58 named geometries. Save three plot-data CSVs, four PDFs and a source note under
`results/manuscript_exhibit_revision_20261002/`; copy the four PDFs to CTJ and the annual heatmap to the chapter
figure directory. **Choice:** common log MSE axes retain the large control-only
differences; a diverging log-ratio heatmap centered at SRS=1 makes lower/higher
MSE visible with exact printed cells. Fixed shapes supplement colors. Map
simplification is for display only. **Inputs/outputs:** completed grid RDS,
authorized annual/cluster aggregates, frozen partition and existing setup geometry;
four vector exhibits and their provenance records. **Gotchas:** requires local
ignored `setup.rds` and completed results; run from the project root. No restricted
county files are read. The export overwrites only its own compact figures/records.
**R notes:** existing `ggplot2` handles plotted facets, `sf` handles geography,
and base `grid` places the two map panels without adding a composition package.

### Tests and source checker

**Purpose:** reject scientifically wrong but plausible displays and prevent layout
drift. **Logic:** fixtures verify exact error decomposition, separated configurations/
regimes and permuted year-matched denominators; deliberately corrupt values or
duplicate denominators must fail. A real-data check independently recomputes all
90 exported means from their 16 original scenarios. The Python checker validates
36 main grid cells in all three sources, 20 annual cells in chapter/SI, 72 plotted
ratios, 360 full annual cells, 72 budget cells, 32 sensitivity cells and all 22
chapter/SI table bodies. It also checks title, author order, assets/citations,
six exhibit definitions/placements, shared reading source and separate legends.
**Choice:** exact numerical fixtures and matched keys test intent, not just output
shape. **Inputs/outputs:** aggregate sources and document text; explicit PASS or
failure. **Gotchas:** source checks do not replace review of scientific prose or
visual rendering. **R notes:** tests use base R; Python uses its standard library.

### Manuscript source changes

**Purpose:** integrate the selected exhibits and make the two versions identical
in scientific content. **Logic:** `CTJ_Exhibits.tex` owns cells/captions/assets;
the master calls each macro once in the selected layout, embedded or at the end.
The reading wrapper contains no manuscript prose. The restored title and new
author order live in shared front matter; legend sheet, cover letter, checklist
and chapter-derived declarations match. The chapter is edited only in its
canonical QMD; the existing post-commit hook renders/syncs Chapter 3. **Choice:**
move detailed graphics rather than delete substantive evidence or impose journal
constraints on the thesis. **Inputs/outputs:** existing scientific text and revised
figures; paired CTJ PDFs, chapter and SI. **Gotchas:** rebuild both CTJ branches;
manual duplicate prose or altered figure numbering can drift. **R notes:** no
outcome simulation or estimator code was edited.

## Review record

- The first export exposed a wrong partition CSV path and floating-point axis
  clipping of two stacked bars. Both were corrected; the final export has no
  warnings or dropped rows. No scientific results were changed to fix presentation.
- Intent-based R checks and manuscript Python checks pass. Caption/legend content,
  author order, title, references, assets and six-exhibit counts are checked.
- Final builds: CTJ review 21 pp, reading 8 pp, chapter 66 pp, SI 44 pp.
  Main prose/heading count **3,215**; with the 38 math-expression units,
  **3,253**. With end declarations/headings, prose totals **3,381** (or **3,419**
  including math units), still within 3,500. The title page and cover letter use
  the conservative 3,381 count. Abstract **283**; references/exhibits excluded.
- Reading has no final LaTeX warning; neither CTJ branch nor SI reports an overfull
  box or unresolved reference/citation. Review retains existing Sage geometry
  override and intentional separate-sheet whitespace notices. SI retains existing
  underfull-spacing/float-only-page notices; these are reviewed visually.
- All 8 reading, 21 submission, 66 chapter and 44 SI pages inspected in rendered
  contact sheets; full-size checks of updated front matter, tables, figure labels/
  captions and the chapter/SI annual heatmap found no clipping or lost exhibits.
  The final author-order changes were re-rendered and inspected separately.
- Existing grid/benchmark/main-application/tail source hashes match their earlier
  provenance record exactly. Current source/figure hashes are in
  `paper/manuscript_provenance_20261002.json`. The approved hook synced
  Chapter 3 from `3fc8385` to bios-dissertation `fa49f47`: 72 pp, 31 verified
  citekeys, one new figure copied, no figure removed, no push. The sync commit
  touches only the authorized Chapter 3 paths.

Remaining work is author review and submission-only metadata, including the
explicitly approved IRB/data-use placeholder. No push or journal submission.

Checkpoint: `3fc8385`. Current aggregate/source checks pass (87 hashes), and the
added effect-size and sensitivity values match completed outputs. The author can
now review the fuller article and updated exhibits; the Claude prompt names the
new source/export checks.
