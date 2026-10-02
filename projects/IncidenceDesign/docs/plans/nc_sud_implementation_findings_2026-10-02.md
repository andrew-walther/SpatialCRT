# NC observed-incidence implementation: verification and walkthrough

Status: approved implementation, full production and same-allocation tail extension
are complete and independently verified. Full chapter/appendix and derived CTJ/SI
drafts are complete. Scientific findings: `nc_sud_application_findings_2026-10-02.md`;
document provenance, review and export/checker teach-back:
`manuscript_completion_findings_2026-10-02.md`. Historical pilot/execution entries
below remain engineering history, not current study status.

## What is implemented and why

The primary outcome is simulated continuous education performance, with
Y = (I − ρW)⁻¹[τZ + S(Z) + ε]. Observed annual SUD rates inform allocation;
they do not enter its baseline. τ = 1 and σ = 1. The explicitly hypothetical
matched baseline/adjustment sensitivity adds X = (average rank − 0.5)/58 to the
mean and fits X. This implements the author's education-outcome clarification.

There are 1,248 production reporting rows: 432 primary queen, 432 matched baseline
queen, 288 primary rook corners and 96 Design 8 summary alternatives. Content
keys reduce these to 936 computational sources. Shared yearly rows
are repeated reporting, not independent yearly evidence. Annual population
shares and incidence contrasts are recalculated even when performance is reused.

The three Design 8 summaries do sometimes change regional ordering. In frozen
region IDs (highest → lowest), the primary mean-rank rule gives 2/1/3/4 in 2018,
4/2/3/1 in 2019, 2/4/1/3 in 2020 and 2/1/4/3 in 2021. Mean raw rates give
2/4/1/3, 4/2/3/1, 2/4/1/3, 2/4/1/3. Pooled regional deaths/population gives
4/3/2/1, 4/3/1/2, 4/3/2/1, 4/3/1/2. These are input-rule diagnostics,
not performance conclusions; keep the approved primary rule regardless of which
alternative performs best.

## Function-by-function teach-back

All vectors/matrices use the cached weight matrix's ordered 58 college names.
Inputs are fixed observations, not synthetic incidence. The following explains
purpose, logic, design choices, inputs/outputs, edge cases and R/package choices.

### Allocation adaptation (`application_designs.R`)

`get_application_designs()` accepts a design ID, allocation count, cluster signal,
neighbors, coordinates, populations, county counts and optional frozen region/
block/summary objects. It returns an integer 58 × J treatment matrix. Existing
rules are retained; supplied block IDs prevent repeated k-means fitting and
supplied regional summaries enable the agreed focused sensitivity. Balanced
Halves starts at floor(size/2), randomly gives the remaining treatment to one
odd-sized half, then permutes within each half: 14 + 15 = 29. Region saturation
retains rounded local budgets, spatial blocks retain their original local rounding,
and the buffer remains a random independent set. `dplyr::ntile` defines equal
rank groups, with random tie ranking preceding it. Designs may therefore have
different treatment counts; do not interpret their comparison as equal-budget
optimization. A graph Checkerboard does not reproduce the rectangular-grid
identity WZ = 1 − Z automatically.

### Frozen setup (`real_sud_setup.R`)

- `rs_seed(...)`: wraps the existing explicit-RNG key seeder with a distinct NC
  namespace. Each allocation index or fixed-size noise batch gets its own seed.
  More replications preserve old prefixes. Numeric fields use the established
  two-decimal convention; do not introduce parameters closer than that resolution
  without revising the seed convention and manifest.
- `rs_hash_files(paths)`: SHA-256 hashes existing inputs/code and stops on missing
  files. It returns a named character vector; `digest` hashes actual contents,
  so filenames or timestamps alone cannot authorize stale reuse.
- `rs_setup(root)`: reads ignored cluster observations and tracked cached weights,
  orders by college name, validates periods/rate formulas/cluster sums and
  nb↔A↔W agreement, dissolves public legal county geometry to clusters, and checks
  its contiguity against cached A. It freezes one geographic/balance-selected
  four-region map and one spatial-block map per year, then saves `setup.rds` and
  `frozen_partitions.csv`. `sf` handles projected geometry, `spdep` checks
  contiguity and `dplyr` dissolves clusters. Geometry is not rebuilt to improve
  performance. Existing setup disagreement stops instead of overwriting it.
- `rs_manifest(root,setup,config)`: saves/checks configuration, input/source hashes,
  design objects, R/packages, RNG and BLAS. It returns the manifest invisibly.
  A mismatch refuses checkpoint reuse. The display renderer is excluded because
  display edits do not change computation; the three computation modules and
  unchanged grid helpers are included.

### Nested simulation (`real_sud_simulation.R`)

- `rs_config(profile)`: returns the approved grid and replication tiers. Smoke
  and pilot shrink replication for engineering verification; production uses
  the approved queen grid and focused sensitivities. These are distinct profiles,
  not progressively pooled scientific evidence.
- `rs_units(config)`: forms the explicit annual/model/neighbor/summary/regime/
  parameter/design reporting grid. It returns a data frame, not outcomes.
  Local `grid()` only constructs its Cartesian product with base `expand.grid`.
- `rs_region_values(a,regions,method)`: returns named means of ranks, means of raw
  rates, or regional deaths/population × 100k, preserving ties. `tapply` aggregates
  by frozen membership; numeric gaps do not matter after the final region ranking.
- `rs_support_key(u,setup)`: hashes every input that changes allocation support.
  Fixed regions/blocks and rank/tie groups enter where used. This permits proven
  equivalence across years or summary rules, never equivalence based on observed
  performance. Extra distinctions in keys cost computing time but do not fabricate
  evidence; omitted distinctions would be a correctness bug.
- `rs_block_key(u,setup,config)`: adds model, X where fitted, actual W, ρ, γ,
  regime, τ and σ to the support key. Identical keys share a simulation source;
  distinct design keys have independent outcome streams.
- `rs_draws(u,setup,J)`: generates 58 × J integer assignments with one keyed seed
  per draw and frozen objects. Extending J appends draws. Equivalent supports
  intentionally share allocation streams across parameter/model settings;
  pooled uncertainty must account for that dependence.
- `rs_model(u,z,setup,config)`: constructs true S(Z), the actual fitted matrix,
  A = (I − ρW)⁻¹ and its mean. It tests treatment identification using residual Z
  after all nuisance columns, and records rank/conditioning. Base `lm.fit`, `qr`
  and `solve` keep those diagnostics tied to the actual matrix. A nuisance alias
  and an unidentified τ are different conditions; finite dropped-column estimates
  do not establish identification.
- `rs_allocation_metrics(fits)`: summarizes independent repeated outcomes at one
  fixed assignment using the existing validated allocation-risk helper. Boundary/
  error messages invalidate fits even if their raw values are finite. Raw cache
  values remain preserved. Counts, bias, MSE, coverage, conditional MC SEs and
  two-half MSEs are returned; any invalid fit marks the allocation incomplete.
- `rs_distribution(a,singleton)`: frequency-weights conditional risks and computes
  allocation risk plus duplicate-aware joint MC uncertainty. For frequency fᵢ,
  conditional MC variance vᵢ and J draws, joint variance equals raw draw variance/J
  + Σfᵢ(fᵢ−1)vᵢ/[J(J−1)]. Singleton repeats reduce to conditional R uncertainty.
  Signed variance corrections are preserved, and q90/worst-10% concern estimated
  risks, not an exactly known tail distribution. Incomplete comparisons are withheld.
- `rs_block(u,setup,config,J,R,old)`: draws assignments, caches unique masks and
  simulates fixed-size keyed Gaussian noise batches. Y = mean + Aε is fit by the
  unchanged lean SAR engine. It preserves raw warnings/aliases/failures, fits,
  draw frequencies, treatment diagnostics and summary. An old matching cache
  supplies previous fits; extending R appends the exact fresh-run suffix.
- `rs_run_block(u,setup,config,checkpoint)`: saves each refinement tier, then stops
  on completion/precision success or explicit unresolved status. It chooses J
  versus R using the dominant component of whichever normalized MSE/coverage
  target has the larger miss. The prescribed tiers bound computing; exhausting
  them does not silently become a precision success. Singleton supports refine R.

### Reporting runner (`run_real_sud.R`)

- `rs_reporting_rows(u,b,setup)`: attaches an annual reporting key to canonical
  performance and recalculates its allocation population shares/incidence
  contrasts from that year's inputs. It returns summary and allocation data frames.
  Base indexing follows college order and preserves every sampled frequency.
- `rs_run(profile,workers,root)`: validates the manifest, groups canonical keys,
  executes distinct blocks, then writes performance/allocation/warning CSVs and
  an RDS bundle. `parallel::mclapply` uses explicit key streams (`mc.set.seed=FALSE`);
  failed workers stop the runner. Full fit caches stay on disk, keeping worker
  return values small. Forking is for the current Mac/Linux environment; Windows
  needs one worker or a separately approved backend. Limit BLAS threads and match
  workers to allocated cores on HPC. Reused Source_ID rows cannot multiply the
  number of independent simulated distributions.

### Companion (`render_real_sud_companion.R`)

- `rs_json(d)`: encodes atomic reporting columns into an embedded JSON array;
  its local scalar converter escapes strings and maps missing/non-finite values
  to null. No new R dependency is introduced. Inputs must remain atomic columns.
- `rs_companion(root,profile)`: reads saved setup and optional performance, draws
  cluster maps with `sf`/`ggplot2`, and writes an offline HTML page with model flow,
  interpretation, methods and interactive native-JavaScript filters/charts/table.
  Display simplification affects maps only, never W. Smoke/pilot get explicit
  preliminary labels. Only authorized cluster aggregates are displayed; no remote
  scripts or requests are needed. Empty filter combinations show no rows because
  the approved sensitivity grid is not a full Cartesian product.

## Verification evidence

- Existing application data suite passed: hand rate fixtures, primary corrected
  numerator, published-count reconciliation and named legal-boundary weights.
- New behavioral suite passed: fixed observed panel, 29-treatment rank/SRS budgets,
  half allocation, buffer independence, regional rule, both outcome specifications,
  both spillover definitions, exact-alias flagging, lean/`lagsarlm` agreement for
  primary and matched sensitivity, exact J/R extensions, annual reporting reuse,
  singleton independent-fit counts, hand-computed duplicate covariance, finite
  boundary invalidation and manifest mismatch refusal.
- Current smoke: 116 reporting rows, no incomplete settings and no warning/failure
  fit rows. Its 115 unresolved precision rows are expected at J3/R8; never cite its
  performance as study evidence. An initial earlier-code smoke is preserved under
  `smoke_initial/`; source manifests distinguish it from the current smoke.
- Larger pilot: 928 reporting rows from 688 computational sources, no incomplete
  rows or warning/failure records; 524 precision misses at deliberately small
  J20/R100 (singletons R200). Production uses the approved refinement tiers.
  A roxygen-only addition followed the pilot; no scientific behavior changed.
- Independent full pilot cache audit and `verify_real_sud.R pilot` reconciled all
  1,116,800 physical outcome fits, every source/draw frequency and annual
  diagnostic. No mismatches. Production's computational manifest matches current
  sources; the pilot's documented comment-only hash difference prevents stale
  reuse under strict manifest checks.
- `tests/verify_real_sud.R` independently reconstructs each conditional summary,
  joint uncertainty and annual diagnostic from every source cache, checks draw/
  fit frequencies and CSV/RDS agreement, and writes `verification.txt`. It checks
  the production manifest and fails on unresolved production precision. It does
  not overwrite simulation outputs or silently drop failed settings.
- Independent read-only implementation reviewer verified model/prefix/covariance/
  budget behavior and requested boundary validity, input consistency and
  coverage-driven refinement safeguards; all three were incorporated.
- Companion maps visually reviewed as local images. In-app browser rejects local
  file URLs; interactive rendering has not yet been visually verified in-browser.
  Continue static/syntax checks and preserve that limitation explicitly.
- HTML static checks parsed the embedded result data, confirmed local image
  assets and passed `node --check` for embedded JavaScript.

### Focused tail refinement (`refine_real_sud_tail.R`)

`rs_cross_tail(b)` expands sampled frequencies, selects the worst 10% using one
outcome half and evaluates those selected allocations using the other half. It
returns both directions and the number of tail draws, withholding incomplete
metrics. Independent-half evaluation helps diagnose noisy selection; it is not
an exact estimate of the population's upper tail. A hand fixture with opposite
half rankings verifies the returned evaluation uses the other half rather than
the selected values themselves.

`rs_refine_tail(workers,root)` requires complete/precision-valid main outputs,
checks the parent computation manifest and its own source manifest, and extends
R to at least 400 at fixed main-study allocations for primary queen ρ = 0/0.5,
γ = 0.8, both regimes, all four years and all nine designs. It returns/saves
separate reporting/allocation/warning outputs under `tail_confirmation/`, verifies
all original Z and fit prefixes, and never overwrites main results. Cached
singletons with R ≥ 1000 need no new fits. This is outcome refinement, not fresh
independent allocation confirmation; conclusions will remain qualified where
sampling or conditional-noise uncertainty persists. Base indexing/order respects
duplicate frequencies; `mclapply` uses the established keyed streams.

## Remaining work

Complete pilot verification, production and focused tail diagnostics/refinement;
check every reporting row, independent-fit count, manifest and warning status;
update this note and companion with verified findings. Then integrate the full
chapter/appendix, render/sync/review, derive the full standalone CTJ/SI, check live
journal rules/verified Project 1 citation, and prepare the comprehensive Claude
review prompt. Never push without asking. Project 1 remains accepted/final/read-only.

## Completed reporting deliverables: teach-back

`real_sud_summary.R` turns independently verified results into manuscript-ready
CSV tables and standard ggplot2 PNG/PDF figures. It validates manifests and gates,
joins SRS within exactly matched settings, joins sensitivity references at the
same year/parameters, and exports annual descriptive averages. Inputs are main
and refined results bundles; outputs are tables, figures, aggregate inputs and
an exhibit manifest. It deliberately omits pooled MC intervals because shared
allocations make independence inappropriate. Singleton risks and main versus
refined outcome replication must remain distinct.

`test_real_sud_summary.R` uses actual and incomplete-reference fixtures to verify
SRS differences/SEs, zero self-comparison error and matched sensitivity contrasts.
`test_companion.mjs` executes the actual saved JavaScript in a small DOM fixture,
changes evidence/year/model/neighbor/summary/regime/parameter controls and metrics,
and checks nonconstant SVG bars, absent unplanned cells and all seven image assets.
Inputs are the HTML and its referenced files; failures throw assertions. Node's
built-in fs/vm/assert modules add no dependency. This is a behavior check, not a
browser visual review. The in-app browser blocks local file URLs, so no actual
browser rendering is claimed. Local scientific/map PNGs were visually reviewed.

Full scientific findings and completed verification counts are in
`nc_sud_application_findings_2026-10-02.md`. Production and tail outputs are
complete; chapter/appendix and CTJ/SI drafting are now the active work.
