# AGENTS.md — application/

**Verified completion checkpoint (2026-10-02):** the observed-incidence application
and focused tail refinement are complete. Production: 1,248 reporting rows,
936 distinct sources, 7,598,000 independent outcome fits; separate tail refinement
adds 2,280,000 outcomes on the same allocations. All settings pass completeness
and the approved mean-MSE/coverage Monte Carlo gates; no warnings, aliases,
boundary fits or failures. Coverage remains below 95% and allocation tails remain
uncertain. See [completed findings](../docs/plans/nc_sud_application_findings_2026-10-02.md).
Chapter/appendix integration and derived CTJ/SI are the active remaining work.

> Loaded only when working under `application/` (via the `@AGENTS.md` stub in `CLAUDE.md`).
> Moved from the project-level [../AGENTS.md](../AGENTS.md) "Current State (2026-09-25)" section.

## Revised application pipeline

- Implementation is approved. `real_sud_setup.R` freezes observed inputs/regions/
  yearly blocks; `real_sud_simulation.R` runs nested allocations and outcomes;
  `run_real_sud.R` writes yearly reports with explicit canonical source reuse.
- Primary education outcomes have β = 0 and fit intercept/Z/Spill. The matched
  β = 1 sensitivity adds rank-scaled X to both DGP and fit. τ = 1 throughout.
- Source manifests cover computation files, not the display renderer. Do not
  alter computation and reuse its checkpoints. Boundary fits invalidate complete
  comparisons; aliases/warnings/failures remain recorded. Mean MSE/coverage
  precision accounts for duplicate frequencies; singleton repeats are not extra
  independent fits. Allocation streams are shared across equivalent supports,
  so cross-setting pooled MC SEs cannot assume independent allocation randomness.
- Outputs: `results/real_sud_rev_20261002/`; ignored fit caches; tracked cluster
  summaries and aggregate maps. `render_real_sud_companion.R` updates the offline
  HTML companion. Run commands are in README; behavioral tests are in
  `tests/test_real_sud.R`. Smoke/pilot are preliminary engineering stages.
- Findings and full function/script teach-back:
  `../docs/plans/nc_sud_implementation_findings_2026-10-02.md`.

## Real NC SUD data (as of 2026-09-25)

- **Sources** (gitignored: restricted death-certificate data, public repo; originals in
  OneDrive `.../Application - Sudden Death/SUD Data - Ashkan/`):
  - **Numerator:** `num_obs` in `application/data/final_county_sudden.csv` (Habib's
    corrected case filtering). 23,523 deaths.
  - **Denominator:** `pop_18_64` in the same file (SEER, year-specific; Σ = 25,594,321
    person-years).
  - **Reconciliation only:** `application/data/sudden_county_year.csv`, Habib's
    `Temporal Trends Data.R` output (CORES = 3-digit county FIPS, DOD_YR, num_obs), the
    counts behind `habib_temporal_2026`. 21,147 deaths. Join: `county_fips = 37000 + CORES`.
- **Decision (author, 2026-09-25; reverses the earlier choice of the 21,147 counts):** the
  numerator is `final_county_sudden$num_obs` (23,523). Reason (Habib, personal communication,
  2026-09-25): that file comes from his corrected case filtering, which no longer excludes
  heart-failure deaths as presumed non-sudden, because adjudicated sudden cardiac death
  overlaps non-negligibly with heart-failure patients. The corrected counts are ≥ the paper's
  in all 400 county-years (equal in 48). Statewide rates: 83.9/85.9/97.1/100.5 by year, 91.9
  pooled (paper: 75.1/77.9/86.7/90.6, 82.6). The 21,147 counts stay as `deaths_habib2026`,
  and the tests check that column still reproduces the paper (including Orange 41.2, Swain
  215.6 vs the paper's 216.0). Cluster ranks barely differ (pooled Spearman 0.982). All
  analysis uses the corrected counts even though they no longer match the paper; the chapter
  text is updated separately. Remaining questions for Ashkan: `application/README.md`.
- **Rates:**
  - yearly: d/P
  - average: Σd/ΣP (per 100k person-years)
  - total: Σd/(ΣP/4) = 4 × average
  - `mean_of_yearly_rates` is diagnostic only
- **Outputs:** `application/data/derived/`, which regenerates via
  `Rscript application/code/run_sud_aggregation.R`. Pooled cluster rates 45.3–194.9
  (median 125.6); min cluster-year count 5; year-to-year cluster rank Spearman 0.72–0.82.
  Maps come from `plot_real_sud_incidence_maps()` in `run_application_profiles.R`.
- **Fixed:** `load_real_sud_data()` guessed columns by regex (it picked `county_name` as the
  count column), and `integrate_real_sud_data()` filled missing values with 0/1. Both now wrap
  the new script. `cc_mapping_data.R`'s directory detection resolved to `getwd()` when sourced
  inside `suppressMessages()`; it now uses the innermost `source()` frame.
- **Contiguity weights:** `application/data/nc_cluster_weights.rds` (tracked; build with
  `build_cluster_weights.R`). Queen degree 1–8 (mean 4.93), rook 1–8 (mean 4.79); 4
  queen-only corner pairs. **Legal boundaries (`cb = FALSE`)** by author decision
  (2026-09-25; nearby counties across water can benefit from an education intervention).
  This adds Albemarle–Martin, Carteret–Pamlico and Carteret–Beaufort County CC relative to
  the shoreline-clipped `cb = TRUE` used by the earlier synthetic runs.
- **Code layout (user, 2026-09-25):** application scripts stay under `application/`, split
  one job per file. Script headers say "Author: Andrew Walther".
- **Habib script notes:** the `hf` (I50) pattern is defined but unused. That's explained:
  under the corrected filtering, heart failure is deliberately not excluded (Habib,
  2026-09-25). The multi-line free-text regex can't match "KIDNEY FAILURE" or "LIVER
  FAILURE".
