# NC SUD Application

This directory adapts the IncidenceDesign simulation framework from a regular 10x10 grid to the North Carolina Community College service-area application. The application uses 58 irregular spatial clusters, Queen-contiguity neighbors, and population-balanced design adaptations to evaluate treatment assignment strategies for a future SUD intervention study.

## Directory Structure

```text
application/
  code/      R scripts for mapping, real-data aggregation, weights, synthetic incidence, designs, and runs
  data/      NC Community College mapping CSVs and nc_cluster_weights.rds (tracked) + restricted
             SUD source data and derived/ outputs (gitignored; see "Real SUD incidence data" below)
  tests/     test_application_data.R
  notes/     implementation notes and full-run instructions
  report/    Quarto report source plus rendered HTML/PDF outputs
  results/   smoke, pilot, full, and cached population outputs
```

## Study Region

Counties are assigned to one of 58 clusters based on primary North Carolina Community College service coverage. The map below shows the county-level assignment used to form the irregular spatial clusters for randomization and outcome simulation.

![NC Community College service-area clusters](report/figures/community_college_service_area_clusters.png)

## Current State

The synthetic-incidence application pipeline is implemented and has been run end to end.

- Synthetic baseline SUD incidence is generated with a Poisson SAR process over the 58 clusters.
- All 8 treatment assignment designs have irregular-map adaptations.
- Smoke, pilot, and full evaluation profiles are available through `code/run_application_profiles.R`.
- The full synthetic run completed 640 scenarios and 160,000 converged MLE fits.
- The current report is available in `report/IncidenceDesign_Application_Report.html` and `report/IncidenceDesign_Application_Report.pdf`.

## Design rules changed 2026-09-24 (not yet re-run)

`code/application_designs.R` now follows the simulation revision's design rules
(`../docs/plans/simulation-revision-spec.md` §3):
- Incidence ties are broken at random for every draw.
- High Incidence Focus treats exactly round(N/2) = 29.
- Balanced Quartiles treats exactly floor(N/2) = 29.
- Only the graph checkerboard is deterministic.

The application has **not** been re-run with these rules (plan decision). The results below
and in `results/` predate them, and `code/14_manuscript_supplement_figures.R` labels its
application table STALE. The application still uses `lagsarlm`; porting the lean engine is
future work.

## Real SUD incidence data (aggregated 2026-09-25)

Real county-level SUD data now replace the synthetic placeholder **as input**. The design
comparison has not yet been re-run on them; that is the next plan.

**Sources** (copied from OneDrive `UNC Dissertation (Lin)/Project 2 - Prior Incidence
Weighting/Application - Sudden Death/SUD Data - Ashkan/data/`; **gitignored**, because they are
restricted death-certificate data and this repo is public):

| File | Role | Content |
|---|---|---|
| `data/sudden_county_year.csv` | **Numerator** | Habib's `Temporal Trends Data.R` output: `CORES` (3-digit county FIPS), `DOD_YR`, `num_obs` = sudden unexpected out-of-hospital deaths, ages 18–64, after the case-definition filters in Habib (2026). 21,147 deaths. |
| `data/final_county_sudden.csv` | **Denominator** + covariates | `pop_18_64` = SEER year-specific population aged 18–64 (Σ 2018–2021 = 25,594,321 person-years); 36 static county covariates plus year-varying `cvd_rate`. Its own `num_obs` (23,523) comes from an undocumented step and is used only for reconciliation. |

The two files join on `county_fips = 37000 + CORES`. With these choices the aggregation
reproduces Habib (2026) exactly: 75.1 / 77.9 / 86.7 / 90.6 per 100,000 by year, 82.6
overall, and county range Orange 41.2 to Swain 215.6 (the paper says 216.0).

**Code** (`code/`; each file has a header saying what it does and where it fits):

| File | Role |
|---|---|
| `run_sud_aggregation.R` | Entry script: runs the three steps below, checks that clusters add up to counties, writes `data/derived/` |
| `sud_load_data.R` | Step 1: `load_sud_county_data()` reads and joins the two sources; stops unless it's a complete 100 × 4 panel |
| `sud_incidence.R` | Step 2: `compute_incidence()` (rate formulas) and `aggregate_sud_to_clusters()` (county → college) |
| `sud_reconcile.R` | Step 3: `reconcile_sud_counts()` writes the QC report against Habib (2026) and the 23,523 counts |
| `build_cluster_weights.R` | Queen and rook contiguity weights for the 58 clusters → `data/nc_cluster_weights.rds` |
| `cc_mapping_data.R` | County → college mapping (strict 1-to-1, 100 counties → 58 colleges) |

**Regenerate** (after placing both sources in `data/`):

```bash
Rscript projects/IncidenceDesign/application/code/run_sud_aggregation.R
Rscript projects/IncidenceDesign/application/code/build_cluster_weights.R
Rscript projects/IncidenceDesign/application/tests/test_application_data.R
```

Outputs in `data/derived/` (gitignored):

| File | Content |
|---|---|
| `sud_cluster_incidence.csv` | Long: 58 colleges × period ∈ {2018, 2019, 2020, 2021, 2018-2021}; deaths, person-years, rates |
| `sud_county_incidence.csv` | The same for the 100 counties |
| `sud_reconciliation_report.txt` | Statewide totals vs the paper, the numerator checks, and the 21,147 vs 23,523 comparison |
| `figures/observed_sud_incidence_<period>.png` | Cluster maps (`plot_real_sud_incidence_maps()` in `run_application_profiles.R`) |

**Rate definitions** (d = deaths, P = population aged 18–64, T = 4 years):
- yearly: d_t / P_t × 100,000
- **average** (`rate_per_100k` in the 2018-2021 rows): Σd / ΣP × 100,000, per 100,000
  person-years
- **total** (`rate_cumulative_per_100k`): Σd / (ΣP / T) × 100,000 = T × average, deaths per
  100,000 over the whole period
- `mean_of_yearly_rates`: diagnostic only. It differs from the average by ≤ 0.75 per 100k.

**Key facts:**
- Pooled cluster rates are 42.0–168.9 per 100k (median 110.2).
- The smallest cluster-year count is 4 deaths.
- Year-to-year cluster rank Spearman is 0.72–0.82.
- Switching to the 23,523 counts barely changes cluster ranks (Spearman 0.98). The code
  switch is `deaths_col = "deaths_final_csv"`.

`load_real_sud_data()` / `integrate_real_sud_data()` in `run_application_profiles.R` now wrap
this code. They used to guess columns by regex, which picked `county_name` as the count
column, and silently filled missing counts with 0 and missing populations with 1.

### Contiguity weights (`data/nc_cluster_weights.rds`, tracked)

Public geography only, so it's committed.
- **Contents:** `queen` and `rook` entries, each holding `nb` (spdep neighbor list),
  `adjacency` (binary, symmetric) and `W` (row-standardized, rows sum to 1), plus a
  `clusters` table. Rows and columns are named by `Primary_College`, so join by name, not
  position.
- **Degree:** queen 1–8 neighbors per cluster (mean 4.93), rook 1–8 (mean 4.79). Tri-County
  CC, the state's southwest corner, has one neighbor.
- **Queen-only pairs (corner touches):** Montgomery–South Piedmont, Richmond–Stanly,
  Johnston–Vance-Granville, Nash–Wake Tech.
- **Boundaries:** built from the full legal county boundaries (`cb = FALSE`), so counties
  that meet across water are neighbors. Author decision, 2026-09-25: an education
  intervention can plausibly benefit nearby counties across a river or sound. Compared with
  the shoreline-clipped cartographic boundaries (`cb = TRUE`, used by the earlier synthetic
  runs), this adds 3 pairs: College of The Albemarle–Martin CC (Chowan River),
  Carteret–Pamlico (Neuse River) and Carteret–Beaufort County CC (Pamlico Sound).

### Source-data notes and open questions for Ashkan Habib

- **Two scripts in `SUD Data - Ashkan/scripts/`:**
  - `Temporal Trends Data.R` builds `sudden_county_year.csv`.
  - `Data Age and County Residence Restriction.R` is a Lenoir-only extract.
  - Neither builds `final_county_sudden.csv`.
- **Two issues in `Temporal Trends Data.R`, neither of which explains the 23,523:**
  - It defines a heart-failure (I50) pattern but never applies it, although the paper
    excludes heart failure.
  - Its multi-line free-text regex can't match "KIDNEY FAILURE" or "LIVER FAILURE".
- **Questions (confirmation items, not blockers):**
  1. Which script built `final_county_sudden.csv`'s `num_obs`, and which filters differ?
  2. Which count should publications cite (21,147 / 82.6 or 23,523 / 91.9)?
  3. Was heart failure (I50) meant to be excluded?
  4. Should the free-text filter match KIDNEY / LIVER FAILURE?
  5. Is `pop_18_64` the SEER mid-year estimate, and the right denominator for both counts?
  6. Does the data-use agreement allow publishing county- or cluster-level counts and maps?

## Current Results

Using synthetic incidence, the best-performing designs by mean MSE were Balanced Quartiles and Balanced Halves, followed closely by Saturation Regions, 2x2 Blocking, and Incidence-Guided Saturation Regions. Block Stratified Sampling, Isolation Buffer, and High Incidence Focus performed substantially worse.

The current working recommendation is to carry Balanced Quartiles, Balanced Halves, Saturation Regions, and Incidence-Guided Saturation Regions forward as the main candidate designs. This recommendation should be finalized only after repeating the analysis with true SUD incidence data.

## Outstanding Work

- **Design comparison on the real surface (next plan).** Starting idea: the four yearly
  surfaces (2018–2021) replace the synthetic surfaces, with design draws × simulated outcomes
  within each year. Port the runner onto the revised engine (`fit_sar_lag()`, key-based
  seeds, 6 manuscript designs, queen/rook). Decide which surface feeds X and the τ/ρ/γ grid.
- Confirm the case definition with Ashkan Habib (questions in `sud_reconciliation_report.txt`).
- Then update the report, the recommendation, and the manuscripts' Application text.
