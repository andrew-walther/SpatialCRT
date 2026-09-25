# Session prompt: real-data cluster-level incidence aggregation

> Paste this into a new session to continue the real-data application work.
> Written 2026-09-25, after the J-drive data pull and the `cc_mapping_data.R`
> full-name crosswalk (see `application/data/cc_name_crosswalk.csv`).

---

I'm working on the `IncidenceDesign` project in `SpatialCRT` (real-data application: NC Sudden Unexpected Death, SUD). This is a multi-step task — **don't start writing code. Start by asking me questions until we're aligned on a plan**, then propose the plan explicitly and wait for my approval before implementing anything.

**End goal:** apply this study's incidence-informed sample design methodology to the real state of North Carolina, clustered by the 58 NC Community College service areas, to determine which of the study's treatment assignment designs (see `projects/IncidenceDesign/03_designs.R`) performs best under the *real* empirical incidence surface — not just the synthetic incidence surfaces used in the main simulation. This means: build real cluster-level incidence → feed it into the existing simulation/estimation machinery → compare designs on MSE/coverage/etc. exactly as the synthetic study does, but on the real data.

**Immediate sub-task (this session's focus):** get from raw county-level source data to cluster-level (58-college) incidence expressed as a rate **per 100,000 population at risk**, in three forms: year-stratified (2018–2021 separately), summed total across years, and averaged across years.

**Read `projects/IncidenceDesign/CLAUDE.md` first for full project context.**

**Critical instruction — investigate before proposing anything:** I do not yet have a confident understanding of the source data, and neither do you. Before proposing any aggregation plan, actually inspect and interrogate:

- `projects/IncidenceDesign/application/data/final_county_sudden.csv` — 400 rows (100 NC counties × 2018–2021), 42 columns. Figure out precisely what `num_obs` and `pop_18_64` represent (is `num_obs` the raw SUD death count for that county-year? is `pop_18_64` the population at risk for that same year, or a static/repeated value across years for a county? check whether it varies by year within a county). Also inventory the other ~38 covariate columns (health/demographic prevalence rates) — we don't need them for incidence-per-100k right now, but note what's there for later.
- Cross-check against what we know from provenance: this file was derived by Ashkan Habib from NC death certificates 2018–2021 (SUDDEN methodology), documented in `paper/dissertation_chapter/Dissertation_Chapter.qmd` around the "County-level SUD incidence for 2018--2021 was derived by Habib" passage — read that passage and reconcile it with the actual CSV (e.g., does summing `num_obs` across all counties/years roughly match the reported "21,147 sudden unexpected out-of-hospital deaths among working-age adults"? If not, figure out why before trusting the column).
- `projects/IncidenceDesign/application/code/cc_mapping_data.R`'s `get_cc_mapping_data()` (100-county → 58-college mapping, strict 1-to-1 partition; Bertie/Northampton/Roanoke-Chowan CC handling is already settled — don't relitigate, just use it) and `projects/IncidenceDesign/application/data/cc_name_crosswalk.csv`.
- The existing (currently generic/placeholder-oriented) pipeline in `projects/IncidenceDesign/application/code/run_application_profiles.R`: `load_real_sud_data()`, `integrate_real_sud_data()`, `build_nc_application_clusters()`. These don't yet match `final_county_sudden.csv`'s actual column names and were written for a single year only — decide with me whether to extend them or write new year-stratified/total/average variants.

**Then ask me directly about anything ambiguous**, for example (don't limit yourself to just these):
- If `pop_18_64` repeats per county across years, is it truly year-specific, or should "total" and "average" population-at-risk be computed differently than naively summing/averaging a possibly-static column?
- Should "average incidence per 100k" be computed as the average of the four yearly rates, or as (summed counts / summed population-at-risk) × 100,000 — these differ when counts or population vary a lot year to year, and I want a rationale before locking one in.
- Whether the per-100k denominator should stay `pop_18_64` (the age group the SUDDEN methodology targets) at every level (county, cluster, year, total, average), for consistency with the rest of the chapter.
- How this real cluster-level incidence surface should plug into the existing simulation/design-comparison machinery (`05_run_simulation.R`, `04_estimation.R`) — same estimator, same 6 designs, τ sweep? Or a fixed real τ if one can be estimated from the data?

Only after we've converged on answers to these should you write up a numbered implementation plan for me to approve.
