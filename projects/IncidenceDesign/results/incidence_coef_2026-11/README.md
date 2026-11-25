# Incidence coefficient sensitivity study

Author: Andrew Walther.

The study varies the incidence coefficient τ on X while preserving the grid's treatment assignments and noise draws. Each spatial lag fit supplies the intervention estimate β̂ on Z and the incidence estimate τ̂ on X. See the [specification](../../docs/plans/2026-11-25-incidence-coefficient-spec.md) for the simulation grid and comparisons.

## Run the verification profiles

From `projects/IncidenceDesign`, set the BLAS thread limits before starting R. The second argument is the number of workers. The verification runs use three workers to leave capacity for concurrent work.

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript code/31_incidence_coefficient.R smoke 3
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript code/31_incidence_coefficient.R pilot 3
```

The smoke run uses the production test profile. The pilot supplies all ten surfaces and 25 draws per surface for the production pilot grid. At τ = 1, the pilot compares its β̂ results with the corresponding September 27 production rows. The `full` profile checks the complete production grid at τ = 1. A mismatch stops the run before the remaining τ values.

The five-value pilot took 6.9 minutes on three workers. At τ = 1, all stored β̂ metrics matched exactly across 1,440 scenario rows and 14,400 surface rows. Across 1.8 million pilot model fits, the message log contains 50,000 rook Checkerboard spillover aliases and zero other warnings or fit errors. Both coefficient estimates and their standard errors were available for every fit. The full grid has 20 times as many fits, so the pilot implies roughly 2 hours 20 minutes with three workers under similar load.

## Read the output

Each profile writes its own subdirectory. `scenarios.rds` and `surfaces.rds` use the production scenario keys, plus `Incidence_Coef` for the true τ and `Coefficient` for the fitted term, Z or X. The production key `True_Tau` denotes the true intervention effect β. `Coverage_Tdf` and `Power_Tdf` use the corrected 95% interval. `Power_Tdf` is the positive-direction rejection rate for Z and the two-sided rejection rate for X. For X at τ = 0 it measures type I error. `Paired_MSD` in the Z rows is the mean squared difference from the matched β̂ at τ = 1, with `N_Paired` giving its denominator.

`fit_messages.csv` records aliases, warnings and errors by scenario, estimator and τ. An aliased fit can still contribute an estimate when its coefficient and standard error are finite. Check the attempted, valid and missing counts alongside metric values. `run_info.rds` records the settings, software versions and timing. The pilot and full profiles write `reproduction_check.txt`. Checkpoints are ignored by Git and are reused only when the manifest matches.

## Build the summary

The full sweep took 145.8 minutes with five workers. To reproduce it, use the same one-thread BLAS settings with `Rscript code/31_incidence_coefficient.R full 5`.

From the project root after the full profile, run `Rscript code/31_incidence_coefficient_summary.R`. The script reads the scenario and surface tables and requires a passing reproduction report. It writes `full/summary/performance.csv` and the oracle figures for queen and rook contiguity. `Intervention_Grid` distinguishes the equally weighted mean over the five β values (`All_Beta`) from the β = 1 subset (`Beta_1`). Each row keeps the incidence configuration, neighbor type and spillover regime. The tables include both estimators. `full/summary/provenance.txt` identifies the summary inputs and scripts by checksum and gives the plotting package versions. The archived full `surfaces.rds` uses xz compression. Reloading verified exact identity with the original data.

## How the runner works

- **`fit_incidence_model()` and `fit_incidence_models()`** take a cluster outcome vector and the production model inputs. The wrapper reuses 04's model matrices and calls the lean fitter once per estimator. It reads Z and X by column name from the fitted coefficients and standard errors. An aliased term can disappear from these vectors, so the wrapper returns a missing value for that term and keeps the fit message. It calculates interval degrees of freedom from the terms the fitter retained.
- **`summarize_coefficient()`** takes one estimate, model-based SE and degrees-of-freedom value per draw. It reuses the production Wald summary and adds the corrected interval with half-width \(t_{0.975, N-p}\,\mathrm{SE}\sqrt{N/(N-p)}\). It counts draws with unavailable estimates or SEs, keeps their surface rows, and uses the valid fits for the reported metrics. Monte Carlo standard errors use variation between surface means.
- **`run_incidence_unit()`** takes one incidence configuration, neighbor type, outcome ρ and true τ. Each assignment uses its own surface \(X_k\). The loop reuses 05's keys and arithmetic for \(Y=(I-\rho W)^{-1}(\beta Z+S(Z)+\tau X+\varepsilon)\). That keeps the input draws paired across τ and preserves the production fit at τ = 1. The unit returns scenario and surface tables, messages, and, at τ = 1, the paired reference estimates in draw order.
- **`match_grid()` and `check_reproduction()`** check the expected keys and compare the stored and rerun Wald metrics. `dplyr::left_join()` keeps every expected key visible so a missing match can stop the run. Comparison includes the stored coverage, rejection and Monte Carlo SE columns and checks the missing-value patterns exactly.
- **`run_incidence_study()`** distributes units with `parallel::mclapply()` and checks that every worker returned a result. `digest::digest()` records the code hashes in the checkpoint manifest. A worker publishes its RDS checkpoint after writing the complete file. The parent checks keys and fit counts before writing the profile outputs.

## How the summary works

- **`summarize_incidence()`** averages the comparable scenario rows within each configuration and design. It uses a `dplyr` join to divide each design's mean MSE by the SRS mean from the same grid slice. The root mean square change in β̂ weights each scenario's squared change by its number of valid matched draws.
- **`summarize_surface_mc()`** takes the surface rows and forms one mean per surface within each slice. It pairs each surface with its τ = 1 counterpart for changes in MSE or model SE. The Monte Carlo SE is the standard deviation between the ten surface means divided by √10. It measures simulation uncertainty within an incidence configuration.
- **`write_incidence_summary()`** reads the verified RDS inputs and writes a CSV and eight vector PDFs. The `ggplot2` figures use the same axis scales across incidence configurations within each neighbor type. In the MSE plot, a ratio below one means lower MSE than SRS. In the rejection plot the 5% reference gives the nominal rejection rate at τ = 0.
