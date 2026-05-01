# Full application run notes

The full profile is now configured but was not executed during the smoke/pilot implementation phase.

Launch command from the `IncidenceDesign` project root:

```bash
APPLICATION_RUN_FULL=TRUE Rscript application/run_application_profiles.R full
```

or:

```bash
bash application/run_full_application.sh
```

The full profile runs:

- `tau = c(0.8, 1.0, 1.5, 2.0, 3.0)`
- `gamma = c(0.5, 0.6, 0.7, 0.8)`
- `rho = c(0, 0.01, 0.20, 0.50)`
- Queen contiguity only
- `spillover_type = "both"`
- all 8 application designs
- 25 design resamples
- 10 outcome resamples
- 160,000 MLE fits total

The runner writes per-scenario checkpoint CSVs to:

```text
application/results/full/scenario_chunks/
```

This makes the run resumable. Re-running the same command skips completed scenario chunks by default and rebuilds the aggregate `iteration_results.csv`, `summary_results.csv`, and `application_full_results.rds` outputs from available chunks plus newly completed chunks.

The full profile requires explicit confirmation through `APPLICATION_RUN_FULL=TRUE` or `run_application_profile("full", allow_full = TRUE)` to avoid accidentally starting the 160,000-fit job.

