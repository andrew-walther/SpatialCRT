# NC SUD application: synthetic incidence and smoke/pilot evaluation

This application layer adapts the completed IncidenceDesign simulation from a 10x10 theoretical grid to 58 irregular North Carolina Community College service-area clusters. The current SUD incidence values are placeholders only. Future true incidence data should be joined using the same schema: `cluster_id`, `population`, `sud_count`, `sud_rate_per_100k`, and `incidence_rank01`.

## Synthetic placeholder incidence

The placeholder generator uses Queen contiguity on dissolved service-area polygons. It first builds a binary adjacency matrix `A`, sets the diagonal to the requested self-weight (`0` by default), and row-standardizes to obtain `W`. For non-isolated clusters, each row of `W` sums to one, so spatial lags such as `WZ` are neighbor averages.

Synthetic log risk is generated with a spatial autoregressive process:

```text
eta = (I - rho_incidence W)^-1 (log(mu) + epsilon)
epsilon_i ~ Normal(0, sigma^2)
```

The Poisson observation layer maps this spatial risk surface to positive event counts:

```text
lambda_i = exp(eta_i) population_i
SUD_count_i ~ Poisson(lambda_i)
SUD_rate_i = 100000 * SUD_count_i / population_i
```

The design-facing incidence covariate is `incidence_rank01`, the rank-normalized SUD rate. This mirrors the original simulation design logic while preserving a clean replacement path for future true incidence.

## Irregular-map design adaptations

The eight designs remain the same conceptual designs as the grid simulation, but grid-specific geometry is replaced with conservative irregular-map analogs:

- Design 1 is displayed as Block Stratified Sampling and uses a graph checkerboard assignment over the Queen-contiguity graph on the irregular NC map.
- Design 2 treats the top 50 percent by incidence rank.
- Design 3 creates four centroid k-means saturation regions and randomly assigns 20, 40, 60, and 80 percent saturation across regions.
- Design 4 greedily assigns treated clusters while blocking adjacent treated neighbors where feasible.
- Design 5 is displayed as 2x2 Blocking and uses local spatial blocking over centroid, incidence, and population features, then randomizes 1:1 within blocks.
- Design 6 randomizes 1:1 within incidence quartiles.
- Design 7 randomizes 1:1 within incidence halves.
- Design 8 uses the same four k-means regions as Design 3 and assigns higher saturation to higher-incidence regions.

For tables and figures, designs use the canonical display order from `IncidenceDesign_ProjectSummary`: Blocking designs first, then Stratified designs, then Saturation designs (`D1, D5, D4, D2, D7, D6, D3, D8`). For Designs 3 and 8, the k-means region selection is repeated across seeded candidates and scored for population balance, cluster-count balance, county-count balance, and compactness. Population defaults to the NC OSBM `county-population-totals` API, using 2024 county total population estimates. If that download is unavailable, the runner falls back to equal county weights and records `population_source = "equal_county_placeholder"` in the incidence output. The output includes `region_balance_diagnostics.csv` with `n_clusters`, `n_counties`, `population`, `population_share`, and `mean_incidence_rank01`.

## Evaluation profiles

Only `smoke` and `pilot` are executed in this phase. The full 160,000-fit grid is intentionally deferred.

| Profile | Tau | Gamma | Rho | Designs | Design resamples | Outcome resamples | Fits |
|---|---:|---:|---:|---:|---:|---:|---:|
| smoke | 1.0 | 0.6 | 0.20 | 8 | 2 | 2 | 32 |
| pilot | 1.0, 2.0 | 0.6, 0.8 | 0, 0.20 | 8 | 5 | 5 | 1,600 |
| full | 0.8, 1.0, 1.5, 2.0, 3.0 | 0.5, 0.6, 0.7, 0.8 | 0, 0.01, 0.20, 0.50 | 8 | 25 | 10 | 160,000 |

The full profile is defined in code and requires explicit confirmation with `APPLICATION_RUN_FULL=TRUE` or `allow_full = TRUE`. It should be run after smoke and pilot outputs have been reviewed.

## Outcome model

The application evaluation keeps the original SDM/SAR parameterization:

```text
Y = (I - rho W)^-1 [tau Z + gamma WZ + beta X + epsilon]
```

Current defaults are Queen neighbors only, `spillover_type = "both"`, `beta = 1`, and `sigma = 1`. Estimation uses `spatialreg::lagsarlm(Y ~ Z + Spill + X)` with the Queen `listw` object.

## Outputs

Each profile writes to `application/results/<profile>/`:

- `synthetic_incidence.csv`
- `region_balance_diagnostics.csv`
- `iteration_results.csv`
- `summary_results.csv`
- `application_<profile>_results.rds`
- `figures/synthetic_incidence_map.png`
- `figures/kmeans_regions_map.png`

The full profile also writes per-scenario checkpoint files to `application/results/full/scenario_chunks/`, allowing interrupted full runs to resume without repeating completed scenarios.

The report clearly labels smoke and pilot outputs as validation runs, not final full-grid evidence.
