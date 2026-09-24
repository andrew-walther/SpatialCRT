# Simulation revision spec (2026-09) — step 0.5, B1

Source of truth for the revised simulation. Written 2026-09-24 from
`simulation-revision-plan.md` (decisions settled there) and the code as it stood
at commit `38a14dc`. `00_mathematical_specification.Rmd` and the chapter Methods are
rewritten from this file in Phases C and D. Fix IDs (M1–M8) follow the plan.

## 1. What is unchanged

- 10×10 grid, N = 100 clusters; rook and queen contiguity; row-standardized W.
- Outcome model (spatial lag with a spillover term, "SDM" in older docs):
  Y = (I − ρW)⁻¹ (τZ + S(Z) + βX + ε), β = 1, σ = 1.
- Spillover term S(Z): `both` = γ·WZ; `control_only` = γ·WZ ∘ (1 − Z).
- Grid: 5 incidence configs {iid; spatial ρX ∈ {0.20, 0.50}; Poisson ρX ∈ {0.20, 0.50}}
  × nb {rook, queen} × ρ {0, 0.01, 0.20, 0.50} × γ {0.5, 0.6, 0.7, 0.8}
  × regime {control_only, both} × 8 designs × τ {0.8, 1.0, 1.5, 2.0, 3.0}
  = 12,800 scenarios.
- Incidence generators: iid U(0,1); spatial = Φ((I − ρX·W_queen)⁻¹u), u ~ N(0, I);
  Poisson = rank-normalized rates (rank − 0.5)/N of C ~ Poisson(λP),
  log λ = log(35/100,000) + (I − ρX·W_queen)⁻¹u. The spatial filters always use queen W.
- Primary estimator: oracle ML spatial lag, Y ~ Z + Spill + X, with Spill = S(Z) (the true
  regime-specific spillover covariate). 95% Wald CI from the ML standard error.

## 2. Data-generating process: before and after

Notation: config c = (mode, ρX); block b = (c, nb, ρ, γ, regime); scenario s = (b, d, τ).

### Before (April 2026 runs)

- X = [X₁ … X₁₀] drawn once per (c, τ), **before** any per-scenario seed.
- Designs saw X₁ only: Z_j = design_d(X₁), j = 1..25. Checkerboard and High Incidence
  Focus were deterministic and copied 25×.
- ε = [ε₁ … ε₁₀] drawn once per (c, nb, ρ), before the per-scenario seed; shared by every
  γ, regime and design.
- Fit (j, k): Y_{jk} = (I − ρW)⁻¹(τZ_j + S(Z_j) + βX_k + ε_k), analyzed with covariate X_k.
  250 fits = 25 designs × 10 (X_k, ε_k) pairs, so each ε_k was reused 25 times, and a
  deterministic design produced only 10 distinct fits.
- cor(X₁, X_k) ≈ 0: the design and the outcome saw unrelated incidence surfaces.
- Poisson P = 1,000 per cluster: λP ≈ 0.35 expected deaths, ~5 distinct X values.

### After (this revision)

For surface k = 1..K (K = 10) and design draw j = 1..J (J = 25):

- **M1, matched surfaces.** Z_{kj} = design_d(X_k; draw j). The outcome and the analysis use
  the same X_k:
  Y_{kj} = (I − ρW)⁻¹(τZ_{kj} + S(Z_{kj}) + βX_k + ε_{kj}).
- Every fit gets its own noise column: ε_{kj}, K·J = 250 columns per block.
- Checkerboard is the only deterministic design: Z_{kj} is the same for all k, j, but each
  fit still has a distinct ε_{kj}.
- **M2.** Poisson P = 100,000 per cluster (λP ≈ 35), `pop_mode = "equal"`.

Fits per scenario: K·J = 250, as before, but now 250 distinct (Z, ε) pairs over 10 surfaces.

## 3. Designs (M3)

Ties in X are broken at random, freshly for every design draw j, with
r = `rank(X_k, ties.method = "random")`. Designs 2, 6, 7 and 8 use r (or random-tie ranks
of quadrant means); the rest are unchanged.

| ID | Design | Rule after revision | Treated |
|---|---|---|---|
| 1 | Checkerboard | (x + y) mod 2; deterministic | 50 |
| 2 | High Incidence Focus | treat r > N/2 | exactly 50 |
| 3 | Saturation Quadrants | random permutation of {0.2, 0.4, 0.6, 0.8} over 5×5 quadrants | 50 |
| 4 | Isolation Buffer | greedy random maximal independent set in nb | rook ≈ 38, queen ≈ 22 |
| 5 | 2×2 Blocking | 2 of 4 within each 2×2 block | 50 |
| 6 | Balanced Quartiles | strata = `ntile(r, 4)` (25 each); floor(25/2) = 12 per stratum plus one more in 2 random strata | exactly 50 |
| 7 | Balanced Halves | strata = `ntile(r, 2)` (50 each); 25 treated per stratum | 50 |
| 8 | Incidence-Guided Saturation Quadrants | quadrant saturations {0.8, 0.6, 0.4, 0.2} by random-tie rank of quadrant mean X_k (highest mean → 0.8) | 50 |

Balanced Quartiles treated 48 before the revision (`round(25/2)` = 12 per stratum); it now
treats exactly N/2 (user decision after the B4 pilot, 2026-09-24). `is_design_deterministic()` is TRUE for Checkerboard only.
The application copy (`application/code/application_designs.R`) gets the same tie rule.

## 4. Seeds (M4)

Every random draw is preceded by
`set.seed(digest::digest2int(key), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")`,
where `key` joins these fields with `"|"` and numeric fields are formatted with `sprintf("%.2f")`:

| Draw | Key fields |
|---|---|
| X_k (one surface) | `"X", mode, ρX, k` |
| Z_{k·} (all J draws for surface k) | `"Z", mode, ρX, nb, ρ, γ, regime, k, d` |
| ε (all 250 columns for a block) | `"eps", mode, ρX, nb, ρ, γ, regime` |

Consequences:
- Each draw depends only on its key, so results are identical sequentially, in
  parallel, and in any order.
- X_k is shared by every block of config c (Invariant 1).
- **Common random numbers:** within a block, all 8 designs share ε and X_k; each τ reuses
  the same Z and ε; the oracle and non-oracle fits use the same Y.
- Z and ε are independent across blocks, but blocks within a config share the X surfaces.

## 5. Estimand and dependence structure

For scenario s, the target is performance **averaged over incidence surfaces**:
MSE(s) = E_X E_{Z|X} E_ε[(τ̂ − τ)²], and likewise for bias, coverage and power. X is drawn
from the config's generator, Z from the design given X, and ε ~ N(0, I).

Monte Carlo structure within a scenario: 10 independent surfaces; given X_k, 25
independent (Z_{kj}, ε_{kj}) pairs. The surface is the primary sampling unit, so Monte
Carlo SEs are computed from the K = 10 surface-level means (§7), not from 250 fits treated
as independent.

Across scenarios:
- **Designs within a block** are paired through common ε and X_k (a within-block
  comparison removes that shared noise).
- **τ levels** share Z and ε. Per-τ Friedman tests are therefore not independent
  confirmations of each other.
- **Blocks within a config** share X_1..X_10. A Friedman test over blocks treats them as
  independent, which overstates the information about surface-to-surface variation. The
  surface-level robustness check (B3) addresses this: for the key pairwise claims it runs a
  paired test over 50 independent (config × k) units, 5 configs × 10 surfaces, with the
  metric for each unit averaged over that surface's blocks.

## 6. Estimators

Both fit to the same Y_{kj}; both include an intercept.
- **Oracle (primary):** Y ~ Z + Spill + X_k, Spill = S(Z_{kj}).
- **Non-oracle (M8, sensitivity):** Y ~ Z + X_k.

Engine: `fit_sar_lag()` (`code/04_estimation.R`) is the same estimator as
`spatialreg::lagsarlm(method = "eigen")`. It maximizes the concentrated log-likelihood over
ρ ∈ 1/range(eig(W)) ± machine ε, fits OLS given ρ̂, and computes SEs from the analytic
(σ², ρ, β) information matrix scaled by s². It was validated on 5,120 fits
(`results/estimator_validation/`): max |Δτ̂| 9.7e-8, relative ΔSE 1.1e-8, identical
aliasing, 0.6 vs 76 ms/fit. Almost all of `lagsarlm`'s time is a `gc()` call on exit.
`engine = "lagsarlm"` is retained for the 1% cross-check.

## 7. Metrics, Monte Carlo SEs, flags (M6/M7)

Per scenario, over the valid fits (τ̂ and SE finite):
- Bias = mean(τ̂ − τ); SD = sd(τ̂); MSE = mean((τ̂ − τ)²).
- Coverage = mean(L ≤ τ ≤ U); Power = mean(L > 0), where L, U = τ̂ ∓ 1.96·SE.
  Power is one-sided, as before.

Surface-level SEs: compute each metric m_k on surface k's valid fits;
SE(m) = sd(m_1..m_K)/√K, and intervals use m ± t_{0.975, K−1}·SE (t₉). Design differences
are paired at the surface level: SE of (m_{d,k} − m_{d′,k}) over k.

Flags:
- `N_Aliased`: fits where the engine warned "Aliased variables found".
- `N_Warn`: all other warnings (every message is written to the warning log).
- `Z_WZ_rank_deficient`: rank([1, Z, WZ]) < 3 for any Z draw of the scenario. Recorded in
  both estimator files; the non-oracle model has no Spill column to alias, but its
  Checkerboard × rook estimates share the same cause (τ is not identified apart from γ).
- Aliased estimates are **kept** (user decision M6) and flagged in exhibits.

## 8. Result schema

`results/sim_data/`, one file per estimator:
`sim_results_MLE_tau_sweep_combined_<ts>.rds` (oracle) and
`sim_results_MLEnonoracle_tau_sweep_combined_<ts>.rds`. The non-oracle filename doesn't
match `load_latest_results("MLE_tau_sweep")`. One row per scenario:

| Column | Meaning |
|---|---|
| Incidence_Mode, Rho_Incidence, Neighbor_Type, Design, Rho, Gamma, Spillover_Type, True_Tau | scenario keys (unchanged) |
| Estimator | "oracle" / "nonoracle" |
| Mean_Estimate, Bias, SD, MSE, Coverage, Power | as §7 |
| Fail_Rate, N_Valid_Est | fraction/count of fits with non-finite τ̂ or SE |
| SE_Bias, SE_MSE, SE_Coverage, SE_Power, N_Surfaces | surface-level MC SEs (§7) |
| Mean_Treated | mean Σ Z over the 250 fits |
| N_Aliased, N_Warn, Z_WZ_rank_deficient | flags (§7) |

Also, per estimator:
- `surface_results_<tag>_<ts>.rds`: one row per (scenario, k) with N_Valid, Bias, MSE,
  Coverage, Power, Mean_Treated.
- `warnings_<tag>_<ts>.csv`: one row per distinct (scenario, message) with a count.

## 9. Runner

- Work unit = (config × nb × ρ): 40 units, each covering 4 γ × 2 regimes × 8 designs
  × 5 τ × 250 fits × 2 estimators.
- Run with `mclapply(mc.preschedule = FALSE)`. If any unit returns NULL or a try-error, the
  run **stops**.
- BLAS threads = 1 in workers. This R links Apple Accelerate (vecLib), whose thread count
  is fixed at process start by `VECLIB_MAXIMUM_THREADS` (also `OPENBLAS_NUM_THREADS` /
  `OMP_NUM_THREADS` for other BLAS builds), so the runner must be launched as
  `VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript 05_run_simulation.R <profile>`
  and stops if these aren't set. The BLAS path goes in the manifest.
- Checkpoints (one per unit) go in `results/checkpoints/rev_2026-09/<profile>/`, with a
  manifest holding the hash of parameters, the hash of code files 01–05 and package
  versions. On any mismatch the runner refuses to load and stops. Pre-revision checkpoints
  move to `results/checkpoints_pre_revision_2026-09/` (git-ignored).
- Profiles: `full` (the grid in §1) and `pilot` (τ = 1; all configs, nb and regimes;
  ρ ∈ {0, 0.5}; γ ∈ {0.5, 0.8}). Pilot output goes to `results/pilot_rev_2026-09/`,
  never `sim_data/`.
