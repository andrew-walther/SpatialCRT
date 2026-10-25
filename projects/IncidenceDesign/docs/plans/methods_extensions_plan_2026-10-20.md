---

editor_options: 
  markdown: 
    wrap: 72
---

# Methods extensions plan — post-advisor-meeting analyses (2026-10-20)

**Status: approved by the author 2026-10-20.** Planning only: no analysis code has been written and nothing has been run. Decisions D1–D9 were agreed with the author in the 2026-10-19/20 planning interview. "Claude's call" marks choices Claude made that the author has seen. "Builder's call" is left to whoever implements the plan.

Sources: [advisor meeting record](advisor_meeting_2026-10-14.md) (Lin's feedback, follow-up reasoning, synthesis, planned extensions), [TODO.md](TODO.md) ("Advisor meeting 2026-10-14 follow-ups"), Section 8 of `application/report/IncidenceDesign_Application_Study_Report.qmd`, `DECISIONS.md`, `application/AGENTS.md`, and the application code (`application/code/real_sud_*.R`, `code/04_estimation.R`).

------------------------------------------------------------------------

## 1. Purpose

The project gives Dr. Simpson a recommendation for how to allocate the SUD-prevention education intervention when it is scaled from a single-county pilot to a statewide community-college trial. The 2026-10-14 meeting raised five gaps between the completed study and that recommendation:

1.  what the trial's outcome is and how it relates to incidence;
2.  not every college will join;
3.  "how many colleges are needed?" (sample size and power);
4.  whether design differences are statistically supported, and how many colleges it takes for a design advantage to matter;
5.  robustness to the estimation method, plus slightly low interval coverage.

This plan closes those gaps before the prelim (Monday 2026-12-07). The author's target is to run the analyses this week; manuscript and slide edits follow.

**Goals** - Power curves, colleges needed for 80% power, and minimum detectable effect, by design and spillover regime, for the NC application. - Whether the regime guidance holds under realistic partial participation. - A targeting measure (share of SUD burden in treated areas) beside accuracy. - Formal design comparisons for NC, reported as a consistency summary. - A small estimation-method sensitivity analysis and a coverage-calibration attempt, enough for chapter appendix, SI and backup slides.

**Non-goals** - No future-SUD-incidence (count / rate-ratio) outcome model. That is future work (D1). - No partial-participation or new power runs on the 10×10 grid (D3). - No change to the completed grid or NC production results. Every new run writes to new output folders with its own manifest. - No Z·X (incidence-dependent) treatment effect.

**Anti-goals** (failure even if it runs) - Presenting Monte Carlo p-values as if they were trial-level evidence. - Pooling neighbor types, regimes or β silently. - Editing the hashed computation files of the completed runs (`application/code/real_sud_{setup,simulation}.R`, `code/01`–`05`), which would break their manifests. New code goes in new files that source existing functions. - Substituting pilot/timing-run numbers for production numbers.

------------------------------------------------------------------------

## 2. Decisions taken (author, 2026-10-19/20)

Recorded in `DECISIONS.md` (2026-10-20 entries); one-line rules in `AGENTS.md`.

**D1. Outcome and β.** - The trial's outcome Y is the measured education endpoint (an assessment of the education received by first responders and healthcare professionals), averaged within each college. It is an intermediate endpoint: we hope it relates to a future reduction in SUD incidence, but the project does not assume it. Calling Y "future incidence" would conflate process with outcome. - β = 0 stays primary. β = 1 stays a hypothetical sensitivity analysis: the education endpoint tracks prior burden in both arms. It is not a claim that the education works better in high-incidence areas. Both are reported side by side. This reaffirms the 2026-10-02 rule and adds the reasoning. - τ is constant wherever the intervention is applied (no Z·X term): there is no information to support anything else. - No count / rate-ratio future-incidence module. The manuscripts say in one sentence that designing for a future-incidence endpoint is future work. Revisit only if the upcoming meeting with Dr. Lin and Dr. Simpson asks for it.

**D2. Targeting.** Incidence guidance is justified by targeting: it sends more treatment to high-burden areas while matching plain saturation's accuracy, which is the ethical rationale. Add a burden-share measure for each design (analysis A1).

**D3. Partial participation (NC only).** - All colleges are invited after the pilot, and each volunteers or not. Full participation (58) is the best case. The full-58 comparison stays primary; partial participation is presented as "planning for realistic participation" and is the base of the power analysis. - Volunteering: random (primary). Sensitivity: the chance of joining rises with incidence rank, the highest-burden college about 3× as likely as the lowest. - Spillover exposure (rule A): a participant's exposure is γ × (treated neighbors) / (all geographic queen neighbors). Non-participants count as untreated neighbors. - Spatial dependence (ρ) links participating colleges only. A participant with no participating neighbor (an island) has a zero row in the lag W. - Designs are rebuilt on the participating set with their current rules (Section 4, A3). - Not run on the grid.

**D4. Power scale.** - Effects are standardized (Cohen's d; 0.2 small, 0.5 medium, 0.8 large). - Main case: a fixed number m of trained professionals per college. Results are translated post hoc into colleges × m × ICC tables using the college-average noise σ_c = √(ICC + (1 − ICC)/m), in individual-SD units. - Sensitivity: m proportional to the college's population (simulated, since unequal noise can change design performance). - The fitter is the fast eigen-based `fit_sar_lag()` / `fit_one_lag_model(engine =   "lean")` (`code/04_estimation.R`), validated against `spatialreg::lagsarlm()` (\|Δτ̂\| \< 1e-7; 0.63 vs 76 ms per fit; `results/estimator_validation/`). Not `lagsarlm()`.

**D5. Formal tests.** - Grid tests (Friedman / Nemenyi / Holm-Wilcoxon) stay in the chapter appendix and SI, with one sentence in the main Results, and get a prelim backup slide. - NC tests are added. The four years are four related observed surfaces, not independent draws, so the p-values are labeled a consistency summary. - The headline comparison is the number of colleges each design needs for 80% power.

**D6. Placement.** See Section 7. CTJ: power and sample size go in the main text as prose, with the figure in the SI (the six main exhibits stay). Anything a committee member could ask about gets a main or backup slide.

------------------------------------------------------------------------

## 3. Decisions on the sensitivity analyses (author, 2026-10-20)

- **D7. Coverage:** try every reasonable correction (t, df-corrected variance, percentile and studentized parametric bootstrap, profile likelihood), with the 94–96% adoption rule. Fallback: explain the issue, what was tried and the result. Backup slide either way.
- **D8. Estimation sensitivity:** CAR and Bayesian CAR, expanded beyond a toy example for robustness. `CARBayes` install approved.
- **D9. Spillover forms:** proportional (main) plus binary and additive-count sensitivities. No second-order (neighbors of neighbors) spillover.
- **Compute:** the A3/A4 power runs use the fast fitter and run overnight locally.

------------------------------------------------------------------------

## 4. Analyses

Shared conventions (all analyses): - Key-seeded draws via the existing seed-key pattern, with new key prefixes (e.g. `"pp-subset"`, `"pp-Z"`, `"pp-eps"`). The effect size is **not** in the noise or allocation keys, so all effects share allocations and noise (common random numbers), which gives smooth power curves. - Queen primary. β, regime, participation mechanism and year are never pooled silently. Results are reported per regime; years are summarized as mean and range. - New outputs go under `application/results/methods_ext_2026-10/<analysis>/`, each with a manifest (parameter hash, code hashes, package versions, BLAS), as for production. - Fail loud: count and report failures, aliased fits, unidentified allocations, islands per subset, and boundary ρ̂ fits. Never drop them silently. - Threads: `VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`, forked workers (10 local cores).

### A0. Quick numbers from stored production fits (no new simulation)

- Power at τ = 1 from the stored NC fits (τ̂, SE), per design / regime / year / β. This is expected to be near 100%, and documents why new runs are needed.
- Outputs: `a0_stored_power.csv`. Compute: minutes.

### A1. Burden share treated (D2; no fitting)

- For each design, year and allocation: share of that year's SUD deaths, and of the person-years, in treated clusters. Summarize the mean and the 10th/90th percentiles across allocations. Use the production allocations (J = 100), N = 58.
- Also at each participation level from A3 (share of the participants' deaths).
- Expected sanity values: SRS mean ≈ treated fraction; High Incidence Focus highest; Incidence-Guided Saturation above Plain Saturation.
- Outputs: `a1_burden_share.csv`, one figure. Compute: minutes.

### A2. NC formal design comparisons (D5; stored results)

- Units: matched settings year × ρ × γ (4 × 3 × 2 = 24 blocks), separately for each regime and each β. Nine designs including SRS.
- Friedman test on setting-level mean MSE; Nemenyi post hoc; Holm-adjusted paired Wilcoxon for all 36 pairs. Adapt `code/12_six_design_statistical_comparisons.R` (and helpers in `code/10_*`) in a new NC script; do not edit 10/12.
- Label the results a consistency summary (shared map, correlated years, allocation streams shared across equivalent supports).
- Also run them per participation level from A3 (N = 20, 29, 58) for the "does the guidance hold" question. Those use A3's 16 blocks (year × ρ ∈ {0, 0.5} × κ), not the production 24.
- Outputs: `a2_tests_*.csv`, critical-difference plots. Compute: minutes.

### A3. Partial participation + power (D3, D4; main new simulation)

**DGP** (queen geography, N participants P ⊂ {1..58}): - Y_P = (I − ρ W_P)⁻¹ [τ Z + S(Z) + β X_P + σ_c ε], with ε \~ N(0, I). - W_P: the queen contiguity among participants, row-standardized; islands get zero rows. - S(Z) = κτ · e(Z), where e_j = (Σ\_{k ∈ nb(j)} Z_k) / \|nb(j)\|, over **all** geographic queen neighbors (rule A; non-participants have Z = 0). Times (1 − Z_j) under control-only spillover. - Spillover scales with the effect, γ = κτ with κ ∈ {0.5, 0.8}. This matches production at τ = 1 (γ = 0.5, 0.8). - Main runs: β = 0, σ_c = 1, so τ is the standardized cluster-level effect.

**Fitted model:** intercept, Z, Spill and X when β = 1. Spill is the *unscaled* exposure e(Z) (× (1 − Z) under control-only), whose coefficient estimates κτ. The analyst can compute it, because which neighbors were treated is known. (Claude's call, from the builder check: a κτ-scaled column would be all zeros at τ = 0 and get dropped as aliased, so the τ = 0 fits would use a different model; τ̂ is identical for τ ≠ 0.) Spatial lag on W_P via `fit_one_lag_model(engine = "lean")`, with the eigen setup recomputed per subset. Wald 95% interval and two-sided test of τ = 0.

**Settings:**

| Factor | Values |
|------------------------------------|------------------------------------|
| Participants N | 15, 20, 25, 29, 30, 35, 40, 45, 50, 58 (29 = half the colleges) |
| Mechanism | random at all N; incidence-weighted (3×) at N = 20, 29 |
| Effect τ (cluster SD) | 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0 |
| κ (spillover / τ) | 0.5, 0.8 |
| ρ | 0, 0.5 |
| Regime | both, control-only |
| Designs | all 8 + SRS |
| Year (planning surface) | 2018–2021 |
| β | 0 (main); 1 at N = 29, 58 with τ ∈ {0, 0.4, 1.0} (on the main grid, so common random numbers apply) |

Incidence-weighted sampling: successive draws without replacement with weights w_j ∝ 1 + 2·(r_j − 1)/57, where r_j is the year's incidence rank, so the top/bottom weight ratio is 3. (Builder's call: an equivalent scheme is fine if the realized inclusion-probability ratio is reported.)

**Replicates per setting:** for N \< 58, 50 participant subsets × 20 allocations × 10 outcomes = 10,000 fits. For N = 58, 1 subset × 100 allocations × 100 outcomes, matching production. Deterministic designs (Checkerboard; High Incidence Focus without ties): 1 allocation per subset × 200 outcomes for N \< 58, and R = 1,000 at N = 58 as in production. Subsets are shared across designs (key excludes design), so designs are compared on the same participant sets. Under random participation the subset key also excludes the year, so year-to-year differences aren't mixed with subset variation. Incidence-weighted subsets are year-specific, since the weights use that year's ranks. (Claude's calls, from the builder check.)

**Designs on a subset** (rebuild with the existing functions in `application/code/application_designs.R`). At N = 58, reuse the frozen production regions and yearly blocks (`setup.rds`), so that V1 is an exact reproduction check. Rebuild them only for N \< 58 (Claude's call): - Saturation regions (Plain and Incidence-Guided): `make_kmeans_regions()` on the participants, 4 regions. Incidence-guided ordering uses the participants' mean incidence rank per region, as in production. - 2×2 Blocking analogue: `make_spatial_blocks()` on the participants, target size 4. - Checkerboard: `graph_checkerboard_assignment()` on the participant neighbor graph. - Isolation Buffer: maximal independent set on the participant graph. - High Incidence Focus, Balanced Quartiles, Balanced Halves, SRS: ranks among participants. The fixed-budget designs all treat ⌊N/2⌋ at odd N. This replaces the existing `round(N/2)` in High Incidence Focus and SRS, which rounds half to even, and is identical at N = 58. The rule-based designs (saturation, Isolation Buffer, Checkerboard) keep their own budgets, as approved on 2026-10-02 (unequal budgets reported). (Claude's call.) High Incidence Focus is deterministic for a subset when the participants' rates have no ties. - Islands: under Isolation Buffer an island is always treated, and under Checkerboard it is always control (unless the colouring flips). Both are properties of those designs, reported, not changed. Other designs randomize islands normally.

**Outputs:** - Per setting: power (Wald), empirical type I error (τ = 0), size-adjusted power (the critical value is the 95th percentile of \|t\| at τ = 0 for that design/setting), MSE, bias, coverage, failure counts, island counts, and MC SEs that treat the subset as the top-level unit. - Do not store all per-fit rows (about 350 million). Store per-setting summaries plus per-(subset, allocation) summaries. Keep the \|t\| values at τ = 0 so size adjustment can be recomputed. - Derived: N₈₀ (colleges needed for 80% power) by design × regime × effect, by monotone interpolation over N; MDE(N) (effect at 80% power) by design × regime × N; bootstrap intervals over subsets. - Translation tables: colleges × m ∈ {10, 25, 50} × ICC ∈ {0.01, 0.05, 0.10} → power for individual effects d ∈ {0.1, 0.2, 0.3, 0.5}. These use σ_c and are exact under β = 0, since ML is scale-equivariant (checked in V3). - "Guidance holds?" table: mean-MSE ratio to SRS by design × regime at N = 20, 29, 40, 58, from the τ = 1 runs, with A2 tests per N.

**Compute:** about 25,900 main settings (10 N values) plus about 5,200 informative and 1,700 β = 1, at about 10,000 fits each, ≈ 330 million fits. The measured fit time is 0.63 ms at n = 100 (equivalence log); assume about 0.8 ms per fit, including allocation and outcome generation. That is ≈ 75 core-hours: **about 7–8 h wall-clock on 10 local cores (one overnight run)**, or about 30–60 min on Longleaf (200 cores, plus queue). Plan: local, with checkpoints by (year, N, mechanism). Longleaf only if the timing pilot exceeds 1.5 ms per fit. A SLURM array over (year, N) is the natural split; memory per worker is small (n ≤ 58).

**Timing and smoke pilot first:** 1 year, N ∈ {20, 58}, τ ∈ {0, 1}, κ = 0.8, ρ = 0.5, reduced replicates. Report ms per fit and project the full wall-clock before launch.

### A4. Unequal trained professionals (D4 sensitivity)

- m_j = m̄ · pop_j / mean(pop) among participants, with m̄ = 25 and ICC = 0.05, so σ_j = √(ICC + (1 − ICC)/m_j). Heteroscedastic noise in the DGP; the fit is unchanged (unweighted SAR), as a naive analyst would fit it. (Builder's call: also report an inverse-variance-weighted note if trivial; not required.)
- N ∈ {20, 29, 40, 58}, random participation, individual d ∈ {0, 0.1, 0.2, 0.3, 0.4}, other factors as A3 main.
- Compare to the equal-m translation at the same m̄: does any design gain or lose?
- Compute: about 5,800 settings, ≈ 58 million fits, ≈ 1.5 h local.

### A5. Coverage calibration (D7)

- **Stage 1 (stored fits, no runs).** From the NC production fits, intervals with:

  (a) a t₅₄ critical value (done 2026-10-14: 93.1% → 93.7%); (b) a df-corrected residual variance, SE × √(n / (n − p)); (c) both. If the grid's per-fit τ̂/SE caches exist, repeat on the grid (queen, τ = 1); otherwise state that the grid was not checked.

- **Stage 2 (small new run).** On the subset N = 58, β = 0, four years, ρ ∈ {0, 0.5}, γ = 0.8, both regimes, nine designs (144 settings) × 500 fits:

  - parametric percentile bootstrap and parametric bootstrap-t, each with B = 499 outcomes simulated from the fitted model and refitted;
  - profile-likelihood (likelihood-ratio) interval for τ, maximizing over ρ and the other coefficients at each fixed τ.

  About 36 million bootstrap fits plus the profile fits, ≈ 1–1.5 h local.

- **Adoption rule (author-approved):** adopt a method if its coverage is within 94.0–96.0% for every design and regime and the design rankings are unchanged. Report it in the appendix/SI and recommend it for the trial analysis. Otherwise keep Wald as the reported interval and write the limitation: close to nominal, about 93%, with what was tried and the results.

- A backup slide either way: the issue, the methods tried, the result.

- Power (A3) is reported as Wald plus size-adjusted power regardless, so A3 does not wait on A5.

- Note: `fit_sar_lag()` reproduces lagsarlm's asymptotic covariance, which already includes ρ̂. The TODO's "SE ignores ρ̂ uncertainty" explanation is untested; Stage 1(b) checks the small-sample variance explanation instead.

### A6. Estimation-method sensitivity (D8)

- Data from the production SAR DGP at N = 58, β = 0. Every fit uses the same mean (intercept + Z + Spill) on the full 58-college W.
- **CAR** (`spatialreg::spautolm(family = "CAR")`) over the full production grid: four years × ρ ∈ {0, 0.2, 0.5} × γ ∈ {0.5, 0.8} × both regimes × nine designs (432 settings) × 500 fits (25 allocations × 20 outcomes), with SAR lag fitted to the same data as the reference. At about 50 ms per CAR fit, 216,000 fits ≈ 20 min wall.
- **Reverse misspecification:** data from a CAR DGP (same mean; the CAR dependence parameter on the same grid as ρ), fitted by SAR lag and by CAR. Same settings and cost.
- **Bayesian CAR** (`CARBayes::S.CARleroux`, Gaussian; install approved): four years × ρ ∈ {0, 0.5} × γ = 0.8 × both regimes × nine designs (144 settings) × 100 fits = 14,400 MCMC fits. Burn-in and iterations are the builder's call (start at 2,000 + 5,000, thin 1). Check convergence (Geweke / effective sample size) on a random 5% of fits and report it. About 5–10 s per fit ≈ 2–3 h wall; run after A3. Posterior-mean τ̂ and 95% credible-interval coverage.
- Outputs per estimator: MSE ratio to SRS by design × regime, Spearman correlation of the design orderings with SAR, bias, coverage, and failure / non-convergence counts.
- Existing related evidence to cite alongside: grid non-oracle and omitted-spillover fits.

### A7. Spillover-form sensitivity: binary and additive count (D9)

Both use N = 58, β = 0, four years, ρ ∈ {0, 0.5}, γ ∈ {0.5, 0.8}, both regimes, nine designs, J = 100 × R = 100. DGP and fit use the same exposure form (correctly specified). Proportional (the main analysis) is the comparison. - **Binary (Project 1-style):** S(Z) = γ · 1{any treated queen neighbor} (× (1 − Z) under control-only). - **Additive count:** S(Z) = (γ / d̄) · (number of treated queen neighbors), where d̄ = 4.93 is the mean queen degree. A college of average degree with every neighbor treated gets γ, so magnitudes are comparable to the proportional form. (Claude's call.) The author notes this form is less plausible for an education intervention; it is included to show what the recommendation would look like. - Report aliased or unidentified designs as such. Expected: under binary, control-only spillover, Isolation Buffer and Checkerboard are not identified. Never drop them. - Outputs: MSE ratio to SRS by design × regime for each form, beside the proportional results. Does the regime guidance change? - Compute: about 5.8 million fits, ≈ 10 min.

## 5. Verification checks (each must pass before results are used)

- **V1 Reproduction:** A3 at N = 58, random participation, τ = 1, κ ∈ {0.5, 0.8}, β = 0 reproduces production mean MSE and coverage per design × regime × ρ × year within 2 MC SEs. Rule A at N = 58 is the production spillover exactly. With frozen partitions and production's R, differences beyond MC error indicate a bug.
- **V2 Size:** at τ = 0, the Wald rejection rate is about 5–8% (consistent with about 93% coverage) for every design. Size-adjusted rejection is 5% by construction; report its MC SE.
- **V3 Scale invariance:** for β = 0, one setting run at (τ, σ_c) = (0.5, 0.5) gives the same τ̂/SE as (1, 1) to floating-point tolerance. This validates the m × ICC translation.
- **V4 Rule A unit test:** a toy map where a control college with 4 neighbors (2 non-participants, 1 treated, 1 control) gets exposure 0.25; an island gets 0.
- **V5 Islands in the fitter:** `fit_sar_lag()` on a W with zero rows matches `lagsarlm(zero.policy = TRUE)` on ≥ 1,000 fits (\|Δτ̂\| \< 1e-6). `sar_lag_setup()` stops when `can.be.simmed(listw)` is true, so check that the subset listw objects (built with `mat2listw(..., zero.policy = TRUE)`) pass that guard, not only that zero rows work.
- **V6 Designs on subsets:** for each design on 1,000 random subsets: budgets as specified; Isolation Buffer has no two adjacent treated participants and is maximal; 4 regions with ≥ 2 members each, recording how often a region falls below 3; the checkerboard is deterministic given the subset.
- **V7 Informative sampling:** the realized inclusion-probability ratio (top vs bottom incidence decile) is about 3.
- **V8 Monotonicity:** power is non-decreasing in τ and in N within MC error; flag violations rather than smoothing them away.
- **V9 Manifests:** the hashes of existing production computation files are unchanged after the work (compare to the production manifest).
- **V10 A6 sanity:** CAR and SAR agree closely at ρ = 0. Bayesian posterior means are close to the ML estimates on a few hand-checked fits. Non-converged MCMC fits are counted and reported.
- **V11 A7 exposure unit tests:** on a toy map, the binary and additive exposures match hand calculations. At N = 58, the proportional form equals production.
- Tests live in new files under `application/tests/` (e.g. `test_methods_extensions.R`), runnable without the restricted data where possible.

------------------------------------------------------------------------

## 6. Order of work

| Step | Work | Depends on | Real time |
|------------------|------------------|------------------|------------------|
| 1 | A0, A1 (N = 58), A2 (stored), A5 Stage 1 | — | ½ day |
| 2 | A3/A4 runner + V3–V7 tests; timing/smoke pilot | — | 1 day |
| 3 | A3 + A4 production (overnight) | 2 | overnight |
| 4 | A5 Stage 2, A6 CAR + reverse, A7 (after 3); A6 Bayesian (next night) | 3 (cores) | 1 day + 1 night |
| 5 | A3 summaries: N₈₀, MDE, translation tables, guidance-holds table, A1/A2 per N; V1, V2, V8 | 3 | ½ day |
| 6 | Findings doc `docs/plans/methods_extensions_findings_2026-10-<dd>.md`, teach-back | 1–5 | ½ day |
| 7 | Manuscript edits (separate approval): chapter → CTJ (both layouts) + SI → verifiers → Chapter 3 sync | 6, author review | 1–2 days |
| 8 | Prelim slides in `bios-dissertation/slides/` (main + backup) | 6 | separate session |

Analyses are realistic within this week (steps 1–6); manuscripts and slides follow, well before 2026-12-07. Commit after each step (checkpoints), with the author's OK.

------------------------------------------------------------------------

## 7. Where each result goes

Manuscript edits route through `paper/thesis_chapter/Thesis_Chapter.qmd` first, then CTJ (`CTJ_Manuscript.tex` / `CTJ_Exhibits.tex`; rebuild `CTJ_Manuscript.pdf` and `CTJ_Reading.pdf`) and the SI. Run `Rscript code/tests/test_manuscript_exhibits.R` and `python3 paper/tools/verify_manuscripts.py` after each, and recount CTJ words (body limit 3,500; currently 3,362).

| Result | Chapter 3 | CTJ | SI | Prelim | Report / companion |
|------------|------------|------------|------------|------------|------------|
| D1 endpoint wording (education endpoint; β = 1 hypothetical; future-incidence endpoint as future work) | Methods + Discussion | Methods + Discussion | — | main | Section 8 |
| A3 colleges needed, MDE, translation tables | main Results | prose in main text | figure + tables | main | new section |
| A3 guidance holds under partial participation | main (short) + appendix | one paragraph | detail | main | new section |
| A1 burden share | main | column/sentence in existing table | detail | main | yes |
| A2 NC tests; grid tests | one sentence + appendix | one sentence | detail | backup | — |
| A4 unequal professionals | appendix | — | brief | backup | — |
| A5 coverage (issue, attempts, result) | appendix + limitation | limitation | detail | backup | — |
| A6 estimation sensitivity | appendix | one sentence | detail | backup | — |
| A7 binary / additive spillover | appendix | — | detail | backup | — |

Also update the application report's Section 8 ("Planned extensions") once results exist, from planned to done.

------------------------------------------------------------------------

## 8. Dependencies outside this plan

- Dr. Lin's sign-off on the 2026-10-02 application design choices and the other open prep questions (meeting record), at the upcoming meeting. If the meeting reverses D1 (e.g. asks for a count endpoint), A3's DGP changes. Step 3 can still run first; its β = 0 results stay valid as the education-endpoint case.
- Habib / Simpson approval of the data use and the IRB determination before any publication of the new NC results. Analysis can proceed now.
- Manuscript edits for Table 1 grouping, the spillover-form sentence, and the Habib citation are tracked in TODO.md, outside this plan.
