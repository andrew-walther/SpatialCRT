# Allocation-risk study — authorized plan, 2026-09-27

The author authorized the necessary pilot/analysis and subagent work after discussion
in `ch2_check_findings_2026-09-27.md`. The aim is to assess the distribution of
allocation-specific MSE for every candidate design and SRS, not to guarantee a
preferred design's ranking. Manuscript framing will follow the results.

## Scope and execution

1. Implement an isolated runner and behavior-based tests using the existing DGP,
   assignment functions, incidence surfaces and validated lean SAR ML estimator.
   Verify fixed allocations/surfaces across repeated outcomes, seed reproducibility,
   MSE and between-allocation variance calculations, and actual estimator behavior.
2. Pilot: oracle, tau=1, queen; all five incidence configurations; the first two
   existing key-seeded incidence surfaces per configuration; rho={0,0.5},
   gamma={0.5,0.8}; both spillover regimes; all eight candidate designs plus SRS.
   Start with 100 allocation draws and 100 independent outcome replicates per draw.
   Checkerboard requires one allocation per fixed surface/scenario. Diagnose
   effectively deterministic High Incidence Focus separately.
3. Assess finite-noise and allocation-sampling precision. Increase outcome or
   allocation replication in targeted follow-ups if necessary. Cross-check nested
   mean MSE against the matching slice of the existing main run. Keep pilot-specific
   estimates clearly separate from the full main-study summaries.
4. Report mean conditional MSE, between-allocation variance/SD, q90, mean of the
   worst 10%, and sampled maximum. Review design trade-offs by spillover regime and
   incidence configuration; do not silently pool those conditions. Save findings
   and limitations for future manuscript sessions.

This initial scope does not repeat the full tau sweep, rook sensitivity, or NC
application. The author authorized expansion if needed to obtain useful findings;
material changes to estimands/design rules remain discussion points.

## Estimand and precision

At fixed scenario and incidence surface X, let

    m_d(a|X) = E_epsilon[(tau_hat - tau)^2 | allocation a, X, design d].

For R independent outcomes, estimate m by the mean squared error and its Monte
Carlo variance by sample variance of squared errors divided by R. For independent
allocation draws with independent outcome batches, estimate allocation variance by

    V_corrected = sample_variance(m_hat_a) - mean(MC_variance(m_hat_a)).

Retain the raw variance, noise correction and corrected variance. Negative corrected
estimates indicate unresolved allocation variability, not proof that it equals zero.
For a known deterministic assignment, allocation variance is structurally zero but
the MSE estimate remains uncertain and may be large. Do not treat outcome replicates
or scenario blocks sharing X as independent incidence surfaces.

Do not cache duplicate allocations with reused noise unless covariance is explicitly
handled. The simple correction above requires independent noise across allocation
draw indices, including duplicate Z. Preserve assignment probabilities/draw frequency.
Check tail stability using independent outcome batches; bootstrap of noisy MSE
estimates alone does not remove ranking noise. Quantiles/maxima are estimates over
sampled allocations, not exhaustive worst-case guarantees.

Implementation clarification: duplicate assignments are cached with their original
draw frequencies. The runner uses the explicit covariance correction
`sum(f_i * (1 - f_i/n) * v_i)/(n-1)`, verified against hand-computed examples.
Early split-half diagnostics showed substantial outcome-noise contamination, so
the authorized precision follow-up extends the same allocations from 100 to 400
outcomes for designs 3, 4, 6, 8 and SRS, first surface of every configuration,
rho=0.5/gamma=0.8, both regimes (50 blocks). This is a precision extension, not an
independent replication or an increase in the number of sampled allocations.

## Files, constraints and delegation

- New runner: `code/16_allocation_risk.R`; tests: `code/tests/test_allocation_risk.R`.
- New outputs only: `results/allocation_risk/`; separate manifest/checkpoints.
- Existing modules 01–05, main-run outputs and validated checkpoints stay unchanged.
- SpillSpatialDepSim and bios-dissertation stay read-only. No manuscript edits.
- Local parallel workers with BLAS pinned to one thread; log failures, warnings,
  aliasing, source hashes and software versions. Never silently omit failed fits.
- Implementation/execution: GPT-6 Sol, high reasoning. Independent methods/code
  review: GPT-6 Sol, high reasoning. Parent handles integration, interpretation and
  documentation. No delegated commits or external messages.
- Commit code/tests, results and final documentation at logical milestones, without
  AI co-author lines. Never push without author permission.
