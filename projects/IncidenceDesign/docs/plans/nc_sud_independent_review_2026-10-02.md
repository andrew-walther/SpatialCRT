# Independent NC SUD study-design review — 2026-10-02

Status: completed read-only by an independent reviewer agent at the author's
explicit request. No code, production outputs or Project 1 files were changed
by the reviewer. This note preserves the findings; the implementation plan
addresses them before production.

## Verdict

The accepted study is ready for implementation planning. No change to the primary
education model or additional parameter sweep is needed. Production requires
the design-freezing, identification, reproducibility and precision safeguards below.

## Alignment with earlier implementations

- SAR construction, constant structural τ, Gaussian innovations, both spillover
  regimes and validated ML are consistent with the revised grid framework.
  `code/05_run_simulation.R:264` specifies the grid DGP; the application primary
  intentionally sets its incidence baseline contribution to zero.
- The hypothetical sensitivity's X = (average rank(rate) − 0.5)/58 follows the
  grid Poisson transformation (`code/02_incidence_generation.R:87`). Mean regional
  ranks follow that version's saturation logic; mean raw rates are an alternative
  allocation rule, not an equivalent implementation.
- Current random tie handling, exact 29-cluster SRS, approved variable budgets,
  named cached queen/rook W, and the approved Balanced Halves 14/15 correction
  are coherent irregular-geography adaptations. The reviewer checked that
  `sar_lag_setup()` accepts both cached matrices through `mat2listw()`.

## Intentional differences to explain in the manuscripts

Primary β = 0 makes observed incidence informative for allocation, without
assuming it predicts baseline education. The β = 1 sensitivity is hypothetical.
Neither analysis estimates an observed education effect or establishes SUD reduction.
Population enters planning diagnostics and some designs, but does not change
Gaussian innovation variance under the accepted model.

Primary `Y ~ Z + Spill` versus sensitivity `Y ~ Z + Spill + X` changes the
baseline DGP and corresponding nuisance adjustment. Compare designs within each
matched specification. Between-specification changes are sensitivity to the
baseline/adjustment specification, not a pure β effect. A β = 0/X-adjusted bridge
is unnecessary unless a later claim specifically tries to separate those effects;
no bridge sweep is recommended.

Graph Checkerboard is a deterministic greedy interspersed assignment, not
necessarily a proper two-coloring (`application/code/application_designs.R:98`).
Spatial blocks use standardized location, incidence and log-population with
unequal block sizes, not literal 2×2 cells (`application_designs.R:134`).
Performance comparisons include their differing treatment budgets and exposure
patterns; report counts and population shares.

With β = 0, fixed W and a fixed region map, Designs 1, 3, 4 and 9 have identical
allocation/outcome distributions across years. Year-specific population/incidence
diagnostics still differ. Verify support equality and reuse primary performance
with explicit source IDs, or explain that independently repeated annual performance
differences are Monte Carlo noise. Reused results must not be treated as independent
yearly confirmations. This equivalence does not generally hold for β = 1.

## Required implementation safeguards

1. Freeze one period-wide region map and one spatial-block map per year outside
   the scenario loop. The existing Design 5 helper rebuilds k-means blocks per
   call (`application_designs.R:235`), so the revised caller/helper must pass saved
   block membership. Fix the canonical named order: graph tie-breaking uses row order.
2. Implement and test randomized 14/15 Balanced Halves; the existing rounding
   produces 28 treated (`application_designs.R:262`).
3. Use `fit_one_lag_model()` with explicit model matrices. The existing
   `fit_tau_models()` always includes X (`code/04_estimation.R:331`); a zero X column
   in the primary model would introduce unnecessary aliasing.
4. Replace the legacy application's both-arms-only, suppressed-warning and
   `na.rm` behavior (`run_application_profiles.R:344`) in a separate revised runner.
5. Prove singleton support where possible. All four observed yearly surfaces have
   58 distinct rates and no tie at the treatment cutoff, so High Incidence Focus
   has one eligible allocation per year despite its general helper metadata.
   Graph Checkerboard also has one implemented assignment per neighbor type.
   One sampled unique allocation alone does not prove singleton support for other designs.
6. Assess structural identification with the actual fitted matrix, including
   residual treatment variation after nuisance adjustment. A finite coefficient
   after spillover is dropped does not prove identification of structural τ.
   Retain errors, failed SEs, warnings, bound hits and alias statuses.

## Monte Carlo precision

For singleton support, increasing allocation draws adds no information. Under an
approximately unbiased normal error benchmark, relative MCSE(MSE) ≈ sqrt(2/R):
R = 100/400/800 gives approximately 14%/7%/5%. At coverage 0.95, R ≈ 475 is needed
for coverage MCSE 0.01. Start singleton production with R = 1,000 and assess
empirical targets; these benchmarks are not guarantees.

For stochastic support, use joint allocation/outcome uncertainty with cached
duplicate covariance: corrected allocation variance/J + sum(w_i² v_i), where w_i
is observed draw frequency/J and v_i is the outcome MC variance of allocation i's
estimated metric. The existing MSE arithmetic is in `code/16_allocation_risk.R:149`.
Apply analogous arithmetic to coverage and bias. Increase J when allocation
variation dominates and R when outcome noise dominates. Check precision in each
reported setting, not only pooled summaries.

Prespecify extension tiers/criteria before interpreting rankings. Descriptive
simulation estimates do not require a fresh full confirmation after every
extension. Firm claims about close superiority or allocation tails should use
a final fixed replication tier chosen from the pilot and fresh seeded confirmation
batches, or remain explicitly qualified. Never extend until a preferred ranking appears.

R400 and split-halves assess outcome noise, not the limited allocation-tail sample.
At J = 100 the worst decile contains ten sampled draws. Report estimated q90 and
worst-decile means, cross-selected diagnostics and unresolved uncertainty;
no exhaustive protection guarantees or exact worst-allocation claims.

## Verification before production

Behavioral tests must cover budgets/ties, frozen partitions, named joins, both
DGP regimes, duplicate-aware uncertainty, seed prefixes/order/worker behavior,
manifest rejection and irregular-map lean/reference agreement. An exact invariant:
Isolation Buffer implies Z ⊙ WZ = 0, so both spillover formulas agree for the same
allocation and noise. Do not alter the completed grid modules/results or reuse
old application caches.

Implementation response: `nc_sud_implementation_plan_2026-10-02.md`.
