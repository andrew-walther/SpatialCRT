# The incidence coefficient τ: grid sensitivity and estimation (plan)

Author: Andrew Walther. Written 2026-11-25. Status: proposed; the spec (step 1) needs the author's
approval before any full run. Phase 4 of `docs/plans/2026-11-25-project2-phased-plan.md`.

## Notation

Documents (since 2026-11-15): β is the intervention effect on **Z**, τ the incidence coefficient
on **X**, ψ the spillover magnitude. Code (unchanged): `tau` is the effect on Z, `beta` the
coefficient on X, `gamma` the spillover magnitude. In this plan, τ always means the incidence
coefficient (`beta` in code).

## Questions

1. **Grid at τ = 0.** The grid study fixed τ = 1, so incidence shifts the outcome and treatment
   that follows incidence is partly confounded with it (the harder case). The NC case study uses
   τ = 0. Do the grid design comparisons change when incidence does not shift the outcome, and at
   an intermediate τ = 0.5? Expected: designs that follow incidence (High Incidence Focus, and to a
   lesser degree Incidence-Guided Saturation) lose less precision at τ = 0.
2. **Estimating τ.** The manuscripts treat τ as pre-specified at the design stage. After a trial,
   an investigator has both prior incidence and outcomes, and the fitted spatial lag model already
   includes incidence as a covariate on the grid, so τ̂ is estimated in every grid fit but has never
   been reported. How well can τ be estimated (bias, MSE, coverage of its interval, power to detect
   τ ≠ 0), and how does the design affect that? Designs that align treatment with incidence leave
   less independent variation in X given Z, so they may estimate τ (and β) less precisely. This
   gives investigators an estimate of the association between prior incidence and later outcomes.

This reopens the 2026-11-14 decision that estimating τ after a trial is future work (author
request, 2026-11-25); record the outcome in DECISIONS.md.

## Approach

- New script `code/31_incidence_coefficient.R` in the pattern of `code/26` and `code/29`: source
  `05` with `sim_define_only <- TRUE`, never edit `01`–`05`, new output folder
  `results/incidence_coef_2026-11/`, manifest-checked checkpoints (gitignored).
- Expectations to check, not assume: (1) β̂ barely moves with τ, because both fits include X, so
  τX is absorbed except through the estimate of ρ; a near-null answer to question 1 is a finding.
  (2) τ̂'s SE is about the same at every τ and depends on the design: roughly 0.35 when treatment
  is unrelated to incidence (X on [0, 1] has variance near 1/12, n = 100, σ = 1) and about twice
  that for High Incidence Focus, whose Z is X above the median.
- Set the code's `beta` (τ) to each value in {0, 0.25, 0.5, 1, 2} (revised 2026-11-25 from
  {0, 0.5, 1}: the extra values give a power curve for detecting τ). Assignment keys and noise keys do not
  include `beta`, so every τ value sees the same assignments and the same noise: comparisons across
  τ are paired, and τ = 1 must reproduce the 2026-09-27 production run exactly (stop otherwise).
- Record τ̂ and its standard error from the oracle fit (and the non-oracle fit) without editing
  `04`: wrap `fit_tau_models()` in the new script so it also returns the coefficient on `X` and its
  standard error from `fit_sar_lag()` (which already returns all coefficients and `rest.se`).
  Summarize τ̂ like β̂: bias, SD, MSE, coverage and power of the t + df interval (df as in `code/26`).
- All nine designs, queen and rook, both regimes, the full grid of ρ, ψ and β. Smoke and pilot
  profiles first.
- Compute: the production run took about 14 minutes on 10 cores; expect about 15 to 20 minutes per
  value of τ, so about 1.5 hours for five values.

## Interpreting τ

τ is a regression coefficient, not a correlation: with X on [0, 1] and residual SD 1, τ = 1 puts
the highest-incidence cluster's outcome one residual SD above the lowest one's before the spatial
multiplier, about the size of β = 1. It is unbounded and could be negative. The findings should
say this plainly, since the documents so far only fix τ in the data generation.

## SI Methods subsection (after the findings, drafted in the main session)

"Estimating the incidence coefficient after a trial" (author request, 2026-11-25): what τ is;
that after a trial both prior incidence and outcomes are known, so τ is estimated in the same
model; what it tells an investigator (the association between prior incidence and the later
outcome, adjusted for treatment and spillover); how the design changes its precision; limits.
Placement (curated or full SI, Chapter 3 mirror) is decided after the findings. Codex does not
draft it.

## Outputs

- `results/incidence_coef_2026-11/`: scenario and surface results with β̂ and τ̂ metrics,
  reproduction check, summaries and figures (no in-figure titles; design names from `00`).
- `docs/findings/2026-11-25-incidence-coefficient-findings.md`: method and findings, then an
  independent review (`docs/reviews/`), as for the sensitivity study of 2026-11-24.
- Rows in the docs READMEs, the AGENTS file map, ROADMAP progress log and TODO.

## Report

- β̂ by design at each τ: MSE relative to SRS by regime (as `code/15`), and whether
  the regime-based conclusions change.
- τ̂ by design: bias, MSE, coverage, power to detect τ at each nonzero value; which designs estimate τ
  well and which poorly; the relation between τ̂ precision and how closely each design's treatment
  follows incidence (for example the correlation between Z and X).
- What an investigator gains by estimating τ after a trial, and the limits (the oracle fit,
  correct spillover form, ten surfaces per configuration).

## Stops

1. After the spec (choices of τ values, metrics, figures). 2. After the smoke and pilot runs.
3. After the findings and review, before anything goes into a manuscript.
