# Methods extensions — findings (2026-10-20)

**Status: complete (steps 1–6), 2026-10-20 12:40.** Implements steps 1–6 of
[methods_extensions_plan_2026-10-20.md](methods_extensions_plan_2026-10-20.md)
(D1–D9). Manuscripts, chapter, CTJ, SI and slides are untouched (steps 7–8 need
separate approval). All outputs are under
`application/results/methods_ext_2026-10/`; fit caches are gitignored.

## Summary (for the author)

- **Ran everything in steps 1–6**: 400.3 million A3/A4 fits, 0.43M A6 data sets fit two
  ways, 14,400 MCMC fits, 36M bootstrap refits, and 6.9M A7 fits. No failed tasks; no
  failed, boundary or error SAR fits in A3/A4/A5/A7. Production files unchanged (V9).
- **Colleges needed (size-adjusted, 80% power, cluster-SD effect τ):** under both-arms
  spillover about 36 colleges for τ = 1 and 19 for τ = 1.5, the same for SRS, the rank
  designs and both saturation designs. Checkerboard and Isolation Buffer need far more
  (> 58 at τ = 1). Under control-only spillover τ = 1 is not reachable with 58 colleges;
  at τ = 1.5 the saturation designs need 28–30 colleges vs 33 for SRS. In individual
  terms (m = 25, ICC = 0.05) d = 0.3 gives 0.70 power at 29 colleges and 0.95 at 58
  (SRS, both-arms).
- **Guidance under partial participation:** Checkerboard and Isolation Buffer stay
  poor; the saturation designs' control-only advantage (22% lower MSE at N = 58)
  shrinks to 5% at N = 29 and vanishes at N = 20 or under incidence-weighted
  volunteering.
- **Robust:** the regime guidance holds under CAR fitting, a CAR DGP, the additive-count
  spillover form and, at the ends, Bayesian CAR. Binary spillover makes Isolation Buffer
  (and Checkerboard under control-only) unidentified.
- **Coverage:** t₅₄ + df-corrected SE reaches 94.1–95.0% by design × regime on the
  production fits, but no method meets 94–96% in every cell in Stage 2 (bootstrap-t 94.7%
  overall).
- **Targeting (A1):** Incidence-Guided Saturation barely out-targets Plain Saturation on
  the NC map. Death share fails the plan's expected sanity ordering (High Incidence
  Focus treats 35% of deaths). This needs an author decision on D2 wording.
- **Failed checks:** V2 (Wald size 8–13% at N ≤ 35; hence size-adjusted headlines), V6
  (tiny rebuilt saturation regions at N = 15–20; immaterial), V8 (drops along N from
  shared subset sets; none along τ), V5 guard (symmetric W_P handled and validated).
  Bayesian CAR convergence fails in 30% of chains (mostly ρ = 0.5).
- **Open author decisions:** `TODO.md` ("Methods-extension decisions"). Steps 7–8 not
  started.


---

## 1. What was run

| Step | Analysis | Output folder | Status |
|---|---|---|---|
| 1 | A0 stored power, A1 burden share (N = 58), A2 NC tests, A5 Stage 1 | `stored/` | done |
| 2 | A3/A4 runner, V3–V7/V11 tests, timing pilot | `participation/pilot/`, `verification/` | done; V6 failed (Section 4) |
| 3 | A3 + A4 production | `participation/{a3,a4}/` | done: 576/576 tasks, 400.3M fits |
| 4 | A7, A5 Stage 2, A6 CAR + reverse, A6 Bayesian CAR | `a7/`, `a5/`, `a6car/`, `a6bayes/` | done 08:09–12:38; all tasks OK |
| 5 | Summaries, V1, V2, V8, V9, V10 | `summaries/`, `sensitivity_summaries/`, `verification/` | done |

Timing pilot (2018, N ∈ {20, 58}, τ ∈ {0, 1}, κ = 0.8, ρ = 0.5, reduced replicates):
18/18 tasks OK, no failures. Production-sized timing (one process): 0.22 ms/fit at
N = 15, 0.26 at N = 29, 0.32 at N = 58, all below the 1.5 ms Longleaf threshold. Under
10 concurrent workers the rate is about 0.57 ms/fit at N = 58 (6 performance + 4
efficiency cores, shared memory bandwidth). Projected: about 5 h for A3 (~308 million
fits) and 1.5 h for A4 (~104 million fits, both A4 arms).

---

## 2. Step 1 results (stored production fits, no new simulation)

### A0. Power at τ = 1 from the stored NC fits
`stored/a0_stored_power.csv`. Not near 100%, contrary to the plan's expectation.
With τ = 1 in cluster-SD units and N = 58, two-sided Wald power is 0.96–0.97 for the
fixed-budget, saturation and blocking designs under both-arms spillover, but only
0.48–0.68 under control-only spillover. Checkerboard (0.60 both-arms / 0.20 control-only)
and Isolation Buffer (0.53; it treats about 27% of colleges) are lowest. High Incidence
Focus drops to 0.52 (both) and 0.33 (control-only) under β = 1. Power is flat across
ρ and γ (SRS 0.96 both, 0.54 control-only). No excluded (boundary/error) fits.
So power is not uninformative at N = 58. The A3 grid is needed for colleges and effect
size, not because τ = 1 is trivially detectable.

### A1. Burden share treated (N = 58, production allocations, J = 100 per year)
`stored/a1_burden_share.csv`, `stored/a1_burden_share.{pdf,png}`. Means over years:

| Design | Death share | Person-year share | Enrichment (deaths/PY) | Rate share | Treated / control mean rate |
|---|---|---|---|---|---|
| High Incidence Focus | 0.354 | 0.229 | 1.55 | 0.622 | 148.6 / 90.3 |
| Incidence-Guided Saturation | 0.492 | 0.495 | 1.01 | 0.511 | 122.0 / 116.8 |
| Plain Saturation | 0.496 | 0.491 | 1.02 | 0.504 | 120.2 / 118.6 |
| Balanced Halves | 0.500 | 0.498 | 1.01 | 0.501 | 119.7 / 119.2 |
| Balanced Quartiles | 0.498 | 0.497 | 1.01 | 0.500 | 119.5 / 119.4 |
| Simple Random Sampling | 0.497 | 0.495 | 1.02 | 0.500 | 119.5 / 119.4 |
| 2x2 Blocking | 0.506 | 0.509 | 0.99 | 0.492 | 118.4 / 120.4 |
| Checkerboard | 0.453 | 0.505 | 0.90 | 0.431 | 110.6 / 127.1 |
| Isolation Buffer | 0.265 | 0.267 | 1.02 | 0.274 | 119.4 / 119.4 |

**The plan's expected sanity ordering fails for the death-count share** and holds only
for the per-college rate share. High-incidence colleges are small: the 29
highest-rate colleges hold about 23% of the person-years, so High Incidence Focus treats
35% of deaths (less than SRS's 50%) while reaching 1.55 times more deaths per
person-year treated. Incidence-Guided Saturation is at most marginally above Plain
Saturation (rate share 0.511 vs 0.504; death share 0.492 vs 0.496) because the frozen
four production regions have similar mean incidence ranks (about 28–34 of 58 in every
year), so ordering saturations by regional mean rank has little to target. This bears
on D2 (the targeting rationale for incidence guidance) and needs the author's
attention before any manuscript wording: on the NC map, the targeting benefit of
Incidence-Guided Saturation over Plain Saturation is small. Sources: death share =
Σ deaths·z / Σ deaths; rate share = Σ rate·z / Σ rate; both per allocation, then mean
and 10th/90th percentiles.

### A2. NC formal design comparisons (consistency summary)
`stored/a2_tests_{friedman,avg_ranks,pairwise}.csv`, `stored/a2_cd_plots.pdf`.
24 blocks (4 years × 3 ρ × 2 γ) per regime and β, nine designs. **These p-values are a
consistency summary, not trial-level evidence**: one map, four correlated years, and
allocation streams shared across equivalent supports.

| Setting | Friedman χ² (df 8) | p | CD | Best three (average rank) |
|---|---|---|---|---|
| β = 0, both-arms | 135.3 | 2e-25 | 2.45 | High Incidence Focus 1.88, Balanced Halves 3.17, Plain Saturation 3.83 |
| β = 0, control-only | 153.8 | 3e-29 | 2.45 | Plain Saturation 1.67, Incidence-Guided Saturation 1.79, Isolation Buffer 3.79 |
| β = 1, both-arms | 150.7 | 1e-28 | 2.45 | Balanced Quartiles 2.17, Balanced Halves 2.46, Plain Saturation 3.58 |
| β = 1, control-only | 179.8 | 1e-34 | 2.45 | Plain Saturation 1.29, Incidence-Guided Saturation 1.71, Isolation Buffer 3.50 |

Holm-adjusted Wilcoxon: 28, 30, 25 and 34 of the 36 pairs have p < 0.05 in the four
settings. Checkerboard ranks last or next to last everywhere; Isolation Buffer is last
under both-arms and third under control-only.

### A5 Stage 1. Coverage corrections on the stored fits
`stored/a5_stage1_coverage.csv`, `stored/a5_stage1_coverage_by_design_regime.csv`.
n = 58, p = mean coefficients + ρ (4 for β = 0, 5 for β = 1).

| Interval | Design × regime mean coverage | Per setting | 94–96% rule |
|---|---|---|---|
| Wald (z) | 92.3–93.6% | 91.5–95.0% | not met |
| (a) t₅₄ critical value | 93.0–94.1% | 92.2–95.4% | not met |
| (b) SE × √(n/(n − p)) | 93.5–94.5% | 92.8–95.9% | not met |
| (b′) as (b) with p = mean coefficients only | 93.3–94.3% | 92.4–95.7% | not met |
| **(c) t₅₄ and SE × √(n/(n − p))** | **94.1–95.0%** | 93.3–96.3% | **met** |

Correction (c) meets the adoption rule at the design × regime level (every design,
both regimes, both β); design rankings are unchanged because τ̂ is unchanged. Per
setting it still ranges 93.3–96.3%. Stage 2 (bootstrap, profile likelihood) results:
Section 4. The **grid was not checked**: `results/checkpoints/rev_2026-09/` keeps
surface-level summaries, not per-fit τ̂/SE.

---

## 3. A3/A4 results (participation and power)

**Run.** A3: 432 tasks, **301,824,000 fits**, 03:17–06:56; A4: 144 tasks,
**98,496,000 fits**, 06:56–08:07 (36.2 + 11.7 core-hours, 0.43 ms/fit under load). All
576 tasks OK. **No missing estimates or intervals, no boundary-ρ fits, no errors, no
unidentified allocations.** 38,160 A3 fits (0.012%) and 2,400 A4 fits had an aliased
Spill column, all at N = 15–20 (subsets/allocations where no control had a treated
neighbor); τ̂ is still estimable there and the fits are kept and counted. Participant
islands: mean 4.1 per subset at N = 15 (all 50 subsets have one), 1.6 at N = 29, 0.5 at
N = 40, 0.2 at N = 50; one N = 15 subset used the symmetric-W_P setup (V5 note).
Tables: `summaries/a3_settings.csv` (every setting), `a3_headline_tables.md`,
`a3_translation_headline.md`, `a3_n80_summary.csv`, `a3_mde_summary.csv`,
`a3_translation_m_icc.csv`, `a3_guidance_ratio_to_srs.csv`, `a2_tests_per_N_*`,
`a4_summary.csv`.

**Size first (V2).** The Wald test is anti-conservative with few colleges: rejection
at τ = 0 averages 11–13% at N = 15, 10% at N = 20, 8.3–8.8% at N = 29, 7.4–7.7% at
N = 40, 6.3–7.0% at N = 58 (by design). This is the small-sample face of the same
problem as the ~93% coverage (A5). **Wald power therefore overstates power at small
N, and the headline colleges-needed numbers below are size-adjusted** (critical value
= 95th percentile of |t| at τ = 0 in the same design/setting; Wald numbers are in the
same files). For example, at τ = 1 under both-arms spillover the Wald N₈₀ is 30–32 colleges and the
size-adjusted N₈₀ 35–37.

### Colleges needed for 80% power (size-adjusted, random volunteering, κ = 0.8)
Effects in cluster-SD units (σ_c = 1). Mean over 2018–2021; years agree within one
college except where a range is given. Results are identical to the unit at ρ = 0 and
0.5, and at κ = 0.5 (`a3_headline_tables.md`).

| Design | Both-arms τ = 1 | τ = 1.5 | τ = 2 | Control-only τ = 1 | τ = 1.5 | τ = 2 |
|---|---|---|---|---|---|---|
| Simple Random Sampling | 36 | 19 | ≤ 15 | > 58 | 33–34 | ≤ 15–16 |
| High Incidence Focus | 35 | 18–19 | ≤ 15 | > 58 | 31–32 (28–34) | ≤ 15 |
| Balanced Halves | 36 | 19 | ≤ 15 | > 58 | 31–32 | 15–16 |
| Balanced Quartiles | 35 | 19 | ≤ 15 | > 58 | 32–34 | ≤ 15–16 |
| Plain Saturation | 36 | 19 | ≤ 15 | > 58 | 29–30 | 16 |
| Incidence-Guided Saturation | 36 | 19 | ≤ 15 | > 58 | 28 | 16–17 |
| 2x2 Blocking | 36–37 | 19 | ≤ 15 | > 58 | 57 to > 58 | 17 |
| Checkerboard | > 58 | 34 | 19 | > 58 | > 58 | 19–20 |
| Isolation Buffer | > 58 | > 58 | 36 | > 58 | > 58 | 36–37 |

MDE at N = 58 (size-adjusted, κ = 0.8): 0.76–0.78 cluster SD for every non-geographic
design under both-arms spillover (Checkerboard 1.30, Isolation Buffer 1.52); under
control-only 1.28–1.29 for the two saturation designs, 1.49–1.59 for SRS and the
rank designs, 1.64 for 2x2 Blocking, > 2 for Checkerboard. At N = 29: 1.13–1.15
(both-arms) and 1.47–1.56 (control-only).

**Reading.** Under both-arms spillover the design barely matters for sample size:
about 36 colleges give 80% power for a one-cluster-SD effect with any
non-geographic design. Checkerboard and Isolation Buffer cost 15–25+ colleges.
Under control-only spillover the effect is harder to detect for every design. Power
almost stops rising beyond about 40 participants for SRS and the rank designs. When
nearly everyone participates, controls' treated-neighbor exposures become uniformly
high, so the Spill column approaches collinearity with (1 − Z). Non-participants
(and saturation designs) restore exposure variation. This is why the saturation
designs save 3–6 colleges at τ = 1.5 and have the lowest control-only MDE at N = 58.

### Translation to trained professionals (m) and ICC
`a3_translation_headline.md` (SRS and Incidence-Guided Saturation, ρ = 0.5, κ = 0.8)
and `a3_translation_m_icc.csv` (all designs/settings). σ_c = √(ICC + (1 − ICC)/m),
τ = d/σ_c; exact by scale invariance under β = 0 (V3). Examples, size-adjusted power,
SRS, both-arms: d = 0.2 with m = 25 and ICC = 0.05 gives 0.39 at N = 29 and 0.69 at
N = 58. d = 0.3 gives 0.70 and 0.95. With ICC = 0.01 and m = 50, d = 0.2 reaches 0.80
at N = 29. Control-only, same m and ICC: d = 0.3 gives 0.46 (N = 29) and 0.50 (N = 58),
and Incidence-Guided Saturation gives 0.49 and 0.60. Only 2 of 36 (m, ICC, d) cells
exceed the simulated τ ≤ 2 (d = 0.5, m = 50, ICC = 0.01/0.05); they are reported as
"≥ power at τ = 2".

### Does the regime guidance hold under partial participation?
MSE ratio to SRS at τ = 1, mean over 16 blocks (year × ρ × κ), random volunteering
(full table with ranges: `a3_headline_tables.md`):

| Design | Both N = 20 | 29 | 40 | 58 | Control-only N = 20 | 29 | 40 | 58 |
|---|---|---|---|---|---|---|---|---|
| High Incidence Focus | 0.98 | 0.97 | 0.96 | 0.92 | 0.99 | 1.00 | 0.99 | 1.06 |
| Balanced Halves | 0.99 | 0.98 | 0.98 | 1.00 | 0.99 | 0.97 | 1.01 | 1.01 |
| Balanced Quartiles | 1.00 | 0.99 | 0.99 | 1.01 | 1.00 | 0.99 | 0.99 | 1.02 |
| Plain Saturation | 1.00 | 0.98 | 1.00 | 1.00 | 0.99 | 0.96 | 0.90 | 0.78 |
| Incidence-Guided Saturation | 0.99 | 0.97 | 1.00 | 1.00 | 0.99 | 0.95 | 0.91 | 0.78 |
| 2x2 Blocking | 1.03 | 1.01 | 1.02 | 1.03 | 1.08 | 1.11 | 1.11 | 1.18 |
| Checkerboard | 1.88 | 2.18 | 2.23 | 2.82 | 1.18 | 1.38 | 1.80 | 3.10 |
| Isolation Buffer | 5.92 | 4.21 | 3.74 | 3.87 | 3.61 | 2.25 | 1.51 | 0.99 |

- **What holds at every N:** Checkerboard and Isolation Buffer are worse than SRS under
  both-arms spillover. 2x2 Blocking is mildly worse under control-only. The fixed-budget
  rank designs are within a few percent of SRS.
- **What weakens:** the saturation designs' control-only advantage (22% lower MSE at
  N = 58) shrinks to 9–10% at N = 40, 4–5% at N = 29 and about 1% at N = 20. With few
  participants, non-participating neighbors already create exposure variation for
  every design. Under incidence-weighted volunteering (N = 20, 29) the saturation
  advantage disappears (0.98–1.01).
- **What changes:** Isolation Buffer under control-only goes from SRS-equivalent at
  N = 58 to 2–4 times SRS's MSE at N ≤ 29. At small N it treats more than half the
  participants (islands are always treated).

A2 per N (Friedman/Nemenyi/Holm-Wilcoxon on the 16 blocks): `a2_tests_per_N_ranks.csv`,
`a2_tests_per_N_pairwise.csv`, `a2_cd_plots_per_N.pdf`. Same consistency-summary caveat.

**Saturation regions at small N (V6).** Restricting to subsets whose four regions all
have ≥ 2 colleges keeps 49/50 subsets at N = 15 and 50/50 at N ≥ 20. It changes the
saturation designs' MSE ratio to SRS by at most 0.005
(`a3_saturation_region_restriction.csv`). The V6 failure has no material effect on
these results.

**β = 1 (hypothetical sensitivity; N = 29, 58).** MSE barely changes for every design
except High Incidence Focus, whose MSE rises about 4× (0.069 → 0.293 at N = 58
both-arms; 0.151 → 0.582 at N = 29). Its treatment is nearly collinear with burden X,
as in production (`a3_beta1_vs_beta0.csv`).

**A1 at participation levels** (`a1_burden_share_per_N.csv`, random volunteering, share
of the participants' deaths): High Incidence Focus treats 33–35% of the participants'
deaths with enrichment 1.50–1.55 at every N. Incidence-Guided Saturation's enrichment
is 1.04 (N = 20), 1.02 (N = 29) and 1.00 (N = 58); the other designs are 0.96–1.02,
except Checkerboard at N = 58 (0.90).

**Incidence-weighted participation** (N = 20, 29): see the guidance table above; power
and MSE per setting in `a3_settings.csv` (Mechanism = incidence3x).

### A4. Unequal trained professionals (m̄ = 25, ICC = 0.05)
Both arms on the same noise (`a4_equal_vs_unequal.csv`, `a4_summary.csv`). With
m_j ∝ population, power falls by about 0.08–0.11 on average (at d = 0.3: 0.59 → 0.46 at
N = 20, 0.88 → 0.77 at N = 58, both-arms). MSE rises by a factor of 1.46–1.63 for every
design, because mean(1/m_j) > 1/m̄ (small colleges are noisy). No design gains: the
MSE ratio is lowest for High Incidence Focus (1.46) and Isolation Buffer (1.51) and
highest for 2x2 Blocking (1.63) and Checkerboard (1.62). The design ordering is
essentially unchanged; unequal participation sizes are a power issue, not a design
issue.

---

## 4. Sensitivity results (N = 58)

### A5 Stage 2. Coverage calibration (144 settings × 500 fits, B = 499)
`a5/fits.rds` (gitignored), `sensitivity_summaries/a5_stage2_coverage_{settings,by_design_regime}.csv`.
Four years × ρ ∈ {0, 0.5} × γ = 0.8 × both regimes × nine designs; 72,000 original
fits, about 36 million bootstrap refits. No failed SAR fits, no failed bootstrap refits,
every profile interval found; the profile log-likelihood at τ̂ equals the full
log-likelihood to 1e-13.

| Interval | Overall coverage | Both-arms | Control-only | Design × regime range | Mean width | Rule |
|---|---|---|---|---|---|---|
| Wald (z) | 92.9% | 93.1% | 92.6% | 91.2–94.3% | 1.66 | not met |
| t₅₄ and df-corrected SE (Stage 1 c) | 94.5% | 94.7% | 94.2% | 92.8–95.8% | 1.75* | not met |
| Parametric percentile bootstrap | 92.3% | 92.7% | 91.8% | 90.2–93.6% | 1.69 | not met |
| Parametric bootstrap-t | **94.7%** | 94.6% | 94.7% | 93.2–95.8% | 1.78 | not met |
| Profile likelihood (LR) | 93.5% | 93.7% | 93.3% | 91.8–94.7% | 1.70 | not met |

\* Wald width × t₀.₉₇₅,₅₄/z₀.₉₇₅ × √(58/54) = 1.66 × 1.054. Per-setting MC SE of
coverage is about 1.0–1.2% (500 fits); a design × regime cell (8 settings) is about
0.4%.

**Result.** No method meets the adoption rule at Stage 2 (94–96% for every design and
regime). Bootstrap-t and the cheap t₅₄ + df-corrected Wald interval come closest. Both
are near nominal overall, with their lowest cell Balanced Halves under control-only
(93.2% and 92.8%). The percentile bootstrap and profile likelihood do not fix the
under-coverage. That points to a small-sample variance problem (SE too small),
consistent with Stage 1, rather than skewness of τ̂. On the larger production fit set
(Stage 1, 864 settings, 7.6M fits) the t₅₄ + df-corrected interval did meet the rule at
the design × regime level (94.1–95.0%).
**Author decision needed:** (i) keep Wald and write the limitation (the rule's
fallback), or (ii) adopt the t₅₄ + df-corrected interval on the strength of Stage 1,
noting that the smaller Stage 2 check put two cells just under 94%. Design rankings
are unchanged by any interval method (τ̂ is unchanged).

### A6. Estimation-method sensitivity (CAR, reverse misspecification)
`a6car/fits.rds` (gitignored), `sensitivity_summaries/a6_car_{settings,ratio_to_srs,ratio_to_srs_by_rho,spearman}.csv`.
432 settings × 500 fits for each DGP (SAR-lag production DGP and a CAR DGP), each fit
by SAR lag and by CAR: 432,000 data sets, 864,000 fits. SAR fits: no failures. CAR
fits: 35 errors (0.008%) and 1,167 warnings (0.27%), counted. No unidentified fits.

MSE ratio to SRS (mean over 24 blocks: year × ρ × γ):

| Design | SAR DGP, SAR fit: both / ctrl | SAR DGP, CAR fit: both / ctrl | CAR DGP, SAR fit: both / ctrl | CAR DGP, CAR fit: both / ctrl |
|---|---|---|---|---|
| High Incidence Focus | 1.06 / 1.08 | 1.03 / 1.03 | 1.07 / 1.10 | 1.07 / 1.09 |
| Balanced Halves | 1.04 / 1.02 | 1.03 / 0.98 | 1.05 / 1.03 | 1.04 / 0.99 |
| Balanced Quartiles | 1.04 / 1.03 | 1.05 / 1.00 | 1.04 / 1.04 | 1.06 / 1.02 |
| Plain Saturation | 1.04 / 0.79 | 1.02 / 0.82 | 1.04 / 0.81 | 1.03 / 0.78 |
| Incidence-Guided Saturation | 1.03 / 0.81 | 1.04 / 0.82 | 1.04 / 0.81 | 1.03 / 0.76 |
| 2x2 Blocking | 0.98 / 1.21 | 0.99 / 1.17 | 1.00 / 1.22 | 1.02 / 1.22 |
| Checkerboard | 2.89 / 2.99 | 2.83 / 2.64 | 2.93 / 3.01 | 3.09 / 3.00 |
| Isolation Buffer | 4.11 / 1.03 | 4.13 / 0.95 | 4.16 / 1.03 | 4.45 / 1.05 |

**The regime guidance is robust to the estimation method and to reverse
misspecification.** Under every DGP and fitter, the saturation designs are about 20%
better than SRS under control-only spillover, Checkerboard and Isolation Buffer are 3–4
times worse under both-arms, and the rank designs are within 10% of SRS. Spearman
correlation of the nine-design MSE ordering, CAR vs SAR on the same data: 0.83 (range
0.52–0.98) both-arms and 0.90 (0.53–1.00) control-only under the SAR DGP; 0.89 and 0.97
under the CAR DGP. The best single design agrees in 11–17 of 24 blocks. Low correlations
come from reshuffling among near-tied designs (ratios within a few percent of 1), not
from changes at the ends.
Bias/coverage: under the SAR DGP at ρ = 0.5 the misspecified CAR fit is biased +0.16
with 90.1% coverage (SAR fit +0.04, 93.2%). Under the CAR DGP both fits are nearly
unbiased (|bias| ≤ 0.05) with 92.6–93.2% coverage, so the SAR fit is robust to a
CAR-type error structure.
Existing related evidence to cite alongside: the grid non-oracle and omitted-spillover
fits.

### A6. Bayesian CAR (CARBayes `S.CARleroux`, Gaussian)
`a6bayes/fits.rds` (gitignored), `sensitivity_summaries/a6_bayes_{settings,ratio_to_srs,spearman}.csv`,
`v10_hand_checked_fits.csv`. 144 settings (four years × ρ ∈ {0, 0.5} × γ = 0.8 × both
regimes × nine designs) × 100 fits = 14,400 MCMC fits, each 5,000 burn-in + 20,000
iterations thinned by 5 (4,000 kept). The SAR ML fit is on the same data. No MCMC errors.
Convergence was checked on every fit: ESS(τ) ≥ 400 and |Geweke z| < 1.96.
**4,378 fits (30.4%) fail at least one** (ESS < 400: 3,637; |Geweke| ≥ 1.96: 1,477,
against about 5% expected by chance; median ESS 3,044, min 12). Failures concentrate at
ρ = 0.5, where 42–50% of chains converge (92–94% at ρ = 0), evenly across designs
(67–72% converged per design). All fits are kept in the tables below.

| | ρ = 0: Bayes / SAR ML | ρ = 0.5: Bayes / SAR ML |
|---|---|---|
| Coverage, both-arms | 94.8% / 93.8% | 90.1% / 93.5% |
| Coverage, control-only | 94.9% / 93.4% | 87.9% / 93.2% |
| Bias, both-arms | −0.001 / +0.008 | +0.137 / +0.023 |
| Bias, control-only | −0.016 / +0.005 | +0.244 / +0.028 |

MSE ratio to SRS, Bayes (SAR ML on the same data): both-arms: High Incidence Focus 0.93
(0.97), Balanced Halves 0.96 (0.96), Balanced Quartiles 1.01 (0.99), Plain Saturation 1.09
(1.10), Incidence-Guided Saturation 1.08 (1.10), 2x2 Blocking 1.08 (1.12), Checkerboard
3.07 (3.11), Isolation Buffer 3.78 (3.85). Control-only: Incidence-Guided Saturation 0.85
(0.79), Plain Saturation 0.90 (0.84), Isolation Buffer 0.95 (1.02), Balanced Halves 1.03
(1.06), Balanced Quartiles 1.10 (1.09), High Incidence Focus 1.15 (1.21), 2x2 Blocking 1.37
(1.49), Checkerboard 2.90 (3.22). Spearman of the design ordering vs SAR ML: 0.90
(both-arms; 0.40–1.00) and 0.83 (control-only; 0.62–1.00); the best design agrees in
7/8 and 5/8 blocks. Converged fits only: the ordering is the same at the ends
(Checkerboard and Isolation Buffer worst under both-arms; Incidence-Guided
Saturation best under control-only, 0.95). Coverage stays 89% (both) and 84%
(control-only) at ρ = 0.5.
**Reading.** At ρ = 0 the Bayesian CAR agrees with ML. Under SAR-lag dependence at
ρ = 0.5 it is misspecified: biased, under-covering, and half its chains do not mix at
this length. The design guidance is still similar at the ends. Bayesian CAR is
therefore a robustness check of the design ordering, not a better trial analysis for
this DGP.


### A7. Spillover-form sensitivity (binary, additive count; proportional as reference)
`a7/fits.rds` (gitignored), `sensitivity_summaries/a7_{forms_settings,identification,ratio_to_srs}.csv`.
288 settings per form, J = 100 × R = 100 (deterministic designs R = 1,000): 6,912,000 fits.
DGP and fit use the same form. The proportional form reproduces production (V11 at
N = 58): 9 of 288 settings with |z| > 2 (13 expected), max |z| 2.46.

MSE ratio to SRS (mean over 16 blocks):

| Design | Proportional both / ctrl | Count both / ctrl | Binary both / ctrl |
|---|---|---|---|
| High Incidence Focus | 0.91 / 1.03 | 0.92 / 1.04 | 0.96 / n.i. (50%) |
| Balanced Halves | 1.02 / 0.98 | 1.02 / 0.99 | 1.02 / 0.95 (59%) |
| Balanced Quartiles | 1.01 / 0.99 | 1.01 / 0.99 | 1.00 / 0.96 (60%) |
| Plain Saturation | 1.01 / 0.76 | 1.02 / 0.81 | 1.03 / 0.70 (88%) |
| Incidence-Guided Saturation | 1.01 / 0.75 | 1.03 / 0.79 | 1.04 / 0.71 (86%) |
| 2x2 Blocking | 1.04 / 1.17 | 1.03 / 1.07 | 1.00 / 1.14 (28%) |
| Checkerboard | 2.97 / 3.09 | 1.78 / 1.41 | 1.07 / **not identified** |
| Isolation Buffer | 3.82 / 0.95 | 3.84 / 1.16 | **not identified** / **not identified** |

(Percentages: minimum share of identified fits per block; the SRS reference is itself
identified in only 66–70% of fits under binary control-only spillover, so those ratios
use identified fits only. "n.i.": identified in about half the fits.)

**Result.** The regime guidance holds under the count form: saturation designs best
under control-only (0.79–0.81), Checkerboard and Isolation Buffer worst under
both-arms. Checkerboard's penalty is smaller (1.78 / 1.41). Under the binary
(Project 1-style) form the design question changes. Isolation Buffer is never
identified (with a maximal independent set every control has a treated neighbor, so
Spill = γ(1 − Z)), and Checkerboard is not identified under control-only, as the plan
expected. Under control-only every design is partly unidentified (in 30–54% of SRS,
rank-design and 2×2 allocations every control has a treated neighbor). The saturation
designs are identified most often (86–90%) and are best when identified (0.70–0.71).
Under binary both-arms spillover every identified design is within 0–7% of SRS
(Checkerboard 1.07).
These designs are reported as not identified, never dropped.



---

## 5. Verification checks

| Check | Result |
|---|---|
| V1 reproduction at N = 58 | PASS: 288 cells (year × design × regime × ρ × κ), 576 MSE + coverage z-scores vs production; 25 with \|z\| > 2 (26.2 expected by chance; binomial p = 0.62), none > 3, max 2.96 (`summaries/v1_reproduction.csv`) |
| V2 size at τ = 0 | **FAIL for N ≤ 35, PASS for N ≥ 40**: design-mean Wald size 6.3–7.0% at N = 58, 6.9–7.7% at N = 40–50, 7.6–8.8% at N = 29–35, 9.5–12.8% at N = 15–20 (settings up to 13.7%). Size-adjusted rejection is 5% by construction (`summaries/v2_size.csv`) |
| V3 scale invariance | PASS: (0.5, 0.5) vs (1, 1), 200 fits, max |Δ| in τ̂, SE and t = 5.6e-8 (the ρ optimizer's tolerance, not exact floating-point equality) |
| V4 rule A toy map | PASS: exposure 0.25; island 0 with a zero W_P row |
| V5 islands in the fitter | PASS with a recorded deviation (below): 1,000 fits vs `lagsarlm(zero.policy = TRUE)`, max |Δτ̂| 3.7e-8, max |ΔSE| 5.2e-9; 300 symmetric-W_P fits vs lagsarlm's symmetric branch, max |Δτ̂| 1.6e-8 |
| V6 designs on subsets | **FAIL** (one part): see below |
| V7 informative sampling | PASS: highest/lowest inclusion ratio 3.05 (N = 20), 2.99 (N = 29); top/bottom decile 2.67 vs design value 2.68 |
| V8 monotonicity | **FAIL as specified (flagged, not smoothed)**: 1,842 of 12,096 curves drop by > 2 MC SEs. None along τ (common random numbers). All are along N, concentrated at 15→20 (Wald only: size distortion) and at 29→30, 45→50, 50→58. Every design and year shares the same 50 subsets at a given N, so one unlucky subset set flags hundreds of curves together. 50→58 also switches from rebuilt to frozen partitions. Mean drop 0.034, max 0.17 (`summaries/v8_monotonicity_violations.csv`). N₈₀ bootstrap intervals over subsets carry this uncertainty |
| V9 manifests | PASS (03:28 and 12:39, after all runs): the SHA-256 hashes of all 10 NC production computation files and the grid manifest hashes of code/01–05 are unchanged (`verification/v9_manifests.txt`) |
| V10 A6 sanity | PASS with a caveat: CAR vs SAR at ρ = 0 on the same data, correlation 0.985, mean \|Δτ̂\| 0.052 (about 0.1 SE), mean Δ −0.015 (71,991 fits). Bayes vs ML: correlation 0.950, mean \|Δ\| 0.12, mean Δ +0.075 (driven by ρ = 0.5); 5 hand-checked fits in `v10_hand_checked_fits.csv` (posterior means within 0.4 of ML, most within 0.1). Non-converged MCMC fits are counted: 4,378 / 14,400 (Section 4) |
| V11 A7 exposures | PASS: toy binary/count/proportional match hand values; N = 58 proportional exposure equals production W z; the N = 58 DGP mean and fitted Spill equal production `rs_model()` (1e-12) |

Other passes: every N = 58 design draw equals production's (same call, same seed) for
all nine designs; the N = 58 participant W equals production's queen W.

**V5 deviation (engineering, not design).** `sar_lag_setup()` refuses listw objects for
which `spatialreg::can.be.simmed()` is TRUE. A participant W_P is itself symmetric
when every non-island participant has the same within-subset degree (e.g. all
components are pairs): 76 of 2,000 random N = 15 subsets. The plan's V5 asks for this
guard to pass; it does not for these subsets. Rather than drop them, `me_engine_setup()`
marks a copy of the listw as non-symmetric and calls the unchanged `sar_lag_setup()`
(similar matrices share eigenvalues, so the likelihood is identical). This was
validated against `lagsarlm()`, which takes its own symmetric branch on the original
listw. Counted per subset (`Guard_Bypass`).

**V6 failure.** Fixed-budget designs treat ⌊N/2⌋; Isolation Buffer is independent,
maximal and treats islands; the saturation/blocking rules hold; the checkerboard is
deterministic. **But the rebuilt k-means saturation regions (`make_kmeans_regions()`
on the participants, as approved) are not always ≥ 2 colleges:** on 1,000 random
subsets, a region of 1 college occurred in 8% of N = 15 subsets and 0.9% at N = 20,
and none at N ≥ 25; a region below 3 occurred in 38% (N = 15), 7% (N = 20), about 1%
(N = 25–30), 0 above. The function scores partitions by balance and compactness but
does not enforce a minimum size. Decision taken unattended: run A3 with all nine
designs (the runner itself is unaffected, and every subset's minimum region size is
stored), and report Plain and Incidence-Guided Saturation at N < 25 with this caveat
plus a restriction to subsets whose regions all have ≥ 2 colleges
(`summaries/a3_saturation_region_restriction.csv`). **Author decision needed:** keep,
restrict, or rebuild with a minimum-size rule (a rerun of designs 3 and 8 only, about
2/9 of the A3 compute).

---

## 6. Teach-back: what was built and why

All new code is in new files; the hashed production files are sourced, never edited
(V9). Every new script has the standard header and roxygen blocks.

### `application/code/methods_ext_core.R` — shared helpers
**Purpose.** One place for everything the participation, power and sensitivity
analyses share, so N = 58 runs reproduce production exactly and N < 58 runs differ
only where the approved plan says.

**Logic.**
- `me_seed(...)` is production's `allocation_seed()` with a new namespace
  (`"nc-methods-ext-v1"`). Key fields never include τ, κ or d, so all effects share
  allocations and noise (common random numbers, smooth power curves).
- `me_load_setup()` reads the frozen `setup.rds` and stops unless its design digest
  equals the production manifest's, so N = 58 reuses exactly the production regions,
  yearly blocks and ranks.
- `me_subset()`: random volunteering = SRS of N colleges keyed by (N, s) only, so the
  same 50 subsets serve every year and design. Incidence-weighted = `me_upsystematic()`
  on πⱼ = N wⱼ / Σw with wⱼ = 1 + 2(rⱼ − 1)/57.
- `me_upsystematic(pi)`: randomized systematic πps sampling. Units are shuffled, then a
  uniform start u is stepped through the cumulative πs; unit k is selected when an
  integer + u falls in its interval, so P(k selected) = πₖ exactly and the sample size
  is exactly N.
- `me_geometry(P)`: the participant adjacency A_P, W_P (row-standardized; zero rows for
  islands), the plain neighbor list for designs, island and component counts, the
  listw (`zero.policy = TRUE`) and the lean-engine setup via `me_engine_setup()`.
- `me_exposure(z, g)`: rule A, eⱼ = (A[P,P] z)ⱼ / deg_full(j). Dividing by the full
  geographic degree is what makes non-participants count as untreated neighbors.
- `me_exposure_form()`: proportional (= rule A), binary 1{any treated neighbor},
  count (treated neighbors)/d̄ with d̄ = 4.931 (A7).
- `me_design_inputs()`: at N = 58 the frozen production objects; at N < 58, ranks among
  participants, X = (r − 0.5)/N, k-means regions on participants (pooled population, as
  production), 2×2-blocking analogue blocks (year's incidence and population), and the
  regional mean participant rank for Incidence-Guided Saturation.
- `me_draw_design()`: designs 1 and 3–8 call the unchanged `get_application_designs()`;
  High Incidence Focus and SRS treat ⌊N/2⌋ at N < 58 (identical to production at 58).
- `me_fit()`: `fit_one_lag_model()` (lean engine) plus classification. Boundary-ρ and
  error fits are invalid (as in production); aliased fits are kept and counted.
- `me_manifest()`: config hash, SHA-256 of the new and sourced production files, setup
  digest, package versions, BLAS. A mismatch stops the run (no stale checkpoint reuse).

**Design choices.** Frozen N = 58 objects make V1 an exact-code reproduction check.
Systematic πps was chosen over the plan's successive-draw scheme (builder's call) because
successive draws compress inclusion ratios as N grows: at N = 29 they would give about
2.2 rather than 3. πps gives an exact 3× ratio between the extreme colleges (V7).
**Gotchas.** `make_kmeans_regions()` calls `set.seed()` internally; it is seeded with
the subset key's integer, so it does not disturb later draws. Regions can be tiny at
small N (V6). Symmetric W_P needs `me_engine_setup()` (V5 deviation).
**R notes.** `spdep::mat2listw(..., zero.policy = TRUE)` warns about sub-graphs; the
warning text is captured per subset, not suppressed silently. `spatialreg::eigenw()` and
`can.be.simmed()` are the functions behind the guard.

### `application/code/methods_ext_participation.R` — A3/A4 runner
**Purpose.** Simulate and fit every A3 (partial participation + power) and A4 (unequal
professionals) setting, checkpointed per (mechanism, year, N, design) task. That is
finer than the plan's (year, N, mechanism), so it is a superset of the planned
checkpointing.
**Logic** (`mp_task`). For each subset s: geometry, design inputs, deterministic flag
(J = 1, R = 200, or R = 1,000 at N = 58), (I − ρW_P)⁻¹ for each ρ. For each allocation j:
seed by (key, year, design, j), draw z, compute e(z), seed the noise, and draw one
n × R standard-normal matrix shared by every cell. For each regime, β, ρ and arm:
Y = A(τz + κτ·s_reg + βX + σ∘ε). That is τAz + κτAs + βAX plus A(σ∘ε), each computed
once per cell family. Fits use intercept, Z and unscaled Spill (and X when β = 1).
Per (cell, subset, allocation) the runner keeps counts and sums (valid, rejections,
coverage, Σe, Σe², alias, boundary, error, warning, unidentified). |t| at τ = 0 is kept
per fit; the size-adjusted critical value is its 95th percentile within the same
(β, ρ, regime, arm).
**I/O.** In: `mp_config(profile)`, the frozen setup. Out: one RDS per task with
`cells`, `stats[cell, subset, allocation, stat]`, `t0`, allocation diagnostics
(treated count, death and person-year shares, treated islands) and subset diagnostics
(islands, components, guard bypass, region and block sizes).
**Choices.** τ = 0 is fitted once and reported under both κ (the DGP is identical). A4
simulates both the equal-m (σ_c) and unequal-m arms on the same noise (builder's call),
so "does any design gain or lose" is a paired comparison rather than an interpolation.
β = 1 uses each college's own rank-scaled burden a$X[P] (not participant ranks): the
outcome tracks the college's burden regardless of who joins.
**Gotchas.** Changing either runner file changes the manifest, so the checkpoints
would be refused and the run would restart. `mc.preschedule = FALSE` balances N = 58
and small-N tasks. Worker errors are logged per task and listed at the end. A subset
with no participant edges would be recorded as fit errors (none expected at N ≥ 15).

### `application/code/methods_ext_stored.R` — A0, A1, A2, A5 Stage 1
- `ms_a0_a5()`: reads each production checkpoint and computes per-allocation means of
  |t| > z, |t| > t₅₄ and the four interval corrections, weighted by allocation draw
  frequency (as production's `rs_distribution()`). It first checks that the recomputed
  Wald coverage equals the stored one (max difference 6e-16).
- `ms_a1()`: production draws for each year and design are shared across ρ, γ and
  regime (support key), and this is asserted. Shares are computed per allocation.
- `me_design_tests()`: Friedman (`stats::friedman.test`), Nemenyi
  (`PMCMRplus::frdAllPairsNemenyiTest`, long format as in code/10), CD =
  q₀.₉₅,ₖ,∞/√2 · √(k(k+1)/6N), and Holm-adjusted paired Wilcoxon for all 36 pairs.
  `me_cd_plot()` evaluates code/10's `plot_cd_diagram()` from its parsed source,
  without sourcing the whole grid visualization stack.
**Gotcha.** A2 uses `Mean_MSE` per setting, so with 24 blocks every Wilcoxon p is
bounded below by the signed-rank distribution. With so many Monte Carlo replicates the
tests are a consistency summary, as the plan says.

### `application/code/methods_ext_sensitivity.R` — A5 Stage 2, A6, A7
- `mx_profile_ll(τ)`: the concentrated SAR-lag likelihood at fixed τ. It uses
  M(y − τz) and MWy (the lag term is not offset), maximized over ρ with the same
  log-determinant as `fit_sar_lag()`. Checked: at τ̂ it equals the full log-likelihood
  (≤ 6e-14 in the smoke run). `mx_profile_ci()` finds the two roots of
  2[ℓ(τ̂) − ℓ_p(τ)] = χ²₁,₀.₉₅ by `uniroot`.
- `mx_boot()`: parametric bootstrap from the fitted model,
  y* = (I − ρ̂W)⁻¹(Xb̂ + ŝε*), B = 499, refit by the lean engine; percentile and
  bootstrap-t intervals; failed refits are counted.
- `mx_fit_car()`: `spautolm(family = "CAR")` with `similar.listw()` weights,
  S = D^(−1/2) A D^(−1/2). S is symmetric with the same eigenvalues as the
  row-standardized W, so λ lives on ρ's scale and spautolm's non-symmetric-weights
  warning (and its NaNs) are avoided.
  Reverse-misspecification DGP: u ~ N(0, (I − λS)⁻¹), i.e. iid N(0, 1) at λ = 0, with
  the same mean τz + Spill.
- `mx_fit_bcar()`: `CARBayes::S.CARleroux` (Gaussian, binary W). It runs 5,000 burn-in
  plus 20,000 iterations thinned by 5 (4,000 kept); the builder's call is longer than
  the plan's 2,000 + 5,000 start because ESS for τ was about 40 at that length. The
  posterior mean and 2.5/97.5% quantiles of the Z coefficient are kept, with
  `coda::effectiveSize` and `coda::geweke.diag` for every fit (a superset of the 5%
  check).
**Gotcha.** A7 binary spillover is not identified for Isolation Buffer (both regimes)
and Checkerboard (control-only), and is partly unidentified for SRS under
control-only. Those fits are kept with `Identified = FALSE`, and MSE is computed over
identified fits with the identified share reported.

### Summaries and checks
- `methods_ext_summaries.R` reduces each checkpoint to setting rows with MC SEs, with
  the subset as top-level unit. It fills the shared τ = 0 rows under both κ and then
  computes:
  - V1 (z-scores vs production with a binomial test of the |z| > 2 count) and V2;
  - V8 (drops beyond 2 MC SEs are flagged, not smoothed);
  - N₈₀ by linear interpolation at the first crossing, with censoring and a subset
    bootstrap (B = 500);
  - MDE(N), and the translation tables (τ_c = d/σ_c, power interpolated over the τ
    grid, flagged beyond τ = 2);
  - the guidance table with per-N A2 tests (16 blocks), the saturation-region
    restriction, A1 per N, β = 1 vs β = 0, and A4 paired comparisons.
- `methods_ext_sensitivity_summary.R`: setting metrics with allocation-level MC SEs,
  the A5 adoption rule per method, A6 ratios to SRS and Spearman agreement of design
  orderings, V10, A7 ratios with identification shares, and an A7 proportional vs
  production check.
- `application/tests/test_methods_extensions.R` (V3–V7, V11, plus N = 58 equivalence)
  and `verify_methods_ext_manifests.R` (V9).
- `application/code/run_methods_ext_queue.sh`: the unattended step-4/5 queue (each step
  logged, failures logged and the queue continues).

