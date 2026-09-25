# A5 technical notes: investigation of four anticipated committee questions

**DRAFT for the coordinator and author. Nothing here has been applied to the chapter.**

Written 2026-09-25. Method authority: `docs/plans/simulation-revision-spec.md`. Data: the full run
`results/sim_data/*_20260924_025509.rds`. Nothing under `results/` was written.

## How the evidence was produced

- **Stored results.** Scenario-level numbers come from
  `results/sim_data/sim_results_MLE_tau_sweep_combined_20260924_025509.rds` (oracle; "stored").
- **Per-fit re-simulation.** The stored files keep neither ρ̂ nor the per-fit SE. We re-ran every
  block of the full run for the six chapter designs, both neighbor types and all five τ (9,600
  scenarios, 2.4 million oracle fits). The re-run sourced `code/05_run_simulation.R` with
  `sim_define_only <- TRUE` and used the same `set_seed_key()` keys, `draw_assignments()` and
  `fit_sar_lag()`. For each fit it recorded τ̂, SE(τ̂), ρ̂ and the decomposition terms below. It
  **reproduces the stored results exactly**: over all 9,600 scenarios, max |ΔBias| = max |ΔMSE| =
  max |ΔSD| = max |ΔCoverage| = 0. Scripts: `resim_a5.R`, `analyze_a5.R` and `analyze_a5_b.R` in the
  session scratchpad
  (`/private/tmp/claude-501/-Users-ajwalther-GithubProjects-SpatialCRT-projects-IncidenceDesign/41fd41d4-9dd6-45d6-89b4-356afa38517d/scratchpad/`).
  Everything labeled "re-sim" below comes from these fits.
- **Null calibration (Q4 only).** For Q4 we added 400 independent noise sets for the 80 rook
  Checkerboard blocks. They use a new seed-key namespace, `("eps_a5rep", mode, ρX, "rook", ρ, γ, r)`,
  so they cannot collide with the run's `"eps"` keys. Scripts: `q4_null_reps.R` and
  `analyze_q4_reps.R` (same scratchpad).

### The error decomposition used in Q1 and Q2 (exact)

Notation follows the chapter. The oracle design matrix is x = [1, Z, S(Z), X_k], with true coefficients θ = (0, τ, 1, β).
Also A = (I − ρW)⁻¹, P = (x′x)⁻¹x′, and M = I − xP. Given ρ̂, the ML coefficients are b̂ = P(I − ρ̂W)Y.
Since (I − ρW)Y = xθ + ε, we have (I − ρ̂W)Y = xθ + ε + (ρ − ρ̂)WY. Writing WY = WA(xθ + ε) gives

  τ̂ − τ = [Pε]_Z + (ρ − ρ̂)·h + **τ·(ρ − ρ̂)·g**,  g = [P·WAZ]_Z,  h = [P·WA(S + βX + ε)]_Z.

Here g is the coefficient on Z when the lagged treatment WAZ = WZ + ρW²Z + ρ²W³Z + … is regressed on
the model's own covariates. [Pε]_Z has mean zero given (Z, X). So a τ-proportional error term exists
exactly when g ≠ 0, and its sign and size are set by the error in ρ̂.

**When ρ̂ does not depend on τ.** The concentrated likelihood depends on Y only through
SSE(ρ̃) = ‖M(I − ρ̃W)Y‖². The τ-dependent part of Y is τAZ, and (I − ρ̃W)A = I + (ρ − ρ̃)WA gives
M(I − ρ̃W)AZ = (ρ − ρ̃)·M·WAZ. If WAZ lies in the column space of x, this is zero for every ρ̃.
Then ρ̂ is the same at every τ. If also g = 0, then τ̂ − τ is exactly the same at every τ.

The identity was checked on every fit. On the 2,199,980 non-aliased fits the largest residual was
4.2 × 10⁻¹³. On the 200,000 rook Checkerboard fits, where the coefficient on Z is τ − γ, it was
2.1 × 10⁻¹⁵ (re-sim).

---

## Q1. Why does the oracle bias of Incidence-Guided Saturation Quadrants grow with τ?

### Proposed appendix paragraph

> The oracle estimate of $\tau$ from Incidence-Guided Saturation Quadrants had a small positive bias that grew with $\tau$, from 0.029 at $\tau = 0.8$ to 0.047 at $\tau = 3$ under queen contiguity, whereas the bias of Isolation Buffer stayed between 0.002 and 0.004. The cause is the estimate of the spatial lag coefficient. Given $\hat\rho$, the coefficients are least squares estimates from $(I - \hat\rho W)Y$, so the error in $\hat\tau$ contains the term $\tau(\rho - \hat\rho)g$. Here $g$ is the coefficient on $Z$ when the lagged treatment $W(I - \rho W)^{-1}Z$ is regressed on the model's covariates $1$, $Z$, $S(Z)$, and $X_k$. In our fits $\hat\rho$ was biased downward by about 0.06 on average, so this term is positive and proportional to $\tau$ for any design with $g > 0$. Under the both-arms regime with $\rho = 0$, the lagged treatment is $WZ = S(Z)/\gamma$, one of the covariates, so $g = 0$. Under the control-only regime, $WZ = S(Z)/\gamma + Z \circ WZ$, and $g$ is approximately the average share of treated neighbors among treated clusters. That share is 0.57 for Incidence-Guided Saturation Quadrants, whose treated clusters are concentrated in high-saturation quadrants, and exactly 0 for Isolation Buffer, whose treated clusters have no treated neighbors. For Isolation Buffer the two regimes define the same spillover covariate, $g$ is zero at $\rho = 0$ and about 0.06 on average at $\rho = 0.5$, and when $\rho = 0$ its estimation error is identical at every value of $\tau$. The same term explains the rise in the control-only MSE of Incidence-Guided Saturation Quadrants, from 0.131 at $\tau = 0.8$ to 0.182 at $\tau = 3$: the variance of $\tau(\rho - \hat\rho)g$ grows with $\tau$ (more slowly than $\tau^2$, because $\hat\rho$ becomes less variable), from 0.008 to 0.061, while the rest of the error variance changes little. The growth in bias slows at larger $\tau$ because the same term strengthens the information about $\rho$ in the likelihood, so the downward bias of $\hat\rho$ shrinks as $\tau$ grows.

(Numbers: see Evidence. The chapter can round 0.566–0.567 to 0.57. If the paragraph is used,
"0.06" should be cited as "between 0.054 and 0.071 for this design across τ" or as the pooled
"about 0.06".)

### Evidence

**Stored, queen, pooled over 160 blocks** (these match the chapter's Stability section):

| τ | Incidence-Guided SQ bias | Incidence-Guided SQ MSE | Isolation Buffer bias | Isolation Buffer MSE |
|---|---|---|---|---|
| 0.8 | 0.0286 | 0.0884 | 0.0025 | 0.1602 |
| 1.0 | 0.0321 | 0.0907 | 0.0026 | 0.1602 |
| 1.5 | 0.0390 | 0.0970 | 0.0030 | 0.1602 |
| 2.0 | 0.0435 | 0.1035 | 0.0034 | 0.1603 |
| 3.0 | 0.0467 | 0.1147 | 0.0042 | 0.1604 |

Squared bias as a share of MSE for Incidence-Guided Saturation Quadrants (stored,
mean(Bias²)/mean(MSE)): 0.016, 0.019, 0.026, 0.030 and 0.030 at the five τ. This matches the
chapter's "at or below 3.0%".

**Exact bias decomposition (re-sim, queen, pooled over 160 blocks):** bias = mean[Pε]_Z +
mean[(ρ − ρ̂)h] + mean[τ(ρ − ρ̂)g].

| Design | τ | Bias | [Pε] term | h term | **τ-term** | mean(ρ̂ − ρ) | mean g |
|---|---|---|---|---|---|---|---|
| Incidence-Guided SQ | 0.8 | 0.0286 | 0.0012 | 0.0070 | **0.0205** | −0.071 | 0.371 |
| Incidence-Guided SQ | 3.0 | 0.0467 | 0.0012 | 0.0024 | **0.0430** | −0.054 | 0.371 |
| Isolation Buffer | 0.8 | 0.0025 | −0.0004 | 0.0019 | **0.0010** | −0.054 | 0.020 |
| Isolation Buffer | 3.0 | 0.0042 | −0.0004 | 0.0010 | **0.0035** | −0.054 | 0.020 |
| Balanced Quartiles | 0.8 → 3.0 | 0.0141 → 0.0371 | −0.0009 | −0.0008 → −0.0039 | 0.0159 → 0.0419 | −0.066 → −0.055 | 0.306 |
| High Incidence Focus | 0.8 → 3.0 | 0.0287 → 0.0423 | 0.0012 | 0.0120 → 0.0061 | 0.0156 → 0.0351 | −0.066 → −0.052 | 0.305 |
| 2×2 Blocking | 0.8 → 3.0 | 0.0073 → 0.0291 | −0.0004 | −0.0031 → −0.0054 | 0.0107 → 0.0349 | −0.057 → −0.053 | 0.235 |
| Checkerboard | 0.8 → 3.0 | −0.0017 → 0.0283 | −0.0079 | −0.0031 → 0.0022 | 0.0093 → 0.0339 | −0.049 | 0.214 |

The growth in every design's bias is the τ-term. The [Pε] term does not depend on τ, and the h term
is small and shrinks. The designs order by g, and Isolation Buffer's g is an order of magnitude
smaller than the others'.

**The design-specific factor g (re-sim, queen, τ = 1; mean over fits):**

| Design | Both arms, ρ = 0 / 0.01 / 0.2 / 0.5 | Control only, ρ = 0 / 0.01 / 0.2 / 0.5 | Mean WZ among treated |
|---|---|---|---|
| Incidence-Guided SQ | ~0 / 0.002 / 0.039 / 0.128 | 0.560 / 0.565 / 0.678 / 0.997 | 0.566–0.567 |
| Balanced Quartiles | ~0 / 0.001 / 0.031 / 0.096 | 0.494 / 0.497 / 0.569 / 0.760 | 0.494–0.495 |
| High Incidence Focus | ~0 / 0.001 / 0.031 / 0.096 | 0.484 / 0.487 / 0.566 / 0.776 | 0.531 |
| Checkerboard | ~0 / −0.0001 / −0.002 / −0.005 | 0.462 / 0.460 / 0.423 / 0.374 | 0.461 |
| 2×2 Blocking | ~0 / 0.001 / 0.018 / 0.054 | 0.421 / 0.421 / 0.445 / 0.518 | 0.421–0.422 |
| Isolation Buffer | ~0 / 0.001 / 0.018 / 0.060 | ~0 / 0.001 / 0.018 / 0.061 | 0.000 |

"~0" means |g| < 10⁻¹⁴. At ρ = 0 under the control-only regime, g is close to the mean treated-neighbor
share among treated clusters, which is the algebraic reason given in the paragraph. It is not exactly
equal to that share, because the regression also adjusts for S and X. For Isolation Buffer the two
regimes give the same g, because Z ∘ WZ = 0.

**ρ̂ is τ-invariant where the algebra says it should be** (re-sim; the range of ρ̂ across the five τ
for the same fit):
- ≤ 1.04 × 10⁻⁷ (optimizer tolerance) for every design under the both-arms regime at ρ = 0, and for
  Isolation Buffer under both regimes at ρ = 0, queen and rook.
- ≤ 7.0 × 10⁻⁸ for every rook Checkerboard fit.
- Up to 0.51 for Incidence-Guided Saturation Quadrants under the control-only regime.

In stored numbers, the queen control-only bias of Isolation Buffer at ρ = 0 is −0.00777 at all five τ,
and its both-arms bias at ρ = 0 is 0.00160 at all five τ. The same holds for Incidence-Guided
Saturation Quadrants under the both-arms regime at ρ = 0 (0.01125 at all five τ).

**ρ̂ is biased downward** (re-sim, queen, τ = 1, all six designs): mean ρ̂ − ρ = −0.056, −0.058,
−0.063 and −0.063 at ρ = 0, 0.01, 0.2 and 0.5. The bias shrinks with τ where WAZ is not in the span of
x. For Incidence-Guided Saturation Quadrants under the control-only regime, mean ρ̂ − ρ goes from −0.068
(τ = 0.8) to −0.038 (τ = 3), and SD(ρ̂ − ρ) from 0.166 to 0.123. The algebraic reason: SSE(ρ̃)
contains τ²(ρ − ρ̃)²‖M·WAZ‖², which sharpens the likelihood around the true ρ as τ grows. This is why
the bias curve flattens (0.0435 at τ = 2, 0.0467 at τ = 3).

**Control-only MSE decomposition (re-sim, queen; within-scenario components averaged over 80 blocks):**

| Design | τ | MSE | Bias² | Var of τ-term | Var of rest | 2·Cov |
|---|---|---|---|---|---|---|
| Incidence-Guided SQ | 0.8 | 0.1311 | 0.0023 | 0.0082 | 0.1190 | 0.0015 |
| Incidence-Guided SQ | 3.0 | 0.1822 | 0.0062 | 0.0612 | 0.1145 | 0.0003 |
| Isolation Buffer | 0.8 | 0.1595 | 0.0006 | 0.00001 | 0.1588 | 0.00002 |
| Isolation Buffer | 3.0 | 0.1596 | 0.0007 | 0.00019 | 0.1588 | 0.00002 |

The control-only MSE values (0.131 → 0.182; 0.1595 → 0.1596) match the chapter.

### Verdict: **SUPPORTED**

The decomposition is an exact identity, verified per fit. The τ-term accounts for the growth in bias
and in control-only variance for every design. The immunity of Isolation Buffer follows algebraically
from Z ∘ WZ = 0, and is exact at ρ = 0. One component is observed rather than derived: the sign and
size of the ρ̂ error (a downward bias of about 0.06). The appendix should present it as observed in
our fits. We did not check it against a published finite-sample result.

---

## Q2. Why does rook Checkerboard coverage rise with τ?

### Proposed appendix paragraph

> Under rook contiguity, Checkerboard's coverage rose from 0.138 at $\tau = 0.8$ to 0.442 at $\tau = 3$, although its coefficient on $Z$ estimates $\tau - \gamma$ at every $\tau$. The rise comes from the standard error, which grew with $\tau$ along with the actual dispersion of $\hat\tau$. For this design $WZ = 1 - Z$, so $W(I - \rho W)^{-1}Z = \{(1 - Z) + \rho Z\}/(1 - \rho^2)$ lies in the span of the intercept and $Z$. Two consequences follow. The estimate $\hat\rho$ is the same at every value of $\tau$, and the error in $\hat\tau$ has the exact form $\hat\tau - \tau = -\gamma + c + (\tau - \gamma)(\hat\rho - \rho)/(1 + \rho)$, where $c$ does not depend on $\tau$. The variability of $\hat\rho$ is therefore multiplied by $\tau - \gamma$: the mean standard deviation of $\hat\tau$ rose from 0.202 to 0.327 over $\tau = 0.8$ to 3. The maximum likelihood standard error accounts for this propagated uncertainty in $\hat\rho$ and grew with it, from 0.199 to 0.322 on average. With the bias fixed near $-\gamma$, wider intervals covered $\tau$ more often. The increase was largest at $\rho = 0$ (coverage 0.145 to 0.541) and smallest at $\rho = 0.5$ (0.130 to 0.265), as the factor $1/(1 + \rho)$ implies. The same term raised the MSE from 0.477 to 0.567: about three quarters of the increase was variance (0.041 to 0.109), and the rest was squared bias (0.436 to 0.458). The squared bias grew because $\hat\rho$ was biased slightly downward, which moved the mean bias from $-0.650$ to $-0.667$. The higher coverage at large $\tau$ reflects less precise estimates of the wrong quantity, not a recovery of $\tau$.

### Evidence

**Algebra (exact).** Under rook contiguity the checkerboard gives WZ = 1 − Z exactly, so WᵐZ
alternates between 1 − Z (m odd) and Z (m even). Then WAZ = Σ ρᵐWᵐ⁺¹Z = [(1 − Z) + ρZ]/(1 − ρ²),
which is in span(1, Z). Hence:
- g = −1/(1 + ρ). Re-sim: max |g + 1/(1 + ρ)| = 8.9 × 10⁻¹⁶ over 200,000 fits.
- ρ̂ is τ-invariant. Re-sim: max range across τ = 7.0 × 10⁻⁸.
- Spill = γ(1 − Z) is aliased and dropped, so the reduced model's coefficient on Z is τ − γ, with
  intercept γ.

The decomposition becomes τ̂ − τ = −γ + c + (τ − γ)k, with k = (ρ̂ − ρ)/(1 + ρ) and
c = [Pε]_Z + (ρ − ρ̂)h′ free of τ. Re-sim: c varies across τ by at most 3.3 × 10⁻⁸. Therefore
Var(τ̂) = Var(c) + 2(τ − γ)Cov(c, k) + (τ − γ)²Var(k), and Bias = −γ + E[c] + (τ − γ)E[k].

**Why the SE tracks this.** In `fit_sar_lag()` (`code/04_estimation.R`) the information matrix
couples ρ with the coefficients through x′AWx·b̂ and ‖AWx·b̂‖². After inversion, the ρ-uncertainty
part of Var(τ̂) is driven by [P·WAx·b̂]_Z. With b̂ ≈ (γ, τ − γ, β), that is (τ − γ)g + β[P·WAX]_Z,
the same quantity that multiplies ρ̂ − ρ in the error decomposition. This is a reading of the
formula. It is supported empirically by the SE/SD agreement below, but we did not derive a closed
form for Var(τ̂).

**Numbers (rook Checkerboard, 160 scenarios per τ):**

| τ | Bias (stored) | mean SD (stored) | mean SE (re-sim) | Coverage (stored) | MSE (stored) | mean Bias² (stored) | mean Var (stored, SD²·249/250) |
|---|---|---|---|---|---|---|---|
| 0.8 | −0.6504 | 0.2024 | 0.1986 | 0.1379 | 0.4767 | 0.436 | 0.0409 |
| 1.0 | −0.6519 | 0.2056 | 0.2016 | 0.1498 | 0.4799 | 0.438 | 0.0422 |
| 1.5 | −0.6556 | 0.2229 | 0.2184 | 0.1979 | 0.4924 | 0.443 | 0.0497 |
| 2.0 | −0.6593 | 0.2511 | 0.2461 | 0.2693 | 0.5110 | 0.448 | 0.0633 |
| 3.0 | −0.6667 | 0.3274 | 0.3215 | 0.4419 | 0.5667 | 0.458 | 0.1087 |

In the re-sim, within-scenario Bias² + Var (denominator n) sums exactly to MSE: 0.4358 + 0.0409 =
0.4767, and 0.4580 + 0.1087 = 0.5667. The MSE increase of 0.0900 is 0.0679 variance (75%) and 0.0222
squared bias (25%). Mean k = −0.0074 (SD 0.111; re-sim, pooled). The predicted change in bias,
(3 − 0.8) × (−0.0074) = −0.0164, matches the observed −0.0163.

**By ρ (stored SD and coverage; re-sim SE), τ = 0.8 → 3:** ρ = 0: SD 0.206 → 0.369, SE 0.199 → 0.359,
coverage 0.145 → 0.541. ρ = 0.01: 0.203 → 0.360, 0.199 → 0.357, 0.138 → 0.533. ρ = 0.2:
0.202 → 0.316, 0.198 → 0.313, 0.138 → 0.430. ρ = 0.5: 0.199 → 0.265, 0.198 → 0.257, 0.130 → 0.265.

**By γ (stored), coverage at τ = 0.8 → 3:** γ = 0.5: 0.298 → 0.633. γ = 0.6: 0.160 → 0.506.
γ = 0.7: 0.070 → 0.377. γ = 0.8: 0.024 → 0.252. A larger bias −γ needs a wider interval to be
covered.

### Verdict: **SUPPORTED**

The hypothesis that the SE and the dispersion of τ̂ grow with τ is confirmed. The mechanism (ρ̂
error multiplied by τ − γ, with ρ̂ τ-invariant) is exact for this design. One nuance for the
coordinator: the rise in MSE is mostly variance, but not entirely. About 25% of the increase is
squared bias, from the (τ − γ)E[k] term.

---

## Q3. Chapter 2's 3×4 control-only block-stratified MSE ≈ 0.0004 vs its own β̂ ≈ β − ψ

Source: `~/GithubProjects/bios-dissertation/prelim/project-proposals/project1-spillover/draft/project1-spillover-draft.qmd`
(last commit touching it: c52a19e, 2026-09-10). The evidence is from the draft alone, as instructed.
We did not inspect Chapter 2's code.

### Finding: a genuine inconsistency within the draft; the definitions do not reconcile it

**The definitions.**
- **L537–540** (and L370–374, L1665): z_i = 1 "if cluster i is a Rook neighbor of at least one
  treated cluster and the spillover type permits it". Under Control Only (TrtNoSpill) only control
  clusters may have z_i = 1.
- **L439–441, L445–447**: block stratification requires that no two same-treatment clusters share a
  Rook edge.
- **L509–511**: on the 3×4 grid exactly 2 of 924 allocations satisfy it, "with a 'checkerboard'
  design". These are the two checkerboards.

**Consequence.** On the 3×4 checkerboard every control cluster has 2–4 Rook neighbors, all treated,
so z_i = 1 for every control. Under TrtNoSpill every treated cluster has z_i = 0. So **z = 1 − x
exactly** for both block-stratified allocations. Under TrtSpill a treated cluster would need a treated
Rook neighbor, and a checkerboard has none, so z = 1 − x in that regime too. With an intercept,
α + βx + ψz = (α + ψ) + (β − ψ)x. The draft's own statements then predict β̂ ≈ β − ψ, with bias ≈ −ψ
and MSE ≈ ψ² = 0.25 at ψ = 0.5, identically in both regimes:
- **L759–766**: "exactly collinear … a systematic negative bias approximately equal to −ψ … visible
  in the control-clusters-only results (Tables 2x4–3x4)".
- **L1351–1356**: "an exactly collinear allocation yields β̂ ≈ β − ψ".

**What the 3×4 table reports.** Block-stratified, Control Only, ρ = 0:
- **L1018**: MSE (×10) = 0.004 (so MSE ≈ 0.0004) and bias = −0.006 at ψ = 0.5.
- **L1019–1021**: the same at ψ = 0.6–0.8.
- **L974–976** (text): "BSS MSE ≈ 0.0004 vs. SRS MSE = 0.034 at ψ = 0.5".

A bias of −0.006 rather than −0.5, and an MSE that does not change with ψ, are incompatible with
z = 1 − x. No scaling of the "MSE (×10)" header (L1011) reconciles them: the bias column carries no
multiplier and is itself −0.006. By the draft's aggregation (**L1671–1673**), the reported quantity is
the MSE of β̂ over replications, averaged over the two block-stratified allocations. That is the same
quantity as in the 2×4 table, so a "different quantity" does not explain it.

**The 2×4 table is consistent with the definitions.** Block-stratified, ρ = 0 (**L887–890** Control
Only; **L897–900** Control & Intervention): bias −0.496, −0.597, −0.697 and −0.798 at ψ = 0.5–0.8 (≈ −ψ),
with identical values under both regimes. That is what z = 1 − x predicts.

**Related inconsistencies found in the same check:**
1. **3×4 regime difference.** Block-stratified at ρ = 0.01: Control Only MSE ≈ 0.035, bias −0.152
   (L1013); TrtSpill MSE ≈ 0.223, bias −0.415 (L1023). Under the draft's definitions the two regimes
   give the same z for a checkerboard, so they define the same model. L409–413 itself says regime
   differences should not appear under block stratification. (That sentence reads "does not satisfy
   block stratification conditions", but its parenthetical describes allocations that do satisfy
   them, so the wording looks inverted.) The draft doesn't define its SD column precisely enough to quantify
   Monte Carlo error. With σ² = 0.01 (L691), a gap of 0.19 in MSE is implausible as noise, but we did
   not verify that.
2. **"Checkerboard paradox" explanation.** At **L988–993** and **L1314–1321**, spillover flows "back
   into neighboring intervention clusters (which are diagonally adjacent but not Rook neighbors)". The
   Rook-only definition of z (L537–540, L1665) gives diagonal neighbors no spillover. The explanation
   therefore relies on a mechanism the model as written excludes.
3. **3×3 grid, same check.** The 6 block-stratified 3×3 allocations (L497–499; L750) include the
   strict one: the 4 edge-middle cells treated, with corners and centre as controls. Every control
   there has only treated Rook neighbors, so z = 1 − x. (We enumerated the other five: 4 of the 5
   corner/centre cells treated. Each leaves one control with no treated neighbor, so z ≠ 1 − x and β
   is identified.) With equal weighting over the 6 allocations, the strict one alone implies an MSE of
   at least ψ²/6 ≈ 0.042 at ψ = 0.5. The draft reports 0.001 (L914, L945), constant in ψ, with bias
   0.003. This is inconsistent in the same way, but less certain, because the draft does not state
   the weighting over allocations explicitly.

**Hypotheses we could not test from the draft** (listed so the Chapter 2 audit knows where to look;
none is asserted):
- The 3×3/3×4 code builds z or W with a different adjacency than stated (the L1316 diagonal wording
  hints at this).
- The allocations flagged as block-stratified in the code are not the checkerboards.
- The column order or aliasing handling in the `lagsarlm()` call differs from the 2×4 run.

### Verdict: **UNRESOLVED as a Chapter 3 note; it is a genuine Chapter 2 inconsistency**

Recommendation: fix it in Chapter 2, not Chapter 3. Audit the 3×3 and 3×4 simulation code: the z
construction, the neighbor definition, the identification of block-stratified allocations, and the
aggregation. Then either regenerate those tables or revise the collinearity text and the
"checkerboard paradox" passages. Chapter 3 needs no note. Its link (Discussion, "Relation to
Chapter 2") cites only Chapter 2's statement β̂ ≈ β − ψ. The 2×4 results support that statement, and
Chapter 3 proves the analogous result algebraically (Methods, Estimation: WZ = 1 − Z, so the
coefficient on Z estimates τ − γ). The same Discussion paragraph also quotes the "checkerboard
paradox" on the 3×4 grid (chapter L454). If the audit changes that Chapter 2 finding, that sentence
needs revisiting. **Author decision.**

---

## Q4. The rook Checkerboard regime gap (0.477 vs 0.483 at τ = 1)

### Proposed appendix paragraph

> Under rook contiguity, Checkerboard had an MSE of 0.477 under the both-arms regime and 0.483 under the control-only regime at $\tau = 1$, and a paired $t$-test over the 80 blocks gives $p = 0.009$. Yet for this design the two regimes define the same outcome model. Because $WZ = 1 - Z$, both spillover covariates equal $\gamma(1 - Z)$. The assignment is fixed and the incidence surfaces are shared, so a block's two regimes differ only in their noise, which is drawn from a seed whose key includes the regime. The two MSEs therefore have the same expectation, and their difference is Monte Carlo error. A paired test over blocks removes what the two members of a pair share, namely the surfaces, the assignment, and the values of $\rho$ and $\gamma$. It does not remove the regime-specific noise, so when there is no regime effect it rejects at its nominal rate, and a small $p$-value is a chance event rather than evidence of an effect. The shared surfaces do not invalidate this particular test. Given the surfaces, the two noise sets of a block are exchangeable, so each block's difference has conditional mean zero and the 80 differences are uncorrelated. To check the calibration, we drew 400 more independent noise sets for these 80 blocks. At $\tau = 1$ the 80-block mean MSE had a standard deviation of 0.0020 across noise sets, so the difference between two independent sets has a standard deviation of 0.0029. The observed gap, 0.0067, is 2.3 of these standard deviations, and paired $t$-tests between 200 disjoint pairs of independent noise sets gave $p \le 0.009$ for 2\% of the pairs. The observed gap is therefore an unusually large, but not implausible, Monte Carlo difference. Most of it traces to the both-arms noise set, in which the average $\hat\rho$ was higher than in 399 of the 400 replicate sets. The gap has the same sign at every value of $\tau$ and grows with $\tau$. This adds no evidence: the five values of $\tau$ share their noise, and the error in $\hat\rho$ is multiplied by $\tau - \gamma$ (Q2).

### Evidence

**Identical models by construction.**
- Spec §2 and chapter Methods: S(Z) = γWZ (both arms) or γWZ ∘ (1 − Z) (control only). Under rook
  WZ = 1 − Z exactly (asserted to 1e-12 in `q4_null_reps.R`), so both equal γ(1 − Z).
- Checkerboard is deterministic (spec §3), so no Z draw is involved.
- X_k is keyed by `("X", mode, ρX, k)` and shared by every block (spec §4; Invariant 1).
- **The noise key includes the regime.** Spec §4 gives `"eps", mode, ρX, nb, ρ, γ, regime`, and
  `code/05_run_simulation.R` `run_unit()` calls
  `set_seed_key("eps", cfg$mode, cfg$rho_x, nb_type, rho, gamma_val, spill_type)`.
- So Y differs between a block's two regimes only through ε, and the two ε matrices are independent
  draws of N(0, I).

**The observed comparison (stored; 80 blocks = 5 configurations × 4 ρ × 4 γ):**

| τ | MSE both arms | MSE control only | mean d | paired t p | Wilcoxon p | sign-flip permutation p | sign test | blocks with d > 0 |
|---|---|---|---|---|---|---|---|---|
| 0.8 | 0.474 | 0.480 | 0.0057 | 0.026 | 0.054 | 0.026 | 0.31 | 45 |
| 1.0 | 0.4766 | 0.4833 | 0.0067 | **0.0090** | 0.0136 | 0.0087 | 0.033 | 50 |
| 1.5 | 0.488 | 0.497 | 0.0092 | 7.2e-4 | 7.9e-4 | 7.4e-4 | 0.0049 | 53 |
| 2.0 | 0.505 | 0.517 | 0.0117 | 1.3e-4 | 2.3e-4 | 1.7e-4 | 0.0097 | 52 |
| 3.0 | 0.5584 | 0.5751 | 0.0167 | 5.6e-5 | 9.0e-5 | 7.0e-5 | 0.0011 | 55 |

The "p = 0.009" in the step-3 draft is the paired t-test at τ = 1. The permutation test uses 200,000
random sign flips with seed 20260925. d = control-only MSE − both-arms MSE.

**Why the gap grows with τ (re-sim).** Using Q2's exact form τ̂ − τ = −γ + c + (τ − γ)k, the
per-block regime differences at τ = 1 break down as follows:

| Quantity | Both arms | Control only | Mean difference | Paired t p (80 blocks) | Blocks with difference > 0 |
|---|---|---|---|---|---|
| Mean k (the ρ̂ error, scaled) | −0.00583 | −0.00905 | −0.00322 | 0.0033 | 30 of 80 |
| Mean c | 0.00249 | −0.00103 | −0.00352 | 0.080 | 38 of 80 |
| Var(k) | — | — | −0.00005 | 0.82 | — |
| Var(c) | — | — | 0.00043 | 0.38 | — |

The gap is a shift in the mean error (both components push the control-only bias more negative), not
a difference in dispersion. The k part is multiplied by (τ − γ), which is why d and its test
statistic grow with τ.

Correlation of the 80 block differences across τ (stored): 0.99 (τ = 0.8 vs 1), 0.64 (1 vs 3),
0.54 (0.8 vs 3). The five tests are not independent, as the chapter already says of the τ levels
(Methods; Limitations, fifth point).

**Shared surfaces.** d_b = f(X_c, ε_b^co) − f(X_c, ε_b^both). The two ε are i.i.d. given X_c, so
E[d_b | X_c] = 0. Distinct blocks use distinct ε keys, so the d_b are conditionally independent given
X_c. Hence Cov(d_b, d_b′) = E[Cov(d_b, d_b′ | X)] + Cov(E[d_b | X], E[d_b′ | X]) = 0. The same argument
makes the sign-flip test exact given X. Empirically, a one-way ANOVA of d on configuration gives
p = 0.948 at τ = 1 (0.79–0.95 across τ). This differs from the design comparisons, where
E[d | X] ≠ 0 and the shared surfaces do induce dependence. The chapter already handles those with the
50 configuration-by-surface units.

**Null calibration (400 new noise sets, key namespace `eps_a5rep`; 80 blocks each).** Each set was fit
once at τ = 1. Other τ were obtained exactly from τ̂ − τ = −γ + c + (τ − γ)k.
- τ = 1:
  - 80-block mean MSE across sets: mean 0.4806, SD 0.0020, range 0.4750–0.4876.
  - The observed both-arms value (0.4766) is at z = −1.96 and the control-only value (0.4833) at
    z = +1.33.
  - SD of the difference between two independent sets: 0.0029. The observed d = 0.0067 is 2.28 SD;
    3.0% of 200 disjoint replicate pairs had |d| at least this large.
  - Paired t-tests over the 200 disjoint pairs: P(p ≤ 0.05) = 0.080, P(p ≤ 0.01) = 0.020,
    P(p ≤ 0.009) = 0.020; SD of t = 1.10; KS test of uniform p-values, p = 0.40. The test is roughly
    calibrated, possibly slightly liberal. With 200 pairs, the binomial SE of 0.08 is about 0.019.
    2.0% of pairs had |t| ≥ 2.68 (the observed value).
  - 80-block mean k across sets: mean −0.00805, SD 0.00076. The both-arms set (−0.00583) is at
    z = +2.93, exceeded by 1 of 400 sets. The control-only set (−0.00905) is at z = −1.32.
- τ = 3:
  - Mean MSE 0.5686, SD 0.0032. Observed both arms 0.5584 (z = −3.17), control only 0.5751
    (z = +1.99).
  - Observed d = 0.0167 is 3.59 SD of the replicate difference; no replicate pair of the 200 reached
    it. P(p ≤ 0.05) = 0.055 and P(p ≤ 0.01) = 0.015 across pairs.
  - This is the same noise as at τ = 1, amplified through (τ − γ)k.

**Negative controls (stored, τ = 1).** Isolation Buffer's two regimes are also identical in expectation
(Z ∘ WZ = 0), although its Z draws also differ by regime key.
- Queen: 0.161 vs 0.159, paired t p = 0.53.
- Rook: 0.131 vs 0.134, p = 0.12.
- For rook Isolation Buffer, an ANOVA of d on configuration gave p = 0.003. We did not investigate
  this. Its four blocks with one aliased fit each fall under different configurations. It doesn't
  bear on the argument above, which holds for Isolation Buffer too.

### Verdict: **SUPPORTED**, with a caveat on how to phrase it

That the gap is Monte Carlo error follows from the construction. It is exact, not statistical: same
model, same assignment, same surfaces, independent noise. The paired test is valid and roughly
calibrated, and the observed p reflects an unusually large draw (about the 2nd percentile of the
null). If the appendix quotes p = 0.009, it should say that the test is calibrated and that the gap
is a roughly 1-in-50 Monte Carlo event. It should not describe the gap as "small". The chapter's
current sentence ("differ only by Monte Carlo error, because the two regimes use separate noise
draws") is accurate as written.

---

## Other likely committee questions (seed for the Q&A companion)

Pointers are to the chapter's section headings in `Dissertation_Chapter.qmd` as of this writing.

1. **Why is queen contiguity primary and rook only a sensitivity?** Under rook, τ is not identified for
   Checkerboard (WZ = 1 − Z), so only queen gives an identified comparison of all six designs; rook is
   kept for continuity with Chapter 2. → Methods, *Spatial structure and outcome model*; Methods,
   *Estimation*; Results, *Rook contiguity*.
2. **Why is τ = 1 the headline?** It is the pre-specified primary scenario. Under queen the design
   order is the same at every τ ∈ {0.8, …, 3}, although the τ levels share draws and are not
   independent replications. → Results, *Stability across treatment effect sizes*; Limitations (fifth
   point).
3. **Would the ranking hold beyond a regular 10×10 grid (for example, the 58 community college
   service areas)?** Untested. Some results depend on geometry (Isolation Buffer's treated count; the
   36 boundary clusters that identify queen Checkerboard under the control-only regime), and applying
   the designs to the real county surfaces is the next step. → Limitations (fourth point); Future work
   (third and fourth); Application.
4. **Why an oracle estimator, which assumes the analyst knows the spillover form?** It isolates the
   design's contribution under a correctly specified model. The non-oracle fit shows the cost of
   omitting the term (bias −0.12 to −0.40, coverage 0.51–0.90 at τ = 1), and Incidence-Guided
   Saturation Quadrants still had the lowest MSE. A misspecified exposure term was not studied.
   → Methods, *Estimation*; Results, *Non-oracle sensitivity analysis*; Limitations (first point).
5. **Why is coverage about 0.94 rather than 0.95?** The intervals are Wald intervals from the ML
   information matrix with N = 100. In our re-simulation, ρ̂ was biased downward by about 0.06 (this
   note, Q1); we have not attributed the undercoverage to that. The Q&A doc should call this
   finite-sample behavior of asymptotic intervals, not a design effect. → Results, *Design performance
   under queen contiguity*. (PARTIAL: a direct link between the ρ̂ bias and coverage is not
   established.)
6. **Does the matched-surface assumption favor the incidence-guided designs?** Possibly. It assumes
   the planning incidence is the incidence that drives the outcome. A noisy-snapshot data-generating
   process is proposed future work. → Methods, *Heterogeneous baseline incidence*; Limitations (second
   point); Future work (first).
7. **Isolation Buffer treats only about 21 clusters under queen: is the comparison fair?** Each design
   is evaluated as it would be used, including its treated fraction. Isolation Buffer's lower power
   (0.731 vs 0.930 for Incidence-Guided Saturation Quadrants at τ = 1) reflects that, while its MSE is
   third. → Table 1; Table 2; Results, *Design performance under queen contiguity*.
8. **Why are Monte Carlo SEs computed from 10 surface means instead of 250 fits?** The surface is the
   primary sampling unit, and fits from the same surface share X_k. The block-level Friedman tests
   overstate the information, which is why adjacent pairs are also tested over 50 independent
   configuration-by-surface units. → Methods, *Simulation design, resampling, and metrics*; Results,
   *Formal statistical comparison*; Appendix A2.
9. **Why were Saturation Quadrants and Balanced Halves dropped?** They were statistically
   indistinguishable from their retained counterparts in the eight-design comparison. → Methods,
   *Treatment assignment designs*; Appendix A1.
10. **Why does Checkerboard swing from near-best (both arms) to worst (control only) under queen?**
    Each interior cluster has WZ = 1/2, so under the control-only regime the spillover covariate is
    collinear with the intercept and Z except at the 36 boundary clusters. The swing is variance, not
    bias. → Methods, *Estimation*; Results, *Sensitivity to … spillover regime*.
11. **Why did the non-oracle model have lower MSE than the oracle for queen Checkerboard?** Omitting
    the nearly collinear spillover column removes the variance inflation at the cost of bias. At
    τ = 1, queen Checkerboard's MSE was 1.087 for the oracle and 0.108 for the non-oracle, with a
    non-oracle bias of −0.206. The non-oracle figures come from a direct computation on
    `sim_results_MLEnonoracle_tau_sweep_combined_20260924_025509.rds`, pooled over 160 scenarios.
    → Results, *Non-oracle sensitivity analysis*.

## Moved to Q&A companion (2026-09-25)

The note below was the third note of the chapter's Technical notes appendix until 2026-09-25, when
the author removed it from the chapter (review of the prelim render). The body keeps the sentence
that rook Checkerboard's two regime MSEs (0.477 and 0.483 at τ = 1) "differ only by Monte Carlo
error", without an appendix pointer. The note is kept here, verbatim (LaTeX as it stood in
`Dissertation_Chapter.qmd`), for the Q&A companion. Its cross-references pointed to Table
`tab:regime` (MSE by spillover regime), Table `tab:seed-keys` (seed keys), the Rook contiguity and
Sensitivity to spatial dependence, spillover magnitude, and spillover regime subsections of the
Results, and the preceding note on rook Checkerboard coverage (still in the appendix). Evidence
and scripts: Q4 above.

**Why Checkerboard's two spillover regimes gave different MSEs under rook contiguity.** Under rook contiguity at $\tau = 1$, Checkerboard had an MSE of 0.477 under the both-arms regime and 0.483 under the control-only regime (Table \ref{tab:regime}), and a paired $t$-test of the difference over the 80 pairs of blocks, matched on incidence configuration, $\rho$, and $\gamma$, gives $p = 0.009$. For this design the two regimes define the same outcome model: because $WZ = 1 - Z$, both spillover covariates equal $\gamma(1 - Z)$ (Rook contiguity). The assignment is fixed and the incidence surfaces are shared, so the two blocks of a pair differ only in their noise, which is drawn from seed keys that include the regime (Table \ref{tab:seed-keys}). The two MSEs therefore have the same expectation, and their difference is Monte Carlo error. The shared surfaces do not invalidate the paired test: given the surfaces, the two noise sets of a pair are exchangeable, so each pair's difference has conditional mean zero, and the 80 differences are uncorrelated. To check the test's calibration, we drew 400 further independent noise sets for these 80 pairs, with seed keys distinct from those of the simulation. At $\tau = 1$, the 80-block mean MSE had a standard deviation of 0.0020 across these noise sets, so the difference between two independent sets has a standard deviation of 0.0029; the observed difference, 0.0067, is 2.3 of these standard deviations. In paired $t$-tests between 200 disjoint pairs of the new noise sets, which by construction have no regime effect, $p$ was at most 0.05 for 8\% of the pairs and at most 0.009 for 2\%. The test is therefore roughly calibrated, perhaps slightly liberal, and the observed difference is an unusually large Monte Carlo difference, of a size that arises about once in 50 such comparisons, not evidence of a regime effect. Most of it traces to the both-arms noise set, whose average of $(\hat\rho - \rho)/(1 + \rho)$ was higher than in 399 of the 400 new sets. The difference had the same sign at every value of $\tau$ and grew with $\tau$. That adds no evidence: the five values of $\tau$ share their noise, and the error in $\hat\rho$ is multiplied by $\tau - \gamma$ (previous note). The argument extends to Isolation Buffer, whose two regimes also define the same outcome model (Sensitivity to spatial dependence, spillover magnitude, and spillover regime). Its assignments, unlike Checkerboard's, are redrawn for each regime, because their seed keys also include the regime (Table \ref{tab:seed-keys}), but given the surfaces the assignments and noise of the two blocks of a pair are both exchangeable, so each pair's difference again has conditional mean zero. Its regime differences at $\tau = 1$ were not significant (paired $t$-test $p = 0.53$ under queen contiguity and $p = 0.12$ under rook contiguity).
