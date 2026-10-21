# Methods extensions: report and proposed manuscript text (2026-10-21)

**Status: draft for the author's review.** No chapter, CTJ, SI or slide file has been
edited. This report combines the 2026-10-20 results
([findings](methods_extensions_findings_2026-10-20.md)), the author's 2026-10-21
decisions (`DECISIONS.md`) and the 2026-10-21 reruns
([overnight plan](overnight_plan_2026-10-21.md)). Every number below was read from
the CSV or log named in Section 6. Section 3 lists three new questions for the author.
Section 4 holds the proposed text, and Section 5 the proposed Part II slides.

---

## 1. What ran on 2026-10-21

| Run | Result |
|---|---|
| Grid burden share (`code/22_grid_burden_share.R`) | Done. 5 incidence configurations × 2 neighbor types × 32 blocks × 9 designs, production surfaces and allocations |
| v2 checks (`test_methods_ext_v2.R`) | 8/8 pass |
| A3 rerun `a3v2`: saturation regions of ≥ 2 colleges, plus the t + df-corrected test | 432/432 tasks, 301,824,000 fits, 0 failed, 03:50–07:30 |
| Reproduction check vs the 2026-10-20 A3 run | **PASS after a correction to the check itself** (below) |
| Bayesian CAR at ρ = 0.5, chains four times longer (`a6bayeslong`) | 36/36 tasks, 7,200 MCMC fits, 0 errors, 07:35–10:40; same data as 2026-10-20 confirmed |
| Summaries (`summaries_v2/`), V1, V2, V8, V9 | V1 PASS (25 of 576 \|z\| > 2, 26.2 expected; none > 3). V9 PASS at 10:41 (queue ended `QUEUEDONE`) |

**The reproduction check failed once, because of an error in the check, not in the
data.** The first version required that no non-saturation task have a subset whose
regions had been rebuilt. But `Regions_Changed` describes the participant subset, so
it is set for every design run on that subset. The corrected check is stricter for
those designs: all seven non-saturation designs must match the 2026-10-20 run on every
subset, including rebuilt ones, in every count and every allocation. They do, bit for
bit, across 336 tasks. Only 10 saturation tasks differ, the ones whose subsets got new
regions. The new rule rebuilt the regions of only two distinct participant subsets:
one random N = 15 subset, which is shared across the four years and so counts four
times, and one incidence-weighted 2021 N = 20 subset (45 subset × design pairs over 9
designs). The queue was resumed at the Bayes step, and the queue
log records the failure, the correction and the resume (`queue_v2.log`,
`a3v2_reproduction_corrected.log`).

---

## 2. Findings, by question

### 2.1 How many colleges does the trial need?

Effects are in cluster-SD units (σ_c = 1), with random volunteering, κ = 0.8, and the
mean over 2018–2021. The three tests are:
- **t + df:** the corrected test, adopted 2026-10-21;
- **size-adj.:** size-adjusted, using a critical value from the simulation;
- **Wald.**

Results at ρ = 0 and ρ = 0.5 agree within one or two colleges
(`summaries_v2/a3_n80_summary.csv`).

| Design | Both-arms τ = 1: t + df / size-adj. / Wald | Control-only τ = 1.5: t + df / size-adj. / Wald |
|---|---|---|
| Simple Random Sampling | 35 / 35–36 / 30 | 28–32 / 33 / 19–22 |
| High Incidence Focus | 34 / 34 / 29–31 | 27–33 / 28–34 / 18–22 |
| Balanced Halves, Balanced Quartiles | 35 / 35 / 30–31 | 28–32 / 29–35 / 20–22 |
| Plain Saturation | 35 / 35 / 30 | 27–28 / 28–31 / 19–20 |
| Incidence-Guided Saturation | 35 / 35 / 30 | 26–27 / 28 / 20 |
| 2×2 Blocking | 36 / 36 / 31 | 44 to > 58 / 57 to > 58 / 25–27 |
| Checkerboard | > 58 / > 58 / > 58 | > 58 |
| Isolation Buffer | > 58 / > 58 / > 58 | 57 to > 58 / > 58 / 52–54 |

Other points:
- **τ = 1 under control-only spillover is not reached with all 58 colleges, for any
  design or test.**
- At τ = 2, every design except Checkerboard (19–20 colleges) and Isolation Buffer (36)
  reaches 80% power with 15–17 colleges.
- **Smallest detectable effect at N = 58** (t + df):
  - both-arms: 0.76–0.78 for every non-geographic design, 1.31–1.32 for Checkerboard,
    1.50 for Isolation Buffer;
  - control-only: 1.24–1.25 for the two saturation designs, 1.45–1.58 for SRS and the
    rank designs, 1.60–1.61 for 2×2 Blocking, above 2 for Checkerboard.
- **The same at N = 29:** 1.12–1.15 both-arms, 1.45–1.51 control-only, excluding the
  geographic designs.

**Why control-only power plateaus.** As nearly every college joins, almost every
control has a treated neighbor. The spillover column then moves almost in lockstep
with (1 − Z), so τ is harder to separate from the spillover effect. Non-participating
neighbors restore some exposure variation, and so do the saturation designs.

**Translation to trained professionals.** SRS, both-arms, ρ = 0.5, m = 25 per college,
ICC = 0.05 (`summaries_v2/a3_translation_m_icc.csv`):

| Individual effect | N = 29 | N = 58 |
|---|---|---|
| d = 0.2 | 0.40 | 0.69 |
| d = 0.3 | 0.71 | 0.96 |

Power is t + df; size-adjusted is within 0.01. Under control-only spillover, d = 0.3
gives 0.49 (N = 29) and 0.53 (N = 58) for SRS, and 0.51 and 0.64 for Incidence-Guided
Saturation.

### 2.2 The t + df-corrected test is correctly sized at every number of colleges

This is new on 2026-10-21 and is the most useful result of the night
(`summaries_v2/summaries_log.txt`). It covers random volunteering, all designs and
settings at each N:

| Colleges N | Wald size | t + df size | Wald coverage at τ = 1 | t + df coverage at τ = 1 (design × regime range) |
|---|---|---|---|---|
| 15 | 11.1–12.8% | 4.4–5.4% | 87.7% | 94.9% (94.5–95.6%) |
| 20 | 9.5–10.4% | 4.7–5.4% | 89.8% | 94.8% (94.4–95.3%) |
| 29 | 8.3–8.8% | 5.0–5.5% | 91.6% | 94.8% (94.4–95.1%) |
| 40 | 7.4–7.7% | 5.2–5.5% | 92.5% | 94.7% (94.2–94.9%) |
| 58 | 6.3–7.0% | 5.0–5.5% | 93.2% | 94.6% (94.0–95.1%) |

- At N = 58 this agrees with Stage 1 on the production fits (94.1–95.0% by design ×
  regime).
- The correction is τ̂ ± t(0.975, N − p) × SE × √(N/(N − p)), with p the number of
  fitted coefficients plus ρ.
- It fixes the small-sample Wald problem at every participation level, so the reported
  interval and test now agree with each other.
- Its colleges-needed numbers match size-adjusted power under both-arms spillover. They
  are 1–5 colleges lower under control-only spillover at τ = 1.5.

### 2.3 Does the design guidance hold when only some colleges join?

The MSE ratio to SRS at τ = 1 is unchanged by the region rule (the largest change is
in one range bound, 0.97 → 0.98). The 2026-10-20 table stands:
- Checkerboard and Isolation Buffer are worse than SRS under both-arms spillover at
  every N.
- The rank designs stay within a few percent of SRS.
- The saturation designs' control-only advantage shrinks with fewer colleges: 22% lower
  MSE than SRS at N = 58, 9–10% at N = 40, 4–5% at N = 29 and about 1% at N = 20. It
  disappears under incidence-weighted volunteering.
- Under control-only spillover, Isolation Buffer goes from SRS-like at N = 58 to 2–4
  times SRS's MSE at N ≤ 29.

### 2.4 Targeting: where does Incidence-Guided Saturation send treatment?

**NC, N = 58** (`stored/a1_burden_share.csv`). Rate share is the share of per-college
SUD rate in treated colleges; death share is the share of deaths.

| Design | Rate share | Death share |
|---|---|---|
| High Incidence Focus | 0.622 | 0.354 |
| Incidence-Guided Saturation | 0.511 | 0.492 |
| Plain Saturation | 0.504 | 0.496 |
| SRS | 0.500 | 0.497 |

The four frozen NC regions have similar mean incidence ranks, 24–34 of 58 in every
year, with a spread of 0.10–0.14 when the gap between the highest and lowest regional
mean rank is divided by 58. So ordering saturations by region has little to work with.
Even regions formed from rank quartiles, ignoring geography, would reach only a rate
share of 0.56–0.57 with a 0.8/0.6/0.4/0.2 saturation, and their death share would fall
to about 0.42. NC's high-rate colleges are small.

**Grid** (`results/grid_burden_share/`: queen, production allocations, mean over 10
surfaces and 32 blocks):

| Incidence | Quadrant spread, mean (range) | Incidence-Guided Saturation | Plain Saturation | High Incidence Focus | Gain vs Plain, mean (range) | corr(spread, gain) |
|---|---|---|---|---|---|---|
| iid | 0.16 (0.08–0.22) | 0.528 | 0.500 | 0.758 | +0.028 (+0.015 to +0.038) | 0.88 |
| spatial ρ_X = 0.2 | 0.17 (0.09–0.26) | 0.528 | 0.500 | 0.766 | +0.029 (+0.018 to +0.040) | 0.84 |
| spatial ρ_X = 0.5 | 0.18 (0.13–0.26) | 0.529 | 0.500 | 0.748 | +0.029 (+0.018 to +0.051) | 0.88 |
| Poisson ρ_X = 0.2 | 0.13 (0.07–0.22) | 0.522 | 0.501 | 0.750 | +0.021 (+0.010 to +0.036) | 0.99 |
| Poisson ρ_X = 0.5 | 0.20 (0.13–0.32) | 0.531 | 0.500 | 0.750 | +0.031 (+0.020 to +0.049) | 0.99 |

**Reading.** Incidence guidance moves 2–3 percentage points of incidence share toward
treated clusters on the grid, three to four times its NC gain (+0.007). Across
surfaces the gain tracks the quadrant spread almost linearly (correlation 0.84–0.99).
This supports the revised wording: incidence guidance pays off in proportion to how
much burden differs between regions, and NC sits at the low end of that range.

Two caveats:
- The gain is modest even on the grid. High Incidence Focus reaches 0.75–0.77 because
  it gives up randomization across the whole map.
- Grid clusters have equal populations, so rate share and death share coincide there.
  The NC tension between them does not arise on the grid.

### 2.5 NC formal design tests (consistency summary)

Unchanged from 2026-10-20 (`stored/a2_tests_*.csv`): Friedman χ² (df 8) 135.3, 153.8,
150.7 and 179.8 (p from 1e-34 to 2e-25) in the four β × regime settings, each over 24
blocks.
- Under control-only spillover, Plain Saturation and Incidence-Guided Saturation have
  the two best average ranks (1.29–1.79).
- One correction to the 2026-10-20 findings doc: Checkerboard is not last or next to
  last everywhere. Under β = 1 both-arms it is seventh of nine, behind High Incidence
  Focus and Isolation Buffer.

### 2.6 Robustness

- **Unequal trained professionals per college** (A4, m̄ = 25, ICC = 0.05, Wald
  power). Power falls by 0.097 on average across designs and settings. At d = 0.3 under
  both-arms spillover, the mean over the nine designs is 0.59 → 0.46 at N = 20 and
  0.88 → 0.77 at N = 58. MSE rises by a factor of 1.46–1.63 for every design. No design
  gains.
- **CAR fitting and CAR data** (A6, MSE ratio to SRS, mean over 24 blocks). The
  saturation designs stay 0.76–0.82 under control-only spillover in all four DGP ×
  fitter combinations. Checkerboard (2.83–3.09) and Isolation Buffer (4.11–4.45) are
  worst under both-arms spillover.
- **Spillover form** (A7):
  - under the count form the saturation designs remain best under control-only
    spillover (0.79–0.81);
  - under the binary form Isolation Buffer is never identified, and Checkerboard is not
    identified under control-only spillover;
  - the saturation designs are identified most often under binary control-only
    spillover (86–88% of fits) and are best when identified (0.70–0.71).
- **Bayesian CAR with longer chains** (ρ = 0.5, 20,000 burn-in + 100,000 iterations,
  same data; `sensitivity_summaries_v2/a6_bayes_long_log.txt`).
  - **Convergence improves.** 83.1% of chains now converge, against 46.2% for the
    same fits with the shorter chains: 87.0% vs 50.4% under both-arms spillover and
    79.1% vs 42.1% under control-only. The median effective sample size is 1,161
    (previously 419).
  - **Convergence still misses the target.** 449 fits have ESS < 400, and 898 (12.5%)
    have \|Geweke z\| ≥ 1.96 against about 5% expected by chance.
  - **The bias and under-coverage remain** even among converged chains:
    - both-arms: bias +0.14 and coverage 90.4% (SAR ML on the same data: +0.02,
      93.5%);
    - control-only: +0.22 and 87.7% (SAR ML: +0.03, 93.2%).

    So the 2026-10-20 problems at ρ = 0.5 come from fitting a CAR model to SAR-lag
    data, not from short chains.
  - **Design ordering, converged chains, MSE ratio to SRS.**
    - Under both-arms spillover the ends agree with SAR ML: Checkerboard 2.90 and
      Isolation Buffer 3.80 are worst, and High Incidence Focus 0.83 is best.
    - Under control-only spillover, Checkerboard (2.31) and 2×2 Blocking (1.34) remain
      worst.
    - The saturation advantage does not carry over at ρ = 0.5: Incidence-Guided
      Saturation 0.95, Saturation Quadrants 1.05 and Isolation Buffer 0.87, against
      SAR ML's 0.81, 0.86 and 1.01 on the same data.
    - Spearman correlation with the SAR ML ordering: 0.87 (both-arms) and 0.77
      (control-only).
  - **Reading.** Under strong SAR-lag dependence, a misspecified Bayesian CAR still
    identifies the poor designs but cannot resolve the control-only advantage of
    saturation. At ρ = 0, where the CAR model is close to correct, it agreed with ML on
    2026-10-20. Slides and text should say this plainly rather than claim that the
    ordering holds everywhere.

### 2.7 Checks

| Check | 2026-10-21 status |
|---|---|
| V1 reproduction at N = 58 (a3v2) | PASS (25 of 576 \|z\| > 2, 26.2 expected; max 2.96) |
| V2 Wald size | Fails for N ≤ 35, as before. The t + df test is 4.4–5.6% at every N |
| V6 regions | Fixed: every a3v2 subset has regions of ≥ 2 colleges. Two distinct subsets were rebuilt: one random N = 15 subset, used in all four years, and one incidence-weighted 2021 N = 20 subset |
| V8 monotonicity | 2,397 of 18,144 curves dip, now counting three tests. Unchanged in kind: all dips are along N. Reported with the subset-bootstrap intervals (author decision) |
| Reproduction of the 2026-10-20 run | PASS (corrected check, Section 1) |
| V9 manifests | PASS (10:41, after all runs) |

---

## 3. Questions for the author

1. **Headline test.** You chose size-adjusted power as the headline, with t + df
   power beside it. The run shows the t + df test holds its 5% size at every N and is
   the test a trial analyst would actually use. **Proposal:** make t + df power the
   headline, use size-adjusted power as a check that agrees with it, and keep Wald in
   the SI. The drafts in Section 4 are written this way. If you keep size-adjusted as
   the headline, only the numbers in the colleges-needed sentences change, by 0–5
   colleges.
2. **Design name.** AGENTS.md lists "Plain Saturation" as the canonical name, but
   Chapter 3 and the SI say "Saturation Quadrants" (22 and 23 uses), and the CTJ main
   text uses neither name. The drafts use "Saturation Quadrants" to match the
   manuscripts. Which name should win everywhere?
3. **Push.** Six commits on main are unpushed (23f1d7a, 72d5041, c63d96b, 40efbf8,
   5922cf8, plus the one that records this report). Push them?

---

## 4. Proposed manuscript text

Conventions: Chapter 3 first, then the CTJ (both layouts) and the SI, with the two
verifiers run after each. Drafts use "we", no contractions, and the manuscripts' design
names. "Replace" quotes the current sentence, which is unchanged until you approve.

### 4.1 Chapter 3

**(a) Application Methods, new subsection after "Outcome model and Monte Carlo
evaluation": "Partial participation and sample size".**

> All 58 colleges will be invited after the pilot, but not every college will join. We
> therefore repeated the application for subsets of $N$ participating colleges,
> $N \in \{15, 20, 25, 29, 30, 35, 40, 45, 50, 58\}$, drawn at random (50 subsets per
> $N$, shared by every design and year). In a sensitivity analysis, a college's chance
> of joining rose linearly with its incidence rank, so the highest-burden college was
> three times as likely to join as the lowest ($N = 20, 29$). Each design was rebuilt on
> the participants with its own rule; the four saturation regions were required to
> contain at least two colleges each. Spatial dependence links participating colleges
> only, so a participant with no participating neighbor has a zero row in $W_P$.
> Spillover exposure counts treated participating neighbors over all geographic
> neighbors, $e_j(Z) = \sum_{k \in \mathcal{N}(j)} Z_k / |\mathcal{N}(j)|$, so a
> non-participating neighbor counts as untreated. Spillover scales with the effect,
> $\gamma = \kappa\tau$ with $\kappa \in \{0.5, 0.8\}$, and the effect $\tau$ is in
> units of the college-level outcome standard deviation, $\tau \in \{0, 0.2, \ldots, 2\}$.
> For each setting we estimated power, and from it the number of colleges needed for
> 80\% power and the smallest effect detectable with 80\% power at each $N$. A
> college-level effect $\tau$ corresponds to an individual-level standardized effect
> $d = \tau\sigma_c$, where $\sigma_c = \sqrt{\rho_I + (1 - \rho_I)/m}$ for $m$ trained
> professionals per college and intraclass correlation $\rho_I$.

The plan uses "ICC". The symbol ρ_I keeps it apart from the spatial ρ. **Decide:** ρ_I,
or spell out "ICC".

**(b) Application Methods, estimation paragraph: add the interval.**

> We report $\hat\tau \pm t_{0.975,\,n-p}\,\widehat{\mathrm{SE}}\sqrt{n/(n-p)}$, where $n$
> is the number of participating colleges and $p$ the number of fitted coefficients
> plus $\rho$, and we test $\tau = 0$ with the same reference distribution. The maximum
> likelihood variance divides by $n$, and the normal reference ignores the estimated
> parameters. With 58 or fewer clusters, both make the Wald interval too narrow.

**(c) Application Results, new subsection "Planning for realistic participation".**

> **Colleges needed.** Under both-arms spillover, the design hardly affected sample
> size: every design except Checkerboard and Isolation Buffer needed 34--36 colleges
> for 80\% power to detect a one-standard-deviation college-level effect and 18--19 for
> 1.5 standard deviations. Checkerboard and Isolation Buffer did not reach 80\% power
> at $\tau = 1$ even with all 58 colleges. Under control-only spillover, no design
> reached 80\% power at $\tau = 1$ with 58 colleges. At $\tau = 1.5$ the two saturation
> designs needed 26--28 colleges, against 28--32 for SRS and the rank-based designs.
> Power under control-only spillover rose slowly beyond about 40 colleges: when nearly
> every college participates, nearly every control has a treated neighbor, so
> spillover exposure varies little and is hard to separate from treatment.
>
> For a trial with $m = 25$ trained professionals per college and $\rho_I = 0.05$, an
> individual-level effect of $d = 0.3$ gives power 0.71 with 29 colleges and 0.96 with
> 58 under both-arms spillover (SRS). Under control-only spillover the same effect
> gives 0.49 and 0.53 with SRS, and 0.51 and 0.64 with Incidence-Guided Saturation.
>
> **Interval calibration.** The Wald test rejected a true null in 11--13\% of trials
> with 15 colleges and 6--7\% with 58. With the corrected interval, the rejection rate
> was 4.4--5.6\% and the coverage at $\tau = 1$ was 94.0--95.6\% for every design and
> regime at every number of colleges.
>
> **Does the guidance hold?** The regime-specific guidance held with fewer
> participants, but the size of the saturation advantage depended on participation.
> Checkerboard and Isolation Buffer remained worse than SRS under both-arms spillover
> at every $N$, and the rank-based designs remained within a few percent of SRS.
> Under control-only spillover, the saturation designs' mean MSE was 22\% below SRS
> with 58 colleges, 9--10\% below with 40, 4--5\% below with 29 and about 1\% below
> with 20. With few participants, non-participating neighbors already create the
> exposure variation that saturation otherwise supplies.

**(d) Application Results, targeting paragraph, after the budgets paragraph.**

> Incidence-Guided Saturation sends more treatment to high-incidence regions when
> incidence differs between regions. On the NC map it differs little: the four
> regions' mean incidence ranks lie between 24 and 34 of 58 in every year. On average,
> Incidence-Guided Saturation treats colleges holding 51.1\% of the summed college
> incidence rates, against 50.4\% for Saturation Quadrants and 62.2\% for High
> Incidence Focus. Shares of deaths run the other way (49.2\%, 49.6\% and 35.4\%),
> because the highest-rate colleges serve small populations. On the simulation grid,
> where quadrant incidence differed more, incidence guidance raised the treated share
> of incidence by 2.1--3.1 percentage points over Saturation Quadrants. Across
> surfaces, the gain rose with the spread of quadrant mean ranks (correlation
> 0.84--0.99).

**(e) Discussion, Summary: replace the incidence-guidance sentence.**

Current: "…incidence guidance lets a trial direct its highest treatment saturations to
the highest-incidence regions with little or no loss of estimation accuracy."

Proposed:

> …incidence guidance lets a trial direct its highest treatment saturations to the
> highest-incidence regions with little or no loss of estimation accuracy. How much
> treatment it redirects depends on how much incidence differs between regions: a few
> percentage points of incidence share on the grid, and less than one on the NC map,
> whose regions have similar incidence.

**(f) Discussion, Limitations.**
- **Third (unequal populations): add a sentence.**
  > In the application, letting the number of trained professionals grow with college
  > population reduced power by about 0.1 for every design and changed no design's
  > relative standing (Appendix A9).
- **Eighth (coverage): replace the first two sentences.**
  > The 95\% Wald intervals under-covered in both the grid and the application, and
  > the Wald test was anti-conservative with few colleges. A small-sample correction,
  > a $t$ reference distribution with a degrees-of-freedom-corrected standard error,
  > restored coverage of 94--96\% for every design in the application at every number
  > of participating colleges (Appendix A9). We did not apply it to the grid study,
  > whose per-fit estimates were not stored.

  Keep the rest of the paragraph.
- **New ninth.**
  > The participation analysis assumes that colleges volunteer at random or with
  > probability rising with incidence, and that non-participating colleges neither
  > transmit nor receive spillover.

**(g) Future work: replace "unequal participant precision" with the work that
remains.**

> …noisy planning incidence, effect heterogeneity, exposure networks beyond contiguity
> and a broader set of graph geometries.

**(h) New appendix section, A9 "Methods extensions for the NC application".**
Subsections, each with one table or figure taken from the CSVs in Section 6:
1. Colleges needed and minimum detectable effect, all designs, both κ, both ρ, under
   all three tests (Wald, size-adjusted, corrected).
2. Translation to m × ICC.
3. Guidance under partial participation, including the incidence-weighted mechanism
   and the region rule.
4. Burden share treated: the NC table and the grid table of Section 2.4.
5. NC Friedman, Nemenyi and Holm–Wilcoxon tests, labeled a consistency summary.
6. Unequal trained professionals (A4).
7. Interval calibration: the Stage 1 and Stage 2 tables, and the coverage-by-N table
   of Section 2.2.
8. Estimation-method sensitivity: CAR, reverse misspecification and Bayesian CAR,
   including the long-chain rerun. State that at ρ = 0.5 the Bayesian CAR is biased,
   under-covers and does not reproduce the control-only saturation advantage.
9. Spillover form: binary and count.

### 4.2 CTJ main text

The body is at 3,362 of 3,500 words, so the additions are about 120 words in total.
Power stays in prose, and its figure goes in the SI (D6).

- **Methods, NC application, one sentence.**
  > We repeated the application for random subsets of 15--58 participating colleges to
  > estimate the number of colleges needed for 80\% power, using a $t$-based interval
  > with a degrees-of-freedom-corrected standard error.
- **Results, NC application, two sentences.**
  > Under both-arms spillover, every design except Checkerboard and Isolation Buffer
  > needed 34--36 colleges for 80\% power at a one-standard-deviation effect. Under
  > control-only spillover, no design reached 80\% power at that effect with all 58
  > colleges, and the saturation designs' advantage over SRS fell from 22\% lower MSE
  > with 58 colleges to 4--5\% with 29.
- **Limitations, replace** "Wald coverage remained below nominal even when Monte Carlo
  precision was adequate. Interval calibration must be addressed before adopting this
  analysis for a final protocol."
  > Wald intervals under-covered, but a $t$-based interval with a
  > degrees-of-freedom-corrected standard error restored 94--96\% coverage at every
  > number of participating colleges (Supplementary Section S11).
- **Implications or Conclusions.** The abstract's "Incidence-guided saturation directs
  the highest treatment saturation to high-burden regions while maintaining accuracy"
  stays true. Add after it, only if the word count allows:
  > , and the redirection grows with regional differences in burden.

### 4.3 SI

A new **Section S11, "Methods extensions for the NC application"**, mirroring
chapter A9 (subsections 1–9), with S-numbered tables. Main figure candidates:
- power curves by N for SRS, both saturation designs and Checkerboard, both regimes;
- the grid quadrant-spread vs targeting-gain scatter.

Both need small new plotting scripts. No figure exists yet.

---

## 5. Proposed Part II slides

Added to the ideas bank in `bios-dissertation/slides/prelim/OUTLINE.md`
(Part II, "Methods extensions"), so they can be slotted later.

**Main**
1. *About 35 colleges detect a one-SD effect under both-arms spillover, whatever the
   design.* Colleges-needed bar chart, t + df, τ = 1 and 1.5. New figure.
2. *Under control-only spillover, adding colleges stops helping.* Power vs N for SRS
   and Incidence-Guided Saturation, both regimes, τ = 1. New figure.
3. *Saturation's advantage shrinks as fewer colleges join.* MSE ratio to SRS by N.
   Table from `a3_guidance_ratio_to_srs.csv`.
4. *Incidence guidance redirects treatment where burden is regional.* Grid
   spread-vs-gain scatter with the NC point marked. New figure.

**Backup**
1. *A t interval with corrected SE restores nominal coverage at every N.* Section 2.2
   table.
2. *Wald overstates power with few colleges.* Size by N, three tests.
3. *Rate share and death share pull apart in NC.* Section 2.4 NC table.
4. *Formal tests agree with the descriptive ranking.* Critical-difference plots
   (`stored/a2_cd_plots.pdf`).
5. *The guidance survives CAR fitting and a CAR outcome model.* A6 table.
6. *A Bayesian CAR fit flags the same poor designs, but under strong dependence it is biased.* A6 Bayes, long chains (Section 2.6).
7. *Binary spillover makes Isolation Buffer unidentifiable.* A7 table.
8. *Unequal college sizes cost power for every design equally.* A4.
9. *From college-level τ to trained professionals.* Translation table.

---

## 6. Sources

All paths are under `application/results/methods_ext_2026-10/` unless they start with
`results/`.

| Numbers | File |
|---|---|
| Colleges needed, MDE (three tests) | `summaries_v2/a3_n80_summary.csv`, `a3_mde_summary.csv`, `a3_headline_tables.md` |
| Size and coverage by N (t + df) | `summaries_v2/summaries_log.txt`, `v2_size.csv`, `a3_settings.csv` |
| Translation | `summaries_v2/a3_translation_m_icc.csv` |
| Guidance under participation | `summaries_v2/a3_guidance_ratio_to_srs.csv` |
| NC burden share | `stored/a1_burden_share.csv`; per N `summaries_v2/a1_burden_share_per_N.csv` |
| Region mean ranks | `application/results/real_sud_rev_20261002/setup.rds` (recomputed 2026-10-21) |
| Grid burden share | `results/grid_burden_share/grid_burden_share_summary.csv`, `grid_quadrant_spread.csv` |
| NC tests | `stored/a2_tests_friedman.csv`, `a2_tests_avg_ranks.csv` |
| Stage 1, Stage 2 coverage | `stored/a5_stage1_coverage_by_design_regime.csv`, `sensitivity_summaries/a5_stage2_coverage_by_design_regime.csv` |
| A4 | `summaries/a4_summary.csv` |
| A6 CAR, A7 | `sensitivity_summaries/a6_car_ratio_to_srs.csv`, `a7_ratio_to_srs.csv` |
| Bayesian CAR, long chains | `sensitivity_summaries_v2/a6_bayes_long_log.txt`, `a6_bayes_long_*.csv` |
| Reproduction | `a3v2_reproduction.log` (first version, failed), `a3v2_reproduction_corrected.log` |
