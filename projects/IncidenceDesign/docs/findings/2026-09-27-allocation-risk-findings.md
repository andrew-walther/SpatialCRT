# Allocation-specific MSE pilot — 2026-09-27

Author: Andrew Walther

This note continues the [authorized plan](allocation_risk_plan_2026-09-27.md)
and [Chapter 2 discussion](ch2_check_findings_2026-09-27.md). It concerns new
IncidenceDesign analyses only. Project 1 and bios-dissertation remain read-only;
no manuscript wording or exhibit placement has been changed.

## Main pilot findings

The saturation designs retain their substantive advantage over SRS under
**control-only spillover**, now on both mean error and estimated upper-tail error.
The advantage does **not** extend to both-arms spillover. This supports a conditional
recommendation, not a universal claim that systematic allocation outperforms SRS.

For incidence-guided saturation quadrants versus SRS under control-only spillover:
mean MSE is 0.14087 versus 0.21717 (35.1% lower), estimated q90 is 0.17475 versus
0.29798 (41.4% lower), and estimated worst-decile mean is 0.19358 versus 0.34856
(44.5% lower). Corrected RMS allocation SD is approximately 0.0166 versus 0.0544,
about 70% lower. This is evidence of a useful combination of low mean error and
reduced allocation risk in the selected settings. It is not merely a low-variance
result divorced from the level of MSE.

Under both-arms spillover, incidence-guided saturation quadrants has mean MSE
0.04575 versus SRS 0.04277 (7.0% higher), with estimated worst-decile mean 0.05823
versus 0.05444 (7.0% higher). Balanced Quartiles and Balanced Halves are close to
SRS. This pilot provides no reason to portray SRS as an unreliable choice in that
regime. Very small differences among those near-tied designs should not become
broad recommendations from only two X surfaces per configuration.

Plain Saturation Quadrants also performs strongly, with lower pilot control-only
mean/tail MSE than the incidence-guided variant. Thus the evidence supports the
saturation approach, and incidence-guided saturation remains strong among the
retained manuscript designs; this analysis does not establish that incidence
guidance itself adds an accuracy benefit or that it is uniquely optimal.

The tables below average conditional metrics across 40 fixed-X/scenario units per
regime. Q90 and worst-decile values are finite-outcome estimates, not deconvolved
true-risk quantiles; read them with the precision check below.

### Control-only spillover

| Design | Mean MSE | Corrected RMS SD | Mean conditional q90 | Mean conditional worst 10% |
|---|---:|---:|---:|---:|
| Checkerboard | 2.1287 | 0.0000 | 2.1287 | 2.1287 |
| High Incidence Focus | 0.2998 | 0.0093 | 0.3038 | 0.3038 |
| Saturation Quadrants | 0.1353 | 0.0168 | 0.1689 | 0.1859 |
| Isolation Buffer | 0.1632 | 0.0156 | 0.1977 | 0.2177 |
| 2x2 Blocking | 0.5812 | 0.1387 | 0.7943 | 0.9174 |
| Balanced Quartiles | 0.2166 | 0.0547 | 0.2958 | 0.3484 |
| Balanced Halves | 0.2168 | 0.0524 | 0.2955 | 0.3432 |
| Incidence-Guided Saturation Quadrants | 0.1409 | 0.0166 | 0.1747 | 0.1936 |
| Simple Random Sampling | 0.2172 | 0.0544 | 0.2980 | 0.3486 |

### Both-arms spillover

| Design | Mean MSE | Corrected RMS SD | Mean conditional q90 | Mean conditional worst 10% |
|---|---:|---:|---:|---:|
| Checkerboard | 0.0602 | 0.0000 | 0.0602 | 0.0602 |
| High Incidence Focus | 0.1673 | 0.0024 | 0.1695 | 0.1695 |
| Saturation Quadrants | 0.0452 | 0.0021 | 0.0538 | 0.0573 |
| Isolation Buffer | 0.1623 | 0.0145 | 0.1970 | 0.2139 |
| 2x2 Blocking | 0.0618 | 0.0083 | 0.0769 | 0.0854 |
| Balanced Quartiles | 0.0422 | 0.0020 | 0.0503 | 0.0543 |
| Balanced Halves | 0.0425 | 0.0016 | 0.0507 | 0.0542 |
| Incidence-Guided Saturation Quadrants | 0.0458 | 0.0020 | 0.0545 | 0.0582 |
| Simple Random Sampling | 0.0428 | 0.0016 | 0.0509 | 0.0544 |

Saturation Quadrants, incidence-guided saturation and Isolation Buffer each have
lower mean and estimated worst-decile MSE than SRS in all 40 control-only fixed
units; none has lower mean MSE than SRS in any of the 40 both-arms units. These are
descriptive shares, not population probabilities.

The control-only incidence-guided/SRS mean-MSE ratio ranges from 0.622 to 0.672
across the five incidence configurations; its estimated worst-decile ratio ranges
from 0.519 to 0.570. Under both-arms spillover these ranges are 1.065–1.078 and
1.051–1.090, respectively. Thus the aggregate regime distinction is not driven by
one incidence configuration in this pilot. Full configuration/corner tables are in
`results/allocation_risk/summary/pilot/`.

Checkerboard is not rescued by considering allocation variability: its implemented
fixed allocation has zero conditional allocation variance but control-only MSE
2.12866, versus SRS's estimated worst-decile mean 0.34856. This concerns the current
Project 2 implementation/settings, not a retraction of Project 1. High Incidence
Focus illustrates a different trade-off: control-only mean MSE exceeds SRS's, while
its estimated worst-decile mean is lower. Neither low variance alone nor low mean
alone summarizes every investigator's design criterion.

## Targeted precision check: 100 versus 400 outcomes

The refinement preserved every original assignment and fitted outcome and added
1,500,000 new fits. All 50 blocks completed: 2,000,000 fits including the original
500,000, with zero warnings/aliases and no failed point estimates or CIs. Combined
with the pilot, the substantive analysis generated **7,117,600 distinct outcome
fits**. Two signed negative corrected-variance estimates remain in the refinement.

The following compares the **same selected subset** at R=100 and R=400. It must
not be compared directly with the full pilot table as if only R changed across
all 720 blocks. These are incidence-guided/SRS ratios of fixed-unit averages.

| Regime | Outcomes per allocation | Mean MSE ratio | q90 ratio | Worst-decile ratio |
|---|---:|---:|---:|---:|
| Control-only | 100 | 0.6643 | 0.5820 | 0.5511 |
| Control-only | 400 | 0.6581 | 0.5549 | 0.5273 |
| Both-arms | 100 | 1.0876 | 1.0905 | 1.0922 |
| Both-arms | 400 | 1.0871 | 1.1041 | 1.0948 |

The control-only direction survives refinement: incidence-guided saturation has
mean MSE 0.14863 versus SRS 0.22586 and estimated worst-decile mean 0.18492 versus
0.35069. Thus the initial allocation-risk benefit is not removed by reducing
outcome Monte Carlo noise. Plain saturation remains similarly strong (mean
0.14380, estimated worst-decile mean 0.18316). Isolation Buffer also improves on
SRS (0.16576 / 0.20839). Balanced Quartiles is closer to SRS (0.22000 / 0.32907).
These results support a saturation-family recommendation under control-only
spillover, with exact tail effect sizes still treated as pilot estimates.

The both-arms comparison remains unfavorable for incidence-guided saturation:
mean MSE 0.04762 versus 0.04381, estimated worst-decile mean 0.05480 versus 0.05005.
Balanced Quartiles' estimated worst-decile ratio moves from slightly below to
slightly above one (R400 ratio 1.0066), reinforcing an approximate-tie interpretation.

Precision is improved, **not fully resolved**. At R400 the median noise fraction
of raw allocation-MSE variance is 0.382 for incidence-guided saturation and 0.083
for SRS under control-only spillover. Under both-arms it remains 0.828 and 0.753,
respectively. Individual upper-tail identities remain unstable: median split-half
top-decile Jaccard overlap is 0.176/0.538 for incidence-guided/SRS under control-only,
and 0.053/0.111 under both-arms. These diagnostics do not invalidate the large
control-only mean/tail contrast, but rule out claiming precise true-tail estimates,
exact worst-allocation identities or exhaustive protection guarantees. No further
expansion was needed to establish the useful regime-specific pilot conclusion;
broader population/tail claims would require a separately scoped larger study.

## Implications to discuss before manuscript edits

1. **Recommendation:** retain a strong, specific recommendation for saturation
   allocation under control-only spillover. The combination of average accuracy
   and allocation-risk evidence is more informative than a pooled rank alone.
   Incidence-guided saturation is the relevant retained manuscript design; the
   plain-saturation findings still need fair acknowledgment in the full comparison.
2. **SRS:** report it as an informative benchmark, with results separated by
   spillover regime. Its good both-arms performance is a finding, not a narrative
   defect to explain away. Do not imply that SRS merely relies on investigators
   getting lucky in every setting.
3. **Project 1 link:** the common question is whether restricting allocations can
   improve accuracy and limit poor-allocation risk. This pilot shows that the
   answer depends on the restriction and spillover regime in Project 2. The two
   projects differ in grid size, outcome/analysis unit, incidence structure,
   exposure specification and allocation support. Their design rankings need not
   transport unchanged. This is not evidence that heterogeneous incidence alone
   caused the difference, and the accepted-exhibit reproduction discrepancy
   remains unresolved. Preserve the accepted study's findings within its scope.
4. **Checkerboard:** low allocation variability does not rescue its current
   control-only performance. Whether to show this in the chapter body or an
   appendix remains an author decision; no design has been removed from analyses.
5. **References:** dissertation prose can say "Chapter 2". CTJ must stand alone,
   using "previous work" with the verified BMC Medical Research Methodology
   citation. Do not insert chapter/project labels in CTJ prose.
6. **Next analysis:** broader tail claims would require more X surfaces and, where
   necessary, more outcomes/allocations. The real-SUD NC application is still a
   separate pending study. Neither is implied complete by this pilot.

Possible interpretation for discussion, **not approved manuscript wording**:

> Under control-only spillover, saturation-based allocation improved average
> accuracy and reduced estimated poor-allocation risk relative to complete
> randomization; under both-arms spillover, complete randomization and balanced
> designs remained competitive.

## Question and interpretation

For fixed incidence X, a scenario, and treatment allocation a, the quantity of
interest is m(a|X) = E[(tau_hat − tau)^2 | a, X]. One fitted dataset produces a
squared error; independent repeated outcomes estimate the allocation's MSE.
Drawing allocations from each design then estimates the distribution of these
conditional MSEs. The objective is to compare both average accuracy and the risk
of drawing an allocation with poor accuracy.

This follows the allocation-risk question raised by the author from Project 1,
but the repeated datasets differ: Project 1 also resampled subject locations;
this pilot holds the 100-cluster grid and X fixed while resampling outcome noise.
It does not reproduce Project 1's design distribution or resolve the historical
accepted-exhibit discrepancy. The previous study's allocation-consistency
argument motivates this analysis without predetermining its result.

## Scope

- Oracle spatial-lag ML, queen weights, tau = 1, existing 10×10 grid and DGP.
- All five incidence configurations: iid; spatial rho_X = 0.2/0.5;
  Poisson rho_X = 0.2/0.5. First two original key-seeded X surfaces per configuration.
- Outcome rho = 0/0.5, spillover gamma = 0.5/0.8, both-arms/control-only regimes.
- All eight candidate designs plus SRS: 80 fixed-X/scenario units × 9 designs.
- Initially 100 assignment draws per stochastic design and 100 outcomes per
  distinct allocation. Checkerboard has one implemented assignment. High Incidence
  Focus is fixed conditional on X unless ties cross the treatment cutoff; such
  ties retain their existing randomization rule.
- Targeted precision check: same 100 assignment draws, increased to 400 outcomes,
  for Saturation Quadrants, Isolation Buffer, Balanced Quartiles,
  Incidence-Guided Saturation Quadrants and SRS; first X surface of each of the five
  configurations, rho = 0.5, gamma = 0.8, both regimes (50 design blocks).

This is a conditional pilot on ten selected incidence surfaces, not a replacement
for the full simulation, a new tau/rook sensitivity study, or the NC application.
The revised design comparison on the real SUD incidence surfaces remains pending.

## Metrics and finite-replication precision

Each allocation has MSE, bias, coverage, Monte Carlo SE and independent split-half
MSE estimates. Equal-weight fixed-unit summaries average the ten selected surfaces
and four rho/gamma corners within each spillover regime. Configuration and corner
tables retain those distinctions. Ratios to SRS are ratios of averages.

Repeated assignment draws that produce the same allocation share cached outcome
fits. Their original frequencies are retained. With frequencies f_i, n = sum f_i,
estimated risk m_i and outcome Monte Carlo variance v_i, the corrected allocation
variance is

    sum[f_i (m_i − mean_m)^2] / (n − 1)
      − sum[f_i (1 − f_i/n) v_i] / (n − 1).

This explicitly accounts for covariance from shared cached estimates. Signed
negative corrected values are retained as evidence of unresolved precision.
The reported RMS within-setting SD is the square root of the average signed
corrected variance, not a pooled SD across scenarios or the average of SDs.
A structurally deterministic design has zero allocation variance conditional on
X; that alone says nothing about whether its MSE is acceptably small.

The q90 and mean of the worst 10% use estimated allocation MSEs. Outcome noise
inflates apparent dispersion and can change which allocations enter the upper
tail. Split-half rank agreement and cross-selected tail evaluation diagnose this:
selecting allocations with one half and evaluating with the other removes that
half's selection noise, but does not directly estimate the true worst-10% risk.
The R=400 check assesses stability of design comparisons on the same allocations.
Reported maxima are sample-dependent, never an exhaustive worst-case guarantee.

Mean-MSE difference MC SEs combine independent design/block simulation streams
conditional on the selected X surfaces. They do not provide uncertainty for an
incidence-surface population or make the shared-X blocks independent population
replicates. No population p-values are reported.

## Reproducibility and code walkthrough

`code/16_allocation_risk.R` sources unchanged modules 01–04. It regenerates the
original X surfaces, draws assignments with a separate key namespace, holds each
assignment fixed, generates repeated Gaussian errors, forms SAR responses, and
fits the existing validated lean ML engine. It records unique allocation metrics,
draw frequencies, block summaries, and individual errors/SEs in local checkpoints.
Fixed-size noise batches preserve the original sequence when outcome replication
increases. Source/package/BLAS manifests prevent incompatible cache reuse.

`code/17_allocation_risk_summary.R` reads completed result objects and writes
stratified comparisons, SRS ratios, conditional Monte Carlo precision, and tail
diagnostics. It refuses incomplete or unmatched design/SRS units. Tables average
conditional metrics instead of mixing different scenario distributions into one
purported allocation-risk distribution.

`code/18_allocation_risk_precision.R` validates the original manifest and selected
scenario keys, loads the existing 100 outcomes, and adds 300 outcomes on exactly
the same assignments. It asserts equality of every original fit and saves the
extension in a separate directory. The pilot is preserved.

Inputs are existing model/design functions and explicit simulation parameters;
outputs are RDS objects and CSV tables under `results/allocation_risk/`. Individual
fit checkpoints, transient logs and redundant allocation/draw CSV exports are
gitignored; tracked `results.rds` objects preserve the full allocation/draw tables,
and compact summary CSVs and manifests are retained. `spdep` supplies existing
weights, `spatialreg` the reference
estimator in validation, `digest` stable keys/manifests, and `parallel::mclapply`
local workers. BLAS must use one thread per worker. This is not a SLURM submission.

From `projects/IncidenceDesign/code`:

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript tests/test_allocation_risk.R
Rscript tests/test_allocation_risk_summary.R
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript 16_allocation_risk.R pilot 8
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript 18_allocation_risk_precision.R
Rscript 17_allocation_risk_summary.R pilot
Rscript 17_allocation_risk_summary.R precision_R400
```

Do not change the locked simulation source and reuse its old checkpoints. More
outcomes improve estimation of each allocation's risk; more allocations improve
sampling of the design distribution; more incidence surfaces are needed for broader
generalization. Those are different precision questions.

## Validation and protected files

- Existing revision tests passed with `SKIP_EQUIVALENCE=1`; the unchanged engine's
  lengthy full equivalence suite was not repeated. The new pilot tests separately
  compare three actual fitted datasets with `spatialreg::lagsarlm`, with tau and
  SE agreeing within 1e-6.
- Pilot tests cover known MSE/variance arithmetic, shared-cache covariance,
  negative corrections, failed CI handling, reproducible assignment/outcome
  prefixes, exact original X surfaces, and actual nonconstant fitted errors.
- Summary tests cover ratio-of-averages versus average ratios, conditional MC SE,
  reporting-stratum isolation, unmatched/incomplete/nonfinite-block rejection,
  frequency weighting and cross-selected tail diagnostics.
- Independent GPT-6 Sol/high review cleared the methods, runner, tests and precision
  extension. Execution used GPT-6 Sol/high; summary work used GPT-6 Sol/medium.
- SHA-256 comparison verified all 18 protected original source/main-result files
  unchanged. Modules 01–05 and the 2026-09-27 full-run result files are preserved.
- No Project 1, bios-dissertation or manuscript files were modified. No pushes.

## Further analysis recommendations — discussion, not approved new runs

The author asked what additional analyses would improve clarity and robustness.
The following priorities are recommendations; no additional simulation or change
to the design rules is authorized by this discussion alone.

1. **Confirm allocation-risk findings across more incidence surfaces.** If downside
   protection will be a substantive manuscript claim, extend the focused comparison
   to all ten existing X surfaces per configuration, retaining both spillover regimes
   and the same selected parameter settings. Include SRS, both saturation designs,
   Balanced Quartiles and Isolation Buffer. Ten surfaces is a practical extension,
   not a guarantee of adequate population precision. Choose further outcome and
   allocation replication using explicit Monte Carlo precision targets, agreed
   before the extension. More outcomes refine each allocation's MSE; more allocation
   draws improve tail sampling; more surfaces address generalization. Do not expand
   every tau/neighbor/scenario combination simply to produce more fits. Exact maxima
   and identification of the individual worst allocations are low priorities.
2. **Explain why the regimes differ.** Use saved assignments to compare treatment
   versus spillover exposure overlap, residual treatment variation after accounting
   for intercept/spillover/incidence, and suitably scaled model-matrix conditioning.
   Relate these diagnostics to allocation MSE, especially for Checkerboard, SRS and
   the saturation designs. This could explain which geometries help distinguish the
   direct effect from spillover. These are diagnostics, not a complete SAR variance
   formula or causal attribution of the Project 1/2 difference to heterogeneity.
3. **Use existing results before adding sensitivity runs.** Present bias, coverage
   and failure/alias rates beside mean MSE, by regime, with Monte Carlo uncertainty.
   Summarize the already available non-oracle sensitivity in the same strata; it
   omits the spillover term and is not a check of every possible misspecification.
   Retain oracle ML as primary. Compare plain and incidence-guided saturation directly
   in the full-study results so the manuscript distinguishes support for saturation
   from evidence for an incremental benefit of incidence guidance. No new estimator
   contest or DIM comparison is proposed.
4. **Clarify implementation and resource constraints.** Report treated-cluster
   counts/fractions and incidence balance with performance. Isolation Buffer does
   not enforce the same 50/100 treatment budget as most grid designs; its risk
   comparison should make that explicit. In the NC application also distinguish
   cluster balance from population balance and define the allocation support,
   spatial grouping, and any deterministic assignment rule. A new budget-matched
   buffer or changed randomization rule would be a separate design decision, not
   a silent modification of the current study.
5. **Make the real-SUD NC application the next major robustness test.** Include SRS,
   report each yearly surface rather than only pooling, and preserve the regime
   distinction. Irregular service-area geometry and unequal populations provide a
   useful complement to the grid results. A focused optional sensitivity would
   construct allocations using an earlier year's incidence and evaluate simulated
   outcomes under a later year's incidence. Distinguish the historical planning
   signal from the covariate driving outcomes and agree analysis adjustment before
   implementation. This tests reliance on perfectly measured/current incidence;
   consecutive year pairs are related settings, not independent population samples.
   Real incidence informs a simulation-based application, not evidence of an
   observed intervention effect. The revised engine/designs still need integration
   and verification for the application before its existing results can be replaced.

Recommended sequence: first extract existing-result/assignment diagnostics (items
2–4); agree whether to promote the allocation-risk pilot to a manuscript analysis
and its precision target (item 1); then run the planned NC application with any
agreed historical-incidence sensitivity (item 5). Avoid a broad new parameter sweep
or reopening Project 1 as a prerequisite to writing. A no-spillover control or
additional grid sizes can remain secondary unless a specific claim requires them.

The emphasis on explicit aims/performance measures and Monte Carlo uncertainty is
consistent with Morris, White and Crowther (2019), *Using simulation studies to
evaluate statistical methods*, Statistics in Medicine,
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086). The specific priorities above
are judgments based on this project's completed analyses, not conclusions supplied
by that methodological reference.

## Checkpoints and continuation

Code/method checkpoints: `455a480` (authorized plan), `ea2452d` (runner/tests),
`ae4f4f6` (summary/precision scripts), `23ccfb2` (completed outputs and reporting
guards). Independent final review checked the findings against both result sets.

Next: agree the SRS narrative, Checkerboard placement and Chapter 2 linking sentence;
then update/render/review the chapter. The revised real-SUD application and the
standalone CTJ draft remain separate subsequent work. Do not extend or revise
Project 1 based on these notes.

Ready-to-paste continuation prompt:

> Continue IncidenceDesign. Read AGENTS.md, docs/plans/TODO.md and
> docs/plans/allocation_risk_findings_2026-09-27.md, with the linked Chapter 2
> findings for context. The authorized nested pilot and R400 precision check are
> complete: saturation designs improve mean and estimated allocation-tail MSE
> versus SRS under control-only spillover, while SRS remains competitive/better
> under both-arms spillover. Tail identities remain noisy; incidence guidance is
> not uniquely optimal. Discuss SRS framing, Checkerboard placement and the neutral
> Chapter 2 link before editing the chapter. Project 1 and bios-dissertation remain
> read-only. CTJ must cite the BMC paper and stand alone; the revised real-SUD
> application remains pending. Never push without asking.
