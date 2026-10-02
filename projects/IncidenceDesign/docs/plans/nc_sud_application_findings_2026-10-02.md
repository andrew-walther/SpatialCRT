# Completed NC SUD education application — 2026-10-02

The observed-incidence application and its focused allocation-risk extension are
complete and independently verified. Thesis/appendix integration and CTJ/SI
derivation remain unfinished. No restricted county data were released or staged.

## Verified study and provenance

The four observed 2018–2021 SUD incidence settings and 58 named college service-area
clusters are fixed. Primary simulated education outcomes use
Y = (I − ρW)⁻¹[τZ + S(Z) + ε], τ = 1, σ = 1, without an incidence baseline term.
Incidence informs allocation. Queen weights are primary; both spillover regimes
and all eight strategies plus SRS are reported. The matched β = 1/X-adjusted,
rook-corner and regional-summary sensitivities follow the approved plan.

Production has **1,248 reporting rows from 936 distinct distributions** and
**7,598,000 independent outcome fits**. Every setting is complete, with no
warnings, aliases, boundary fits or failures. All allocation draws use J100.
Outcomes use R100 for stochastic assignments, R1000 for proven singletons, and
R2000 in two matched-sensitivity HIF sources. Maximum relative MCSE(mean MSE)
is **0.0487115**; maximum MCSE(coverage) is **0.00891687**. Both approved gates
pass everywhere. Precision success does not establish nominal CI coverage.

The separate tail extension has **144 reporting rows from 96 sources** at
primary queen ρ ∈ {0,0.5}, γ = 0.8, both regimes/all years/designs. It retains
the same J100 allocations and extends outcomes to R400 or retains existing
singleton R1000. All tail settings also pass completeness and mean/coverage
precision gates, with no warnings or failures. The 3,060,000 tail outcomes include
780,000 copied main outcomes and **2,280,000 additional independent outcomes**.
Combined distinct outcome effort is therefore **9,878,000**, not 10,658,000.

Sources under `application/results/real_sud_rev_20261002/`:

- `production/{performance.csv,allocation_metrics.csv,results.rds,manifest.rds,verification.txt}`:
  complete parameter-grid performance and diagnostics.
- `tail_confirmation/`: separate refined corner risks, cross-selected halves,
  source manifest and independent verification. Main outputs are preserved.
- `exhibits/yearly_primary_design_means.csv`: six equally weighted ρ/γ pairs per
  year/regime/design; no pooled MC intervals because allocation streams are shared.
- `exhibits/primary_setting_srs_comparison.csv`: within-setting SRS matches,
  mean-MSE differences/ratios and approximate Monte Carlo error intervals.
- `exhibits/sensitivity_setting_matched_comparison.csv` and
  `yearly_sensitivity_means.csv`: sensitivity comparisons with primary references
  matched at the exact same year/ρ/γ/regime/design. Rook's four corners must not
  be compared with unconditioned six-setting queen means.
- `exhibits/refined_tail_comparison.csv`: refined q90/worst-decile estimates,
  signed variance correction, conditional-noise fractions and half diagnostics.
- `exhibits/observed_cluster_inputs.csv`, `spatial_block_diagnostics.csv`,
  `exhibit_inventory.csv`, `exhibit_manifest.rds`, `manuscript_extract.txt`:
  authorized aggregate inputs, mechanical details and number/figure provenance.

## Mean accuracy, bias and coverage

The table below reports primary queen MSE averaged equally over the six ρ/γ
pairs, **separately for each observed planning year**. Values are descriptive.
SRS and plain saturation reuse identical primary performance across years because
their supports and the education outcome do not depend on incidence. Those
repeated values are not four independent annual confirmations.

| Year | Both arms: SRS | Both arms: guided saturation | Control only: SRS | Control only: plain saturation | Control only: guided saturation |
|---|---:|---:|---:|---:|---:|
| 2018 | 0.07638 | 0.07428 | 0.30761 | 0.22329 | 0.24285 |
| 2019 | 0.07638 | 0.07628 | 0.30761 | 0.22329 | 0.21565 |
| 2020 | 0.07638 | 0.07682 | 0.30761 | 0.22329 | 0.22029 |
| 2021 | 0.07638 | 0.07450 | 0.30761 | 0.22329 | 0.23732 |

Under control-only spillover, plain saturation reduces descriptive mean MSE by
**27.4%** relative to SRS; guided saturation's yearly reductions are
**21.1–29.9%**. Both have lower MSE in every primary setting, including the
setting-specific approximate Monte Carlo uncertainty. Neither is uniquely best;
the study does not establish that incidence guidance itself improves upon plain
saturation. HIF's control-only ratio to SRS varies substantially (about 0.697 in
2018, 1.398 in 2019 and 1.148 in 2021), precluding a uniform recommendation.

Under both arms, SRS is competitive with balanced, spatial-blocking and saturation
designs. Small point differences and few clearly separated within-setting MC
intervals do not establish equivalence or universal optimality. HIF's low primary
both-arms MSE is conditional on a deterministic, incidence-selected assignment
and the assumed education model, and is not evidence of randomization-based
causal validity for a future trial.

Primary yearly mean coverage is approximately **92.6–93.6%**, below nominal 95%.
Do not call it nominal coverage or equate precision-gate success with valid
95% interval coverage. Bias and coverage are retained for every setting in the
source tables; do not silently replace finite-sample undercoverage with a preferred
estimator or an unapproved additional sweep.

## Allocation-specific risk and its uncertainty

Risk is m(a) = Eε[(τ̂ − τ)² | a, observed inputs], evaluated across eligible
assignment masks. Its q90 is a quantile of **conditional MSE across allocations**,
not a quantile of realized treatment-effect errors within one allocation.

At the two prespecified queen ρ corners with γ = 0.8 under control-only spillover,
plain saturation's estimated q90 is **31.8–32.9% below SRS** and its estimated
worst-decile mean is **34.6–37.8% below**. Guided saturation's corresponding
reductions are **20.8–38.6%** and **20.7–43.4%**. Cross-selected estimates retain
the pattern when one outcome half selects poor allocations and the other
evaluates them. These findings support a qualified descriptive downside benefit.

They do not establish exact population tails, exhaustive protection or fresh
allocation confirmation. J100 supplies only ten sampled worst-decile draws.
All 112 non-singleton reporting rows flag uncertain membership at the tail
boundary. Both-arms identities are particularly unstable: guided saturation's
half Jaccard overlap is 0–0.25, with roughly 61–86% of raw estimated-risk variance
attributable to outcome MC noise. Plain saturation has overlap 0.11–0.25 and
about 69–71% noise. Small both-arms tail differences cannot support superiority.

Graph Checkerboard and HIF have proven singleton supports for these fixed inputs.
Their allocation q90 and worst-decile mean equal conditional mean MSE; Jaccard=1
is vacuous. Zero allocation variation does not imply zero estimation error,
good average performance or no outcome uncertainty. Signed variance corrections
are preserved; none is negative in the refined outputs.

## Budgets, adaptation and sensitivity

- HIF treats 29 clusters but approximately **22–24% of working-age population**.
  Buffer treats a mean **15.89 clusters** and about **27%** of that population.
  Shares measure geographic reach, not the number/proportion of healthcare
  professionals educated. Equal treatment-cluster counts need not imply equal
  participant or implementation budgets.
- Graph Checkerboard treats 27 under queen and 26 under rook. On this irregular
  graph τ is identified; its poor precision reflects weak treatment/spillover
  separation, not the rectangular-rook grid's exact non-identification.
- Frozen spatial-block budgets are **28/31/27/29** by year. Each year includes a
  singleton block assigned control under retained round(1/2)=0. Disclose that
  support restriction and avoid describing these irregular groups as literal 2×2
  squares. Other current rules/budgets were retained; Balanced Halves was corrected
  to 14/15 treatments with a random choice of the extra-treatment half.
- The matched β=1 baseline/X-adjusted sensitivity substantially increases HIF MSE
  (roughly **3.68–4.44×** under both arms; **1.42–2.52×** control-only). In a checked
  2018/ρ0.5/γ0.8/both setting, residual treatment SS declines from 14.50 to 3.60
  after adjustment. This is consistent with treatment–X alignment reducing
  information, but changes baseline and adjustment together; it is not an isolated
  β effect or a claim that SUD incidence predicts education ability.
- Matched rook/queen corner comparisons preserve broad regime-specific conclusions.
  Graph MSE increases roughly **11–25%** both-arms and **18–29%** control-only under
  rook. Match parameters before interpreting any neighbor difference.
- Mean raw rates and pooled regional deaths/population sometimes change the regional
  ordering. Mean-rate rules reuse primary sources with identical order/tie support;
  population-weighted rules add two distinct orderings. Preserve the prespecified
  mean-rank rule; do not select summaries retrospectively by performance.

## Review and manuscript implications

The author-requested independent reviewer audited all 688 pilot sources, reviewed
finished production/exhibits and matched sensitivities, and checked all 96 tail
sources against main caches with no allocation/draw/fit-prefix mismatch. Actual
graph/HIF queen/rook fits matched `lagsarlm` within numerical tolerance. Existing
data, allocation-risk/summary and simulation-revision behavior suites pass;
the latter reused its previously validated engine equivalence, while new actual
NC equivalence checks were performed. Application summary/tail tests also pass.

The companion at `application/report/real_sud_companion.html` now displays all
maps and traced figures, full-study/refined-evidence selectors and key findings.
Actual JavaScript passes a DOM-fixture behavioral check and static/syntax checks.
Local map/publication images were visually reviewed. In-app browser blocks file
URLs; real browser visual rendering remains unverified, not silently claimed.

Keep the selected-surface grid allocation-risk pilot distinct from these completed
NC results and from the full grid study. Represent Chapter 2's practical BSS
argument fairly within its studied settings; do not reopen Project 1's audit or
claim its historical discrepancy resolved. Recommendations remain regime-specific,
with budgets, model assumptions, undercoverage and tail uncertainty explicit.

Next: integrate into the single thesis chapter/appendix source, trace all numbers
and exhibits and render/sync/review; then derive complete standalone CTJ/SI with
verified Project 1 citation and live journal requirements. Author confirmed the
existing author/corresponding list, no current funding and no conflicts. The
IRB/data-use statement remains an explicitly authorized placeholder before
submission. Prepare the comprehensive ready-to-paste Claude review prompt last.
