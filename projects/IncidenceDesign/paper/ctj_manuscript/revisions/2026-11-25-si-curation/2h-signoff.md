# Phase 2h sign-off: plain rewrites, short list titles and light SI trim

Author: Andrew Walther. Drafted 2026-11-25. Nothing in the SI or the article has changed yet. Approve by ID.

Three parts, each measured in a real build of the curated SI (49 pages today):

- **NY and NYC items** rewrite every "X, not Y" contrast (and the few "rather than" contrasts that answer an unasked question) as a plain statement of what the thing is: 25 in the SI and 6 in the article. Where a rewrite also fixes a mid-sentence "therefore" or a ", so" join, the item says so. Several touch wording approved in 2g, which each item names.
- **SC items** shorten the List of Figures and List of Tables entries (the captions themselves do not change) and put both lists on one page: 49 to 46 pages.
- **PC items** are the low-hanging cuts the author asked for: text whose numbers a table or figure already gives, and a table that repeats the article's Table 1. With SC, PC-2, 4, 7, 8, 12, 13 and 14 give 44 pages, and adding PC-15 and PC-11 (author request) gives 43. Everything moved stays in the unabridged build.

Page counts move in whole pages because of float placement, so most PC items save nothing alone and work as a set.

## Summary

| ID | Recommendation | What |
|---|---|---|
| NY-1 | Accept | Overview: what the case study is |
| NY-2 | Accept | Budget table caption: what the shares measure |
| NY-3 | Accept | Population denominators |
| NY-4 | Accept | Row-standardization |
| NY-5 | Accept | The 100/0 boundary setting |
| NY-6 | Accept | Adjusted associations (Section 3.2) |
| NY-7 | Accept | What the service-area map is built from |
| NY-8 | Accept | The weights matrix as a proxy |
| NY-9 | Accept | Prespecified regional summaries |
| NY-10 | Accept | Mean ranks as an adaptation of the grid rule |
| NY-11 | Accept | Unequal treatment budgets |
| NY-12 | Accept | Target of the case-study analysis |
| NY-13 | Accept | The hypothetical incidence association |
| NY-14 | Accept | Weights of repeated allocations |
| NY-15 | Accept | What the Jaccard index describes |
| NY-16 | Accept | Case-study rank-test p-values |
| NY-17 | Accept | Flagged rook Checkerboard estimates (grid table text) |
| NY-18 | Accept | Source of the control-only increase |
| NY-19 | Accept | Shared draws across effect sizes |
| NY-20 | Accept | Rook Checkerboard estimates (rook section) |
| NY-21 | Accept | High Incidence Focus and randomization |
| NY-22 | Accept | Tail summaries and the reused assignments |
| NY-23 | Accept | Prespecified mean-rank rule |
| NY-24 | Accept | Matched sensitivity ratios |
| NY-25 | Accept | Case-study formal tests |
| NY-26 | Accept | Mid-sentence therefore (Wilcoxon tests) |
| NY-27 | Accept | Mid-sentence therefore (Isolation Buffer regimes) |
| NY-28 | Accept | Mid-sentence therefore (powers of W) |
| NY-29 | Accept | Mid-sentence therefore (error in the spatial coefficient) |
| NYC-1 | Accept | Abstract scope sentence |
| NYC-2 | Accept | Direct effect definition |
| NYC-3 | Accept | Units for Monte Carlo standard errors |
| NYC-4 | Accept | Planning years |
| NYC-5 | Accept | Checkerboard and Isolation Buffer accuracy |
| NYC-6 | Accept | Conclusions scope sentence |
| SC-1 | Accept | Short titles in the List of Figures and List of Tables |
| SC-2 | Accept | Both lists on one page |
| PC-2 | Accept | Table S1 (six design rules) moves; the article's Table 1 gives the same rules with example maps |
| PC-4 | Accept | Section 5.4 second paragraph (refinement schedule and frequency-weight bookkeeping) moves; Section 5.3 states the precision targets |
| PC-7 | Accept | Section 10 rook rank tests move; Table S6 and Section 12.2 keep the rook results and the comparison with SRS |
| PC-8 | Accept | Section 9 last paragraph (why rook Checkerboard coverage rises with β) moves; the coverage values stay in Section 8 |
| PC-12 | Accept | Section 6.2, two sentences that the article's Figure 2 shows |
| PC-13 | Accept | Section 11, the bias and coverage ranges and a closing sentence that Table S7 prints |
| PC-14 | Accept | Section 12.2 last paragraph; Table S9 and Sections 12.3 and 15 give its worst-decile ratios |
| PC-15 | Accept (author) | Section 8 second paragraph (six-design ranks and MSE by β); Figure S4 shows the curves, and the CTJ's β claim rests on the paragraph that stays |
| PC-11 | Accept (author) | Case-study Friedman/Nemenyi tests (Section 5.6 first paragraph and the 16.4 block); Tables S8, S14 and S15 show the ordering |

## Batch NY · SI rewrites

### NY-1. Overview: what the case study is

**Recommendation:** Accept.

**Current**

```latex
The case study is a simulation-based design evaluation, not observed intervention efficacy.
```

**Proposed**

```latex
The case study is a simulation-based evaluation of designs for a planned trial.
```

### NY-2. Budget table caption: what the shares measure

**Recommendation:** Accept.

**Current**

```latex
Shares are geographic reach, not healthcare-professional enrollment.}
```

**Proposed**

```latex
Shares measure geographic reach.}
```

### NY-3. Population denominators

**Recommendation:** Accept.

**Current**

```latex
These population
denominators measure geographic reach, not the number of healthcare
professionals educated, intervention cost or participant-level
precision.
```

**Proposed**

```latex
These population denominators measure geographic reach only.
```

### NY-4. Row-standardization

**Recommendation:** Accept.

**Current**

```latex
Row-standardization
expresses each cluster's spillover exposure as an average over its own
neighborhood rather than a raw sum, which would otherwise depend on how
many neighbors the cluster has.
```

**Proposed**

```latex
Row-standardization makes each cluster's spillover exposure an average over its own neighborhood, which keeps exposure comparable between clusters with different numbers of neighbors.
```

### NY-5. The 100/0 boundary setting

**Recommendation:** Accept. Touches Section 2.2, approved 2026-11-25; changed only for the rule.

**Current**

```latex
The 100/0 setting is a boundary case, included to show whether going beyond 90/10 helps or hurts, not as a design to use.
```

**Proposed**

```latex
The 100/0 setting is a boundary case, included only to show whether going beyond 90/10 helps or hurts.
```

### NY-6. Adjusted associations (Section 3.2)

**Recommendation:** Accept. Touches SR-2 wording approved in 2g. The sentence before it says that incidence is not randomized.

**Current**

```latex
Therefore, these coefficients describe
adjusted associations rather than effects of incidence.
```

**Proposed**

```latex
Therefore, these coefficients describe adjusted associations between incidence and the outcome.
```

### NY-7. What the service-area map is built from

**Recommendation:** Accept.

**Current**

```latex
The map is built from county boundaries and the county-to-college assignment alone, not from the SUD data.
```

**Proposed**

```latex
The map uses only county boundaries and the county-to-college assignment.
```

### NY-8. The weights matrix as a proxy

**Recommendation:** Accept.

**Current**

```latex
Legal boundaries shared across water
can create adjacency, so $\mathbf{W}$ is a prespecified proxy for geographic
exposure, not a measured network of communication among professionals.
```

**Proposed**

```latex
Legal boundaries shared across water can make two service areas adjacent. $\mathbf{W}$ is a prespecified geographic proxy for exposure between neighboring service areas.
```

### NY-9. Prespecified regional summaries

**Recommendation:** Accept.

**Current**

```latex
Mean raw rates and pooled regional rates were
prespecified as sensitivity analyses rather than chosen by performance,
and
```

**Proposed**

```latex
Mean raw rates and pooled regional rates were prespecified as sensitivity analyses, and
```

### NY-10. Mean ranks as an adaptation of the grid rule

**Recommendation:** Accept.

**Current**

```latex
Mean ranks agree with it in the rank-scaled
Poisson configurations but are a deliberate adaptation for the
case study, not the exact grid rule for continuous incidence.
```

**Proposed**

```latex
Mean ranks agree with it in the rank-scaled Poisson configurations and can order regions differently when incidence is continuous. They are a deliberate adaptation for the case study.
```

### NY-11. Unequal treatment budgets

**Recommendation:** Accept. Also moves a mid-sentence "therefore" to the start of its sentence.

**Current**

```latex
Treatment budgets (numbers of treated clusters) therefore differ across
designs. We report them in Table \ref{tab:nc-budgets} rather than modify
the designs to force equal budgets.
```

**Proposed**

```latex
Therefore, treatment budgets (numbers of treated clusters) differ across designs. We kept each design's own rule and report its budget in Table \ref{tab:nc-budgets}.
```

### NY-12. Target of the case-study analysis

**Recommendation:** Accept. Touches HW-7 wording approved in 2g.

**Current**

```latex
The target is the direct coefficient \(\beta\)
before spatial propagation, not an average total effect after applying
\((\mathbf{I}-\rho \mathbf{W})^{-1}\).
```

**Proposed**

```latex
The target is the direct coefficient \(\beta\), the effect of a college's own treatment before outcomes spread through \((\mathbf{I}-\rho \mathbf{W})^{-1}\).
```

### NY-13. The hypothetical incidence association

**Recommendation:** Accept. Touches HW-8 wording approved in 2g.

**Current**

```latex
The association it builds in is hypothetical, not an
estimated relationship between incidence and education.
```

**Proposed**

```latex
The association it builds in is hypothetical.
```

### NY-14. Weights of repeated allocations

**Recommendation:** Accept. Also moves a mid-sentence "therefore" to the start of its sentence.

**Current**

```latex
The variance due to outcome noise therefore involves squared weights \((n_a/J)^2\), not \(n_a/J^2\) as if
duplicate estimates were independent.
```

**Proposed**

```latex
Therefore, the variance due to outcome noise involves the squared weights \((n_a/J)^2\), because repeated draws of one allocation share one risk estimate.
```

### NY-15. What the Jaccard index describes

**Recommendation:** Accept. Touches HW-9 wording approved in 2g.

**Current**

```latex
It describes the stability of sampled
tail membership, not a confidence interval or a population tail
probability.
```

**Proposed**

```latex
It describes how stable sampled tail membership is.
```

### NY-16. Case-study rank-test p-values

**Recommendation:** Accept. Touches HW-11 wording approved in 2g. The sentence before it says the years share one map and some simulations.

**Current**

```latex
Therefore, these $p$-values describe how consistently the designs ranked, not independent evidence.
```

**Proposed**

```latex
Therefore, these $p$-values describe how consistently the designs ranked across blocks that share data.
```

### NY-17. Flagged rook Checkerboard estimates (grid table text)

**Recommendation:** Accept.

**Current**

```latex
The table flags Checkerboard's rook-contiguity estimates, which target $\beta-\psi$ rather than $\beta$
```

**Proposed**

```latex
The table flags Checkerboard's rook-contiguity estimates, which target $\beta-\psi$
```

### NY-18. Source of the control-only increase

**Recommendation:** Accept. The next sentence explains why the oracle model has no omitted-variable bias.

**Current**

```latex
The increase under the control-only regime was variance, not bias.
```

**Proposed**

```latex
The control-only increase came almost entirely from variance.
```

### NY-19. Shared draws across effect sizes

**Recommendation:** Accept. Touches HW-13 wording approved in 2g; also removes its ", so" join.

**Current**

```latex
The five values of
\(\beta\) share their assignments and noise (Section
\ref{simulation-design-resampling-and-metrics}), so the stable order is
one finding on the same draws, not five independent confirmations.
```

**Proposed**

```latex
The five values of \(\beta\) share their assignments and noise (Section \ref{simulation-design-resampling-and-metrics}). Therefore, the stable order is one finding on shared draws.
```

### NY-20. Rook Checkerboard estimates (rook section)

**Recommendation:** Accept. Touches HW-16 wording approved in 2g.

**Current**

```latex
Under rook contiguity its estimates were precise but centered on
\(\beta-\psi\) rather than \(\beta\)
```

**Proposed**

```latex
Under rook contiguity its estimates were precise but centered on \(\beta-\psi\)
```

### NY-21. High Incidence Focus and randomization

**Recommendation:** Accept.

**Current**

```latex
Its low MSE reflects model-based
accuracy for that one allocation, not randomization-based causal
validity.
```

**Proposed**

```latex
Its low MSE describes model-based accuracy for that one allocation. With only one possible allocation, its validity rests on the outcome model alone.
```

### NY-22. Tail summaries and the reused assignments

**Recommendation:** Accept. Touches HW-24 wording approved in 2g.

**Current**

```latex
These summaries describe the spread of conditional MSE
across sampled assignments, not quantiles of realized estimation errors.
The refinement reused the main run's assignments rather than drawing a
fresh sample.
```

**Proposed**

```latex
These summaries describe the spread of conditional MSE across sampled assignments. The refinement reused the main run's assignments.
```

### NY-23. Prespecified mean-rank rule

**Recommendation:** Accept.

**Current**

```latex
We retain the
prespecified mean-rank rule rather than select a summary
retrospectively.
```

**Proposed**

```latex
We retain the prespecified mean-rank rule.
```

### NY-24. Matched sensitivity ratios

**Recommendation:** Accept. The sentence before it lists which settings each sensitivity mean uses.

**Current**

```latex
Each ratio is a matched contrast, not a ratio of
unmatched pooled analyses.
```

**Proposed**

```latex
Each ratio compares the sensitivity analysis with the primary analysis on the same settings.
```

### NY-25. Case-study formal tests

**Recommendation:** Accept. Touches HW-30 wording approved in 2g; also removes its ", so" join.

**Current**

```latex
(Section \ref{further-sensitivity-analyses}), so they describe the consistency of the ordering rather than independent evidence.
```

**Proposed**

```latex
(Section \ref{further-sensitivity-analyses}). Therefore, they describe how consistent the ordering is.
```

### NY-26. Mid-sentence therefore (Wilcoxon tests)

**Recommendation:** Accept. Style rule only.

**Current**

```latex
We
therefore also ran paired Wilcoxon signed-rank tests
```

**Proposed**

```latex
Therefore, we also ran paired Wilcoxon signed-rank tests
```

### NY-27. Mid-sentence therefore (Isolation Buffer regimes)

**Recommendation:** Accept. Style rule only.

**Current**

```latex
The two regimes in Equation \ref{eq:spillover} therefore define the same
outcome model
```

**Proposed**

```latex
Therefore, the two regimes in Equation \ref{eq:spillover} define the same outcome model
```

### NY-28. Mid-sentence therefore (powers of W)

**Recommendation:** Accept. Style rule only. The sentence prints only in the unabridged build if PC-8 is accepted.

**Current**

```latex
The powers of \(\mathbf{W}\) applied
to \(\mathbf{Z}\) therefore alternate
```

**Proposed**

```latex
Therefore, the powers of \(\mathbf{W}\) applied to \(\mathbf{Z}\) alternate
```

### NY-29. Mid-sentence therefore (error in the spatial coefficient)

**Recommendation:** Accept. Style rule only.

**Current**

```latex
The error in
\(\hat\rho\) is therefore multiplied by
```

**Proposed**

```latex
Therefore, the error in \(\hat\rho\) is multiplied by
```

## Batch SC · List titles (SI)

### SC-1. Short titles in the List of Figures and List of Tables

**Recommendation:** Accept. Each caption gets a short title for the lists only (`\caption[short]{full}`); the caption under each exhibit is unchanged. The lists stay in both builds (2e decision E2). Together with SC-2: 49 to 46 pages.

| Label | Short title |
|---|---|
| `tab:designs` | The six treatment assignment designs |
| `fig:design-samples` | Example grid assignment from each design and SRS |
| `tab:excluded-designs` | Saturation Quadrants and Balanced Halves, the two related variants |
| `fig:incidence-rejection` | Rejection rate of the interval for the incidence coefficient $\tau$ |
| `tab:parameter-grid` | Parameter grid of the grid simulation |
| `tab:seed-keys` | Seed keys |
| `fig:service-areas` | Counties assigned to the 58 community college service areas |
| `fig:nc-observed` | Observed annual SUD incidence by service-area cluster |
| `tab:nc-budgets` | Treatment budgets in the case study |
| `tab:nc-order` | Regional ordering under each incidence summary |
| `fig:nc-regions` | Fixed four-region partition for the saturation designs |
| `fig:nc-examples` | Example case-study assignment from each design |
| `tab:performance` | Grid performance of the six designs at $\beta = 1$ |
| `tab:performance` | Grid performance of the eight designs and SRS at $\beta = 1$ |
| `fig:mse` | Grid MSE by design under queen contiguity |
| `fig:coverage` | Grid interval coverage by design under queen contiguity |
| `tab:by-config` | MSE by design within each incidence configuration |
| `tab:rho` | MSE by design across levels of $\rho$ |
| `tab:psi` | MSE by design across levels of $\psi$ |
| `tab:regime` | MSE by design under each spillover regime |
| `fig:biasvar` | Squared bias and variance by design under queen contiguity |
| `tab:eight-design` | Grid performance of the eight designs and SRS at $\beta = 1$ |
| `fig:cd` | Nemenyi critical-difference diagrams |
| `fig:tau-mse` | MSE by design across effect sizes $\beta$ |
| `fig:tau-coverage` | Interval coverage by design across effect sizes $\beta$ |
| `tab:checkerboard-coverage` | Checkerboard coverage across levels of $\rho$ |
| `fig:mse-rook` | Grid MSE by design under rook contiguity |
| `fig:coverage-rook` | Grid interval coverage by design under rook contiguity |
| `fig:biasvar-rook` | Squared bias and variance by design under rook contiguity |
| `tab:oracle-nonoracle` | Oracle and non-oracle estimates of $\beta$ |
| `tab:robustness` | Distribution of block-level MSE by design |
| `fig:vs-srs` | Design MSE relative to SRS by spillover regime |
| `tab:grid-srs` | Full grid benchmark against SRS |
| `fig:igs-sensitivity` | MSE of the variants of incidence guidance relative to SRS |
| `fig:nc-mse` | Case-study yearly MSE relative to SRS |
| `tab:nc-year` | Case-study yearly mean MSE |
| `fig:nc-pop` | Population shares treated in the case study |
| `fig:nc-setting-ratios` | Case-study MSE ratios to SRS by setting |
| `fig:nc-mse-detail` | Case-study MSE by design, year and regime |
| `tab:nc-annual` | Case-study performance by design and planning year |
| `fig:spread-gain` | Gain in treated incidence share against regional spread |
| `tab:downside` | Mean and worst-decile conditional MSE relative to SRS |
| `fig:alloc-dist` | Distribution of allocation-specific MSE relative to SRS |
| `fig:nc-risk` | Allocation-specific conditional MSE in the case study |
| `tab:nc-summary-sens` | Incidence-Guided Saturation sensitivity analyses |
| `fig:colleges-needed` | Colleges needed for 80\% power by design |
| `fig:power-n` | Power by number of participating colleges |
| `tab:nc-power` | Colleges needed for 80\% power and minimum detectable effect |
| `tab:nc-translation` | Power for an individual-level standardized effect |
| `tab:nc-calibration` | Size and coverage of the Wald and corrected $t$ tests |
| `tab:interval-comparison` | Coverage of the Wald and corrected $t$ intervals |
| `tab:nc-guidance` | MSE relative to SRS by number of participating colleges |
| `tab:nc-car` | MSE relative to SRS under spatial lag and CAR models |
| `tab:nc-bayes` | MSE relative to SRS under maximum likelihood and Bayesian fits |
| `tab:nc-forms` | MSE relative to SRS under three spillover forms |

The four subfigure captions ("Queen contiguity (primary)", "Rook contiguity (sensitivity)") are already short and stay as they are.

### SC-2. List of Figures and List of Tables on one page

**Recommendation:** Accept. With short titles both lists fit on one page.

**Current**

```latex
\listoffigures
\clearpage
```

**Proposed**

```latex
\listoffigures
\bigskip
```

## Batch PC · Light trim of the curated SI

Each item moves text into the unabridged build by wrapping it in `\iffullsi ... \fi`. A moved block is named by its first and last words. Pointer sentences that would otherwise print a missing number change as shown. The CTJ cites none of the moved labels.

### PC-2. Table S1 (six design rules) moves

**Recommendation:** Accept. Table S1 (six design rules) moves; the article's Table 1 gives the same rules with example maps.


**Current (edit 1)**

```latex
Table \ref{tab:designs} describes the six retained designs.
\iffullsi Figure \ref{fig:design-samples}\else
Table~\ref{M-tab:allocation-rules} of the article\fi\ shows one assignment
```

**Proposed (edit 1)**

```latex
\iffullsi Table \ref{tab:designs} describes the six retained designs. Figure \ref{fig:design-samples}\else
Table~\ref{M-tab:allocation-rules} of the article describes the six retained designs and\fi\ shows one assignment
```

- **Moves to the unabridged build:** from "\begin{table}[!htb] \small \centering \caption{The six treatment assignment designs" to "\end{tabular} \end{table} \iffullsi \begin{figure}".

**Current (edit 2)**

```latex
\caption{The two related variants outside the detailed six-design grid ranking, in the format of Table \ref{tab:designs}. Saturation Quadrants shares its mechanism with Incidence-Guided Saturation, and Balanced Halves with Balanced Quartiles. Incidence rank and clusters treated are defined as in Table \ref{tab:designs}.}
```

**Proposed (edit 2)**

```latex
\iffullsi\caption{The two related variants outside the detailed six-design grid ranking, in the format of Table \ref{tab:designs}. Saturation Quadrants shares its mechanism with Incidence-Guided Saturation, and Balanced Halves with Balanced Quartiles. Incidence rank and clusters treated are defined as in Table \ref{tab:designs}.}
\else\caption{The two related variants outside the detailed six-design grid ranking. Saturation Quadrants shares its mechanism with Incidence-Guided Saturation, and Balanced Halves with Balanced Quartiles. Incidence rank is a cluster's rank by the baseline incidence surface $\mathbf{X}_k$, with ties broken at random in every draw.}
\fi
```


**Current (edit 3)**

```latex
Balanced Halves in the format of Table \ref{tab:designs}. Each shares
```

**Proposed (edit 3)**

```latex
Balanced Halves\iffullsi\ in the format of Table \ref{tab:designs}\fi. Each shares
```


**Current (edit 4)**

```latex
Design & The six designs of Table \ref{tab:designs} and the two of Table \ref{tab:excluded-designs} & 8 \\
```

**Proposed (edit 4)**

```latex
Design & The six designs of \iffullsi Table \ref{tab:designs}\else Table~\ref{M-tab:allocation-rules} of the article\fi{} and the two of Table \ref{tab:excluded-designs} & 8 \\
```


### PC-4. Section 5.4 second paragraph (refinement schedule and frequency-weight bookkeeping) moves

**Recommendation:** Accept. Section 5.4 second paragraph (refinement schedule and frequency-weight bookkeeping) moves; Section 5.3 states the precision targets.

- **Moves to the unabridged build:** from "For randomized designs, the prespecified refinement schedule allowed" to "from one with only one possible assignment.".

### PC-7. Section 10 rook rank tests move

**Recommendation:** Accept. Section 10 rook rank tests move; Table S6 and Section 12.2 keep the rook results and the comparison with SRS.

- **Moves to the unabridged build:** from "With SRS as a seventh member, the Friedman test rejected equal ranks (\(\chi^2 = 795.77\)" to "including SRS against \(2\times2\) Blocking, differed in all 50 units.".

### PC-8. Section 9 last paragraph (why rook Checkerboard coverage rises with β) moves

**Recommendation:** Accept. Section 9 last paragraph (why rook Checkerboard coverage rises with β) moves; the coverage values stay in Section 8.


**Current (edit 1)**

```latex
Section \ref{stability-across-treatment-effect-sizes} reports that
Incidence-Guided Saturation's bias and rook Checkerboard's coverage grow
with \(\beta\).
```

**Proposed (edit 1)**

```latex
Section \ref{stability-across-treatment-effect-sizes} reports that
Incidence-Guided Saturation's bias\iffullsi\ and rook Checkerboard's coverage grow\else\ grows\fi\
with \(\beta\).
```


**Current (edit 2)**

```latex
re-run fit, with a largest residual of \(4.2\times10^{-13}\); for
Checkerboard under rook contiguity it held in the reduced form given
below.
```

**Proposed (edit 2)**

```latex
re-run fit, with a largest residual of \(4.2\times10^{-13}\)\iffullsi; for
Checkerboard under rook contiguity it held in the reduced form given
below\else. For Checkerboard under rook contiguity, whose spillover column is aliased, it held with that column dropped\fi.
```

- **Moves to the unabridged build:** from "The decomposition also explains why the coverage of Checkerboard under rook" to "\(-0.667\).".

**Current (edit 3)**

```latex
0.286, and 0.464 at \(\beta = 1\), 1.5, 2, and 3 (Section
\ref{app:technical-notes}).
```

**Proposed (edit 3)**

```latex
0.286, and 0.464 at \(\beta = 1\), 1.5, 2, and 3\iffullsi\ (Section
\ref{app:technical-notes})\fi.
```


### PC-12. Section 6.2, two sentences that the article's Figure 2 shows

**Recommendation:** Accept. Section 6.2, two sentences that the article's Figure 2 shows.


**Current (edit 1)**

```latex
The configuration changed each design's MSE only slightly. In this
order, each design's highest MSE across the five configurations was
lower than the next design's lowest. Incidence-Guided Saturation had the lowest MSE in 18 to 21 of the 32 blocks in each
configuration. Coverage
```

**Proposed (edit 1)**

```latex
The configuration changed each design's MSE only slightly.\iffullsi\ In this
order, each design's highest MSE across the five configurations was
lower than the next design's lowest. Incidence-Guided Saturation had the lowest MSE in 18 to 21 of the 32 blocks in each
configuration.\fi\ Coverage
```


### PC-13. Section 11, the bias and coverage ranges and a closing sentence that Table S7 prints

**Recommendation:** Accept. Section 11, the bias and coverage ranges and a closing sentence that Table S7 prints.


**Current (edit 1)**

```latex
coverage from 0.54 to 0.91 across designs and neighbor types. Under
queen contiguity its bias ranged
```

**Proposed (edit 1)**

```latex
coverage from 0.54 to 0.91 across designs and neighbor types.\iffullsi\ Under
queen contiguity its bias ranged
```


**Current (edit 2)**

```latex
Checkerboard was Isolation Buffer (bias \(-0.400\), coverage 0.536). For
```

**Proposed (edit 2)**

```latex
Checkerboard was Isolation Buffer (bias \(-0.400\), coverage 0.536).\fi\ For
```

- **Moves to the unabridged build:** from "Table \ref{tab:oracle-nonoracle} gives the mean squared error, bias, and coverage of the oracle and" to "as summarized above.".

### PC-14. Section 12.2 last paragraph

**Recommendation:** Accept. Section 12.2 last paragraph; Table S9 and Sections 12.3 and 15 give its worst-decile ratios.

- **Moves to the unabridged build:** from "The mean is only part of the comparison, because a trial realizes a single allocation." to "gives the precision of the case study tail estimates.".

### PC-15. Section 8 second paragraph (six-design ranks and MSE by β)

**Recommendation:** Accept (author). Section 8 second paragraph (six-design ranks and MSE by β); Figure S4 shows the curves, and the CTJ's β claim rests on the paragraph that stays.

- **Moves to the unabridged build:** from "The gaps between designs changed as the treatment effect increased." to "between 0.037 and 0.041.".

**Current (edit 1)**

```latex
Isolation Buffer gained under the
control-only regime, where Incidence-Guided Saturation had a
lower MSE
```

**Proposed (edit 1)**

```latex
\iffullsi Isolation Buffer gained under the
control-only regime, where\else Under the control-only regime,\fi\ Incidence-Guided Saturation had a
lower MSE
```


### PC-11. Case-study Friedman/Nemenyi tests (Section 5.6 first paragraph and the 16.4 block)

**Recommendation:** Accept (author). Case-study Friedman/Nemenyi tests (Section 5.6 first paragraph and the 16.4 block); Tables S8, S14 and S15 show the ordering.

- **Moves to the unabridged build:** from "We compared the eight designs and SRS with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid," to "across blocks that share data. To check that".

**Current (edit 1)**

```latex
whether $\beta$ was identifiable depended on the assignment.
Formal rank tests summarize the ordering of the eight designs and SRS
within each regime.
```

**Proposed (edit 1)**

```latex
whether $\beta$ was identifiable depended on the assignment.\iffullsi\
Formal rank tests summarize the ordering of the eight designs and SRS
within each regime.\fi
```

- **Moves to the unabridged build:** from "\sihead{Formal tests on the case study results}" to "\iffullsi \subsection{Implications for trial planning}".

## Batch E (NYC) · Article rewrites (CTJ_Manuscript.tex)

These change the article body, so the title-page word count (3,431) is recomputed when they are applied. Content edits to the CTJ usually need the same change in Chapter 3; they join the 2i list.

### NYC-1. Abstract scope sentence

**Recommendation:** Accept.

**Current**

```latex
This study evaluates trial designs, not the effectiveness of an implemented intervention.
```

**Proposed**

```latex
This study uses simulation to evaluate designs for a planned trial.
```

### NYC-2. Direct effect definition

**Recommendation:** Accept.

**Current**

```latex
before outcomes spread between neighbors through $(\mathbf{I}-\rho \mathbf{W})^{-1}$, rather than a total effect.
```

**Proposed**

```latex
before outcomes spread between neighbors through $(\mathbf{I}-\rho \mathbf{W})^{-1}$.
```

### NYC-3. Units for Monte Carlo standard errors

**Recommendation:** Accept.

**Current**

```latex
(five configurations by ten surfaces), rather than treating settings that share a surface as independent replications.
```

**Proposed**

```latex
(five configurations by ten surfaces), since results from settings that share a surface are correlated.
```

### NYC-4. Planning years

**Recommendation:** Accept.

**Current**

```latex
The four years were analyzed as separate planning scenarios rather than as follow-up of an intervention.
```

**Proposed**

```latex
The four years were analyzed as separate planning scenarios.
```

### NYC-5. Checkerboard and Isolation Buffer accuracy

**Recommendation:** Accept.

**Current**

```latex
Their worse control-only accuracy reflects imprecision rather than bias
```

**Proposed**

```latex
Their worse control-only accuracy comes almost entirely from variance
```

### NYC-6. Conclusions scope sentence

**Recommendation:** Accept.

**Current**

```latex
The education case study is a design evaluation, not an implemented intervention or evidence that education reduces SUD.
```

**Proposed**

```latex
The education case study evaluates designs for a planned trial of education for healthcare professionals. Its outcomes are simulated and give no evidence about whether education reduces SUD.
```

## Chapter 3 follow-ups for 2i

- NYC-1 to NYC-6 and the NY items with chapter counterparts: carry the same rewrites to `Thesis_Chapter.qmd`.
- SC and PC items are SI-only.

## Independent review (2026-11-25)

A fresh Opus reviewer checked every item against the sources. Its fixes are already in the items above: NY-5 keeps "included only to show", NY-8 puts the cause first, NY-10 keeps that mean ranks can order regions differently for continuous incidence, NY-21 says the design's validity rests on the outcome model alone, NY-18 and NYC-5 say "almost entirely from variance" (squared bias is at most 2.2% of MSE), NYC-1 no longer implies no intervention has run (a pilot has), NYC-6 keeps that the simulated outcomes give no evidence about SUD, NYC-3 says the results (not the settings) are correlated, PC-4 keeps the curated sentence on how the MCSEs of setting-specific differences are computed (the article relies on it), PC-8 keeps a self-contained rook Checkerboard qualifier, and three short titles are tightened. It also found four remaining mid-sentence "therefore"s, now NY-26 to NY-29. Apply order: NY items first, then PC (PC-11 is anchored on NY-16's new wording), then SC. Noted, not changed: PC-2 removes the only curated statement of Balanced Quartiles' 12/12/13/13 split (the article's Table 1 says half of each quartile).

## Author decisions

- **PC-11 and PC-15 (2026-11-25):** the author asked to include both cuts.
- Everything else pending.

