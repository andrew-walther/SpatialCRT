# Phase 2g merged sign-off: SI prose, statistics, new items and article

Author: Andrew Walther. Drafted 2026-11-25.

This file merges the Codex proposals ([2g-signoff.md](2g-signoff.md), HW-1 to HW-30, and [2g-review-signoff.md](2g-review-signoff.md), SR-1 to SR-18) with four independent reviews: a biostatistical read, a trial-planner read, a human-writer read of the SI and a human-writer read of the article. It replaces those two files for approval. No manuscript file has been changed.

How to read an item. The verdict line gives each reader's call. Where any reader revised a Codex item, the proposed block below is the single merged wording, which restores each dropped fact or caveat and keeps Codex's prose gains. Where every reader confirmed, Codex's blocks appear unchanged. An SR item that builds on an HW item quotes the final HW wording in its current block, so decide both IDs together. SN items are new SI items, and HC items concern `CTJ_Manuscript.tex`. Every current block was checked by script (see the end of the file). Approve by ID.

## Summary

| ID | Recommendation | Reason |
|---|---|---|
| HW-17 | Accept revised | Restores the section's main finding, that no single design was best under both spillover regimes, and the regime-specific guidance. |
| HW-3 | Accept revised | Restores which model is primary and which is a sensitivity analysis, and that the oracle assumes a known spillover form. |
| HW-7 | Accept revised | Restores "but not X", the direct-versus-total distinction and the one-sided grid power at level 0.025, and adds the planner's plain reading of the two regimes and of \(\psi\). |
| HW-23 | Accept revised | Restores the regions-differ-little framing, the rate-versus-death-share "because" and the grid comparison, and moves the High Incidence Focus share from the caption back to the text with its unit fixed. |
| HW-25 | Accept revised | Restores that the added outcomes did not enlarge the sample of assignments, the reason each worst decile stays at ten draws. |
| HW-13 | Accept revised | Restores that the five effect sizes on shared draws give one finding, not five confirmations, and makes the caption say what rook Checkerboard estimates. |
| HW-11 | Accept revised | Block 1 restores "not independent evidence". |
| HW-30 | Accept revised | Block 3 replaces "within a few percent" with the verified "within 8%" (C-6). |
| HW-6 | Accept revised (prose), reject table block | The prose states round half to even plainly. |
| SR-7 | Accept revised | Restores the takeaway that the design hardly affected the number of colleges, labels the example as an illustration and states which number of colleges each control-only power belongs to. |
| HW-24 | Accept revised | Restores that the tail summaries are not quantiles of realized errors and not a fresh sample of assignments. |
| HW-9 | Accept revised | Restores the Jaccard caveat (not a confidence interval or tail probability) and fixes two unclear sentences. |
| HW-20 | Accept revised | The caption restores that worst deciles are selected and evaluated with the same outcomes. |
| HW-22 | Accept revised | Block 1 keeps the claim that the shared saturation structure drives the advantage. |
| HW-16 | Accept revised | Block 1 says what rook estimates target and keeps the queen coverage contrast. |
| HW-14 | Accept revised | Labels the average ranks as ranks among the six designs without SRS, which differ from Table S6, and splits the opening sentence. |
| HW-15 | Accept revised | Every term of the decomposition, not only those involving \(\hat\rho\), comes from the rerun. |
| SR-1 | Accept revised | Block 1 names the comparisons that did not hold (CAR fits to spatial lag data). |
| SR-2 | Accept revised | Block 1 says a one-unit change spans the whole 0-to-1 scale. |
| SR-6 | Accept revised | Says that averaging divides the variance, states the assumption behind \(d=\beta\sigma_c\) and writes \(d\) in one plain sentence. |
| SR-10 | Accept revised | Restores "classical contamination" and that neighboring treatment adds to a treated cluster's own effect. |
| SR-8 | Accept revised | Restores the run-in heading that the replacement dropped and that the other paragraphs of this section use. |
| HW-1 | Accept revised | Restores that the earlier grids used rook neighbors and that 35 per 100,000 is not the case-study rate, and uses the SI's name for the earlier study. |
| HW-2 | Accept revised | Keeps "the only deterministic design", which the next sentences rely on, and the SI's name "the earlier small-grid study". |
| HW-4 | Accept revised | Keeps the technical term "common random numbers" and fixes the tense. |
| HW-8 | Accept revised | Says which Monte Carlo standard error each target applies to and restores "not independent confirmations". |
| HW-10 | Accept revised | Restores "only" (nonparticipants do not enter the spatial dependence), uses "that" for a college and removes the colon before the display. |
| HW-18 | Accept revised | Block 1 keeps the subject of the equivalence caveat and orders the two caveats. |
| HW-12 | Accept revised | Uses "unit", as the paragraph does, and says that the losses outweigh the wins in the pooled mean. |
| HW-5 | Accept as proposed | Removes a repeated access statement. |
| HW-19 | Accept as proposed | Keeps the numbers, the extreme-setting scope and the noise share. |
| HW-21 | Accept as proposed | Removes pipeline narration. 7,598,000 plus 2,280,000 equals 9,878,000, so "giving" is exact. |
| HW-26 | Accept as proposed | Clearer row label. |
| HW-27 | Accept as proposed | Removes the manifest jargon and leads with the conditions. |
| HW-28 | Accept as proposed | Removes "narrowly" and separates the two runs. |
| HW-29 | Accept as proposed | Three paragraphs, with the diagnostics kept. |
| SR-3 | Accept as proposed | With \(\rho\) fixed at its true value, \(\hat\beta\) does not depend on \(\tau\). |
| SR-4 | Accept as proposed | States the two conditional means correctly. |
| SR-5 | Accept as proposed | Replaces the sign argument with the product \((\rho-\hat\rho)g\). |
| SR-9 | Accept as proposed | Observed rejection rates of 4.4% to 5.6% support "close to nominal". |
| SR-11 | Accept as proposed | Scope verified. 95 + 58 + 7 = 160, so SRS is excluded. |
| SR-12 | Accept as proposed | The non-oracle order differs, so dropping "do not change the comparison" is right. |
| SR-13 | Accept as proposed | Describes the pooled standard error and the approximate ratio interval correctly. |
| SR-14 | Accept as proposed | Matches the downside table and removes the bare "NC". |
| SR-15 | Accept as proposed | Its sentence is applied inside the final HW-23 caption, so it has no separate edit. |
| SR-16 | Accept as proposed | Wording unchanged. |
| SR-17 | Accept as proposed | Accurate diagnostic language. |
| SR-18 | Accept as proposed | Two-decimal ratios follow the author's rule. |
| SN-1 | Accept as proposed | The 2.700 for Incidence-Guided Saturation is a rank among the eight designs and SRS, while Table S6 prints 1.775 among six designs and SRS. |
| SN-2 | Accept revised | Each setting used 500 fits, and the per-setting MSE is far less precise than in the main case study, which explains why, for example, High Incidence Focus reads 1.06 here against about 0.94 in the main analysis. |
| SN-3 | Accept as proposed | The original design's control-only ratio is 0.67 in Section 13 and 0.631 in Section 12.2. |
| SN-4 | Accept as proposed | In the curated build these two sentences follow no figure and repeat the 0.162 coverage and the squared-bias share already reported above them. |
| SN-5 | Accept as proposed | The opening paragraph uses "queen and rook results" before either term is defined, and "labeled separately" is pipeline language. |
| SN-6 | Accept as proposed | "Education ability" is not the study's outcome. |
| SN-7 | Accept as proposed | The term is used twice in Section 16 and never defined. |
| SN-8 | Accept revised | The first case-study MSE values (0.223 against 0.308) appear with no unit. |
| SN-9 | Accept as proposed | Author rule against a bare "NC" descriptor. |
| SN-10 | Accept as proposed | Author rule puts "However" at the start of the sentence, and "which matters because" announces rather than states. |
| SN-11 | Accept revised | As written, "significant" could attach to either regime. |
| SN-12 | Accept as proposed | Approve with HW-26. "Baseline-adjusted" would be a third name for the analysis that HW-26 and SR-14 call the incidence sensitivity, and the first sentence repeats the pointer that ends the paragraph before it. |
| SN-13 | Accept as proposed | The caveat matters to planners, and "That is model-based accuracy" has no clear antecedent. |
| SN-14 | Accept as proposed | HW-7 now points readers to this definition, and the 70-word sentence hides the values of \(p\) in a parenthesis split by a semicolon. |
| SN-15 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-16 | Accept as proposed | Semicolon rule. "Are noted here" is dropped because the sentence itself is the note. |
| SN-17 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-18 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-19 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-20 | Accept as proposed | Semicolon rule. |
| SN-21 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-22 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-23 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| SN-24 | Accept as proposed | Semicolon rule. |
| SN-25 | Accept as proposed | Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. |
| HC-1 | Accept as proposed | Abbreviation used before it is defined. "SI" first appears in Methods (outcome model paragraph) but is defined only later in "Treatment assignment designs". |
| HC-2 | Accept as proposed | Undefined terms. "Rook" and "queen" contiguity are used in the identifiability paragraph before the Simulation study defines them, and a clinical reader may not know them. |
| HC-3 | Accept as proposed | Pivot sentence (AI trope). "The two spillover regimes gave opposite results" announces a contrast that the next two sentences and the subsection heading already show. |
| HC-4 | Accept as proposed | Author rule: never a bare "NC" as a descriptor. "On the NC map" is the only body use of the abbreviation, so the definition in the case-study paragraph can go too (Table 1 caption defines NC on its own). |
| HC-5 | Accept as proposed | Agentless passive with the verb repeated, and MSE is defined a second time (it was defined in the Introduction). |
| HC-6 | Accept as proposed | Sentence of about 50 words (rank-scaling of Poisson counts). |
| HC-7 | Decline | Replaces a sentence-initial "Therefore" with "Because" or ", so", against the author's style. |
| HC-8 | Decline | Replaces a sentence-initial "Therefore" with "Because" or ", so", against the author's style. |
| HC-9 | Decline | Replaces a sentence-initial "Therefore" with "Because" or ", so", against the author's style. |
| HC-10 | Decline | Replaces a sentence-initial "Therefore" with "Because" or ", so", against the author's style. |
| HC-11 | Accept as proposed | "In this paper," opens two consecutive Introduction paragraphs. |
| HC-12 | Accept as proposed | Two unrelated facts yoked by "while" (coverage across designs, then Checkerboard's rook bias), and "Checkerboard's was" leaves the noun implied. |
| HC-13 | Accept as proposed (optional) | Repeated caveat. |
| HC-14 | Accept as proposed | Repeated caveat. |
| HC-15 | Accept as proposed (optional) | Redundant gloss. "Elementwise product" already means multiplying cluster by cluster. |
| HC-16 | Accept as proposed | Stacked "and ... and" in a fronted phrase delays the subject. |
| HC-17 | Accept as proposed | Sentence of 33 words joining two separate points with "and". |
| HC-18 | Accept as proposed | Abbreviation defined but never used again in the body. |
| HC-19 | Accept as proposed | Abstract (not counted in the body limit). |
| HC-20 | Accept as proposed | The title page says 3,435 but the body count is 3,473 now. |

## Conflicts resolved

- **HW-6 table (Codex vs Writer).** Rejected. `verify_manuscripts.py` formats the budget column of Chapter 3's `tab:nc-budgets` with two decimals and requires the SI table body to equal the chapter's, so "27" would fail both checks. The prose half of HW-6 is accepted in revised form.
- **HW-26 row labels (Writer condition).** No chapter change is needed. `tab:nc-summary-sens` is not in Chapter 3 (held out 2026-11-15), and the verifier checks only the SI table's numbers.
- **HW-23 High Incidence Focus share (Codex, Biostat, Planner).** Codex moved "0.75 to 0.77" into the spread-gain caption. Biostat and the planner objected because the figure does not plot it and the unit did not match. It returns to the text as "75% to 77% of total incidence", checked against the saved burden shares.
- **HW-13 caption (Planner vs Writer).** The planner's "plotted MSE is not that of an estimate of \(\beta\)" and the writer's "plotted estimates target \(\beta-\psi\)" are merged into one sentence that names what the coefficient estimates.
- **HW-7 regime description (Planner only).** The planner's sentences translating the two regimes and \(\psi\) into trial terms are added. Biostat and the writer did not object to the content, and the planner flagged the gap as the main one for a trial reader.
- **HW-22 block 1 (Writer only).** The writer restored "their advantage over SRS comes from the saturation structure they share". Section 2.1 says the same, but this is where the case-study result is reported, so it stays.
- **SN-8 (Planner).** The units sentence is added. The planner's deletion of "descriptive" is not, because that word marks results reported without pooled Monte Carlo intervals.
- **SR-8 (all confirmed).** Accepted with the run-in heading restored, which the replacement had dropped.
- **Spill caption (C-5, N6, HN-1).** Merged into HW-3 as its third block instead of a separate item, so HW-3 cannot be approved without it.
- **C-6 (Biostat).** Merged into HW-30 block 3.

## Numbers verified against saved results

- **SR-7.** In `application/results/methods_ext_2026-10/summaries_v2/a3_translation_m_icc.csv`, Incidence-Guided Saturation has control-only power 0.508 at 29 colleges and 0.641 at 58. SRS has 0.488 and 0.525. The planner's assumed order is correct.
- **SN-11 (HN-5).** In `results/srs_benchmark/srs_benchmark_summary.txt` (oracle, rook, \(\beta=1\)), the both-arms Holm-adjusted p-values are 0.0075 for Incidence-Guided Saturation and 0.150 for Saturation Quadrants. Both control-only tests have p = 6.2e-09. "Significant" therefore refers to the both-arms tests.
- **HW-23.** In `results/grid_burden_share/grid_burden_share_summary.csv`, High Incidence Focus treats 0.748 to 0.766 of total grid incidence (queen, five configurations). That share belongs in the text, not in the caption of a figure that does not plot it.
- **SN-2 (C-2).** In `a6_car_settings.csv`, all 1,728 settings have 500 fits, and the Monte Carlo standard error of MSE is 4.3% to 11.1% of MSE.
- **HW-30 block 3 (C-6) and SN-7 (N3).** Checked against the printed `tab:nc-guidance` and `tab:nc-power` bodies.

## Batch A1 · Restored facts and resolved conflicts (highest value)

### HW-17. Robustness across parameter blocks

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores the section's main finding, that no single design was best under both spillover regimes, and the regime-specific guidance. All three reviewers found both dropped.  The session lead also restored the approved framing that guided saturation matches Saturation Quadrants while targeting high-incidence regions.

**Current**

```latex
No single design was best under both spillover regimes. Under
control-only spillover, Incidence-Guided Saturation had the lowest MSE
of the six designs, and it had the lowest worst-case block MSE\iffullsi\ (Table
\ref{tab:robustness})\fi. Under both-arms spillover, SRS and Balanced
Quartiles were slightly more accurate (Section
\ref{srs-benchmark-and-allocation-specific-risk}). A
saturation design is a reasonable choice when control-only spillover is
plausible, and SRS remains competitive under both-arms spillover. Saturation Quadrants performed similarly, so incidence guidance lets a trial direct its highest treatment saturations to the highest-incidence regions with little or no loss of estimation accuracy. \iffullsi The\else These\fi\ scenario-level worst cases\iffullsi\ in Table \ref{tab:robustness}\fi\ summarize variation across
parameter settings, which is different from the conditional risk across
assignments within one setting (Section \ref{allocation-risk-results}).
```

**Proposed**

```latex
No single design was best under both spillover regimes. Under
control-only spillover, Incidence-Guided Saturation had the lowest
mean and worst-case block MSE of the six designs\iffullsi\ (Table
\ref{tab:robustness})\fi. Under both-arms spillover, SRS and Balanced
Quartiles had slightly lower mean MSE (Section
\ref{srs-benchmark-and-allocation-specific-risk}). A saturation design
is a reasonable choice when control-only spillover is plausible, and SRS
remains competitive under both-arms spillover. Saturation Quadrants
performed similarly, so incidence guidance lets a trial direct its
highest treatment saturations to the highest-incidence regions with
little or no loss of estimation accuracy (Section \ref{app:variants}).
Worst-case block MSE describes variation across parameter settings.
Section \ref{allocation-risk-results} examines variation across
assignments within a setting.
```

### HW-3. Oracle model terminology and the identification pointer

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores which model is primary and which is a sensitivity analysis, and that the oracle assumes a known spillover form. Drops the unused "ML" and fixes the Table S7 caption, which would otherwise cite an undefined "Spill" (C-5, N6, HN-1).

The third block is added from C-5 (Biostat), N6 (Planner) and HN-1 (Writer), which all flagged it. Approve it with blocks 1 and 2. Block 2 drops the identity \(\mathbf{W}\mathbf{Z}=\mathbf{1}-\mathbf{Z}\) here because Section `identification-diagnostics` states it in full.

**Current (block 1 of 3)**

```latex
We estimated \(\beta\) with two maximum likelihood spatial lag models,
both fit to every simulated outcome and both including an intercept. The
primary, oracle model has the spatial lag \(\rho \mathbf{W}\mathbf{Y}\) and the covariates
\(\mathbf{Z}\), \(\text{Spill} = \mathbf{S}(\mathbf{Z})\), and \(\mathbf{X}_k\), so it matches the outcome
model in Equation \ref{eq:outcome}, including the true regime-specific
spillover covariate from Equation \ref{eq:spillover}. The oracle model
assumes that the analyst knows the form of the spillover mechanism and
has the exact exposure term, so its results are a best case for each
design; we return to this assumption as a limitation. The non-oracle
model, a sensitivity analysis, has the covariates \(\mathbf{Z}\) and \(\mathbf{X}_k\) only
and omits the spillover term. Both models are fit to the same outcomes.
```

**Proposed (block 1 of 3)**

```latex
We fit two maximum likelihood spatial lag models with intercepts to each
simulated outcome. The primary, oracle model includes the spatial lag
\(\rho \mathbf{W}\mathbf{Y}\) and the covariates \(\mathbf{Z}\),
\(\mathbf{S}(\mathbf{Z})\) and \(\mathbf{X}_k\). It matches the outcome
model in Equation \ref{eq:outcome}, using the true regime-specific
spillover term in Equation \ref{eq:spillover}. The oracle model assumes
that the analyst knows the form of spillover and has the exact exposure
term, so its results are a best case for each design. The non-oracle
model, a sensitivity analysis, includes \(\mathbf{Z}\) and
\(\mathbf{X}_k\) and omits the spillover term.
```

**Current (block 2 of 3)**

```latex
For Checkerboard, strict interspersion makes treatment and spillover hard or impossible to separate. Under rook contiguity on the rectangular grid, \(\mathbf{W}\mathbf{Z}=\mathbf{1}-\mathbf{Z}\), so
treatment and spillover are not separately identifiable; the fitted
coefficient on treatment estimates \(\beta-\psi\). Under queen
contiguity, treatment is identifiable but control-only spillover is weakly
separated from treatment. We retain Checkerboard in the main design
comparison and flag aliased estimates. Exact matrix identities, boundary
information and the rare Isolation Buffer aliases are detailed in
Section \ref{identification-diagnostics}. These diagnostics concern the models
implemented here.
```

**Proposed (block 2 of 3)**

```latex
We retain Checkerboard in the design comparison and flag aliased estimates.
Section \ref{identification-diagnostics} explains why rook Checkerboard
cannot separate treatment from spillover and why queen Checkerboard
estimates the treatment effect imprecisely under control-only spillover.
It also describes the rare Isolation Buffer aliases.
```

**Current (block 3 of 3)**

```latex
The non-oracle model omits $\text{Spill}$.
```

**Proposed (block 3 of 3)**

```latex
The non-oracle model omits the spillover term $\mathbf{S}(\mathbf{Z})$.
```

### HW-7. Case study interval definition

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores "but not X", the direct-versus-total distinction and the one-sided grid power at level 0.025, and adds the planner's plain reading of the two regimes and of \(\psi\).

The explicit interval formula and its 0.975 quantile move to Section `estimation` (see SN-14), which HW-7 now cites. The added regime sentences come from the planner's revision and use the case study's values \(\psi=0.5, 0.8\) with \(\beta=1\).

**Current**

```latex
For the primary education outcome, we set \(\tau=0\) in Equation
\ref{eq:outcome}, use the 58-cluster \(\mathbf{W}\), and fix \(\beta=1\),
\(\sigma=1\): \begin{equation}
\mathbf{Y}=(\mathbf{I}-\rho \mathbf{W})^{-1}[\beta \mathbf{Z}+\psi\,\mathbf{S}(\mathbf{Z})+\boldsymbol{\varepsilon}],\qquad \boldsymbol{\varepsilon}\sim N(\mathbf{0},\mathbf{I}).
\label{eq:nc-education}
\end{equation} Observed incidence affects assignments but neither
baseline education nor the constant structural intervention effect. We estimate $\beta$ with the same validated maximum likelihood implementation as in the grid study, with an intercept, \(\mathbf{Z}\)
and the true regime-specific \(\mathbf{S}(\mathbf{Z})\), without \(\mathbf{X}\). As in the grid
model, \(\beta\) is the coefficient before spatial propagation, not an
average total impact after applying \((\mathbf{I}-\rho \mathbf{W})^{-1}\). For the main case study and for the analyses of partial participation and sample size below, we report $\hat\beta \pm t_{0.975,\,N-p}\,\widehat{\mathrm{SE}}\sqrt{N/(N-p)}$, where $N$ is the number of participating colleges ($N = 58$ in the main case study) and $p$ the number of fitted coefficients plus $\rho$, and we test $\beta = 0$ two-sided at level 0.05, rejecting when the 95\% interval excludes zero. Grid power (Section \ref{simulation-design-resampling-and-metrics}) is one-sided at level 0.025. The maximum likelihood variance divides by $N$, and the normal reference distribution ignores the estimated parameters. With 58 or fewer clusters, both make the Wald interval too narrow.
```

**Proposed**

```latex
For the primary education outcome we set \(\tau=0\) in Equation
\ref{eq:outcome}, use the 58-cluster \(\mathbf{W}\), and fix \(\beta=1\)
and \(\sigma=1\).
\begin{equation}
\mathbf{Y}=(\mathbf{I}-\rho \mathbf{W})^{-1}[\beta \mathbf{Z}+\psi\,\mathbf{S}(\mathbf{Z})+\boldsymbol{\varepsilon}],\qquad \boldsymbol{\varepsilon}\sim N(\mathbf{0},\mathbf{I}).
\label{eq:nc-education}
\end{equation}
Observed incidence guides assignment only. We use the same validated
maximum likelihood estimator as on the grid and include an intercept,
\(\mathbf{Z}\) and the true regime-specific \(\mathbf{S}(\mathbf{Z})\),
but not \(\mathbf{X}\). The target is the direct coefficient \(\beta\)
before spatial propagation, not an average total effect after applying
\((\mathbf{I}-\rho \mathbf{W})^{-1}\). Under control-only spillover, a
control college's outcome rises with the share of its neighbors that are
treated. Under both-arms spillover, treated colleges also gain from
treated neighbors. When every neighbor is treated, the spillover term
adds \(\psi=0.5\) or \(0.8\), half or four-fifths of the direct effect.
The case study and the partial-participation analyses use the corrected
$t$ interval in Section \ref{estimation}, with $N$ equal to the number
of participating colleges and $N=58$ in the main case study. For that
interval, $p$ counts the regression coefficients plus $\rho$. We test
$\beta=0$ two-sided at level 0.05 by checking whether the 95\% interval
excludes zero, whereas grid power (Section
\ref{simulation-design-resampling-and-metrics}) is one-sided at level
0.025. In the Wald interval the maximum likelihood variance divides by
$N$, and the normal reference distribution ignores estimation of the
parameters. Both make the Wald interval too narrow with 58 or fewer
colleges.
```

### HW-23. Geographic reach and spread of regional incidence

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores the regions-differ-little framing, the rate-versus-death-share "because" and the grid comparison, and moves the High Incidence Focus share from the caption back to the text with its unit fixed.

Checked in `results/grid_burden_share/grid_burden_share_summary.csv`. The High Incidence Focus treated share of total grid incidence is 0.748 to 0.766 across the five queen configurations, so "0.75 to 0.77" becomes "75% to 77%" (same values, unit fixed). It stays in the text because the figure does not plot it. The gain of 2.1 to 3.1 points and the 0.84 to 0.99 correlations are Biostat's check of `grid_quadrant_spread.csv`. SR-15's sentence is the caption's second sentence.

**Current (block 1 of 2)**

```latex
Incidence-Guided Saturation sends more treatment to high-incidence regions when incidence differs between regions. On the NC map it differs little, because the four regions' mean incidence ranks lie between 24 and 34 of 58 in every year. With all 58 colleges, Incidence-Guided Saturation's treated colleges held on average 51.1\% of the summed college incidence rates, against 50.4\% for Saturation Quadrants, 50.0\% for SRS and 62.2\% for High Incidence Focus. Shares of deaths run the other way (49.2\%, 49.6\%, 49.7\% and 35.4\%), because the highest-rate colleges serve small populations. On the grid, where regional incidence differed more and equal cluster populations make rate share and death share coincide, incidence guidance raised the treated share of incidence over Saturation Quadrants by 2.1 to 3.1 percentage points on average across the five incidence configurations, while High Incidence Focus reached 0.75 to 0.77 by giving up randomization over the whole map. Across surfaces the gain rose almost linearly with the spread of the regions' mean incidence ranks (correlation 0.84 to 0.99), and the NC regions lie at the low end of that spread (Figure \ref{fig:spread-gain}).
```

**Proposed (block 1 of 2)**

```latex
Incidence-Guided Saturation sends more treatment to high-incidence
regions when incidence differs between regions. In North Carolina it
differs little, because the four regions' mean incidence ranks lay
between 24 and 34 of 58 in every year. With all 58 colleges, the treated
colleges held on average 51.1\% of the summed college incidence rates
under Incidence-Guided Saturation, 50.4\% under Saturation Quadrants,
50.0\% under SRS and 62.2\% under High Incidence Focus. Their shares of
deaths ran the other way, at 49.2\%, 49.6\%, 49.7\% and 35.4\%, because
the highest-rate colleges serve small populations. On the grid, where
regional incidence differed more and equal cluster populations make rate
share and death share coincide, incidence guidance raised the treated
share of incidence over Saturation Quadrants by 2.1 to 3.1 percentage
points on average across the five incidence configurations, and High
Incidence Focus treated 75\% to 77\% of total incidence by giving up
randomization over the whole map. Figure \ref{fig:spread-gain} shows that
the gain grows with the spread of regional mean ranks and places the
North Carolina regions at the low end of that spread. Section
\ref{igs-sensitivity-results} shows how the number of regions and the
saturation targets affect targeting and MSE on the grid.
```

**Current (block 2 of 2)**

```latex
\caption{Gain in treated incidence share from incidence guidance (Incidence-Guided Saturation minus Saturation Quadrants) against the spread of the regions' mean incidence ranks, the largest minus the smallest regional mean rank divided by the number of clusters. Filled symbols are the 50 grid surfaces (queen, production allocations) by incidence configuration, with a least-squares line fitted to them. Open triangles are the fixed NC regions in each of the four planning years, whose spread divides by 58 clusters rather than 100.}
```

**Proposed (block 2 of 2)**

```latex
\caption{Gain in treated incidence share under Incidence-Guided Saturation compared with Saturation Quadrants, by spread in regional mean incidence ranks. The share uses sums of simulated incidence values \(\mathbf{X}\) on the grid and sums of observed college rates in the case study. Spread is the highest minus the lowest regional mean rank divided by the number of clusters, 100 on the grid and 58 in the case study. Filled symbols show the 50 grid surfaces under queen contiguity by incidence configuration, with a least-squares line. Within configurations, the correlation between gain and spread ranged from 0.84 to 0.99. Open triangles show the fixed North Carolina regions in the four planning years.}
```

### HW-25. Tail membership and Monte Carlo error

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores that the added outcomes did not enlarge the sample of assignments, the reason each worst decile stays at ten draws.

**Current**

```latex
The tail estimates remain uncertain. With \(J=100\), the worst decile
contains only ten sampled draws, and extra outcomes reduce noise in
\(\hat m(a)\) but do not enlarge the allocation sample. All 112 reported settings in the tail refinement for designs with more than one possible assignment flag uncertain worst-decile membership. Under both arms, the overlap between the worst deciles identified by the two outcome halves (the Jaccard index: allocations in both divided by allocations in either) was 0--0.25 for Incidence-Guided Saturation, and approximately 61--86\% of the raw variance of estimated risk across allocations was attributable to outcome Monte Carlo noise. Saturation Quadrants' worst-decile membership was similarly unstable. Small both-arms tail differences therefore do not establish superior downside protection. Checkerboard and High Incidence Focus each have only one possible assignment for these inputs, so their 90th percentile and worst-decile mean equal their conditional mean MSE. That equality does not imply low error or absence of outcome uncertainty.
```

**Proposed**

```latex
With \(J=100\), each worst decile contains ten sampled draws. The
refinement added outcomes to the same assignments, which reduced noise
in \(\hat m(a)\) but did not enlarge the sample of assignments (Section
\ref{outcome-model-and-monte-carlo-evaluation}). Worst-decile membership
was flagged as uncertain in all 112 reported refinement settings for
designs with more than one possible assignment. Under both-arms
spillover, Incidence-Guided Saturation had a Jaccard index of 0--0.25
(Section \ref{app:nc}). Outcome Monte Carlo noise accounted for
approximately 61--86\% of its raw estimated variance in assignment risk.
Saturation Quadrants had similarly unstable worst-decile membership. The
both-arms tail differences are too small to establish better downside
protection.

Checkerboard and High Incidence Focus each have one assignment for
these inputs. Their 90th percentile and worst-decile mean therefore
equal their conditional mean MSE. That equality does not imply low
error, and their conditional MSE is still estimated from simulated
outcomes.
```

### HW-13. Shared draws and uncertainty bands in the effect-size sweep

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores that the five effect sizes on shared draws give one finding, not five confirmations, and makes the caption say what rook Checkerboard estimates.

The caption keeps Codex's band explanation. Its last sentence merges the planner's and the writer's wording.

**Current (block 1 of 2)**

```latex
Isolation Buffer ahead only at \(\beta = 3\). The five values of
\(\beta\) share their assignments and noise, so these five comparisons
are not independent confirmations of one another. They show only that, under
queen contiguity, the order did not change across the five values of
\(\beta\) we tested, given the same draws. In Figure \ref{fig:tau-mse}, the shaded bands describe the precision of a single scenario's MSE, not of the plotted average, whose standard error is smaller (Table \ref{tab:performance}).
```

**Proposed (block 1 of 2)**

```latex
Isolation Buffer ahead only at \(\beta = 3\). The five values of
\(\beta\) share their assignments and noise (Section
\ref{simulation-design-resampling-and-metrics}), so the stable order is
one finding on the same draws, not five independent confirmations.
```

**Current (block 2 of 2)**

```latex
\caption{Mean squared error by design across $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$ under (a) queen and (b) rook contiguity, pooled over the five incidence configurations. Each point averages the 160 scenarios of one design at one value of $\beta$. Shaded bands show $\pm 1$ times the average scenario-level Monte Carlo standard error of MSE. The legend follows queen-contiguity MSE at $\beta = 1$, vertical scales differ, and the gray dashed line marks SRS. Under rook contiguity, $\beta$ is not identifiable for Checkerboard. Its estimates are kept and flagged.}
```

**Proposed (block 2 of 2)**

```latex
\caption{Mean squared error by design across $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$ under (a) queen and (b) rook contiguity, pooled over the five incidence configurations. Each point averages the 160 scenarios of one design at one value of $\beta$. Shaded bands extend one average scenario-level Monte Carlo standard error of MSE above and below each point. The Monte Carlo standard error of each plotted average is smaller (Table \ref{tab:performance}). The legend follows queen-contiguity MSE at $\beta = 1$, vertical scales differ, and the gray dashed line marks SRS. For rook Checkerboard, $\beta$ is not identifiable and the coefficient on $\mathbf{Z}$ estimates $\beta-\psi$ (Section \ref{identification-diagnostics}).}
```

### HW-11. Sensitivity methods prose

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Block 1 restores "not independent evidence". Blocks 2 to 4 are unchanged.

**Current (block 1 of 4)**

```latex
We compared the nine strategies with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid, treating the 24 combinations of year, $\rho$ and $\psi$ as blocks, separately for each spillover regime and each value of $\tau$. The four years share one map and the strategies whose allocations ignore annual incidence share their simulations across years, so these $p$-values summarize consistency rather than independent evidence.
```

**Proposed (block 1 of 4)**

```latex
We compared the eight designs and SRS with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid, treating the 24 combinations of year, $\rho$ and $\psi$ as blocks, separately for each spillover regime and each value of $\tau$. The four years share one map, and some designs and SRS reuse their simulations across years (Section \ref{outcome-model-and-monte-carlo-evaluation}). These $p$-values therefore describe how consistently the designs ranked, not independent evidence.
```

**Current (block 2 of 4)**

```latex
To check that the design comparison does not depend on the spatial model, we generated data at $N = 58$ with $\tau = 0$ from the spatial lag model and fitted them by a conditional autoregressive (CAR) model \citep{besag_spatial_1974, cressie_statistics_1993} as well as by the spatial lag model. In the reverse direction, we generated data from a CAR model with the same mean and fitted them by both models. We also fitted two Bayesian models to the spatial lag data, a CAR model with the prior of Leroux and colleagues \citep{leroux_estimation_2000, lee_carbayes_2013} and a Bayesian spatial lag model \citep{lesage_introduction_2009}. The paragraphs below give the four models, their priors, the Markov chain Monte Carlo settings and the convergence checks. Finally, we replaced the proportional spillover exposure by a binary indicator of any treated neighbor and by a count of treated neighbors scaled by the mean number of neighbors, with the same form in the outcome model and the fit.
```

**Proposed (block 2 of 4)**

```latex
To check that the design comparison does not depend on the spatial model, we generated data at $N = 58$ with $\tau = 0$ from the spatial lag model and fit a conditional autoregressive (CAR) model \citep{besag_spatial_1974, cressie_statistics_1993} and the spatial lag model to those data. In the reverse direction, we generated data from a CAR model with the same mean and fit both models. We also fit two Bayesian models to the spatial lag data, a CAR model with the prior of Leroux and colleagues \citep{leroux_estimation_2000, lee_carbayes_2013} and a Bayesian spatial lag model \citep{lesage_introduction_2009}. Finally, we replaced the proportional spillover exposure by a binary indicator of any treated neighbor and by a count of treated neighbors scaled by the mean number of neighbors, with the same form in the outcome model and the fit.
```

**Current (block 3 of 4)**

```latex
\textbf{Estimation models.} Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fitted by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes and the nine strategies, with 500 fits per setting.
```

**Proposed (block 3 of 4)**

```latex
\textbf{Estimation models.} Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fit by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```

**Current (block 4 of 4)**

```latex
where $a_{jk}$ is the binary queen adjacency. We used the defaults of the \texttt{S.CARleroux} function of the R package \texttt{CARBayes} \citep{lee_carbayes_2013}, version 6.1.1: independent $N(0, 10^5)$ priors for the regression coefficients, inverse-gamma(1, 0.01) priors for $\nu^2$ and $\sigma_\phi^2$, and a uniform(0, 1) prior for $\lambda$. Each chain ran 25,000 iterations, of which the first 5,000 were discarded, thinned by 5. At $\rho = 0.5$, where fewer than half the chains converged, we reran the same data with 100,000 iterations, 20,000 discarded and thinning by 20. The Bayesian spatial lag model \citep{lesage_introduction_2009} has the same form as the maximum likelihood model, with the defaults of the \texttt{spBreg\_lag} function of \texttt{spatialreg}, version 1.4.3: a $N(\mathbf{0}, 10^{12} \mathbf{I})$ prior for the regression coefficients, the prior $p(\sigma^2) \propto 1/\sigma^2$, a Beta(1.01, 1.01) prior for $\rho$ rescaled to its feasible interval, and 10,000 draws of which the first 2,000 were discarded. Both Bayesian models used four years, $\rho \in \{0, 0.5\}$, $\psi = 0.8$, both regimes and the nine strategies, with 100 fits per setting (14,400 fits per model), and both were fitted to the same data as the maximum likelihood model. We report the posterior mean of $\beta$ and the 95\% equal-tailed credible interval. A chain was counted as converged when the effective sample size of $\beta$ was at least 400 and the absolute Geweke $z$ statistic \citep{geweke_evaluating_1992} was below 1.96.
```

**Proposed (block 4 of 4)**

```latex
where $a_{jk}$ is the binary queen adjacency. We used the defaults of the \texttt{S.CARleroux} function of the R package \texttt{CARBayes} \citep{lee_carbayes_2013}, version 6.1.1: independent $N(0, 10^5)$ priors for the regression coefficients, inverse-gamma(1, 0.01) priors for $\nu^2$ and $\sigma_\phi^2$, and a uniform(0, 1) prior for $\lambda$. Each chain ran 25,000 iterations, of which the first 5,000 were discarded, thinned by 5. At $\rho = 0.5$, where fewer than half the chains converged, we reran the same data with 100,000 iterations, 20,000 discarded and thinning by 20. The Bayesian spatial lag model \citep{lesage_introduction_2009} has the same form as the maximum likelihood model, with the defaults of the \texttt{spBreg\_lag} function of \texttt{spatialreg}, version 1.4.3: a $N(\mathbf{0}, 10^{12} \mathbf{I})$ prior for the regression coefficients, the prior $p(\sigma^2) \propto 1/\sigma^2$, a Beta(1.01, 1.01) prior for $\rho$ rescaled to its feasible interval, and 10,000 draws of which the first 2,000 were discarded. Both Bayesian models used four years, $\rho \in \{0, 0.5\}$, $\psi = 0.8$, both regimes, the eight designs and SRS, with 100 fits per setting (14,400 fits per model), and both used the same data as the maximum likelihood model. We report the posterior mean of $\beta$ and the 95\% equal-tailed credible interval. A chain was counted as converged when the effective sample size of $\beta$ was at least 400 and the absolute Geweke $z$ statistic \citep{geweke_evaluating_1992} was below 1.96.
```

### HW-30. Case study design terminology and block dependence

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Block 3 replaces "within a few percent" with the verified "within 8%" (C-6). Block 4 removes the bare "NC" and the imperative ending and keeps the independence caveat. Blocks 1 and 2 are unchanged.

C-6 checked against `tab:nc-guidance`. Under both-arms spillover the rank-based designs range from 0.92 (High Incidence Focus, 58 colleges) to 1.01, so "within 8%" is exact and "a few percent" understated it. SN-7 defines "rank-based designs" at its first use, before this paragraph.

**Current (block 1 of 4)**

```latex
When the number of trained professionals varies, with $n_c$ proportional to college population and averaging 25, and ICC $= 0.05$, power fell by 0.097 on average across designs, settings and effect sizes. For $d = 0.3$ under both-arms spillover, the mean power over the nine strategies fell from 0.59 to 0.46 with 20 colleges and from 0.88 to 0.77 with 58. MSE rose by a factor of 1.24 to 1.84 across designs and settings, and no design gained relative to the others.
```

**Proposed (block 1 of 4)**

```latex
When the number of trained professionals varies, with $n_c$ proportional to college population and averaging 25, and ICC $= 0.05$, power fell by 0.097 on average across designs, settings and effect sizes. For $d = 0.3$ under both-arms spillover, the mean power over the eight designs and SRS fell from 0.59 to 0.46 with 20 colleges and from 0.88 to 0.77 with 58. MSE rose by a factor of 1.24 to 1.84 across designs and settings, and no design gained relative to the others.
```

**Current (block 2 of 4)**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives, for each strategy, the number of colleges needed for 80\% power and the smallest college-level effect detectable with 80\% power at 29 and 58 colleges. The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not perfectly monotone in $N$ within Monte Carlo error, because the participant subsets differ between values of $N$, and we report the dips rather than smooth them.
```

**Proposed (block 2 of 4)**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives, for each design and for SRS, the number of colleges needed for 80\% power and the smallest college-level effect detectable with 80\% power at 29 and 58 colleges. The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not perfectly monotone in $N$ within Monte Carlo error, because the participant subsets differ between values of $N$, and we report the dips rather than smooth them.
```

**Current (block 3 of 4)**

```latex
Table \ref{tab:nc-guidance} gives each strategy's mean MSE relative to SRS at $\beta = 1$ for 20, 29, 40 and 58 participating colleges. The regime-specific guidance held with fewer participants, but the size of the saturation advantage depended on participation. Checkerboard and Isolation Buffer remained worse than SRS under both-arms spillover at every number of colleges, and the rank-based designs remained within a few percent of SRS. Under control-only spillover, the saturation designs' mean MSE was 22\% below that of SRS with 58 colleges, 9\% to 10\% below with 40, 4\% to 5\% below with 29 and about 1\% below with 20. With few participants, non-participating neighbors already create the exposure variation that saturation otherwise supplies. Isolation Buffer went from matching SRS under control-only spillover with 58 colleges to 2.3 times SRS's MSE with 29 colleges and 6.2 times with 15.
```

**Proposed (block 3 of 4)**

```latex
Table \ref{tab:nc-guidance} reports mean MSE for each design relative to SRS at $\beta = 1$ for 20, 29, 40 and 58 participating colleges. The regime-specific guidance held with fewer participants, but the size of the saturation advantage depended on participation. Checkerboard and Isolation Buffer remained worse than SRS under both-arms spillover at every number of colleges, and the rank-based designs remained within 8\% of SRS. Under control-only spillover, the saturation designs' mean MSE was 22\% below that of SRS with 58 colleges, 9\% to 10\% below with 40, 4\% to 5\% below with 29 and about 1\% below with 20. With few participants, non-participating neighbors already create the exposure variation that saturation otherwise supplies. Isolation Buffer went from matching SRS under control-only spillover with 58 colleges to 2.3 times SRS's MSE with 29 colleges and 6.2 times with 15.
```

**Current (block 4 of 4)**

```latex
\textbf{Formal tests on the NC results.} Friedman tests over the 24 combinations of year, $\rho$ and $\psi$ rejected equal ranks of the nine strategies in all four combinations of regime and $\tau$ ($\chi^2 = 135.3$, 153.8, 150.7 and 179.8 on 8 degrees of freedom, $p \le 2.3\times10^{-25}$; Nemenyi critical difference 2.45). Under control-only spillover, Saturation Quadrants and Incidence-Guided Saturation had the two best average ranks (1.67 and 1.79 with $\tau = 0$; 1.29 and 1.71 with $\tau = 1$), and SRS ranked seventh with $\tau = 0$ and sixth with $\tau = 1$. Under both-arms spillover with $\tau = 0$, High Incidence Focus ranked first (1.88) and SRS sixth (5.17), with the balancing and saturation designs between them. These tests treat the 24 blocks as independent, although the years share one map and some strategies share their simulations across years, so they describe the consistency of the ordering rather than independent evidence.
```

**Proposed (block 4 of 4)**

```latex
\textbf{Formal tests on the case study results.} Friedman tests over the 24 combinations of year, $\rho$ and $\psi$ rejected equal ranks of the eight designs and SRS in all four combinations of regime and $\tau$ ($\chi^2 = 135.3$, 153.8, 150.7 and 179.8 on 8 degrees of freedom, $p \le 2.3\times10^{-25}$; Nemenyi critical difference 2.45). Under control-only spillover, Saturation Quadrants and Incidence-Guided Saturation had the two best average ranks (1.67 and 1.79 with $\tau = 0$; 1.29 and 1.71 with $\tau = 1$), and SRS ranked seventh with $\tau = 0$ and sixth with $\tau = 1$. Under both-arms spillover with $\tau = 0$, High Incidence Focus ranked first (1.88) and SRS sixth (5.17), with the balancing and saturation designs between them. These tests treat the 24 blocks as independent, although the four years share one map and some designs and SRS reuse their simulations across years (Section \ref{further-sensitivity-analyses}), so they describe the consistency of the ordering rather than independent evidence.
```

### HW-6. Spatial blocking and integer treatment budgets

Codex: Accept · Biostat: CONFIRM · Planner: REVISE · Writer: REVISE (prose), REJECT (table)

**Recommendation:** Accept revised (prose), reject table block. The prose states round half to even plainly. The table change would fail `verify_manuscripts.py`, which requires two-decimal budgets in Chapter 3's `tab:nc-budgets` and an identical SI table body.

Conflict. Codex proposed "27" for "27.00". The writer rejected it because `verify_manuscripts.py` formats the budgets column of Chapter 3's `tab:nc-budgets` with two decimals and requires every SI table body to equal the chapter's. Biostat and planner did not object to the table. Decision: keep "27.00" (reject block 2), because changing it would also need a verifier change and would mix integers with 29.03 and 15.89 in one column. Chapter 3 (`Thesis_Chapter.qmd`, line 416) repeats the misstated rule "rounded to the nearest even whole number". Carry the prose fix to the chapter in 2i.

**Current (block 1 of 2)**

```latex
On the NC map, $2\times2$ Blocking uses 14 irregular groups of nearby
clusters rather than square blocks. The groups are formed once per year
by \(k\)-means on standardized projected coordinates, incidence ranks and
log(1+population), then held fixed, and contain one to eight clusters.
Half of each group, rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), is treated at random, so a one-cluster group is never treated. Each year has one
single-cluster group, which stays control because half of one cluster
rounds to zero, so the yearly totals are 28, 31, 27 and 29.
```

**Proposed (block 1 of 2)**

```latex
On the North Carolina map, $2\times2$ Blocking is the spatial-blocking
adaptation described in Table~\ref{M-tab:allocation-rules} of the article.
The rule uses 14 groups of nearby clusters. Each year, \(k\)-means
on standardized projected coordinates, incidence ranks and
$\log(1+\text{population})$ forms groups of one to eight clusters. We
hold the groups fixed within each year. The rule randomly treats half
of each group, rounded to the nearest whole number, with an exact half
rounded to the even number (1.5 to 2 and 0.5 to 0). Each year has one
single-cluster group, which therefore stays in control, and the yearly
totals are 28, 31, 27 and 29.
```

**Current (block 2 of 2)**

```latex
Checkerboard & 27.00 / 50.4 & 27.00 / 50.4 & 27.00 / 50.4 & 27.00 / 50.7 \\
High Incidence Focus & 29.00 / 22.2 & 29.00 / 23.5 & 29.00 / 23.3 & 29.00 / 22.4 \\
Saturation Quadrants & 29.03 / 49.2 & 29.03 / 49.2 & 29.03 / 49.1 & 29.03 / 49.1 \\
Isolation Buffer & 15.89 / 26.7 & 15.89 / 26.7 & 15.89 / 26.7 & 15.89 / 26.7 \\
$2\times2$ Blocking & 28.00 / 50.0 & 31.00 / 55.1 & 27.00 / 48.0 & 29.00 / 50.3 \\
Balanced Quartiles & 29.00 / 48.5 & 29.00 / 49.7 & 29.00 / 50.6 & 29.00 / 49.8 \\
Balanced Halves & 29.00 / 50.0 & 29.00 / 50.5 & 29.00 / 49.2 & 29.00 / 49.4 \\
Incidence-Guided Saturation & 29.00 / 53.2 & 29.00 / 45.1 & 29.00 / 47.9 & 29.00 / 51.8 \\
Simple Random Sampling & 29.00 / 49.5 & 29.00 / 49.5 & 29.00 / 49.5 & 29.00 / 49.5 \\
```

**Proposed (block 2 of 2):** no change. This block is rejected.

### SR-7. Power assumptions and the participation example

Codex: Accept · Biostat: CONFIRM · Planner: REVISE · Writer: not reviewed

**Recommendation:** Accept revised. Restores the takeaway that the design hardly affected the number of colleges, labels the example as an illustration and states which number of colleges each control-only power belongs to.

Checked in `application/results/methods_ext_2026-10/summaries_v2/a3_translation_m_icc.csv` (Power_Tdf, \(\rho=0.5\), \(\kappa=0.8\), \(n_c=25\), ICC 0.05, \(d=0.3\)). Control-only power is 0.488 at 29 and 0.525 at 58 colleges for SRS, and 0.508 at 29 and 0.641 at 58 for Incidence-Guided Saturation, so 0.51 belongs to 29 and 0.64 to 58. Both-arms power is 0.711 to 0.718 at 29 and 0.957 to 0.958 at 58. The file's Tau_c is 1.011, so \(d=0.3\) is about one college-level residual standard deviation.

**Current (block 1 of 2)**

```latex
Under both-arms spillover, the design hardly affected the number of colleges a trial needs. Every design except Checkerboard and Isolation Buffer needed 35 to 37 colleges for 80\% power to detect a college-level effect of one standard deviation, and about 19 colleges for 1.5 standard deviations (Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). Checkerboard and Isolation Buffer did not reach 80\% power at $\beta = 1$ even with all 58 colleges.
```

**Proposed (block 1 of 2)**

```latex
Under both-arms spillover, the design hardly affected the number of
colleges needed. With $\rho=0.5$ and $\psi=0.8\beta$, the estimated
numbers for 80\% power were 35 to 37 colleges at an effect of one
college-level residual standard deviation and about 19 colleges at 1.5
standard deviations, except for Checkerboard and Isolation Buffer
(Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). At
$\beta=1$, those two designs did not reach 80\% power with all 58
colleges.
```

**Current (block 2 of 2)**

```latex
For a trial with $n_c = 25$ trained professionals per college and an ICC of 0.05, so that $\sigma_c = 0.30$, an individual-level effect of $d = 0.3$ gives power 0.71 to 0.72 with 29 colleges and 0.96 with 58 under both-arms spillover, for SRS and for both saturation designs. Under control-only spillover the same effect gives power 0.49 and 0.53 with SRS, and 0.51 and 0.64 with Incidence-Guided Saturation\iffullsi\ (Table \ref{tab:nc-translation})\fi.
```

**Proposed (block 2 of 2)**

```latex
As an illustration, suppose $n_c = 25$ trained professionals per college and an ICC of 0.05. Then $\sigma_c = 0.30$, and an individual-level effect of $d = 0.3$ is close to one college-level residual standard deviation. With $\rho=0.5$, this effect gives power 0.71 to 0.72 with 29 colleges and 0.96 with 58 under both-arms spillover, for SRS and for both saturation designs. Under control-only spillover, the same effect gives power 0.49 with 29 colleges and 0.53 with 58 under SRS, and 0.51 and 0.64 under Incidence-Guided Saturation\iffullsi\ (Table \ref{tab:nc-translation})\fi.
```

## Batch A2 · Restored caveats and clarified scope

### HW-24. Control-only downside and the article figure

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept revised. Restores that the tail summaries are not quantiles of realized errors and not a fresh sample of assignments.

**Current (block 1 of 2)**

```latex
At the two refined queen settings, $\rho=0$ and $\rho=0.5$ with $\psi=0.8$, under control-only spillover, Saturation Quadrants' estimated 90th percentile of allocation-specific MSE is 31.8--32.9\%
below SRS and its estimated worst-decile mean 34.6--37.8\% below. Incidence-Guided
Saturation's corresponding reductions are 20.8--38.6\% and 20.7--43.4\%.
The pattern persists when the worst decile is identified and evaluated on separate halves of the simulated outcomes (Section \ref{app:nc}). \iffullsi Figure \ref{fig:nc-risk} shows the $\rho=0.5$ setting. \fi The tail estimates rest on a finite number of outcomes per allocation, so worst-decile membership is uncertain. They describe the spread of conditional MSE across the sampled allocations, not quantiles of realized errors, and because the refinement reused the main run's allocations they are not a fresh allocation sample. Read this way, the results show a descriptive downside benefit
alongside lower average error, conditional on the studied model and
sampled assignments.
```

**Proposed (block 1 of 2)**

```latex
At the two refined queen settings, $\rho=0$ and $\rho=0.5$ with
$\psi=0.8$, the control-only 90th percentile of assignment-specific MSE
was 31.8--32.9\% below SRS for Saturation Quadrants, and the worst-decile
mean was 34.6--37.8\% below. The corresponding reductions for
Incidence-Guided Saturation were 20.8--38.6\% and 20.7--43.4\%.
The split-half check gave a similar pattern (Section \ref{app:nc}).
\iffullsi Figure \ref{fig:nc-risk} shows the $\rho=0.5$ setting. \fi
Under the study model, the saturation designs had a lower worst-decile
mean of conditional MSE across the sampled assignments, as well as lower
average MSE. These summaries describe the spread of conditional MSE
across sampled assignments, not quantiles of realized estimation errors.
The refinement reused the main run's assignments rather than drawing a
fresh sample.
```

**Current (block 2 of 2)**

```latex
\iffullsi Figure \ref{fig:alloc-dist}\else Figure~\ref{M-fig:alloc-dist} of the article\fi\ shows the full allocation-to-allocation distribution for the six designs and SRS. For each allocation in \iffullsi Figure \ref{fig:alloc-dist}\else that figure\fi, MSE is estimated from repeated simulated outcomes (100 in the simulation study, at least 400 in the case study). Percentiles are computed within each setting and then averaged, so the bars show allocation-to-allocation spread, and outcome Monte Carlo noise in the per-allocation estimates widens them slightly. Under control-only spillover, most Incidence-Guided Saturation allocations had lower MSE than the average SRS allocation in both studies, whereas the SRS distribution extended well above its mean.
```

**Proposed (block 2 of 2)**

```latex
\iffullsi Figure \ref{fig:alloc-dist}\else Figure~\ref{M-fig:alloc-dist} of the article\fi\
summarizes assignment-specific MSE for the six designs and SRS. Each
estimate uses 100 outcomes per grid assignment and at least 400 per
case study assignment. Under control-only spillover, most
Incidence-Guided Saturation assignments had lower MSE than the
average SRS assignment in both studies. The SRS distribution extended
well above its mean.
```

### HW-9. Case study reproducibility and tail overlap

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Restores the Jaccard caveat (not a confidence interval or tail probability) and fixes two unclear sentences.

Block 1 drops the repeated seed and data-access sentences (with their "58-cluster"). HW-8 keeps the seed statement and the overview keeps the access statement.

**Current (block 1 of 2)**

```latex
The primary case-study numbers come from the aggregate performance and allocation-risk summaries of the main run and of the separate tail refinement. Annual tables average the six queen parameter pairs equally, without pooled Monte Carlo intervals. Each random draw used a fixed, reproducible seed. An independent check recomputed performance from the stored fits, checked annual treatment budgets and regional orderings, and checked the allocations and outcomes reused from the main run in the tail refinement. The county death-certificate source data are restricted and not publicly available; only 58-cluster summaries are reported.
```

**Proposed (block 1 of 2)**

```latex
The annual results give equal weight to the six queen parameter pairs
and are reported without pooled Monte Carlo intervals. An independent
check recomputed performance from the simulated trial fits, verified
annual treatment budgets and regional orderings, and confirmed that the
tail refinement reused the main run's assignments and outcomes.
```

**Current (block 2 of 2)**

```latex
Estimated risk variance is reported raw and with its outcome Monte Carlo
noise contribution subtracted. Signed corrected estimates are retained
without truncation at zero. None is negative in the NC refinement.
Quantiles and worst-decile means use noisy estimated conditional risks.
To check for selection bias, worst-decile assignments were also identified with one half of the simulated outcomes and evaluated with the other half, and the roles of the two halves were then reversed. The cross-selected reductions relative to SRS were similar to or slightly larger than the reported ones. The Jaccard overlap of the two worst deciles describes sampled tail membership, not a confidence interval or population tail probability. For a design with only one possible assignment, an overlap of one is uninformative.
```

**Proposed (block 2 of 2)**

```latex
We report the variance of the estimated allocation risks before and after
subtracting the contribution of outcome Monte Carlo noise. We did not
truncate the corrected estimates at zero, and none was negative in the
North Carolina refinement. The estimated risks used for quantiles and
worst-decile means still contain outcome noise.

We checked for selection bias by selecting the worst decile with one half
of the simulated outcomes, evaluating it with the other half, then
exchanging the roles of the halves. The resulting reductions relative to
SRS were similar to or slightly larger than the reductions reported from
all outcomes. The Jaccard index measures agreement between the two
sampled worst deciles. It divides the number of assignments in both
deciles by the number in either. It describes the stability of sampled
tail membership, not a confidence interval or a population tail
probability. With only one possible assignment, the index equals one by
construction and carries no information.
```

### HW-20. Case study headings and the downside table caption

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept revised. The caption restores that worst deciles are selected and evaluated with the same outcomes. The three headings are unchanged.

**Current (block 1 of 4)**

```latex
\section{Case study: accuracy, inference and reach}\label{completed-results-accuracy-inference-and-reach}
```

**Proposed (block 1 of 4)**

```latex
\section{Accuracy, inference and reach in the case study}\label{completed-results-accuracy-inference-and-reach}
```

**Current (block 2 of 4)**

```latex
\section{Case study: allocation downside}\label{allocation-downside-and-sensitivity}
```

**Proposed (block 2 of 4)**

```latex
\section{Allocation downside in the case study}\label{allocation-downside-and-sensitivity}
```

**Current (block 3 of 4)**

```latex
\section{Case study: participation, power and sensitivity}\label{app:nc-ext}
```

**Proposed (block 3 of 4)**

```latex
\section{Participation, power and sensitivity in the case study}\label{app:nc-ext}
```

**Current (block 4 of 4)**

```latex
\caption{Mean conditional MSE and worst-decile mean conditional MSE of each design divided by those of SRS in the same setting and spillover regime (below 1 is lower error than SRS). Grid: allocation-risk study, queen contiguity, $\beta=1$, averaged over 200 units, each a fixed incidence surface with one extreme combination of $\rho$ and $\psi$. NC: refined queen settings $\rho\in\{0,0.5\}$, $\psi=0.8$, averaged equally over 2018--2021 and both values of $\rho$, with the map adaptations of the designs. A design with one allocation has worst-decile mean equal to its mean. Worst-decile means are estimated from all outcomes of each assignment, without split-half selection.}
```

**Proposed (block 4 of 4)**

```latex
\caption{Mean and worst-decile mean conditional MSE divided by the corresponding SRS value in the same setting and spillover regime. Ratios below 1 favor the design over SRS. The queen grid study fixes $\beta=1$ and averages over 200 units, each an incidence surface at one extreme combination of $\rho$ and $\psi$. The North Carolina analysis uses the map adaptations and averages the refined queen settings equally over 2018--2021 and both values of $\rho\in\{0,0.5\}$, at $\psi=0.8$. For a design with one assignment, the worst-decile mean equals its mean. Worst deciles are selected and evaluated with the same simulated outcomes, without split-half selection.}
```

### HW-22. Saturation accuracy and interval coverage

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Block 1 keeps the claim that the shared saturation structure drives the advantage. Block 2 restores the contrast with exact rook non-identifiability.

**Current (block 1 of 2)**

```latex
\ref{tab:nc-year}\else Table \ref{tab:nc-annual}\fi). Both saturation designs have lower mean MSE than SRS in every primary parameter setting, and approximate Monte Carlo intervals for each setting-specific difference excluded zero. Neither is uniquely best. The two performed similarly, so incidence guidance lets a trial direct its highest treatment saturations to the highest-incidence regions with little or no loss of estimation accuracy, and their advantage over SRS comes from the saturation structure both designs share. Under control-only spillover, the accuracy of High Incidence Focus
depended on the planning year.
```

**Proposed (block 1 of 2)**

```latex
\ref{tab:nc-year}\else Table \ref{tab:nc-annual}\fi). Both saturation designs have lower mean MSE than SRS in every primary parameter setting, and approximate Monte Carlo intervals for each setting-specific difference excluded zero. The two had similar accuracy, and their advantage over SRS comes from the saturation structure they share (Section \ref{app:variants}). Under control-only spillover, the accuracy of High Incidence Focus
depended on the planning year.
```

**Current (block 2 of 2)**

```latex
Mean coverage of the 95\% intervals across the four years was between 94.1\% and 94.9\% for every design and regime, and between 94.1\% and 95.1\% in single years. The Wald interval covered only 92.6\% to 93.6\% (Table \ref{tab:interval-comparison}). Meeting the Monte Carlo precision targets did not repair that finite-sample undercoverage, but the corrected $t$ interval described in Section \ref{outcome-model-and-monte-carlo-evaluation} did (Section \ref{interval-calibration}). Table \ref{tab:nc-annual} gives bias and coverage for every design, year and regime. Checkerboard's treatment effect is identifiable on this irregular graph, but weak separation of treatment from spillover makes its estimates imprecise. This differs from the exact non-identifiability under rook contiguity on the rectangular grid.
```

**Proposed (block 2 of 2)**

```latex
Mean coverage of the 95\% intervals across the four years was between 94.1\% and 94.9\% across the eight designs and SRS in both regimes, and between 94.1\% and 95.1\% in single years. The Wald interval covered only 92.6\% to 93.6\% (Table \ref{tab:interval-comparison}). Section \ref{interval-calibration} compares the corrected $t$ interval with the Wald interval and other alternatives. Table \ref{tab:nc-annual} gives bias and coverage for every design, year and regime. Checkerboard's treatment effect is identifiable on this irregular graph, but weak separation of treatment from spillover makes its estimates imprecise. On the rook grid, by contrast, Checkerboard's treatment effect is not identifiable (Section \ref{identification-diagnostics}).
```

### HW-16. Rook Checkerboard reporting

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Block 1 says what rook estimates target and keeps the queen coverage contrast. Block 2 says what the flag marks.

**Current (block 1 of 2)**

```latex
regimes use separate noise draws. Under queen contiguity, by contrast,
Checkerboard's coverage stayed between 0.949 and 0.952 at every level of \(\rho\); under rook contiguity it was between 0.146 and 0.171\iffullsi\ (Table
\ref{tab:checkerboard-coverage})\fi. Rook Checkerboard estimates were
precise but centered on the wrong quantity, the opposite of its
imprecise but valid queen estimates. We keep these estimates in every
ranking and flag them.
```

**Proposed (block 1 of 2)**

```latex
regimes use separate noise draws. Checkerboard's coverage by level of
\(\rho\) was 0.949 to 0.952 under queen contiguity and 0.146 to 0.171
under rook contiguity\iffullsi\ (Table \ref{tab:checkerboard-coverage})\fi.
Under rook contiguity its estimates were precise but centered on
\(\beta-\psi\) rather than \(\beta\) (Section \ref{identification-diagnostics}).
Under queen contiguity they were imprecise under control-only spillover
but kept nominal coverage.
```

**Current (block 2 of 2)**

```latex
The table keeps Checkerboard's rook-contiguity estimates and flags them.
```

**Proposed (block 2 of 2)**

```latex
The table flags Checkerboard's rook-contiguity estimates, which target $\beta-\psi$ rather than $\beta$ (Section \ref{identification-diagnostics}).
```

### HW-14. Treatment-effect sweep results by regime

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Labels the average ranks as ranks among the six designs without SRS, which differ from Table S6, and splits the opening sentence.

Scope checked by Biostat against `dissertation_results_extract.txt`. These ranks exclude SRS, while Table S6 ranks six designs and SRS.

**Current**

```latex
Within that stable order, the margins moved. Under queen contiguity, the
average rank of Incidence-Guided Saturation rose from 1.41 at
\(\beta = 0.8\) to 1.85 at \(\beta = 3\), but Balanced Quartiles did not
gain on it: Incidence-Guided Saturation had the lower MSE of
the two in 103, 102, 102, 102, and 101 of the 160 blocks at
\(\beta = 0.8\), 1, 1.5, 2, and 3, and the average rank of Balanced
Quartiles stayed between 2.13 and 2.14. Isolation Buffer gained instead,
and its average rank fell from 3.71 to 3.31. The MSE of Isolation Buffer
stayed at 0.160 at every value of \(\beta\), while that of
Incidence-Guided Saturation rose from 0.088 to 0.115 and that
of Balanced Quartiles from 0.130 to 0.152, so the gap between
Incidence-Guided Saturation and Balanced Quartiles stayed
between 0.037 and 0.041. Isolation Buffer's gain came under the
control-only regime, where Incidence-Guided Saturation had a
lower MSE than Isolation Buffer in 75, 73, 56, 40, and 12 of the 80
blocks at the five values of \(\beta\). The mean control-only MSE of
Incidence-Guided Saturation rose from 0.131 at \(\beta = 0.8\)
to 0.182 at \(\beta = 3\), while that of Isolation Buffer stayed between
0.159 and 0.160 (Section \ref{app:technical-notes}). At \(\beta = 2\) the
two were level (0.160 for Isolation Buffer against 0.160), and at
\(\beta = 3\) Isolation Buffer had the lower control-only MSE (0.160
against 0.182). Incidence-Guided Saturation remained second
under the control-only regime at those two values, second to Balanced
Quartiles under the both-arms regime at every value of \(\beta\), and
lowest in pooled MSE at every value of \(\beta\). Its worst case across
blocks was the lowest of the six designs at \(\beta \le 1.5\), but at
\(\beta = 2\) and \(\beta = 3\) that of Isolation Buffer was lower (0.207
against 0.210, and 0.208 against 0.238). The mean bias of
Incidence-Guided Saturation also grew with \(\beta\), from
0.029 at \(\beta = 0.8\) to 0.047 at \(\beta = 3\), but squared bias,
averaged over scenarios, stayed at or below 3.0\% of its MSE, so the
rise in its MSE was mostly variance.
```

**Proposed**

```latex
The gaps between designs changed as the treatment effect increased.
Under queen contiguity, the average rank of Incidence-Guided Saturation
among the six designs, excluding SRS, rose from 1.41 at \(\beta = 0.8\)
to 1.85 at \(\beta = 3\). Balanced Quartiles did not gain on it.
Incidence-Guided Saturation had the lower MSE of the two in 103, 102,
102, 102, and 101 of the 160 blocks at \(\beta = 0.8\), 1, 1.5, 2, and
3, and the average rank of Balanced Quartiles stayed between 2.13 and
2.14. Isolation Buffer gained instead, and its average rank fell from
3.71 to 3.31. The MSE of Isolation Buffer stayed at 0.160 at every value
of \(\beta\), while that of Incidence-Guided Saturation rose from 0.088
to 0.115 and that of Balanced Quartiles from 0.130 to 0.152, so the gap
between Incidence-Guided Saturation and Balanced Quartiles stayed
between 0.037 and 0.041.

Isolation Buffer gained under the
control-only regime, where Incidence-Guided Saturation had a
lower MSE than Isolation Buffer in 75, 73, 56, 40, and 12 of the 80
blocks at the five values of \(\beta\). The mean control-only MSE of
Incidence-Guided Saturation rose from 0.131 at \(\beta = 0.8\)
to 0.182 at \(\beta = 3\), while that of Isolation Buffer stayed between
0.159 and 0.160 (Section \ref{app:technical-notes}). At \(\beta = 2\)
both rounded to 0.160. At
\(\beta = 3\) Isolation Buffer had the lower control-only MSE (0.160
against 0.182). Incidence-Guided Saturation remained second
under the control-only regime at those two values, second to Balanced
Quartiles under the both-arms regime at every value of \(\beta\), and
lowest in pooled MSE at every value of \(\beta\) among the six designs.

Its worst case across
blocks was the lowest of the six designs at \(\beta \le 1.5\), but at
\(\beta = 2\) and \(\beta = 3\) that of Isolation Buffer was lower (0.207
against 0.210, and 0.208 against 0.238). The mean bias of
Incidence-Guided Saturation also grew with \(\beta\), from
0.029 at \(\beta = 0.8\) to 0.047 at \(\beta = 3\), but squared bias,
averaged over scenarios, stayed at or below 3.0\% of its MSE, so the
rise in its MSE was mostly variance.
```

### HW-15. Error decomposition simulations

Codex: Accept · Biostat: REVISE · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Every term of the decomposition, not only those involving \(\hat\rho\), comes from the rerun.

**Current**

```latex
This section explains why the bias of Incidence-Guided Saturation
and the coverage of Checkerboard under rook contiguity both grow with
\(\beta\), as reported in Section \ref{stability-across-treatment-effect-sizes}. The stored results keep neither \(\hat\rho\) nor
the standard error of each fit, so we re-ran every block for the six
designs, both neighbor types, and all five values of \(\beta\), 2.4
million fits of the oracle model, recording \(\hat\rho\) and the
standard error, and the re-run reproduced the stored results in all
9,600 scenarios. Numbers below that involve \(\hat\rho\), the standard
error, or the terms of Equation \ref{eq:error-decomposition} come from
this re-run; the others come from the stored results.
```

**Proposed**

```latex
Section \ref{stability-across-treatment-effect-sizes} reports that
Incidence-Guided Saturation's bias and rook Checkerboard's coverage grow
with \(\beta\). To examine why, we reran 2.4 million oracle fits for
the six designs, both neighbor types and all five values of \(\beta\).
We recorded \(\hat\rho\) and each fit's standard error. The rerun
reproduced all 9,600 scenario summaries. Every number below that involves
\(\hat\rho\), the standard error of \(\hat\beta\) or a term of Equation
\ref{eq:error-decomposition} comes from the rerun. The other numbers come
from the original simulations.
```

### SR-1. Binary spillover and identifiable treatment effects

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: not reviewed

**Recommendation:** Accept revised. Block 1 names the comparisons that did not hold (CAR fits to spatial lag data). Block 4 drops the confusing "if exposed". Blocks 2 and 3 are unchanged.

**Current (block 1 of 4)**

```latex
The design comparison did not depend on the spatial model used to generate or fit the data, except for a CAR model fit to spatial lag data under strong spatial dependence. It also held under two other forms of spillover exposure, and formal rank tests found a consistent ordering of the strategies within each regime.
```

**Proposed (block 1 of 4)**

```latex
The design comparison held when spillover depended on the number of
treated neighbors rather than their share, and in most spatial-model
comparisons. The exceptions were CAR
models fit to spatial lag data. The maximum likelihood CAR fit was
biased under strong spatial dependence, and the Bayesian CAR fit did not
reproduce the saturation designs' control-only advantage. When spillover depended only on whether any neighbor was treated,
whether $\beta$ was identifiable depended on the assignment.
Formal rank tests summarize the ordering of the eight designs and SRS
within each regime.
```

**Current (block 2 of 4)**

```latex
\textbf{Spillover form.} Table \ref{tab:nc-forms} compares the three forms. Under the binary form, Isolation Buffer is never identifiable, because every control has a treated neighbor, and under control-only spillover no strategy was identifiable in every fit. The saturation designs were identifiable most often there (in 89.7\% and 91.7\% of fits) and, among identifiable fits, had MSE 0.70 to 0.71 times that of SRS. Under the count form, the saturation designs still had the lowest MSE under control-only spillover.
```

**Proposed (block 2 of 4)**

```latex
\textbf{Spillover form.} Table \ref{tab:nc-forms} compares the three
forms. Under binary exposure, Isolation Buffer assigns every treated
cluster zero exposure and every control cluster exposure one. The
model therefore cannot separate $\beta$ from spillover in either regime. Under
control-only spillover, $\beta$ was unidentified in some fits under
every assignment rule, including SRS. The two saturation designs
identified $\beta$ most often, in 89.7\% and 91.7\% of fits. Among
the identifiable fits for each design and for SRS, their MSEs were
0.70 to 0.71 times that of SRS. These ratios describe the subsets that
were identifiable under each rule. Under the count form, the
saturation designs still had the lowest MSE under control-only spillover.
```

**Current (block 3 of 4)**

```latex
\caption{Mean MSE relative to simple random sampling (SRS) under three spillover forms, 58 colleges, $\tau = 0$, mean over 16 blocks (year $\times$ $\rho$ $\times$ $\psi$). n.i.: not identifiable in every fit (the spillover covariate is aliased with the intercept and treatment in some or all allocations).}
```

**Proposed (block 3 of 4)**

```latex
\caption{Mean MSE relative to simple random sampling (SRS) under three spillover forms, 58 colleges, $\tau = 0$, mean over 16 blocks (year $\times$ $\rho$ $\times$ $\psi$). An entry of n.i.\ means that the treatment coefficient $\beta$ was unidentifiable in at least one sampled assignment.}
```

**Current (block 4 of 4)**

```latex
\textbf{Spillover form.} Under the binary form, $\mathbf{S}(\mathbf{Z}) = \mathbb{I}\{\text{any treated queen neighbor}\}$, times $\mathbf{1} - \mathbf{Z}$ under control-only spillover. Under the count form, $\mathbf{S}(\mathbf{Z}) = (1/\bar d) \times$ (number of treated queen neighbors), where $\bar d = 4.93$ is the mean number of queen neighbors, so a college of average degree with every neighbor treated has $\mathbf{S}(\mathbf{Z}) = 1$ and receives $\psi$, as under the proportional form.
```

**Proposed (block 4 of 4)**

```latex
\textbf{Spillover form.} Under the binary form, $\mathbf{S}(\mathbf{Z}) = \mathbb{I}\{\text{any treated queen neighbor}\}$. Under the count form, $\mathbf{S}(\mathbf{Z}) = (1/\bar d) \times$ (number of treated queen neighbors), where $\bar d = 4.93$ is the mean number of queen neighbors, so a college of average degree with every neighbor treated has $\mathbf{S}(\mathbf{Z}) = 1$ and receives $\psi$, as under the proportional form. Under control-only spillover, both the binary and the count exposure are multiplied by $\mathbf{1}-\mathbf{Z}$.
```

### SR-2. Meaning of the incidence coefficient in the new Section 3.2

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: not reviewed

**Recommendation:** Accept revised. Block 1 says a one-unit change spans the whole 0-to-1 scale. Block 2 gives the right reason for "association" (incidence is not randomized).

**Current (block 1 of 2)**

```latex
The incidence coefficient \(\tau\) in Equation~\ref{eq:outcome} is a regression coefficient, not a correlation. Incidence \(\mathbf{X}\) lies between 0 and 1 and \(\sigma = 1\), so \(\tau\) is the shift in the outcome, in residual standard deviations, from the lowest to the highest incidence, holding treatment and spillover fixed and before the spatial multiplier spreads it to neighbors. At \(\tau = 1\) that shift equals the primary intervention effect \(\beta = 1\). The value \(\tau = 0\) means that incidence has no association with the outcome, and \(\tau\) can exceed 1 or be negative.
```

**Proposed (block 1 of 2)**

```latex
The incidence coefficient \(\tau\) in Equation~\ref{eq:outcome} gives
the change in the outcome's linear predictor for a one-unit increase in
scaled incidence \(\mathbf{X}\), the full width of its 0-to-1 range,
holding treatment and spillover fixed. The change is in residual
standard deviations because \(\sigma=1\), and is measured before
spatial propagation. At \(\tau=1\), a one-unit incidence contrast
contributes as much to the linear predictor as the primary intervention
effect \(\beta=1\). At \(\tau=0\) incidence makes no direct
contribution to that predictor. The coefficient can exceed 1 or be
negative.
```

**Current (block 2 of 2)**

```latex
Before a trial, the investigator chooses \(\tau\) for planning, because there are no outcomes yet. After the trial, the investigator has both the prior incidence used for allocation and the measured outcomes, so \(\tau\) can be estimated in the same spatial lag model that estimates \(\beta\). Both grid models already include \(\mathbf{X}_k\) (Section~\ref{estimation}), so every grid fit gives \(\hat\tau\) beside \(\hat\beta\). The estimate describes how prior incidence is associated with the later outcome after adjustment for treatment and spillover. Incidence is not randomized, so \(\hat\tau\) is an adjusted association rather than an effect of incidence.
```

**Proposed (block 2 of 2)**

```latex
Before outcomes are observed, investigators specify \(\tau\) for the
planning model. After a trial, the same spatial lag model can estimate
\(\tau\) and \(\beta\) from the observed outcomes and the incidence
used for allocation. Both grid models include \(\mathbf{X}_k\)
(Section~\ref{estimation}), and every fit provides both estimates.
The oracle estimate describes the association of prior incidence
with the outcome after adjustment for treatment and spillover. The
non-oracle model adjusts for treatment and omits spillover.
Incidence is not randomized. These coefficients therefore describe
adjusted associations rather than effects of incidence.
```

### SR-6. Hypothetical education endpoint and the power scale

Codex: Accept · Biostat: REVISE · Planner: REVISE · Writer: not reviewed

**Recommendation:** Accept revised. Says that averaging divides the variance, states the assumption behind \(d=\beta\sigma_c\) and writes \(d\) in one plain sentence. Block 1 now quotes only the changed sentence of the final HW-10.

Builds on HW-10 and HW-27. Current blocks that come from those items quote their final merged wording.

**Current (block 1 of 4)**

```latex
values of $\beta$ from 0 to 2 in units of the college-level outcome
standard deviation.
```

**Proposed (block 1 of 4)**

```latex
values of $\beta$ from 0 to 2 in units of the college-level residual
standard deviation before spatial propagation.
```

**Current (block 2 of 4)**

```latex
The trial outcome is a score measured on each trained professional and averaged within each college. We write a professional's score as a college effect plus individual variation, with total variance 1, so that effects are in individual standard-deviation units. The intraclass correlation coefficient (ICC) is the share of that variance that lies between colleges, so the between-college variance is the ICC and the within-college variance is $1 - \text{ICC}$. Averaging the scores of $n_c$ professionals leaves the college part unchanged and divides the individual part by $n_c$, so the college mean has standard deviation $\sigma_c = \sqrt{\text{ICC} + (1 - \text{ICC})/n_c}$ \citep{donner_design_2000, hayes_cluster_2017}. A college-level effect of $\beta$ standard deviations is therefore an individual-level standardized effect $d = \beta\sigma_c$. With $\tau = 0$, multiplying the effect and the noise by the same constant leaves $\hat\beta/\widehat{\mathrm{SE}}$ unchanged, so power for any $n_c$, ICC and $d$ follows from the simulated college-level results. Rerunning one setting with the effect and the noise both halved reproduced its estimates to within $6\times10^{-8}$. In a sensitivity analysis, the number of trained professionals grew with college population ($n_c$ averaging 25 over participants, ICC $= 0.05$), so the noise variance differed between colleges while the fitted model stayed unweighted, as an analyst unaware of the difference would fit it.
```

**Proposed (block 2 of 4)**

```latex
For planning, suppose the education outcome is a score measured on each trained professional and averaged within each college. We represent a professional's score as a college effect plus individual variation, with total variance 1, so that effects are in individual standard-deviation units. The intraclass correlation coefficient (ICC) is the share of that variance that lies between colleges, so the between-college variance is the ICC and the within-college variance is $1 - \text{ICC}$. Averaging the scores of $n_c$ professionals leaves the college part unchanged and divides the variance of the individual part by $n_c$, so the college mean has standard deviation $\sigma_c = \sqrt{\text{ICC} + (1 - \text{ICC})/n_c}$ \citep{donner_design_2000, hayes_cluster_2017}. A structural effect of $\beta$ college-level residual standard deviations, before spatial propagation, therefore corresponds to an individual-level standardized effect $d = \beta\sigma_c$. This treats the sampling error of each college mean as the residual $\boldsymbol{\varepsilon}$ in Equation~\ref{eq:nc-education}, which the spatial multiplier also propagates when $\rho > 0$. With $\tau = 0$, multiplying the effect and the noise by the same constant leaves $\hat\beta/\widehat{\mathrm{SE}}$ unchanged, so power for any $n_c$, ICC and $d$ follows from the simulated college-level results. Rerunning one setting with the effect and the noise both halved reproduced its estimates to within $6\times10^{-8}$. In a sensitivity analysis, the number of trained professionals grew with college population ($n_c$ averaging 25 over participants, ICC $= 0.05$), so the noise variance differed between colleges while the fitted model stayed unweighted, as an analyst unaware of the difference would fit it.
```

**Current (block 3 of 4)**

```latex
college-level standard deviations under both-arms spillover for every design except Checkerboard (1.32) and Isolation Buffer (1.50)
```

**Proposed (block 3 of 4)**

```latex
college-level residual standard deviations under both-arms spillover for every design except Checkerboard (1.32) and Isolation Buffer (1.50)
```

**Current (block 4 of 4)**

```latex
minimum detectable effect (MDE) in college-level standard deviations
```

**Proposed (block 4 of 4)**

```latex
minimum detectable effect (MDE) in college-level residual standard deviations
```

## Batch A3 · Smaller restorations and wording fixes

### SR-10. Structural spillover and spatial propagation

Codex: Accept · Biostat: CONFIRM · Planner: REVISE · Writer: not reviewed

**Recommendation:** Accept revised. Restores "classical contamination" and that neighboring treatment adds to a treated cluster's own effect.

**Current (block 1 of 2)**

```latex
\end{equation} where \(\mathbf{Z} \in \{0,1\}^{100}\) is the treatment assignment
produced by the design under evaluation, \(\beta\) is the direct
treatment effect of interest, \(\mathbf{S}(\mathbf{Z})\) is the spillover term defined in
```

**Proposed (block 1 of 2)**

```latex
\end{equation} where \(\mathbf{Z} \in \{0,1\}^{100}\) is the treatment assignment
produced by the design under evaluation, \(\beta\) is the direct
treatment coefficient before spatial propagation, \(\mathbf{S}(\mathbf{Z})\) is the spillover term defined in
```

**Current (block 2 of 2)**

```latex
\end{equation} where \(\circ\) denotes the elementwise product, which multiplies the two vectors cluster by cluster. The share of treated neighbors \(\mathbf{W}\mathbf{Z}\) is an exposure mapping, a function that assigns each cluster an exposure from the full assignment vector \citep{aronow_estimating_2017,forastiere_identification_2021}. Under
the control-only regime, only untreated clusters respond to the
treatment status of their neighbors, which represents classical
contamination: an untreated cluster's outcome rises with the share of
its neighbors that are treated. Under the both-arms regime, every
cluster responds to its neighbors' treatment regardless of its own
assignment, so neighboring treatment adds to a treated cluster's own
effect. We tested four spillover magnitudes,
```

**Proposed (block 2 of 2)**

```latex
\end{equation} where \(\circ\) multiplies the two vectors cluster by cluster. The share of treated neighbors \(\mathbf{W}\mathbf{Z}\) maps the full assignment to each cluster's exposure \citep{aronow_estimating_2017,forastiere_identification_2021}.
Control-only spillover represents classical contamination. The exposure
term enters the outcome equation only for untreated clusters. An
untreated cluster's outcome rises with the share of its neighbors that
are treated. The spatial multiplier can still propagate that contribution
to treated clusters. Under both-arms spillover, the exposure term
enters the equation for treated and untreated clusters. Neighboring
treatment then adds to a treated cluster's own effect. We tested four spillover magnitudes,
```

### SR-8. How the number of colleges was interpolated

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept revised. Restores the run-in heading that the replacement dropped and that the other paragraphs of this section use. Both reviewers confirmed the content.

Builds on HW-30. Current blocks that come from those items quote their final merged wording.

Departure from all reviewers, who confirmed SR-8. The departure only restores the heading.

**Current**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives, for each design and for SRS, the number of colleges needed for 80\% power and the smallest college-level effect detectable with 80\% power at 29 and 58 colleges. The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not perfectly monotone in $N$ within Monte Carlo error, because the participant subsets differ between values of $N$, and we report the dips rather than smooth them.
```

**Proposed**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives the estimated number of colleges
needed for 80\% power and the MDE at 29 and 58 colleges. Within each
year, we used linear interpolation between the simulated sample sizes
where estimated power first reached 80\%. When we obtained an estimate
for every year, we averaged those estimates and rounded up. Individual
years can require more colleges than the rounded average. Estimates at
$\rho=0$ differed from those at $\rho=0.5$ by about two colleges or
fewer. We retain the small dips in the simulated power curves. Each
sample size uses its own participant subsets.
```

### HW-1. Neighbor definitions and Poisson incidence

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Restores that the earlier grids used rook neighbors and that 35 per 100,000 is not the case-study rate, and uses the SI's name for the earlier study.

**Current (block 1 of 2)**

```latex
We simulated a \(10\times10\) regular grid of 100 clusters under two
definitions of adjacency. Under queen contiguity, two clusters are
neighbors if they share an edge or a corner; under rook contiguity, they
must share an edge. Queen contiguity is the primary analysis. Rook contiguity is a sensitivity analysis that also keeps
continuity with earlier work, whose grids used rook neighbors. For
```

**Proposed (block 1 of 2)**

```latex
We simulated a \(10\times10\) regular grid of 100 clusters. Under queen
contiguity, two clusters are neighbors if they share an edge or a corner.
Under rook contiguity, they must share an edge. Queen contiguity is
primary. Rook contiguity is a sensitivity analysis and matches the
neighbors used on the grids of the earlier small-grid study. For
```

**Current (block 2 of 2)**

```latex
their average rank. Ties remain at this population size: the Poisson
surfaces used here have between 60 and 76 distinct values among the 100
clusters. The base rate of 35 per 100,000 was chosen to be
representative of published sudden unexpected death (SUD) incidence: it lies between the SUD
incidence rates that Mirzaei and colleagues estimated for US men (49.3
per 100,000) and women (21.7 per 100,000) aged 20--64, close to their
average \citep{mirzaei_years_2019}. It is not the observed rate of the
county data described in the case study section, which come from a
different source, age range, period, and case definition. The incidence
```

**Proposed (block 2 of 2)**

```latex
their average rank. Ties remain at this population size, and each
Poisson surface has only 60 to 76 distinct values among its 100
clusters. We set the base rate at 35 per 100,000, near the average of
the sudden unexpected death (SUD) rates that Mirzaei and colleagues
estimated for US men (49.3 per 100,000) and women (21.7 per 100,000)
aged 20--64 \citep{mirzaei_years_2019}. It is not the observed rate in
the case study, whose data come from a different source, age range,
period, and case definition. The incidence
```

### HW-2. Define SRS and shorten the Checkerboard history

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Keeps "the only deterministic design", which the next sentences rely on, and the SI's name "the earlier small-grid study". Block 2 is unchanged.

**Current (block 1 of 2)**

```latex
Table \ref{tab:designs} describes the six retained designs, and \iffullsi Figure
\ref{fig:design-samples}\else Table~\ref{M-tab:allocation-rules} of the article\fi\ shows one assignment from each, and from SRS, on a single
illustrative incidence surface. That surface, spatially autocorrelated with $\rho_X = 0.5$ and generated with a fixed seed, is not one of the simulation's incidence surfaces. Checkerboard is the strict
block-stratified allocation of the earlier small-grid study, in which no
two rook neighbors share an arm. The earlier small-grid study's Block
Stratified Sampling design also admitted relaxed allocations on its
\(3\times3\) grid, where the arms cannot be balanced; on its
\(2\times4\) and \(3\times4\) grids, as on the \(10\times10\) grid here,
the strict allocations are the two checkerboard patterns. We use the two
names interchangeably for this allocation and use Checkerboard in the
rest of the detailed analysis. Checkerboard is the only deterministic
design: it gives the same assignment in every draw, although each of its
fits still receives its own outcome noise. Every other design was drawn
independently 25 times for each incidence surface. Designs that rank
clusters by incidence break ties at random, afresh for each draw. High
Incidence Focus is random only through that tie-breaking; under the iid
and spatial mechanisms incidence values are distinct, so it gives the
same assignment in all 25 draws on a surface, and its assignment changes
only from one surface to the next. Isolation Buffer builds its
assignment from the neighbor type of the scenario, and it treats 38
clusters on average under rook contiguity and 21 under queen contiguity.
Every other design treats exactly 50.
```

**Proposed (block 1 of 2)**

```latex
Table \ref{tab:designs} describes the six retained designs.
\iffullsi Figure \ref{fig:design-samples}\else
Table~\ref{M-tab:allocation-rules} of the article\fi\ shows one assignment
from each design and one from simple random sampling (SRS). The examples
use a separate spatially autocorrelated incidence surface with
$\rho_X = 0.5$, generated with a fixed seed.

Checkerboard is the strict Block Stratified Sampling allocation of the
earlier small-grid study, in which no two rook neighbors share an arm.
That design allowed relaxed allocations on its \(3\times3\) grid because
the arms could not be balanced. On its \(2\times4\) and \(3\times4\)
grids, as on the \(10\times10\) grid here, the strict allocations are the
two checkerboard patterns. We call this allocation Checkerboard hereafter.

Checkerboard is the only deterministic design. Its assignment is fixed,
and every fit gets its own outcome noise. We drew every other design
independently 25 times per incidence surface. Designs that rank clusters
by incidence break ties at random for each draw. High Incidence Focus
varies only through tie-breaking. Under the iid and spatial mechanisms,
incidence values have no ties, and its assignment is fixed within each
surface and varies across surfaces. Isolation Buffer uses the neighbor
type of the scenario and treats an average of 38 clusters with rook
contiguity and 21 with queen contiguity. Every other design treats exactly 50.
```

**Current (block 2 of 2)**

```latex
and simple random sampling (SRS) at \(\beta = 1\), neither could be
```

**Proposed (block 2 of 2)**

```latex
and SRS at \(\beta = 1\), neither could be
```

### HW-4. Grid simulation size and shared random draws

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Keeps the technical term "common random numbers" and fixes the tense.

**Current**

```latex
Within each scenario, we drew \(J = 25\) assignments from each of the
\(K = 10\) incidence surfaces, for 250 fits per scenario, and gave every
fit its own noise vector, so that fit \((k, j)\) has outcome
\(\mathbf{Y}_{kj} = (\mathbf{I} - \rho \mathbf{W})^{-1}[\beta \mathbf{Z}_{kj} + \psi\,\mathbf{S}(\mathbf{Z}_{kj}) + \tau \mathbf{X}_k + \boldsymbol{\varepsilon}_{kj}]\).
Checkerboard's 250 fits share one assignment but have 250 distinct noise
vectors. Given a surface, its 25 pairs of assignment and noise are
independent, and the ten surfaces are independent of one another, so the
surface is the primary sampling unit. Each estimator was fit 3.2 million
times. SRS assigns exactly 50 of the 100 clusters to treatment by
complete randomization, ignoring incidence and geography. It was run
afterward with the same matched surfaces, noise keys and 250-fit
structure, bringing the complete benchmark to 14,400 scenarios and 3.6
million fits per estimator, and the original eight-design results were
unchanged.

Every random draw was seeded from a hash of a key that names the draw:
the configuration and surface index for an incidence surface; the block,
surface, and design for a set of 25 assignments; and the block for the
noise. Each draw therefore depends only on its key, and the results are
identical whether scenarios are run sequentially, in parallel, or in any
order. The keys also give common random numbers. Within a block, all
designs share the same noise vectors and incidence surfaces; the five
values of \(\beta\) reuse the same assignments and noise; and the oracle
and non-oracle models are fit to the same outcomes. Assignments and
noise are independent across blocks, but blocks within an incidence
configuration share its ten surfaces.
```

**Proposed**

```latex
We drew \(J = 25\) assignments for each of the \(K = 10\) incidence
surfaces, for 250 fits per scenario, with a separate noise vector for
each fit. Fit \((k,j)\) uses outcome
\(\mathbf{Y}_{kj} = (\mathbf{I} - \rho \mathbf{W})^{-1}[\beta \mathbf{Z}_{kj} + \psi\,\mathbf{S}(\mathbf{Z}_{kj}) + \tau \mathbf{X}_k + \boldsymbol{\varepsilon}_{kj}]\).
The assignment is fixed for all 250 Checkerboard fits. Within each
surface the 25 pairs of assignments and noise vectors are independent.
The ten surfaces are also independent, making the surface the primary
sampling unit. Each model was fit 3.2 million times for the eight
designs. SRS randomly treats 50 of the 100 clusters with no restrictions
on incidence or geography. We ran it afterward using the same surfaces,
noise and 250 fits per scenario. This brought the benchmark to 14,400
scenarios and 3.6 million fits per model and left the eight-design
results unchanged.

We seeded each draw by hashing its identifiers. A surface uses its
configuration and index. A set of 25 assignments uses its block,
surface and design. The outcome noise uses its block. This makes the
draws reproducible in any run order, and the shared keys give common
random numbers. Within each block, all designs share the same surfaces
and noise, the five values of \(\beta\) share the same assignments and
noise, and both models use the same outcomes. Assignments and noise are
independent across blocks. Blocks in the same incidence configuration
share its ten surfaces.
```

### HW-8. Simulation counts and shared case study runs

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Says which Monte Carlo standard error each target applies to and restores "not independent confirmations".

**Current**

```latex
Primary queen settings cross \(\rho=0,0.2,0.5\), \(\psi=0.5,0.8\),
both regimes, four years and nine designs (432 reported settings). A
matched sensitivity analysis lets incidence shift the outcome. It adds
\(\mathbf{X}=(\operatorname{rank}_{\rm average}(\text{rate})-0.5)/58\) with \(\tau=1\)
to both the simulated outcome and the fitted covariates (432 settings). It
changes the outcome and the adjustment together and is not an empirical
incidence--education association. A rook-contiguity sensitivity analysis of the primary model uses the extreme values $\rho=0,0.5$ (288 settings), and the two alternative regional summaries add 96 Incidence-Guided Saturation settings. The full case study reports 1,248 settings, computed from 936 distinct simulation settings, because reported settings with the same allocations and outcome model share one simulation.

We initially sampled $J=100$ assignments and simulated $R=100$ independent outcomes per distinct assignment. A design whose set of possible assignments contains only one assignment instead received $R=1,000$ outcomes for that assignment, rather than 100 repeated copies of it counted as separate allocations. An adaptive, prespecified replication schedule extended allocations or outcomes until the Monte Carlo standard error of mean MSE was at most 5\% of mean MSE and that of coverage at most 0.01. Two single-assignment settings in the matched sensitivity analysis required $R=2,000$; all randomized designs met these targets at $J=100$ and $R=100$. Each random draw used a fixed, reproducible seed. Settings with the same set of possible assignments shared their simulations, with annual treatment counts and population shares recalculated. In the primary model, Checkerboard, Saturation Quadrants, Isolation Buffer and SRS reuse the same simulations across years, because annual incidence enters neither their allocation nor the outcome model, so their repeated values are not independent
annual confirmations.
```

**Proposed**

```latex
The primary queen comparison includes 432 reported settings, crossing
\(\rho=0,0.2,0.5\), \(\psi=0.5,0.8\), both regimes, four years,
and the eight designs and SRS. In a matched sensitivity analysis,
we add \(\mathbf{X}=(\operatorname{rank}_{\rm average}(\text{rate})-0.5)/58\)
with \(\tau=1\) to both the outcome model and the fitted covariates
(432 settings). The association it builds in is hypothetical, not an
estimated relationship between incidence and education. The rook
sensitivity uses \(\rho=0,0.5\) (288 settings). Two alternative regional
summaries add 96 Incidence-Guided Saturation settings. Settings with the
same possible assignments and outcome model share a simulation, giving
936 distinct simulations for 1,248 reported settings.

We initially sampled $J=100$ assignments and simulated $R=100$
independent outcomes per distinct assignment. A design with one possible
assignment received $R=1,000$ outcomes for that assignment. The
prespecified schedule extended the number of assignments or outcomes
until the Monte Carlo standard error of mean MSE was at most 5\% of mean
MSE and that of coverage was at most 0.01. Two single-assignment settings
in the matched sensitivity analysis required $R=2,000$. All randomized
designs met both targets at $J=100$ and $R=100$. Each random draw used a
fixed seed.

In the primary model, Checkerboard, Saturation Quadrants, Isolation
Buffer and SRS share simulations across years because annual incidence
affects neither their assignments nor their outcomes. We recalculate
treatment counts and population shares for each year. Their results in
different years are therefore not independent confirmations.
```

### HW-10. Partial participation methods

Codex: Accept · Biostat: CONFIRM · Planner: REVISE · Writer: REVISE

**Recommendation:** Accept revised. Restores "only" (nonparticipants do not enter the spatial dependence), uses "that" for a college and removes the colon before the display.

**Current**

```latex
All 58 colleges will be invited to take part after the pilot, but not every college will join. We therefore repeated the case study for subsets of $N$ participating colleges, $N \in \{15, 20, 25, 29, 30, 35, 40, 45, 50, 58\}$, under the education outcome model with $\tau = 0$. Under random volunteering, 50 subsets were drawn for each $N$ and shared by every design and year, so designs are compared on the same participants. In a sensitivity analysis at $N = 20$ and $N = 29$, a college's chance of joining rose linearly with its incidence rank, so that the highest-burden college was three times as likely to join as the lowest. Each design was rebuilt on the participants with its own rule, and the four saturation regions were required to contain at least two colleges each. Spatial dependence links participating colleges only, so a participant with no participating neighbor has a zero row in the participants' weights matrix $\mathbf{W}_P$. Spillover exposure counts treated participating neighbors over all geographic neighbors, $e_j(\mathbf{Z}) = \sum_{k \in \mathcal{N}(j)} Z_k / |\mathcal{N}(j)|$, so a non-participating neighbor counts as untreated. Spillover scales with the effect, $\psi = \kappa\beta$ with $\kappa \in \{0.5, 0.8\}$, and the effect $\beta$ took nine values from 0 to 2 in units of the college-level outcome standard deviation. Settings crossed these with $\rho \in \{0, 0.5\}$, both spillover regimes, the nine strategies and the four planning years. With fewer than 58 participants, each setting used 20 allocations and 10 outcomes on each of the 50 subsets (10,000 fits), and a design with a single possible allocation on a subset used 200 outcomes instead. With all 58 colleges, each setting used 100 allocations and 100 outcomes, as in the main case study, and it reproduced the main case study's MSE and coverage within Monte Carlo error. All effect sizes share the same subsets, allocations and noise, so power curves are smooth in $\beta$. For each setting we estimated power with the corrected $t$ test above, the number of colleges needed for 80\% power, and the smallest effect detectable with 80\% power at each $N$.
```

**Proposed**

```latex
All 58 colleges will be invited to the proposed extension. To allow for
incomplete participation, we repeated the case study for
$N \in \{15,20,25,29,30,35,40,45,50,58\}$ participating colleges
under the education outcome model with $\tau=0$. Under random
volunteering, we drew 50 subsets for each $N$ and used those same
subsets for all designs, SRS and planning years. At $N=20$ and $N=29$
we also examined volunteering probabilities that increased linearly
with incidence rank. The highest-burden college was three times as
likely to join as the lowest. We reapplied each assignment rule on the
participants and required at least two colleges in each saturation region.

Spatial dependence links participating colleges only. A participant with
no participating neighbors has a zero row in $\mathbf{W}_P$. Spillover
exposure is the number of treated participating neighbors divided by the
number of geographic neighbors,
$e_j(\mathbf{Z})=\sum_{k\in\mathcal{N}(j)}Z_k/|\mathcal{N}(j)|$.
A neighbor that does not participate counts as untreated.

We set $\psi=\kappa\beta$ with $\kappa \in \{0.5,0.8\}$ and used nine
values of $\beta$ from 0 to 2 in units of the college-level outcome
standard deviation. We crossed these settings with $\rho \in \{0,0.5\}$,
both regimes, the eight designs and SRS, and four planning years. With
fewer than 58 colleges, each of the 50 subsets used 20 assignments and
10 outcomes per assignment (10,000 fits per setting). A design with a
single possible assignment on a subset used 200 outcomes for that
assignment. With 58 colleges, each setting used 100 assignments and 100
outcomes as in the main analysis and reproduced its MSE and coverage
within Monte Carlo error.

Every effect size used the same subsets, assignments and noise, which
keeps the simulated power curves smooth in $\beta$. For each setting we
estimated power with the corrected $t$ test. From the power curves we
obtained the number of colleges needed for 80\% power and the minimum
detectable effect (MDE), the smallest effect detectable with 80\% power
at each $N$.
```

### HW-18. Benchmark comparison and allocation downside

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Block 1 keeps the subject of the equivalence caveat and orders the two caveats. Block 2 is unchanged.

**Current (block 1 of 2)**

```latex
At \(\beta=1\) under queen contiguity, pooled over the five incidence
configurations and all settings of \(\rho\), \(\psi\) and spillover regime, Incidence-Guided Saturation has
mean MSE 0.091 against 0.128 for SRS, and Saturation Quadrants has 0.090. \iffullsi Table
\ref{tab:grid-srs} separates the regimes and includes all eight
strategies.\else Table~\ref{M-tab:grid-srs} of the article separates the regimes for the six designs and SRS, and Figure \ref{fig:vs-srs} adds the two variants.\fi\ Incidence-Guided Saturation's control-only ratio to SRS is 0.631,
whereas its both-arms ratio is 1.083. SRS, Balanced Quartiles and
Balanced Halves have similar both-arms accuracy. A small mean difference
does not show equivalence, and pooling the regimes would hide that the
two regimes favor different designs.
```

**Proposed (block 1 of 2)**

```latex
At \(\beta=1\) under queen contiguity, pooled over the five incidence
configurations and all settings of \(\rho\), \(\psi\) and spillover regime, Incidence-Guided Saturation has
mean MSE 0.091 against 0.128 for SRS, and Saturation Quadrants has 0.090. \iffullsi Table
\ref{tab:grid-srs} separates the regimes and includes the eight designs and SRS.\else Table~\ref{M-tab:grid-srs} of the article separates the regimes for the six designs and SRS, and Figure \ref{fig:vs-srs} adds the two variants.\fi\ Incidence-Guided Saturation's ratio to SRS was 0.631 under
control-only spillover and 1.083 under both-arms spillover, so the
pooled mean hides that the two regimes favor different designs.
SRS, Balanced Quartiles and Balanced Halves had similar both-arms MSE,
but small mean differences do not establish equivalence.
```

**Current (block 2 of 2)**

```latex
The mean is only part of the comparison, because a trial realizes a single allocation. Under control-only spillover, the worst 10\% of SRS allocations on the grid had mean conditional MSE 60\% above that of the average SRS allocation (0.346 against 0.217), and the saturation designs' worst-decile MSE was 0.54 to 0.55 times SRS's on the grid and 0.64 to 0.68 times in NC. Under both-arms spillover, their worst-decile MSE was 1.06 to 1.07 times SRS's on the grid and 0.93 to 0.94 times in NC (Table \ref{tab:downside}), so the small mean disadvantage there comes with no worse downside in NC.
```

**Proposed (block 2 of 2)**

```latex
The mean is only part of the comparison, because a trial realizes a single allocation. Under control-only spillover, the worst 10\% of SRS allocations on the grid had mean conditional MSE 60\% above that of the average SRS allocation (0.346 against 0.217), and the saturation designs' worst-decile MSE was 0.54 to 0.55 times SRS's on the grid and 0.64 to 0.68 times in the North Carolina case study. Under both-arms spillover, their worst-decile MSE was 1.06 to 1.07 times SRS's on the grid and 0.93 to 0.94 times in the North Carolina case study (Table \ref{tab:downside}). Section \ref{allocation-downside-and-sensitivity} gives the precision of the case study tail estimates.
```

### HW-12. Queen grid rank comparisons

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: REVISE

**Recommendation:** Accept revised. Uses "unit", as the paragraph does, and says that the losses outweigh the wins in the pooled mean.

**Current**

```latex
48 (96\%). Because this test follows the mean-MSE order, it separated High
Incidence Focus from \(2\times2\) Blocking, which the Nemenyi comparison
of within-block ranks did not, because each unit's MSE averages over
both spillover regimes, and there the control-only losses of
\(2\times2\) Blocking outweigh its both-arms wins.
```

**Proposed**

```latex
48 (96\%). The adjacent-pair test separated High Incidence Focus
from \(2\times2\) Blocking, which the Nemenyi comparison of within-block
ranks did not. Each unit's MSE pools the two regimes, and in that pooled
mean the control-only losses of \(2\times2\) Blocking outweigh its
both-arms wins.
```

## Batch B · Remaining prose items, accepted as proposed

### HW-5. Restriction on the death-certificate data

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Removes a repeated access statement. The overview keeps it.

**Current**

```latex
Annual incidence is \(100,000D_{iy}/P_{iy}\), and pooled incidence is
\(100,000\sum_yD_{iy}/\sum_yP_{iy}\), with the latter population sum
interpreted as person-years. All 58 clusters are present in each year,
and counts and populations sum back to statewide totals. We hold these
observed annual surfaces fixed, generate no synthetic incidence surface
and report years separately (\iffullsi Figure \ref{fig:nc-observed}\else Figure~\ref{M-fig:nc-map}C of the article\fi). The county death-certificate source data are restricted and not publicly available; only 58-cluster summaries are reported.
```

**Proposed**

```latex
Annual incidence is \(100,000D_{iy}/P_{iy}\), and pooled incidence is
\(100,000\sum_yD_{iy}/\sum_yP_{iy}\), with the latter population sum
interpreted as person-years. All 58 clusters are present in each year,
and counts and populations sum back to statewide totals. We hold these
observed annual surfaces fixed, generate no synthetic incidence surface
and report years separately (\iffullsi Figure \ref{fig:nc-observed}\else Figure~\ref{M-fig:nc-map}C of the article\fi).
```

### HW-19. Grid allocation-risk results

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Keeps the numbers, the extreme-setting scope and the noise share.

**Current**

```latex
Averaged over the 200 units per regime (fixed surfaces crossed with the extreme combinations of $\rho$ and $\psi$), Incidence-Guided Saturation's control-only mean MSE is 0.140 against SRS 0.217, the estimated 90th percentile of allocation-specific MSE is 0.173 against 0.296, and the estimated worst-decile mean is 0.191 against 0.346, lower than SRS in all 200 units on both measures. Under both arms, the
corresponding mean MSE is 0.046 against 0.043 and worst-decile mean
0.058 against 0.054. Saturation Quadrants also performs favorably under
control-only spillover. Table \ref{tab:downside} (Section \ref{allocation-downside-and-sensitivity}) compares every design with SRS on both criteria. These results come from the extreme parameter settings only and give qualified support to saturation on downside risk. The refined setting kept the control-only direction but
showed that outcome noise accounts for much of the estimated variation
between allocations (a median 63\% of the variance of Incidence-Guided Saturation's estimated risk under control-only spillover and 92\% under both arms), especially under both arms. Quantiles and worst-decile membership estimated from a finite number of outcomes remain uncertain. We therefore treat
the NC case study (Sections \ref{completed-results-accuracy-inference-and-reach} and \ref{allocation-downside-and-sensitivity}) as a separate, setting-specific
assessment rather than extrapolate exact tail protection from these
settings.
```

**Proposed**

```latex
Under control-only spillover, Incidence-Guided Saturation's mean MSE
was 0.140 against 0.217 for SRS, its estimated 90th percentile of
assignment-specific MSE was 0.173 against 0.296, and its estimated
worst-decile mean was 0.191 against 0.346. The last two were below SRS
in all 200 units, which cross fixed surfaces with the extreme values of
\(\rho\) and \(\psi\). Saturation Quadrants also performed favorably.
Under both-arms spillover, the corresponding mean MSE was 0.046 against
0.043 and the worst-decile mean was 0.058 against 0.054 (Table
\ref{tab:downside}).

In the refined grid setting, the control-only advantage persisted.
Outcome Monte Carlo noise accounted for a median 63\% of the estimated
variance in Incidence-Guided Saturation's assignment risks under
control-only spillover and 92\% under both arms. That noise makes the
estimated quantiles and worst-decile membership uncertain. Section
\ref{allocation-downside-and-sensitivity} gives a separate assessment
on the North Carolina map.
```

### HW-21. Accuracy, precision and records of fits

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Removes pipeline narration. 7,598,000 plus 2,280,000 equals 9,878,000, so "giving" is exact.

**Current (block 1 of 2)**

```latex
The main run generated 7,598,000 independent outcome fits. Every setting is complete and meets its prespecified precision targets for mean MSE and coverage. The largest relative Monte Carlo standard error of mean MSE is 0.049 and the largest Monte Carlo standard error of coverage is 0.00892. There are no warnings, aliased columns, estimates of $\rho$ at the boundary of its parameter space or estimation failures. The tail refinement adds 2,280,000 independent outcomes, and counting the outcomes reused from the main run only once gives 9,878,000 distinct outcomes across the two stages.
```

**Proposed (block 1 of 2)**

```latex
The main analysis used 7,598,000 independent outcome fits. Every setting
met its precision targets for mean MSE and coverage. The largest relative
Monte Carlo standard error of MSE was 0.049, and the largest Monte Carlo
standard error of coverage was 0.00892. We recorded no warnings,
aliased columns, boundary estimates of $\rho$ or estimation failures.
The tail refinement added 2,280,000 independent outcomes, giving
9,878,000 distinct outcomes across the two stages.
```

**Current (block 2 of 2)**

```latex
Table \ref{tab:nc-annual}\iffullsi\ and Figure \ref{fig:nc-mse-detail} give\else\ gives\fi\ MSE, bias and coverage for every design in every year, and Table \ref{tab:nc-budgets} gives the budgets. More granular setting-specific results,
precision components, residual treatment information and warning and failure records are kept as aggregate output tables. These records were checked, so the absence of warnings and failures is a verified result.
```

**Proposed (block 2 of 2)**

```latex
Table \ref{tab:nc-annual}\iffullsi\ and Figure \ref{fig:nc-mse-detail} give\else\ gives\fi\ MSE, bias and coverage for every design in every year, and Table \ref{tab:nc-budgets} gives the budgets.
```

### HW-26. Incidence sensitivity row label

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Clearer row label. It does not need a Chapter 3 change (see the note).

Conflict. The writer approved this only if Chapter 3's `tab:nc-summary-sens` changed too. Checked: that table was held out of the chapter on 2026-11-15, and `verify_manuscripts.py` reads it from the SI and compares only its numbers. The row label change passes the verifier without a chapter edit.

**Current (block 1 of 2)**

```latex
\caption{Incidence-Guided Saturation sensitivity MSE divided by matched primary MSE, by year and regime. Baseline adds rank-scaled $\mathbf{X}$ to both generation and fit; rook uses the matched extreme combinations of $\rho$ and $\psi$. Regional-summary alternatives preserve fixed regions.}
```

**Proposed (block 1 of 2)**

```latex
\caption{Incidence-Guided Saturation sensitivity MSE divided by matched primary MSE, by year and regime. The incidence sensitivity adds rank-scaled $\mathbf{X}$ to the simulated outcome and the fitted model. Rook contiguity uses the matched extreme combinations of $\rho$ and $\psi$. Regional-summary alternatives preserve fixed regions.}
```

**Current (block 2 of 2)**

```latex
Baseline / adjusted & Both & 1.016 & 1.030 & 1.015 & 1.035 \\
Baseline / adjusted & Control & 1.026 & 1.023 & 1.042 & 1.040 \\
```

**Proposed (block 2 of 2)**

```latex
Incidence in outcome and model & Both & 1.016 & 1.030 & 1.015 & 1.035 \\
Incidence in outcome and model & Control & 1.026 & 1.023 & 1.042 & 1.040 \\
```

### HW-27. Power results opening and use of MDE

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Removes the manifest jargon and leads with the conditions.

**Current (block 1 of 2)**

```latex
This section reports the participation, power and sensitivity analyses of the case study. All numbers come from aggregate summaries of runs that wrote to their own output folders, each with a manifest of parameters, code hashes and package versions. These runs did not change the main case-study results. Unless stated otherwise, results use random volunteering, $\tau = 0$, $\kappa = 0.8$ and the corrected $t$ test, and they average the four planning years.
```

**Proposed (block 1 of 2)**

```latex
Unless stated otherwise, these results use random volunteering,
$\tau=0$, $\kappa=0.8$ and the corrected $t$ test, averaged over the
four planning years.
```

**Current (block 2 of 2)**

```latex
\caption{Colleges needed for 80\% power and smallest effect detectable with 80\% power (MDE, college-level standard deviations), corrected $t$ test, random volunteering, $\rho = 0.5$, spillover $\psi = 0.8\beta$, mean over the 2018--2021 planning years. $>58$: 80\% power not reached with all 58 colleges in at least one year. $>2$: 80\% power not reached at the largest simulated effect, $\beta = 2$.}
```

**Proposed (block 2 of 2)**

```latex
\caption{Colleges needed for 80\% power and minimum detectable effect (MDE) in college-level standard deviations, corrected $t$ test, random volunteering, $\rho = 0.5$, spillover $\psi = 0.8\beta$, mean over the 2018--2021 planning years. An entry of $>58$ means estimated power stayed below 80\% with all 58 colleges in at least one year. An entry of $>2$ means estimated power stayed below 80\% at the largest simulated effect, $\beta=2$.}
```

### HW-28. Interval calibration

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Removes "narrowly" and separates the two runs. The ranges match the saved coverage summary.

**Current**

```latex
We compared several interval corrections and adopted a method only if its coverage lay between 94.0\% and 96.0\% for every design and regime. On the stored fits of the main case study ($\tau = 0$), the Wald interval covered 92.6\% to 93.6\% by design and regime. A $t$ reference distribution alone gave 93.2\% to 94.1\%, a degrees-of-freedom-corrected standard error alone 93.5\% to 94.4\%, and the two together 94.1\% to 94.9\%. In a separate run of 500 fits per setting at 58 colleges, parametric percentile bootstrap intervals with 499 replicates covered 90.1\% to 93.6\%, parametric bootstrap-$t$ intervals 93.2\% to 95.8\% and profile-likelihood intervals 91.8\% to 94.7\%. Only the $t$ interval with the corrected standard error met the rule on the stored fits. In the separate 500-fit run it covered 92.8\% to 95.8\%, narrowly missing the rule there, as every other method did. Table \ref{tab:nc-calibration} shows its size and coverage at each number of participating colleges. Table \ref{tab:interval-comparison} compares the Wald and corrected intervals for each design in the grid study and the main case study, computed from the same fits.
```

**Proposed**

```latex
We compared interval corrections against a coverage target of 94.0\%
to 96.0\% for every design and regime. For the main case study
($\tau=0$), coverage across the eight designs and SRS was 92.6\% to
93.6\% with the Wald interval, 93.2\% to 94.1\% with a $t$ reference
distribution alone, and 93.5\% to 94.4\% with a degrees-of-freedom
correction to the standard error alone. The two corrections together
gave 94.1\% to 94.9\%, the only interval that met the target.

In a separate analysis with 500 fits per setting at 58 colleges, coverage
was 90.1\% to 93.6\% for parametric percentile bootstrap intervals
with 499 replicates, 93.2\% to 95.8\% for parametric bootstrap-$t$
intervals, and 91.8\% to 94.7\% for profile-likelihood intervals. The
corrected $t$ interval covered 92.8\% to 95.8\%. Every method in that
analysis had coverage below 94.0\% in at least one design and regime.
Table \ref{tab:nc-calibration} gives size and coverage by number of
colleges. Table \ref{tab:interval-comparison} compares the Wald and
corrected intervals using the same fits within each study.
```

### HW-29. Spatial sensitivity results in three paragraphs

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: CONFIRM

**Recommendation:** Accept as proposed. Three paragraphs, with the diagnostics kept.

**Current**

```latex
\textbf{Spatial model.} Fitting a CAR model, generating the data from a CAR model, or both, left the saturation designs with 0.76 to 0.82 times SRS's MSE under control-only spillover, and Checkerboard and Isolation Buffer the worst designs under both-arms spillover. When the CAR model was fit to spatial lag data by maximum likelihood, its estimates were biased by 0.035 and 0.085 at $\rho = 0.2$ and by 0.110 and 0.205 at $\rho = 0.5$ under both-arms and control-only spillover. When it was fit to CAR data, its bias was at most 0.005. Under the spatial lag model the outcome mean includes the spatial multiplier, and the CAR coefficient on treatment therefore estimates a different quantity from $\beta$.\citep{wall_close_2004} With the longer chains, 83.1\% of the Bayesian CAR chains at $\rho = 0.5$ converged, against 46.2\% for the shorter chains, but 12.5\% still had an absolute Geweke $z$ of at least 1.96, against about 5\% expected by chance. Among converged chains, the Bayesian CAR estimate was biased by $+0.138$ under both-arms spillover and $+0.219$ under control-only spillover, with credible-interval coverage of 90.4\% and 87.7\%. For maximum likelihood on the same data, the biases were $+0.023$ and $+0.028$ and the corrected $t$ coverage 95.0\% and 94.6\%. The Bayesian CAR model also did not reproduce the saturation designs' control-only advantage. Its bias has the same source as that of the maximum likelihood CAR fit to spatial lag data, and the correctly specified Bayesian spatial lag model does not share the problem. Of the Bayesian spatial lag chains, 94.9\% converged. The median effective sample size of $\beta$ was 8,000 (minimum 5,897) and 5.1\% had an absolute Geweke $z$ of at least 1.96. Its posterior means correlated 0.999 with the maximum likelihood estimates (mean absolute difference 0.017). Its bias was at most 0.039 in absolute value in every combination of $\rho$ and regime, and its credible intervals covered 94.6\% to 95.3\% over all chains (94.5\% to 95.3\% over converged chains). The Spearman correlation of its design ordering with that of maximum likelihood averaged 0.96 to 1.00 across the four combinations of $\rho$ and regime (0.83 in the lowest single block), and the best design agreed in every block. \iffullsi Tables \ref{tab:nc-car} and \ref{tab:nc-bayes} give each strategy's MSE relative to SRS under each model.\else Table \ref{tab:nc-car} gives each strategy's MSE relative to SRS under each combination of data-generating and fitted model.\fi
```

**Proposed**

```latex
\textbf{Spatial model.} Fitting a CAR model, generating the data from a CAR model, or both, left the saturation designs with 0.76 to 0.82 times SRS's MSE under control-only spillover, and Checkerboard and Isolation Buffer the worst designs under both-arms spillover. When the CAR model was fit to spatial lag data by maximum likelihood, its estimates were biased by 0.035 and 0.085 at $\rho = 0.2$ and by 0.110 and 0.205 at $\rho = 0.5$ under both-arms and control-only spillover. When it was fit to CAR data, its bias was at most 0.005. Under the spatial lag model the outcome mean includes the spatial multiplier, and the CAR coefficient on treatment therefore estimates a different quantity from $\beta$ \citep{wall_close_2004}.

\textbf{Bayesian CAR model.} With the longer chains, 83.1\% of the Bayesian CAR chains at $\rho = 0.5$ converged, against 46.2\% for the shorter chains, but 12.5\% still had an absolute Geweke $z$ of at least 1.96, against about 5\% expected by chance. Among converged chains, the Bayesian CAR estimate was biased by $+0.138$ under both-arms spillover and $+0.219$ under control-only spillover, with credible-interval coverage of 90.4\% and 87.7\%. For maximum likelihood on the same data, the biases were $+0.023$ and $+0.028$ and the corrected $t$ coverage 95.0\% and 94.6\%. The Bayesian CAR model also did not reproduce the saturation designs' control-only advantage. Its bias has the same source as that of the maximum likelihood CAR fit to spatial lag data, and the Bayesian spatial lag model uses the same spatial multiplier as the simulated outcome.

\textbf{Bayesian spatial lag model.} Of its chains, 94.9\% converged. The median effective sample size of $\beta$ was 8,000 (minimum 5,897) and 5.1\% had an absolute Geweke $z$ of at least 1.96. Its posterior means correlated 0.999 with the maximum likelihood estimates (mean absolute difference 0.017). Its bias was at most 0.039 in absolute value in every combination of $\rho$ and regime, and its credible intervals covered 94.6\% to 95.3\% over all chains (94.5\% to 95.3\% over converged chains). The Spearman correlation of its design ordering with that of maximum likelihood averaged 0.96 to 1.00 across the four combinations of $\rho$ and regime (0.83 in the lowest single block), and the best design agreed in every block. \iffullsi Tables \ref{tab:nc-car} and \ref{tab:nc-bayes} give each strategy's MSE relative to SRS under each model.\else Table \ref{tab:nc-car} gives the MSE of each design relative to SRS under each combination of data-generating and fitted model.\fi
```

## Batch C · Remaining statistical items, accepted as proposed

### SR-3. Incidence sensitivity and error in the spatial lag coefficient

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. With \(\rho\) fixed at its true value, \(\hat\beta\) does not depend on \(\tau\).

**Current**

```latex
Changing \(\tau\) barely moved \(\hat\beta\). Over the five values of \(\tau\), no design's mean MSE for \(\beta\) changed by more than 3.4\% from its value at \(\tau = 1\) in any incidence configuration, neighbor type or regime. Both fitted models include \(\mathbf{X}\), so a change in \(\tau\) is absorbed by the coefficient on \(\mathbf{X}\) and reaches \(\hat\beta\) only through the estimate of \(\rho\). The spillover regime still decided which designs led. Under control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest MSE in every comparison. Under both-arms spillover, the best design's MSE was at most 5\% below that of SRS. The lowest-MSE design changed with \(\tau\) in 6 of the 20 comparisons, and in 5 of these the two saturation designs, whose MSEs were nearly equal, exchanged first place.
```

**Proposed**

```latex
Changing \(\tau\) barely moved \(\hat\beta\). Over the five values of \(\tau\), no design's mean MSE for \(\beta\) changed by more than 3.4\% from its value at \(\tau = 1\) in any incidence configuration, neighbor type or regime. Both fitted models include \(\mathbf{X}\). Changing \(\tau\) leaves the fitted intervention coefficient unchanged when \(\rho\) is held at its true value. Error in estimating \(\rho\) can pass part of the change into \(\hat\beta\) (Section~\ref{app:technical-notes}). The spillover regime still decided which designs led. Under control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest MSE in every comparison. Under both-arms spillover, the best design's MSE was at most 5\% below that of SRS. The lowest-MSE design changed with \(\tau\) in 6 of the 20 comparisons, and in 5 of these the two saturation designs, whose MSEs were nearly equal, exchanged first place.
```

### SR-4. Spatial lag and CAR mean structures

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. States the two conditional means correctly.

Builds on HW-11. Current blocks that come from those items quote their final merged wording.

**Current (block 1 of 2)**

```latex
In the reverse direction, we generated data from a CAR model with the same mean and fit both models.
```

**Proposed (block 1 of 2)**

```latex
We also generated data from a CAR model with the same linear predictor and fit both models.
```

**Current (block 2 of 2)**

```latex
\textbf{Estimation models.} Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fit by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```

**Proposed (block 2 of 2)**

```latex
\textbf{Estimation models.} The two spatial model families use the linear predictor $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and use the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fit by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. Thus the spatial lag outcome has conditional mean $(\mathbf{I}-\rho\mathbf{W})^{-1}\boldsymbol{\mu}$ given the assignment, whereas the CAR outcome has conditional mean $\boldsymbol{\mu}$. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```

### SR-5. Observed contribution of spatial lag bias

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Replaces the sign argument with the product \((\rho-\hat\rho)g\).

The proposal drops the old "\(g > 0\)" sign argument on purpose, as its reason states.

**Current**

```latex
For designs whose treated clusters have treated neighbors, the
last term gives a small positive bias that grows with \(\beta\). Under
queen contiguity, pooled over the 160 blocks, the bias of
Incidence-Guided Saturation rose from 0.029 at \(\beta = 0.8\)
to 0.047 at \(\beta = 3\) (Section \ref{stability-across-treatment-effect-sizes}),
while that of Isolation Buffer stayed between 0.002 and 0.004. Almost
all of the rise is the last term of Equation
\ref{eq:error-decomposition}. In our fits, \(\hat\rho\) was biased
downward by about 0.06 on average (for Incidence-Guided Saturation, by 0.071 at \(\beta = 0.8\) and 0.054 at \(\beta = 3\)). A finite-sample bias of this estimator in spatial autoregressive models is known,\citep{bao_finite-sample_2013} and we did not derive its size for this design. With
\(\rho - \hat\rho\) positive on average, the last term is positive and
grows with \(\beta\) for any design with \(g > 0\). For Incidence-Guided
Saturation, \(g\) averaged 0.371 over the queen fits, and the
last term contributed 0.021 of the bias at \(\beta = 0.8\) and 0.043 at
\(\beta = 3\); the first term contributed 0.001 at both values, and the
second 0.007 and 0.002. For Isolation Buffer, \(g\) averaged 0.020, and
the last term contributed 0.001 and 0.004.
```

**Proposed**

```latex
Under queen contiguity, pooled over the 160 blocks, the bias of
Incidence-Guided Saturation rose from 0.029 at \(\beta=0.8\) to 0.047
at \(\beta=3\) (Section \ref{stability-across-treatment-effect-sizes}).
The bias of Isolation Buffer stayed between 0.002 and 0.004. The last
term of Equation \ref{eq:error-decomposition} accounted for almost all
of the rise for Incidence-Guided Saturation.

The estimate \(\hat\rho\) was biased downward by about 0.06 on average.
For Incidence-Guided Saturation, the downward bias was 0.071 at
\(\beta=0.8\) and 0.054 at \(\beta=3\). Finite-sample bias of
\(\hat\rho\) is known in spatial autoregressive models
\citep{bao_finite-sample_2013}. The last term depends on the product
\((\rho-\hat\rho)g\) as well as on \(\beta\). For Incidence-Guided
Saturation, \(g\) averaged 0.371 over the queen fits. The last term
contributed 0.021 and 0.043 to bias at the two effect sizes. The first
term contributed 0.001 at both values, and the second contributed
0.007 and 0.002. For Isolation Buffer, \(g\) averaged 0.020 and the
last term contributed 0.001 and 0.004.
```

### SR-9. Size and coverage of the corrected test

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Observed rejection rates of 4.4% to 5.6% support "close to nominal".

Builds on HW-22 and HW-28. Current blocks that come from those items quote their final merged wording.

**Current**

```latex
The colleges-needed numbers above therefore come from a test with the nominal size.
```

**Proposed**

```latex
The sample-size estimates therefore use a test whose rejection rate was close to the nominal 5\% level in these simulations.
```

### SR-11. Scope of the grid wins and robustness results

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Scope verified. 95 + 58 + 7 = 160, so SRS is excluded.

**Current (block 1 of 4)**

```latex
Incidence-Guided Saturation also led block by block. Of the
160 queen blocks at \(\beta = 1\), Incidence-Guided Saturation
had the lowest MSE in 95 (59\%), Balanced Quartiles in 58 (36\%), and
Isolation Buffer in 7 (4\%); no other design had the lowest MSE in any
block.
```

**Proposed (block 1 of 4)**

```latex
Among the six detailed designs, Incidence-Guided Saturation had the
lowest MSE in 95 (59\%) of the 160 queen blocks at \(\beta=1\).
Balanced Quartiles led in 58 (36\%) and Isolation Buffer in 7 (4\%).
```

**Current (block 2 of 4)**

```latex
Incidence-Guided Saturation had the
lowest MSE in 97 of the 160 rook blocks (61\%) and Balanced Quartiles in
60 (38\%).
```

**Proposed (block 2 of 4)**

```latex
Among the six detailed designs, Incidence-Guided Saturation had the
lowest MSE in 97 of the 160 rook blocks (61\%) and Balanced Quartiles in
60 (38\%).
```

**Current (block 3 of 4)**

```latex
the lowest median (0.080) and the lowest worst-case MSE of any design
```

**Proposed (block 3 of 4)**

```latex
the lowest median (0.080) and worst-case MSE among the six detailed designs
```

**Current (block 4 of 4)**

```latex
Incidence-Guided Saturation also had the lowest worst-case MSE
(0.138).
```

**Proposed (block 4 of 4)**

```latex
Incidence-Guided Saturation also had the lowest worst-case MSE
among the six detailed designs (0.138).
```

### SR-12. Non-oracle rankings

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. The non-oracle order differs, so dropping "do not change the comparison" is right.

**Current**

```latex
and Checkerboard (identical). These results do not change the
comparison of designs. Incidence-Guided Saturation also had the
lowest non-oracle MSE under both neighbor types (0.087 under queen
contiguity and 0.086 under rook contiguity), and where a design's
non-oracle MSE was lower, the estimate was biased and its intervals
under-covered, so the lower MSE does not support valid inference.
```

**Proposed**

```latex
and Checkerboard (identical). Among the six detailed designs,
Incidence-Guided Saturation had the lowest pooled MSE under both
models and neighbor types. Its non-oracle MSE was 0.087 under queen
contiguity and 0.086 under rook contiguity. Designs with lower MSE
under the non-oracle model still had bias and intervals that covered
less often than nominal.
```

### SR-13. Monte Carlo standard errors and paired comparisons

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Describes the pooled standard error and the approximate ratio interval correctly.

**Current (block 1 of 3)**

```latex
Each metric is estimated from a finite number of simulated datasets, so it carries Monte Carlo error, which its Monte Carlo standard error
measures and which shrinks as more datasets are simulated. We computed
these standard errors from the surface-level metrics \citep{morris_using_2019}. For a
metric \(u\), we computed \(u_k\) on the 25 fits from surface \(k\) and
set \(\text{SE}(u) = \text{sd}(u_1, \dots, u_{10})/\sqrt{10}\), with
intervals \(u \pm t_{0.975,\,9}\,\text{SE}(u)\). A difference between
two designs is paired at the surface level: its standard error is that
of \(u_{d,k} - u_{d',k}\) over \(k\). For a metric pooled over the five
incidence configurations, we apply the same formula to the 50
configuration-by-surface means, each averaged over that surface's 32
blocks; these 50 means are independent of one another. \iffullsi Section
```

**Proposed (block 1 of 3)**

```latex
We computed Monte Carlo standard errors for bias, MSE, coverage and
power from the surface-level metrics \citep{morris_using_2019}. For a
metric \(u\), we computed \(u_k\) on the 25 fits from surface \(k\) and
set \(\text{SE}(u) = \text{sd}(u_1,\dots,u_{10})/\sqrt{10}\), with
intervals \(u \pm t_{0.975,\,9}\,\text{SE}(u)\). The standard error of
a paired difference uses \(u_{d,k}-u_{d',k}\) over \(k\). For pooled
results we apply the same formula to the 50 independent
configuration-by-surface means. When pooling both regimes, each mean
averages that surface's 32 blocks. The pooled standard deviation
includes variation across configurations. \iffullsi Section
```

**Current (block 2 of 3)**

```latex
that surface's blocks, and \(t_{49}\) intervals for the mean
difference.
```

**Proposed (block 2 of 3)**

```latex
that surface's blocks, and approximate \(t_{49}\) intervals for the mean
difference.
```

**Current (block 3 of 3)**

```latex
\caption{Each design's mean squared error (MSE) relative to simple random sampling (SRS) at $\beta = 1$ under queen contiguity, within each spillover regime. Each point is one plus the mean paired MSE difference over the 50 configuration-by-surface units divided by SRS's mean MSE, with its $t_{49}$ 95\% interval; 1 marks equal error. Filled points are significant after Holm's adjustment over the eight comparisons with SRS (Wilcoxon signed-rank test). The axis is logarithmic, ratios above 2 are printed, and rows are sorted by the control-only ratio.}
```

**Proposed (block 3 of 3)**

```latex
\caption{Each design's mean squared error (MSE) relative to simple random sampling (SRS) at $\beta = 1$ under queen contiguity, within each spillover regime. Each point is one plus the mean paired MSE difference over the 50 configuration-by-surface units divided by SRS's mean MSE, with an approximate 95\% interval obtained by transforming the paired $t_{49}$ interval for the mean difference with the SRS mean held fixed. The value 1 marks equal error. Filled points are significant after Holm's adjustment over the eight comparisons with SRS (Wilcoxon signed-rank test). The axis is logarithmic, ratios above 2 are printed, and rows are sorted by the control-only ratio.}
```

### SR-14. Study averages in the downside comparison

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Matches the downside table and removes the bare "NC".

**Current**

```latex
Table \ref{tab:downside} compares every design with SRS on both criteria, mean conditional MSE and the mean of the worst 10\% of sampled allocations, on the grid and at the refined NC settings. Under control-only spillover, both saturation designs had lower mean and worst-decile MSE than SRS in both settings. Isolation Buffer also had a lower worst decile than SRS under control-only spillover (0.62 of SRS on the grid and 0.91 in NC), but under both-arms spillover its mean was 3.80 and 3.69 times that of SRS, and it treated fewer clusters. High Incidence Focus has a single allocation in NC. That allocation's MSE was below the SRS worst-decile mean under both regimes (0.76 under both arms and 0.64 under control-only spillover), with mean MSE close to SRS (0.93 and 1.03). On the grid, however, its both-arms mean was 3.81 times that of SRS, and its NC result depends on the incidence model (see the matched baseline sensitivity below). Balanced Quartiles and Balanced Halves matched SRS on both criteria. Checkerboard and $2\times2$ Blocking did not improve on SRS on either criterion in any setting.
```

**Proposed**

```latex
Table \ref{tab:downside} compares every design with SRS on both criteria, mean conditional MSE and the mean of the worst 10\% of sampled allocations, on the grid and at the refined North Carolina settings. Under control-only spillover, both saturation designs had lower mean and worst-decile MSE than SRS in both settings. Isolation Buffer also had a lower worst decile than SRS under control-only spillover (0.62 of SRS on the grid and 0.91 in the case study), but under both-arms spillover its mean was 3.80 and 3.69 times that of SRS, and it treated fewer clusters. High Incidence Focus has a single allocation in the case study. That allocation's MSE was below the SRS worst-decile mean under both regimes (0.76 under both arms and 0.64 under control-only spillover), with mean MSE close to SRS (0.93 and 1.03). In the grid allocation-risk study, its both-arms mean MSE was 3.81 times that of SRS, and its case study result depended on the incidence model (see the matched incidence sensitivity below). Balanced Quartiles and Balanced Halves matched SRS on both criteria. For Checkerboard and $2\times2$ Blocking, mean and worst-decile MSE were at least as large as those for SRS in each study-by-regime average.
```

### SR-15. Incidence shares use different inputs on the two maps

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Its sentence is applied inside the final HW-23 caption, so it has no separate edit.

Builds on HW-23. Current blocks that come from those items quote their final merged wording.

Codex's proposed caption sentence ("The share uses sums of simulated incidence values...") is already the second sentence of the final HW-23 caption. Approving HW-23 applies SR-15.

### SR-16. What the case study replication counts count

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Wording unchanged. Block 1 now quotes the final HW-10 text.

Builds on HW-10 and HW-21. Current blocks that come from those items quote their final merged wording.

**Current (block 1 of 3)**

```latex
With 58 colleges, each setting used 100 assignments and 100
outcomes as in the main analysis and reproduced its MSE and coverage
within Monte Carlo error.

Every effect size used the same subsets, assignments and noise, which
keeps the simulated power curves smooth in $\beta$.
```

**Proposed (block 1 of 3)**

```latex
With 58 colleges, randomized rules used 100 assignments and 100
outcomes per assignment. A rule with one possible assignment used
1,000 outcomes for that assignment. The resulting MSE and coverage
matched the main analysis within Monte Carlo error.

We used the same subsets, assignments and noise across effect sizes
to reduce Monte Carlo noise in the comparisons.
```

**Current (block 2 of 3)**

```latex
The main analysis used 7,598,000 independent outcome fits.
```

**Proposed (block 2 of 3)**

```latex
The main analysis used 7,598,000 distinct simulated outcomes.
```

**Current (block 3 of 3)**

```latex
The tail refinement added 2,280,000 independent outcomes, giving
```

**Proposed (block 3 of 3)**

```latex
The tail refinement added 2,280,000 distinct outcomes, giving
```

### SR-17. Bayesian chain diagnostics

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Accurate diagnostic language.

Builds on HW-11 and HW-29. Current blocks that come from those items quote their final merged wording.

**Current (block 1 of 7)**

```latex
where fewer than half the chains converged,
```

**Proposed (block 1 of 7)**

```latex
where fewer than half the chains met the diagnostics for $\beta$,
```

**Current (block 2 of 7)**

```latex
A chain was counted as converged when the effective sample size
```

**Proposed (block 2 of 7)**

```latex
A chain met the diagnostics for $\beta$ when its effective sample size
```

**Current (block 3 of 7)**

```latex
83.1\% of the Bayesian CAR chains at $\rho = 0.5$ converged, against 46.2\% for the shorter chains,
```

**Proposed (block 3 of 7)**

```latex
83.1\% of the Bayesian CAR chains at $\rho=0.5$ met the diagnostics for $\beta$, compared with 46.2\% of the shorter chains,
```

**Current (block 4 of 7)**

```latex
Among converged chains, the Bayesian CAR estimate
```

**Proposed (block 4 of 7)**

```latex
Among chains that met those diagnostics, the Bayesian CAR estimate
```

**Current (block 5 of 7)**

```latex
For maximum likelihood on the same data,
```

**Proposed (block 5 of 7)**

```latex
Across all maximum likelihood fits to the same datasets,
```

**Current (block 6 of 7)**

```latex
Of its chains, 94.9\% converged.
```

**Proposed (block 6 of 7)**

```latex
Of its chains, 94.9\% met the diagnostics for $\beta$.
```

**Current (block 7 of 7)**

```latex
over converged chains
```

**Proposed (block 7 of 7)**

```latex
over chains meeting the diagnostics
```

### SR-18. Display precision for SRS ratios

Codex: Accept · Biostat: CONFIRM · Planner: CONFIRM · Writer: not reviewed

**Recommendation:** Accept as proposed. Two-decimal ratios follow the author's rule. Block 1 now quotes the final HW-18 text.

Builds on HW-18. Current blocks that come from those items quote their final merged wording.

Rounding changes are the point of this item. 0.631 and 1.083 print as 0.63 and 1.08, and "about 17.3" becomes the saved 17.25 and 17.26 (`igs_contrasts.csv`, verified by Biostat).

**Current (block 1 of 2)**

```latex
ratio to SRS was 0.631 under
control-only spillover and 1.083 under both-arms spillover
```

**Proposed (block 1 of 2)**

```latex
ratio to SRS was 0.63 under
control-only spillover and 1.08 under both-arms spillover
```

**Current (block 2 of 2)**

```latex
Its halves at 100/0 had slightly higher MSE than SRS under control-only spillover and about 17.3 times that of SRS under both-arms spillover.
```

**Proposed (block 2 of 2)**

```latex
Both versions of halves at 100/0 had slightly higher MSE than SRS under control-only spillover. Under both-arms spillover, their MSEs were 17.25 and 17.26 times that of SRS for the unguided and guided versions, respectively.
```

## Batch D1 · New SI items: errors, scope and planner clarity

### SN-1. Scope of the average ranks in Section 2.1

Sources: C-1 (Biostat)

**Recommendation:** Accept as proposed. The 2.700 for Incidence-Guided Saturation is a rank among the eight designs and SRS, while Table S6 prints 1.775 among six designs and SRS. Naming the set removes the apparent conflict.

**Current**

```latex
the better average rank (2.394 against 2.700 under queen contiguity and
```

**Proposed**

```latex
the better average rank among the eight designs and SRS (2.394 against 2.700 under queen contiguity and
```

### SN-2. Monte Carlo precision of the spatial-model table

Sources: C-2 (Biostat)

**Recommendation:** Accept revised. Each setting used 500 fits, and the per-setting MSE is far less precise than in the main case study, which explains why, for example, High Incidence Focus reads 1.06 here against about 0.94 in the main analysis. The comparative clause in the source wording is dropped because the two ranges overlap.

Checked against `application/results/methods_ext_2026-10/sensitivity_summaries/a6_car_settings.csv` (the table's source family): all 1,728 settings have 500 fits, and MCSE_MSE/MSE runs from 0.043 to 0.111 (median 0.068).

**Current**

```latex
\caption{Mean MSE relative to simple random sampling (SRS) under four combinations of data-generating model and fitted model: spatial lag (SAR) or conditional autoregressive (CAR). 58 colleges, $\tau = 0$, mean over 24 blocks (year $\times$ $\rho$ $\times$ $\psi$). Below 1 is lower error than SRS.}
```

**Proposed**

```latex
\caption{Mean MSE relative to simple random sampling (SRS) under four combinations of data-generating model and fitted model, each a spatial lag (SAR) or conditional autoregressive (CAR) model. 58 colleges, $\tau = 0$, mean over 24 blocks (year $\times$ $\rho$ $\times$ $\psi$). Below 1 is lower error than SRS. Each setting used 500 fits, and the Monte Carlo standard error of a setting's MSE was 4\% to 11\% of that MSE.}
```

### SN-3. Section 13 pools over the effect size

Sources: C-3 (Biostat)

**Recommendation:** Accept as proposed. The original design's control-only ratio is 0.67 in Section 13 and 0.631 in Section 12.2. Without saying that Section 13 pools over the five values of \(\beta\), the two read as a contradiction.

Checked: `docs/findings/2026-11-24-igs-sensitivity-findings.md` states that the results pool over the five effect sizes and five incidence configurations.

**Current**

```latex
SRS had MSE 0.227 under control-only spillover and 0.043 under both-arms spillover.
```

**Proposed**

```latex
Pooled over the five values of \(\beta\), SRS had MSE 0.227 under control-only spillover and 0.043 under both-arms spillover. These pooled ratios differ from those at \(\beta = 1\) in Section~\ref{srs-benchmark-and-allocation-specific-risk}.
```

### SN-4. Orphan restatements at the end of the rook section

Sources: C-4 (Biostat), HN-2 (Writer)

**Recommendation:** Accept as proposed. In the curated build these two sentences follow no figure and repeat the 0.162 coverage and the squared-bias share already reported above them. Moving them inside the existing full-only blocks keeps them in the full build.

**Current (block 1 of 2)**

```latex
\iffullsi Figure \ref{fig:coverage-rook} corresponds to Figure \ref{fig:coverage}.
\fi Every design except Checkerboard has coverage close to the nominal
level, whereas Checkerboard's coverage, 0.162 pooled over the five incidence
configurations (Table \ref{tab:performance}), is far below it.
```

**Proposed (block 1 of 2)**

```latex
\iffullsi Figure \ref{fig:coverage-rook} corresponds to Figure \ref{fig:coverage}.
Every design except Checkerboard has coverage close to the nominal
level, whereas Checkerboard's coverage, 0.162 pooled over the five incidence
configurations (Table \ref{tab:performance}), is far below it.\fi
```

**Current (block 2 of 2)**

```latex
\iffullsi Figure \ref{fig:biasvar-rook} corresponds to Figure \ref{fig:biasvar}.
\fi Nearly all of Checkerboard's MSE is squared bias, because its
coefficient on \(\mathbf{Z}\) estimates \(\beta - \psi\) (Section \ref{identification-diagnostics}).
```

**Proposed (block 2 of 2)**

```latex
\iffullsi Figure \ref{fig:biasvar-rook} corresponds to Figure \ref{fig:biasvar}.
Nearly all of Checkerboard's MSE is squared bias, because its
coefficient on \(\mathbf{Z}\) estimates \(\beta - \psi\) (Section \ref{identification-diagnostics}).\fi
```

### SN-5. Queen and rook named before they are defined

Sources: N1 (Planner)

**Recommendation:** Accept as proposed. The opening paragraph uses "queen and rook results" before either term is defined, and "labeled separately" is pipeline language.

**Current**

```latex
Queen and rook results, planning years and spillover regimes are labeled separately throughout.
```

**Proposed**

```latex
Results for the two neighbor definitions (queen and rook contiguity), the four planning years of the case study and the two spillover regimes are reported separately throughout.
```

### SN-6. "Education ability" and what is observed versus assumed

Sources: N2 (Planner)

**Recommendation:** Accept as proposed. "Education ability" is not the study's outcome. A planner also needs one statement separating the observed inputs from the planning assumptions, which are otherwise scattered over three subsections.

**Current**

```latex
estimated relationship between local SUD incidence and education
ability. Any eventual effect on SUD incidence requires a separate
clinical evaluation.
```

**Proposed**

```latex
estimated relationship between local SUD incidence and the education
outcome. Any eventual effect on SUD incidence requires a separate
clinical evaluation. The observed inputs are the service-area map, the
county-to-college assignment, and annual deaths and populations. The
outcome model, the effect size, the spatial dependence $\rho$, the
spillover regime and its magnitude $\psi$ are planning assumptions, as
are the volunteering mechanism, the number of professionals per college
and the intraclass correlation in the participation analyses (Section
\ref{partial-participation-methods}).
```

### SN-7. Define "rank-based designs" at first use

Sources: N3 (Planner)

**Recommendation:** Accept as proposed. The term is used twice in Section 16 and never defined. This is its first use. The second use is in HW-30.

Checked against Table `tab:nc-power`: control-only colleges needed at \(\beta=1.5\) are 30, 31, 29 and 32 for High Incidence Focus, Balanced Quartiles, Balanced Halves and SRS ("29 to 32"), and their MDEs at 58 colleges are 1.57, 1.46, 1.46 and 1.45 ("1.45 to 1.57").

**Current**

```latex
for SRS and the rank-based designs, while
```

**Proposed**

```latex
for SRS and the rank-based designs (High Incidence Focus, Balanced Quartiles and Balanced Halves), while
```

### SN-8. Units of the case-study MSE

Sources: N4 (Planner)

**Recommendation:** Accept revised. The first case-study MSE values (0.223 against 0.308) appear with no unit. The planner also proposed deleting "descriptive". That word is kept because it marks results reported without pooled Monte Carlo intervals.

**Current**

```latex
Under control-only spillover, Saturation Quadrants has descriptive mean MSE
```

**Proposed**

```latex
With $\sigma=1$, MSE is in squared college-level residual standard deviations. Under control-only spillover, Saturation Quadrants has descriptive mean MSE
```

### SN-9. Bare "NC" as a descriptor

Sources: N5 (Planner), HN-20 (Writer)

**Recommendation:** Accept as proposed. Author rule against a bare "NC" descriptor. HW and SR items cover the other prose uses, and the narrow table headers keep the abbreviation.

**Current (block 1 of 3)**

```latex
saturation structure they share. On the NC map, both designs use four
```

**Proposed (block 1 of 3)**

```latex
saturation structure they share. On the North Carolina map, both designs use four
```

**Current (block 2 of 3)**

```latex
draws a greedy random maximal independent set on the NC adjacency graph.
```

**Proposed (block 2 of 3)**

```latex
draws a greedy random maximal independent set on the North Carolina adjacency graph.
```

**Current (block 3 of 3)**

```latex
The NC map was not part of this study.
```

**Proposed (block 3 of 3)**

```latex
The North Carolina map was not part of this study.
```

## Batch D2 · New SI items: rule fixes and ambiguity

### SN-10. Mid-sentence "however" and an importance announcement

Sources: HN-4 (Writer)

**Recommendation:** Accept as proposed. Author rule puts "However" at the start of the sentence, and "which matters because" announces rather than states.

**Current**

```latex
Among the six detailed designs, Incidence-Guided Saturation was first or second under both regimes, at this and every other value of $\beta$, which matters because a trial team may not know in advance which regime applies. In the full benchmark that includes SRS (\iffullsi Table \ref{tab:grid-srs}\else Table~\ref{M-tab:grid-srs} of the article\fi), however, its both-arms MSE was 8\% above that of SRS, a small but statistically detectable difference (Figure \ref{fig:vs-srs}).
```

**Proposed**

```latex
A trial team may not know in advance which regime applies. Among the six detailed designs, Incidence-Guided Saturation was first or second under both regimes, at this and every other value of $\beta$. However, in the full benchmark that includes SRS (\iffullsi Table \ref{tab:grid-srs}\else Table~\ref{M-tab:grid-srs} of the article\fi), its both-arms MSE was 8\% above that of SRS, a small but statistically detectable difference (Figure \ref{fig:vs-srs}).
```

### SN-11. Which rook difference was significant

Sources: HN-5 (Writer)

**Recommendation:** Accept revised. As written, "significant" could attach to either regime. It refers to the both-arms differences, and the p-values are Holm-adjusted, as in the queen paragraph before it.

Checked in `results/srs_benchmark/srs_benchmark_summary.txt` (oracle, rook, section 3b, \(\beta=1\)): both-arms Holm p = 7.53e-03 for Incidence-Guided Saturation and 1.50e-01 for Saturation Quadrants. Control-only differences are -22.7% and -23.6% (Holm p = 6.2e-09), and Isolation Buffer is +5.8%. The p-values are the both-arms tests.

**Current**

```latex
Under rook contiguity the saturation designs had 23\% to 24\% lower control-only MSE than SRS and 3\% to 5\% higher both-arms MSE, significant for Incidence-Guided Saturation ($p = 0.0075$) but not for Saturation Quadrants ($p = 0.15$), and Isolation Buffer had 5.8\% higher control-only MSE than SRS.
```

**Proposed**

```latex
Under rook contiguity, the saturation designs had 23\% to 24\% lower control-only MSE than SRS and 3\% to 5\% higher both-arms MSE. The both-arms difference was significant after Holm's adjustment for Incidence-Guided Saturation ($p = 0.0075$) but not for Saturation Quadrants ($p = 0.15$). Isolation Buffer had 5.8\% higher control-only MSE than SRS.
```

### SN-12. Name of the incidence sensitivity and a doubled pointer

Sources: HN-6 (Writer)

**Recommendation:** Accept as proposed. Approve with HW-26. "Baseline-adjusted" would be a third name for the analysis that HW-26 and SR-14 call the incidence sensitivity, and the first sentence repeats the pointer that ends the paragraph before it.

**Current**

```latex
Table \ref{tab:nc-summary-sens} reports Incidence-Guided Saturation's matched
annual mean-MSE ratios. Rook means use only the four extreme combinations of $\rho$ and $\psi$ shared with the primary analysis. Baseline-adjusted and regional-summary means use the six
queen pairs. These are matched contrasts, not ratios of unmatched pooled
analyses.
```

**Proposed**

```latex
In Table \ref{tab:nc-summary-sens}, the rook means use only the four
extreme combinations of $\rho$ and $\psi$ shared with the primary
analysis, and the incidence-sensitivity and regional-summary means use
the six queen pairs. Each ratio is a matched contrast, not a ratio of
unmatched pooled analyses.
```

### SN-13. Vague "That is" in the High Incidence Focus caveat

Sources: HN-9 (Writer)

**Recommendation:** Accept as proposed. The caveat matters to planners, and "That is model-based accuracy" has no clear antecedent.

**Current**

```latex
years. SRS remains competitive with balancing, blocking and saturation.
Small point differences do not establish equivalence. High Incidence Focus has low both-arms MSE in the primary
model, but its allocation is deterministic for these fixed rates. That
is model-based accuracy conditional on one allocation, not evidence of
randomization-based causal validity.
```

**Proposed**

```latex
years. SRS remains competitive with balancing, blocking and saturation,
and small point differences do not establish equivalence. High Incidence
Focus has low both-arms MSE in the primary model, but its allocation is
deterministic for these fixed rates. Its low MSE reflects model-based
accuracy for that one allocation, not randomization-based causal
validity.
```

### SN-14. Semicolon and parenthetical in the interval definition

Sources: HN-8 (Writer)

**Recommendation:** Accept as proposed. HW-7 now points readers to this definition, and the 70-word sentence hides the values of \(p\) in a parenthesis split by a semicolon.

**Current**

```latex
roughly 1\% of the computing time; we also refit a 1\% sample of the
full simulation with \texttt{lagsarlm} itself as a check. Confidence intervals are 95\% intervals with a $t$ reference distribution and a degrees-of-freedom-corrected standard error, $\hat\beta \pm t_{0.975,\,N-p}\,\widehat{\text{SE}}\sqrt{N/(N-p)}$, where $\widehat{\text{SE}}$ is the estimated standard error, $N = 100$ is the number of clusters and $p$ is the number of regression coefficients plus $\rho$ (5 for the oracle model; 4 for the non-oracle model and for rook Checkerboard, whose aliased spillover column is dropped).
```

**Proposed**

```latex
roughly 1\% of the computing time. We also refit a 1\% sample of the
full simulation with \texttt{lagsarlm} itself as a check. Confidence intervals are 95\% intervals with a $t$ reference distribution and a degrees-of-freedom-corrected standard error, $\hat\beta \pm t_{0.975,\,N-p}\,\widehat{\text{SE}}\sqrt{N/(N-p)}$. Here $\widehat{\text{SE}}$ is the estimated standard error, $N = 100$ is the number of clusters, and $p$ is the number of regression coefficients plus $\rho$. The oracle model has $p = 5$, and the non-oracle model and rook Checkerboard, whose aliased spillover column is dropped, have $p = 4$.
```

### SN-15. Colon and semicolon in the identification diagnostics

Sources: HN-11 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current (block 1 of 2)**

```latex
\(\beta\) is not identifiable: the fitting routine drops the aliased
spillover column, and the coefficient on \(\mathbf{Z}\) estimates
```

**Proposed (block 1 of 2)**

```latex
\(\beta\) is not identifiable. The fitting routine drops the aliased
spillover column, and the coefficient on \(\mathbf{Z}\) estimates
```

**Current (block 2 of 2)**

```latex
information that separates it from the intercept and \(\mathbf{Z}\); \(\beta\) is
identifiable but imprecisely estimated.
```

**Proposed (block 2 of 2)**

```latex
information that separates it from the intercept and \(\mathbf{Z}\).
Therefore, \(\beta\) is identifiable but imprecisely estimated.
```

### SN-16. Semicolon in the flagging rule

Sources: HN-12 (Writer)

**Recommendation:** Accept as proposed. Semicolon rule. "Are noted here" is dropped because the sentence itself is the note.

**Current**

```latex
kept in all summaries. In the exhibits, we flag Checkerboard under rook
contiguity, where every fit is aliased; the Isolation Buffer scenarios,
each with one aliased fit among its 250, are noted here and carry no
flag in the exhibits.
```

**Proposed**

```latex
kept in all summaries. In the exhibits, we flag Checkerboard under rook
contiguity, where every fit is aliased. The Isolation Buffer scenarios,
each with one aliased fit among its 250, carry no flag in the exhibits.
```

### SN-17. Semicolons in the incidence-surface paragraph

Sources: HN-13 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current (block 1 of 2)**

```latex
are used by every scenario in that configuration, across neighbor types,
outcome-model parameters, designs, and values of \(\beta\); the design
draws and the outcome noise are shared more selectively, as described
in Section \ref{simulation-design-resampling-and-metrics}.
```

**Proposed (block 1 of 2)**

```latex
are used by every scenario in that configuration, across neighbor types,
outcome-model parameters, designs, and values of \(\beta\). Section
\ref{simulation-design-resampling-and-metrics} describes how the design
draws and the outcome noise are shared.
```

**Current (block 2 of 2)**

```latex
performance averaged over incidence surfaces; for example, the mean squared error (MSE) is
```

**Proposed (block 2 of 2)**

```latex
performance averaged over incidence surfaces. For example, the mean squared error (MSE) is
```

## Batch D3 · New SI items: remaining semicolons

### SN-18. Semicolons in the queen regime comparison

Sources: HN-14 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current (block 1 of 2)**

```latex
Blocking (0.576), and Checkerboard (2.111); at larger values of
\(\beta\), Isolation Buffer matched or passed Incidence-Guided Saturation under this regime
```

**Proposed (block 1 of 2)**

```latex
Blocking (0.576), and Checkerboard (2.111). At larger values of
\(\beta\), Isolation Buffer matched or passed Incidence-Guided Saturation under this regime
```

**Current (block 2 of 2)**

```latex
outcome model, and so the same expected MSE; the observed values, 0.161
and 0.159, differ only by Monte Carlo error, because the two regimes use
separate draws.
```

**Proposed (block 2 of 2)**

```latex
outcome model, and so the same expected MSE. The observed values, 0.161
and 0.159, differ only by Monte Carlo error because the two regimes use
separate draws.
```

### SN-19. Semicolons in the rank and mean reversals

Sources: HN-15 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current (block 1 of 2)**

```latex
\(2\times2\) Blocking (0.319), and Checkerboard (1.087); by average rank
it was \(2\times2\) Blocking (5.24), Checkerboard (5.74), and High
```

**Proposed (block 1 of 2)**

```latex
\(2\times2\) Blocking (0.319), and Checkerboard (1.087). By average rank
it was \(2\times2\) Blocking (5.24), Checkerboard (5.74), and High
```

**Current (block 2 of 2)**

```latex
under the control-only regime; mean MSE weights those losses by their
size.
```

**Proposed (block 2 of 2)**

```latex
under the control-only regime. Mean MSE weights those losses by their
size.
```

### SN-20. Semicolons and a repeated fact in the case-study precision paragraph

Sources: HN-16 (Writer)

**Recommendation:** Accept as proposed. Semicolon rule. The removed sentence ("All randomized designs met the targets at J=100 and R=100") stays in Section 5.3 through HW-8.

**Current (block 1 of 2)**

```latex
For randomized designs, the prespecified refinement schedule allowed $J=100,200,400$ and $R=100,200,400,1000,2000,4000$; whichever source of Monte Carlo error, allocations or outcomes, contributed more to the shortfall from the precision targets was extended. All randomized designs met the targets at $J=100$ and $R=100$. For designs with a single possible assignment, only outcomes were extended.
```

**Proposed (block 1 of 2)**

```latex
For randomized designs, the prespecified refinement schedule allowed $J=100,200,400$ and $R=100,200,400,1000,2000,4000$. It extended whichever source of Monte Carlo error, allocations or outcomes, contributed more to the shortfall from the precision targets. For designs with a single possible assignment, only outcomes were extended.
```

**Current (block 2 of 2)**

```latex
The overall precision includes both variation across assignments and outcome noise given an assignment; a design whose sampled assignments happened to coincide is distinguished from one with only one possible assignment.
```

**Proposed (block 2 of 2)**

```latex
The overall precision includes both variation across assignments and outcome noise given an assignment. It distinguishes a design whose sampled assignments happened to coincide from one with only one possible assignment.
```

### SN-21. Semicolons in the rook results

Sources: HN-17 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current (block 1 of 4)**

```latex
pairs; the four they did not separate were Incidence-Guided Saturation against
```

**Proposed (block 1 of 4)**

```latex
pairs. The four they did not separate were Incidence-Guided Saturation against
```

**Current (block 2 of 4)**

```latex
standard deviation stayed between 0.201 and 0.208; squared bias, averaged
```

**Proposed (block 2 of 4)**

```latex
standard deviation stayed between 0.201 and 0.208. Squared bias, averaged
```

**Current (block 3 of 4)**

```latex
regimes define the same outcome model and the same expected MSE; the
observed values, 0.477 under the both-arms regime and 0.483 under the
```

**Proposed (block 3 of 4)**

```latex
regimes define the same outcome model and the same expected MSE. The
observed values, 0.477 under the both-arms regime and 0.483 under the
```

**Current (block 4 of 4)**

```latex
Checkerboard had the highest MSE in all 80 control-only blocks; in the
80 both-arms blocks, the highest MSE belonged to High Incidence Focus in
```

**Proposed (block 4 of 4)**

```latex
Checkerboard had the highest MSE in all 80 control-only blocks. In the
80 both-arms blocks, the highest MSE belonged to High Incidence Focus in
```

### SN-22. Semicolon in the annual-table caption

Sources: HN-18 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current**

```latex
Each entry averages the six $\rho/\psi$ pairs of one spillover regime with equal weight; coverage is that of the corrected $t$ interval. Checkerboard, Saturation Quadrants, Isolation Buffer and SRS use the same simulations in every year, because annual incidence enters neither their allocation nor the outcome model, so their rows repeat and are not independent annual confirmation.}
```

**Proposed**

```latex
Each entry averages the six $\rho/\psi$ pairs of one spillover regime with equal weight, and coverage is that of the corrected $t$ interval. Checkerboard, Saturation Quadrants, Isolation Buffer and SRS use the same simulations in every year, because annual incidence enters neither their allocation nor the outcome model. Their rows therefore repeat and are not independent annual confirmations.}
```

### SN-23. Semicolon in the regional-summary sensitivity

Sources: HN-19 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current**

```latex
Alternative regional summaries change the regional ordering in some years; where an alternative summary gives the same ordering and ties, it gives the same set of possible assignments and so the same results.
```

**Proposed**

```latex
Alternative regional summaries change the regional ordering in some years. Where an alternative summary gives the same ordering and ties, it gives the same set of possible assignments and so the same results.
```

### SN-24. Semicolon in the overview's data-access sentence

Sources: HN-21 (Writer)

**Recommendation:** Accept as proposed. Semicolon rule. This is the only copy of the data-access statement left after HW-5 and HW-9.

**Current**

```latex
The North Carolina death-certificate data used in the case study are restricted and are not publicly available; only 58-cluster summaries are reported.
```

**Proposed**

```latex
The North Carolina death-certificate data used in the case study are restricted and not publicly available, so we report only 58-cluster summaries.
```

### SN-25. Remaining prose semicolons in the methods

Sources: HN-22 (Writer)

**Recommendation:** Accept as proposed. Semicolons in prose are a fixed rule of the human-writer pass, not a judgment call. No number or claim changes.

**Current (block 1 of 4)**

```latex
enters only through the spillover term; larger values of \(\rho\) induce
```

**Proposed (block 1 of 4)**

```latex
enters only through the spillover term. Larger values of \(\rho\) induce
```

**Current (block 2 of 4)**

```latex
neighbor type, \(\rho\), \(\psi\), and spillover regime a block; there
are
```

**Proposed (block 2 of 4)**

```latex
neighbor type, \(\rho\), \(\psi\), and spillover regime a block. There
are
```

**Current (block 3 of 4)**

```latex
was refined to 400 outcomes on the same assignments as a precision check; a design with only one possible assignment received its own outcome replication.
```

**Proposed (block 3 of 4)**

```latex
was refined to 400 outcomes on the same assignments as a precision check. A design with only one possible assignment received its own outcome replication.
```

**Current (block 4 of 4)**

```latex
Roanoke-Chowan Community College; these are analysis conventions, not
exclusive service commitments.
```

**Proposed (block 4 of 4)**

```latex
Roanoke-Chowan Community College. These are analysis conventions, not
exclusive service commitments.
```

## Batch E1 · Article (CTJ_Manuscript.tex)

### HC-1. Define SI at first use

Source: Writer (article) HC-1 (+2 words)

**Recommendation:** Accept as proposed.

Abbreviation used before it is defined. "SI" first appears in Methods (outcome model paragraph) but is defined only later in "Treatment assignment designs". Recommend: define it at first use and use the bare abbreviation afterward.

**Current (block 1 of 2)**

```latex
SI Section~\ref{SI-nc-sensitivity-results} compares proportional, count and binary forms in the case study.
```

**Proposed (block 1 of 2)**

```latex
Section~\ref{SI-nc-sensitivity-results} of the supplementary information (SI) compares proportional, count and binary forms in the case study.
```

**Current (block 2 of 2)**

```latex
The supplementary information (SI) gives full specifications.
```

**Proposed (block 2 of 2)**

```latex
The SI gives full specifications.
```

### HC-2. Define rook and queen contiguity at first use

Source: Writer (article) HC-2 (-3 words)

**Recommendation:** Accept as proposed.

Undefined terms. "Rook" and "queen" contiguity are used in the identifiability paragraph before the Simulation study defines them, and a clinical reader may not know them. Recommend: define them at first use and shorten the later sentence to a plain statement of which was primary. The citation stays in place.

**Current (block 1 of 2)**

```latex
With rook contiguity, Checkerboard gives $\mathbf{W}\mathbf{Z}=\mathbf{1}-\mathbf{Z}$, which means $\beta$ is not identifiable and the treatment coefficient estimates $\beta-\psi$. With queen contiguity, $\beta$ is identifiable under Checkerboard but can be estimated imprecisely.
```

**Proposed (block 1 of 2)**

```latex
When neighbors share an edge (rook contiguity), Checkerboard gives $\mathbf{W}\mathbf{Z}=\mathbf{1}-\mathbf{Z}$, which means $\beta$ is not identifiable and the treatment coefficient estimates $\beta-\psi$. When neighbors share an edge or a corner (queen contiguity), $\beta$ is identifiable under Checkerboard but can be estimated imprecisely.
```

**Current (block 2 of 2)**

```latex
Neighbors were clusters sharing an edge or a corner (queen contiguity) in the primary analysis and clusters sharing an edge (rook contiguity) in a sensitivity analysis.\cite{bivand_applied_2008}
```

**Proposed (block 2 of 2)**

```latex
The primary analysis used queen contiguity and a sensitivity analysis used rook contiguity.\cite{bivand_applied_2008}
```

### HC-3. Pivot sentence before the grid benchmark

Source: Writer (article) HC-3 (-7 words)

**Recommendation:** Accept as proposed.

Pivot sentence (AI trope). "The two spillover regimes gave opposite results" announces a contrast that the next two sentences and the subsection heading already show. Recommend: cut.

**Current**

```latex
Table~\ref{tab:grid-srs} reports the queen grid benchmark at $\beta=1$. The two spillover regimes gave opposite results. Under
```

**Proposed**

```latex
Table~\ref{tab:grid-srs} reports the queen grid benchmark at $\beta=1$. Under
```

### HC-4. Bare "NC" in the body

Source: Writer (article) HC-4 (+0 words)

**Recommendation:** Accept as proposed.

Author rule: never a bare "NC" as a descriptor. "On the NC map" is the only body use of the abbreviation, so the definition in the case-study paragraph can go too (Table 1 caption defines NC on its own). Recommend: apply both parts.

**Current (block 1 of 2)**

```latex
On the NC map the best design
```

**Proposed (block 1 of 2)**

```latex
On the North Carolina map the best design
```

**Current (block 2 of 2)**

```latex
partition the 100 counties of North Carolina (NC).
```

**Proposed (block 2 of 2)**

```latex
partition the 100 counties of North Carolina.
```

### HC-5. Active voice for the accuracy measures

Source: Writer (article) HC-5 (-5 words)

**Recommendation:** Accept as proposed.

Agentless passive with the verb repeated, and MSE is defined a second time (it was defined in the Introduction). Recommend: make "we" the subject and use the abbreviation.

**Current**

```latex
Accuracy was measured by the mean squared error (MSE), the average of $(\hat\beta-\beta)^2$ over simulated trials, and coverage was measured as the proportion of 95\% confidence intervals containing $\beta$.
```

**Proposed**

```latex
We measured accuracy by the MSE, the average of $(\hat\beta-\beta)^2$ over simulated trials, and coverage by the proportion of 95\% confidence intervals containing $\beta$.
```

### HC-6. Split the rank-scaling sentence

Source: Writer (article) HC-6 (-1 words)

**Recommendation:** Accept as proposed.

Sentence of about 50 words (rank-scaling of Poisson counts). Recommend: split after the definition of the rank. Revised 2026-11-25 so the split does not leave a ", so" join.

**Current**

```latex
where $r_i$ is cluster $i$'s rank with average ranks for ties, to put them on the same $(0,1)$ scale as the other configurations, so that one value
```

**Proposed**

```latex
where $r_i$ is cluster $i$'s rank with average ranks for ties. This put them on the same $(0,1)$ scale as the other configurations. As a result, one value
```

### HC-7. Fold a "Therefore" sentence into "Because"

Source: Writer (article) HC-7 (+0 words)

**Recommendation:** Decline (revised 2026-11-25). The author keeps sentence-initial "Therefore" and avoids sentence-initial "Because" and ", so" joins, so the current sentence stays.

Same template ending three Results paragraphs and others ("... . Therefore, X does not show / need not ..."). The body has seven sentence-initial "Therefore"s. Recommend: fold this one into a "Because" clause. HC-8 to HC-10 trim three more, leaving three (the Introduction, the low-MSE-and-valid-inference sentence and the grid-ranking sentence).

**Current**

```latex
The two changes were made together. Therefore, their separate contributions cannot be distinguished.
```

**Proposed**

```latex
Because the two changes were made together, their separate contributions cannot be distinguished.
```

### HC-8. Join a "Therefore" sentence (tail noise)

Source: Writer (article) HC-8 (+0 words)

**Recommendation:** Decline (revised 2026-11-25). The author keeps sentence-initial "Therefore" and avoids sentence-initial "Because" and ", so" joins, so the current sentence stays.

Same "Therefore" template (see HC-7). Recommend: join with ", so".

**Current**

```latex
was Monte Carlo error (SI Section~\ref{SI-allocation-downside-and-sensitivity}). Therefore, small both-arms differences
```

**Proposed**

```latex
was Monte Carlo error (SI Section~\ref{SI-allocation-downside-and-sensitivity}), so small both-arms differences
```

### HC-9. Join a "Therefore" sentence (High Incidence Focus)

Source: Writer (article) HC-9 (+0 words)

**Recommendation:** Decline (revised 2026-11-25). The author keeps sentence-initial "Therefore" and avoids sentence-initial "Because" and ", so" joins, so the current sentence stays.

Same "Therefore" template (see HC-7). Recommend: join with ", so".

**Current**

```latex
was sensitive to incidence adjustment. Therefore, it does not support a general recommendation.
```

**Proposed**

```latex
was sensitive to incidence adjustment, so it does not support a general recommendation.
```

### HC-10. Join a "Therefore" sentence (noise)

Source: Writer (article) HC-10 (+0 words)

**Recommendation:** Decline (revised 2026-11-25). The author keeps sentence-initial "Therefore" and avoids sentence-initial "Because" and ", so" joins, so the current sentence stays.

Same "Therefore" template (see HC-7). A two-sentence chain for a one-step inference. Recommend: join with ", so".

**Current**

```latex
The noise $\boldsymbol{\varepsilon}$ is independent standard normal. Therefore, the intervention effect
```

**Proposed**

```latex
The noise $\boldsymbol{\varepsilon}$ is independent standard normal, so the intervention effect
```

## Batch E2 · Article (CTJ_Manuscript.tex)

### HC-11. Repeated "In this paper"

Source: Writer (article) HC-11 (-3 words)

**Recommendation:** Accept as proposed.

"In this paper," opens two consecutive Introduction paragraphs. Recommend: drop it from the second.

**Current**

```latex
In this paper, we also use a case study,
```

**Proposed**

```latex
We also use a case study,
```

### HC-12. Split the coverage and Checkerboard sentence

Source: Writer (article) HC-12 (+0 words)

**Recommendation:** Accept as proposed.

Two unrelated facts yoked by "while" (coverage across designs, then Checkerboard's rook bias), and "Checkerboard's was" leaves the noun implied. Recommend: split and name "effect".

**Current**

```latex
Coverage was 94--95\% for every design whose effect was identifiable, while under rook contiguity Checkerboard's was not identifiable and its estimate was substantially biased.
```

**Proposed**

```latex
Coverage was 94--95\% for every design whose effect was identifiable. Under rook contiguity, Checkerboard's effect was not identifiable and its estimate was substantially biased.
```

### HC-13. Repeated separation caveat (optional)

Source: Writer (article) HC-13 (-13 words)

**Recommendation:** Accept as proposed (optional).

Repeated caveat. The Introduction's second paragraph already says balancing incidence does not ensure that treatment can be separated from neighboring exposure, with the Checkerboard example. Recommend: cut the second statement. Optional: it is broader (it covers targeting and saturation too), so keep it if the author wants that generalization stated.

**Current**

```latex
(saturation). None of these goals ensures that treatment can be separated from neighboring exposure. A trial draws only one assignment.
```

**Proposed**

```latex
(saturation). A trial draws only one assignment.
```

### HC-14. Repeated incidence-guides-allocation caveat

Source: Writer (article) HC-14 (-4 words)

**Recommendation:** Accept as proposed.

Repeated caveat. That incidence guides only the allocation when tau = 0 is said in the Introduction, in the outcome-model paragraph ("With tau=0, incidence guides only the allocation, as in the case study") and in the Figure 2 caption. Recommend: keep the case-specific reason and drop the restatement.

**Current**

```latex
The education outcome is not expected to depend on SUD incidence. Incidence guided only the allocation, and $\tau=1$ served as a hypothetical sensitivity analysis.
```

**Proposed**

```latex
The education outcome is not expected to depend on SUD incidence, so $\tau=1$ served only as a hypothetical sensitivity analysis.
```

### HC-15. Redundant gloss on the elementwise product (optional)

Source: Writer (article) HC-15 (-8 words)

**Recommendation:** Accept as proposed (optional).

Redundant gloss. "Elementwise product" already means multiplying cluster by cluster. Recommend: cut. Keep it if a coauthor asked for the plain-language gloss for clinical readers.

**Current**

```latex
Here $\psi$ is the spillover magnitude and $\circ$ is the elementwise product, which multiplies the two vectors cluster by cluster.
```

**Proposed**

```latex
Here $\psi$ is the spillover magnitude and $\circ$ is the elementwise product.
```

### HC-16. Reorder the coverage sentence

Source: Writer (article) HC-16 (-1 words)

**Recommendation:** Accept as proposed.

Stacked "and ... and" in a fronted phrase delays the subject. Recommend: lead with the averaging and move the scope after the subject.

**Current**

```latex
Across the six designs and SRS and both regimes, averaged over years, mean coverage of the $t$-based interval was 94.2--94.9\%,
```

**Proposed**

```latex
Averaged over years, mean coverage of the $t$-based interval across the six designs, SRS and both regimes was 94.2--94.9\%,
```

### HC-17. Split the SRS downside sentence

Source: Writer (article) HC-17 (-1 words)

**Recommendation:** Accept as proposed.

Sentence of 33 words joining two separate points with "and". Recommend: split.

**Current**

```latex
is not attractive if its error is high for every assignment, and small, noisy differences among the worst assignments are not enough to dismiss SRS.
```

**Proposed**

```latex
is not attractive if its error is high for every assignment. Small, noisy differences among the worst assignments are not enough to dismiss SRS.
```

### HC-18. Unused SEER abbreviation

Source: Writer (article) HC-18 (-1 words)

**Recommendation:** Accept as proposed.

Abbreviation defined but never used again in the body. Recommend: drop "(SEER)". Low priority.

**Current**

```latex
Surveillance, Epidemiology, and End Results (SEER) population estimates
```

**Proposed**

```latex
Surveillance, Epidemiology, and End Results population estimates
```

### HC-19. Split the abstract Results sentence

Source: Writer (article) HC-19 (+0 words)

**Recommendation:** Accept as proposed.

Abstract (not counted in the body limit). The Results sentence runs to about 45 words and joins the two regimes with "while". Recommend: split at the regime change.

**Current**

```latex
had a control-only mean MSE of 0.135, compared with 0.215 for SRS, while under both-arms spillover no design was detectably more accurate than SRS.
```

**Proposed**

```latex
had a control-only mean MSE of 0.135, compared with 0.215 for SRS. Under both-arms spillover, no design was detectably more accurate than SRS.
```

### HC-20. Stale title-page word count

Source: lead instruction (hw-ctj "Not proposed, but noticed")

**Recommendation:** Accept as proposed. The title page says 3,435, but `python3 paper/tools/ctj_wordcount.py` reports 3,473 now (text and headers). If HC-1 to HC-19 are all accepted, the writer measured 3,428 (3,449 without the optional HC-13 and HC-15). Recompute with the script after applying the accepted HC items and enter that count. The abstract count (422) is unchanged by HC-19, which only splits a sentence.

**Current**

```latex
Main-text word count: 3,435.
```

**Proposed**

```latex
Main-text word count: 3,428.
```

## New items not carried forward

- C-5 (Biostat), N6 (Planner) and HN-1 (Writer) became the third block of HW-3.
- C-6 (Biostat) became part of HW-30 block 3.
- HN-20 (Writer) duplicates N5 and is merged into SN-9. HN-2 duplicates C-4 and is merged into SN-4.
- HN-3 (a second pointer to the non-oracle table), HN-7 (a delayed verb) and HN-10 (a pivot sentence and a repeated budget list) are taste edits with no error or material clarity problem. They are dropped and can be revisited in the prose pass.
- Biostat's minor point that Section 4.2 ("surface 1") and Section 12.3 ("median 63%") do not name what the median is taken over was not itemized by any reader.

## Chapter 3 follow-ups for 2i

- HW-6: Chapter 3 line 416 has the same misstated rounding rule ("nearest even whole number"). Fix it with the approved HW-6 prose.
- HW-26 and SN-12 need no chapter change, because `tab:nc-summary-sens` is SI-only.
- Any approved HW, SR or SN wording that has a counterpart in the chapter carries across as listed in the Chapter 3 table of [2g-signoff.md](2g-signoff.md).
- HC items are article-only. HC-2 (rook and queen defined at first use), HC-5 and HC-12 have chapter counterparts worth checking.

## Checks run on this file

A Python script parsed this file and checked every item (93 items, 30 HW, 18 SR, 25 SN, 20 HC).

- **Current blocks.** Every current block of an HW, SR or SN item occurs exactly once in `Supplementary_Information.tex`, or exactly once in the final proposed text of the HW item it builds on (SR-4, SR-6, SR-8, SR-16, SR-17, SR-18). Every HC current block occurs exactly once in `CTJ_Manuscript.tex`. Failures: none.
- **Overlap.** No SN current block overlaps the passage of any HW or SR current block in the SI.
- **Numbers.** No proposed block drops a number from its current block except where the item says why: HW-3 (identity moved to the diagnostics), HW-7 (0.975 moves with the interval formula to Section `estimation`), HW-9 (repeated "58-cluster" access sentence), HW-23 (0.75 to 0.77 restated as 75% to 77%), SR-5 (old sign argument), SR-18 (rounding), SN-20 (repeated sentence, kept in Section 5.3) and HC-20 (word count). No number, citation or reference was invented. The added sentences in HW-7, SR-7 and SN-6 use only values already printed in the SI.
- **Style.** The prose outside the LaTeX blocks has no em-dashes or semicolons.
- **Verifier.** Read, not run: `verify_manuscripts.py` checks SI table bodies against Chapter 3 tables and the budget and sensitivity cells against saved CSVs. The only table-body edits here are the HW-6 table block (rejected) and the HW-26 row labels (SI-only table, numbers unchanged).

## Author decisions

Recorded as the author replies, batch by batch. Nothing is applied until every batch is decided.

- **Style revision (2026-11-25, after the author's A2 reply):** every proposed sentence that newly opened with "Because" or joined two clauses with ", so" was rewritten: SR-7 block 2 (already accepted in A1), HW-24 block 1, SR-1 blocks 1 and 2 (block 1 also explains "count form" in plain words), SR-2 block 2, SR-10 block 2, HW-2 block 1, HW-10 (and the SR-16 current block that quotes it), SN-3, SN-8, HC-2 block 1 and HC-6. HC-7 to HC-10, which turned "Therefore" sentences into "Because" or ", so", are now recommended for decline. Such joins already in the current SI or article text are left alone.
- **Batch A1 (2026-11-25):** accepted as recommended, all ten: HW-17, HW-3 (blocks 1 to 3), HW-7, HW-23 (both blocks), HW-25, HW-13 (both blocks), HW-11, HW-30, HW-6 (prose accepted, table block rejected as recommended) and SR-7.
- **Batch A2 (2026-11-25):** HW-9, HW-20, HW-22, HW-16, HW-14, HW-15 and SR-6 accepted as recommended. HW-24 accepted with block 1 revised to drop the sentence-initial "Because". SR-1 accepted with block 1 revised to explain "count form" in plain words (block 2's ", so" also removed). SR-2 accepted with its ", so" joins removed.
