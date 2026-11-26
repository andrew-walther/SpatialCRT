# Phase 2i sign-off: carrying the 2a–2h, 4e and 3c revisions to Chapter 3

Author: Andrew Walther. Drafted 2026-11-26 against `paper/thesis_chapter/Thesis_Chapter.qmd` as it stands after commit 9bea05f (chapter last aligned at 3e62909/8926344). Nothing is applied yet. Approve, revise or decline each item by ID.

The chapter is the long form and follows the unabridged SI wording (the `\iffullsi` A branch). Curated-only cuts, "of the article" pointers and short list titles stay out. Each text item carries wording the author already approved in 2c, 2d, 2g, 2h or 4e/3c, adapted only for chapter syntax (`@sec-` and `Appendix \ref{}` references, `$...$` math, raw LaTeX blocks) and for the standing style rules. Where approved SI text had a ", so" join, the chapter version splits it and the item says so. Batch 1 lists the places where the chapter now states something the revised SI shows to be wrong or stale, and Batch 2 brings in the new 4e and 3c material.

## Summary

| ID | Recommendation | Reason |
|---|---|---|
| C3T-1 | Accept | Control-only spillover also reaches treated clusters through the multiplier (SR-10) |
| C3T-2 | Accept | Ratios to SRS print with two decimals (SR-18) |
| C3T-3 | Accept | The `fig:vs-srs` interval is a transformed paired interval (SR-13) |
| C3T-4 | Accept | Non-oracle "lowest MSE" holds among the six designs only, and drops "recommendation" (R2, SR-12) |
| C3T-5 | Accept | Rounding rule misstated as "nearest even whole number" (HW-6, NY-11) |
| C3T-6 | Accept | Partial-participation units, 58-college run sizes and MDE (SR-6, HW-10) |
| C3T-7 | Accept | "Did not improve in any setting" overstates the table (SR-14) |
| C3T-8 | Accept | "Nominal size" overstates the calibration (SR-9) |
| C3T-9 | Accept | Sign argument replaced by the product $(\rho-\hat\rho)g$ (SR-5) |
| C3T-10 | Accept | Spatial lag and CAR share a linear predictor, not a mean (SR-4) |
| C3T-11 | Accept | The college mean divides the variance by $n_c$, and effects are in residual SD (SR-6) |
| C3T-12 | Accept | The colleges-needed rule is linear interpolation where power first reaches 80% (SR-8) |
| C3T-13 | Accept | Binary-exposure identification explained correctly (SR-1) |
| C3N-1 | Accept | Point readers to the new regions-and-spread study (3c) |
| C3N-2 | Accept | "Estimating it from trial data is left for future work" is now stale (C1, 4e) |
| C3N-3 | Accept | Saturation-spread result in the Discussion (C2) |
| C3N-4 | Accept | Limitation says other saturation sets "were not studied" (3c) |
| C3N-5 | Accept | Future work says $\tau$ "can be estimated" as if untested (4e) |
| C3N-6 | Accept | Future work promises the study 3c has now done |
| C3N-7 | Accept | Pointer from the reach appendix to 3c, with "North Carolina" (HW-23) |
| C3N-8 | Accept | New Appendix A13 from SI Section 3.2 (N1, SR-2) |
| C3N-9 | Accept | New Appendix A14 from SI Sections 2.2 and 13 (N4, N5, SR-18) |
| C3X-1 | Accept | Parameter-grid caption is stale about SRS (S3) |
| C3X-2 | Accept | SRS rows in `tab:adjacent-pairs` and its text (T2, D2, C3-5) |
| C3F-1 | Accept option A | Remove the verifier's `PENDING_SI` entry, which needs the table in the unabridged SI |
| C3F-2 | Accept | Rebuild the chapter copy of the rejection figure to the thesis rule |
| C3F-3 | Accept | Track the two new chapter figures and protect the rebuilt one in `copy_figures.sh` |
| C3F-4 | Accept | `fig_tau_sensitivity_6design_combined.pdf` still prints the band note (F5, F8) |
| C3F-5 | Accept | `fig_power_vs_n.pdf` prints minus signs in "college-level" (F8) |
| C3F-6 | Accept | `refined_allocation_risk.pdf` prints "Estimated q90" (F3) |
| C3F-7 | Accept | The chapter's NC example map still has numbered panel titles (F4) |
| C3T-14 | Accept | Define queen and rook at first use (HC-2) |
| C3T-15 | Accept | Row-standardization without "rather than" (NY-4) |
| C3T-16 | Accept | Direct effect without "not an average total effect" (NYC-2) |
| C3T-17 | Accept | Poisson ties sentence without a colon (HW-1) |
| C3T-18 | Accept | Two colons removed in the surface paragraph (SN-17) |
| C3T-19 | Accept | Deterministic and tie-breaking sentences (HW-2) |
| C3T-20 | Accept | Oracle model paragraph, the `lagsarlm` refit and the "Spill" label (HW-3, SN-14) |
| C3T-21 | Accept | Grid simulation size, "for the eight designs" (HW-4) |
| C3T-22 | Accept | Monte Carlo standard errors without a colon, with the pooling detail (SR-13) |
| C3T-23 | Accept | Source of the control-only increase (NY-18) |
| C3T-24 | Accept | Rook comparison names the Holm adjustment (SN-11) |
| C3T-25 | Accept | Worst-case MSE scoped to the six designs (SR-11, HW-17) |
| C3T-26 | Accept | Shared draws across effect sizes (NY-19) |
| C3T-27 | Accept | Uncertainty bands and rook Checkerboard in the caption (HW-13) |
| C3T-67 | Accept | Treatment-effect sweep: "both rounded to 0.160" and the six-design scope (HW-14) |
| C3T-28 | Accept | Grid noise-share sentence without a colon (HW-19) |
| C3T-68 | Accept | Grid allocation-risk numbers and their scope (HW-19 block 1) |
| C3T-29 | Accept | Observed inputs versus planning assumptions (SN-6) |
| C3T-30 | Accept | County-to-college assignments (SN-25) |
| C3T-31 | Accept | Weights matrix as a proxy (NY-8) |
| C3T-32 | Accept | Budget caption (NY-2) |
| C3T-71 | Accept | Colon before the case-study outcome equation (HW-7) |
| C3T-33 | Accept | Case-study target and spillover regimes (NY-12, SR-10) |
| C3T-34 | Accept | Case-study settings and the hypothetical association (HW-8, NY-13) |
| C3T-35 | Accept | Case-study precision paragraph (HW-8, SN-20) |
| C3T-36 | Accept | Case-study rank tests (NY-16, HW-30) |
| C3T-70 | Accept | Units of the case-study MSE (SN-8) |
| C3T-37 | Accept | High Incidence Focus by year (R4, HW-22) |
| C3T-38 | Accept | High Incidence Focus and randomization (NY-21, SN-13) |
| C3T-39 | Accept | "Eight designs and SRS" in two results paragraphs (HW-30) |
| C3T-40 | Accept | Population denominators (NY-3) |
| C3T-41 | Accept | "In North Carolina it differs little" (HW-23) |
| C3T-42 | Accept | Tail summaries and reused assignments (HW-24, NY-22) |
| C3T-43 | Accept | Downside table caption (HW-20) |
| C3T-44 | Accept | Tail uncertainty paragraph (HW-25) |
| C3T-45 | Accept | Prespecified mean-rank rule (NY-23) |
| C3T-46 | Accept | Bayesian CAR diagnostics and the spatial lag model (SR-17, HW-29) |
| C3T-47 | Accept | Colleges needed, with the conditions and residual units (SR-7, SR-6) |
| C3T-48 | Accept | High Incidence Focus in the Discussion (NY-21) |
| C3T-49 | Accept | "Identifiable but imprecisely estimated" (I4, C3-2) |
| C3T-50 | Accept | What the rerun supplies (HW-15) |
| C3T-51 | Accept | Case-study reproducibility paragraph (HW-9) |
| C3T-52 | Accept | Regional summaries (NY-9, NY-10) |
| C3T-53 | Accept | Bayesian chain diagnostics and "eight designs and SRS" (SR-17, HW-30) |
| C3T-54 | Accept | Refinement bookkeeping (SN-20, NY-14) |
| C3T-55 | Accept | Selection-bias check and the Jaccard index (HW-9, NY-15) |
| C3T-56 | Accept | Stale "given here to three decimal places" (S2) |
| C3T-57 | Accept | SRS placed in the adjacent-pair order (T1) |
| C3T-58 | Accept | Rook block counts scoped to the six designs (SR-11) |
| C3T-69 | Accept | Rook 50-unit tests with SRS in the adjacent-pair order (2c T3) |
| C3T-59 | Accept | Rook Checkerboard estimates (NY-20, HW-16) |
| C3T-60 | Accept | Run-sizes head on its own line (SH-1, HW-21) |
| C3T-61 | Accept | "Annual confirmations" (SN-22) |
| C3T-62 | Accept | Illustration of the individual-level effect (SR-7) |
| C3T-63 | Accept | Colleges-needed caption (HW-27, SR-6) |
| C3T-64 | Accept | Interval calibration paragraph (HW-28) |
| C3T-65 | Accept | Case-study Friedman tests (NY-25, HW-30) |
| C3T-66 | Accept | Spatial-model results get their own heading and three paragraphs (HW-29, SR-17) |

Item IDs follow the order of the batches. Within each batch, items run in chapter order. C3T-67 to C3T-71 were added by the independent review and sit in chapter order inside their batches.

## Batch 1 · Stale or wrong statements in the chapter (priority)

### C3T-1. Control-only spillover and the spatial multiplier (line 87)

**Recommendation:** Accept. The chapter says only untreated clusters respond under control-only spillover, but with $\rho > 0$ the multiplier carries that contribution to treated clusters too. From SR-10 and HC-15 (the elementwise gloss).

**Current**

```qmd
where $\circ$ denotes the elementwise product, which multiplies the two vectors cluster by cluster. The share of treated neighbors $\mathbf{W}\mathbf{Z}$ is an exposure mapping, a function that assigns each cluster an exposure from the full assignment vector [@aronow_estimating_2017; @forastiere_identification_2021]. Under the control-only regime, only untreated clusters respond to the treatment status of their neighbors. This is classical contamination: an untreated cluster's outcome rises with the share of its neighbors that are treated. Under the both-arms regime, every cluster responds to its neighbors' treatment regardless of its own assignment. Therefore, neighboring treatment adds to a treated cluster's own effect.
```

**Proposed**

```qmd
where $\circ$ multiplies the two vectors cluster by cluster. The share of treated neighbors $\mathbf{W}\mathbf{Z}$ maps the full assignment to each cluster's exposure [@aronow_estimating_2017; @forastiere_identification_2021]. Control-only spillover represents classical contamination. The exposure term enters the outcome equation only for untreated clusters. An untreated cluster's outcome rises with the share of its neighbors that are treated. The spatial multiplier can still propagate that contribution to treated clusters. Under both-arms spillover, the exposure term enters the equation for treated and untreated clusters. Neighboring treatment then adds to a treated cluster's own effect.
```

### C3T-2. Ratios to SRS with two decimals (line 283)

**Recommendation:** Accept. The chapter prints 0.631 and 1.083, against the two-decimal rule for ratios to SRS. SR-18 rounded the same numbers in the SI. This is the only number change in the list.

**Current**

```qmd
Incidence-Guided Saturation's MSE was 0.631 times that of SRS under control-only spillover and 1.083 times under both arms.
```

**Proposed**

```qmd
Incidence-Guided Saturation's MSE was 0.63 times that of SRS under control-only spillover and 1.08 times under both arms.
```

### C3T-3. The interval in `fig:vs-srs` (line 293)

**Recommendation:** Accept. The plotted interval is the paired $t_{49}$ interval for the difference, transformed with the SRS mean held fixed. Calling it "the $t_{49}$ 95% interval" misdescribes it. From SR-13 block 3.

**Current**

```qmd
The $t_{49}$ 95\% interval is shown for each point.
```

**Proposed**

```qmd
Each point has an approximate 95\% interval, obtained by transforming the paired $t_{49}$ interval for the mean difference with the SRS mean held fixed.
```

### C3T-4. Non-oracle results: scope and the word "recommendation" (lines 344, 613, 1183)

**Recommendation:** Accept. The block counts and rankings behind "lowest MSE" cover the six detailed designs, not SRS, and line 1183 still speaks of "the recommendation", which the no-single-design rule removed. From R2 and SR-12.

**Current (edit 1)**

```qmd
although Incidence-Guided Saturation still had the lowest MSE under both neighbor types (Appendix \ref{app:oracle-nonoracle}).
```

**Proposed (edit 1)**

```qmd
although Incidence-Guided Saturation still had the lowest MSE of the six detailed designs under both neighbor types (Appendix \ref{app:oracle-nonoracle}).
```

**Current (edit 2)**

```qmd
but Incidence-Guided Saturation still had the lowest MSE under both neighbor types.
```

**Proposed (edit 2)**

```qmd
but Incidence-Guided Saturation still had the lowest MSE of the six detailed designs under both neighbor types.
```

**Current (edit 3)**

```qmd
However, where a design's non-oracle MSE was lower, its estimate was biased and its intervals under-covered, which means the lower MSE does not support valid inference. These results do not change the recommendation: Incidence-Guided Saturation also had the lowest non-oracle MSE under both neighbor types (0.087 under queen contiguity and 0.086 under rook contiguity).
```

**Proposed (edit 3)**

```qmd
Among the six detailed designs, Incidence-Guided Saturation had the lowest pooled MSE under both models and neighbor types. Its non-oracle MSE was 0.087 under queen contiguity and 0.086 under rook contiguity. Designs with lower MSE under the non-oracle model still had bias and intervals that covered less often than nominal.
```

### C3T-5. Spatial-blocking rounding rule and budgets (line 416)

**Recommendation:** Accept. "Rounded to the nearest even whole number" is wrong (2g HW-6 follow-up). The rule rounds to the nearest whole number and sends an exact half to the even one. Also from NY-11 ("rather than modify the designs") and SN-9 ("On the NC map").

**Current (edit 1)**

```qmd
On the NC map, $2\times2$ Blocking becomes spatial blocking.
```

**Proposed (edit 1)**

```qmd
On the North Carolina map, $2\times2$ Blocking becomes spatial blocking.
```

**Current (edit 2)**

```qmd
Half of each group, rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), is treated at random. Therefore, treatment budgets, the numbers of treated clusters, differ across designs. We report them (Table \ref{tab:nc-budgets}) rather than modify the designs to force equal budgets.
```

**Proposed (edit 2)**

```qmd
The rule randomly treats half of each group, rounded to the nearest whole number, with an exact half rounded to the even number (1.5 to 2 and 0.5 to 0). Therefore, treatment budgets, the numbers of treated clusters, differ across designs. We kept each design's own rule and report its budget in Table \ref{tab:nc-budgets}.
```

### C3T-6. Partial-participation units and run sizes (lines 470, 472)

**Recommendation:** Accept. Three corrections. Effects are in college-level residual standard deviations before spatial propagation (SR-6), single-assignment designs used 1,000 outcomes at 58 colleges (the chapter says every setting used 100), and the shared draws serve to reduce Monte Carlo noise in the comparisons (HW-10). Also defines MDE where it is first needed and uses "the eight designs and SRS" (HW-30).

**Current (edit 1)**

```qmd
and the effect $\beta$ took nine values from 0 to 2 in units of the college-level outcome standard deviation.
```

**Proposed (edit 1)**

```qmd
and the effect $\beta$ took nine values from 0 to 2 in units of the college-level residual standard deviation before spatial propagation.
```

**Current (edit 2)**

```qmd
Settings crossed these with $\rho \in \{0, 0.5\}$, both spillover regimes, the nine designs and the four planning years.
```

**Proposed (edit 2)**

```qmd
Settings crossed these with $\rho \in \{0, 0.5\}$, both spillover regimes, the eight designs and SRS, and the four planning years.
```

**Current (edit 3)**

```qmd
With all 58 colleges, each setting used 100 allocations and 100 outcomes, as in the main case study, and reproduced that case study's MSE and coverage within Monte Carlo error. All effect sizes share the same subsets, allocations and noise, which makes power curves smooth in $\beta$. For each setting we estimated three quantities with the corrected $t$ test above: power, the number of colleges needed for 80\% power, and the smallest effect detectable with 80\% power at each $N$.
```

**Proposed (edit 3)**

```qmd
With all 58 colleges, randomized designs used 100 assignments and 100 outcomes per assignment, as in the main case study, and a design with one possible assignment used 1,000 outcomes for that assignment. The resulting MSE and coverage matched the main case study within Monte Carlo error. We used the same subsets, assignments and noise across effect sizes to reduce Monte Carlo noise in the comparisons. For each setting we estimated power with the corrected $t$ test. From the power curves we obtained the number of colleges needed for 80\% power and the minimum detectable effect (MDE), the smallest effect detectable with 80\% power at each $N$.
```

### C3T-7. Checkerboard and 2×2 Blocking against SRS on the downside criteria (line 503)

**Recommendation:** Accept. "Did not improve on SRS on either criterion in any setting" overstates what the table shows, which is study-by-regime averages (and $2\times2$ Blocking's case-study both-arms mean ratio is 1.00). From SR-14, with "North Carolina" and "case study" for bare "NC" (SN-9).

**Current (edit 1)**

```qmd
on the grid and at the refined NC settings.
```

**Proposed (edit 1)**

```qmd
on the grid and at the refined North Carolina settings.
```

**Current (edit 2)**

```qmd
(0.62 of SRS on the grid and 0.91 in NC). Under both-arms spillover, its mean was 3.80 and 3.69 times that of SRS, and it treated fewer clusters. High Incidence Focus has a single allocation in NC.
```

**Proposed (edit 2)**

```qmd
(0.62 of SRS on the grid and 0.91 in the case study). Under both-arms spillover, its mean was 3.80 and 3.69 times that of SRS, and it treated fewer clusters. High Incidence Focus has a single allocation in the case study.
```

**Current (edit 3)**

```qmd
However, on the grid its both-arms mean was 3.81 times that of SRS, and its NC result depends on the incidence model, as the matched sensitivity analysis below shows. Balanced Quartiles and Balanced Halves matched SRS on both criteria. Checkerboard and $2\times2$ Blocking did not improve on SRS on either criterion in any setting.
```

**Proposed (edit 3)**

```qmd
However, in the grid allocation-risk study its both-arms mean MSE was 3.81 times that of SRS, and its case study result depended on the incidence model, as the matched sensitivity analysis below shows. Balanced Quartiles and Balanced Halves matched SRS on both criteria. For Checkerboard and $2\times2$ Blocking, mean and worst-decile MSE were at least as large as those for SRS in each study-by-regime average.
```

### C3T-8. What the calibration supports (line 581)

**Recommendation:** Accept. The simulations show a rejection rate close to 5%, not an exact nominal size. From SR-9 (as revised by the author).

**Current**

```qmd
Therefore, the colleges-needed numbers above come from a test with the nominal size.
```

**Proposed**

```qmd
Therefore, the sample-size estimates use a test with a rejection rate close to the nominal 5\% level in these simulations.
```

### C3T-9. The term that grows with $\beta$ (line 776)

**Recommendation:** Accept. A positive $\rho-\hat\rho$ on average and $g > 0$ do not by themselves fix the sign of the term's mean, which depends on the product $(\rho-\hat\rho)g$. From SR-5. The variance figures in the next sentences stay.

**Current**

```qmd
With $\rho - \hat\rho$ positive on average, the last term adds a positive bias and extra variance that both grow with $\beta$ for any design with $g > 0$.
```

**Proposed**

```qmd
The last term depends on the product $(\rho-\hat\rho)g$ as well as on $\beta$.
```

### C3T-10. Spatial lag and CAR models share a linear predictor (lines 851, 859)

**Recommendation:** Accept. The two models do not have the same mean. The spatial lag outcome has mean $(\mathbf{I}-\rho\mathbf{W})^{-1}\boldsymbol{\mu}$. From SR-4, with "the eight designs and SRS" (HW-30). Edit 2 also replaces the "Spill" label, as C3T-20 does in the main text.

**Current (edit 1)**

```qmd
In the reverse direction, we generated data from a CAR model with the same mean and fit them with both models.
```

**Proposed (edit 1)**

```qmd
We also generated data from a CAR model with the same linear predictor and fit both models to those data.
```

**Current (edit 2)**

```qmd
Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\text{Spill}$, where $\text{Spill} = \mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography.
```

**Proposed (edit 2)**

```qmd
The two spatial families use the linear predictor $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and use the 58-college queen geography.
```

**Current (edit 3)**

```qmd
We fit it by maximum likelihood with the `spautolm` function of the R package `spatialreg` [@bivand_r_2022].
```

**Proposed (edit 3)**

```qmd
We fit it by maximum likelihood with the `spautolm` function of the R package `spatialreg` [@bivand_r_2022]. Thus the spatial lag outcome has conditional mean $(\mathbf{I}-\rho\mathbf{W})^{-1}\boldsymbol{\mu}$ given the assignment, whereas the CAR outcome has conditional mean $\boldsymbol{\mu}$.
```

**Current (edit 4)**

```qmd
These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes and the nine designs, with 500 fits per setting.
```

**Proposed (edit 4)**

```qmd
These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```

### C3T-11. From college-level to individual-level effects (line 855)

**Recommendation:** Accept. Averaging divides the variance of the individual part by $n_c$, not the part itself, and the college-level unit is the residual standard deviation before spatial propagation. From SR-6, with the approved sentence on what the translation assumes.

**Current (edit 1)**

```qmd
The trial outcome is a score measured on each trained professional and averaged within each college.
```

**Proposed (edit 1)**

```qmd
For planning, suppose the education outcome is a score measured on each trained professional and averaged within each college.
```

**Current (edit 2)**

```qmd
leaves the college part unchanged and divides the individual part by $n_c$.
```

**Proposed (edit 2)**

```qmd
leaves the college part unchanged and divides the variance of the individual part by $n_c$.
```

**Current (edit 3)**

```qmd
A college-level effect of $\beta$ standard deviations then corresponds to an individual-level standardized effect $d = \beta\sigma_c$.
```

**Proposed (edit 3)**

```qmd
A structural effect of $\beta$ college-level residual standard deviations, before spatial propagation, then corresponds to an individual-level standardized effect $d = \beta\sigma_c$. This treats the sampling error of each college mean as the residual $\boldsymbol{\varepsilon}$ in the case-study outcome model (@sec-nc-outcome), which the spatial multiplier also propagates when $\rho > 0$.
```

### C3T-12. How the number of colleges needed was computed (line 1361)

**Recommendation:** Accept. "Interpolated monotonically" is not what the code does. It interpolates linearly where estimated power first reaches 80%, and a year-averaged count can be exceeded in single years. From SR-8. The chapter's own statement that some dips exceed Monte Carlo error stays, without the colon and "rather than".

**Current**

```qmd
The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not always monotone in $N$: some drops in power between successive values of $N$ are larger than Monte Carlo error, because the participant subsets differ between values of $N$. We report the dips rather than smooth them.
```

**Proposed**

```qmd
Within each year, we used linear interpolation between the simulated sample sizes where estimated power first reached 80\%. When we obtained an estimate for every year, we averaged those estimates and rounded up. Individual years can require more colleges than the rounded average. Estimates at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not always monotone in $N$. Some drops in power between successive values of $N$ are larger than Monte Carlo error, because each sample size uses its own participant subsets. We retain these dips in the simulated power curves.
```

### C3T-13. Identification under the binary spillover form (lines 1519, 1525)

**Recommendation:** Accept. The old reason for Isolation Buffer's failure is incomplete (its binary exposure equals $\mathbf{1}-\mathbf{Z}$), and the ratio among identifiable fits needs its scope. The table note now says what "n.i." marks. From SR-1.

**Current (edit 1)**

```qmd
With the binary form, Isolation Buffer is never identifiable, because every control has a treated neighbor, and under control-only spillover no design was identifiable in every fit. The saturation designs were identifiable most often there (in 89.7\% and 91.7\% of fits) and, among identifiable fits, had MSE 0.70 to 0.71 times that of SRS.
```

**Proposed (edit 1)**

```qmd
Under binary exposure, Isolation Buffer assigns every treated cluster zero exposure and every control cluster exposure one. Therefore, the model cannot separate $\beta$ from spillover in either regime. Under control-only spillover, $\beta$ was unidentified in some fits under every assignment rule, including SRS. The two saturation designs identified $\beta$ most often, in 89.7\% and 91.7\% of fits. Among the identifiable fits for each design and for SRS, their MSEs were 0.70 to 0.71 times that of SRS. These ratios describe the subsets that were identifiable under each rule.
```

**Current (edit 2)**

```qmd
n.i.: not identifiable in every fit (the spillover covariate is aliased with the intercept and treatment in some or all allocations).}
```

**Proposed (edit 2)**

```qmd
An entry of n.i.\ means that the treatment coefficient $\beta$ was unidentifiable in at least one sampled assignment.}
```

## Batch 2 · New 4e and 3c material

Placement decision inside this batch. The new material becomes two appendices at the end, A13 (incidence coefficient, from SI Section 3.2) and A14 (variants of incidence guidance, from SI Sections 2.2 and 13). Appending keeps A1–A12 and every existing `\ref` unchanged. The alternative is to insert A13 after A2 and A14 after A6 and renumber A3–A12 by hand. Say so if you prefer it. Every sentence below comes from the applied SI text (N1, N4, N5 as revised by SR-2 and SR-18), with these chapter adaptations. The math uses `$...$`, references and figure paths follow the chapter (no `../`), short captions follow the chapter's style, and the SI's ", so" joins become two sentences.

### C3N-1. Pointer to the regions-and-spread study (line 101)

**Recommendation:** Accept. The design section is where a reader asks whether four quadrants and 80/20 were tuned. From N4.

**Current**

```qmd
Appendix \ref{app:eight-design} describes both and gives the comparison of all eight designs with SRS.
```

**Proposed**

```qmd
Appendix \ref{app:eight-design} describes both and gives the comparison of all eight designs with SRS. Appendix \ref{app:igs-sensitivity} varies the number of regions and the saturation spread of Incidence-Guided Saturation on the grid.
```

### C3N-2. Estimating $\tau$ is no longer future work (line 162)

**Recommendation:** Accept (priority). The chapter still says estimating $\tau$ "is left for future work", which 4e has done. Carries the approved C1 sentence.

**Current**

```qmd
The value may need adjusting for a particular outcome, and estimating it from trial data is left for future work.
```

**Proposed**

```qmd
The value may need adjusting for a particular outcome. After a trial, the same model estimates it, less precisely when treatment follows incidence (Appendix \ref{app:incidence-coefficient}).
```

### C3N-3. Saturation spread in the Discussion (line 607)

**Recommendation:** Accept. Carries the approved C2 sentence to the paragraph on the shared saturation structure, and replaces "the NC map" (SN-9).

**Current**

```qmd
a few percentage points of incidence share on the grid and less than one on the NC map, whose regions have similar incidence.
```

**Proposed**

```qmd
a few percentage points of incidence share on the grid and less than one on the North Carolina map, whose regions have similar incidence. On the grid, widening the quadrants' saturation spread beyond 80/20 lowered control-only MSE but raised both-arms MSE (Appendix \ref{app:igs-sensitivity}).
```

### C3N-4. Limitation on saturation sets (line 619)

**Recommendation:** Accept (priority). "Other saturation sets and numbers of regions were not studied" is stale after 3c. The replacement keeps the scope caveats N5 states (equal-size partitions, linearly spaced targets, grid only).

**Current**

```qmd
The saturation designs used one set of saturations, 80\%, 60\%, 40\% and 20\%, and four regions. Other saturation sets and numbers of regions [@baird_optimal_2018] were not studied, and the saturation advantage under control-only spillover may depend on them.
```

**Proposed**

```qmd
The main comparison used one set of saturations, 80\%, 60\%, 40\% and 20\%, and four regions. A grid sensitivity analysis varied the number of equal-size regions and the spread of the saturations [@baird_optimal_2018] (Appendix \ref{app:igs-sensitivity}). Widening the spread up to 90/10 lowered control-only MSE for every number of regions and raised both-arms MSE for halves and quadrants. That analysis covered equal-size partitions of the grid with linearly spaced targets, and the North Carolina map was not part of it.
```

### C3N-5. Future work on $\tau$ (line 635)

**Recommendation:** Accept (priority). The chapter describes estimating $\tau$ after a trial as an untested possibility. From N1 and C1.

**Current**

```qmd
After a trial, it can be estimated from the trial data by including incidence as a covariate in the fitted model, and earlier trials or observational data on the same outcome could inform its value at the design stage.
```

**Proposed**

```qmd
After a trial, the fitted model estimates it, with precision that depends on the design (Appendix \ref{app:incidence-coefficient}). Earlier trials or observational data on the same outcome could inform its value at the design stage.
```

### C3N-6. Future work on saturations and regions (line 639)

**Recommendation:** Accept (priority). The paragraph promises the study that 3c has now run. The replacement names what 3c left open, using only its stated scope.

**Current**

```qmd
Incidence-Guided Saturation outperformed SRS under control-only spillover with one set of saturations and four regions. Future work will vary the saturations and the number of regions to find whether another configuration estimates the intervention effect more accurately.
```

**Proposed**

```qmd
Incidence-Guided Saturation outperformed SRS under control-only spillover with four regions and saturations from 80\% to 20\%. The grid sensitivity analysis of regions and saturation spread used equal-size partitions with linearly spaced targets (Appendix \ref{app:igs-sensitivity}). Irregular regions, such as those on the North Carolina map, and targets that are not linearly spaced remain to be studied.
```

### C3N-7. Reach appendix: pointer and "North Carolina" (line 1505)

**Recommendation:** Accept. The SI's reach paragraph now ends with this pointer to 3c and says "North Carolina regions" (HW-23).

**Current**

```qmd
and the NC regions lie at the low end of that spread.
```

**Proposed**

```qmd
and the North Carolina regions lie at the low end of that spread. Appendix \ref{app:igs-sensitivity} shows how the number of regions and the saturation targets affect targeting and MSE on the grid.
```

### C3N-8. New Appendix A13: estimating the incidence coefficient (after A12)

**Recommendation:** Accept. SI Section 3.2 as applied (N1 with SR-2's definition of $\tau$). Four ", so" joins are split. The figure needs the rebuilt chapter copy (C3F-2) before the chapter renders to the thesis rule.

**Current**

````qmd
Simple Random Sampling & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & n.i. \\
\bottomrule
\end{tabular}
\end{table}
```
````

**Proposed**

````qmd
Simple Random Sampling & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & n.i. \\
\bottomrule
\end{tabular}
\end{table}
```

```{=latex}
\FloatBarrier
```

## A13. Estimating the Incidence Coefficient After a Trial {.unnumbered .unlisted}

\applabel{A13}{app:incidence-coefficient}The incidence coefficient $\tau$ in Equation \ref{eq:outcome} gives the change in the outcome's linear predictor for a one-unit increase in scaled incidence $\mathbf{X}$, the full width of its 0-to-1 range, holding treatment and spillover fixed. The change is in residual standard deviations because $\sigma = 1$, and is measured before spatial propagation. At $\tau = 1$, a one-unit incidence contrast contributes as much to the linear predictor as the primary intervention effect $\beta = 1$. At $\tau = 0$ incidence makes no direct contribution to that predictor. The coefficient can exceed 1 or be negative.

Before outcomes are observed, investigators specify $\tau$ for the planning model. After a trial, the same spatial lag model can estimate $\tau$ and $\beta$ from the observed outcomes and the incidence used for allocation. Both grid models include $\mathbf{X}_k$ (@sec-estimation), and every fit provides both estimates. The oracle estimate describes the association of prior incidence with the outcome after adjustment for treatment and spillover. The non-oracle model adjusts for treatment and omits spillover. Incidence is not randomized. Therefore, these coefficients describe adjusted associations between incidence and the outcome.

To see how well $\tau$ can be estimated, and whether its value changes the design comparison, we repeated the full grid benchmark at $\tau = 0$, 0.25, 0.5, 1 and 2. The incidence surfaces, assignments and outcome noise were the same at every value of $\tau$. Only the incidence term of the outcome changed. The 95\% interval for $\tau$ has the same $t$-based, degrees-of-freedom-corrected form as the interval for $\beta$, and power is the proportion of intervals that exclude zero, a two-sided test at level 0.05. At $\tau = 1$ the rerun reproduced every stored intervention result of the benchmark exactly.

With the oracle model, the largest absolute bias of $\hat\tau$ across designs, incidence configurations, neighbor types, regimes and values of $\tau$ was 0.023. The intervals covered $\tau$ in 93.8--95.8\% of fits, and at $\tau = 0$ they excluded zero in 4.3--6.2\% of fits. Without the spillover term, bias reached 0.137, because omitting spillover can shift the fitted association with incidence.

The design affected the precision of $\hat\tau$ far more than the value of $\tau$ did. On the queen grid with independent uniform incidence, the mean model standard error of $\hat\tau$ was 0.667--0.670 for High Incidence Focus and 0.336--0.337 for SRS, and every other design was between 0.334 and 0.341. At $\tau = 1$, the test rejected in 29.3--30.3\% of High Incidence Focus trials and in 81.1--81.9\% of SRS trials across the two spillover regimes (Figure \ref{fig:incidence-rejection}). High Incidence Focus reached that level of power only at $\tau = 2$, where it rejected in 81.3--81.5\% of trials. It had the largest standard error of any design in every incidence configuration, neighbor type and regime. High Incidence Focus treats the clusters above the median incidence. Its mean absolute correlation between $\mathbf{Z}$ and $\mathbf{X}$ was 0.866, against 0.079--0.080 for SRS, and little variation in $\mathbf{X}$ remains once treatment is in the model. For uniform $\mathbf{X}$, splitting treatment at the median leaves one quarter of the variance of $\mathbf{X}$ within the arms. In a simple regression with residual standard deviation 1, this doubles the standard error of the incidence coefficient from $\sqrt{12/100} = 0.35$ to 0.69. Incidence-Guided Saturation also uses incidence, but by region and with every region partly treated. Its correlation was 0.109--0.110 and its standard error 0.339--0.341, close to that of SRS.

Changing $\tau$ barely moved $\hat\beta$. Over the five values of $\tau$, no design's mean MSE for $\beta$ changed by more than 3.4\% from its value at $\tau = 1$ in any incidence configuration, neighbor type or regime. Both fitted models include $\mathbf{X}$. Changing $\tau$ leaves the fitted intervention coefficient unchanged when $\rho$ is held at its true value. Error in estimating $\rho$ can pass part of the change into $\hat\beta$ (Appendix \ref{app:technical-notes}). The spillover regime still decided which designs led. Under control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest MSE in every comparison. Under both-arms spillover, the best design's MSE was at most 5\% below that of SRS. The lowest-MSE design changed with $\tau$ in 6 of the 20 comparisons, and in 5 of these the two saturation designs, whose MSEs were nearly equal, exchanged first place.

The oracle results assume the correct form of spillover, and the non-oracle bias shows that $\hat\tau$ depends on it. The grid used one number of clusters, ten incidence surfaces per configuration and nonnegative values of $\tau$, and $\hat\tau$ is on the scale of incidence after it was scaled to the unit interval. The primary case-study model omits incidence because the education outcome is not expected to depend on it. Therefore, these results give no estimate of $\tau$ for that trial.

```{=latex}
\begin{figure}[!htb]
\centering
\includegraphics[width=\linewidth]{figures/fig_incidence_rejection_queen.pdf}
\caption[Two-sided rejection rate of the interval for the incidence coefficient $\tau$ under queen contiguity]{Two-sided rejection rate of the 95\% interval for the incidence coefficient $\tau$ under queen contiguity with the oracle model, by incidence configuration (rows) and spillover regime (columns), averaged over the five values of $\beta$. At $\tau = 0$ the rate is the test size, and the dashed line marks 0.05. High Incidence Focus has the lowest rejection rate at every $\tau > 0$, and the other designs and SRS overlap.}
\label{fig:incidence-rejection}
\end{figure}

\FloatBarrier
```
````

### C3N-9. New Appendix A14: variants of incidence guidance (end of the appendix)

**Recommendation:** Accept. SI Sections 2.2 (N4) and 13 (N5 with SR-18's rook numbers and D3's two-decimal ratios), as applied. Three ", so" joins in the methods are split. The SI's cross-reference to its own Section 2.2 becomes "described above".

**Current**

````qmd
```{=latex}
% Restore the table-of-contents depth for the chapters that follow.
````

**Proposed**

````qmd
## A14. Variants of Incidence Guidance on the Grid {.unnumbered .unlisted}

\applabel{A14}{app:igs-sensitivity}Incidence-Guided Saturation fixes two choices, four quadrants and saturations from 80\% down to 20\%. To see how much each matters, we varied the number of regions and the spread of the saturations on the grid. The grid was divided into 2 halves, 4 quadrants, 10 rectangles or 25 blocks of equal size. With $R$ regions and spread $s$, the region ranked $r$ by mean incidence received the target saturation $0.5 + s\{1 - 2(r-1)/(R-1)\}$. The highest-ranked region was then treated at $50\% + s$ and the lowest at $50\% - s$. We used $s = 0.1$, 0.2, 0.3, 0.4 and 0.5, which give highest and lowest saturations of 60/40, 70/30, 80/20, 90/10 and 100/0 percent. Four quadrants at 80/20 is the original design.

The 100/0 setting is a boundary case, included only to show whether going beyond 90/10 helps or hurts. It treats one region fully and leaves another untreated. For quadrants its targets are 100\%, 67\%, 33\% and 0\%, or 25, 17, 8 and 0 treated clusters. Targets were rounded by largest remainders so that every assignment treated exactly 50 of the 100 clusters, and ties in regional mean incidence were broken at random. Each guided variant has an unguided twin that gives the same targets to the regions in random order. The pair isolates the effect of incidence guidance. The unguided twin of the original design is Saturation Quadrants. With 25 blocks at 60/40, every block has two treated clusters. Both versions are then $2\times2$ Blocking, and the pair serves as a null control.

The 38 new variants, the original design, Saturation Quadrants and SRS were run on the full grid of @sec-simulation, with the same incidence surfaces and outcome noise as the other designs, 25 assignments per surface and both fitted models. We report queen contiguity and the oracle model, pooled over the five incidence configurations and the five values of $\beta$ within each spillover regime. Each comparison uses the 50 paired configuration-by-surface differences in MSE, with a two-sided Wilcoxon signed-rank test and Holm adjustment within each family. The families compare each design with SRS (40 comparisons), each guided variant with the original design (19) and each guided variant with its twin (20). For the figure, the paired 95\% $t_{49}$ interval for the mean difference is divided by the comparator's MSE and shifted by 1, which puts it on the ratio scale. Every 100/0 assignment gave full-rank design matrices, and no fit of a new variant was aliased or failed.

### Results {.unnumbered .unlisted}

Figure \ref{fig:igs-sensitivity} shows the MSE of each variant described above relative to SRS. Pooled over the five values of $\beta$, SRS had MSE 0.227 under control-only spillover and 0.043 under both-arms spillover. These pooled ratios differ from those at $\beta = 1$ in @sec-vs-srs. Under control-only spillover, widening the guided saturations from 60/40 to 90/10 lowered MSE for every number of regions, and at each of these spreads fewer, larger regions had lower MSE. At 90/10 the guided ratios to SRS were 0.55 for halves, 0.59 for quadrants, 0.77 for ten regions and 1.05 for 25 blocks, against 0.67 for the original design. Under both-arms spillover, wider spreads raised MSE for halves and quadrants. Guided quadrants rose from 1.00 times SRS at 60/40 to 1.22 at 90/10, and halves from 0.98 to 1.94. With ten regions or 25 blocks, MSE stayed between 1.00 and 1.22 times that of SRS from 70/30 to 100/0. Some variants traded less between the regimes than the original design. Quadrants at 70/30 and ten regions at 80/20 had ratios of 0.84 and 0.89 under control-only spillover and 1.02 under both-arms spillover, where neither differed detectably from SRS.

Moving from 90/10 to the 100/0 boundary had different effects by number of regions under control-only spillover. MSE rose for halves, from 0.55 to 0.88 times that of SRS for the guided version, and the bias of $\hat\beta$ rose from 0.086 to 0.161. Quadrants changed little (0.59 at both spreads), while MSE fell for ten regions (0.77 to 0.68) and 25 blocks (1.05 to 0.86). At 100/0, both 25-block versions had detectably lower MSE than SRS for the first time. Under both-arms spillover, MSE rose at 100/0 for every number of regions in both versions. The largest rise was for halves, whose MSE reached 8.13 times that of SRS. When one half is fully treated, $\mathbf{W}\mathbf{Z}$ differs from $\mathbf{Z}$ only at the 20 clusters along the dividing line, and the two columns have correlation 0.941 under queen contiguity. The model then has full rank but separates treatment from spillover poorly. At 100/0, coverage of the $t$-based interval was between 0.932 and 0.953, lowest for halves under control-only spillover. The changes from 90/10 to 100/0 are descriptive, because the tests compare each design with SRS, the original design or its twin.

Incidence guidance changed where treatment went more than how accurate the estimate was. At the original 80/20 quadrants, guided and unguided MSE were close (ratios 1.00 under control-only and 1.04 under both-arms spillover). No guided variant had detectably lower MSE than its unguided twin, and every detectable difference favored the unguided version. At 100/0 under control-only spillover, guided MSE was 1.04 to 1.10 times that of the twin, and all four differences remained after Holm adjustment. In the targeting check (queen contiguity, $\rho = 0$, $\psi = 0.5$, control-only spillover, averaged over the five incidence configurations), the original design treated 52.7\% of total incidence. The share grew with the spread and the number of regions, to 57.7\% for 25 blocks at 90/10 and 53.6--60.1\% at 100/0, while the unguided versions and SRS treated about 50\%.

Rook contiguity gave the same directions from 90/10 to 100/0. Both versions of halves at 100/0 had slightly higher MSE than SRS under control-only spillover. Under both-arms spillover, their MSEs were 17.25 and 17.26 times that of SRS for the unguided and guided versions, respectively. With the non-oracle model, MSE fell from 90/10 to 100/0 for quadrants under control-only spillover with both neighbor types, and the both-arms penalty for halves remained, with coverage as low as 0.607 under rook contiguity. Adding the 100/0 variants enlarged the Holm families, and 3 of the 756 earlier comparisons rose above 0.05, all in rook sensitivity analyses. No queen oracle conclusion changed. These results cover equal-size partitions of the grid with linearly spaced targets. The North Carolina map was not part of this study.

```{=latex}
\begin{figure}[!htb]
\centering
\includegraphics[width=\linewidth]{figures/fig_igs_re_vs_srs.pdf}
\caption[MSE of the variants of incidence guidance relative to SRS under queen contiguity]{MSE of the variants of incidence guidance relative to SRS on the grid under queen contiguity with the oracle model, pooled over the five incidence configurations and five values of $\beta$. The horizontal axis gives the highest and lowest target saturations in percent, and colors give the number of equal-size regions. The left panels give the targets in order of regional mean incidence, and the right panels give the same targets in random order. Rows are the spillover regimes and have different vertical scales. Bars are paired 95\% $t_{49}$ intervals on the ratio scale, and the line at 1 marks SRS. Circles mark Incidence-Guided Saturation (left) and Saturation Quadrants (right).}
\label{fig:igs-sensitivity}
\end{figure}

\FloatBarrier
```

```{=latex}
% Restore the table-of-contents depth for the chapters that follow.
````

## Batch 3 · Tables, figures and the verifier

### C3X-1. Parameter-grid caption (line 673)

**Recommendation:** Accept. The caption calls the grid "the original eight-design simulation" and omits the SRS scenarios. From S3, adapted to the chapter ("compared in the body of this chapter").

**Current**

```qmd
\caption[Parameter grid of the original eight-design simulation]{Parameter grid of the original eight-design simulation. The crossed parameters form $5 \times 5 \times 2 \times 4 \times 4 \times 2 \times 8 = 12{,}800$ scenarios, 9,600 of them for the six designs compared in the body of this chapter. The settings in the lower panel are the same in every scenario.}
```

**Proposed**

```qmd
\caption[Parameter grid of the grid simulation]{Parameter grid of the grid simulation. The crossed parameters form $5 \times 5 \times 2 \times 4 \times 4 \times 2 \times 8 = 12{,}800$ scenarios for the eight designs, 9,600 of them for the six designs compared in the body of this chapter. SRS, run afterward on the same grid, adds 1,600 scenarios (14,400 in all). The settings in the lower panel are the same in every scenario.}
```

### C3X-2. SRS in the adjacent-pair table (line 965 and `tab:adjacent-pairs`)

**Recommendation:** Accept. The rule "SRS in every design test" (DECISIONS 2026-10-25) and 2d D2 put SRS in the chain. Under queen contiguity SRS falls between Incidence-Guided Saturation (0.091) and Balanced Quartiles (0.131). The Incidence-Guided Saturation against Balanced Quartiles row ($-0.040$, $(-0.043, -0.037)$, 50) leaves the table because the two are no longer adjacent, as in T2. New numbers come from `results/srs_benchmark/srs_benchmark_summary.txt` (written by `code/15_srs_benchmark.R`), queen oracle section 3, "Paired surface-level tests vs SRS (tau = 1)": Incidence-Guided Saturation minus SRS $-0.03781$ $(-0.04046, -0.03516)$, frac_better 1.00. Balanced Quartiles minus SRS $+0.00233$ $(-0.00073, 0.00540)$, Wilcoxon $p = 0.153$, frac_better 0.42 (21 of 50 for Balanced Quartiles and 29 for SRS). They match the T2 prose already in the SI. Apply together with C3F-1.

**Current (edit 1)**

```qmd
We also compared designs adjacent in the six-design mean-MSE order over the 50 independent configuration-by-surface units. Table \ref{tab:adjacent-pairs} gives the mean paired differences in MSE with $t_{49}$ 95\% intervals, and the Wilcoxon $p$-value was at most $4.9\times10^{-9}$ in every case. The better design of each pair had the lower MSE in all 50 units,
```

**Proposed (edit 1)**

```qmd
We also compared designs adjacent in the mean-MSE order, with SRS in its place, over the 50 independent configuration-by-surface units. Table \ref{tab:adjacent-pairs} gives the mean paired differences in MSE with $t_{49}$ 95\% intervals. Every adjacent pair differed significantly except SRS and Balanced Quartiles (Wilcoxon $p = 0.15$, SRS lower in 29 of the 50 units), and the Wilcoxon $p$-value was at most $4.9\times10^{-9}$ in every other case. The better design of each of those pairs had the lower MSE in all 50 units,
```

**Current (edit 2)**

```qmd
{Paired differences in MSE between designs adjacent in the six-design mean-MSE order under queen contiguity
```

**Proposed (edit 2)**

```qmd
{Paired differences in MSE between designs adjacent in the mean-MSE order of the six designs and SRS under queen contiguity
```

**Current (edit 3)**

```qmd
Incidence-Guided Saturation against Balanced Quartiles & $-0.040$ & $(-0.043, -0.037)$ & 50 \\
```

**Proposed (edit 3)**

```qmd
Incidence-Guided Saturation against Simple Random Sampling & $-0.038$ & $(-0.040, -0.035)$ & 50 \\
Simple Random Sampling against Balanced Quartiles & $-0.002$ & $(-0.005, 0.001)$ & 29 \\
```

### C3F-1. Verifier: remove the `PENDING_SI` entry for `tab:adjacent-pairs`

**Recommendation:** Accept option A. `paper/tools/verify_manuscripts.py` (lines 26–28) holds `tab:adjacent-pairs` as pending, and its supplement-derivation loop (from line 284) requires every chapter table to appear in the SI with an identical body. Removing the entry alone therefore makes the verifier fail, because the SI reports these tests in prose.

- **Option A (recommended).** Add `tab:adjacent-pairs` to the unabridged SI only, inside `\iffullsi ... \fi` next to the T2 paragraph, with the body of the C3X-2 table, and delete lines 26–28 of the verifier. The verifier already reads the unabridged branch, the curated SI keeps its prose, and both documents then report the same SRS rows.
- **Option B.** Keep the SI as it is, delete the `PENDING_SI` entry and add a named chapter-only exception for this label. This leaves the chapter table without an SI counterpart.

Either option is an SI or tool edit outside the chapter and needs your choice. Run `python3 paper/tools/verify_manuscripts.py` after C3X-2.

### C3F-2. Rebuild the chapter copy of `fig_incidence_rejection_queen.pdf` to the thesis rule

**Recommendation:** Accept. The untracked chapter copy is the SI file: 8.5 × 10 in, text 6–10 pt at that size (about 4.6–7.6 pt when printed at 6.5 in), and Helvetica and Symbol are not embedded (`pdffonts`). It fails all three parts of the thesis rule. Producer: `code/31_incidence_coefficient_summary.R`. Proposed: draw the chapter version in `code/28_thesis_chapter_figures.R` from the saved summary (`results/incidence_coef_2026-11/full/summary/`), at 6.5 in wide with `thesis_text_sizes()` and `save_thesis_pdf()`, as code/28 does for the other chapter figures. Nothing is simulated. The SI copy stays as it is.

### C3F-3. Track the two new chapter figures

**Recommendation:** Accept. `paper/thesis_chapter/figures/fig_igs_re_vs_srs.pdf` and `fig_incidence_rejection_queen.pdf` are untracked. Commit them with C3N-8 and C3N-9. If C3F-2 builds the rejection figure in code/28, add `fig_incidence_rejection_queen.pdf` to the `THESIS_BUILT` list in `paper/tools/copy_figures.sh` (line 55), or the next run of that script overwrites the chapter copy with the SI version. `fig_igs_re_vs_srs.pdf` already meets the rule (6.5 × 5.5 in, 8–9 pt, ArialMT embedded), and the shared copy can stay.

### C3F-4. `fig_tau_sensitivity_6design_combined.pdf`: band note and hyphen (F5, F8)

**Recommendation:** Accept. The chapter copy still prints the y-axis title "Mean MSE (band: ± average scenario−level Monte Carlo SE)", with a minus glyph for the hyphen. The caption already explains the band (and C3T-27 rewrites that sentence). Producer: `code/28_thesis_chapter_figures.R`, line 248. Proposed: axis title "MSE", as F5 set for the SI figure. Rerun code/28.

### C3F-5. `fig_power_vs_n.pdf`: "college-level" prints with minus signs (F8)

**Recommendation:** Accept. The panel labels read "Effect of 1 college−level SD". Producer: `code/28_thesis_chapter_figures.R`, lines 326–327, which build the label without the `hy()` soft-hyphen helper used elsewhere in the script. Proposed: wrap both `sprintf()` calls in `hy()`. Rerun code/28.

### C3F-6. `real_sud/refined_allocation_risk.pdf`: legend labels (F3)

**Recommendation:** Accept. The chapter copy (Figure `fig:nc-risk`) prints "Estimated q90" and "Estimated worst decile". Producer: `code/28_thesis_chapter_figures.R`, lines 490–493. Proposed: "Estimated 90th percentile" and "Estimated worst-decile mean", as F3 set for the SI. Rerun code/28. The same script's `real_sud/primary_setting_ratios.pdf` prints "rho=0.5, psi=0.8" (F1, line 469), but the chapter does not include that file (see "Checked, no change").

### C3F-7. `real_sud/design_examples.png`: numbered panel titles (F4)

**Recommendation:** Accept. The chapter copy (Figure `fig:nc-examples`) still has "1. Checkerboard (27/58)" through "9. Simple Random Sampling (29/58)". The CTJ folder copy was fixed in 2c. Producer: `application/code/manuscript_maps.R`, line 56, which writes both folders from `rs_maps()`. The F4 fix is already in `render_real_sud_companion.R`. Rerunning it needs the restricted `setup.rds`. Separate pre-existing issue for your judgment. The map is saved 14 in wide and printed at 6.5 in, which puts its titles well under 8 pt. A thesis-width export (about 6.5 × 4 in with 8 pt text) would fix both.

## Batch 4 · Methods text (grid study)

### C3T-14. Define queen and rook contiguity at first use (line 61)

**Recommendation:** Accept. The aims paragraph uses both terms before Methods defines them. From HC-2.

**Current**

```qmd
with queen contiguity as the primary neighbor definition and rook contiguity as a sensitivity analysis.
```

**Proposed**

```qmd
with neighbors that share an edge or a corner (queen contiguity) as the primary definition and neighbors that share an edge (rook contiguity) as a sensitivity analysis.
```

### C3T-15. Row-standardization (line 71)

**Recommendation:** Accept. From NY-4.

**Current**

```qmd
Row-standardization expresses each cluster's spillover exposure as an average over its own neighborhood rather than a raw sum, which would otherwise depend on how many neighbors the cluster has.
```

**Proposed**

```qmd
Row-standardization makes each cluster's spillover exposure an average over its own neighborhood, which keeps exposure comparable between clusters with different numbers of neighbors.
```

### C3T-16. Direct effect before propagation (line 78)

**Recommendation:** Accept. From NYC-2, the article's wording of the same sentence.

**Current**

```qmd
The effect $\beta$ is the coefficient of a cluster's own treatment before spatial propagation, not an average total effect after applying $(\mathbf{I}-\rho \mathbf{W})^{-1}$.
```

**Proposed**

```qmd
The effect $\beta$ is the coefficient of a cluster's own treatment before outcomes spread between neighbors through $(\mathbf{I}-\rho \mathbf{W})^{-1}$.
```

### C3T-17. Ties in the Poisson surfaces (line 91)

**Recommendation:** Accept. Removes a prose colon. From HW-1.

**Current**

```qmd
Ties remain at this population size: the Poisson surfaces used here have between 60 and 76 distinct values among the 100 clusters.
```

**Proposed**

```qmd
Ties remain at this population size, and each Poisson surface has only 60 to 76 distinct values among its 100 clusters.
```

### C3T-18. Shared draws and the matched-surface assumption (line 97)

**Recommendation:** Accept. From SN-17 block 1 and HW-1. The parenthetical name becomes a closing clause.

**Current (edit 1)**

```qmd
The design draws and the outcome noise are shared more selectively, as described in @sec-simulation.
```

**Proposed (edit 1)**

```qmd
@sec-simulation describes how the design draws and the outcome noise are shared.
```

**Current (edit 2)**

```qmd
We call this the matched-surface assumption: the incidence observed when a trial is planned is the incidence that drives its outcome.
```

**Proposed (edit 2)**

```qmd
This assumes that the incidence observed when a trial is planned is the incidence that drives its outcome, which we call the matched-surface assumption.
```

### C3T-19. Deterministic and tie-broken designs (line 103)

**Recommendation:** Accept. From HW-2. Edit 2 drops "all 25 draws on a surface", and the 25 draws per surface stay in edit 1.

**Current (edit 1)**

```qmd
Checkerboard is the only deterministic design: it gives the same assignment in every draw, although each of its fits still receives its own outcome noise. Every other design was drawn independently 25 times for each incidence surface.
```

**Proposed (edit 1)**

```qmd
Checkerboard is the only deterministic design. Its assignment is fixed, and every fit gets its own outcome noise. We drew every other design independently 25 times per incidence surface.
```

**Current (edit 2)**

```qmd
High Incidence Focus is random only through that tie-breaking. Under the iid and spatial mechanisms, incidence values are distinct, which means High Incidence Focus gives the same assignment in all 25 draws on a surface, and that assignment changes only from one surface to the next. Isolation Buffer builds its assignment from the neighbor type of the scenario, and it treats 38 clusters on average under rook contiguity and 21 under queen contiguity.
```

**Proposed (edit 2)**

```qmd
High Incidence Focus varies only through tie-breaking. Under the iid and spatial mechanisms, incidence values have no ties, and its assignment is fixed within each surface and varies across surfaces. Isolation Buffer uses the neighbor type of the scenario and treats an average of 38 clusters with rook contiguity and 21 with queen contiguity.
```

### C3T-20. Oracle model, the `lagsarlm` refit and the "Spill" label (lines 146, 148, 1191)

**Recommendation:** Accept. From HW-3 and SN-14. The SI's ", so its results" stays as the chapter's "Therefore," sentence. The 1% `lagsarlm` refit is a fact the chapter lacks. Edit 3 retires the "Spill" label in the A10 caption, which edit 1 and C3T-10 remove elsewhere.

**Current (edit 1)**

```qmd
We estimated $\beta$ with two maximum likelihood spatial lag models, both fit to every simulated outcome and both including an intercept. The primary, oracle model has the spatial lag $\rho \mathbf{W}\mathbf{Y}$ and the covariates $\mathbf{Z}$, $\text{Spill} = \mathbf{S}(\mathbf{Z})$, and $\mathbf{X}_k$, which matches the outcome model in Equation \ref{eq:outcome}, including the true regime-specific spillover covariate from Equation \ref{eq:spillover}. The oracle model assumes that the analyst knows the form of the spillover mechanism and has the exact exposure term. Therefore, its results are a best case for each design. We return to this assumption as a limitation. The non-oracle model, a sensitivity analysis, has the covariates $\mathbf{Z}$ and $\mathbf{X}_k$ only and omits the spillover term. Both models are fit to the same outcomes.
```

**Proposed (edit 1)**

```qmd
We fit two maximum likelihood spatial lag models with intercepts to each simulated outcome. The primary, oracle model includes the spatial lag $\rho \mathbf{W}\mathbf{Y}$ and the covariates $\mathbf{Z}$, $\mathbf{S}(\mathbf{Z})$ and $\mathbf{X}_k$. It matches the outcome model in Equation \ref{eq:outcome}, using the true regime-specific spillover term in Equation \ref{eq:spillover}. The oracle model assumes that the analyst knows the form of spillover and has the exact exposure term. Therefore, its results are a best case for each design. We return to this assumption as a limitation. The non-oracle model, a sensitivity analysis, includes $\mathbf{Z}$ and $\mathbf{X}_k$ and omits the spillover term.
```

**Current (edit 2)**

```qmd
at roughly 1\% of the computing time.
```

**Proposed (edit 2)**

```qmd
at roughly 1\% of the computing time. We also refit a 1\% sample of the full simulation with `lagsarlm` itself as a check.
```

**Current (edit 3)**

```qmd
The oracle model has the covariates $\mathbf{Z}$, $\text{Spill}$, and $\mathbf{X}_k$, and the non-oracle model omits $\text{Spill}$.
```

**Proposed (edit 3)**

```qmd
The oracle model has the covariates $\mathbf{Z}$, $\mathbf{S}(\mathbf{Z})$, and $\mathbf{X}_k$, and the non-oracle model omits the spillover term $\mathbf{S}(\mathbf{Z})$.
```

### C3T-21. Grid simulation size (line 164)

**Recommendation:** Accept. "Each estimator was fit 3.2 million times" counts the eight designs only, and the SRS run is described in plain words. From HW-4. The chapter keeps "estimator" elsewhere, and this paragraph follows the SI's "model".

**Current**

```qmd
Within each scenario, we drew $J = 25$ assignments from each of the $K = 10$ incidence surfaces, for 250 fits per scenario. Every fit had its own noise vector, and fit $(k, j)$ has outcome $\mathbf{Y}_{kj} = (\mathbf{I} - \rho \mathbf{W})^{-1}[\beta \mathbf{Z}_{kj} + \psi\,\mathbf{S}(\mathbf{Z}_{kj}) + \tau \mathbf{X}_k + \boldsymbol{\varepsilon}_{kj}]$. Checkerboard's 250 fits share one assignment but have 250 distinct noise vectors. The ten surfaces are independent, which makes the surface the primary sampling unit. Each estimator was fit 3.2 million times. SRS assigns exactly 50 of the 100 clusters to treatment by complete randomization, ignoring incidence and geography. We ran it afterward with the same matched surfaces, noise keys and 250-fit structure, which brought the complete benchmark to 14,400 scenarios and 3.6 million fits per estimator. Adding SRS left the original eight-design results unchanged.
```

**Proposed**

```qmd
We drew $J = 25$ assignments for each of the $K = 10$ incidence surfaces, for 250 fits per scenario, with a separate noise vector for each fit. Fit $(k, j)$ uses outcome $\mathbf{Y}_{kj} = (\mathbf{I} - \rho \mathbf{W})^{-1}[\beta \mathbf{Z}_{kj} + \psi\,\mathbf{S}(\mathbf{Z}_{kj}) + \tau \mathbf{X}_k + \boldsymbol{\varepsilon}_{kj}]$. The assignment is fixed for all 250 Checkerboard fits. Within each surface the 25 pairs of assignments and noise vectors are independent. The ten surfaces are also independent, making the surface the primary sampling unit. Each model was fit 3.2 million times for the eight designs. SRS randomly treats 50 of the 100 clusters with no restrictions on incidence or geography. We ran it afterward using the same surfaces, noise and 250 fits per scenario. This brought the benchmark to 14,400 scenarios and 3.6 million fits per model and left the eight-design results unchanged.
```

### C3T-22. Monte Carlo standard errors (line 168)

**Recommendation:** Accept. Removes a prose colon and adds what the pooled standard error averages over. From SR-13 block 1.

**Current**

```qmd
Monte Carlo standard errors come from the surface-level metrics: with $u_k$ a metric computed on the 25 fits from surface $k$, $\text{SE}(u) = \text{sd}(u_1, \dots, u_{10})/\sqrt{10}$, and, for metrics pooled over the five incidence configurations, from the 50 independent configuration-by-surface means.
```

**Proposed**

```qmd
We computed Monte Carlo standard errors from the surface-level metrics. For a metric $u$, we computed $u_k$ on the 25 fits from surface $k$ and set $\text{SE}(u) = \text{sd}(u_1, \dots, u_{10})/\sqrt{10}$. For pooled results we apply the same formula to the 50 independent configuration-by-surface means. When pooling both regimes, each mean averages that surface's 32 blocks. The pooled standard deviation includes variation across configurations.
```

### C3T-23. Source of the control-only increase (line 279)

**Recommendation:** Accept. From NY-18.

**Current**

```qmd
The increase under the control-only regime was variance, not bias.
```

**Proposed**

```qmd
The control-only increase came almost entirely from variance.
```

## Batch 5 · Grid results and case-study methods

### C3T-24. Rook comparison with SRS (line 287)

**Recommendation:** Accept. Says which difference was tested and that the p-value is Holm-adjusted. From SN-11.

**Current**

```qmd
Under rook contiguity the saturation designs had 23\% to 24\% lower control-only MSE than SRS and 3\% to 5\% higher both-arms MSE, significant for Incidence-Guided Saturation ($p = 0.0075$) but not for Saturation Quadrants ($p = 0.15$).
```

**Proposed**

```qmd
Under rook contiguity, the saturation designs had 23\% to 24\% lower control-only MSE than SRS and 3\% to 5\% higher both-arms MSE. The both-arms difference was significant after Holm's adjustment for Incidence-Guided Saturation ($p = 0.0075$) but not for Saturation Quadrants ($p = 0.15$).
```

### C3T-25. Worst-case block MSE (line 313)

**Recommendation:** Accept. `tab:robustness` has no SRS row, and "of any design" claims more than the table shows. From SR-11 blocks 3 and 4 and HW-17. The pointer names the grid allocation-risk subsection, which has no label.

**Current (edit 1)**

```qmd
Incidence-Guided Saturation had the lowest median (0.080) and the lowest worst-case MSE of any design (0.169).
```

**Proposed (edit 1)**

```qmd
Incidence-Guided Saturation had the lowest median (0.080) and worst-case MSE among the six detailed designs (0.169).
```

**Current (edit 2)**

```qmd
Under rook contiguity, Incidence-Guided Saturation also had the lowest worst-case MSE (0.138). These worst cases summarize variation across parameter settings, not the conditional risk across assignments within one setting.
```

**Proposed (edit 2)**

```qmd
Under rook contiguity, Incidence-Guided Saturation also had the lowest worst-case MSE among the six detailed designs (0.138). Worst-case block MSE describes variation across parameter settings. The allocation-specific risk results below examine variation across assignments within a setting.
```

### C3T-26. Shared draws across effect sizes (line 348)

**Recommendation:** Accept. From NY-19.

**Current**

```qmd
Therefore, these comparisons are not independent confirmations of one another.
```

**Proposed**

```qmd
Therefore, the stable order is one finding on shared draws.
```

### C3T-27. Bands and rook Checkerboard in `fig:tau-mse` (line 356)

**Recommendation:** Accept. From HW-13 block 2. Pairs with C3F-4, which removes the in-figure band note.

**Current (edit 1)**

```qmd
Shaded bands show $\pm 1$ times the average of the scenario-level Monte Carlo standard errors of MSE. They describe the precision of a single scenario's MSE, not of the plotted average, whose standard error is smaller (Table \ref{tab:performance}).
```

**Proposed (edit 1)**

```qmd
Shaded bands extend one average scenario-level Monte Carlo standard error of MSE above and below each point. The Monte Carlo standard error of each plotted average is smaller (Table \ref{tab:performance}).
```

**Current (edit 2)**

```qmd
Checkerboard under rook contiguity is flagged because $\beta$ is not identifiable (Section~\ref{sec-estimation}). The five values of $\beta$ share their assignments and noise.}
```

**Proposed (edit 2)**

```qmd
For rook Checkerboard, $\beta$ is not identifiable and the coefficient on $\mathbf{Z}$ estimates $\beta-\psi$ (Section~\ref{sec-estimation}). The five values of $\beta$ share their assignments and noise.}
```

### C3T-67. Treatment-effect sweep: the $\beta = 2$ tie and the six-design scope (line 362)

**Recommendation:** Accept. Added by the independent review. HW-14 (2g) replaced "the two were level (0.160 ... against 0.160)" with "both rounded to 0.160" and scoped "lowest in pooled MSE at every value of $\beta$" to the six designs. The chapter still has the old wording.

**Current**

```qmd
The two were level at $\beta = 2$ (0.160 each), and at $\beta = 3$ Isolation Buffer had the lower control-only MSE. Incidence-Guided Saturation remained lowest in pooled MSE at every value of $\beta$,
```

**Proposed**

```qmd
At $\beta = 2$ both rounded to 0.160, and at $\beta = 3$ Isolation Buffer had the lower control-only MSE. Among the six designs, Incidence-Guided Saturation remained lowest in pooled MSE at every value of $\beta$,
```

### C3T-28. Noise share of the grid risk estimates (line 370)

**Recommendation:** Accept. Removes the colon join. From HW-19. The chapter's split-half and upward-bias sentences that follow stay.

**Current**

```qmd
The precision check preserved the control-only direction, but outcome noise contributed substantially to the estimated variation across allocations, especially under both arms: it accounted for a median 63\% of the variance of Incidence-Guided Saturation's estimated risk under control-only spillover and 92\% under both arms.
```

**Proposed**

```qmd
The precision check preserved the control-only direction. Outcome Monte Carlo noise accounted for a median 63\% of the estimated variance in Incidence-Guided Saturation's assignment risks under control-only spillover and 92\% under both arms.
```

### C3T-68. Grid allocation-risk numbers and their scope (line 370)

**Recommendation:** Accept. Added by the independent review. C3T-28 carries only the second half of HW-19. Block 1 of HW-19 rewrote the first half: the all-200-units claim applies to the 90th percentile and the worst-decile mean, "allocation-specific" became "assignment-specific", and "both arms" became "both-arms spillover". The chapter keeps its own Saturation Quadrants sentence (more specific than the SI's "performed favorably") and its table pointer. The chapter's next sentence ("These results come from the extreme parameter settings only and give qualified support to saturation on downside risk.") stays, although HW-19 cut it from the SI. See the review notes.

**Current**

```qmd
Averaged over the 200 units per regime (fixed surfaces crossed with the extreme combinations of $\rho$ and $\psi$), Incidence-Guided Saturation's control-only mean MSE was 0.140 against 0.217 for SRS. Its estimated 90th percentile of allocation-specific MSE was 0.173 against 0.296, and its estimated worst-decile mean was 0.191 against 0.346, lower than SRS in all 200 units on both measures. Under both arms, the corresponding mean MSE was 0.046 against 0.043 and the worst-decile mean 0.058 against 0.054.
```

**Proposed**

```qmd
Under control-only spillover, Incidence-Guided Saturation's mean MSE was 0.140 against 0.217 for SRS, its estimated 90th percentile of assignment-specific MSE was 0.173 against 0.296, and its estimated worst-decile mean was 0.191 against 0.346. The last two were below SRS in all 200 units, which cross fixed surfaces with the extreme values of $\rho$ and $\psi$. Under both-arms spillover, the corresponding mean MSE was 0.046 against 0.043 and the worst-decile mean was 0.058 against 0.054.
```

### C3T-29. Observed inputs and planning assumptions (line 378)

**Recommendation:** Accept. Tells the reader which inputs are observed and which are assumed. From SN-6. The SI's pointer to its participation section becomes Appendix A5.

**Current**

```qmd
The case study uses no observed intervention outcomes and no estimated relationship between local SUD incidence and education performance. Any eventual effect on SUD incidence requires a separate clinical evaluation.
```

**Proposed**

```qmd
The case study uses no observed intervention outcomes and no estimated relationship between local SUD incidence and the education outcome. Any eventual effect on SUD incidence requires a separate clinical evaluation. The observed inputs are the service-area map, the county-to-college assignment, and annual deaths and populations. The outcome model, the effect size, the spatial dependence $\rho$, the spillover regime and its magnitude $\psi$ are planning assumptions, as are the volunteering mechanism, the number of professionals per college and the intraclass correlation in the participation analyses (Appendix \ref{app:nc-ext}).
```

### C3T-30. County-to-college assignments (lines 380, 790)

**Recommendation:** Accept. From SN-25 block 4, as revised by the author.

**Current (edit 1)**

```qmd
The county-to-college assignments are analysis conventions, not exclusive service commitments (Appendix \ref{app:nc}).
```

**Proposed (edit 1)**

```qmd
We made these county-to-college assignments for the analysis only, and they do not limit which counties a college serves (Appendix \ref{app:nc}).
```

**Current (edit 2)**

```qmd
These county-to-college assignments are analysis conventions, not exclusive service commitments.
```

**Proposed (edit 2)**

```qmd
We made these assignments for the analysis only, and they do not limit which counties a college serves.
```

### C3T-31. Adjacency as a proxy (line 410)

**Recommendation:** Accept. From NY-8.

**Current**

```qmd
Adjacency here is a prespecified geographic proxy for exposure, not a measured network of professional communication.
```

**Proposed**

```qmd
Adjacency here is a prespecified geographic proxy for exposure between neighboring service areas.
```

### C3T-32. Budget table caption (line 423)

**Recommendation:** Accept. From NY-2.

**Current**

```qmd
Shares are geographic reach, not healthcare-professional enrollment.}
```

**Proposed**

```qmd
Shares measure geographic reach.}
```

### C3T-71. Colon before the case-study outcome equation (line 448)

**Recommendation:** Accept. Added by the independent review. HW-7 (2g) ended this sentence with a period and "and" in place of the colon. The rest of HW-7 is already in the chapter in its own form (interval formula, two-sided test, one-sided grid power).

**Current**

```qmd
and fix $\beta=1$, $\sigma=1$:
```

**Proposed**

```qmd
and fix $\beta=1$ and $\sigma=1$.
```

### C3T-33. Target of the case-study analysis (line 452)

**Recommendation:** Accept. From NY-12, with the approved SR-10 sentences on what each regime means for a college.

**Current**

```qmd
As in the grid model, $\beta$ is the coefficient before spatial propagation, not an average total effect after applying $(\mathbf{I}-\rho \mathbf{W})^{-1}$.
```

**Proposed**

```qmd
As in the grid model, the target is the direct coefficient $\beta$, the effect of a college's own treatment before outcomes spread through $(\mathbf{I}-\rho \mathbf{W})^{-1}$. Under control-only spillover, a control college's outcome rises with the share of its neighbors that are treated. Under both-arms spillover, treated colleges also gain from treated neighbors. When every neighbor is treated, the spillover term adds $\psi = 0.5$ or $0.8$, half or four-fifths of the direct effect.
```

## Batch 6 · Case-study methods and results

### C3T-34. Case-study settings (line 454)

**Recommendation:** Accept. From HW-8 and NY-13. The sentence on the hypothetical association keeps the chapter's "changes the outcome and the adjustment together".

**Current (edit 1)**

```qmd
four years and nine designs (432 reported settings).
```

**Proposed (edit 1)**

```qmd
four years, and the eight designs and SRS (432 reported settings).
```

**Current (edit 2)**

```qmd
It changes the outcome and the adjustment together and is not an empirical association between incidence and education.
```

**Proposed (edit 2)**

```qmd
It changes the outcome and the adjustment together, and the association it builds in is hypothetical.
```

**Current (edit 3)**

```qmd
The full case study reports 1,248 settings, computed from 936 distinct simulation settings, because reported settings with the same allocations and outcome model share one simulation.
```

**Proposed (edit 3)**

```qmd
Settings with the same possible assignments and outcome model share a simulation, giving 936 distinct simulations for 1,248 reported settings.
```

### C3T-35. Case-study precision and shared simulations (line 456)

**Recommendation:** Accept. From HW-8 and SN-20. Edit 2 takes the seed sentence that C3T-51 removes from A4. Edit 3 takes the recalculation sentence that C3T-54 removes from A5.

**Current (edit 1)**

```qmd
A design with only one possible assignment instead received $R=1,000$ outcomes for that assignment, rather than 100 repeated copies of it counted as separate allocations. A prespecified adaptive schedule extended allocations or outcomes where needed
```

**Proposed (edit 1)**

```qmd
A design with one possible assignment received $R=1,000$ outcomes for that assignment. A prespecified adaptive schedule extended the number of assignments or outcomes where needed
```

**Current (edit 2)**

```qmd
required $R=2,000$, and all randomized designs met these targets at $J=100$ and $R=100$ (Appendix \ref{app:nc-ext}).
```

**Proposed (edit 2)**

```qmd
required $R=2,000$. All randomized designs met both targets at $J=100$ and $R=100$ (Appendix \ref{app:nc-ext}). Each random draw used a fixed seed.
```

**Current (edit 3)**

```qmd
reuse the same simulations across years because annual incidence enters neither their allocation nor the outcome model. As a result, their repeated values are not independent annual confirmations.
```

**Proposed (edit 3)**

```qmd
share simulations across years because annual incidence affects neither their assignments nor their outcomes. We recalculate treatment counts and population shares for each year. Therefore, their results in different years are not independent confirmations.
```

### C3T-36. Case-study rank tests (line 460)

**Recommendation:** Accept. From NY-16 and HW-30.

**Current (edit 1)**

```qmd
We compared the nine designs with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid,
```

**Proposed (edit 1)**

```qmd
We compared the eight designs and SRS with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid,
```

**Current (edit 2)**

```qmd
The four years share one map and the designs whose allocations ignore annual incidence share their simulations across years, which means these $p$-values summarize consistency rather than independent evidence.
```

**Proposed (edit 2)**

```qmd
The four years share one map, and some designs and SRS reuse their simulations across years. Therefore, these $p$-values describe how consistently the designs ranked across blocks that share data.
```

### C3T-70. Units of the case-study MSE (line 480)

**Recommendation:** Accept. Added by the independent review. SN-8 (2g) states the units where the first case-study MSE values appear. The chapter has no such sentence.

**Current**

```qmd
Figure \ref{fig:nc-mse} compares the eight designs with SRS by year and regime. Under control-only spillover, Saturation Quadrants had mean MSE 0.223 against 0.308 for SRS, a 27.4\% reduction.
```

**Proposed**

```qmd
Figure \ref{fig:nc-mse} compares the eight designs with SRS by year and regime. With $\sigma = 1$, MSE is in squared college-level residual standard deviations. Under control-only spillover, Saturation Quadrants had mean MSE 0.223 against 0.308 for SRS, a 27.4\% reduction.
```

### C3T-37. High Incidence Focus by planning year (line 480)

**Recommendation:** Accept. "We cannot recommend it for every year" sits awkwardly with the no-single-design rule. From R4 and HW-22 block 1.

**Current**

```qmd
The performance of High Incidence Focus under control-only spillover depended on the year. As a result, we cannot recommend it for every year.
```

**Proposed**

```qmd
Under control-only spillover, the accuracy of High Incidence Focus depended on the planning year.
```

### C3T-38. High Incidence Focus and randomization (line 482)

**Recommendation:** Accept. From NY-21 and SN-13.

**Current**

```qmd
SRS remained competitive with balancing, blocking and saturation. High Incidence Focus had low both-arms MSE in the primary model, but its allocation is deterministic for these fixed rates. Its low MSE reflects model-based accuracy conditional on that allocation and does not show that the design has the protection that randomization provides.
```

**Proposed**

```qmd
SRS remained competitive with balancing, blocking and saturation, and small point differences do not establish equivalence. High Incidence Focus had low both-arms MSE in the primary model, but its allocation is deterministic for these fixed rates. Its low MSE describes model-based accuracy for that one allocation. With only one possible allocation, its validity rests on the outcome model alone.
```

### C3T-39. "Eight designs and SRS" (lines 484, 1513)

**Recommendation:** Accept. From HW-30 (block 1 for line 1513).

**Current (edit 1)**

```qmd
was between 0.941 and 0.949 for all nine designs and both regimes, and between 0.941 and 0.951 in single years. The Wald interval covered only 0.926 to 0.936 for the same nine designs.
```

**Proposed (edit 1)**

```qmd
was between 0.941 and 0.949 for the eight designs and SRS in both regimes, and between 0.941 and 0.951 in single years. The Wald interval covered only 0.926 to 0.936 for the same eight designs and SRS.
```

**Current (edit 2)**

```qmd
the mean power over the nine designs fell from 0.59 to 0.46 with 20 colleges
```

**Proposed (edit 2)**

```qmd
the mean power over the eight designs and SRS fell from 0.59 to 0.46 with 20 colleges
```

### C3T-40. Population denominators (line 486)

**Recommendation:** Accept. From NY-3.

**Current**

```qmd
These population denominators measure geographic reach, not the number of healthcare professionals educated, intervention cost or participant-level precision.
```

**Proposed**

```qmd
These population denominators measure geographic reach only.
```

### C3T-41. Regional incidence on the North Carolina map (line 488)

**Recommendation:** Accept. From HW-23 block 1, which also replaces the bare "NC" (SN-9).

**Current**

```qmd
On the NC map incidence differs little between regions, because
```

**Proposed**

```qmd
In North Carolina it differs little, because
```

### C3T-42. Tail summaries and the reused assignments (line 501)

**Recommendation:** Accept. Says what the tail summaries describe and that the refinement reused the main run's assignments. From HW-24 block 1 and NY-22.

**Current**

```qmd
These reductions suggest a benefit for downside risk as well as for average error, under the studied model and sampled assignments.
```

**Proposed**

```qmd
Under the study model, the saturation designs had a lower worst-decile mean of conditional MSE across the sampled assignments, as well as lower average MSE. These summaries describe the spread of conditional MSE across sampled assignments. The refinement reused the main run's assignments.
```

### C3T-43. Downside table caption (line 512)

**Recommendation:** Accept. From HW-20 block 4. The chapter's list of map adaptations stays as its own sentence, and the short caption is unchanged.

**Current**

```qmd
{Mean conditional MSE and worst-decile mean conditional MSE of each design divided by those of SRS in the same setting and spillover regime (below 1 is lower error than SRS). Worst-decile means are estimated from all outcomes of each assignment, without split-half selection. Grid: allocation-risk study, queen contiguity, $\beta=1$, averaged over 200 units, each a fixed incidence surface with one extreme combination of $\rho$ and $\psi$. NC: refined queen settings $\rho\in\{0,0.5\}$, $\psi=0.8$, averaged equally over 2018--2021 and both values of $\rho$. NC designs are the map adaptations (Checkerboard as a greedy graph-coloring assignment, $2\times2$ Blocking on spatial blocks of nearby clusters, and both saturation designs on four fixed regions). A design with one allocation has worst-decile mean equal to its mean.}
```

**Proposed**

```qmd
{Mean and worst-decile mean conditional MSE divided by the corresponding SRS value in the same setting and spillover regime. Ratios below 1 favor the design over SRS. The queen grid study fixes $\beta=1$ and averages over 200 units, each an incidence surface at one extreme combination of $\rho$ and $\psi$. The North Carolina analysis uses the map adaptations and averages the refined queen settings equally over 2018--2021 and both values of $\rho\in\{0,0.5\}$, at $\psi=0.8$. In the map adaptations, Checkerboard is a greedy graph-coloring assignment, $2\times2$ Blocking uses spatial blocks of nearby clusters, and both saturation designs use four fixed regions. For a design with one assignment, the worst-decile mean equals its mean. Worst deciles are selected and evaluated with the same simulated outcomes, without split-half selection.}
```

## Batch 7 · Case-study results, Discussion and appendix methods

### C3T-44. Uncertainty of the tail estimates (line 552)

**Recommendation:** Accept. Adds that the refinement did not enlarge the sample of assignments, the reason each worst decile stays at ten draws, and says what the single-assignment equality means. From HW-25. The chapter's in-line definition of the Jaccard index stays. "All 112" comes from the SI.

**Current (edit 1)**

```qmd
The tail estimates remain uncertain. With $J=100$ the worst decile contains only ten sampled draws, and every reported setting in the tail refinement for designs with more than one possible assignment flags uncertain worst-decile membership.
```

**Proposed (edit 1)**

```qmd
The tail estimates remain uncertain. With $J=100$, each worst decile contains ten sampled draws. The refinement added outcomes to the same assignments, which reduced noise in $\hat m(a)$ but did not enlarge the sample of assignments. Worst-decile membership was flagged as uncertain in all 112 reported refinement settings for designs with more than one possible assignment.
```

**Current (edit 2)**

```qmd
came from outcome Monte Carlo noise, which means small both-arms tail differences do not establish better downside protection. Checkerboard and High Incidence Focus each have only one possible assignment in this setting. They show no spread across assignments, but that does not mean their error is low.
```

**Proposed (edit 2)**

```qmd
came from outcome Monte Carlo noise. The both-arms tail differences are too small to establish better downside protection. Checkerboard and High Incidence Focus each have only one possible assignment in this setting. Therefore, their 90th percentile and worst-decile mean equal their conditional mean MSE. That equality does not imply low error, and their conditional MSE is still estimated from simulated outcomes.
```

### C3T-45. Prespecified mean-rank rule (line 554)

**Recommendation:** Accept. From NY-23.

**Current**

```qmd
We retain the prespecified mean-rank rule rather than select a summary retrospectively.
```

**Proposed**

```qmd
We retain the prespecified mean-rank rule.
```

### C3T-46. Bayesian CAR chains and the spatial lag model (line 560)

**Recommendation:** Accept. "Converged" becomes "met the diagnostics for $\beta$" (SR-17), and the last sentence states what the Bayesian spatial lag model shares with the data and its bias (HW-29). The 0.039 bound is the chapter's own number from line 1515.

**Current (edit 1)**

```qmd
even when the chains converged, with coverage of 0.88 to 0.90,
```

**Proposed (edit 1)**

```qmd
even in chains that met the convergence diagnostics for $\beta$, with coverage of 0.88 to 0.90,
```

**Current (edit 2)**

```qmd
The correctly specified Bayesian spatial lag model does not share the problem.
```

**Proposed (edit 2)**

```qmd
The Bayesian spatial lag model uses the same spatial multiplier as the simulated outcome, and its bias was at most 0.039 (Appendix \ref{app:nc-ext-results}).
```

### C3T-47. Colleges needed: conditions and units (lines 566, 577, 579)

**Recommendation:** Accept. States the $\rho$ and $\psi$ behind the 35 to 37 colleges (SR-7) and puts every college-level effect in residual standard deviations before spatial propagation (SR-6).

**Current (edit 1)**

```qmd
Every other design needed 35 to 37 colleges for 80\% power to detect a college-level effect of one standard deviation, and about 19 colleges for 1.5 standard deviations (Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). Checkerboard and Isolation Buffer did not reach 80\% power at $\beta = 1$ even with all 58 colleges.
```

**Proposed (edit 1)**

```qmd
With $\rho = 0.5$ and $\psi = 0.8\beta$, every other design needed an estimated 35 to 37 colleges for 80\% power at an effect of one college-level residual standard deviation, and about 19 colleges at 1.5 standard deviations (Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). At $\beta = 1$, Checkerboard and Isolation Buffer did not reach 80\% power with all 58 colleges.
```

**Current (edit 2)**

```qmd
was 0.76 to 0.78 college-level standard deviations for every design
```

**Proposed (edit 2)**

```qmd
was 0.76 to 0.78 college-level residual standard deviations for every design
```

**Current (edit 3)**

```qmd
A college-level effect of $\beta$ standard deviations corresponds to an individual-level standardized effect
```

**Proposed (edit 3)**

```qmd
A structural college-level effect of $\beta$ residual standard deviations, before spatial propagation, corresponds to an individual-level standardized effect
```

### C3T-48. High Incidence Focus in the Discussion (line 605)

**Recommendation:** Accept. From NY-21, with "North Carolina" for "NC" (SN-9).

**Current**

```qmd
High Incidence Focus gave one fixed allocation whose MSE in NC was below the SRS worst-decile mean under both regimes. However, its accuracy is model-based rather than randomization-based.
```

**Proposed**

```qmd
High Incidence Focus gave one fixed allocation whose MSE in North Carolina was below the SRS worst-decile mean under both regimes. However, with only one possible allocation, its validity rests on the outcome model alone.
```

### C3T-49. Queen Checkerboard identification (line 763)

**Recommendation:** Accept. The chapter already says "identifiable" everywhere else (C3-2). This is the one sentence where I4 also changed the wording. Also moves "As a result" to the house "Therefore".

**Current**

```qmd
As a result, $\beta$ is identifiable, but imprecisely.
```

**Proposed**

```qmd
Therefore, $\beta$ is identifiable but imprecisely estimated.
```

### C3T-50. What the rerun supplies (line 769)

**Recommendation:** Accept. Every term of the decomposition, as well as $\hat\rho$, comes from the rerun. From HW-15.

**Current**

```qmd
The numbers that involve $\hat\rho$ come from a re-run of every block for the six designs, both neighbor types and all five values of $\beta$ (2.4 million fits of the oracle model), which reproduced the stored results in all 9,600 scenarios.
```

**Proposed**

```qmd
The numbers that involve $\hat\rho$, the standard error of $\hat\beta$ or a term of Equation \ref{eq:error-decomposition} come from a rerun of every block for the six designs, both neighbor types and all five values of $\beta$ (2.4 million fits of the oracle model), which reproduced the stored results in all 9,600 scenarios. The other numbers come from the original simulations.
```

### C3T-51. Case-study reproducibility (line 786)

**Recommendation:** Accept. From HW-9 block 1. The seed sentence moves to the precision paragraph in the body (C3T-35, edit 2).

**Current**

```qmd
The primary case-study numbers come from the aggregate performance and allocation-risk summaries of the main run and of the separate tail refinement. Annual tables average the six queen parameter pairs equally, without pooled Monte Carlo intervals. Each random draw used a fixed, reproducible seed. An independent check recomputed performance from the stored fits and verified the annual treatment budgets, the regional orderings, and the allocations and outcomes reused from the main run in the tail refinement.
```

**Proposed**

```qmd
The annual results give equal weight to the six queen parameter pairs and are reported without pooled Monte Carlo intervals. An independent check recomputed performance from the simulated trial fits, verified annual treatment budgets and regional orderings, and confirmed that the tail refinement reused the main run's assignments and outcomes.
```

### C3T-52. Regional summaries (line 797)

**Recommendation:** Accept. From NY-9 and NY-10.

**Current (edit 1)**

```qmd
Mean raw rates and pooled regional rates are evaluated as prespecified sensitivities rather than chosen by performance.
```

**Proposed (edit 1)**

```qmd
Mean raw rates and pooled regional rates are evaluated as prespecified sensitivity analyses.
```

**Current (edit 2)**

```qmd
Mean ranks align with its rank-scaled Poisson cases but are a deliberate adaptation for the case study, not the exact rule used for continuous grid incidence.
```

**Proposed (edit 2)**

```qmd
Mean ranks agree with it in the rank-scaled Poisson configurations and can order regions differently when incidence is continuous. They are a deliberate adaptation for the case study.
```

### C3T-53. Bayesian chain diagnostics (line 865)

**Recommendation:** Accept. From SR-17 and HW-30.

**Current (edit 1)**

```qmd
where fewer than half the chains converged,
```

**Proposed (edit 1)**

```qmd
where fewer than half the chains met the diagnostics for $\beta$,
```

**Current (edit 2)**

```qmd
both regimes and the nine designs, with 100 fits per setting
```

**Proposed (edit 2)**

```qmd
both regimes, the eight designs and SRS, with 100 fits per setting
```

**Current (edit 3)**

```qmd
A chain was counted as converged when the effective sample size of $\beta$
```

**Proposed (edit 3)**

```qmd
A chain met the diagnostics for $\beta$ when the effective sample size of $\beta$
```

## Batch 8 · Appendix methods and results

### C3T-54. Refinement bookkeeping (line 869)

**Recommendation:** Accept. From SN-20 (the repeated "met the targets" sentence stays in the body only) and NY-14. The last sentence moves to the body (C3T-35, edit 3).

**Current (edit 1)**

```qmd
Whichever source of Monte Carlo error, allocations or outcomes, contributed more to the shortfall from the precision targets was extended. All randomized designs met the targets at $J=100$ and $R=100$.
```

**Proposed (edit 1)**

```qmd
It extended whichever source of Monte Carlo error, allocations or outcomes, contributed more to the shortfall from the precision targets.
```

**Current (edit 2)**

```qmd
Therefore, the variance due to outcome noise involves squared weights $(n_a/J)^2$, not $n_a/J^2$ as if duplicate estimates were independent.
```

**Proposed (edit 2)**

```qmd
Therefore, the variance due to outcome noise involves the squared weights $(n_a/J)^2$, because repeated draws of one allocation share one risk estimate.
```

**Current (edit 3)**

```qmd
A design whose sampled assignments happened to coincide is distinguished from one with only one possible assignment. Approximate Monte Carlo standard errors of setting-specific differences between designs treat the designs' random draws as independent, but averages across settings or models do not assume that shared allocations are independent. Settings with the same set of possible assignments shared their simulations, with annual treatment counts and population shares recalculated.
```

**Proposed (edit 3)**

```qmd
It distinguishes a design whose sampled assignments happened to coincide from one with only one possible assignment. Approximate Monte Carlo standard errors of setting-specific differences between designs treat the designs' random draws as independent, whereas averages across settings or models do not assume that shared allocations are independent.
```

### C3T-55. Selection-bias check and the Jaccard index (line 871)

**Recommendation:** Accept. From HW-9 block 2 and NY-15, with "North Carolina" (SN-9).

**Current**

```qmd
We report the estimated variance of risk both raw and with the contribution of outcome Monte Carlo noise subtracted. Signed corrected estimates are retained without truncation at zero. None is negative in the NC refinement. Quantiles and worst-decile means use noisy estimated conditional risks. To check for selection bias, worst-decile assignments were also identified with one half of the simulated outcomes and evaluated with the other half, and the roles of the two halves were then reversed. The cross-selected reductions relative to SRS were similar to or slightly larger than the reported ones. The Jaccard overlap of the two worst deciles describes sampled tail membership, not a confidence interval or population tail probability. For a design with only one possible assignment, an overlap of one is uninformative.
```

**Proposed**

```qmd
We report the variance of the estimated allocation risks before and after subtracting the contribution of outcome Monte Carlo noise. We did not truncate the corrected estimates at zero, and none was negative in the North Carolina refinement. The estimated risks used for quantiles and worst-decile means still contain outcome noise. We checked for selection bias by selecting the worst decile with one half of the simulated outcomes, evaluating it with the other half, then exchanging the roles of the halves. The resulting reductions relative to SRS were similar to or slightly larger than the reductions reported from all outcomes. The Jaccard index measures agreement between the two sampled worst deciles. It divides the number of assignments in both deciles by the number in either. It describes how stable sampled tail membership is. With only one possible assignment, the index equals one by construction and carries no information.
```

### C3T-56. Stale precision note (line 904)

**Recommendation:** Accept. `tab:performance` also prints three decimals, and the phrase says nothing. From S2.

**Current**

```qmd
The MSEs and coverages of the six retained designs are those of Table \ref{tab:performance}, given here to three decimal places.
```

**Proposed**

```qmd
The MSEs and coverages of the six retained designs are those of Table \ref{tab:performance}.
```

### C3T-57. SRS in the adjacent-pair order (line 961)

**Recommendation:** Accept. From T1. Pairs with C3X-2.

**Current**

```qmd
These tests compared each pair of designs adjacent in the mean-MSE order, and each design with SRS within each spillover regime,
```

**Proposed**

```qmd
These tests compared each pair adjacent in the mean-MSE order, with SRS in its place in that order, and each design with SRS within each spillover regime,
```

### C3T-58. Rook block counts (line 1130)

**Recommendation:** Accept. The counts (97 + 60 + 2 + 1 = 160) rank the six designs without SRS. From SR-11 block 2, keeping the chapter's two smaller counts.

**Current**

```qmd
Incidence-Guided Saturation had the lowest MSE in 97 of the 160 rook blocks (61\%),
```

**Proposed**

```qmd
Among the six detailed designs, Incidence-Guided Saturation had the lowest MSE in 97 of the 160 rook blocks (61\%),
```

### C3T-69. Rook 50-unit tests with SRS in the order (line 1132)

**Recommendation:** Accept. Added by the independent review. 2c T3 (approved with T1 and T2, which C3T-57 and C3X-2 carry) added the Balanced Quartiles against SRS comparison to the rook 50-unit tests. Under rook contiguity SRS (0.085) falls between Balanced Quartiles (0.084) and $2\times2$ Blocking (0.128). Numbers from `results/srs_benchmark/srs_benchmark_summary.txt`, rook oracle section 3: Balanced Quartiles minus SRS $-0.001329$ $(-0.00294, 0.000278)$, Wilcoxon $p = 0.105$, frac_better 0.64 (32 of 50), and $2\times2$ Blocking frac_better 0.00 (SRS lower in all 50). They match the SI (unabridged branch) word for word except for the chapter's interval format. The paragraph already gives the $2\times2$ Blocking and Isolation Buffer result (36 of 50), so "every other" has the SI's scope.

**Current**

```qmd
Incidence-Guided Saturation nevertheless had the lower MSE in 48 of the 50 independent configuration-by-surface units against Balanced Quartiles (test methods in Appendix \ref{app:supp-comparisons}), and its paired comparisons with SRS within each regime appear in @sec-vs-srs.
```

**Proposed**

```qmd
Incidence-Guided Saturation nevertheless had the lower MSE in 48 of the 50 independent configuration-by-surface units against Balanced Quartiles (test methods in Appendix \ref{app:supp-comparisons}), and its paired comparisons with SRS within each regime appear in @sec-vs-srs. Balanced Quartiles and SRS were not separated. Balanced Quartiles had the lower MSE in 32 of the 50 units (mean difference $-0.001$, $t_{49}$ 95\% interval $(-0.003, 0.000)$, Wilcoxon $p = 0.11$). Every other adjacent pair, including SRS against $2\times2$ Blocking, differed in all 50 units.
```

### C3T-59. Rook Checkerboard estimates (line 1134)

**Recommendation:** Accept. From NY-20 and HW-16 block 1, which also says what queen Checkerboard kept.

**Current**

```qmd
Under rook contiguity, Checkerboard's estimates were precise but centered on $\beta - \psi$ rather than $\beta$. Under queen contiguity they were centered on $\beta$ but imprecise.
```

**Proposed**

```qmd
Under rook contiguity, Checkerboard's estimates were precise but centered on $\beta - \psi$. Under queen contiguity they were imprecise under control-only spillover but kept nominal coverage.
```

### C3T-60. Run sizes (line 1232)

**Recommendation:** Accept. SH-1 puts the bold run-in head on its own line (it is the chapter appendix's only bold run-in head), and HW-21 removes the pipeline narration. Every number stays.

**Current**

```qmd
\applabel{A11}{app:nc-results}**Run sizes.** The main run generated 7,598,000 independent outcome fits. Every setting was complete and met its prespecified precision targets for mean MSE and coverage. The largest relative Monte Carlo standard error of mean MSE was 0.049 and the largest Monte Carlo standard error of coverage was 0.00892. No setting produced warnings, aliased columns, estimates of $\rho$ at the boundary of its parameter space or estimation failures. The tail refinement added 2,280,000 independent outcomes. Counting the outcomes reused from the main run only once gives 9,878,000 distinct outcomes across the two stages.
```

**Proposed**

```qmd
\applabel{A11}{app:nc-results}**Run sizes**

The main analysis used 7,598,000 independent outcome fits. Every setting met its precision targets for mean MSE and coverage. The largest relative Monte Carlo standard error of MSE was 0.049, and the largest Monte Carlo standard error of coverage was 0.00892. We recorded no warnings, aliased columns, boundary estimates of $\rho$ or estimation failures. The tail refinement added 2,280,000 independent outcomes, giving 9,878,000 distinct outcomes across the two stages.
```

### C3T-61. Annual table caption (line 1279)

**Recommendation:** Accept. Plural, as in SN-22.

**Current**

```qmd
Therefore, their rows repeat and are not independent annual confirmation.}
```

**Proposed**

```qmd
Therefore, their rows repeat and are not independent annual confirmations.}
```

## Batch 9 · Participation, power and sensitivity appendix

### C3T-62. Illustration of the individual-level effect (line 1359)

**Recommendation:** Accept. From SR-7, which states that $d = 0.3$ is about one college-level residual standard deviation and that the powers assume $\rho = 0.5$.

**Current**

```qmd
With $n_c = 25$ and an ICC of 0.05, values chosen for illustration rather than estimated for this trial, $\sigma_c = 0.30$. Table \ref{tab:nc-translation} gives the results for other values. Under both-arms spillover, an individual-level effect of $d = 0.3$ then gives power 0.71 to 0.72 with 29 colleges and 0.96 with 58, for SRS and for both saturation designs.
```

**Proposed**

```qmd
As an illustration, suppose $n_c = 25$ trained professionals per college and an ICC of 0.05. Then $\sigma_c = 0.30$, and an individual-level effect of $d = 0.3$ is close to one college-level residual standard deviation. Table \ref{tab:nc-translation} gives the results for other values. With $\rho = 0.5$, this effect gives power 0.71 to 0.72 with 29 colleges and 0.96 with 58 under both-arms spillover, for SRS and for both saturation designs.
```

### C3T-63. Colleges-needed caption (line 1368)

**Recommendation:** Accept. From HW-27 block 2 (no colon labels) and SR-6 (residual units).

**Current**

```qmd
(minimum detectable effect, MDE, in college-level standard deviations), corrected $t$ test, random volunteering, $\rho = 0.5$, spillover $\psi = 0.8\beta$, mean over the 2018--2021 planning years. $>58$: 80\% power not reached with all 58 colleges in at least one year. $>2$: 80\% power not reached at the largest simulated effect, $\beta = 2$.}
```

**Proposed**

```qmd
(minimum detectable effect, MDE) in college-level residual standard deviations, corrected $t$ test, random volunteering, $\rho = 0.5$, spillover $\psi = 0.8\beta$, mean over the 2018--2021 planning years. An entry of $>58$ means estimated power stayed below 80\% with all 58 colleges in at least one year. An entry of $>2$ means estimated power stayed below 80\% at the largest simulated effect, $\beta = 2$.}
```

### C3T-64. Interval calibration (line 1419)

**Recommendation:** Accept. From HW-28. It states the target once, lists each interval's coverage, and replaces "narrowly missed the rule there, which every other method also missed" with what happened. The chapter's four-decimal proportions stay (the SI prints percentages), and its two sentences on how the method was chosen stay.

**Current**

```qmd
We tried several corrections and adopted a method only if its coverage lay between 0.9400 and 0.9600 for every design and regime. On the stored fits of the main case study ($\tau = 0$), the Wald interval covered 0.9260 to 0.9360 by design and regime. A $t$ reference distribution alone gave 0.9320 to 0.9410, a degrees-of-freedom-corrected standard error alone 0.9350 to 0.9440, and the two together 0.9410 to 0.9490. In a separate run of 500 fits per setting at 58 colleges, parametric percentile bootstrap intervals with 499 replicates covered 0.9010 to 0.9360, parametric bootstrap-$t$ intervals 0.9320 to 0.9580 and profile-likelihood intervals 0.9180 to 0.9470. Only the $t$ interval with the corrected standard error met the rule on the stored fits. In the separate 500-fit run it covered 0.9280 to 0.9580 and narrowly missed the rule there, which every other method also missed.
```

**Proposed**

```qmd
We compared interval corrections against a coverage target of 0.9400 to 0.9600 for every design and regime. For the main case study ($\tau = 0$), coverage across the eight designs and SRS was 0.9260 to 0.9360 with the Wald interval, 0.9320 to 0.9410 with a $t$ reference distribution alone, and 0.9350 to 0.9440 with a degrees-of-freedom correction to the standard error alone. The two corrections together gave 0.9410 to 0.9490, the only interval that met the target. In a separate analysis with 500 fits per setting at 58 colleges, coverage was 0.9010 to 0.9360 for parametric percentile bootstrap intervals with 499 replicates, 0.9320 to 0.9580 for parametric bootstrap-$t$ intervals, and 0.9180 to 0.9470 for profile-likelihood intervals. The corrected $t$ interval covered 0.9280 to 0.9580. Every method in that analysis had coverage below 0.9400 in at least one design and regime.
```

### C3T-65. Case-study Friedman tests (line 1509)

**Recommendation:** Accept. From NY-25 and HW-30 block 4.

**Current (edit 1)**

```qmd
rejected equal ranks of the nine designs in all four combinations of regime and $\tau$
```

**Proposed (edit 1)**

```qmd
rejected equal ranks of the eight designs and SRS in all four combinations of regime and $\tau$
```

**Current (edit 2)**

```qmd
These tests treat the 24 blocks as independent, although the years share one map and some designs share their simulations across years. Therefore, they describe the consistency of the ordering rather than independent evidence.
```

**Proposed (edit 2)**

```qmd
These tests treat the 24 blocks as independent, although the four years share one map and some designs and SRS reuse their simulations across years (@sec-nc-outcome). Therefore, they describe how consistent the ordering is.
```

### C3T-66. Spatial-model results under their own heading (line 1515)

**Recommendation:** Accept. The spatial-model and Bayesian results now sit under "Unequal Numbers of Trained Professionals". HW-29 gave the SI three headed paragraphs for them, and SR-17 replaced "converged". The chapter version adds one heading and splits the paragraph in the SI's order (maximum likelihood CAR, Bayesian CAR, Bayesian spatial lag). Every number stays, in the chapter's format.

**Current**

```qmd
With the longer chains, 83.1\% of the Bayesian CAR chains at $\rho = 0.5$ converged, against 46.2\% for the shorter chains, but 12.5\% still had an absolute Geweke $z$ of at least 1.96, against about 5\% expected by chance. Among converged chains, the Bayesian CAR estimate was biased by $+0.138$ under both-arms spillover and $+0.219$ under control-only spillover, with credible-interval coverage of 0.904 and 0.877. For maximum likelihood on the same data, the biases were $+0.023$ and $+0.028$ and the corrected $t$ coverage 0.950 and 0.946. Of the Bayesian spatial lag chains, 94.9\% converged. The median effective sample size of $\beta$ was 8,000 (minimum 5,897), and 5.1\% of chains had an absolute Geweke $z$ of at least 1.96. The model's posterior means correlated 0.999 with the maximum likelihood estimates (mean absolute difference 0.017). Its bias was at most 0.039 in absolute value in every combination of $\rho$ and regime, and its credible intervals covered 0.946 to 0.953 over all chains (0.945 to 0.953 over converged chains). The Spearman correlation of its design ordering with that of maximum likelihood averaged 0.96 to 1.00 across the four combinations of $\rho$ and regime (0.83 in the lowest single block), and the best design agreed in every block. In the maximum likelihood comparison of the spatial lag and CAR models, the saturation designs' MSE under control-only spillover was 0.76 to 0.82 times that of SRS in every combination of generating model and fitted model. When the CAR model was fit to spatial lag data, its estimates were biased by 0.035 and 0.085 at $\rho = 0.2$ and by 0.110 and 0.205 at $\rho = 0.5$ under both-arms and control-only spillover. When it was fit to CAR data, its bias was at most 0.005.
```

**Proposed**

```qmd
### Spatial Model and Bayesian Fits {.unnumbered .unlisted}

In the maximum likelihood comparison of the spatial lag and CAR models, the saturation designs' MSE under control-only spillover was 0.76 to 0.82 times that of SRS in every combination of generating model and fitted model. When the CAR model was fit to spatial lag data, its estimates were biased by 0.035 and 0.085 at $\rho = 0.2$ and by 0.110 and 0.205 at $\rho = 0.5$ under both-arms and control-only spillover. When it was fit to CAR data, its bias was at most 0.005.

With the longer chains, 83.1\% of the Bayesian CAR chains at $\rho = 0.5$ met the diagnostics for $\beta$, compared with 46.2\% of the shorter chains, but 12.5\% still had an absolute Geweke $z$ of at least 1.96, against about 5\% expected by chance. Among chains that met those diagnostics, the Bayesian CAR estimate was biased by $+0.138$ under both-arms spillover and $+0.219$ under control-only spillover, with credible-interval coverage of 0.904 and 0.877. Across all maximum likelihood fits to the same datasets, the biases were $+0.023$ and $+0.028$ and the corrected $t$ coverage 0.950 and 0.946.

Of the Bayesian spatial lag chains, 94.9\% met the diagnostics for $\beta$. The median effective sample size of $\beta$ was 8,000 (minimum 5,897), and 5.1\% of chains had an absolute Geweke $z$ of at least 1.96. The model's posterior means correlated 0.999 with the maximum likelihood estimates (mean absolute difference 0.017). Its bias was at most 0.039 in absolute value in every combination of $\rho$ and regime, and its credible intervals covered 0.946 to 0.953 over all chains (0.945 to 0.953 over chains meeting the diagnostics). The Spearman correlation of its design ordering with that of maximum likelihood averaged 0.96 to 1.00 across the four combinations of $\rho$ and regime (0.83 in the lowest single block), and the best design agreed in every block.
```

## Checked, no change

- Line 78, "$\beta$ is the direct treatment effect of interest": the same paragraph defines $\beta$ before propagation (C3T-16), and H1's wording adds nothing.
- Line 78, N3 pointer where $\tau$ is first set: the chapter sets $\tau$ in @sec-simulation, where C3N-2 adds the pointer.
- Line 93, base rate (HW-1): the chapter's fuller wording (between the two Mirzaei rates, close to their average) is accurate.
- Line 138, caption "It is not one of the simulation's incidence surfaces" (HW-2): a plain fact the long form keeps.
- Line 148, interval definition (SN-14): already split into short sentences with $p$ stated plainly.
- Line 150, Checkerboard alias detail (HW-3): the long form keeps the identities the SI moved to its diagnostics.
- Lines 158, 160, 172, 221, 239–243, 382 (HW-1, SN-25, SR-11 block 1, SN-18, SN-10, HW-5): already match the approved wording, and line 221 already says "of the six designs".
- Line 166, seed scheme (HW-4): the chapter's short pointer to A1 is clear.
- Line 390, map caption (NY-7): already plain ("shows county boundaries ... only").
- Lines 466 and 468 (HW-10): already sentence-level, and only lines 470 and 472 change (C3T-6).
- Line 484, Checkerboard identifiability contrast (HW-22 block 2): already aligned.
- Lines 541 and 547, `fig:alloc-dist` text and caption (HW-24 block 2): the chapter caption keeps the outcome counts and the percentile method that the SI moved or cut.
- Line 583, rank-based designs "within 2% ... except High Incidence Focus (4% and 8%)" (HW-30 block 3): more precise than the SI's "within 8%".
- Line 593 (HC-12): already split.
- Remaining "NC" in captions, the table header at line 1456, and lines 370, 412, 547, 599, 795: the chapter defines "North Carolina (NC)" at line 55, and every edited sentence above uses "North Carolina" or "case study".
- Line 753, adjacent-pair formula: unchanged by SRS joining the order.
- Line 799, spatial-blocking groups (HW-6 counterpart): says "rounded" and sits after the rule C3T-5 states.
- Line 851, "This section gives the four models, their priors, ..." (SI cut): a useful roadmap in the long form.
- Line 875 (S1): the chapter body's benchmark table includes both variants, and "included in the nine-design benchmark" is accurate here.
- Line 885 and the two variant tables (2d, PC-2): already aligned in 2d.
- Captions at lines 1001, 1024, 1064, 1089, 1159, 1171 and 1451 (I5–I9): already say "identifiable".
- Line 1272, checks of warning and failure records (HW-21 block 2): the chapter's own sentence is accurate.
- Line 1344, `fig:power-n` caption: keeps detail the SI moved into prose.
- Lines 1474 and 1480 (H84, H85): no wording change applies.
- Line 1505, High Incidence Focus share "0.75 to 0.77": same values as the SI's 75% to 77%, in the chapter's format.
- Line 621, limitation on shared draws: chapter-only text with no approved counterpart.
- HW-26, SN-12 (`tab:nc-summary-sens`): SI-only table.
- `real_sud/primary_setting_ratios.pdf` (F1), `primary_yearly_mse.pdf` and `population_shares.pdf`: copies sit in the chapter's figure folder, but the chapter includes none of them.

## Checks

A Python script parsed this file and tested every item against `Thesis_Chapter.qmd` as of 9bea05f.

- **Current blocks.** All 115 Current blocks in the 77 text and table items (66 C3T, 9 C3N, 2 C3X) occur exactly once in `Thesis_Chapter.qmd`, and no two of them overlap. Every edit can therefore be applied on its own. The 7 C3F items have no chapter text. Failures: none.
- **Numbers.** A script compared the numbers in each Current block with those in its Proposed block. The only drops are intended: C3T-2 (0.631 and 1.083 rounded to 0.63 and 1.08, SR-18), C3X-2 edit 3 (the Incidence-Guided Saturation against Balanced Quartiles row, $-0.040$ $(-0.043, -0.037)$, leaves because the two are no longer adjacent, T2), C3T-9 ("$g > 0$" goes with the sign argument, SR-5), C3T-19 edit 2 ("all 25 draws on a surface", kept as "25 times per incidence surface" in edit 1), C3T-21 (a repeated "250 distinct noise vectors"), C3T-35 edit 1 ("100 repeated copies", NY-style cut in HW-8), C3T-54 edits 1 and 2 (the "met the targets at $J=100$ and $R=100$" sentence, kept in the body at line 456 by C3T-35, and the $n_a/J^2$ contrast removed by NY-14), C3T-64 (a repeated "500-fit run") and C3T-27 ("$\pm 1$" written as "one"). Two flags were tokenization only (C3T-5 "2," and C3T-62 "58,"). New numbers come from the applied SI text or, for C3X-2, from `results/srs_benchmark/srs_benchmark_summary.txt`; none was typed from memory.
- **Style.** Every Proposed block was scanned for semicolons, em dashes, ", so" joins, "rather than", "not only", "is fitted", sentence-initial "Because", mid-sentence "however" or "therefore", prose colons and bare "NC". The one semicolon is inside the existing citation `[@aronow_estimating_2017; @forastiere_identification_2021]`. "Fitted" appears only before a noun. Design names match the canonical set.
- **New-figure check.** `fig_igs_re_vs_srs.pdf`: 6.5 × 5.5 in, text 8–9 pt, ArialMT embedded. It passes. `fig_incidence_rejection_queen.pdf`: 8.5 × 10 in, text 6–10 pt (about 4.6–7.6 pt at 6.5 in), Helvetica and Symbol not embedded. It fails (C3F-2). Measured with `pdfinfo`, `pdffonts` and PyMuPDF span sizes.
- **Other chapter figures.** Text extraction of every PDF in `paper/thesis_chapter/figures/` found the four problems listed as C3F-4 to C3F-7 and no others ("q90", "rho=", "gamma", band notes, minus-glyph hyphens). The PNG example map was inspected by eye.
- **Not run.** No manuscript, figure script or verifier was run or edited. Quarto rendering and `verify_manuscripts.py` follow application.

## Independent review (2026-11-26)

A fresh reviewer checked every item against the applied SI (unabridged branch), the CTJ, the 2c, 2d, 2g-merged, 2h and 4e/3c sign-offs, the hunk list and the saved results. It did not edit the chapter, SI, CTJ, code or figures.

### How it was checked

- **Fidelity.** Each Proposed sentence was compared, after normalizing references and math delimiters, with the SI, the CTJ and the approved sign-off text. Every sentence that is not verbatim is a chapter adaptation (chapter labels, `@sec-` references, ", so" splits), a chapter sentence the item says it keeps, or a numeric format change (proportions for the SI's percentages). No Proposed block contains text from a curated-only `\else` branch.
- **Numbers.** C3X-2's new rows match `results/srs_benchmark/srs_benchmark_summary.txt` (queen oracle, section 3: $-0.03781$ $(-0.04046, -0.03516)$, and Balanced Quartiles minus SRS $+0.00233$ $(-0.00073, 0.00540)$, Wilcoxon $p = 0.153$, SRS lower in 29 units). The A13 and A14 numbers match SI Sections 3.2, 2.2 and 13 as applied.
- **Coverage.** Every old SI sentence that the hunks changed was matched against the chapter. Each approved 2g, 2h and 2c ID was traced to an item or to "Checked, no change". The gaps found are the five items added below.

### Changes made to this file

- **C3T-67 (new, Batch 5, line 362).** HW-14: "The two were level at $\beta = 2$ (0.160 each)" becomes "At $\beta = 2$ both rounded to 0.160", and "lowest in pooled MSE at every value of $\beta$" is scoped to the six designs.
- **C3T-68 (new, Batch 5, line 370).** HW-19 block 1, which C3T-28 left out. The all-200-units claim is limited to the 90th percentile and worst-decile mean, with "assignment-specific" and "both-arms spillover".
- **C3T-71 (new, Batch 5, line 448).** HW-7 removed the colon before the case-study outcome equation.
- **C3T-70 (new, Batch 6, line 480).** SN-8: MSE units at the first case-study MSE values.
- **C3T-69 (new, Batch 8, line 1132).** 2c T3: the rook 50-unit result for Balanced Quartiles against SRS (32 of 50, $-0.001$, $(-0.003, 0.000)$, $p = 0.11$, checked against the rook section of `srs_benchmark_summary.txt`) and "Every other adjacent pair, including SRS against $2\times2$ Blocking, differed in all 50 units." T1 and T2 were already carried (C3T-57, C3X-2), T3 was not.
- **C3N-8.** Added `\FloatBarrier` after the A13 figure, as the other appendix sections end, so the rejection figure cannot float into A14.
- **Summary table and the ID note** list the five new items. The "Checks" section above predates them. The current counts are 120 Current blocks in 82 text and table items (71 C3T, 9 C3N, 2 C3X).

### Flagged for the author, not changed

- **C3T-68.** The chapter sentence after it ("These results come from the extreme parameter settings only and give qualified support to saturation on downside risk.") stays. HW-19 cut it from the SI. It is accurate, but you may prefer to drop it for alignment.
- **Chapter-only contrasts with no approved counterpart.** Line 412 ("an approximation of strict interspersion rather than an exact checkerboard"), line 623 ("a general design study rather than a calibration") and line 625 ("not an observed education intervention or demonstrated reduction in SUD", "rather than a sample of healthcare professionals"). NYC-6 rewrote the article's matching scope sentence ("Its outcomes are simulated and give no evidence about whether education reduces SUD."), which could serve line 625. These sentences are left for your human-writer pass.
- **Line 158** defines queen and rook again after C3T-14 defines them at line 61. The repetition is harmless in the long form.
- **C3T-12 and C3T-66** keep chapter-only detail on purpose (the dips larger than Monte Carlo error, and the CAR bias explanation, which the chapter gives in the body at lines 556 to 560). No content is lost.

### Mechanics check (after these edits)

A script parsed all 89 items and confirmed the following. Every one of the 120 Current blocks occurs exactly once in `Thesis_Chapter.qmd`, no two overlap, and each Current has one Proposed. Applying all 120 edits in sequence succeeds, and the result has `\applabel` A11 to A14. Every `\ref` and `@sec-` in the Proposed blocks resolves to an existing label or to one the items add. The scan of the Proposed blocks for semicolons outside citations, em dashes, ", so", "rather than", "not only", ", not", "fitted" as a verb, sentence-initial "Because", mid-sentence "however" or "therefore", prose colons and bare "NC" found nothing.

## Author decisions

Pending.
