# Phase 2g sign-off for the statistical and trial-planning reviews

Author: Andrew Walther. Drafted 2026-11-25.

Two readers independently examined the current curated SI before any prose
proposals were applied. The statistical read compared the methods and claims
with the exhibits, saved results and analysis code. The trial-planning read
followed the case study from the geographic inputs through recruitment,
assignment, spillover assumptions and power.

Each item includes its source passage and exact proposed replacement. When an
item builds on an HW item, its current block is the text that would result from
approving that HW item. Decide both IDs together. The SI, article and chapter
remain at their current versions while these proposals are under review.


## Batch 1 · Interpretation and trial-planning decisions


### SR-1. Binary spillover and identifiable treatment effects

**Accept.** The binary control-only ratios condition on identifiable fits, and the set of identifiable assignments can differ by design. See the identification check in application/code/methods_ext_sensitivity.R and the two saved a7 summaries cited by the table.


**Current**

```latex
The design comparison did not depend on the spatial model used to generate or fit the data, except for a CAR model fit to spatial lag data under strong spatial dependence. It also held under two other forms of spillover exposure, and formal rank tests found a consistent ordering of the strategies within each regime.
```

**Proposed**

```latex
The broad comparison under proportional spillover held under the count
form and across most spatial-model comparisons. Under binary spillover,
identification depended on the assignment. Formal rank tests summarize
the ordering of the eight designs and SRS within each regime.
```


**Current**

```latex
\textbf{Spillover form.} Table \ref{tab:nc-forms} compares the three forms. Under the binary form, Isolation Buffer is never identifiable, because every control has a treated neighbor, and under control-only spillover no strategy was identifiable in every fit. The saturation designs were identifiable most often there (in 89.7\% and 91.7\% of fits) and, among identifiable fits, had MSE 0.70 to 0.71 times that of SRS. Under the count form, the saturation designs still had the lowest MSE under control-only spillover.
```

**Proposed**

```latex
\textbf{Spillover form.} Table \ref{tab:nc-forms} compares the three
forms. Under binary exposure, Isolation Buffer assigns every treated
cluster zero exposure and every control cluster exposure one, so the
model cannot separate $\beta$ from spillover in either regime. Under
control-only spillover, $\beta$ was unidentified in some fits under
every assignment rule, including SRS. The two saturation designs
identified $\beta$ most often, in 89.7\% and 91.7\% of fits. Among
the identifiable fits for each design and for SRS, their MSEs were
0.70 to 0.71 times that of SRS. These ratios describe the subsets that
were identifiable under each rule. Under the count form, the
saturation designs still had the lowest MSE under control-only spillover.
```


**Current**

```latex
\caption{Mean MSE relative to simple random sampling (SRS) under three spillover forms, 58 colleges, $\tau = 0$, mean over 16 blocks (year $\times$ $\rho$ $\times$ $\psi$). n.i.: not identifiable in every fit (the spillover covariate is aliased with the intercept and treatment in some or all allocations).}
```

**Proposed**

```latex
\caption{Mean MSE relative to simple random sampling (SRS) under three spillover forms, 58 colleges, $\tau = 0$, mean over 16 blocks (year $\times$ $\rho$ $\times$ $\psi$). An entry of n.i.\ means that the treatment coefficient $\beta$ was unidentifiable in at least one sampled assignment.}
```


**Current**

```latex
\textbf{Spillover form.} Under the binary form, $\mathbf{S}(\mathbf{Z}) = \mathbb{I}\{\text{any treated queen neighbor}\}$, times $\mathbf{1} - \mathbf{Z}$ under control-only spillover. Under the count form, $\mathbf{S}(\mathbf{Z}) = (1/\bar d) \times$ (number of treated queen neighbors), where $\bar d = 4.93$ is the mean number of queen neighbors, so a college of average degree with every neighbor treated has $\mathbf{S}(\mathbf{Z}) = 1$ and receives $\psi$, as under the proportional form.
```

**Proposed**

```latex
\textbf{Spillover form.} Under the binary form, $\mathbf{S}(\mathbf{Z}) = \mathbb{I}\{\text{any treated queen neighbor}\}$, times $\mathbf{1} - \mathbf{Z}$ under control-only spillover. Under the count form, $\mathbf{S}(\mathbf{Z}) = (1/\bar d) \times$ (number of treated queen neighbors), where $\bar d = 4.93$ is the mean number of queen neighbors, so a college of average degree with every neighbor treated has $\mathbf{S}(\mathbf{Z}) = 1$ and receives $\psi$ if exposed, as under the proportional form. Under control-only spillover, both the binary and the count exposure are multiplied by $\mathbf{1}-\mathbf{Z}$.
```


### SR-2. Meaning of the incidence coefficient in the new Section 3.2

**Accept.** This touches newly approved text. A one-unit contrast in scaled incidence differs from the realized range, and omitting the structural incidence term does not remove associations through assignment. Spillover adjustment applies to the oracle model.


**Current**

```latex
The incidence coefficient \(\tau\) in Equation~\ref{eq:outcome} is a regression coefficient, not a correlation. Incidence \(\mathbf{X}\) lies between 0 and 1 and \(\sigma = 1\), so \(\tau\) is the shift in the outcome, in residual standard deviations, from the lowest to the highest incidence, holding treatment and spillover fixed and before the spatial multiplier spreads it to neighbors. At \(\tau = 1\) that shift equals the primary intervention effect \(\beta = 1\). The value \(\tau = 0\) means that incidence has no association with the outcome, and \(\tau\) can exceed 1 or be negative.
```

**Proposed**

```latex
The incidence coefficient \(\tau\) in Equation~\ref{eq:outcome} gives
the change in the outcome's linear predictor for a one-unit increase in
scaled incidence \(\mathbf{X}\), holding treatment and spillover fixed.
The change is in residual standard deviations because \(\sigma=1\), and
is measured before spatial propagation. At \(\tau=1\), a one-unit
incidence contrast contributes as much to the linear predictor as the
primary intervention effect \(\beta=1\). At \(\tau=0\) incidence makes
no direct contribution to that predictor. The coefficient can exceed
1 or be negative.
```


**Current**

```latex
Before a trial, the investigator chooses \(\tau\) for planning, because there are no outcomes yet. After the trial, the investigator has both the prior incidence used for allocation and the measured outcomes, so \(\tau\) can be estimated in the same spatial lag model that estimates \(\beta\). Both grid models already include \(\mathbf{X}_k\) (Section~\ref{estimation}), so every grid fit gives \(\hat\tau\) beside \(\hat\beta\). The estimate describes how prior incidence is associated with the later outcome after adjustment for treatment and spillover. Incidence is not randomized, so \(\hat\tau\) is an adjusted association rather than an effect of incidence.
```

**Proposed**

```latex
Before outcomes are observed, investigators specify \(\tau\) for the
planning model. After a trial, the same spatial lag model can estimate
\(\tau\) and \(\beta\) from the observed outcomes and the incidence
used for allocation. Both grid models include \(\mathbf{X}_k\)
(Section~\ref{estimation}), so every fit provides both estimates.
The oracle estimate describes the association of prior incidence
with the outcome after adjustment for treatment and spillover. The
non-oracle model adjusts for treatment and omits spillover.
Incidence is observed, so these coefficients describe associations.
```


### SR-3. Incidence sensitivity and error in the spatial lag coefficient

**Accept.** This touches newly approved Section 3.2. With the spatial lag coefficient held at its true value, changing the incidence coefficient leaves the treatment estimate unchanged.


**Current**

```latex
Changing \(\tau\) barely moved \(\hat\beta\). Over the five values of \(\tau\), no design's mean MSE for \(\beta\) changed by more than 3.4\% from its value at \(\tau = 1\) in any incidence configuration, neighbor type or regime. Both fitted models include \(\mathbf{X}\), so a change in \(\tau\) is absorbed by the coefficient on \(\mathbf{X}\) and reaches \(\hat\beta\) only through the estimate of \(\rho\). The spillover regime still decided which designs led. Under control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest MSE in every comparison. Under both-arms spillover, the best design's MSE was at most 5\% below that of SRS. The lowest-MSE design changed with \(\tau\) in 6 of the 20 comparisons, and in 5 of these the two saturation designs, whose MSEs were nearly equal, exchanged first place.
```

**Proposed**

```latex
Changing \(\tau\) barely moved \(\hat\beta\). Over the five values of \(\tau\), no design's mean MSE for \(\beta\) changed by more than 3.4\% from its value at \(\tau = 1\) in any incidence configuration, neighbor type or regime. Both fitted models include \(\mathbf{X}\). Changing \(\tau\) leaves the fitted intervention coefficient unchanged when \(\rho\) is held at its true value. Error in estimating \(\rho\) can pass part of the change into \(\hat\beta\) (Section~\ref{app:technical-notes}). The spillover regime still decided which designs led. Under control-only spillover, Incidence-Guided Saturation or Saturation Quadrants had the lowest MSE in every comparison. Under both-arms spillover, the best design's MSE was at most 5\% below that of SRS. The lowest-MSE design changed with \(\tau\) in 6 of the 20 comparisons, and in 5 of these the two saturation designs, whose MSEs were nearly equal, exchanged first place.
```


### SR-4. Spatial lag and CAR mean structures

**Accept.** Builds on HW-11. The two models share a linear predictor. The spatial lag model propagates it through the weights matrix, as the simulation code does.


**Current**

```latex
In the reverse direction, we generated data from a CAR model with the same mean and fit both models.
```

**Proposed**

```latex
We also generated data from a CAR model with the same linear predictor and fit both models.
```


**Current**

```latex
\textbf{Estimation models.} Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fit by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```

**Proposed**

```latex
\textbf{Estimation models.} The two spatial model families use the linear predictor $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and use the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fit by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. Thus the spatial lag outcome has conditional mean $(\mathbf{I}-\rho\mathbf{W})^{-1}\boldsymbol{\mu}$ given the assignment, whereas the CAR outcome has conditional mean $\boldsymbol{\mu}$. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```


### SR-5. Observed contribution of spatial lag bias

**Accept.** A positive average of the spatial lag error alone does not determine the sign of its average product with g or how the product changes with the effect size.


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


### SR-6. Hypothetical education endpoint and the power scale

**Accept.** Builds on HW-10 and HW-27. The planning model uses residual standard deviations before spatial propagation, consistent with the CTJ. The example assumes a professional-level score, whose final definition is a trial decision.


**Current**

```latex
All 58 colleges will be invited to the proposed extension. To allow for
incomplete participation, we repeated the case study for
$N \in \{15,20,25,29,30,35,40,45,50,58\}$ participating colleges
under the education outcome model with $\tau=0$. Under random
volunteering, we drew 50 subsets for each $N$ and used those same
subsets for all designs, SRS and planning years. At $N=20$ and $N=29$
we also examined volunteering probabilities that increased linearly
with incidence rank, so the highest-burden college was three times as
likely to join as the lowest. We reapplied each assignment rule on the
participants and required at least two colleges in each saturation region.

Spatial dependence links participating colleges. A participant with no
participating neighbors has a zero row in $\mathbf{W}_P$. For spillover,
we divided the number of treated participating neighbors by the number
of geographic neighbors:
$e_j(\mathbf{Z})=\sum_{k\in\mathcal{N}(j)}Z_k/|\mathcal{N}(j)|$.
A neighbor who does not participate counts as untreated.

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

We used the same subsets, assignments and noise at all effect sizes to
keep the simulated power curves smooth in $\beta$. We used the
corrected $t$ test to estimate power at each setting, the number of
colleges needed for 80\% power, and the minimum detectable effect
(MDE), the smallest effect detectable with 80\% power at each $N$.
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
with incidence rank, so the highest-burden college was three times as
likely to join as the lowest. We reapplied each assignment rule on the
participants and required at least two colleges in each saturation region.

Spatial dependence links participating colleges. A participant with no
participating neighbors has a zero row in $\mathbf{W}_P$. For spillover,
we divided the number of treated participating neighbors by the number
of geographic neighbors:
$e_j(\mathbf{Z})=\sum_{k\in\mathcal{N}(j)}Z_k/|\mathcal{N}(j)|$.
A neighbor who does not participate counts as untreated.

We set $\psi=\kappa\beta$ with $\kappa \in \{0.5,0.8\}$ and used nine
values of $\beta$ from 0 to 2 in units of the college-level residual
standard deviation before spatial propagation. We crossed these settings with $\rho \in \{0,0.5\}$,
both regimes, the eight designs and SRS, and four planning years. With
fewer than 58 colleges, each of the 50 subsets used 20 assignments and
10 outcomes per assignment (10,000 fits per setting). A design with a
single possible assignment on a subset used 200 outcomes for that
assignment. With 58 colleges, each setting used 100 assignments and 100
outcomes as in the main analysis and reproduced its MSE and coverage
within Monte Carlo error.

We used the same subsets, assignments and noise at all effect sizes to
keep the simulated power curves smooth in $\beta$. We used the
corrected $t$ test to estimate power at each setting, the number of
colleges needed for 80\% power, and the minimum detectable effect
(MDE), the smallest effect detectable with 80\% power at each $N$.
```


**Current**

```latex
The trial outcome is a score measured on each trained professional and averaged within each college. We write a professional's score as a college effect plus individual variation, with total variance 1, so that effects are in individual standard-deviation units. The intraclass correlation coefficient (ICC) is the share of that variance that lies between colleges, so the between-college variance is the ICC and the within-college variance is $1 - \text{ICC}$. Averaging the scores of $n_c$ professionals leaves the college part unchanged and divides the individual part by $n_c$, so the college mean has standard deviation $\sigma_c = \sqrt{\text{ICC} + (1 - \text{ICC})/n_c}$ \citep{donner_design_2000, hayes_cluster_2017}. A college-level effect of $\beta$ standard deviations is therefore an individual-level standardized effect $d = \beta\sigma_c$. With $\tau = 0$, multiplying the effect and the noise by the same constant leaves $\hat\beta/\widehat{\mathrm{SE}}$ unchanged, so power for any $n_c$, ICC and $d$ follows from the simulated college-level results. Rerunning one setting with the effect and the noise both halved reproduced its estimates to within $6\times10^{-8}$. In a sensitivity analysis, the number of trained professionals grew with college population ($n_c$ averaging 25 over participants, ICC $= 0.05$), so the noise variance differed between colleges while the fitted model stayed unweighted, as an analyst unaware of the difference would fit it.
```

**Proposed**

```latex
For planning, suppose the education outcome is a score measured on each trained professional and averaged within each college. We represent a professional's score as a college effect plus individual variation, with total variance 1, so that effects are in individual standard-deviation units. The intraclass correlation coefficient (ICC) is the share of that variance that lies between colleges, so the between-college variance is the ICC and the within-college variance is $1 - \text{ICC}$. Averaging the scores of $n_c$ professionals leaves the college part unchanged and divides the individual part by $n_c$, so the college mean has standard deviation $\sigma_c = \sqrt{\text{ICC} + (1 - \text{ICC})/n_c}$ \citep{donner_design_2000, hayes_cluster_2017}. To express the structural effect before spatial propagation on the individual standard-deviation scale, we use $d = \beta\sigma_c$. With $\tau = 0$, multiplying the effect and the noise by the same constant leaves $\hat\beta/\widehat{\mathrm{SE}}$ unchanged, so power for any $n_c$, ICC and $d$ follows from the simulated college-level results. Rerunning one setting with the effect and the noise both halved reproduced its estimates to within $6\times10^{-8}$. In a sensitivity analysis, the number of trained professionals grew with college population ($n_c$ averaging 25 over participants, ICC $= 0.05$), so the noise variance differed between colleges while the fitted model stayed unweighted, as an analyst unaware of the difference would fit it.
```


**Current**

```latex
college-level standard deviations under both-arms spillover for every design except Checkerboard (1.32) and Isolation Buffer (1.50)
```

**Proposed**

```latex
college-level residual standard deviations under both-arms spillover for every design except Checkerboard (1.32) and Isolation Buffer (1.50)
```


**Current**

```latex
minimum detectable effect (MDE) in college-level standard deviations
```

**Proposed**

```latex
minimum detectable effect (MDE) in college-level residual standard deviations
```


### SR-7. Power assumptions and the participation example

**Accept.** State the spatial dependence and spillover assumptions beside the figures a planner would use. The saved corrected-t sample-size and translation summaries support both examples.


**Current**

```latex
Under both-arms spillover, the design hardly affected the number of colleges a trial needs. Every design except Checkerboard and Isolation Buffer needed 35 to 37 colleges for 80\% power to detect a college-level effect of one standard deviation, and about 19 colleges for 1.5 standard deviations (Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). Checkerboard and Isolation Buffer did not reach 80\% power at $\beta = 1$ even with all 58 colleges.
```

**Proposed**

```latex
With both-arms spillover, $\rho=0.5$ and $\psi=0.8\beta$,
the estimated numbers for 80\% power were 35 to 37 colleges at an
effect of one college-level residual standard deviation and about 19
colleges at 1.5 standard deviations, except for Checkerboard and
Isolation Buffer (Figure \ref{fig:colleges-needed}, Table
\ref{tab:nc-power}). At $\beta=1$, those two designs did not reach
80\% power with all 58 colleges.
```


**Current**

```latex
For a trial with $n_c = 25$ trained professionals per college and an ICC of 0.05, so that $\sigma_c = 0.30$, an individual-level effect of $d = 0.3$ gives power 0.71 to 0.72 with 29 colleges and 0.96 with 58 under both-arms spillover, for SRS and for both saturation designs. Under control-only spillover the same effect gives power 0.49 and 0.53 with SRS, and 0.51 and 0.64 with Incidence-Guided Saturation\iffullsi\ (Table \ref{tab:nc-translation})\fi.
```

**Proposed**

```latex
With $\rho=0.5$ and $n_c = 25$ trained professionals per college and an ICC of 0.05, so that $\sigma_c = 0.30$, an individual-level effect of $d = 0.3$ gives power 0.71 to 0.72 with 29 colleges and 0.96 with 58 under both-arms spillover, for SRS and for both saturation designs. Under control-only spillover the same effect gives power 0.49 and 0.53 with SRS, and 0.51 and 0.64 with Incidence-Guided Saturation\iffullsi\ (Table \ref{tab:nc-translation})\fi.
```


### SR-8. How the number of colleges was interpolated

**Accept.** Builds on HW-30. The saved summaries use the first crossing on the original simulated curve. The mean across years describes this planning scenario and need not attain 80 percent in each year.


**Current**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives, for each design and for SRS, the number of colleges needed for 80\% power and the smallest college-level effect detectable with 80\% power at 29 and 58 colleges. The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not perfectly monotone in $N$ within Monte Carlo error, because the participant subsets differ between values of $N$, and we report the dips rather than smooth them.
```

**Proposed**

```latex
Table \ref{tab:nc-power} gives the estimated number of colleges
needed for 80\% power and the MDE at 29 and 58 colleges. Within each
year, we used linear interpolation between the simulated sample sizes
where estimated power first reached 80\%. When we obtained an estimate
for every year, we averaged those estimates and rounded up. Individual
years can require more colleges than the rounded average. Estimates at
$\rho=0$ differed from those at $\rho=0.5$ by about two colleges or
fewer. We retain the small dips in the simulated power curves. Each
sample size uses its own participant subsets.
```


### SR-9. Size and coverage of the corrected test

**Accept.** Builds on HW-22 and HW-28. The final sentence about nominal size should reflect the observed rejection rates. The two HW items already scope the coverage ranges and distinguish the calibration run.


**Current**

```latex
The colleges-needed numbers above therefore come from a test with the nominal size.
```

**Proposed**

```latex
The sample-size estimates therefore use a test whose rejection rate was close to the nominal 5\% level in these simulations.
```


## Batch 2 · Scope, simulation precision and reporting


### SR-10. Structural spillover and spatial propagation

**Accept.** The spillover term may enter only the controls’ equations and still propagate to treated clusters through the spatial multiplier.


**Current**

```latex
\end{equation} where \(\mathbf{Z} \in \{0,1\}^{100}\) is the treatment assignment
produced by the design under evaluation, \(\beta\) is the direct
treatment effect of interest, \(\mathbf{S}(\mathbf{Z})\) is the spillover term defined in
```

**Proposed**

```latex
\end{equation} where \(\mathbf{Z} \in \{0,1\}^{100}\) is the treatment assignment
produced by the design under evaluation, \(\beta\) is the direct
treatment coefficient before spatial propagation, \(\mathbf{S}(\mathbf{Z})\) is the spillover term defined in
```


**Current**

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

**Proposed**

```latex
\end{equation} where \(\circ\) multiplies the two vectors cluster by cluster. The share of treated neighbors \(\mathbf{W}\mathbf{Z}\) maps the full assignment to each cluster's exposure \citep{aronow_estimating_2017,forastiere_identification_2021}.
Under control-only spillover, the exposure term enters the outcome
equation only for untreated clusters. A larger share of treated
neighbors raises that term for an untreated cluster. The spatial
multiplier can propagate the contribution to treated clusters.
Under both-arms spillover, the exposure term enters the equation for
treated and untreated clusters. We tested four spillover magnitudes,
```


### SR-11. Scope of the grid wins and robustness results

**Accept.** The saved grid scenarios support these numbers among the six detailed designs. Adding SRS or the two variants changes who wins some blocks.


**Current**

```latex
Incidence-Guided Saturation also led block by block. Of the
160 queen blocks at \(\beta = 1\), Incidence-Guided Saturation
had the lowest MSE in 95 (59\%), Balanced Quartiles in 58 (36\%), and
Isolation Buffer in 7 (4\%); no other design had the lowest MSE in any
block.
```

**Proposed**

```latex
Among the six detailed designs, Incidence-Guided Saturation had the
lowest MSE in 95 (59\%) of the 160 queen blocks at \(\beta=1\).
Balanced Quartiles led in 58 (36\%) and Isolation Buffer in 7 (4\%).
```


**Current**

```latex
Incidence-Guided Saturation had the
lowest MSE in 97 of the 160 rook blocks (61\%) and Balanced Quartiles in
60 (38\%).
```

**Proposed**

```latex
Among the six detailed designs, Incidence-Guided Saturation had the
lowest MSE in 97 of the 160 rook blocks (61\%) and Balanced Quartiles in
60 (38\%).
```


**Current**

```latex
the lowest median (0.080) and the lowest worst-case MSE of any design
```

**Proposed**

```latex
the lowest median (0.080) and worst-case MSE among the six detailed designs
```


**Current**

```latex
Incidence-Guided Saturation also had the lowest worst-case MSE
(0.138).
```

**Proposed**

```latex
Incidence-Guided Saturation also had the lowest worst-case MSE
among the six detailed designs (0.138).
```


### SR-12. Non-oracle rankings

**Accept.** The leading pooled design stayed the same in the detailed comparison, while other positions changed.


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

**Accept.** The simulation reports surface-level Monte Carlo standard errors for four metrics. The paired intervals address mean differences and use the pooled variation across the independent surfaces.


**Current**

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

**Proposed**

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


**Current**

```latex
that surface's blocks, and \(t_{49}\) intervals for the mean
difference.
```

**Proposed**

```latex
that surface's blocks, and approximate \(t_{49}\) intervals for the mean
difference.
```


**Current**

```latex
\caption{Each design's mean squared error (MSE) relative to simple random sampling (SRS) at $\beta = 1$ under queen contiguity, within each spillover regime. Each point is one plus the mean paired MSE difference over the 50 configuration-by-surface units divided by SRS's mean MSE, with its $t_{49}$ 95\% interval; 1 marks equal error. Filled points are significant after Holm's adjustment over the eight comparisons with SRS (Wilcoxon signed-rank test). The axis is logarithmic, ratios above 2 are printed, and rows are sorted by the control-only ratio.}
```

**Proposed**

```latex
\caption{Each design's mean squared error (MSE) relative to simple random sampling (SRS) at $\beta = 1$ under queen contiguity, within each spillover regime. Each point is one plus the mean paired MSE difference over the 50 configuration-by-surface units divided by SRS's mean MSE, with an approximate 95\% interval obtained by transforming the paired $t_{49}$ interval for the mean difference with the SRS mean held fixed. The value 1 marks equal error. Filled points are significant after Holm's adjustment over the eight comparisons with SRS (Wilcoxon signed-rank test). The axis is logarithmic, ratios above 2 are printed, and rows are sorted by the control-only ratio.}
```


### SR-14. Study averages in the downside comparison

**Accept.** The tail table averages by study and regime. Saved refinement rows include individual settings where 2×2 Blocking is below SRS, and the 3.81 ratio comes from the grid risk study.


**Current**

```latex
Table \ref{tab:downside} compares every design with SRS on both criteria, mean conditional MSE and the mean of the worst 10\% of sampled allocations, on the grid and at the refined NC settings. Under control-only spillover, both saturation designs had lower mean and worst-decile MSE than SRS in both settings. Isolation Buffer also had a lower worst decile than SRS under control-only spillover (0.62 of SRS on the grid and 0.91 in NC), but under both-arms spillover its mean was 3.80 and 3.69 times that of SRS, and it treated fewer clusters. High Incidence Focus has a single allocation in NC. That allocation's MSE was below the SRS worst-decile mean under both regimes (0.76 under both arms and 0.64 under control-only spillover), with mean MSE close to SRS (0.93 and 1.03). On the grid, however, its both-arms mean was 3.81 times that of SRS, and its NC result depends on the incidence model (see the matched baseline sensitivity below). Balanced Quartiles and Balanced Halves matched SRS on both criteria. Checkerboard and $2\times2$ Blocking did not improve on SRS on either criterion in any setting.
```

**Proposed**

```latex
Table \ref{tab:downside} compares every design with SRS on both criteria, mean conditional MSE and the mean of the worst 10\% of sampled allocations, on the grid and at the refined North Carolina settings. Under control-only spillover, both saturation designs had lower mean and worst-decile MSE than SRS in both settings. Isolation Buffer also had a lower worst decile than SRS under control-only spillover (0.62 of SRS on the grid and 0.91 in the case study), but under both-arms spillover its mean was 3.80 and 3.69 times that of SRS, and it treated fewer clusters. High Incidence Focus has a single allocation in the case study. That allocation's MSE was below the SRS worst-decile mean under both regimes (0.76 under both arms and 0.64 under control-only spillover), with mean MSE close to SRS (0.93 and 1.03). In the grid allocation-risk study, its both-arms mean MSE was 3.81 times that of SRS, and its case study result depended on the incidence model (see the matched incidence sensitivity below). Balanced Quartiles and Balanced Halves matched SRS on both criteria. For Checkerboard and $2\times2$ Blocking, mean and worst-decile MSE were at least as large as those for SRS in each study-by-regime average.
```


### SR-15. Incidence shares use different inputs on the two maps

**Accept.** Builds on HW-23. The grid calculation sums the generated incidence covariate. The case study sums observed rates and also reports death shares.


**Current**

```latex
Gain in the share of incidence assigned to treatment under Incidence-Guided Saturation compared with Saturation Quadrants, by spread in regional mean incidence ranks.
```

**Proposed**

```latex
Gain in treated incidence share under Incidence-Guided Saturation compared with Saturation Quadrants, by spread in regional mean incidence ranks. The share uses sums of simulated incidence values \(\mathbf{X}\) on the grid and sums of observed college rates in the case study.
```


### SR-16. What the case study replication counts count

**Accept.** Builds on HW-10 and HW-21. Replications share assignments, and the full-participation runs used a separate schedule for deterministic assignments.


**Current**

```latex
With 58 colleges, each setting used 100 assignments and 100
outcomes as in the main analysis and reproduced its MSE and coverage
within Monte Carlo error.

We used the same subsets, assignments and noise at all effect sizes to
keep the simulated power curves smooth in $\beta$.
```

**Proposed**

```latex
With 58 colleges, randomized rules used 100 assignments and 100
outcomes per assignment. A rule with one possible assignment used
1,000 outcomes for that assignment. The resulting MSE and coverage
matched the main analysis within Monte Carlo error.

We used the same subsets, assignments and noise across effect sizes
to reduce Monte Carlo noise in the comparisons.
```


**Current**

```latex
The main analysis used 7,598,000 independent outcome fits.
```

**Proposed**

```latex
The main analysis used 7,598,000 distinct simulated outcomes.
```


**Current**

```latex
The tail refinement added 2,280,000 independent outcomes, giving
```

**Proposed**

```latex
The tail refinement added 2,280,000 distinct outcomes, giving
```


### SR-17. Bayesian chain diagnostics

**Accept.** Builds on HW-11 and HW-29. The effective sample size and Geweke checks directly concern the treatment coefficient. Say what the diagnostics measure and identify the maximum likelihood comparison set.


**Current**

```latex
where fewer than half the chains converged,
```

**Proposed**

```latex
where fewer than half the chains met the diagnostics for $\beta$,
```


**Current**

```latex
A chain was counted as converged when the effective sample size
```

**Proposed**

```latex
A chain met the diagnostics for $\beta$ when its effective sample size
```


**Current**

```latex
83.1\% of the Bayesian CAR chains at $\rho = 0.5$ converged, against 46.2\% for the shorter chains,
```

**Proposed**

```latex
83.1\% of the Bayesian CAR chains at $\rho=0.5$ met the diagnostics for $\beta$, compared with 46.2\% of the shorter chains,
```


**Current**

```latex
Among converged chains, the Bayesian CAR estimate
```

**Proposed**

```latex
Among chains that met those diagnostics, the Bayesian CAR estimate
```


**Current**

```latex
For maximum likelihood on the same data,
```

**Proposed**

```latex
Across all maximum likelihood fits to the same datasets,
```


**Current**

```latex
Of its chains, 94.9\% converged.
```

**Proposed**

```latex
Of its chains, 94.9\% met the diagnostics for $\beta$.
```


**Current**

```latex
over converged chains
```

**Proposed**

```latex
over chains meeting the diagnostics
```


### SR-18. Display precision for SRS ratios

**Approve display changes before applying.** Builds on HW-18 and touches newly approved Section 13. The two-decimal display follows the author’s rule. The two versions of rook halves have distinct rounded values, ordered below from unguided to guided.


**Current**

```latex
Incidence-Guided Saturation's ratio to SRS was 0.631 under
control-only spillover and 1.083 under both-arms spillover.
```

**Proposed**

```latex
Incidence-Guided Saturation's ratio to SRS was 0.63 under
control-only spillover and 1.08 under both-arms spillover.
```


**Current**

```latex
Its halves at 100/0 had slightly higher MSE than SRS under control-only spillover and about 17.3 times that of SRS under both-arms spillover.
```

**Proposed**

```latex
Both versions of halves at 100/0 had slightly higher MSE than SRS under control-only spillover. Under both-arms spillover, their MSEs were 17.25 and 17.26 times that of SRS for the unguided and guided versions, respectively.
```


## Read-only companion checks

- The article’s 3.74 both-arms MSE ratio for Isolation Buffer agrees with
  the ratio of the unrounded annual means in
  `application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv`.
  Taking a ratio of the already rounded SI table entries gives a different
  rounded value. The estimates require no revision.
- The article reports coverage for the six detailed designs and SRS. The
  case study results and calibration passages in the SI include Saturation
  Quadrants and Balanced Halves as well. HW-22 and HW-28 name the SI set.
- The article’s grid methods give the lowest-to-highest interpretation of
  the incidence coefficient and use the weights symbol without the queen
  subscript in their description of incidence generation. Review both
  when the article is next open for exact-wording approval.
- The comparisons for 70/30 quadrants and 80/20 ten-region partitions
  agree with the saved Holm results. Keep that observation in the SI.
- The data-generating and fitting assumptions here concern the structural
  education outcome. Any change to the planned endpoint, participant
  population or geographic exposure proxy belongs to trial planning.

## Author decisions

Pending. Record accepted, declined or revised IDs here before applying edits.
