# Phase 2g sign-off for the SI prose

Author: Andrew Walther. Drafted 2026-11-25.

The proposals refer to the curated SI after Phases 4e and 3c. Each pair gives the
LaTeX source to replace and its proposed replacement. Where the source selects
text by build, the replacement keeps the necessary branch and the references it
uses. A blank proposed block deletes the quoted passage. Approve by ID. The
paired statistical and trial-planning proposals are in [the independent review
sign-off](2g-review-signoff.md).


## Batch 1 · Grid and case study methods


### HW-1. Neighbor definitions and Poisson incidence

**Accept.** Split the definitions and put the base rate beside the rates that inform it.


**Current**

```latex
We simulated a \(10\times10\) regular grid of 100 clusters under two
definitions of adjacency. Under queen contiguity, two clusters are
neighbors if they share an edge or a corner; under rook contiguity, they
must share an edge. Queen contiguity is the primary analysis. Rook contiguity is a sensitivity analysis that also keeps
continuity with earlier work, whose grids used rook neighbors. For
```

**Proposed**

```latex
We simulated a \(10\times10\) regular grid of 100 clusters. Under queen
contiguity, two clusters are neighbors if they share an edge or a corner.
Under rook contiguity, they must share an edge. Queen contiguity is
primary. We use rook contiguity, as in our earlier study, for sensitivity
analysis. For
```


**Current**

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

**Proposed**

```latex
their average rank. Each Poisson surface still has only 60 to 76 distinct
values among its 100 clusters. We set the base rate at 35 per 100,000,
near the average of the sudden unexpected death (SUD) rates that Mirzaei
and colleagues estimated for US men (49.3 per 100,000) and women
(21.7 per 100,000) aged 20--64 \citep{mirzaei_years_2019}. The observed
rates in the case study come from a different source, age range, period,
and case definition. The incidence
```


### HW-2. Define SRS and shorten the Checkerboard history

**Accept.** Define the reference at first use and keep the historical grid dimensions together.


**Current**

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

**Proposed**

```latex
Table \ref{tab:designs} describes the six retained designs.
\iffullsi Figure \ref{fig:design-samples}\else
Table~\ref{M-tab:allocation-rules} of the article\fi\ shows one assignment
from each design and one from simple random sampling (SRS). The examples
use a separate spatially autocorrelated incidence surface with
$\rho_X = 0.5$, generated with a fixed seed.

Checkerboard is the strict Block Stratified Sampling allocation from our
earlier small-grid study. No two rook neighbors share an arm. The earlier
design allowed relaxed allocations on its \(3\times3\) grid because the arms
could not be balanced. On its \(2\times4\) and \(3\times4\) grids, as on
the \(10\times10\) grid here, the strict allocations are the two
checkerboard patterns. We call this allocation Checkerboard hereafter.

Checkerboard is deterministic. Its assignment is fixed, and every fit
gets its own outcome noise. We drew every other design independently
25 times per incidence surface. Designs that rank clusters by incidence
break ties at random for each draw. High Incidence Focus varies only
through tie-breaking. Under the independent and spatial mechanisms,
incidence values are distinct, so the assignment is fixed within each
surface and varies across surfaces. Isolation Buffer uses the neighbor
type of the scenario and treats an average of 38 clusters with rook
contiguity and 21 with queen contiguity. Every other design treats exactly 50.
```


**Current**

```latex
and simple random sampling (SRS) at \(\beta = 1\), neither could be
```

**Proposed**

```latex
and SRS at \(\beta = 1\), neither could be
```


### HW-3. Oracle model terminology and the identification pointer

**Accept.** Use the spillover notation directly and leave the Checkerboard derivation in the diagnostics.


**Current**

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

**Proposed**

```latex
We fit two maximum likelihood (ML) spatial lag models with intercepts to
each simulated outcome. The oracle model includes the spatial lag
\(\rho \mathbf{W}\mathbf{Y}\) and the covariates \(\mathbf{Z}\),
\(\mathbf{S}(\mathbf{Z})\) and \(\mathbf{X}_k\). It matches the outcome
model in Equation \ref{eq:outcome}, using the true regime-specific
spillover term in Equation \ref{eq:spillover}. Knowing that term gives a
best-case analysis for each design. The non-oracle model includes
\(\mathbf{Z}\) and \(\mathbf{X}_k\) and omits the spillover term.
```


**Current**

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

**Proposed**

```latex
We retain Checkerboard in the design comparison and flag aliased estimates.
Section \ref{identification-diagnostics} explains why rook Checkerboard
cannot separate treatment from spillover and why queen Checkerboard
estimates the treatment effect imprecisely under control-only spillover.
It also describes the rare Isolation Buffer aliases.
```


### HW-4. Grid simulation size and shared random draws

**Accept.** Keep the run totals and the dependence between comparisons in two short paragraphs.


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
scenarios and 3.6 million fits per model and preserved the eight-design
results.

We seed each draw by hashing its identifiers. A surface uses its
configuration and index. A set of 25 assignments uses its block,
surface and design. The outcome noise uses its block. This makes the
draws reproducible in any run order. Within each block, all designs share
the same surfaces and noise, the five values of \(\beta\) share the same
assignments and noise, and both models use the same outcomes. Assignments
and noise are independent across blocks. Blocks in the same incidence
configuration share its ten surfaces.
```


### HW-5. Restriction on the death-certificate data

**Accept.** Keep the access statement in the overview and retain the details about the observed rates here.


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


### HW-6. Spatial blocking and integer treatment budgets

**Accept.** Connect the case study rule to the article and show fixed whole counts as integers.


**Current**

```latex
On the NC map, $2\times2$ Blocking uses 14 irregular groups of nearby
clusters rather than square blocks. The groups are formed once per year
by \(k\)-means on standardized projected coordinates, incidence ranks and
log(1+population), then held fixed, and contain one to eight clusters.
Half of each group, rounded to the nearest even whole number (1.5 to 2, 0.5 to 0), is treated at random, so a one-cluster group is never treated. Each year has one
single-cluster group, which stays control because half of one cluster
rounds to zero, so the yearly totals are 28, 31, 27 and 29.
```

**Proposed**

```latex
On the North Carolina map, $2\times2$ Blocking is the spatial-blocking
adaptation described in Table~\ref{M-tab:allocation-rules} of the article.
The rule uses 14 groups of nearby clusters. Each year, \(k\)-means
on standardized projected coordinates, incidence ranks and
$\log(1+\text{population})$ forms groups of one to eight clusters. We
hold the groups fixed within each year. The rule randomly treats half
of each group, rounded to the nearest even integer in a tie (1.5 to 2
and 0.5 to 0). Each year has one single-cluster group that stays in
control, and the yearly totals are 28, 31, 27 and 29.
```


**Current**

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

**Proposed**

```latex
Checkerboard & 27 / 50.4 & 27 / 50.4 & 27 / 50.4 & 27 / 50.7 \\
High Incidence Focus & 29 / 22.2 & 29 / 23.5 & 29 / 23.3 & 29 / 22.4 \\
Saturation Quadrants & 29.03 / 49.2 & 29.03 / 49.2 & 29.03 / 49.1 & 29.03 / 49.1 \\
Isolation Buffer & 15.89 / 26.7 & 15.89 / 26.7 & 15.89 / 26.7 & 15.89 / 26.7 \\
$2\times2$ Blocking & 28 / 50.0 & 31 / 55.1 & 27 / 48.0 & 29 / 50.3 \\
Balanced Quartiles & 29 / 48.5 & 29 / 49.7 & 29 / 50.6 & 29 / 49.8 \\
Balanced Halves & 29 / 50.0 & 29 / 50.5 & 29 / 49.2 & 29 / 49.4 \\
Incidence-Guided Saturation & 29 / 53.2 & 29 / 45.1 & 29 / 47.9 & 29 / 51.8 \\
Simple Random Sampling & 29 / 49.5 & 29 / 49.5 & 29 / 49.5 & 29 / 49.5 \\
```


### HW-7. Case study interval definition

**Accept.** Refer to the common interval formula and keep the case study test and cluster count here.


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
\(\mathbf{Z}\) and the true regime-specific \(\mathbf{S}(\mathbf{Z})\).
The target is the direct coefficient \(\beta\) before spatial
propagation. The case study and the partial-participation analyses use
the corrected $t$ interval in Section \ref{estimation}, with $N$ equal
to the number of participating colleges and $N=58$ in the main case
study. We test $\beta=0$ two-sided at level 0.05 by checking whether
the 95\% interval excludes zero. For the corrected interval, $p$ counts the regression coefficients
plus $\rho$. In the Wald interval the maximum likelihood variance
divides by $N$, and the normal reference
distribution ignores estimation of the parameters. Both make the
Wald interval too narrow with 58 or fewer colleges.
```


### HW-8. Simulation counts and shared case study runs

**Accept.** State what was simulated and explain shared years once, alongside the precision targets.


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
(432 settings). This is a hypothetical incidence association and
adjustment. The rook sensitivity uses \(\rho=0,0.5\) (288 settings).
Two alternative regional summaries add 96 Incidence-Guided Saturation
settings. Settings with the same possible assignments and outcome
model share a simulation, giving 936 distinct simulations for
1,248 reported settings.

We initially sampled $J=100$ assignments and simulated $R=100$
independent outcomes per distinct assignment. A design with one possible
assignment received $R=1,000$ outcomes for that assignment. The
prespecified schedule extended the number of assignments or outcomes
until the Monte Carlo standard error was at most 5\% of mean MSE and
at most 0.01 for coverage. Two single-assignment settings in the matched
sensitivity analysis required $R=2,000$. All randomized designs met
both targets at $J=100$ and $R=100$. Each random draw used a fixed seed.

In the primary model, Checkerboard, Saturation Quadrants, Isolation
Buffer and SRS share simulations across years because annual incidence
affects neither their assignments nor their outcomes. We recalculate
treatment counts and population shares for each year. Repeated values
come from the same simulated trials.
```


### HW-9. Case study reproducibility and tail overlap

**Accept.** Define the Jaccard index at first use and keep the full split-half method in one place.


**Current**

```latex
The primary case-study numbers come from the aggregate performance and allocation-risk summaries of the main run and of the separate tail refinement. Annual tables average the six queen parameter pairs equally, without pooled Monte Carlo intervals. Each random draw used a fixed, reproducible seed. An independent check recomputed performance from the stored fits, checked annual treatment budgets and regional orderings, and checked the allocations and outcomes reused from the main run in the tail refinement. The county death-certificate source data are restricted and not publicly available; only 58-cluster summaries are reported.
```

**Proposed**

```latex
The annual results give equal weight to the six queen parameter pairs
and are reported without pooled Monte Carlo intervals. An independent
check recomputed performance from the simulated trial fits, verified
annual treatment budgets and regional orderings, and checked which
assignments and outcomes the refinement reused.
```


**Current**

```latex
Estimated risk variance is reported raw and with its outcome Monte Carlo
noise contribution subtracted. Signed corrected estimates are retained
without truncation at zero. None is negative in the NC refinement.
Quantiles and worst-decile means use noisy estimated conditional risks.
To check for selection bias, worst-decile assignments were also identified with one half of the simulated outcomes and evaluated with the other half, and the roles of the two halves were then reversed. The cross-selected reductions relative to SRS were similar to or slightly larger than the reported ones. The Jaccard overlap of the two worst deciles describes sampled tail membership, not a confidence interval or population tail probability. For a design with only one possible assignment, an overlap of one is uninformative.
```

**Proposed**

```latex
We report the variance of the estimated allocation risks before and after
subtracting the contribution of outcome Monte Carlo noise. The correction
retains signed estimates. All estimates were nonnegative in the North
Carolina refinement. The estimated risks used for quantiles and
worst-decile means still contain outcome noise.

We checked for selection bias by selecting the worst decile with one half
of the simulated outcomes, evaluating it with the other half, then
exchanging the roles of the halves. The resulting reductions relative to
SRS were similar to or slightly larger than the reductions reported from
all outcomes. The Jaccard index measures agreement between the two
sampled worst deciles. It divides the number of assignments in both
deciles by the number in either. With only one possible assignment,
the index equals one by construction.
```


### HW-10. Partial participation methods

**Accept.** Give recruitment, exposure and simulation settings their own paragraphs.


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


## Batch 2 · Grid results and comparisons


### HW-11. Sensitivity methods prose

**Accept.** Use fit as the verb and name the assignment designs and SRS explicitly.


**Current**

```latex
We compared the nine strategies with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid, treating the 24 combinations of year, $\rho$ and $\psi$ as blocks, separately for each spillover regime and each value of $\tau$. The four years share one map and the strategies whose allocations ignore annual incidence share their simulations across years, so these $p$-values summarize consistency rather than independent evidence.
```

**Proposed**

```latex
We compared the eight designs and SRS with the Friedman, Nemenyi and Holm-adjusted Wilcoxon tests used on the grid, treating the 24 combinations of year, $\rho$ and $\psi$ as blocks, separately for each spillover regime and each value of $\tau$. These $p$-values summarize how consistently the designs ranked over the settings. The four years share one map and reuse some simulations (Section \ref{outcome-model-and-monte-carlo-evaluation}).
```


**Current**

```latex
To check that the design comparison does not depend on the spatial model, we generated data at $N = 58$ with $\tau = 0$ from the spatial lag model and fitted them by a conditional autoregressive (CAR) model \citep{besag_spatial_1974, cressie_statistics_1993} as well as by the spatial lag model. In the reverse direction, we generated data from a CAR model with the same mean and fitted them by both models. We also fitted two Bayesian models to the spatial lag data, a CAR model with the prior of Leroux and colleagues \citep{leroux_estimation_2000, lee_carbayes_2013} and a Bayesian spatial lag model \citep{lesage_introduction_2009}. The paragraphs below give the four models, their priors, the Markov chain Monte Carlo settings and the convergence checks. Finally, we replaced the proportional spillover exposure by a binary indicator of any treated neighbor and by a count of treated neighbors scaled by the mean number of neighbors, with the same form in the outcome model and the fit.
```

**Proposed**

```latex
To check that the design comparison does not depend on the spatial model, we generated data at $N = 58$ with $\tau = 0$ from the spatial lag model and fit a conditional autoregressive (CAR) model \citep{besag_spatial_1974, cressie_statistics_1993} and the spatial lag model to those data. In the reverse direction, we generated data from a CAR model with the same mean and fit both models. We also fit two Bayesian models to the spatial lag data, a CAR model with the prior of Leroux and colleagues \citep{leroux_estimation_2000, lee_carbayes_2013} and a Bayesian spatial lag model \citep{lesage_introduction_2009}. Finally, we replaced the proportional spillover exposure by a binary indicator of any treated neighbor and by a count of treated neighbors scaled by the mean number of neighbors, with the same form in the outcome model and the fit.
```


**Current**

```latex
\textbf{Estimation models.} Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fitted by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes and the nine strategies, with 500 fits per setting.
```

**Proposed**

```latex
\textbf{Estimation models.} Every model in the sensitivity analysis has the mean $\boldsymbol{\mu} = \alpha \mathbf{1} + \beta \mathbf{Z} + \psi\,\mathbf{S}(\mathbf{Z})$, where $\mathbf{S}(\mathbf{Z})$ is the proportional spillover covariate of the main analysis, and uses the 58-college queen geography. The spatial lag model of the main analysis is $\mathbf{Y} = \rho \mathbf{W} \mathbf{Y} + \boldsymbol{\mu} + \boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$. The CAR model is $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{u}$, $\mathbf{u} \sim N\{\mathbf{0}, \sigma^2 (\mathbf{I} - \lambda \mathbf{W}_{\mathrm{sym}})^{-1}\}$, where $\mathbf{W}_{\mathrm{sym}}$ is the symmetric matrix similar to the row-standardized $\mathbf{W}$ \citep{ord_estimation_1975}, fit by maximum likelihood with the \texttt{spautolm} function of the R package \texttt{spatialreg} \citep{bivand_r_2022}. For the reverse comparison, data were generated as $\mathbf{Y} = \boldsymbol{\mu} + \mathbf{L}\boldsymbol{\varepsilon}$ with $\mathbf{L}\mathbf{L}^\top = (\mathbf{I} - \rho \mathbf{W}_{\mathrm{sym}})^{-1}$ and the same values of $\rho$ as the spatial lag data. These runs used four years, $\rho \in \{0, 0.2, 0.5\}$, $\psi \in \{0.5, 0.8\}$, both regimes, the eight designs and SRS, with 500 fits per setting.
```


**Current**

```latex
where $a_{jk}$ is the binary queen adjacency. We used the defaults of the \texttt{S.CARleroux} function of the R package \texttt{CARBayes} \citep{lee_carbayes_2013}, version 6.1.1: independent $N(0, 10^5)$ priors for the regression coefficients, inverse-gamma(1, 0.01) priors for $\nu^2$ and $\sigma_\phi^2$, and a uniform(0, 1) prior for $\lambda$. Each chain ran 25,000 iterations, of which the first 5,000 were discarded, thinned by 5. At $\rho = 0.5$, where fewer than half the chains converged, we reran the same data with 100,000 iterations, 20,000 discarded and thinning by 20. The Bayesian spatial lag model \citep{lesage_introduction_2009} has the same form as the maximum likelihood model, with the defaults of the \texttt{spBreg\_lag} function of \texttt{spatialreg}, version 1.4.3: a $N(\mathbf{0}, 10^{12} \mathbf{I})$ prior for the regression coefficients, the prior $p(\sigma^2) \propto 1/\sigma^2$, a Beta(1.01, 1.01) prior for $\rho$ rescaled to its feasible interval, and 10,000 draws of which the first 2,000 were discarded. Both Bayesian models used four years, $\rho \in \{0, 0.5\}$, $\psi = 0.8$, both regimes and the nine strategies, with 100 fits per setting (14,400 fits per model), and both were fitted to the same data as the maximum likelihood model. We report the posterior mean of $\beta$ and the 95\% equal-tailed credible interval. A chain was counted as converged when the effective sample size of $\beta$ was at least 400 and the absolute Geweke $z$ statistic \citep{geweke_evaluating_1992} was below 1.96.
```

**Proposed**

```latex
where $a_{jk}$ is the binary queen adjacency. We used the defaults of the \texttt{S.CARleroux} function of the R package \texttt{CARBayes} \citep{lee_carbayes_2013}, version 6.1.1: independent $N(0, 10^5)$ priors for the regression coefficients, inverse-gamma(1, 0.01) priors for $\nu^2$ and $\sigma_\phi^2$, and a uniform(0, 1) prior for $\lambda$. Each chain ran 25,000 iterations, of which the first 5,000 were discarded, thinned by 5. At $\rho = 0.5$, where fewer than half the chains converged, we reran the same data with 100,000 iterations, 20,000 discarded and thinning by 20. The Bayesian spatial lag model \citep{lesage_introduction_2009} has the same form as the maximum likelihood model, with the defaults of the \texttt{spBreg\_lag} function of \texttt{spatialreg}, version 1.4.3: a $N(\mathbf{0}, 10^{12} \mathbf{I})$ prior for the regression coefficients, the prior $p(\sigma^2) \propto 1/\sigma^2$, a Beta(1.01, 1.01) prior for $\rho$ rescaled to its feasible interval, and 10,000 draws of which the first 2,000 were discarded. Both Bayesian models used four years, $\rho \in \{0, 0.5\}$, $\psi = 0.8$, both regimes, the eight designs and SRS, with 100 fits per setting (14,400 fits per model), and both used the same data as the maximum likelihood model. We report the posterior mean of $\beta$ and the 95\% equal-tailed credible interval. A chain was counted as converged when the effective sample size of $\beta$ was at least 400 and the absolute Geweke $z$ statistic \citep{geweke_evaluating_1992} was below 1.96.
```


### HW-12. Queen grid rank comparisons

**Accept.** Tie the different test outcomes to how the tests summarize the two regimes.


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
ranks did not. Each surface-level mean pools the two regimes, so
\(2\times2\) Blocking's control-only losses outweigh its both-arms wins.
```


### HW-13. Shared draws and uncertainty bands in the effect-size sweep

**Accept.** Point to the shared-draw methods and explain the bands once in the caption.


**Current**

```latex
Isolation Buffer ahead only at \(\beta = 3\). The five values of
\(\beta\) share their assignments and noise, so these five comparisons
are not independent confirmations of one another. They show only that, under
queen contiguity, the order did not change across the five values of
\(\beta\) we tested, given the same draws. In Figure \ref{fig:tau-mse}, the shaded bands describe the precision of a single scenario's MSE, not of the plotted average, whose standard error is smaller (Table \ref{tab:performance}).
```

**Proposed**

```latex
Isolation Buffer ahead only at \(\beta = 3\). The comparisons hold the
assignments and noise fixed across the five effects, as described in
Section \ref{simulation-design-resampling-and-metrics}.
```


**Current**

```latex
\caption{Mean squared error by design across $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$ under (a) queen and (b) rook contiguity, pooled over the five incidence configurations. Each point averages the 160 scenarios of one design at one value of $\beta$. Shaded bands show $\pm 1$ times the average scenario-level Monte Carlo standard error of MSE. The legend follows queen-contiguity MSE at $\beta = 1$, vertical scales differ, and the gray dashed line marks SRS. Under rook contiguity, $\beta$ is not identifiable for Checkerboard. Its estimates are kept and flagged.}
```

**Proposed**

```latex
\caption{Mean squared error by design across $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$ under (a) queen and (b) rook contiguity, pooled over the five incidence configurations. Each point averages the 160 scenarios of one design at one value of $\beta$. Shaded bands extend one average scenario-level Monte Carlo standard error of MSE above and below each point. The Monte Carlo standard error of each plotted average is smaller (Table \ref{tab:performance}). The legend follows queen-contiguity MSE at $\beta = 1$, vertical scales differ, and the gray dashed line marks SRS. For rook Checkerboard, the plotted coefficient estimates $\beta-\psi$ (Section \ref{identification-diagnostics}).}
```


### HW-14. Treatment-effect sweep results by regime

**Accept.** Split the long paragraph at the shift from pooled error to control-only results.


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
The gaps between designs changed as the treatment effect increased. Under queen contiguity, the
average rank of Incidence-Guided Saturation rose from 1.41 at
\(\beta = 0.8\) to 1.85 at \(\beta = 3\), while its MSE stayed below that of Balanced Quartiles in
103, 102, 102, 102, and 101 of the 160 blocks at
\(\beta = 0.8\), 1, 1.5, 2, and 3, and the average rank of Balanced
Quartiles stayed between 2.13 and 2.14. Isolation Buffer gained instead,
and its average rank fell from 3.71 to 3.31. The MSE of Isolation Buffer
stayed at 0.160 at every value of \(\beta\), while that of
Incidence-Guided Saturation rose from 0.088 to 0.115 and that
of Balanced Quartiles from 0.130 to 0.152, so the gap between
Incidence-Guided Saturation and Balanced Quartiles stayed
between 0.037 and 0.041.

Isolation Buffer gained under the
control-only regime, where Incidence-Guided Saturation had a
lower MSE than Isolation Buffer in 75, 73, 56, 40, and 12 of the 80
blocks at the five values of \(\beta\). The mean control-only MSE of
Incidence-Guided Saturation rose from 0.131 at \(\beta = 0.8\)
to 0.182 at \(\beta = 3\), while that of Isolation Buffer stayed between
0.159 and 0.160 (Section \ref{app:technical-notes}). At \(\beta = 2\) the
two both rounded to 0.160. At
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

**Accept.** Lead with the rerun and identify which estimates use it.


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
reproduced all 9,600 scenario summaries. It supplies the results below
for \(\hat\rho\), its associated terms in Equation
\ref{eq:error-decomposition}, and the standard errors of \(\hat\beta\).
The other numbers come from the original simulations.
```


### HW-16. Rook Checkerboard reporting

**Accept.** Keep the aliasing policy in the identification diagnostics and give readers a local pointer.


**Current**

```latex
regimes use separate noise draws. Under queen contiguity, by contrast,
Checkerboard's coverage stayed between 0.949 and 0.952 at every level of \(\rho\); under rook contiguity it was between 0.146 and 0.171\iffullsi\ (Table
\ref{tab:checkerboard-coverage})\fi. Rook Checkerboard estimates were
precise but centered on the wrong quantity, the opposite of its
imprecise but valid queen estimates. We keep these estimates in every
ranking and flag them.
```

**Proposed**

```latex
regimes use separate noise draws. Checkerboard's coverage by level of
\(\rho\) was 0.949 to 0.952 under queen contiguity and 0.146 to 0.171
under rook contiguity\iffullsi\ (Table \ref{tab:checkerboard-coverage})\fi.
Under rook contiguity the estimates were precise for the combined
coefficient \(\beta-\psi\) (Section \ref{identification-diagnostics}).
Under queen contiguity \(\beta\) was identifiable but imprecisely
estimated under control-only spillover.
```


**Current**

```latex
The table keeps Checkerboard's rook-contiguity estimates and flags them.
```

**Proposed**

```latex
The rook Checkerboard flag refers to Section \ref{identification-diagnostics}.
```


### HW-17. Robustness across parameter blocks

**Accept.** Keep the result and point to the saturation and allocation-risk comparisons for interpretation.


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
Under control-only spillover, Incidence-Guided Saturation had the lowest
mean and worst-case block MSE of the six designs\iffullsi\ (Table
\ref{tab:robustness})\fi. Under both-arms spillover, SRS and Balanced
Quartiles had slightly lower mean MSE (Section
\ref{srs-benchmark-and-allocation-specific-risk}). Section
\ref{app:variants} compares the guided and unguided saturation designs.
Worst-case block MSE describes variation across parameter settings.
Section \ref{allocation-risk-results} examines variation across
assignments within a setting.
```


### HW-18. Benchmark comparison and allocation downside

**Accept.** Distinguish the grid benchmark from the allocation-risk study and qualify the both-arms tail estimates.


**Current**

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

**Proposed**

```latex
At \(\beta=1\) under queen contiguity, pooled over the five incidence
configurations and all settings of \(\rho\), \(\psi\) and spillover regime, Incidence-Guided Saturation has
mean MSE 0.091 against 0.128 for SRS, and Saturation Quadrants has 0.090. \iffullsi Table
\ref{tab:grid-srs} separates the regimes and includes the eight designs and SRS.\else Table~\ref{M-tab:grid-srs} of the article separates the regimes for the six designs and SRS, and Figure \ref{fig:vs-srs} adds the two variants.\fi\ Incidence-Guided Saturation's ratio to SRS was 0.631 under
control-only spillover and 1.083 under both-arms spillover.
SRS, Balanced Quartiles and Balanced Halves had similar both-arms MSE.
The pooled results combine regimes that favor different designs. The comparisons do not establish equivalence.
```


**Current**

```latex
The mean is only part of the comparison, because a trial realizes a single allocation. Under control-only spillover, the worst 10\% of SRS allocations on the grid had mean conditional MSE 60\% above that of the average SRS allocation (0.346 against 0.217), and the saturation designs' worst-decile MSE was 0.54 to 0.55 times SRS's on the grid and 0.64 to 0.68 times in NC. Under both-arms spillover, their worst-decile MSE was 1.06 to 1.07 times SRS's on the grid and 0.93 to 0.94 times in NC (Table \ref{tab:downside}), so the small mean disadvantage there comes with no worse downside in NC.
```

**Proposed**

```latex
The mean is only part of the comparison, because a trial realizes a single allocation. Under control-only spillover, the worst 10\% of SRS allocations on the grid had mean conditional MSE 60\% above that of the average SRS allocation (0.346 against 0.217), and the saturation designs' worst-decile MSE was 0.54 to 0.55 times SRS's on the grid and 0.64 to 0.68 times in the North Carolina case study. Under both-arms spillover, their worst-decile MSE was 1.06 to 1.07 times SRS's on the grid and 0.93 to 0.94 times in the North Carolina case study (Table \ref{tab:downside}). Section \ref{allocation-downside-and-sensitivity} gives the precision of the case study tail estimates.
```


### HW-19. Grid allocation-risk results

**Accept.** Keep the grid numbers and the distinction between the grid study and the case study.


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


### HW-20. Case study headings and the downside table caption

**Accept.** Use plain headings, keep the section labels, and put the simulation settings in sentences.


**Current**

```latex
\section{Case study: accuracy, inference and reach}\label{completed-results-accuracy-inference-and-reach}
```

**Proposed**

```latex
\section{Accuracy, inference and reach in the case study}\label{completed-results-accuracy-inference-and-reach}
```


**Current**

```latex
\section{Case study: allocation downside}\label{allocation-downside-and-sensitivity}
```

**Proposed**

```latex
\section{Allocation downside in the case study}\label{allocation-downside-and-sensitivity}
```


**Current**

```latex
\section{Case study: participation, power and sensitivity}\label{app:nc-ext}
```

**Proposed**

```latex
\section{Participation, power and sensitivity in the case study}\label{app:nc-ext}
```


**Current**

```latex
\caption{Mean conditional MSE and worst-decile mean conditional MSE of each design divided by those of SRS in the same setting and spillover regime (below 1 is lower error than SRS). Grid: allocation-risk study, queen contiguity, $\beta=1$, averaged over 200 units, each a fixed incidence surface with one extreme combination of $\rho$ and $\psi$. NC: refined queen settings $\rho\in\{0,0.5\}$, $\psi=0.8$, averaged equally over 2018--2021 and both values of $\rho$, with the map adaptations of the designs. A design with one allocation has worst-decile mean equal to its mean. Worst-decile means are estimated from all outcomes of each assignment, without split-half selection.}
```

**Proposed**

```latex
\caption{Mean and worst-decile mean conditional MSE divided by the corresponding SRS value in the same setting and spillover regime. Ratios below 1 favor the design over SRS. The queen grid study fixes $\beta=1$ and averages over 200 units, each an incidence surface at one extreme combination of $\rho$ and $\psi$. The North Carolina analysis uses the map adaptations and averages the refined queen settings equally over 2018--2021 and both values of $\rho\in\{0,0.5\}$, at $\psi=0.8$. For a design with one assignment, the worst-decile mean equals its mean. Each worst-decile estimate uses all simulated outcomes for the assignment.}
```


## Batch 3 · Case study results


### HW-21. Accuracy, precision and records of fits

**Accept.** Keep the diagnostics next to the fit counts and drop the later account of file storage.


**Current**

```latex
The main run generated 7,598,000 independent outcome fits. Every setting is complete and meets its prespecified precision targets for mean MSE and coverage. The largest relative Monte Carlo standard error of mean MSE is 0.049 and the largest Monte Carlo standard error of coverage is 0.00892. There are no warnings, aliased columns, estimates of $\rho$ at the boundary of its parameter space or estimation failures. The tail refinement adds 2,280,000 independent outcomes, and counting the outcomes reused from the main run only once gives 9,878,000 distinct outcomes across the two stages.
```

**Proposed**

```latex
The main analysis used 7,598,000 independent outcome fits. Every setting
met its precision targets for mean MSE and coverage. The largest relative
Monte Carlo standard error of MSE was 0.049, and the largest Monte Carlo
standard error of coverage was 0.00892. We recorded no warnings,
aliased columns, boundary estimates of $\rho$ or estimation failures.
The tail refinement added 2,280,000 independent outcomes, giving
9,878,000 distinct outcomes across the two stages.
```


**Current**

```latex
Table \ref{tab:nc-annual}\iffullsi\ and Figure \ref{fig:nc-mse-detail} give\else\ gives\fi\ MSE, bias and coverage for every design in every year, and Table \ref{tab:nc-budgets} gives the budgets. More granular setting-specific results,
precision components, residual treatment information and warning and failure records are kept as aggregate output tables. These records were checked, so the absence of warnings and failures is a verified result.
```

**Proposed**

```latex
Table \ref{tab:nc-annual}\iffullsi\ and Figure \ref{fig:nc-mse-detail} give\else\ gives\fi\ MSE, bias and coverage for every design in every year, and Table \ref{tab:nc-budgets} gives the budgets.
```


### HW-22. Saturation accuracy and interval coverage

**Accept.** Keep the case study findings here and point to the comparison that explains incidence guidance.


**Current**

```latex
\ref{tab:nc-year}\else Table \ref{tab:nc-annual}\fi). Both saturation designs have lower mean MSE than SRS in every primary parameter setting, and approximate Monte Carlo intervals for each setting-specific difference excluded zero. Neither is uniquely best. The two performed similarly, so incidence guidance lets a trial direct its highest treatment saturations to the highest-incidence regions with little or no loss of estimation accuracy, and their advantage over SRS comes from the saturation structure both designs share. Under control-only spillover, the accuracy of High Incidence Focus
depended on the planning year.
```

**Proposed**

```latex
\ref{tab:nc-year}\else Table \ref{tab:nc-annual}\fi). Both saturation designs have lower mean MSE than SRS in every primary parameter setting, and approximate Monte Carlo intervals for each setting-specific difference excluded zero. Both had similar accuracy. Section \ref{app:variants} describes how their shared saturation structure and different regional ordering affect accuracy and targeting. Under control-only spillover, the accuracy of High Incidence Focus
depended on the planning year.
```


**Current**

```latex
Mean coverage of the 95\% intervals across the four years was between 94.1\% and 94.9\% for every design and regime, and between 94.1\% and 95.1\% in single years. The Wald interval covered only 92.6\% to 93.6\% (Table \ref{tab:interval-comparison}). Meeting the Monte Carlo precision targets did not repair that finite-sample undercoverage, but the corrected $t$ interval described in Section \ref{outcome-model-and-monte-carlo-evaluation} did (Section \ref{interval-calibration}). Table \ref{tab:nc-annual} gives bias and coverage for every design, year and regime. Checkerboard's treatment effect is identifiable on this irregular graph, but weak separation of treatment from spillover makes its estimates imprecise. This differs from the exact non-identifiability under rook contiguity on the rectangular grid.
```

**Proposed**

```latex
Mean coverage of the 95\% intervals across the four years was between 94.1\% and 94.9\% across the eight designs and SRS in both regimes, and between 94.1\% and 95.1\% in single years. The Wald interval covered only 92.6\% to 93.6\% (Table \ref{tab:interval-comparison}). Section \ref{interval-calibration} compares the corrected $t$ interval with the Wald interval and other alternatives. Table \ref{tab:nc-annual} gives bias and coverage for every design, year and regime. Checkerboard's treatment effect is identifiable on this irregular graph, but weak separation of treatment from spillover makes its estimates imprecise. Section \ref{identification-diagnostics} gives the corresponding grid results.
```


### HW-23. Geographic reach and spread of regional incidence

**Accept.** Keep the case study figures in the case study and use the spread figure caption for the grid comparisons.


**Current**

```latex
Incidence-Guided Saturation sends more treatment to high-incidence regions when incidence differs between regions. On the NC map it differs little, because the four regions' mean incidence ranks lie between 24 and 34 of 58 in every year. With all 58 colleges, Incidence-Guided Saturation's treated colleges held on average 51.1\% of the summed college incidence rates, against 50.4\% for Saturation Quadrants, 50.0\% for SRS and 62.2\% for High Incidence Focus. Shares of deaths run the other way (49.2\%, 49.6\%, 49.7\% and 35.4\%), because the highest-rate colleges serve small populations. On the grid, where regional incidence differed more and equal cluster populations make rate share and death share coincide, incidence guidance raised the treated share of incidence over Saturation Quadrants by 2.1 to 3.1 percentage points on average across the five incidence configurations, while High Incidence Focus reached 0.75 to 0.77 by giving up randomization over the whole map. Across surfaces the gain rose almost linearly with the spread of the regions' mean incidence ranks (correlation 0.84 to 0.99), and the NC regions lie at the low end of that spread (Figure \ref{fig:spread-gain}).
```

**Proposed**

```latex
Regional mean incidence ranks in North Carolina ranged from 24 to 34
out of 58 in every year. With all 58 colleges, the treated colleges held
51.1\% of the sum of college incidence rates under Incidence-Guided
Saturation, 50.4\% under Saturation Quadrants, 50.0\% under SRS and
62.2\% under High Incidence Focus. The corresponding shares of deaths
were 49.2\%, 49.6\%, 49.7\% and 35.4\%. Colleges with the highest
rates serve small populations. Figure \ref{fig:spread-gain} compares
the spread of regional ranks on the two maps. Section
\ref{igs-sensitivity-results} shows how changing regions and saturation
targets affects targeting and MSE on the grid.
```


**Current**

```latex
\caption{Gain in treated incidence share from incidence guidance (Incidence-Guided Saturation minus Saturation Quadrants) against the spread of the regions' mean incidence ranks, the largest minus the smallest regional mean rank divided by the number of clusters. Filled symbols are the 50 grid surfaces (queen, production allocations) by incidence configuration, with a least-squares line fitted to them. Open triangles are the fixed NC regions in each of the four planning years, whose spread divides by 58 clusters rather than 100.}
```

**Proposed**

```latex
\caption{Gain in the share of incidence assigned to treatment under Incidence-Guided Saturation compared with Saturation Quadrants, by spread in regional mean incidence ranks. Spread is the highest minus the lowest regional mean rank divided by the number of clusters. Filled symbols show the 50 grid surfaces under queen contiguity by incidence configuration, with a least-squares line. Incidence guidance added 2.1 to 3.1 percentage points on average across the five configurations, and the correlation between gain and spread ranged from 0.84 to 0.99. High Incidence Focus reached 0.75 to 0.77 by assigning treatment to the highest-incidence half of the grid, with randomness only in ties. Open triangles show North Carolina's four planning years, at the low end of the spread. Spread divides by 58 for the case study and by 100 for the grid.}
```


### HW-24. Control-only downside and the article figure

**Accept.** Keep the tail findings here and refer to the methods and the article caption for how the estimates were made.


**Current**

```latex
At the two refined queen settings, $\rho=0$ and $\rho=0.5$ with $\psi=0.8$, under control-only spillover, Saturation Quadrants' estimated 90th percentile of allocation-specific MSE is 31.8--32.9\%
below SRS and its estimated worst-decile mean 34.6--37.8\% below. Incidence-Guided
Saturation's corresponding reductions are 20.8--38.6\% and 20.7--43.4\%.
The pattern persists when the worst decile is identified and evaluated on separate halves of the simulated outcomes (Section \ref{app:nc}). \iffullsi Figure \ref{fig:nc-risk} shows the $\rho=0.5$ setting. \fi The tail estimates rest on a finite number of outcomes per allocation, so worst-decile membership is uncertain. They describe the spread of conditional MSE across the sampled allocations, not quantiles of realized errors, and because the refinement reused the main run's allocations they are not a fresh allocation sample. Read this way, the results show a descriptive downside benefit
alongside lower average error, conditional on the studied model and
sampled assignments.
```

**Proposed**

```latex
At the two refined queen settings, $\rho=0$ and $\rho=0.5$ with
$\psi=0.8$, the control-only 90th percentile of assignment-specific MSE
was 31.8--32.9\% below SRS for Saturation Quadrants, and the worst-decile
mean was 34.6--37.8\% below. The corresponding reductions for
Incidence-Guided Saturation were 20.8--38.6\% and 20.7--43.4\%.
The split-half check gave a similar pattern (Section \ref{app:nc}).
\iffullsi Figure \ref{fig:nc-risk} shows the $\rho=0.5$ setting. \fi
Under the study model, the saturation designs had lower worst-decile
error across the sampled assignments and lower average MSE.
```


**Current**

```latex
\iffullsi Figure \ref{fig:alloc-dist}\else Figure~\ref{M-fig:alloc-dist} of the article\fi\ shows the full allocation-to-allocation distribution for the six designs and SRS. For each allocation in \iffullsi Figure \ref{fig:alloc-dist}\else that figure\fi, MSE is estimated from repeated simulated outcomes (100 in the simulation study, at least 400 in the case study). Percentiles are computed within each setting and then averaged, so the bars show allocation-to-allocation spread, and outcome Monte Carlo noise in the per-allocation estimates widens them slightly. Under control-only spillover, most Incidence-Guided Saturation allocations had lower MSE than the average SRS allocation in both studies, whereas the SRS distribution extended well above its mean.
```

**Proposed**

```latex
\iffullsi Figure \ref{fig:alloc-dist}\else Figure~\ref{M-fig:alloc-dist} of the article\fi\
summarizes assignment-specific MSE for the six designs and SRS. Each
estimate uses 100 outcomes per grid assignment and at least 400 per
case study assignment. Under control-only spillover, most
Incidence-Guided Saturation assignments had lower MSE than the
average SRS assignment in both studies. The SRS distribution extended
well above its mean.
```


### HW-25. Tail membership and Monte Carlo error

**Accept.** Use the methods definition of the tail overlap and keep the numerical limits on interpretation here.


**Current**

```latex
The tail estimates remain uncertain. With \(J=100\), the worst decile
contains only ten sampled draws, and extra outcomes reduce noise in
\(\hat m(a)\) but do not enlarge the allocation sample. All 112 reported settings in the tail refinement for designs with more than one possible assignment flag uncertain worst-decile membership. Under both arms, the overlap between the worst deciles identified by the two outcome halves (the Jaccard index: allocations in both divided by allocations in either) was 0--0.25 for Incidence-Guided Saturation, and approximately 61--86\% of the raw variance of estimated risk across allocations was attributable to outcome Monte Carlo noise. Saturation Quadrants' worst-decile membership was similarly unstable. Small both-arms tail differences therefore do not establish superior downside protection. Checkerboard and High Incidence Focus each have only one possible assignment for these inputs, so their 90th percentile and worst-decile mean equal their conditional mean MSE. That equality does not imply low error or absence of outcome uncertainty.
```

**Proposed**

```latex
With \(J=100\), each worst decile contains ten sampled draws. The
refinement added outcomes to the same assignments, reducing noise in
\(\hat m(a)\) (Section \ref{outcome-model-and-monte-carlo-evaluation}).
Worst-decile membership was flagged as uncertain in all 112 reported
refinement settings for designs with more than one possible assignment.
Under both-arms spillover, Incidence-Guided Saturation had a Jaccard
index of 0--0.25 (Section \ref{app:nc}). Outcome Monte Carlo noise
accounted for approximately 61--86\% of its raw estimated variance in
assignment risk. Saturation Quadrants had similarly unstable worst-decile
membership. The both-arms tail differences are too small to establish
better downside protection.

Checkerboard and High Incidence Focus each have one assignment for
these inputs. Their 90th percentile and worst-decile mean therefore
equal their conditional mean MSE, which still depends on the
simulated outcomes.
```


### HW-26. Incidence sensitivity row label

**Accept.** Name the covariate and explain that it enters both simulation and estimation.


**Current**

```latex
\caption{Incidence-Guided Saturation sensitivity MSE divided by matched primary MSE, by year and regime. Baseline adds rank-scaled $\mathbf{X}$ to both generation and fit; rook uses the matched extreme combinations of $\rho$ and $\psi$. Regional-summary alternatives preserve fixed regions.}
```

**Proposed**

```latex
\caption{Incidence-Guided Saturation sensitivity MSE divided by matched primary MSE, by year and regime. The incidence sensitivity adds rank-scaled $\mathbf{X}$ to the simulated outcome and the fitted model. Rook contiguity uses the matched extreme combinations of $\rho$ and $\psi$. Regional-summary alternatives preserve fixed regions.}
```


**Current**

```latex
Baseline / adjusted & Both & 1.016 & 1.030 & 1.015 & 1.035 \\
Baseline / adjusted & Control & 1.026 & 1.023 & 1.042 & 1.040 \\
```

**Proposed**

```latex
Incidence in outcome and model & Both & 1.016 & 1.030 & 1.015 & 1.035 \\
Incidence in outcome and model & Control & 1.026 & 1.023 & 1.042 & 1.040 \\
```


### HW-27. Power results opening and use of MDE

**Accept.** Lead with the conditions for the results and use the abbreviation already defined in the methods if HW-10 is approved.


**Current**

```latex
This section reports the participation, power and sensitivity analyses of the case study. All numbers come from aggregate summaries of runs that wrote to their own output folders, each with a manifest of parameters, code hashes and package versions. These runs did not change the main case-study results. Unless stated otherwise, results use random volunteering, $\tau = 0$, $\kappa = 0.8$ and the corrected $t$ test, and they average the four planning years.
```

**Proposed**

```latex
Unless stated otherwise, these results use random volunteering,
$\tau=0$, $\kappa=0.8$ and the corrected $t$ test, averaged over the
four planning years.
```


**Current**

```latex
\caption{Colleges needed for 80\% power and smallest effect detectable with 80\% power (MDE, college-level standard deviations), corrected $t$ test, random volunteering, $\rho = 0.5$, spillover $\psi = 0.8\beta$, mean over the 2018--2021 planning years. $>58$: 80\% power not reached with all 58 colleges in at least one year. $>2$: 80\% power not reached at the largest simulated effect, $\beta = 2$.}
```

**Proposed**

```latex
\caption{Colleges needed for 80\% power and minimum detectable effect (MDE) in college-level standard deviations, corrected $t$ test, random volunteering, $\rho = 0.5$, spillover $\psi = 0.8\beta$, mean over the 2018--2021 planning years. An entry of $>58$ means estimated power stayed below 80\% with all 58 colleges in at least one year. An entry of $>2$ means estimated power stayed below 80\% at the largest simulated effect, $\beta=2$.}
```


### HW-28. Interval calibration

**Accept.** Name the nine-comparator scope, distinguish the two runs and describe the shortfall directly.


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

**Accept.** Separate maximum likelihood CAR, Bayesian CAR and Bayesian spatial lag while retaining the reported diagnostics.


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


### HW-30. Case study design terminology and block dependence

**Accept.** Use designs and SRS consistently and refer to the methods for the dependence between years.


**Current**

```latex
When the number of trained professionals varies, with $n_c$ proportional to college population and averaging 25, and ICC $= 0.05$, power fell by 0.097 on average across designs, settings and effect sizes. For $d = 0.3$ under both-arms spillover, the mean power over the nine strategies fell from 0.59 to 0.46 with 20 colleges and from 0.88 to 0.77 with 58. MSE rose by a factor of 1.24 to 1.84 across designs and settings, and no design gained relative to the others.
```

**Proposed**

```latex
When the number of trained professionals varies, with $n_c$ proportional to college population and averaging 25, and ICC $= 0.05$, power fell by 0.097 on average across designs, settings and effect sizes. For $d = 0.3$ under both-arms spillover, the mean power over the eight designs and SRS fell from 0.59 to 0.46 with 20 colleges and from 0.88 to 0.77 with 58. MSE rose by a factor of 1.24 to 1.84 across designs and settings, and no design gained relative to the others.
```


**Current**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives, for each strategy, the number of colleges needed for 80\% power and the smallest college-level effect detectable with 80\% power at 29 and 58 colleges. The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not perfectly monotone in $N$ within Monte Carlo error, because the participant subsets differ between values of $N$, and we report the dips rather than smooth them.
```

**Proposed**

```latex
\textbf{Colleges needed and smallest detectable effect.} Table \ref{tab:nc-power} gives, for each design and for SRS, the number of colleges needed for 80\% power and the smallest college-level effect detectable with 80\% power at 29 and 58 colleges. The number of colleges needed is interpolated monotonically between the simulated values of $N$ within each year, averaged over years and rounded up to a whole college. Results at $\rho = 0$ differed from those at $\rho = 0.5$ by about two colleges or fewer. Power curves are not perfectly monotone in $N$ within Monte Carlo error, because the participant subsets differ between values of $N$, and we report the dips rather than smooth them.
```


**Current**

```latex
Table \ref{tab:nc-guidance} gives each strategy's mean MSE relative to SRS at $\beta = 1$ for 20, 29, 40 and 58 participating colleges. The regime-specific guidance held with fewer participants, but the size of the saturation advantage depended on participation. Checkerboard and Isolation Buffer remained worse than SRS under both-arms spillover at every number of colleges, and the rank-based designs remained within a few percent of SRS. Under control-only spillover, the saturation designs' mean MSE was 22\% below that of SRS with 58 colleges, 9\% to 10\% below with 40, 4\% to 5\% below with 29 and about 1\% below with 20. With few participants, non-participating neighbors already create the exposure variation that saturation otherwise supplies. Isolation Buffer went from matching SRS under control-only spillover with 58 colleges to 2.3 times SRS's MSE with 29 colleges and 6.2 times with 15.
```

**Proposed**

```latex
Table \ref{tab:nc-guidance} reports mean MSE for each design relative to SRS at $\beta = 1$ for 20, 29, 40 and 58 participating colleges. The regime-specific guidance held with fewer participants, but the size of the saturation advantage depended on participation. Checkerboard and Isolation Buffer remained worse than SRS under both-arms spillover at every number of colleges, and the rank-based designs remained within a few percent of SRS. Under control-only spillover, the saturation designs' mean MSE was 22\% below that of SRS with 58 colleges, 9\% to 10\% below with 40, 4\% to 5\% below with 29 and about 1\% below with 20. With few participants, non-participating neighbors already create the exposure variation that saturation otherwise supplies. Isolation Buffer went from matching SRS under control-only spillover with 58 colleges to 2.3 times SRS's MSE with 29 colleges and 6.2 times with 15.
```


**Current**

```latex
\textbf{Formal tests on the NC results.} Friedman tests over the 24 combinations of year, $\rho$ and $\psi$ rejected equal ranks of the nine strategies in all four combinations of regime and $\tau$ ($\chi^2 = 135.3$, 153.8, 150.7 and 179.8 on 8 degrees of freedom, $p \le 2.3\times10^{-25}$; Nemenyi critical difference 2.45). Under control-only spillover, Saturation Quadrants and Incidence-Guided Saturation had the two best average ranks (1.67 and 1.79 with $\tau = 0$; 1.29 and 1.71 with $\tau = 1$), and SRS ranked seventh with $\tau = 0$ and sixth with $\tau = 1$. Under both-arms spillover with $\tau = 0$, High Incidence Focus ranked first (1.88) and SRS sixth (5.17), with the balancing and saturation designs between them. These tests treat the 24 blocks as independent, although the years share one map and some strategies share their simulations across years, so they describe the consistency of the ordering rather than independent evidence.
```

**Proposed**

```latex
\textbf{Formal tests on the NC results.} Friedman tests over the 24 combinations of year, $\rho$ and $\psi$ rejected equal ranks of the eight designs and SRS in all four combinations of regime and $\tau$ ($\chi^2 = 135.3$, 153.8, 150.7 and 179.8 on 8 degrees of freedom, $p \le 2.3\times10^{-25}$; Nemenyi critical difference 2.45). Under control-only spillover, Saturation Quadrants and Incidence-Guided Saturation had the two best average ranks (1.67 and 1.79 with $\tau = 0$; 1.29 and 1.71 with $\tau = 1$), and SRS ranked seventh with $\tau = 0$ and sixth with $\tau = 1$. Under both-arms spillover with $\tau = 0$, High Incidence Focus ranked first (1.88) and SRS sixth (5.17), with the balancing and saturation designs between them. The four years share one map and reuse some simulations, as described in Section \ref{further-sensitivity-analyses}. Interpret the ordering descriptively.
```


## 2g changes for Chapter 3

These entries extend the Phase 2i hand-over after C3-1 to C3-6 and the
4e/3c sign-off. Carry across only approved wording where the same material
appears. Use the chapter's own labels and table layout.

| Area | Proposals to carry forward if approved |
|---|---|
| Grid methods | HW-1 to HW-4. Align incidence and neighbor definitions, the SRS introduction, the treatment and spillover notation, and the explanation of shared draws. Keep the Checkerboard alias derivation in the appendix diagnostics. |
| Case study methods | HW-5 to HW-11. Consolidate data access, participation and precision details. Use the spatial-blocking name, integer fixed budgets and the common interval definition. Define the Jaccard index and MDE on first use. |
| Grid results | HW-12 to HW-19. Match comparisons, the uncertainty-band explanation, the rerun provenance and the location of the allocation-risk caveats. |
| Case study results and exhibits | HW-20 to HW-30. Align headings where useful, reach versus incidence targeting, tail interpretation and diagnostics, the incidence sensitivity row, and calibration and model sensitivity. Retain the chapter's own figure and section numbering. |
| Independent review | Carry across the approved statistical and trial-planning corrections listed in [the independent review sign-off](2g-review-signoff.md), including corrections to the new sections before transferring them. |

The curated overview keeps “the full methods.” The two “Completed primary
case study” captions flagged in the older review now belong solely to
unabridged figures. The curated build cites the article's exhibits instead.

## Author decisions

Pending. Record accepted, declined or revised IDs here before applying edits.
