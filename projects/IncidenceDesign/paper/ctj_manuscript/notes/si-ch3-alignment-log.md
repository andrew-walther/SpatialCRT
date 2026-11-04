# SI edits with counterparts in Chapter 3 or the CTJ manuscript (2026-11-04)

The SI was reorganized, its captions trimmed, and its prose revised (commits 582562b to c347f42).
Chapter 3 and the CTJ manuscript were deliberately left unchanged in that pass. Each entry below is an SI
caption or passage whose old text also appears in the file named (at least half of its 12-word windows
verbatim, or for short passages the first 60 characters), so the two now differ. Use this list in the
later SI/Chapter 3/CTJ alignment pass. Blocks moved between SI sections without rewording are not
listed (see si-reorg-map.md). Entry headings use the section numbers
of the reorganized SI with an S prefix (S7.1 = Section 7.1); the prefix was later dropped from
section numbers, which are otherwise unchanged.

42 entries (25 captions, 17 prose passages).

## 1. Caption tab:designs
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

The six treatment assignment designs compared in the detailed analysis, ordered by mean squared error under queen contiguity at $\tau = 1$, pooled over the five incidence configurations (Table \ref{tab:performance}). A cluster's incidence rank is its rank by the baseline incidence surface $X_k$, with ties broken at random in every draw. Clusters treated is the number of clusters treated per assignment; for Isolation Buffer, which depends on the neighbor type, it is the mean over its assignments under each neighbor type.

**New SI text**

The six treatment assignment designs compared in the detailed analysis, ordered by mean squared error under queen contiguity at $\tau = 1$, pooled over the five incidence configurations (Table \ref{tab:performance}). Incidence rank is a cluster's rank by the baseline incidence surface $X_k$, with ties broken at random in every draw. Clusters treated is per assignment; for Isolation Buffer it is the mean under each neighbor type.

## 2. Caption fig:design-samples
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

One example assignment from each of the six designs and SRS on the $10\times10$ grid, with the designs in the order of Table \ref{tab:designs} followed by SRS. Every example uses the illustrative incidence surface in the first panel (darker is higher; dashed lines mark the quadrants), a spatially autocorrelated surface with $\rho_X = 0.5$ generated with a fixed seed; it is not one of the simulation's incidence surfaces. Blue clusters are treated and white clusters are controls. Lines mark the $2\times2$ blocks and the quadrants, numbers in the Balanced Quartiles panel give each cluster's incidence quartile, and percentages in the Incidence-Guided Saturation panel give each quadrant's treated share. Isolation Buffer is shown under queen contiguity.

**New SI text**

One example assignment from each of the six designs and SRS on the $10\times10$ grid, in the order of Table \ref{tab:designs}, all on the illustrative incidence surface in the first panel (darker is higher). Blue clusters are treated and white clusters are controls. Lines mark the $2\times2$ blocks and the quadrants, numbers in the Balanced Quartiles panel give incidence quartiles, and percentages in the Incidence-Guided Saturation panel give treated shares. Isolation Buffer is shown under queen contiguity.

## 3. Caption tab:excluded-designs
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

The two related variants outside the detailed six-design grid ranking, in the format of Table \ref{tab:designs}. Saturation Quadrants shares its mechanism with Incidence-Guided Saturation, and Balanced Halves with Balanced Quartiles (see Treatment assignment designs). A cluster's incidence rank is its rank by the baseline incidence surface $X_k$, with ties broken at random in every draw. Clusters treated is the number of clusters treated per assignment.

**New SI text**

The two related variants outside the detailed six-design grid ranking, in the format of Table \ref{tab:designs}. Saturation Quadrants shares its mechanism with Incidence-Guided Saturation, and Balanced Halves with Balanced Quartiles. Incidence rank and clusters treated are defined as in Table \ref{tab:designs}.

## 4. Caption fig:service-areas
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

The assignment of North Carolina's 100 counties to the 58 community college service areas used as clusters in the planned application. Each county is colored by its primary serving college, and each service area carries a single number, keyed to the legend. Hatching marks the two counties that also lie within Roanoke-Chowan Community College's broader service network: Northampton County, assigned to Halifax Community College (24), and Bertie County, assigned to Martin Community College (30). Roanoke-Chowan Community College (41) serves Hertford County alone in this analysis. The map shows county boundaries and the county-to-college assignment only; it does not depend on the SUD data.

**New SI text**

Assignment of North Carolina's 100 counties to the 58 community college service areas used as clusters in the planned application. Each county is colored by its primary serving college, and each service area carries a number keyed to the legend. Hatching marks Northampton County, assigned to Halifax Community College (24), and Bertie County, assigned to Martin Community College (30), which also lie within the broader service network of Roanoke-Chowan Community College (41).

## 5. Caption tab:performance
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Performance of the six designs at the primary scenario, $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages the 160 scenarios of one design and neighbor type (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes), each summarizing 250 fits of the oracle model. Rows in both panels are in order of mean squared error (MSE) under queen contiguity. Monte Carlo standard errors of MSE, in parentheses, are computed from the 50 configuration-by-surface means (Section \ref{app:mechanics}); those of bias, coverage, and power are all below 0.007. Average rank is the mean of each design's within-block MSE rank among the six designs and SRS over the 160 blocks (1 = lowest MSE). SRS's average rank is 2.600 under queen contiguity and 2.450 under rook contiguity. Clusters treated is the mean number of clusters treated per assignment. Under rook contiguity, $\tau$ is not identified for Checkerboard ($WZ = 1 - Z$); its estimates are kept and flagged, and its coefficient on $Z$ estimates $\tau - \gamma$.

**New SI text**

Performance of the six designs at $\tau = 1$, pooled over the five incidence configurations, under queen and rook contiguity. Each entry averages 160 scenarios of 250 oracle fits each. Rows follow queen-contiguity mean squared error (MSE); SE is its Monte Carlo standard error. Average rank is the within-block MSE rank among the six designs and SRS (1 = lowest); clusters treated is a mean per assignment. Under rook contiguity Checkerboard is not identified ($WZ = 1 - Z$), and its coefficient on $Z$ estimates $\tau - \gamma$.

## 6. Caption fig:mse
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Mean squared error (MSE) of the oracle estimate of $\tau$ by design at $\tau = 1$, under queen contiguity. Columns of panels correspond to the four levels of outcome spatial dependence $\rho$, rows to the both-arms (top) and control-only (bottom) spillover regimes, and bar colors to the spillover magnitude $\gamma$. Each bar overlays the MSEs of the five incidence configurations for one combination of design, $\rho$, $\gamma$, and regime: the bar's height is the largest of the five, and the horizontal lines within it mark the other four. Designs are ordered by pooled mean MSE (Table \ref{tab:performance}). The gray dashed line in each panel marks simple random sampling (SRS), the largest of its five configuration MSEs averaged over the four values of $\gamma$. The vertical scale is set separately for each row. Figure \ref{fig:mse-rook} in Section \ref{app:rook-figures} shows rook contiguity.

**New SI text**

Oracle mean squared error (MSE) by design at $\tau = 1$ under queen contiguity, with designs in the order of Table \ref{tab:performance}. Columns are the levels of $\rho$, rows the both-arms (top) and control-only (bottom) regimes, and colors the values of $\gamma$. Each bar's height is the largest MSE of the five incidence configurations, and lines within it mark the other four. The gray dashed line marks SRS, the largest of its five configuration MSEs averaged over $\gamma$. Vertical scales differ by row.

## 7. Caption fig:coverage
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Coverage of the 95\% confidence intervals by design at $\tau = 1$, under queen contiguity. Each box summarizes the coverage of the 160 scenarios of one design (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes), each estimated from 250 fits. Designs are ordered by median coverage, so the order differs from the MSE order used elsewhere. The dashed red line marks the nominal 0.95 level, and the gray dotted line the median coverage of simple random sampling. Figure \ref{fig:coverage-rook} in Section \ref{app:rook-figures} shows rook contiguity.

**New SI text**

Coverage of the 95\% confidence intervals by design at $\tau = 1$ under queen contiguity. Each box summarizes the 160 scenarios of one design, each estimated from 250 fits. Designs are ordered by median coverage, not by MSE. The dashed red line marks the nominal 0.95 level, and the gray dotted line the median coverage of simple random sampling.

## 8. Caption tab:by-config
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Mean squared error by design within each of the five incidence configurations at $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity). Each entry averages the 32 scenarios of one design, neighbor type, and configuration (4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes), each summarizing 250 fits of the oracle model. The spatially autocorrelated and Poisson configurations are labeled by the strength of dependence in the incidence surface, $\rho_X$. Rows in both panels are ordered by pooled mean MSE under queen contiguity (Table \ref{tab:performance}). Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

**New SI text**

Mean squared error by design within each of the five incidence configurations at $\tau = 1$, under queen and rook contiguity. Each entry averages the 32 scenarios of one design, neighbor type, and configuration (4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes), each from 250 oracle fits. Spatial and Poisson configurations are labeled by $\rho_X$. Rows follow Table \ref{tab:performance}. Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

## 9. Caption tab:gamma
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Mean squared error by design across the four levels of spillover magnitude $\gamma$ at $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages 40 scenarios (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 2 spillover regimes). Rows in both panels are ordered by pooled mean MSE under queen contiguity (Table \ref{tab:performance}). Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

**New SI text**

Mean squared error by design across the four levels of spillover magnitude $\gamma$ at $\tau = 1$, under queen and rook contiguity, pooled over the five incidence configurations. Each entry averages 40 scenarios (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 2 spillover regimes). Rows follow Table \ref{tab:performance}. Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

## 10. Caption tab:regime
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Mean squared error by design under each spillover regime at $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages 80 scenarios (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\gamma$). Rows are ordered by pooled mean MSE under queen contiguity (Table \ref{tab:performance}). Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

**New SI text**

Mean squared error by design under each spillover regime at $\tau = 1$, under queen and rook contiguity, pooled over the five incidence configurations. Each entry averages 80 scenarios (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\gamma$). Rows follow Table \ref{tab:performance}. Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

## 11. Caption fig:biasvar
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Decomposition of mean squared error into squared bias and variance by design at $\tau = 1$ for one incidence configuration, Poisson with $\rho_X = 0.20$, under queen contiguity. For each design, the two stacked components are the squared bias and the variance ($\text{SD}^2$) of $\hat\tau$, each averaged over the 32 scenarios of that configuration (4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes). Designs are ordered by pooled mean MSE (Table \ref{tab:performance}). The gray dashed line marks the MSE of simple random sampling (SRS) in the same configuration. Figure \ref{fig:biasvar-rook} in Section \ref{app:rook-figures} shows rook contiguity.

**New SI text**

Squared bias and variance ($\text{SD}^2$) of $\hat\tau$ by design at $\tau = 1$ under queen contiguity, for the Poisson configuration with $\rho_X = 0.20$. Each stacked component is averaged over the 32 scenarios of that configuration (4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes). Designs follow Table \ref{tab:performance}, and the gray dashed line marks the MSE of simple random sampling (SRS) in the same configuration.

## 12. Caption tab:eight-design
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Performance of all eight designs and SRS at the primary scenario, $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages the 160 scenarios of one design and neighbor type (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes), each summarizing 250 fits of the oracle model. Rows in both panels are in order of mean squared error (MSE) under queen contiguity. Retained marks the six designs of the detailed comparison, and SRS is the benchmark. Average rank is the mean of each design's within-block MSE rank among the eight designs and SRS over the 160 blocks (1 = lowest MSE), so it differs from the average rank in Table \ref{tab:performance}. Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

**New SI text**

Performance of all eight designs and SRS at $\tau = 1$, pooled over the five incidence configurations, under queen and rook contiguity. Each entry averages the 160 scenarios of one design and neighbor type, each from 250 oracle fits. Rows follow queen-contiguity mean squared error (MSE), and Retained marks the six designs of the detailed comparison. Average rank is the mean within-block MSE rank among the eight designs and SRS (1 = lowest). Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

## 13. Caption fig:tau-mse
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Mean squared error by design across the five values of the treatment effect, $\tau \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$, under (a) queen contiguity (primary) and (b) rook contiguity (sensitivity). Each point averages the 160 scenarios of one design at one value of $\tau$, pooled over the five incidence configurations. Shaded bands show $\pm 1$ times the average of the scenario-level Monte Carlo standard errors of MSE; they describe the precision of a single scenario's MSE, not of the plotted average, whose standard error is smaller (Table \ref{tab:performance}). The shared legend lists the designs in order of pooled mean MSE under queen contiguity at $\tau = 1$. The vertical scales of the two panels differ. The gray dashed line marks simple random sampling (SRS). Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged. The five values of $\tau$ share their assignments and noise.

**New SI text**

Mean squared error by design across $\tau \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$ under (a) queen and (b) rook contiguity, pooled over the five incidence configurations. Each point averages the 160 scenarios of one design at one value of $\tau$. Shaded bands show $\pm 1$ times the average scenario-level Monte Carlo standard error of MSE. The legend follows queen-contiguity MSE at $\tau = 1$, vertical scales differ, and the gray dashed line marks SRS. Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged.

## 14. Caption tab:oracle-nonoracle
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Oracle and non-oracle estimates of $\tau$ by design at $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages the 160 scenarios of one design and neighbor type. The oracle model has the covariates $Z$, $\text{Spill}$, and $X_k$; the non-oracle model omits $\text{Spill}$. Rows in both panels are in order of oracle MSE under queen contiguity. Under rook contiguity, $\tau$ is not identified for Checkerboard ($WZ = 1 - Z$); its estimates are kept and flagged, and the two models coincide.

**New SI text**

Oracle and non-oracle estimates of $\tau$ by design at $\tau = 1$, under queen and rook contiguity, pooled over the five incidence configurations. Each entry averages the 160 scenarios of one design and neighbor type. The non-oracle model omits $\text{Spill}$. Rows follow oracle MSE under queen contiguity. Under rook contiguity, $\tau$ is not identified for Checkerboard ($WZ = 1 - Z$); its estimates are kept and flagged, and the two models coincide.

## 15. Caption fig:vs-srs
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Each design's mean squared error (MSE) compared with simple random sampling (SRS) at $\tau = 1$ under queen contiguity, within each spillover regime. Each point is one plus the mean paired difference in MSE over the 50 independent configuration-by-surface units divided by SRS's mean MSE in that regime, so that 1 marks equal error and values below 1 mark lower error than SRS, with its $t_{49}$ 95\% interval. Filled points are significant after Holm's adjustment over the eight designs compared with SRS within the regime (Wilcoxon signed-rank test). The horizontal axis is logarithmic, ratios above 2 are printed, and rows are sorted by the control-only ratio.

**New SI text**

Each design's mean squared error (MSE) relative to simple random sampling (SRS) at $\tau = 1$ under queen contiguity, within each spillover regime. Each point is one plus the mean paired MSE difference over the 50 configuration-by-surface units divided by SRS's mean MSE, with its $t_{49}$ 95\% interval; 1 marks equal error. Filled points are significant after Holm's adjustment over the eight comparisons with SRS (Wilcoxon signed-rank test). The axis is logarithmic, ratios above 2 are printed, and rows are sorted by the control-only ratio.

## 16. Caption tab:grid-srs
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Full grid benchmark at $\tau=1$, queen contiguity. Each regime-specific entry averages 80 scenarios (five incidence configurations, four $\rho$ and four $\gamma$ values). MSE is the mean squared error of the estimate of $\tau$; coverage is the proportion of 95\% $t$-based confidence intervals, with a degrees-of-freedom-corrected standard error, that contain $\tau$. Both use the oracle fit, the spatial lag model that includes the true spillover term. Simple random sampling randomly selects 50 of the 100 clusters for treatment. No value pools queen and rook.

**New SI text**

Full grid benchmark at $\tau=1$, queen contiguity. Each regime-specific entry averages 80 scenarios (five incidence configurations, four $\rho$ and four $\gamma$ values). MSE and coverage use the oracle fit; coverage is that of 95\% $t$-based intervals with a degrees-of-freedom-corrected standard error. Simple random sampling treats 50 of the 100 clusters at random. No value pools queen and rook.

## 17. Caption fig:nc-mse
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Completed primary application: yearly mean MSE divided by year- and regime-matched SRS MSE, for the eight other designs with queen weights and education outcomes. Each cell averages six $\rho/\gamma$ settings equally before taking the ratio. Blue indicates lower error than SRS, orange higher error; colors use a logarithmic ratio scale and cells print the ratios. SRS mean MSE is 0.07638 under both arms and 0.30761 control-only. Years that share the same simulations are not independent replications. Absolute MSE, bias and coverage are in Table \ref{tab:nc-annual}, budgets in Table \ref{tab:nc-budgets} and setting-specific comparisons in Figure \ref{fig:nc-setting-ratios}.

**New SI text**

Completed primary application: yearly mean MSE divided by year- and regime-matched SRS MSE for the eight other designs, with queen weights and education outcomes. Each cell averages six $\rho/\gamma$ settings equally before taking the ratio. Blue marks lower error than SRS and orange higher; colors use a logarithmic ratio scale and cells print the ratios. Years that share simulations are not independent replications.

## 18. Caption tab:nc-annual
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Primary queen application performance by design and planning year. Each entry averages the six $\rho/\gamma$ pairs of one spillover regime with equal weight; MSE, bias and coverage of the corrected $t$ interval are Monte Carlo estimates. Checkerboard, Saturation Quadrants, Isolation Buffer and SRS use the same simulations in every year, because annual incidence enters neither their allocation nor the outcome model, so their rows repeat; repeated values are not independent annual confirmation. Treatment counts and population shares are in Table \ref{tab:nc-budgets}.

**New SI text**

Primary queen application performance by design and planning year. Each entry averages the six $\rho/\gamma$ pairs of one spillover regime with equal weight; coverage is that of the corrected $t$ interval. Checkerboard, Saturation Quadrants, Isolation Buffer and SRS use the same simulations in every year, because annual incidence enters neither their allocation nor the outcome model, so their rows repeat and are not independent annual confirmation.

## 19. Caption fig:spread-gain
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Gain in treated incidence share from incidence guidance, Incidence-Guided Saturation minus Saturation Quadrants, against the spread of the regions' mean incidence ranks, (largest minus smallest regional mean rank) divided by the number of clusters. Filled symbols are the 50 grid surfaces (queen, production allocations) by incidence configuration, with a least-squares line fitted to the grid surfaces only. Open triangles are the four fixed NC regions in each planning year, whose spread divides by 58 clusters rather than 100.

**New SI text**

Gain in treated incidence share from incidence guidance (Incidence-Guided Saturation minus Saturation Quadrants) against the spread of the regions' mean incidence ranks, the largest minus the smallest regional mean rank divided by the number of clusters. Filled symbols are the 50 grid surfaces (queen, production allocations) by incidence configuration, with a least-squares line fitted to them. Open triangles are the fixed NC regions in each of the four planning years, whose spread divides by 58 clusters rather than 100.

## 20. Caption tab:downside
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Mean conditional MSE and worst-decile mean conditional MSE of each design relative to SRS: each entry is the design's average divided by SRS's average in the same setting and spillover regime (below 1 is lower error than SRS). Grid: allocation-risk study, queen contiguity, $\tau=1$, averaged over 200 units, each a fixed incidence surface with one of the extreme combinations of $\rho$ and $\gamma$. NC: refined queen settings $\rho\in\{0,0.5\}$, $\gamma=0.8$, averaged with equal weight over 2018--2021 and both values of $\rho$; NC designs are the map adaptations (Checkerboard as a greedy graph-coloring assignment, $2\times2$ Blocking on spatial blocks of nearby clusters, and both saturation designs on four fixed regions). A design with one allocation has worst-decile mean equal to its mean.

**New SI text**

Mean conditional MSE and worst-decile mean conditional MSE of each design divided by those of SRS in the same setting and spillover regime (below 1 is lower error than SRS). Grid: allocation-risk study, queen contiguity, $\tau=1$, averaged over 200 units, each a fixed incidence surface with one extreme combination of $\rho$ and $\gamma$. NC: refined queen settings $\rho\in\{0,0.5\}$, $\gamma=0.8$, averaged equally over 2018--2021 and both values of $\rho$, with the map adaptations of the designs. A design with one allocation has worst-decile mean equal to its mean.

## 21. Caption fig:alloc-dist
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Distribution of allocation-specific MSE for the six designs and simple random sampling (SRS), relative to the mean MSE of SRS in the same setting (queen weights). Simulation study: grid allocation-risk study at $\tau=1$, 200 fixed surface/parameter settings per spillover regime. Application study: refined NC settings ($\rho=0, 0.5$; $\gamma=0.8$) in each year 2018--2021. For each allocation, MSE is estimated from repeated simulated outcomes (100 in the simulation study, at least 400 in the application study). Thin lines span the 5th to 95th percentiles of allocations and thick bars the 25th to 75th; circles mark the mean. Percentiles are computed within each setting and then averaged, so bars show allocation-to-allocation spread; a design with one allocation per setting appears as a circle only. Outcome Monte Carlo noise in the per-allocation estimates widens the bars slightly.

**New SI text**

Distribution of allocation-specific MSE for the six designs and simple random sampling (SRS), relative to the mean MSE of SRS in the same setting (queen weights). Simulation study: grid allocation-risk study at $\tau=1$, 200 fixed surface/parameter settings per spillover regime. Application study: refined NC settings ($\rho=0, 0.5$; $\gamma=0.8$) in each year 2018--2021. Thin lines span the 5th to 95th percentiles of allocations, thick bars the 25th to 75th, and circles the mean; a design with one allocation per setting appears as a circle only.

## 22. Caption fig:nc-risk
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Allocation-specific conditional MSE in the refined primary queen setting $\rho=0.5,\gamma=0.8$, by year and spillover regime. Mean MSE, the estimated 90th percentile of allocation-specific MSE and the estimated worst-decile mean summarize the same $J=100$ sampled assignments, with at least 400 independent outcomes per assignment of a randomized design or 1,000 for designs with a single possible assignment. These tail estimates, based on a finite number of outcomes, are descriptive, with uncertain membership; they are not quantiles of realized errors or a fresh allocation confirmation.

**New SI text**

Allocation-specific conditional MSE in the refined primary queen setting $\rho=0.5,\gamma=0.8$, by year and spillover regime. Mean MSE, the estimated 90th percentile of allocation-specific MSE and the estimated worst-decile mean summarize the same $J=100$ sampled assignments, with at least 400 independent outcomes per assignment of a randomized design or 1,000 for designs with a single possible assignment.

## 23. Caption fig:colleges-needed
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Number of participating colleges needed for 80\% power with the corrected $t$ test, by design, for a college-level effect of $\tau = 1$ standard deviation under both-arms spillover and $\tau = 1.5$ under control-only spillover. Random volunteering, $\rho = 0.5$ and spillover $\gamma = 0.8\tau$. Dark bars show the mean over the 2018--2021 planning years and whiskers the range. Light bars mark designs that did not reach 80\% power with all 58 colleges in the stated number of years, with the smallest number needed in the other years.

**New SI text**

Number of participating colleges needed for 80\% power with the corrected $t$ test, by design, for a college-level effect of $\tau = 1$ standard deviation under both-arms spillover and $\tau = 1.5$ under control-only spillover, with random volunteering, $\rho = 0.5$ and $\gamma = 0.8\tau$. Dark bars show the mean over the 2018--2021 planning years and whiskers the range. Light bars mark designs that did not reach 80\% power with all 58 colleges in the stated number of years, with the smallest number needed in the other years.

## 24. Caption fig:power-n
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Power of the corrected $t$ test against the number of participating colleges for simple random sampling (SRS), the two saturation designs and Checkerboard, by spillover regime (columns) and college-level effect (rows). Random volunteering, $\rho = 0.5$ and spillover $\gamma = 0.8\tau$. Lines show the mean over the 2018--2021 planning years and bands the range, which is narrower than the lines except near 58 colleges, and the dashed line marks 80\% power. The two saturation designs nearly coincide. Each number of colleges uses its own 50 random participant subsets, so the small drop from 29 to 30 colleges reflects which subsets were drawn.

**New SI text**

Power of the corrected $t$ test against the number of participating colleges for simple random sampling (SRS), the two saturation designs and Checkerboard, by spillover regime (columns) and college-level effect (rows), with random volunteering, $\rho = 0.5$ and $\gamma = 0.8\tau$. Lines show the mean over the 2018--2021 planning years and bands the range, and the dashed line marks 80\% power. Each number of colleges uses its own 50 random participant subsets.

## 25. Caption tab:interval-comparison
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Coverage of the 95\% Wald interval and the corrected $t$ interval ($t$ reference distribution with a degrees-of-freedom-corrected standard error), the interval reported throughout, for the six designs and SRS. Grid: oracle model at $\tau = 1$, pooled over the five incidence configurations and both spillover regimes (160 scenarios per entry). Application: queen weights, $\beta = 0$, averaged over the four planning years and the six $\rho/\gamma$ pairs of each regime. Both intervals were computed from the same fits. \textsuperscript{a}Under rook contiguity on the grid, $\tau$ is not identified for Checkerboard.

**New SI text**

Coverage of the 95\% Wald interval and the corrected $t$ interval ($t$ reference distribution with a degrees-of-freedom-corrected standard error) for the six designs and SRS. Grid: oracle model at $\tau = 1$, pooled over the five incidence configurations and both spillover regimes (160 scenarios per entry). Application: queen weights, $\beta = 0$, averaged over the four planning years and the six $\rho/\gamma$ pairs of each regime. \textsuperscript{a}Under rook contiguity on the grid, $\tau$ is not identified for Checkerboard.

## 26. S1.2
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

draws and the outcome noise are shared more selectively, as described
under Simulation design, resampling, and metrics.

**New SI text**

draws and the outcome noise are shared more selectively, as described
in Section \ref{simulation-design-resampling-and-metrics}.

## 27. S5.1
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

includes heart-failure deaths formerly excluded as presumed non-sudden.
It totals 23,523
deaths over 25,594,321 person-years. The earlier numerator reproduces
the totals behind Habib's unpublished report
\citep{habib_temporal_2026}; the corrected numerator is the application
input.

**New SI text**

includes heart-failure deaths formerly excluded as presumed non-sudden.
It totals 23,523
deaths over 25,594,321 person-years. The earlier numerator, which
excluded them, reproduces the totals behind Habib's unpublished report
\citep{habib_temporal_2026}. The application uses the corrected
numerator.

## 28. S5.2
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Candidates use standardized projected coordinates and are selected
without outcome performance.

**New SI text**

Candidates use standardized projected coordinates and were selected
without reference to outcome performance.

## 29. S5.3
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

annual incidence enters neither their allocation nor the outcome model. These repeated values are not independent
annual confirmation.

**New SI text**

annual incidence enters neither their allocation nor the outcome model, so their repeated values are not independent
annual confirmations.

## 30. S7.2
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

The configuration changed each design's MSE only slightly: in this
order, each design's highest MSE

**New SI text**

The configuration changed each design's MSE only slightly. In this
order, each design's highest MSE

## 31. S7.3
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Isolation Buffer matched or passed Incidence-Guided Saturation under this regime (see Stability across treatment effect
sizes).

**New SI text**

Isolation Buffer matched or passed Incidence-Guided Saturation under this regime (Section \ref{stability-across-treatment-effect-sizes}).

## 32. S7.3
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Among the six detailed designs, Incidence-Guided Saturation was first or second under both regimes, at this and every other value of $\tau$. In the full benchmark that includes SRS (Table \ref{tab:grid-srs}), its both-arms MSE was 8\% above that of SRS, a small but statistically detectable difference (Figure \ref{fig:vs-srs}). This matters because a trial team may not know in advance which regime applies.

**New SI text**

Among the six detailed designs, Incidence-Guided Saturation was first or second under both regimes, at this and every other value of $\tau$, which matters because a trial team may not know in advance which regime applies. In the full benchmark that includes SRS (Table \ref{tab:grid-srs}), however, its both-arms MSE was 8\% above that of SRS, a small but statistically detectable difference (Figure \ref{fig:vs-srs}).

## 33. S7.3
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

(0.004) and its coverage stayed at 0.951: its intervals widened with its
variance. Its pooled MSE of 1.087 was almost entirely variance: squared
bias, averaged over scenarios, was 0.0045, less than 0.5\% of MSE

**New SI text**

(0.004) and its coverage stayed at 0.951, because its intervals widened with its
variance. Its pooled MSE of 1.087 was almost entirely variance, and squared
bias, averaged over scenarios, was 0.0045, less than 0.5\% of MSE

## 34. S12.2
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

their small both-arms disadvantage is statistically detectable. Under rook contiguity

**New SI text**

their small both-arms disadvantage is statistically detectable.

Under rook contiguity

## 35. S13
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Every setting is complete and meets its prespecified precision targets for mean MSE and coverage: the largest relative

**New SI text**

Every setting is complete and meets its prespecified precision targets for mean MSE and coverage. The largest relative

## 36. S13
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

The tail refinement adds 2,280,000 independent outcomes; counting the outcomes reused from the main run only once gives 9,878,000 distinct outcomes across the two stages. Monte Carlo precision describes
uncertainty in simulated performance; it does not certify interval
calibration.

**New SI text**

The tail refinement adds 2,280,000 independent outcomes, and counting the outcomes reused from the main run only once gives 9,878,000 distinct outcomes across the two stages. Monte Carlo precision describes
uncertainty in simulated performance, not interval
calibration.

## 37. S13
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Table \ref{tab:nc-annual} retains MSE, bias and coverage for every design in every year; budgets are in Table \ref{tab:nc-budgets}. More granular setting-specific results,
precision components, residual treatment information and warning and failure records are kept as aggregate output tables. The absence of warnings and failures is a verified result, not a sign that diagnostics were skipped.

**New SI text**

Table \ref{tab:nc-annual} and Figure \ref{fig:nc-mse-detail} give MSE, bias and coverage for every design in every year, and Table \ref{tab:nc-budgets} gives the budgets. More granular setting-specific results,
precision components, residual treatment information and warning and failure records are kept as aggregate output tables. These records were checked, so the absence of warnings and failures is a verified result.

## 38. S14
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

Treatment--incidence alignment leaves less residual treatment
information after adjustment; this sensitivity changes baseline and
adjustment jointly and does not imply that SUD incidence predicts
education ability. Rook and queen comparisons matched at the same extreme combinations of $\rho$ and $\gamma$ preserve broad regime-specific conclusions.

**New SI text**

When treatment follows incidence, adjusting for incidence leaves less
residual information about treatment. This sensitivity analysis changes baseline and
adjustment together and does not imply that SUD incidence predicts
education ability. Comparing rook and queen contiguity at the same extreme combinations of $\rho$ and $\gamma$ leaves the regime-specific conclusions broadly unchanged.

## 39. S15.1
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

and about 19 colleges for 1.5 standard deviations (Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). Checkerboard and Isolation Buffer did not reach 80\% power at $\tau = 1$ even with all 58 colleges. Under control-only

**New SI text**

and about 19 colleges for 1.5 standard deviations (Figure \ref{fig:colleges-needed}, Table \ref{tab:nc-power}). Checkerboard and Isolation Buffer did not reach 80\% power at $\tau = 1$ even with all 58 colleges.

Under control-only

## 40. S15.4
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

and under a count exposure the saturation designs remained the best under control-only spillover.

**New SI text**

and under a count exposure the saturation designs still had the lowest MSE under control-only spillover.

## 41. S15.4
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

With the longer chains, 83.1\% of the Bayesian CAR chains

**New SI text**

\textbf{Bayesian models.} With the longer chains, 83.1\% of the Bayesian CAR chains

## 42. S15.5
Also in: Chapter 3 (paper/thesis_chapter/Thesis_Chapter.qmd)

**Old SI text**

and SRS alongside other randomized
strategies when both-arm exposure is expected.

**New SI text**

and SRS alongside other randomized
strategies when both-arms exposure is expected.

## Later structural cleanup (same day, after the entries above)

Three SI passages were rewritten more heavily after the human-writer pass. Chapter 3
has the old versions of each. Compare section by section rather than sentence by sentence.

- **SI Section 2 opening and 2.1 (excluded variants).** The case for leaving out
  Saturation Quadrants and Balanced Halves (MSEs, Nemenyi/Wilcoxon p-values, average
  ranks) moved from the start of Section 2 into 2.1. Section 2 now opens with a short
  pointer. The restated quadrant saturations (20/40/60/80%) were dropped from the prose
  because Table S2 gives them. Chapter 3 counterpart: "Treatment assignment designs" in
  Methods.
- **SI Section 5.2 (geographic adaptations).** Rewritten as one sequence: adjacency,
  rank and interspersion designs, saturation regions with the criterion C, blocking,
  budgets. The criterion paragraph that sat after the budgets table was folded
  in. No numbers dropped. Chapter 3 counterpart: "Geographic adaptations and treatment
  budgets".
- **SI Sections 15.1 and 15.4 (participation and sensitivity).** The unequal-professionals
  result moved from 15.4 into 15.1 and replaced the rougher "about 0.1" sentence. In
  15.4, the summary paragraph that repeated the detailed CAR/Bayesian, spillover-form
  and formal-test numbers is now a three-sentence lead-in, followed by Spatial model,
  Spillover form and Formal tests paragraphs. Unique facts from the summary (CAR ratios
  0.76–0.82, the count-form result, why the Bayesian CAR fails) were kept in those
  paragraphs. Chapter 3 counterpart: "Colleges needed under partial participation",
  "Sensitivity analyses" and appendix A9.
- **SI Sections 13, 14, 15.3 and 5.3 (repeated passages).** Section 13's two targeting
  paragraphs (incidence and death shares, grid gain of 2.1–3.1 points) are merged into
  one. The two 15.3 paragraphs that both reported the guidance ratios are now one paragraph
  for random volunteering and one for incidence-weighted volunteering. The split-half
  worst-decile sentence is kept in full only in 5.4: it is dropped from 5.3 and shortened
  to a pointer in 14. Chapter 3 counterparts: "Accuracy, coverage and reach",
  "Allocation-specific risk in North Carolina", "Colleges needed under partial
  participation" and appendix A9.
