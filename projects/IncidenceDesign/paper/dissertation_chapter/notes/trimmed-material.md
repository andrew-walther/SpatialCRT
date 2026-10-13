# Chapter 3 trimmed material

Text removed or relocated while preparing the concise Chapter 3 copy
(`Dissertation_Chapter_concise.qmd`, branch `ch3-concise`, started 2026-10-12;
plan: `docs/plans/chapter3_length_review_plan_2026-10-12.md`). Nothing here is
deleted from the canonical `Dissertation_Chapter.qmd`, which stays unchanged
until the author approves the concise copy as its replacement.

Each entry records the date, the source section, the reason, and the verbatim
removed text. Text that was moved (not cut) is noted with its new location only.

## 2026-10-12 — Batch 1 (Introduction and Methods)

### 3.1.1 The design-stage gap — light trim
Reason: the same literature is reviewed in Chapter 1; the gap argument and every
design-stage citation except the two analysis-stage statistics are kept. Dropped
citations (`moran_notes_1950`, `geary_contiguity_1954`) remain cited in Chapter 1.

Removed wording (original two paragraphs): see canonical
`Dissertation_Chapter.qmd`, section "The design-stage gap in the spatial interference
literature"; the concise copy merges the restricted/covariate-constrained randomization
sentences and lists the four recent theoretical papers in one citation group instead
of describing each.

Original paragraphs:

> The methodological literature addressing spatial interference in CRTs has concentrated on analysis-stage correction. Spatial autocorrelation statistics such as Moran's $I$ [@moran_interpretation_1948; @moran_notes_1950] and Geary's $C$ [@geary_contiguity_1954] provide tools for detecting residual spatial dependence once outcomes have been observed, and spatial regression models, including spatial lag and spatial Durbin specifications [@lesage_introduction_2009], estimated by maximum likelihood [@ord_estimation_1975] and implemented in modern software [@bivand_r_2022], allow an analyst to separate a direct treatment effect from a spillover effect after the fact. In their systematic review of spatial analysis methods for cluster randomized trials, Jarvis and colleagues placed trials that used spatial separation or buffers outside the primary review of statistical methods [@jarvis_spatial_2017]. The literature therefore offers much more guidance on analyzing spatial effects than on how the choice of treatment assignment itself, made before any data are collected, might reduce the problem.
>
> Design-stage tools do exist, but they were largely developed for a different purpose. Restricted randomization addresses the baseline imbalance that a limited number of clusters can produce [@raab_balance_2001]. Covariate-constrained randomization generates candidate assignments, retains those meeting a prespecified covariate-balance criterion, and randomizes among them [@moulton_covariatebased_2004]; Crisp and colleagues extended it to the selection of clusters from a larger eligible set and noted that their procedure also increases spatial dispersion among the selected clusters [@crisp_covariate-constrained_2023]. Their dispersion criterion concerns which clusters are selected rather than how the selected clusters are assigned to arms, and none of these methods targets the spatial pattern of arm assignment under spillover: balancing a covariate across arms does not by itself determine whether treated and control clusters are spatially interspersed in a way that induces contamination. A more recent body of work addresses the problem studied in this chapter more directly. Randomized saturation designs assign different saturations, the fractions of units treated, to different groups and randomize treatment within each group, which allows direct and spillover effects to be jointly identified rather than merely controlled for [@baird_optimal_2018]. A conventional "fried-egg" design omits each cluster's outer area from analysis to limit contamination; McCann and colleagues proposed instead selecting whole clusters for inclusion or exclusion while preserving randomization [@mccann_spatial_2018]. Most recently, a group of theoretical papers published or posted between 2022 and 2025 has derived rate-optimal cluster-randomized designs under spatial interference [@leung_rateoptimal_2022], independent-set designs under network interference [@cai_independentset_2023], designs that explicitly account for cross-cluster interference [@leung_clusterrandomized_2025], and designs for spatio-temporally heterogeneous treatment effects [@watson_design_2025]. This chapter complements this work with a simulation comparison of concrete assignment strategies that a trial team could implement, under heterogeneous baseline incidence and a range of spatial dependence and spillover conditions.

### 3.2.3 Treatment assignment designs — moved to Appendix B.1
Reason: the Nemenyi/Wilcoxon p-values, MSEs and average ranks that justify keeping
one design from each related pair were stated in Methods before the tests are defined,
and B.1 only pointed back to them. All numbers now appear in B.1. Original paragraph:

> We implemented eight treatment assignment designs and present six of them as the primary comparison set in this chapter. The other two, Saturation Quadrants and Balanced Halves, each share the mechanism of a retained design, and in the eight-design comparison at $\tau = 1$ neither could be distinguished from its counterpart under either neighbor type, by the Nemenyi and Holm-adjusted Wilcoxon signed-rank tests described under Simulation design, resampling, and metrics below, applied here to all eight designs. Saturation Quadrants assigns the saturations 20\%, 40\%, 60\%, and 80\% to the four quadrants in random order, without reference to incidence; its retained counterpart, Incidence-Guided Saturation Quadrants, assigns the higher saturations to the quadrants with higher mean incidence (Table \ref{tab:designs}). We shorten the second name to Incidence-Guided Saturation in the rest of this chapter, as in the figures. Saturation Quadrants's mean squared error was 0.0900, against 0.0907 for Incidence-Guided Saturation, under queen contiguity (Nemenyi $p = 0.98$; Holm-adjusted Wilcoxon signed-rank $p = 0.41$), and 0.0708 against 0.0717 under rook contiguity (Nemenyi $p = 1.00$; Wilcoxon $p = 0.82$). We retain one representative of the saturation mechanism, Incidence-Guided Saturation. Saturation Quadrants had the marginally lower mean under both neighbor types and the better average rank in the eight-design comparison (2.100 against 2.350 under queen contiguity and 2.281 against 2.306 under rook contiguity), so this choice does not rest on performance. The two saturation designs performed similarly, so incidence guidance lets a trial direct its highest treatment saturations to the highest-incidence quadrants with little or no loss of estimation accuracy; the accuracy advantage of both designs over SRS, which arises under control-only spillover (see SRS benchmark and allocation-specific risk), comes from the saturation structure they share. Balanced Halves divides clusters into two halves by incidence rank and treats 25 at random in each. Its mean squared error was 0.1299, against 0.1308 for Balanced Quartiles, under queen contiguity (Nemenyi $p = 1.00$; Wilcoxon $p = 0.88$), and 0.0842 against 0.0837 under rook contiguity (Nemenyi $p = 1.00$; Wilcoxon $p = 0.82$). The two designs differ only in how finely they stratify on incidence, and we retain the finer, four-stratum version. Balanced Halves had the marginally lower mean and the better average rank under queen contiguity (3.312 against 3.356), and Balanced Quartiles had both under rook contiguity (2.919 against 3.038), so this choice does not rest on performance either. Both related variants, and the full eight-design comparison, are reported in Appendix \ref{app:eight-design}.

### 3.2.5 Simulation design, resampling, and metrics — condensed
Reason: the seed-key scheme, metric formulas, interval definitions and the pooled /
paired Monte Carlo standard-error rules are given in full in Appendix B.2 (Simulation
mechanics, Tables B.3–B.4 and the formulas). The concise text keeps the fit structure,
benchmark size, common-random-number consequences, metric list and SE rule.
Original paragraphs:

> Within each scenario, we drew $J = 25$ assignments from each of the $K = 10$ incidence surfaces, for 250 fits per scenario, and gave every fit its own noise vector, so that fit $(k, j)$ has outcome $Y_{kj} = (I - \rho W)^{-1}[\tau Z_{kj} + S(Z_{kj}) + \beta X_k + \epsilon_{kj}]$. Checkerboard's 250 fits share one assignment but have 250 distinct noise vectors. Given a surface, its 25 pairs of assignment and noise are independent, and the ten surfaces are independent of one another, so the surface is the primary sampling unit. Each estimator was fit 3.2 million times. SRS was subsequently run with the same matched surfaces, noise keys and 250-fit structure, bringing the complete benchmark to 14,400 scenarios and 3.6 million fits per estimator; the original eight-design results were unchanged. SRS assigns exactly 50 of the 100 clusters to treatment by complete randomization, ignoring incidence and geography.
>
> Every random draw was seeded from a hash of a key that names the draw: the configuration and surface index for an incidence surface; the block, surface, and design for a set of 25 assignments; and the block for the noise. Each draw therefore depends only on its key, and the results are identical whether scenarios are run sequentially, in parallel, or in any order. The keys also give common random numbers. Within a block, all designs share the same noise vectors and incidence surfaces; the five values of $\tau$ reuse the same assignments and noise; and the oracle and non-oracle models are fit to the same outcomes. Assignments and noise are independent across blocks, but blocks within an incidence configuration share its ten surfaces.
>
> For each scenario and estimator, we computed five metrics over the valid fits, those with a finite estimate and standard error: bias, the mean of $\hat\tau - \tau$; the standard deviation (SD) of $\hat\tau$; mean squared error, $\text{MSE} = \text{mean}\{(\hat\tau - \tau)^2\}$; coverage, the fraction of fits with $L \le \tau \le U$, where $L, U = \hat\tau \mp 1.96\,\widehat{\text{SE}}$; and power, the fraction of fits with $L > 0$, a one-sided test of $\tau = 0$ at level 0.025. For both estimators, every one of the 250 fits in every scenario was valid, and the fitting routine raised no warning other than the aliasing flags described above.
>
> Monte Carlo standard errors come from the surface-level metrics. For a metric $m$, we computed $m_k$ on the 25 fits from surface $k$ and set $\text{SE}(m) = \text{sd}(m_1, \dots, m_{10})/\sqrt{10}$, with intervals $m \pm t_{0.975,\,9}\,\text{SE}(m)$. A difference between two designs is paired at the surface level: its standard error is that of $m_{d,k} - m_{d',k}$ over $k$. For a metric pooled over the five incidence configurations, we apply the same formula to the 50 configuration-by-surface means, each averaged over that surface's 32 blocks; these 50 means are independent of one another. Appendix \ref{app:mechanics} gives the parameter grid, the seed keys, and the exact formulas for these metrics and standard errors.

## 2026-10-12 — Batch 2 (grid results)

### 3.3.2 Formal statistical comparison — merged into 3.3.1
Reason: the average ranks repeat Table 3.2; the Nemenyi p-values for unseparated pairs move to Appendix B.3. The subsection heading is removed and the condensed test summary follows the block-level paragraph of 3.3.1.

Original text:

> ## Formal statistical comparison
>
> The Friedman test rejected equal MSE ranks across the six designs under queen contiguity at $\tau = 1$ ($\chi^2 = 451.57$, $\text{df} = 5$, 160 blocks, $p < 2.2\times10^{-16}$). The average within-block ranks were 1.43 for Incidence-Guided Saturation, 2.14 for Balanced Quartiles, 3.69 for Isolation Buffer, 4.24 for $2\times2$ Blocking, 4.74 for Checkerboard, and 4.76 for High Incidence Focus. Nemenyi comparisons (critical difference 0.596) separated 11 of the 15 pairs (Appendix \ref{app:supp-comparisons} shows the critical-difference diagrams). Incidence-Guided Saturation ranked significantly better than every other design, including Balanced Quartiles ($p = 0.0086$), and Balanced Quartiles ranked significantly better than the remaining four. The four pairs that Nemenyi comparisons did not separate all lay among the bottom four designs: $2\times2$ Blocking against Isolation Buffer ($p = 0.090$), High Incidence Focus ($p = 0.130$), and Checkerboard ($p = 0.170$), and Checkerboard against High Incidence Focus ($p = 1.00$). Holm-adjusted Wilcoxon signed-rank tests separated all 15 pairs.

### 3.3.2 rank-versus-mean and surface-level paragraphs — condensed
Reason: the reversal explanation is kept in shorter form; the paired differences and t49 intervals over the 50 units move to Appendix B.3.

Original text:

> The bottom three designs were ordered differently by average rank and by mean MSE. By mean MSE the order was High Incidence Focus (0.245), $2\times2$ Blocking (0.319), and Checkerboard (1.087); by average rank it was $2\times2$ Blocking (4.24), Checkerboard (4.74), and High Incidence Focus (4.76). Both reversals have the same source: the MSEs of $2\times2$ Blocking and Checkerboard depended strongly on the spillover regime. Each had a lower MSE than High Incidence Focus in all 80 both-arms blocks (mean MSE 0.062 for each, against 0.162) and a higher MSE in all 80 control-only blocks (0.576 and 2.111, against 0.329; see the spillover-regime results below). Average rank counts each block once, so their wins under the both-arms regime offset their losses under the control-only regime; mean MSE weights those losses by their size.
>
> Because blocks within an incidence configuration share their ten surfaces, we also compared each pair of designs adjacent in the mean-MSE order over the 50 independent configuration-by-surface units. Every adjacent pair differed significantly. The mean paired differences in MSE, with $t_{49}$ 95\% intervals, were $-0.0401$ ($-0.0431$, $-0.0372$) for Incidence-Guided Saturation against Balanced Quartiles, $-0.0293$ ($-0.0332$, $-0.0255$) for Balanced Quartiles against Isolation Buffer, $-0.0853$ ($-0.0954$, $-0.0752$) for Isolation Buffer against High Incidence Focus, $-0.0738$ ($-0.0858$, $-0.0617$) for High Incidence Focus against $2\times2$ Blocking, and $-0.7677$ ($-0.7893$, $-0.7460$) for $2\times2$ Blocking against Checkerboard, with Wilcoxon $p \le 4.9\times10^{-9}$ in every case. The better design of each pair had the lower MSE in all 50 units, except High Incidence Focus against $2\times2$ Blocking, where it did so in 48 (96\%). This test follows the mean-MSE order. It separated High Incidence Focus from $2\times2$ Blocking, which the Nemenyi comparison of within-block ranks did not, because each unit's MSE averages over both spillover regimes, and there the control-only losses of $2\times2$ Blocking outweigh its both-arms wins.

### 3.3.4 Tables 3.4 (rho) and 3.5 (gamma) — moved to new Appendix B.4
Reason: both tables show that neither factor changes the ranking, which the text states with the ranges; the tables are unchanged in the new Appendix section "Additional sensitivity tables".

### 3.3.5 Rook contiguity — condensed
Reason: the rook test statistics and within-configuration coverage ranges move to Appendix B.3; Table 3.7 (Checkerboard coverage by rho) moves to the new Appendix B.4 and its range is stated in the text. The Checkerboard tau - gamma explanation is kept in short form (it is derived in Methods and B.6.1).

Original text:

> Under rook contiguity at $\tau = 1$, pooled over the five incidence configurations, the mean-MSE order was Incidence-Guided Saturation (0.072), Balanced Quartiles (0.084), $2\times2$ Blocking (0.128), Isolation Buffer (0.133), High Incidence Focus (0.206), and Checkerboard (0.480) (Table \ref{tab:performance}; Figure \ref{fig:mse-rook} in Appendix \ref{app:rook-figures}), and the same order held within each of the five incidence configurations (Table \ref{tab:by-config}). Incidence-Guided Saturation had the lowest MSE in 97 of the 160 rook blocks (61\%) and Balanced Quartiles in 60 (38\%). Coverage was between 0.938 and 0.942 for every design except Checkerboard, whose coverage was 0.150 (Figure \ref{fig:coverage-rook} in Appendix \ref{app:rook-figures}); within single incidence configurations, the corresponding ranges were 0.935 to 0.946, and 0.14 to 0.15 for Checkerboard. The Friedman test rejected equal ranks ($\chi^2 = 703.76$, $\text{df} = 5$, 160 blocks, $p < 2.2\times10^{-16}$). Nemenyi comparisons separated 13 of 15 pairs; the two they did not separate were Incidence-Guided Saturation against Balanced Quartiles ($p = 0.74$) and $2\times2$ Blocking against Isolation Buffer ($p = 1.00$). Holm-adjusted Wilcoxon tests separated every pair except $2\times2$ Blocking against Isolation Buffer. In the surface-level comparison, Incidence-Guided Saturation had a lower MSE than Balanced Quartiles in 48 of the 50 units (mean difference $-0.0121$; $-0.0139$, $-0.0102$; Wilcoxon $p = 9.9\times10^{-10}$). $2\times2$ Blocking and Isolation Buffer were close: $2\times2$ Blocking had the lower MSE in 36 of the 50 units (72\%; mean difference $-0.0048$; $-0.0078$, $-0.0017$; Wilcoxon $p = 0.0027$), but Isolation Buffer had the slightly better average rank (3.38 against 3.40), and the two designs traded places at some values of $\tau$ (see the results for $\tau$ below). Every other adjacent pair differed in all 50 units.
>
> Checkerboard under rook contiguity is a different case: $\tau$ is not identified, and the coefficient on $Z$ estimates $\tau - \gamma$ (Methods). The estimates behaved accordingly. At $\tau = 1$, Checkerboard's mean bias was $-0.502$, $-0.600$, $-0.704$, and $-0.802$ at $\gamma = 0.5$, 0.6, 0.7, and 0.8, while its mean standard deviation stayed between 0.20 and 0.21; squared bias, averaged over scenarios, accounted for 0.438 of its MSE of 0.480 (Figure \ref{fig:biasvar-rook} in Appendix \ref{app:rook-figures} shows the decomposition for the Poisson configuration with $\rho_X = 0.20$). Its coverage fell from 0.316 at $\gamma = 0.5$ to 0.027 at $\gamma = 0.8$. Because $WZ = 1 - Z$, the two spillover covariates in Equation \ref{eq:spillover} coincide, $\gamma\, WZ \circ (1 - Z) = \gamma(1 - Z) = \gamma\, WZ$, so the two regimes define the same outcome model and the same expected MSE; the observed values, 0.477 under the both-arms regime and 0.483 under the control-only regime, differ only by Monte Carlo error, because the two regimes use separate noise draws. Under queen contiguity, by contrast, Checkerboard's coverage stayed between 0.939 and 0.943 at every level of $\rho$; under rook contiguity it was between 0.134 and 0.161 (Table \ref{tab:checkerboard-coverage}). Rook Checkerboard estimates were precise but centered on the wrong quantity, the opposite of its imprecise but valid queen estimates. We keep these estimates in every ranking and flag them.
>
> ```{=latex}
> \begin{table}[!htb]
>
> \small
> \centering
>
> \caption{Checkerboard's 95\% confidence interval coverage across the four levels of outcome spatial dependence $\rho$ at $\tau = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages 40 scenarios (5 incidence configurations $\times$ 4 values of $\gamma$ $\times$ 2 spillover regimes). Under rook contiguity, $\tau$ is not identified for Checkerboard; its estimates are kept and flagged. The nominal level is 0.95.}
> \label{tab:checkerboard-coverage}
>
> \begin{tabular}{lcccc}
>
> \toprule
> Neighbor type & $\rho=0$ & $\rho=0.01$ & $\rho=0.20$ & $\rho=0.50$ \\
>
> \midrule
> Queen & 0.941 & 0.943 & 0.940 & 0.939 \\
> Rook & 0.161 & 0.156 & 0.148 & 0.134 \\
>
> \bottomrule
>
> \end{tabular}
>
> \end{table}
> ```

### 3.3.6 Non-oracle sensitivity analysis — condensed to one paragraph
Reason: design-by-design bias and coverage values are all in Table B.5 (Appendix B.6 in the concise copy).

Original text:

> The non-oracle model, which omits the spillover term, was biased for every design (Table \ref{tab:oracle-nonoracle} in Appendix \ref{app:oracle-nonoracle} gives both models' MSE, bias, and coverage for each design). At $\tau = 1$, and setting aside Checkerboard under rook contiguity, its bias ranged from $-0.12$ to $-0.40$ and its coverage from 0.51 to 0.90 across designs and neighbor types. Under queen contiguity its bias ranged from $-0.123$ (Incidence-Guided Saturation) to $-0.244$ ($2\times2$ Blocking) and its coverage from 0.709 ($2\times2$ Blocking) to 0.896 (High Incidence Focus); under rook contiguity the most biased design other than Checkerboard was Isolation Buffer (bias $-0.400$, coverage 0.513). For Checkerboard under rook contiguity, the non-oracle estimates were identical to the oracle estimates (MSE 0.480, bias $-0.652$, coverage 0.150), because once the aliased spillover column is dropped the two models are the same.
>
> Omitting the spillover term also removes a column that can be nearly collinear with $Z$, and the non-oracle estimates were accordingly less variable for some designs. Under queen contiguity, the non-oracle MSE was lower than the oracle MSE for every design, most of all for Checkerboard (oracle 1.087, non-oracle 0.108) and $2\times2$ Blocking (oracle 0.319, non-oracle 0.121). Under rook contiguity it was higher for every design except $2\times2$ Blocking (0.128 under both models) and Checkerboard (identical). These results do not change the recommendation. Incidence-Guided Saturation also had the lowest non-oracle MSE under both neighbor types (0.087 under queen contiguity and 0.086 under rook contiguity), and where a design's non-oracle MSE was lower, the estimate was biased and its intervals under-covered, so the lower MSE does not support valid inference.

### 3.3.7 Stability across treatment effect sizes — margins and coverage condensed
Reason: block counts by tau, average-rank trajectories, the tau = 2 worst case and the rook margin detail move to Appendix B.6.2 (Block-level detail); the intermediate rook-Checkerboard coverage values move there too.

Original text:

> Within that stable order, the margins moved. Under queen contiguity, the average rank of Incidence-Guided Saturation rose from 1.41 at $\tau = 0.8$ to 1.85 at $\tau = 3$, but Balanced Quartiles did not gain on it: Incidence-Guided Saturation had the lower MSE of the two in 103, 102, 102, 102, and 101 of the 160 blocks at $\tau = 0.8$, 1, 1.5, 2, and 3, and the average rank of Balanced Quartiles stayed between 2.13 and 2.14. Isolation Buffer gained instead, and its average rank fell from 3.71 to 3.31. The MSE of Isolation Buffer stayed at 0.160 at every value of $\tau$, while that of Incidence-Guided Saturation rose from 0.088 to 0.115 and that of Balanced Quartiles from 0.130 to 0.152, so the gap between Incidence-Guided Saturation and Balanced Quartiles stayed between 0.037 and 0.041. Isolation Buffer's gain came under the control-only regime, where Incidence-Guided Saturation had a lower MSE than Isolation Buffer in 75, 73, 56, 40, and 12 of the 80 blocks at the five values of $\tau$. The mean control-only MSE of Incidence-Guided Saturation rose from 0.131 at $\tau = 0.8$ to 0.182 at $\tau = 3$, while that of Isolation Buffer stayed between 0.159 and 0.160 (Appendix \ref{app:technical-notes}). At $\tau = 2$ the two were level (0.1596 for Isolation Buffer against 0.1604), and at $\tau = 3$ Isolation Buffer had the lower control-only MSE (0.160 against 0.182). Incidence-Guided Saturation remained second under the control-only regime at those two values, second to Balanced Quartiles under the both-arms regime at every value of $\tau$, and lowest in pooled MSE at every value of $\tau$. Its worst case across blocks was the lowest of the six designs at $\tau \le 1.5$, but at $\tau = 2$ and $\tau = 3$ that of Isolation Buffer was lower (0.207 against 0.210, and 0.208 against 0.238). The mean bias of Incidence-Guided Saturation also grew with $\tau$, from 0.029 at $\tau = 0.8$ to 0.047 at $\tau = 3$, but squared bias, averaged over scenarios, stayed at or below 3.0\% of its MSE, so the rise in its MSE was mostly variance.
>
> Under rook contiguity, Incidence-Guided Saturation kept the lowest control-only MSE at every value of $\tau$, although its margin over Isolation Buffer narrowed, from 0.095 against 0.134 at $\tau = 0.8$ to 0.131 against 0.136 at $\tau = 3$, where it had the lower MSE in 47 of the 80 control-only blocks. Its worst case across rook blocks was the lowest of the six designs at every value of $\tau$ except $\tau = 3$, where that of Isolation Buffer was lower (0.173 against 0.174).
>
> Except for Checkerboard under rook contiguity, coverage stayed close to 0.94 at every value of $\tau$. Pooled over the five incidence configurations, every design's coverage was between 0.936 and 0.941 at every value of $\tau$ under queen contiguity, and every design except Checkerboard was between 0.937 and 0.942 at every value of $\tau$ under rook contiguity. The coverage of rook Checkerboard, for which $\tau$ is not identified, rose from 0.138 at $\tau = 0.8$ to 0.150, 0.198, 0.269, and 0.442 at $\tau = 1$, 1.5, 2, and 3 (Appendix \ref{app:technical-notes}). Appendix \ref{app:supp-comparisons} shows coverage across $\tau$ for each design.

### 3.3.8 recommendation paragraph — cut (heading shortened to "Robustness across blocks")
Reason: the recommendation and the Saturation Quadrants statement are given in the SRS benchmark subsection, the application and the Discussion; the closing distinction between scenario-level and allocation-level worst cases is kept at the end of the robustness paragraph.

Original text:

> Among the six designs in the detailed comparison, Incidence-Guided Saturation combines low pooled MSE with comparatively favorable control-only performance. The SRS benchmark below changes the practical interpretation: recommend the saturation approach when control-only spillover is plausible, while retaining SRS as a competitive option under both-arms spillover. Saturation Quadrants performed similarly, so incidence guidance lets a trial direct its highest treatment saturations to the highest-incidence quadrants with little or no loss of estimation accuracy. Scenario-level worst cases summarize variation across parameter settings; they are distinct from conditional risk across assignments within one setting.

### New Appendix B.4 "Additional sensitivity tables"
Reason: holds the three moved tables unchanged; later appendix sections renumber (B.4-B.7 become B.5-B.8; labels unchanged).

Note: Table 3.4/3.5/3.7 order in B.4 is rho, Checkerboard coverage, gamma, with a float barrier so all three stay in B.4.

## 2026-10-12 — Batch 3 (application)

### 3.4.2 Geographic adaptations and treatment budgets — procedural detail moved to Appendix B.8
Reason: the step-by-step adaptation rules (graph-coloring order, stratum rounding, the 250-candidate regional search, regional ranking rationale and the 2x2 group construction) now appear in B.8 under "Design adaptations"; the body keeps one line per design and the budget table.

Original text:

> Queen contiguity between service areas, computed from legal county boundaries, defines primary adjacency, and rook contiguity is a sensitivity analysis. Before simulation, we checked the 58-cluster row-standardized contiguity matrix $W$, its binary adjacency and its row order against the boundary geometry. Shared legal boundaries across water can create adjacency; this is a prespecified geographic exposure proxy, not a measured professional communication network.
>
> All eight designs and SRS are evaluated under the names used on the grid, adapted to the map as follows. On the NC map, Checkerboard becomes a deterministic greedy graph-coloring assignment: visit clusters by decreasing degree, choose the arm that conflicts least with assigned neighbors, break equal conflicts toward control, and complement the assignment if more than half are treated. This is an interspersion analog, not a literal rectangular checkerboard or an exact bipartite coloring. High Incidence Focus treats the 29 highest annual-rate clusters. Isolation Buffer draws a greedy random maximal independent set on the NC adjacency graph. SRS treats a random 29 of 58 clusters. Balanced Quartiles and Balanced Halves stratify annual rate ranks and treat 29 overall; odd-size strata receive their necessary extra treatments at random. In Balanced Halves, one half receives 14 treatments and the other 15, with the extra-treatment half chosen randomly.
>
> The two saturation designs use one fixed four-region partition for all years. A geographic $k$-means procedure generates 250 candidates, selecting the partition by compactness and balance in reference population, cluster count and county count, without consulting outcomes. Reference population is pooled person-years divided by four. Saturation Quadrants randomly assigns 20/40/60/80\% to these regions each draw. Incidence-Guided Saturation first ranks clusters by annual per-capita rates, averages their ranks equally within each region, ranks those four means, and assigns 80/60/40/20\% from highest to lowest. Regional ties are broken randomly per allocation. Local treatment counts are rounded under the retained design rules. Per-capita cluster rates already account for population size; regional mean ranks answer an equally weighted cluster question, whereas pooled regional deaths/population weights larger clusters more heavily. Mean raw rates and pooled regional rates are evaluated as prespecified sensitivities rather than chosen by performance. The grid rule averages $X$ directly: mean ranks align with its rank-scaled Poisson cases but are a deliberate adaptation for the application, rather than the exact rule used for continuous grid incidence.
>
> On the NC map, $2\times2$ Blocking uses spatial blocks of nearby clusters: 14 irregular groups, formed once per year by $k$-means on standardized projected coordinates, incidence ranks and log(1+population), then held fixed. Half of each group, rounded, is treated randomly. This is not a set of literal square blocks. In each year one single-cluster group remains control, because half of one cluster rounds to zero. Budgets differ across designs; we report them rather than modify the designs to force equal budgets (Table \ref{tab:nc-budgets}). Appendix \ref{app:nc} maps the fixed regions and example assignments.

### 3.4.3 Outcome model and Monte Carlo evaluation — replication and tail-method detail condensed
Reason: the adaptive schedule, the R = 2,000 exception, shared-simulation accounting and the duplicate-allocation covariance are described in B.8; the body keeps J, R, the precision targets, the shared-years caveat, the definition of m(a) and the tail-refinement design.

Original text:

> We initially sampled $J=100$ assignments and simulated $R=100$ independent outcomes per distinct assignment. A design whose set of possible assignments contains only one assignment instead received $R=1,000$ outcomes for that assignment, rather than 100 repeated copies of it counted as separate allocations. An adaptive, prespecified replication schedule extended allocations or outcomes until the Monte Carlo standard error of mean MSE was at most 5\% of mean MSE and that of coverage at most 0.01. Two single-assignment settings in the matched sensitivity analysis required $R=2,000$; all randomized designs met these targets at $J=100$ and $R=100$. Each random draw used a fixed, reproducible seed. Settings with the same set of possible assignments shared their simulations, with annual treatment counts and population shares recalculated. In the primary model, Checkerboard, Saturation Quadrants, Isolation Buffer and SRS reuse the same simulations across years, because annual incidence enters neither their allocation nor the outcome model. These repeated values are not independent annual confirmation.
>
> For each allocation $a$, the repeated outcomes estimate the conditional MSE $m(a)=E_\epsilon[(\hat\tau-\tau)^2\mid a,\text{observed inputs}]$. We report its distribution across sampled allocations, its estimated 90th percentile and the estimated mean among the worst 10\% of sampled allocations, alongside unconditional mean MSE, bias and coverage. An allocation drawn more than once keeps its sampled frequency but uses one conditional estimate, and the Monte Carlo uncertainty accounts for the resulting covariance. A separate tail refinement at queen $\rho=0,0.5$, $\gamma=0.8$, both regimes and all years and designs kept the same allocations and extended the randomized designs to $R=400$ outcomes per allocation, keeping $R=1,000$ for single-assignment designs. Worst-decile assignments were identified with one half of the simulated outcomes and evaluated with the other half, to avoid selection bias. The refinement is not a new allocation sample. The rank of the fitted design matrix, residual treatment information, aliased columns, warnings, invalid estimates or intervals, and estimates of $\rho$ at the boundary of its parameter space were recorded. Appendix \ref{app:nc} gives the reproducibility and precision details.

### 3.4.4 repeated saturation sentence — shortened
Reason: the "little or no loss" statement is made in Methods and the Discussion.

Original text:

> Neither is uniquely best: the two performed similarly, so incidence guidance lets a trial direct its highest treatment saturations to the highest-incidence regions with little or no loss of estimation accuracy, and their advantage over SRS comes from the saturation structure both designs share.

### 3.4.4 population-share paragraph — condensed
Reason: the budgets are in Table 3.10; the 2x2 yearly counts are in B.8.

Original text:

> Population shares contextualize performance. High Incidence Focus treats 29 clusters but only about 22--24\% of working-age population, while Isolation Buffer treats a mean of 15.89 clusters and about 27\%. These population denominators measure geographic reach, not the number of healthcare professionals educated, intervention cost or participant-level precision. $2\times2$ Blocking treatment counts are 28/31/27/29 by year; complete budgets and shares are reported in Table \ref{tab:nc-budgets} and the appendix.

### 3.4.5 split-half sentence — shortened
Reason: the split-half method is defined in 3.4.3.

Original text:

> The pattern persists when worst-decile assignments are identified with one half of the simulated outcomes and evaluated with the other half, which avoids selection bias.

### B.8 Figure "Annual population shares treated" (fig:nc-pop) — cut from the chapter
Reason: it repeats the population shares of Table 3.10 (body); it remains in the Supplementary Information.

Original text:

> \begin{figure}[!htb]
> \centering
> \includegraphics[width=\linewidth]{figures/real_sud/population_shares.pdf}
>
> \caption{Annual population shares treated under primary queen assignments. Designs whose assignments do not depend on annual incidence share performance across years, but their annual population shares are recalculated. Population aged 18--64 is geographic reach rather than education-participant enrollment.}
> \label{fig:nc-pop}
>
> \end{figure}

### B.8 Figure "mean MSE by design, planning year and spillover regime" (fig:nc-mse-detail) — cut from the chapter
Reason: it plots the absolute MSEs listed in the merged annual table and summarized as SRS ratios in Figure 3.8; it remains in the Supplementary Information.

Original text:

> \begin{figure}[!htb]
>
> \centering
> \includegraphics[width=\linewidth]{figures/real_sud/primary_yearly_mse.pdf}
>
> \caption{Completed primary application: mean MSE by design, planning year and spillover regime, queen contiguity and education outcome. Each entry equally averages the six $\rho/\gamma$ pairs. Both regime columns use a common logarithmic MSE scale. Years that share the same simulations are repeated reporting, not independent replications; pooled Monte Carlo intervals are not supplied.}
> \label{fig:nc-mse-detail}
>
> \end{figure}

### B.8 four annual tables (tab:nc-full-2018 ... 2021) — merged into tab:nc-annual
Reason: author-approved merge: one table lists MSE, bias and coverage for every design in every year (36 rows); the Treated and Pop.% columns repeated Table 3.10 and are dropped here. Values regenerated from yearly_primary_design_means.csv with the original rounding; the verifier checks them.

Original text:

> ```{=latex}
> \begin{table}[!htb]
> \small
> \centering
>
> \caption{Primary queen 2018 application performance. Equal weights over six parameter pairs within each regime; bias, MSE and coverage are Monte Carlo estimates. Treatment count and population share describe the allocation budget.}
> \label{tab:nc-full-2018}
> % SOURCE: application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv
>
> \begin{tabular}{p{3.5cm}rrrrr}
> \toprule
> Design & MSE & Bias & Coverage & Treated & Pop.\% \\
> \midrule
> \multicolumn{6}{l}{\textit{Both arms}} \\
> Checkerboard & 0.21591 & -0.0074 & 0.936 & 27.00 & 50.4 \\
> High Incidence Focus & 0.07172 & 0.0160 & 0.934 & 29.00 & 22.2 \\
> Saturation Quadrants & 0.07518 & 0.0208 & 0.935 & 29.03 & 49.2 \\
> Isolation Buffer & 0.28591 & 0.0088 & 0.932 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.07609 & 0.0099 & 0.933 & 28.00 & 50.0 \\
> Balanced Quartiles & 0.07610 & 0.0170 & 0.933 & 29.00 & 48.5 \\
> Balanced Halves & 0.07458 & 0.0158 & 0.934 & 29.00 & 50.0 \\
> Incidence-Guided Saturation & 0.07428 & 0.0174 & 0.935 & 29.00 & 53.2 \\
> Simple Random Sampling & 0.07638 & 0.0135 & 0.932 & 29.00 & 49.5 \\
> \multicolumn{6}{l}{\textit{Control only}} \\
> Checkerboard & 0.92036 & 0.0100 & 0.933 & 27.00 & 50.4 \\
> High Incidence Focus & 0.21426 & 0.0237 & 0.932 & 29.00 & 22.2 \\
> Saturation Quadrants & 0.22329 & 0.0464 & 0.926 & 29.03 & 49.2 \\
> Isolation Buffer & 0.28575 & 0.0106 & 0.931 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.30622 & 0.0281 & 0.929 & 28.00 & 50.0 \\
> Balanced Quartiles & 0.29352 & 0.0315 & 0.929 & 29.00 & 48.5 \\
> Balanced Halves & 0.29101 & 0.0395 & 0.930 & 29.00 & 50.0 \\
> Incidence-Guided Saturation & 0.24285 & 0.0464 & 0.928 & 29.00 & 53.2 \\
> Simple Random Sampling & 0.30761 & 0.0378 & 0.929 & 29.00 & 49.5 \\
>
> \bottomrule
> \end{tabular}
> \end{table}
> ```
>
> ```{=latex}
> \begin{table}[!htb]
> \small
> \centering
>
> \caption{Primary queen 2019 application performance. Equal weights over six parameter pairs within each regime; bias, MSE and coverage are Monte Carlo estimates. Treatment count and population share describe the allocation budget.}
> \label{tab:nc-full-2019}
> % SOURCE: application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv
>
> \begin{tabular}{p{3.5cm}rrrrr}
> \toprule
> Design & MSE & Bias & Coverage & Treated & Pop.\% \\
> \midrule
> \multicolumn{6}{l}{\textit{Both arms}} \\
> Checkerboard & 0.21591 & -0.0074 & 0.936 & 27.00 & 50.4 \\
> High Incidence Focus & 0.07180 & 0.0063 & 0.936 & 29.00 & 23.5 \\
> Saturation Quadrants & 0.07518 & 0.0208 & 0.935 & 29.03 & 49.2 \\
> Isolation Buffer & 0.28591 & 0.0088 & 0.932 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.07808 & 0.0143 & 0.930 & 31.00 & 55.1 \\
> Balanced Quartiles & 0.07624 & 0.0175 & 0.933 & 29.00 & 49.7 \\
> Balanced Halves & 0.07517 & 0.0167 & 0.934 & 29.00 & 50.5 \\
> Incidence-Guided Saturation & 0.07628 & 0.0233 & 0.935 & 29.00 & 45.1 \\
> Simple Random Sampling & 0.07638 & 0.0135 & 0.932 & 29.00 & 49.5 \\
> \multicolumn{6}{l}{\textit{Control only}} \\
> Checkerboard & 0.92036 & 0.0100 & 0.933 & 27.00 & 50.4 \\
> High Incidence Focus & 0.43019 & 0.0227 & 0.931 & 29.00 & 23.5 \\
> Saturation Quadrants & 0.22329 & 0.0464 & 0.926 & 29.03 & 49.2 \\
> Isolation Buffer & 0.28575 & 0.0106 & 0.931 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.43630 & 0.0314 & 0.928 & 31.00 & 55.1 \\
> Balanced Quartiles & 0.30616 & 0.0367 & 0.929 & 29.00 & 49.7 \\
> Balanced Halves & 0.29955 & 0.0360 & 0.930 & 29.00 & 50.5 \\
> Incidence-Guided Saturation & 0.21565 & 0.0522 & 0.926 & 29.00 & 45.1 \\
> Simple Random Sampling & 0.30761 & 0.0378 & 0.929 & 29.00 & 49.5 \\
>
> \bottomrule
> \end{tabular}
> \end{table}
> ```
>
> ```{=latex}
> \begin{table}[!htb]
> \small
> \centering
>
> \caption{Primary queen 2020 application performance. Equal weights over six parameter pairs within each regime; bias, MSE and coverage are Monte Carlo estimates. Treatment count and population share describe the allocation budget.}
> \label{tab:nc-full-2020}
> % SOURCE: application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv
>
> \begin{tabular}{p{3.5cm}rrrrr}
> \toprule
> Design & MSE & Bias & Coverage & Treated & Pop.\% \\
> \midrule
> \multicolumn{6}{l}{\textit{Both arms}} \\
> Checkerboard & 0.21591 & -0.0074 & 0.936 & 27.00 & 50.4 \\
> High Incidence Focus & 0.07292 & 0.0094 & 0.931 & 29.00 & 23.3 \\
> Saturation Quadrants & 0.07518 & 0.0208 & 0.935 & 29.03 & 49.1 \\
> Isolation Buffer & 0.28591 & 0.0088 & 0.932 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.07764 & 0.0122 & 0.933 & 27.00 & 48.0 \\
> Balanced Quartiles & 0.07530 & 0.0158 & 0.933 & 29.00 & 50.6 \\
> Balanced Halves & 0.07466 & 0.0150 & 0.934 & 29.00 & 49.2 \\
> Incidence-Guided Saturation & 0.07682 & 0.0209 & 0.933 & 29.00 & 47.9 \\
> Simple Random Sampling & 0.07638 & 0.0135 & 0.932 & 29.00 & 49.5 \\
> \multicolumn{6}{l}{\textit{Control only}} \\
> Checkerboard & 0.92036 & 0.0100 & 0.933 & 27.00 & 50.4 \\
> High Incidence Focus & 0.25918 & 0.0327 & 0.930 & 29.00 & 23.3 \\
> Saturation Quadrants & 0.22329 & 0.0464 & 0.926 & 29.03 & 49.1 \\
> Isolation Buffer & 0.28575 & 0.0106 & 0.931 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.31848 & 0.0281 & 0.930 & 27.00 & 48.0 \\
> Balanced Quartiles & 0.28932 & 0.0356 & 0.927 & 29.00 & 50.6 \\
> Balanced Halves & 0.30610 & 0.0372 & 0.930 & 29.00 & 49.2 \\
> Incidence-Guided Saturation & 0.22029 & 0.0474 & 0.926 & 29.00 & 47.9 \\
> Simple Random Sampling & 0.30761 & 0.0378 & 0.929 & 29.00 & 49.5 \\
>
> \bottomrule
> \end{tabular}
> \end{table}
> ```
>
> ```{=latex}
> \begin{table}[!htb]
> \small
> \centering
>
> \caption{Primary queen 2021 application performance. Equal weights over six parameter pairs within each regime; bias, MSE and coverage are Monte Carlo estimates. Treatment count and population share describe the allocation budget.}
> \label{tab:nc-full-2021}
> % SOURCE: application/results/real_sud_rev_20261002/exhibits/yearly_primary_design_means.csv
>
> \begin{tabular}{p{3.5cm}rrrrr}
> \toprule
> Design & MSE & Bias & Coverage & Treated & Pop.\% \\
> \midrule
> \multicolumn{6}{l}{\textit{Both arms}} \\
> Checkerboard & 0.21591 & -0.0074 & 0.936 & 27.00 & 50.7 \\
> High Incidence Focus & 0.07031 & 0.0104 & 0.934 & 29.00 & 22.4 \\
> Saturation Quadrants & 0.07518 & 0.0208 & 0.935 & 29.03 & 49.1 \\
> Isolation Buffer & 0.28591 & 0.0088 & 0.932 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.07736 & 0.0107 & 0.931 & 29.00 & 50.3 \\
> Balanced Quartiles & 0.07453 & 0.0175 & 0.934 & 29.00 & 49.8 \\
> Balanced Halves & 0.07584 & 0.0160 & 0.932 & 29.00 & 49.4 \\
> Incidence-Guided Saturation & 0.07450 & 0.0197 & 0.936 & 29.00 & 51.8 \\
> Simple Random Sampling & 0.07638 & 0.0135 & 0.932 & 29.00 & 49.5 \\
> \multicolumn{6}{l}{\textit{Control only}} \\
> Checkerboard & 0.92036 & 0.0100 & 0.933 & 27.00 & 50.7 \\
> High Incidence Focus & 0.35318 & 0.0296 & 0.929 & 29.00 & 22.4 \\
> Saturation Quadrants & 0.22329 & 0.0464 & 0.926 & 29.03 & 49.1 \\
> Isolation Buffer & 0.28575 & 0.0106 & 0.931 & 15.89 & 26.7 \\
> $2\times2$ Blocking & 0.36168 & 0.0289 & 0.928 & 29.00 & 50.3 \\
> Balanced Quartiles & 0.28439 & 0.0354 & 0.929 & 29.00 & 49.8 \\
> Balanced Halves & 0.29089 & 0.0340 & 0.929 & 29.00 & 49.4 \\
> Incidence-Guided Saturation & 0.23732 & 0.0451 & 0.926 & 29.00 & 51.8 \\
> Simple Random Sampling & 0.30761 & 0.0378 & 0.929 & 29.00 & 49.5 \\
>
> \bottomrule
> \end{tabular}
> \end{table}
> ```
