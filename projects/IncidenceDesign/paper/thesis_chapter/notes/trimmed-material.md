# Chapter 3 trimmed material

Text removed or relocated while preparing the concise Chapter 3 copy
(`Dissertation_Chapter_concise.qmd`, branch `ch3-concise`, started 2026-10-12;
plan: `docs/plans/2026-10-12-chapter3-length-review-plan.md`). Nothing here is
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
Reason: the rook test statistics and within-configuration coverage ranges move to Appendix B.3; Table 3.7 (Checkerboard coverage by rho) moves to the new Appendix B.4 and its range is stated in the text. The Checkerboard tau - gamma explanation is kept in short form (it is derived in Methods and concise B.7.1, canonical B.6.1).

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
Reason: design-by-design bias and coverage values are all in the oracle/non-oracle table (tab:oracle-nonoracle; canonical B.5, concise B.6).

Original text:

> The non-oracle model, which omits the spillover term, was biased for every design (Table \ref{tab:oracle-nonoracle} in Appendix \ref{app:oracle-nonoracle} gives both models' MSE, bias, and coverage for each design). At $\tau = 1$, and setting aside Checkerboard under rook contiguity, its bias ranged from $-0.12$ to $-0.40$ and its coverage from 0.51 to 0.90 across designs and neighbor types. Under queen contiguity its bias ranged from $-0.123$ (Incidence-Guided Saturation) to $-0.244$ ($2\times2$ Blocking) and its coverage from 0.709 ($2\times2$ Blocking) to 0.896 (High Incidence Focus); under rook contiguity the most biased design other than Checkerboard was Isolation Buffer (bias $-0.400$, coverage 0.513). For Checkerboard under rook contiguity, the non-oracle estimates were identical to the oracle estimates (MSE 0.480, bias $-0.652$, coverage 0.150), because once the aliased spillover column is dropped the two models are the same.
>
> Omitting the spillover term also removes a column that can be nearly collinear with $Z$, and the non-oracle estimates were accordingly less variable for some designs. Under queen contiguity, the non-oracle MSE was lower than the oracle MSE for every design, most of all for Checkerboard (oracle 1.087, non-oracle 0.108) and $2\times2$ Blocking (oracle 0.319, non-oracle 0.121). Under rook contiguity it was higher for every design except $2\times2$ Blocking (0.128 under both models) and Checkerboard (identical). These results do not change the recommendation. Incidence-Guided Saturation also had the lowest non-oracle MSE under both neighbor types (0.087 under queen contiguity and 0.086 under rook contiguity), and where a design's non-oracle MSE was lower, the estimate was biased and its intervals under-covered, so the lower MSE does not support valid inference.

### 3.3.7 Stability across treatment effect sizes — margins and coverage condensed
Reason: block counts by tau, average-rank trajectories, the tau = 2 worst case and the rook margin detail move to the technical-notes appendix (concise B.7.2, canonical B.6.2; Block-level detail); the intermediate rook-Checkerboard coverage values move there too.

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
Reason: the budgets are in the body budget table (tab:nc-budgets; canonical Table 3.10); the 2x2 yearly counts are in B.8.

Original text:

> Population shares contextualize performance. High Incidence Focus treats 29 clusters but only about 22--24\% of working-age population, while Isolation Buffer treats a mean of 15.89 clusters and about 27\%. These population denominators measure geographic reach, not the number of healthcare professionals educated, intervention cost or participant-level precision. $2\times2$ Blocking treatment counts are 28/31/27/29 by year; complete budgets and shares are reported in Table \ref{tab:nc-budgets} and the appendix.

### 3.4.5 split-half sentence — shortened
Reason: the split-half method is defined in 3.4.3.

Original text:

> The pattern persists when worst-decile assignments are identified with one half of the simulated outcomes and evaluated with the other half, which avoids selection bias.

### B.8 Figure "Annual population shares treated" (fig:nc-pop) — cut from the chapter
Reason: it repeats the population shares of the body budget table (tab:nc-budgets); it remains in the Supplementary Information.

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
Reason: author-approved merge: one table lists MSE, bias and coverage for every design in every year (36 rows); the Treated and Pop.% columns repeated the body budget table (tab:nc-budgets) and are dropped here. Values regenerated from yearly_primary_design_means.csv with the original rounding; the verifier checks them.

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

## 2026-10-12 — Batch 4 (discussion)

### 3.5.1 Summary and relevance — restated results condensed
Reason: the queen/rook, regime and Checkerboard paragraphs restated Results numbers at length; the condensed version keeps each finding with headline numbers, but drops the Checkerboard detail values (bias 0.004, coverage 0.942, rook coverage 0.150, MSE growing with gamma), all reported in the Results.

Original text:

> Under queen contiguity, Incidence-Guided Saturation had the lowest mean squared error (0.091) and Balanced Quartiles the second lowest (0.131). The order of the six designs was the same in each of the five incidence configurations and at every value of $\tau$, and neither the strength of spatial dependence in the outcome nor the spillover magnitude changed it. The mechanism that generated the incidence surface changed each design's MSE modestly (by under 10\%, except High Incidence Focus). The same two designs also led under rook contiguity (0.072 and 0.084). Coverage of the 95\% intervals was close to 0.94 for every design under both neighbor types, except Checkerboard under rook contiguity, so the designs differed in precision and coverage did not separate them.
>
> Under queen contiguity, Checkerboard and $2\times2$ Blocking were close to the best designs when spillover reached both arms and were the two worst when it reached only control clusters: at $\tau = 1$, Checkerboard's MSE was 0.062 under the both-arms regime and 2.111 under the control-only regime, and that of $2\times2$ Blocking was 0.062 and 0.576. Balanced Quartiles also lost ground under the control-only regime (0.042 and 0.219), while Incidence-Guided Saturation was first or second among the six detailed designs under both regimes at every value of $\tau$ and remained competitive with SRS under both-arms spillover. Its control-only MSE did rise with $\tau$, and Isolation Buffer, whose MSE did not, matched it under that regime at $\tau = 2$ and passed it at $\tau = 3$. A design whose performance depends on the regime carries a risk that its average MSE does not show.
>
> Checkerboard failed for different reasons under the two neighbor types. Under queen contiguity, its high MSE was variance, not bias: under the control-only regime, only the boundary clusters separate the spillover effect from the treatment effect (Methods). The estimates stayed centered on $\tau$ (bias 0.004 under the control-only regime at $\tau = 1$) and their intervals kept close to nominal coverage (0.942), but they were imprecise. Under rook contiguity, $\tau$ is not identified (Methods): the estimates were centered on $\tau - \gamma$, so rook Checkerboard's MSE grew with $\gamma$, and its coverage at $\tau = 1$ was 0.150. Checkerboard intersperses treated and control clusters as finely as the grid allows; both failures come from the dependence this pattern creates between $Z$ and the spillover covariate, which is a linear combination of the intercept and $Z$ at every cluster under rook contiguity and at every interior cluster under queen contiguity and the control-only regime.

### 3.5.1 Saturation Quadrants paragraph — shortened
Reason: same statement as in Methods and the application; kept once in short form here.

Original text:

> Saturation Quadrants, which assigns the same saturations without reference to incidence, performed similarly to Incidence-Guided Saturation (see Treatment assignment designs). Both designs treat the four quadrants at saturations from 20\% to 80\%, and their accuracy advantage over SRS under control-only spillover comes from this shared saturation structure. Incidence guidance therefore lets a trial direct its highest treatment saturations to the highest-incidence regions with little or no loss of estimation accuracy.

### 3.5.3 Limitation 1 — numbers trimmed
Reason: the non-oracle bias and coverage ranges are reported in the Results and the oracle/non-oracle table (tab:oracle-nonoracle).

Original text:

> We fit a non-oracle model that omits the spillover term as a sensitivity analysis. At $\tau = 1$, setting aside Checkerboard under rook contiguity, it was biased for every design, with bias between $-0.12$ and $-0.40$ and coverage between 0.51 and 0.90, and Incidence-Guided Saturation still had the lowest MSE under both neighbor types.

### 3.5.3 Limitation 6 — shortened
Reason: the reason for keeping the base rate is given in Heterogeneous baseline incidence.

Original text:

> Sixth, the simulation's base rate of 35 per 100,000 is not the rate of the county data described in the Application (see Heterogeneous baseline incidence). We kept it deliberately, because the simulation is a general design study, not a calibration to one data set. In the Poisson mechanism, incidence enters the designs and the fitted model only through the ranks of the observed rates, so the base rate affects the simulation only through the spatial pattern of those ranks and the number of tied rates.

## 2026-10-12 — Batch 5 (figure placement)

### Queen coverage figure (fig:coverage) and bias-variance figure (fig:biasvar) — moved to Appendix B.5
Reason: both support statements made in the text (coverage 0.936-0.941 for every design; MSE increase is variance, not bias) and sit next to their rook counterparts; B.5 retitled "Coverage, bias-variance and rook contiguity figures" (label app:rook-figures unchanged). Captions unchanged except the pointer to the rook figure.

## 2026-10-13 — Review fixes and grid figure

### Review fixes (independent read-only review, 2026-10-13)
Reason: three meaning-drift errors (queen-only order claim in the Discussion, replication schedule wording, antecedent for the alternative regional summaries), worst-case/regime and tau qualifiers, the non-oracle rook contrast, the bivand citation and intro wording restored; block-level detail moved to the end of B.7; duplicate 2x2 statement removed. Grid-accuracy figure (CTJ Figure 2, grid_accuracy_srs.pdf) added to 3.3.9 at the author's request.

Numbering note (2026-10-13): appendix and table numbers above are labelled canonical or concise where they differ; LaTeX labels are given where numbering shifts.

## 2026-11-05 — Round 3 human-writer review (AI tropes)

Reason: Round 3 pass on Chapter 3 (edit lists and reviews in bios-dissertation `prelim/full-draft/notes/human-writer-review/round3/ch3/`). Each entry is a sentence or clause cut as a repeated caveat, recap, pointer, counted package or announcement; the point it made is still stated once where it applies. Matching SI edits are in `ch3-mirror-final-edits.md`.

### E3.tb.meth.2
Reason: cut: "How much treatment it redirects depends on how much incidence differs between regions (see Application)." The Discussion (line 696) says the same with the numbers, and line 548 shows it; in Methods it anticipates a result.

Original text:

> to the quadrants or regions with the highest incidence. How much treatment it redirects depends on how much incidence differs between regions (see Application). The accuracy advantage

Replaced by:

> to the quadrants or regions with the highest incidence. The accuracy advantage

### E3.tb.meth.4
Reason: cut: "Because the values of $\tau$ share assignments and noise, Friedman tests at different values of $\tau$ are not independent confirmations of one another." Line 151 already states the sharing and line 372, where the tau results are read, gives the same caveat with its consequence; Limitations repeats it again. Removes two $\tau$ tokens (intended).

Original text:

> with Holm's adjustment over the eight designs compared with SRS. Because the values of $\tau$ share assignments and noise, Friedman tests at different values of $\tau$ are not independent confirmations of one another.

Replaced by:

> with Holm's adjustment over the eight designs compared with SRS.

### E3.tb.res.2
Reason: cut: "The mechanism that generated the incidence surface therefore did not change the ranking." It restates the paragraph's opening sentence; the High Incidence Focus point is kept. SHARED: SI (M3.tb.2).

Original text:

> The mechanism that generated the incidence surface therefore did not change the ranking. The same held for High Incidence Focus, the design that uses the incidence ranks most directly: its MSE was no higher

Replaced by:

> High Incidence Focus, the design that uses the incidence ranks most directly, also had no higher MSE

### E3.tb.res.3
Reason: cut: "Its stable rank across regimes matters because a trial team may not know in advance which regime applies." Announces importance; the Discussion (lines 680 and 686) makes the unknown-regime point where it belongs.

Original text:

> Its stable rank across regimes matters because a trial team may not know in advance which regime applies. In the full benchmark

Replaced by:

> In the full benchmark

### E3.tb.res.4
Reason: cut: "The advantage of the saturation designs over SRS is therefore specific to control-only spillover, and under queen contiguity their small both-arms disadvantage is statistically detectable." Restates the numbers two sentences up and line 274's "small but statistically detectable difference". SHARED: SI (M3.tb.4).

Original text:

> and Checkerboard under control-only spillover. The advantage of the saturation designs over SRS is therefore specific to control-only spillover, and under queen contiguity their small both-arms disadvantage is statistically detectable.

Replaced by:

> and Checkerboard under control-only spillover.

### E3.tb.res.6
Reason: cut: "The mean is only part of the comparison, because a trial realizes a single allocation." Methods line 141 makes this the stated reason for the second criterion; the subsection can open with the result. Not mirrored: the SI has no earlier statement of this reason for the grid.

Original text:

> The mean is only part of the comparison, because a trial realizes a single allocation. Under control-only spillover, the worst 10\% of SRS

Replaced by:

> Under control-only spillover, the worst 10\% of SRS

### E3.tb.app.1
Reason: cut: "Monte Carlo precision describes uncertainty in simulated performance, not interval calibration." Line 544, two paragraphs later, makes the same point with the coverage numbers.

Original text:

> gives 9,878,000 distinct outcomes across the two stages. Monte Carlo precision describes uncertainty in simulated performance, not interval calibration.

Replaced by:

> gives 9,878,000 distinct outcomes across the two stages.

### E3.tb.app.2
Reason: The shared-structure explanation appears in Methods (line 90) and the Discussion (line 696); cut the clause from the middle instance.

Original text:

> as on the grid, the two performed similarly, and their advantage over SRS comes from the saturation structure both designs share.

Replaced by:

> as on the grid, the two performed similarly.

### E3.tb.app.5
Reason: cut: "This sensitivity analysis changes baseline and adjustment together and does not imply that SUD incidence predicts education ability." Near-verbatim repeat of line 520, where the analysis is set up. SHARED: SI (M3.tb.6; the SI has the setup sentence at its line 829).

Original text:

> leaves less residual treatment information. This sensitivity analysis changes baseline and adjustment together and does not imply that SUD incidence predicts education ability. Comparing rook

Replaced by:

> leaves less residual treatment information. Comparing rook

### E3.tb.app.7
Reason: cut: "The findings do not establish intervention effectiveness or a final trial protocol." Effectiveness is disclaimed at lines 58 and 447 and in Limitations (716); the protocol point is the preceding sentence. Not mirrored: the SI states the effectiveness caveat only once nearby (its line 625 is 1,800 lines earlier).

Original text:

> before selecting one realized allocation. The findings do not establish intervention effectiveness or a final trial protocol.

Replaced by:

> before selecting one realized allocation.

### E3.tb.disc.3
Reason: cut: "Together, the grid benchmark, the allocation-risk results (...) and the NC application give guidance by scenario." The heading already says this; the table reference moves to the claim it supports.

Original text:

> Together, the grid benchmark, the allocation-risk results (Table \ref{tab:downside}) and the NC application give guidance by scenario. When spillover is expected to reach mainly control clusters, Incidence-Guided Saturation and Saturation Quadrants had the lowest mean MSE and among the lowest worst-decile MSE, so

Replaced by:

> When spillover is expected to reach mainly control clusters, Incidence-Guided Saturation and Saturation Quadrants had the lowest mean MSE and among the lowest worst-decile MSE (Table \ref{tab:downside}), so

### E3.tb.disc.5
Reason: cut: "The NC model assesses education outcomes rather than assuming that historical SUD incidence predicts education ability or a heterogeneous intervention effect." The incidence-education caveat is stated at lines 54, 447, 518 and 520, and the constant effect at 518 and 716.

Original text:

> did not improve on SRS on either criterion. The NC model assesses education outcomes rather than assuming that historical SUD incidence predicts education ability or a heterogeneous intervention effect.

Replaced by:

> did not improve on SRS on either criterion.

### E3.tb.disc.6
Reason: cut: "No direct cross-chapter result comparison is needed to motivate these extensions." Repeats the third sentence of the same paragraph ("do not compare numerical results across chapters directly") and line 54.

Original text:

> from assumptions about education outcomes. No direct cross-chapter result comparison is needed to motivate these extensions.

Replaced by:

> from assumptions about education outcomes.

### E3.tb.disc.7
Reason: cut: "These conclusions have several limitations." The heading says it.

Original text:

> These conclusions have several limitations. First, the primary estimator is an oracle.

Replaced by:

> First, the primary estimator is an oracle.

### E3.tb.disc.8
Reason: cut: "These would extend the completed application and are not needed to characterize the prespecified settings reported here." Defends the completeness of the reported work against no stated objection; ends at the last concrete point.

Original text:

> if that precision is needed for a final choice. These would extend the completed application and are not needed to characterize the prespecified settings reported here.

Replaced by:

> if that precision is needed for a final choice.

### E3.ta.A3.2
Reason: cut: "A gray bar joins designs that no pairwise Nemenyi comparison separated." It repeats the caption of Figure fig:cd almost word for word. SHARED: SI (M3.ta.2).

Original text:

> The red arrow shows the critical difference, 0.712 at $\alpha = 0.05$. A gray bar joins designs that no pairwise Nemenyi comparison separated. Under queen contiguity, the bars join

Replaced by:

> The red arrow shows the critical difference, 0.712 at $\alpha = 0.05$. Under queen contiguity, gray bars join

### E3.ta.A3.3
Reason: cut: "Its estimates are kept and flagged." Already said at the end of the previous text paragraph (line 960). Chapter only: in the SI the two passages sit in different sections.

Original text:

> is Checkerboard, for which $\tau$ is not identified. Its estimates are kept and flagged. The five values

Replaced by:

> is Checkerboard, for which $\tau$ is not identified. The five values

### E3.ta.A7.3
Reason: cut: "Two consequences of Equation \ref{eq:error-decomposition} explain the results for $\tau$ reported in the body." It repeats the section's opening sentence 20 lines later; "First/Second" are separated by an intervening paragraph and do not organize the text. Removes one `\ref` and one `$\tau$` (intended invariant change). SHARED: SI (M3.ta.9).

Original text:

> Two consequences of Equation \ref{eq:error-decomposition} explain the results for $\tau$ reported in the body.  First, for designs whose treated clusters have treated neighbors,

Replaced by:

> For designs whose treated clusters have treated neighbors,

### E3.ta.A7.5
Reason: cut: "The higher coverage at large $\tau$ reflects wider intervals around the wrong quantity, not a recovery of $\tau$." The same paragraph already says the coefficient estimates $\tau - \gamma$ at every $\tau$ and that, "with the bias held near $-\gamma$, wider intervals covered $\tau$ more often". Removes two `$\tau$` (intended invariant change). SHARED: SI (M3.ta.11).

Original text:

> moving the mean bias from $-0.650$ to $-0.667$. The higher coverage at large $\tau$ reflects wider intervals around the wrong quantity, not a recovery of $\tau$.

Replaced by:

> moving the mean bias from $-0.650$ to $-0.667$.

### E3.tb.res.7 (reworded)
Reason: "They do not replace the full grid benchmark" answers a claim nobody made; the limitation it guards (extreme settings only) is now stated directly. SHARED: SI (M3.tb.5).

Original text:

> These results at the extreme parameter settings support a qualified downside argument for saturation. They do not replace the full grid benchmark.

Replaced by:

> These results come from the extreme parameter settings only and give qualified support to saturation on downside risk.

### E3.tb.app.3 (reworded)
Reason: No one claims universal optimality; the equivalence caveat is the one the reader needs.

Original text:

> Small point differences do not establish equivalence or universal optimality.

Replaced by:

> Small point differences do not establish equivalence.

### E3.tb.disc.2 (reworded)
Reason: "complement ... rather than replace it" rebuts an unstated claim; the paragraph's actual point is the dependence on the spillover model. Same strength as CTJ line 246 ("A good allocation does not remove the need for an analysis that models spillover").

Original text:

> These are design-stage results, and they complement analysis-stage modeling rather than replace it. Every design was analyzed with a spatial lag model that included the spillover term.

Replaced by:

> These design-stage results assume an analysis that models spillover: every design was analyzed with a spatial lag model that included the spillover term.

## 2026-11-11 — Introduction regrouped (Background & Motivation, Aims & Contributions)

Andrew asked for the Chapter 2 structure: drop the "Extending Chapter 2" and "The Design-Stage Gap..." headings, drop the "(Chapter 1)" callback, and drop the sentence reconciling beta/tau notation with Chapter 2 (beta becomes the intervention effect in both chapters after the reviewer-requested reparameterization). The opening and design-stage-gap paragraphs were kept; the Chapter 2, aims and outline paragraphs were rewritten. The full introduction as it stood before:

> # Introduction
> 
> Cluster randomized trials (CRTs) randomize intact groups, such as clinics, schools, or geographic areas, rather than individuals. They are the standard design when an intervention is delivered at the group level or when individual randomization would create unacceptable contamination between arms [@donner_design_2000; @hayes_cluster_2017]. Standard CRT analysis relies on the no-interference component of the stable unit treatment value assumption (SUTVA): a cluster's outcome must not depend on which other clusters were assigned to treatment [@Rubin1986]. This assumption comes from the theory of individually randomized trials. It can fail once clusters are spatially arranged, as can the weaker assumption of partial interference, under which outcomes may depend on assignments within prespecified groups but not across them [@hudgens_toward_2008]. Tobler's first law of geography states that "everything is related to everything else, but near things are more related than distant things" [@tobler_computer_1970]. By this principle, an intervention delivered in one location can plausibly influence outcomes in adjacent locations through shared populations, service networks, information diffusion, or environmental pathways. Documented examples of spatial spillover come from public health and development economics. In a trial of insecticide-treated bed nets in northern Ghana, child mortality was lower in control villages close to treated villages than in control villages farther away [@binka_impact_1998]. A school-based deworming trial in Kenya found health and attendance benefits among untreated children near treated schools [@miguel_worms_2004]. When control clusters benefit from nearby treatment, an unadjusted comparison of treated and control clusters may understate the direct effect. How much contamination occurs can depend on the spatial arrangement of the two arms, not only on the number of clusters allocated to each.
> 
> ## The Design-Stage Gap in the Spatial Interference Literature
> 
> Methods for spatial interference in CRTs have concentrated on correction at the analysis stage (Chapter 1). Spatial autocorrelation statistics such as Moran's $I$ [@moran_interpretation_1948] detect residual spatial dependence once outcomes have been observed. Spatial regression models, such as the spatial lag model [@lesage_introduction_2009] estimated by maximum likelihood [@ord_estimation_1975] and implemented in software [@bivand_r_2022], allow an analyst to separate a direct treatment effect from a spillover effect after the fact. In their systematic review of spatial analysis methods for cluster randomized trials, Jarvis and colleagues placed trials that used spatial separation or buffers outside the primary review of statistical methods [@jarvis_spatial_2017]. The literature therefore offers much more guidance on analyzing spatial effects than on how the choice of treatment assignment itself, made before any data are collected, might reduce the problem.
> 
> Design-stage tools do exist, but they were largely developed for a different purpose. Restricted and covariate-constrained randomization address the baseline imbalance that a limited number of clusters can produce [@raab_balance_2001; @moulton_covariatebased_2004]. Crisp and colleagues noted that their extension, which selects clusters from a larger eligible set, also increases spatial dispersion among the selected clusters [@crisp_covariate-constrained_2023]. Their dispersion criterion concerns which clusters are selected, not how they are assigned to arms. Balancing a covariate across arms also does not by itself determine whether treated and control clusters are spatially interspersed in a way that induces contamination. Other work addresses the problem studied here more directly. Randomized saturation designs assign different fractions of units treated to different groups so that direct and spillover effects can be jointly identified [@baird_optimal_2018]. Another approach selects whole clusters for inclusion or exclusion in place of a "fried-egg" design that omits each cluster's outer area from analysis [@mccann_spatial_2018]. Recent theoretical work develops designs under spatial, network and cross-cluster interference and spatio-temporal heterogeneity [@leung_rateoptimal_2022; @cai_independentset_2023; @leung_clusterrandomized_2025; @watson_design_2025]. This chapter complements that work with a simulation comparison of concrete assignment strategies that a trial team could implement, under heterogeneous baseline incidence and a range of spatial dependence and spillover conditions.
> 
> ## Extending Chapter 2
> 
> Chapter 2 examined block-stratified allocation as a practical way to obtain reasonably good estimation while limiting poor-allocation risk in the settings it studied. This chapter evaluates both criteria across a broader set of designs under heterogeneous incidence. We retain Chapter 2's simple random assignment benchmark and strict alternating allocation, called Checkerboard here, while adding incidence balancing, targeting, saturation and isolation strategies. The grid comparison increases the number of clusters from 8--12 to 100, considers queen as well as rook neighbors, and includes heterogeneous baseline incidence in the outcome model. These changes extend the design question, so the results are not a direct numerical replication of Chapter 2. In Chapter 2, $\beta$ denotes the intervention effect. Here that effect is $\tau$, and $\beta$ denotes the coefficient of the incidence covariate. The observed-incidence application later in the chapter uses the North Carolina (NC) community college geography and a continuous education outcome, without assuming that historical sudden unexpected death (SUD) incidence predicts baseline education performance.
> 
> ## Aims and Outline
> 
> This chapter has three aims: (1) compare assignment designs by mean squared error (MSE), bias and confidence interval coverage, (2) assess how conclusions depend on incidence heterogeneity, spatial dependence, spillover regime and treatment-effect magnitude, and (3) assess both average accuracy and allocation-specific downside for planning a SUD-prevention education pilot on an observed North Carolina geography. The detailed grid results present six designs, and the appendix documents the two related variants. A main-text benchmark includes all eight strategies and simple random sampling (SRS). The completed application evaluates all nine strategies separately for each observed planning year. The application simulates outcomes and assignments and does not report outcomes from an implemented intervention.
> 
> The Methods set out the outcome model, the baseline incidence mechanisms used on the grid, the treatment assignment designs, and the estimator and performance criteria, and then describe the grid simulation. The Results compare the designs on the 100-cluster grid, overall and by incidence configuration, spatial dependence and spillover regime, and then compare each design with SRS. Later subsections report how the comparison changes with the treatment effect, across parameter settings and across allocations, and the last two give sensitivity analyses for rook contiguity and for a model without the spillover term. The application then evaluates the same designs on the observed NC geography, including participation and power. The Discussion draws out implications for trial design and for the SUD education pilot.
> 

## 2026-11-11 — BSS/Checkerboard link stated once; two citation sentences tightened

Reason: Andrew asked that the Checkerboard = Chapter 2 BSS link be made once,
at first mention in the Introduction, and not repeated in Methods, the design
table or the Discussion. The Binka and Jarvis sentences were reworded after a
source check (PubMed abstract for Binka 1998; full text, exclusion criterion 5,
for Jarvis 2017).

Removed from Methods (Treatment Assignment Designs paragraph):

> Checkerboard is the strict block-stratified allocation of Chapter 2, in which no two rook neighbors share an arm. Chapter 2's Block Stratified Sampling design also admitted relaxed allocations on its $3\times3$ grid, where the arms cannot be balanced. On its $2\times4$ and $3\times4$ grids, as on the $10\times10$ grid here, the strict allocations are the two checkerboard patterns. We use the two names interchangeably for this allocation and use Checkerboard in the rest of this chapter.

Design table row label, old: "Checkerboard (strict Block Stratified Sampling allocation in Chapter 2)".
Discussion, old: "Checkerboard, the strict allocation of Chapter 2, failed ...".

Introduction, old sentences:

> It compared simple random assignment with block stratified sampling, which in its strict form lets no two rook neighbors share a treatment arm, on small grids of 8 to 12 clusters.

> The strategies include Chapter 2's block stratified allocation, called Checkerboard here, along with incidence balancing, incidence targeting, saturation and isolation designs.

> In a trial of insecticide-treated bed nets in northern Ghana, child mortality was lower in control villages close to treated villages than in control villages farther away [@binka_impact_1998].

> In their systematic review of spatial analysis methods for cluster randomized trials, Jarvis and colleagues placed trials that used spatial separation or buffers outside the primary review of statistical methods [@jarvis_spatial_2017].

## 2026-11-15 — Round 5 carry-over, Results (Thesis_Chapter.qmd)

Carry-over list: bios-dissertation `prelim/full-draft/notes/human-writer-review/round5/ch3/carry-over-edits.md`, items 3-01 to 3-07.

### Results, Design Performance Under Queen Contiguity, paragraph 3 (test detail)
Reason: Moved in shortened form to Appendix A3; the full paragraph is kept here.

> The Friedman test rejected equal MSE ranks across the six designs and SRS under queen contiguity at $\beta = 1$ ($\chi^2 = 560.33$, $\text{df} = 6$, 160 blocks, $p < 2.2\times10^{-16}$). Nemenyi comparisons (critical difference 0.712) separated 17 of the 21 pairs. Incidence-Guided Saturation ranked significantly better than every other design and than SRS ($p = 0.011$). SRS and Balanced Quartiles, whose average ranks were 2.60 and 2.69, could not be distinguished from each other ($p = 1.00$), and both ranked significantly better than the remaining four designs. The other three pairs not separated all lay among the bottom three designs. Holm-adjusted Wilcoxon signed-rank tests separated 20 of the 21 pairs, all except Balanced Quartiles against SRS ($p = 0.15$). We also compared each pair of designs adjacent in the mean-MSE order over the 50 independent configuration-by-surface units. Every adjacent pair differed significantly (Wilcoxon $p \le 4.9\times10^{-9}$), and the better design had the lower MSE in all 50 units for every pair except High Incidence Focus against $2\times2$ Blocking (48 of 50). Appendix \ref{app:supp-comparisons} gives the pairwise results and the critical-difference diagrams.

### Results, Design Performance Under Queen Contiguity, rank-versus-mean paragraph
Reason: Shortened to two sentences; the numbers are in tab:performance and tab:grid-srs.

> The bottom three designs were ordered differently by average rank and by mean MSE: by mean MSE, High Incidence Focus (0.245), $2\times2$ Blocking (0.319), and Checkerboard (1.087). By average rank, the order was $2\times2$ Blocking (5.24), Checkerboard (5.74), and High Incidence Focus (5.76). Both reversals arise because the MSEs of $2\times2$ Blocking and Checkerboard depended strongly on the spillover regime. Each had a lower MSE than High Incidence Focus in all 80 both-arms blocks (mean MSE 0.062 for each, against 0.162) and a higher MSE in all 80 control-only blocks (0.576 and 2.111, against 0.329). Average rank counts each block once, so their wins under the both-arms regime offset their losses under the control-only regime. Mean MSE weights those losses by their size. The surface-level comparison, which averages each unit over both regimes, also separated High Incidence Focus from $2\times2$ Blocking.

### Results, Spillover Regime, Table tab:regime
Reason: Merged into tab:grid-srs (queen values identical). The rook columns are summarized in one sentence in Appendix A5.

> ```{=latex}
> \begin{table}[!htb]
>
> \small
> \centering
> \setlength{\tabcolsep}{4pt}
>
> \caption[Mean squared error by design under each spillover regime at $\beta = 1$]{Mean squared error by design under each spillover regime at $\beta = 1$, under queen contiguity (primary) and rook contiguity (sensitivity), pooled over the five incidence configurations. Each entry averages 80 scenarios (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\psi$). Rows are ordered by pooled mean MSE under queen contiguity (Table \ref{tab:performance}). Under rook contiguity, $\beta$ is not identifiable for Checkerboard. Its estimates are kept and flagged.}
> \label{tab:regime}
>
> \begin{tabular}{lcccc}
>
> \toprule
> & \multicolumn{2}{c}{Queen} & \multicolumn{2}{c}{Rook} \\
> Design & Both arms & Control only & Both arms & Control only \\
>
> \midrule
> Incidence-Guided Saturation & 0.046 & 0.135 & 0.045 & 0.098 \\
> Balanced Quartiles & 0.042 & 0.219 & 0.042 & 0.125 \\
> Isolation Buffer & 0.161 & 0.159 & 0.131 & 0.134 \\
> High Incidence Focus & 0.162 & 0.329 & 0.164 & 0.249 \\
> 2$\times$2 Blocking & 0.062 & 0.576 & 0.051 & 0.206 \\
> Checkerboard (flagged under rook) & 0.062 & 2.111 & 0.477 & 0.483 \\
>
> \bottomrule
>
> \end{tabular}
>
> \end{table}
> ```

### Results, Stability Across Treatment Effect Sizes (full text before shortening)
Reason: Shortened to about 300 words (carry-over 3-05). Block-level counts, worst cases, Friedman statistics across beta and the Balanced Quartiles gap are kept here; they are also in SI Section 6.

> ## Stability Across Treatment Effect Sizes {#sec-tau}
>
> Under queen contiguity, the order by average within-block rank and the order by mean MSE were each the same at every $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$. The two orders differ from each other in the bottom three designs (see @sec-queen; Figure \ref{fig:tau-mse}a shows mean MSE). With SRS included, the order by average rank placed SRS second at every value of $\beta$, and the Friedman statistic was 564.09, 560.33, 540.92, 526.55, and 497.55 at these five values ($p < 2.2\times10^{-16}$ at each). Under rook contiguity the order was also stable (Figure \ref{fig:tau-mse}b). By average rank among the six designs and SRS, it was the same at every value of $\beta$, with SRS third behind Incidence-Guided Saturation and Balanced Quartiles. By mean MSE, the close pair of $2\times2$ Blocking and Isolation Buffer traded places, with Isolation Buffer ahead only at $\beta = 3$. The five values of $\beta$ share their assignments and noise, so these five comparisons are not independent confirmations of one another. They show only that, under queen contiguity, the order did not change across the five values of $\beta$ we tested, given the same draws.
>
> [figure fig:tau-mse, unchanged]
> Within that stable order, the gap in MSE between Incidence-Guided Saturation and Balanced Quartiles under queen contiguity stayed between 0.037 and 0.041, and Incidence-Guided Saturation had the lower MSE of the two in 101 to 103 of the 160 blocks at every value of $\beta$. Isolation Buffer closed the gap instead: its MSE stayed at 0.160 at every value of $\beta$, while that of Incidence-Guided Saturation rose from 0.088 to 0.115 and that of Balanced Quartiles from 0.130 to 0.152. Isolation Buffer's gain came under the control-only regime, where the mean MSE of Incidence-Guided Saturation rose from 0.131 at $\beta = 0.8$ to 0.182 at $\beta = 3$, while that of Isolation Buffer stayed between 0.159 and 0.160. At $\beta = 2$ the two were level (0.160 for Isolation Buffer against 0.160), and at $\beta = 3$ Isolation Buffer had the lower control-only MSE (0.160 against 0.182). At $\beta = 2$ and $\beta = 3$, Isolation Buffer also had the lower worst case across all 160 blocks (0.207 against 0.210, and 0.208 against 0.238). Incidence-Guided Saturation remained second under the control-only regime at those two values, second to Balanced Quartiles under the both-arms regime at every value of $\beta$, and lowest in pooled MSE at every value of $\beta$. The rise in its MSE was mostly variance: its mean bias grew from 0.029 at $\beta = 0.8$ to 0.047 at $\beta = 3$, but squared bias stayed at or below 3.0\% of its MSE. Appendix \ref{app:technical-notes} traces both rises to a term of the estimation error that grows with $\beta$ when treated clusters have treated neighbors, and gives the block-level detail.
>
> Under rook contiguity, Incidence-Guided Saturation kept the lowest control-only MSE at every value of $\beta$, although its margin over Isolation Buffer narrowed, from 0.095 against 0.134 at $\beta = 0.8$ to 0.131 against 0.136 at $\beta = 3$.
>
> Except for Checkerboard under rook contiguity, coverage stayed close to 0.95 at every value of $\beta$: between 0.946 and 0.950 for every design under queen contiguity, and between 0.947 and 0.952 for every other design under rook contiguity. The coverage of rook Checkerboard, for which $\beta$ is not identifiable, rose from 0.150 at $\beta = 0.8$ to 0.464 at $\beta = 3$, because its intervals widened around the wrong quantity (Appendix \ref{app:technical-notes}). Appendix \ref{app:supp-comparisons} shows coverage across $\beta$ for each design.

## 2026-11-15 — Round 5 carry-over, Case Study (Thesis_Chapter.qmd)

Items 4-01 to 4-08. Moved (not cut): county conventions, numerator history and run sizes to A8; tab:nc-year to A8; CAR/Bayesian/spillover-form methods, the ICC translation and fig:power-n to A9.

### Case Study, Allocation-Specific Risk in North Carolina, tail-uncertainty paragraph (full)
Reason: Shortened to four sentences (carry-over 4-07).

> The tail estimates remain uncertain. With $J=100$, the worst decile contains only ten sampled draws. Extra outcomes reduce noise in $\hat m(a)$ but do not enlarge the allocation sample. All 112 reported settings in the tail refinement for designs with more than one possible assignment flag uncertain worst-decile membership. Under both arms, the Jaccard index of the worst deciles identified by the two outcome halves (allocations in both divided by allocations in either) was 0--0.25 for Incidence-Guided Saturation, and approximately 61--86\% of the raw variance of estimated risk across allocations was attributable to outcome Monte Carlo noise. Saturation Quadrants' worst-decile membership was similarly unstable. Small both-arms tail differences therefore do not establish superior downside protection. Checkerboard and High Incidence Focus each have only one possible assignment for these inputs. Their 90th percentile and worst-decile mean equal their conditional mean MSE, which does not imply low error or absence of outcome uncertainty.

### Case Study, closing paragraph of the power subsection
Reason: Removed: the Discussion's Design Guidance by Spillover Regime says the same (carry-over 4-08); its practical list is carried into that subsection.

> The case study supports considering either saturation strategy when control-only spillover is plausible, and SRS alongside other randomized strategies when both-arms exposure is expected. A trial team should consider likely spillover, recruitment, cluster and participant budgets, achievable measurement precision and interval calibration before selecting one realized allocation.

## 2026-11-15 — Round 5 carry-over, Discussion (Thesis_Chapter.qmd)

Items 5-01 to 5-05. The numbers dropped from the Summary are all stated in the Results.

### Summary of Findings (full, before shortening)

> We compared six treatment assignment designs in the detailed grid ranking, plus two related variants and SRS, for a spatial cluster randomized trial on a grid of 100 clusters. The comparison covered three mechanisms of heterogeneous incidence, two spillover regimes, and a range of values for spatial dependence, spillover magnitude, and the treatment effect. Some findings held across the incidence configurations and the values of $\rho$, $\psi$, and $\beta$ under queen contiguity. Others depended on the spillover regime, which a trial team may not know in advance.
>
> We judged each design on two criteria: its average MSE, and the spread of MSE across the allocations it can produce. A trial realizes one allocation, so a design with a low average MSE can still give a poor estimate if some of its allocations perform badly. Chapter 2 compared SRS and BSS on the same two criteria. For SRS under control-only spillover, the spread was large: on the grid, the worst 10\% of SRS allocations had mean conditional MSE 0.346, against 0.217 for SRS allocations overall.
>
> Under queen contiguity, Incidence-Guided Saturation had the lowest mean squared error (0.091) and Balanced Quartiles the second lowest (0.131), and the same two designs led under rook contiguity (0.072 and 0.084). Under queen contiguity, the order of the six designs was the same in each of the five incidence configurations and at every value of $\beta$, and neither the strength of spatial dependence in the outcome nor the spillover magnitude changed it. The mechanism that generated the incidence surface changed each design's MSE by under 10\%, except for High Incidence Focus. Coverage of the 95\% intervals was close to 0.95 for every design under both neighbor types, except Checkerboard under rook contiguity, so the designs differed in precision and coverage did not separate them.
>
> The spillover regime was the one factor that changed the queen ranking. Under queen contiguity at $\beta = 1$, Checkerboard and $2\times2$ Blocking were close to the best designs when spillover reached both arms and were the two worst when it reached only control clusters (MSE 0.062 and 2.111 for Checkerboard, 0.062 and 0.576 for $2\times2$ Blocking). Balanced Quartiles also lost ground under the control-only regime (0.042 and 0.219). Incidence-Guided Saturation, by contrast, was first or second among the six detailed designs under both regimes at every value of $\beta$, with a both-arms MSE 8\% above that of SRS. Its control-only MSE did rise with $\beta$, and Isolation Buffer, whose MSE did not, matched or passed it under that regime at $\beta = 2$ and passed it at $\beta = 3$. For a design whose performance depends on the regime, average MSE does not show the risk of choosing it when the regime is unknown.
>
> Checkerboard failed for different reasons under the two neighbor types. Under queen contiguity its estimates stayed centered on $\beta$ with close to nominal coverage but were imprecise, because under the control-only regime only the boundary clusters separate the spillover effect from the treatment effect. Under rook contiguity, $\beta$ is not identifiable, and its estimates were centered on $\beta - \psi$. Checkerboard intersperses treated and control clusters as finely as the grid allows, and both failures come from the dependence this pattern creates between $\mathbf{Z}$ and the spillover covariate (@sec-methods).
>
> These design-stage results assume an analysis that models spillover: every design was analyzed with a spatial lag model that included the spillover term. In the grid study, the estimates were close to unbiased for every design whose $\beta$ was identifiable: at $\beta = 1$, the absolute bias of each design under each regime was at most 0.048 under queen contiguity and 0.027 under rook contiguity. Given that model, the choice of design determined how precisely $\beta$ was estimated.

### Relation to prior Spatial CRT work (before rewrite)

> ### Relation to prior Spatial CRT work
>
> Chapter 2's practical argument considered reasonably good estimation together with protection against poor allocation in the small-grid settings it studied. We preserve both criteria while extending the designs, geography and incidence conditions. Our findings concern these new models and settings and do not compare numerical results across chapters directly. The rectangular Checkerboard and its graph adaptation illustrate why a restricted allocation's low assignment variation must be interpreted together with its level of estimation error. The case study also separates incidence-guided assignment from assumptions about education outcomes.

### Limitations, Fifth (before)

> Fifth, the simulation shares random draws by design. Blocks within an incidence configuration share its ten incidence surfaces, so the block-level tests treat dependent blocks as independent, which is why we also compared adjacent designs over the 50 independent configuration-by-surface units. The five values of $\beta$ also share their draws, so the stability of the design order across $\beta$ is not an independent replication at each value (see @sec-tau). Settings and years that share the same simulations cannot be treated as independent confirmation.

### Limitations, Sixth (before)

> Sixth, the simulation's base rate of 35 per 100,000 is not the rate of the county data described in @sec-application (see @sec-incidence). We kept it deliberately, because the simulation is a general design study rather than a calibration to one data set. In the Poisson mechanism, incidence enters the designs and the fitted model only through the ranks of the observed rates, so the base rate affects the simulation only through the spatial pattern of those ranks and the number of tied rates.

### Limitations, Eighth (before)

> Eighth, the 95\% Wald intervals under-covered in both the grid and the case study, and the Wald test was anti-conservative with few colleges. A small-sample correction, a $t$ reference distribution with a degrees-of-freedom-corrected standard error, restored coverage to between 94\% and 96\% for every design in the case study at every number of participating colleges (Appendix \ref{app:nc-ext}) and to between 0.946 and 0.952 for every identifiable design on the grid. We therefore report it in both studies. The finite-outcome allocation-tail estimates do not remove noise from individual risks or establish population worst cases. Budgets also differ, so these are evaluations of the retained strategies, not an equal-cost optimization.

### Future Considerations (before expansion)

> The next trial-design steps are to specify the education endpoint, eligible professionals, implementation budget and plausible spillover network, then assess interval calibration for that protocol before randomization. Additional methodological work could study noisy planning incidence, effect heterogeneity, exposure networks beyond contiguity and a broader set of graph geometries. The trial endpoint here is a measured education score, an intermediate outcome. Designing for a later reduction in SUD incidence would need a count outcome model and is left for future work. A larger allocation sample could make the downside estimates more precise if that precision is needed for a final choice.


## 2026-11-15 — Round 5 carry-over, Appendix holds (H-01 to H-12)

Held out of the appendix when it was reorganized into A1–A12. Each block is verbatim from the chapter at the commit before this step. All numbers appear in the SI and in the sentences named in the check table of the session report. Reasons are those of Section H of `bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/carry-over-edits.md`.

### H-01 `tab:seed-keys` (old A2, now A1)

**Reason:** Reproducibility detail; the SI and the code have it.

*Sentence removed from the Random-number streams paragraph: *

~~~~latex
Table \ref{tab:seed-keys} lists the three kinds of draw, the parameters in their keys, and the scenarios that therefore share them.
~~~~

*Table*

~~~~latex
```{=latex}
\begin{table}[!htb]

\small
\centering

\caption[Seed keys]{Seed keys. Each key names the kind of draw and the parameters listed here. The incidence mechanism is iid, spatial, or Poisson, and the iid configuration, which has no $\rho_X$, uses 0 in its place. $k$ indexes the incidence surface and $d$ the design. A draw is shared by every scenario whose key parameters match its own.}
\label{tab:seed-keys}

\begin{tabular}{p{4.2cm}p{5.0cm}p{5.2cm}}

\toprule
Draw & Key parameters & Shared by \\

\midrule
Incidence surface $\mathbf{X}_k$ & Incidence mechanism, $\rho_X$, $k$ & Every block, design, and value of $\beta$ of the configuration \\
The $J = 25$ assignments of design $d$ on surface $k$ & Incidence mechanism, $\rho_X$, neighbor type, $\rho$, $\psi$, spillover regime, $k$, $d$ & All five values of $\beta$ and both estimators \\
The 250 noise vectors of a block & Incidence mechanism, $\rho_X$, neighbor type, $\rho$, $\psi$, spillover regime & All eight designs, all five values of $\beta$, and both estimators \\

\bottomrule

\end{tabular}

\end{table}
```
~~~~

### H-02 `fig:cd-rook` and the rook half of the pairwise text (old A3, now A7)

**Reason:** The rook ranking is stated in A9; the queen diagram stays (single figure, label `fig:cd-queen`).

*Paragraph*

~~~~latex
**Pairwise tests under rook contiguity.** Nemenyi comparisons, with SRS as a seventh member, separated 17 of 21 pairs. The four they did not separate were Incidence-Guided Saturation against Balanced Quartiles ($p = 0.66$) and SRS ($p = 0.14$), Balanced Quartiles against SRS ($p = 0.97$), and $2\times2$ Blocking against Isolation Buffer ($p = 1.00$). Holm-adjusted Wilcoxon tests separated every pair except $2\times2$ Blocking against Isolation Buffer and Balanced Quartiles against SRS (both $p = 0.21$). In the surface-level comparison, Incidence-Guided Saturation had a lower MSE than Balanced Quartiles in 48 of the 50 units (mean difference $-0.012$; $-0.014$, $-0.010$; Wilcoxon $p = 9.9\times10^{-10}$). $2\times2$ Blocking and Isolation Buffer were close: $2\times2$ Blocking had the lower MSE in 36 of the 50 units (72\%; mean difference $-0.005$; $-0.008$, $-0.002$; Wilcoxon $p = 0.0027$). Isolation Buffer, however, had the slightly better average rank (4.21 against 4.33), and the two designs traded places in mean MSE at $\beta = 3$. Every other adjacent pair differed in all 50 units. Within single incidence configurations, rook coverage was between 0.944 and 0.954 for every design except Checkerboard, whose coverage was between 0.155 and 0.165.
~~~~

*Original two-panel figure (queen panel kept, rook panel held)*

~~~~latex
```{=latex}
\begin{figure}[!htb]

\centering
\begin{subfigure}{\linewidth}

\centering
\includegraphics[width=\linewidth]{figures/app_cd_diagram_rank_queen.pdf}

\caption{Queen contiguity (primary)}
\label{fig:cd-queen}
\end{subfigure}

\vspace{0.5em}
\begin{subfigure}{\linewidth}

\centering
\includegraphics[width=\linewidth]{figures/app_cd_diagram_rank_rook.pdf}

\caption{Rook contiguity (sensitivity)}
\label{fig:cd-rook}
\end{subfigure}

\caption[Nemenyi critical-difference diagrams for the six designs and simple random sampling (SRS) at $\beta = 1$]{Nemenyi critical-difference diagrams for the six designs and simple random sampling (SRS) at $\beta = 1$, pooled over the five incidence configurations. Designs are placed at their average within-block MSE rank, given in parentheses (1 = lowest MSE). Gray bars join designs that the Nemenyi comparisons did not separate.}
\label{fig:cd}

\end{figure}
~~~~

*Sentences removed from the figure description*

~~~~latex
Figure \ref{fig:cd} gives the Nemenyi critical-difference diagrams for the comparison of the six designs and SRS at $\beta = 1$, whose results are reported in @sec-queen and Appendix \ref{app:rook-figures}. Each design is placed at its average within-block MSE rank among the six designs and SRS over the 160 blocks of its neighbor type. The red arrow shows the critical difference, 0.712 at $\alpha = 0.05$. Under queen contiguity, gray bars join Balanced Quartiles with SRS, and $2\times2$ Blocking with Checkerboard and High Incidence Focus. Under rook contiguity, they join Incidence-Guided Saturation, Balanced Quartiles and SRS, and Isolation Buffer with $2\times2$ Blocking. The rook diagram includes Checkerboard, for which $\beta$ is not identifiable, with its estimates kept and flagged.
~~~~

### H-03 `fig:coverage-rook` (old A5, now A9)

**Reason:** Only Checkerboard differs, and its 0.162 is in the A9 text; `fig:mse-rook` stays.

*Lead-in sentence*

~~~~latex
In Figure \ref{fig:coverage-rook}, the rook version of Figure \ref{fig:coverage}, every design except Checkerboard has coverage close to the nominal level. Checkerboard's coverage, 0.162 pooled over the five incidence configurations (Table \ref{tab:performance}), is far below it.
~~~~

*Figure*

~~~~latex
```{=latex}
\begin{figure}[!htb]

\centering
\includegraphics[width=0.9\linewidth]{figures/fig_coverage_by_design_6design_rook.pdf}

\caption[Coverage of the 95\% confidence intervals by design at $\beta = 1$]{Coverage of the 95\% confidence intervals by design at $\beta = 1$, under rook contiguity, with the layout of Figure \ref{fig:coverage}. Checkerboard is flagged ($\beta$ not identifiable).}
\label{fig:coverage-rook}

\end{figure}
~~~~

### H-04 `fig:tau-coverage-queen` and `fig:tau-coverage-rook` (old A3, now A7)

**Reason:** Coverage is 0.946–0.952 at every β, stated in one sentence in the body (Results across β).

*Figure*

~~~~latex
```{=latex}
\begin{figure}[!htb]

\centering
\begin{subfigure}{\linewidth}

\centering
\includegraphics[width=\linewidth]{figures/app_tau_coverage_queen.pdf}

\caption{Queen contiguity (primary)}
\label{fig:tau-coverage-queen}
\end{subfigure}

\vspace{0.5em}
\begin{subfigure}{\linewidth}

\centering
\includegraphics[width=\linewidth]{figures/app_tau_coverage_rook.pdf}

\caption{Rook contiguity (sensitivity)}
\label{fig:tau-coverage-rook}
\end{subfigure}

\caption[Coverage of the 95\% confidence intervals by design across $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$]{Coverage of the 95\% confidence intervals by design across $\beta \in \{0.8, 1.0, 1.5, 2.0, 3.0\}$, pooled over the five incidence configurations. Each point averages the 160 scenarios of one design at one value of $\beta$. The dashed red line marks the nominal 0.95 level, and the gray dashed line marks the mean coverage of simple random sampling.}
\label{fig:tau-coverage}

\end{figure}
~~~~

*Paragraph*

~~~~latex
Across the five values of $\beta$, Figure \ref{fig:tau-coverage} gives the coverage of each design, as summarized in @sec-tau. The shaded bands show $\pm 1$ times the average of the scenario-level Monte Carlo standard errors of coverage. They describe the precision of a single scenario's coverage, not of the plotted average. The legend lists the designs in order of pooled mean MSE at $\beta = 1$ for the neighbor type of the panel (Table \ref{tab:performance}). The vertical axis runs from 0 to 1, so lines close to 0.95 overlap: all six under queen contiguity, between 0.946 and 0.950, and all but Checkerboard under rook contiguity, between 0.947 and 0.952. The lowest line under rook contiguity, which the legend does not mark, is Checkerboard, for which $\beta$ is not identifiable. The five values of $\beta$ share their assignments and noise.
~~~~

### H-05 Estimation-error subsection of old A7 (now A3, replaced by the shortened text)

**Reason:** Long derivation for a secondary result; the full version stays in SI Section 6. Kept in A3: the one-line error decomposition, the downward bias of ρ̂ (about 0.06), g for Incidence-Guided Saturation (0.57) and Isolation Buffer (0), the variance of the last term, and the rook Checkerboard coverage (0.150 to 0.464).

*Old A7 introduction (after the applabel)*

~~~~latex
\applabel{A7}{app:technical-notes}This section explains why the bias of Incidence-Guided Saturation and the coverage of Checkerboard under rook contiguity both grow with $\beta$, as reported in the body of this chapter. It closes with block-level detail for the results across $\beta$. The stored results keep neither $\hat\rho$ nor the standard error of each fit, so we re-ran every block for the six designs, both neighbor types, and all five values of $\beta$, 2.4 million fits of the oracle model, recording $\hat\rho$ and the standard error. The re-run reproduced the stored results in all 9,600 scenarios. Numbers below that involve $\hat\rho$, the standard error, or the terms of Equation \ref{eq:error-decomposition} come from this re-run. The others come from the stored results.
~~~~

*Old subsection 'Estimation Error and Its Dependence on β', including the block-level detail paragraph*

~~~~latex
### Estimation Error and Its Dependence on $\beta$ {.unnumbered .unlisted}

For one fit $(k, j)$, write $\mathbf{Y}$, $\mathbf{Z}$, and $\boldsymbol{\varepsilon}$ for $\mathbf{Y}_{kj}$, $\mathbf{Z}_{kj}$, and $\boldsymbol{\varepsilon}_{kj}$ (Simulation design, resampling, and metrics). Let $\mathbf{V} = [\mathbf{1}, \mathbf{Z}, \mathbf{S}(\mathbf{Z}), \mathbf{X}_k]$ be the $100 \times 4$ covariate matrix of the oracle model, so that Equation \ref{eq:outcome} reads $(\mathbf{I} - \rho \mathbf{W})\mathbf{Y} = \mathbf{V}\boldsymbol{\theta} + \boldsymbol{\varepsilon}$ with $\boldsymbol{\theta} = (0, \beta, \psi, \tau)^\top$. Let $\mathbf{A} = (\mathbf{I} - \rho \mathbf{W})^{-1}$ and $\mathbf{V}^{+} = (\mathbf{V}^\top \mathbf{V})^{-1}\mathbf{V}^\top$, and for a vector $\mathbf{v}$ write $[\mathbf{V}^{+}\mathbf{v}]_Z$ for the coefficient on $\mathbf{Z}$ when $\mathbf{v}$ is regressed on the columns of $\mathbf{V}$ by least squares. Given $\hat\rho$, the estimated coefficients are $\mathbf{V}^{+}(\mathbf{I} - \hat\rho \mathbf{W})\mathbf{Y}$ (@sec-estimation). Because $(\mathbf{I} - \hat\rho \mathbf{W})\mathbf{Y} = \mathbf{V}\boldsymbol{\theta} + \boldsymbol{\varepsilon} + (\rho - \hat\rho)\mathbf{W}\mathbf{Y}$ and $\mathbf{W}\mathbf{Y} = \mathbf{W}\mathbf{A}(\mathbf{V}\boldsymbol{\theta} + \boldsymbol{\varepsilon})$,
\begin{equation}
\hat\beta - \beta = [\mathbf{V}^{+}\boldsymbol{\varepsilon}]_Z + (\rho - \hat\rho)\,h + \beta(\rho - \hat\rho)\,g,
\label{eq:error-decomposition}
\end{equation}
where
$$
g = [\mathbf{V}^{+}\mathbf{W}\mathbf{A}\mathbf{Z}]_Z, \qquad h = [\mathbf{V}^{+}\mathbf{W}\mathbf{A}\{\psi\,\mathbf{S}(\mathbf{Z}) + \tau \mathbf{X}_k + \boldsymbol{\varepsilon}\}]_Z.
$$
The first term has mean zero given $\mathbf{Z}$ and $\mathbf{X}_k$. The last term has a factor of $\beta$. It is present exactly when $g \neq 0$, where $g$ is the coefficient on $\mathbf{Z}$ when the lagged treatment $\mathbf{W}\mathbf{A}\mathbf{Z} = \mathbf{W}\mathbf{Z} + \rho \mathbf{W}^2 \mathbf{Z} + \rho^2 \mathbf{W}^3 \mathbf{Z} + \cdots$ is regressed on the model's own covariates. The sign and size of the last term are set by the error in $\hat\rho$. The concentrated likelihood depends on $\mathbf{Y}$ only through the residual sum of squares of $(\mathbf{I} - \tilde\rho \mathbf{W})\mathbf{Y}$ on $\mathbf{V}$, and the part of that residual that involves $\beta$ is $\beta(\rho - \tilde\rho)$ times the residual of $\mathbf{W}\mathbf{A}\mathbf{Z}$ on $\mathbf{V}$. If $\mathbf{W}\mathbf{A}\mathbf{Z}$ lies in the column space of $\mathbf{V}$, that part vanishes for every $\tilde\rho$, so $\hat\rho$ is the same at every value of $\beta$. If, in addition, $g = 0$, then $\hat\beta - \beta$ is the same at every value of $\beta$. Equation \ref{eq:error-decomposition} held on every re-run fit, with a largest residual of $4.2\times10^{-13}$. For Checkerboard under rook contiguity, it held in the reduced form given below.

For designs whose treated clusters have treated neighbors, the last term gives a small positive bias that grows with $\beta$. Under queen contiguity, pooled over the 160 blocks, the bias of Incidence-Guided Saturation rose from 0.029 at $\beta = 0.8$ to 0.047 at $\beta = 3$ (@sec-tau), while that of Isolation Buffer stayed between 0.002 and 0.004. Almost all of the rise comes from the last term of Equation \ref{eq:error-decomposition}. In our fits, $\hat\rho$ was biased downward by about 0.06 on average (for Incidence-Guided Saturation, by 0.071 at $\beta = 0.8$ and 0.054 at $\beta = 3$). We observed this downward bias in these simulations but did not derive it or check it against a published finite-sample result. With $\rho - \hat\rho$ positive on average, the last term is positive and grows with $\beta$ for any design with $g > 0$. For Incidence-Guided Saturation, $g$ averaged 0.371 over the queen fits, and the last term contributed 0.021 of the bias at $\beta = 0.8$ and 0.043 at $\beta = 3$. The first term contributed 0.001 at both values, and the second 0.007 and 0.002. For Isolation Buffer, $g$ averaged 0.020, and the last term contributed 0.001 and 0.004.

Most of the difference in $g$ comes from the control-only regime. Under the both-arms regime with $\rho = 0$, the lagged treatment is $\mathbf{W}\mathbf{Z} = \mathbf{S}(\mathbf{Z})$, one of the covariates, so $g = 0$ for every design. Under the control-only regime, $\mathbf{W}\mathbf{Z} = \mathbf{S}(\mathbf{Z}) + \mathbf{Z} \circ \mathbf{W}\mathbf{Z}$, and at $\rho = 0$, $g$ is close to the average share of treated neighbors among treated clusters. That share is 0.57 for Incidence-Guided Saturation, whose treated clusters are concentrated in the high-saturation quadrants, and exactly 0 for Isolation Buffer, whose treated clusters have no treated neighbors. For Isolation Buffer, $\mathbf{Z} \circ \mathbf{W}\mathbf{Z} = \mathbf{0}$, so the two regimes define the same spillover covariate and the same $g$, which is 0 at $\rho = 0$ and about 0.06 at $\rho = 0.5$. At $\rho = 0$, its $\hat\rho$ and its estimation error were the same at every value of $\beta$, up to the optimizer's tolerance.

The same term explains why the control-only MSE of Incidence-Guided Saturation rose from 0.131 at $\beta = 0.8$ to 0.182 at $\beta = 3$ while that of Isolation Buffer did not. Under that regime, averaged over the 80 queen blocks, the variance of the last term of Equation \ref{eq:error-decomposition} rose from 0.008 to 0.061 for Incidence-Guided Saturation, while the variance of the rest of the error changed little (0.119 and 0.115). The variance of the last term grew more slowly than $\beta^2$ because $\hat\rho$ became less variable as $\beta$ grew. For Isolation Buffer, the variance of the last term stayed below 0.0002. The rise in bias slowed at larger values of $\beta$ for a related reason: the part of the likelihood that involves $\beta$ sharpens the likelihood around the true $\rho$ as $\beta$ grows, so the downward bias of $\hat\rho$ shrank.

The decomposition also explains why the coverage of Checkerboard under rook contiguity rose with $\beta$, from 0.150 at $\beta = 0.8$ to 0.464 at $\beta = 3$ (@sec-tau), although its coefficient on $\mathbf{Z}$ estimates $\beta - \psi$ at every value of $\beta$. For this design $\mathbf{W}\mathbf{Z} = \mathbf{1} - \mathbf{Z}$, and because $\mathbf{W}$ is row-standardized, $\mathbf{W}^2 \mathbf{Z} = \mathbf{W}(\mathbf{1} - \mathbf{Z}) = \mathbf{Z}$. The powers of $\mathbf{W}$ applied to $\mathbf{Z}$ therefore alternate between $\mathbf{1} - \mathbf{Z}$ and $\mathbf{Z}$, and $\mathbf{W}\mathbf{A}\mathbf{Z} = \{(\mathbf{1} - \mathbf{Z}) + \rho \mathbf{Z}\}/(1 - \rho^2)$ lies in the span of the intercept and $\mathbf{Z}$. The derivation of Equation \ref{eq:error-decomposition} applies with the aliased spillover column dropped, $\mathbf{V} = [\mathbf{1}, \mathbf{Z}, \mathbf{X}_k]$ and $\boldsymbol{\theta} = (\psi, \beta - \psi, \tau)^\top$, so that $\beta - \psi$ takes the place of $\beta$, and $g = -1/(1 + \rho)$. Two consequences follow: $\hat\rho$ is the same at every value of $\beta$, and
$$
\hat\beta - \beta = -\psi + \delta + (\beta - \psi)\,\frac{\hat\rho - \rho}{1 + \rho},
$$
where $\delta$ collects the first two terms of that decomposition and, like $\hat\rho$, does not depend on $\beta$. The error in $\hat\rho$ is therefore multiplied by $\beta - \psi$, which grows with $\beta$. The mean standard deviation of $\hat\beta$ rose from 0.202 at $\beta = 0.8$ to 0.327 at $\beta = 3$, and the standard error tracked it, rising on average from 0.199 to 0.322. With the bias held near $-\psi$, wider intervals covered $\beta$ more often. The rise was largest at $\rho = 0$ (coverage 0.157 to 0.562) and smallest at $\rho = 0.5$ (0.142 to 0.284), consistent with the factor $1/(1 + \rho)$. The same term raised the MSE from 0.477 to 0.567. About three quarters of that increase was variance, $(n-1)\text{SD}^2/n$ averaged over scenarios (Appendix \ref{app:mechanics}), which rose from 0.041 to 0.109. The other quarter was squared bias, which rose from 0.436 to 0.458 because $\hat\rho$ was biased slightly downward, moving the mean bias from $-0.650$ to $-0.667$.

**Block-level detail.** Under queen contiguity, the average rank of Incidence-Guided Saturation among the six designs rose from 1.41 at $\beta = 0.8$ to 1.85 at $\beta = 3$, while that of Balanced Quartiles stayed between 2.13 and 2.14 and that of Isolation Buffer fell from 3.71 to 3.31. Incidence-Guided Saturation had a lower MSE than Balanced Quartiles in 103, 102, 102, 102, and 101 of the 160 blocks at $\beta = 0.8$, 1, 1.5, 2, and 3, and a lower MSE than Isolation Buffer in 75, 73, 56, 40, and 12 of the 80 control-only blocks. Its worst case across blocks was the lowest of the six designs at $\beta \le 1.5$, but at $\beta = 2$ and $\beta = 3$ that of Isolation Buffer was lower (0.207 against 0.210, and 0.208 against 0.238). Under rook contiguity, at $\beta = 3$, Incidence-Guided Saturation had a lower MSE than Isolation Buffer in 47 of the 80 control-only blocks, and its worst case across rook blocks was the lowest of the six designs at every value of $\beta$ except $\beta = 3$, where that of Isolation Buffer was lower (0.173 against 0.174). The coverage of rook Checkerboard was 0.150, 0.162, 0.212, 0.286, and 0.464 at $\beta = 0.8$, 1, 1.5, 2, and 3.
~~~~

### H-06 `tab:nc-summary-sens` and `fig:nc-setting-ratios` (old A8)

**Reason:** Overlap with `fig:nc-mse` and `tab:nc-annual`.

*Figure*

~~~~latex
\begin{figure}[!htb]
\centering
\includegraphics[width=\linewidth]{figures/real_sud/primary_setting_ratios.pdf}

\caption[Primary queen setting-specific mean MSE ratios to year-matched SRS]{Primary queen setting-specific mean MSE ratios to year-matched SRS. Teal cells have lower MSE than SRS and red cells higher. Setting-specific Monte Carlo errors are kept in the aggregate output tables and not shown here. Parameters and spillover regimes are shown separately, and neighbor definitions are not pooled.}
\label{fig:nc-setting-ratios}

\end{figure}
~~~~

*Paragraph*

~~~~latex
For Incidence-Guided Saturation, Table \ref{tab:nc-summary-sens} reports the matched annual mean-MSE ratios. Rook means use only the four extreme combinations of $\rho$ and $\psi$ shared with the primary analysis. Means for the matched sensitivity analysis in which incidence shifts the outcome and for the regional summaries use the six queen pairs. These are matched contrasts, not ratios of unmatched pooled analyses.
~~~~

*Table*

~~~~latex
```{=latex}
\begin{table}[!htb]
\small
\centering

\caption[Incidence-Guided Saturation sensitivity MSE divided by matched primary MSE]{Incidence-Guided Saturation sensitivity MSE divided by matched primary MSE, by year and regime. Incidence shift adds rank-scaled $\mathbf{X}$ to both generation and fit. Rook uses the matched extreme combinations of $\rho$ and $\psi$. Regional-summary alternatives preserve fixed regions.}
\label{tab:nc-summary-sens}
% SOURCE: application/results/real_sud_rev_20261002/exhibits/yearly_sensitivity_means.csv

\begin{tabular}{llrrrr}
\toprule
Sensitivity & Regime & 2018 & 2019 & 2020 & 2021 \\
\midrule
Incidence shift / adjusted & Both & 1.016 & 1.030 & 1.015 & 1.035 \\
Incidence shift / adjusted & Control & 1.026 & 1.023 & 1.042 & 1.040 \\
Rook contiguity & Both & 1.005 & 0.998 & 0.985 & 0.998 \\
Rook contiguity & Control & 0.967 & 0.969 & 0.971 & 1.036 \\
Mean regional rate & Both & 1.034 & 1.000 & 1.000 & 1.031 \\
Mean regional rate & Control & 0.907 & 1.000 & 1.000 & 0.928 \\
Pooled regional rate & Both & 1.015 & 0.979 & 0.981 & 1.003 \\
Pooled regional rate & Control & 0.962 & 1.133 & 1.061 & 1.030 \\

\bottomrule
\end{tabular}
\end{table}
```
~~~~

### H-07 `tab:nc-bayes` (old A9)

**Reason:** The body and A5/A12 say the Bayesian lag model reproduced the ML results and the Bayesian CAR model did not.

*Table*

~~~~latex
\begin{table}[!htb]
\small
\centering
\setlength{\tabcolsep}{4pt}
\caption[Mean MSE relative to simple random sampling (SRS) at $\rho = 0.5$]{Mean MSE relative to simple random sampling (SRS) at $\rho = 0.5$ for maximum likelihood (ML), the Bayesian spatial lag model (Bayes SAR) and the Bayesian CAR model with converged chains only (Bayes CAR), all fitted to the same spatial lag data. 58 colleges, $\tau = 0$, $\psi = 0.8$, mean over the four planning years.}
\label{tab:nc-bayes}
% SOURCE: application/results/methods_ext_2026-10/sensitivity_summaries_v2/a6_bayes_lag_ratio_to_srs.csv (code/24_methods_extension_tables.R)
\begin{tabular}{lcccccc}
\toprule
& \multicolumn{3}{c}{Both arms} & \multicolumn{3}{c}{Control only} \\
\cmidrule(lr){2-4}\cmidrule(lr){5-7}
Design & ML & Bayes SAR & Bayes CAR & ML & Bayes SAR & Bayes CAR \\
\midrule
Checkerboard & 3.14 & 3.18 & 2.90 & 3.16 & 3.16 & 2.31 \\
High Incidence Focus & 0.96 & 0.96 & 0.83 & 1.19 & 1.19 & 1.04 \\
Saturation Quadrants & 1.12 & 1.10 & 1.11 & 0.86 & 0.86 & 1.05 \\
Isolation Buffer & 3.89 & 3.94 & 3.80 & 1.01 & 1.02 & 0.87 \\
$2\times2$ Blocking & 1.14 & 1.14 & 1.11 & 1.48 & 1.48 & 1.34 \\
Balanced Quartiles & 1.00 & 1.00 & 1.05 & 1.10 & 1.10 & 1.15 \\
Balanced Halves & 0.97 & 0.97 & 0.99 & 1.05 & 1.06 & 1.01 \\
Incidence-Guided Saturation & 1.12 & 1.10 & 1.13 & 0.81 & 0.81 & 0.95 \\
Simple Random Sampling & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 \\
\bottomrule
\end{tabular}
\end{table}
~~~~

*Sentence removed from the convergence paragraph (shared with H-08)*

~~~~latex
Tables \ref{tab:nc-car} and \ref{tab:nc-bayes} give each strategy's MSE relative to SRS under each model.
~~~~

### H-08 `tab:nc-car` (old A9)

**Reason:** Same reason; one sentence of its result is kept in A12.

*Table*

~~~~latex
\begin{table}[!htb]
\small
\centering
\caption[Mean MSE relative to simple random sampling (SRS) under four combinations of data-generating model and fitted model]{Mean MSE relative to simple random sampling (SRS) under four combinations of data-generating model and fitted model: spatial lag (SAR) or conditional autoregressive (CAR). 58 colleges, $\tau = 0$, mean over 24 blocks (year $\times$ $\rho$ $\times$ $\psi$). Below 1 is lower error than SRS.}
\label{tab:nc-car}
% SOURCE: application/results/methods_ext_2026-10/sensitivity_summaries/a6_car_ratio_to_srs.csv (code/24_methods_extension_tables.R)
\begin{tabular}{lcccccccc}
\toprule
& \multicolumn{4}{c}{Both arms} & \multicolumn{4}{c}{Control only} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-9}
Data & SAR & SAR & CAR & CAR & SAR & SAR & CAR & CAR \\
Fit & SAR & CAR & SAR & CAR & SAR & CAR & SAR & CAR \\
\midrule
Checkerboard & 2.88 & 2.83 & 2.93 & 3.09 & 2.99 & 2.64 & 3.01 & 3.00 \\
High Incidence Focus & 1.06 & 1.02 & 1.07 & 1.07 & 1.08 & 1.03 & 1.10 & 1.09 \\
Saturation Quadrants & 1.04 & 1.02 & 1.04 & 1.03 & 0.79 & 0.81 & 0.81 & 0.78 \\
Isolation Buffer & 4.11 & 4.13 & 4.16 & 4.45 & 1.03 & 0.95 & 1.03 & 1.05 \\
$2\times2$ Blocking & 0.98 & 0.99 & 1.00 & 1.02 & 1.21 & 1.17 & 1.22 & 1.22 \\
Balanced Quartiles & 1.04 & 1.05 & 1.04 & 1.06 & 1.03 & 1.00 & 1.03 & 1.02 \\
Balanced Halves & 1.04 & 1.03 & 1.05 & 1.04 & 1.02 & 0.98 & 1.03 & 0.99 \\
Incidence-Guided Saturation & 1.03 & 1.04 & 1.04 & 1.03 & 0.81 & 0.82 & 0.81 & 0.76 \\
Simple Random Sampling & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 \\
\bottomrule
\end{tabular}
\end{table}
~~~~

### H-09 `tab:rho` (old A4; `tab:regime` was already removed)

**Reason:** ρ changed no MSE by more than 0.020; stated in the body (Sensitivity subsection).

*Table*

~~~~latex
```{=latex}
\begin{table}[!htb]

\small
\centering

\caption[Mean squared error by design across the four levels of outcome spatial dependence $\rho$]{Mean squared error by design across the four levels of outcome spatial dependence $\rho$, under queen contiguity at $\beta = 1$, pooled over the five incidence configurations. Each entry averages 40 scenarios (5 incidence configurations $\times$ 4 values of $\psi$ $\times$ 2 spillover regimes). Rows are ordered by pooled mean MSE under queen contiguity (Table \ref{tab:performance}).}
\label{tab:rho}

\begin{tabular}{lcccc}

\toprule
Design & $\rho=0$ & $\rho=0.01$ & $\rho=0.20$ & $\rho=0.50$ \\

\midrule
Incidence-Guided Saturation & 0.090 & 0.087 & 0.093 & 0.093 \\
Balanced Quartiles & 0.129 & 0.129 & 0.132 & 0.133 \\
Isolation Buffer & 0.161 & 0.157 & 0.158 & 0.165 \\
High Incidence Focus & 0.243 & 0.242 & 0.249 & 0.248 \\
2$\times$2 Blocking & 0.318 & 0.320 & 0.314 & 0.324 \\
Checkerboard & 1.095 & 1.077 & 1.084 & 1.092 \\

\bottomrule

\end{tabular}

\end{table}
~~~~

### H-10 `fig:spread-gain` (old A9)

**Reason:** Its point (gain rises with the spread of regional incidence, correlation 0.84 to 0.99) is one sentence in the body.

*Figure*

~~~~latex
\begin{figure}[!htb]
\centering
\includegraphics[width=0.9\linewidth]{figures/fig_spread_vs_gain.pdf}
\caption[Gain in treated incidence share from incidence guidance against the spread of the regions' mean incidence ranks]{Gain in treated incidence share from incidence guidance (Incidence-Guided Saturation minus Saturation Quadrants) against the spread of the regions' mean incidence ranks, the largest minus the smallest regional mean rank divided by the number of clusters. Filled symbols are the 50 grid surfaces (queen, production allocations) by incidence configuration, with a least-squares line fitted to them. Open triangles are the fixed NC regions in each of the four planning years, whose spread divides by 58 clusters rather than 100.}
\label{fig:spread-gain}
\end{figure}
~~~~

*Sentence removed from the Targeting paragraph*

~~~~latex
Figure \ref{fig:spread-gain} shows that the gain rose almost linearly with the spread of the regions' mean incidence ranks, and that the NC regions lie at the low end of that spread.
~~~~

### H-11 `fig:coverage` (old A5)

**Reason:** Every design is at 0.946–0.950, and `tab:performance` already has a coverage column.

*Figure*

~~~~latex
\begin{figure}[!htb]

\centering
\includegraphics[width=0.9\linewidth]{figures/fig_coverage_by_design_6design_queen.pdf}

\caption[Coverage of the 95\% confidence intervals by design at $\beta = 1$]{Coverage of the 95\% confidence intervals by design at $\beta = 1$, under queen contiguity. Each box summarizes the coverage of the 160 scenarios of one design (5 incidence configurations $\times$ 4 values of $\rho$ $\times$ 4 values of $\psi$ $\times$ 2 spillover regimes), each estimated from 250 fits. Designs are ordered by median coverage, not by MSE. The dashed red line marks the nominal 0.95 level, and the gray dotted line the median coverage of simple random sampling. Figure \ref{fig:coverage-rook} shows rook contiguity.}
\label{fig:coverage}

\end{figure}
~~~~

### H-12 `fig:biasvar` and `fig:biasvar-rook` (old A5), combined into one two-panel figure (not held; originals kept for the record)

**Reason:** One figure instead of two; the same layout and configuration.

*Original `fig:biasvar`*

~~~~latex
```{=latex}
\begin{figure}[!htb]

\centering
\includegraphics[width=0.9\linewidth]{figures/fig_biasvar_6design_queen.pdf}

\caption[Decomposition of mean squared error into squared bias and variance by design at $\beta = 1$ for one incidence configuration]{Decomposition of mean squared error into squared bias and variance by design at $\beta = 1$ for one incidence configuration, Poisson with $\rho_X = 0.20$, under queen contiguity. For each design, the two stacked components are the squared bias and the variance ($\text{SD}^2$) of $\hat\beta$, each averaged over the 32 scenarios of that configuration (4 values of $\rho$ $\times$ 4 values of $\psi$ $\times$ 2 spillover regimes). Designs are ordered by pooled mean MSE (Table \ref{tab:performance}). The gray dashed line marks the MSE of simple random sampling (SRS) in the same configuration. Figure \ref{fig:biasvar-rook} shows rook contiguity.}
\label{fig:biasvar}

\end{figure}
```
~~~~

*Original `fig:biasvar-rook` and its lead-in sentence*

~~~~latex
Nearly all of Checkerboard's MSE in Figure \ref{fig:biasvar-rook}, the rook version of Figure \ref{fig:biasvar}, is squared bias, because its coefficient on $\mathbf{Z}$ estimates $\beta - \psi$ (@sec-methods).

```{=latex}
\begin{figure}[!htb]

\centering
\includegraphics[width=0.9\linewidth]{figures/fig_biasvar_6design_rook.pdf}

\caption[Squared bias and variance of $\hat\beta$ by design at $\beta = 1$ for the Poisson configuration with $\rho_X = 0.20$]{Squared bias and variance of $\hat\beta$ by design at $\beta = 1$ for the Poisson configuration with $\rho_X = 0.20$, under rook contiguity, with the layout of Figure \ref{fig:biasvar}. Checkerboard is flagged ($\beta$ not identifiable).}
\label{fig:biasvar-rook}

\end{figure}
~~~~

## 2026-11-15 (round 5 writing review)

Sentences removed from `Thesis_Chapter.qmd` when the approved round 5 writing edits (W-01 to W-117, `bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch3/proposed-edits.md`) were applied. Each entry gives the item, the former location (section), the verbatim text and the reason. Sentences that were only reworded are not listed.

*W-83. Appendix A5 (opening paragraph). Reason: describes file handling, not the study.*

~~~~
These runs left the production case-study files unchanged.
~~~~

*W-112. Methods, Estimation and Performance Criteria. Reason: no outcome of the check was reported, and the sentence is not critical.*

~~~~
As a further check, we refit a 1\% sample of the full simulation with `lagsarlm` itself.
~~~~

*W-59. Results, Case Study, Accuracy, Coverage and Share of Population Treated (both-arms paragraph). Reason: repeats the caveat given in the grid results (W-65).*

~~~~
Small point differences do not establish equivalence.
~~~~

*W-59. Same paragraph. Replaced by "...does not show that the design has the protection that randomization provides." Reason: "randomization-based causal validity" is jargon.*

~~~~
Its low MSE is model-based accuracy conditional on that allocation, not evidence of randomization-based causal validity.
~~~~

*W-65. Results, Sampling Designs Compared to Simple Random Sampling. Merged into one sentence with "although" and a rewritten last sentence. Reason: a stand-alone caveat sentence, and "hides" gave the pooled average an action.*

~~~~
Small mean differences do not establish equivalence. Pooling over the regimes also hides how the design guidance differs between them.
~~~~

*W-16. Case Study, Accuracy paragraph. Reason: duplicate figure reference (the figure is cited in the first sentence of the paragraph).*

~~~~
(Figure \ref{fig:nc-mse}; Table \ref{tab:nc-year} in Appendix \ref{app:nc-results})
~~~~

*W-99. Appendix A4 (opening paragraph). Reason: repeated from the body (Case Study, Observed Inputs).*

~~~~
The county death-certificate source data are restricted and not publicly available. We report only 58-cluster summaries.
~~~~

*W-99. Appendix A4, Data Conventions. Reason: repeated from the body; the remaining sentence points to the body for the corrected numerator.*

~~~~
The corrected numerator supplied by Habib includes heart-failure deaths formerly excluded as presumed non-sudden. It totals 23,523 deaths over 25,594,321 person-years. The earlier numerator, which excluded them, reproduces the totals behind Habib's unpublished report [@habib_temporal_2026]. The case study uses the corrected numerator.
~~~~

*W-100. Appendix A6 (variants paragraph). Reason: self-reference to its own appendix.*

~~~~
(Appendix \ref{app:eight-design})
~~~~

*W-100. Appendix A6 (variants paragraph). Reason: repeats the Discussion (Design Guidance by Spillover Regime).*

~~~~
The accuracy advantage of both designs over SRS, which arises under control-only spillover (see @sec-vs-srs), comes from the saturation structure they share.
~~~~

*W-88. Appendix A9 (rook results). Reason: repeats the sentence that squared bias accounted for 0.438 of Checkerboard's rook MSE of 0.480, which already cites the figure.*

~~~~
Nearly all of Checkerboard's MSE under rook contiguity (Figure \ref{fig:biasvar-rook}) is squared bias, because its coefficient on $\mathbf{Z}$ estimates $\beta - \psi$ (@sec-methods).
~~~~

*W-88. Appendix A9 (before the rook MSE figure). Reason: repeats the opening paragraph of A9, which says Figure \ref{fig:mse-rook} gives the rook version of Figure \ref{fig:mse}.*

~~~~
Figure \ref{fig:mse-rook} corresponds to Figure \ref{fig:mse}.
~~~~

*W-29. Appendix A11 (opening text before Table `tab:nc-annual`). Replaced by a statement that the records were checked. Reason: file-management language.*

~~~~
More granular setting-specific results, precision components, residual treatment information and warning and failure records are kept as aggregate output tables. These records were checked, so the absence of warnings and failures is a verified result.
~~~~

*W-82. Appendix A11, caption of `fig:nc-risk`. Replaced by the statement of what the tail estimates are. Reason: "a fresh allocation confirmation" was unclear.*

~~~~
They are not quantiles of realized errors or a fresh allocation confirmation.
~~~~
