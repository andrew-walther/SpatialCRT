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
