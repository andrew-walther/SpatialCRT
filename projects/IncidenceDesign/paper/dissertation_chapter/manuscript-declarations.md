# Manuscript Abstract and Declarations — Preserved, Not in the Chapter

Preserved per `docs/plans/manuscript-unification-plan.md` step 0 (component-fate
table: no standalone abstract in the chapter body; the dissertation front matter
carries its own). Kept verbatim, with LaTeX escaping cleaned up only (`\noindent`
dropped, `\%` → `%`; no wording changed), so re-inserting it later is a paste,
not a rewrite. Source: `Dissertation_Chapter.qmd` as of 2026-09-24 (identical to
the 2026-09-23 snapshot in
`archive/IncidenceDesign_Manuscripts_PreUnification_2026-09-23/`, lines 26–28).

**Numbers superseded:** the figures below (MSE 0.079/0.802, coverage 94%/53%,
Friedman χ² 557.3–800.9, 9,600 scenarios) come from the April 2026 simulation and
are replaced by the 2026-09 revision (`docs/plans/simulation-revision-plan.md`).

CTJ's Declarations (conflicting interests, funding) are added here at step 7.

## Abstract

Cluster randomized trials (CRTs) are the standard design when an intervention must be delivered at the group level, but they are vulnerable to two closely related threats when clusters are spatially arranged: spatial spillover, in which a treatment effect diffuses from treated to nearby control clusters, and spatial dependence, in which a cluster's outcome is correlated with its neighbors' outcomes independent of treatment. The existing methodological literature on spatial interference is concentrated almost entirely on analysis-stage correction, detecting and modeling spillover after outcomes are observed, and has largely neglected design-stage mitigation: choosing a treatment assignment strategy, before randomization, that limits spillover's ability to bias estimation in the first place. This chapter extends an earlier proof-of-concept study of a single block-stratified design on small, exhaustively enumerable grids to a substantially larger and more realistic setting: a 100-cluster spatial grid, six conceptually distinct treatment assignment strategies, and three mechanisms of heterogeneous baseline outcome incidence, so that some designs can exploit prior knowledge of spatial risk while others cannot. Outcomes were generated from a spatial Durbin model and estimated with maximum likelihood spatial regression across 9,600 simulated scenarios spanning five magnitudes of the true treatment effect. An incidence-guided saturation design, which varies treatment intensity across spatial regions in proportion to historical incidence, achieved the lowest mean squared error (MSE = 0.079) with near-nominal 95% confidence interval coverage (94%), while a design that maximizes spatial separation between treated and control clusters (checkerboard-style alternation) performed worst by a wide margin (MSE = 0.802) with coverage collapsing to 53%. This ranking was stable across all five tested effect magnitudes (Friedman $\chi^2$ range 557.3--800.9, all $p < 2.2\times10^{-16}$) and across all three incidence-heterogeneity mechanisms, though the second-best design's advantage depended heavily on how informative the observed incidence surface was. A preliminary, explicitly placeholder application to a planned sudden unexpected death prevention pilot in North Carolina, delivered through the state's community college network, illustrated the same qualitative pattern on a real, irregular geography. These findings support incorporating spatial spillover structure into the design stage of cluster randomized trials as an actionable complement to, rather than a substitute for, analysis-stage correction.
